"""Coverage as a facet (RM278): depth tools' outputs, polars-bio depth over a BAM, callable positions."""

import gzip
import re
from datetime import timedelta
from pathlib import Path

import polars as pl
import pytest

from just_vcf import (
    Coverage,
    CoverageFacet,
    CoverageFormat,
    Resolution,
    UnsupportedCoverageError,
    alignment_coverage,
    callable_at,
    callable_sites,
    clear_cache,
    read_coverage,
)
from just_vcf import cache
from just_vcf import coverage as coverage_module
from just_vcf.assembly import ASSEMBLY_LENGTHS
from just_vcf.coverage import CALLABLE_REASONS, COVERAGE_INTERVALS

# --- An independent oracle: per-base depth from read definitions ----------------------------------

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def depth_of(reads: list[tuple[str, int, str]], *, contig: str) -> dict[int, int]:
    """1-based position → depth on `contig`, from `(contig, 0-based start, cigar)` reads. Only M, =
    and X put a read base on the reference; D and N move along it without one (mosdepth's rule)."""
    depth: dict[int, int] = {}
    for name, start, cigar in reads:
        if name != contig:
            continue
        at = start
        for length, op in ((int(n), o) for n, o in _CIGAR.findall(cigar)):
            if op in "M=X":
                for base in range(at + 1, at + length + 1):
                    depth[base] = depth.get(base, 0) + 1
            if op in "M=XDN":
                at += length
    return depth


def windows_of(depth: dict[int, int], length: int, window: int) -> list[tuple[int, int, float, int, int]]:
    """`(pos, end_pos, mean, min, max)` per window over a contig of `length`, zeros filled in."""
    out = []
    for start in range(1, length + 1, window):
        end = min(start + window - 1, length)
        values = [depth.get(base, 0) for base in range(start, end + 1)]
        out.append((start, end, sum(values) / len(values), min(values), max(values)))
    return out


# --- A tiny BAM written with pysam --------------------------------------------------------------

CONTIGS = [("chr1", 100), ("chrX", 50), ("chrY", 40)]
READS = [("chr1", 10, "20M"), ("chr1", 15, "5M3D10M"), ("chr1", 15, "20M"), ("chr1", 60, "4M2N4M"), ("chrX", 0, "10M")]


def write_bam(path: Path, reads: list[tuple[str, int, str]], *, flags: list[int] | None = None, mapq: list[int] | None = None) -> Path:
    pysam = pytest.importorskip("pysam")
    names = [name for name, _ in CONTIGS]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS]}
    order = sorted(range(len(reads)), key=lambda i: (names.index(reads[i][0]), reads[i][1]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for i in order:
            name, start, cigar = reads[i]
            record = pysam.AlignedSegment()
            record.query_name = f"r{i}"
            record.flag = flags[i] if flags else 0
            record.reference_id = names.index(name)
            record.reference_start = start
            record.mapping_quality = mapq[i] if mapq else 60
            record.cigarstring = cigar
            n = record.infer_query_length()
            record.query_sequence = "A" * n
            record.query_qualities = pysam.qualitystring_to_array("I" * n)
            handle.write(record)
    pysam.index(str(path))
    return path


def rows_of(coverage: Coverage, chrom: str) -> list[tuple[int, int, float, int, int]]:
    frame = coverage.intervals.filter(pl.col("chrom") == chrom).collect()
    return [(r["pos"], r["end_pos"], pytest.approx(r["depth_mean"], abs=1e-6), r["depth_min"], r["depth_max"]) for r in frame.iter_rows(named=True)]


def test_alignment_windows_equal_the_reads_depth(tmp_path: Path) -> None:
    coverage = alignment_coverage(write_bam(tmp_path / "a.bam", READS), window=10)
    for source, length in CONTIGS:
        chrom = source.removeprefix("chr")
        assert rows_of(coverage, chrom) == windows_of(depth_of(READS, contig=source), length, 10)
    frame = coverage.intervals.collect()
    assert frame.schema == pl.Schema(COVERAGE_INTERVALS)
    # Every window of every header contig is written and described at every base, zeros included.
    assert frame.height == sum(-(-length // 10) for _, length in CONTIGS)
    assert (frame["bases"] == frame["end_pos"] - frame["pos"] + 1).all()
    facet = coverage.facet
    assert (facet.format_id, facet.resolution, facet.window, facet.zero_filled) == (CoverageFormat.ALIGNMENT_DEPTH, Resolution.WINDOW, 10, True)
    assert facet.options == {"filter_flag": 1796, "min_mapping_quality": 0}
    assert [c["chrom"] for c in facet.contig_lengths] == ["1", "X", "Y"]
    assert facet.source_sha256 is None and facet.assembly_label is None


def test_a_read_on_a_contigs_last_base_is_counted(tmp_path: Path) -> None:
    # polars-bio 0.36's block mode emits nothing for a contig whose reads reach its last base (on a
    # real BAM: the last window of 60 contigs read a zero there); per-base mode counts it.
    reads = [("chrX", 30, "10M"), ("chrX", 40, "10M"), ("chrY", 0, "40M")]
    coverage = alignment_coverage(write_bam(tmp_path / "end.bam", reads), window=10)
    for source, length in CONTIGS:
        assert rows_of(coverage, source.removeprefix("chr")) == windows_of(depth_of(reads, contig=source), length, 10)


def test_alignment_depth_is_deterministic_and_its_id_follows_the_options(tmp_path: Path) -> None:
    path = write_bam(tmp_path / "a.bam", READS)
    first, second = alignment_coverage(path, window=10), alignment_coverage(path, window=10)
    assert first.facet.facet_id == second.facet.facet_id
    assert first.intervals.collect().equals(second.intervals.collect())
    assert alignment_coverage(path, window=20).facet.facet_id != first.facet.facet_id
    assert alignment_coverage(path, window=10, min_mapping_quality=20).facet.facet_id != first.facet.facet_id


def test_the_read_filters_are_the_substrates(tmp_path: Path) -> None:
    # A duplicate (0x400) and a MAPQ-10 read: the default flag mask drops the first, a MAPQ floor the second.
    path = write_bam(tmp_path / "f.bam", READS, flags=[0, 0, 0x400, 0, 0], mapq=[60, 10, 60, 60, 60])
    kept = [read for i, read in enumerate(READS) if i != 2]
    assert rows_of(alignment_coverage(path, window=10), "1") == windows_of(depth_of(kept, contig="chr1"), 100, 10)
    floored = [read for i, read in enumerate(kept) if i != 1]
    assert rows_of(alignment_coverage(path, window=10, min_mapping_quality=20), "1") == windows_of(depth_of(floored, contig="chr1"), 100, 10)


def test_a_depth_past_int16_is_left_out_and_counted(tmp_path: Path) -> None:
    # polars-bio 0.36 writes depth as Int16 by a wrapping cast: 32,770 reads over one base come back negative.
    deep = [("chr1", 40, "5M")] * 32_770
    coverage = alignment_coverage(write_bam(tmp_path / "deep.bam", [*READS, *deep]), window=10)
    window = coverage.intervals.filter((pl.col("chrom") == "1") & (pl.col("pos") == 41)).collect().row(0, named=True)
    assert window["bases"] == 5  # the five wrapped bases are not described
    assert (window["depth_min"], window["depth_max"]) == (0, 0)
    (warning,) = coverage.facet.warnings
    assert warning.startswith("5 base(s) in 1 window(s) came back with a negative depth")


def test_an_alignment_is_told_by_its_magic_bytes(tmp_path: Path) -> None:
    fake = tmp_path / "fake.bam"
    fake.write_bytes(b"not a bam at all")
    with pytest.raises(UnsupportedCoverageError, match="magic bytes"):
        alignment_coverage(fake)


def test_a_refusal_names_the_kind_of_input_never_its_file(tmp_path: Path) -> None:
    """A path is a label and can carry a person's name: the message says what the input is instead."""
    planted = tmp_path / "PlantedPerson"
    planted.mkdir()
    (planted / "PlantedPerson.bam").write_bytes(b"not a bam at all")
    (planted / "PlantedPerson.bed").write_bytes(b"")
    headerless = write_bam(planted / "PlantedPerson-noref.bam", [])
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(headerless), "wb", header={"HD": {"VN": "1.6"}}):
        pass
    cases = [
        (lambda: alignment_coverage(planted / "PlantedPerson.bam"), "the alignment input is not a BAM or CRAM by its magic bytes"),
        (lambda: read_coverage(planted / "PlantedPerson.bed"), "the coverage input is an empty file"),
        (lambda: alignment_coverage(headerless), "the BAM header lists no reference sequence, so no window can be placed"),
    ]
    for call, message in cases:
        with pytest.raises(UnsupportedCoverageError) as caught:
            call()
        assert str(caught.value) == message and "Planted" not in str(caught.value)


# --- Depth tools' outputs -----------------------------------------------------------------------


def write(path: Path, lines: list[str], *, compress: bool = False) -> Path:
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(gzip.compress(data) if compress else data)
    return path


def per_base_lines(reads: list[tuple[str, int, str]]) -> list[str]:
    """mosdepth `.per-base.bed` runs (0-based, half-open) over CONTIGS, zero runs included."""
    out = []
    for source, length in CONTIGS:
        depth = depth_of(reads, contig=source)
        values = [depth.get(base, 0) for base in range(1, length + 1)]
        start = 0
        for i in range(1, length + 1):
            if i == length or values[i] != values[start]:
                out.append(f"{source}\t{start}\t{i}\t{values[start]}")
                start = i
    return out


@pytest.mark.parametrize("compress", [False, True])
def test_per_base_runs_fold_into_the_same_windows_as_the_alignment(tmp_path: Path, compress: bool) -> None:
    # The BED is named like a regions file: the kind comes from the content, compression from the magic bytes.
    path = write(tmp_path / "x.regions.bed", per_base_lines(READS), compress=compress)
    coverage = read_coverage(path, window=10)
    assert coverage.facet.format_id is CoverageFormat.PER_BASE_BED
    assert (coverage.facet.resolution, coverage.facet.window, coverage.facet.zero_filled) == (Resolution.WINDOW, 10, False)
    for source, length in CONTIGS:
        assert rows_of(coverage, source.removeprefix("chr")) == windows_of(depth_of(READS, contig=source), length, 10)


def test_folding_across_chunks_gives_the_same_windows(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = write(tmp_path / "p.bed", per_base_lines(READS))
    whole = read_coverage(path, window=7).intervals.collect()
    monkeypatch.setattr(coverage_module, "_CHUNK_BYTES", 23)  # cuts lines mid-way, many times
    assert read_coverage(path, window=7).intervals.collect().equals(whole)


def test_without_a_window_each_run_is_an_interval(tmp_path: Path) -> None:
    lines = per_base_lines(READS)
    coverage = read_coverage(write(tmp_path / "p.bed", lines), window=None)
    frame = coverage.intervals.collect()
    assert coverage.facet.resolution is Resolution.BASE and coverage.facet.window is None
    assert frame.height == len(lines)
    first = lines[0].split("\t")
    assert frame.row(0) == ("1", int(first[1]) + 1, int(first[2]), float(first[3]), int(first[3]), int(first[3]), int(first[2]) - int(first[1]), None)


def test_samtools_depth_without_zero_rows_leaves_windows_partly_described(tmp_path: Path) -> None:
    depth = depth_of(READS, contig="chr1")
    lines = ["#CHROM\tPOS\ta.bam", *(f"chr1\t{pos}\t{value}" for pos, value in sorted(depth.items()))]
    coverage = read_coverage(write(tmp_path / "d.txt", lines), window=10)
    assert coverage.facet.format_id is CoverageFormat.SAMTOOLS_DEPTH
    frame = coverage.intervals.collect()
    # Only the windows a row reached, each describing just the bases samtools wrote.
    expected = {w: [p for p in depth if (p - 1) // 10 == w] for w in {(p - 1) // 10 for p in depth}}
    assert sorted(zip(frame["pos"], frame["bases"])) == sorted((w * 10 + 1, len(ps)) for w, ps in expected.items())
    row = frame.filter(pl.col("pos") == 31).row(0, named=True)
    assert row["depth_mean"] == pytest.approx(sum(depth[p] for p in expected[3]) / len(expected[3]))


def test_a_headerless_three_column_file_of_intervals_needs_its_kind_stated(tmp_path: Path) -> None:
    # A BED3 and a samtools depth near a contig's start (chrM at thousands deep) look alike.
    path = write(tmp_path / "x.txt", ["chrM\t1\t3000", "chrM\t2\t3010", "chrM\t3\t3021"])
    with pytest.raises(UnsupportedCoverageError, match="kind=CoverageFormat.SAMTOOLS_DEPTH"):
        read_coverage(path)
    stated = read_coverage(path, window=None, kind=CoverageFormat.SAMTOOLS_DEPTH)
    assert stated.intervals.select("chrom", "pos", "depth_min").collect().rows() == [("MT", 1, 3000), ("MT", 2, 3010), ("MT", 3, 3021)]
    # Further in, depth stays under the position: no ambiguity.
    assert read_coverage(write(tmp_path / "y.txt", ["chr1\t10001\t3", "chr1\t10002\t4"])).facet.format_id is CoverageFormat.SAMTOOLS_DEPTH


def test_samtools_depth_of_several_files_is_refused(tmp_path: Path) -> None:
    path = write(tmp_path / "d.txt", ["#CHROM\tPOS\ta.bam\tb.bam", "chr1\t11\t1\t2"])
    with pytest.raises(UnsupportedCoverageError, match="2 depth columns"):
        read_coverage(path)


@pytest.mark.parametrize("named", [False, True])
def test_regions_keep_their_mean_and_bound_nothing(tmp_path: Path, named: bool) -> None:
    lines = ["chr1\t0\t500\tgeneA\t12.50", "chrX\t500\t1000\tgeneB\t3.25"] if named else ["chr1\t0\t500\t12.50", "chrX\t500\t1000\t3.25"]
    coverage = read_coverage(write(tmp_path / "r.bed.gz", lines, compress=True))
    assert (coverage.facet.format_id, coverage.facet.resolution) == (CoverageFormat.REGIONS_BED, Resolution.REGION)
    assert coverage.intervals.collect().rows() == [
        ("1", 1, 500, 12.5, None, None, 500, "geneA" if named else None),
        ("X", 501, 1000, 3.25, None, None, 500, "geneB" if named else None),
    ]


def test_quantized_bins_become_depth_bounds(tmp_path: Path) -> None:
    lines = ["chr1\t0\t10\t0:1", "chr1\t10\t30\t1:4", "chr1\t30\t60\t4:inf", "chr1\t60\t70\tCALLABLE"]
    coverage = read_coverage(write(tmp_path / "q.bed", lines))
    assert (coverage.facet.format_id, coverage.facet.resolution) == (CoverageFormat.QUANTIZED_BED, Resolution.BIN)
    assert coverage.intervals.select("pos", "end_pos", "depth_min", "depth_max", "depth_mean").collect().rows() == [
        (1, 10, 0, 0, None), (11, 30, 1, 3, None), (31, 60, 4, None, None), (61, 70, None, None, None),
    ]
    assert coverage.facet.warnings == ["1 quantized interval(s) carry a bin label that names no depth range (a renamed bin); their depth is unknown"]


def test_mosdepth_summary_gives_contig_bounds_and_the_build(tmp_path: Path) -> None:
    grch38 = {chrom: next(iter(lengths)) if isinstance(lengths, (set, frozenset)) else lengths for chrom, lengths in ASSEMBLY_LENGTHS["GRCh38"].items()}
    names = ["1", "2", "3", "20", "X"]
    lines = ["chrom\tlength\tbases\tmean\tmin\tmax", *(f"chr{c}\t{grch38[c]}\t{grch38[c] * 30}\t30.00\t0\t900" for c in names)]
    lines += ["chr1_region\t1000\t30000\t30.00\t0\t50", f"total\t{sum(grch38[c] for c in names)}\t0\t30.00\t0\t900"]
    coverage = read_coverage(write(tmp_path / "s.mosdepth.summary.txt", lines))
    facet = coverage.facet
    assert (facet.format_id, facet.assembly_label, facet.assembly_label_basis.value) == (CoverageFormat.MOSDEPTH_SUMMARY, "GRCh38", "inferred")
    assert coverage.intervals.select("chrom", "pos", "end_pos", "depth_min", "depth_max").collect().rows() == [
        (c, 1, grch38[c], 0, 900) for c in names
    ]
    assert facet.warnings == ["2 summary row(s) for the total or for --by regions carry no interval and were left out"]
    # A stated build that the lengths contradict is kept, with a warning.
    stated = read_coverage(tmp_path / "s.mosdepth.summary.txt", assembly_label="GRCh37").facet
    assert (stated.assembly_label, stated.assembly_label_basis.value) == ("GRCh37", "stated")
    assert "disagrees with the contig lengths" in stated.warnings[-1]


def test_samtools_coverage_keeps_the_mean_and_names_what_it_drops(tmp_path: Path) -> None:
    lines = [
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq",
        "chrX\t1\t156040895\t1000\t150000000\t96.1\t15.2\t35\t58",
    ]
    coverage = read_coverage(write(tmp_path / "c.txt", lines))
    assert coverage.facet.format_id is CoverageFormat.SAMTOOLS_COVERAGE
    assert coverage.intervals.collect().rows() == [("X", 1, 156_040_895, pytest.approx(15.2), None, None, 156_040_895, None)]
    assert coverage.facet.dropped_fields == ["numreads", "covbases", "coverage", "meanbaseq", "meanmapq"]


@pytest.mark.parametrize(
    ("lines", "match"),
    [
        (["#chrom\tstart\tend\tregion\t1X\t10X", "chr1\t0\t10\tunknown\t10\t5"], "thresholds"),
        (["chr1\tfoo\t10\t3"], "not a BED of depths"),
        ([""], "empty"),
    ],
)
def test_what_is_not_a_coverage_output_is_refused(tmp_path: Path, lines: list[str], match: str) -> None:
    with pytest.raises(UnsupportedCoverageError, match=match):
        read_coverage(write(tmp_path / "bad.txt", lines))


# --- Callable --------------------------------------------------------------------------------------


def synthetic(rows: list[tuple]) -> Coverage:
    intervals = pl.DataFrame(rows, schema=list(COVERAGE_INTERVALS), orient="row").cast(COVERAGE_INTERVALS).lazy()
    facet = CoverageFacet(
        facet_id="c", format_id=CoverageFormat.REGIONS_BED, source_path="synthetic", source_sha256=None, assembly_label=None,
        assembly_label_basis=None, contig_lengths=None, resolution=Resolution.WINDOW, window=10, options={}, zero_filled=False,
        dropped_fields=[], warnings=[],
    )
    return Coverage(facet=facet, intervals=intervals)


COV = synthetic([
    ("1", 1, 10, 12.0, 10, 14, 10, None),     # every base at 10 or more
    ("1", 11, 20, 1.5, 0, 3, 10, None),       # every base under 4
    ("1", 21, 30, 6.0, 2, 9, 10, None),       # straddles 4
    ("1", 31, 40, 20.0, None, None, 10, None),  # a mean only
    ("1", 41, 50, 15.0, 10, 20, 6, None),     # only 6 bases described
    ("1", 61, 70, None, 4, None, 10, "4:inf"),  # an open bin
])


@pytest.mark.parametrize(
    ("pos", "expected", "reason"),
    [(5, True, None), (15, False, None), (25, None, "straddles"), (35, None, "no_bounds"), (45, None, "partly_described"),
     (55, None, "no_interval"), (65, True, None), (100, None, "no_interval")],
)
def test_callable_at_reads_the_bounds_never_the_mean(pos: int, expected: bool | None, reason: str | None) -> None:
    found = callable_at(COV, "chr1", pos, min_depth=4)
    assert (found.callable, found.reason) == (expected, reason)
    assert bool(found) is (expected is True)
    assert reason is None or reason in CALLABLE_REASONS


def test_an_open_bin_under_the_threshold_straddles() -> None:
    assert callable_at(COV, "1", 65, min_depth=10).reason == "straddles"


def test_callable_sites_keeps_the_sites_order_and_columns() -> None:
    sites = pl.DataFrame({"chrom": ["1", "1", "2", "1"], "pos": [25, 5, 5, 15], "tag": ["a", "b", "c", "d"]})
    out = callable_sites(COV, sites, min_depth=4).collect()
    assert out["tag"].to_list() == ["a", "b", "c", "d"]
    assert out["callable"].to_list() == [None, True, None, False]
    assert out["callable_reason"].to_list() == ["straddles", None, "no_interval", None]
    assert out["coverage_depth_min"].to_list() == [2, 10, None, 0]
    assert out.columns == ["chrom", "pos", "tag", "coverage_depth_min", "coverage_depth_max", "coverage_depth_mean", "callable", "callable_reason"]


# --- The cache spill: a BAM is read once (RM278) -----------------------------------------------------


@pytest.fixture
def depth_calls(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Each path polars-bio's `depth()` is run on, recorded; the real `depth()` still runs."""
    seen: list[str] = []
    real = coverage_module.pb.depth

    def counted(path, *args, **kwargs):
        seen.append(path)
        return real(path, *args, **kwargs)

    monkeypatch.setattr(coverage_module.pb, "depth", counted)
    return seen


def coverage_manifests() -> dict[str, cache.CacheEntry]:
    base = cache.cache_dir() / "coverage"
    return {path.stem: cache.CacheEntry.model_validate_json(path.read_bytes()) for path in sorted(base.glob("*.json"))} if base.is_dir() else {}


def test_a_second_call_reads_the_spill_not_the_bam(tmp_path: Path, depth_calls: list[str]) -> None:
    path = write_bam(tmp_path / "a.bam", [*READS, *[("chr1", 40, "5M")] * 32_770])
    first = alignment_coverage(path, window=10)
    second = alignment_coverage(path, window=10)
    assert depth_calls == [str(path)]
    assert second.facet == first.facet  # the Int16 warning included, given back from the spill
    assert second.facet.warnings and second.facet.warnings[0].startswith("5 base(s) in 1 window(s)")
    assert second.intervals.collect().equals(first.intervals.collect())
    for source, length in CONTIGS[1:]:
        assert rows_of(second, source.removeprefix("chr")) == windows_of(depth_of(READS, contig=source), length, 10)
    ((name, manifest),) = coverage_manifests().items()
    assert manifest.kind == cache.COVERAGE_SPILL
    assert manifest.versions == {
        "COVERAGE_VERSION": str(coverage_module.COVERAGE_VERSION),
        "POLARS_BIO_VERSION": coverage_module.POLARS_BIO_VERSION,
        "SPILL_LAYOUT": str(cache.SPILL_LAYOUT),
    }
    assert (manifest.facet_id, manifest.source_sha256, manifest.source_path) == (first.facet.facet_id, None, str(path.resolve()))
    assert manifest.size_bytes == sum(f.stat().st_size for f in (cache.cache_dir() / "coverage" / name).rglob("*") if f.is_file())


def test_a_spill_is_never_served_for_other_parameters_or_versions(tmp_path: Path, depth_calls: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_bam(tmp_path / "a.bam", READS, mapq=[60, 10, 60, 60, 60])
    base = alignment_coverage(path, window=10).intervals.collect()
    runs = [
        (lambda: alignment_coverage(path, window=20), 20, 0),
        (lambda: alignment_coverage(path, window=10, min_mapping_quality=20), 10, 20),
        (lambda: alignment_coverage(path, window=10, filter_flag=0), 10, 0),
    ]
    for call, window, mapq in runs:
        before = len(depth_calls)
        found = call()
        assert len(depth_calls) == before + 1
        kept = [read for i, read in enumerate(READS) if not (mapq and i == 1)]
        assert rows_of(found, "1") == windows_of(depth_of(kept, contig="chr1"), 100, window)
    # The same parameters under another COVERAGE_VERSION or another polars-bio: computed again.
    for name, value in (("COVERAGE_VERSION", coverage_module.COVERAGE_VERSION + 1), ("POLARS_BIO_VERSION", "0.0.1")):
        with monkeypatch.context() as patch:
            patch.setattr(coverage_module, name, value)
            before = len(depth_calls)
            assert alignment_coverage(path, window=10).intervals.collect().equals(base)
            assert len(depth_calls) == before + 1
    assert len(coverage_manifests()) == 6
    # A stated build changes no interval: the spill is reused, under the facet id the build gives.
    before = len(depth_calls)
    stated = alignment_coverage(path, window=10, assembly_label="GRCh38")
    assert len(depth_calls) == before
    assert stated.intervals.collect().equals(base)


def test_the_spill_is_swept_and_cleared_like_the_loaders(tmp_path: Path, depth_calls: list[str]) -> None:
    path = write_bam(tmp_path / "a.bam", READS)
    other = write_bam(tmp_path / "b.bam", READS[:2])
    alignment_coverage(other, window=10)
    old = {"COVERAGE_VERSION": coverage_module.COVERAGE_VERSION - 1, "POLARS_BIO_VERSION": coverage_module.POLARS_BIO_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    newer = {**old, "COVERAGE_VERSION": coverage_module.COVERAGE_VERSION + 1}
    root = cache.cache_dir() / "coverage"
    entries = {}
    for name, versions, days in (("a", old, cache.GRACE_DAYS + 1), ("b", old, 0), ("c", newer, cache.GRACE_DAYS + 1)):
        entry = root / (name * 32)
        with cache.writing(entry, cache.COVERAGE_SPILL, versions, cache.Origin(source_sha256=None, source_path=None)):
            entry.mkdir(parents=True)
            (entry / "intervals.parquet").write_bytes(b"x" * 100)
        manifest_path = root / f"{entry.name}.json"
        then = cache._now() - timedelta(days=days)
        manifest = cache.CacheEntry.model_validate_json(manifest_path.read_bytes())
        cache._write_manifest(manifest_path, manifest.model_copy(update={"created": then, "last_used": then}))
        entries[name] = entry
    alignment_coverage(path, window=10)
    # The older version's entry unused past the grace period goes; a recent one and a newer install's stay.
    assert not entries["a"].exists() and entries["b"].exists() and entries["c"].exists()
    # clear_cache by the BAM removes its spill and leaves the other BAM's.
    report = clear_cache(source=path)
    assert [info.kind for info in report.removed] == [cache.COVERAGE_SPILL]
    assert sorted(m.source_path or "" for m in coverage_manifests().values()) == ["", "", str(other.resolve())]
    alignment_coverage(path, window=10)
    assert depth_calls.count(str(path)) == 2
