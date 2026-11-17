"""Coverage as a facet (RM278): depth tools' outputs, polars-bio depth over a BAM, callable positions."""

import gzip
import hashlib
import re
from datetime import timedelta
from itertools import pairwise
from pathlib import Path

import polars as pl
import pytest
from polars_bio import polars_bio as pb_native

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
from just_vcf.coverage import CALLABLE_REASONS, COVERAGE_INTERVALS, COVERAGE_THRESHOLDS

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


def write_bam(
    path: Path, reads: list[tuple[str, int, str]], *, flags: list[int] | None = None, mapq: list[int] | None = None, index: bool = True
) -> Path:
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
    if index:
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


def test_a_quantized_label_written_empty_is_unbounded_and_a_missing_one_is_short(tmp_path: Path) -> None:
    """A label written empty names no depth range, as a renamed bin does; a row short of the label field
    is counted with the short rows instead, never twice (RM291)."""
    lines = ["chr1\t0\t10\t0:1", "chr1\t10\t20\t", "chr1\t20\t30\tCALLABLE", "chr1\t30\t40", "chr1\t40\t50\t1:4"]
    warnings = read_coverage(write(tmp_path / "q.bed", lines), kind=CoverageFormat.QUANTIZED_BED).facet.warnings
    fields = [line.split("\t") for line in lines]
    unbounded = sum(len(f) == 4 and not re.fullmatch(r"\d+:(\d+|inf)", f[3]) for f in fields)
    assert warnings[0].startswith(f"{unbounded} quantized interval(s) carry a bin label that names no depth range")
    assert counted(warnings) == {"short": sum(len(f) < 4 for f in fields)}
    assert len(warnings) == 2


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


THRESHOLDS = [0, 1, 2, 3]
REGIONS = [
    # `--by 10` on chr1 (mosdepth names such regions `unknown`), a named BED region on X, and Y, where no read lands.
    *[("chr1", start, min(start + 10, 100), "unknown") for start in range(0, 100, 10)],
    ("chrX", 5, 30, "geneX"),
    ("chrY", 0, 40, "unknown"),
]


def thresholds_lines(reads: list[tuple[str, int, str]], regions: list[tuple[str, int, int, str]], thresholds: list[int]) -> list[str]:
    """mosdepth's `.thresholds.bed` (`write_header`, `write_thresholds`): per region, the bases at a depth of each threshold or more."""
    out = ["#chrom\tstart\tend\tregion\t" + "\t".join(f"{t}X" for t in thresholds)]
    for source, start, end, name in regions:
        depth = depth_of(reads, contig=source)
        counts = [sum(1 for base in range(start + 1, end + 1) if depth.get(base, 0) >= t) for t in thresholds]
        out.append("\t".join([source, str(start), str(end), name, *map(str, counts)]))
    return out


def test_a_thresholds_file_keeps_its_counts_apart_from_every_depth(tmp_path: Path) -> None:
    coverage = read_coverage(write(tmp_path / "t.thresholds.bed.gz", thresholds_lines(READS, REGIONS, THRESHOLDS), compress=True))
    facet = coverage.facet
    assert (facet.format_id, facet.resolution, facet.thresholds) == (CoverageFormat.THRESHOLDS_BED, Resolution.THRESHOLDS, THRESHOLDS)
    assert (facet.window, facet.zero_filled, facet.contig_lengths, facet.dropped_fields, facet.options) == (None, False, None, [], {})
    intervals = coverage.intervals.collect()
    assert intervals.schema == pl.Schema(COVERAGE_INTERVALS)
    # The regions, and not one depth: a count of bases at a threshold is no depth of any base.
    assert intervals.select("chrom", "pos", "end_pos", "bases", "label").rows() == [
        (source.removeprefix("chr"), start + 1, end, end - start, name) for source, start, end, name in REGIONS
    ]
    assert intervals.select(pl.col("depth_mean", "depth_min", "depth_max").null_count()).row(0) == (len(REGIONS),) * 3
    counts = coverage.thresholds.collect()
    assert counts.schema == pl.Schema(COVERAGE_THRESHOLDS)
    expected = [
        (source.removeprefix("chr"), start + 1, end, name, t, sum(1 for b in range(start + 1, end + 1) if depth_of(READS, contig=source).get(b, 0) >= t))
        for source, start, end, name in REGIONS
        for t in THRESHOLDS
    ]
    assert counts.rows() == expected
    # mosdepth's own invariants hold: at 0X every base, no count above its region, none rising with the threshold.
    by_region = counts.group_by("chrom", "pos", maintain_order=True).agg("bases_at_or_above", length=(pl.col("end_pos") - pl.col("pos") + 1).first())
    for row in by_region.iter_rows(named=True):
        assert row["bases_at_or_above"][0] == row["length"]
        assert row["bases_at_or_above"] == sorted(row["bases_at_or_above"], reverse=True)
    # What is callable cannot be read from it: the depth columns bound nothing.
    assert callable_at(coverage, "1", 25, min_depth=1).reason == "no_bounds"


def test_a_thresholds_file_with_no_reads_reads_as_written(tmp_path: Path) -> None:
    # A BAM with no read at all (mosdepth's `tid == -2`) writes 0 at every threshold, 0X included.
    lines = ["#chrom\tstart\tend\tregion\t0X\t10X", "chr1\t0\t500\tunknown\t0\t0"]
    counts = read_coverage(write(tmp_path / "t.bed", lines)).thresholds.collect()
    assert counts.select("threshold", "bases_at_or_above").rows() == [(0, 0), (10, 0)]


@pytest.mark.parametrize(
    ("lines", "match"),
    [
        (["#chrom\tstart\tend\tregion\t1X\tdeep", "chr1\t0\t10\tunknown\t10\t5"], "not all `<n>X` thresholds"),
        (["#chrom\tstart\tend\tregion\t10X\t1X", "chr1\t0\t10\tunknown\t5\t10"], "do not rise strictly"),
        (["#chrom\tstart\tend\tregion\t1X\t10X", "chr1\t0\t10\tunknown\t11\t5"], "1 row(s) counting more bases than their region"),
        (["#chrom\tstart\tend\tregion\t1X\t10X", "chr1\t0\t10\tunknown\t5\t6"], "1 counting more bases at a higher threshold"),
    ],
)
def test_a_thresholds_file_mosdepth_could_not_have_written_is_refused(tmp_path: Path, lines: list[str], match: str) -> None:
    with pytest.raises(UnsupportedCoverageError, match=re.escape(match)):
        read_coverage(write(tmp_path / "t.bed", lines))


def test_a_thresholds_file_without_its_header_is_never_read_as_depths(tmp_path: Path) -> None:
    """tabix drops a file's header unless run with -h. Without it one threshold's counts look like a
    regions file's means, and several look like nothing: both are refused rather than read."""
    one = thresholds_lines(READS, REGIONS, [1])[1:]
    with pytest.raises(UnsupportedCoverageError, match="kind=CoverageFormat.REGIONS_BED"):
        read_coverage(write(tmp_path / "one.bed", one))
    with pytest.raises(UnsupportedCoverageError, match="without its header line, which names the thresholds"):
        read_coverage(write(tmp_path / "many.bed", thresholds_lines(READS, REGIONS, THRESHOLDS)[1:]))
    with pytest.raises(UnsupportedCoverageError, match="must start with mosdepth's header line"):
        read_coverage(write(tmp_path / "stated.bed", one), kind=CoverageFormat.THRESHOLDS_BED)
    # A regions file of whole-number means, stated, still reads.
    stated = read_coverage(write(tmp_path / "r.bed", ["chr1\t0\t500\tgeneA\t12"]), kind=CoverageFormat.REGIONS_BED)
    assert stated.intervals.select("depth_mean", "label").collect().rows() == [(12.0, "geneA")]


@pytest.mark.parametrize(
    ("lines", "match"),
    [
        (["chr1\tfoo\t10\t3"], "not a BED of depths"),
        ([""], "empty"),
    ],
)
def test_what_is_not_a_coverage_output_is_refused(tmp_path: Path, lines: list[str], match: str) -> None:
    with pytest.raises(UnsupportedCoverageError, match=match):
        read_coverage(write(tmp_path / "bad.txt", lines))


def widened(lines: list[str], rows: list[int]) -> list[str]:
    """`lines` with two trailing fields added to the lines at `rows`."""
    return [f"{line}\textra\t9" if i in rows else line for i, line in enumerate(lines)]


def fixed_width_files() -> dict[CoverageFormat, list[str]]:
    """One file of each format read with a fixed set of columns."""
    return {
        CoverageFormat.PER_BASE_BED: per_base_lines(READS),
        CoverageFormat.SAMTOOLS_DEPTH: [f"chr1\t{pos}\t{value}" for pos, value in sorted(depth_of(READS, contig="chr1").items())],
        CoverageFormat.QUANTIZED_BED: ["chr1\t0\t10\t0:1", "chr1\t10\t30\t1:4", "chr1\t30\t60\t4:inf", "chr1\t60\t70\t1:4"],
        CoverageFormat.THRESHOLDS_BED: thresholds_lines(READS, REGIONS, THRESHOLDS),
    }


@pytest.mark.parametrize("window", [None, 7])
@pytest.mark.parametrize("wide", ["first", "later"])
@pytest.mark.parametrize("kind", list(fixed_width_files()))
def test_rows_wider_than_the_format_read_as_their_leading_fields(tmp_path: Path, kind: CoverageFormat, wide: str, window: int | None) -> None:
    """A row with fields past the format's reads as the same row without them, whether the first data row
    is the wide one or later ones are (RM289): polars 2 refuses a headerless read whose first row is wider
    than the names, where 1.x read the extra fields into columns of their own."""
    lines = fixed_width_files()[kind]
    data = [i for i, line in enumerate(lines) if not line.startswith("#")]
    mixed = widened(lines, data[:1] if wide == "first" else data[1::2])
    plain = read_coverage(write(tmp_path / "plain.txt", lines), window=window, kind=kind)
    read = read_coverage(write(tmp_path / "mixed.txt", mixed), window=window, kind=kind)
    assert read.intervals.collect().equals(plain.intervals.collect())
    assert (read.thresholds is None) == (plain.thresholds is None)
    if plain.thresholds is not None:
        assert read.thresholds.collect().equals(plain.thresholds.collect())


@pytest.mark.parametrize("kind", list(fixed_width_files()))
def test_a_first_row_short_of_a_field_reads_it_as_empty(tmp_path: Path, kind: CoverageFormat) -> None:
    """A row lacking the format's last field reads as the row with that field written empty, the first
    data row included, which polars 2 checks the columns against (RM289)."""
    lines = fixed_width_files()[kind]
    first = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    cut = lines[first].rsplit("\t", 1)[0]
    empty = read_coverage(write(tmp_path / "empty.txt", [*lines[:first], cut + "\t", *lines[first + 1:]]), window=None, kind=kind)
    short = read_coverage(write(tmp_path / "short.txt", [*lines[:first], cut, *lines[first + 1:]]), window=None, kind=kind)
    assert short.intervals.collect().equals(empty.intervals.collect())


@pytest.mark.parametrize(
    ("lines", "kind"),
    [
        (["\t".join(["chrom", "length", "bases", "mean", "min", "max"])], None),
        (["#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"], None),
        (["#chrom\tstart\tend\tregion\t1X\t10X"], CoverageFormat.THRESHOLDS_BED),
        (["# a comment, no row"], CoverageFormat.PER_BASE_BED),
        (["# a comment, no row"], CoverageFormat.QUANTIZED_BED),
    ],
)
@pytest.mark.parametrize("window", [None, 10])
def test_a_file_with_no_rows_reads_as_no_intervals(tmp_path: Path, lines: list[str], kind: CoverageFormat | None, window: int | None) -> None:
    """A header or comments and no row read as a facet with no intervals on polars 1.x and 2.x alike: 2.x
    reads an empty frame where 1.x raises `NoDataError` (RM289)."""
    coverage = read_coverage(write(tmp_path / "empty.txt", lines), kind=kind, window=window)
    intervals = coverage.intervals.collect()
    assert intervals.schema == pl.Schema(COVERAGE_INTERVALS) and intervals.height == 0
    if coverage.thresholds is not None:
        assert coverage.thresholds.collect().height == 0


REASONS = {"wide": "more than the format's", "short": "fewer than the format's", "blank": "blank line", "depths": "depth columns"}
"""The phrase that names each reason in a facet's warnings (RM291)."""


def counted(warnings: list[str]) -> dict[str, int]:
    """The rows each reason's warning counts, read from its leading number; one line per reason."""
    out: dict[str, int] = {}
    for line in warnings:
        for reason, phrase in REASONS.items():
            if phrase in line:
                assert reason not in out, f"two lines for {reason}"
                out[reason] = int(line.split(" ", 1)[0].replace(",", ""))
    return out


def expected_counts(lines: list[str], columns: int, kind: CoverageFormat) -> dict[str, int]:
    """The rows of `lines` wider than, short of and blank against `columns`, from the lines as written."""
    widths = [len(line.split("\t")) if line else 0 for line in lines if not line.startswith("#")]
    counts = {
        "depths" if kind is CoverageFormat.SAMTOOLS_DEPTH else "wide": sum(w > columns for w in widths),
        "short": sum(0 < w < columns for w in widths),
        "blank": widths.count(0),
    }
    return {reason: rows for reason, rows in counts.items() if rows}


def shortened(lines: list[str], rows: list[int]) -> list[str]:
    """`lines` with the last field taken off the lines at `rows`."""
    return [line.rsplit("\t", 1)[0] if i in rows else line for i, line in enumerate(lines)]


RAGGED_VARIANTS = ("none", "wide_first", "wide_later", "short", "mixed", "blank")


def ragged_variants(lines: list[str]) -> dict[str, list[str]]:
    """A fixture's lines with rows cut or padded each way RM291 counts. `mixed` widens as many rows by one
    field as it shortens by one, so the file's field total is the format's: a count from totals sees nothing."""
    data = [i for i, line in enumerate(lines) if not line.startswith("#")]
    pairs = min(len(data[0::4]), len(data[2::4]))
    mixed = [f"{line}\textra" if i in data[0::4][:pairs] else line for i, line in enumerate(shortened(lines, data[2::4][:pairs]))]
    return {
        "none": lines,
        "wide_first": widened(lines, data[:1]),
        "wide_later": widened(lines, data[1::2]),
        "short": shortened(lines, data[2::3]),
        "mixed": mixed,
        "blank": [*lines[: data[1]], "", *lines[data[1] :], ""],
    }


@pytest.mark.parametrize("window", [None, 7])
@pytest.mark.parametrize("variant", RAGGED_VARIANTS)
@pytest.mark.parametrize("kind", list(fixed_width_files()))
def test_rows_cut_or_padded_are_counted_by_reason(tmp_path: Path, kind: CoverageFormat, variant: str, window: int | None) -> None:
    """Each row read as its leading fields, or with missing fields as null, is counted under its reason in
    the facet's warnings, one line per reason (RM291); a file whose rows all hold the format's fields
    warns of nothing."""
    plain = fixed_width_files()[kind]
    columns = len(next(line for line in plain if not line.startswith("#")).split("\t"))
    variants = ragged_variants(plain)
    assert tuple(variants) == RAGGED_VARIANTS
    lines = variants[variant]
    warnings = read_coverage(write(tmp_path / "rows.txt", lines), window=window, kind=kind).facet.warnings
    expected = expected_counts(lines, columns, kind)
    assert (variant == "none") == (not expected)
    assert counted(warnings) == expected
    assert len(warnings) == len(expected)


@pytest.mark.parametrize("window", [None, 7])
@pytest.mark.parametrize("kind", [CoverageFormat.PER_BASE_BED, CoverageFormat.SAMTOOLS_DEPTH])
def test_rows_are_counted_across_the_read_s_blocks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: CoverageFormat, window: int | None) -> None:
    """A per-base source read in blocks counts its cut and padded rows over all of them, as one read does."""
    path = write(tmp_path / "rows.txt", ragged_variants(fixed_width_files()[kind])["mixed"])
    whole = read_coverage(path, window=window, kind=kind)
    monkeypatch.setattr(coverage_module, "_CHUNK_BYTES", 64)
    blocks = read_coverage(path, window=window, kind=kind)
    assert blocks.facet.warnings == whole.facet.warnings and counted(whole.facet.warnings)
    assert blocks.intervals.collect().equals(whole.intervals.collect())


def test_samtools_depth_of_several_files_stated_reads_the_first_column_and_says_so(tmp_path: Path) -> None:
    """`samtools depth` over several alignment files, read with `kind` stated, reads the first depth column
    and says how many there were; the header's file names stay out of the warning."""
    depth = sorted(depth_of(READS, contig="chr1").items())
    rows = [f"chr1\t{pos}\t{value}\t{value + 1}" + ("\t7" if i % 3 == 0 else "") for i, (pos, value) in enumerate(depth)]
    several = write(tmp_path / "several.txt", ["#CHROM\tPOS\tperson_a.bam\tperson_b.bam", *rows])
    single = write(tmp_path / "single.txt", [f"chr1\t{pos}\t{value}" for pos, value in depth])
    read = read_coverage(several, window=None, kind=CoverageFormat.SAMTOOLS_DEPTH)
    assert read.intervals.collect().equals(read_coverage(single, window=None, kind=CoverageFormat.SAMTOOLS_DEPTH).intervals.collect())
    columns = sorted({len(row.split("\t")) - 2 for row in rows})
    (warning,) = read.facet.warnings
    assert counted([warning]) == {"depths": len(rows)}
    assert f"carry {columns[0]} to {columns[-1]} depth columns" in warning and "the first was read" in warning
    assert "person" not in warning and ".bam" not in warning
    two = read_coverage(write(tmp_path / "two.txt", [row.rsplit("\t", 1)[0] if row.count("\t") == 4 else row for row in rows]), window=None, kind=CoverageFormat.SAMTOOLS_DEPTH)
    assert "carry 2 depth columns" in two.facet.warnings[0]


def test_wide_rows_past_the_sniffed_lines_leave_the_format_and_the_runs_as_they_are(tmp_path: Path) -> None:
    # The format is told from the first lines alone, so a file whose wide rows come later is read unstated.
    # Rows past the hundredth are past the lines polars infers a width from, which polars 1.x refused too.
    lines = [f"chr1\t{i}\t{i + 1}\t{i % 4}" for i in range(150)]
    mixed = widened(lines, list(range(120, 150, 3)))
    plain = read_coverage(write(tmp_path / "plain.bed", lines), window=None)
    read = read_coverage(write(tmp_path / "mixed.bed", mixed), window=None)
    assert read.facet.format_id is plain.facet.format_id is CoverageFormat.PER_BASE_BED
    assert read.intervals.collect().equals(plain.intervals.collect())


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


def test_an_unindexed_bam_rewritten_with_other_reads_is_not_served_the_old_spill(tmp_path: Path, depth_calls: list[str]) -> None:
    """Without an index, a BAM's header and byte size are all that named it, so a BAM rewritten in
    place with the same header and size but other reads was served the first one's windows, under its
    facet id. An unindexed BAM is now known by its bytes' sha256: one more read of the file, against a
    depth pass on every call had it not been spilled at all. An indexed one keeps its header, index and
    size."""
    path = tmp_path / "a.bam"
    first = alignment_coverage(write_bam(path, READS, index=False), window=10)
    size = path.stat().st_size
    assert first.facet.source_sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    moved = None
    for start in range(1, 30):  # the X read moved along until the rewritten file has the same size
        reads = [*READS[:-1], ("chrX", start, "10M")]
        write_bam(path, reads, index=False)
        if path.stat().st_size == size:
            moved = reads
            break
    assert moved is not None, "no X read start gives a BAM of the first one's size"
    second = alignment_coverage(path, window=10)
    assert depth_calls == [str(path), str(path)]
    assert second.facet.facet_id != first.facet.facet_id
    assert rows_of(second, "X") == windows_of(depth_of(moved, contig="chrX"), 50, 10) != rows_of(first, "X")
    # Indexed, the file is not hashed: the index names it.
    indexed = alignment_coverage(write_bam(tmp_path / "b.bam", READS), window=10)
    assert indexed.facet.source_sha256 is None


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


# --- Quantized runs from the same depth pass (RM278) ---------------------------------------------

# chrX's last base (50) is reached by a read of its own: its depth (1) differs from the base before it (0).
QREADS = [*READS, ("chrX", 40, "5M"), ("chrX", 49, "1M")]


def quantized_lines(reads: list[tuple[str, int, str]], edges: list[int | None]) -> list[str]:
    """mosdepth's `.quantized.bed` for `edges` by what its bins mean (`linear_search`, `make_lookup` in
    mosdepth.nim): each base in the bin `[lower, upper)` its depth falls in, runs of one bin joined, a
    base outside every bin in no run, bins labelled `lower:upper` or `lower:inf`."""
    pairs = list(pairwise(edges))
    labels = [f"{low}:{'inf' if high is None else high}" for low, high in pairs]

    def bin_of(depth: int) -> int | None:
        return next((i for i, (low, high) in enumerate(pairs) if depth >= low and (high is None or depth < high)), None)

    out = []
    for source, length in CONTIGS:
        depth = depth_of(reads, contig=source)
        current, start = None, 0
        for base in range(1, length + 2):
            found = bin_of(depth.get(base, 0)) if base <= length else None
            if found == current:
                continue
            if current is not None:
                out.append(f"{source}\t{start - 1}\t{base - 1}\t{labels[current]}")
            current, start = found, base
    return out


@pytest.mark.parametrize("edges", [(0, 1, 3, None), (1, 3), (0, 2), (2, None)])
def test_quantized_runs_are_what_read_coverage_gives_for_mosdepths_file(tmp_path: Path, edges: tuple[int | None, ...]) -> None:
    bam = write_bam(tmp_path / "q.bam", QREADS)
    found = alignment_coverage(bam, window=10, quantize=edges)
    mosdepth = read_coverage(write(tmp_path / "q.quantized.bed.gz", quantized_lines(QREADS, list(edges)), compress=True))
    assert (mosdepth.facet.format_id, mosdepth.facet.resolution) == (CoverageFormat.QUANTIZED_BED, Resolution.BIN)
    assert found.quantized.intervals.collect().equals(mosdepth.intervals.collect())
    facet = found.quantized.facet
    assert (facet.format_id, facet.resolution, facet.window, facet.zero_filled, facet.quantize) == (
        CoverageFormat.ALIGNMENT_DEPTH, Resolution.BIN, None, False, list(edges),
    )
    assert (facet.options, facet.contig_lengths, facet.source_sha256) == (found.facet.options, found.facet.contig_lengths, found.facet.source_sha256)
    assert facet.facet_id != found.facet.facet_id
    # The windows of the same pass are those of a call without bins, under the same facet id.
    plain = alignment_coverage(bam, window=10)
    assert found.intervals.collect().equals(plain.intervals.collect())
    assert found.facet == plain.facet and plain.quantized is None


def test_batch_edges_change_no_window_and_no_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, depth_calls: list[str]) -> None:
    bam = write_bam(tmp_path / "q.bam", QREADS)
    whole = alignment_coverage(bam, window=10, quantize=(0, 1, 3, None))
    clear_cache(source=bam)
    monkeypatch.setattr(coverage_module, "_BATCH_ROWS", 7)  # cuts runs, windows and contigs, many times
    cut = alignment_coverage(bam, window=10, quantize=(0, 1, 3, None))
    assert depth_calls == [str(bam)] * 2
    assert cut.intervals.collect().equals(whole.intervals.collect())
    assert cut.quantized.intervals.collect().equals(whole.quantized.intervals.collect())


def test_a_contigs_last_base_takes_its_own_bin(tmp_path: Path) -> None:
    """mosdepth's `gen_quantized` loops to the second-to-last base and ends its last run at the
    contig's end, so a last base whose bin differs from the base before it goes into that run. The bin
    here is the base's own (docs/SCHEMA.md, the lead's call)."""
    found = alignment_coverage(write_bam(tmp_path / "q.bam", QREADS), window=10, quantize=(0, 1, None)).quantized
    rows = found.intervals.filter(pl.col("chrom") == "X").select("pos", "end_pos", "label").collect().rows()
    assert rows[-2:] == [(46, 49, "0:1"), (50, 50, "1:inf")]


def test_the_quantized_runs_come_from_the_one_depth_pass(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`depth_calls` counts `depth()` calls, not reads of the file: each execution of polars-bio's
    source reads it (`py_read_table`), so that is what is counted here."""
    executions: list[str] = []
    real = pb_native.py_read_table

    def counted(*args, **kwargs):
        executions.append(args[1])
        return real(*args, **kwargs)

    monkeypatch.setattr(pb_native, "py_read_table", counted)
    bam = write_bam(tmp_path / "q.bam", QREADS)
    alignment_coverage(bam, window=10)
    plain = len(executions)
    assert plain == 1
    alignment_coverage(write_bam(tmp_path / "r.bam", QREADS), window=10, quantize=(0, 1, 3, None))
    assert len(executions) == plain + 1
    # Read again from the spill: no execution at all.
    alignment_coverage(tmp_path / "r.bam", window=10, quantize=(0, 1, 3, None))
    assert len(executions) == plain + 1


def test_quantized_runs_are_spilled_under_their_bins(tmp_path: Path, depth_calls: list[str]) -> None:
    bam = write_bam(tmp_path / "q.bam", QREADS)
    first = alignment_coverage(bam, window=10, quantize=(0, 1, None))
    again = alignment_coverage(bam, window=10, quantize=[0, 1, None])
    assert depth_calls == [str(bam)]
    assert again.quantized.facet == first.quantized.facet
    assert again.quantized.intervals.collect().equals(first.quantized.intervals.collect())
    # Other bins, and no bins, are other spills: the file is read for each (the bins key the spill).
    other = alignment_coverage(bam, window=10, quantize=(0, 2, None))
    alignment_coverage(bam, window=10)
    assert depth_calls == [str(bam)] * 3
    assert other.quantized.facet.facet_id != first.quantized.facet.facet_id
    assert len(coverage_manifests()) == 3


def test_a_wrapped_depth_is_in_no_quantized_run(tmp_path: Path) -> None:
    deep = [("chr1", 40, "5M")] * 32_770
    found = alignment_coverage(write_bam(tmp_path / "deep.bam", [*READS, *deep]), window=10, quantize=(0, 1, None))
    runs = found.quantized.intervals.filter(pl.col("chrom") == "1").collect()
    covered = {base for row in runs.iter_rows(named=True) for base in range(row["pos"], row["end_pos"] + 1)}
    assert covered == set(range(1, 101)) - set(range(41, 46))
    (warning,) = found.quantized.facet.warnings
    assert warning == found.facet.warnings[0]
    assert warning.startswith("5 base(s) in 1 window(s) came back with a negative depth")
    assert warning.endswith("no quantized run covers them, so their bin is unknown")


@pytest.mark.parametrize(
    ("edges", "match"),
    [((), "two or more"), ((0,), "two or more"), ((1, 1), "rise strictly"), ((2, 1), "rise strictly"), ((0, None, 5), "None only as the last"),
     ((-1, 2), "not negative"), ((True, 2), "whole numbers"), ((0, 1.5), "whole numbers")],
)
def test_quantize_edges_are_checked(tmp_path: Path, edges: tuple, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        alignment_coverage(tmp_path / "never-read.bam", quantize=edges)
