import gzip
import hashlib
from importlib.metadata import version
from pathlib import Path

import polars as pl
import pytest

from just_vcf import (
    ContigNames,
    Loaded,
    UnsortedRowsError,
    UnwritableRowsError,
    load,
    rejoin_split,
    vcf_samples,
    write_vcf,
)
from just_vcf.bgzf import EOF_BLOCK
from just_vcf.sniff import sniff

FORMAT_LINES = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    "PL": '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "GL": '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods">',
    "GP": '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posteriors">',
    "DS": '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">',
    "PS": '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "FT": '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">',
    "MIN_DP": '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Block minimum depth">',
}

CONTIGS = ("##contig=<ID=1,length=249250621>", "##contig=<ID=2,length=243199373>", "##contig=<ID=MT,length=16569>")

# Every shape the writer has to carry back: a multi-allelic het with a six-entry PL, an indel, a
# monomorphic record, haploid and half calls, a whole no-call, a phased pair with PS, a phased row in
# the implicit set, GP posteriors, a dosage, ids, filters, a sample filter, and a row with fields left out.
RICH = [
    "1\t100\trs1;rs2\tA\tG,T\t50\tq10;s50\t.\tGT:GQ:DP:AD:PL:FT\t1/2:30:7:0,3,4:90,40,0,50,20,80:LowGQ",
    "1\t200\t.\tACG\tA\t12.5\tPASS\t.\tGT:GQ:DP:AD:PL\t0/1:20:9:5,4:30,0,60",
    "1\t300\t.\tT\t.\t.\t.\t.\tGT:GQ:DP\t0/0:.:11",
    "1\t400\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:400",
    "1\t500\t.\tC\tT\t.\tPASS\t.\tGT:PS\t1|0:.",
    "1\t600\t.\tG\tC\t.\t.\t.\tGT:GP:DS\t0/1:0.1,0.8,0.1:1",
    "1\t700\t.\tG\tC,A\t.\t.\t.\tGT\t1/.",
    "1\t800\t.\tG\tC\t.\t.\t.\tGT\t./.",
    "2\t100\t.\tA\tAT\t3\t.\t.\tGT:GL\t1/1:-3.2,-0.3,0",
    "2\t200\t.\tT\tC\t.\t.\t.\tGQ:DP\t40:15",
    "MT\t150\t.\tC\tT\t99\tPASS\t.\tGT:AD\t1:0,40",
]
RICH_FORMAT = ("GT", "GQ", "DP", "AD", "PL", "GL", "GP", "DS", "PS", "FT")

KEY = ["chrom", "pos", "source_ref", "alts", "gt", "gt_alleles"]
VALUES = ["ids", "end_pos", "qual", "filters", "gq", "dp", "ad", "gl", "gp", "dosage", "sample_filters", "phase_sets", "source_alleles"]


def vcf_text(records: list[str], fmt: tuple[str, ...] = ("GT",), info: tuple[str, ...] = (), contigs: tuple[str, ...] = CONTIGS,
             version_line: str = "VCFv4.2", sample: str = "SYNTHETIC_SAMPLE") -> str:
    meta = [f"##fileformat={version_line}", *contigs]
    meta += ['##INFO=<ID=END,Number=1,Type=Integer,Description="End">'] if "END" in info else []
    meta += [FORMAT_LINES[key] for key in fmt]
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
    return "\n".join([*meta, "\t".join(columns), *records]) + "\n"


def loaded_from(tmp_path: Path, records: list[str], name: str = "in.vcf", **kwargs) -> Loaded:
    path = tmp_path / name
    path.write_text(vcf_text(records, **kwargs))
    return load(path)


def as_set(frame: pl.DataFrame, columns: list[str]) -> set[tuple]:
    """Rows as hashable tuples; lists become tuples."""
    return {tuple(tuple(v) if isinstance(v, list) else v for v in row) for row in frame.select(columns).iter_rows()}


def plain_text(path: Path) -> str:
    return gzip.decompress(path.read_bytes()).decode()


def test_round_trip_keeps_sites_genotypes_and_format_values(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    written = write_vcf(loaded, tmp_path / "out.vcf.gz")
    original = loaded.calls.collect()
    back = load(written.path).calls.collect()

    assert written.rows_written == original.height == back.height
    assert written.rows_skipped == {} and written.blocks_not_written == 0 and written.warnings == ()
    assert as_set(back, KEY) == as_set(original, KEY)
    assert as_set(back, KEY + VALUES) == as_set(original, KEY + VALUES)
    # Phase: VCF 4.3 has no leading indicator, so the first haplotype's phase is the implicit one on both sides.
    assert as_set(back, KEY + ["phased", "phase_set_kinds"]) == as_set(original, KEY + ["phased", "phase_set_kinds"])
    # The GT text itself comes back as the source wrote it, for every row that had one.
    assert as_set(back.filter(pl.col("source_gt").is_not_null()), KEY + ["source_gt"]) == as_set(
        original.filter(pl.col("source_gt").is_not_null()), KEY + ["source_gt"]
    )


def test_two_writes_of_one_frame_are_byte_identical(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    first = write_vcf(loaded, tmp_path / "a.vcf.gz")
    second = write_vcf(loaded, tmp_path / "b.vcf.gz")
    assert first.sha256 == second.sha256 == hashlib.sha256((tmp_path / "a.vcf.gz").read_bytes()).hexdigest()
    assert (tmp_path / "a.vcf.gz").read_bytes() == (tmp_path / "b.vcf.gz").read_bytes()
    # A reload of the same file is the same frames, so the same bytes again.
    third = write_vcf(load(tmp_path / "in.vcf"), tmp_path / "c.vcf.gz")
    assert third.sha256 == first.sha256


def test_output_is_bgzf_unless_asked_for_plain_text(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    packed = write_vcf(loaded, tmp_path / "out.vcf.gz").path.read_bytes()
    assert packed[:4] == b"\x1f\x8b\x08\x04" and packed[12:14] == b"BC" and packed.endswith(EOF_BLOCK)
    plain = write_vcf(loaded, tmp_path / "out.vcf", compress=False)
    assert plain.path.read_text() == gzip.decompress(packed).decode()
    assert plain.sha256 == hashlib.sha256(plain.path.read_bytes()).hexdigest()


def test_header_declares_every_key_written_and_names_the_facet(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    text = plain_text(write_vcf(loaded, tmp_path / "out.vcf.gz").path)
    meta = [line for line in text.splitlines() if line.startswith("##")]
    body = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
    facet = loaded.facets.collect().row(0, named=True)

    assert meta[0] == "##fileformat=VCFv4.3"
    assert f"##source=just-vcf {version('just-vcf')}" in meta
    assert f"##just_vcf_facet_id={facet['facet_id']}" in meta
    assert f"##just_vcf_source_sha256={facet['source_sha256']}" in meta
    assert str(tmp_path) not in text, "the source path is a label that can carry a name; never written"
    assert text.splitlines()[len(meta)].split("\t")[9] == facet["sample"]

    used_format = {key for row in body for key in row[8].split(":")}
    declared_format = {line.split("ID=")[1].split(",")[0] for line in meta if line.startswith("##FORMAT=")}
    assert used_format == declared_format == {"GT", "GQ", "DP", "AD", "GL", "GP", "DS", "FT", "PS"}
    used_filters = {f for row in body for f in row[6].split(";") if f != "."}
    used_filters |= {row[9].split(":")[row[8].split(":").index("FT")] for row in body if "FT" in row[8].split(":")}
    declared_filters = {line.split("ID=")[1].split(",")[0] for line in meta if line.startswith("##FILTER=")}
    assert used_filters <= declared_filters and declared_filters - used_filters == {"PASS"} - used_filters
    declared_contigs = [line.split("ID=")[1].split(",")[0].rstrip(">") for line in meta if line.startswith("##contig=")]
    assert declared_contigs == ["1", "2", "MT"]
    assert [row[0] for row in body] == sorted([row[0] for row in body], key=declared_contigs.index), "rows follow the header's contig order"
    assert not any(line.startswith("##INFO=") for line in meta), "no INFO written, none declared"


def test_unresolved_ref_is_refused_with_its_count(tmp_path: Path) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(
        "# This data file generated by 23andMe\n# rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t100\tAG\nrs2\t1\t200\tCC\nrs3\tMT\t300\tT\n"
    )
    loaded = load(path)
    calls = loaded.calls.collect()
    expected = calls.filter(pl.col("ref").is_null() & pl.col("source_ref").is_null()).height
    assert expected == calls.height > 0
    with pytest.raises(UnwritableRowsError) as refused:
        write_vcf(loaded, tmp_path / "out.vcf.gz")
    assert refused.value.counts == {"ref_unresolved": expected}
    assert not (tmp_path / "out.vcf.gz").exists() and not (tmp_path / "out.vcf.gz.partial").exists()

    written = write_vcf(loaded, tmp_path / "out.vcf.gz", resolved_only=True)
    assert written.rows_written == 0 and written.rows_skipped == {"ref_unresolved": expected}
    assert any(f"{expected} row(s) not written" in warning for warning in written.warnings)
    header = plain_text(written.path).splitlines()
    assert header[-1].split("\t")[9] == loaded.facets.collect()["facet_id"][0], "an array names no sample: the facet id stands in"


def test_resolved_only_writes_the_rest_and_counts_each_reason(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    calls = loaded.calls.collect()
    unresolved = pl.col("pos") == 100  # on chromosomes 1 and 2
    mismatch = (pl.col("chrom") == "1") & (pl.col("pos") == 200)
    no_alts = (pl.col("chrom") == "1") & (pl.col("pos") == 300)
    edited = calls.with_columns(
        pl.when(unresolved).then(None).otherwise(pl.col("source_ref")).alias("source_ref"),
        pl.when(mismatch).then(pl.lit("mismatch")).otherwise(pl.col("ref_check")).cast(calls.schema["ref_check"]).alias("ref_check"),
        pl.when(no_alts).then(None).otherwise(pl.col("alts")).alias("alts"),
    )
    damaged = Loaded(facets=loaded.facets, calls=edited.lazy(), blocks=loaded.blocks)
    expected = {
        "ref_unresolved": calls.filter(unresolved).height,
        "ref_mismatch": calls.filter(mismatch).height,
        "alts_unresolved": calls.filter(no_alts).height,
    }
    assert expected == {"ref_unresolved": 2, "ref_mismatch": 1, "alts_unresolved": 1}
    with pytest.raises(UnwritableRowsError) as refused:
        write_vcf(damaged, tmp_path / "out.vcf.gz")
    assert refused.value.counts == expected
    written = write_vcf(damaged, tmp_path / "out.vcf.gz", resolved_only=True)
    assert written.rows_skipped == expected
    assert written.rows_written == calls.height - sum(expected.values())
    back = load(written.path).calls.collect()
    assert as_set(back, KEY) == as_set(calls.filter(~(unresolved | mismatch | no_alts)), KEY)


def test_vectors_indexed_over_other_alleles_are_left_out(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    calls = loaded.calls.collect()
    moved = calls.with_columns(
        pl.when(pl.col("pos") == 200).then(pl.lit(["ACG", "C"])).otherwise(pl.col("source_alleles")).alias("source_alleles")
    )
    written = write_vcf(Loaded(facets=loaded.facets, calls=moved.lazy(), blocks=loaded.blocks), tmp_path / "out.vcf.gz")
    assert any(w.startswith("1 row(s) had GL/GP/AD") for w in written.warnings)
    row = load(written.path).calls.filter(pl.col("pos") == 200).collect().row(0, named=True)
    assert row["gl"] is None and row["ad"] is None and row["gq"] == 20, "only the vectors go; the scalars stay"


@pytest.mark.parametrize(
    ("style", "expected"),
    [
        (ContigNames.STORED, ["1", "2", "MT"]),
        (ContigNames.CHR, ["chr1", "chr2", "chrM"]),
        (ContigNames.SOURCE, ["chr1", "chr2", "chrM"]),
    ],
)
def test_contig_names(tmp_path: Path, style: ContigNames, expected: list[str]) -> None:
    # chr2 is declared and has no rows, as in an exome or a region extract.
    contigs = ("##contig=<ID=chr1,length=249250621>", "##contig=<ID=chr2,length=243199373>", "##contig=<ID=chrM,length=16571>")
    records = ["chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1", "chrM\t200\t.\tC\tT\t.\t.\t.\tGT\t1"]
    loaded = loaded_from(tmp_path, records, contigs=contigs)
    assert loaded.calls.collect()["chrom"].to_list() == ["1", "MT"]
    text = plain_text(write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=style).path)
    assert [line.split("\t")[0] for line in text.splitlines() if not line.startswith("#")] == [expected[0], expected[2]]
    assert [line for line in text.splitlines() if line.startswith("##contig=")] == [
        f"##contig=<ID={expected[0]},length=249250621>",
        f"##contig=<ID={expected[1]},length=243199373>",
        f"##contig=<ID={expected[2]},length=16571>",
    ]
    assert load(tmp_path / "out.vcf.gz").calls.collect()["chrom"].to_list() == ["1", "MT"]


def test_chr_style_from_plain_source_maps_the_mitochondrion_by_name(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1", "MT\t200\t.\tC\tT\t.\t.\t.\tGT\t1"]
    contigs = ("##contig=<ID=1,length=249250621>", "##contig=<ID=MT,length=16569>", "##contig=<ID=GL000192.1,length=547496>")
    loaded = loaded_from(tmp_path, records, contigs=contigs)
    text = plain_text(write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=ContigNames.CHR).path)
    ids = [line.split("ID=")[1].split(",")[0] for line in text.splitlines() if line.startswith("##contig=")]
    assert ids == ["chr1", "chrM", "GL000192.1"], "MT is chrM, never chrMT; an unplaced contig keeps its name"


def test_split_rows_are_written_as_held_and_tagged(tmp_path: Path) -> None:
    records = [
        "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/0",
        "1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1",
        "1\t300\t.\tC\tG,T\t50\tPASS\t.\tGT\t1/2",
    ]
    loaded = loaded_from(tmp_path, records)
    flagged = loaded.calls.filter(pl.col("split_signal").is_not_null()).collect()
    written = write_vcf(loaded, tmp_path / "out.vcf.gz")
    assert any(w.startswith(f"{flagged.height} row(s) look split") for w in written.warnings)
    body = [line.split("\t") for line in plain_text(written.path).splitlines() if not line.startswith("#")]
    tagged = {(row[1], row[4]) for row in body if "JUST_VCF_SPLIT=" in row[7]}
    assert tagged == {(str(r["pos"]), ",".join(r["alts"])) for r in flagged.iter_rows(named=True)}

    rejoined = write_vcf(rejoin_split(loaded), tmp_path / "joined.vcf.gz")
    assert not any("look split" in w for w in rejoined.warnings)
    assert rejoined.rows_written == written.rows_written - 1


def test_blocks_are_counted_not_dropped_silently(tmp_path: Path) -> None:
    records = [
        "1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ:MIN_DP:PL\t0/0:60:23:0,60,900",
        "1\t114\t.\tT\tTC,<NON_REF>\t213\t.\t.\tGT:GQ:MIN_DP:PL\t0/1:99:.:51,0,36,93,92,86",
        "1\t116\t.\tC\t<*>\t.\t.\tEND=120\tGT:GQ:MIN_DP:PL\t./.:.:.:.",
    ]
    loaded = loaded_from(tmp_path, records, fmt=("GT", "GQ", "MIN_DP", "PL"), info=("END",))
    blocks = loaded.blocks.collect().height
    written = write_vcf(loaded, tmp_path / "out.vcf.gz")
    assert written.blocks_not_written == blocks == 2
    assert any(w.startswith(f"{blocks} reference/no-call block(s) not written") for w in written.warnings)
    assert written.rows_written == loaded.calls.collect().height == 1


def test_symbolic_extent_is_written_as_end(tmp_path: Path) -> None:
    records = ["1\t121\t.\tA\t<*>\t.\t.\tEND=125\tGT:PL\t0/1:40,0,10"]
    loaded = loaded_from(tmp_path, records, fmt=("GT", "PL"), info=("END",))
    text = plain_text(write_vcf(loaded, tmp_path / "out.vcf.gz").path)
    assert text.splitlines()[-1].split("\t")[7] == "END=125"
    assert '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the record">' in text
    assert as_set(load(tmp_path / "out.vcf.gz").calls.collect(), KEY + ["end_pos"]) == as_set(loaded.calls.collect(), KEY + ["end_pos"])


def test_rows_out_of_order_are_refused_and_leave_nothing_behind(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    shuffled = Loaded(facets=loaded.facets, calls=loaded.calls.collect().reverse().lazy(), blocks=loaded.blocks)
    with pytest.raises(UnsortedRowsError):
        write_vcf(shuffled, tmp_path / "out.vcf.gz")
    assert not (tmp_path / "out.vcf.gz").exists() and not (tmp_path / "out.vcf.gz.partial").exists()
    contigs_interleaved = pl.concat([loaded.calls.collect().filter(pl.col("chrom") == c) for c in ("1", "2", "1")])
    with pytest.raises(UnsortedRowsError, match="comes back"):
        write_vcf(Loaded(facets=loaded.facets, calls=contigs_interleaved.lazy(), blocks=loaded.blocks), tmp_path / "out.vcf.gz")


def test_htslib_reads_the_output(tmp_path: Path) -> None:
    """Opt-in external validation: htslib (through pysam, the `[reference]` extra) parses every record."""
    pysam = pytest.importorskip("pysam")
    bcftools = pytest.importorskip("pysam.bcftools")
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    written = write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=ContigNames.CHR)
    calls = loaded.calls.collect()
    with pysam.VariantFile(str(written.path)) as reader:
        records = list(reader)
        sample = list(reader.header.samples)[0]
    assert len(records) == calls.height
    got = {(r.chrom, r.pos, r.ref, r.alts or (), r.samples[sample].get("GT")) for r in records}
    names = {"1": "chr1", "2": "chr2", "MT": "chrM"}
    expected = {
        (names[row["chrom"]], row["pos"], row["source_ref"], tuple(row["alts"]), tuple(row["gt"]) if row["gt"] is not None else None)
        for row in calls.iter_rows(named=True)
    }
    assert got == expected
    viewed = bcftools.view("--no-header", str(written.path), catch_stdout=True)
    assert len(viewed.splitlines()) == calls.height


def test_one_facet_only(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    two = Loaded(facets=pl.concat([loaded.facets, loaded.facets]), calls=loaded.calls, blocks=loaded.blocks)
    with pytest.raises(ValueError, match="one facet"):
        write_vcf(two, tmp_path / "out.vcf.gz")


def test_real_vcf_upload_round_trips(pgp_participants: Path, tmp_path: Path) -> None:
    """Opt-in: a small PGP-HMS VCF upload comes back from its own export with the same sites and genotypes."""
    candidates = sorted(pgp_participants.glob("*/upload-*/*.vcf*"), key=lambda p: p.stat().st_size)
    found = next((p for p in candidates[:200] if p.stat().st_size > 10_000 and sniff(p).format == "vcf"), None)
    if found is None:
        pytest.skip("no plain VCF upload found")
    samples = vcf_samples(found)
    loaded = load(found, sample=samples[0] if len(samples) > 1 else None)
    written = write_vcf(loaded, tmp_path / "out.vcf.gz", resolved_only=True)
    original, back = loaded.calls.collect(), load(written.path).calls.collect()
    skipped = sum(written.rows_skipped.values())
    assert written.rows_written == back.height == original.height - skipped
    if skipped:
        assert as_set(back, KEY) <= as_set(original, KEY)
    else:
        assert as_set(back, KEY) == as_set(original, KEY)
