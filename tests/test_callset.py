"""Callset type and what absence means, from the spread of a VCF facet's calls (RM4)."""

import json
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff
from just_vcf.callset import (
    FLANK_BP,
    MAX_PANEL_SITES,
    MAX_TARGETED_BREADTH,
    MAX_WES_SITES,
    MIN_PANEL_SITES,
    MIN_WES_SITES,
    MIN_WGS_BREADTH,
    MIN_WGS_SITES,
    absence_means,
    classify,
    imputation_marker,
    spread,
)
from just_vcf.schema import AbsenceMeans, CallsetType

AUTOSOMES = [str(n) for n in range(1, 23)]
READS = ["GT", "GQ", "DP", "AD", "PL"]


def sites(rows: list[tuple[str, int]], alts: list[str] | None = None) -> pl.LazyFrame:
    return pl.LazyFrame(
        {"chrom": [c for c, _ in rows], "pos": [p for _, p in rows], "alts": [alts if alts is not None else ["G"]] * len(rows)},
        schema={"chrom": pl.String, "pos": pl.UInt32, "alts": pl.List(pl.String)},
    )


def genome_shaped(total: int) -> pl.LazyFrame:
    """`total` sites (rounded down to a multiple of 22) spread evenly over 22 autosomes of 100 Mb each."""
    per = total // len(AUTOSOMES)
    step = 100_000_000 // per
    positions = pl.int_range(0, per, eager=True).cast(pl.UInt32) * step + 1
    return pl.concat(
        pl.LazyFrame({"chrom": pl.repeat(chrom, per, eager=True), "pos": positions, "alts": pl.repeat(["G"], per, eager=True)})
        for chrom in AUTOSOMES
    )


def clustered(total_per_cluster: int, clusters: int = 25) -> pl.LazyFrame:
    """`clusters` dense clusters (sites 100 bp apart) 4 Mb apart on each autosome, built in Polars."""
    within = pl.int_range(0, total_per_cluster, eager=True).cast(pl.UInt32) * 100
    frames = [
        pl.LazyFrame({"chrom": pl.repeat(chrom, total_per_cluster, eager=True), "pos": within + 1_000_000 + g * 4_000_000})
        for chrom in AUTOSOMES
        for g in range(clusters)
    ]
    return pl.concat(frames).with_columns(pl.lit(["G"]).alias("alts"))


def exome_shaped(genes: int = 25, calls: int = 20, chromosomes: list[str] = AUTOSOMES) -> list[tuple[str, int]]:
    """Clusters of `calls` sites 100 bp apart, one per gene, genes 4 Mb apart on each chromosome."""
    return [(chrom, 1_000_000 + g * 4_000_000 + c * 100) for chrom in chromosomes for g in range(genes) for c in range(calls)]


def test_spread_is_the_ported_breadth() -> None:
    # Two contigs: gaps of 5 kb (counted whole) and 50 kb (counted as 2 x flank), spans 10 kb and 50 kb.
    measured = spread(sites([("1", 1000), ("1", 6000), ("1", 11000), ("2", 100), ("2", 50100)]))
    assert measured["sites"] == 5
    assert measured["reach_bp"] == 5000 + 5000 + 2 * FLANK_BP
    assert measured["breadth"] == pytest.approx((10000 + 2 * FLANK_BP) / (10000 + 50000))
    assert measured["autosomes"] == 2


def test_spread_counts_distinct_primary_small_variant_sites() -> None:
    rows = [("1", 100), ("1", 100), ("chrUn_x", 5), ("1", 300)]
    frame = pl.concat([sites(rows), sites([("1", 500)], alts=["<DEL>"]), sites([("1", 700)], alts=[])])
    measured = spread(frame)
    # The duplicate position counts once, the scaffold and the symbolic-only row not at all, and the
    # monomorphic record (ALT `.`) does.
    assert measured["sites"] == 3


def test_a_genome_shaped_callset_is_wgs() -> None:
    inferred = classify(genome_shaped(MIN_WGS_SITES + 22), READS, None)
    assert inferred and inferred.callset_type is CallsetType.WGS
    assert inferred.sites >= MIN_WGS_SITES and inferred.breadth >= MIN_WGS_BREADTH


def test_genome_shaped_but_short_of_the_wgs_floor_is_unknown() -> None:
    inferred = classify(genome_shaped(600_000), READS, None)
    assert not inferred and inferred.callset_type is None
    assert inferred.sites == 600_000 - 600_000 % 22 and inferred.breadth > MIN_WGS_BREADTH


def test_millions_of_clustered_sites_are_not_wgs() -> None:
    inferred = classify(clustered(2000), READS, None)
    assert inferred.sites > MIN_WGS_SITES and inferred.breadth < MAX_TARGETED_BREADTH
    assert inferred.callset_type is None


def test_an_exome_shaped_callset_is_wes() -> None:
    inferred = classify(sites(exome_shaped()), READS, None)
    assert inferred.callset_type is CallsetType.WES
    assert MIN_WES_SITES <= inferred.sites <= MAX_WES_SITES and inferred.breadth <= MAX_TARGETED_BREADTH


def test_an_exome_missing_an_autosome_is_unknown() -> None:
    inferred = classify(sites(exome_shaped(genes=26, chromosomes=AUTOSOMES[:-1])), READS, None)
    assert inferred.sites >= MIN_WES_SITES and inferred.autosomes == 21
    assert inferred.callset_type is None


def test_a_few_genes_are_a_panel_and_a_handful_of_sites_is_nothing() -> None:
    panel = [(chrom, 5_000_000 + c * 150) for chrom in ("7", "13", "17") for c in range(40)]
    inferred = classify(sites(panel), READS, None)
    assert inferred.callset_type is CallsetType.PANEL and MIN_PANEL_SITES <= inferred.sites <= MAX_PANEL_SITES
    few = classify(sites(panel[: MIN_PANEL_SITES - 1]), READS, None)
    assert few.callset_type is None and few.sites == MIN_PANEL_SITES - 1
    assert classify(sites([]), READS, None).callset_type is None


def test_a_panel_too_large_for_a_panel_and_too_small_for_an_exome_is_unknown() -> None:
    between = exome_shaped(genes=12, calls=30)
    assert MAX_PANEL_SITES < len(between) < MIN_WES_SITES
    assert classify(sites(between), READS, None).callset_type is None


def test_without_read_evidence_no_sequencing_class_is_named() -> None:
    inferred = classify(genome_shaped(MIN_WGS_SITES + 22), ["GT", "GP", "DS"], None)
    assert inferred.callset_type is None
    assert "read-evidence" in inferred.reason


def test_an_imputation_marker_wins_over_the_spread() -> None:
    inferred = classify(genome_shaped(MIN_WGS_SITES + 22), READS, "the header names an imputation tool")
    assert inferred.callset_type is CallsetType.IMPUTED


@pytest.mark.parametrize(
    ("meta", "info", "found"),
    [
        (["##source=Minimac4.v1.0.2"], [], True),
        (['##source="beagle.27Jan18.7e1.jar"'], [], True),
        (["##source=loimpute"], [], True),
        (["##source=GLIMPSE_phase v1.1.1"], [], True),
        ([], ["AF", "R2", "ER2"], True),
        ([], ["DR2", "AF", "IMP"], True),
        # DRAGEN declares R2_5P_bias, and a caller's ##source is no imputation tool.
        (["##source=DRAGEN", "##source=SelectVariants"], ["R2_5P_bias", "DP"], False),
        (["##imputation_note=minimac"], [], False),
    ],
)
def test_imputation_marker(meta: list[str], info: list[str], found: bool) -> None:
    assert (imputation_marker(meta, info) is not None) is found


@pytest.mark.parametrize(
    ("callset_type", "has_blocks", "expected"),
    [
        (CallsetType.WGS, False, AbsenceMeans.WEAK_REFERENCE),
        (CallsetType.WGS, True, AbsenceMeans.BLOCKS),
        (CallsetType.WES, True, AbsenceMeans.BLOCKS),
        (CallsetType.WES, False, AbsenceMeans.UNKNOWN),
        (CallsetType.PANEL, False, AbsenceMeans.UNKNOWN),
        (CallsetType.IMPUTED, False, AbsenceMeans.UNKNOWN),
        (CallsetType.ARRAY, False, AbsenceMeans.UNKNOWN),
        (None, False, AbsenceMeans.UNKNOWN),
    ],
)
def test_absence_means_follows_the_type_and_blocks(callset_type, has_blocks: bool, expected: AbsenceMeans) -> None:
    assert absence_means(callset_type, has_blocks) is expected


@pytest.mark.parametrize(
    ("callset_type", "has_blocks", "expected"),
    [
        # CGAPipeline's VCF: what it omits it called reference, beside its no-call blocks.
        (CallsetType.WGS, True, AbsenceMeans.REFERENCE),
        (CallsetType.WGS, False, AbsenceMeans.REFERENCE),
        # Rows that do not show a whole genome (a region cut, or too few to tell) keep the other rules.
        (CallsetType.PANEL, True, AbsenceMeans.BLOCKS),
        (None, True, AbsenceMeans.BLOCKS),
        (CallsetType.WES, False, AbsenceMeans.UNKNOWN),
    ],
)
def test_absence_is_a_reference_call_where_the_producer_omits_only_reference(callset_type, has_blocks: bool, expected: AbsenceMeans) -> None:
    assert absence_means(callset_type, has_blocks, omitted_called_reference=True) is expected


HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={chrom},length=200000000>" for chrom in AUTOSOMES],
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
]
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"


def write_vcf(tmp_path: Path, positions: list[tuple[str, int]], alt: str = "G", meta: tuple[str, ...] = (), info: str = ".") -> Path:
    rows = [f"{chrom}\t{pos}\t.\tA\t{alt}\t50\tPASS\t{info}\tGT:DP\t0/1:30" for chrom, pos in positions]
    path = tmp_path / "calls.vcf"
    path.write_text("\n".join([*HEADER, *meta, COLUMNS, *rows]) + "\n")
    return path


def loaded_facet(path: Path) -> dict:
    return load(path).facets.collect().row(0, named=True)


def test_load_types_a_genome_vcf_and_reads_its_absences_as_weak_reference(tmp_path: Path) -> None:
    rows = genome_shaped(MIN_WGS_SITES + 22).collect()
    path = tmp_path / "calls.vcf"
    path.write_text("\n".join([*HEADER, COLUMNS]) + "\n")
    with path.open("a") as handle:
        rows.select(
            "chrom",
            "pos",
            pl.lit(".").alias("id"),
            pl.lit("A").alias("ref"),
            pl.lit("G").alias("alt"),
            pl.lit("50").alias("qual"),
            pl.lit("PASS").alias("filter"),
            pl.lit(".").alias("info"),
            pl.lit("GT:DP").alias("format"),
            pl.lit("0/1:30").alias("S1"),
        ).write_csv(handle, separator="\t", include_header=False)
    row = loaded_facet(path)
    assert (row["callset_type"], row["absence_means"]) == ("wgs", "weak_reference")


def test_load_types_an_exome_vcf(tmp_path: Path) -> None:
    row = loaded_facet(write_vcf(tmp_path, exome_shaped()))
    assert (row["callset_type"], row["absence_means"]) == ("wes", "unknown")


def test_load_types_a_panel_vcf_and_leaves_a_tiny_one_unknown(tmp_path: Path) -> None:
    (tmp_path / "panel").mkdir()
    (tmp_path / "tiny").mkdir()
    panel = [(chrom, 5_000_000 + c * 150) for chrom in ("7", "13", "17") for c in range(40)]
    assert loaded_facet(write_vcf(tmp_path / "panel", panel))["callset_type"] == "panel"
    tiny = loaded_facet(write_vcf(tmp_path / "tiny", panel[:3]))
    assert (tiny["callset_type"], tiny["absence_means"]) == (None, "unknown")


def test_load_does_not_read_a_structural_variant_callset_as_any_class(tmp_path: Path) -> None:
    meta = ('##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">', '##ALT=<ID=DEL,Description="Deletion">')
    row = loaded_facet(write_vcf(tmp_path, exome_shaped(), alt="<DEL>", meta=meta, info="SVLEN=-50"))
    assert row["callset_type"] is None


def test_load_types_imputation_output(tmp_path: Path) -> None:
    row = loaded_facet(write_vcf(tmp_path, exome_shaped(), meta=("##source=Minimac4.v1.0.2",)))
    assert (row["callset_type"], row["absence_means"]) == ("imputed", "unknown")


def test_load_keeps_blocks_deciding_absence_on_a_typed_gvcf(tmp_path: Path) -> None:
    meta = ('##INFO=<ID=END,Number=1,Type=Integer,Description="End">', '##ALT=<ID=*,Description="Any other allele">')
    calls = [f"{chrom}\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30" for chrom, pos in exome_shaped()]
    block = "1\t10\t.\tA\t<*>\t.\t.\tEND=500\tGT:DP\t0/0:20"
    path = tmp_path / "calls.g.vcf"
    path.write_text("\n".join([*HEADER, *meta, COLUMNS, block, *calls]) + "\n")
    row = loaded_facet(path)
    assert (row["callset_type"], row["absence_means"]) == ("wes", "blocks")


def test_real_exome_uploads_are_wes(pgp_participants: Path) -> None:
    """PGP-HMS uploads labelled Genos exome (consented-public): the spread alone names them `wes`."""
    manifest = pgp_participants / "manifest.jsonl"
    if not manifest.exists():
        pytest.skip("no participants manifest")
    exomes = []
    for line in manifest.read_text().splitlines():
        record = json.loads(line)
        if record["item"]["kind"] == "Genos Exome" and record["path"] and record["outcome"] != "duplicate":
            path = pgp_participants / record["path"].split("/participants/", 1)[1]
            if path.exists() and path.stat().st_size > 1_000_000:
                exomes.append(path)
    if not exomes:
        pytest.skip("no Genos exome upload on disk")
    checked = 0
    for path in exomes:
        if sniff(path).format != "vcf":
            continue
        assert loaded_facet(path)["callset_type"] == "wes"
        checked += 1
    assert checked > 0
