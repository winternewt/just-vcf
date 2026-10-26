"""What a VCF facet carries by variant class, the content filter it was cut by, and what that does to
absence (SCHEMA.md Decision 4; RM4's filtered-genome question)."""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import AbsenceReason, ContentFilter, ContentFilterBasis, VariantClass, absence_is_reference, absence_support, load
from just_vcf.callset import (
    MAX_FILTER_RESIDUE,
    MIN_FILTER_ALLELES,
    MIN_UNFILTERED_RESIDUE,
    MIN_WGS_SITES,
    READS_REFERENCE,
    infer_content_filter,
    named_content_filter,
    rows_content_filter,
    stated_content_filters,
    variant_class_counts,
)
from just_vcf.formats import get_spec

SNV, MNV, INDEL, SV = VariantClass.SNV, VariantClass.MNV, VariantClass.INDEL, VariantClass.SV
AUTOSOMES = [str(n) for n in range(1, 23)]


def counts(snv: int = 0, mnv: int = 0, indel: int = 0, sv: int = 0) -> dict[VariantClass, int]:
    return {SNV: snv, MNV: mnv, INDEL: indel, SV: sv}


# --- counting ---------------------------------------------------------------------------------


def test_each_alt_is_counted_in_its_class() -> None:
    rows = [
        ("A", ["G"]),  # snv
        ("A", ["G", "AT"]),  # snv + insertion
        ("AC", ["GT"]),  # mnv
        ("ACGT", ["A"]),  # deletion
        ("A", [""]),  # an unanchored deletion (F12's `-`, read as the empty allele)
        ("", ["T"]),  # an insertion point with no anchor base (CGI's native files)
        ("A" * 400, ["A"]),  # a 399 bp deletion written out base by base is still an indel
        ("N", ["<DEL>"]),  # symbolic
        ("G", ["G]17:198982]"]),  # breakend
        ("G", [".G"]),  # single breakend
        ("A", ["<*>"]),  # gVCF reference allele: no class
        ("A", ["G", "*"]),  # spanning deletion: no class
        ("A", ["<CGA_NOCALL>"]),  # no-call: no class
        ("A", []),  # monomorphic record (ALT .)
    ]
    frame = pl.LazyFrame(
        {"source_ref": [r for r, _ in rows], "alts": [a for _, a in rows]},
        schema={"source_ref": pl.String, "alts": pl.List(pl.String)},
    )
    assert variant_class_counts(frame) == counts(snv=3, mnv=1, indel=5, sv=3)


# --- what the rows show -----------------------------------------------------------------------


def test_rows_show_one_class_with_a_small_residue() -> None:
    kept = MIN_FILTER_ALLELES * 10
    residue = int(kept * MAX_FILTER_RESIDUE * 0.9)
    assert rows_content_filter(counts(snv=kept, indel=residue))[:2] == (ContentFilter.SNV_ONLY, True)
    assert rows_content_filter(counts(indel=kept, snv=residue))[:2] == (ContentFilter.INDEL_ONLY, True)
    assert rows_content_filter(counts(sv=kept))[:2] == (ContentFilter.SV_ONLY, True)


def test_rows_of_a_natural_mix_show_no_filter() -> None:
    # A whole genome: roughly one indel to seven SNVs.
    assert rows_content_filter(counts(snv=4_000_000, indel=600_000, mnv=0))[:2] == (None, True)


def test_rows_too_few_or_between_the_bands_settle_nothing() -> None:
    assert rows_content_filter(counts(snv=MIN_FILTER_ALLELES - 1))[:2] == (None, False)
    kept = MIN_FILTER_ALLELES * 100
    between = int(kept * (MAX_FILTER_RESIDUE + MIN_UNFILTERED_RESIDUE) / 2)
    assert rows_content_filter(counts(snv=kept, indel=between))[:2] == (None, False)
    # MNVs alone name no filter in the vocabulary.
    assert rows_content_filter(counts(mnv=kept))[:2] == (None, False)
    assert rows_content_filter(counts())[:2] == (None, False)


# --- what the header states -------------------------------------------------------------------

DATE = "; Date=Tue Oct 15 10:00:00 2019"


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (f"##bcftools_viewCommand=view -v snps -Oz -o out.vcf.gz in.vcf.gz{DATE}", ContentFilter.SNV_ONLY),
        (f"##bcftools_viewCommand=view -v indels in.vcf.gz{DATE}", ContentFilter.INDEL_ONLY),
        (f"##bcftools_viewCommand=view --types=snps in.vcf.gz{DATE}", ContentFilter.SNV_ONLY),
        (f"##bcftools_viewCommand=view -V indels,mnps,other in.vcf.gz{DATE}", ContentFilter.SNV_ONLY),
        (f"##bcftools_viewCommand=view -v snps,indels in.vcf.gz{DATE}", None),
        (f"##bcftools_viewCommand=view -V indels in.vcf.gz{DATE}", None),
        (f"##bcftools_viewCommand=view -i 'TYPE==\"indel\"' in.vcf.gz{DATE}", ContentFilter.INDEL_ONLY),
        (f"##bcftools_filterCommand=filter -i 'TYPE=\"snp\" ' -Oz in.vcf.gz{DATE}", ContentFilter.SNV_ONLY),
        (f"##bcftools_filterCommand=filter -e 'TYPE!=\"indel\"' in.vcf.gz{DATE}", ContentFilter.INDEL_ONLY),
        (f"##bcftools_filterCommand=filter -i 'TYPE=\"snp\" && QUAL>30' in.vcf.gz{DATE}", None),
        # A soft filter flags rows and removes none.
        (f"##bcftools_filterCommand=filter -s NotSNP -i 'TYPE=\"snp\"' in.vcf.gz{DATE}", None),
        (f"##bcftools_viewCommand=view -f PASS in.vcf.gz{DATE}", None),
        (
            '##GATKCommandLine.SelectVariants=<ID=SelectVariants,Version=3.4-46,CommandLineOptions="analysis_type='
            'SelectVariants select=false selectTypeToInclude=[SNP] selectTypeToExclude=[]">',
            ContentFilter.SNV_ONLY,
        ),
        (
            '##GATKCommandLine.SelectVariants.2=<ID=SelectVariants,CommandLineOptions="select=false '
            'selectTypeToInclude=[] selectTypeToExclude=[]">',
            None,
        ),
        ('##GATKCommandLine.SelectVariants=<ID=SelectVariants,CommandLineOptions="selectTypeToInclude=[SNP, MNP]">', None),
        (
            '##GATKCommandLine=<ID=SelectVariants,CommandLine="SelectVariants --select-type-to-include INDEL '
            '--variant in.vcf.gz",Version="4.2.0.0">',
            ContentFilter.INDEL_ONLY,
        ),
        ('##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --select-type-to-include SNP">', None),
        # DRAGEN's small-variant caller is named DRAGEN_SNV, and its unsplit output carries indels.
        ("##source=DRAGEN_SNV", None),
    ],
)
def test_stated_content_filters(line: str, expected: ContentFilter | None) -> None:
    found = stated_content_filters([line])
    assert [kept for kept, _ in found] == ([expected] if expected else [])


# --- what the file name says ------------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("sample.snp.vcf.gz", ContentFilter.SNV_ONLY),
        ("sample.filtered.snp.vcf.gz", ContentFilter.SNV_ONLY),
        ("sample.hard-filtered.snv.vcf.gz", ContentFilter.SNV_ONLY),
        ("genome.snps.vcf.gz", ContentFilter.SNV_ONLY),
        ("sample_indels.vcf", ContentFilter.INDEL_ONLY),
        ("sample.sv.vcf.gz", ContentFilter.SV_ONLY),
        ("sample.cnv.vcf.gz", ContentFilter.SV_ONLY),
        ("dir/SAMPLE.INDEL.VCF", ContentFilter.INDEL_ONLY),
        ("sample.hard-filtered.vcf.gz", None),
        ("sample.vcf.gz", None),
        ("snps_and_indels.vcf", None),
        ("sample.svm.vcf", None),
        ("synopsis.vcf", None),
    ],
)
def test_named_content_filter(name: str, expected: ContentFilter | None) -> None:
    found = named_content_filter(name)
    assert (found[0] if found else None) is expected


# --- weighing the evidence --------------------------------------------------------------------

GENOME_SNVS = counts(snv=4_000_000, indel=12)
GENOME_MIX = counts(snv=4_000_000, indel=600_000)
SMALL = counts(snv=500, indel=60)
VIEW_SNPS = [(ContentFilter.SNV_ONLY, "##bcftools_viewCommand")]
INDEL_NAME = (ContentFilter.INDEL_ONLY, ".indel.")
SNP_NAME = (ContentFilter.SNV_ONLY, ".snp.")


def test_rows_win_and_agreement_is_silent() -> None:
    inferred = infer_content_filter(GENOME_SNVS, VIEW_SNPS, SNP_NAME)
    assert (inferred.content_filter, inferred.basis, inferred.warnings) == (ContentFilter.SNV_ONLY, ContentFilterBasis.ROWS, ())


def test_rows_win_over_a_contradicting_header_and_name_with_a_warning_each() -> None:
    # The `.indel.` name and the header's `view -v snps` both lose to rows holding both classes.
    inferred = infer_content_filter(GENOME_MIX, VIEW_SNPS, INDEL_NAME)
    assert (inferred.content_filter, inferred.basis) == (None, None)
    assert len(inferred.warnings) == 2
    assert all("the rows are kept" in w for w in inferred.warnings)
    assert any(".indel." in w for w in inferred.warnings)


def test_the_header_settles_what_the_rows_cannot_and_beats_the_name() -> None:
    inferred = infer_content_filter(SMALL, VIEW_SNPS, INDEL_NAME)
    assert (inferred.content_filter, inferred.basis) == (ContentFilter.SNV_ONLY, ContentFilterBasis.STATED)
    assert len(inferred.warnings) == 1 and "the header is kept" in inferred.warnings[0]


def test_the_name_alone_is_the_basis_only_when_nothing_else_speaks() -> None:
    inferred = infer_content_filter(SMALL, [], SNP_NAME)
    assert (inferred.content_filter, inferred.basis, inferred.warnings) == (ContentFilter.SNV_ONLY, ContentFilterBasis.FILE_NAME, ())
    assert ".snp." in inferred.reason


def test_header_lines_that_disagree_settle_nothing() -> None:
    stated = [(ContentFilter.SNV_ONLY, "##bcftools_viewCommand"), (ContentFilter.INDEL_ONLY, "##GATKCommandLine.SelectVariants")]
    inferred = infer_content_filter(SMALL, stated, SNP_NAME)
    assert (inferred.content_filter, inferred.basis) == (None, None)
    assert len(inferred.warnings) == 1 and "settles none" in inferred.warnings[0]


def test_no_evidence_is_no_filter() -> None:
    inferred = infer_content_filter(SMALL, [], None)
    assert (inferred.content_filter, inferred.basis, inferred.warnings) == (None, None, ())


# --- what absence means per class -------------------------------------------------------------


def facet(absence: str, classes: list[tuple[str, int]] | None, content_filter: str | None = None) -> dict:
    listed = None if classes is None else [{"variant_class": c, "alleles": n} for c, n in classes]
    return {"absence_means": absence, "variant_classes": listed, "content_filter": content_filter}


WGS_FULL = facet("weak_reference", [("snv", 4_000_000), ("indel", 600_000)])
WGS_SNP_HALF = facet("weak_reference", [("snv", 4_000_000), ("indel", 12)], "snv_only")


@pytest.mark.parametrize(
    ("row", "asked", "force", "expected"),
    [
        (WGS_FULL, "snv", False, True),
        (WGS_FULL, VariantClass.INDEL, False, True),
        # The genome wrote no SVs and no MNVs: their absence says nothing.
        (WGS_FULL, "sv", False, None),
        (WGS_FULL, "mnv", False, None),
        # Cut to SNVs: held back without force, and the cut class stays unknown even with it,
        # though a residue of indels is present.
        (WGS_SNP_HALF, "snv", False, None),
        (WGS_SNP_HALF, "snv", True, True),
        (WGS_SNP_HALF, "indel", False, None),
        (WGS_SNP_HALF, "indel", True, None),
        # A producer that called every omitted position reference: the same gates by class and cut.
        (facet("reference", [("snv", 3_000_000), ("indel", 400_000)]), "indel", False, True),
        (facet("reference", [("snv", 3_000_000), ("indel", 400_000)]), "sv", False, None),
        (facet("reference", [("snv", 3_000_000), ("indel", 12)], "snv_only"), "snv", False, None),
        (facet("reference", [("snv", 3_000_000), ("indel", 12)], "snv_only"), "snv", True, True),
        # Blocks decide per position, which a facet-level answer cannot see.
        (facet("blocks", [("snv", 10)]), "snv", True, None),
        (facet("unknown", [("snv", 10)]), "snv", True, None),
        # Classes not counted: nothing to go on.
        (facet("weak_reference", None), "snv", True, None),
        ({"absence_means": "weak_reference"}, "snv", True, None),
    ],
)
def test_absence_is_reference_is_true_or_unknown(row: dict, asked: str, force: bool, expected: bool | None) -> None:
    assert absence_is_reference(row, asked, force=force) is expected


def test_absence_is_reference_refuses_an_unknown_class() -> None:
    with pytest.raises(ValueError):
        absence_is_reference(WGS_FULL, "snp")
    with pytest.raises(ValueError):
        absence_support(WGS_FULL, "snp")


def typed(absence: str, callset_type: str | None) -> dict:
    return {**facet(absence, [("snv", 10)]), "callset_type": callset_type}


@pytest.mark.parametrize(
    ("row", "asked", "force", "reason"),
    [
        (WGS_FULL, "snv", False, AbsenceReason.WHOLE_GENOME_SILENCE),
        (WGS_FULL, "sv", False, AbsenceReason.CLASS_NOT_CALLED),
        (WGS_SNP_HALF, "snv", False, AbsenceReason.CONTENT_FILTER_HELD),
        (WGS_SNP_HALF, "snv", True, AbsenceReason.CONTENT_FILTER_FORCED),
        (WGS_SNP_HALF, "indel", False, AbsenceReason.CLASS_CUT_BY_FILTER),
        (WGS_SNP_HALF, "indel", True, AbsenceReason.CLASS_CUT_BY_FILTER),
        # A pure cut wrote none of the class: the cut, not the count, is the reason.
        (facet("weak_reference", [("snv", 4_000_000)], "snv_only"), "indel", True, AbsenceReason.CLASS_CUT_BY_FILTER),
        (facet("weak_reference", [("snv", 4_000_000)], "snv_only"), "sv", False, AbsenceReason.CLASS_CUT_BY_FILTER),
        (facet("reference", [("snv", 3_000_000), ("indel", 400_000)]), "indel", False, AbsenceReason.OMITTED_CALLED_REFERENCE),
        (facet("reference", [("snv", 3_000_000), ("indel", 12)], "snv_only"), "snv", True, AbsenceReason.CONTENT_FILTER_FORCED),
        (facet("blocks", [("snv", 10)]), "snv", True, AbsenceReason.OUTSIDE_BLOCKS),
        (typed("unknown", "array"), "snv", True, AbsenceReason.ARRAY),
        (typed("unknown", "wes"), "snv", True, AbsenceReason.TARGETED),
        (typed("unknown", "panel"), "snv", True, AbsenceReason.TARGETED),
        (typed("unknown", "imputed"), "snv", True, AbsenceReason.IMPUTED),
        (typed("unknown", None), "snv", True, AbsenceReason.CALLSET_TYPE_UNKNOWN),
        (facet("weak_reference", None), "snv", True, AbsenceReason.CLASS_NOT_COUNTED),
        ({"absence_means": "weak_reference"}, "snv", True, AbsenceReason.CLASS_NOT_COUNTED),
    ],
)
def test_absence_support_gives_the_gate_that_decided(row: dict, asked: str, force: bool, reason: AbsenceReason) -> None:
    """Each gate of `absence_is_reference` answers with its own reason, and the answer is the same one."""
    support = absence_support(row, asked, force=force)
    assert support.reason is reason
    assert support.reads_reference is absence_is_reference(row, asked, force=force)
    assert support.reads_reference is (True if reason in READS_REFERENCE else None)
    assert bool(support) is (support.reads_reference is True)


# --- at load ----------------------------------------------------------------------------------

HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={chrom},length=200000000>" for chrom in AUTOSOMES],
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
]
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"


def genome_vcf(path: Path, indel_every: int | None) -> Path:
    """A whole-genome-shaped VCF of SNVs, with an insertion every `indel_every` rows (None: none)."""
    per = (MIN_WGS_SITES + 22) // len(AUTOSOMES)
    step = 100_000_000 // per
    rows = pl.concat(
        pl.DataFrame({"chrom": pl.repeat(chrom, per, eager=True), "pos": pl.int_range(0, per, eager=True) * step + 1})
        for chrom in AUTOSOMES
    ).with_row_index()
    alt = pl.lit("G") if indel_every is None else pl.when(pl.col("index") % indel_every == 0).then(pl.lit("AT")).otherwise(pl.lit("G"))
    path.write_text("\n".join([*HEADER, COLUMNS]) + "\n")
    with path.open("a") as handle:
        rows.select(
            "chrom",
            "pos",
            pl.lit(".").alias("id"),
            pl.lit("A").alias("ref"),
            alt.alias("alt"),
            pl.lit("50").alias("qual"),
            pl.lit("PASS").alias("filter"),
            pl.lit(".").alias("info"),
            pl.lit("GT:DP").alias("format"),
            pl.lit("0/1:30").alias("S1"),
        ).write_csv(handle, separator="\t", include_header=False)
    return path


def test_load_keeps_an_snv_only_genome_wgs_and_holds_its_absences_back(tmp_path: Path) -> None:
    # Named `.indel.` to show that the rows, not the name, decide, and that the conflict is reported.
    loaded = load(genome_vcf(tmp_path / "sample.indel.vcf", indel_every=None))
    row = loaded.facets.collect().row(0, named=True)
    assert (row["callset_type"], row["absence_means"]) == ("wgs", "weak_reference")
    assert row["variant_classes"] == [{"variant_class": "snv", "alleles": (MIN_WGS_SITES + 22) // 22 * 22}]
    assert (row["content_filter"], row["content_filter_basis"]) == ("snv_only", "rows")
    assert [w for w in loaded.warnings if "content filter" in w and ".indel." in w and "sample" not in w]
    assert absence_is_reference(row, "snv") is None
    assert absence_is_reference(row, "snv", force=True) is True
    assert absence_is_reference(row, "indel", force=True) is None


def test_load_reads_a_mixed_genome_as_unfiltered(tmp_path: Path) -> None:
    loaded = load(genome_vcf(tmp_path / "sample.vcf", indel_every=8))
    row = loaded.facets.collect().row(0, named=True)
    classes = {entry["variant_class"]: entry["alleles"] for entry in row["variant_classes"]}
    assert set(classes) == {"snv", "indel"} and sum(classes.values()) == (MIN_WGS_SITES + 22) // 22 * 22
    assert (row["content_filter"], row["content_filter_basis"]) == (None, None)
    assert absence_is_reference(row, "snv") is True and absence_is_reference(row, "indel") is True
    assert not [w for w in loaded.warnings if "content filter" in w]


def small_vcf(path: Path, meta: tuple[str, ...] = ()) -> Path:
    rows = [f"1\t{100 + i * 1000}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30" for i in range(60)]
    path.write_text("\n".join([*HEADER, *meta, COLUMNS, *rows]) + "\n")
    return path


def test_load_takes_a_stated_filter_and_then_the_name(tmp_path: Path) -> None:
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    stated = load(small_vcf(tmp_path / "a" / "x.vcf", (f"##bcftools_viewCommand=view -v snps in.vcf.gz{DATE}",))).facets.collect()
    assert (stated["content_filter"][0], stated["content_filter_basis"][0]) == ("snv_only", "stated")
    named = load(small_vcf(tmp_path / "b" / "x.snps.vcf")).facets.collect()
    assert (named["content_filter"][0], named["content_filter_basis"][0]) == ("snv_only", "file_name")


def test_a_vcf_of_blocks_only_counts_an_empty_list(tmp_path: Path) -> None:
    """Counted and found nothing (`[]`) is not the same as not counted (null)."""
    meta = ('##INFO=<ID=END,Number=1,Type=Integer,Description="End">', '##ALT=<ID=*,Description="Any other allele">')
    blocks = [f"1\t{100 + i * 1000}\t.\tA\t<*>\t.\t.\tEND={600 + i * 1000}\tGT:DP\t0/0:20" for i in range(5)]
    path = tmp_path / "x.g.vcf"
    path.write_text("\n".join([*HEADER, *meta, COLUMNS, *blocks]) + "\n")
    row = load(path).facets.collect().row(0, named=True)
    assert (row["variant_classes"], row["content_filter"], row["absence_means"]) == ([], None, "blocks")
    assert absence_is_reference(row, "snv", force=True) is None


def test_a_genotype_table_counts_no_classes(tmp_path: Path) -> None:
    """An array file writes no REF/ALT, so nothing is counted: null, never an empty list."""
    path = tmp_path / "x.snps.txt"
    path.write_bytes(get_spec("23andme").example.encode())
    row = load(path).facets.collect().row(0, named=True)
    assert row["absence_means"] == "unknown"
    assert (row["variant_classes"], row["content_filter"], row["content_filter_basis"]) == (None, None, None)
    assert absence_is_reference(row, "snv", force=True) is None
    assert absence_support(row, "snv", force=True).reason is AbsenceReason.ARRAY
