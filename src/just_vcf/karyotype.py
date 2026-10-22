"""Sex-chromosome evidence for one facet, and the karyotype estimate it supports (RM269).

`sex_chromosome_evidence(loaded)` measures what a facet's own calls say about its X and Y
chromosomes and turns that into an estimate with a reason. Nothing is stored: like `vrs_withheld()`,
the result is computed from the frames on request, so it never drifts from the rows it summarises and
its thresholds can move without a schema change. There is no person in the schema yet, so the result
is per facet; liftovers of one callset are one piece of evidence, not several, and weighing facets
against each other belongs to a later person-level layer.

What is measured, over the rows a caller did not filter out (`filters` is `["PASS"]` or null, which
is VCF's "filters not applied"; FORMAT `FT` likewise; `evidence_basis = assay`):

1. **GT ploidy as written** per group (autosomes, non-PAR X, non-PAR Y, PAR, MT): rows whose genotype
   holds one, two or more haplotypes. A one-allele genotype on X is the caller's (or vendor's) decision
   that the sample has one X. It is recorded as that assertion and never counted as evidence.
2. **Heterozygosity on calls written diploid**: the share of heterozygous genotypes among the
   non-PAR X calls that carry an ALT, and the median allele balance of those hets from AD (a
   contamination or mosaic signal; reported, not used in the estimate). The same for the
   X-degenerate Y, reported only: it does not separate samples consistently across callers.
3. **Depth from the calls' own FORMAT DP**: the median DP of ALT-carrying calls on non-PAR X and on
   the X-degenerate male-specific Y (`X_DEGENERATE_Y`), each divided by the same median on chr20.

The estimate:

- **X copies** (one, or two or more) from two independent lines, the het share and the X depth ratio.
  Where both are measured they must agree. The het share counts only where the caller wrote most X
  ALT calls diploid: a sex-aware caller writes a one-X sample's X haploid and keeps a few diploid
  calls it chose itself (all hets, in PGP-HMS huBE0518), which are no sample of anything. Then depth
  alone decides. A facet with neither line answers unknown.
- **Y present** from the X-degenerate Y depth ratio, and only ever `True`. Y calls, even in that
  region, occur in samples with and without a Y (mismapped reads), so their presence says nothing. A
  high depth there shows a Y; a low one does not show its absence: one man's GATK genome read 0.433,
  above the 0.400 the highest two-X genome read (PGP-HMS calibration, docs/ROADMAP.md RM269).
- **Karyotype** from the two: one X and Y present is `XY`, two X and Y present `XXY`. `XX` and `XO`
  need a Y shown absent, which call-site depth cannot do, so they wait for a coverage facet with
  XXX, XYY and the rarer classes (a Y ratio above the usual range still reads as one Y present). A
  value in an unknown band, lines that disagree, too few sites, a facet that is not a whole genome or a
  build without a region set all answer `UNKNOWN_KARYOTYPE` with the reason, never the nearest class.
  `x_copies` is still answered wherever the X lines settle it.

The estimate names a karyotype, not a sex: a person with an XX karyotype and a male phenotype has
`x_copies = two_or_more`, and would be `XX` once Y absence can be shown. It is read from variant
calls; it is not a cytogenetic result and not a diagnosis.

The X thresholds come from the lead's measurement of 2026-10-22 on 16 local single-sample whole-genome
VCFs (docs/ROADMAP.md, RM269): each line separated them into two clusters with a gap, and each unknown
band below is the central half of its gap, so every class edge keeps a quarter of the gap as margin.
The Y threshold was checked against 27 PGP-HMS VCFs with self-reported sex (same section).
"""

from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict, field_validator

from just_vcf.loading import Loaded
from just_vcf.schema import CallsetType

# --- Regions -----------------------------------------------------------------------------------
# Derived from GIAB genome stratifications v3.6, `XY` group (docs/PROVENANCE.md), converted from BED
# (0-based, half-open) to 1-based closed `pos` ranges; tests/test_karyotype.py re-derives them from the
# BED files when they are on disk. The derivation is just-vcf's own, not a GIAB product.

PAR_BOUNDS: dict[str, dict[str, tuple[tuple[int, int], tuple[int, int]]]] = {
    "GRCh37": {
        "X": ((60_002, 2_699_520), (154_931_045, 155_260_560)),
        "Y": ((10_001, 2_749_806), (59_133_471, 59_363_566)),
    },
    "GRCh38": {
        "X": ((10_001, 2_781_479), (155_701_384, 156_030_895)),
        "Y": ((10_002, 2_781_479), (56_887_903, 57_217_415)),
    },
}
"""PAR1 and PAR2 per build and chromosome, 1-based closed: from the first to the last base of GIAB's
`chr{X,Y}_PAR` intervals on each arm (GIAB leaves assembly gaps out, so its PAR file is fragmented)."""

X_DEGENERATE_Y: dict[str, tuple[tuple[int, int], ...]] = {
    "GRCh37": (
        (2_749_807, 2_918_085), (6_616_755, 7_442_522), (13_104_554, 13_143_954), (13_193_955, 13_748_578),
        (13_798_579, 16_096_353), (16_170_614, 17_986_973), (18_017_095, 18_271_675), (18_537_846, 19_568_146),
        (21_032_221, 22_369_679), (22_419_680, 23_467_839), (28_561_583, 28_819_361), (58_819_362, 58_917_656),
        (58_967_657, 59_133_470),
    ),
    "GRCh38": (
        (2_781_480, 3_050_044), (6_748_714, 7_574_481), (10_316_945, 10_544_039), (10_594_040, 10_633_440),
        (10_645_834, 10_649_989), (10_651_422, 10_669_737), (10_670_733, 10_674_058), (10_676_545, 10_679_715),
        (10_682_443, 10_691_573), (10_691_903, 10_694_192), (10_744_193, 10_747_280), (10_747_306, 10_809_671),
        (10_810_480, 10_816_758), (10_817_293, 10_852_900), (10_855_136, 10_871_218), (10_871_898, 10_890_419),
        (10_891_943, 10_896_525), (10_898_185, 10_908_519), (10_909_102, 10_922_486), (10_923_565, 10_956_868),
        (10_957_768, 10_961_625), (10_962_534, 10_965_694), (10_967_285, 10_969_992), (10_970_621, 10_986_652),
        (10_986_678, 11_002_276), (11_003_060, 11_012_528), (11_013_047, 11_016_992), (11_017_248, 11_023_374),
        (11_024_267, 11_028_342), (11_029_161, 11_036_252), (11_037_033, 11_592_902), (11_642_903, 11_647_442),
        (11_647_710, 11_653_064), (11_653_103, 11_660_374), (11_662_182, 11_663_986), (11_664_007, 11_669_948),
        (11_670_089, 11_671_800), (11_671_821, 11_673_795), (11_674_124, 13_984_473), (14_058_734, 15_875_093),
        (15_905_215, 16_159_795), (16_425_966, 17_456_266), (18_870_335, 20_207_793), (20_257_794, 21_305_953),
        (26_415_436, 26_673_214), (56_673_215, 56_771_509), (56_821_510, 56_887_902),
    ),
}
"""The X-degenerate male-specific Y, 1-based closed: GIAB's `chrY_nonPAR` minus `chrY_XTR` minus
`chrY_ampliconic`, kept between PAR1 and PAR2 (GRCh38's `chrY_nonPAR` also lists the single base before
PAR1, which is dropped). Sequence here has no close copy on X, so reads that land in it come from a Y."""

REGION_SET = "GIAB genome stratifications v3.6 XY, derived by just-vcf"

DEPTH_REFERENCE_CHROM = "20"
"""The autosome whose median DP is the denominator of both depth ratios. The thresholds were measured
against it, so they inherit it."""

# --- Thresholds (lead, 2026-10-22, 16 local WGS VCFs; the unknown band is the gap's central half) ---

ONE_X_MAX_HET_SHARE = 0.23
TWO_X_MIN_HET_SHARE = 0.46
"""Measured clusters 0.009–0.115 and 0.580–0.626; gap 0.115–0.580."""

ONE_X_MAX_DEPTH_RATIO = 0.60
TWO_X_MIN_DEPTH_RATIO = 0.79
"""Measured clusters 0.42–0.50 and 0.89–1.00; gap 0.50–0.89."""

TWO_X_MAX_DEPTH_RATIO = 1.25
"""Above it the depth fits more than two X better than two (three copies would sit near 1.5): unknown,
not two. The midpoint between two and three copies, from copy number alone."""

Y_PRESENT_MIN_DEPTH_RATIO = 0.52
"""Local genomes: 0.35–0.40 with two X, 0.56–0.82 with one. PGP-HMS: self-reported women 0.125 and
0.222, men 0.433–1.970. 0.52 keeps 0.12 above every two-X genome; the two male GATK files below it
answer unknown. There is no absent threshold: men reach the two-X range."""

MIN_X_SITES = 1_000
"""Diploid-written non-PAR X ALT calls needed for a het share, and X ALT calls with DP for a depth
ratio. A whole genome has tens of thousands; at 1,000 a share's standard error is under 0.016."""

MIN_REFERENCE_SITES = 1_000
"""chr20 ALT calls with DP needed for the denominator."""

MIN_Y_SITES = 300
"""X-degenerate Y ALT calls with DP needed for the Y ratio; the measured files had thousands each."""

ESTIMATE_BUILDS = frozenset(PAR_BOUNDS)

_AUTOSOMES = [str(n) for n in range(1, 23)]


class KaryotypicSex(str, Enum):
    """GA4GH Phenopackets v2 `KaryotypicSex`, values verbatim. Only XY, XXY and UNKNOWN_KARYOTYPE are
    estimated here so far; the rest are listed so the vocabulary is the standard's."""

    UNKNOWN_KARYOTYPE = "UNKNOWN_KARYOTYPE"
    XX = "XX"
    XY = "XY"
    XO = "XO"
    XXY = "XXY"
    XXX = "XXX"
    XXYY = "XXYY"
    XXXY = "XXXY"
    XXXX = "XXXX"
    XYY = "XYY"
    OTHER_KARYOTYPE = "OTHER_KARYOTYPE"


class XCopies(str, Enum):
    """How many X chromosomes the evidence supports. It cannot tell two from three."""

    ONE = "one"
    TWO_OR_MORE = "two_or_more"


UNKNOWN_REASONS: frozenset[str] = frozenset(
    {
        "no_build",
        "unsupported_build",
        "not_wgs",
        "too_few_sites",
        "intermediate",
        "lines_disagree",
        "caller_asserted_only",
        "no_depth",
        "y_unusable",
        "y_not_shown",
    }
)
"""Why an answer is unknown. An open vocabulary: a new kind of evidence brings new reasons.

- `no_build`: the facet has no `assembly_label`, so the PARs cannot be placed.
- `unsupported_build`: a label with no region set here (T2T-CHM13).
- `not_wgs`: `callset_type` is not `wgs`; the thresholds were measured on whole genomes.
- `too_few_sites`: a line had fewer sites than its minimum.
- `intermediate`: a value fell in the unknown band between two classes, or outside every band.
- `lines_disagree`: the het share and the X depth ratio point to different X counts.
- `caller_asserted_only`: the only X evidence is the caller's own haploid genotypes.
- `no_depth`: the calls carry no DP (a liftover that dropped FORMAT, for one), so no depth line can be
  measured; the karyotype then waits even where the het share gives X copies, since Y needs depth.
- `y_unusable`: X copies are known but the Y line does not settle (the karyotype only; `y_reason`
  says why).
- `y_not_shown`: the Y depth ratio is below the present band. It does not mean no Y: call-site depth
  cannot show a Y absent (`y_reason` only).
"""


class GtPloidyCounts(BaseModel):
    """Rows by the number of haplotypes their genotype was written with. This is what the source wrote,
    not a ploidy anyone measured."""

    model_config = ConfigDict(frozen=True)

    one: int
    two: int
    more: int


class Heterozygosity(BaseModel):
    """Genotypes carrying an ALT in one region. Only calls written diploid are evidence of copies."""

    model_config = ConfigDict(frozen=True)

    diploid_alt_sites: int
    """Fully called two-allele genotypes with at least one ALT (rows flagged by `split_signal` left out)."""
    het_sites: int
    het_share: float | None
    """`het_sites / diploid_alt_sites`; None with no site."""
    haploid_alt_sites: int
    """One-allele genotypes with an ALT: the caller's assertion of one copy, not evidence of it."""
    balance_sites: int
    """`0/k` hets with an AD covering the site's alleles and some depth on the two of them."""
    median_het_allele_balance: float | None
    """Median of `AD[k] / (AD[0] + AD[k])` over those hets; near 0.5 for two X, lower with
    contamination or a mosaic loss of one X. Reported only."""


class DepthRatios(BaseModel):
    """Median FORMAT DP of ALT-carrying calls per region, and the two ratios against chr20."""

    model_config = ConfigDict(frozen=True)

    reference_chrom: str = DEPTH_REFERENCE_CHROM
    reference_sites: int
    reference_median_dp: float | None
    x_sites: int
    x_median_dp: float | None
    y_sites: int
    y_median_dp: float | None
    x_ratio: float | None
    """Non-PAR X median over chr20 median; None when either side is under its minimum."""
    y_ratio: float | None
    """X-degenerate Y median over chr20 median; None when either side is under its minimum."""


class SexChromosomeEvidence(BaseModel):
    """What one facet's calls say about its sex chromosomes, and the estimate it supports.

    Truthy only when a karyotype was estimated. `x_copies` can be known when the karyotype is not:
    a facet whose Y line is unusable still says one X or two.
    """

    model_config = ConfigDict(frozen=True)

    facet_id: str
    assembly_label: str | None
    callset_type: CallsetType | None
    caller: str | None
    """`facets.caller`, carried for context."""
    ploidy_statements: list[str] | None
    """`facets.ploidy_statements`, verbatim: the header's own ploidy lines (DRAGEN's
    `##referenceSexKaryotype`, GATK's `--sample-ploidy`). With `x_heterozygosity.haploid_alt_sites`
    this is what the caller asserted; neither is evidence of the karyotype."""
    region_set: str | None
    """The source of the PAR bounds and the X-degenerate Y; None without a supported build."""
    rows_counted: int
    """Rows the measurements read: assay rows no filter failed."""
    ploidy: dict[str, GtPloidyCounts]
    """By group: `autosomes`, `x_nonpar`, `y_nonpar`, `par`, `mt` (the sex-chromosome groups only with
    a supported build)."""
    x_heterozygosity: Heterozygosity | None
    """Non-PAR X."""
    y_heterozygosity: Heterozygosity | None
    """The X-degenerate Y. Reported, not used in the estimate."""
    depth: DepthRatios | None
    x_copies: XCopies | None
    x_copies_reason: str | None
    y_present: bool | None
    """`True` when the Y depth shows a Y; otherwise None, never `False` (module docstring)."""
    y_reason: str | None
    inferred_karyotypic_sex: KaryotypicSex
    """A Phenopackets `KaryotypicSex` value, estimated from calls: not a cytogenetic result."""
    reason: str | None
    """Why `inferred_karyotypic_sex` is unknown; None when it is not."""
    detail: str
    """The same in words, with the numbers."""

    @field_validator("x_copies_reason", "y_reason", "reason")
    @classmethod
    def _known_reason(cls, value: str | None) -> str | None:
        if value is not None and value not in UNKNOWN_REASONS:
            raise ValueError(f"unknown reason {value!r}")
        return value

    def __bool__(self) -> bool:
        return self.inferred_karyotypic_sex is not KaryotypicSex.UNKNOWN_KARYOTYPE


def sex_chromosome_evidence(loaded: Loaded) -> tuple[SexChromosomeEvidence, ...]:
    """Measure each loaded facet's sex-chromosome evidence and estimate its karyotype (module docstring).

    `loaded` is what `just_vcf.load` (or `rejoin_split`) returns. The result has one entry per row of
    `loaded.facets`, in that order: a freeBayes+VarScan2 merge loads as two facets and gets two
    answers, never one blended from both callers. The calls are read lazily, a few columns at a time.

    The answer is a karyotype, not a sex: an XX male reads two X copies (and `XX` once Y absence can
    be shown), never `XY` for his phenotype. It is an estimate from variant calls,
    never a cytogenetic or clinical result.
    """
    facets = loaded.facets.select("facet_id", "assembly_label", "callset_type", "caller", "ploidy_statements").collect()
    return tuple(
        estimate(
            measure(loaded.calls.filter(pl.col("facet_id") == facet["facet_id"]), facet["assembly_label"]),
            facet["facet_id"],
            facet["assembly_label"],
            CallsetType(facet["callset_type"]) if facet["callset_type"] is not None else None,
            facet["caller"],
            facet["ploidy_statements"],
        )
        for facet in facets.iter_rows(named=True)
    )


def measure(calls: pl.LazyFrame, assembly_label: str | None) -> dict:
    """The raw measurements of a `calls` frame, per group, as plain numbers (one small collect)."""
    supported = assembly_label in ESTIMATE_BUILDS
    counted = (
        calls.filter(
            (pl.col("evidence_basis") == "assay")
            & _not_failed("filters")
            & _not_failed("sample_filters")
        )
        .select("chrom", "pos", "gt", "gt_alleles", "dp", "ad", "alts", "source_alleles", "split_signal")
        .with_columns(_group(assembly_label).alias("group"))
        .filter(pl.col("group").is_not_null())
    )
    n_haplotypes = pl.col("gt_alleles").list.len()
    gt0, gt1 = pl.col("gt").list.get(0, null_on_oob=True), pl.col("gt").list.get(1, null_on_oob=True)
    has_alt = pl.col("gt").list.eval(pl.element() >= 1).list.any().fill_null(False)
    diploid_called = (pl.col("gt").list.len() == 2) & gt0.is_not_null() & gt1.is_not_null()
    diploid_alt = diploid_called & has_alt & pl.col("split_signal").is_null()
    het = diploid_alt & (gt0 != gt1)
    k = pl.max_horizontal(gt0, gt1)
    ad_ok = (
        pl.col("ad").is_not_null()
        & (pl.col("ad").list.len() == pl.col("alts").list.len() + 1)
        & (pl.col("source_alleles").list.len() == pl.col("alts").list.len() + 1)
    )
    ad_ref = pl.col("ad").list.get(0, null_on_oob=True).cast(pl.Float64)
    ad_alt = pl.col("ad").list.get(k.cast(pl.Int64), null_on_oob=True).cast(pl.Float64)
    balance_row = het & (pl.min_horizontal(gt0, gt1) == 0) & ad_ok & ((ad_ref + ad_alt) > 0)
    depth_row = has_alt & pl.col("dp").is_not_null()
    dp = pl.col("dp").cast(pl.Float64)
    xdeg = _in_intervals(X_DEGENERATE_Y.get(assembly_label, ())) if supported else pl.lit(False)
    reference = pl.col("chrom") == DEPTH_REFERENCE_CHROM

    per_group = (
        counted.group_by("group")
        .agg(
            pl.len().alias("rows"),
            (n_haplotypes == 1).sum().alias("one"),
            (n_haplotypes == 2).sum().alias("two"),
            (n_haplotypes > 2).sum().alias("more"),
            diploid_alt.sum().alias("diploid_alt"),
            het.sum().alias("het"),
            ((n_haplotypes == 1) & has_alt).sum().alias("haploid_alt"),
            balance_row.sum().alias("balance_sites"),
            (ad_alt / (ad_ref + ad_alt)).filter(balance_row).median().alias("balance"),
            depth_row.sum().alias("depth_sites"),
            dp.filter(depth_row).median().alias("median_dp"),
            (depth_row & reference).sum().alias("reference_sites"),
            dp.filter(depth_row & reference).median().alias("reference_median_dp"),
            (depth_row & xdeg).sum().alias("xdeg_sites"),
            dp.filter(depth_row & xdeg).median().alias("xdeg_median_dp"),
            (diploid_alt & xdeg).sum().alias("xdeg_diploid_alt"),
            (het & xdeg).sum().alias("xdeg_het"),
            ((n_haplotypes == 1) & has_alt & xdeg).sum().alias("xdeg_haploid_alt"),
            (balance_row & xdeg).sum().alias("xdeg_balance_sites"),
            (ad_alt / (ad_ref + ad_alt)).filter(balance_row & xdeg).median().alias("xdeg_balance"),
        )
        .collect()
    )
    return {"supported": supported, "groups": {row.pop("group"): row for row in per_group.to_dicts()}}


def _not_failed(column: str) -> pl.Expr:
    """No filter failed: the column is null (VCF `.`, filters not applied) or exactly `["PASS"]`."""
    values = pl.col(column)
    return values.is_null() | ((values.list.len() == 1) & (values.list.first() == "PASS"))


def _group(assembly_label: str | None) -> pl.Expr:
    """Each row's group; null for contigs outside the primary assembly. Without a supported build the
    sex chromosomes have no group, since their PARs cannot be placed."""
    chrom = pl.col("chrom")
    expr = pl.when(chrom.is_in(_AUTOSOMES)).then(pl.lit("autosomes")).when(chrom == "MT").then(pl.lit("mt"))
    if assembly_label in ESTIMATE_BUILDS:
        bounds = PAR_BOUNDS[assembly_label]
        expr = (
            expr.when(chrom == "XY")
            .then(pl.lit("par"))
            .when((chrom == "X") & _in_intervals(bounds["X"]))
            .then(pl.lit("par"))
            .when(chrom == "X")
            .then(pl.lit("x_nonpar"))
            .when((chrom == "Y") & _in_intervals(bounds["Y"]))
            .then(pl.lit("par"))
            .when(chrom == "Y")
            .then(pl.lit("y_nonpar"))
        )
    return expr.otherwise(pl.lit(None, dtype=pl.String))


def _in_intervals(intervals: tuple[tuple[int, int], ...]) -> pl.Expr:
    """`pos` inside any of the 1-based closed intervals."""
    return pl.any_horizontal([pl.col("pos").is_between(start, end) for start, end in intervals]) if intervals else pl.lit(False)


def estimate(
    measured: dict,
    facet_id: str,
    assembly_label: str | None,
    callset_type: CallsetType | None,
    caller: str | None = None,
    ploidy_statements: list[str] | None = None,
) -> SexChromosomeEvidence:
    """Turn `measure`'s numbers into the evidence and the estimate (rules in the module docstring).

    `caller` and `ploidy_statements` are the facet's own columns, carried through unread: what the
    caller said about ploidy is its assertion, kept beside the evidence and never counted in it."""
    groups = measured["groups"]
    empty = {"rows": 0, "one": 0, "two": 0, "more": 0}
    names = ["autosomes", "x_nonpar", "y_nonpar", "par", "mt"] if measured["supported"] else ["autosomes", "mt"]
    ploidy = {name: GtPloidyCounts(**{key: int(groups.get(name, empty)[key]) for key in ("one", "two", "more")}) for name in names}
    base = {
        "facet_id": facet_id,
        "assembly_label": assembly_label,
        "callset_type": callset_type,
        "caller": caller,
        "ploidy_statements": ploidy_statements,
        "rows_counted": int(sum(row["rows"] for row in groups.values())),
        "ploidy": ploidy,
    }

    def unknown(reason: str, detail: str, **known) -> SexChromosomeEvidence:
        fields = {"region_set": None, "x_heterozygosity": None, "y_heterozygosity": None, "depth": None, "x_copies": None, "x_copies_reason": reason, "y_present": None, "y_reason": reason}
        return SexChromosomeEvidence(**base, **{**fields, **known}, inferred_karyotypic_sex=KaryotypicSex.UNKNOWN_KARYOTYPE, reason=reason, detail=detail)

    if assembly_label is None:
        return unknown("no_build", "the facet has no assembly label, so the pseudo-autosomal regions cannot be placed")
    if not measured["supported"]:
        return unknown("unsupported_build", f"no PAR bounds or Y region set for {assembly_label} (only {', '.join(sorted(ESTIMATE_BUILDS))})")

    x, y, auto = groups.get("x_nonpar", {}), groups.get("y_nonpar", {}), groups.get("autosomes", {})
    x_het = _heterozygosity(x, "")
    y_het = _heterozygosity(y, "xdeg_")
    ref_sites, ref_dp = int(auto.get("reference_sites", 0)), auto.get("reference_median_dp")
    x_sites, x_dp = int(x.get("depth_sites", 0)), x.get("median_dp")
    y_sites, y_dp = int(y.get("xdeg_sites", 0)), y.get("xdeg_median_dp")
    reference_ok = ref_sites >= MIN_REFERENCE_SITES and bool(ref_dp)
    depth = DepthRatios(
        reference_sites=ref_sites,
        reference_median_dp=ref_dp,
        x_sites=x_sites,
        x_median_dp=x_dp,
        y_sites=y_sites,
        y_median_dp=y_dp,
        x_ratio=x_dp / ref_dp if reference_ok and x_sites >= MIN_X_SITES and x_dp is not None else None,
        y_ratio=y_dp / ref_dp if reference_ok and y_sites >= MIN_Y_SITES and y_dp is not None else None,
    )
    known = {"region_set": REGION_SET, "x_heterozygosity": x_het, "y_heterozygosity": y_het, "depth": depth}
    if callset_type is not CallsetType.WGS:
        kind = callset_type.value if callset_type else "unknown"
        return unknown("not_wgs", f"callset type {kind}; the thresholds were measured on whole genomes", **known)

    no_dp = int(auto.get("depth_sites", 0)) + x_sites + int(y.get("depth_sites", 0)) == 0
    x_copies, x_reason, x_detail = _x_copies(x_het, depth, no_dp)
    y_present, y_reason, y_detail = _y_present(depth, no_dp)
    detail = f"X: {x_detail}; Y: {y_detail}"
    lines = {"x_copies": x_copies, "x_copies_reason": x_reason, "y_present": y_present, "y_reason": y_reason}
    if x_copies is None:
        return unknown(x_reason, detail, **known, **{**lines})
    if y_present is None:
        return unknown("no_depth" if y_reason == "no_depth" else "y_unusable", detail, **known, **lines)
    karyotype = {XCopies.ONE: KaryotypicSex.XY, XCopies.TWO_OR_MORE: KaryotypicSex.XXY}[x_copies]
    return SexChromosomeEvidence(**base, **known, **lines, inferred_karyotypic_sex=karyotype, reason=None, detail=detail)


def _heterozygosity(group: dict, prefix: str) -> Heterozygosity:
    diploid, het = int(group.get(f"{prefix}diploid_alt", 0)), int(group.get(f"{prefix}het", 0))
    return Heterozygosity(
        diploid_alt_sites=diploid,
        het_sites=het,
        het_share=het / diploid if diploid else None,
        haploid_alt_sites=int(group.get(f"{prefix}haploid_alt", 0)),
        balance_sites=int(group.get(f"{prefix}balance_sites", 0)),
        median_het_allele_balance=group.get(f"{prefix}balance"),
    )


def _x_copies(x_het: Heterozygosity, depth: DepthRatios, no_dp: bool) -> tuple[XCopies | None, str | None, str]:
    """X copies from the het share and the X depth ratio; both must agree where both are measured."""
    written_diploid = x_het.diploid_alt_sites >= MIN_X_SITES and x_het.diploid_alt_sites > x_het.haploid_alt_sites
    share = x_het.het_share if written_diploid else None
    by_share = _band(share, ONE_X_MAX_HET_SHARE, TWO_X_MIN_HET_SHARE, None)
    by_depth = _band(depth.x_ratio, ONE_X_MAX_DEPTH_RATIO, TWO_X_MIN_DEPTH_RATIO, TWO_X_MAX_DEPTH_RATIO)
    said = (
        f"het share {_fmt(x_het.het_share)} over {x_het.diploid_alt_sites:,} diploid ALT call(s), "
        f"{x_het.haploid_alt_sites:,} haploid; depth ratio {_fmt(depth.x_ratio)} over {depth.x_sites:,} call(s) with DP"
    )
    if by_share == "intermediate" or by_depth == "intermediate":
        return None, "intermediate", said
    if by_share and by_depth and by_share != by_depth:
        return None, "lines_disagree", said
    found = by_share or by_depth
    if found:
        return XCopies(found), None, said
    if x_het.haploid_alt_sites >= MIN_X_SITES:
        return None, "caller_asserted_only", said
    return None, ("no_depth" if no_dp else "too_few_sites"), said


def _y_present(depth: DepthRatios, no_dp: bool) -> tuple[bool | None, str | None, str]:
    """Y present from the X-degenerate Y depth ratio: `True` or unknown, never `False`."""
    said = f"X-degenerate depth ratio {_fmt(depth.y_ratio)} over {depth.y_sites:,} call(s) with DP"
    if no_dp:
        return None, "no_depth", said
    if depth.y_ratio is None:
        return None, "too_few_sites", said
    if depth.y_ratio >= Y_PRESENT_MIN_DEPTH_RATIO:
        return True, None, said
    return None, "y_not_shown", said


def _band(value: float | None, one_max: float, two_min: float, two_max: float | None) -> str | None:
    """`one`, `two_or_more`, `intermediate`, or None when there is no value."""
    if value is None:
        return None
    if value <= one_max:
        return XCopies.ONE.value
    if value >= two_min and (two_max is None or value <= two_max):
        return XCopies.TWO_OR_MORE.value
    return "intermediate"


def _fmt(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.3f}"
