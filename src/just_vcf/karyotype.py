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
- **Y present** from the X-degenerate Y depth ratio. From the calls' DP it is only ever `True`. Y
  calls, even in that region, occur in samples with and without a Y (mismapped reads), so their
  presence says nothing. A high depth there shows a Y; a low one does not show its absence: one man's
  GATK genome read 0.433, above the 0.400 the highest two-X genome read (PGP-HMS calibration,
  docs/ROADMAP.md RM269). Depth at called sites is depth where a caller found something to call.
- **A coverage facet** (`just_vcf.coverage`, RM278), when one is given, adds a third X line and a
  second Y line: the median depth of the reads themselves over non-PAR X and over the X-degenerate Y,
  each over the median on chr20, read at every base and not only at called sites. There a Y can be
  shown absent (`y_present = False`): with no Y, reads reach that region only by mismapping. Where the
  two Y lines disagree (the calls show a Y, the coverage shows none) the answer is unknown.
- **Counts past two** come from the coverage line only. Copy number puts each copy at half an
  autosome's depth, so a coverage ratio counts X copies (one 0.5, two 1.0, three 1.5) and Y copies
  (none 0, one 0.5, two 1.0), each in a band of 0.125 either side: the central half of every gap
  between two counts is unknown, as for the other lines. The call-based lines cannot count: the het
  share says one X or two or more, the calls' Y depth says a Y is present or nothing. The calls' X
  depth says one, two, or (above 1.25, the midpoint between two and three copies) more than two,
  which a coverage three confirms and nothing else turns into a count. Each line allows a range of
  counts and the answer is the range every line allows (`x_copy_count`, `y_copy_count` where it is
  one number). No XXX or XYY genome was measured: these bands are copy number alone.
- **Karyotype** from the two counts: `XO`, `XY`, `XYY`; `XX`, `XXY`, `XXYY`; `XXX`, `XXXY`; three X
  with two Y is `OTHER_KARYOTYPE` (Phenopackets lists no XXXYY). Where a line could not count, the
  reading from before counts were measured stands: a het share's two or more reads two, a Y shown
  only by the calls reads one (so `XY` and `XXY` from calls alone, which cannot tell XY from XYY);
  more than two X with no coverage line to count them is unknown (`copies_uncounted`). A value in an
  unknown band, lines that disagree, too few sites, a facet that is not a whole genome, a build
  without a region set, or one whose thresholds were never measured (T2T-CHM13, measured and
  reported only) all answer `UNKNOWN_KARYOTYPE` with the reason, never the nearest class.
  `x_copies` is still answered wherever the X lines settle it.
- **An array facet** (`callset_type = array`) has no depth; its lines are read from the letters as
  written (`array`, `ArrayLines`), so no REF is needed. X: the het share of two-letter base calls on
  non-PAR X over the autosomes' het share, one X near 0 (genotyping errors only), two X near 3/4 (X
  diversity under the neutral expectation), the central half of the gap unknown. An X written with one
  letter is the vendor's decision (23andMe does so for a sample it calls male), counted and never read
  (`caller_asserted_only`). Y: the called share of the X-degenerate Y probes over the autosomes'
  called share, near 1 with a Y and near 0 without. A vendor may no-call every Y probe because it
  decided the sample is female, and nothing in the file tells that mask from a chip that saw no Y:
  on PGP-HMS, 23andMe exports with a two-letter X read exactly 0, AncestryDNA women 0.008-0.098
  (docs/ROADMAP.md RM269). So the Y line is read only with `array_y=True`; otherwise `y_reason` is
  `vendor_asserted` and the karyotype unknown. Arrays cannot count copies.

**The build** places the PARs and the X-degenerate Y, and it is the facet's `assembly_label`, whatever
in the library filled it: stated by the source at load, inferred at load from the header's `##contig`
lengths (RM4), inferred by `resolve_ref` from an array's own rows scored against a reference (F9), or
a liftover's target collection. The result says which (`assembly_label_basis`, the facet's column,
and `assembly_label_source`, the finer source read off the facet's other columns by SCHEMA.md
Decision 5's rule), and carries `ref_seqcol_digest`, the collection REF was read from, where it was.
Nothing here infers a build or falls back to one: an array table that states no build answers
`no_build` as loaded, and is estimated on the `resolve_ref` output, which carries the build the rows
gave it.

The estimate names a karyotype, not a sex: a person with an XX karyotype and a male phenotype has
`x_copies = two_or_more`, and `XX` where a coverage facet shows no Y. It is read from variant calls
and read depth; it is not a cytogenetic result and not a diagnosis.

The X thresholds come from the lead's measurement of 2026-10-22 on 16 local single-sample whole-genome
VCFs (docs/ROADMAP.md, RM269): each line separated them into two clusters with a gap, and each unknown
band below is the central half of its gap, so every class edge keeps a quarter of the gap as margin.
The Y threshold was checked against 27 PGP-HMS VCFs with self-reported sex (same section).
"""

from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict, field_validator

from just_vcf.coverage import Coverage
from just_vcf.links import LIFTOVER_TO
from just_vcf.loading import Loaded
from just_vcf.schema import AssemblyLabelBasis, CallsetType

# --- Regions -----------------------------------------------------------------------------------
# Derived from GIAB genome stratifications v3.6, `XY` group, and for T2T-CHM13's X-degenerate Y from the
# T2T consortium's sequence classes (docs/PROVENANCE.md), converted from BED (0-based, half-open) to
# 1-based closed `pos` ranges; tests/test_karyotype.py re-derives them from the BED files carried in
# assets/. The derivation is just-vcf's own, not a GIAB or T2T product.

PAR_BOUNDS: dict[str, dict[str, tuple[tuple[int, int], tuple[int, int]]]] = {
    "GRCh37": {
        "X": ((60_002, 2_699_520), (154_931_045, 155_260_560)),
        "Y": ((10_001, 2_749_806), (59_133_471, 59_363_566)),
    },
    "GRCh38": {
        "X": ((10_001, 2_781_479), (155_701_384, 156_030_895)),
        "Y": ((10_002, 2_781_479), (56_887_903, 57_217_415)),
    },
    "T2T-CHM13": {
        "X": ((1, 2_394_410), (153_925_835, 154_259_566)),
        "Y": ((1, 2_458_320), (62_122_810, 62_460_029)),
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
    "T2T-CHM13": (
        (2_458_321, 2_727_072), (6_400_876, 7_234_999), (12_656_464, 14_891_106), (14_964_928, 16_781_335),
        (16_811_394, 17_066_094), (17_332_347, 18_362_331), (19_776_672, 20_961_203), (21_226_264, 22_270_252),
    ),
}
"""The X-degenerate male-specific Y, 1-based closed. Sequence here has no close copy on X, so reads that
land in it come from a Y.

GRCh37 and GRCh38: GIAB's `chrY_nonPAR` minus `chrY_XTR` minus `chrY_ampliconic`, kept between PAR1
and PAR2 (GRCh38's `chrY_nonPAR` also lists the single base before PAR1, which is dropped). On those
builds the rest of the male-specific Y (heterochromatin, centromere) is an assembly gap, which GIAB's
nonPAR leaves out.

T2T-CHM13: the T2T consortium's `X-DEG` sequence class (`chm13v2.0_chrXY_sequence_class_v1.bed`, on
the complete Y of HG002), less GIAB's `chrY_XTR` and `chrY_PAR`. The class file's second X-DEG interval
starts two bases inside XTR2, and those two bases are dropped. GIAB's CHM13 set has no ampliconic file,
and on a gapless Y its nonPAR minus XTR would keep the ampliconic sequence, the heterochromatin, the
centromere and the satellites, so the GRCh37/38 derivation does not carry over."""

REGION_SET: dict[str, str] = {
    "GRCh37": "GIAB genome stratifications v3.6 XY, derived by just-vcf",
    "GRCh38": "GIAB genome stratifications v3.6 XY, derived by just-vcf",
    "T2T-CHM13": (
        "GIAB genome stratifications v3.6 XY (PAR, XTR) and T2T chm13v2.0 chrXY sequence classes v1 "
        "(X-degenerate Y), derived by just-vcf"
    ),
}
"""Where each build's PAR bounds and X-degenerate Y come from, as `region_set` reports it."""

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
"""Above it the depth fits more than two X better than two (three copies would sit near 1.5): the
line then allows three or more, never two, and names no count itself. The midpoint between two and
three copies, from copy number alone."""

Y_PRESENT_MIN_DEPTH_RATIO = 0.52
"""Local genomes: 0.35–0.40 with two X, 0.56–0.82 with one. PGP-HMS: self-reported women 0.125 and
0.222, men 0.433–1.970. 0.52 keeps 0.12 above every two-X genome; the two male GATK files below it
answer unknown. There is no absent threshold: men reach the two-X range."""

# --- Coverage thresholds (RM278): ratios of read-depth medians, not of call-site DP ---------------

COVERAGE_RATIO_PER_COPY = 0.5
"""An autosome holds two copies, so each copy of X or Y reads half its depth: one copy 0.5, two 1.0,
three 1.5, none 0 (plus mismapped reads)."""

COVERAGE_BAND_HALF_WIDTH = 0.125
"""A count's band reaches this far either side of its copy-number value: a quarter of the 0.5 between
two counts, so the unknown band between them is the gap's central half, as for the call-based lines.
A mosaic (a loss of Y or of one X in part of the cells) sits between two counts and reads unknown
there, not as either. On four BAMs every ratio lay within 0.034 of its copy-number value
(docs/ROADMAP.md RM269)."""

COVERAGE_X_COUNTS = (1, 2, 3)
COVERAGE_Y_COUNTS = (0, 1, 2)
"""The counts a coverage line reads. Outside their bands (no X, four X, three Y) the line is unknown,
not the nearest count: none of those was asked for or measured."""

COVERAGE_MIN_EVEN_SHARE = 0.75
"""Share of chr20's described bases (gaps left out) whose interval mean lies between half and twice
the chr20 median. A whole genome sits near 1; an exome or panel, deep on target and shallow off it,
sits far under, and its medians say nothing about copies. The 0.75 echoes RM4's whole-genome breadth
rule; the measured BAMs are in docs/ROADMAP.md RM278."""

MIN_COVERAGE_BASES = 1_000_000
"""Bases a coverage line needs in each region (chr20, non-PAR X, and the X-degenerate Y, which holds
about 10 Mb on each build), counted over intervals described at every base."""

MIN_X_SITES = 1_000
"""Diploid-written non-PAR X ALT calls needed for a het share, and X ALT calls with DP for a depth
ratio. A whole genome has tens of thousands; at 1,000 a share's standard error is under 0.016."""

MIN_REFERENCE_SITES = 1_000
"""chr20 ALT calls with DP needed for the denominator."""

MIN_Y_SITES = 300
"""X-degenerate Y ALT calls with DP needed for the Y ratio; the measured files had thousands each."""

REGION_BUILDS = frozenset(PAR_BOUNDS)
"""Builds whose sex chromosomes can be grouped and measured."""

ESTIMATE_BUILDS = frozenset({"GRCh37", "GRCh38"})
"""Builds the thresholds were measured on, and so the only ones estimated. T2T-CHM13 is measured and
reported, but its estimate waits for T2T callsets of known karyotype: its X-degenerate Y is another
region set, and no T2T genome was among the calibration files. Rejected: applying the GRCh37/38 bands
to it unmeasured."""

_AUTOSOMES = [str(n) for n in range(1, 23)]
_BASES = ["A", "C", "G", "T"]


class KaryotypicSex(str, Enum):
    """GA4GH Phenopackets v2 `KaryotypicSex`, values verbatim. Every value but XXXX is estimated here;
    all but XY and XXY only with a coverage facet, and OTHER_KARYOTYPE for three X with two Y."""

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
    """How many X chromosomes the evidence supports, one or more than one. The count itself, where a
    coverage line settles it, is `SexChromosomeEvidence.x_copy_count`."""

    ONE = "one"
    TWO_OR_MORE = "two_or_more"


UNKNOWN_REASONS: frozenset[str] = frozenset(
    {
        "no_build",
        "unsupported_build",
        "uncalibrated_build",
        "not_wgs",
        "too_few_sites",
        "intermediate",
        "lines_disagree",
        "caller_asserted_only",
        "no_depth",
        "y_unusable",
        "y_not_shown",
        "coverage_no_build",
        "coverage_build_mismatch",
        "coverage_uneven",
        "copies_uncounted",
        "vendor_asserted",
        "no_calls_not_written",
    }
)
"""Why an answer is unknown. An open vocabulary: a new kind of evidence brings new reasons.

- `no_build`: the facet has no `assembly_label`, so the PARs cannot be placed. An array table that
  states none gets one from `resolve_ref`, which infers it from the rows (F9); pass its output.
- `unsupported_build`: a label with no region set here (GRCh36).
- `uncalibrated_build`: a build with a region set but no measured thresholds (T2T-CHM13): the
  measurements are reported, the estimate is withheld.
- `not_wgs`: `callset_type` is neither `wgs` nor `array` (an array has its own lines); the depth
  thresholds were measured on whole genomes.
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
- `coverage_no_build`, `coverage_build_mismatch`: a coverage facet was given, but it has no assembly
  label, or another one than the facet's, so its regions cannot be placed; it is set aside
  (`coverage.reason` only).
- `coverage_uneven`: the coverage facet's chr20 depth is not even enough for a whole genome
  (`COVERAGE_MIN_EVEN_SHARE`): an exome's or a panel's reads; it is set aside (`coverage.reason` only).
- `vendor_asserted`: an array's Y probes were measured but not read: whether they are called follows
  the vendor's own sex call where it masks Y for the samples it calls female (`array_y`; `y_reason`
  only).
- `no_calls_not_written`: an array export with no no-call row on the autosomes (a genuine probe
  export calls 97-99.9% of them): it drops its no-calls, so its Y rows are only the called ones and
  their share says nothing (`y_reason` only).
- `copies_uncounted`: the X lines agree on more than two copies, but none counts them (the calls'
  X depth above two copies' band, with no coverage line): `x_copies` is `two_or_more`, the karyotype
  unknown.
"""


ASSEMBLY_LABEL_SOURCES: frozenset[str] = frozenset({"stated", "contig_lengths", "rows", "liftover"})
"""Where a facet's `assembly_label` came from (`assembly_label_source`). An open vocabulary: a new way
of filling the label brings a new value.

- `stated`: the source stated it, read at load (`assembly_label_basis = stated`).
- `contig_lengths`: inferred at load from the header's `##contig` lengths (RM4).
- `rows`: inferred by `resolve_ref` from an array's own rows, scored against each reference collection
  (F9); `ref_seqcol_digest` names the collection that fitted.
- `liftover`: a lifted facet's target collection (`relation = liftover_to`).
"""


def assembly_label_source(facet: dict) -> str | None:
    """Which of `ASSEMBLY_LABEL_SOURCES` filled the facet's `assembly_label`, from its own columns.

    SCHEMA.md Decision 5's rule: `inferred` is a liftover's where `relation = liftover_to`, else the
    contig lengths' where `header_contigs` gives a length (only a VCF read at load has them), else the
    rows' where REF was read (`ref_seqcol_digest`) for a facet that is no VCF, the only kind
    `resolve_ref` infers a build for. None without a label, and for a label no rule accounts for (one a
    caller wrote into the frame): the label is still used, its source is unknown.
    """
    if facet["assembly_label"] is None:
        return None
    if facet["assembly_label_basis"] == AssemblyLabelBasis.STATED.value:
        return "stated"
    if facet["assembly_label_basis"] != AssemblyLabelBasis.INFERRED.value:
        return None
    if facet["relation"] == LIFTOVER_TO:
        return "liftover"
    if any(entry["length"] is not None for entry in facet["header_contigs"] or []):
        return "contig_lengths"
    if facet["ref_seqcol_digest"] is not None and facet["vcf_version"] is None:
        return "rows"
    return None


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


class CoverageRatios(BaseModel):
    """Median read depth per region from a coverage facet, and the two ratios against chr20.

    The medians weight each interval by its length and take only intervals that lie wholly inside the
    region, carry a mean and were described at every base. Zero-depth intervals are left out on chr20
    and non-PAR X, where they are assembly gaps (each span holds a centromere), and kept on the
    X-degenerate Y, whose GIAB intervals hold no gap and where zero is the signal."""

    model_config = ConfigDict(frozen=True)

    coverage_facet_id: str
    format_id: str
    reason: str | None
    """Why the coverage facet was set aside for this facet; None when it was used."""
    reference_bases: int
    reference_median: float | None
    reference_even_share: float | None
    """Share of chr20's bases whose interval mean is within half to twice its median."""
    x_bases: int
    x_median: float | None
    y_bases: int
    y_median: float | None
    x_ratio: float | None
    """Non-PAR X median over chr20 median; None when either side has under `MIN_COVERAGE_BASES`."""
    y_ratio: float | None
    """X-degenerate Y median over chr20 median; None likewise."""

    @field_validator("reason")
    @classmethod
    def _known_reason(cls, value: str | None) -> str | None:
        if value is not None and value not in UNKNOWN_REASONS:
            raise ValueError(f"unknown reason {value!r}")
        return value


class ArrayLines(BaseModel):
    """What an array export's letters say about its X and Y, read from `gt_alleles` as written (no REF
    is needed). Rows no filter failed; a row is called when it lists alleles and none is null."""

    model_config = ConfigDict(frozen=True)

    autosome_rows: int
    autosome_called: int
    autosome_no_calls: int
    """Rows written as no-calls (`--`, a letter null): what says the export writes its no-calls. A row
    whose genotype could not be read at all (`gt_alleles` null) is neither called nor one."""
    autosome_base_pairs: int
    """Called two-letter rows whose letters are both bases (A, C, G, T)."""
    autosome_hets: int
    x_rows: int
    """Non-PAR X rows."""
    x_called: int
    x_called_one_letter: int
    """Called non-PAR X rows written with one letter: the vendor's decision that the sample has one X,
    never evidence of it."""
    x_base_pairs: int
    x_hets: int
    y_rows: int
    """Rows inside the X-degenerate Y."""
    y_called: int
    x_het_ratio: float | None
    """The het share of non-PAR X base pairs over the autosomes' het share. None when X is not written
    with two letters on most called rows, or either side is under its minimum."""
    y_call_ratio: float | None
    """The called share of X-degenerate Y rows over the autosomes' called share; None under the minimum,
    or when the export wrote no no-call on the autosomes (it drops them, so rows are not probes)."""


class SexChromosomeEvidence(BaseModel):
    """What one facet's calls say about its sex chromosomes, and the estimate it supports.

    Truthy only when a karyotype was estimated. `x_copies` can be known when the karyotype is not:
    a facet whose Y line is unusable still says one X or two.
    """

    model_config = ConfigDict(frozen=True)

    facet_id: str
    assembly_label: str | None
    """The facet's build, as the frames carry it (module docstring): never inferred or defaulted here."""
    assembly_label_basis: AssemblyLabelBasis | None = None
    """`facets.assembly_label_basis`, verbatim: `stated` or `inferred`."""
    assembly_label_source: str | None = None
    """Which of `ASSEMBLY_LABEL_SOURCES` filled the label (`assembly_label_source`); None without one."""
    ref_seqcol_digest: str | None = None
    """`facets.ref_seqcol_digest`, verbatim: the collection `resolve_ref` read REF from, which the
    facet's contigs were checked against (by length, or by label); None before REF is read."""
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
    """`True` when a Y depth line shows a Y; `False` only when a coverage facet shows none and no line
    shows one; otherwise None (module docstring)."""
    y_reason: str | None
    x_copy_count: int | None = None
    """The number of X copies every X line allows, where that is one number: a coverage line counts one
    to three, the calls' X depth one or two (above two copies' band it allows three or more), the het
    share alone only one. None where `x_copies` is None or more than one count fits."""
    y_copy_count: int | None = None
    """The same for Y, zero to two, from a coverage line only: the calls' Y depth shows a Y, never how
    many."""
    inferred_karyotypic_sex: KaryotypicSex
    """A Phenopackets `KaryotypicSex` value, estimated from calls: not a cytogenetic result."""
    reason: str | None
    """Why `inferred_karyotypic_sex` is unknown; None when it is not."""
    detail: str
    """The same in words, with the numbers."""
    coverage: CoverageRatios | None = None
    """The coverage facet's line, when `sex_chromosome_evidence` was given one."""
    array: ArrayLines | None = None
    """An array facet's letter counts (`callset_type = array`); None for every other facet."""

    @field_validator("x_copies_reason", "y_reason", "reason")
    @classmethod
    def _known_reason(cls, value: str | None) -> str | None:
        if value is not None and value not in UNKNOWN_REASONS:
            raise ValueError(f"unknown reason {value!r}")
        return value

    @field_validator("assembly_label_source")
    @classmethod
    def _known_source(cls, value: str | None) -> str | None:
        if value is not None and value not in ASSEMBLY_LABEL_SOURCES:
            raise ValueError(f"unknown assembly label source {value!r}")
        return value

    def __bool__(self) -> bool:
        return self.inferred_karyotypic_sex is not KaryotypicSex.UNKNOWN_KARYOTYPE


def sex_chromosome_evidence(loaded: Loaded, coverage: Coverage | None = None, *, array_y: bool = False) -> tuple[SexChromosomeEvidence, ...]:
    """Measure each loaded facet's sex-chromosome evidence and estimate its karyotype (module docstring).

    `loaded` is what `just_vcf.load` (or `rejoin_split`, `resolve_ref`, `liftover`) returns. The build
    is the facet's `assembly_label` as those wrote it; an array table that states none answers
    `no_build` until `resolve_ref` has inferred one from its rows, so pass that output.

    The result has one entry per row of `loaded.facets`, in that order: a freeBayes+VarScan2 merge loads as two facets and gets two
    answers, never one blended from both callers. The calls are read lazily, a few columns at a time.

    `coverage` is a coverage facet of the same person (`read_coverage`, `alignment_coverage`); pairing
    the two is the caller's call, since the schema has no person. It adds a read-depth line to X and
    Y for each facet on its own build, and is set aside for any other facet, with the reason in
    `coverage.reason`.

    `array_y` reads an array facet's X-degenerate Y probes as a Y line (module docstring). It is off by
    default: a vendor that no-calls Y for the samples it calls female makes the line its own sex call.
    Turn it on for an export whose vendor does not (on PGP-HMS, AncestryDNA's women called a few Y
    probes; 23andMe's none). It changes nothing for other facets.

    The answer is a karyotype, not a sex: an XX male reads two X copies (and `XX` where a coverage facet
    shows no Y), never `XY` for his phenotype. It is an estimate from variant calls and read depth,
    never a cytogenetic or clinical result.
    """
    facets = loaded.facets.select(
        "facet_id", "assembly_label", "assembly_label_basis", "relation", "header_contigs", "vcf_version", "ref_seqcol_digest",
        "callset_type", "caller", "ploidy_statements",
    ).collect()
    measured = measure_coverage(coverage) if coverage is not None else None
    return tuple(
        estimate(
            measure(loaded.calls.filter(pl.col("facet_id") == facet["facet_id"]), facet["assembly_label"]),
            facet["facet_id"],
            facet["assembly_label"],
            CallsetType(facet["callset_type"]) if facet["callset_type"] is not None else None,
            facet["caller"],
            facet["ploidy_statements"],
            measured if measured is None else _paired(measured, coverage.facet.assembly_label, facet["assembly_label"]),
            array_y,
            assembly_label_basis=AssemblyLabelBasis(facet["assembly_label_basis"]) if facet["assembly_label_basis"] is not None else None,
            assembly_label_source=assembly_label_source(facet),
            ref_seqcol_digest=facet["ref_seqcol_digest"],
        )
        for facet in facets.iter_rows(named=True)
    )


def measure_coverage(coverage: Coverage) -> CoverageRatios:
    """A coverage facet's region medians and ratios, on the coverage's own build (`CoverageRatios`)."""
    label = coverage.facet.assembly_label
    base = {"coverage_facet_id": coverage.facet.facet_id, "format_id": coverage.facet.format_id.value}
    if label not in ESTIMATE_BUILDS:
        return CoverageRatios(
            **base, reason="coverage_no_build" if label is None else "unsupported_build", reference_bases=0, reference_median=None,
            reference_even_share=None, x_bases=0, x_median=None, y_bases=0, y_median=None, x_ratio=None, y_ratio=None,
        )
    (_, par1_end), (par2_start, _) = PAR_BOUNDS[label]["X"]
    start, end = pl.col("pos").cast(pl.Int64), pl.col("end_pos").cast(pl.Int64)
    region = (
        pl.when(pl.col("chrom") == DEPTH_REFERENCE_CHROM)
        .then(pl.lit("reference"))
        .when((pl.col("chrom") == "X") & (start > par1_end) & (end < par2_start))
        .then(pl.lit("x"))
        .when((pl.col("chrom") == "Y") & pl.any_horizontal([(start >= s) & (end <= e) for s, e in X_DEGENERATE_Y[label]]))
        .then(pl.lit("y"))
    )
    histogram = (
        coverage.intervals.filter(pl.col("depth_mean").is_not_null() & (pl.col("bases").cast(pl.Int64) == end - start + 1))
        .with_columns(region=region)
        .filter(pl.col("region").is_not_null() & ((pl.col("region") == "y") | (pl.col("depth_mean") > 0)))
        .group_by("region", "depth_mean")
        .agg(pl.col("bases").cast(pl.Int64).sum())
        .collect()
    )
    (ref_bases, ref_median), (x_bases, x_median), (y_bases, y_median) = (
        _weighted_median(histogram.filter(pl.col("region") == name)) for name in ("reference", "x", "y")
    )
    reference = histogram.filter(pl.col("region") == "reference")
    even = (
        float(reference.filter(pl.col("depth_mean").is_between(ref_median / 2, ref_median * 2))["bases"].sum()) / ref_bases
        if ref_bases and ref_median
        else None
    )
    uneven = even is not None and even < COVERAGE_MIN_EVEN_SHARE
    reference_ok = ref_bases >= MIN_COVERAGE_BASES and bool(ref_median) and not uneven

    def ratio(bases: int, median: float | None) -> float | None:
        return median / ref_median if reference_ok and bases >= MIN_COVERAGE_BASES and median is not None else None

    return CoverageRatios(
        **base, reason="coverage_uneven" if uneven else None, reference_bases=ref_bases, reference_median=ref_median,
        reference_even_share=even, x_bases=x_bases, x_median=x_median,
        y_bases=y_bases, y_median=y_median, x_ratio=ratio(x_bases, x_median), y_ratio=ratio(y_bases, y_median),
    )


def _weighted_median(histogram: pl.DataFrame) -> tuple[int, float | None]:
    """The total bases and the length-weighted median of `depth_mean, bases` rows."""
    total = int(histogram["bases"].sum()) if histogram.height else 0
    if total == 0:
        return 0, None
    ordered = histogram.sort("depth_mean").with_columns(cumulative=pl.col("bases").cum_sum())
    return total, float(ordered.filter(pl.col("cumulative") * 2 >= total)["depth_mean"][0])


def _paired(measured: CoverageRatios, coverage_label: str | None, facet_label: str | None) -> CoverageRatios:
    """The coverage line as one facet may use it: set aside unless both builds are known and equal."""
    if measured.reason is not None or coverage_label == facet_label:
        return measured
    return measured.model_copy(update={"reason": "coverage_build_mismatch", "x_ratio": None, "y_ratio": None})


def measure(calls: pl.LazyFrame, assembly_label: str | None) -> dict:
    """The raw measurements of a `calls` frame, per group, as plain numbers (one small collect)."""
    supported = assembly_label in REGION_BUILDS
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
    # Letters as written, for arrays (no REF before `resolve_ref`, so no `gt`): a row is called when it
    # lists alleles and none is null; a two-letter base call is het when its letters differ.
    letters = pl.col("gt_alleles")
    called = (n_haplotypes > 0) & letters.list.eval(pl.element().is_not_null()).list.all().fill_null(False)
    no_call = ((n_haplotypes > 0) & letters.list.eval(pl.element().is_null()).list.any()).fill_null(False)
    base_pair = called & (n_haplotypes == 2) & letters.list.eval(pl.element().is_in(_BASES)).list.all().fill_null(False)
    letter_het = base_pair & (letters.list.get(0, null_on_oob=True) != letters.list.get(1, null_on_oob=True))

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
            called.sum().alias("called"),
            no_call.sum().alias("no_calls"),
            (called & (n_haplotypes == 1)).sum().alias("called_one"),
            base_pair.sum().alias("base_pairs"),
            letter_het.sum().alias("letter_hets"),
            xdeg.sum().alias("xdeg_rows"),
            (called & xdeg).sum().alias("xdeg_called"),
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
    if assembly_label in REGION_BUILDS:
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
    coverage: CoverageRatios | None = None,
    array_y: bool = False,
    *,
    assembly_label_basis: AssemblyLabelBasis | None = None,
    assembly_label_source: str | None = None,
    ref_seqcol_digest: str | None = None,
) -> SexChromosomeEvidence:
    """Turn `measure`'s numbers into the evidence and the estimate (rules in the module docstring).

    `caller` and `ploidy_statements` are the facet's own columns, carried through unread: what the
    caller said about ploidy is its assertion, kept beside the evidence and never counted in it.
    `coverage` is a coverage facet's line (`measure_coverage`); it counts only without a `reason`.
    `assembly_label_basis`, `assembly_label_source` and `ref_seqcol_digest` say where `assembly_label`
    came from; they are carried, and decide nothing."""
    groups = measured["groups"]
    empty = {"rows": 0, "one": 0, "two": 0, "more": 0}
    names = ["autosomes", "x_nonpar", "y_nonpar", "par", "mt"] if measured["supported"] else ["autosomes", "mt"]
    ploidy = {name: GtPloidyCounts(**{key: int(groups.get(name, empty)[key]) for key in ("one", "two", "more")}) for name in names}
    base = {
        "facet_id": facet_id,
        "assembly_label": assembly_label,
        "assembly_label_basis": assembly_label_basis,
        "assembly_label_source": assembly_label_source,
        "ref_seqcol_digest": ref_seqcol_digest,
        "callset_type": callset_type,
        "caller": caller,
        "ploidy_statements": ploidy_statements,
        "rows_counted": int(sum(row["rows"] for row in groups.values())),
        "ploidy": ploidy,
        "coverage": coverage,
    }
    usable = coverage if coverage is not None and coverage.reason is None else None

    def unknown(reason: str, detail: str, **known) -> SexChromosomeEvidence:
        fields = {"region_set": None, "x_heterozygosity": None, "y_heterozygosity": None, "depth": None, "x_copies": None, "x_copies_reason": reason, "y_present": None, "y_reason": reason}
        return SexChromosomeEvidence(**base, **{**fields, **known}, inferred_karyotypic_sex=KaryotypicSex.UNKNOWN_KARYOTYPE, reason=reason, detail=detail)

    def named(x_range: Copies | None, x_reason: str | None, y_range: Copies | None, y_reason: str | None, detail: str) -> SexChromosomeEvidence:
        lines = {
            "x_copies": None if x_range is None else (XCopies.ONE if x_range[0] == 1 else XCopies.TWO_OR_MORE),
            "x_copies_reason": x_reason,
            "y_present": None if y_range is None else y_range[0] >= 1,
            "y_reason": y_reason,
            "x_copy_count": _count(x_range),
            "y_copy_count": _count(y_range),
        }
        if x_range is None:
            return unknown(x_reason, detail, **known, **lines)
        if y_range is None:
            return unknown("no_depth" if y_reason == "no_depth" else "y_unusable", detail, **known, **lines)
        x_read = _count(x_range) or _UNCOUNTED_READING.get(x_range)
        y_read = _count(y_range) if _count(y_range) is not None else _UNCOUNTED_READING.get(y_range)
        if x_read is None or y_read is None:
            return unknown("copies_uncounted", detail, **known, **lines)
        return SexChromosomeEvidence(**base, **known, **lines, inferred_karyotypic_sex=KARYOTYPES[(x_read, y_read)], reason=None, detail=detail)

    if assembly_label is None:
        return unknown(
            "no_build",
            "the facet has no assembly label, so the pseudo-autosomal regions cannot be placed "
            "(resolve_ref infers an array table's build from its rows; estimate its output)",
        )
    if not measured["supported"]:
        return unknown("unsupported_build", f"no PAR bounds or Y region set for {assembly_label} (only {', '.join(sorted(REGION_BUILDS))})")

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
    known = {"region_set": REGION_SET[assembly_label], "x_heterozygosity": x_het, "y_heterozygosity": y_het, "depth": depth}
    if callset_type is CallsetType.ARRAY:
        lines = _array_lines(groups)
        known["array"] = lines
        if assembly_label not in ESTIMATE_BUILDS:
            return unknown("uncalibrated_build", f"no array rule was checked on {assembly_label}", **known)
        x_range, x_reason, x_detail = _array_x(lines)
        y_range, y_reason, y_detail = _array_y(lines, array_y)
        return named(x_range, x_reason, y_range, y_reason, f"X: {x_detail}; Y: {y_detail}")
    if callset_type is not CallsetType.WGS:
        kind = callset_type.value if callset_type else "unknown"
        return unknown("not_wgs", f"callset type {kind}; the thresholds were measured on whole genomes", **known)
    if assembly_label not in ESTIMATE_BUILDS:
        measured_here = (
            f"X het share {_fmt(x_het.het_share)}, X depth ratio {_fmt(depth.x_ratio)}, "
            f"X-degenerate Y depth ratio {_fmt(depth.y_ratio)}"
        )
        calibrated = " and ".join(sorted(ESTIMATE_BUILDS))
        return unknown("uncalibrated_build", f"the thresholds were measured on {calibrated} genomes, none on {assembly_label}; measured: {measured_here}", **known)

    no_dp = int(auto.get("depth_sites", 0)) + x_sites + int(y.get("depth_sites", 0)) == 0
    x_range, x_reason, x_detail = _x_copies(x_het, depth, no_dp, usable)
    y_range, y_reason, y_detail = _y_present(depth, no_dp, usable)
    return named(x_range, x_reason, y_range, y_reason, f"X: {x_detail}; Y: {y_detail}")



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


Copies = tuple[int, int | None]
"""The copy counts a line allows, lowest and highest; None as the highest is no upper bound."""

Band = tuple[float | None, float | None, Copies]
"""A line's value range (inclusive, None open) and the counts it allows there."""

HET_SHARE_BANDS: tuple[Band, ...] = ((None, ONE_X_MAX_HET_SHARE, (1, 1)), (TWO_X_MIN_HET_SHARE, None, (2, None)))
DEPTH_X_BANDS: tuple[Band, ...] = (
    (None, ONE_X_MAX_DEPTH_RATIO, (1, 1)),
    (TWO_X_MIN_DEPTH_RATIO, TWO_X_MAX_DEPTH_RATIO, (2, 2)),
    (TWO_X_MAX_DEPTH_RATIO, None, (3, None)),
)
"""The calls' X depth: one, two, or more than two (the first band holding a value wins, so 1.25 is two)."""


def _coverage_bands(counts: tuple[int, ...]) -> tuple[Band, ...]:
    """A band per count around its copy-number value; a count of zero has no lower edge."""
    return tuple(
        (
            None if n == 0 else n * COVERAGE_RATIO_PER_COPY - COVERAGE_BAND_HALF_WIDTH,
            n * COVERAGE_RATIO_PER_COPY + COVERAGE_BAND_HALF_WIDTH,
            (n, n),
        )
        for n in counts
    )


COVERAGE_X_BANDS = _coverage_bands(COVERAGE_X_COUNTS)
COVERAGE_Y_BANDS = _coverage_bands(COVERAGE_Y_COUNTS)

TWO_X_EXPECTED_HET_RATIO = 0.75
"""What two X should read against the autosomes: X diversity is about 3/4 of an autosome's, since a
population with as many men as women carries three X for every four copies of an autosome
(Ne_X = 3/4 Ne_A, the neutral expectation). One X reads 0, genotyping errors aside."""

ARRAY_ONE_X_MAX_HET_RATIO = TWO_X_EXPECTED_HET_RATIO / 4
ARRAY_TWO_X_MIN_HET_RATIO = TWO_X_EXPECTED_HET_RATIO * 3 / 4
"""An array's X het share over its autosomal het share, banded at the central half of the gap
between 0 and `TWO_X_EXPECTED_HET_RATIO` (0.1875 and 0.5625), the rule the other lines follow.
Rejected: an expectation of 1 (X as diverse as an autosome), which gives 0.25 and 0.75 like PLINK
`--check-sex`'s F cut-offs of 0.2 and 0.8 (F from expected het by allele frequency, which a single
export has not). Two X has no upper edge: hets cannot count copies. PGP-HMS, docs/ROADMAP.md RM269."""

ARRAY_Y_ABSENT_MAX_CALL_RATIO = 0.25
ARRAY_Y_PRESENT_MIN_CALL_RATIO = 0.75
"""The X-degenerate Y probes' called share over the autosomes'. A Y calls like an autosome (near 1),
none leaves its probes uncalled (near 0); the unknown band is the gap's central half. Read only with
`array_y` (the vendor's no-calls on Y follow its own sex call)."""

ARRAY_X_BANDS: tuple[Band, ...] = ((None, ARRAY_ONE_X_MAX_HET_RATIO, (1, 1)), (ARRAY_TWO_X_MIN_HET_RATIO, None, (2, None)))
ARRAY_Y_BANDS: tuple[Band, ...] = ((None, ARRAY_Y_ABSENT_MAX_CALL_RATIO, (0, 0)), (ARRAY_Y_PRESENT_MIN_CALL_RATIO, None, (1, None)))

KARYOTYPES: dict[tuple[int, int], KaryotypicSex] = {
    (1, 0): KaryotypicSex.XO, (1, 1): KaryotypicSex.XY, (1, 2): KaryotypicSex.XYY,
    (2, 0): KaryotypicSex.XX, (2, 1): KaryotypicSex.XXY, (2, 2): KaryotypicSex.XXYY,
    (3, 0): KaryotypicSex.XXX, (3, 1): KaryotypicSex.XXXY, (3, 2): KaryotypicSex.OTHER_KARYOTYPE,
}  # fmt: skip
"""(X copies, Y copies) to the Phenopackets value. Phenopackets lists no XXXYY, so three X with two Y
is its `OTHER_KARYOTYPE`."""

_UNCOUNTED_READING: dict[Copies, int] = {(2, None): 2, (1, None): 1}
"""How an uncounted range reads in the karyotype, as it did before counts were measured: a het share's
two or more as two, a Y the calls show as one. Three or more X with nothing counting them reads
nothing (`copies_uncounted`)."""


def _read(value: float | None, bands: tuple[Band, ...]) -> Copies | str | None:
    """The counts of the first band holding `value`; `intermediate` in none; None with no value."""
    if value is None:
        return None
    for low, high, copies in bands:
        if (low is None or value >= low) and (high is None or value <= high):
            return copies
    return "intermediate"


def _agree(found: list[Copies | str]) -> Copies | str | None:
    """The counts every line allows: `intermediate` if one fell between bands, `lines_disagree` if
    they share none, None without a line."""
    if "intermediate" in found:
        return "intermediate"
    ranges = [line for line in found if not isinstance(line, str)]
    if not ranges:
        return None
    low = max(lo for lo, _ in ranges)
    highs = [hi for _, hi in ranges if hi is not None]
    high = min(highs) if highs else None
    return "lines_disagree" if high is not None and low > high else (low, high)


def _count(copies: Copies | None) -> int | None:
    """The count a range holds when it is one number."""
    return copies[0] if copies is not None and copies[0] == copies[1] else None


def _x_copies(
    x_het: Heterozygosity, depth: DepthRatios, no_dp: bool, coverage: CoverageRatios | None
) -> tuple[Copies | None, str | None, str]:
    """The X counts the het share, the X depth ratio and the coverage ratio all allow."""
    written_diploid = x_het.diploid_alt_sites >= MIN_X_SITES and x_het.diploid_alt_sites > x_het.haploid_alt_sites
    share = x_het.het_share if written_diploid else None
    by_share = _read(share, HET_SHARE_BANDS)
    by_depth = _read(depth.x_ratio, DEPTH_X_BANDS)
    by_coverage = _read(coverage.x_ratio if coverage else None, COVERAGE_X_BANDS)
    said = (
        f"het share {_fmt(x_het.het_share)} over {x_het.diploid_alt_sites:,} diploid ALT call(s), "
        f"{x_het.haploid_alt_sites:,} haploid; depth ratio {_fmt(depth.x_ratio)} over {depth.x_sites:,} call(s) with DP"
    )
    if coverage is not None:
        said += f"; coverage ratio {_fmt(coverage.x_ratio)} over {coverage.x_bases:,} base(s)"
    agreed = _agree([line for line in (by_share, by_depth, by_coverage) if line is not None])
    if isinstance(agreed, str):
        return None, agreed, said
    if agreed is not None:
        return agreed, None, said
    if x_het.haploid_alt_sites >= MIN_X_SITES:
        return None, "caller_asserted_only", said
    return None, ("no_depth" if no_dp and coverage is None else "too_few_sites"), said


def _array_lines(groups: dict) -> ArrayLines:
    """An array facet's letter counts and the two ratios (`ArrayLines`)."""
    auto, x, y = groups.get("autosomes", {}), groups.get("x_nonpar", {}), groups.get("y_nonpar", {})
    count = {
        "autosome_rows": int(auto.get("rows", 0)), "autosome_called": int(auto.get("called", 0)),
        "autosome_no_calls": int(auto.get("no_calls", 0)),
        "autosome_base_pairs": int(auto.get("base_pairs", 0)), "autosome_hets": int(auto.get("letter_hets", 0)),
        "x_rows": int(x.get("rows", 0)), "x_called": int(x.get("called", 0)), "x_called_one_letter": int(x.get("called_one", 0)),
        "x_base_pairs": int(x.get("base_pairs", 0)), "x_hets": int(x.get("letter_hets", 0)),
        "y_rows": int(y.get("xdeg_rows", 0)), "y_called": int(y.get("xdeg_called", 0)),
    }  # fmt: skip
    auto_share = count["autosome_hets"] / count["autosome_base_pairs"] if count["autosome_base_pairs"] >= MIN_REFERENCE_SITES else None
    auto_called = count["autosome_called"] / count["autosome_rows"] if count["autosome_rows"] >= MIN_REFERENCE_SITES else None
    written_two = count["x_base_pairs"] >= MIN_X_SITES and count["x_base_pairs"] > count["x_called_one_letter"]
    x_ratio = count["x_hets"] / count["x_base_pairs"] / auto_share if written_two and auto_share else None
    # Rows count the probes only where the export writes its no-calls; one that drops them has nothing
    # but called Y rows, and a woman's cross-hybridising calls would read as a Y. A written no-call says
    # so; a row that is merely not called (its genotype unreadable, `gt_alleles` null) does not.
    writes_no_calls = count["autosome_no_calls"] > 0
    y_ratio = count["y_called"] / count["y_rows"] / auto_called if count["y_rows"] >= MIN_Y_SITES and auto_called and writes_no_calls else None
    return ArrayLines(**count, x_het_ratio=x_ratio, y_call_ratio=y_ratio)


def _array_x(lines: ArrayLines) -> tuple[Copies | None, str | None, str]:
    """X copies from the het share of two-letter X calls over the autosomes'."""
    said = (
        f"het ratio {_fmt(lines.x_het_ratio)} ({lines.x_hets:,} het of {lines.x_base_pairs:,} two-letter calls on non-PAR X, "
        f"{lines.x_called_one_letter:,} written with one letter)"
    )
    found = _read(lines.x_het_ratio, ARRAY_X_BANDS)
    if isinstance(found, str):
        return None, found, said
    if found is not None:
        return found, None, said
    if lines.x_called_one_letter >= MIN_X_SITES:
        return None, "caller_asserted_only", said
    return None, "too_few_sites", said


def _array_y(lines: ArrayLines, array_y: bool) -> tuple[Copies | None, str | None, str]:
    """Y from the called share of X-degenerate Y probes, read only when `array_y` is set."""
    said = f"call ratio {_fmt(lines.y_call_ratio)} ({lines.y_called:,} of {lines.y_rows:,} X-degenerate Y rows called)"
    if lines.y_call_ratio is None:
        if lines.y_rows >= MIN_Y_SITES and lines.autosome_rows >= MIN_REFERENCE_SITES and lines.autosome_no_calls == 0:
            return None, "no_calls_not_written", said + "; the export writes no no-call, so its rows are not its probes"
        return None, "too_few_sites", said
    if not array_y:
        return None, "vendor_asserted", said + "; not read: the vendor's Y calls and no-calls follow its own sex call"
    found = _read(lines.y_call_ratio, ARRAY_Y_BANDS)
    if isinstance(found, str):
        return None, found, said
    return found, None, said + "; read as the vendor's calls (array_y)"


def _y_present(depth: DepthRatios, no_dp: bool, coverage: CoverageRatios | None) -> tuple[Copies | None, str | None, str]:
    """The Y counts the X-degenerate Y lines allow: the calls' DP shows a Y (one or more) or nothing,
    a coverage line counts none, one or two."""
    said = f"X-degenerate depth ratio {_fmt(depth.y_ratio)} over {depth.y_sites:,} call(s) with DP"
    if no_dp:
        by_calls, calls_reason = None, "no_depth"
    elif depth.y_ratio is None:
        by_calls, calls_reason = None, "too_few_sites"
    elif depth.y_ratio >= Y_PRESENT_MIN_DEPTH_RATIO:
        by_calls, calls_reason = (1, None), None
    else:
        by_calls, calls_reason = None, "y_not_shown"
    by_coverage = None
    if coverage is not None:
        said += f"; coverage ratio {_fmt(coverage.y_ratio)} over {coverage.y_bases:,} base(s)"
        by_coverage = _read(coverage.y_ratio, COVERAGE_Y_BANDS)
    agreed = _agree([line for line in (by_calls, by_coverage) if line is not None])
    if isinstance(agreed, str):
        return None, agreed, said
    if agreed is None:
        return None, calls_reason, said
    return agreed, None, said


def _fmt(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.3f}"
