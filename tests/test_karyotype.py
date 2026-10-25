"""Sex-chromosome evidence and the karyotype estimate it supports (RM269)."""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import (
    BLOCKS,
    FACETS,
    Coverage,
    CoverageFacet,
    CoverageFormat,
    KaryotypicSex,
    Loaded,
    Resolution,
    XCopies,
    load,
    sex_chromosome_evidence,
)
from just_vcf.coverage import COVERAGE_INTERVALS
from just_vcf.karyotype import (
    COVERAGE_X_BANDS,
    COVERAGE_Y_BANDS,
    MIN_X_SITES,
    MIN_Y_SITES,
    ESTIMATE_BUILDS,
    ONE_X_MAX_HET_SHARE,
    PAR_BOUNDS,
    REGION_BUILDS,
    REGION_SET,
    TWO_X_MAX_DEPTH_RATIO,
    TWO_X_MIN_HET_SHARE,
    UNKNOWN_REASONS,
    X_DEGENERATE_Y,
    SexChromosomeEvidence,
    estimate,
    measure,
)
from just_vcf.schema import CALLS, AssemblyLabelBasis, CallsetType, conform, empty

ASSETS = Path(__file__).resolve().parent.parent / "assets"
GIAB = ASSETS / "giab_stratifications_v3.6_xy"
T2T_CLASSES = ASSETS / "t2t_chm13v2.0_annotation" / "chm13v2.0_chrXY_sequence_class_v1.bed"
GIAB_PREFIX = {"GRCh37": "GRCh37", "GRCh38": "GRCh38", "T2T-CHM13": "CHM13"}
"""The assembly label → the prefix GIAB names its files with."""

# --- The regions, re-derived from the GIAB files ---------------------------------------------------


def bed(build: str, name: str) -> list[tuple[int, int]]:
    """A GIAB BED as 0-based half-open `(start, end)` pairs, sorted."""
    lines = (GIAB / f"{GIAB_PREFIX[build]}_{name}.bed").read_text().splitlines()
    return sorted((int(start), int(end)) for _, start, end, *_ in (line.split("\t") for line in lines if line.strip()))


def subtract(keep: list[tuple[int, int]], cut: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Half-open interval difference."""
    out = []
    for start, end in keep:
        pieces = [(start, end)]
        for cut_start, cut_end in cut:
            pieces = [
                part
                for piece_start, piece_end in pieces
                for part in (
                    [(piece_start, piece_end)]
                    if cut_end <= piece_start or cut_start >= piece_end
                    else [(piece_start, cut_start)] * (cut_start > piece_start) + [(cut_end, piece_end)] * (cut_end < piece_end)
                )
            ]
        out += pieces
    return sorted(out)


def to_closed(start: int, end: int) -> tuple[int, int]:
    """BED (0-based, half-open) to 1-based closed `pos` bounds: the one conversion, tested below."""
    return start + 1, end


def test_bed_conversion_keeps_the_bases() -> None:
    # BED `chrY 10000 10001` is one base, POS 10001.
    assert to_closed(10000, 10001) == (10001, 10001)
    start, end = to_closed(2_781_479, 3_050_044)
    assert end - start + 1 == 3_050_044 - 2_781_479


@pytest.mark.parametrize("build", sorted(PAR_BOUNDS))
def test_par_bounds_are_the_giab_par_arms(build: str) -> None:
    for chrom in ("X", "Y"):
        intervals = bed(build, f"chr{chrom}_PAR")
        middle = 30_000_000
        par1 = [i for i in intervals if i[0] < middle]
        par2 = [i for i in intervals if i[0] >= middle]
        expected = (to_closed(par1[0][0], max(e for _, e in par1)), to_closed(par2[0][0], max(e for _, e in par2)))
        assert PAR_BOUNDS[build][chrom] == expected


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_x_degenerate_y_is_nonpar_minus_xtr_minus_ampliconic(build: str) -> None:
    derived = subtract(bed(build, "chrY_nonPAR"), bed(build, "chrY_XTR") + bed(build, "chrY_ampliconic"))
    (_, par1_end), (par2_start, _) = PAR_BOUNDS[build]["Y"]
    closed = tuple(to_closed(s, e) for s, e in derived)
    kept = tuple(i for i in closed if i[0] > par1_end and i[1] < par2_start)
    assert X_DEGENERATE_Y[build] == kept
    # Only GRCh38's single base before PAR1 falls outside the PARs' span.
    assert len(closed) - len(kept) == (1 if build == "GRCh38" else 0)
    for start, end in kept:
        assert not any(s < end and start <= e for s, e in (to_closed(*i) for i in bed(build, "chrY_XTR") + bed(build, "chrY_ampliconic")))


def t2t_classes(chrom: str) -> dict[str, list[tuple[int, int]]]:
    """The T2T sequence-class BED for one chromosome: class name → 0-based half-open intervals."""
    found: dict[str, list[tuple[int, int]]] = {}
    for line in T2T_CLASSES.read_text().splitlines():
        name, start, end, kind, *_ = line.split("\t")
        if name == f"chr{chrom}":
            found.setdefault(kind, []).append((int(start), int(end)))
    return found


def test_t2t_x_degenerate_y_is_the_t2t_class_less_giabs_xtr_and_par() -> None:
    classes = t2t_classes("Y")
    derived = subtract(sorted(classes["X-DEG"]), bed("T2T-CHM13", "chrY_XTR") + bed("T2T-CHM13", "chrY_PAR"))
    closed = tuple(to_closed(s, e) for s, e in derived)
    assert X_DEGENERATE_Y["T2T-CHM13"] == closed
    # Only the two bases the class file's second X-DEG interval shares with XTR2 were dropped.
    assert sum(e - s for s, e in classes["X-DEG"]) - sum(e - s for s, e in derived) == 2
    (_, par1_end), (par2_start, _) = PAR_BOUNDS["T2T-CHM13"]["Y"]
    assert all(par1_end < start and end < par2_start for start, end in closed)
    nonpar = bed("T2T-CHM13", "chrY_nonPAR")
    assert all(any(s <= start - 1 and end <= e for s, e in nonpar) for start, end in closed)


def test_t2t_classes_agree_with_giab_on_the_par_and_the_xtr() -> None:
    for chrom in ("X", "Y"):
        classes = t2t_classes(chrom)
        assert sorted(classes["PAR1"] + classes["PAR2"]) == bed("T2T-CHM13", f"chr{chrom}_PAR")
    y = t2t_classes("Y")
    assert sorted(y["XTR1"] + y["XTR2"]) == bed("T2T-CHM13", "chrY_XTR")


def test_giabs_t2t_nonpar_less_xtr_is_not_x_degenerate() -> None:
    """The rejected derivation: on a gapless Y, nonPAR minus XTR keeps the ampliconic, heterochromatic,
    centromeric and satellite classes, most of the male-specific Y."""
    rest = subtract(bed("T2T-CHM13", "chrY_nonPAR"), bed("T2T-CHM13", "chrY_XTR"))
    classes = t2t_classes("Y")
    xdeg = sum(e - s for s, e in classes["X-DEG"])
    assert sum(e - s for s, e in rest) > 4 * xdeg
    assert all(any(s <= start and end <= e for s, e in rest) for start, end in classes["HET"] + classes["CEN"])


# --- Synthetic calls ----------------------------------------------------------------------------

ALLELES = ["A", "G", "T"]


def rows(chrom: str, start: int, n: int, gt: list[int | None], dp: int | None, *, step: int = 10, ad: list[int] | None = None,
         filters: list[str] | None = ["PASS"], split: str | None = None) -> pl.DataFrame:  # noqa: B006 - read only
    return pl.DataFrame(
        {
            "facet_id": "f",
            "chrom": chrom,
            "pos": pl.int_range(0, n, eager=True) * step + start,
            "source_ref": "A",
            "alts": [ALLELES[1:]] * n,
            "source_alleles": [ALLELES] * n,
            "gt": [gt] * n,
            "gt_alleles": [[None if i is None else ALLELES[i] for i in gt]] * n,
            "dp": [dp] * n,
            "ad": [ad] * n,
            "filters": [filters] * n,
            "evidence_basis": "assay",
            "split_signal": split,
        },
        schema_overrides={"gt": pl.List(pl.UInt8), "dp": pl.UInt32, "ad": pl.List(pl.UInt32), "pos": pl.UInt32, "filters": pl.List(pl.String)},
    )


def frame(*parts: pl.DataFrame) -> pl.LazyFrame:
    return conform(pl.concat(parts, how="diagonal_relaxed").lazy(), CALLS)


REF_DP = 30
X_START = 5_000_000  # non-PAR X on both builds
Y_START = 2_790_000  # inside GRCh38's first X-degenerate interval (2,781,480–3,050,044)


def chr20(n: int = 2000, dp: int | None = REF_DP) -> pl.DataFrame:
    return rows("20", 1_000_000, n, [0, 1], dp)


def x_calls(het: int, hom: int, dp: int | None) -> pl.DataFrame:
    return pl.concat([rows("X", X_START, het, [0, 1], dp, ad=[15, 15, 0] if dp else None), rows("X", X_START + 1_000_000, hom, [1, 1], dp)])


def y_calls(n: int, dp: int | None, gt: list[int | None] = [1, 1]) -> pl.DataFrame:  # noqa: B006 - read only
    return rows("Y", Y_START, n, gt, dp, step=20)


def evidence(calls: pl.LazyFrame, build: str | None = "GRCh38", kind: CallsetType | None = CallsetType.WGS) -> SexChromosomeEvidence:
    return estimate(measure(calls, build), "f", build, kind)


XX = (x_calls(1200, 800, REF_DP), y_calls(1000, 6))
XY = (x_calls(40, 1960, REF_DP // 2), y_calls(1000, 18))


@pytest.mark.parametrize(
    ("x", "y", "expected", "copies", "reason"),
    [
        (XY[0], XY[1], KaryotypicSex.XY, XCopies.ONE, None),
        (XX[0], XY[1], KaryotypicSex.XXY, XCopies.TWO_OR_MORE, None),
        # A low Y depth does not show a Y absent, so XX and XO are not named; X copies still are.
        (XX[0], XX[1], KaryotypicSex.UNKNOWN_KARYOTYPE, XCopies.TWO_OR_MORE, "y_unusable"),
        (XY[0], XX[1], KaryotypicSex.UNKNOWN_KARYOTYPE, XCopies.ONE, "y_unusable"),
    ],
)
def test_the_classes(x: pl.DataFrame, y: pl.DataFrame, expected: KaryotypicSex, copies: XCopies, reason: str | None) -> None:
    found = evidence(frame(chr20(), x, y))
    assert (found.inferred_karyotypic_sex, found.x_copies, found.reason) == (expected, copies, reason)
    assert bool(found) is (reason is None)
    assert found.y_present is (True if reason is None else None)
    assert found.y_reason == (None if reason is None else "y_not_shown")


def test_measurements_are_what_the_rows_hold() -> None:
    found = evidence(frame(chr20(), *XX))
    het = found.x_heterozygosity
    assert (het.diploid_alt_sites, het.het_sites, het.het_share) == (2000, 1200, 0.6)
    assert (het.balance_sites, het.median_het_allele_balance) == (1200, 0.5)
    assert (found.depth.reference_sites, found.depth.x_sites, found.depth.y_sites) == (2000, 2000, 1000)
    assert (found.depth.x_ratio, found.depth.y_ratio) == (1.0, 6 / REF_DP)
    assert found.ploidy["x_nonpar"].two == 2000 and found.ploidy["y_nonpar"].two == 1000
    assert (found.y_heterozygosity.diploid_alt_sites, found.y_heterozygosity.het_share) == (1000, 0.0)
    assert found.rows_counted == 5000


def test_par_rows_count_as_par_not_as_x() -> None:
    # Hets across PAR1 on both arms would push a one-X sample's het share into the two-X band.
    (par1_start, _), (par2_start, _) = PAR_BOUNDS["GRCh38"]["X"]
    par = pl.concat([rows("X", par1_start, 3000, [0, 1], REF_DP), rows("X", par2_start, 3000, [0, 1], REF_DP), rows("XY", 100, 50, [0, 1], REF_DP)])
    found = evidence(frame(chr20(), *XY, par))
    assert found.inferred_karyotypic_sex is KaryotypicSex.XY
    assert found.ploidy["par"].two == 6050
    assert found.x_heterozygosity.diploid_alt_sites == 2000


def test_the_par_bounds_follow_the_build() -> None:
    # X:2,700,000 is PAR1 on GRCh38 (it ends at 2,781,479) and non-PAR on GRCh37 (it ends at 2,699,520).
    calls = frame(rows("X", 2_700_000, 10, [0, 1], REF_DP))
    assert evidence(calls, "GRCh38").ploidy["par"].two == 10
    assert evidence(calls, "GRCh37").ploidy["x_nonpar"].two == 10


def test_y_depth_counts_only_the_x_degenerate_region() -> None:
    # The same Y depth in the X-transposed region (GRCh38 3,050,045–6,235,111) is not counted.
    xtr = rows("Y", 4_000_000, 1000, [1, 1], 18, step=20)
    found = evidence(frame(chr20(), XY[0], xtr))
    assert found.depth.y_sites == 0
    assert found.y_heterozygosity.diploid_alt_sites == 0
    assert found.ploidy["y_nonpar"].two == 1000
    assert (found.y_present, found.y_reason, found.reason) == (None, "too_few_sites", "y_unusable")
    assert found.x_copies is XCopies.ONE


def test_a_sex_aware_callers_own_diploid_minority_is_not_a_het_share() -> None:
    # PGP-HMS huBE0518 (DRAGEN): 103,154 haploid X ALT calls and 1,315 diploid ones, all hets. The
    # caller picked those; read as a share they say two X against a one-X depth.
    x = pl.concat([rows("X", X_START, 20_000, [1], REF_DP // 2), rows("X", X_START + 1_000_000, 1500, [0, 1], REF_DP // 2)])
    found = evidence(frame(chr20(), x, y_calls(1000, 18, gt=[1])))
    assert found.x_heterozygosity.het_share == 1.0
    assert (found.x_copies, found.inferred_karyotypic_sex) == (XCopies.ONE, KaryotypicSex.XY)


def test_a_sex_aware_callers_haploid_x_is_decided_by_depth_alone() -> None:
    x = rows("X", X_START, 2000, [1], REF_DP // 2)
    found = evidence(frame(chr20(), x, y_calls(1000, 18, gt=[1])))
    assert found.x_heterozygosity.haploid_alt_sites == 2000
    assert found.x_heterozygosity.diploid_alt_sites == 0
    assert found.ploidy["x_nonpar"].one == 2000
    assert found.inferred_karyotypic_sex is KaryotypicSex.XY


def test_haploid_genotypes_alone_are_the_callers_assertion() -> None:
    x = rows("X", X_START, 2000, [1], None)
    found = evidence(frame(chr20(dp=None), x))
    assert (found.x_copies, found.x_copies_reason) == (None, "caller_asserted_only")
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE


def test_without_dp_the_het_share_still_gives_x_copies_and_the_karyotype_waits() -> None:
    found = evidence(frame(chr20(dp=None), x_calls(1200, 800, None), y_calls(1000, None)))
    assert found.x_copies is XCopies.TWO_OR_MORE
    assert (found.y_reason, found.reason, found.inferred_karyotypic_sex) == ("no_depth", "no_depth", KaryotypicSex.UNKNOWN_KARYOTYPE)


def test_too_few_x_sites_is_unknown() -> None:
    few = MIN_X_SITES - 1
    found = evidence(frame(chr20(), rows("X", X_START, few, [0, 1], REF_DP), y_calls(1000, 6)))
    assert (found.x_copies, found.reason) == (None, "too_few_sites")
    assert found.depth.x_ratio is None


def test_too_few_y_sites_is_unknown() -> None:
    found = evidence(frame(chr20(), XY[0], y_calls(MIN_Y_SITES - 1, 18)))
    assert (found.y_reason, found.reason) == ("too_few_sites", "y_unusable")
    assert found.x_copies is XCopies.ONE


def test_a_share_between_the_bands_is_intermediate_not_the_nearest_class() -> None:
    share = (ONE_X_MAX_HET_SHARE + TWO_X_MIN_HET_SHARE) / 2
    het = round(2000 * share)
    found = evidence(frame(chr20(), x_calls(het, 2000 - het, REF_DP), y_calls(1000, 6)))
    assert (found.x_copies_reason, found.reason) == ("intermediate", "intermediate")


def test_lines_that_disagree_are_unknown() -> None:
    # A two-X het share at one-X depth.
    found = evidence(frame(chr20(), x_calls(1200, 800, REF_DP // 2), y_calls(1000, 6)))
    assert (found.x_copies, found.reason) == (None, "lines_disagree")


def test_depth_above_two_copies_is_not_called_xx() -> None:
    # Three X copies would put the X depth near 1.5 times the autosomes', with a two-X het share. The
    # calls' depth allows three or more there, never two, and counts nothing: no XX, no XXY.
    dp = round(REF_DP * 1.5)
    assert dp / REF_DP > TWO_X_MAX_DEPTH_RATIO
    found = evidence(frame(chr20(), x_calls(1200, 800, dp), y_calls(1000, 6)))
    assert (found.x_copies, found.x_copy_count, found.reason) == (XCopies.TWO_OR_MORE, None, "y_unusable")
    with_y = evidence(frame(chr20(), x_calls(1200, 800, dp), y_calls(1000, 18)))
    assert (with_y.inferred_karyotypic_sex, with_y.reason, with_y.y_present) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "copies_uncounted", True)


def test_calls_alone_count_x_but_never_y() -> None:
    # The calls' X depth counts one or two; their Y depth shows a Y, never how many (men read 0.43-1.97).
    xy, xxy = evidence(frame(chr20(), *XY)), evidence(frame(chr20(), XX[0], XY[1]))
    assert [(e.inferred_karyotypic_sex, e.x_copy_count, e.y_copy_count) for e in (xy, xxy)] == [
        (KaryotypicSex.XY, 1, None), (KaryotypicSex.XXY, 2, None),
    ]  # fmt: skip
    # Without DP only the het share speaks: two or more X, no count.
    assert evidence(frame(chr20(dp=None), x_calls(1200, 800, None))).x_copy_count is None


def test_a_low_y_depth_never_reads_as_no_y() -> None:
    # PGP-HMS hu3BDC4B, self-reported male, GATK on bowtie2 (GRCh37): one X, X-degenerate Y ratio
    # 0.433, above every two-X genome measured (at most 0.400). Depth at called sites cannot say a Y
    # is absent, so this must not read XO, and an XX facet must not read XX either.
    dp = 13  # 13/30 = 0.433
    found = evidence(frame(chr20(), XY[0], y_calls(1000, dp)))
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE
    assert (found.x_copies, found.y_present, found.y_reason, found.reason) == (XCopies.ONE, None, "y_not_shown", "y_unusable")
    assert evidence(frame(chr20(), *XX)).inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE


def test_the_y_present_band_starts_at_its_threshold() -> None:
    # 16/30 = 0.533 is present, 15/30 = 0.5 is not shown (Y_PRESENT_MIN_DEPTH_RATIO = 0.52).
    assert evidence(frame(chr20(), XX[0], y_calls(1000, 16))).inferred_karyotypic_sex is KaryotypicSex.XXY
    found = evidence(frame(chr20(), XX[0], y_calls(1000, 15)))
    assert (found.y_reason, found.reason, found.x_copies) == ("y_not_shown", "y_unusable", XCopies.TWO_OR_MORE)


def test_failed_filters_are_left_out_and_unapplied_ones_kept() -> None:
    noise = rows("X", X_START + 3_000_000, 5000, [0, 1], REF_DP // 2, filters=["LowQual"])
    found = evidence(frame(chr20(), *XY, noise))
    assert found.inferred_karyotypic_sex is KaryotypicSex.XY
    unapplied = [part.with_columns(pl.lit(None, dtype=pl.List(pl.String)).alias("filters")) for part in (chr20(), *XY)]
    assert evidence(frame(*unapplied)).inferred_karyotypic_sex is KaryotypicSex.XY


def test_split_rows_do_not_count_as_hets() -> None:
    # A `1/2` written one row per ALT reads `0/1` twice; flagged rows stay out of the het share.
    split = rows("X", X_START + 3_000_000, 5000, [0, 1], REF_DP // 2, split="slot_pattern")
    found = evidence(frame(chr20(), *XY, split))
    assert found.x_heterozygosity.diploid_alt_sites == 2000
    assert found.x_copies is XCopies.ONE


def test_allele_balance_reads_ad_at_the_called_alt() -> None:
    x = pl.concat([rows("X", X_START, 1200, [0, 2], REF_DP, ad=[20, 0, 10]), rows("X", X_START + 1_000_000, 800, [1, 1], REF_DP)])
    found = evidence(frame(chr20(), x))
    assert found.x_heterozygosity.median_het_allele_balance == pytest.approx(10 / 30)


@pytest.mark.parametrize(
    ("build", "kind", "reason"),
    [
        (None, CallsetType.WGS, "no_build"),
        ("GRCh36", CallsetType.WGS, "unsupported_build"),
        ("GRCh38", CallsetType.WES, "not_wgs"),
        ("GRCh38", None, "not_wgs"),
    ],
)
def test_no_estimate_without_a_build_or_a_genome(build: str | None, kind: CallsetType | None, reason: str) -> None:
    found = evidence(frame(chr20(), *XX), build, kind)
    assert (found.inferred_karyotypic_sex, found.reason, found.x_copies) == (KaryotypicSex.UNKNOWN_KARYOTYPE, reason, None)
    assert found.ploidy["autosomes"].two == 2000
    assert ("x_nonpar" in found.ploidy) is (build == "GRCh38")
    if reason == "not_wgs":
        assert found.x_heterozygosity.het_share == 0.6


T2T_Y_START = 2_500_000  # inside T2T-CHM13's first X-degenerate interval (2,458,321–2,727,072)


def test_t2t_is_measured_on_its_own_regions_but_not_estimated() -> None:
    """RM286: T2T-CHM13 has region sets, but the thresholds were measured on GRCh37/38 only."""
    assert ESTIMATE_BUILDS < REGION_BUILDS and "T2T-CHM13" in REGION_BUILDS - ESTIMATE_BUILDS
    y = rows("Y", T2T_Y_START, 1000, [1, 1], 18, step=20)
    par = rows("X", 100_000, 500, [0, 1], REF_DP)  # PAR1 on T2T-CHM13 (to 2,394,410)
    found = evidence(frame(chr20(), XY[0], y, par), "T2T-CHM13")
    assert (found.inferred_karyotypic_sex, found.reason, found.x_copies, found.y_present) == (
        KaryotypicSex.UNKNOWN_KARYOTYPE, "uncalibrated_build", None, None,
    )  # fmt: skip
    assert found.region_set == REGION_SET["T2T-CHM13"] and "T2T chm13v2.0 chrXY sequence classes" in found.region_set
    assert found.ploidy["par"].two == 500 and found.ploidy["x_nonpar"].two == 2000
    assert (found.depth.y_sites, found.depth.y_ratio, found.depth.x_ratio) == (1000, 18 / REF_DP, 0.5)
    assert "X-degenerate Y depth ratio 0.600" in found.detail
    # The same evidence on a calibrated build (its Y rows placed in GRCh38's X-degenerate Y) is estimated.
    grch38 = evidence(frame(chr20(), *XY), "GRCh38")
    assert grch38.inferred_karyotypic_sex is KaryotypicSex.XY


def test_reasons_are_a_closed_list_checked_on_the_model() -> None:
    found = evidence(frame(chr20(), *XX), None)
    with pytest.raises(ValueError, match="unknown reason"):
        SexChromosomeEvidence(**{**found.model_dump(), "reason": "nearest_class"})
    assert {"no_build", "y_unusable", "caller_asserted_only", "too_few_sites", "intermediate", "not_wgs"} <= UNKNOWN_REASONS


def test_karyotypic_sex_is_the_phenopackets_vocabulary() -> None:
    assert [member.value for member in KaryotypicSex] == [
        "UNKNOWN_KARYOTYPE", "XX", "XY", "XO", "XXY", "XXX", "XXYY", "XXXY", "XXXX", "XYY", "OTHER_KARYOTYPE",
    ]


# --- Through `load` -----------------------------------------------------------------------------

GRCH38_LENGTHS = {"20": 64_444_167, "X": 156_040_895, "Y": 57_227_415, "1": 248_956_422, "2": 242_193_529}


def write_vcf(path: Path, parts: list[pl.DataFrame], prefix: str = "chr", meta: tuple[str, ...] = ()) -> Path:
    header = [
        "##fileformat=VCFv4.2",
        *meta,
        *[f"##contig=<ID={prefix}{name},length={length}>" for name, length in GRCH38_LENGTHS.items()],
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    body = (
        pl.concat(parts, how="diagonal_relaxed")
        .select(
            (pl.lit(prefix) + pl.col("chrom")).alias("chrom"),
            pl.col("pos").cast(pl.String),
            pl.lit(".").alias("id"),
            pl.lit("A").alias("ref"),
            pl.lit("G,T").alias("alt"),
            pl.lit("50").alias("qual"),
            pl.col("filters").list.join(";").fill_null(".").alias("filter"),
            pl.lit(".").alias("info"),
            pl.lit("GT:DP:AD").alias("format"),
            pl.concat_str(
                pl.col("gt").list.eval(pl.element().cast(pl.String)).list.join("/"),
                pl.col("dp").cast(pl.String).fill_null("."),
                pl.col("ad").list.eval(pl.element().cast(pl.String)).list.join(",").fill_null("."),
                separator=":",
            ).alias("S1"),
        )
        .sort("chrom", pl.col("pos").cast(pl.Int64))
    )
    path.write_text("\n".join(header) + "\n")
    with path.open("a") as handle:
        body.write_csv(handle, separator="\t", include_header=False, quote_style="never")
    return path


def test_through_load_a_small_vcf_is_measured_but_not_estimated(tmp_path: Path) -> None:
    path = write_vcf(tmp_path / "small.vcf", [chr20(), *XY])
    loaded = load(path)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["assembly_label"] == "GRCh38"
    (found,) = sex_chromosome_evidence(loaded)
    assert (found.facet_id, found.reason) == (facet["facet_id"], "not_wgs")
    # `chrX`, `chrY` and `chr20` reach the groups through the loader's canonical names.
    assert found.x_heterozygosity.diploid_alt_sites == 2000
    assert (found.depth.x_ratio, found.depth.y_ratio) == (0.5, 0.6)
    assert found.x_heterozygosity.median_het_allele_balance == 0.5


def test_through_load_a_genome_shaped_vcf_is_estimated(tmp_path: Path) -> None:
    # 1,000,000+ autosomal sites spread over two chromosomes make the facet `wgs` (RM4).
    spread = [rows(chrom, 1, 520_000, [0, 1], REF_DP, step=400) for chrom in ("1", "2")]
    # The header's karyotype line is the caller's statement: carried verbatim, and it does not decide.
    path = write_vcf(tmp_path / "genome.vcf", [*spread, chr20(), XX[0], XY[1]], prefix="", meta=("##referenceSexKaryotype=XY",))
    (found,) = sex_chromosome_evidence(load(path))
    assert found.callset_type is CallsetType.WGS
    assert found.inferred_karyotypic_sex is KaryotypicSex.XXY
    assert found.ploidy_statements == ["##referenceSexKaryotype=XY"]


def test_each_facet_of_a_load_gets_its_own_answer() -> None:
    # Two facets in one `Loaded` (a freeBayes+VarScan2 merge loads so): one XXY, one XY, never blended.
    first, second = frame(chr20(), XX[0], XY[1]), frame(chr20(), *XY).with_columns(pl.lit("g").alias("facet_id"))
    facets = conform(
        pl.LazyFrame({"facet_id": ["f", "g"], "assembly_label": ["GRCh38"] * 2, "callset_type": ["wgs"] * 2, "caller": ["freeBayes", "VarScan2"]}),
        FACETS,
    )
    loaded = Loaded(facets=facets, calls=pl.concat([first, second]), blocks=empty(BLOCKS))
    found = sex_chromosome_evidence(loaded)
    assert [(e.facet_id, e.caller, e.inferred_karyotypic_sex) for e in found] == [
        ("f", "freeBayes", KaryotypicSex.XXY),
        ("g", "VarScan2", KaryotypicSex.XY),
    ]


# --- With a coverage facet (RM278) -----------------------------------------------------------------

COVERAGE_WINDOW = 10_000
X_DEGENERATE_SPAN = (11_674_124, 13_984_473)  # GRCh38's longest X-degenerate interval, 2.3 Mb


def windows(chrom: str, start: int, count: int, depth: float) -> list[tuple]:
    return [
        (chrom, start + i * COVERAGE_WINDOW, start + (i + 1) * COVERAGE_WINDOW - 1, depth, int(depth), int(depth) + 1, COVERAGE_WINDOW, None)
        for i in range(count)
    ]


def coverage_facet(x_depth: float, y_depth: float, *, build: str | None = "GRCh38", chr20_gap: int = 0) -> Coverage:
    """30x on chr20 (with `chr20_gap` zero windows, an assembly gap), `x_depth` over 2 Mb of non-PAR X,
    `y_depth` over 2 Mb inside the X-degenerate Y."""
    rows_ = windows("20", 1_000_001, 200, 30.0) + windows("20", 3_000_001, chr20_gap, 0.0)
    rows_ += windows("X", X_START, 200, x_depth) + windows("Y", X_DEGENERATE_SPAN[0], 200, y_depth)
    intervals = pl.DataFrame(rows_, schema=list(COVERAGE_INTERVALS), orient="row").cast(COVERAGE_INTERVALS).lazy()
    facet = CoverageFacet(
        facet_id="cov", format_id=CoverageFormat.ALIGNMENT_DEPTH, source_path="synthetic", source_sha256=None, assembly_label=build,
        assembly_label_basis=None if build is None else AssemblyLabelBasis.STATED, contig_lengths=None, resolution=Resolution.WINDOW,
        window=COVERAGE_WINDOW, options={}, zero_filled=True, dropped_fields=[], warnings=[],
    )
    return Coverage(facet=facet, intervals=intervals)


def with_coverage(calls: pl.LazyFrame, coverage: Coverage, build: str = "GRCh38") -> SexChromosomeEvidence:
    facets = conform(pl.LazyFrame({"facet_id": ["f"], "assembly_label": [build], "callset_type": ["wgs"]}), FACETS)
    (found,) = sex_chromosome_evidence(Loaded(facets=facets, calls=calls, blocks=empty(BLOCKS)), coverage=coverage)
    return found


@pytest.mark.parametrize(
    ("x", "y", "x_depth", "y_depth", "expected"),
    [
        # Coverage shows a Y absent, which call-site depth never could: XX and XO are named.
        (XX[0], XX[1], 30.0, 0.6, KaryotypicSex.XX),
        (XY[0], XX[1], 15.0, 0.0, KaryotypicSex.XO),
        # Coverage shows a Y the calls did not: XY.
        (XY[0], XX[1], 15.0, 15.0, KaryotypicSex.XY),
        (XX[0], XY[1], 30.0, 15.0, KaryotypicSex.XXY),
    ],
)
def test_coverage_names_xx_and_xo(x: pl.DataFrame, y: pl.DataFrame, x_depth: float, y_depth: float, expected: KaryotypicSex) -> None:
    found = with_coverage(frame(chr20(), x, y), coverage_facet(x_depth, y_depth))
    assert (found.inferred_karyotypic_sex, found.reason) == (expected, None)
    assert found.y_present is (expected in (KaryotypicSex.XY, KaryotypicSex.XXY))
    assert found.coverage.reason is None
    assert (found.coverage.x_ratio, found.coverage.y_ratio) == (pytest.approx(x_depth / 30), pytest.approx(y_depth / 30))


def test_coverage_absent_against_a_call_depth_present_is_unknown() -> None:
    found = with_coverage(frame(chr20(), *XY), coverage_facet(15.0, 0.0))
    assert (found.inferred_karyotypic_sex, found.reason, found.y_reason) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "y_unusable", "lines_disagree")


def test_a_coverage_y_between_the_bands_is_intermediate() -> None:
    found = with_coverage(frame(chr20(), XX[0], XX[1]), coverage_facet(30.0, 7.5))
    assert (found.inferred_karyotypic_sex, found.y_present, found.y_reason) == (KaryotypicSex.UNKNOWN_KARYOTYPE, None, "intermediate")


def test_a_coverage_x_that_disagrees_with_the_het_share_is_unknown() -> None:
    found = with_coverage(frame(chr20(), XX[0], XX[1]), coverage_facet(15.0, 0.0))
    assert (found.x_copies, found.x_copies_reason, found.reason) == (None, "lines_disagree", "lines_disagree")


def test_coverage_alone_decides_y_where_the_calls_carry_no_depth() -> None:
    # A liftover that dropped FORMAT/DP: the het share gives two X, coverage shows no Y.
    found = with_coverage(frame(chr20(dp=None), x_calls(1200, 800, None), y_calls(1000, None)), coverage_facet(30.0, 0.0))
    assert found.inferred_karyotypic_sex is KaryotypicSex.XX


def test_assembly_gaps_on_chr20_do_not_pull_the_reference_down() -> None:
    # 300 zero windows against 200 at 30x: counted, the median would be zero and no ratio measured.
    found = with_coverage(frame(chr20(), XX[0], XX[1]), coverage_facet(30.0, 0.0, chr20_gap=300))
    assert found.coverage.reference_median == 30.0
    assert found.inferred_karyotypic_sex is KaryotypicSex.XX


@pytest.mark.parametrize(("coverage_build", "reason"), [("GRCh37", "coverage_build_mismatch"), (None, "coverage_no_build")])
def test_coverage_on_another_or_no_build_is_set_aside(coverage_build: str | None, reason: str) -> None:
    found = with_coverage(frame(chr20(), XX[0], XX[1]), coverage_facet(30.0, 0.0, build=coverage_build))
    assert found.coverage.reason == reason
    assert (found.coverage.x_ratio, found.coverage.y_ratio) == (None, None)
    # Without its line the answer is the calls' alone: two X, Y not shown.
    assert (found.inferred_karyotypic_sex, found.x_copies, found.y_reason) == (KaryotypicSex.UNKNOWN_KARYOTYPE, XCopies.TWO_OR_MORE, "y_not_shown")


def test_uneven_coverage_is_set_aside() -> None:
    # An exome's shape: chr20 deep on target, shallow off it. Its medians say nothing about copies.
    even = coverage_facet(30.0, 0.0)
    exome_like = pl.concat([even.intervals, pl.DataFrame(windows("20", 5_000_001, 300, 2.0), schema=list(COVERAGE_INTERVALS), orient="row").cast(COVERAGE_INTERVALS).lazy()])
    found = with_coverage(frame(chr20(), XX[0], XX[1]), Coverage(facet=even.facet, intervals=exome_like))
    assert found.coverage.reason == "coverage_uneven"
    assert found.coverage.reference_even_share == pytest.approx(0.6)  # the median falls on the shallow 300
    assert (found.coverage.x_ratio, found.coverage.y_ratio, found.inferred_karyotypic_sex) == (None, None, KaryotypicSex.UNKNOWN_KARYOTYPE)
    assert with_coverage(frame(chr20(), XX[0], XX[1]), even).coverage.reference_even_share == 1.0


def test_coverage_with_too_few_bases_gives_no_ratio() -> None:
    full = coverage_facet(30.0, 0.0)
    kept = (pl.col("chrom") != "Y") | (pl.col("pos") < X_DEGENERATE_SPAN[0] + 50 * COVERAGE_WINDOW)
    found = with_coverage(frame(chr20(), XX[0], XX[1]), Coverage(facet=full.facet, intervals=full.intervals.filter(kept)))
    assert (found.coverage.y_bases, found.coverage.y_ratio) == (500_000, None)
    assert found.y_reason == "y_not_shown"


# --- Counts past two, from coverage only ---------------------------------------------------------

THREE_X_DP = round(REF_DP * 1.5)  # calls' X depth 1.5: above two copies' band, three or more


@pytest.mark.parametrize(
    ("x", "y", "x_depth", "y_depth", "expected", "counts"),
    [
        # Copy number: X at 1.5 of chr20 is three copies, Y at 1.0 two.
        (x_calls(1200, 800, THREE_X_DP), XX[1], 45.0, 0.0, KaryotypicSex.XXX, (3, 0)),
        (XY[0], XY[1], 15.0, 30.0, KaryotypicSex.XYY, (1, 2)),
        (XX[0], XY[1], 30.0, 30.0, KaryotypicSex.XXYY, (2, 2)),
        (x_calls(1200, 800, THREE_X_DP), XY[1], 45.0, 15.0, KaryotypicSex.XXXY, (3, 1)),
        # Phenopackets has no XXXYY: the standard's own OTHER_KARYOTYPE.
        (x_calls(1200, 800, THREE_X_DP), XY[1], 45.0, 30.0, KaryotypicSex.OTHER_KARYOTYPE, (3, 2)),
        (XY[0], XX[1], 15.0, 15.0, KaryotypicSex.XY, (1, 1)),
    ],
)
def test_coverage_counts_three_x_and_two_y(
    x: pl.DataFrame, y: pl.DataFrame, x_depth: float, y_depth: float, expected: KaryotypicSex, counts: tuple[int, int]
) -> None:
    found = with_coverage(frame(chr20(), x, y), coverage_facet(x_depth, y_depth))
    assert (found.inferred_karyotypic_sex, found.reason) == (expected, None)
    assert (found.x_copy_count, found.y_copy_count) == counts
    assert found.x_copies is (XCopies.ONE if counts[0] == 1 else XCopies.TWO_OR_MORE)
    assert found.y_present is (counts[1] > 0)


@pytest.mark.parametrize(
    ("x", "x_depth", "y_depth", "x_reason", "y_reason"),
    [
        # Between two counts' bands: the central half of each gap (0.5 per copy, 0.125 either side).
        (XX[0], 37.5, 0.0, "intermediate", None),  # X 1.25, between two and three
        (XX[0], 30.0, 22.5, None, "intermediate"),  # Y 0.75, between one and two
        (XY[0], 7.5, 15.0, "intermediate", None),  # X 0.25, under one copy's band
        (XX[0], 30.0, 37.5, None, "intermediate"),  # Y 1.25: three Y are not read
    ],
)
def test_a_coverage_count_between_bands_is_unknown(x: pl.DataFrame, x_depth: float, y_depth: float, x_reason: str | None, y_reason: str | None) -> None:
    found = with_coverage(frame(chr20(), x, XX[1]), coverage_facet(x_depth, y_depth))
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE
    assert (found.x_copies_reason, found.y_reason) == (x_reason, y_reason)
    assert found.reason == (x_reason or "y_unusable")


def test_a_coverage_three_against_a_calls_two_disagrees() -> None:
    # The calls' X depth 1.0 is two copies; coverage at 1.5 says three. Neither wins.
    found = with_coverage(frame(chr20(), *XX), coverage_facet(45.0, 0.0))
    assert (found.x_copies, found.x_copies_reason, found.inferred_karyotypic_sex) == (None, "lines_disagree", KaryotypicSex.UNKNOWN_KARYOTYPE)


def test_the_coverage_bands_are_copy_number() -> None:
    # Each count k at k/2 of chr20's depth, 0.125 either side: the edges hold, the next step out does not.
    for count, (low, high, copies) in zip((1, 2, 3), COVERAGE_X_BANDS, strict=True):
        assert (low, high, copies) == (count / 2 - 0.125, count / 2 + 0.125, (count, count))
    assert [band[:2] for band in COVERAGE_Y_BANDS] == [(None, 0.125), (0.375, 0.625), (0.875, 1.125)]
    edge = with_coverage(frame(chr20(), x_calls(1200, 800, THREE_X_DP), XX[1]), coverage_facet(41.25, 3.75))  # 1.375, 0.125
    assert (edge.x_copy_count, edge.y_copy_count) == (3, 0)


# --- Array exports: letters as written, no REF ---------------------------------------------------


def letters(chrom: str, start: int, n: int, alleles: list[str | None], *, step: int = 10) -> pl.DataFrame:
    """Array rows as a loader writes them before `resolve_ref`: letters in `gt_alleles`, no `gt`."""
    return pl.DataFrame(
        {
            "facet_id": "f",
            "chrom": chrom,
            "pos": pl.int_range(0, n, eager=True) * step + start,
            "gt_alleles": [alleles] * n,
            "evidence_basis": "assay",
        },
        schema_overrides={"facet_id": pl.String, "chrom": pl.String, "pos": pl.UInt32, "gt_alleles": pl.List(pl.String), "evidence_basis": pl.String},
    )


def autosomes(hets: int = 600, homs: int = 1400, no_calls: int = 0) -> pl.DataFrame:
    """Autosomal het share 0.3 by default (600 of 2,000 two-letter calls)."""
    return pl.concat([
        letters("1", 1_000_000, hets, ["A", "G"]), letters("2", 1_000_000, homs, ["C", "C"]),
        letters("3", 1_000_000, no_calls, [None, None]),
    ])  # fmt: skip


def array_x(hets: int, homs: int, *, one_letter: int = 0) -> pl.DataFrame:
    return pl.concat([
        letters("X", X_START, hets, ["A", "G"]), letters("X", X_START + 1_000_000, homs, ["T", "T"]),
        letters("X", X_START + 2_000_000, one_letter, ["T"]),
    ])  # fmt: skip


def array_y(called: int, no_calls: int) -> pl.DataFrame:
    return pl.concat([letters("Y", Y_START, called, ["G"], step=20), letters("Y", Y_START + 20 * called, no_calls, [None], step=20)])


def on_array(*parts: pl.DataFrame, read_y: bool = False, build: str = "GRCh38") -> SexChromosomeEvidence:
    return estimate(measure(frame(*parts), build), "f", build, CallsetType.ARRAY, array_y=read_y)


def test_an_array_with_two_letters_on_x_reads_two_x_from_its_hets() -> None:
    # X het share 500/2,000 = 0.25 over the autosomes' 0.3: ratio 0.833, two X. Y no-called throughout.
    found = on_array(autosomes(), array_x(500, 1500), array_y(0, 700))
    lines = found.array
    assert (lines.x_base_pairs, lines.x_hets, lines.autosome_base_pairs, lines.autosome_hets) == (2000, 500, 2000, 600)
    assert lines.x_het_ratio == pytest.approx((500 / 2000) / (600 / 2000))
    assert (lines.y_rows, lines.y_called, lines.y_call_ratio) == (700, 0, 0.0)
    assert (found.x_copies, found.x_copy_count) == (XCopies.TWO_OR_MORE, None)
    # The vendor's Y no-calls follow its own sex call: by default they are not read.
    assert (found.y_present, found.y_reason, found.reason) == (None, "vendor_asserted", "y_unusable")
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE
    read = on_array(autosomes(), array_x(500, 1500), array_y(0, 700), read_y=True)
    assert (read.y_present, read.y_copy_count, read.inferred_karyotypic_sex) == (False, 0, KaryotypicSex.XX)


def test_an_arrays_one_letter_x_is_the_vendors_assertion() -> None:
    # 23andMe writes a one-X sample's X with one letter: its call, not evidence. Y called (one letter).
    found = on_array(autosomes(), array_x(0, 0, one_letter=2000), array_y(650, 50), read_y=True)
    assert found.array.x_called_one_letter == 2000 and found.array.x_het_ratio is None
    assert (found.x_copies, found.x_copies_reason) == (None, "caller_asserted_only")
    assert found.y_present is True
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE


def test_an_array_writing_one_x_with_two_letters_reads_one_x() -> None:
    # Doubled letters with a genotyping error's few hets: ratio (4/2000)/0.3, about 0.007.
    found = on_array(autosomes(), array_x(4, 1996), array_y(650, 50), read_y=True)
    assert (found.x_copies, found.x_copy_count) == (XCopies.ONE, 1)
    assert found.array.y_call_ratio == pytest.approx(650 / 700)
    assert (found.inferred_karyotypic_sex, found.reason) == (KaryotypicSex.XY, None)


@pytest.mark.parametrize(
    ("x", "y", "x_reason", "y_reason"),
    [
        (array_x(300, 1700), array_y(0, 700), "intermediate", None),  # ratio 0.5, between one and two X
        (array_x(500, 1500), array_y(350, 350), None, "intermediate"),  # Y called half as often as the autosomes
        (array_x(500, 1500), array_y(100, 100), None, "too_few_sites"),  # 200 X-degenerate Y rows
    ],
)
def test_weak_array_lines_are_unknown(x: pl.DataFrame, y: pl.DataFrame, x_reason: str | None, y_reason: str | None) -> None:
    found = on_array(autosomes(), x, y, read_y=True)
    assert (found.x_copies_reason, found.y_reason) == (x_reason, y_reason)
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE


def test_array_rows_without_two_bases_are_not_het_evidence() -> None:
    # 23andMe's I/D letters and written no-calls are neither base pairs nor hets; no-calls lower the
    # autosomes' called share, which the Y ratio divides by.
    indels = letters("X", X_START + 3_000_000, 3000, ["I", "D"])
    found = on_array(autosomes(no_calls=500), array_x(500, 1500), indels, array_y(560, 140), read_y=True)
    assert (found.array.x_base_pairs, found.array.x_hets) == (2000, 500)
    assert (found.array.autosome_called, found.array.autosome_rows) == (2000, 2500)
    assert found.array.y_call_ratio == pytest.approx((560 / 700) / (2000 / 2500))
    assert found.inferred_karyotypic_sex is KaryotypicSex.XXY


def test_a_genome_carries_no_array_lines() -> None:
    assert evidence(frame(chr20(), *XY)).array is None
