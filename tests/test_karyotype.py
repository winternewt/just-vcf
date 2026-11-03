"""Sex-chromosome evidence and the karyotype estimate it supports (RM269)."""

import random
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
    Person,
    PersonKaryotype,
    ReferenceStore,
    Resolution,
    XCopies,
    load,
    person_karyotype,
    person_sex_chromosome_evidence,
    resolve_ref,
    sex_chromosome_evidence,
)
from just_vcf.coverage import COVERAGE_INTERVALS
from just_vcf.formats import REGISTRY
from just_vcf.karyotype import (
    ARRAY_ESTIMATE_BUILDS,
    ARRAY_ONE_X_MAX_HET_RATIO,
    ARRAY_TWO_X_MIN_HET_RATIO,
    ARRAY_Y_ABSENT_MAX_CALL_RATIO,
    ARRAY_Y_PRESENT_MIN_CALL_RATIO,
    ARRAY_Y_VENDORS,
    ASSEMBLY_LABEL_SOURCES,
    CALLSET_BASES,
    COVERAGE_PAIRINGS,
    COVERAGE_X_BANDS,
    COVERAGE_Y_BANDS,
    MIN_X_SITES,
    MIN_Y_SITES,
    ESTIMATE_BUILDS,
    ONE_X_MAX_HET_SHARE,
    PAR_BOUNDS,
    PERSON_KARYOTYPE_REASONS,
    REGION_BUILDS,
    REGION_SET,
    SET_ASIDE_REASONS,
    TWO_X_MAX_DEPTH_RATIO,
    TWO_X_MIN_HET_SHARE,
    UNKNOWN_REASONS,
    X_DEGENERATE_Y,
    CallsetKaryotype,
    CoveragePairing,
    SetAsideFacet,
    SexChromosomeEvidence,
    assembly_label_source,
    estimate,
    measure,
)
from just_vcf.liftover import read_chain
from just_vcf.links import (
    CALLER_SPLIT,
    CONSENSUS,
    LIFTOVER_TO,
    MERGE,
    NORMALIZED,
    derived_links,
    kin_links,
    link_disagreements,
    links_frame,
)
from just_vcf.schema import CALLS, AssemblyLabelBasis, CallsetType, conform, empty
from test_reference import write_fasta

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


@pytest.mark.parametrize("build", sorted(GIAB_PREFIX))
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


# --- GRCh36 (NCBI36): Ensembl 54's PARs, and GRCh37's X-degenerate Y lifted through UCSC's chain ---

ENSEMBL_54 = ASSETS / "ensembl_54_ncbi36"
HG19_TO_HG18_Y = ASSETS / "ucsc_liftover" / "hg19ToHg18.over.chain.chrY.chain"


def ensembl_table(name: str) -> list[list[str]]:
    return [line.split("\t") for line in (ENSEMBL_54 / name).read_text().splitlines() if line]


def test_grch36_par_bounds_are_ensembl_54s_par_exceptions() -> None:
    """Ensembl stores each PAR once, on X, and Y's copy as an `assembly_exception` of type PAR pointing at
    it. Coordinates are 1-based closed, as `pos` is; the seq_regions are NCBI36 chromosomes."""
    (ncbi36,) = [row[0] for row in ensembl_table("coord_system.txt") if row[2:4] == ["chromosome", "NCBI36"]]
    names = {row[0]: row[1] for row in ensembl_table("seq_region_XY.txt") if row[2] == ncbi36}
    pars = sorted(
        ((int(y_start), int(y_end)), (int(x_start), int(x_end)))
        for _, y, y_start, y_end, kind, x, x_start, x_end, _ in ensembl_table("assembly_exception.txt")
        if kind == "PAR" and (names.get(y), names.get(x)) == ("Y", "X")
    )
    assert PAR_BOUNDS["GRCh36"] == {"X": tuple(x for _, x in pars), "Y": tuple(y for y, _ in pars)}
    # Each PAR is one sequence on both chromosomes: the same length on X and Y.
    assert [e - s for s, e in PAR_BOUNDS["GRCh36"]["X"]] == [e - s for s, e in PAR_BOUNDS["GRCh36"]["Y"]]


def lifted_to_grch36(chain_blocks: pl.DataFrame, intervals: tuple[tuple[int, int], ...]) -> tuple[list[tuple[int, int]], int, int]:
    """Each 1-based closed GRCh37 Y interval through the chain's blocks: the span from its first to its
    last base landing on chrY. Also the bases landing elsewhere and the bases in no block."""
    blocks = chain_blocks.filter(pl.col("t_chrom") == "Y").sort("t_start")
    spans, elsewhere, unaligned = [], 0, 0
    for start, end in intervals:
        hit = blocks.filter((pl.col("t_end") > start - 1) & (pl.col("t_start") < end))
        pieces = []
        for block in hit.iter_rows(named=True):
            assert not block["q_minus"]
            t_start, t_end = max(start - 1, block["t_start"]), min(end, block["t_end"])
            q_start = block["q_start"] + t_start - block["t_start"]
            pieces.append((block["q_chrom"], block["chain"], q_start + 1, q_start + t_end - t_start))
        on_y = [piece for piece in pieces if piece[0] == "chrY"]
        elsewhere += sum(e - s + 1 for chrom, _, s, e in pieces if chrom != "chrY")
        unaligned += end - start + 1 - sum(e - s + 1 for _, _, s, e in pieces)
        assert len({chain for _, chain, _, _ in on_y}) == 1  # one chain, in order: a span means something
        assert all(a[3] < b[2] for a, b in zip(on_y, on_y[1:]))
        spans.append((on_y[0][2], on_y[-1][3]))
    return spans, elsewhere, unaligned


def test_grch36_x_degenerate_y_is_grch37s_lifted_through_hg19_to_hg18() -> None:
    chain = read_chain(HG19_TO_HG18_Y)
    spans, elsewhere, unaligned = lifted_to_grch36(chain.blocks, X_DEGENERATE_Y["GRCh37"])
    (_, par1_end), (par2_start, _) = PAR_BOUNDS["GRCh36"]["Y"]
    kept = tuple((max(s, par1_end + 1), min(e, par2_start - 1)) for s, e in spans if s < par2_start and e > par1_end)
    assert X_DEGENERATE_Y["GRCh36"] == kept
    # One interval to every GRCh37 one; only the last crossed into PAR2 (NCBI's, sequence-checked).
    assert len(kept) == len(X_DEGENERATE_Y["GRCh37"])
    assert [span != k for span, k in zip(spans, kept)] == [False] * (len(kept) - 1) + [True]
    assert (elsewhere, unaligned) == (117, 536)  # 117 bases go to hg18's chrX; 536 are in no block
    assert all(a[1] < b[0] for a, b in zip(kept, kept[1:]))


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
        ("NCBI35", CallsetType.WGS, "unsupported_build"),  # hg17: a real build with no region set here
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


GRCH36_Y_START = 11_220_000  # in GRCh36's third X-degenerate interval (11,214,554–11,253,954); not one on GRCh37


def test_grch36_is_measured_on_its_own_regions_but_not_estimated_from_calls() -> None:
    """GRCh36 has region sets (Ensembl 54's PARs, GRCh37's X-degenerate Y lifted), but no NCBI36 genome
    was among the calibration files: a whole genome on it is measured and reported, not estimated."""
    assert "GRCh36" in REGION_BUILDS - ESTIMATE_BUILDS
    y = rows("Y", GRCH36_Y_START, 1000, [1, 1], 18, step=20)
    par1 = rows("X", 2_700_001, 500, [0, 1], REF_DP)  # PAR1 on GRCh36 (to 2,709,520), not on GRCh37 (to 2,699,520)
    par2 = rows("Y", 57_450_000, 300, [0, 1], REF_DP)  # Y's PAR2 on GRCh36 (from 57,443,438)
    parts = (chr20(), XY[0], y, par1, par2)
    found = evidence(frame(*parts), "GRCh36")
    assert (found.inferred_karyotypic_sex, found.reason, found.x_copies, found.y_present) == (
        KaryotypicSex.UNKNOWN_KARYOTYPE, "uncalibrated_build", None, None,
    )  # fmt: skip
    assert found.region_set == REGION_SET["GRCh36"] and "hg19ToHg18" in found.region_set
    assert (found.ploidy["par"].two, found.ploidy["x_nonpar"].two, found.ploidy["y_nonpar"].two) == (800, 2000, 1000)
    assert (found.depth.y_sites, found.depth.y_ratio, found.depth.x_ratio) == (1000, 18 / REF_DP, 0.5)
    # On GRCh37's regions the same rows place otherwise: PAR1 rows are non-PAR X, the Y rows no X-degenerate Y.
    grch37 = evidence(frame(*parts), "GRCh37")
    assert (grch37.ploidy["x_nonpar"].two, grch37.depth.y_sites) == (2500, 0)


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


@pytest.mark.parametrize(("build", "reason"), [("GRCh36", "uncalibrated_build"), ("T2T-CHM13", "uncalibrated_build"), ("NCBI35", "unsupported_build")])
def test_coverage_on_a_build_without_thresholds_says_whether_it_has_regions(build: str, reason: str) -> None:
    found = with_coverage(frame(chr20(), XX[0], XX[1]), coverage_facet(30.0, 0.0, build=build), build=build)
    assert (found.coverage.reason, found.coverage.x_ratio, found.coverage.y_ratio) == (reason, None, None)
    assert found.reason == reason


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


AUTOSOME_CALL_RATE = 2000 / 2040


def autosomes(hets: int = 600, homs: int = 1400, no_calls: int = 40) -> pl.DataFrame:
    """Autosomal het share 0.3 by default (600 of 2,000 two-letter calls), 40 written no-calls."""
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
    assert (lines.autosome_rows, lines.autosome_called) == (2040, 2000)
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
    assert found.array.y_call_ratio == pytest.approx((650 / 700) / AUTOSOME_CALL_RATE)
    assert (found.inferred_karyotypic_sex, found.reason) == (KaryotypicSex.XY, None)


def test_an_array_that_drops_its_no_calls_gives_no_y_line() -> None:
    # Without no-call rows the Y rows are only the called ones: 700 of 700 says nothing about a Y.
    found = on_array(autosomes(no_calls=0), array_x(500, 1500), array_y(700, 0), read_y=True)
    assert (found.array.autosome_called, found.array.autosome_rows, found.array.y_call_ratio) == (2000, 2000, None)
    assert (found.y_present, found.y_reason, found.inferred_karyotypic_sex) == (None, "no_calls_not_written", KaryotypicSex.UNKNOWN_KARYOTYPE)
    assert found.x_copies is XCopies.TWO_OR_MORE


def test_an_unreadable_autosomal_row_is_not_a_written_no_call() -> None:
    # A row whose genotype no loader could read (`gt_alleles` null) is not called, but it is no no-call
    # the export wrote: one such row must not make a no-call-dropping export's Y rows read as its probes.
    unreadable = letters("4", 1_000_000, 1, None)
    found = on_array(autosomes(no_calls=0), unreadable, array_x(500, 1500), array_y(700, 0), read_y=True)
    assert (found.array.y_call_ratio, found.y_present, found.y_reason) == (None, None, "no_calls_not_written")
    assert (found.array.autosome_rows, found.array.autosome_called, found.array.autosome_no_calls) == (2001, 2000, 0)
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE
    # With the export's own no-calls written beside it, the Y line reads.
    written = on_array(autosomes(no_calls=40), unreadable, array_x(500, 1500), array_y(700, 0), read_y=True)
    assert written.array.autosome_no_calls == 40 and written.y_present is True


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


def test_an_array_on_grch36_is_estimated_on_grch36s_regions() -> None:
    """A GRCh36 export (23andMe's v1 to v3) is read as a GRCh37 one is, on NCBI36's PARs and X-degenerate Y."""
    assert "GRCh36" in ARRAY_ESTIMATE_BUILDS - ESTIMATE_BUILDS
    y = pl.concat([letters("Y", GRCH36_Y_START, 650, ["G"], step=20), letters("Y", GRCH36_Y_START + 20 * 650, 50, [None], step=20)])
    par1 = letters("X", 2_700_001, 500, ["A", "A"])  # PAR1 on GRCh36 (to 2,709,520): no non-PAR X call
    parts = (autosomes(), array_x(4, 1996), par1, y)
    found = on_array(*parts, read_y=True, build="GRCh36")
    assert (found.array.x_base_pairs, found.array.y_rows, found.array.y_called) == (2000, 700, 650)
    assert found.array.y_call_ratio == pytest.approx((650 / 700) / AUTOSOME_CALL_RATE)
    assert (found.inferred_karyotypic_sex, found.reason, found.x_copy_count) == (KaryotypicSex.XY, None, 1)
    assert found.region_set == REGION_SET["GRCh36"]
    off = on_array(*parts, build="GRCh36")
    assert (off.x_copies, off.y_reason, off.reason) == (XCopies.ONE, "vendor_asserted", "y_unusable")
    # On GRCh37's regions the same rows find no X-degenerate Y probe, and the PAR1 rows are non-PAR X.
    grch37 = on_array(*parts, read_y=True, build="GRCh37")
    assert (grch37.array.y_rows, grch37.array.x_base_pairs, grch37.y_reason) == (0, 2500, "too_few_sites")


# --- Whether an array's Y is read: per vendor by default (desk #136, 2026-10-26) -----------------

TWO_X_Y_UNCALLED = (autosomes(), array_x(500, 1500), array_y(0, 700))
"""Hets on X at 0.833 of the autosomes' share, every X-degenerate Y probe no-called."""
ONE_X_Y_CALLED = (autosomes(), array_x(4, 1996), array_y(650, 50))
"""A genotyping error's few hets on X, the Y probes called like the autosomes."""
ARRAY_VENDORS = ("AncestryDNA", "23andMe", "FamilyTreeDNA", None)
"""Two vendors checked against self-reported sex on PGP-HMS, one never checked, and none stated."""


def by_the_bands(found: SexChromosomeEvidence) -> KaryotypicSex:
    """The karyotype the array bands give the result's own ratios."""
    x, y = found.array.x_het_ratio, found.array.y_call_ratio
    x_count = 1 if x <= ARRAY_ONE_X_MAX_HET_RATIO else 2 if x >= ARRAY_TWO_X_MIN_HET_RATIO else None
    y_count = 0 if y <= ARRAY_Y_ABSENT_MAX_CALL_RATIO else 1 if y >= ARRAY_Y_PRESENT_MIN_CALL_RATIO else None
    return {(2, 0): KaryotypicSex.XX, (1, 1): KaryotypicSex.XY}[(x_count, y_count)]


def test_the_vendors_whose_y_is_read_by_default_are_ancestrydnas_alone() -> None:
    # The person's decision: on for AncestryDNA, off for 23andMe and every other vendor; spelled as the
    # format registry writes `facets.vendor`, or the default would never match a loaded facet.
    assert ARRAY_Y_VENDORS == frozenset({"AncestryDNA"})
    assert ARRAY_Y_VENDORS <= {spec.vendor for spec in REGISTRY}


@pytest.mark.parametrize("parts", [TWO_X_Y_UNCALLED, ONE_X_Y_CALLED], ids=["two_x_y_uncalled", "one_x_y_called"])
@pytest.mark.parametrize("array_y", [None, True, False])
@pytest.mark.parametrize("vendor", ARRAY_VENDORS)
def test_an_arrays_y_is_read_by_its_vendors_default_unless_the_caller_says(
    parts: tuple[pl.DataFrame, ...], array_y: bool | None, vendor: str | None
) -> None:
    found = estimate(measure(frame(*parts), "GRCh38"), "f", "GRCh38", CallsetType.ARRAY, array_y=array_y, vendor=vendor)
    assert found.x_copies is not None, "X is read from the hets whatever the Y"
    read = array_y if array_y is not None else vendor == "AncestryDNA"
    if read:
        assert (found.y_reason, found.reason) == (None, None)
        assert found.inferred_karyotypic_sex is by_the_bands(found)
        assert ("array_y=True" if array_y else "on by default for AncestryDNA") in found.detail
        return
    assert (found.y_present, found.y_reason, found.reason) == (None, "vendor_asserted", "y_unusable")
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE
    if array_y is False:
        assert "array_y=False" in found.detail and "by default" not in found.detail
    elif vendor is None:
        assert "off by default" in found.detail and "vendor is unknown" in found.detail
    else:
        assert f"off by default for {vendor}" in found.detail


def array_facets(vendors: dict[str, str | None], *, subject: str | None = None) -> Loaded:
    """One array facet per id, each with `vendor` as given and the same two-X letters."""
    columns = {
        "facet_id": list(vendors), "vendor": list(vendors.values()), "assembly_label": ["GRCh38"] * len(vendors),
        "callset_type": ["array"] * len(vendors), "source_sha256": [facet_id * 64 for facet_id in vendors],
    }  # fmt: skip
    if subject is not None:
        columns["subject"] = [subject] * len(vendors)
    calls = pl.concat([frame(*TWO_X_Y_UNCALLED).with_columns(facet_id=pl.lit(facet_id)) for facet_id in vendors])
    return Loaded(facets=conform(pl.LazyFrame(columns), FACETS), calls=calls, blocks=empty(BLOCKS))


def test_each_array_facet_of_a_load_takes_its_own_vendors_default() -> None:
    loaded = array_facets({"a": "AncestryDNA", "t": "23andMe", "n": None})
    default = by_facet(sex_chromosome_evidence(loaded))
    assert default["a"].inferred_karyotypic_sex is by_the_bands(default["a"]) is KaryotypicSex.XX
    assert [default[f].y_reason for f in ("t", "n")] == ["vendor_asserted"] * 2
    # An explicit value applies to every array facet, whatever its vendor.
    assert {f.inferred_karyotypic_sex for f in sex_chromosome_evidence(loaded, array_y=True)} == {KaryotypicSex.XX}
    assert {f.y_reason for f in sex_chromosome_evidence(loaded, array_y=False)} == {"vendor_asserted"}


def test_a_person_with_an_ancestrydna_and_a_23andme_facet_reads_y_on_the_first_only() -> None:
    ancestry, twenty_three = array_facets({"a": "AncestryDNA"}, subject="p"), array_facets({"t": "23andMe"}, subject="p")
    found = by_facet(person_sex_chromosome_evidence(Person.of(ancestry, twenty_three, check=None)))
    assert (found["a"].y_present, found["a"].inferred_karyotypic_sex) == (False, by_the_bands(found["a"]))
    assert (found["t"].y_present, found["t"].y_reason, found["t"].inferred_karyotypic_sex) == (
        None, "vendor_asserted", KaryotypicSex.UNKNOWN_KARYOTYPE,
    )  # fmt: skip
    assert found["a"].x_copies is found["t"].x_copies is XCopies.TWO_OR_MORE


# --- The build: whatever in the library filled `assembly_label`, and which source that was ----------

NON_PAR_X_GRCH37 = 2_700_001  # just past GRCh37's PAR1 (ends 2,699,520)


def drawn(rng: random.Random, chrom: str, seq: str, start: int, n: int, het_share: float) -> list[tuple[str, str, int, str]]:
    """FTDNA rows around `seq`'s bases, every 10th base from `start`: hets at `het_share`, the rest hom-ref."""
    out = []
    for i in range(n):
        pos = start + 10 * i
        ref = seq[pos - 1]
        alt = rng.choice([b for b in "ACGT" if b != ref])
        out.append((f"rs{chrom}{pos}", chrom, pos, "".join(sorted(ref + alt)) if rng.random() < het_share else ref * 2))
    return out


def test_an_array_table_without_a_build_is_estimated_on_the_build_resolve_ref_inferred(tmp_path: Path) -> None:
    """An FTDNA table states no build: as loaded it is `no_build`; `resolve_ref` infers GRCh37 from its rows
    (F9), and the estimate reads that build and says it came from the rows."""
    rng = random.Random(5)
    lengths = {"1": 25_000, "X": NON_PAR_X_GRCH37 + 25_000}
    grch37 = {name: "".join(rng.choices("ACGT", k=length)) for name, length in lengths.items()}
    other = {name: "".join(rng.choices("ACGT", k=length)) for name, length in lengths.items()}
    store = ReferenceStore.open(tmp_path / "store")
    digest = store.add_fasta(write_fasta(tmp_path / "a.fa", grch37), assembly_label="GRCh37")
    store.add_fasta(write_fasta(tmp_path / "b.fa", other), assembly_label="GRCh38")
    table = drawn(rng, "1", grch37["1"], 1, 2000, 0.3) + drawn(rng, "X", grch37["X"], NON_PAR_X_GRCH37, 2000, 0.25)
    path = tmp_path / "ftdna.csv"
    path.write_text("RSID,CHROMOSOME,POSITION,RESULT\n" + "".join(f'"{r}","{c}","{p}","{g}"\n' for r, c, p, g in table))
    loaded = load(path)

    (as_loaded,) = sex_chromosome_evidence(loaded)
    assert (as_loaded.assembly_label, as_loaded.assembly_label_basis, as_loaded.assembly_label_source) == (None, None, None)
    assert as_loaded.reason == "no_build" and "resolve_ref" in as_loaded.detail

    resolved = resolve_ref(loaded, store)
    (found,) = sex_chromosome_evidence(resolved)
    assert (found.assembly_label, found.assembly_label_basis, found.assembly_label_source) == ("GRCh37", AssemblyLabelBasis.INFERRED, "rows")
    assert found.ref_seqcol_digest == digest
    # The letters are read as written: the same lines as on the loaded rows placed on GRCh37 by hand.
    by_hand = estimate(measure(loaded.calls, "GRCh37"), found.facet_id, "GRCh37", CallsetType.ARRAY)
    assert found.array == by_hand.array
    hets = {chrom: sum(g[0] != g[1] for _, c, _, g in table if c == chrom) for chrom in ("1", "X")}
    assert found.array.x_het_ratio == pytest.approx((hets["X"] / 2000) / (hets["1"] / 2000))
    assert found.x_copies is XCopies.TWO_OR_MORE  # ratio about 0.83
    assert (found.y_reason, found.inferred_karyotypic_sex) == ("too_few_sites", KaryotypicSex.UNKNOWN_KARYOTYPE)


def test_a_vcfs_build_says_whether_the_header_stated_it_or_its_contig_lengths_gave_it(tmp_path: Path) -> None:
    parts = [chr20(), *XY]
    (from_lengths,) = sex_chromosome_evidence(load(write_vcf(tmp_path / "lengths.vcf", parts)))
    assert (from_lengths.assembly_label, from_lengths.assembly_label_basis, from_lengths.assembly_label_source) == (
        "GRCh38", AssemblyLabelBasis.INFERRED, "contig_lengths",
    )  # fmt: skip
    (stated,) = sex_chromosome_evidence(load(write_vcf(tmp_path / "stated.vcf", parts, meta=("##reference=hg38.fa",))))
    assert (stated.assembly_label, stated.assembly_label_basis, stated.assembly_label_source) == ("GRCh38", AssemblyLabelBasis.STATED, "stated")
    assert from_lengths.ref_seqcol_digest is None and stated.ref_seqcol_digest is None  # REF not read


def facet_row(**columns: object) -> dict:
    """One `facets` row with every column null but those given."""
    return conform(pl.LazyFrame({name: [value] for name, value in columns.items()}), FACETS).collect().row(0, named=True)


LENGTH_37 = [{"name": "1", "length": 249_250_621, "md5": None}]


@pytest.mark.parametrize(
    ("columns", "source"),
    [
        ({"assembly_label_basis": "stated", "header_contigs": LENGTH_37}, "stated"),
        ({"assembly_label_basis": "inferred", "header_contigs": LENGTH_37, "vcf_version": "VCFv4.2"}, "contig_lengths"),
        ({"assembly_label_basis": "inferred", "ref_seqcol_digest": "d"}, "rows"),
        ({"assembly_label_basis": "inferred", "relation": "liftover_to", "ref_seqcol_digest": "d"}, "liftover"),
        # A VCF's REF read on a label its rows did not give: resolve_ref infers no build for a VCF.
        ({"assembly_label_basis": "inferred", "ref_seqcol_digest": "d", "vcf_version": "VCFv4.2"}, None),
        # A label a caller wrote into the frame, with no basis: used, its source unknown.
        ({}, None),
    ],
)
def test_the_label_source_follows_schema_decision_5(columns: dict, source: str | None) -> None:
    facet = facet_row(facet_id="f", assembly_label="GRCh37", **columns)
    assert assembly_label_source(facet) == source
    assert assembly_label_source({**facet, "assembly_label": None}) is None
    assert source is None or source in ASSEMBLY_LABEL_SOURCES


def test_the_label_source_is_a_closed_list_checked_on_the_model() -> None:
    found = evidence(frame(chr20(), *XY))
    with pytest.raises(ValueError, match="unknown assembly label source"):
        SexChromosomeEvidence(**{**found.model_dump(), "assembly_label_source": "guessed"})


# --- A person's coverage facets paired with its callsets (RM278) -----------------------------------


def callset(facet_id: str, x: pl.DataFrame, y: pl.DataFrame, *, build: str = "GRCh38", subject: str = "person-a") -> Loaded:
    """One whole-genome facet under `subject`, a source of its own (its own `source_sha256`)."""
    facets = conform(
        pl.LazyFrame({"facet_id": [facet_id], "assembly_label": [build], "callset_type": ["wgs"], "subject": [subject], "source_sha256": [facet_id * 64]}),
        FACETS,
    )
    calls = frame(chr20(), x, y).with_columns(facet_id=pl.lit(facet_id))
    return Loaded(facets=facets, calls=calls, blocks=empty(BLOCKS))


def covered(facet_id: str, x_depth: float, y_depth: float, *, build: str = "GRCh38", subject: str | None = "person-a", uneven: bool = False) -> Coverage:
    base = coverage_facet(x_depth, y_depth, build=build)
    intervals = base.intervals
    if uneven:  # an exome's shape on chr20 (`test_uneven_coverage_is_set_aside`)
        intervals = pl.concat([intervals, pl.DataFrame(windows("20", 5_000_001, 300, 2.0), schema=list(COVERAGE_INTERVALS), orient="row").cast(COVERAGE_INTERVALS).lazy()])
    return Coverage(facet=base.facet.model_copy(update={"facet_id": facet_id, "subject": subject}), intervals=intervals)


def by_facet(found: tuple[SexChromosomeEvidence, ...]) -> dict[str, SexChromosomeEvidence]:
    return {item.facet_id: item for item in found}


def without_pairing(found: SexChromosomeEvidence) -> SexChromosomeEvidence:
    return found.model_copy(update={"coverage_pairing": None})


def test_a_persons_one_coverage_pairs_with_its_callsets_on_its_build() -> None:
    f, g = callset("f", *XX), callset("g", *XX, build="GRCh37")
    cov = covered("cov", 30.0, 0.0)
    found = by_facet(person_sex_chromosome_evidence(Person.of(g, f, coverage=[cov], check=None)))
    assert list(found) == ["f", "g"], "one answer per facet, in the person's facet order"
    # The facet on the coverage's build gets exactly what the hand-paired facet-level call gives.
    assert without_pairing(found["f"]) == sex_chromosome_evidence(f, coverage=cov)[0]
    assert found["f"].inferred_karyotypic_sex is KaryotypicSex.XX
    assert found["f"].coverage_pairing.model_dump(exclude={"detail"}) == {"coverage_facet_id": "cov", "basis": "only_one_on_build", "candidates": ("cov",)}
    # The other build's facet is answered from its calls alone, and says why.
    assert without_pairing(found["g"]) == sex_chromosome_evidence(g)[0]
    assert found["g"].coverage is None
    assert found["g"].coverage_pairing.model_dump(exclude={"detail"}) == {"coverage_facet_id": None, "basis": "none_on_build", "candidates": ()}
    assert "GRCh38" in found["g"].coverage_pairing.detail


def test_without_coverage_the_person_level_answer_is_the_facet_levels() -> None:
    f = callset("f", *XY)
    (found,) = person_sex_chromosome_evidence(Person.of(f, check=None))
    assert without_pairing(found) == sex_chromosome_evidence(f)[0]
    assert (found.coverage_pairing.basis, found.coverage_pairing.coverage_facet_id) == ("no_coverage", None)


def test_two_usable_coverages_on_one_build_are_named_and_neither_used() -> None:
    f = callset("f", *XX)
    person = Person.of(f, coverage=[covered("c2", 30.0, 15.0), covered("c1", 30.0, 0.0)], check=None)
    (found,) = person_sex_chromosome_evidence(person)
    assert found.coverage_pairing.model_dump(exclude={"detail"}) == {"coverage_facet_id": None, "basis": "several_on_build", "candidates": ("c1", "c2")}
    assert without_pairing(found) == sex_chromosome_evidence(f)[0]
    assert found.inferred_karyotypic_sex is KaryotypicSex.UNKNOWN_KARYOTYPE, "the two would answer XX and XXY: no guess between them"
    # The caller states the pairing; it is used and says it was stated.
    (stated,) = person_sex_chromosome_evidence(person, pairs={"f": "c2"})
    assert without_pairing(stated) == sex_chromosome_evidence(f, coverage=covered("c2", 30.0, 15.0))[0]
    assert (stated.inferred_karyotypic_sex, stated.coverage_pairing.basis, stated.coverage_pairing.coverage_facet_id) == (KaryotypicSex.XXY, "stated", "c2")
    (none,) = person_sex_chromosome_evidence(person, pairs={"f": None})
    assert (none.coverage, none.coverage_pairing.basis, none.coverage_pairing.coverage_facet_id) == (None, "stated", None)


def test_of_several_on_one_build_the_only_usable_one_is_used_and_the_rest_named() -> None:
    f = callset("f", *XX)
    even, exome = covered("c1", 30.0, 0.0), covered("c2", 30.0, 0.0, uneven=True)
    (found,) = person_sex_chromosome_evidence(Person.of(f, coverage=[exome, even], check=None))
    assert found.coverage_pairing.model_dump(exclude={"detail"}) == {"coverage_facet_id": "c1", "basis": "only_usable_on_build", "candidates": ("c1", "c2")}
    assert "c2: coverage_uneven" in found.coverage_pairing.detail
    assert without_pairing(found) == sex_chromosome_evidence(f, coverage=even)[0]
    # None usable: none is used, and each one's reason is named.
    (neither,) = person_sex_chromosome_evidence(Person.of(f, coverage=[exome, covered("c3", 30.0, 0.0, uneven=True)], check=None))
    assert (neither.coverage_pairing.basis, neither.coverage_pairing.coverage_facet_id, neither.coverage) == ("none_usable_on_build", None, None)
    assert "c2: coverage_uneven" in neither.coverage_pairing.detail and "c3: coverage_uneven" in neither.coverage_pairing.detail


def test_a_stated_pairing_must_name_the_persons_own_facets() -> None:
    person = Person.of(callset("f", *XX), coverage=[covered("c1", 30.0, 0.0)], check=None)
    with pytest.raises(ValueError, match="not a facet of this person: x"):
        person_sex_chromosome_evidence(person, pairs={"x": "c1"})
    with pytest.raises(ValueError, match="not a coverage facet of this person: c9"):
        person_sex_chromosome_evidence(person, pairs={"f": "c9"})


def test_the_pairing_basis_is_a_closed_list_checked_on_the_model() -> None:
    assert COVERAGE_PAIRINGS == {"stated", "only_one_on_build", "only_usable_on_build", "several_on_build", "none_usable_on_build", "none_on_build", "no_coverage"}
    with pytest.raises(ValueError, match="unknown coverage pairing"):
        CoveragePairing(coverage_facet_id=None, basis="guessed", candidates=(), detail="")


# --- The person level: callsets weighed against each other (RM269) ---------------------------------

WITHHELD = (pl.DataFrame(), pl.DataFrame())
"""No X or Y rows: a facet that answers nothing (`too_few_sites` on both lines)."""
NO_DP_XY = (x_calls(40, 1960, None), y_calls(1000, None))
"""A one-X het share with no DP: `x_copies = one`, the karyotype withheld (`no_depth`)."""


def made_from(parent: Loaded, facet_id: str, x: pl.DataFrame, y: pl.DataFrame, *, relation: str = LIFTOVER_TO, build: str = "GRCh37") -> Loaded:
    """A facet one step made from `parent`: the parent's row with its own id, `derived_from` and
    `relation` (the steps copy their parent's row, `source_sha256` and subject included), and the
    parent's links plus the one this step adds. Y_START lies in the X-degenerate Y on both builds."""
    (parent_id,) = parent.facets.select("facet_id").collect()["facet_id"].to_list()
    facets = parent.facets.with_columns(
        facet_id=pl.lit(facet_id), derived_from=pl.lit(parent_id), relation=pl.lit(relation), assembly_label=pl.lit(build)
    )
    links = derived_links(parent.links, parent_id, facet_id, relation, step="liftover", step_version="test")
    calls = frame(chr20(), x, y).with_columns(facet_id=pl.lit(facet_id))
    return Loaded(facets=conform(facets, FACETS), calls=calls, blocks=empty(BLOCKS), links=links)


def pooled(facet_id: str, parents: list[str], x: pl.DataFrame, y: pl.DataFrame, relation: str = CONSENSUS) -> Loaded:
    """A facet made from several parents (a consensus or a merge), a source of its own, with its links."""
    facets = conform(
        pl.LazyFrame({"facet_id": [facet_id], "assembly_label": ["GRCh38"], "callset_type": ["wgs"], "subject": ["person-a"], "relation": [relation]}),
        FACETS,
    )
    provenance = {"step": "consensus", "step_version": "test", "derivation_inputs": []}
    links = links_frame([{"from_facet_id": p, "to_facet_id": facet_id, "relation": relation, "provenance": provenance} for p in parents])
    calls = frame(chr20(), x, y).with_columns(facet_id=pl.lit(facet_id))
    return Loaded(facets=facets, calls=calls, blocks=empty(BLOCKS), links=links)


def weighed(*loaded: Loaded, **kwargs: object) -> tuple[PersonKaryotype, Person]:
    person = Person.of(*loaded, check=None, coverage=kwargs.pop("coverage", ()))
    assert link_disagreements(person.facets.collect(), person.links.collect()) == [], "the fixtures' lineage is well formed"
    found = person_karyotype(person, **kwargs)
    # The facet-level answers are carried exactly as the facet-level function gives them.
    assert found.evidence == person_sex_chromosome_evidence(person, **kwargs)
    return found, person


def callsets(found: PersonKaryotype) -> dict[str, CallsetKaryotype]:
    return {item.root_facet_id: item for item in found.callsets}


def test_a_liftover_of_one_callset_answers_once_with_its_root() -> None:
    root = callset("r", *XY)
    found, _ = weighed(made_from(root, "l", *XY), root)
    (one,) = found.callsets
    assert (one.root_facet_id, one.root_loaded, one.members, one.standing_facet_id, one.basis) == ("r", True, ("r", "l"), "r", "root")
    assert (found.inferred_karyotypic_sex, found.reason, found.supported_by) == (KaryotypicSex.XY, None, ("r",))
    assert (found.x_copies, found.x_copies_supported_by) == (XCopies.ONE, ("r",))
    assert bool(found) and found.set_aside == ()
    # Two facets answered XY; the person-level answer counts one witness, not two.
    assert [item.inferred_karyotypic_sex for item in found.evidence] == [KaryotypicSex.XY, KaryotypicSex.XY]


def test_a_withheld_root_gives_way_to_the_member_that_answers() -> None:
    root = callset("r", *NO_DP_XY)
    found, _ = weighed(root, made_from(root, "l", *XY))
    (one,) = found.callsets
    r = by_facet(found.evidence)["r"]
    assert (r.inferred_karyotypic_sex, r.reason, r.x_copies) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "y_unusable", XCopies.ONE)
    assert (one.standing_facet_id, one.basis, one.inferred_karyotypic_sex) == ("l", "nearest_answering", KaryotypicSex.XY)
    assert found.inferred_karyotypic_sex is KaryotypicSex.XY and found.supported_by == ("r",)
    # A member that answers nothing is neither agreement nor disagreement: the root stands for itself.
    found, _ = weighed(root, made_from(root, "l", *WITHHELD))
    (one,) = found.callsets
    assert (one.standing_facet_id, one.basis, one.x_copies, one.reason) == ("r", "root", XCopies.ONE, "y_unusable")
    assert (found.reason, found.x_copies, found.x_copies_supported_by) == ("all_withheld", XCopies.ONE, ("r",))
    # No member says anything: the root's withheld answer stands, with its own reason.
    blank = callset("r", *WITHHELD)
    found, _ = weighed(made_from(blank, "l", *WITHHELD), blank)
    (one,) = found.callsets
    assert (one.standing_facet_id, one.basis, one.reason) == ("r", "none_answered", "too_few_sites")


def test_a_lineage_joins_through_a_middle_facet_the_person_does_not_hold() -> None:
    root = callset("r", *XY)
    middle = made_from(root, "n", *XY, relation=NORMALIZED, build="GRCh38")
    found, _ = weighed(made_from(middle, "l", *XY), root)
    (one,) = found.callsets
    assert (one.root_facet_id, one.members) == ("r", ("r", "l"))
    # Without the root, its descendants still form one callset named for it.
    found, _ = weighed(made_from(middle, "l", *XY), made_from(root, "g", *XY))
    (one,) = found.callsets
    assert (one.root_facet_id, one.root_loaded, one.members, one.basis, one.standing_facet_id) == ("r", False, ("g", "l"), "nearest_answering", "g")


def test_a_callset_whose_members_disagree_is_withheld_with_its_own_reason() -> None:
    root, other = callset("r", *XY), callset("o", *XY)
    found, _ = weighed(root, made_from(root, "l", XX[0], XY[1]), other)
    torn = callsets(found)["r"]
    assert (torn.basis, torn.standing_facet_id, torn.reason, torn.inferred_karyotypic_sex) == (
        "members_disagree", None, "members_disagree", KaryotypicSex.UNKNOWN_KARYOTYPE,
    )  # fmt: skip
    assert "r: XY" in torn.detail and "l: XXY" in torn.detail
    # The other callset answers XY, but one piece of evidence contradicting itself is not set aside for it.
    assert callsets(found)["o"].inferred_karyotypic_sex is KaryotypicSex.XY
    assert (found.inferred_karyotypic_sex, found.reason, found.supported_by) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "members_disagree", ())
    assert (found.x_copies, found.x_copies_reason) == (None, "members_disagree")


def test_independent_callsets_that_agree_support_the_answer_together() -> None:
    found, _ = weighed(callset("b", *XY, build="GRCh37"), callset("a", *XY))
    assert [item.basis for item in found.callsets] == ["only_member", "only_member"]
    assert (found.inferred_karyotypic_sex, found.supported_by, found.x_copies_supported_by) == (KaryotypicSex.XY, ("a", "b"), ("a", "b"))


def test_independent_callsets_that_disagree_withhold_and_list_each_answer() -> None:
    found, _ = weighed(callset("a", *XY), callset("b", XX[0], XY[1]))
    assert (found.inferred_karyotypic_sex, found.reason, found.supported_by) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "facets_disagree", ())
    assert {root: item.inferred_karyotypic_sex for root, item in callsets(found).items()} == {"a": KaryotypicSex.XY, "b": KaryotypicSex.XXY}
    assert "a: XY" in found.detail and "b: XXY" in found.detail
    assert (found.x_copies, found.x_copies_reason) == (None, "facets_disagree")
    assert not found


def test_a_measured_genome_does_not_outrank_an_estimated_array() -> None:
    """The lead's pin: no ranking between independent callsets. An AncestryDNA export (Y read by
    default) answering XX against a genome's XY is a disagreement, not the genome's answer."""
    found, _ = weighed(callset("g", *XY), array_facets({"a": "AncestryDNA"}, subject="person-a"))
    assert callsets(found)["a"].inferred_karyotypic_sex is KaryotypicSex.XX
    assert (found.reason, found.inferred_karyotypic_sex) == ("facets_disagree", KaryotypicSex.UNKNOWN_KARYOTYPE)


def test_a_withheld_karyotype_still_contradicts_with_the_x_it_states() -> None:
    """A 23andMe export's two-letter X: Y not read (`vendor_asserted`), so its karyotype is withheld,
    but it states two or more X, which a genome's XY contradicts and a coverage-measured XX does not."""
    array = array_facets({"t": "23andMe"}, subject="person-a")
    found, _ = weighed(callset("g", *XY), array)
    t = callsets(found)["t"]
    assert (t.inferred_karyotypic_sex, t.reason, t.x_copies) == (KaryotypicSex.UNKNOWN_KARYOTYPE, "y_unusable", XCopies.TWO_OR_MORE)
    assert (found.reason, found.x_copies_reason) == ("facets_disagree", "facets_disagree")
    found, _ = weighed(callset("g", *XX), array, coverage=[covered("cov", 30.0, 0.0)])
    assert (found.inferred_karyotypic_sex, found.supported_by) == (KaryotypicSex.XX, ("g",)), "a fitting line is listed, not support"
    assert (found.x_copies, found.x_copies_supported_by) == (XCopies.TWO_OR_MORE, ("g", "t"))


def test_one_withheld_callset_beside_one_answering_leaves_the_answer() -> None:
    found, _ = weighed(callset("a", *XY), callset("b", *WITHHELD))
    assert callsets(found)["b"].reason == "too_few_sites"
    assert (found.inferred_karyotypic_sex, found.reason, found.supported_by) == (KaryotypicSex.XY, None, ("a",))
    assert found.x_copies_supported_by == ("a",)


def test_only_withheld_callsets_withhold_with_their_own_reason() -> None:
    found, _ = weighed(callset("a", *WITHHELD), callset("b", *WITHHELD))
    assert (found.inferred_karyotypic_sex, found.reason, found.x_copies, found.x_copies_reason) == (
        KaryotypicSex.UNKNOWN_KARYOTYPE, "all_withheld", None, "all_withheld",
    )  # fmt: skip
    assert [item.reason for item in found.callsets] == ["too_few_sites", "too_few_sites"]


def test_a_consensus_and_what_is_made_from_it_are_set_aside_and_its_parents_still_count() -> None:
    a, b = callset("a", *XY), callset("b", *XY)
    consensus = pooled("c", ["a", "b"], XX[0], XY[1])  # would answer XXY, against both parents
    found, _ = weighed(a, b, consensus, made_from(consensus, "d", XX[0], XY[1], build="GRCh38"))
    assert [(item.facet_id, item.reason) for item in found.set_aside] == [("c", "multi_parent"), ("d", "descends_from_multi_parent")]
    assert by_facet(found.evidence)["c"].inferred_karyotypic_sex is KaryotypicSex.XXY, "measured at the facet level as always"
    assert [item.root_facet_id for item in found.callsets] == ["a", "b"]
    assert (found.inferred_karyotypic_sex, found.supported_by) == (KaryotypicSex.XY, ("a", "b"))
    # A merge is set aside the same way; with nothing else held, no callset is left.
    found, _ = weighed(pooled("m", ["a", "b"], *XY, relation=MERGE))
    assert (found.callsets, found.reason, found.x_copies_reason) == ((), "no_callset", "no_callset")
    assert [(item.facet_id, item.reason) for item in found.set_aside] == [("m", "multi_parent")]


def test_facets_split_by_caller_from_one_file_are_separate_callsets() -> None:
    """Kin links (`caller_split`) join facets of one source neither of which was made from the other."""
    two = callset("v", *XY)
    facets = pl.concat([two.facets, two.facets.with_columns(facet_id=pl.lit("f"))])
    calls = pl.concat([two.calls, frame(chr20(), *XY).with_columns(facet_id=pl.lit("f"))])
    found, _ = weighed(Loaded(facets=facets, calls=calls, blocks=empty(BLOCKS), links=kin_links(["v", "f"], CALLER_SPLIT)))
    assert [item.root_facet_id for item in found.callsets] == ["f", "v"]
    assert found.supported_by == ("f", "v")


def test_the_person_level_answer_does_not_depend_on_the_order_facets_are_given() -> None:
    z = callset("z", *XY)
    loads = [z, made_from(z, "a", *XY), made_from(made_from(z, "m", *XY, relation=NORMALIZED, build="GRCh38"), "b", *XY), callset("k", *WITHHELD)]
    first, _ = weighed(*loads)
    for seed in range(5):
        shuffled = loads[:]
        random.Random(seed).shuffle(shuffled)
        assert weighed(*shuffled)[0] == first
    assert [item.root_facet_id for item in first.callsets] == ["k", "z"]
    # Members run nearest the root first, then by id: the root, its one-step liftover, then the two-step one.
    assert callsets(first)["z"].members == ("z", "a", "b")


def test_the_person_level_vocabularies_are_closed_lists_checked_on_the_models() -> None:
    assert PERSON_KARYOTYPE_REASONS == {"facets_disagree", "members_disagree", "all_withheld", "no_callset"}
    assert CALLSET_BASES == {"only_member", "root", "nearest_answering", "none_answered", "members_disagree"}
    assert SET_ASIDE_REASONS == {"multi_parent", "descends_from_multi_parent"}
    found, _ = weighed(callset("a", *XY))
    with pytest.raises(ValueError, match="unknown reason"):
        PersonKaryotype(**{**found.model_dump(), "reason": "nearest_class"})
    with pytest.raises(ValueError, match="unknown callset basis"):
        CallsetKaryotype(**{**found.callsets[0].model_dump(), "basis": "deepest"})
    with pytest.raises(ValueError, match="unknown set-aside reason"):
        SetAsideFacet(facet_id="c", reason="ranked_lower", detail="")
