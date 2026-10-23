"""Graded genotype evidence across sources (layer 7; docs/SCHEMA.md, "Graded genotype evidence").

Expected values are computed here from the fixture's own text, by the VCF spec's definitions, never
read back from the function under test.
"""

import math
from pathlib import Path

import polars as pl
import pytest

from just_vcf import AbsenceMeans, Loaded, load, read_coverage
from just_vcf.evidence import (
    GENOTYPE_EVIDENCE,
    GENOTYPE_EVIDENCE_AT,
    ConfidenceKind,
    GradeBasis,
    VectorPrior,
    genotype_evidence,
    genotype_evidence_at,
)
from test_load_vcf import vcf_text

DS_LINE = '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">'
GL_LINE = '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods">'
FMT = ("GT", "GQ", "PL", "GP")


def _vcf(tmp_path: Path, records: list[str], *, fmt: tuple[str, ...] = FMT, version: str = "VCFv4.3", extra: tuple[str, ...] = (),
         info: tuple[str, ...] = ("END",), name: str = "calls.vcf") -> Loaded:  # fmt: skip
    path = tmp_path / name
    path.write_text(vcf_text(records, fmt=fmt, info=info, version=version, extra_meta=extra))
    return load(path)


def _graded(loaded: Loaded, **kwargs) -> pl.DataFrame:
    return genotype_evidence(loaded, **kwargs).collect()


def _row(graded: pl.DataFrame, pos: int, origin: str = "assay") -> dict:
    rows = graded.filter((pl.col("pos") == pos) & (pl.col("evidence_basis") == origin))
    assert rows.height == 1, rows
    return rows.row(0, named=True)


def _vcf_genotype_order(alleles: int, ploidy: int) -> list[tuple[int, ...]]:
    """The VCF spec's genotype ordering, written out by its loop (VCF 4.3 §1.6.2, "Genotype ordering")."""
    if ploidy == 1:
        return [(a,) for a in range(alleles)]
    return [(j, k) for k in range(alleles) for j in range(k + 1)]


def _phred_top(probabilities: list[float]) -> float:
    total = sum(probabilities)
    top = max(probabilities)
    rest = total - top
    return math.inf if rest == 0 else -10 * math.log10(rest / total)


def test_the_frame_is_graded_never_boolean() -> None:
    """Charter 7 and CLAUDE.md §2: no column of either answer is a boolean "called"."""
    for schema in (GENOTYPE_EVIDENCE, GENOTYPE_EVIDENCE_AT):
        assert pl.Boolean not in schema.dtypes()
    assert {m.value for m in GradeBasis} >= {"gp", "gl", "gq", "gt_only", "block", "absent_reference", "absent_unknown"}


def test_likelihoods_are_normalized_with_a_flat_prior(tmp_path: Path) -> None:
    pl_text = "51,0,36,93,92,86"
    loaded = _vcf(tmp_path, [f"1\t114\t.\tT\tTC,TA\t213\t.\t.\tGT:GQ:PL:GP\t0/1:36:{pl_text}:."])
    row = _row(_graded(loaded), 114)
    weights = [10 ** (-int(v) / 10) for v in pl_text.split(",")]
    expected = [w / sum(weights) for w in weights]
    assert row["grade_basis"] == GradeBasis.GL.value and row["prior"] == VectorPrior.FLAT.value
    assert row["confidence_kind"] == ConfidenceKind.QUANTIFIED.value
    assert row["genotype_alleles"] == ["T", "TC", "TA"]
    assert row["genotype_probabilities"] == pytest.approx(expected, rel=1e-5)
    assert sum(row["genotype_probabilities"]) == pytest.approx(1.0)
    assert row["confidence"] == pytest.approx(_phred_top(expected), rel=1e-4)
    order = _vcf_genotype_order(3, 2)
    assert order[row["genotype_index"]] == (0, 1), "the top genotype, in VCF order, is the PL 0 one"


def test_gl_and_pl_stating_one_likelihood_grade_alike(tmp_path: Path) -> None:
    """The frame folds GL and PL into `gl` (Decision 3), so the basis is `gl` for both and the values agree."""
    by_pl = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:PL\t0/1:40,0,300"], fmt=("GT", "PL"), name="pl.vcf")
    by_gl = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GL\t0/1:-4,0,-30"], fmt=("GT",), extra=(GL_LINE,), name="gl.vcf")
    a, b = _row(_graded(by_pl), 100), _row(_graded(by_gl), 100)
    assert a["grade_basis"] == b["grade_basis"] == GradeBasis.GL.value
    assert a["genotype_probabilities"] == pytest.approx(b["genotype_probabilities"], rel=1e-6)
    assert a["confidence"] == pytest.approx(b["confidence"], rel=1e-6)


def test_extreme_likelihoods_keep_a_finite_confidence(tmp_path: Path) -> None:
    """A PL of thousands is 1 - p of about 10^-300: computed in log space, never rounded to infinity."""
    loaded = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:PL\t1/1:9000,3000,0"], fmt=("GT", "PL"))
    row = _row(_graded(loaded), 100)
    assert row["confidence"] == pytest.approx(3000.0, rel=1e-6)
    assert row["genotype_index"] == 2


def test_posteriors_are_taken_as_given_and_win_over_likelihoods(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, [
        "1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/1:20:30,0,40:0.02,0.95,0.03",
        "1\t101\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t1/1:.:.:0,0,1",
    ])  # fmt: skip
    graded = _graded(loaded)
    first, certain = _row(graded, 100), _row(graded, 101)
    assert first["grade_basis"] == GradeBasis.GP.value and first["prior"] == VectorPrior.SOURCE.value
    assert first["genotype_probabilities"] == pytest.approx([0.02, 0.95, 0.03], rel=1e-6)
    assert first["confidence"] == pytest.approx(_phred_top([0.02, 0.95, 0.03]), rel=1e-5)
    assert first["genotype_index"] == 1
    assert certain["confidence"] == math.inf, "every other genotype has probability zero as written"
    assert certain["confidence_kind"] == ConfidenceKind.QUANTIFIED.value


def test_gq_alone_puts_its_confidence_on_the_called_genotype(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, [
        "1\t100\t.\tA\tG,T\t50\t.\t.\tGT:GQ:PL:GP\t1/2:27:.:.",
        "1\t101\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/1:0:.:.",
        "X\t102\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t1:13:.:.",
    ])  # fmt: skip
    graded = _graded(loaded)
    row = _row(graded, 100)
    order = _vcf_genotype_order(3, 2)
    called = order.index((1, 2))
    error = 10 ** (-27 / 10)
    expected = [1 - error if i == called else error / (len(order) - 1) for i in range(len(order))]
    assert row["grade_basis"] == GradeBasis.GQ.value and row["prior"] == VectorPrior.EVEN_REMAINDER.value
    assert row["genotype_index"] == called
    assert row["genotype_probabilities"] == pytest.approx(expected, rel=1e-9)
    assert row["confidence"] == 27.0
    zero = _row(graded, 101)
    assert zero["confidence"] == 0.0 and zero["confidence_kind"] == ConfidenceKind.QUANTIFIED.value, "GQ 0 is a number: no confidence"
    haploid = _row(graded, 102)
    assert haploid["genotype_index"] == 1 and len(haploid["genotype_probabilities"]) == 2


def test_gq_grades_a_reference_genotype_only_when_asked(tmp_path: Path) -> None:
    """VCF defines GQ as conditional on the site being variant (SCHEMA, VCF 4.5 points)."""
    loaded = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/0:45:.:."])
    plain, asked = _row(_graded(loaded), 100), _row(_graded(loaded, gq_on_reference=True), 100)
    assert plain["grade_basis"] == GradeBasis.GT_ONLY.value and plain["confidence"] is None
    assert plain["genotype_index"] == 0
    assert asked["grade_basis"] == GradeBasis.GQ.value and asked["confidence"] == 45.0
    assert asked["genotype_probabilities"][0] == pytest.approx(1 - 10 ** -4.5)


def test_genotypes_without_a_quantity_are_asserted_or_unknown_never_zero(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, [
        "1\t100\t.\tA\tG,T\t50\t.\t.\tGT:GQ:PL:GP\t1/2:.:.:.",
        "1\t101\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t./.:.:.:.",
        "1\t102\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/.:30:.:.",
        "1\t103\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t./.:0:.:.",
    ])  # fmt: skip
    graded = _graded(loaded)
    hard = _row(graded, 100)
    assert hard["grade_basis"] == GradeBasis.GT_ONLY.value and hard["confidence_kind"] == ConfidenceKind.ASSERTED.value
    assert hard["confidence"] is None and hard["genotype_probabilities"] is None
    assert _vcf_genotype_order(3, 2)[hard["genotype_index"]] == (1, 2)
    for pos in (101, 103):
        none = _row(graded, pos)
        assert none["grade_basis"] == GradeBasis.NO_CALL.value and none["confidence_kind"] == ConfidenceKind.UNKNOWN.value
        assert none["confidence"] is None and none["genotype_index"] is None
    half = _row(graded, 102)
    assert half["grade_basis"] == GradeBasis.GT_ONLY.value, "a half-call has no one genotype for GQ to grade"
    assert half["genotype_index"] is None


def test_a_dosage_alone_has_its_own_basis(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:DS\t0/1:0.98"], fmt=("GT",), extra=(DS_LINE,))
    row = _row(_graded(loaded), 100)
    assert row["grade_basis"] == GradeBasis.DOSAGE.value and row["confidence"] is None and row["genotype_probabilities"] is None


def test_a_vector_that_does_not_fit_the_ploidy_falls_to_the_next_basis(tmp_path: Path) -> None:
    triploid = [10 ** (-v / 10) for v in (60, 0, 20, 90)]
    loaded = _vcf(tmp_path, [
        "1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/0/1:25:60,0,20,90:.",
        "1\t101\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/0/1:25:60,0,20:.",
    ])  # fmt: skip
    graded = _graded(loaded)
    fits = _row(graded, 100)
    assert fits["grade_basis"] == GradeBasis.GL.value
    assert fits["genotype_probabilities"] == pytest.approx([w / sum(triploid) for w in triploid], rel=1e-5)
    short = _row(graded, 101)
    assert short["grade_basis"] == GradeBasis.GQ.value and short["genotype_probabilities"] is None, "no vector beyond diploid from GQ"
    assert short["confidence"] == 25.0


GVCF = [
    "1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ:PL:GP\t0/0:60:0,60,900:.",
    "1\t114\t.\tTAC\tT,<NON_REF>\t213\t.\t.\tGT:GQ:PL:GP\t0/1:99:51,0,36,93,92,86:.",
    "1\t120\t.\tC\t<*>\t.\t.\tEND=125\tGT:GQ:PL:GP\t./.:.:.:.",
    "1\t126\t.\tC\t<*>\t.\t.\tEND=130\tGT:GQ:PL:GP\t0/0:33:.:.",
]


def test_blocks_are_graded_by_their_likelihoods_or_asserted(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, GVCF)
    graded = _graded(loaded)
    ref = _row(graded, 100, "ref_block")
    weights = [10 ** (-v / 10) for v in (0, 60, 900)]
    assert ref["grade_basis"] == GradeBasis.GL.value and ref["block_kind"] == "ref"
    assert ref["genotype_alleles"] == [None, "<*>"], "a block keeps no REF base"
    assert ref["genotype_probabilities"] == pytest.approx([w / sum(weights) for w in weights], rel=1e-6)
    assert ref["genotype_index"] == 0
    no_call = _row(graded, 120, "ref_block")
    assert no_call["grade_basis"] == GradeBasis.NO_CALL.value and no_call["confidence_kind"] == ConfidenceKind.UNKNOWN.value
    bare = _row(graded, 126, "ref_block")
    assert bare["grade_basis"] == GradeBasis.BLOCK.value and bare["confidence"] is None and bare["genotype_index"] == 0
    with_gq = _row(_graded(loaded, gq_on_reference=True), 126, "ref_block")
    assert with_gq["grade_basis"] == GradeBasis.GQ.value and with_gq["confidence"] == 33.0
    assert with_gq["genotype_probabilities"] is None, "a block's ploidy is not in the frame without its likelihoods"
    keys = graded.select("pos", "evidence_basis").rows()
    assert keys == sorted(keys, key=lambda k: (k[0], k[1])), "canonical order, calls before blocks at one position"
    assert graded.height == pl.collect_all([loaded.calls.select(pl.len()), loaded.blocks.select(pl.len())])[0].item() + 3


def test_an_array_genotype_is_asserted_and_its_no_call_unknown(tmp_path: Path) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(
        "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n"
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
        "# rsid\tchromosome\tposition\tgenotype\n"
        "rs4477212\t1\t82154\tAA\nrs3094315\t1\t752566\tAG\nrs3131972\t1\t752721\t--\n"
    )
    loaded = load(path)
    graded = _graded(loaded)
    assert graded["grade_basis"].to_list() == ["gt_only", "gt_only", "no_call"]
    assert graded["confidence"].null_count() == graded.height
    assert graded["confidence_kind"].to_list() == ["asserted", "asserted", "unknown"]


def _weak_reference(loaded: Loaded) -> Loaded:
    """The same facet read as a whole-genome callset without blocks (Decision 4's `weak_reference`)."""
    enum = pl.Enum([m.value for m in AbsenceMeans])
    return Loaded(facets=loaded.facets.with_columns(pl.lit("weak_reference").cast(enum).alias("absence_means")),
                  calls=loaded.calls, blocks=loaded.blocks.clear())  # fmt: skip


def test_asked_positions_read_records_blocks_and_absence(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, GVCF)
    sites = pl.DataFrame({"chrom": ["1", "1", "1", "1", "1", "1"], "pos": [105, 115, 116, 140, 114, 105]})
    at = genotype_evidence_at(loaded, sites, variant_class="snv")
    assert at.schema == GENOTYPE_EVIDENCE_AT
    assert at["pos"].to_list() == [105, 115, 116, 140, 114], "sites in the order given, a duplicate once"
    by_pos = {r["pos"]: r for r in at.iter_rows(named=True)}
    assert by_pos[105]["evidence_basis"] == "ref_block" and (by_pos[105]["record_pos"], by_pos[105]["record_end_pos"]) == (100, 113)
    assert by_pos[115]["record_pos"] == 114 and by_pos[115]["grade_basis"] == "gl", "inside a called deletion is not absent"
    assert by_pos[116]["record_pos"] == 114
    assert by_pos[140]["grade_basis"] == GradeBasis.ABSENT_UNKNOWN.value, "a gVCF's gap says nothing (absence_means = blocks)"
    assert by_pos[140]["record_pos"] is None and by_pos[140]["absence_means"] == "blocks"
    assert by_pos[140]["confidence"] is None and by_pos[140]["confidence_kind"] == ConfidenceKind.UNKNOWN.value


def test_a_whole_genome_silence_is_asserted_reference_for_the_classes_it_called(tmp_path: Path) -> None:
    loaded = _weak_reference(_vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/1:30:.:."]))
    sites = pl.DataFrame({"chrom": ["1"], "pos": [500]})
    snv = genotype_evidence_at(loaded, sites, variant_class="snv").row(0, named=True)
    assert snv["grade_basis"] == GradeBasis.ABSENT_REFERENCE.value and snv["confidence_kind"] == ConfidenceKind.ASSERTED.value
    assert snv["confidence"] is None and snv["genotype_index"] == 0 and snv["absence_means"] == "weak_reference"
    sv = genotype_evidence_at(loaded, sites, variant_class="sv").row(0, named=True)
    assert sv["grade_basis"] == GradeBasis.ABSENT_UNKNOWN.value, "a caller that wrote no SVs says nothing about SVs"


def test_coverage_rides_beside_the_answer_and_changes_no_basis(tmp_path: Path) -> None:
    loaded = _weak_reference(_vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/1:30:.:."]))
    bed = tmp_path / "depth.bed"
    bed.write_text("1\t399\t450\t0\n1\t450\t600\t31\n")
    coverage = read_coverage(bed, window=None)
    sites = pl.DataFrame({"chrom": ["1", "1", "1"], "pos": [420, 500, 900]})
    plain = genotype_evidence_at(loaded, sites, variant_class="snv")
    covered = genotype_evidence_at(loaded, sites, variant_class="snv", coverage=coverage, min_depth=10)
    assert covered.select(GENOTYPE_EVIDENCE_AT.names()).equals(plain)
    assert covered["callable"].to_list() == [False, True, None]
    assert covered["callable_reason"].to_list() == [None, None, "no_interval"]
    with pytest.raises(ValueError, match="min_depth"):
        genotype_evidence_at(loaded, sites, variant_class="snv", coverage=coverage)


def test_asking_positions_of_several_facets_is_refused(tmp_path: Path) -> None:
    a = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/1:30:.:."], name="a.vcf")
    b = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t1/1:30:.:."], name="b.vcf")
    both = Loaded(facets=pl.concat([a.facets, b.facets]), calls=pl.concat([a.calls, b.calls]), blocks=pl.concat([a.blocks, b.blocks]))
    assert genotype_evidence(both).collect().height == 2, "the row-wise frame takes any number of facets"
    with pytest.raises(ValueError, match="one facet"):
        genotype_evidence_at(both, pl.DataFrame({"chrom": ["1"], "pos": [100]}), variant_class="snv")
