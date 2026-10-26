"""Reference-supported positions (RM3): what one facet says about the reference allele at asked positions.

Expected grades are computed here from the fixture's own text (a PL normalized under a flat prior, by
the VCF spec), never read back from the function under test. Facets whose type a small fixture cannot
earn (an exome, a CGI-style `reference` genome) are the loaded facet with its `facets` columns set, as
`tests/test_evidence.py` does; the whole-genome gates run on whole-genome-shaped files loaded for real.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import (
    REFERENCE_SUPPORT,
    AbsenceMeans,
    AbsenceReason,
    BlockGrade,
    CallsetType,
    ConfidenceKind,
    ContentFilter,
    GradeBasis,
    Loaded,
    RecordReason,
    SupportKind,
    load,
    reference_support,
)
from just_vcf.schema import enum_dtype
from test_content_filter import genome_vcf
from test_evidence import BLOCKS, GVCF, _block_records, _expected_block, _flat_pl, _phred_top, _vcf

SNV_HET = "0/1:99:51,0,36:."
"""A het SNV with a PL, at positions the tests choose."""


def _support(loaded: Loaded, positions: list[tuple[str, int]], variant_class: str = "snv", **kwargs) -> dict[tuple[str, int], dict]:
    sites = pl.DataFrame({"chrom": [c for c, _ in positions], "pos": [p for _, p in positions]})
    frame = reference_support(loaded, sites, variant_class=variant_class, **kwargs)
    assert isinstance(frame, pl.LazyFrame), "lazy over the facet's frames"
    answer = frame.collect()
    assert answer.schema == REFERENCE_SUPPORT
    assert answer.height == len(set(positions)), "one row per distinct position"
    return {(row["chrom"], row["pos"]): row for row in answer.iter_rows(named=True)}


def _as(loaded: Loaded, **columns: object) -> Loaded:
    """The loaded facet with some `facets` columns set, each cast to the column's own dtype."""
    facets = loaded.facets
    dtypes = facets.collect_schema()
    return Loaded(
        facets=facets.with_columns(pl.lit(value).cast(dtypes[name]).alias(name) for name, value in columns.items()),
        calls=loaded.calls,
        blocks=loaded.blocks,
    )


def _assert_none(row: dict, reason: AbsenceReason) -> None:
    assert (row["support"], row["reason"]) == (SupportKind.NONE.value, reason.value), row
    assert row["grade_basis"] == GradeBasis.ABSENT_UNKNOWN.value and row["confidence_kind"] == ConfidenceKind.UNKNOWN.value
    assert row["genotype_index"] is None and row["confidence"] is None and row["records"] == 0


def _assert_absence(row: dict, reason: AbsenceReason, absence_means: AbsenceMeans) -> None:
    assert (row["support"], row["reason"], row["absence_means"]) == (SupportKind.ABSENCE.value, reason.value, absence_means.value), row
    assert row["grade_basis"] == GradeBasis.ABSENT_REFERENCE.value and row["confidence_kind"] == ConfidenceKind.ASSERTED.value
    assert row["genotype_index"] == 0, "absence that reads as reference states hom-ref"
    assert row["confidence"] is None, "no number is made up for a silence (Decision 3)"
    assert row["record_pos"] is None and row["records"] == 0


def test_the_answer_is_graded_never_boolean() -> None:
    assert pl.Boolean not in REFERENCE_SUPPORT.dtypes()
    assert REFERENCE_SUPPORT["support"] == enum_dtype(SupportKind)
    reasons = set(REFERENCE_SUPPORT["reason"].categories.to_list())
    assert reasons == {m.value for m in RecordReason} | {m.value for m in AbsenceReason}


def test_a_gvcf_answers_by_its_blocks_with_inclusive_edges(tmp_path: Path) -> None:
    """GVCF: a `ref` block 100–113 (PL 0,60,900), a het deletion TAC>T at 114 (span 114–116), a `no_call`
    block 120–125, a bare `ref` block 126–130 (GQ 33, no PL). Its absence means `blocks`."""
    loaded = _vcf(tmp_path, GVCF)
    asked = [("1", p) for p in (99, 100, 113, 114, 115, 116, 117, 119, 120, 125, 126, 130, 131)]
    answer = _support(loaded, asked)
    for pos in (100, 113):
        row = answer[("1", pos)]
        assert (row["support"], row["reason"], row["records"]) == ("block", "ref_block", 1)
        assert (row["record_pos"], row["record_end_pos"]) == (100, 113), "both ends of a block are inside it"
        assert row["grade_basis"] == "gl" and row["genotype_index"] == 0
        assert row["confidence"] == pytest.approx(_phred_top(_flat_pl((0, 60, 900))), rel=1e-6)
        assert row["genotype_probabilities"] == pytest.approx(_flat_pl((0, 60, 900)), rel=1e-6)
    deletion_pl = (51, 0, 36, 93, 92, 86)
    for pos, reason in ((114, "call_at_position"), (115, "call_spans_position"), (116, "call_spans_position")):
        row = answer[("1", pos)]
        assert (row["support"], row["reason"], row["record_pos"], row["record_end_pos"]) == ("call", reason, 114, 116)
        assert row["genotype_index"] == _flat_pl(deletion_pl).index(max(_flat_pl(deletion_pl))) == 1, "het: one haplotype deleted"
        assert row["confidence"] == pytest.approx(_phred_top(_flat_pl(deletion_pl)), rel=1e-6)
    no_call = answer[("1", 125)]
    assert (no_call["support"], no_call["reason"], no_call["record_pos"]) == ("block", "no_call_block", 120)
    assert no_call["grade_basis"] == "no_call" and no_call["genotype_index"] is None
    bare = answer[("1", 130)]
    assert (bare["reason"], bare["grade_basis"], bare["confidence"], bare["genotype_index"]) == ("ref_block", "block", None, 0)
    assert _support(loaded, [("1", 130)], gq_on_reference=True)[("1", 130)]["confidence"] == 33.0
    for pos in (99, 117, 119, 131):
        _assert_none(answer[("1", pos)], AbsenceReason.OUTSIDE_BLOCKS)
        assert answer[("1", pos)]["absence_means"] == "blocks"


@pytest.mark.parametrize("rule", list(BlockGrade))
def test_each_block_rule_grades_the_reference_it_supports(tmp_path: Path, rule: BlockGrade) -> None:
    loaded = _vcf(tmp_path, _block_records(), fmt=("GT", "GQ", "PL"))
    answer = _support(loaded, [("1", p + 9) for p in BLOCKS], block_grade=rule)
    for pos, (gq, pls) in BLOCKS.items():
        row = answer[("1", pos + 9)]
        basis, confidence, vector = _expected_block(rule, gq, pls)
        assert (row["support"], row["reason"], row["record_end_pos"]) == ("block", "ref_block", pos + 9), (rule, pos)
        assert row["grade_basis"] == basis, (rule, pos)
        assert row["confidence"] == (pytest.approx(confidence, rel=1e-6) if confidence is not None else None), (rule, pos)
        assert row["genotype_probabilities"] == (pytest.approx(vector, rel=1e-6) if vector is not None else None), (rule, pos)


GVCF_OVERLAPS = [
    "1\t200\t.\tA\t<NON_REF>\t.\t.\tEND=210\tGT:GQ:PL:GP\t0/0:60:0,60,900:.",
    "1\t205\t.\tC\tT,<NON_REF>\t50\t.\t.\tGT:GQ:PL:GP\t0/1:99:51,0,36,90,90,90:.",
    "1\t300\t.\tACGT\tA,<NON_REF>\t50\t.\t.\tGT:GQ:PL:GP\t0/1:40:40,0,50,90,90,90:.",
    "1\t302\t.\tG\t<NON_REF>\t.\t.\tEND=310\tGT:GQ:PL:GP\t0/0:30:0,30,450:.",
    "1\t400\t.\tA\tG\t50\t.\t.\tGT:GQ:PL:GP\t0/0:30:0,30,300:.",
    "1\t400\t.\tA\tT\t50\t.\t.\tGT:GQ:PL:GP\t0/1:30:30,0,300:.",
]


def test_a_call_answers_before_a_block_and_a_variant_before_a_reference(tmp_path: Path) -> None:
    """Where records overlap, one answers and `records` counts them: a call starting at the position, then
    the call spanning it, then the block; among calls of one span, the one stating an ALT."""
    loaded = _vcf(tmp_path, GVCF_OVERLAPS)
    answer = _support(loaded, [("1", 204), ("1", 205), ("1", 301), ("1", 302), ("1", 303), ("1", 304), ("1", 400)])
    inside = answer[("1", 205)]
    assert (inside["support"], inside["reason"], inside["records"], inside["record_pos"]) == ("call", "call_at_position", 2, 205)
    assert inside["genotype_index"] == 1 and inside["confidence"] == pytest.approx(_phred_top(_flat_pl((51, 0, 36, 90, 90, 90))), rel=1e-6)
    assert (answer[("1", 204)]["support"], answer[("1", 204)]["records"]) == ("block", 1)
    for pos in (302, 303):
        row = answer[("1", pos)]
        assert (row["support"], row["reason"], row["records"], row["record_pos"], row["record_end_pos"]) == ("call", "call_spans_position", 2, 300, 303), pos
        assert row["genotype_index"] == 1, "a het deletion over the block states an ALT there"
    assert (answer[("1", 301)]["reason"], answer[("1", 301)]["records"]) == ("call_spans_position", 1)
    after = answer[("1", 304)]
    assert (after["support"], after["reason"], after["record_pos"], after["records"]) == ("block", "ref_block", 302, 1)
    assert after["confidence"] == pytest.approx(_phred_top(_flat_pl((0, 30, 450))), rel=1e-6)
    both = answer[("1", 400)]
    assert (both["support"], both["records"], both["genotype_index"]) == ("call", 2, 1), "a hom-ref row never hides the het row at the site"
    assert both["genotype_alleles"] == ["A", "T"]


def test_positions_come_out_in_canonical_order_once(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, GVCF)
    sites = pl.DataFrame({"chrom": ["2", "1", "1", "2", "1"], "pos": [5, 900, 100, 5, 113]})
    answer = reference_support(loaded, sites, variant_class="snv").collect()
    assert answer.select("chrom", "pos").rows() == [("1", 100), ("1", 113), ("1", 900), ("2", 5)]
    again = reference_support(loaded, sites.reverse(), variant_class="snv").collect()
    assert again.equals(answer), "the same positions in any order give the same frame"


def test_an_array_absence_is_never_reference(tmp_path: Path) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(
        "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n"
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
        "# rsid\tchromosome\tposition\tgenotype\n"
        "rs4477212\t1\t82154\tAA\nrs3094315\t1\t752566\tAG\nrs3131972\t1\t752721\t--\n"
    )
    loaded = load(path)
    assert loaded.facets.select("callset_type").collect().item() == CallsetType.ARRAY.value
    answer = _support(loaded, [("1", 82154), ("1", 752721), ("1", 752600)], force=True)
    called = answer[("1", 82154)]
    assert (called["support"], called["reason"], called["grade_basis"], called["confidence_kind"]) == ("call", "call_at_position", "gt_only", "asserted")
    assert called["genotype_index"] is None, "the array's REF is not read, so its genotype indexes no REF"
    no_call = answer[("1", 752721)]
    assert (no_call["support"], no_call["grade_basis"]) == ("call", "no_call")
    _assert_none(answer[("1", 752600)], AbsenceReason.ARRAY)


def test_a_vcf_of_unknown_type_or_an_exome_says_nothing_off_its_rows(tmp_path: Path) -> None:
    loaded = _vcf(tmp_path, [f"1\t100\t.\tC\tT\t50\t.\t.\tGT:GQ:PL:GP\t{SNV_HET}"])
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["callset_type"], facet["absence_means"]) == (None, "unknown"), "two rows earn no callset type"
    _assert_none(_support(loaded, [("1", 500)], force=True)[("1", 500)], AbsenceReason.CALLSET_TYPE_UNKNOWN)
    exome = _as(loaded, callset_type=CallsetType.WES.value, variant_classes=[{"variant_class": "snv", "alleles": 1}])
    answer = _support(exome, [("1", 100), ("1", 500)], force=True)
    _assert_none(answer[("1", 500)], AbsenceReason.TARGETED)
    assert answer[("1", 100)]["support"] == "call" and answer[("1", 100)]["genotype_index"] == 1


def test_a_producers_reference_call_by_omission_is_its_own_absence(tmp_path: Path) -> None:
    """CGI's VCF (`absence_means = reference`): the producer called the omitted position reference."""
    loaded = _vcf(tmp_path, [f"1\t100\t.\tC\tT\t50\t.\t.\tGT:GQ:PL:GP\t{SNV_HET}"])
    cgi = _as(loaded, callset_type="wgs", absence_means="reference", variant_classes=[{"variant_class": "snv", "alleles": 1}])
    answer = _support(cgi, [("1", 500)])
    _assert_absence(answer[("1", 500)], AbsenceReason.OMITTED_CALLED_REFERENCE, AbsenceMeans.REFERENCE)
    _assert_none(_support(cgi, [("1", 500)], variant_class="indel")[("1", 500)], AbsenceReason.CLASS_NOT_CALLED)


@pytest.fixture(scope="module")
def genomes(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Loaded]:
    """Two whole-genome-shaped VCFs loaded for real: one cut to SNVs, one mixed (RM4 types them `wgs`)."""
    root = tmp_path_factory.mktemp("genomes")
    return {"snv_only": load(genome_vcf(root / "sample.vcf", indel_every=None)), "mixed": load(genome_vcf(root / "mixed.vcf", indel_every=8))}


def test_a_whole_genome_silence_reads_as_reference_for_the_classes_it_called(genomes: dict[str, Loaded]) -> None:
    mixed = genomes["mixed"]
    facet = mixed.facets.collect().row(0, named=True)
    assert (facet["callset_type"], facet["absence_means"], facet["content_filter"]) == ("wgs", "weak_reference", None)
    called, absent = ("1", 1), ("1", 2)
    for asked in ("snv", "indel"):
        answer = _support(mixed, [called, absent], variant_class=asked)
        _assert_absence(answer[absent], AbsenceReason.WHOLE_GENOME_SILENCE, AbsenceMeans.WEAK_REFERENCE)
        row = answer[called]
        assert (row["support"], row["grade_basis"], row["confidence_kind"]) == ("call", "gt_only", "asserted")
        assert row["genotype_index"] == 1, "the file's 0/1"
    _assert_none(_support(mixed, [absent], variant_class="sv")[absent], AbsenceReason.CLASS_NOT_CALLED)


def test_a_filtered_genome_needs_force_and_only_for_the_class_it_kept(genomes: dict[str, Loaded]) -> None:
    """The person's 2026-10-22 ruling: a whole genome cut to one class stays `wgs`, and reading its
    absences as reference is additionally gated by `force=True`, for the kept class only."""
    cut = genomes["snv_only"]
    facet = cut.facets.collect().row(0, named=True)
    assert (facet["callset_type"], facet["absence_means"], facet["content_filter"]) == ("wgs", "weak_reference", ContentFilter.SNV_ONLY.value)
    absent = ("1", 2)
    _assert_none(_support(cut, [absent], variant_class="snv")[absent], AbsenceReason.CONTENT_FILTER_HELD)
    _assert_absence(_support(cut, [absent], variant_class="snv", force=True)[absent], AbsenceReason.CONTENT_FILTER_FORCED, AbsenceMeans.WEAK_REFERENCE)
    for force in (False, True):
        _assert_none(_support(cut, [absent], variant_class="indel", force=force)[absent], AbsenceReason.CLASS_CUT_BY_FILTER)
    assert _support(cut, [("1", 1)], variant_class="snv")[("1", 1)]["support"] == "call", "a call answers whatever the filter"
