"""The aggregate-only query helper (RM283): answers that carry no row, no sample name and no path."""

import json
import statistics
from collections import Counter
from pathlib import Path

import polars as pl
import pytest
from pydantic import BaseModel, ValidationError

from just_vcf import load, mint_vrs_ids, resolve_ref
from just_vcf.formats import get_spec
from just_vcf.loading import Loaded
from just_vcf.query import (
    CountGroup,
    CountTable,
    CountUnit,
    Dimension,
    FacetProvenance,
    Frame,
    GenotypeClass,
    NumericSummary,
    Presence,
    PresenceBasis,
    Query,
    StatTable,
    WarningGroup,
    query,
)
from just_vcf.schema import AbsenceMeans, BlockKind, VariantClass
from test_load_vcf import vcf_text
from test_reference import SEQ, base, contig_lines, other, store, vcf  # noqa: F401  (store is a fixture)

SAMPLE = "XQZSAMPLE77"
DIRECTORY = "XQZDIR77"
FILE = "XQZFILE77.vcf"

# Positions picked to be distinctive strings, so an answer that leaks one is found by a substring search.
SNV_HET = 1234567
SNV_HOM = 2345678
SNV_HOMREF = 3456789
BLOCK_REF = 4567890
BLOCK_NOCALL = 5678901
DELETION = 6789012
MULTI = 7890123
NOCALL = 8901234
PARTIAL = 9012345
POSITIONS = [SNV_HET, SNV_HOM, SNV_HOMREF, BLOCK_REF, BLOCK_NOCALL, DELETION, MULTI, NOCALL, PARTIAL]

RECORDS = [
    f"1\t{SNV_HET}\trs7001\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:50:30",
    f"1\t{SNV_HOM}\trs7002\tC\tT\t60\tPASS\t.\tGT:GQ:DP\t1/1:99:40",
    f"1\t{SNV_HOMREF}\t.\tG\tA\t10\tLowQual\t.\tGT:GQ:DP\t0/0:12:8",
    f"1\t{BLOCK_REF}\t.\tA\t<*>\t.\t.\tEND={BLOCK_REF + 100}\tGT:GQ:DP\t0/0:30:20",
    f"1\t{BLOCK_NOCALL}\t.\tC\t<*>\t.\t.\tEND={BLOCK_NOCALL + 100}\tGT:GQ:DP\t./.:0:0",
    f"2\t{DELETION}\t.\tAT\tA\t40\tPASS\t.\tGT:GQ:DP\t0/1:45:25",
    f"2\t{MULTI}\t.\tA\tC,T\t70\tPASS\t.\tGT:GQ:DP\t1/2:80:35",
    f"2\t{NOCALL}\t.\tC\tG\t.\t.\t.\tGT:GQ:DP\t./.:.:.",
    f"2\t{PARTIAL}\t.\tT\tC\t20\tPASS;LowQual\t.\tGT:GQ:DP\t0/.:5:3",
]


@pytest.fixture
def path(tmp_path: Path) -> Path:
    directory = tmp_path / DIRECTORY
    directory.mkdir()
    target = directory / FILE
    target.write_text(vcf_text(RECORDS, samples=(SAMPLE,), fmt=("GT", "GQ", "DP"), info=("END",)))
    return target


@pytest.fixture
def loaded(path: Path) -> Loaded:
    return load(path)


def everything(q: Query) -> list[BaseModel]:
    """Every answer the helper gives on the fixture: each dimension, each statistic, each presence path."""
    out: list[BaseModel] = [*q.facets(), *q.warnings()]
    for dimension in Dimension:
        frame = Frame.BLOCKS if dimension in (Dimension.BLOCK_KIND,) else Frame.CALLS
        out += q.counts(by=[dimension], frame=frame)
    out += q.counts(by=[Dimension.CHROM, Dimension.GENOTYPE_CLASS, Dimension.FILTER])
    out += q.counts(by=[Dimension.CHROM, Dimension.BLOCK_KIND], frame=Frame.BLOCKS)
    out += q.counts(by=[])
    out += q.stats("gq") + q.stats("dp", by=["chrom"]) + q.stats("qual") + q.stats("dp_min") + q.stats("gq", frame="blocks")
    return out


def dumped(answers: list[BaseModel]) -> str:
    return "\n".join(a.model_dump_json() for a in answers)


# --- what never comes out -------------------------------------------------------------------------


def test_no_answer_carries_the_sample_the_path_or_a_position(loaded: Loaded) -> None:
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["sample"] == SAMPLE and DIRECTORY in facet["source_path"]  # the fixture does hold them
    text = dumped(everything(query(loaded)))
    assert SAMPLE not in text
    assert DIRECTORY not in text and FILE not in text
    for position in POSITIONS:
        assert str(position) not in text


def test_a_presence_answer_names_no_position_but_the_one_asked(loaded: Loaded) -> None:
    q = query(loaded)
    for asked in POSITIONS:
        answers = q.carries(chrom="1" if asked < DELETION else "2", pos=asked, ref="A", alt="G")
        text = dumped(list(answers))
        assert SAMPLE not in text and DIRECTORY not in text
        assert [p for p in POSITIONS if str(p) in text] == []  # not even the asked one is echoed back


def test_warnings_lose_the_sample_and_any_path(loaded: Loaded, path: Path) -> None:
    # The loaders name facets, never samples; this checks the helper's own guard, on lines that break the rule.
    lines = (
        f"3 row(s) in {path} dropped",
        f"2 row(s) in {path} dropped",
        f"the header names {SAMPLE} twice; sample_meta left as written",
        f"read from {path.name}, see /home/someone/genomes/raw.vcf",
    )
    q = query(Loaded(facets=loaded.facets, calls=loaded.calls, blocks=loaded.blocks, warnings=lines))
    groups = q.warnings()
    text = dumped(list(groups))
    assert SAMPLE not in text and DIRECTORY not in text and FILE not in text and "someone" not in text
    assert groups[0] == WarningGroup(reason="row(s) in <path> dropped", count=5, lines=2)
    assert groups[1].count is None and "<sample>" in groups[1].reason
    assert groups[2].reason == "read from <file name>, see <path>"


FORBIDDEN_FIELDS = {"sample", "source_path", "member", "sample_meta", "source_sha256", "pos", "end_pos", "chrom", "ref", "alt", "alts", "ids", "gt", "gt_alleles", "source_gt", "source_ref", "ploidy_statements", "header_contigs"}


def _field_names(schema: dict) -> set[str]:
    names: set[str] = set()
    for definition in [schema, *schema.get("$defs", {}).values()]:
        names |= set(definition.get("properties", {}))
    return names


@pytest.mark.parametrize("model", [CountGroup, CountTable, NumericSummary, StatTable, FacetProvenance, WarningGroup, Presence])
def test_no_answer_type_has_a_field_for_a_row_or_a_person(model: type[BaseModel]) -> None:
    schema = model.model_json_schema()
    assert _field_names(schema) & FORBIDDEN_FIELDS == set()
    assert model.model_config.get("extra") == "forbid" and model.model_config.get("frozen") is True
    # No field accepts an arbitrary object, where a row could be put.
    assert '"additionalProperties": true' not in json.dumps(schema)
    assert "{}" not in json.dumps(schema)


def test_an_answer_refuses_a_field_it_does_not_declare() -> None:
    with pytest.raises(ValidationError):
        Presence(facet_id="f", carried=None, basis=PresenceBasis.ABSENCE, sample=SAMPLE)  # type: ignore[call-arg]
    with pytest.raises(ValidationError):
        CountGroup(values=("x" * 200,), count=1)


# --- counts and statistics ------------------------------------------------------------------------


def _expected_genotype_class(gt: list[int | None] | None, alleles: list[str | None] | None) -> str:
    """The test's own reading of a genotype, slot by slot."""
    if gt is None:
        return "not_stated" if alleles is None else "unindexed"
    called = [g for g in gt if g is not None]
    if not called:
        return "no_call"
    if len(called) < len(gt):
        return "partial"
    if len(gt) == 1:
        return "haploid_ref" if gt[0] == 0 else "haploid_alt"
    if all(g == 0 for g in gt):
        return "hom_ref"
    if len({g for g in gt if g != 0}) >= 2:
        return "het_alt"
    return "het" if 0 in gt else "hom_alt"


def test_counts_by_contig_and_genotype_class_equal_the_rows(loaded: Loaded) -> None:
    calls = loaded.calls.collect()
    expected = Counter((row["chrom"], _expected_genotype_class(row["gt"], row["gt_alleles"])) for row in calls.iter_rows(named=True))
    (table,) = query(loaded).counts(by=["chrom", "genotype_class"])
    assert table.unit is CountUnit.ROWS and table.total == calls.height
    assert {g.values: g.count for g in table.groups} == dict(expected)
    assert {g.values[1] for g in table.groups} >= {"het", "hom_alt", "hom_ref", "het_alt", "no_call", "partial"}
    # Fixed order: contigs by rank, then the enum's own order.
    order = [(g.values[0], list(GenotypeClass).index(GenotypeClass(g.values[1]))) for g in table.groups]
    assert order == sorted(order)


def test_counts_by_variant_class_count_alt_alleles(loaded: Loaded) -> None:
    calls = loaded.calls.collect()
    (table,) = query(loaded).counts(by="variant_class")
    alleles = sum(len(a) for a in calls["alts"].to_list() if a is not None)
    assert table.unit is CountUnit.ALT_ALLELES and table.total == alleles
    facet = loaded.facets.collect().row(0, named=True)
    stated = {e["variant_class"]: e["alleles"] for e in facet["variant_classes"]}
    assert {g.values[0]: g.count for g in table.groups if g.values[0] is not None} == stated
    assert table.rows_without_alts == calls.filter(pl.col("alts").is_null()).height


def test_counts_by_filter_keep_the_filter_set_and_the_unstated(loaded: Loaded) -> None:
    calls = loaded.calls.collect()
    expected = Counter(None if f is None else ";".join(sorted(f)) for f in calls["filters"].to_list())
    (table,) = query(loaded).counts(by="filter")
    assert {g.values[0]: g.count for g in table.groups} == dict(expected)
    assert "LowQual;PASS" in expected and None in expected  # the fixture has both cases


def test_block_counts_report_their_bases(loaded: Loaded) -> None:
    blocks = loaded.blocks.collect()
    (table,) = query(loaded).counts(by=["block_kind"], frame="blocks")
    assert table.unit is CountUnit.BLOCKS
    for group in table.groups:
        rows = blocks.filter(pl.col("kind") == group.values[0])
        assert group.count == rows.height
        assert group.bases == sum(e - p + 1 for p, e in zip(rows["pos"], rows["end_pos"], strict=True))
    assert {g.values[0] for g in table.groups} == {BlockKind.REF.value, BlockKind.NO_CALL.value}


def test_a_dimension_the_frame_lacks_is_refused(loaded: Loaded) -> None:
    q = query(loaded)
    with pytest.raises(ValueError, match="cannot be grouped"):
        q.counts(by=["genotype_class"], frame="blocks")
    with pytest.raises(ValueError):
        q.counts(by=["pos"])
    with pytest.raises(ValueError, match="has no"):
        q.stats("dp_min", frame="calls")


def test_stats_summarise_the_column(loaded: Loaded) -> None:
    calls = loaded.calls.collect()
    gq = [v for v in calls["gq"].to_list() if v is not None]
    ((summary,),) = [t.groups for t in query(loaded).stats("gq")]
    assert summary.present == len(gq) and summary.missing == calls.height - len(gq)
    assert summary.minimum == min(gq) and summary.maximum == max(gq)
    assert summary.median == statistics.median(gq)
    assert summary.mean == pytest.approx(statistics.mean(gq))


def test_stats_of_an_empty_column_are_null_not_zero(tmp_path: Path) -> None:
    path = tmp_path / "nogq.vcf"
    path.write_text(vcf_text(["1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"]))
    ((summary,),) = [t.groups for t in query(path).stats("gq")]
    assert summary.present == 0 and summary.missing == 1
    assert (summary.minimum, summary.median, summary.maximum, summary.mean) == (None, None, None, None)


def test_provenance_keeps_the_facet_and_drops_the_person(loaded: Loaded) -> None:
    facet = loaded.facets.collect().row(0, named=True)
    (provenance,) = query(loaded).facets()
    assert provenance.facet_id == facet["facet_id"]
    assert provenance.format_id == facet["format_id"]
    assert provenance.absence_means == AbsenceMeans(facet["absence_means"])
    assert provenance.contigs == tuple(sorted({e["chrom"] for e in facet["contig_map"]}))
    assert {c.variant_class for c in provenance.variant_classes or ()} == {VariantClass(e["variant_class"]) for e in facet["variant_classes"]}


# --- presence -------------------------------------------------------------------------------------


def test_a_stated_allele_answers_with_its_genotype(loaded: Loaded) -> None:
    q = query(loaded)
    (het,) = q.carries(chrom="chr1", pos=SNV_HET, ref="a", alt="g")
    assert het and het.carried is True and het.basis is PresenceBasis.STATED
    assert (het.copies, het.ploidy, het.genotype_class) == (1, 2, GenotypeClass.HET)
    assert (het.gq, het.dp, het.filters_passed) == (50, 30, True)
    (hom,) = q.carries(chrom="1", pos=SNV_HOM, ref="C", alt="T")
    assert (hom.carried, hom.copies, hom.genotype_class) == (True, 2, GenotypeClass.HOM_ALT)
    (homref,) = q.carries(chrom="1", pos=SNV_HOMREF, ref="G", alt="A")
    assert not homref and homref.carried is False and homref.basis is PresenceBasis.STATED
    assert (homref.genotype_class, homref.filters_passed) == (GenotypeClass.HOM_REF, False)
    (second,) = q.carries(chrom="2", pos=MULTI, ref="A", alt="T")
    assert (second.carried, second.copies, second.genotype_class) == (True, 1, GenotypeClass.HET_ALT)
    (indel,) = q.carries(chrom="2", pos=DELETION, ref="AT", alt="A")
    assert indel.carried is True and indel.copies == 1


def test_a_no_call_or_half_call_is_unknown_not_absent(loaded: Loaded) -> None:
    q = query(loaded)
    (nocall,) = q.carries(chrom="2", pos=NOCALL, ref="C", alt="G")
    assert nocall.carried is None and nocall.basis is PresenceBasis.STATED and nocall.genotype_class is GenotypeClass.NO_CALL
    (partial,) = q.carries(chrom="2", pos=PARTIAL, ref="T", alt="C")
    assert partial.carried is None and partial.genotype_class is GenotypeClass.PARTIAL
    # The site genotyped without the asked allele: not carried only where every haplotype was called.
    (other_snv,) = q.carries(chrom="1", pos=SNV_HET, ref="A", alt="T")
    assert other_snv.carried is False and other_snv.basis is PresenceBasis.SITE_GENOTYPED
    (half,) = q.carries(chrom="2", pos=PARTIAL, ref="T", alt="G")
    assert half.carried is None and half.basis is PresenceBasis.SITE_GENOTYPED


def test_blocks_answer_for_a_position_they_cover(loaded: Loaded) -> None:
    q = query(loaded)
    (in_ref,) = q.carries(chrom="1", pos=BLOCK_REF + 50, ref="A", alt="G")
    assert in_ref.carried is False and in_ref.basis is PresenceBasis.REF_BLOCK and in_ref.gq == 30
    (indel_in_ref,) = q.carries(chrom="1", pos=BLOCK_REF + 50, ref="A", alt="AT")
    assert indel_in_ref.carried is None and indel_in_ref.basis is PresenceBasis.REF_BLOCK
    (in_nocall,) = q.carries(chrom="1", pos=BLOCK_NOCALL + 50, ref="A", alt="G")
    assert in_nocall.carried is None and in_nocall.basis is PresenceBasis.NO_CALL_BLOCK


def test_an_snv_inside_a_deletion_is_unknown(loaded: Loaded) -> None:
    (inside,) = query(loaded).carries(chrom="2", pos=DELETION + 1, ref="T", alt="G")
    assert inside.carried is None and inside.basis is PresenceBasis.OVERLAPPING_CALL


def test_absence_is_read_through_absence_means(loaded: Loaded) -> None:
    (gvcf,) = query(loaded).carries(chrom="1", pos=100, ref="A", alt="G")
    assert gvcf.carried is None and gvcf.basis is PresenceBasis.ABSENCE and gvcf.absence_means is AbsenceMeans.BLOCKS
    # The same rows as a whole-genome callset without blocks: its silence is weak reference for an SNV.
    wgs = Loaded(facets=loaded.facets.with_columns(pl.lit("weak_reference").cast(pl.Enum([m.value for m in AbsenceMeans])).alias("absence_means")), calls=loaded.calls, blocks=loaded.blocks.clear())
    q = query(wgs)
    (snv,) = q.carries(chrom="1", pos=100, ref="A", alt="G")
    assert snv.carried is False and snv.basis is PresenceBasis.ABSENCE and snv.absence_means is AbsenceMeans.WEAK_REFERENCE
    # An indel the facet does not state may be written elsewhere in another spelling: unknown.
    (indel,) = q.carries(chrom="1", pos=100, ref="A", alt="AC")
    assert indel.carried is None and indel.basis is PresenceBasis.ABSENCE


def test_rows_that_disagree_are_not_resolved(tmp_path: Path) -> None:
    path = tmp_path / "twice.vcf"
    path.write_text(vcf_text(["1\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", "1\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0/0"]))
    loaded = load(path)
    assert loaded.calls.collect().height == 2
    (answer,) = query(loaded).carries(chrom="1", pos=500, ref="A", alt="G")
    assert answer.carried is None and answer.basis is PresenceBasis.ROWS_DISAGREE and answer.rows == 2


def test_an_array_before_resolve_ref_is_unresolved(tmp_path: Path) -> None:
    spec = get_spec("23andme")
    path = tmp_path / "genome.txt"
    path.write_text(spec.example)
    loaded = load(path)
    calls = loaded.calls.collect()
    row = calls.filter(pl.col("gt_alleles").is_not_null() & (pl.col("pos") > 0)).row(0, named=True)
    (answer,) = query(loaded).carries(chrom=row["chrom"], pos=row["pos"], ref="A", alt="G")
    assert answer.carried is None and answer.basis is PresenceBasis.ALLELES_UNRESOLVED
    (table,) = query(loaded).counts(by="genotype_class")
    expected = Counter(_expected_genotype_class(r["gt"], r["gt_alleles"]) for r in calls.iter_rows(named=True))
    assert {g.values[0]: g.count for g in table.groups} == dict(expected)
    assert "unindexed" in expected


def test_a_vrs_id_finds_its_row_and_an_unknown_id_says_why(tmp_path: Path, store) -> None:  # noqa: F811
    het, hom = 200, 900
    records = [("1", het, base("1", het), other(base("1", het)), "0/1"), ("1", hom, base("1", hom), other(base("1", hom)), "1/1")]
    minted = mint_vrs_ids(resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": len(SEQ["1"])}))), store), store)
    ids = minted.calls.collect()["vrs_allele_ids"].to_list()
    q = query(minted)
    (by_id,) = q.carries(vrs_id=ids[1][1])
    assert by_id.carried is True and by_id.copies == 2 and by_id.basis is PresenceBasis.STATED
    (both,) = q.carries(vrs_id=ids[0][1], chrom="1", pos=het, ref=base("1", het), alt=other(base("1", het)))
    assert both.carried is True and both.copies == 1 and both.rows == 1
    (unknown,) = q.carries(vrs_id="ga4gh:VA." + "A" * 32)
    assert unknown.carried is None and unknown.basis is PresenceBasis.NOT_FOUND_BY_ID
    # The REF slot's id is not an ALT allele: it is not looked up.
    (ref_slot,) = q.carries(vrs_id=ids[0][0])
    assert ref_slot.basis is PresenceBasis.NOT_FOUND_BY_ID


def test_a_malformed_question_is_refused(loaded: Loaded) -> None:
    q = query(loaded)
    with pytest.raises(ValueError):
        q.carries()
    with pytest.raises(ValueError):
        q.carries(chrom="1", pos=100)
    with pytest.raises(ValueError):
        q.carries(vrs_id="rs123")
    with pytest.raises(ValueError):
        q.carries(chrom="1", pos=100, ref="A", alt="A")
    with pytest.raises(ValueError):
        q.carries(chrom="1", pos=100, ref="A", alt="<DEL>")


def test_query_reads_a_file_itself_and_shows_no_rows_in_its_repr(path: Path) -> None:
    q = query(path)
    assert SAMPLE not in repr(q) and DIRECTORY not in repr(q)
    assert q.facet_ids == tuple(load(path).facets.collect()["facet_id"].to_list())
