"""A symbolic `<DEL>`, `<DUP>` or `<INV>` meeting the breakends another caller writes for it (RM282).

truvari 5.5's default `decompose` (`VariantRecord.decompose`, `cpx_match`) turns a symbolic call into
the breakend records it implies and meets them with `bnd_match`. `implied` below restates that
decomposition; `meets` restates `bnd_match` (as tests/test_concordance_bnd.py does). The junctions
themselves are pinned to the VCF spec's own records, independently of truvari:

- the deletion: VCF 4.5 "Encoding Structural Variants", `DEL_split_bp_cn` (chrA, `<DEL>` at 2 with
  SVLEN 2, so END 4; `delbp1` is `T[chrA:5[` at 2, `delbp2` `]chrA:2]A`, whose POS the spec prints as
  2 though its own mate, delbp1's ALT, places it at 5);
- the inversion: VCF 4.3 section 5.4.3 (`<INV>` at 2:321681 with END 421681, and `bnd_W`, `bnd_V`,
  `bnd_U`, `bnd_X`);
- the tandem duplication: no spec writes one as breakends, so its records are derived from VCF 4.5's
  prose ("a breakpoint based duplication call requires that the end of the duplicated region be
  adjacent to the start of the duplicated region") and VCF 4.3's sentence for `t[p[` ("piece extending
  to the right of p is joined after t"), over VCF 4.5's `homology_dup` (`<DUP>` at 5, SVLEN 3: the
  duplicated bases 6 to 8): the base at 8 joined to the base at 6.

The spec's bases are not this fixture's (its sequence is random), so each record is written with the
fixture's base at its position; the forms and the mate positions are the spec's.
"""

import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, SvMatching, load, resolve_ref
from just_vcf.concordance import SITES, concordance
from just_vcf.loading import Loaded

RNG = random.Random(2822)
SEQ = {"1": "".join(RNG.choice("ACGT") for _ in range(30_000)), "2": "".join(RNG.choice("ACGT") for _ in range(430_000))}

META = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
)

BNDDIST = SvMatching().bnddist

SPEC_DEL = ("1", 2, 4)
SPEC_DEL_RECORDS = [("1", 2, "t[p[", "1", 5), ("1", 5, "]p]t", "1", 2)]
"""VCF 4.5's `delbp1` and `delbp2` (the latter at the POS its mate states)."""

SPEC_INV = ("2", 321681, 421681)
SPEC_INV_RECORDS = [
    ("2", 321681, "t]p]", "2", 421681),  # bnd_W
    ("2", 321682, "[p[t", "2", 421682),  # bnd_V
    ("2", 421681, "t]p]", "2", 321681),  # bnd_U
    ("2", 421682, "[p[t", "2", 321682),  # bnd_X
]
"""VCF 4.3 section 5.4.3's breakends for `<INV>` 2:321681-421681."""


def join_after(left: int, right: int, chrom: str = "1") -> list[tuple]:
    """The two records VCF 4.3 writes for the base at `left` followed by the base at `right`: `t[p[` at
    `left` (the piece right of `right` joined after `t`) and its mate `]p]t` at `right`."""
    return [(chrom, left, "t[p[", chrom, right), (chrom, right, "]p]t", chrom, left)]


SPEC_DUP = ("1", 5, 8)
SPEC_DUP_RECORDS = join_after(8, 6)
"""VCF 4.5's `homology_dup` as breakends: the end of the duplicated bases (8) adjacent to their start (6)."""


def base(chrom: str, pos: int) -> str:
    return SEQ[chrom][pos - 1]


def alt(form: str, t: str, mate_chrom: str, mate_pos: int) -> str:
    return form.replace("t", t).replace("p", f"{mate_chrom}:{mate_pos}")


def info_text(fields: dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in fields.items()) or "."


def bnd(record: tuple, gt: str = "0/1", **info: str) -> str:
    chrom, pos, form, mate_chrom, mate_pos = record
    t = base(chrom, pos)
    return f"{chrom}\t{pos}\t.\t{t}\t{alt(form, t, mate_chrom, mate_pos)}\t.\tPASS\t{info_text({'SVTYPE': 'BND', **info})}\tGT\t{gt}"


def sv(call: tuple, symbol: str, gt: str = "0/1", svtype: str | None = None, **info: str) -> str:
    chrom, pos, end = call
    fields = {"SVTYPE": svtype or symbol.strip("<>").split(":")[0], "END": str(end), **info}
    return f"{chrom}\t{pos}\t.\t{base(chrom, pos)}\t{symbol}\t.\tPASS\t{info_text(fields)}\tGT\t{gt}"


def implied(kind: str, call: tuple, cipos_hi: int = 0, ciend_hi: int = 0) -> list[tuple]:
    """truvari 5.5 `VariantRecord.decompose`, restated: the breakend records a symbolic call implies,
    as `(chrom, pos, form, mate_chrom, mate_pos, cipos_hi, ciend_hi)`. Each is a copy of the symbolic
    record, so each carries its CIPOS and CIEND, whichever end it stands at."""
    chrom, p, e = call
    shapes = {
        "INV": [(p, "t]p]", e), (p + 1, "[p[t", e + 1), (e, "t]p]", p), (e + 1, "[p[t", p + 1)],
        "DEL": [(p, "t[p[", e + 1), (e + 1, "]p]t", p)],
        "DUP": [(p + 1, "]p]t", e), (e, "t[p[", p + 1)],
    }
    return [(chrom, pos, form, chrom, mate, cipos_hi, ciend_hi) for pos, form, mate in shapes[kind]]


def meets(a: tuple, b: tuple, bnddist: int = BNDDIST) -> bool:
    """truvari 5.5 `bnd_match`, restated over `(chrom, pos, form, mate_chrom, mate_pos, cipos_hi, ciend_hi)`."""

    def overlaps(p1: int, c1: int, p2: int, c2: int) -> bool:
        return max(p1 - bnddist - c1, p2 - bnddist - c2) < min(p1 + bnddist + c1, p2 + bnddist + c2)

    return a[0] == b[0] and overlaps(a[1], a[5], b[1], b[5]) and a[3] == b[3] and overlaps(a[4], a[6], b[4], b[6]) and a[2] == b[2]


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n" for name, seq in SEQ.items()))
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(fasta)
    return s


def facet(tmp_path: Path, name: str, records: list[str], store: ReferenceStore) -> Loaded:
    meta = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items()), *META]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    path = tmp_path / f"{name}.vcf"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([*meta, header, *records]) + "\n")
    return resolve_ref(load(path), store)


def compare(tmp_path: Path, store: ReferenceStore, a: list[str], b: list[str], **kwargs):
    return concordance(facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store), **kwargs)


def sv_rows(tmp_path: Path, store: ReferenceStore, a: list[str], b: list[str], **kwargs) -> pl.DataFrame:
    return compare(tmp_path, store, a, b, **kwargs).sites.collect().filter(pl.col("variant_class") == "sv")


SPEC_CASES = {
    "DEL": ("<DEL>", SPEC_DEL, SPEC_DEL_RECORDS),
    "DUP": ("<DUP>", SPEC_DUP, SPEC_DUP_RECORDS),
    "INV": ("<INV>", SPEC_INV, SPEC_INV_RECORDS),
}


SHOWN = {"DEL": SPEC_DEL_RECORDS[0], "DUP": SPEC_DUP_RECORDS[1], "INV": SPEC_INV_RECORDS[0]}
"""The breakend record a pair shows: of the records that met best (each spec record meets an implied
one exactly), the one at its adjacency's first end; of an inversion's two adjacencies, the one its
first implied breakend (at POS) meets."""


@pytest.mark.parametrize("kind", sorted(SPEC_CASES))
def test_truvaris_decomposition_is_the_specs_breakends(kind: str) -> None:
    # The restatement of truvari is checked against the spec's records before anything leans on it.
    _, call, records = SPEC_CASES[kind]
    assert sorted(r[:5] for r in implied(kind, call)) == sorted(records)


@pytest.mark.parametrize("kind", sorted(SPEC_CASES))
def test_a_symbolic_call_meets_the_breakends_the_spec_writes_for_it(tmp_path: Path, store: ReferenceStore, kind: str) -> None:
    symbol, call, records = SPEC_CASES[kind]
    rows = sv_rows(tmp_path, store, [sv(call, symbol)], [bnd(r) for r in records])
    met = rows.filter(pl.col("key") == "decomposed")
    # Taken once (`--pick single`): one row for the symbolic call, the breakend shown is the record that
    # met best, every record of the spec meeting one of the implied ones exactly.
    assert met.height == 1, rows
    row = met.row(0, named=True)
    first = SHOWN[kind]
    assert (row["chrom"], row["pos"], row["alt"], row["sv_kind"], row["agreement"], row["absence"]) == (call[0], call[1], symbol, kind, "match", None)
    assert (row["b_pos"], row["b_alt"]) == (first[1], alt(first[2], base(first[0], first[1]), first[3], first[4]))
    assert (row["a_end_pos"], row["b_end_pos"], row["ci_needed"]) == (call[2], first[1], False)
    assert row["size_similarity"] is None and row["reciprocal_overlap"] is None
    # An inversion is two adjacencies: the one it did not take is near one it implies, so `b`'s silence
    # reads `other_alleles`, never `unknown`.
    alone = rows.filter(pl.col("key").is_null())
    adjacencies = 2 if kind == "INV" else 1
    assert alone.select("agreement", "absence").rows() == [("only_b", "other_alleles")] * (adjacencies - 1)


@pytest.mark.parametrize("kind", sorted(SPEC_CASES))
def test_each_implied_breakend_lies_exactly_where_the_spec_puts_it(tmp_path: Path, store: ReferenceStore, kind: str) -> None:
    # Each spec record alone, beside four of its neighbours one base off (its position or its mate's):
    # all meet within `bnddist`, and the exact one scores best (truvari's score falls with each base of
    # distance), so only an implied breakend at the spec's coordinates shows it.
    symbol, call, records = SPEC_CASES[kind]
    for index, record in enumerate(records):
        chrom, pos, form, mate_chrom, mate_pos = record
        off = [(chrom, pos + d, form, mate_chrom, mate_pos) for d in (-1, 1)] + [(chrom, pos, form, mate_chrom, mate_pos + d) for d in (-1, 1)]
        rows = sv_rows(tmp_path / f"{index}", store, [sv(call, symbol)], [bnd(r) for r in [*off[:2], record, *off[2:]]])
        met = rows.filter(pl.col("key") == "decomposed")
        assert met["b_alt"].to_list() == [alt(form, base(chrom, pos), mate_chrom, mate_pos)], (record, met)
        assert met["b_pos"].to_list() == [pos]


def test_a_tandem_duplication_subtype_is_decomposed_and_another_is_not(tmp_path: Path, store: ReferenceStore) -> None:
    records = [bnd(r) for r in SPEC_DUP_RECORDS]
    tandem = sv_rows(tmp_path / "tandem", store, [sv(SPEC_DUP, "<DUP:TANDEM>")], records)
    assert tandem["key"].to_list() == ["decomposed"]
    assert tandem["sv_kind"].to_list() == ["DUP"]
    # A dispersed duplication (VCF 4.5) puts its copy elsewhere: its end is not joined to its start.
    dispersed = sv_rows(tmp_path / "dispersed", store, [sv(SPEC_DUP, "<DUP:DISPERSED>", svtype="DUP")], records)
    assert dispersed["key"].is_null().all() and sorted(dispersed["agreement"].to_list()) == ["only_a", "only_b"]


def test_a_copy_number_call_is_not_decomposed(tmp_path: Path, store: ReferenceStore) -> None:
    # DRAGEN's copy-number caller writes `<DEL>` with SVTYPE=CNV: an abundance claim, no junction (VCF
    # 4.5's SVCLAIM=D). truvari reads SVTYPE first and does not decompose it either.
    result = compare(tmp_path, store, [sv(SPEC_DEL, "<DEL>", svtype="CNV")], [bnd(r) for r in SPEC_DEL_RECORDS])
    rows = result.sites.collect()
    assert rows["key"].is_null().all() and sorted(rows["agreement"].to_list()) == ["only_a", "only_b"]
    assert any("not decomposed into breakends" in line for line in result.warnings), result.warnings


@pytest.mark.parametrize(
    "shift",
    [0, 2 * BNDDIST - 1, 2 * BNDDIST, -(2 * BNDDIST - 1), -(2 * BNDDIST)],
)
def test_implied_breakends_meet_within_truvaris_breakend_distance(tmp_path: Path, store: ReferenceStore, shift: int) -> None:
    call = ("1", 5000, 6000)
    # b's junction of the deletion, its position moved by `shift`, its mate where the deletion ends.
    record = ("1", 5000 + shift, "t[p[", "1", 6001)
    rows = sv_rows(tmp_path, store, [sv(call, "<DEL>")], [bnd(record)])
    expected = any(meets(r, (*record, 0, 0)) for r in implied("DEL", call))
    if expected:
        assert rows.select("key", "a_pos", "b_pos").rows() == [("decomposed", 5000, 5000 + shift)]
    else:
        assert sorted(rows["agreement"].to_list()) == ["only_a", "only_b"]


def test_another_join_at_the_same_ends_does_not_meet_and_reads_other_alleles(tmp_path: Path, store: ReferenceStore) -> None:
    # The ends a deletion joins, written in a duplication's forms (the piece right of 5000 after the
    # piece left of 6001): same positions, another adjacency.
    call = ("1", 5000, 6000)
    records = [("1", 5000, "]p]t", "1", 6001), ("1", 6001, "t[p[", "1", 5000)]
    assert not any(meets(r, (*s, 0, 0)) for r in implied("DEL", call) for s in records)
    rows = sv_rows(tmp_path, store, [sv(call, "<DEL>")], [bnd(r) for r in records])
    # Neither meets the other, and each facet's silence is read beside the other's SVs: other_alleles.
    assert sorted(rows.select("agreement", "absence").rows()) == [("only_a", "other_alleles"), ("only_b", "other_alleles")]


def test_confidence_intervals_widen_as_truvaris_copies_carry_them(tmp_path: Path, store: ReferenceStore) -> None:
    call = ("1", 5000, 6000)
    far = 5000 + 2 * BNDDIST + 5
    cases = {
        # The symbolic record's CIPOS widens each implied breakend's own position.
        "cipos": ({"CIPOS": "-10,10"}, (10, 0), ("1", far, "t[p[", "1", 6001), True),
        # Its CIEND widens the mate's, at index 1 as truvari reads it; the record at END+1 is a copy too,
        # so it widens its own position by CIPOS, not CIEND, and its mate (POS) by CIEND.
        "ciend": ({"CIEND": "-10,10"}, (0, 10), ("1", 5000, "t[p[", "1", 6001 + 2 * BNDDIST + 5), True),
        "none": ({}, (0, 0), ("1", far, "t[p[", "1", 6001), False),
    }
    for name, (info, (ci_pos, ci_end), record, want) in cases.items():
        rows = sv_rows(tmp_path / name, store, [sv(call, "<DEL>", **info)], [bnd(record)])
        expected = any(meets(r, (*record, 0, 0)) for r in implied("DEL", call, ci_pos, ci_end))
        assert expected is want, name
        assert (rows["key"].to_list() == ["decomposed"]) is expected, (name, rows)
        if expected:
            assert rows["ci_needed"].to_list() == [True]


def test_a_pair_of_one_representation_is_taken_before_decomposition(tmp_path: Path, store: ReferenceStore) -> None:
    call = ("1", 5000, 6000)
    records = [bnd(r[:5]) for r in implied("DEL", call)]
    # b wrote the deletion twice, as a symbolic call and as breakends: a's `<DEL>` meets b's `<DEL>` by
    # extent, and b's breakends, left alone, are near a's implied ones.
    rows = sv_rows(tmp_path / "both", store, [sv(call, "<DEL>")], [sv(call, "<DEL>"), *records])
    assert rows.select("key", "agreement", "absence").sort("key", nulls_last=True).rows() == [("extent", "match", None), (None, "only_b", "other_alleles")]
    # Likewise a's breakends meet b's breakends, and b's `<DEL>` is left alone beside them.
    rows = sv_rows(tmp_path / "bnd", store, records, [sv(call, "<DEL>"), *records])
    assert rows.select("key", "agreement", "absence").sort("key", nulls_last=True).rows() == [("breakend", "match", None), (None, "only_b", "other_alleles")]


def test_decomposed_pairs_compare_genotypes_and_swap_with_the_facets(tmp_path: Path, store: ReferenceStore) -> None:
    call = ("1", 5000, 6000)
    records = [bnd(r[:5], gt="1/1") for r in implied("DEL", call)]
    rows = sv_rows(tmp_path / "ab", store, [sv(call, "<DEL>", gt="0/1")], records)
    assert rows.select("key", "sv_kind", "agreement", "a_copies", "b_copies").rows() == [("decomposed", "DEL", "mismatch", 1, 2)]
    # Swapped, `a` is the breakend side: the row stands at its record, the `<DEL>` is `b_alt`, and the
    # type matched is still the deletion's.
    swapped = sv_rows(tmp_path / "ba", store, records, [sv(call, "<DEL>", gt="0/1")])
    row = swapped.row(0, named=True)
    t = base("1", 5000)
    assert swapped.height == 1
    assert (row["key"], row["sv_kind"], row["alt"], row["b_alt"], row["a_pos"], row["b_pos"]) == ("decomposed", "DEL", alt("t[p[", t, "1", 6001), "<DEL>", 5000, 5000)
    assert (row["agreement"], row["a_copies"], row["b_copies"], row["a_end_pos"], row["b_end_pos"]) == ("mismatch", 2, 1, 5000, 6000)


def test_decomposed_rows_join_the_sites_frame_and_its_summary(tmp_path: Path, store: ReferenceStore) -> None:
    calls = [("1", 5000, 6000), ("1", 15_000, 15_400)]
    a = [sv(calls[0], "<DEL>"), sv(calls[1], "<INV>")]
    b = [*(bnd(r[:5]) for r in implied("DEL", calls[0])), *(bnd(r[:5]) for r in implied("INV", calls[1]))]
    result = compare(tmp_path / "x", store, a, b)
    sites = result.sites.collect()
    assert sites.schema == SITES
    counts = {(row["key"], row["sv_kind"]): row["alleles"] for row in result.summary.iter_rows(named=True)}
    assert counts == {("decomposed", "DEL"): 1, ("decomposed", "INV"): 1, (None, "BND"): 1}
    assert sum(counts.values()) == sites.height
    # The same inputs in another order give the same rows.
    again = compare(tmp_path / "y", store, list(reversed(a)), list(reversed(b))).sites.collect()
    assert again.equals(sites)


def test_a_lone_breakend_near_no_implied_breakend_keeps_the_absence_rule(tmp_path: Path, store: ReferenceStore) -> None:
    # Far from the deletion's junctions (more than 2 * bnddist from both), and an `<INS>`, which implies
    # no breakend, beside it: b's silence is not `other_alleles`.
    lone = ("1", 20_000, "t[p[", "2", 9000)
    rows = sv_rows(tmp_path, store, [bnd(lone)], [sv(("1", 5000, 6000), "<DEL>"), sv(("1", 20_000, 20_000), "<INS>")])
    assert rows.filter(pl.col("agreement") == "only_a")["absence"].to_list() != ["other_alleles"]
