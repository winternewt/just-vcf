"""Concordance of breakends between two facets (RM282): truvari bench's `bnddist` matching.

Each case writes two one-sample VCFs of breakend records over a small refget store, reads REF and
compares. The rule mirrored is truvari 5.5's `VariantRecord.bnd_match`; `meets` below restates it from
truvari's source, independently of the library, and every expectation is computed from it and the
coordinates the fixture writes:

- the two records lie on one contig, and their positions' windows `[pos - bnddist - c, pos + bnddist
  + c]` overlap strictly (`c` the record's `|CIPOS[1]|`), so two positions meet when they lie less than
  `2 * bnddist + ca + cb` apart;
- the mates parsed from ALT lie on one contig, their windows likewise (`|CIEND[1]|`);
- the ALT has the same one of the four bracket forms (truvari's direction and strand).

A pair of mates written as two records is one adjacency, counted once (docs/SCHEMA.md "Breakends").
"""

import random
import threading
import time
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, SvMatching, load, resolve_ref
from just_vcf.concordance import EXCLUDED_REASONS, SITES, concordance
from just_vcf.loading import Loaded

RNG = random.Random(2820)
SEQ = {name: "".join(RNG.choice("ACGT") for _ in range(30_000)) for name in ("1", "2", "X")}

META = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
)

FORMS = ("t[p[", "t]p]", "]p]t", "[p[t")
"""VCF 4.2-4.5's four breakend forms, `t` the base(s) at POS and `p` the mate position."""

MATE_FORM = {"t[p[": "]p]t", "]p]t": "t[p[", "t]p]": "t]p]", "[p[t": "[p[t"}
"""The form the mate record writes, as in VCF 4.3 section 5.4's examples (`bnd_W`/`bnd_Y`,
`bnd_V`/`bnd_U`, `bnd_X`/`bnd_Z`)."""

BNDDIST = SvMatching().bnddist


def base(chrom: str, pos: int) -> str:
    return SEQ[chrom.removeprefix("chr")][pos - 1]


def alt(form: str, t: str, mate_chrom: str, mate_pos: int) -> str:
    return form.replace("t", t).replace("p", f"{mate_chrom}:{mate_pos}")


def bnd(chrom: str, pos: int, form: str, mate_chrom: str, mate_pos: int, gt: str = "0/1", ident: str = ".", **info: str) -> str:
    fields = {"SVTYPE": "BND", **info}
    pl_value = fields.pop("PL", None)
    fmt, sample = ("GT:PL", f"{gt}:{pl_value}") if pl_value else ("GT", gt)
    t = base(chrom, pos)
    return f"{chrom}\t{pos}\t{ident}\t{t}\t{alt(form, t, mate_chrom, mate_pos)}\t.\tPASS\t{';'.join(f'{k}={v}' for k, v in fields.items())}\t{fmt}\t{sample}"


def raw(chrom: str, pos: int, alt_text: str, gt: str = "0/1") -> str:
    return f"{chrom}\t{pos}\t.\t{base(chrom, pos)}\t{alt_text}\t.\tPASS\tSVTYPE=BND\tGT\t{gt}"


def mates(chrom: str, pos: int, form: str, mate_chrom: str, mate_pos: int, gt: str = "0/1", name: str = "bnd") -> list[str]:
    """One adjacency written as two records, each naming the other as MATEID (Manta's shape)."""
    return [
        bnd(chrom, pos, form, mate_chrom, mate_pos, gt, f"{name}_0", MATEID=f"{name}_1"),
        bnd(mate_chrom, mate_pos, MATE_FORM[form], chrom, pos, gt, f"{name}_1", MATEID=f"{name}_0"),
    ]


def meets(a: tuple, b: tuple, bnddist: int = BNDDIST) -> bool:
    """truvari 5.5 `bnd_match`, restated: `a` and `b` are `(chrom, pos, form, mate_chrom, mate_pos,
    cipos_hi, ciend_hi)`."""

    def overlaps(p1: int, c1: int, p2: int, c2: int) -> bool:
        return max(p1 - bnddist - c1, p2 - bnddist - c2) < min(p1 + bnddist + c1, p2 + bnddist + c2)

    return (
        a[0] == b[0]
        and overlaps(a[1], abs(a[5]), b[1], abs(b[5]))
        and a[3] == b[3]
        and overlaps(a[4], abs(a[6]), b[4], abs(b[6]))
        and a[2] == b[2]
    )


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa", SEQ))
    return s


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n" for name, seq in seqs.items()))
    return path


def write(path: Path, records: list[str], contigs: dict[str, int] | None = None) -> Path:
    contigs = contigs or {c: len(s) for c, s in SEQ.items()}
    meta = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={n}>" for c, n in contigs.items()), *META]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([*meta, header, *records]) + "\n")
    return path


def facet(tmp_path: Path, name: str, records: list[str], store: ReferenceStore) -> Loaded:
    return resolve_ref(load(write(tmp_path / f"{name}.vcf", records)), store)


def bnd_rows(tmp_path: Path, store: ReferenceStore, a: list[str], b: list[str], **kwargs) -> pl.DataFrame:
    result = concordance(facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store), **kwargs)
    return result.sites.collect().filter(pl.col("sv_kind") == "BND")


@pytest.mark.parametrize("form_b", FORMS)
@pytest.mark.parametrize("form_a", FORMS)
def test_each_form_meets_only_its_own_form(tmp_path: Path, store: ReferenceStore, form_a: str, form_b: str) -> None:
    rows = bnd_rows(tmp_path, store, [bnd("1", 5000, form_a, "2", 8000)], [bnd("1", 5000, form_b, "2", 8000)])
    if meets(("1", 5000, form_a, "2", 8000, 0, 0), ("1", 5000, form_b, "2", 8000, 0, 0)):
        assert rows.height == 1, rows
        row = rows.row(0, named=True)
        assert (row["key"], row["agreement"], row["variant_class"], row["absence"]) == ("breakend", "match", "sv", None)
        t = base("1", 5000)
        assert (row["alt"], row["b_alt"], row["a_pos"], row["b_pos"]) == (alt(form_a, t, "2", 8000), alt(form_b, t, "2", 8000), 5000, 5000)
        assert row["ci_needed"] is False and row["size_similarity"] is None and row["reciprocal_overlap"] is None
    else:
        # Another join at the same position: the other facet stated a breakend there, so its silence on
        # this one is `other_alleles`, never a disagreement.
        assert sorted(rows.select("agreement", "absence").rows()) == [("only_a", "other_alleles"), ("only_b", "other_alleles")]
        assert rows["key"].is_null().all()


@pytest.mark.parametrize(
    ("b_pos", "b_mate"),
    [
        (5000, 8000),
        (5000 + 2 * BNDDIST - 1, 8000),  # the position one base inside truvari's strict window overlap
        (5000 + 2 * BNDDIST, 8000),  # on its edge: outside
        (5000 - 2 * BNDDIST + 1, 8000),
        (5000 - 2 * BNDDIST, 8000),
        (5000, 8000 + 2 * BNDDIST - 1),  # the mate likewise
        (5000, 8000 + 2 * BNDDIST),
        (5000 + BNDDIST, 8000 - BNDDIST),
    ],
)
def test_both_ends_meet_within_truvaris_breakend_distance(tmp_path: Path, store: ReferenceStore, b_pos: int, b_mate: int) -> None:
    rows = bnd_rows(tmp_path, store, [bnd("1", 5000, "t[p[", "2", 8000)], [bnd("1", b_pos, "t[p[", "2", b_mate)])
    a_call, b_call = ("1", 5000, "t[p[", "2", 8000, 0, 0), ("1", b_pos, "t[p[", "2", b_mate, 0, 0)
    if meets(a_call, b_call):
        assert rows.select("key", "a_pos", "b_pos").rows() == [("breakend", 5000, b_pos)]
    else:
        assert sorted(rows["agreement"].to_list()) == ["only_a", "only_b"]
        # Near when the positions alone meet (truvari's first test), whatever the mate.
        near = meets(a_call, (*b_call[:3], *a_call[3:]))
        assert set(rows["absence"].to_list()) == {"other_alleles" if near else "unknown"}


def test_a_mate_on_another_contig_does_not_meet(tmp_path: Path, store: ReferenceStore) -> None:
    rows = bnd_rows(tmp_path, store, [bnd("1", 5000, "t[p[", "2", 8000)], [bnd("1", 5000, "t[p[", "X", 8000)])
    assert sorted(rows["agreement"].to_list()) == ["only_a", "only_b"]


def test_confidence_intervals_widen_as_truvari_reads_them(tmp_path: Path, store: ReferenceStore) -> None:
    # 205 bases apart: 5 beyond the window overlap without an interval.
    far = 5000 + 2 * BNDDIST + 5
    far_mate = 8000 + 2 * BNDDIST + 5
    cases = {
        "cipos": ({"CIPOS": "-10,10"}, (10, 0), far, 8000, True),
        # truvari's `bounds` reads index 1 of CIPOS (the upper offset) for the position, whatever the
        # lower one says: `key = 'CI' + key` runs before `idx = 0 if key == 'POS' else 1`.
        "cipos_lower_only": ({"CIPOS": "-30,5"}, (5, 0), far, 8000, False),
        "ciend_for_the_mate": ({"CIEND": "-10,10"}, (0, 10), 5000, far_mate, True),
        # The mate's distance is widened by CIEND only: CIPOS does not stand in for it (truvari reads no fallback).
        "cipos_not_for_the_mate": ({"CIPOS": "-10,10"}, (10, 0), 5000, far_mate, False),
    }
    for name, (info, (ci_pos, ci_end), b_pos, b_mate, want) in cases.items():
        rows = bnd_rows(tmp_path / name, store, [bnd("1", 5000, "t[p[", "2", 8000, **info)], [bnd("1", b_pos, "t[p[", "2", b_mate)])
        expected = meets(("1", 5000, "t[p[", "2", 8000, ci_pos, ci_end), ("1", b_pos, "t[p[", "2", b_mate, 0, 0))
        assert expected is want, name
        assert (rows["key"].to_list() == ["breakend"]) is expected, (name, rows)
        if expected:
            assert rows["ci_needed"].to_list() == [True]


def test_mates_on_two_records_count_once(tmp_path: Path, store: ReferenceStore) -> None:
    pair = mates("1", 5000, "t[p[", "2", 8000)
    # Both facets wrote both records: one adjacency, one row, shown at its first end.
    rows = bnd_rows(tmp_path / "both", store, pair, pair)
    assert rows.select("chrom", "pos", "key", "agreement").rows() == [("1", 5000, "breakend", "match")]
    # One facet wrote only the mate: still one row, the pair of records on one contig shown.
    rows = bnd_rows(tmp_path / "one", store, pair, pair[1:])
    assert rows.select("chrom", "a_pos", "b_pos", "key").rows() == [("2", 8000, 8000, "breakend")]
    # Alone, the two records are one adjacency only `a` states, never two.
    rows = bnd_rows(tmp_path / "alone", store, pair, [bnd("1", 20_000, "t]p]", "X", 100)])
    assert rows.filter(pl.col("agreement") == "only_a").height == 1
    # Each mate form the spec pairs folds; a record whose form is not its partner's mate is another adjacency.
    for form in FORMS:
        assert bnd_rows(tmp_path / f"fold{FORMS.index(form)}", store, mates("1", 5000, form, "2", 8000), []).height == 1
    wrong = [bnd("1", 5000, "t[p[", "2", 8000), bnd("2", 8000, "t[p[", "1", 5000)]
    assert bnd_rows(tmp_path / "wrong", store, wrong, []).height == 2
    # A foldback joins an end to itself: its two mates stand at one end and are one adjacency, not a
    # repeat; a third such record is one.
    fold = bnd("1", 5000, "t]p]", "1", 5000)
    for copies, repeated in ((2, False), (3, True)):
        result = concordance(facet(tmp_path / f"self{copies}", "a", [fold] * copies, store), facet(tmp_path / f"self{copies}", "b", [], store))
        assert result.sites.collect().height == 1
        assert any("repeat an end" in line for line in result.warnings) is repeated, result.warnings
    # Two foldbacks a few bases apart within their intervals are one adjacency, each record at its own end.
    near = [bnd("1", 5000, "t]p]", "1", 5000, CIPOS="-3,3"), bnd("1", 5002, "t]p]", "1", 5002, CIPOS="-3,3")]
    result = concordance(facet(tmp_path / "near", "a", near, store), facet(tmp_path / "near", "b", [], store))
    assert result.sites.collect().height == 1
    assert not any("repeat an end" in line for line in result.warnings), result.warnings


def test_mates_written_apart_within_their_intervals_count_once(tmp_path: Path, store: ReferenceStore) -> None:
    # Manta's shape over 3 bases of microhomology: each record's POS at the left of its homology
    # (CIPOS 0,3) and its ALT pointing at the right end of the other's.
    def apart(cipos: str | None, offset: int = 3) -> list[str]:
        info = {"CIPOS": cipos} if cipos else {}
        return [bnd("1", 5000, "t[p[", "2", 8000 + offset, **info), bnd("2", 8000, "]p]t", "1", 5000 + offset, **info)]

    # Each mate position lies in the other record's [POS + CIPOS[0], POS + CIPOS[1]]: one adjacency.
    assert bnd_rows(tmp_path / "within", store, apart("0,3"), []).height == 1
    # One base beyond the interval, or no interval at all: two adjacencies, as stated.
    assert bnd_rows(tmp_path / "beyond", store, apart("0,3", offset=4), []).height == 2
    assert bnd_rows(tmp_path / "none", store, apart(None), []).height == 2
    # Two records that could each be the mate: none folds, and the count says why.
    third = bnd("2", 8001, "]p]t", "1", 5002, CIPOS="0,3")
    result = concordance(facet(tmp_path / "several", "a", [*apart("0,3"), third], store), facet(tmp_path / "several", "b", [], store))
    assert result.sites.collect().height == 3
    assert any("could each be the mate of several" in line for line in result.warnings), result.warnings
    # Against a facet that wrote the mates exactly, the folded adjacency meets once.
    rows = bnd_rows(tmp_path / "against", store, apart("0,3"), mates("1", 5000, "t[p[", "2", 8000))
    assert rows.select("key", "agreement").rows() == [("breakend", "match")]


def test_a_mate_record_that_states_another_genotype_is_counted(tmp_path: Path, store: ReferenceStore) -> None:
    pair = mates("1", 5000, "t[p[", "2", 8000)
    split = [pair[0], pair[1].replace("\t0/1", "\t1/1")]
    result = concordance(facet(tmp_path, "a", split, store), facet(tmp_path, "b", pair, store))
    assert result.sites.collect().filter(pl.col("sv_kind") == "BND").height == 1
    assert any("mate records state different genotypes" in line for line in result.warnings), result.warnings


def test_single_and_unreadable_breakends_are_counted_apart(tmp_path: Path, store: ReferenceStore) -> None:
    t = base("1", 5000)
    records = [
        raw("1", 5000, f".{t}"),  # single breakends: truvari's `no_single_bnd` default leaves them out
        raw("1", 6000, f"{base('1', 6000)}."),
        raw("1", 7000, f"{base('1', 7000)}[2:abc["),  # no position
        raw("1", 8000, f"{base('1', 8000)}[2:100]"),  # brackets that do not pair
        raw("1", 9000, f"[2:100[{base('1', 9000)}["),
        bnd("1", 12_000, "t[p[", "2", 8000),
    ]
    result = concordance(facet(tmp_path, "a", records, store), facet(tmp_path, "b", records, store))
    excluded = {(row["facet"], row["reason"]): row["alleles"] for row in result.excluded.iter_rows(named=True)}
    assert excluded == {(s, r): n for s in ("a", "b") for r, n in (("single_breakend", 2), ("breakend_unparsed", 3))}
    assert {"single_breakend", "breakend_unparsed"} <= EXCLUDED_REASONS and "breakend" not in EXCLUDED_REASONS
    sites = result.sites.collect()
    assert sites.select("pos", "key", "agreement").rows() == [(12_000, "breakend", "match")]


def test_genotypes_of_met_breakends_compare_as_copies(tmp_path: Path, store: ReferenceStore) -> None:
    cases = [("0/1", "0/1", "match"), ("0/1", "1/1", "mismatch"), ("1/1", "./.", "uncalled")]
    a = [bnd("1", 2000 + 5000 * i, "]p]t", "X", 1000 + 5000 * i, gt=ga) for i, (ga, _, _) in enumerate(cases)]
    b = [bnd("1", 2010 + 5000 * i, "]p]t", "X", 1000 + 5000 * i, gt=gb) for i, (_, gb, _) in enumerate(cases)]
    rows = bnd_rows(tmp_path, store, a, b).sort("a_pos")
    assert rows["agreement"].to_list() == [want for *_, want in cases]
    assert rows["a_copies"].to_list() == [1, 1, 2] and rows["b_copies"].to_list() == [1, 2, None]


def test_met_breakends_with_likelihoods_are_graded(tmp_path: Path, store: ReferenceStore) -> None:
    pa, pb = [0, 30, 300], [20, 0, 40]
    a = [bnd("1", 5000, "t]p]", "2", 8000, PL=",".join(map(str, pa)))]
    b = [bnd("1", 5000, "t]p]", "2", 8000, PL=",".join(map(str, pb)))]

    def probs(values: list[int]) -> list[float]:
        weights = [10 ** (-v / 10) for v in values]
        return [w / sum(weights) for w in weights]

    rows = bnd_rows(tmp_path, store, a, b)
    assert rows["basis"].to_list() == ["gl_flat_prior"]
    assert rows["p_same"][0] == pytest.approx(sum(x * y for x, y in zip(probs(pa), probs(pb))))


def test_each_adjacency_is_taken_once_by_truvaris_score(tmp_path: Path, store: ReferenceStore) -> None:
    # Two of b's breakends meet a's; truvari's score, `1 - mean(|start|, |end| distance) / (2 * bnddist)`,
    # picks the nearer in sum, then the nearer start.
    a = [bnd("1", 5000, "t[p[", "2", 8000)]
    b = [bnd("1", 5050, "t[p[", "2", 8000), bnd("1", 5010, "t[p[", "2", 8030)]

    def score(pos: int, mate: int) -> float:
        return max(0.0, 1 - ((abs(5000 - pos) + abs(8000 - mate)) / 2) / (2 * BNDDIST))

    best = max(((5050, 8000), (5010, 8030)), key=lambda c: (score(*c), -abs(5000 - c[0])))
    rows = bnd_rows(tmp_path / "x", store, a, b)
    assert rows.filter(pl.col("key") == "breakend")["b_pos"].to_list() == [best[0]]
    assert rows.filter(pl.col("key").is_null()).select("agreement", "absence").rows() == [("only_b", "other_alleles")]
    again = bnd_rows(tmp_path / "y", store, a, list(reversed(b)))
    assert again.equals(rows)
    # Swapping the facets swaps the sides and keeps the pair.
    swapped = bnd_rows(tmp_path / "z", store, b, a).filter(pl.col("key") == "breakend")
    assert swapped.select("a_pos", "b_pos").rows() == [(best[0], 5000)]


def test_bnddist_is_a_parameter(tmp_path: Path, store: ReferenceStore) -> None:
    a, b = [bnd("1", 5000, "t[p[", "2", 8000)], [bnd("1", 5300, "t[p[", "2", 8000)]
    for bnddist in (100, 151, 150):
        rows = bnd_rows(tmp_path / str(bnddist), store, a, b, sv=SvMatching(bnddist=bnddist))
        expected = meets(("1", 5000, "t[p[", "2", 8000, 0, 0), ("1", 5300, "t[p[", "2", 8000, 0, 0), bnddist)
        assert (rows["key"].to_list() == ["breakend"]) is expected, bnddist


def test_breakend_rows_join_the_sites_frame_in_order(tmp_path: Path, store: ReferenceStore) -> None:
    a = [*mates("2", 300, "[p[t", "1", 9000), *mates("1", 400, "t]p]", "X", 50, name="w")]
    b = [*reversed(a)]
    result = concordance(facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store))
    sites = result.sites.collect()
    assert sites.schema == SITES
    # Each adjacency once, at the end the canonical order puts first, whatever order the file wrote.
    assert sites.select("chrom", "pos", "key").rows() == [("1", 400, "breakend"), ("1", 9000, "breakend")]
    by_kind = dict(result.summary.group_by("sv_kind").agg(pl.col("alleles").sum()).iter_rows())
    assert by_kind == {"BND": 2}
    assert concordance(facet(tmp_path, "a2", a, store), facet(tmp_path, "b2", b, store)).sites.collect().equals(sites)


def with_facet(loaded: Loaded, **values: object) -> Loaded:
    """`loaded` with facet columns overwritten, typed as the facets frame types them."""
    schema = loaded.facets.collect_schema()
    return Loaded(
        facets=loaded.facets.with_columns(*[pl.lit(v).cast(schema[k]).alias(k) for k, v in values.items()]), calls=loaded.calls, blocks=loaded.blocks
    )


def test_a_lone_adjacency_reads_the_other_facets_absence_rule_at_both_ends(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [bnd("1", 5000, "t[p[", "2", 8000), bnd("1", 15_000, "t[p[", "X", 300)], store)
    # A whole-genome callset that called SVs (one breakend far away) and declared contigs 1 and 2, not X.
    b_path = write(tmp_path / "b.vcf", [bnd("1", 25_000, "t]p]", "2", 20_000)], contigs={c: len(SEQ[c]) for c in ("1", "2")})
    b = with_facet(resolve_ref(load(b_path), store), absence_means="weak_reference")
    rows = concordance(a, b).sites.collect().filter(pl.col("agreement") == "only_a")
    # Absence reads as reference only where b could have spoken for both ends.
    assert rows.select("pos", "absence").sort("pos").rows() == [(5000, "reference"), (15_000, "unknown")]


def test_mate_contigs_are_read_by_their_canonical_name_and_must_be_shared(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [bnd("1", 5000, "t[p[", "X", 300), bnd("1", 9000, "[p[t", "2", 500)], store)
    # Another collection: UCSC names over the same contigs 1 and X, contig 2 left out.
    ucsc = {"chr1": SEQ["1"], "chrX": SEQ["X"]}
    digest = store.add_fasta(write_fasta(tmp_path / "ucsc.fa", ucsc))
    b_path = write(tmp_path / "b.vcf", [bnd("chr1", 5000, "t[p[", "chrX", 300)], contigs={k: len(v) for k, v in ucsc.items()})
    b = resolve_ref(load(b_path), store, collection=digest)
    result = concordance(a, b, reference=store)
    assert result.sites.collect().select("chrom", "pos", "key").rows() == [("1", 5000, "breakend")]
    # a's join into contig 2, which b's collection does not hold, is not compared.
    assert {(r["facet"], r["reason"]): r["alleles"] for r in result.excluded.iter_rows(named=True)} == {("a", "contig_not_shared"): 1}


def test_an_assembly_mate_is_compared_by_its_name_as_written(tmp_path: Path, store: ReferenceStore) -> None:
    """VCF 4.3 §5.4.2's large insertion names its mate on a contig of the `##assembly` file, `<ctg1>`,
    not of the reference (§1.6.1). truvari 5.5's `bnd_position` reads it as the ALT spells it (its
    `[^\\[\\]:]+` takes `<ctg1>`, brackets and all), and so it is here: on facets of two collections, as
    on one, it is no contig either collection must share, and two such mates meet only by that name."""
    t5, t9 = base("1", 5000), base("1", 9000)
    a = facet(tmp_path, "a", [raw("1", 5000, f"{t5}[<ctg1>:7["), raw("1", 9000, f"]<ctg2>:214]{t9}")], store)
    ucsc = {"chr1": SEQ["1"], "chrX": SEQ["X"]}
    digest = store.add_fasta(write_fasta(tmp_path / "ucsc.fa", ucsc))
    b_records = [raw("chr1", 5000, f"{t5}[<ctg1>:7["), raw("chr1", 9000, f"]<ctg3>:214]{t9}")]
    b = resolve_ref(load(write(tmp_path / "b.vcf", b_records, contigs={k: len(v) for k, v in ucsc.items()})), store, collection=digest)
    stated_a = [("1", 5000, "t[p[", "<ctg1>", 7, 0, 0), ("1", 9000, "]p]t", "<ctg2>", 214, 0, 0)]
    stated_b = [("1", 5000, "t[p[", "<ctg1>", 7, 0, 0), ("1", 9000, "]p]t", "<ctg3>", 214, 0, 0)]
    expected = [(x[1], "match") for x in stated_a if any(meets(x, y) for y in stated_b)]
    expected += [(x[1], "only_a") for x in stated_a if not any(meets(x, y) for y in stated_b)]
    expected += [(y[1], "only_b") for y in stated_b if not any(meets(x, y) for x in stated_a)]
    result = concordance(a, b, reference=store)
    assert result.excluded.is_empty()
    rows = result.sites.collect()
    assert sorted(rows.select("pos", "agreement").rows()) == sorted(expected) == [(5000, "match"), (9000, "only_a"), (9000, "only_b")]
    assert rows.filter(pl.col("agreement") == "match")["key"].to_list() == ["breakend"]


def test_a_lone_adjacency_into_the_assembly_file_reads_unknown(tmp_path: Path, store: ReferenceStore) -> None:
    """A whole-genome callset that called SVs reads as reference where it declared both ends' contigs.
    No caller declares a contig of another file's `##assembly`: it could state a join into one only
    with an assembly file of its own, so its silence on such an adjacency is unknown, not reference."""
    a = facet(tmp_path, "a", [bnd("1", 5000, "t[p[", "2", 8000), raw("1", 15_000, f"{base('1', 15_000)}[<ctg1>:7[")], store)
    b = with_facet(resolve_ref(load(write(tmp_path / "b.vcf", [bnd("1", 25_000, "t]p]", "2", 20_000)])), store), absence_means="weak_reference")
    rows = concordance(a, b).sites.collect().filter(pl.col("agreement") == "only_a")
    assert rows.select("pos", "absence").sort("pos").rows() == [(5000, "reference"), (15_000, "unknown")]


def _rss() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    raise AssertionError("/proc/self/status carries no VmRSS line, so the peak cannot be measured")


@pytest.mark.skipif(not Path("/proc/self/status").exists(), reason="the peak is read from /proc/self/status (Linux)")
def test_breakends_crowded_on_one_contig_compare_within_a_bounded_peak(tmp_path: Path, store: ReferenceStore) -> None:
    """1,500 breakend records on one contig compared with themselves stay within 1 GB above the start.
    The near pairs are a range join per contig; a contig equality beside the range made polars hold
    every pair of records on the contig at once. Measured 2026-10-26 with 2,000 records on one contig:
    5.7 GB that way, 29 MB as a range join; this test's 1,500 took 3 GB that way."""
    rng = random.Random(1500)
    positions = sorted(rng.sample(range(100, len(SEQ["1"]) - 100), 1500))
    records = [raw("1", pos, f"{base('1', pos)}[1:{rng.randrange(100, len(SEQ['1']) - 100)}[") for pos in positions]
    one = facet(tmp_path, "a", records, store)
    peak = [0]
    done = threading.Event()

    def sample() -> None:
        while not done.is_set():
            peak[0] = max(peak[0], _rss())
            time.sleep(0.02)

    before = _rss()
    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    sites = concordance(one, one).sites.collect()
    done.set()
    sampler.join()
    # Each record meets itself and is taken once.
    assert sites.select("pos", "agreement").unique()["pos"].sort().to_list() == positions
    assert peak[0] - before < 2**30, f"peak {(peak[0] - before) / 2**20:.0f} MB above the start"


def per_alt_facet(tmp_path: Path, name: str, records: list[str], store: ReferenceStore) -> Loaded:
    """A facet whose header declares CIPOS and CIEND with any count (VCF 4.5's two per ALT)."""
    meta = [line.replace("ID=CIPOS,Number=2", "ID=CIPOS,Number=.").replace("ID=CIEND,Number=2", "ID=CIEND,Number=.") for line in META]
    contigs = [f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items()]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    path = tmp_path / f"{name}.vcf"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(["##fileformat=VCFv4.3", *contigs, *meta, header, *records]) + "\n")
    return resolve_ref(load(path), store)


def test_each_alt_of_a_record_reads_its_own_interval_pair(tmp_path: Path, store: ReferenceStore) -> None:
    """Two bracket ALTs in one record, each with its own CIPOS and CIEND pair (VCF 4.5): ALT `k`'s
    windows read pair `k`'s upper offsets, and so do the position intervals mates fold within."""
    t = base("1", 5000)
    record = f"1\t5000\t.\t{t}\t{t}[2:8000[,{t}]X:300]\t.\tPASS\tSVTYPE=BND;CIPOS=0,0,0,300;CIEND=0,400,0,0\tGT\t1/2"
    a = per_alt_facet(tmp_path, "a", [record], store)
    assert a.calls.collect().select("cipos", "ciend").rows() == [([0, 0, 0, 300], [0, 400, 0, 0])]
    b_records = [bnd("1", 5000, "t[p[", "2", 8450), bnd("1", 5350, "t]p]", "X", 300)]
    alts = [("1", 5000, "t[p[", "2", 8000, 0, 400), ("1", 5000, "t]p]", "X", 300, 300, 0)]
    stated = [("1", 5000, "t[p[", "2", 8450, 0, 0), ("1", 5350, "t]p]", "X", 300, 0, 0)]
    # Each ALT meets its own b record with its own pair, and with the other ALT's pair would not.
    assert [meets(alt, other) for alt, other in zip(alts, stated, strict=True)] == [True, True]
    swapped = [(*alts[0][:5], *alts[1][5:]), (*alts[1][:5], *alts[0][5:])]
    assert [meets(alt, other) for alt, other in zip(swapped, stated, strict=True)] == [False, False]
    result = concordance(a, per_alt_facet(tmp_path, "b", b_records, store))
    rows = result.sites.collect().filter(pl.col("sv_kind") == "BND")
    assert sorted(rows.filter(pl.col("key") == "breakend").select("a_pos", "b_pos").rows()) == [(5000, 5000), (5000, 5350)]
    assert rows.height == 2

    # Folding: the second ALT and a mate record written 3 bases apart, within ALT 1's CIPOS pair only.
    mate_record = bnd("2", 8000, "]p]t", "1", 5003, CIPOS="0,3")
    two = f"1\t5000\t.\t{t}\t{t}]X:300],{t}[2:8003[\t.\tPASS\tSVTYPE=BND;CIPOS=0,0,0,3\tGT\t1/2"
    folded = concordance(per_alt_facet(tmp_path / "fold", "a", [two, mate_record], store), per_alt_facet(tmp_path / "fold", "b", [], store))
    assert folded.sites.collect().height == 2  # the adjacency to X, and the folded pair of mates
    wide_first = two.replace("CIPOS=0,0,0,3", "CIPOS=0,3,0,0")
    unfolded = concordance(per_alt_facet(tmp_path / "nofold", "a", [wide_first, mate_record], store), per_alt_facet(tmp_path / "nofold", "b", [], store))
    assert unfolded.sites.collect().height == 3
