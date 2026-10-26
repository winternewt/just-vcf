"""Concordance of structural-variant alleles between two facets (RM282): matching by type and extent.

Each case writes two one-sample SV VCFs over a small refget store, reads REF and compares. The
thresholds are truvari bench's defaults (`SvMatching`); every expected similarity, overlap and
distance is recomputed here from the coordinates the fixture writes, never copied from a run. Each
call's extent in truvari's terms is the half-open `[POS - 1, END)`: the anchor base before the event
and the event itself (docs/SCHEMA.md "Concordance", "Structural variants").
"""

import itertools
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, SvMatching, load, resolve_ref
from just_vcf.concordance import SITES, SUMMARY, concordance, summarize
from just_vcf.loading import Loaded

RNG = random.Random(282)
SEQ = {name: "".join(RNG.choice("ACGT") for _ in range(30_000)) for name in ("1", "2")}

SV_META = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
)


def base(chrom: str, pos: int) -> str:
    return SEQ[chrom][pos - 1]


def sv(pos: int, alt: str, end: int | None = None, gt: str = "0/1", chrom: str = "1", **info: str) -> str:
    """One SV record; `end` is INFO END (an `<INS>` gets END = POS, as Manta writes it)."""
    fields = {"SVTYPE": alt.strip("<>").split(":")[0] if alt.startswith("<") else "BND", "END": str(end if end is not None else pos)}
    fields.update(info)
    pl_value = fields.pop("PL", None)
    fmt, sample = ("GT:PL", f"{gt}:{pl_value}") if pl_value else ("GT", gt)
    return f"{chrom}\t{pos}\t.\t{base(chrom, pos)}\t{alt}\t.\tPASS\t{';'.join(f'{k}={v}' for k, v in fields.items())}\t{fmt}\t{sample}"


def snv(pos: int, gt: str = "0/1", chrom: str = "1") -> str:
    ref = base(chrom, pos)
    return f"{chrom}\t{pos}\t.\t{ref}\t{next(x for x in 'ACGT' if x != ref)}\t.\tPASS\t.\tGT\t{gt}"


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n" for name, seq in SEQ.items()))
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(fasta)
    return s


def facet(tmp_path: Path, name: str, records: list[str], store: ReferenceStore) -> Loaded:
    meta = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items()), *SV_META]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    path = tmp_path / f"{name}.vcf"  # not `.sv.vcf`: that name marks a file cut to SVs (content_filter)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([*meta, header, *records]) + "\n")
    return resolve_ref(load(path), store)


def with_facet(loaded: Loaded, **values: object) -> Loaded:
    """`loaded` with facet columns overwritten, typed as the facets frame types them."""
    schema = loaded.facets.collect_schema()
    return Loaded(
        facets=loaded.facets.with_columns(*[pl.lit(v).cast(schema[k]).alias(k) for k, v in values.items()]), calls=loaded.calls, blocks=loaded.blocks
    )


def sv_rows(tmp_path: Path, store: ReferenceStore, a: list[str], b: list[str], **kwargs) -> pl.DataFrame:
    result = concordance(facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store), **kwargs)
    return result.sites.collect().filter(pl.col("variant_class") == "sv")


def span(pos: int, end: int) -> tuple[int, int]:
    """truvari's extent of a call: half-open, from the anchor base."""
    return pos - 1, end


def size_similarity(a: int, b: int) -> float:
    if a == 0 or b == 0:
        return 1.0 if a == b else min(max(a, 1), max(b, 1)) / max(max(a, 1), max(b, 1))
    return min(a, b) / max(a, b)


def reciprocal_overlap(a: tuple[int, int], b: tuple[int, int]) -> float:
    start, end = max(a[0], b[0]), min(a[1], b[1])
    return (end - start) / max(a[1] - a[0], b[1] - b[0]) if start < end else 0.0


def within_refdist(a: tuple[int, int], b: tuple[int, int], refdist: int = 500) -> bool:
    return max(a[0] - refdist, b[0]) < min(a[1] + refdist, b[1])


A_DEL = (1000, 1300)


@pytest.mark.parametrize(
    ("b_pos", "b_end"),
    [
        (1800, 2100),  # b's extent starts 500 bases past a's end: inside refdist (a strict overlap of a's span widened by 500)
        (1801, 2101),  # one base further: outside
        (1000, 1210),  # size 210 against 300: similarity 0.7, the threshold
        (1000, 1209),  # 209 against 300: below it
        (1000, 1300),  # the same call
        (400, 700),  # before a, inside refdist
    ],
)
def test_deletions_match_by_reference_distance_and_size(tmp_path: Path, store: ReferenceStore, b_pos: int, b_end: int) -> None:
    rows = sv_rows(tmp_path, store, [sv(*A_DEL[:1], "<DEL>", A_DEL[1])], [sv(b_pos, "<DEL>", b_end)])
    a_span, b_span = span(*A_DEL), span(b_pos, b_end)
    size = size_similarity(A_DEL[1] - A_DEL[0], b_end - b_pos)
    expected = within_refdist(a_span, b_span) and size >= 0.7
    if expected:
        assert rows.height == 1, rows
        row = rows.row(0, named=True)
        assert (row["key"], row["agreement"], row["sv_kind"], row["absence"]) == ("extent", "match", "DEL", None)
        assert (row["a_pos"], row["a_end_pos"], row["b_pos"], row["b_end_pos"]) == (A_DEL[0], A_DEL[1], b_pos, b_end)
        assert row["size_similarity"] == pytest.approx(size)
        assert row["reciprocal_overlap"] == pytest.approx(reciprocal_overlap(a_span, b_span))
        assert row["ci_needed"] is False
    else:
        assert sorted(rows["agreement"].to_list()) == ["only_a", "only_b"]
        # b states a deletion within refdist of a's (only too small): that is what b's silence on a's call is.
        near = within_refdist(a_span, b_span)
        assert set(rows["absence"].to_list()) == {"other_alleles" if near else "unknown"}
        assert rows["key"].is_null().all()


def test_reciprocal_overlap_threshold(tmp_path: Path, store: ReferenceStore) -> None:
    for b_pos in (1150, 1152):
        a_span, b_span = span(*A_DEL), span(b_pos, b_pos + 300)
        overlap = reciprocal_overlap(a_span, b_span)
        rows = sv_rows(tmp_path / str(b_pos), store, [sv(1000, "<DEL>", 1300)], [sv(b_pos, "<DEL>", b_pos + 300)], sv=SvMatching(pctovl=0.5))
        assert (rows["key"].to_list() == ["extent"]) == (overlap >= 0.5), (b_pos, overlap, rows)
        # truvari's default, 0.0, takes either.
        assert sv_rows(tmp_path / f"{b_pos}d", store, [sv(1000, "<DEL>", 1300)], [sv(b_pos, "<DEL>", b_pos + 300)]).height == 1


def test_types_match_by_the_alt_and_subtypes_are_ignored(tmp_path: Path, store: ReferenceStore) -> None:
    # A duplication over a deletion's extent is another event; `<DUP:TANDEM>` is a `<DUP>`.
    a = [sv(1000, "<DEL>", 1300), sv(5000, "<DUP:TANDEM>", 5400), sv(9000, "<INV>", 9500)]
    b = [sv(1000, "<DUP>", 1300), sv(5010, "<DUP>", 5400), sv(9400, "<INV>", 9900)]
    rows = sv_rows(tmp_path / "typed", store, a, b)
    by = {(row["a_pos"], row["b_pos"]): row for row in rows.iter_rows(named=True)}
    assert by[(5000, 5010)]["agreement"] == "match" and by[(5000, 5010)]["b_alt"] == "<DUP>" and by[(5000, 5010)]["sv_kind"] == "DUP"
    # The inversions: shifted 400 bases, the same size, within refdist.
    assert by[(9000, 9400)]["size_similarity"] == pytest.approx(1.0)
    assert by[(9000, 9400)]["reciprocal_overlap"] == pytest.approx(reciprocal_overlap(span(9000, 9500), span(9400, 9900)))
    lone = {(row["agreement"], row["alt"], row["absence"]) for row in rows.filter(pl.col("key").is_null()).iter_rows(named=True)}
    assert lone == {("only_a", "<DEL>", "other_alleles"), ("only_b", "<DUP>", "other_alleles")}
    # truvari's --typeignore pairs the deletion with the duplication.
    ignored = sv_rows(tmp_path / "ignored", store, a, b, sv=SvMatching(typeignore=True))
    assert ignored.filter(pl.col("a_pos") == 1000)["b_alt"].to_list() == ["<DUP>"]
    # SVTYPE is the record's claim, not the allele's type: DRAGEN's CNV caller writes SVTYPE=CNV on a
    # `<DEL>`, which meets another caller's `<DEL>`.
    dragen = sv_rows(tmp_path / "dragen", store, [sv(1000, "<DEL>", 1300, SVTYPE="CNV")], [sv(1000, "<DEL>", 1300)])
    assert dragen["key"].to_list() == ["extent"]


def test_insertions_meet_by_distance_and_their_size_is_unknown(tmp_path: Path, store: ReferenceStore) -> None:
    # The frame keeps no SVLEN, so an `<INS>` has no size: its similarity is withheld, never 1.
    for b_pos, meets in ((5300, True), (5500, True), (5501, False)):
        assert within_refdist(span(5000, 5000), span(b_pos, b_pos)) is meets
        result = concordance(facet(tmp_path / str(b_pos), "a", [sv(5000, "<INS>")], store), facet(tmp_path / str(b_pos), "b", [sv(b_pos, "<INS>")], store))
        rows = result.sites.collect()
        if meets:
            assert rows["key"].to_list() == ["extent"] and rows["size_similarity"].is_null().all()
            assert any("without a size" in line for line in result.warnings), result.warnings
        else:
            assert rows["agreement"].sort().to_list() == ["only_a", "only_b"]


def test_confidence_intervals_widen_the_reference_distance(tmp_path: Path, store: ReferenceStore) -> None:
    far = sv(1840, "<DEL>", 2140)  # its extent starts at 1839, 39 bases beyond a's span widened by 500
    assert not within_refdist(span(1000, 1300), span(1840, 2140))
    cases = {
        "none": (sv(1000, "<DEL>", 1300), False),
        "ciend": (sv(1000, "<DEL>", 1300, CIPOS="-5,5", CIEND="-50,50"), True),  # the end may lie 50 bases later
        "cipos_for_ciend": (sv(1000, "<DEL>", 1300, CIPOS="-50,50"), True),  # no CIEND: VCF 4.5 says assume CIPOS
        "too_narrow": (sv(1000, "<DEL>", 1300, CIEND="-30,30"), False),
    }
    for name, (record, meets) in cases.items():
        rows = sv_rows(tmp_path / name, store, [record], [far])
        assert (rows["key"].to_list() == ["extent"]) is meets, (name, rows)
        if meets:
            assert rows["ci_needed"].to_list() == [True]
    # b's own interval widens its start toward a.
    rows = sv_rows(tmp_path / "b_cipos", store, [sv(1000, "<DEL>", 1300)], [sv(1840, "<DEL>", 2140, CIPOS="-40,40")])
    assert rows["key"].to_list() == ["extent"] and rows["ci_needed"].to_list() == [True]


def test_genotypes_of_matched_calls_compare_as_copies(tmp_path: Path, store: ReferenceStore) -> None:
    cases = [("0/1", "0/1", "match"), ("0/1", "1/1", "mismatch"), ("1/1", "./.", "uncalled"), ("1", "1/1", "match")]
    a = [sv(1000 + 3000 * i, "<DEL>", 1300 + 3000 * i, gt=ga) for i, (ga, _, _) in enumerate(cases)]
    b = [sv(1050 + 3000 * i, "<DEL>", 1300 + 3000 * i, gt=gb) for i, (_, gb, _) in enumerate(cases)]
    rows = sv_rows(tmp_path, store, a, b).sort("a_pos")
    assert rows["agreement"].to_list() == [want for *_, want in cases]
    assert rows["a_copies"].to_list() == [1, 1, 2, 1] and rows["b_copies"].to_list() == [1, 2, None, 2]
    assert rows["ploidy_differs"].to_list() == [False, False, False, True]


def test_matched_calls_with_likelihoods_are_graded(tmp_path: Path, store: ReferenceStore) -> None:
    pa, pb = [0, 30, 300], [20, 0, 40]
    rows = sv_rows(tmp_path, store, [sv(1000, "<DEL>", 1300, PL=",".join(map(str, pa)))], [sv(1020, "<DEL>", 1300, PL=",".join(map(str, pb)))])

    def probs(values: list[int]) -> list[float]:
        weights = [10 ** (-v / 10) for v in values]
        return [w / sum(weights) for w in weights]

    expected = sum(x * y for x, y in zip(probs(pa), probs(pb)))
    assert rows["basis"].to_list() == ["gl_flat_prior"]
    assert rows["p_same"][0] == pytest.approx(expected)


def test_each_call_is_matched_once_to_its_closest(tmp_path: Path, store: ReferenceStore) -> None:
    # Two of b's deletions both meet a's; the one with the better score (size similarity plus
    # overlap), then the nearer start, is taken; the other is b's alone, near a call of a's.
    b = [sv(1100, "<DEL>", 1400), sv(1000, "<DEL>", 1290)]
    rows = sv_rows(tmp_path / "x", store, [sv(1000, "<DEL>", 1300)], b)
    scores = {
        p: size_similarity(300, e - p) + reciprocal_overlap(span(1000, 1300), span(p, e))
        for p, e in ((1100, 1400), (1000, 1290))
    }
    best = max(scores, key=scores.get)
    assert rows.filter(pl.col("key") == "extent")["b_pos"].to_list() == [best]
    assert rows.filter(pl.col("key").is_null()).select("agreement", "absence").rows() == [("only_b", "other_alleles")]
    # The same inputs give the same rows, whichever order b wrote them in.
    again = sv_rows(tmp_path / "y", store, [sv(1000, "<DEL>", 1300)], list(reversed(b)))
    assert again.equals(rows)


def test_a_facet_that_cannot_see_svs_leaves_them_unknown(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [sv(1000, "<DEL>", 1300), snv(20_000)], store)
    small = facet(tmp_path, "small", [snv(20_000), snv(25_000)], store)
    # A whole-genome small-variant callset reads absence as reference for the classes it called, never
    # for SVs it did not; the same callset cut to SNVs says nothing; one that called SVs says reference.
    wgs_small = with_facet(small, absence_means="weak_reference")
    sv_caller = with_facet(facet(tmp_path, "caller", [sv(9000, "<DEL>", 9300), snv(20_000), snv(25_000)], store), absence_means="weak_reference")
    cut = with_facet(sv_caller, content_filter="snv_only")
    for other, reading in ((small, "unknown"), (wgs_small, "unknown"), (cut, "unknown"), (sv_caller, "reference")):
        rows = concordance(a, other).sites.collect().filter(pl.col("a_pos") == 1000)
        assert rows.select("agreement", "absence").rows() == [("only_a", reading)]
    # Never a disagreement: a one-sided SV is only ever `only_a` / `only_b`.
    assert "mismatch" not in concordance(a, wgs_small).sites.collect()["agreement"].to_list()


def test_copy_number_alleles_are_read_per_alt(tmp_path: Path, store: ReferenceStore) -> None:
    # Canvas: `<CN0>,<CN2>` with GT 1/2 holds one copy of each; a `<CN0>` in b meets the first.
    a = [sv(8000, "<CN0>,<CN2>", 8999, gt="1/2", SVTYPE="LOH")]
    b = [sv(8000, "<CN0>", 8999, gt="0/1", SVTYPE="CNV")]
    rows = sv_rows(tmp_path, store, a, b)
    met = rows.filter(pl.col("key") == "extent").row(0, named=True)
    assert (met["alt"], met["sv_kind"], met["a_copies"], met["b_copies"], met["agreement"]) == ("<CN0>", "CN0", 1, 1, "match")
    assert rows.filter(pl.col("alt") == "<CN2>").select("agreement", "a_copies").rows() == [("only_a", 1)]


def test_sv_rows_join_the_sites_frame_its_summary_and_excluded_counts(tmp_path: Path, store: ReferenceStore) -> None:
    a = [sv(1000, "<DEL>", 1300), sv(5000, "<INS>"), "1\t7000\t.\t" + base("1", 7000) + "\t" + base("1", 7000) + "]2:7000]\t.\tPASS\tSVTYPE=BND\tGT\t0/1", snv(20_000)]
    b = [sv(1010, "<DEL>", 1300), snv(20_000)]
    result = concordance(facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store))
    sites = result.sites.collect()
    assert sites.schema == SITES and result.summary.schema == SUMMARY
    assert result.summary.equals(summarize(sites)) and int(result.summary["alleles"].sum()) == sites.height
    by_kind = dict(result.summary.group_by("sv_kind").agg(pl.col("alleles").sum()).iter_rows())
    assert by_kind == {"DEL": 1, "INS": 1, None: 1}
    excluded = {(row["facet"], row["reason"]): row["alleles"] for row in result.excluded.iter_rows(named=True)}
    assert excluded == {("a", "breakend"): 1}
    # Deterministic: the same comparison twice gives identical frames, SV rows in canonical order.
    assert concordance(facet(tmp_path, "a2", a, store), facet(tmp_path, "b2", b, store)).sites.collect().equals(sites)
    assert sites["pos"].to_list() == sorted(sites["pos"].to_list())


def test_symmetric_comparison(tmp_path: Path, store: ReferenceStore) -> None:
    # Swapping the facets swaps the sides and keeps the pairs.
    a = [sv(1000 + 2000 * i, kind, 1000 + 2000 * i + size) for i, (kind, size) in enumerate(itertools.product(("<DEL>", "<DUP>"), (200, 400)))]
    b = [sv(1100 + 2000 * i, kind, 1100 + 2000 * i + 300) for i, kind in enumerate(("<DEL>", "<DEL>", "<DUP>", "<DUP>"))]
    fa, fb = facet(tmp_path, "a", a, store), facet(tmp_path, "b", b, store)
    forward = concordance(fa, fb).sites.collect().filter(pl.col("key") == "extent").select("a_pos", "b_pos")
    backward = concordance(fb, fa).sites.collect().filter(pl.col("key") == "extent").select(pl.col("b_pos").alias("a_pos"), pl.col("a_pos").alias("b_pos"))
    assert forward.sort("a_pos").equals(backward.sort("a_pos"))
    assert forward.height == sum(size_similarity(s, 300) >= 0.7 for s in (200, 400, 200, 400))
