"""RM3: the anchor base of an insertion point (CGI, Genos, CASAVA) and REF over a Genos row's stated span,
read by `resolve_ref` from a small synthetic reference; expectations computed from the sequence."""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, genotype_evidence, load, mint_vrs_ids, normalize, resolve_ref, write_vcf
from just_vcf.loading import Loaded
from test_load_cgi import Locus, Seg, master_var_text, var_text
from test_load_cgi import write as write_cgi
from test_reference import SEQ, base, frames, other, write_fasta


@pytest.fixture
def grch37_store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "grch37")
    s.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    return s


# --- CGI ---------------------------------------------------------------------------------------------


def cgi_point(chrom: str, begin: int, *seqs: str) -> Locus:
    """A zero-length CGI locus (an insertion point at 0-based `begin`), one haplotype per sequence; `""`
    is the reference haplotype there."""
    haps = tuple((Seg(begin, begin, "ref" if seq == "" else "ins", "", seq),) for seq in seqs)
    return Locus(chrom, begin, begin, haps=haps)


P_HET, P_HOM, P_NONE, P_DEL = 100, 200, 300, 400


def point_loci() -> list[Locus]:
    return [
        cgi_point("chr1", P_HET, "", "GG"),
        cgi_point("chr1", P_HOM, "TCA", "TCA"),
        Locus("chr1", P_NONE, P_NONE, uniform="no-call"),
        Locus(
            "chr1", P_DEL - 1, P_DEL,
            haps=((Seg(P_DEL - 1, P_DEL, "ref", base("1", P_DEL), base("1", P_DEL)),), (Seg(P_DEL - 1, P_DEL, "del", base("1", P_DEL), ""),)),
        ),
        cgi_point("chr2", 0, "", "GG"),  # before the contig's first base
        cgi_point("chr3", 50, "", "A"),  # a contig the reference lacks
    ]  # fmt: skip


def resolved_var(tmp_path: Path, store: ReferenceStore) -> Loaded:
    return resolve_ref(load(write_cgi(tmp_path, var_text(point_loci()), "var-GS000000001-ASM.tsv.bz2")), store)


def test_a_cgi_insertion_point_gets_its_anchor_base(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """CGI writes an insertion at a zero-length locus, `begin == end`, between 0-based `begin - 1` and
    `begin`, which the loader places at 1-based `pos = begin`: the base before it, VCF's anchor. REF is
    that base read from the reference, prefixed onto every allele; there is nothing to compare it with,
    since the source wrote an empty REF."""
    resolved = resolved_var(tmp_path, grch37_store)
    _, calls = frames(resolved)
    rows = {(r["chrom"], r["source_ref"], r["pos"]): r for r in calls.iter_rows(named=True)}
    fields = ("ref", "alts", "gt", "gt_alleles", "end_pos", "ref_check")

    def row(key: tuple) -> tuple:
        return tuple(rows[key][name] for name in fields)

    x = base("1", P_HET)
    assert row(("1", "", P_HET)) == (x, [x + "GG"], [0, 1], [x, x + "GG"], P_HET, None)
    y = base("1", P_HOM)
    assert row(("1", "", P_HOM)) == (y, [y + "TCA"], [1, 1], [y + "TCA"] * 2, P_HOM, None)
    # VCF's rule before a contig's first base: the base after is the anchor, at POS 1.
    z = base("2", 1)
    assert row(("2", "", 1)) == (z, ["GG" + z], [0, 1], [z, "GG" + z], 1, None)
    # Nothing called at the point: no allele to anchor, and REF X with no ALT would read as a no-call
    # at the base X. Left as loaded, and so is a point on a contig the reference lacks.
    assert row(("1", "", P_NONE)) == (None, [], [None, None], [None, None], None, "not_run")
    assert row(("3", "", 50)) == (None, ["A"], [0, 1], ["", "A"], None, "not_run")
    # A deletion's REF is read and checked as before; `normalize` anchors it.
    d = base("1", P_DEL)
    assert row(("1", d, P_DEL)) == (d, [""], [0, 1], [d, ""], P_DEL, "match")

    # Graded evidence names the alleles `gt` indexes: the anchored ones, REF the anchor base.
    evidence = genotype_evidence(resolved).filter((pl.col("chrom") == "1") & (pl.col("pos") == P_HET)).collect()
    assert evidence["genotype_alleles"].to_list() == [[x, x + "GG"]]

    warnings = "\n".join(resolved.warnings)
    assert "3 insertion point(s) given their anchor base from the reference (1 at a contig's start, anchored on the base after)" in warnings, warnings
    assert "1 insertion point(s) with no haplotype called: REF left unread" in warnings, warnings
    assert "reaching past their contig" not in warnings, "an insertion point is not a row past its contig"


def test_masterVarBeta_insertion_points_are_anchored_as_var_ones(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    var = resolved_var(tmp_path, grch37_store).calls.collect()
    master = resolve_ref(load(write_cgi(tmp_path, master_var_text(point_loci()), "masterVarBeta-GS000000001-ASM.tsv.bz2")), grch37_store)
    calls = master.calls.collect()
    key = ["chrom", "pos", "ref", "alts", "end_pos"]
    assert calls.select(key).sort(key).equals(var.select(key).sort(key))
    # The read counts stay indexed over the alleles as the row now spells them.
    counted = calls.filter(pl.col("source_alleles").is_not_null() & (pl.col("source_ref") == "") & pl.col("ref").is_not_null())
    assert counted.height > 0
    for row in counted.iter_rows(named=True):
        assert row["source_alleles"] == [row["ref"], *row["alts"]]


def test_a_split_anchored_point_keeps_its_read_counts(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """Two insertions at one point (masterVarBeta, with read counts): `normalize` splits them, and each
    row keeps `ad` sliced to its REF and ALT, the counts being indexed over the anchored alleles."""
    p = next(q for q in range(100, 2900) if base("1", q) in "CG")  # neither insertion can then move left
    path = write_cgi(tmp_path, master_var_text([cgi_point("chr1", p, "A", "AT")]), "masterVarBeta-GS000000001-ASM.tsv.bz2")
    resolved = resolve_ref(load(path), grch37_store)
    parent = resolved.calls.collect().row(0, named=True)
    x = base("1", p)
    assert (parent["ref"], parent["alts"], parent["source_alleles"]) == (x, [x + "A", x + "AT"], [x, x + "A", x + "AT"])
    assert len(parent["ad"]) == 3
    split = normalize(resolved, grch37_store).calls.collect()
    got = {row["alts"][0]: row for row in split.iter_rows(named=True)}
    assert set(got) == {x + "A", x + "AT"}
    for k, alt in enumerate(parent["alts"], 1):
        assert got[alt]["ad"] == [parent["ad"][0], parent["ad"][k]], alt
        assert got[alt]["source_alleles"] == [x, alt]


def test_an_anchored_insertion_point_is_minted_and_written(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    resolved = resolved_var(tmp_path, grch37_store)
    anchored = (pl.col("source_ref") == "") & pl.col("ref").is_not_null()
    minted = mint_vrs_ids(resolved, grch37_store).calls.filter(anchored).collect()
    assert minted.height == 3
    assert all(len(ids) == 2 and None not in ids for ids in minted["vrs_allele_ids"].to_list())
    normalized = normalize(resolved, grch37_store)
    written = write_vcf(normalized, tmp_path / "points.vcf.gz", resolved_only=True)
    back = {(r["chrom"], r["pos"], r["source_ref"], tuple(r["alts"])) for r in load(written.path).calls.collect().iter_rows(named=True)}
    read = normalized.calls.filter(pl.col("ref").is_not_null()).collect()
    assert back == {(r["chrom"], r["pos"], r["ref"], tuple(r["alts"])) for r in read.iter_rows(named=True)}
    z = base("2", 1)
    assert ("2", 1, z, ("GG" + z,)) in back
    assert len(back) == 4, "three insertion points and the deletion"


def test_an_anchored_cgi_insertion_moves_left_through_a_repeat(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """`normalize` aligns the anchored row as it does any insertion whose REF was read."""
    run = next(p for p in range(500, 2900) if base("1", p) == base("1", p + 1) != base("1", p - 1))
    b = base("1", run)
    end = run + 1  # the run's last base; the inserted copy follows it, 0-based `begin = end`
    while base("1", end + 1) == b:
        end += 1
    resolved = resolve_ref(load(write_cgi(tmp_path, var_text([cgi_point("chr1", end, "", b)]), "var-GS000000001-ASM.tsv.bz2")), grch37_store)
    row = resolved.calls.collect().row(0, named=True)
    assert (row["pos"], row["ref"], row["alts"]) == (end, b, [b + b])
    out = normalize(resolved, grch37_store).calls.collect().row(0, named=True)
    a = base("1", run - 1)
    assert (out["pos"], out["ref"], out["alts"], out["gt_alleles"]) == (run - 1, a, [a + b], [a, a + b])


# --- Genos -------------------------------------------------------------------------------------------

GENOS_HEAD = "##FileFormat=Genos\n##ReferenceGenome=UCSC.hg19\n##Columns=CHROM\tBEGIN\tEND\tID\tGENOTYPE\n"


def genos(tmp_path: Path, body: str) -> Loaded:
    path = tmp_path / "genos.tsv"
    path.write_text(GENOS_HEAD + body)
    return load(path)


def shifted(bases: str) -> str:
    """Every base replaced by another: never the reference over the same span."""
    return "".join(other(b) for b in bases)


def third(bases: str) -> str:
    """Every base replaced by one that is neither it nor `shifted`'s."""
    return "".join(next(c for c in "ACGT" if c not in (b, other(b))) for b in bases)


def test_genos_rows_get_ref_over_their_span_and_insertion_sites_their_anchor(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """Genos states each row's span (0-based half-open BEGIN/END) and its alleles on the plus strand,
    and no REF. REF is read over the span; ALT are the called alleles other than REF, sorted. An
    insertion site (`BEGIN == END`) is anchored on the base before, or after at a contig's start."""
    two = base("1", 200, 2)
    x = shifted(two)
    x3 = shifted(base("1", 300, 2))
    lo, hi = sorted([shifted(base("1", 700, 2)), third(base("1", 700, 2))])
    p_ins = next(p for p in range(500, 600) if base("1", p) != "A")  # an inserted A cannot move left past it
    p_del = next(p for p in range(400, 500) if base("1", p - 1) != base("1", p))  # nor can this deletion
    body = (
        f"chr1\t99\t100\trs1\t{base('1', 100)}/{other(base('1', 100))}\n"
        f"chr1\t199\t201\trs2\t{two}/{x}\n"
        f"chr1\t299\t301\trs3\t{x3}/{x3}\n"
        f"chr1\t{p_del - 1}\t{p_del}\trs4\t{base('1', p_del)}/-\n"
        f"chr1\t{p_ins}\t{p_ins}\trs5\t-/A\n"
        f"chr1\t600\t600\trs6\tAC/AC\n"
        f"chr1\t699\t701\trs7\t{hi}/{lo}\n"
        f"chr2\t0\t0\trs8\t-/G\n"
        f"chr2\t0\t1\trs9\t{base('2', 1)}/-\n"
        "chr3\t50\t50\trs10\t-/A\n"  # a contig the reference lacks
        "chr3\t59\t61\trs11\tAC/GT\n"
        f"chr1\t800\t800\trs12\t{base('1', 800)}/{base('1', 800)}\n"  # one base inserted, the anchor's own
    )
    resolved = resolve_ref(genos(tmp_path, body), grch37_store)
    _, calls = frames(resolved)
    got = {r["ids"][0]: r for r in calls.iter_rows(named=True)}
    fields = ("pos", "end_pos", "ref", "alts", "gt", "gt_alleles")

    def row(rsid: str) -> tuple:
        return tuple(got[rsid][name] for name in fields)

    b1 = base("1", 100)
    assert row("rs1") == (100, 100, b1, [other(b1)], [0, 1], [b1, other(b1)])
    assert row("rs2") == (200, 201, two, [x], [0, 1], [two, x])
    assert row("rs3") == (300, 301, base("1", 300, 2), [x3], [1, 1], [x3, x3])
    d = base("1", p_del)
    assert row("rs4") == (p_del, p_del, d, [""], [0, 1], [d, ""])
    a = base("1", p_ins)
    assert row("rs5") == (p_ins, p_ins, a, [a + "A"], [0, 1], [a, a + "A"])
    c = base("1", 600)
    assert row("rs6") == (600, 600, c, [c + "AC"], [1, 1], [c + "AC"] * 2)
    assert row("rs7") == (700, 701, base("1", 700, 2), [lo, hi], [2, 1], [hi, lo])
    z = base("2", 1)
    assert row("rs8") == (1, 1, z, ["G" + z], [0, 1], [z, "G" + z])
    assert row("rs9") == (1, 1, z, [""], [0, 1], [z, ""])
    assert row("rs10") == (50, None, None, None, None, ["", "A"])
    assert row("rs11") == (60, 61, None, None, None, ["AC", "GT"])
    # Single bases at an insertion site are inserted bases, not a call at the anchor: never hom-ref.
    e = base("1", 800)
    assert row("rs12") == (800, 800, e, [e + e], [1, 1], [e + e] * 2)
    assert "2 row(s) left with REF unread, on contig(s) the reference does not have: 3" in "\n".join(resolved.warnings)
    assert all(r["source_ref"] is None and r["ref_check"] is None for r in got.values())
    assert got["rs5"]["source_gt"] == "-/A"
    warnings = "\n".join(resolved.warnings)
    assert "4 insertion point(s) given their anchor base from the reference (1 at a contig's start" in warnings, warnings
    assert "array row(s) where the reference base is not A, C, G or T" not in warnings, warnings

    # `normalize` reaches them: a deletion is anchored on the base before (after, at a contig's start),
    # the insertions stay put.
    out = {r["ids"][0]: r for r in normalize(resolved, grch37_store).calls.collect().iter_rows(named=True)}
    assert (out["rs4"]["pos"], out["rs4"]["ref"], out["rs4"]["alts"]) == (p_del - 1, base("1", p_del - 1, 2), [base("1", p_del - 1)])
    assert (out["rs5"]["pos"], out["rs5"]["ref"], out["rs5"]["alts"]) == (p_ins, a, [a + "A"])
    assert (out["rs9"]["pos"], out["rs9"]["ref"], out["rs9"]["alts"]) == (1, base("2", 1, 2), [base("2", 2)])
    assert all("" not in r["alts"] for r in out.values() if r["alts"] is not None)


# --- "No insertion here" ----------------------------------------------------------------------------


def vcf_records(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines() if not line.startswith("#")]


def test_a_point_that_says_no_insertion_gets_no_genotype_at_the_base_before(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """A point whose every called haplotype is the empty allele (CGI `""/""` or `""/?`, Genos `-/-`)
    says nothing is inserted between two bases. Anchored, it became `REF=X, ALT=.`, `0/0`: a hom-ref
    call at the base before, which the source never made, written beside a real SNV call at that base.
    No frame can carry a statement about a gap without a genotype at a base (a `calls` row at `pos` is
    one; a block spans bases), so such a point is left as loaded, counted under its own reason. A point
    with an inserted allele on one haplotype keeps its anchor: its `0` is the haplotype without it."""
    s, h = 200, 300
    b = base("1", s)
    loci = [
        Locus("chr1", s - 1, s, haps=((Seg(s - 1, s, "snp", b, other(b)),), (Seg(s - 1, s, "snp", b, other(b)),))),
        cgi_point("chr1", s, "", ""),
        Locus("chr1", h, h, haps=((Seg(h, h, "ref", "", ""),), (Seg(h, h, "no-call", "", "?"),))),
        cgi_point("chr1", P_HET, "", "GG"),
    ]
    resolved = resolve_ref(load(write_cgi(tmp_path, var_text(loci), "var-GS000000001-ASM.tsv.bz2")), grch37_store)
    _, calls = frames(resolved)
    points = calls.filter(pl.col("source_ref") == "")
    fields = ("pos", "ref", "alts", "gt", "gt_alleles", "end_pos", "ref_check")
    got = {row["pos"]: tuple(row[name] for name in fields) for row in points.iter_rows(named=True)}
    assert got[s] == (s, None, [], [0, 0], ["", ""], None, "not_run")
    assert got[h] == (h, None, [], [0, None], ["", None], None, "not_run")
    x = base("1", P_HET)
    assert got[P_HET] == (P_HET, x, [x + "GG"], [0, 1], [x, x + "GG"], P_HET, None)
    warnings = "\n".join(resolved.warnings)
    assert "1 insertion point(s) given their anchor base" in warnings, warnings
    assert "2 insertion point(s) where every called haplotype is the empty allele" in warnings, warnings
    assert "no haplotype called" not in warnings, warnings

    written = write_vcf(normalize(resolved, grch37_store), tmp_path / "points.vcf", compress=False, resolved_only=True)
    records = vcf_records(written.path)
    at = [(r[3], r[4], r[9].split(":")[0]) for r in records if int(r[1]) == s]
    assert at == [(b, other(b), "1/1")], "only the SNV is written at the base before the point"
    assert not [r for r in records if r[4] == "."], "no record states a hom-ref call the source never made"
    minted = mint_vrs_ids(resolved, grch37_store).calls.filter((pl.col("source_ref") == "") & (pl.col("pos") == s)).collect()
    assert minted["vrs_allele_ids"].to_list() == [None]


def test_a_genos_insertion_site_that_says_no_insertion_gets_no_genotype(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """Genos `-/-` at an insertion site, beside a hom-alt SNV at the base before: the SNV is read, the
    site is left as loaded, counted once and not as an array row, and nothing reads hom-ref at that base."""
    b = base("1", 500)
    body = f"chr1\t499\t500\trsA\t{other(b)}/{other(b)}\nchr1\t500\t500\trsB\t-/-\n"
    resolved = resolve_ref(genos(tmp_path, body), grch37_store)
    got = {r["ids"][0]: r for r in resolved.calls.collect().iter_rows(named=True)}
    assert (got["rsA"]["ref"], got["rsA"]["alts"], got["rsA"]["gt"]) == (b, [other(b)], [1, 1])
    fields = ("pos", "end_pos", "ref", "alts", "gt", "gt_alleles")
    assert tuple(got["rsB"][name] for name in fields) == (500, None, None, None, None, ["", ""])
    warnings = "\n".join(resolved.warnings)
    assert "1 insertion point(s) where every called haplotype is the empty allele" in warnings, warnings
    assert "for 1 of 1 array row(s)" in warnings, warnings
    written = write_vcf(normalize(resolved, grch37_store), tmp_path / "genos.vcf", compress=False, resolved_only=True)
    assert [(r[2], r[3], r[4], r[9]) for r in vcf_records(written.path)] == [("rsA", b, other(b), "1/1")]
