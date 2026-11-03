"""F23: a source REF made only of `N` states no base (`ref_check = unstated`).

DRAGEN's and Canvas's copy-number VCFs write REF `N` on every record: the caller never reads the
reference. `resolve_ref` reads the base and records that the source stated none, distinct from `match`
(a check that ran and agreed) and `mismatch` (the source stated another base). A partial `N` (`AN`)
states a base, so it is compared as before. Every consumer of `ref_check` is exercised here on
synthetic inline VCFs, and the opt-in test at the end counts the PGP-HMS DRAGEN CNV files.
"""

import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import (
    ReferenceStore,
    cis_phased_blocks,
    concordance,
    consensus,
    genotype_evidence,
    genotype_evidence_at,
    liftover,
    load,
    mint_vrs_ids,
    normalize,
    query,
    reference_support,
    rejoin_split,
    resolve_ref,
    vrs_withheld,
    write_vcf,
)
from just_vcf.concordance import EXCLUDED_REASONS
from just_vcf.evidence import GradeBasis
from just_vcf.karyotype import measure
from just_vcf.liftover import UNLIFTED_REASONS
from just_vcf.loading import Loaded
from just_vcf.schema import CALLS
from just_vcf.site_genotypes import SiteBasis, panel_sites, read_site_genotypes
from just_vcf.site_genotypes import facet_of as site_facet
from just_vcf.vcf_zarr import write_vcf_zarr
from just_vcf.vrs import WITHHELD_REASONS
from just_vcf.writing import UnwritableRowsError
from test_cli import report, run
from test_phasing import Record, phased_vcf
from test_liftover import MULTI_BASE, MULTI_SWAP_POS, TARGET, by_case, expected_sv, s1, setup, sv  # noqa: F401 (setup is a fixture)
from test_liftover import SV_HEADER, other as lift_other
from test_vcf_zarr import read_back

_rng = random.Random(23)
N_RUN = 1201
"""1-based start of ten `N` on contig 1: a gap, where the reference states no base either."""
SEQ = {
    "1": "".join(_rng.choice("ACGT") for _ in range(N_RUN - 1)) + "N" * 10 + "".join(_rng.choice("ACGT") for _ in range(3000 - N_RUN + 1 - 10)),
}


def base(pos: int, length: int = 1) -> str:
    return SEQ["1"][pos - 1 : pos - 1 + length]


def other(b: str) -> str:
    return next(x for x in "ACGT" if x != b)


def write_fasta(path: Path) -> Path:
    path.write_text(">1\n" + "\n".join(SEQ["1"][i : i + 60] for i in range(0, len(SEQ["1"]), 60)) + "\n")
    return path


def vcf(tmp_path: Path, records: list[tuple[int, str, str, str, str]], name: str = "calls.vcf") -> Path:
    """A one-sample VCF on contig 1 of `(pos, id, ref, alt, gt)` records; the id names the case."""
    meta = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={len(SEQ['1'])}>",
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"1\t{p}\t{i}\t{r}\t{a}\t50\tPASS\t{'END=' + str(p + 20) if a.startswith('<') else '.'}\tGT\t{g}" for p, i, r, a, g in records]
    path = tmp_path / name
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa"))
    return s


# The id names the case and what `ref_check` should read after `resolve_ref`.
RECORDS = [
    (10, "match", base(10), other(base(10)), "0/1"),
    (20, "mismatch", other(base(20)), base(20), "0/1"),
    (30, "unstated_del", "N", "<DEL>", "0/1"),
    (60, "unstated_lower_dup", "n", "<DUP>", "0/1"),
    (90, "unstated_two_bases", "NN", ".", "0/0"),
    (120, "unstated_snv", "N", other(base(120)), "0/1"),
    (150, "unstated_insertion", "N", "NAC", "0/1"),
    (240, "unstated_snv_hom_ref", "N", other(base(240)), "0/0"),
    (180, "partial_after", base(180) + "N", "<DEL>", "0/1"),
    (210, "partial_before", "N" + base(211), "<DEL>", "0/1"),
    (N_RUN + 2, "gap", "N", "A", "0/1"),
]
EXPECTED = {
    "match": "match",
    "mismatch": "mismatch",
    "unstated_del": "unstated",
    "unstated_lower_dup": "unstated",
    "unstated_two_bases": "unstated",
    "unstated_snv": "unstated",
    "unstated_insertion": "unstated",
    "unstated_snv_hom_ref": "unstated",
    # A real base beside an N states that base: compared as before, and the reference has no N there.
    "partial_after": "mismatch",
    "partial_before": "mismatch",
    # Over a gap the reference reads N too: the letters agree, so the comparison ran and matched.
    "gap": "match",
}


UNSTATED_SEQUENCE_ALT = {case for _, case, _, alt, _ in RECORDS if EXPECTED[case] == "unstated" and alt != "." and not alt.startswith("<")}
"""Unstated rows with a sequence ALT, spelled against the `N`: no reader takes that ALT as stated."""


def resolved(tmp_path: Path, store: ReferenceStore) -> Loaded:
    return resolve_ref(load(vcf(tmp_path, RECORDS)), store)


def by_id(frame: pl.DataFrame) -> dict[str, dict]:
    return {r["ids"][0]: r for r in frame.iter_rows(named=True)}


def test_a_source_ref_of_only_n_is_unstated_and_counted_apart(tmp_path: Path, store: ReferenceStore) -> None:
    result = resolved(tmp_path, store)
    calls = result.calls.collect()
    rows = by_id(calls)
    assert {case: r["ref_check"] for case, r in rows.items()} == EXPECTED
    assert calls.schema == CALLS and "unstated" in CALLS["ref_check"].categories.to_list()
    # REF is always the reference's bases over the span the source wrote; the source REF is kept as written.
    for pos, case, ref, _, _ in RECORDS:
        assert rows[case]["ref"] == base(pos, len(ref)) and rows[case]["source_ref"] == ref, case
    counts = {value: sum(1 for v in EXPECTED.values() if v == value) for value in ("match", "mismatch", "unstated")}
    line = next(w for w in result.warnings if w.startswith("REF read from collection"))
    assert f"for {len(RECORDS)} of {len(RECORDS)} row(s)" in line
    assert f"{counts['match']} match, {counts['mismatch']} mismatch, {counts['unstated']} with a source REF of only N" in line
    assert "IUPAC" not in line, "a source N is not a reference ambiguity code"


def test_minting_reads_the_ref_slot_and_withholds_an_alt_spelled_against_n(tmp_path: Path, store: ReferenceStore) -> None:
    """The REF slot is the reference's allele, read from the sequence: minted. A sequence ALT on an
    unstated row is spelled against a base the source did not state (`N` → `NAC` read as `C` → `NAC`
    would be a delins), so it is withheld as `ref_unstated`; a symbolic ALT keeps its own reason."""
    assert "ref_unstated" in WITHHELD_REASONS
    rows = by_id(vrs_withheld(mint_vrs_ids(resolved(tmp_path, store), store).calls).collect())
    assert rows["unstated_del"]["vrs_withheld"] == [None, "sv_no_vrs_model"]
    assert rows["unstated_two_bases"]["vrs_withheld"] == [None]
    assert rows["unstated_snv"]["vrs_withheld"] == [None, "ref_unstated"]
    assert rows["unstated_insertion"]["vrs_withheld"] == [None, "ref_unstated"]
    assert rows["partial_after"]["vrs_withheld"] == ["ref_mismatch", "sv_no_vrs_model"]
    assert rows["gap"]["vrs_withheld"] == [None, None]
    # The REF slot's id is the one a source that wrote the reference's base gets at the same place.
    stated = [(120, "stated_snv", base(120), other(base(120)), "0/1"), (150, "stated_ins", base(150), base(150) + "AC", "0/1")]
    twin = by_id(mint_vrs_ids(resolve_ref(load(vcf(tmp_path, stated, "stated.vcf")), store), store).calls.collect())
    assert rows["unstated_snv"]["vrs_allele_ids"][0] == twin["stated_snv"]["vrs_allele_ids"][0]
    assert rows["unstated_insertion"]["vrs_allele_ids"][0] == twin["stated_ins"]["vrs_allele_ids"][0]


def test_normalize_leaves_an_unstated_row_as_written_and_says_why(tmp_path: Path, store: ReferenceStore) -> None:
    source = resolved(tmp_path, store)
    out = normalize(source, store)
    before, after = by_id(source.calls.collect()), by_id(out.calls.collect())
    for case in ("unstated_snv", "unstated_insertion"):
        assert (after[case]["pos"], after[case]["ref"], after[case]["alts"], after[case]["source_ref"]) == (
            before[case]["pos"],
            before[case]["ref"],
            before[case]["alts"],
            before[case]["source_ref"],
        )
    # The insertion only: an SNV has nothing to align, and the symbolic and ALT-less rows are no shape
    # normalize aligns.
    assert any(w == "1 row(s) not left-aligned or trimmed: the source REF states no base (ref_check unstated)" for w in out.warnings)


def unwritable() -> dict[str, int]:
    """What the writers leave out of `RECORDS`: a mismatch, and an unstated row's sequence ALT, which
    written beside the reference's base would read back as a stated `match` (`N` → `NAC` as `C` → `NAC`)."""
    return {"ref_mismatch": sum(1 for v in EXPECTED.values() if v == "mismatch"), "ref_unstated": len(UNSTATED_SEQUENCE_ALT)}


def written_cases() -> set[str]:
    return {case for case, value in EXPECTED.items() if value != "mismatch" and case not in UNSTATED_SEQUENCE_ALT}


def test_write_vcf_writes_an_unstated_row_with_the_reference_base_unless_it_spells_a_sequence_alt(tmp_path: Path, store: ReferenceStore) -> None:
    source = resolved(tmp_path, store)
    with pytest.raises(UnwritableRowsError) as refused:
        write_vcf(source, tmp_path / "refused.vcf", compress=False)
    assert refused.value.counts == unwritable()
    written = write_vcf(source, tmp_path / "out.vcf", compress=False, resolved_only=True)
    assert written.rows_skipped == unwritable()
    back = by_id(load(written.path).calls.collect())
    assert set(back) == written_cases()
    for pos, case, ref, alt, _ in RECORDS:
        if case in back:
            assert back[case]["source_ref"] == base(pos, len(ref)), case
            assert back[case]["alts"] == ([] if alt == "." else [alt]), case


def test_write_vcf_zarr_keeps_the_unstated_rows_without_a_sequence_alt(tmp_path: Path, store: ReferenceStore) -> None:
    source = resolved(tmp_path, store)
    with pytest.raises(UnwritableRowsError) as refused:
        write_vcf_zarr(source, tmp_path / "refused.vcz")
    assert refused.value.counts == unwritable()
    written = write_vcf_zarr(source, tmp_path / "out.vcz", resolved_only=True)
    assert written.rows_skipped == unwritable()
    sites = read_back(tmp_path / "out.vcz")
    kept = sorted(pos for pos, case, *_ in RECORDS if case in written_cases())
    assert sites["start"].to_list() == kept
    refs = dict(zip(sites["start"].to_list(), sites["ref"].to_list()))
    assert all(refs[pos] == base(pos, len(ref)) for pos, case, ref, *_ in RECORDS if pos in refs)


def test_query_counts_unstated_and_reads_neither_site_nor_alt_from_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """An unstated row never says which base the site holds, so an SNV at its POS that it does not spell
    is not answered `site_genotyped`; and its sequence ALT was spelled against the `N`, not against the
    reference's base, so the same SNV that it does spell is not `stated` either (a `0/0` there would read
    `hom_ref`, a base the source never named). Both are `overlapping_call`, carried unknown."""
    q = query(resolved(tmp_path, store))
    (table,) = q.counts(by=["ref_check"])
    assert {group.values[0]: group.count for group in table.groups} == {
        value: sum(1 for v in EXPECTED.values() if v == value) for value in set(EXPECTED.values())
    }
    (spelled,) = q.carries(chrom="1", pos=120, ref=base(120), alt=other(base(120)))
    assert (spelled.basis.value, spelled.carried, spelled.genotype_class) == ("overlapping_call", None, None)
    (hom_ref,) = q.carries(chrom="1", pos=240, ref=base(240), alt=other(base(240)))
    assert (hom_ref.basis.value, hom_ref.carried, hom_ref.genotype_class) == ("overlapping_call", None, None)
    elsewhere = next(b for b in "ACGT" if b not in (base(120), other(base(120))))
    (unstated,) = q.carries(chrom="1", pos=120, ref=base(120), alt=elsewhere)
    assert (unstated.basis.value, unstated.carried) == ("overlapping_call", None)
    (matched,) = q.carries(chrom="1", pos=10, ref=base(10), alt=next(b for b in "ACGT" if b not in (base(10), other(base(10)))))
    assert (matched.basis.value, matched.carried) == ("site_genotyped", False), "a stated REF still genotypes its site"


def test_the_cli_counts_unstated_rows(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "ref.fa")
    [facet] = report(run("check", vcf(tmp_path, RECORDS), "--reference", fasta, "--json"))["facets"]
    assert facet["ref_check"] == {value: sum(1 for v in EXPECTED.values() if v == value) for value in set(EXPECTED.values())}


def third(b: str) -> str:
    """A base that is neither `b` nor `other(b)`."""
    return next(x for x in "ACGT" if x not in (b, other(b)))


def facet_of(tmp_path: Path, store: ReferenceStore, name: str, records: list[tuple[int, str, str, str, str]], **columns: object) -> Loaded:
    """`records` loaded with REF read, and any `facets` column set as given."""
    return with_facet(resolve_ref(load(vcf(tmp_path, records, f"{name}.vcf")), store), **columns)


def test_concordance_compares_no_alt_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """`a` writes the reference's base, `b` writes `N` at the same three sites. `b`'s sequence ALTs were
    spelled against the `N`, so none is compared (`ref_unstated`, counted), and `b`'s row there is not
    read for copies of `a`'s allele either: `a`'s allele is `only_a`, `b`'s silence `other_alleles`.
    Read beside the reference's base, `N` → `X` would match `a`'s `T` → `X`, a concordance the source
    never stated; the twin `b` that writes the base shows what the `N` alone changes."""
    same, hom, differ = 300, 400, 500
    a = facet_of(tmp_path, store, "a", [(p, f"a{p}", base(p), other(base(p)), g) for p, g in ((same, "1/1"), (hom, "0/0"), (differ, "1/1"))])
    alts = {same: other(base(same)), hom: other(base(hom)), differ: third(base(differ))}
    gts = {same: "1/1", hom: "0/0", differ: "1/1"}
    b = facet_of(tmp_path, store, "b", [(p, f"b{p}", "N", alts[p], gts[p]) for p in (same, hom, differ)])
    twin = facet_of(tmp_path, store, "twin", [(p, f"t{p}", base(p), alts[p], gts[p]) for p in (same, hom, differ)])

    # The twin's REF states the base: the same records compare (the `differ` one by reading its row, `site`).
    a_alleles = {(r["pos"], r["agreement"]) for r in concordance(a, twin).sites.collect().iter_rows(named=True) if r["alt"] == other(base(r["pos"]))}
    assert a_alleles == {(same, "match"), (hom, "match"), (differ, "mismatch")}
    result = concordance(a, b)
    sites = result.sites.collect()
    assert sorted(sites["pos"].to_list()) == [same, hom, differ], "only a's alleles: none of b's is compared"
    for r in sites.iter_rows(named=True):
        assert (r["alt"], r["agreement"], r["key"], r["absence"], r["b_gt_alleles"], r["b_copies"]) == (other(base(r["pos"])), "only_a", None, "other_alleles", None, None), r["pos"]
    assert "ref_unstated" in EXCLUDED_REASONS
    assert {(r["facet"], r["reason"]): r["alleles"] for r in result.excluded.iter_rows(named=True)} == {("b", "ref_unstated"): 3}


def test_concordance_matches_a_symbolic_alt_on_an_unstated_row_as_before(tmp_path: Path, store: ReferenceStore) -> None:
    """A symbolic ALT on an unstated row (a copy-number call) meets by type and extent as it does over a
    stated base: its sites are the ones the same records with the reference's base give. On a row with
    a sequence ALT beside it (`N` → `<DEL>,X`) only the sequence ALT is left out."""
    del_pos, dup_pos, mixed = 600, 700, 800
    a = facet_of(tmp_path, store, "a", [(del_pos, "d", base(del_pos), "<DEL>", "0/1"), (dup_pos, "u", base(dup_pos), "<DUP>", "1/1"), (mixed, "m", base(mixed), "<DEL>", "0/1")])

    def b_with(ref: dict[int, str], name: str) -> Loaded:
        return facet_of(
            tmp_path, store, name,
            [(del_pos, "d", ref[del_pos], "<DEL>", "0/1"), (dup_pos, "u", ref[dup_pos], "<DUP>", "1/1"), (mixed, "m", ref[mixed], f"<DEL>,{other(base(mixed))}", "0/1")],
        )  # fmt: skip

    unstated = concordance(a, b_with({p: "N" for p in (del_pos, dup_pos, mixed)}, "b_n"))
    stated = concordance(a, b_with({p: base(p) for p in (del_pos, dup_pos, mixed)}, "b_base"))
    sv = pl.col("variant_class") == "sv"
    sv_sites = unstated.sites.collect().filter(sv)
    # The genotype as written is the one difference: its REF slot is the source's `N`.
    assert sv_sites.drop("b_gt_alleles").equals(stated.sites.collect().filter(sv).drop("b_gt_alleles"))
    assert [gt[0] for gt in sv_sites["b_gt_alleles"].to_list()] == ["N", "<DUP>", "N"]
    assert sv_sites.select("pos", "key", "agreement").rows() == [(del_pos, "extent", "match"), (dup_pos, "extent", "match"), (mixed, "extent", "match")]
    assert unstated.sites.collect().filter(~sv).height == 0, "the mixed row's sequence ALT is not compared"
    assert {(r["facet"], r["reason"]): r["alleles"] for r in unstated.excluded.iter_rows(named=True)} == {("b", "ref_unstated"): 1}
    assert stated.excluded.filter(pl.col("reason") == "ref_unstated").height == 0


def test_consensus_states_no_alt_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """`y` writes `N` where `x` writes the reference's base. `y`'s sequence ALTs are no allele of the
    consensus (none stated by it alone becomes a row), and at a site `x` states, `y`'s row is not read
    for copies: `y` is silent there as `ref_unstated`. It is still a row at the position, so a
    whole-genome `y` is not read as reference there (its absence rule would give zero copies, a
    contradiction from a base `y` never named)."""
    e = 0.01
    snv_x = other(base(300)), other(base(400))
    x = facet_of(tmp_path, store, "x", [(300, "x300", base(300), snv_x[0], "1/1"), (400, "x400", base(400), snv_x[1], "0/1")], caller="X", caller_version="1")
    wgs = {"absence_means": "weak_reference", "variant_classes": [{"variant_class": "snv", "alleles": 3}]}
    y_records = [(300, "y300", "N", snv_x[0], "1/1"), (400, "y400", "N", third(base(400)), "0/0"), (500, "y500", "N", other(base(500)), "0/1")]
    y = facet_of(tmp_path, store, "y", y_records, caller="Y", caller_version="1", **wgs)
    result = consensus([x, y], atlas=None, default_rate=e)
    support = result.support.collect()
    assert support.select("pos", "alt").rows() == [(300, snv_x[0]), (400, snv_x[1])], "no allele y spelled against its N"
    assert support.select("supported", "contradicted", "silent", "facets").rows() == [(1, 0, 1, 2), (1, 0, 1, 2)]
    y_id = y.facets.collect()["facet_id"].item()
    assert f"facet {y_id} silent on alleles others stated: ref_unstated 2" in result.warnings
    # x alone speaks: a hard call under the default rate, flat prior.
    for row, copies in zip(result.calls.collect().iter_rows(named=True), (2, 1)):
        alone = [1 - e if k == copies else e / 2 for k in range(3)]
        assert row["gp"] == pytest.approx([v / sum(alone) for v in alone], rel=1e-5)


# Lifted through test_liftover's construction (block A: source 101..400 at -50; 701..800 reversed onto
# target 2; SWAP_POS 251 is a base the target changes).
LIFT_CASES = {
    "cnv_del_plus": (sv("cnv_del_plus", 201, "<DEL>", 230, ref="N"), None),
    "cnv_dup_minus": (sv("cnv_dup_minus", 720, "<DUP>", 760, ref="N"), None),
    "cnv_segment": (sv("cnv_segment", 460, ".", 520, sample="0/0:2", ref="N"), None),
    # The lift's own check still runs: on a plus strand the target's base must be the base read at the source.
    "cnv_anchor_changed": (sv("cnv_anchor_changed", 251, "<DEL>", 280, ref="N"), "ref_mismatch_after_lift"),
    "cnv_partial": (sv("cnv_partial", 340, "<DEL>", 360, ref=s1(340) + "N"), "source_ref_mismatch"),
    # A sequence ALT on an unstated row was spelled against the `N`: the lift would compare the target's
    # base with it, swap, and re-index the genotype on an allele the source never set against a base.
    "unstated_snv": (f"1\t205\tunstated_snv\tN\t{lift_other(s1(205))}\t50\tPASS\t.\tGT\t0/1", "source_ref_unstated"),
    # REF N, ALT G, GT 0/0 over source base A where the target reads G (a base the target changes;
    # `cnv_anchor_changed` holds the other such position): lifted, it would be REF G ALT A GT 1/1, an
    # allele the source never called.
    "unstated_swap": (f"1\t{MULTI_SWAP_POS}\tunstated_swap\tN\t{MULTI_BASE}\t50\tPASS\t.\tGT\t0/0", "source_ref_unstated"),
}


def test_rejoin_leaves_a_split_site_spelled_against_an_unstated_ref_as_written(tmp_path: Path, store: ReferenceStore) -> None:
    """A split site whose rows' source REF is `N` (`N` → `X` and `N` → `Y`, one row per ALT) joined into
    one row with `ref_check = not_run`, which every reader reads as written: its sequence ALTs, spelled
    against the `N`, became readable again, and `resolve_ref` refuses a facet already read, so nothing
    could set the row back to `unstated`. The group is now left as written, counted under its own
    reason, as the other groups that cannot be read as one site are; a stated split site beside it
    still joins."""
    tag = "OLD_MULTIALLELIC=1:{}:A/G/T"
    x, y = other(base(300)), third(base(300))
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={len(SEQ['1'])}>",
        '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        f"1\t300\t.\tN\t{x}\t50\tPASS\t{tag.format(300)}\tGT\t1/0",
        f"1\t300\t.\tN\t{y}\t50\tPASS\t{tag.format(300)}\tGT\t0/1",
        f"1\t400\t.\t{base(400)}\t{other(base(400))}\t50\tPASS\t{tag.format(400)}\tGT\t1/0",
        f"1\t400\t.\t{base(400)}\t{third(base(400))}\t50\tPASS\t{tag.format(400)}\tGT\t0/1",
    ]
    path = tmp_path / "split.vcf"
    path.write_text("\n".join(lines) + "\n")
    rejoined = rejoin_split(resolve_ref(load(path), store))
    calls = {pos: group for (pos,), group in rejoined.calls.collect().group_by(["pos"])}
    assert calls[300].select("source_ref", "ref_check", "alts", "split_signal").sort("alts").rows() == [
        ("N", "unstated", [x], "marker"),
        ("N", "unstated", [y], "marker"),
    ]
    assert calls[400]["alts"].to_list() == [sorted([other(base(400)), third(base(400))])]
    assert "1 split site(s) left as written: a source REF of only N, which states no base (ref_check unstated)" in rejoined.warnings


def test_liftover_lifts_unstated_rows_like_matching_ones_unless_they_spell_a_sequence_alt(tmp_path: Path, setup) -> None:  # noqa: F811
    assert "source_ref_unstated" in UNLIFTED_REASONS
    store, _, _, chain = setup
    lines = sorted((line for line, _ in LIFT_CASES.values()), key=lambda line: int(line.split("\t")[1]))
    path = tmp_path / "cnv.vcf"
    path.write_text("\n".join([*SV_HEADER, *lines]) + "\n")
    source = resolve_ref(load(path), store)
    parent = by_case(source.calls.collect())
    assert {case: r["ref_check"] for case, r in parent.items()} == {
        case: ("mismatch" if case == "cnv_partial" else "unstated") for case in LIFT_CASES
    }
    lifted = liftover(source, chain, store)
    report_rows = by_case(lifted.report.collect())
    assert {case: r["lift_reason"] for case, r in report_rows.items()} == {case: reason for case, (_, reason) in LIFT_CASES.items()}
    assert report_rows["unstated_swap"]["lift_swapped"] is None, "never compared with the target, so never swapped"
    unstated_alts = sum(1 for _, reason in LIFT_CASES.values() if reason == "source_ref_unstated")
    assert lifted.reasons["source_ref_unstated"] == unstated_alts
    assert any(w.startswith(f"{unstated_alts} row(s) not lifted:") and w.endswith("[source_ref_unstated]") for w in lifted.loaded.warnings)
    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == {case for case, (_, reason) in LIFT_CASES.items() if reason is None}
    for case, r in calls.items():
        p = parent[case]
        chrom, pos, end, _ = expected_sv(p["pos"], p["end_pos"], padded=bool(p["alts"]))
        assert r["end_pos"] == end, case
        # REF is the target's base at the new POS; what the source wrote is not carried (a lifted row's
        # source_ref and ref_check are null); the ALTs stay as written.
        assert (r["chrom"], r["pos"]) == (chrom, pos), case
        assert r["ref"] == TARGET[chrom][pos - 1] and r["source_ref"] is None and r["ref_check"] is None, case
        assert r["alts"] == p["alts"], case


def test_cnv_files_on_pgp_hms_resolve_unstated_and_lift(pgp_participants: Path, test_reference_dir: Path, test_chain_dir: Path) -> None:
    """Opt-in: every `.cnv.` VCF upload on disk (F23: DRAGEN writes REF `N` on every record). Every
    all-`N` row over bases the reference reads is `unstated`, none is `mismatch`, and the lift leaves
    behind as `source_ref_mismatch` exactly the rows whose stated REF disagrees, none of the unstated.
    Aggregates only. Skips without the PGP-HMS tree, a real reference store or the UCSC chains."""
    store = ReferenceStore.open(test_reference_dir)
    files = sorted(pgp_participants.glob("*/upload-*/*.cnv.vcf.gz"))
    if not files:
        pytest.skip("no .cnv. VCF upload under the PGP-HMS tree")
    chains = {"GRCh37": test_chain_dir / "hg19ToHg38.over.chain.gz", "GRCh38": test_chain_dir / "hg38ToHg19.over.chain.gz"}
    all_n = pl.col("source_ref").str.contains(r"^[Nn]+$").fill_null(False)
    over_bases = all_n & pl.col("ref").is_not_null() & (pl.col("source_ref").str.to_uppercase() != pl.col("ref"))
    unstated_rows = 0
    lifted_unstated = 0
    for path in files:
        source = resolve_ref(load(path), store)
        counts = source.calls.select(
            pl.len().alias("rows"),
            over_bases.sum().alias("over_bases"),
            (pl.col("ref_check") == "unstated").sum().alias("unstated"),
            (all_n & (pl.col("ref_check") == "mismatch")).sum().alias("n_mismatch"),
            (pl.col("ref_check") == "mismatch").sum().alias("mismatch"),
        ).collect().row(0, named=True)
        assert counts["unstated"] == counts["over_bases"] and counts["n_mismatch"] == 0, path.name
        unstated_rows += counts["unstated"]
        label = store.assembly_label(source.facets.collect()["ref_seqcol_digest"].item())
        if label not in chains or not chains[label].is_file():
            continue
        lifted = liftover(source, chains[label], store)
        lifted_rows = lifted.loaded.calls.select(pl.len()).collect().item()
        assert lifted.reasons.get("source_ref_mismatch", 0) == counts["mismatch"], path.name
        assert lifted_rows + sum(lifted.reasons.values()) == counts["rows"], path.name
        if counts["unstated"] == counts["rows"]:
            lifted_unstated += lifted_rows
    assert unstated_rows > 0, "no CNV file with a source REF of only N"
    assert lifted_unstated > 0, "no row of an all-N CNV file lifted"


def vcf_with(tmp_path: Path, name: str, keys: str, records: list[tuple[int, str, str, str, str]]) -> Path:
    """`vcf` with FORMAT `keys` (`GT:GP`, `GT:DP`) and each record's sample column written as given."""
    declared = {
        "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "GP": '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posteriors">',
        "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    }
    meta = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={len(SEQ['1'])}>",
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        *(declared[k] for k in keys.split(":")),
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"1\t{p}\t{i}\t{r}\t{a}\t50\tPASS\t{'END=' + str(p + 20) if a.startswith('<') else '.'}\t{keys}\t{s}" for p, i, r, a, s in records]
    path = tmp_path / f"{name}.vcf"
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def with_facet(loaded: Loaded, **columns: object) -> Loaded:
    """`loaded` with any `facets` column set as given."""
    if not columns:
        return loaded
    schema = loaded.facets.collect_schema()
    facets = loaded.facets.with_columns(*[pl.lit(value, dtype=schema[key]).alias(key) for key, value in columns.items()])
    return Loaded(facets=facets, calls=loaded.calls, blocks=loaded.blocks, links=loaded.links)


WGS_SNV = {"absence_means": "weak_reference", "variant_classes": [{"variant_class": "snv", "alleles": 4}, {"variant_class": "sv", "alleles": 1}]}
"""A whole-genome facet's absence columns: its silence reads as reference for an SNV."""

EVIDENCE_RECORDS = [
    (300, "unstated_hom_ref", "N", other(base(300)), "0/0:0.9,0.08,0.02"),
    (400, "unstated_het", "N", other(base(400)), "0/1:."),
    (500, "unstated_del", "N", "<DEL>", "0/1:."),
    (600, "stated_hom_ref", base(600), other(base(600)), "0/0:."),
]


def test_genotype_evidence_grades_no_genotype_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """`N` → `X` with GT `0/0` names no base at the site: its genotype (and its GP vector, which is over
    `[N, X]`) was spelled against the `N`. Read as written it would grade hom-ref at a base the source
    never named. It is `ref_unstated`: no genotype, no vector, unknown, and still a row at its position,
    so `reference_support` answers `call` there and never reads the whole-genome silence as reference.
    A symbolic ALT of such a row is graded as before; a row loaded but not resolved carries
    `ref_check = not_run`, so the rule cannot apply to it yet."""
    path = vcf_with(tmp_path, "evidence", "GT:GP", EVIDENCE_RECORDS)
    loaded = with_facet(resolve_ref(load(path), store), **WGS_SNV)
    pos_of = {case: pos for pos, case, *_ in EVIDENCE_RECORDS}
    graded = {r["pos"]: r for r in genotype_evidence(loaded).collect().iter_rows(named=True)}
    ungraded = {"grade_basis": "ref_unstated", "confidence_kind": "unknown", "genotype_alleles": None, "genotype_probabilities": None,
                "prior": None, "genotype_index": None, "confidence": None}  # fmt: skip
    for case in ("unstated_hom_ref", "unstated_het"):
        assert {k: graded[pos_of[case]][k] for k in ungraded} == ungraded, case
    assert (graded[500]["grade_basis"], graded[500]["genotype_alleles"], graded[500]["genotype_index"]) == ("gt_only", ["N", "<DEL>"], 1)
    assert (graded[600]["grade_basis"], graded[600]["genotype_index"]) == ("gt_only", 0)

    asked = pl.DataFrame({"chrom": ["1"] * 5, "pos": [300, 400, 510, 600, 700]})
    support = {r["pos"]: r for r in reference_support(loaded, asked, variant_class="snv").collect().iter_rows(named=True)}
    assert {p: (r["support"], r["reason"], r["grade_basis"], r["genotype_index"]) for p, r in support.items()} == {
        300: ("call", "call_at_position", "ref_unstated", None),
        400: ("call", "call_at_position", "ref_unstated", None),
        510: ("call", "call_spans_position", "gt_only", 1),
        600: ("call", "call_at_position", "gt_only", 0),
        700: ("absence", "whole_genome_silence", "absent_reference", 0),
    }
    at = genotype_evidence_at(loaded, asked.filter(pl.col("pos") == 300), variant_class="snv")
    assert at.select("evidence_basis", "grade_basis", "genotype_index").rows() == [("assay", "ref_unstated", None)], "a row there, never absence"
    assert GradeBasis("ref_unstated") is GradeBasis.REF_UNSTATED

    # Before `resolve_ref` no REF was compared (`not_run`): the rule cannot apply, and the row is graded as written.
    raw = {r["pos"]: r for r in genotype_evidence(with_facet(load(path), **WGS_SNV)).collect().iter_rows(named=True)}
    assert (raw[300]["grade_basis"], raw[300]["genotype_alleles"], raw[300]["genotype_index"]) == ("gp", ["N", other(base(300))], 0)


KARYOTYPE_RECORDS = [
    (300, "unstated_het", "N", other(base(300)), "0/1:30"),
    (400, "unstated_hom", "N", other(base(400)), "1/1:30"),
    (500, "unstated_haploid", "N", other(base(500)), "1:30"),
    (600, "unstated_hom_ref", "N", other(base(600)), "0/0:30"),
    (700, "unstated_del", "N", "<DEL>", "0/1:30"),
    (800, "stated_het", base(800), other(base(800)), "0/1:30"),
]
_ALT_COUNTS = ("diploid_alt", "het", "haploid_alt", "depth_sites", "median_dp", "reference_sites")
_PLOIDY_COUNTS = ("rows", "one", "two", "more")


def test_sex_chromosome_measures_read_no_alt_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """`measure` counts the ALT-carrying calls its het share and depth lines read. An unstated row's
    sequence ALT was spelled against the `N` (`N` → `X`, `0/1`, says nothing of which two bases the
    site holds), so it carries no ALT there: its counts are those of the rows without it. The ploidy as
    written is the caller's and counts as before. A symbolic ALT of such a row (`<DEL>`, `0/1`) counts
    as it always has, as any symbolic call does; before `resolve_ref` the rule cannot apply."""

    def counts(records: list[tuple[int, str, str, str, str]], name: str, resolved_ref: bool = True) -> dict:
        loaded = load(vcf_with(tmp_path, name, "GT:DP", records))
        return measure((resolve_ref(loaded, store) if resolved_ref else loaded).calls, None)["groups"]["autosomes"]

    unstated = counts(KARYOTYPE_RECORDS, "unstated")
    read = counts([r for r in KARYOTYPE_RECORDS if r[1] in ("unstated_del", "stated_het", "unstated_hom_ref")], "read")
    twin = counts([(p, i, base(p) if ref == "N" else ref, alt, s) for p, i, ref, alt, s in KARYOTYPE_RECORDS], "twin")
    assert {k: unstated[k] for k in _ALT_COUNTS} == {k: read[k] for k in _ALT_COUNTS}
    assert {k: unstated[k] for k in _PLOIDY_COUNTS} == {k: twin[k] for k in _PLOIDY_COUNTS}
    assert (unstated["diploid_alt"], unstated["het"], unstated["haploid_alt"]) == (2, 2, 0), "the <DEL> and the stated het"
    assert {k: counts(KARYOTYPE_RECORDS, "raw", resolved_ref=False)[k] for k in _ALT_COUNTS} == {k: twin[k] for k in _ALT_COUNTS}


def test_sex_chromosome_measures_read_a_symbolic_alt_of_an_unstated_row_by_its_genotype_index(tmp_path: Path, store: ReferenceStore) -> None:
    """On a row `N` → `<DEL>,X` only the sequence ALT was spelled against the `N`: a GT index pointing at
    `<DEL>` is a symbolic call, read as any is (SCHEMA Decision 2), and one pointing at `X` is not read.
    So `0/1` counts as a diploid ALT call, a het and a depth site, as `N` → `<DEL>` with `0/1` does,
    while `0/2` and `1/2` (a slot on `X`) count as no ALT call, as a row with an unread allele. The whole
    row was dropped from these measures before, `0/1` included."""
    records = [
        (300, "del_het", "N", f"<DEL>,{other(base(300))}", "0/1:30"),
        (400, "sequence_het", "N", f"<DEL>,{other(base(400))}", "0/2:30"),
        (500, "del_and_sequence", "N", f"<DEL>,{other(base(500))}", "1/2:30"),
    ]
    got = measure(resolve_ref(load(vcf_with(tmp_path, "mixed", "GT:DP", records)), store).calls, None)["groups"]["autosomes"]
    only_del = measure(resolve_ref(load(vcf_with(tmp_path, "del", "GT:DP", [(300, "del_het", "N", "<DEL>", "0/1:30")])), store).calls, None)
    assert {k: got[k] for k in _ALT_COUNTS} == {k: only_del["groups"]["autosomes"][k] for k in _ALT_COUNTS}
    assert (got["diploid_alt"], got["het"], got["depth_sites"]) == (1, 1, 1)
    assert {k: got[k] for k in _PLOIDY_COUNTS} == {"rows": 3, "one": 0, "two": 3, "more": 0}


def test_array_letter_measures_read_no_pair_spelled_against_an_unstated_ref(tmp_path: Path) -> None:
    """An array's X is read from its letters (`base_pairs`, `letter_hets`): a chip VCF row `N` → `C,G`
    with GT `1/2` has letters `C` and `G`, both spelled against the `N`, so it is neither a base pair nor
    a letter het; its stated twin is both. Its call as written still counts (`called`, `no_calls`), as
    its ploidy does: the chip said it called the probe, whatever the letters name. Read on an autosome
    and on X outside the PARs."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n" + "\n".join(SEQ["1"][i : i + 60] for i in range(0, len(SEQ["1"]), 60)) + "\n" for name in ("1", "X")))
    two = ReferenceStore.open(tmp_path / "two")
    two.add_fasta(fasta)

    def letters(ref_of: object, name: str) -> dict[str, dict]:
        lines = [
            "##fileformat=VCFv4.2",
            *(f"##contig=<ID={c},length={len(SEQ['1'])}>" for c in ("1", "X")),
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]),
        ]
        for chrom in ("1", "X"):
            for pos in (300, 400):
                b = base(pos)
                alts = ",".join(x for x in "ACGT" if x != b)[:3]
                lines.append(f"{chrom}\t{pos}\t.\t{ref_of(b)}\t{alts}\t50\tPASS\t.\tGT\t1/2")
            lines.append(f"{chrom}\t500\t.\t{ref_of(base(500))}\t{other(base(500))}\t50\tPASS\t.\tGT\t./.")
        path = tmp_path / f"{name}.vcf"
        path.write_text("\n".join(lines) + "\n")
        return measure(resolve_ref(load(path), two).calls, "GRCh37")["groups"]

    unstated, stated = letters(lambda b: "N", "unstated"), letters(lambda b: b, "stated")
    for group in ("autosomes", "x_nonpar"):
        assert (stated[group]["base_pairs"], stated[group]["letter_hets"]) == (2, 2), group
        assert (unstated[group]["base_pairs"], unstated[group]["letter_hets"]) == (0, 0), group
        assert {k: unstated[group][k] for k in ("rows", "called", "no_calls", "called_one")} == {
            k: stated[group][k] for k in ("rows", "called", "no_calls", "called_one")
        }, group


def test_variant_classes_are_counted_at_load_before_any_ref_is_read(tmp_path: Path, store: ReferenceStore) -> None:
    """`variant_classes` is counted at load, where every row's `ref_check` is `not_run`: an unstated
    row's sequence ALT is not known to be one yet, so `N` → `X` counts as the SNV its lengths make it,
    and `resolve_ref` keeps the facet's count (Decision 2: the rule cannot apply to a reader that runs
    before the REF is read). Pinned, so a recount is a decision and not an accident."""
    records = [(300, "unstated_snv", "N", other(base(300)), "0/1"), (500, "unstated_del", "N", "<DEL>", "0/1")]
    loaded = load(vcf(tmp_path, records, "classes.vcf"))
    at_load = loaded.facets.collect()["variant_classes"].item().to_list()
    assert at_load == [{"variant_class": "snv", "alleles": 1}, {"variant_class": "sv", "alleles": 1}]
    resolved_facet = resolve_ref(loaded, store)
    assert resolved_facet.calls.collect()["ref_check"].to_list() == ["unstated", "unstated"]
    assert resolved_facet.facets.collect()["variant_classes"].item().to_list() == at_load


SITE_RECORDS = [
    (300, "unstated_hom_alt", "N", other(base(300)), "1/1"),
    (400, "unstated_het", "N", other(base(400)), "0/1"),
    (500, "unstated_hom_ref", "N", other(base(500)), "0/0"),
    (600, "stated_hom_alt", base(600), other(base(600)), "1/1"),
    (700, "unstated_del", "N", "<DEL>", "0/1"),
    (800, "beside_stated_unstated", "N", other(base(800)), "1/1"),
    (800, "beside_stated_stated", base(800), other(base(800)), "0/1"),
]
SITE_POSITIONS = [300, 400, 500, 600, 700, 800, 900]
"""Panel sites: each record's position, and 900, where the facet wrote nothing."""


def site_panel() -> pl.DataFrame:
    """A panel over `SITE_POSITIONS`: the reference's base, and `other` of it as the effect allele."""
    positions = SITE_POSITIONS
    return panel_sites(
        pl.DataFrame({"chrom": ["1"] * len(positions), "pos": positions, "ref_allele": [base(p) for p in positions], "effect_allele": [other(base(p)) for p in positions]})
    )


def read_sites(loaded: Loaded) -> dict[int, tuple[int | None, str]]:
    """`(count, basis)` per panel position, as relatedness and ancestry read the facet."""
    got = read_site_genotypes(loaded, site_facet(loaded), site_panel())
    return {SITE_POSITIONS[r["site"]]: (r["count"], r["basis"]) for r in got.iter_rows(named=True)}


def test_site_genotypes_read_no_genotype_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """Relatedness and ancestry read a facet at panel sites from its rows' literal alleles. `N` → `X`
    with GT `1/1` names `X` twice and no base beside it: read as written, it is two copies of the
    panel's effect allele `X`, the same as a stated `A` → `X`, at a site whose base the source never
    named. Its `0/1` and `0/0` were unknown before only because the letter `N` is neither panel allele.
    All three are `ref_unstated` now: unknown, and still a row at the position, so a whole-genome
    facet's silence is not read as reference there. A stated row beside one at the same position is
    the one usable row. A symbolic ALT reads as before (here `<DEL>` is no panel allele: unusable);
    before `resolve_ref` the rule cannot apply (`not_run`)."""
    path = vcf(tmp_path, SITE_RECORDS, "sites.vcf")
    loaded = with_facet(resolve_ref(load(path), store), absence_means="weak_reference")
    assert read_sites(loaded) == {
        300: (None, "ref_unstated"),
        400: (None, "ref_unstated"),
        500: (None, "ref_unstated"),
        600: (2, "observed"),
        700: (None, "unusable"),
        800: (1, "observed"),
        900: (0, "absent_as_reference"),
    }
    assert SiteBasis("ref_unstated") is SiteBasis.REF_UNSTATED
    # The twin writes the reference's base where the source wrote N: the same genotypes are read.
    twin = [(p, i, base(p), alt, g) for p, i, ref, alt, g in SITE_RECORDS if p in (300, 400, 500)]
    stated = with_facet(resolve_ref(load(vcf(tmp_path, twin, "twin.vcf")), store), absence_means="weak_reference")
    assert {p: read_sites(stated)[p] for p in (300, 400, 500)} == {300: (2, "observed"), 400: (1, "observed"), 500: (0, "observed")}
    # Before `resolve_ref` no REF was compared: the row is read as written.
    assert read_sites(with_facet(load(path), absence_means="weak_reference"))[300] == (2, "observed")


PHASED_RECORDS: list[Record] = [
    # PS 100: an unstated SNV, heterozygous, among stated rows.
    ("1", 700, base(700), other(base(700)), "0|1", 100),
    ("1", 710, base(710), other(base(710)), "1|0", 100),
    ("1", 720, base(720), other(base(720)), "1|1", 100),
    ("1", 730, "N", other(base(730)), "0|1", 100),
    # PS 200: an unstated SNV called hom-ref.
    ("1", 800, base(800), other(base(800)), "1|0", 200),
    ("1", 810, base(810), other(base(810)), "0|1", 200),
    ("1", 820, base(820), other(base(820)), "1|1", 200),
    ("1", 830, "N", other(base(830)), "0|0", 200),
    # PS 300, the control: an unstated row whose only ALT is symbolic.
    ("1", 900, base(900), other(base(900)), "1|1", 300),
    ("1", 910, base(910), other(base(910)), "1|1", 300),
    ("1", 920, "N", "<DEL>", "0|1", 300),
]
"""Three PS sets on contig 1, each with one row whose source REF is `N`."""


def phased_blocks(tmp_path: Path, store: ReferenceStore, records: list[Record], name: str, reference_alleles: bool) -> dict[tuple[int, int], dict]:
    loaded = mint_vrs_ids(resolve_ref(load(phased_vcf(tmp_path, records, {"1": len(SEQ["1"])}, name=name)), store), store)
    frame = cis_phased_blocks(loaded, reference_alleles=reference_alleles).collect()
    return {(int(r["phase_set"]), r["haplotype"]): r for r in frame.iter_rows(named=True)}


def withheld_by(row: dict) -> tuple:
    return (row["withheld"], row["withheld_chrom"], row["withheld_pos"], row["withheld_allele"])


@pytest.mark.parametrize("reference_alleles", [False, True])
def test_phase_blocks_take_no_haplotype_of_a_row_spelled_against_an_unstated_ref(tmp_path: Path, store: ReferenceStore, reference_alleles: bool) -> None:
    """A genotype on `N` → `X` indexes `[N, X]`: its slot 0 names a base the source never stated, as
    much as its slot 1 does, which is why `genotype_evidence` grades none of it. A phase block follows:
    no haplotype of such a row joins one, REF slot included, and a block it would have joined is
    withheld `ref_unstated`, named by that slot (its allele as `gt_alleles` holds it, the source's
    `N` on a REF slot). Read as written, `reference_alleles=True` took the REF slot's minted id as a
    member, a block stating the reference's base on that haplotype. By default a REF slot joins no
    block anyway, so nothing changes there. An unstated row with only a symbolic ALT is read as before:
    its REF slot joins with `reference_alleles=True`, its `<DEL>` withholds its own reason. The twin
    writes the reference's base where the source wrote `N`: its blocks are what the stated rows give."""
    got = phased_blocks(tmp_path, store, PHASED_RECORDS, "unstated.vcf", reference_alleles)
    twin = phased_blocks(tmp_path, store, [(c, p, base(p) if r == "N" else r, a, g, s) for c, p, r, a, g, s in PHASED_RECORDS], "twin.vcf", reference_alleles)
    assert set(got) == set(twin) == {(s, h) for s in (100, 200, 300) for h in (0, 1)}
    assert all(row["cpb_id"] is not None for key, row in twin.items() if key != (300, 1))
    assert all(got[key]["sites"] == twin[key]["sites"] for key in got)
    # The ALT slot of the heterozygous row withholds its block either way (`mint_vrs_ids` withholds the ALT).
    assert withheld_by(got[(100, 1)]) == ("ref_unstated", "1", 730, other(base(730)))
    # The control: `<DEL>` has no VRS model, and the REF slot reads as on the twin.
    assert withheld_by(got[(300, 1)]) == withheld_by(twin[(300, 1)]) == ("sv_no_vrs_model", "1", 920, "<DEL>")
    assert got[(300, 0)]["cpb_id"] == twin[(300, 0)]["cpb_id"] is not None
    if reference_alleles:
        assert withheld_by(got[(100, 0)]) == ("ref_unstated", "1", 730, "N")
        assert withheld_by(got[(200, 0)]) == withheld_by(got[(200, 1)]) == ("ref_unstated", "1", 830, "N")
        assert all(got[key]["members"] is None and got[key]["cpb_id"] is None for key in [(100, 0), (200, 0), (200, 1)])
        assert len(got[(300, 0)]["members"]) == 3  # 900's and 910's ALT, 920's REF slot
    else:
        # A REF slot joins no block by default: the hom-ref row and the het row's slot 0 change nothing.
        stated = [(100, 0), (200, 0), (200, 1)]
        assert {key: got[key]["cpb_id"] for key in stated} == {key: twin[key]["cpb_id"] for key in stated}
        assert len(got[(300, 0)]["members"]) == 2
