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

from just_vcf import ReferenceStore, liftover, load, mint_vrs_ids, normalize, query, resolve_ref, vrs_withheld, write_vcf
from just_vcf.loading import Loaded
from just_vcf.schema import CALLS
from just_vcf.vcf_zarr import write_vcf_zarr
from just_vcf.vrs import WITHHELD_REASONS
from test_cli import report, run
from test_liftover import TARGET, by_case, expected_sv, expected_target, s1, setup, sv  # noqa: F401 (setup is a fixture)
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
    # A real base beside an N states that base: compared as before, and the reference has no N there.
    "partial_after": "mismatch",
    "partial_before": "mismatch",
    # Over a gap the reference reads N too: the letters agree, so the comparison ran and matched.
    "gap": "match",
}


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


def test_write_vcf_writes_an_unstated_row_with_the_reference_base(tmp_path: Path, store: ReferenceStore) -> None:
    source = resolved(tmp_path, store)
    mismatched = sum(1 for v in EXPECTED.values() if v == "mismatch")
    written = write_vcf(source, tmp_path / "out.vcf", compress=False, resolved_only=True)
    assert written.rows_skipped == {"ref_mismatch": mismatched}
    back = by_id(load(written.path).calls.collect())
    assert set(back) == {case for case, value in EXPECTED.items() if value != "mismatch"}
    for pos, case, ref, alt, _ in RECORDS:
        if case in back:
            assert back[case]["source_ref"] == base(pos, len(ref)), case
            assert back[case]["alts"] == ([] if alt == "." else [alt]), case


def test_write_vcf_zarr_keeps_the_unstated_rows(tmp_path: Path, store: ReferenceStore) -> None:
    source = resolved(tmp_path, store)
    mismatched = sum(1 for v in EXPECTED.values() if v == "mismatch")
    written = write_vcf_zarr(source, tmp_path / "out.vcz", resolved_only=True)
    assert written.rows_skipped == {"ref_mismatch": mismatched}
    sites = read_back(tmp_path / "out.vcz")
    kept = sorted(pos for pos, case, *_ in RECORDS if EXPECTED[case] != "mismatch")
    assert sites["start"].to_list() == kept
    refs = dict(zip(sites["start"].to_list(), sites["ref"].to_list()))
    assert all(refs[pos] == base(pos, len(ref)) for pos, case, ref, *_ in RECORDS if pos in refs)


def test_query_counts_unstated_and_reads_no_site_from_an_unstated_ref(tmp_path: Path, store: ReferenceStore) -> None:
    """An unstated row states its ALT as written, so it answers for that allele; it never says which
    base the site holds, so an SNV at its POS that it does not spell is not answered `site_genotyped`."""
    q = query(resolved(tmp_path, store))
    (table,) = q.counts(by=["ref_check"])
    assert {group.values[0]: group.count for group in table.groups} == {
        value: sum(1 for v in EXPECTED.values() if v == value) for value in set(EXPECTED.values())
    }
    (stated,) = q.carries(chrom="1", pos=120, ref=base(120), alt=other(base(120)))
    assert (stated.basis.value, stated.carried) == ("stated", True)
    elsewhere = next(b for b in "ACGT" if b not in (base(120), other(base(120))))
    (unstated,) = q.carries(chrom="1", pos=120, ref=base(120), alt=elsewhere)
    assert (unstated.basis.value, unstated.carried) == ("overlapping_call", None)
    (matched,) = q.carries(chrom="1", pos=10, ref=base(10), alt=next(b for b in "ACGT" if b not in (base(10), other(base(10)))))
    assert (matched.basis.value, matched.carried) == ("site_genotyped", False), "a stated REF still genotypes its site"


def test_the_cli_counts_unstated_rows(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "ref.fa")
    [facet] = report(run("check", vcf(tmp_path, RECORDS), "--reference", fasta, "--json"))["facets"]
    assert facet["ref_check"] == {value: sum(1 for v in EXPECTED.values() if v == value) for value in set(EXPECTED.values())}


# Lifted through test_liftover's construction (block A: source 101..400 at -50; 701..800 reversed onto
# target 2; SWAP_POS 251 is a base the target changes).
LIFT_CASES = {
    "cnv_del_plus": (sv("cnv_del_plus", 201, "<DEL>", 230, ref="N"), None),
    "cnv_dup_minus": (sv("cnv_dup_minus", 720, "<DUP>", 760, ref="N"), None),
    "cnv_segment": (sv("cnv_segment", 460, ".", 520, sample="0/0:2", ref="N"), None),
    # The lift's own check still runs: on a plus strand the target's base must be the base read at the source.
    "cnv_anchor_changed": (sv("cnv_anchor_changed", 251, "<DEL>", 280, ref="N"), "ref_mismatch_after_lift"),
    "cnv_partial": (sv("cnv_partial", 340, "<DEL>", 360, ref=s1(340) + "N"), "source_ref_mismatch"),
    "unstated_snv": (f"1\t205\tunstated_snv\tN\t{lift_other(s1(205))}\t50\tPASS\t.\tGT\t0/1", None),
}


def test_liftover_lifts_unstated_rows_like_matching_ones(tmp_path: Path, setup) -> None:  # noqa: F811
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
    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == {case for case, (_, reason) in LIFT_CASES.items() if reason is None}
    for case, r in calls.items():
        p = parent[case]
        if case == "unstated_snv":
            chrom, pos, _ = expected_target(p["pos"])
        else:
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
