"""VRS allele ids (RM271): `mint_vrs_ids` against vrs-python, the oracle, on the same refget store.

Every expected id is computed at runtime by vrs-python (`scripts/vrs_oracle.py`) from the row's own
position, REF and allele, so the tests compare two implementations of the VRS 2.x normalization and
digest over identical sequence, never an id copied from a run.
"""

import random
from pathlib import Path

import polars as pl
import pytest
from vrs_oracle import StoreProxy, oracle_id

from just_vcf import ReferenceMismatchError, ReferenceStore, load, mint_vrs_ids, resolve_ref, vrs_withheld
from just_vcf.loading import Loaded
from just_vcf.schema import CALLS
from just_vcf.vrs import WITHHELD_REASONS, _ContigSequence, _justify, _Tally


def _contig_one(rng: random.Random) -> str:
    """3000 bases with repeats where justification has to roll, and a run of N."""
    def rnd(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    seq = rnd(500) + "CACACACACA" + rnd(290) + "AAAAAAAA" + rnd(192) + "ACGACGACG" + rnd(191) + "N" * 10 + rnd(290)
    seq += "GTTGTTGTTGTT" + rnd(3000 - len(seq) - 12)
    return seq


def sequences(seed: int = 11) -> dict[str, str]:
    rng = random.Random(seed)
    return {"1": _contig_one(rng), "2": "".join(rng.choice("ACGT") for _ in range(2000)), "MT": "".join(rng.choice("ACGT") for _ in range(400))}


SEQ = sequences()
LENGTHS = {name: len(seq) for name, seq in SEQ.items()}
CA_REPEAT, POLY_A, TRIPLET, N_RUN = 501, 801, 1001, 1201
"""1-based starts of the features `_contig_one` lays out."""


def base(chrom: str, pos: int, length: int = 1) -> str:
    return SEQ[chrom][pos - 1 : pos - 1 + length]


def other(b: str) -> str:
    return next(x for x in "ACGT" if x != b)


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    lines = []
    for name, seq in seqs.items():
        lines.append(f">{name}")
        lines += [seq[i : i + 60] for i in range(0, len(seq), 60)]
    path.write_text("\n".join(lines) + "\n")
    return path


def vcf(tmp_path: Path, records: list[tuple[str, int, str, str]], contigs: dict[str, int], name: str = "calls.vcf") -> Path:
    """A one-sample VCF of `(chrom, pos, ref, alt)` records, each genotyped `0/1`."""
    meta = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={n}>" for c, n in contigs.items())]
    meta += [
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1" for c, p, r, a in records]
    path = tmp_path / name
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa", SEQ))
    return s


def minted(path: Path, store: ReferenceStore) -> tuple[Loaded, pl.DataFrame]:
    result = mint_vrs_ids(resolve_ref(load(path), store), store)
    return result, result.calls.collect()


def accessions(store: ReferenceStore) -> dict[str, str]:
    (digest,) = store.collections()
    return {contig.name: f"SQ.{contig.digest}" for contig in store.contigs(digest)}


def expected_ids(calls: pl.DataFrame, store: ReferenceStore, reference_name: dict[str, str] | None = None) -> list[list[str]]:
    """vrs-python's id for every slot of `[ref] + alts`, from each row's own `pos`, `ref` and allele."""
    proxy = StoreProxy(store)
    sq = accessions(store)
    names = reference_name or {}
    return [
        [oracle_id(proxy, sq[names.get(row["chrom"], row["chrom"])], row["pos"], row["ref"], allele) for allele in [row["ref"], *row["alts"]]]
        for row in calls.iter_rows(named=True)
    ]


def test_ids_equal_vrs_python_for_every_kind_of_allele(tmp_path: Path, store: ReferenceStore) -> None:
    ca, poly, trip = CA_REPEAT, POLY_A, TRIPLET
    mnv_ref = base("1", 60, 3)
    records = [
        ("1", 10, base("1", 10), other(base("1", 10))),  # SNV
        ("1", 20, base("1", 20, 2), other(base("1", 20)) + other(base("1", 21))),  # MNV
        ("1", 60, mnv_ref, mnv_ref[0] + other(mnv_ref[1]) + mnv_ref[2]),  # MNV whose flanks trim to an SNV
        ("1", 100, base("1", 100), base("1", 100) + "TTGCA"),  # insertion the reference does not repeat
        ("1", 150, base("1", 150, 4), base("1", 150)),  # deletion
        ("1", 200, base("1", 200, 3), other(base("1", 200))),  # deletion-insertion
        ("1", ca + 5, base("1", ca + 5), base("1", ca + 5) + "AC"),  # insertion inside the CA repeat
        ("1", ca + 7, base("1", ca + 7, 3), base("1", ca + 7)),  # CA deletion written at the repeat's right end
        ("1", poly + 6, base("1", poly + 6, 2), base("1", poly + 6)),  # A deletion that has to roll left
        ("1", poly + 2, base("1", poly + 2), base("1", poly + 2) + "AAG"),  # insertion in a repeat, not a cycle of it
        ("1", trip + 2, base("1", trip + 2), base("1", trip + 2) + "ACGACG"),  # two units into a triplet repeat
        ("1", N_RUN - 1, base("1", N_RUN - 1, 2), base("1", N_RUN - 1)),  # deletion of the first N
        ("1", 1400, base("1", 1400), ",".join(sorted({other(base("1", 1400)), base("1", 1400) + "G"}))),  # multi-allelic
        ("1", 1450, base("1", 1450), other(base("1", 1450)).lower()),  # lowercase ALT
        ("1", 1480, base("1", 1480), "."),  # monomorphic: only the REF slot
        ("chrM", 30, base("MT", 30, 2), base("MT", 30)),  # chrM in the VCF, MT in the reference
    ]
    lengths = {"1": LENGTHS["1"], "chrM": LENGTHS["MT"]}
    result, calls = minted(vcf(tmp_path, records, lengths), store)

    assert calls["ref_check"].cast(pl.String).unique().to_list() == ["match"]
    assert calls["vrs_allele_ids"].to_list() == expected_ids(calls, store, {"MT": "MT"})
    assert calls["vrs_errors"].null_count() == calls.height
    assert calls["vrs_allele_ids"].list.len().to_list() == (1 + calls["alts"].list.len()).to_list()
    assert calls.schema == CALLS
    slots = int((1 + calls["alts"].list.len()).sum())
    assert any(f"for {slots} of {slots} allele(s): {calls.height} REF, {slots - calls.height} ALT" in w for w in result.warnings)


def test_withheld_alleles_are_null_with_a_reason(tmp_path: Path, store: ReferenceStore) -> None:
    snv_alt = other(base("1", 300))
    records = [
        ("1", 300, base("1", 300), f"{snv_alt},*"),  # the SNV mints, the spanning deletion does not
        ("1", 310, base("1", 310), "<DEL>"),
        ("1", 320, base("1", 320), f"{base('1', 320)}[2:100["),  # breakend
        ("1", 330, base("1", 330), "<NON_REF>"),  # read as <*>
        ("1", 340, base("1", 340), "R"),  # an ambiguity code is not a sequence the id can carry
        ("1", 350, other(base("1", 350)), base("1", 350)),  # REF the reference contradicts
        ("1", 2999, base("1", 2999) + "AAAA", base("1", 2999)),  # REF reaching past the contig's end
        ("3", 10, "A", "G"),  # a contig the reference does not have
    ]
    result, calls = minted(vcf(tmp_path, records, {"1": LENGTHS["1"], "3": 500}), store)
    reasons = vrs_withheld(calls.lazy()).collect()
    by_pos = {(row["chrom"], row["pos"]): row for row in reasons.iter_rows(named=True)}

    def withheld(chrom: str, pos: int) -> list[str | None]:
        return by_pos[(chrom, pos)]["vrs_withheld"]

    proxy, sq = StoreProxy(store), accessions(store)
    row = by_pos[("1", 300)]
    assert row["vrs_allele_ids"] == [oracle_id(proxy, sq["1"], 300, row["ref"], a) for a in (row["ref"], snv_alt)] + [None]
    assert withheld("1", 300) == [None, None, "spanning_deletion"]
    assert withheld("1", 310) == [None, "sv_no_vrs_model"]
    assert withheld("1", 320) == [None, "sv_no_vrs_model"]
    assert withheld("1", 330) == [None, "unspecified"]
    assert withheld("1", 340) == [None, "not_sequence"]
    assert withheld("1", 350) == ["ref_mismatch", "ref_mismatch"]
    assert by_pos[("1", 350)]["vrs_allele_ids"] is None  # nothing minted: the whole list is null
    assert withheld("1", 2999) == ["ref_unchecked", "ref_unchecked"]
    assert withheld("3", 10) == ["ref_unchecked", "ref_unchecked"]
    assert calls["vrs_errors"].null_count() == calls.height  # minting refused nothing; every reason is the frame's

    # Every slot has an id or a reason, never both, and every reason is in the vocabulary.
    for row in reasons.iter_rows(named=True):
        ids = row["vrs_allele_ids"] or [None] * len(row["vrs_withheld"])
        assert [(i is None) != (r is None) for i, r in zip(ids, row["vrs_withheld"])] == [True] * len(ids)
        assert {r for r in row["vrs_withheld"] if r is not None} <= WITHHELD_REASONS
    warnings = "\n".join(result.warnings)
    assert "4 allele(s) without a VRS id: REF not read from the reference" in warnings
    assert "2 allele(s) without a VRS id: the source REF disagrees with the reference (ref_mismatch)" in warnings
    assert "2 allele(s) without a VRS id: a structural variant (a symbolic allele or a breakend)" in warnings


def test_before_minting_the_reason_is_the_pipeline_state(tmp_path: Path, store: ReferenceStore) -> None:
    records = [("1", 10, base("1", 10), other(base("1", 10))), ("1", 20, base("1", 20), "*")]
    loaded = load(vcf(tmp_path, records, {"1": LENGTHS["1"]}))
    unread = vrs_withheld(loaded.calls).collect()["vrs_withheld"].to_list()
    assert unread == [["ref_unchecked", "ref_unchecked"], ["ref_unchecked", "spanning_deletion"]]
    resolved = vrs_withheld(resolve_ref(loaded, store).calls).collect()["vrs_withheld"].to_list()
    assert resolved == [["not_run", "not_run"], ["not_run", "spanning_deletion"]]


def test_random_alleles_equal_vrs_python(tmp_path: Path, store: ReferenceStore) -> None:
    """Hundreds of alleles of every kind at random places, repeats and the N run included."""
    rng = random.Random(3)
    records = []
    for pos in sorted(rng.sample(range(2, LENGTHS["1"] - 10), 600)):
        kind = rng.choice(["snv", "mnv", "ins", "rep_ins", "del", "delins"])
        if kind == "snv":
            ref = base("1", pos)
            alt = other(ref) if ref != "N" else "A"
        elif kind == "mnv":
            ref = base("1", pos, 3)
            alt = "".join(rng.choice("ACGT") for _ in range(3))
        elif kind == "ins":
            ref = base("1", pos)
            alt = ref + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
        elif kind == "rep_ins":
            ref = base("1", pos)
            alt = ref + base("1", pos + 1, rng.randint(1, 6))
        elif kind == "del":
            ref = base("1", pos, rng.randint(2, 7))
            alt = ref[0]
        else:
            ref = base("1", pos, rng.randint(2, 4))
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
        if ref != alt:
            records.append(("1", pos, ref, alt))
    _, calls = minted(vcf(tmp_path, records, {"1": LENGTHS["1"]}), store)
    assert calls.height == len(records)
    assert calls["vrs_allele_ids"].to_list() == expected_ids(calls, store)


def test_a_deletion_without_its_anchor_base_mints_the_anchored_id(tmp_path: Path, store: ReferenceStore) -> None:
    """The empty allele is VRS's zero-length literal: `CT>""` at 31 (an upload's ALT `-`, F12) is the
    deletion `GCT>G` at 30, and gets its id, the one vrs-python gives either spelling. One inside the
    poly-A run rolls over the run like any deletion."""
    poly = POLY_A
    records = [
        ("1", 31, base("1", 31, 2), "-"),
        ("1", 400, base("1", 400, 2), base("1", 400)),  # the same kind of event, anchored, elsewhere
        ("1", poly + 4, base("1", poly + 4), "-"),
    ]
    anchored = [("1", 30, base("1", 30, 3), base("1", 30)), ("1", poly + 3, base("1", poly + 3, 2), base("1", poly + 3))]
    _, calls = minted(vcf(tmp_path, records, {"1": LENGTHS["1"]}), store)
    _, twins = minted(vcf(tmp_path, anchored, {"1": LENGTHS["1"]}, name="anchored.vcf"), store)
    rows = {row["pos"]: row for row in calls.iter_rows(named=True)}
    assert rows[31]["alts"] == [""] and rows[31]["ref_check"] == "match"
    assert calls["vrs_allele_ids"].to_list() == expected_ids(calls, store)
    by_pos = {row["pos"]: row["vrs_allele_ids"][1] for row in twins.iter_rows(named=True)}
    assert rows[31]["vrs_allele_ids"][1] == by_pos[30]
    assert rows[poly + 4]["vrs_allele_ids"][1] == by_pos[poly + 3]
    assert vrs_withheld(calls.lazy()).collect()["vrs_withheld"].to_list() == [[None, None]] * calls.height


PREAMBLE_23 = (
    "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n"
    "#\n"
    "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
    "#\n"
    "# rsid\tchromosome\tposition\tgenotype\n"
)


def test_array_rows_mint_once_their_alleles_are_derived(tmp_path: Path) -> None:
    labelled = ReferenceStore.open(tmp_path / "grch37")
    labelled.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    r10, r20 = base("1", 10), base("1", 20)
    rows = [
        ("rs1", "1", 10, r10 * 2),  # hom-ref: alts are [], so only the REF slot
        ("rs2", "1", 20, other(r20) + r20),  # het
        ("i3", "1", 30, "DI"),  # indel codes: nothing to place
        ("rs4", "1", 40, "--"),  # no-call: REF span unknown
        ("rs5", "25", 6, "AA"),  # PAR code, positioned on X, which this reference lacks
    ]
    path = tmp_path / "genome.txt"
    path.write_text(PREAMBLE_23 + "".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows))
    result = mint_vrs_ids(resolve_ref(load(path), labelled), labelled)
    calls = vrs_withheld(result.calls).collect()
    by_id = {row["ids"][0]: row for row in calls.iter_rows(named=True)}

    proxy, sq = StoreProxy(labelled), accessions(labelled)
    assert by_id["rs1"]["vrs_allele_ids"] == [oracle_id(proxy, sq["1"], 10, r10, r10)]
    assert by_id["rs2"]["vrs_allele_ids"] == [oracle_id(proxy, sq["1"], 20, r20, a) for a in (r20, other(r20))]
    assert by_id["i3"]["vrs_allele_ids"] is None and by_id["i3"]["vrs_withheld"] == ["array_code"]
    assert by_id["rs4"]["vrs_withheld"] == ["ref_unresolved"]
    assert by_id["rs5"]["vrs_withheld"] == ["ref_unresolved"]


def test_minting_is_deterministic_and_order_independent(tmp_path: Path, store: ReferenceStore) -> None:
    records = [
        ("1", 10, base("1", 10), other(base("1", 10))),
        ("1", POLY_A + 6, base("1", POLY_A + 6, 2), base("1", POLY_A + 6)),
        ("1", CA_REPEAT + 5, base("1", CA_REPEAT + 5), base("1", CA_REPEAT + 5) + "AC"),
        ("2", 50, base("2", 50, 3), base("2", 50)),
    ]
    _, first = minted(vcf(tmp_path, records, {"1": LENGTHS["1"], "2": LENGTHS["2"]}), store)
    _, again = minted(vcf(tmp_path, records, {"1": LENGTHS["1"], "2": LENGTHS["2"]}), store)
    shuffled = records[::-1]
    _, reordered = minted(vcf(tmp_path, shuffled, {"1": LENGTHS["1"], "2": LENGTHS["2"]}, name="shuffled.vcf"), store)
    assert first.equals(again)
    assert first.drop("facet_id").equals(reordered.drop("facet_id"))


def test_an_allele_gtars_refuses_is_withheld_and_counted_by_kind(tmp_path: Path, store: ReferenceStore) -> None:
    """After `resolve_ref` gtars cannot see a REF that contradicts the sequence; driven directly here."""
    (digest,) = store.collections()
    contig = store.contigs_by_chrom(digest)["1"][0]
    tally = _Tally()
    sequence = _ContigSequence(store, contig)
    wrong = [(start, other(base("1", start + 1)) + base("1", start + 2), base("1", start + 2)) for start in (100, 200)]
    results = [_justify(start, ref, alt, sequence, tally) for start, ref, alt in wrong]
    assert [result[-1] for result in results] == ["refused", "refused"]
    assert all(value is None for result in results for value in result[:-1])
    assert tally.warnings(digest)[-1] == "2 allele(s) refused by gtars: ref allele mismatch at interbase #"


def test_a_facet_whose_ref_was_not_read_is_refused(tmp_path: Path, store: ReferenceStore) -> None:
    loaded = load(vcf(tmp_path, [("1", 10, base("1", 10), "A")], {"1": LENGTHS["1"]}))
    with pytest.raises(ValueError, match="run resolve_ref first"):
        mint_vrs_ids(loaded, store)


def test_a_store_without_the_facets_collection_is_refused(tmp_path: Path, store: ReferenceStore) -> None:
    resolved = resolve_ref(load(vcf(tmp_path, [("1", 10, base("1", 10), "A")], {"1": LENGTHS["1"]})), store)
    other_store = ReferenceStore.open(tmp_path / "other")
    other_store.add_fasta(write_fasta(tmp_path / "other.fa", sequences(seed=99)))
    with pytest.raises(ReferenceMismatchError, match="does not hold collection"):
        mint_vrs_ids(resolved, other_store)
