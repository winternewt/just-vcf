"""VRS allele ids (RM271): `mint_vrs_ids` against vrs-python, the oracle, on the same refget store.

Every expected id is computed at runtime by vrs-python (`scripts/vrs_oracle.py`) from the row's own
position, REF and allele, so the tests compare two implementations of the VRS 2.x normalization and
digest over identical sequence, never an id copied from a run.
"""

import logging
import random
import re
from collections import Counter
from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("ga4gh.vrs", reason="vrs-python, the oracle, is a dev dependency")

from vrs_oracle import StoreProxy, oracle_id  # noqa: E402  (after the skip: it imports vrs-python)

from just_vcf import ReferenceMismatchError, ReferenceStore, load, mint_vrs_ids, resolve_ref, sniff, vrs_withheld  # noqa: E402
from just_vcf.assembly import ASSEMBLY_LENGTHS  # noqa: E402
from just_vcf.loading import Loaded  # noqa: E402
from just_vcf.schema import BLOCKS, CALLS  # noqa: E402
from just_vcf.vrs import WITHHELD_REASONS, _ContigSequence, _justify, _Tally  # noqa: E402

logger = logging.getLogger(__name__)


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


# Opt-in: the real refget store `JUST_VCF_TEST_REFERENCE_DIR` names (GRCh37 and GRCh38), the same
# comparison on real sequence. Inputs are synthetic alleles placed on that sequence, and slices of
# consented-public PGP-HMS uploads; aggregates only are logged, never a row.

REAL_ASSEMBLIES = ("GRCh37", "GRCh38")
WINDOW = (30_000_000, 30_400_000)
"""0-based half-open: a stretch of chromosome 22 with no N on either assembly, where the repeats the
synthetic alleles need are found at test time."""
ABSENT_CONTIG = {"GRCh37": "hs37d5", "GRCh38": "chrEBV"}
"""A real contig of each assembly's analysis set that the primary-assembly collection lacks
(1000 Genomes' decoy, the Epstein-Barr virus sequence), checked absent at test time."""

Real = tuple[str, int, str, str, "list[str | None] | None"]
"""`(chrom, pos, ref, alt text, expected reasons)`: reasons per slot of `[ref] + alts` as SCHEMA's
"Why an id is missing" gives them, or None where every slot must mint."""


def real_collection(store: ReferenceStore, label: str) -> str:
    digests = [digest for digest in store.collections() if store.assembly_label(digest) == label]
    if len(digests) != 1:
        pytest.skip(f"the store holds {len(digests)} {label} collection(s)")
    return digests[0]


def real_sequence(store: ReferenceStore, digest: str, chrom: str) -> str:
    (contig,) = store.contigs_by_chrom(digest)[chrom]
    return store.fetch(contig, [(0, contig.length)])[0]


def oracle_slots(proxy: StoreProxy, store: ReferenceStore, digest: str, row: dict, slots: list[int] | None = None) -> list[str]:
    """vrs-python's id for each of `slots` (every one by default) of the row's `[ref] + alts`, from its
    own pos, REF and allele."""
    (contig,) = store.contigs_by_chrom(digest)[row["chrom"]]
    alleles = [row["ref"], *(row["alts"] or [])]
    return [oracle_id(proxy, f"SQ.{contig.digest}", row["pos"], row["ref"], alleles[slot]) for slot in (slots if slots is not None else range(len(alleles)))]


def repeats(window: str, offset: int) -> list[tuple[int, str, int]]:
    """`(1-based start, unit, length)` of runs in `window` (which starts at 0-based `offset`):
    homopolymers of 10 or more, and di-, tri- and tetranucleotide repeats of five units or more."""
    found = []
    for size, floor in ((1, 10), (2, 6), (3, 5), (4, 5)):
        for match in re.finditer(rf"(?=(([ACGT]{{{size}}})\2{{{floor - 1},}}))", window):
            unit = match.group(2)
            if size > 1 and len(set(unit)) == 1:
                continue
            found.append((offset + match.start() + 1, unit, len(match.group(1))))
    # A match inside a longer run of the same unit is the same repeat; keep the run's first start.
    kept: list[tuple[int, str, int]] = []
    for start, unit, length in sorted(found):
        if kept and start < kept[-1][0] + kept[-1][2]:
            continue
        kept.append((start, unit, length))
    return kept


def real_records(label: str, seq: str, mt: str, rng: random.Random) -> list[Real]:
    """Every kind of allele on chromosome 22 and MT, the withheld cases, and a random batch."""
    prefix = "chr" if label == "GRCh38" else ""
    c22, cmt = f"{prefix}22", f"{prefix}M" if prefix else "MT"
    used: set[tuple[str, int]] = set()
    records: list[Real] = []

    def b(pos: int, n: int = 1, on: str = seq) -> str:
        return on[pos - 1 : pos - 1 + n]

    def add(chrom: str, pos: int, ref: str, alt: str, reasons: "list[str | None] | None" = None) -> None:
        """One record, unless it would overlap another (each row is then its own site)."""
        span = {(chrom, p) for p in range(pos - 2, pos + len(ref) + 2)}
        if span & used:
            return
        used.update(span)
        records.append((chrom, pos, ref, alt, reasons))

    # Fixed kinds at spaced places in the window, as the synthetic test lays them out.
    start = WINDOW[0] + 1
    p = start + 100
    mnv = b(p + 40, 3)
    add(c22, p, b(p), other(b(p)))  # SNV
    add(c22, p + 20, b(p + 20, 2), other(b(p + 20)) + other(b(p + 21)))  # MNV
    add(c22, p + 40, mnv, mnv[0] + other(mnv[1]) + mnv[2])  # MNV whose flanks trim to an SNV
    add(c22, p + 60, b(p + 60), b(p + 60) + "TTGCA")  # insertion
    add(c22, p + 80, b(p + 80, 5), b(p + 80))  # deletion
    add(c22, p + 100, b(p + 100, 3), other(b(p + 100)))  # deletion-insertion
    multi = b(p + 120, 3)
    alts = sorted({multi[0], other(multi[0]) + multi[1:], multi + "G", multi[0] + "AT" + multi[1:]})
    add(c22, p + 120, multi, ",".join(alts))  # multi-allelic: deletion, MNV, insertions
    add(c22, p + 140, b(p + 140), f"{other(b(p + 140))},{b(p + 140)}AC")  # multi-allelic: SNV and insertion
    add(c22, p + 160, b(p + 160), other(b(p + 160)).lower())  # lowercase ALT
    add(c22, p + 180, b(p + 180), ".")  # monomorphic: the REF slot only

    # Alleles inside the repeats the window holds: VOCA's expansion has to roll each one.
    runs = repeats(seq[WINDOW[0] : WINDOW[1]], WINDOW[0])
    by_size: dict[int, list[tuple[int, str, int]]] = {}
    for run in runs:
        by_size.setdefault(len(run[1]), []).append(run)
    for _, found in sorted(by_size.items()):
        for run_start, unit, length in found[:4]:
            k = len(unit)
            inside = run_start + k  # an anchor inside the run
            add(c22, inside, b(inside), b(inside) + b(inside + 1, k))  # one unit inserted: a cycle of the repeat
            two = run_start + length - 2  # two units inserted near the run's right end
            add(c22, two, b(two), b(two) + b(two - 2 * k + 1, 2 * k))
            right = run_start + length - k - 1  # one unit deleted, written at the run's right end
            add(c22, right, b(right, k + 1), b(right))
            if k > 1:  # part of a unit deleted: inside the run, but no cycle of it
                add(c22, run_start + length // 2, b(run_start + length // 2, 2), b(run_start + length // 2))
            else:  # inserted bases that continue the run, then leave it
                add(c22, run_start + 3, b(run_start + 3), b(run_start + 3) + unit * 2 + other(unit))
    # The first bases after chromosome 22's leading Ns.
    first = re.search(r"N[ACGT]", seq).start() + 2
    add(c22, first - 1, b(first - 1, 3), b(first - 1))  # a deletion anchored on the last N
    add(c22, first + 5, b(first + 5), other(b(first + 5)))  # an SNV a few bases in

    # MT, the rCRS: its D-loop C run, a deletion, a multi-allelic row, and REF reaching past its end.
    def m(pos: int, n: int = 1) -> str:
        return b(pos, n, mt)

    poly_c = re.search(r"C{7,}", mt).start() + 1
    add(cmt, 73, m(73), other(m(73)))
    add(cmt, poly_c + 2, m(poly_c + 2), m(poly_c + 2) + "C")
    add(cmt, poly_c + 30, m(poly_c + 30, 3), m(poly_c + 30))
    add(cmt, 750, m(750, 2), f"{m(750)},{other(m(750))}{m(751)}")
    add(cmt, len(mt) - 1, m(len(mt) - 1, 2) + "AC", m(len(mt) - 1), ["ref_unchecked", "ref_unchecked"])

    # The withheld cases, each with the reason SCHEMA's ordered rules give it.
    q = start + 1000
    snv = other(b(q))
    add(c22, q, b(q), f"{snv},*", [None, None, "spanning_deletion"])
    add(c22, q + 20, b(q + 20), "<DEL>", [None, "sv_no_vrs_model"])
    add(c22, q + 40, b(q + 40), f"{b(q + 40)}[{c22}:{q + 5000}[", [None, "sv_no_vrs_model"])
    add(c22, q + 60, b(q + 60), "<NON_REF>", [None, "unspecified"])
    add(c22, q + 80, b(q + 80), "R", [None, "not_sequence"])
    add(c22, q + 100, other(b(q + 100)), b(q + 100), ["ref_mismatch", "ref_mismatch"])
    add(ABSENT_CONTIG[label], 10, "A", "G", ["ref_unchecked", "ref_unchecked"])

    # A random batch over the window and inside its repeats.
    places = rng.sample(range(start + 2000, WINDOW[1] - 20), 500)
    places += [run_start + rng.randrange(0, length) for run_start, _, length in rng.sample(runs, min(150, len(runs)))]
    for pos in sorted(places):
        kind = rng.choice(["snv", "mnv", "ins", "rep_ins", "del", "delins"])
        if kind == "snv":
            ref, alt = b(pos), other(b(pos))
        elif kind == "mnv":
            ref, alt = b(pos, 3), "".join(rng.choice("ACGT") for _ in range(3))
        elif kind == "ins":
            ref = b(pos)
            alt = ref + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
        elif kind == "rep_ins":
            ref = b(pos)
            alt = ref + b(pos + 1, rng.randint(1, 8))
        elif kind == "del":
            ref = b(pos, rng.randint(2, 9))
            alt = ref[0]
        else:
            ref = b(pos, rng.randint(2, 4))
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
        if ref != alt:
            add(c22, pos, ref, alt)
    order = {c22: 0, cmt: 1, ABSENT_CONTIG[label]: 2}
    return sorted(records, key=lambda record: (order[record[0]], record[1]))


def real_vcf(tmp_path: Path, label: str, records: list[Real]) -> Path:
    """A one-sample VCF whose header states the assembly's primary-chromosome lengths (domain
    constants, `assembly.ASSEMBLY_LENGTHS`, never read from the store that checks them)."""
    prefix = "chr" if label == "GRCh38" else ""
    contigs = []
    for chrom, lengths in ASSEMBLY_LENGTHS[label].items():
        name = (f"{prefix}M" if prefix else "MT") if chrom == "MT" else f"{prefix}{chrom}"
        contigs.append(f"##contig=<ID={name},length={min(lengths)}>")  # GRCh37's MT: the rCRS, 16569
    contigs.append(f"##contig=<ID={ABSENT_CONTIG[label]}>")
    meta = ["##fileformat=VCFv4.2", *contigs, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1" for c, p, r, a, _ in records]
    path = tmp_path / f"{label}.vcf"
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def allele_kind(ref: str, allele: str, slot: int) -> str:
    """The path an allele takes through the normalization, for the counts."""
    if slot == 0:
        return "ref"
    if len(allele) == len(ref):
        return "snv" if len(ref) == 1 else "mnv"
    return "insertion" if len(allele) > len(ref) else "deletion"


@pytest.mark.parametrize("label", REAL_ASSEMBLIES)
def test_real_store_ids_equal_vrs_python(label: str, test_reference_dir: Path, tmp_path: Path) -> None:
    """Opt-in: every kind of allele, on real chromosome 22 and MT of GRCh37 and GRCh38, minted and
    compared slot by slot with vrs-python's id; every withheld slot has SCHEMA's reason."""
    store = ReferenceStore.open(test_reference_dir)
    digest = real_collection(store, label)
    assert ABSENT_CONTIG[label].removeprefix("chr") not in store.contigs_by_chrom(digest)
    assert ABSENT_CONTIG[label] not in store.contigs_by_chrom(digest)
    seq, mt = real_sequence(store, digest, "22"), real_sequence(store, digest, "MT")
    assert "N" not in seq[WINDOW[0] : WINDOW[1]]
    records = real_records(label, seq, mt, random.Random(271))
    by_pos = {record[1]: record for record in records}  # chromosome 22, MT and the absent contig never share a position here
    assert len(by_pos) == len(records)

    result = mint_vrs_ids(resolve_ref(load(real_vcf(tmp_path, label, records)), store), store)
    assert result.facets.collect()["ref_seqcol_digest"].item() == digest
    calls = vrs_withheld(result.calls).collect()
    assert calls.height == len(records)
    assert calls["vrs_errors"].null_count() == calls.height  # minting itself refused nothing

    proxy = StoreProxy(store, [digest])
    compared: Counter[str] = Counter()
    withheld: Counter[str] = Counter()
    for row in calls.iter_rows(named=True):
        expected = by_pos[row["pos"]][4]
        slots = 1 + len(row["alts"] or [])
        ids = row["vrs_allele_ids"] or [None] * slots
        assert len(ids) == len(row["vrs_withheld"]) == slots
        if expected is None:
            assert row["ref_check"] == "match"
            expected = [None] * slots
        assert row["vrs_withheld"] == expected, f"{label}, a {row['chrom']} row: reasons"
        mintable = [slot for slot, reason in enumerate(expected) if reason is None]
        if mintable:
            oracle = oracle_slots(proxy, store, digest, row, mintable)
            assert [ids[slot] for slot in mintable] == oracle, f"{label}, a {row['chrom']} row: ids"
        alleles = [row["ref"] or row["source_ref"], *(row["alts"] or [])]
        compared.update(allele_kind(alleles[0], alleles[slot], slot) for slot in mintable)
        withheld.update(reason for reason in expected if reason is not None)

    assert sum(withheld.values()) == sum(reason is not None for record in records for reason in record[4] or [])
    minted = sum(compared.values())
    assert any(f"for {minted} of {minted + sum(withheld.values())} allele(s)" in warning for warning in result.warnings)
    logger.info(
        "%s: %d rows, %d alleles equal to vrs-python %s; withheld %s",
        label, calls.height, minted, dict(sorted(compared.items())), dict(sorted(withheld.items())),
    )


def test_real_store_array_rows(test_reference_dir: Path, tmp_path: Path) -> None:
    """Opt-in: a 23andMe-shaped build 37 export placed on real chromosome 22. Called rows mint the
    ids vrs-python gives; an I/D code is `array_code` and a no-call `ref_unresolved`."""
    store = ReferenceStore.open(test_reference_dir)
    digest = real_collection(store, "GRCh37")
    seq = real_sequence(store, digest, "22")
    start = WINDOW[0] + 1
    rows = [
        ("rs1", start + 10, seq[start + 9] * 2),  # hom-ref: the REF slot only
        ("rs2", start + 20, other(seq[start + 19]) + seq[start + 19]),  # het
        ("i3", start + 30, "DI"),  # indel codes
        ("rs4", start + 40, "--"),  # no-call
    ]
    path = tmp_path / "genome.txt"
    path.write_text(PREAMBLE_23 + "".join(f"{rsid}\t22\t{pos}\t{gt}\n" for rsid, pos, gt in rows))
    result = mint_vrs_ids(resolve_ref(load(path), store), store)
    assert result.facets.collect()["ref_seqcol_digest"].item() == digest
    calls = vrs_withheld(result.calls).collect()
    by_id = {row["ids"][0]: row for row in calls.iter_rows(named=True)}

    proxy = StoreProxy(store, [digest])
    for rsid, slots in (("rs1", 1), ("rs2", 2)):
        row = by_id[rsid]
        assert row["vrs_allele_ids"] == oracle_slots(proxy, store, digest, row)
        assert len(row["vrs_allele_ids"]) == slots and row["vrs_withheld"] == [None] * slots
    assert by_id["i3"]["vrs_allele_ids"] is None and by_id["i3"]["vrs_withheld"] == ["array_code"]
    assert by_id["rs4"]["vrs_allele_ids"] is None and by_id["rs4"]["vrs_withheld"] == ["ref_unresolved"]


PGP_SLICE = {"22": 3000, "MT": 500}
"""Rows taken from each PGP-HMS upload, per contig, in file order."""


def test_pgp_hms_slices_equal_vrs_python(pgp_participants: Path, test_reference_dir: Path) -> None:
    """Opt-in: the smallest PGP-HMS VCF uploads the store can place, one per assembly it holds; the
    first rows of chromosome 22 and MT of each, every minted allele compared with vrs-python's id and
    every withheld one checked to carry a reason (those the row alone decides, checked by name).
    Aggregates only.

    Skips when the PGP-HMS tree is not under `JUST_VCF_SAMPLES_DIR`."""
    store = ReferenceStore.open(test_reference_dir)
    wanted = {store.assembly_label(digest) for digest in store.collections()} & set(REAL_ASSEMBLIES)
    sizes = sorted(pgp_participants.glob("*/upload-*/*.vcf*"), key=lambda p: p.stat().st_size)
    found: dict[str, Loaded] = {}
    for path in [p for p in sizes if p.stat().st_size > 5_000_000]:
        if found.keys() == wanted:
            break
        sniffed = sniff(path)
        if sniffed.format != "vcf" or sniffed.build not in wanted - found.keys():
            continue  # the build the header states picks the candidates; resolve_ref still checks it
        try:
            loaded = load(path)
        except ValueError:  # several samples, or a shape the loader refuses: not this test's subject
            continue
        rows = pl.concat([loaded.calls.filter(pl.col("chrom") == chrom).head(n) for chrom, n in PGP_SLICE.items()])
        if rows.filter(pl.col("chrom") == "22").select(pl.len()).collect().item() < 1000:
            continue
        window = Loaded(facets=loaded.facets, calls=rows, blocks=pl.LazyFrame(schema=BLOCKS))
        try:
            resolved = resolve_ref(window, store)
        except ReferenceMismatchError:  # not an assembly the store holds
            continue
        label = store.assembly_label(resolved.facets.collect()["ref_seqcol_digest"].item())
        if label in wanted and label not in found:
            found[label] = mint_vrs_ids(resolved, store)
    if not found:
        pytest.skip("no PGP-HMS VCF upload the store can place")

    for label, result in sorted(found.items()):
        digest = result.facets.collect()["ref_seqcol_digest"].item()
        proxy = StoreProxy(store, [digest])
        calls = vrs_withheld(result.calls).collect()
        compared: Counter[str] = Counter()
        withheld: Counter[str] = Counter()
        for row in calls.iter_rows(named=True):
            alleles = [row["ref"] or row["source_ref"], *(row["alts"] or [])]
            ids = row["vrs_allele_ids"] or [None] * len(row["vrs_withheld"])
            assert len(ids) == len(row["vrs_withheld"])
            for slot, (minted, reason) in enumerate(zip(ids, row["vrs_withheld"])):
                assert (minted is None) != (reason is None)
                if minted is not None:
                    assert [minted] == oracle_slots(proxy, store, digest, row, [slot]), f"{label}, a {row['chrom']} row: slot {slot}"
                    compared[allele_kind(alleles[0], alleles[slot], slot)] += 1
                    continue
                withheld[reason] += 1
                assert reason in WITHHELD_REASONS
                allele = alleles[slot]
                if allele == "*":
                    assert reason == "spanning_deletion"
                elif allele is not None and allele.startswith("<"):
                    assert reason in {"unspecified", "sv_no_vrs_model", "symbolic"}
                elif row["ref_check"] == "mismatch" and allele is not None and re.fullmatch(r"[ACGTNacgtn]*", allele):
                    assert reason == "ref_mismatch"
        assert compared["snv"] + compared["insertion"] + compared["deletion"] > 0
        logger.info(
            "PGP-HMS %s slice: %d rows, %d alleles equal to vrs-python %s; withheld %s",
            label, calls.height, sum(compared.values()), dict(sorted(compared.items())), dict(sorted(withheld.items())),
        )
