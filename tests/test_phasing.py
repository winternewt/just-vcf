"""VRS CisPhasedBlock export (RM284): `cis_phased_blocks` against vrs-python, the oracle.

The oracle works from the records as the test wrote them: it splits each GT text itself, assigns the
alleles to their phase set and haplotype slot, builds vrs-python's own normalized Alleles from the
row's position, REF and allele, and identifies a `CisPhasedBlock` of those inline Alleles. Nothing in
it reads `phased`, `phase_sets` or `vrs_allele_ids`, so an id that agrees was computed twice.
"""

import gzip
import json
import random
import re
from collections import defaultdict
from pathlib import Path

import polars as pl
import pytest
from ga4gh.core import ga4gh_identify
from test_load_cgi import Locus, Seg, master_var_text, snp_locus, var_text
from test_load_cgi import write as write_cgi
from test_vrs import CA_REPEAT, LENGTHS, POLY_A, accessions, base, other, write_fasta
from test_vrs import SEQ as SEQUENCES
from vrs_oracle import StoreProxy, oracle_allele, oracle_block_id

from just_vcf import ReferenceStore, cis_phased_blocks, load, mint_vrs_ids, resolve_ref, sniff
from just_vcf.loading import Loaded
from just_vcf.phasing import CIS_PHASED_BLOCKS, CPB_WITHHELD_REASONS
from just_vcf.schema import BLOCKS
from just_vcf.vrs import WITHHELD_REASONS

Record = tuple[str, int, str, str, str, int | None]
"""`(chrom, pos, ref, alts, GT text, PS)`; a PS of None is written `.` (the chromosome's implicit set)."""


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa", SEQUENCES))
    return s


def phased_vcf(tmp_path: Path, records: list[Record], contigs: dict[str, int], version: str = "4.2", name: str = "phased.vcf") -> Path:
    meta = [f"##fileformat=VCFv{version}", *(f"##contig=<ID={c},length={n}>" for c, n in contigs.items())]
    meta += [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT:PS\t{gt}:{'.' if ps is None else ps}" for c, p, r, a, gt, ps in records]
    path = tmp_path / name
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def snv(chrom: str, pos: int, gt: str, ps: int | None) -> Record:
    return (chrom, pos, base(chrom, pos), other(base(chrom, pos)), gt, ps)


def expected_blocks(records: list[Record], store: ReferenceStore, reference_alleles: bool = False) -> dict[tuple, str | None]:
    """`(chrom, PS, slot)` → vrs-python's block id, or None where fewer than two alleles qualify.

    VCF 4.2 diploid only: a GT without `/` is phased on both slots. `*` is never a member; REF is one
    only with `reference_alleles`. A set holding an uncalled phased slot is left out (the caller
    asserts those separately).
    """
    proxy, sq = StoreProxy(store), accessions(store)
    members: dict[tuple, dict[str, object]] = defaultdict(dict)
    uncalled: set[tuple] = set()
    for chrom, pos, ref, alts, gt, ps in records:
        if "/" in gt:
            continue
        alleles = [ref, *alts.split(",")]
        for slot, index in enumerate(re.split(r"[|/]", gt)):
            key = (chrom, ps, slot)
            if index == ".":
                uncalled.add(key)
                continue
            allele = alleles[int(index)]
            if allele == "*" or (int(index) == 0 and not reference_alleles):
                members.setdefault(key, {})
                continue
            vrs = oracle_allele(proxy, sq[chrom], pos, ref, allele)
            members[key][ga4gh_identify(vrs)] = vrs  # one member per allele, as a set
    return {key: (oracle_block_id(list(alleles.values())) if len(alleles) >= 2 else None) for key, alleles in members.items() if key not in uncalled}


def blocks_by_key(frame: pl.DataFrame) -> dict[tuple, dict]:
    return {(row["chrom"], None if row["phase_set"] is None else int(row["phase_set"]), row["haplotype"]): row for row in frame.iter_rows(named=True)}


def minted(path: Path, store: ReferenceStore) -> Loaded:
    return mint_vrs_ids(resolve_ref(load(path), store), store)


SNV_PS = 100
RECORDS: list[Record] = [
    snv("1", 10, "0|1", SNV_PS),
    snv("1", 20, "1|0", SNV_PS),
    snv("1", 20, "1|0", SNV_PS),  # the same record twice: the haplotype carries the allele once
    ("1", 30, base("1", 30), ",".join(sorted({other(base("1", 30)), base("1", 30) + "G"})), "1|2", SNV_PS),  # multi-allelic
    ("1", CA_REPEAT + 5, base("1", CA_REPEAT + 5), base("1", CA_REPEAT + 5) + "AC", "1|1", SNV_PS),  # insertion in a repeat
    ("1", POLY_A + 6, base("1", POLY_A + 6, 2), base("1", POLY_A + 6), "0|1", SNV_PS),  # a deletion that rolls left
    ("1", 150, base("1", 150, 4), base("1", 150), "1|0", SNV_PS),  # a deletion over 151..153
    ("1", 152, base("1", 152), f"*,{other(base('1', 152))}", "1|2", SNV_PS),  # slot 0 is inside the deletion
    snv("1", 400, "1|0", 300),  # set 300: one ALT per slot, too few for a block
    snv("1", 410, "0|1", 300),
    snv("1", 500, "1|0", None),  # the chromosome's implicit set
    snv("1", 700, "1|0", None),
    snv("1", 800, "0/1", None),  # unphased: in no set
    snv("2", 50, "0|1", SNV_PS),  # PS 100 again, on another chromosome: another set
    ("2", 60, base("2", 60, 3), base("2", 60), "1|1", SNV_PS),
]
CONTIGS = {"1": LENGTHS["1"], "2": LENGTHS["2"]}


def test_block_ids_equal_vrs_python(tmp_path: Path, store: ReferenceStore) -> None:
    blocks = cis_phased_blocks(minted(phased_vcf(tmp_path, RECORDS, CONTIGS), store)).collect()
    expected = expected_blocks(RECORDS, store)
    got = blocks_by_key(blocks)

    assert blocks.schema == CIS_PHASED_BLOCKS
    assert set(got) == set(expected)
    assert {key: row["cpb_id"] for key, row in got.items()} == expected
    assert sum(cpb is not None for cpb in expected.values()) >= 4  # both slots of chr1 PS 100, chr2's, the implicit set's slot 0
    for row in got.values():
        assert (row["cpb_id"] is None) == (row["withheld"] is not None) == (row["members"] is None)
        if row["members"] is not None:
            assert row["members"] == sorted(set(row["members"])) and all(m.startswith("ga4gh:VA.") for m in row["members"])
    # chr1 PS 100 slot 0 covers eight rows (the duplicate and the row inside the deletion included), slot 1 the same rows.
    assert got[("1", SNV_PS, 0)]["sites"] == got[("1", SNV_PS, 1)]["sites"] == 8
    assert {row["withheld"] for key, row in got.items() if key[1] == 300} == {"too_few_members"}
    assert got[("1", None, 1)]["withheld"] == "too_few_members"  # only REF on that slot
    assert all(row["phase_set_kind"] == "ps" for row in got.values())


def test_reference_alleles_join_when_asked(tmp_path: Path, store: ReferenceStore) -> None:
    blocks = cis_phased_blocks(minted(phased_vcf(tmp_path, RECORDS, CONTIGS), store), reference_alleles=True).collect()
    expected = expected_blocks(RECORDS, store, reference_alleles=True)
    assert {key: row["cpb_id"] for key, row in blocks_by_key(blocks).items()} == expected
    assert all(cpb is not None for key, cpb in expected.items())  # REF makes every slot at least two members
    without = blocks_by_key(cis_phased_blocks(minted(phased_vcf(tmp_path, RECORDS, CONTIGS, name="again.vcf"), store)).collect())
    assert blocks_by_key(blocks)[("1", SNV_PS, 0)]["cpb_id"] != without[("1", SNV_PS, 0)]["cpb_id"]


def test_withheld_blocks_say_why(tmp_path: Path, store: ReferenceStore) -> None:
    records: list[Record] = [
        snv("1", 10, ".|1", 7),  # set 7: slot 0 uncalled at one site
        snv("1", 20, "1|1", 7),
        snv("1", 30, "1|1", 7),
        ("1", 310, base("1", 310), "<DEL>", "1|0", 8),  # set 8: slot 0 carries a symbolic allele
        snv("1", 320, "1|1", 8),
        snv("1", 330, "0|1", 8),
    ]
    path = phased_vcf(tmp_path, records, {"1": LENGTHS["1"]})
    got = blocks_by_key(cis_phased_blocks(minted(path, store)).collect())
    assert got[("1", 7, 0)]["withheld"] == "uncalled"
    assert got[("1", 7, 1)]["cpb_id"] == expected_blocks([r for r in records if r[5] == 7], store)[("1", 7, 1)]
    assert got[("1", 8, 0)]["withheld"] == "sv_no_vrs_model"  # RM282: an SV allele has no VRS 2.x Allele model
    assert got[("1", 8, 1)]["cpb_id"] is not None

    # Before minting, every block that needs an ALT's id says minting has not run.
    unminted = blocks_by_key(cis_phased_blocks(resolve_ref(load(path), store)).collect())
    assert unminted[("1", 7, 1)]["withheld"] == unminted[("1", 8, 1)]["withheld"] == "not_run"
    assert unminted[("1", 7, 0)]["withheld"] == "uncalled"  # the uncalled site comes first in row order

    every = [row["withheld"] for row in [*got.values(), *unminted.values()] if row["withheld"] is not None]
    assert set(every) <= CPB_WITHHELD_REASONS and CPB_WITHHELD_REASONS >= WITHHELD_REASONS


def test_a_vcf_44_first_slot_is_phase_unknown_not_absent(tmp_path: Path, store: ReferenceStore) -> None:
    """polars-bio drops a 4.4 GT's leading indicator (P4), so slot 0's phase is unknown: its block is
    reported withheld, never left out as if the set had one haplotype."""
    records = [snv("1", 10, "1|1", 5), snv("1", 20, "1|1", 5)]
    loaded = minted(phased_vcf(tmp_path, records, {"1": LENGTHS["1"]}, version="4.4"), store)
    assert loaded.calls.collect()["phased"].to_list() == [[None, True], [None, True]]
    got = blocks_by_key(cis_phased_blocks(loaded).collect())
    assert set(got) == {("1", 5, 0), ("1", 5, 1)}
    assert got[("1", 5, 0)]["withheld"] == "phase_unknown" and got[("1", 5, 0)]["sites"] == 0
    assert got[("1", 5, 1)]["cpb_id"] == expected_blocks(records, store)[("1", 5, 1)]

    # A slot that carries REF would not be a member whatever its phase: it withholds nothing.
    reference_first = [snv("1", 10, "0|1", 5), snv("1", 20, "0|1", 5)]
    loaded = minted(phased_vcf(tmp_path, reference_first, {"1": LENGTHS["1"]}, version="4.4", name="ref.vcf"), store)
    got = blocks_by_key(cis_phased_blocks(loaded).collect())
    assert got[("1", 5, 0)]["withheld"] == "too_few_members"
    assert got[("1", 5, 1)]["cpb_id"] == expected_blocks(reference_first, store)[("1", 5, 1)]
    assert blocks_by_key(cis_phased_blocks(loaded, reference_alleles=True).collect())[("1", 5, 0)]["withheld"] == "phase_unknown"


def test_a_psl_set_is_one_haplotype_and_spans_chromosomes(tmp_path: Path, store: ReferenceStore) -> None:
    """VCF 4.5 PSL names a molecule, not a slot: `chr19*5*1` can hold slot 0 of one record and slot 2
    of another, on another chromosome. No loader reads PSL yet, so the frame is rewritten by hand."""
    records = [snv("1", 10, "1|1", 1), snv("2", 50, "1|1", 1), snv("1", 20, "1|1", 1)]
    loaded = minted(phased_vcf(tmp_path, records, CONTIGS), store)
    names = {("1", 10): ["a", "b"], ("1", 20): ["b", "a"], ("2", 50): [None, "a"]}
    table = loaded.calls.collect()
    table = table.with_columns(
        pl.Series("phase_sets", [names[(c, p)] for c, p in zip(table["chrom"], table["pos"])], dtype=table.schema["phase_sets"]),
        pl.Series("phase_set_kinds", [[None if n is None else "psl" for n in names[(c, p)]] for c, p in zip(table["chrom"], table["pos"])], dtype=table.schema["phase_set_kinds"]),
    )
    table = table.with_columns(pl.when(pl.col("chrom") == "2").then(pl.lit([False, True])).otherwise(pl.col("phased")).alias("phased"))
    blocks = cis_phased_blocks(Loaded(facets=loaded.facets, calls=table.lazy(), blocks=pl.LazyFrame(schema=BLOCKS))).collect()

    assert blocks["phase_set"].to_list() == ["a", "b"]  # first row starts both; the tie goes to the name
    assert blocks["chrom"].to_list() == [None, None] and blocks["haplotype"].to_list() == [None, None]
    proxy, sq = StoreProxy(store), accessions(store)

    def allele(chrom: str, pos: int):
        return oracle_allele(proxy, sq[chrom], pos, base(chrom, pos), other(base(chrom, pos)))

    a = blocks.filter(pl.col("phase_set") == "a").row(0, named=True)
    assert a["sites"] == 3 and a["cpb_id"] == oracle_block_id([allele("1", 10), allele("1", 20), allele("2", 50)])
    b = blocks.filter(pl.col("phase_set") == "b").row(0, named=True)
    assert b["sites"] == 2 and b["cpb_id"] == oracle_block_id([allele("1", 10), allele("1", 20)])


def test_blocks_are_ordered_by_first_row_and_deterministic(tmp_path: Path, store: ReferenceStore) -> None:
    first = cis_phased_blocks(minted(phased_vcf(tmp_path, RECORDS, CONTIGS), store)).collect()
    again = cis_phased_blocks(minted(phased_vcf(tmp_path, RECORDS, CONTIGS), store)).collect()
    shuffled = RECORDS[:]
    random.Random(5).shuffle(shuffled)
    reordered = cis_phased_blocks(minted(phased_vcf(tmp_path, shuffled, CONTIGS, name="shuffled.vcf"), store)).collect()
    assert first.equals(again)
    assert first.drop("facet_id").equals(reordered.drop("facet_id"))
    # File order: chr1 PS 100, set 300, the implicit set, then chr2 PS 100; slot 0 before slot 1.
    assert [(r["chrom"], r["phase_set"], r["haplotype"]) for r in first.iter_rows(named=True)] == [
        ("1", "100", 0), ("1", "100", 1), ("1", "300", 0), ("1", "300", 1), ("1", None, 0), ("1", None, 1), ("2", "100", 0), ("2", "100", 1),
    ]


def cgi_linked_loci() -> list[Locus]:
    """CGI loci over the synthetic reference: two molecules of one 2.5 phase block with the ALTs on
    either slot; an integer id over an unanchored deletion and an SNV; one over a haplotype CGI could
    not call; and the same integer on chr2, another molecule."""
    def snp(chrom: str, pos: int, links: tuple[str, str], alt_slot: int = 0) -> Locus:
        return snp_locus(f"chr{chrom}", pos - 1, base(chrom, pos), other(base(chrom, pos)), links, alt_slot)

    return [
        snp("1", 10, ("Phased_1_1_1", "Phased_1_1_2")),
        snp("1", 20, ("Phased_1_1_2", "Phased_1_1_1"), alt_slot=1),
        snp("1", 30, ("Phased_1_1_2", "Phased_1_1_1")),
        snp("1", 40, ("Phased_1_1_1", "Phased_1_1_2"), alt_slot=1),
        Locus("chr1", 149, 150, haps=((Seg(149, 150, "del", base("1", 150), "", link="7"),), (Seg(149, 150, "ref", base("1", 150), base("1", 150)),))),
        snp("1", 160, ("7", "")),
        Locus("chr1", 199, 200, haps=((Seg(199, 200, "no-call", base("1", 200), "?", filter="", link="9"),), (Seg(199, 200, "snp", base("1", 200), other(base("1", 200))),))),
        snp("1", 210, ("9", "")),
        snp("2", 50, ("7", "")),
        snp("2", 60, ("7", "")),
    ]


def test_cgi_hap_link_sets_equal_vrs_python(tmp_path: Path) -> None:
    """A hapLink set is one molecule on one chromosome: one block per id and chromosome, `haplotype`
    null, whichever slot each ALT sits in. `var` and `masterVarBeta` (which writes the ALT slot first)
    give the same blocks, and the ids are vrs-python's for the same alleles."""
    # CGI's files state GRCh37 and no contig lengths: the store says which assembly it holds.
    store = ReferenceStore.open(tmp_path / "grch37")
    store.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQUENCES), assembly_label="GRCh37")
    proxy, sq = StoreProxy(store), accessions(store)

    def snv_allele(chrom: str, pos: int):
        return oracle_allele(proxy, sq[chrom], pos, base(chrom, pos), other(base(chrom, pos)))

    deletion = oracle_allele(proxy, sq["1"], 149, base("1", 149, 2), base("1", 149))  # CGI's `""` at 150, anchored
    expected = {
        ("1", "Phased_1_1_1"): oracle_block_id([snv_allele("1", 10), snv_allele("1", 20)]),
        ("1", "Phased_1_1_2"): oracle_block_id([snv_allele("1", 30), snv_allele("1", 40)]),
        ("1", "7"): oracle_block_id([deletion, snv_allele("1", 160)]),
        ("2", "7"): oracle_block_id([snv_allele("2", 50), snv_allele("2", 60)]),
    }
    seen = []
    for name, text in (("var-GS000000005-ASM.tsv.bz2", var_text(cgi_linked_loci())), ("masterVarBeta-GS000000005-ASM.tsv.bz2", master_var_text(cgi_linked_loci()))):
        blocks = cis_phased_blocks(minted(write_cgi(tmp_path, text, name), store)).collect()
        assert set(blocks["phase_set_kind"].to_list()) == {"cgi_haplink"} and blocks["haplotype"].null_count() == blocks.height
        got = {(row["chrom"], row["phase_set"]): row for row in blocks.iter_rows(named=True)}
        assert set(got) == {*expected, ("1", "9")}, name
        for key, cpb_id in expected.items():
            # Each of the four block loci has a slot on each molecule; the REF ones are sites, not members.
            sites = 4 if key[1].startswith("Phased") else 2
            assert (got[key]["cpb_id"], got[key]["withheld"], got[key]["sites"]) == (cpb_id, None, sites), (name, key)
        assert (got[("1", "9")]["cpb_id"], got[("1", "9")]["withheld"]) == (None, "uncalled"), "a molecule CGI could not call is not shortened"
        seen.append(blocks.drop("facet_id"))
    assert seen[0].equals(seen[1])


def test_unphased_facets_give_no_blocks(tmp_path: Path, store: ReferenceStore) -> None:
    records = [snv("1", 10, "0/1", None), snv("1", 20, "1/1", None)]
    assert cis_phased_blocks(minted(phased_vcf(tmp_path, records, {"1": LENGTHS["1"]}), store)).collect().height == 0


def test_haploid_calls_are_one_molecule(tmp_path: Path, store: ReferenceStore) -> None:
    """A haploid GT has no separator, so VCF reads it phased: on a chromosome with one copy (male X
    outside the PARs, Y, MT) its ALTs are in cis by ploidy, and the implicit set holds them, beside
    the unphased diploid rows, which are in no set."""
    records = [snv("1", 10, "1", None), snv("1", 20, "0/1", None), snv("1", 30, "1", None), snv("1", 40, "0", None)]
    got = blocks_by_key(cis_phased_blocks(minted(phased_vcf(tmp_path, records, {"1": LENGTHS["1"]}), store)).collect())
    assert set(got) == {("1", None, 0)}
    assert got[("1", None, 0)]["sites"] == 3
    assert got[("1", None, 0)]["cpb_id"] == expected_blocks([r for r in records if "/" not in r[4]], store)[("1", None, 0)]


def first_gt_phased(path: Path) -> bool:
    """Whether the first record's genotype is phased: a cheap look before a whole load. A zipped
    upload is passed over."""
    with path.open("rb") as raw:
        magic = raw.read(2)
    if magic == b"PK":
        return False
    opener = gzip.open if magic == b"\x1f\x8b" else open
    with opener(path, "rt", errors="replace") as text:
        for line in text:
            if not line.startswith("#"):
                fields = line.split("\t")
                return len(fields) > 9 and "|" in fields[9].split(":")[0]
    return False


def test_a_phased_pgp_hms_facet_equals_vrs_python(pgp_participants: Path, test_reference_dir: Path) -> None:
    """Opt-in: the smallest Genes for Good VCF upload whose chromosome 22 is phased (statistically
    phased imputation, no PS: one implicit set per chromosome), the first 1500 phased rows of
    chromosome 22. The oracle reads each row's GT text and alleles itself. Aggregates only.

    Skips when the PGP-HMS tree is not under `JUST_VCF_SAMPLES_DIR` (it may live on another copy of
    the samples directory; point the variable there)."""
    store = ReferenceStore.open(test_reference_dir)
    manifest = pgp_participants / "manifest.jsonl"
    if not manifest.is_file():
        pytest.skip(f"{manifest} is absent")
    candidates = []
    for line in manifest.read_text().splitlines():
        entry = json.loads(line)
        if entry["item"].get("kind") != "Genes for Good" or not entry.get("path"):
            continue
        path = pgp_participants / entry["path"].split("/participants/", 1)[-1]
        if path.is_file() and sniff(path).format == "vcf" and first_gt_phased(path):
            candidates.append(path)
    window = None
    for path in sorted(candidates, key=lambda p: p.stat().st_size):
        loaded = load(path)
        phased = pl.col("phased").list.eval(pl.element().fill_null(False)).list.all().fill_null(False)
        rows = loaded.calls.filter(pl.col("chrom").is_in(["22", "chr22"]) & phased).head(1500)
        if rows.select(pl.len()).collect().item() == 1500:
            window = Loaded(facets=loaded.facets, calls=rows, blocks=pl.LazyFrame(schema=BLOCKS))
            break
    if window is None:
        pytest.skip("no Genes for Good VCF upload with a phased chromosome 22")
    result = mint_vrs_ids(resolve_ref(window, store), store)
    blocks = cis_phased_blocks(result).collect()
    calls = result.calls.collect()

    digest = result.facets.collect()["ref_seqcol_digest"].item()
    contig = next(c for c in store.contigs(digest) if c.name in ("22", "chr22"))
    proxy = StoreProxy(store, [digest])
    members: dict[int, dict[str, object]] = {0: {}, 1: {}}
    for row in calls.iter_rows(named=True):
        alleles = [row["source_ref"], *row["alts"]]
        for slot, index in enumerate(row["source_gt"].split("|")):
            allele = alleles[int(index)]
            if int(index) == 0 or allele == "*":
                continue
            if allele.startswith("<"):
                members[slot] = None  # a symbolic allele: the block must be withheld
            if members[slot] is not None:
                vrs = oracle_allele(proxy, f"SQ.{contig.digest}", row["pos"], row["source_ref"], allele)
                members[slot][ga4gh_identify(vrs)] = vrs

    assert blocks.height == 2 and blocks["phase_set"].null_count() == 2 and blocks["haplotype"].to_list() == [0, 1]
    assert blocks["sites"].to_list() == [1500, 1500]
    for slot in (0, 1):
        row = blocks.row(slot, named=True)
        if members[slot] is None:
            assert row["withheld"] in {"symbolic", "sv_no_vrs_model"}
        else:
            assert row["cpb_id"] == oracle_block_id(list(members[slot].values()))
            assert len(row["members"]) == len(members[slot]) > 100
    assert blocks["cpb_id"].null_count() < 2


PslRecord = tuple[str, int, str, str, str, str]
"""`(chrom, pos, ref, alts, GT text, PSL text)` of a VCF 4.5 record."""


def psl_vcf(tmp_path: Path, records: list[PslRecord], name: str = "psl.vcf") -> Path:
    meta = ["##fileformat=VCFv4.5", *(f"##contig=<ID={c},length={n}>" for c, n in CONTIGS.items())]
    meta += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PSL,Number=P,Type=String,Description="Phase set list">',
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT:PSL\t{gt}:{names}" for c, p, r, a, gt, names in records]
    path = tmp_path / name
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def psl(chrom: str, pos: int, gt: str, names: str) -> PslRecord:
    return (chrom, pos, base(chrom, pos), other(base(chrom, pos)), gt, names)


def withheld_at(row: dict) -> tuple:
    return tuple(row[f"withheld_{name}"] for name in ("chrom", "pos", "allele"))


def test_a_loaded_psl_facet_equals_vrs_python(tmp_path: Path, store: ReferenceStore) -> None:
    """PSL sets read by the VCF loader, not rewritten by hand: each named slot is in its set wherever
    the record sits, a first slot named by PSL is phased (a VCF 4.5 unphased allele has `.`), and a
    phased slot with `.` is in the chromosome's implicit set. The oracle reads the GT and PSL text."""
    records = [
        psl("1", 10, "|1/0", "a,."),  # slot 0 in `a`; slot 1 unphased
        psl("1", 20, "1|1", ".,a"),  # slot 0 says nothing; slot 1 in `a`
        psl("1", 30, "1|1", "b,."),  # slot 1 phased with `.`: the implicit set
        psl("1", 40, "/1|1", ".,b"),
        psl("2", 50, "1|1", "b,a"),  # both sets cross to chromosome 2
    ]
    loaded = minted(psl_vcf(tmp_path, records), store)
    blocks = cis_phased_blocks(loaded).collect()

    proxy, sq = StoreProxy(store), accessions(store)
    expected: dict[str, dict[str, object]] = defaultdict(dict)
    for chrom, pos, ref, alt, gt, names in records:
        for index, name in zip(re.split(r"[|/]", gt.lstrip("|/")), names.split(",")):
            if name != "." and index != "0":
                vrs = oracle_allele(proxy, sq[chrom], pos, ref, alt)
                expected[name][ga4gh_identify(vrs)] = vrs
    got = {row["phase_set"]: row for row in blocks.filter(pl.col("phase_set_kind") == "psl").iter_rows(named=True)}
    assert set(got) == set(expected) == {"a", "b"}
    for name, alleles in expected.items():
        assert got[name]["cpb_id"] == oracle_block_id(list(alleles.values())), name
        assert got[name]["chrom"] is None and got[name]["haplotype"] is None
    assert got["a"]["sites"] == got["b"]["sites"] == 3
    implicit = blocks.filter(pl.col("phase_set_kind") == "ps").rows(named=True)
    assert [(r["chrom"], r["phase_set"], r["haplotype"], r["sites"], r["withheld"]) for r in implicit] == [("1", None, 1, 1, "too_few_members")]


def test_a_withheld_block_names_the_first_allele_that_withheld_it(tmp_path: Path, store: ReferenceStore) -> None:
    """On a whole-chromosome implicit set one symbolic ALT withholds the haplotype's block (RM284);
    the block names that record and allele, the first in row order when there are several."""
    records: list[Record] = [
        snv("1", 10, "1|0", None),
        ("1", 310, base("1", 310), "<DEL>", "1|0", None),
        snv("1", 320, "1|1", None),
        ("1", 500, base("1", 500), "<DUP>", "1|1", None),
        snv("1", 600, ".|1", 9),  # set 9: slot 0 uncalled
        snv("1", 610, "1|1", 9),
    ]
    got = blocks_by_key(cis_phased_blocks(minted(phased_vcf(tmp_path, records, {"1": LENGTHS["1"]}), store)).collect())
    assert got[("1", None, 0)]["withheld"] == "sv_no_vrs_model" and withheld_at(got[("1", None, 0)]) == ("1", 310, "<DEL>")
    assert got[("1", None, 1)]["withheld"] == "sv_no_vrs_model" and withheld_at(got[("1", None, 1)]) == ("1", 500, "<DUP>")
    assert got[("1", 9, 0)]["withheld"] == "uncalled" and withheld_at(got[("1", 9, 0)]) == ("1", 600, None)
    assert got[("1", 9, 1)]["cpb_id"] is not None and withheld_at(got[("1", 9, 1)]) == (None, None, None)
    one = blocks_by_key(cis_phased_blocks(minted(phased_vcf(tmp_path, records[:1], {"1": LENGTHS["1"]}, name="one.vcf"), store)).collect())
    assert one[("1", None, 0)]["withheld"] == "too_few_members" and withheld_at(one[("1", None, 0)]) == (None, None, None)
