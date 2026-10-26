"""RM6: `liftover` lifts a facet to another assembly through a chain file, as a derived facet.

A synthetic source assembly, a target built from it by known edits, and the chain file describing
those edits, all under tmp_path. Expected positions and bases are computed from the construction;
pyliftover (the `[liftover]` extra) is a second, independent reading of the chain.
"""

import gzip
import hashlib
import random
from collections import Counter
from pathlib import Path

import polars as pl
import pytest

from just_vcf import Loaded, ReferenceMismatchError, ReferenceStore, liftover, load, mint_vrs_ids, normalize, read_chain, resolve_ref, vrs_withheld, write_vcf
from just_vcf.liftover import RELATION, REPORT, UNLIFTED_REASONS, _segments
from just_vcf.schema import BLOCKS, CALLS, FACETS, canonical_order

rng = random.Random(7)


def random_bases(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT", "TGCA"))


def comp(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))


def other(*bases: str) -> str:
    return next(b for b in "ACGT" if b not in bases)


S1 = random_bases(1000)
S2 = random_bases(500)
SOURCE = {"1": S1, "2": S2}

# Positions on source contig 1 (1-based) that the target edits, and what it writes there.
SWAP_POS, MISMATCH_POS, MULTI_SWAP_POS = 251, 301, 271
SWAP_BASE = other(S1[SWAP_POS - 1])
MISMATCH_ALT = other(S1[MISMATCH_POS - 1])
MISMATCH_BASE = other(S1[MISMATCH_POS - 1], MISMATCH_ALT)
MULTI_ALT_1 = other(S1[MULTI_SWAP_POS - 1])
MULTI_BASE = other(S1[MULTI_SWAP_POS - 1], MULTI_ALT_1)


def _edited(segment: str, offset: int) -> str:
    """Source bases from 0-based `offset`, with the target's base changes applied."""
    chars = list(segment)
    for pos, base in ((SWAP_POS, SWAP_BASE), (MISMATCH_POS, MISMATCH_BASE), (MULTI_SWAP_POS, MULTI_BASE)):
        if offset <= pos - 1 < offset + len(chars):
            chars[pos - 1 - offset] = base
    return "".join(chars)


# Target contig 1: 50 bases of filler, source 100..400 (block A), source 410..600 (block B; source
# 400..410 is deleted in the target), filler. Target contig 2: filler, source 700..800 reverse-
# complemented (a minus-strand chain), filler, source 150..170 again (a second chain over block A).
F1_HEAD, F1_TAIL = random_bases(50), random_bases(60)
T1 = F1_HEAD + _edited(S1[100:400], 100) + S1[410:600] + F1_TAIL
F2_HEAD, F2_MID, F2_TAIL = random_bases(30), random_bases(40), random_bases(25)
T2 = F2_HEAD + revcomp(S1[700:800]) + F2_MID + S1[150:170] + F2_TAIL
TARGET = {"1": T1, "2": T2}
MINUS_F0 = len(F2_HEAD)  # 0-based forward start of the reverse-complemented segment on target 2
DUP_F0 = len(F2_HEAD) + 100 + len(F2_MID)


def expected_target(pos: int) -> tuple[str, int, str]:
    """The construction's own mapping of a 1-based source position on contig 1: (chrom, pos, strand)."""
    i = pos - 1
    if 100 <= i < 400:
        return "1", i - 100 + 50 + 1, "+"
    if 410 <= i < 600:
        return "1", i - 410 + 350 + 1, "+"
    if 700 <= i < 800:
        return "2", MINUS_F0 + (799 - i) + 1, "-"
    raise AssertionError(pos)


def chain_text(source_1_size: int = len(S1)) -> str:
    q2 = len(T2)
    return "\n".join(
        [
            f"chain 1000 1 {source_1_size} + 100 600 1 {len(T1)} + 50 540 1",
            "300 10 0",
            "190",
            "",
            f"chain 900 1 {source_1_size} + 700 800 2 {q2} - {q2 - MINUS_F0 - 100} {q2 - MINUS_F0} 2",
            "100",
            "",
            f"chain 800 1 {source_1_size} + 150 170 2 {q2} + {DUP_F0} {DUP_F0 + 20} 3",
            "20",
            "",
            f"chain 700 1 {source_1_size} + 850 860 chrUn_q 100 + 0 10 4",
            "10",
            "",
        ]
    )


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()))
    return path


HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SOURCE.items()],
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
    '##INFO=<ID=CN,Number=A,Type=Integer,Description="Copy number per ALT">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def s1(pos: int, length: int = 1) -> str:
    return S1[pos - 1 : pos - 1 + length]


def row(pos: int, ref: str, alt: str, gt: str = "0/1", chrom: str = "1", rsid: str = ".", info: str = ".") -> str:
    return f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t50\tPASS\t{info}\tGT\t{gt}"


# One row per case; the rsid names the case so the report can be read by it.
CASES = {
    "plus_snv": row(201, s1(201), other(s1(201)), "0/1", rsid="plus_snv"),
    "swap": f"1\t{SWAP_POS}\tswap\t{s1(SWAP_POS)}\t{SWAP_BASE}\t50\tPASS\t.\tGT:PL:DS\t0/1:10,0,40:0.9",
    "multi_swap": row(MULTI_SWAP_POS, s1(MULTI_SWAP_POS), f"{MULTI_ALT_1},{MULTI_BASE}", "1/2", rsid="multi_swap"),
    "mismatch": row(MISMATCH_POS, s1(MISMATCH_POS), MISMATCH_ALT, rsid="mismatch"),
    "plus_deletion": row(450, s1(450, 3), s1(450), rsid="plus_deletion"),
    "plus_insertion": row(500, s1(500), s1(500) + "GT", "1/1", rsid="plus_insertion"),
    "across_edge": row(398, s1(398, 5), s1(398), rsid="across_edge"),
    "unmapped": row(50, s1(50), other(s1(50)), rsid="unmapped"),
    "between_blocks": row(650, s1(650), other(s1(650)), rsid="between_blocks"),
    "ambiguous": row(160, s1(160), other(s1(160)), rsid="ambiguous"),
    "minus_snv": row(751, s1(751), other(s1(751)), "0|1", rsid="minus_snv"),
    "minus_mnv": row(761, s1(761, 2), other(s1(761)) + s1(762), rsid="minus_mnv"),
    "minus_indel": row(770, s1(770), s1(770) + "A", rsid="minus_indel"),
    "minus_deletion": row(780, s1(780, 3), s1(780), "1/1", rsid="minus_deletion"),
    "minus_mixed": row(790, s1(790), f"{other(s1(790))},{s1(790)}G", "1/2", rsid="minus_mixed"),
    "minus_edge": row(800, s1(800), s1(800) + "A", rsid="minus_edge"),
    "source_mismatch": row(210, other(s1(210)), s1(210), rsid="source_mismatch"),
    "breakend": row(220, s1(220), s1(220) + "[2:20[", rsid="breakend"),
    "star": row(230, s1(230), f"{other(s1(230))},*", "0/1", rsid="star"),
    "unknown_target": row(855, s1(855), other(s1(855)), rsid="unknown_target"),
    "contig_2": row(20, S2[19], other(S2[19]), chrom="2", rsid="contig_2"),
}
REANCHORED = {"minus_indel", "minus_deletion", "minus_edge"}
"""Minus-strand indels whose alleles share their first base: after reverse-complementing, that base is
at their end, so the row is re-anchored on the target's base before the lifted span."""

EXPECTED_REASON = {
    "plus_snv": None,
    "swap": None,
    "multi_swap": None,
    "mismatch": "ref_mismatch_after_lift",
    "plus_deletion": None,
    "plus_insertion": None,
    "across_edge": "split_across_blocks",
    "unmapped": "no_chain_block",
    "between_blocks": "no_chain_block",
    "ambiguous": "multiple_targets",
    "minus_snv": None,
    "minus_mnv": None,
    "minus_indel": None,
    "minus_deletion": None,
    "minus_mixed": None,
    "minus_edge": None,
    "source_mismatch": "source_ref_mismatch",
    "breakend": "breakend_mate_not_lifted",  # its mate, 2:20, is on a contig the chain does not cover
    "star": "allele_not_sequence",
    "unknown_target": "target_contig_unmatched",
    "contig_2": "contig_not_in_chain",
}


@pytest.fixture
def setup(tmp_path: Path) -> tuple[ReferenceStore, str, str, Path]:
    store = ReferenceStore.open(tmp_path / "store")
    source = store.add_fasta(write_fasta(tmp_path / "source.fa", SOURCE), assembly_label="GRCh37")
    target = store.add_fasta(write_fasta(tmp_path / "target.fa", TARGET), assembly_label="GRCh38")
    chain = tmp_path / "source_to_target.chain.gz"
    chain.write_bytes(gzip.compress(chain_text().encode()))
    return store, source, target, chain


def lifted_from_cases(tmp_path: Path, store: ReferenceStore, chain: Path, cases: dict[str, str] = CASES):
    order = sorted(cases.values(), key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    vcf = tmp_path / "source.vcf"
    vcf.write_text("\n".join([*HEADER, *order]) + "\n")
    resolved = resolve_ref(load(vcf), store)
    return resolved, liftover(resolved, chain, store)


def by_case(frame: pl.DataFrame) -> dict[str, dict]:
    return {r["ids"][0]: r for r in frame.iter_rows(named=True)}


def test_every_row_gets_its_reason_and_lifted_rows_land_where_the_construction_puts_them(tmp_path: Path, setup) -> None:
    store, _, target, chain = setup
    resolved, lifted = lifted_from_cases(tmp_path, store, chain)
    report = by_case(lifted.report.collect())
    assert {case: r["lift_reason"] for case, r in report.items()} == EXPECTED_REASON
    assert set(lifted.reasons) <= UNLIFTED_REASONS
    assert sum(lifted.reasons.values()) == sum(reason is not None for reason in EXPECTED_REASON.values())

    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == {case for case, reason in EXPECTED_REASON.items() if reason is None}
    for case, r in calls.items():
        source_pos, source_len = report[case]["pos"], len(report[case]["ref"])
        if case == "minus_edge":
            chrom, pos, strand = "2", MINUS_F0, "-"  # the base before the span is filler, past the block
        elif case in REANCHORED:
            # Re-anchored: the new first base is the image of the source base after the span.
            chrom, pos, strand = expected_target(source_pos + source_len)
        else:
            chrom, pos, strand = expected_target(source_pos)
            if strand == "-":
                pos -= source_len - 1  # a span's forward start is its last source base's image
        assert (r["chrom"], r["pos"], report[case]["lift_strand"]) == (chrom, pos, strand), case
        assert (report[case]["lift_chrom"], report[case]["lift_pos"]) == (chrom, pos)
        # REF is read from the target over the lifted span, never carried over.
        assert r["ref"] == TARGET[chrom][pos - 1 : pos - 1 + len(r["ref"])]
        assert r["end_pos"] == pos + len(r["ref"]) - 1
        assert r["source_ref"] is None and r["ref_check"] is None and r["source_gt"] is None
        assert r["vrs_allele_ids"] is None
        # gt_alleles stay the person's alleles, on the target strand, and gt indexes [ref] + alts.
        alleles = [r["ref"], *r["alts"]]
        assert [alleles[g] for g in r["gt"]] == r["gt_alleles"]


def test_strand_and_swaps(tmp_path: Path, setup) -> None:
    store, _, _, chain = setup
    _, lifted = lifted_from_cases(tmp_path, store, chain)
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())

    plus = calls["plus_snv"]
    assert plus["gt"] == [0, 1] and plus["alts"] == [CASES["plus_snv"].split("\t")[4]]
    assert report["plus_snv"]["lift_swapped"] is False

    swap = calls["swap"]
    assert swap["ref"] == SWAP_BASE and swap["alts"] == [s1(SWAP_POS)]
    assert swap["gt"] == [1, 0] and swap["gt_alleles"] == [s1(SWAP_POS), SWAP_BASE]
    assert report["swap"]["lift_swapped"] is True
    # Vectors stay indexed over source_alleles, never reordered; the dosage counted the other allele.
    assert swap["gl"] == report["swap"]["gl"] and swap["source_alleles"] == report["swap"]["source_alleles"]
    assert report["swap"]["dosage"] is not None and swap["dosage"] is None

    multi = calls["multi_swap"]
    assert multi["ref"] == MULTI_BASE and multi["alts"] == [MULTI_ALT_1, s1(MULTI_SWAP_POS)]
    assert multi["gt"] == [1, 0]

    minus = calls["minus_snv"]
    assert minus["ref"] == revcomp(s1(751)) and minus["alts"] == [revcomp(report["minus_snv"]["alts"][0])]
    assert minus["gt"] == [0, 1]
    # It landed on another contig (target 2): a PS set is per contig, so its phase is cleared, not merged.
    assert report["minus_snv"]["phased"] is not None and minus["phased"] is None
    mnv = calls["minus_mnv"]
    assert mnv["ref"] == revcomp(s1(761, 2)) and mnv["alts"] == [revcomp(report["minus_mnv"]["alts"][0])]

    indel = calls["plus_deletion"]
    assert (indel["ref"], indel["alts"]) == (s1(450, 3), [s1(450)])


def test_haplotype_scores_ride_along_per_slot(tmp_path: Path, setup) -> None:
    """`hq` (a CGI facet's allele scores) is per GT slot, so it is carried as it was, swapped or not:
    a swap renumbers the alleles, not the slots."""
    store, _, _, chain = setup
    cases = {
        "plus_snv": row(201, s1(201), other(s1(201)), "0/1", rsid="plus_snv"),
        "swap": row(SWAP_POS, s1(SWAP_POS), SWAP_BASE, "0/1", rsid="swap"),
    }
    resolved, _ = lifted_from_cases(tmp_path, store, chain, cases)
    scores = pl.int_ranges(30, 32, dtype=pl.UInt16).alias("hq")
    scored = Loaded(facets=resolved.facets, calls=resolved.calls.with_columns(scores + pl.int_range(pl.len(), dtype=pl.UInt16)), blocks=resolved.blocks)
    source = by_case(scored.calls.collect())
    lifted = liftover(scored, chain, store)
    calls, report = by_case(lifted.loaded.calls.collect()), by_case(lifted.report.collect())
    assert report["swap"]["lift_swapped"] and not report["plus_snv"]["lift_swapped"]
    assert {case: calls[case]["hq"] for case in cases} == {case: source[case]["hq"] for case in cases}
    assert all(source[case]["hq"] is not None for case in cases)


def test_a_swapped_row_withholds_its_per_alt_copy_numbers(tmp_path: Path, setup) -> None:
    """INFO CN is one per ALT: after a swap the ALT slot holds the old REF, so the numbers are withheld."""
    store, _, _, chain = setup
    cases = {
        "plus_snv": row(201, s1(201), other(s1(201)), "0/1", rsid="plus_snv", info="CN=2"),
        "swap": row(SWAP_POS, s1(SWAP_POS), SWAP_BASE, "0/1", rsid="swap", info="CN=3"),
        "multi_swap": row(MULTI_SWAP_POS, s1(MULTI_SWAP_POS), f"{MULTI_ALT_1},{MULTI_BASE}", "1/2", rsid="multi_swap", info="CN=1,4"),
    }
    resolved, lifted = lifted_from_cases(tmp_path, store, chain, cases)
    source = by_case(resolved.calls.collect())
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())
    assert all(source[case]["alt_copy_numbers"] is not None for case in cases)
    for case in cases:
        swapped = report[case]["lift_swapped"]
        assert swapped is (case != "plus_snv")
        assert calls[case]["alt_copy_numbers"] == (None if swapped else source[case]["alt_copy_numbers"])
    assert any(line.startswith("2 swapped row(s): alt_copy_numbers withheld") for line in lifted.loaded.warnings)


def test_minus_strand_indels_are_reanchored_on_the_base_before(tmp_path: Path, setup) -> None:
    """Hand-computed from the construction: source 701..800 is target 2's MINUS_F0 + 1..100, reversed."""
    store, _, _, chain = setup
    _, lifted = lifted_from_cases(tmp_path, store, chain)
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())

    # Insertion of A after source base 770: on the target, T inserted after the image of source 771.
    ins = calls["minus_indel"]
    assert (ins["chrom"], ins["pos"]) == ("2", MINUS_F0 + 30)
    assert (ins["ref"], ins["alts"], ins["gt"]) == (comp(s1(771)), [comp(s1(771)) + "T"], [0, 1])
    assert ins["end_pos"] == ins["pos"]
    assert (report["minus_indel"]["lift_pos"], report["minus_indel"]["lift_strand"]) == (MINUS_F0 + 30, "-")

    # Deletion of source 781..782 (anchored on 780): the target deletes their images, anchored on 783's.
    deletion = calls["minus_deletion"]
    assert (deletion["chrom"], deletion["pos"]) == ("2", MINUS_F0 + 18)
    assert (deletion["ref"], deletion["alts"]) == (revcomp(s1(781, 3)), [comp(s1(783))])
    assert deletion["gt_alleles"] == [comp(s1(783))] * 2 and deletion["end_pos"] == MINUS_F0 + 20

    # An SNV beside an insertion shares no first base: nothing to re-anchor, every allele stays a
    # sequence, so it lifts reverse-complemented in place like an MNV.
    mixed = calls["minus_mixed"]
    assert (mixed["chrom"], mixed["pos"]) == expected_target(790)[:2]
    assert (mixed["ref"], mixed["alts"]) == (comp(s1(790)), [comp(other(s1(790))), "C" + comp(s1(790))])
    assert mixed["gt_alleles"] == mixed["alts"]

    # At the block's edge the base before the lifted span is outside the block (target filler). REF is
    # read from the target, so the record states the target's own base there.
    edge = calls["minus_edge"]
    assert (edge["chrom"], edge["pos"], edge["ref"], edge["alts"]) == ("2", MINUS_F0, F2_HEAD[-1], [F2_HEAD[-1] + "T"])
    assert any("re-anchored" in w for w in lifted.loaded.warnings)


def test_a_minus_strand_indel_at_the_target_contigs_first_base_stays_behind(tmp_path: Path) -> None:
    """Target 1 is source 1 reverse-complemented whole: an insertion after the last source base lands
    at the target's first base, with nothing before it to anchor on."""
    source_seq = random_bases(200)
    store = ReferenceStore.open(tmp_path / "store")
    source = store.add_fasta(write_fasta(tmp_path / "s.fa", {"1": source_seq}), assembly_label="GRCh37")
    target = store.add_fasta(write_fasta(tmp_path / "t.fa", {"1": revcomp(source_seq)}), assembly_label="GRCh38")
    chain = tmp_path / "flip.chain"
    chain.write_text("chain 100 1 200 + 0 200 1 200 - 0 200 1\n200\n\n")
    vcf = tmp_path / "flip.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=200>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        f"1\t100\tinner\t{source_seq[99]}\t{source_seq[99]}C\t50\tPASS\t.\tGT\t0/1",
        f"1\t200\tlast\t{source_seq[199]}\t{source_seq[199]}C\t50\tPASS\t.\tGT\t0/1",
    ]
    vcf.write_text("\n".join(lines) + "\n")
    lifted = liftover(resolve_ref(load(vcf), store, collection=source), chain, store, target=target)
    report = by_case(lifted.report.collect())
    assert report["last"]["lift_reason"] == "indel_minus_strand"
    assert report["inner"]["lift_reason"] is None and lifted.reasons == {"indel_minus_strand": 1}
    inner = lifted.loaded.calls.collect().row(0, named=True)
    # Source 100 maps to target 101; its insertion is anchored on target 100, the image of source 101.
    assert (inner["pos"], inner["ref"], inner["alts"]) == (100, comp(source_seq[100]), [comp(source_seq[100]) + "G"])


def test_the_derived_facet(tmp_path: Path, setup) -> None:
    store, source, target, chain = setup
    resolved, lifted = lifted_from_cases(tmp_path, store, chain)
    parent = resolved.facets.collect().row(0, named=True)
    facet = lifted.loaded.facets.collect()
    assert facet.schema == FACETS
    facet = facet.row(0, named=True)
    assert facet["relation"] == RELATION == "liftover_to"
    assert facet["derived_from"] == parent["facet_id"] != facet["facet_id"]
    assert facet["ref_seqcol_digest"] == target and parent["ref_seqcol_digest"] == source
    assert (facet["assembly_label"], facet["assembly_label_basis"]) == ("GRCh38", "inferred")
    sha = hashlib.sha256(chain.read_bytes()).hexdigest()
    assert facet["derivation_inputs"] == [{"role": "chain", "name": chain.name, "sha256": sha}]
    assert facet["absence_means"] == "unknown"
    assert facet["header_contigs"] is None and facet["seqcol_digest"] is None
    assert {entry["chrom"] for entry in facet["contig_map"]} == set(lifted.loaded.calls.collect()["chrom"])
    assert lifted.loaded.calls.collect().schema == CALLS
    assert lifted.report.collect_schema() == REPORT
    # The parent is untouched, and every parent row is in the report once, in the parent's order.
    assert lifted.report.collect().select(CALLS.names()).equals(resolved.calls.collect())
    # Determinism: the same inputs give the same facet and rows.
    again = liftover(resolved, chain, store)
    assert again.loaded.facets.collect().equals(lifted.loaded.facets.collect())
    assert again.loaded.calls.collect().equals(lifted.loaded.calls.collect())
    assert any("not lifted" in w and "[multiple_targets]" in w for w in lifted.loaded.warnings)


def test_a_lifted_facet_names_no_reference_and_keeps_the_aligner(tmp_path: Path, setup) -> None:
    """The parent's `reference` named the source assembly's file; the lifted rows are against the
    target, which `ref_seqcol_digest` names. The reads were still aligned by the parent's aligner."""
    store, _, _, chain = setup
    dragen = (
        '##DRAGENCommandLine=<ID=dragen,Version="SW: 07.021.624.3.10.4, HW: 07.021.624",Date="Wed Apr 01 09:02:15 UTC 2020",'
        'CommandLineOptions="--ref-dir /staging/ref/hg19 --enable-map-align true --enable-variant-caller true">'
    )
    vcf = tmp_path / "source.vcf"
    vcf.write_text("\n".join([*HEADER[:1], "##reference=file:///staging/ref/hg19/reference.bin", dragen, *HEADER[1:], CASES["plus_snv"]]) + "\n")
    resolved = resolve_ref(load(vcf), store)
    parent = resolved.facets.collect().row(0, named=True)
    facet = liftover(resolved, chain, store).loaded.facets.collect().row(0, named=True)
    assert parent["reference"] == "hg19" and parent["aligner"] == "DRAGEN"
    assert facet["reference"] is None
    for name in ("aligner", "aligner_version", "method_steps"):
        assert facet[name] == parent[name], name


def test_lifted_rows_mint_on_the_target(tmp_path: Path, setup) -> None:
    store, _, _, chain = setup
    _, lifted = lifted_from_cases(tmp_path, store, chain)
    reasons = vrs_withheld(lifted.loaded.calls).collect()["vrs_withheld"].explode(empty_as_null=False).drop_nulls().unique().to_list()
    assert reasons == ["not_run"]
    minted = mint_vrs_ids(lifted.loaded, store).calls.collect()
    assert minted["vrs_allele_ids"].list.drop_nulls().list.len().to_list() == (minted["alts"].list.len() + 1).to_list()


def test_positions_agree_with_pyliftover(tmp_path: Path, setup) -> None:
    pyliftover = pytest.importorskip("pyliftover")
    _, _, _, chain = setup
    lifter = pyliftover.LiftOver(str(chain))
    plain = read_chain(chain)
    for pos in [*range(101, 401), *range(411, 601), *range(701, 801)]:
        # pyliftover takes 0-based positions and names; ambiguous ones give two hits.
        hits = lifter.convert_coordinate("1", pos - 1)
        if 151 <= pos <= 170:
            assert len(hits) == 2
            continue
        assert len(hits) == 1
        chrom, zero, strand, _ = hits[0]
        assert (chrom, zero + 1, strand) == expected_target(pos)
    assert plain.blocks.height == 5  # two blocks in the first chain, one in each other

    # A re-anchored minus-strand indel starts at the image of the source base after its span, where
    # pyliftover puts that base (the edge case's base is past the block, so pyliftover has none).
    store, _, _, _ = setup
    _, lifted = lifted_from_cases(tmp_path, store, chain)
    report = by_case(lifted.report.collect())
    for case in REANCHORED - {"minus_edge"}:
        after = report[case]["pos"] - 1 + len(report[case]["ref"])
        [(chrom, zero, strand, _)] = lifter.convert_coordinate("1", after)
        assert (report[case]["lift_chrom"], report[case]["lift_pos"], strand) == (chrom, zero + 1, "-"), case
    assert lifter.convert_coordinate("1", 800) == []


def test_segments_count_targets_not_blocks() -> None:
    # Block 1 overlaps block 0 with another target; block 3 repeats block 2's mapping over part of it.
    blocks = pl.DataFrame(
        {
            "block": [0, 1, 2, 3],
            "t_start": [0, 10, 30, 35],
            "t_end": [20, 15, 40, 38],
            "q_chrom": ["1", "2", "1", "1"],
            "q_minus": [False, False, False, False],
            "q_start": [100, 0, 200, 205],
        }
    ).with_columns(pl.col("block").cast(pl.UInt32))
    segments = _segments(blocks)
    assert segments.select("seg_start", "depth", "seg_block").rows() == [
        (0, 1, 0),
        (10, 2, None),
        (15, 1, 0),
        (20, 0, None),
        (30, 1, 2),
        (40, 0, None),
    ]


def test_a_chain_listing_its_blocks_twice_lifts_as_once(tmp_path: Path, setup) -> None:
    """Ensembl's GRCh38_to_GRCh37 chain repeats every chain under a second id (F14)."""
    store, _, _, chain = setup
    _, once = lifted_from_cases(tmp_path, store, chain)
    doubled = tmp_path / "doubled.chain"
    doubled.write_text(chain_text() + chain_text())
    _, twice = lifted_from_cases(tmp_path, store, doubled)
    assert twice.reasons == once.reasons
    assert twice.loaded.calls.drop("facet_id").collect().equals(once.loaded.calls.drop("facet_id").collect())


def test_array_rows_rederive_alts_on_a_minus_block(tmp_path: Path, setup) -> None:
    store, _, _, chain = setup
    path = tmp_path / "genome.txt"
    het = (751, s1(751) + other(s1(751)))
    hom_alt = (752, other(s1(752)) * 2)
    swap = (SWAP_POS, s1(SWAP_POS) + SWAP_BASE)
    lines = "".join(f"rs{p}\t1\t{p}\t{gt}\n" for p, gt in sorted([het, hom_alt, swap]))
    path.write_text(
        "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n"
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
        "# rsid\tchromosome\tposition\tgenotype\n" + lines
    )
    lifted = liftover(resolve_ref(load(path), store), chain, store)
    calls = {r["ids"][0]: r for r in lifted.loaded.calls.collect().iter_rows(named=True)}
    for rsid, r in calls.items():
        assert r["alts"] == sorted({a for a in r["gt_alleles"] if a != r["ref"]})
        assert [([r["ref"], *r["alts"]])[g] for g in r["gt"]] == r["gt_alleles"]
    assert calls[f"rs{SWAP_POS}"]["ref"] == SWAP_BASE
    assert sorted(calls[f"rs{het[0]}"]["gt_alleles"]) == sorted(revcomp(a) for a in het[1])


def test_a_chain_for_another_assembly_is_refused(tmp_path: Path, setup) -> None:
    store, _, _, _ = setup
    wrong = tmp_path / "wrong.chain"
    wrong.write_text(chain_text(source_1_size=len(S1) - 1))
    vcf = tmp_path / "source.vcf"
    vcf.write_text("\n".join([*HEADER, CASES["plus_snv"]]) + "\n")
    resolved = resolve_ref(load(vcf), store)
    with pytest.raises(ReferenceMismatchError, match="another assembly"):
        liftover(resolved, wrong, store)


def test_refuses_a_facet_whose_ref_was_not_read(tmp_path: Path, setup) -> None:
    store, _, _, chain = setup
    vcf = tmp_path / "source.vcf"
    vcf.write_text("\n".join([*HEADER, CASES["plus_snv"]]) + "\n")
    with pytest.raises(ValueError, match="resolve_ref first"):
        liftover(load(vcf), chain, store)


def test_normalize_left_aligns_lifted_indels_on_the_target(tmp_path: Path) -> None:
    """A flip chain reverses the source: an indel left-aligned in a repeat on the source lands at the
    repeat's right end on the target. `normalize` on the lifted facet moves it to the left end, read
    from the target's own sequence; bcftools norm on the target agrees."""
    local = random.Random(23)
    flank, middle, tail = ("".join(local.choice("ACGT") for _ in range(n)) for n in (80, 40, 80))
    # A T run of six between G and C, and four CA copies between T and G.
    source_seq = flank + "G" + "T" * 6 + "C" + middle + "T" + "CA" * 4 + "G" + tail
    target_seq = revcomp(source_seq)
    g = len(flank) + 1  # 1-based position of the G before the T run
    t = g + 8 + len(middle)  # 1-based position of the T before the CA copies
    store = ReferenceStore.open(tmp_path / "store")
    source = store.add_fasta(write_fasta(tmp_path / "s.fa", {"1": source_seq}), assembly_label="GRCh37")
    target_fasta = write_fasta(tmp_path / "t.fa", {"1": target_seq})
    target = store.add_fasta(target_fasta, assembly_label="GRCh38")
    chain = tmp_path / "flip.chain"
    size = len(source_seq)
    chain.write_text(f"chain 100 1 {size} + 0 {size} 1 {size} - 0 {size} 1\n{size}\n\n")
    vcf = tmp_path / "repeats.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={size}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        f"1\t{g}\tdel\tGT\tG\t50\tPASS\t.\tGT\t0/1",
        f"1\t{t}\tins\tT\tTCA\t50\tPASS\t.\tGT\t1/1",
    ]
    vcf.write_text("\n".join(lines) + "\n")
    lifted = liftover(resolve_ref(load(vcf), store, collection=source), chain, store, target=target)
    # On the target the T run reads as an A run between G and C, the CA copies as TG copies between C and A.
    a_run, tg_run = target_seq.find("G" + "A" * 6 + "C"), target_seq.find("C" + "TG" * 4 + "A")
    assert a_run == target_seq.rfind("G" + "A" * 6 + "C") >= 0 and tg_run == target_seq.rfind("C" + "TG" * 4 + "A") >= 0
    leftmost = {"del": (a_run + 1, "GA", ["G"]), "ins": (tg_run + 1, "C", ["CTG"])}
    as_lifted = {r["ids"][0]: r["pos"] for r in lifted.loaded.calls.collect().iter_rows(named=True)}
    assert all(as_lifted[case] > leftmost[case][0] for case in leftmost)  # at the run's right end after the lift

    normalized = normalize(lifted.loaded, store)
    rows = {r["ids"][0]: r for r in normalized.calls.collect().iter_rows(named=True)}
    assert {case: (r["pos"], r["ref"], r["alts"]) for case, r in rows.items()} == leftmost
    assert rows["ins"]["gt_alleles"] == ["CTG", "CTG"] and sorted(rows["del"]["gt_alleles"]) == ["G", "GA"]
    assert not any("not left-aligned or trimmed" in w for w in normalized.warnings)
    assert any("2 row(s) left-aligned to another position" in w for w in normalized.warnings)

    pysam = pytest.importorskip("pysam")
    bcftools = pytest.importorskip("pysam.bcftools")
    pysam.faidx(str(target_fasta))
    exported = write_vcf(lifted.loaded, tmp_path / "lifted.vcf", compress=False)
    theirs = tmp_path / "bcftools.vcf"
    bcftools.norm("-f", str(target_fasta), "--no-version", "-o", str(theirs), str(exported.path), catch_stdout=False)
    by_bcftools = {r["ids"][0]: (r["pos"], r["source_ref"], r["alts"]) for r in load(theirs).calls.collect().iter_rows(named=True)}
    assert by_bcftools == leftmost


# Blocks (gVCF reference and no-call intervals) on source contig 1, each named by where it lies in the
# construction above; and one on contig 2, which the chain does not cover.
GVCF_BLOCKS = {
    "whole": ("1", 201, 240, "0/0:30:12:0,30,450"),  # inside block A
    "gap": ("1", 381, 430, "0/0:40:15:0,40,600"),  # across the deleted source 401..410, A then B
    "ambiguous": ("1", 141, 180, "0/0:20:9:0,20,300"),  # over source 151..170, which a second chain also maps
    "minus": ("1", 711, 740, "0/0:50:20:0,50,800"),  # on the minus-strand chain to target 2
    "unmapped": ("1", 21, 60, "0/0:10:4:0,10,100"),  # before the first chain block
    "unknown_target": ("1", 851, 858, "0/0:11:4:0,11,100"),  # on the chain to chrUn_q, which the target lacks
    "no_call": ("1", 581, 620, "./.:.:.:."),  # the end of block B, then no chain
    "other_contig": ("2", 1, 50, "0/0:60:23:0,60,900"),
}
GVCF_HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SOURCE.items()],
    '##ALT=<ID=NON_REF,Description="Any other allele">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Block minimum depth">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def base_fate(chrom: str, pos: int) -> tuple[str, int, str] | str:
    """Where the construction sends one source base: (chrom, 1-based pos, strand), or why nowhere."""
    i = pos - 1
    if chrom != "1":
        return "contig_not_in_chain"
    if 150 <= i < 170:
        return "multiple_targets"
    if 850 <= i < 860:
        return "target_contig_unmatched"
    if 100 <= i < 400 or 410 <= i < 600 or 700 <= i < 800:
        return expected_target(pos)
    return "no_chain_block"


def expected_pieces(chrom: str, start: int, end: int) -> tuple[list[tuple[str, int, int]], Counter]:
    """A block's lifted pieces, base by base from the construction: runs of source bases that are
    consecutive on the source and on the target alike; and the bases left behind per reason."""
    pieces: list[list] = []
    behind: Counter = Counter()
    previous = None
    for pos in range(start, end + 1):
        fate = base_fate(chrom, pos)
        if isinstance(fate, str):
            behind[fate] += 1
            previous = None
            continue
        target_chrom, target_pos, strand = fate
        step = 1 if strand == "+" else -1
        if previous is not None and (previous[0], previous[1] + step, previous[2]) == fate:
            pieces[-1][2 if strand == "+" else 1] = target_pos
        else:
            pieces.append([target_chrom, target_pos, target_pos])
        previous = fate
    return [tuple(piece) for piece in pieces], behind


def gvcf_lifted(tmp_path: Path, store: ReferenceStore, chain: Path):
    records = [
        f"{chrom}\t{start}\t{name}\t{SOURCE[chrom][start - 1]}\t<NON_REF>\t.\t.\tEND={end}\tGT:GQ:MIN_DP:PL\t{sample}"
        for name, (chrom, start, end, sample) in GVCF_BLOCKS.items()
    ]
    call = f"1\t250\tcall\t{s1(250)}\t{other(s1(250))},<NON_REF>\t40\tPASS\t.\tGT:GQ:PL\t0/1:40:40,0,50,90,95,99"
    lines = sorted([*records, call], key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    vcf = tmp_path / "source.g.vcf"
    vcf.write_text("\n".join([*GVCF_HEADER, *lines]) + "\n")
    resolved = resolve_ref(load(vcf), store)
    return resolved, liftover(resolved, chain, store)


def block_key(r: dict) -> tuple:
    return (r["chrom"], r["pos"], r["end_pos"], r["kind"], r["gq"], r["dp_min"], tuple(r["gl"] or ()))


def test_blocks_lift_cut_at_every_chain_gap_and_edge(tmp_path: Path, setup) -> None:
    """Each block's lifted pieces are the runs the construction maps base by base; nothing is
    stretched over a gap or merged across one, and what stays behind is counted in bases per reason."""
    store, _, _, chain = setup
    resolved, lifted = gvcf_lifted(tmp_path, store, chain)
    source = resolved.blocks.collect()
    assert source.height == len(GVCF_BLOCKS)
    blocks = lifted.loaded.blocks.collect()
    assert blocks.schema == BLOCKS

    expected_rows = []
    behind: Counter = Counter()
    for evidence in source.iter_rows(named=True):
        pieces, left = expected_pieces(evidence["chrom"], evidence["pos"], evidence["end_pos"])
        behind += left
        # Every piece carries its block's evidence unchanged.
        expected_rows += [block_key({**evidence, "chrom": chrom, "pos": pos, "end_pos": end}) for chrom, pos, end in pieces]
    got = [block_key(r) for r in blocks.iter_rows(named=True)]
    assert sorted(got) == sorted(expected_rows)
    assert lifted.block_reasons == dict(sorted(behind.items()))
    assert set(lifted.block_reasons) <= UNLIFTED_REASONS
    # Every source base is either in a lifted block or counted behind.
    lifted_bases = sum(end - pos + 1 for _, pos, end, *_ in got)
    assert lifted_bases + sum(lifted.block_reasons.values()) == sum(r["end_pos"] - r["pos"] + 1 for r in source.iter_rows(named=True))

    # The block across the deleted source bases lands as two blocks that touch on the target: kept
    # apart, since the source between them is not there.
    assert [(r[1], r[2]) for r in got if r[4] == 40] == [(331, 350), (351, 370)]
    # The minus-strand block lands on target 2 reversed, as one block.
    assert [r[:3] for r in got if r[4] == 50] == [("2", expected_target(740)[1], expected_target(711)[1])]
    # Canonical order on the target, and the facet's contig map names the contigs the blocks use.
    assert blocks.equals(canonical_order(blocks.lazy()).collect())
    facet = lifted.loaded.facets.collect().row(0, named=True)
    assert set(blocks["chrom"]) <= {entry["chrom"] for entry in facet["contig_map"]}
    assert facet["absence_means"] == "unknown"
    assert set(blocks["facet_id"]) == {facet["facet_id"]}
    assert any("block(s):" in w and "left behind" in w for w in lifted.loaded.warnings)
    assert any("[multiple_targets]" in w and "block base(s)" in w for w in lifted.loaded.warnings)

    # Determinism, and a reused generation answers as the run that wrote it.
    again = liftover(resolved, chain, store)
    assert again.loaded.blocks.collect().equals(blocks)
    assert again.block_reasons == lifted.block_reasons


def test_lifted_blocks_take_their_ref_from_the_target_on_export(tmp_path: Path, setup) -> None:
    """The frame keeps no block REF: a lifted gVCF written out reads each block's first base from the
    target collection, and `load` reads the file back to the lifted blocks."""
    store, _, _, chain = setup
    _, lifted = gvcf_lifted(tmp_path, store, chain)
    written = write_vcf(lifted.loaded, tmp_path / "lifted.g.vcf", blocks=True, reference=store, compress=False)
    records = [line.split("\t") for line in written.path.read_text().splitlines() if not line.startswith("#")]
    block_records = [fields for fields in records if fields[4] == "<*>"]
    assert len(block_records) == lifted.loaded.blocks.select(pl.len()).collect().item() > 0
    for fields in block_records:
        assert fields[3] == TARGET[fields[0]][int(fields[1]) - 1], fields[:2]
    reloaded = load(written.path)
    assert reloaded.blocks.drop("facet_id").collect().equals(lifted.loaded.blocks.drop("facet_id").collect())


SV_HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SOURCE.items()],
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CI around POS">',
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CI around END">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def sv(name: str, pos: int, alt: str, end: int | None = None, info: str = "", sample: str = "0/1", ref: str | None = None) -> str:
    fields = [f"END={end}"] if end is not None else []
    fields += [info] if info else []
    fmt, value = ("GT:CN", sample) if ":" in sample else ("GT", sample)
    return f"1\t{pos}\t{name}\t{ref or s1(pos)}\t{alt}\t50\tPASS\t{';'.join(fields) or '.'}\t{fmt}\t{value}"


# Symbolic SVs and copy-number segments on source contig 1, by where their two ends fall in the
# construction (block A: source 101..400 at -50; source 401..410 deleted; block B: 411..600 at -60;
# 701..800 reversed onto target 2; 151..170 also mapped by a second chain).
SV_CASES = {
    "del_plus": (sv("del_plus", 201, "<DEL>", 230, "SVTYPE=DEL;CIPOS=-5,10;CIEND=-3,4"), None),
    # Across the 10 deleted bases: 190 against 200, exactly SV_LENGTH_TOLERANCE, kept.
    "dup_across_gap": (sv("dup_across_gap", 211, "<DUP:TANDEM>", 411), None),
    # Across them too: 20 against 30, too far.
    "del_shrunk": (sv("del_shrunk", 390, "<DEL>", 420), "sv_length_changed"),
    "inv_end_unmapped": (sv("inv_end_unmapped", 550, "<INV>", 650), "sv_end_unlifted"),
    "cnv_apart": (sv("cnv_apart", 500, "<CNV>", 750), "sv_ends_apart"),
    "del_minus": (sv("del_minus", 720, "<DEL>", 760, "CIPOS=-6,2;CIEND=-1,7"), None),
    "ins_plus": (sv("ins_plus", 300, "<INS>", info="SVTYPE=INS"), None),
    "ins_minus": (sv("ins_minus", 750, "<INS>"), None),
    "segment": (sv("segment", 460, ".", 520, sample="0/0:2"), None),
    "segment_minus": (sv("segment_minus", 705, ".", 795, "CIPOS=-4,6", sample="0/0:3"), None),
    # The target changes the padding base at SWAP_POS: a plus-strand anchor must be the source's.
    "del_anchor_changed": (sv("del_anchor_changed", SWAP_POS, "<DEL>", 280), "ref_mismatch_after_lift"),
    "del_ambiguous": (sv("del_ambiguous", 155, "<DEL>", 190), "multiple_targets"),
    "del_beside_snv": (sv("del_beside_snv", 330, f"{other(s1(330))},<DEL>", 360), "allele_not_sequence"),
    "del_long_ref": (sv("del_long_ref", 340, "<DEL>", 360, ref=s1(340, 2)), "allele_not_sequence"),
}


def expected_sv(pos: int, end: int, padded: bool) -> tuple[str, int, int, str]:
    """Where the construction puts an SV lifted by both ends: (chrom, pos, end_pos, strand). On the
    minus strand the event's image is reversed; a padded SV is anchored on the base before it."""
    (chrom, a, strand), (_, b, _) = expected_target(pos), expected_target(end)
    if strand == "+":
        return chrom, a, b, strand
    return (chrom, b - 1, a - 1, strand) if padded else (chrom, b, a, strand)


def mirrored(pair: list[int]) -> list[int]:
    return [-pair[1], -pair[0]]


def sv_lifted(tmp_path: Path, store: ReferenceStore, chain: Path):
    lines = sorted((line for line, _ in SV_CASES.values()), key=lambda line: int(line.split("\t")[1]))
    vcf = tmp_path / "sv.vcf"
    vcf.write_text("\n".join([*SV_HEADER, *lines]) + "\n")
    resolved = resolve_ref(load(vcf), store)
    return resolved, liftover(resolved, chain, store)


def test_symbolic_svs_lift_by_both_ends(tmp_path: Path, setup) -> None:
    store, _, _, chain = setup
    resolved, lifted = sv_lifted(tmp_path, store, chain)
    source = by_case(resolved.calls.collect())
    report = by_case(lifted.report.collect())
    assert {case: r["lift_reason"] for case, r in report.items()} == {case: reason for case, (_, reason) in SV_CASES.items()}
    assert set(lifted.reasons) <= UNLIFTED_REASONS

    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == {case for case, (_, reason) in SV_CASES.items() if reason is None}
    for case, r in calls.items():
        parent = source[case]
        chrom, pos, end, strand = expected_sv(parent["pos"], parent["end_pos"], padded=bool(parent["alts"]))
        assert (r["chrom"], r["pos"], r["end_pos"]) == (chrom, pos, end), case
        assert (report[case]["lift_chrom"], report[case]["lift_pos"], report[case]["lift_strand"]) == (chrom, pos, strand)
        # The length is the source's wherever both ends sit in one chain block.
        if case != "dup_across_gap":
            assert r["end_pos"] - r["pos"] == parent["end_pos"] - parent["pos"], case
        # REF is the target's base at the new POS; the ALTs stay as written; gt indexes [ref] + alts.
        assert r["ref"] == TARGET[chrom][pos - 1] and r["source_ref"] is None
        assert r["alts"] == parent["alts"]
        assert [[r["ref"], *r["alts"]][g] for g in r["gt"]] == r["gt_alleles"]
        assert (r["sv_type"], r["copy_number"]) == (parent["sv_type"], parent["copy_number"])
    assert calls["dup_across_gap"]["end_pos"] - calls["dup_across_gap"]["pos"] == 190

    # CIPOS and CIEND carry on the plus strand; on the minus strand the new POS is the old END's
    # image, so they swap and mirror, and a missing CIEND reads as CIPOS (VCF 4.5).
    assert (calls["del_plus"]["cipos"], calls["del_plus"]["ciend"]) == (source["del_plus"]["cipos"], source["del_plus"]["ciend"])
    assert (calls["del_minus"]["cipos"], calls["del_minus"]["ciend"]) == (mirrored(source["del_minus"]["ciend"]), mirrored(source["del_minus"]["cipos"]))
    assert source["segment_minus"]["ciend"] is None
    assert calls["segment_minus"]["cipos"] == calls["segment_minus"]["ciend"] == mirrored(source["segment_minus"]["cipos"])
    assert any("symbolic SVs or copy-number segments" in w for w in lifted.loaded.warnings)


def test_lifted_svs_round_trip_through_the_writer(tmp_path: Path, setup) -> None:
    """A lifted SV's END is its lifted `end_pos`: written out and loaded back, the extents are the same."""
    store, _, _, chain = setup
    _, lifted = sv_lifted(tmp_path, store, chain)
    written = write_vcf(lifted.loaded, tmp_path / "lifted.vcf", compress=False, resolved_only=True)
    back = by_case(load(written.path).calls.collect())
    calls = by_case(lifted.loaded.calls.collect())
    assert {case for case, (line, reason) in SV_CASES.items() if reason is None and "<" in line} <= set(back)
    for case in back:
        assert (back[case]["chrom"], back[case]["pos"], back[case]["end_pos"]) == (calls[case]["chrom"], calls[case]["pos"], calls[case]["end_pos"])


def one_contig_lift(tmp_path: Path, source_seq: str, target_seq: str, chain: str, records: list[str]):
    store = ReferenceStore.open(tmp_path / "store")
    source = store.add_fasta(write_fasta(tmp_path / "s.fa", {"1": source_seq}), assembly_label="GRCh37")
    target = store.add_fasta(write_fasta(tmp_path / "t.fa", {"1": target_seq}), assembly_label="GRCh38")
    path = tmp_path / "one.chain"
    path.write_text(chain)
    vcf = tmp_path / "one.vcf"
    header = [line for line in SV_HEADER if not line.startswith("##contig")]
    vcf.write_text("\n".join([header[0], f"##contig=<ID=1,length={len(source_seq)}>", *header[1:], *records]) + "\n")
    return liftover(resolve_ref(load(vcf), store, collection=source), path, store, target=target)


def test_a_minus_strand_sv_keeps_a_precise_end_precise(tmp_path: Path) -> None:
    """On a minus-strand block the new POS is the old END's image and the new END the old POS's, so the
    confidence intervals swap and mirror. A source with CIEND and no CIPOS has an imprecise END and a
    precise POS; lifted, the new END is that precise POS. Writing it null made VCF 4.5's fallback (a
    missing CIEND reads as CIPOS) give the new END the imprecise interval. It is written `[0, 0]`."""
    source_seq = random_bases(200)
    chain = "chain 100 1 200 + 0 200 1 200 - 0 200 1\n200\n\n"
    end_only = f"1\t100\tend_only\t{source_seq[99]}\t<DEL>\t50\tPASS\tEND=150;CIEND=-5,7\tGT\t0/1"
    pos_only = f"1\t110\tpos_only\t{source_seq[109]}\t<DEL>\t50\tPASS\tEND=160;CIPOS=-3,4\tGT\t0/1"
    precise = f"1\t120\tprecise\t{source_seq[119]}\t<DEL>\t50\tPASS\tEND=170\tGT\t0/1"
    lifted = one_contig_lift(tmp_path, source_seq, revcomp(source_seq), chain, [end_only, pos_only, precise])
    calls = by_case(lifted.loaded.calls.collect())
    assert (calls["end_only"]["cipos"], calls["end_only"]["ciend"]) == (mirrored([-5, 7]), [0, 0])
    # A missing CIEND is CIPOS at the source, so both lifted ends get its mirror; nothing stated stays nothing.
    assert calls["pos_only"]["cipos"] == calls["pos_only"]["ciend"] == mirrored([-3, 4])
    assert (calls["precise"]["cipos"], calls["precise"]["ciend"]) == (None, None)


def test_an_sv_whose_ends_swap_order_stays_behind(tmp_path: Path) -> None:
    """Two plus-strand chains swap the source's halves: an SV across the middle has its END land before its POS."""
    source_seq = random_bases(300)
    target_seq = source_seq[150:] + source_seq[:150]
    chain = "chain 100 1 300 + 0 150 1 300 + 150 300 1\n150\n\nchain 90 1 300 + 150 300 1 300 + 0 150 2\n150\n\n"
    inside = f"1\t20\tinside\t{source_seq[19]}\t<DEL>\t50\tPASS\tEND=60\tGT\t0/1"
    across = f"1\t100\tacross\t{source_seq[99]}\t<DEL>\t50\tPASS\tEND=200\tGT\t0/1"
    lifted = one_contig_lift(tmp_path, source_seq, target_seq, chain, [inside, across])
    report = by_case(lifted.report.collect())
    assert report["across"]["lift_reason"] == "sv_ends_reversed"
    assert (report["inside"]["lift_reason"], report["inside"]["lift_pos"]) == (None, 170)


def test_a_minus_strand_sv_at_the_target_contigs_first_base_stays_behind(tmp_path: Path) -> None:
    """Target 1 is source 1 reverse-complemented whole: an insertion after the last source base would be
    anchored before the target's first base. A deletion inside is anchored on the base before its image."""
    source_seq = random_bases(200)
    chain = "chain 100 1 200 + 0 200 1 200 - 0 200 1\n200\n\n"
    last = f"1\t200\tlast\t{source_seq[199]}\t<INS>\t50\tPASS\t.\tGT\t0/1"
    inner = f"1\t100\tinner\t{source_seq[99]}\t<DEL>\t50\tPASS\tEND=150\tGT\t0/1"
    lifted = one_contig_lift(tmp_path, source_seq, revcomp(source_seq), chain, [inner, last])
    report = by_case(lifted.report.collect())
    assert report["last"]["lift_reason"] == "indel_minus_strand"
    inner_row = by_case(lifted.loaded.calls.collect())["inner"]
    # Source 101..150 deleted: their images are 51..100, so the deletion is anchored on 50 and ends at 100.
    assert (inner_row["pos"], inner_row["end_pos"], inner_row["ref"]) == (50, 100, revcomp(source_seq)[49])
