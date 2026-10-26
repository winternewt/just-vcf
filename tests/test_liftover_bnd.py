"""RM6: `liftover` lifts breakends: the record's own position and its mate's, each through the chain.

A synthetic source assembly (contigs `1`, `2`, `3`), a target built from it by known moves (contigs
`chr1`, `chr2`, `chr3`, so a lifted ALT must name the mate as the target collection does), and the chain
describing those moves, written both ways so a lift can be undone. Expected positions come from the
construction's own mapping (`image`), the expected bracket forms from the strand table written out
below by hand, and every lifted record is checked a third way: the novel sequence the ALT asserts, read
from the source by VCF 4.3 §5.4's four sentences, must be the sequence the lifted ALT asserts on the
target, read either way along the molecule.
"""

import gzip
import json
import random
import re
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ContigNames, ReferenceStore, liftover, load, resolve_ref, write_vcf
from just_vcf.cache import cache_dir
from just_vcf.liftover import UNLIFTED_REASONS
from just_vcf.schema import canonical_chrom

rng = random.Random(11)


def random_bases(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT[]", "TGCA]["))


def comp(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))


def other(*bases: str) -> str:
    return next(b for b in "ACGT" if b not in bases)


S1, S2, S3 = random_bases(1000), random_bases(600), random_bases(300)
SOURCE = {"1": S1, "2": S2, "3": S3}

# Blocks: (source contig, 0-based source start, end, target contig, 0-based forward target start, strand).
# Source 3 is in no block, as are source 1's 0..20, 80..100, 400..500, 800..850, 950..1000 and 2's 500..600.
BLOCKS = [
    ("1", 100, 400, "chr1", 50, "+"),  # P1
    ("2", 0, 300, "chr1", 400, "+"),  # P2: a second source contig onto the same target contig
    ("1", 850, 950, "chr1", 800, "-"),  # M3: source 1 onto chr1 again, reversed
    ("1", 20, 80, "chr1", 920, "+"),  # P3: source 1's head after the rest, so order flips
    ("1", 500, 800, "chr2", 30, "-"),  # M1
    ("2", 300, 500, "chr3", 20, "-"),  # M2
]
FILL = {name: random_bases(n) for name, n in (("a", 50), ("b", 50), ("c", 100), ("d", 20), ("e", 20), ("f", 30), ("g", 30), ("h", 20), ("i", 30))}
T1 = FILL["a"] + S1[100:400] + FILL["b"] + S2[0:300] + FILL["c"] + revcomp(S1[850:950]) + FILL["d"] + S1[20:80] + FILL["e"]
T2 = FILL["f"] + revcomp(S1[500:800]) + FILL["g"]
T3 = FILL["h"] + revcomp(S2[300:500]) + FILL["i"]
# One base the target changes: the image of source 1:251 (target chr1:201).
EDIT_SOURCE, EDIT_TARGET = 251, 201
T1 = T1[: EDIT_TARGET - 1] + other(T1[EDIT_TARGET - 1]) + T1[EDIT_TARGET:]
TARGET = {"chr1": T1, "chr2": T2, "chr3": T3}


def image(chrom: str, pos: int) -> tuple[str, int, str] | None:
    """The construction's mapping of a 1-based source base: (target contig, 1-based pos, strand)."""
    i = pos - 1
    for s_chrom, start, end, t_chrom, t_start, strand in BLOCKS:
        if s_chrom == chrom and start <= i < end:
            return t_chrom, (t_start + (i - start) if strand == "+" else t_start + (end - 1 - i)) + 1, strand
    return None


def chain_text(inverse: bool = False) -> str:
    """The chain of `BLOCKS`, one chain per block; `inverse` writes target to source."""
    lines = []
    for index, (s_chrom, start, end, t_chrom, t_start, strand) in enumerate(BLOCKS):
        length = end - start
        s_size, t_size = len(SOURCE[s_chrom]), len(TARGET[t_chrom])
        if not inverse:
            q0 = t_start if strand == "+" else t_size - (t_start + length)
            lines.append(f"chain 100 {s_chrom} {s_size} + {start} {end} {t_chrom} {t_size} {strand} {q0} {q0 + length} {index + 1}")
        else:
            q0 = start if strand == "+" else s_size - end
            lines.append(f"chain 100 {t_chrom} {t_size} + {t_start} {t_start + length} {s_chrom} {s_size} {strand} {q0} {q0 + length} {index + 1}")
        lines += [str(length), ""]
    return "\n".join(lines)


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()))
    return path


HEADER = [
    "##fileformat=VCFv4.3",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SOURCE.items()],
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="CI around POS">',
    '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="CI around the mate">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def base(chrom: str, pos: int) -> str:
    return SOURCE[chrom][pos - 1]


def bnd(case: str, chrom: str, pos: int, alt: str, gt: str = "0/1", info: str = "SVTYPE=BND") -> str:
    """A record; `{r}` in `alt` is the source's base at POS."""
    return f"{chrom}\t{pos}\t{case}\t{base(chrom, pos)}\t{alt.format(r=base(chrom, pos))}\t50\tPASS\t{info}\tGT\t{gt}"


def four(prefix: str, chrom: str, pos: int, mate: str, info: str = "SVTYPE=BND") -> dict[str, str]:
    """The four forms at one position with one mate: `_lr` is `t[p[` (own sequence left of POS, the
    mate's right of p), `_ll` `t]p]`, `_rl` `]p]t`, `_rr` `[p[t`."""
    alts = {"lr": "{r}[" + mate + "[", "ll": "{r}]" + mate + "]", "rl": "]" + mate + "]{r}", "rr": "[" + mate + "[{r}"}
    return {f"{prefix}_{key}": bnd(f"{prefix}_{key}", chrom, pos, alt, info=info) for key, alt in alts.items()}


CI = "SVTYPE=BND;CIPOS=-3,5;CIEND=-2,7"
CASES = {
    # Own end and mate on plus-strand blocks of one source contig, landing on one target contig.
    **four("pp", "1", 200, "1:300", info=CI),
    # Own end plus (P1 to chr1), mate minus on another source contig (M2 to chr3).
    **four("pm", "1", 200, "2:400", info=CI),
    # Own end minus (M1 to chr2), mate plus on another source contig (P2 to chr1).
    **four("mp", "1", 600, "2:100", info=CI),
    # Both on one minus block (M1): their order reverses on the target, as the strand does.
    **four("mm", "1", 600, "1:700", info=CI),
    # One source contig, the ends on two strands of one target contig: no order to keep, distance 700 -> 701.
    "strands": bnd("strands", "1", 200, "{r}[1:900["),
    # One source contig, ends on two target contigs (liftoverSV's case 2).
    "apart": bnd("apart", "1", 200, "{r}[1:600["),
    # One source contig, both ends plus on chr1 but the order flipped (liftoverSV's case 3).
    "reversed": bnd("reversed", "1", 50, "{r}[1:200["),
    # One source contig, ends on two strands, distance 820 -> 741 (liftoverSV's case 4).
    "distance": bnd("distance", "1", 120, "{r}[1:940["),
    # The mate in a chain gap, the mate on a contig the chain does not cover, the own end in a gap.
    "mate_gap": bnd("mate_gap", "1", 210, "{r}[1:450["),
    "mate_no_chain": bnd("mate_no_chain", "1", 215, "{r}]3:100]"),
    "own_gap": bnd("own_gap", "1", 450, "{r}[1:200["),
    # Single breakends: one position, the dot moving side on a minus block.
    "single_plus_r": bnd("single_plus_r", "1", 220, "{r}."),
    "single_plus_l": bnd("single_plus_l", "1", 225, ".{r}"),
    "single_minus_r": bnd("single_minus_r", "1", 610, "{r}."),
    "single_minus_l": bnd("single_minus_l", "1", 620, ".{r}"),
    # Inserted bases: reverse-complemented to the other side of REF when the own end flips; kept when only the mate does.
    "ins_minus": bnd("ins_minus", "1", 630, "{r}GTCA[2:120["),
    "ins_mate_minus": bnd("ins_mate_minus", "1", 230, "]2:350]AC{r}"),
    # Two breakend ALTs in one record, the second's mate on a minus block.
    "two_alts": bnd("two_alts", "1", 270, "{r}[1:300[,{r}]2:400]", gt="1/2"),
    # The target changed the base at POS.
    "ref_mismatch": bnd("ref_mismatch", "1", EDIT_SOURCE, "{r}[1:300["),
    # Brackets that do not pair, and a breakend beside a sequence ALT: still not lifted.
    "unparsed": bnd("unparsed", "1", 260, "{r}[1:300]"),
    "mixed": bnd("mixed", "1", 265, "{r}[1:300[," + other(base("1", 265))),
}

FORM_AFTER = {
    # (source form, own end strand, mate strand) -> lifted form, from VCF 4.3 §5.4 by hand: a minus
    # block at the own end moves t across the brackets (the sequence kept beside POS changes side);
    # a minus block at the mate turns the brackets (the piece runs the other way from p).
    ("lr", "+", "+"): "lr", ("lr", "+", "-"): "ll", ("lr", "-", "+"): "rr", ("lr", "-", "-"): "rl",
    ("ll", "+", "+"): "ll", ("ll", "+", "-"): "lr", ("ll", "-", "+"): "rl", ("ll", "-", "-"): "rr",
    ("rl", "+", "+"): "rl", ("rl", "+", "-"): "rr", ("rl", "-", "+"): "ll", ("rl", "-", "-"): "lr",
    ("rr", "+", "+"): "rr", ("rr", "+", "-"): "rl", ("rr", "-", "+"): "lr", ("rr", "-", "-"): "ll",
}  # fmt: skip

EXPECTED_REASON = {
    **{case: None for case in CASES},
    "apart": "breakend_ends_apart",
    "reversed": "breakend_ends_reversed",
    "distance": "breakend_distance_changed",
    "mate_gap": "breakend_mate_not_lifted",
    "mate_no_chain": "breakend_mate_not_lifted",
    "own_gap": "no_chain_block",
    "ref_mismatch": "ref_mismatch_after_lift",
    "unparsed": "breakend_not_lifted",
    "mixed": "breakend_not_lifted",
}

_PAIRED = re.compile(r"^(?:([A-Z]+)([\[\]])(.+):(\d+)([\[\]])|([\[\]])(.+):(\d+)([\[\]])([A-Z]+))$")


def parse(alt: str) -> tuple[str, str | None, int | None, str]:
    """(form, mate contig, mate pos, t) of a breakend ALT: forms `lr`, `ll`, `rl`, `rr`, or `r.`/`.l`
    for single breakends (`t.`, `.t`)."""
    if alt.endswith("."):
        return "r.", None, None, alt[:-1]
    if alt.startswith("."):
        return ".l", None, None, alt[1:]
    match = _PAIRED.match(alt)
    assert match, alt
    if match.group(1):
        t, bracket, contig, pos = match.group(1), match.group(2), match.group(3), int(match.group(4))
        return ("lr" if bracket == "[" else "ll"), contig, pos, t
    bracket, contig, pos, t = match.group(6), match.group(7), int(match.group(8)), match.group(10)
    return ("rr" if bracket == "[" else "rl"), contig, pos, t


def junction(seqs: dict[str, str], chrom: str, pos: int, alt: str, k: int = 6) -> str:
    """The novel sequence an ALT asserts around its join, read from `seqs` by VCF 4.3 §5.4: s is
    replaced by t, then `t[p[` joins the piece right of p after t, `t]p]` the reverse complement of
    the piece left of p after t, `]p]t` the piece left of p before t, `[p[t` the reverse complement of
    the piece right of p before t. A single breakend joins t to unknown sequence (`.`)."""
    form, mate_chrom, p, t = parse(alt)
    own = seqs[chrom]
    before, after = own[pos - 1 - k : pos - 1], own[pos : pos + k]
    if form == "r.":
        return before + t + "."
    if form == ".l":
        return "." + t + after
    mate = seqs[mate_chrom]
    right_of_p, left_of_p = mate[p - 1 : p - 1 + k], mate[p - k : p]
    return {
        "lr": before + t + right_of_p,
        "ll": before + t + revcomp(left_of_p),
        "rl": left_of_p + t + after,
        "rr": revcomp(right_of_p) + t + after,
    }[form]


@pytest.fixture
def setup(tmp_path: Path) -> tuple[ReferenceStore, str, str, Path, Path]:
    store = ReferenceStore.open(tmp_path / "store")
    source = store.add_fasta(write_fasta(tmp_path / "source.fa", SOURCE), assembly_label="GRCh37")
    target = store.add_fasta(write_fasta(tmp_path / "target.fa", TARGET), assembly_label="GRCh38")
    forward, back = tmp_path / "forward.chain", tmp_path / "back.chain"
    forward.write_text(chain_text())
    back.write_text(chain_text(inverse=True))
    return store, source, target, forward, back


def lift_cases(tmp_path: Path, store: ReferenceStore, source: str, chain: Path, cases: dict[str, str] = CASES):
    order = sorted(cases.values(), key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    vcf = tmp_path / "bnd.vcf"
    vcf.write_text("\n".join([*HEADER, *order]) + "\n")
    resolved = resolve_ref(load(vcf), store, collection=source)
    return resolved, liftover(resolved, chain, store)


def by_case(frame: pl.DataFrame) -> dict[str, dict]:
    return {r["ids"][0]: r for r in frame.iter_rows(named=True)}


def test_every_breakend_gets_its_reason(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    report = by_case(lifted.report.collect())
    assert {case: r["lift_reason"] for case, r in report.items()} == EXPECTED_REASON
    assert set(lifted.reasons) <= UNLIFTED_REASONS
    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == {case for case, reason in EXPECTED_REASON.items() if reason is None}


def test_each_form_on_each_strand_turns_as_the_spec_says(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())
    for prefix in ("pp", "pm", "mp", "mm"):
        for key in ("lr", "ll", "rl", "rr"):
            case = f"{prefix}_{key}"
            r, source_row = calls[case], report[case]
            _, mate_chrom, mate_pos, _ = parse(source_row["alts"][0])
            own = image("1", source_row["pos"])
            mate = image(mate_chrom, mate_pos)
            assert (r["chrom"], r["pos"], report[case]["lift_strand"]) == (own[0].removeprefix("chr"), own[1], own[2]), case
            assert r["end_pos"] == r["pos"]
            form, lifted_mate_chrom, lifted_mate_pos, t = parse(r["alts"][0])
            assert form == FORM_AFTER[(key, own[2], mate[2])], case
            # The mate is named as the target collection names it, at its image.
            assert (lifted_mate_chrom, lifted_mate_pos) == (mate[0], mate[1]), case
            # The anchor is the target's base, the source's complemented on a minus block.
            assert r["ref"] == TARGET[own[0]][own[1] - 1] == (comp(source_row["ref"]) if own[2] == "-" else source_row["ref"])
            assert t == r["ref"]
            assert r["gt_alleles"] == [r["ref"], r["alts"][0]] and r["gt"] == [0, 1]
            assert r["sv_type"] == source_row["sv_type"]


def test_a_lifted_breakend_asserts_the_same_novel_sequence(tmp_path: Path, setup) -> None:
    """Read through VCF's four sentences, the lifted ALT on the target and the source ALT on the
    source describe one molecule, read in one direction or the other."""
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())
    target_by_chrom = {name.removeprefix("chr"): seq for name, seq in TARGET.items()}
    names = {**TARGET, **target_by_chrom}
    checked = 0
    for case, r in calls.items():
        source_row = report[case]
        for source_alt, lifted_alt in zip(source_row["alts"], r["alts"], strict=True):
            before = junction(SOURCE, source_row["chrom"], source_row["pos"], source_alt)
            after = junction(names, r["chrom"], r["pos"], lifted_alt)
            assert after in (before, revcomp(before)), (case, source_alt, lifted_alt)
            checked += 1
    assert checked == sum(len(r["alts"]) for r in calls.values()) > 20


def test_singles_inserted_bases_and_two_alts(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = by_case(lifted.loaded.calls.collect())
    r = {case: calls[case] for case in calls}

    def at(chrom: str, pos: int) -> str:
        return TARGET[image(chrom, pos)[0]][image(chrom, pos)[1] - 1]

    assert r["single_plus_r"]["alts"] == [f"{base('1', 220)}."] and r["single_plus_r"]["pos"] == image("1", 220)[1]
    assert r["single_plus_l"]["alts"] == [f".{base('1', 225)}"]
    assert r["single_minus_r"]["alts"] == [f".{comp(base('1', 610))}"] and r["single_minus_r"]["pos"] == image("1", 610)[1]
    assert r["single_minus_l"]["alts"] == [f"{comp(base('1', 620))}."]
    # t = REF + GTCA on a minus own end becomes TGAC + REF', with the brackets on the other side.
    mate = image("2", 120)
    assert r["ins_minus"]["alts"] == [f"[{mate[0]}:{mate[1]}[TGAC{at('1', 630)}"]
    mate = image("2", 350)
    assert r["ins_mate_minus"]["alts"] == [f"[{mate[0]}:{mate[1]}[AC{at('1', 230)}"]
    one, two = image("1", 300), image("2", 400)
    ref = at("1", 270)
    assert r["two_alts"]["alts"] == [f"{ref}[{one[0]}:{one[1]}[", f"{ref}[{two[0]}:{two[1]}["]
    assert r["two_alts"]["gt"] == [1, 2] and r["two_alts"]["gt_alleles"] == r["two_alts"]["alts"]


def test_intervals_mirror_per_end_and_never_swap(tmp_path: Path, setup) -> None:
    """CIPOS is about POS, CIEND (truvari's reading) about the mate: each mirrors when its own end
    lands on a minus block. Unlike a symbolic SV, POS stays POS, so they never trade places."""
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = by_case(lifted.loaded.calls.collect())
    assert (calls["pp_lr"]["cipos"], calls["pp_lr"]["ciend"]) == ([-3, 5], [-2, 7])
    assert (calls["pm_lr"]["cipos"], calls["pm_lr"]["ciend"]) == ([-3, 5], [-7, 2])
    assert (calls["mp_lr"]["cipos"], calls["mp_lr"]["ciend"]) == ([-5, 3], [-2, 7])
    assert (calls["mm_lr"]["cipos"], calls["mm_lr"]["ciend"]) == ([-5, 3], [-7, 2])
    assert (calls["strands"]["cipos"], calls["strands"]["ciend"]) == (None, None)


def test_the_lifted_facet_names_its_mate_contigs(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    facet = lifted.loaded.facets.collect().row(0, named=True)
    # No row lands on chr3, but pm_* point at it: every contig an ALT names is in the map.
    assert {entry["source"] for entry in facet["contig_map"]} == {"chr1", "chr2", "chr3"}
    assert "3" not in set(lifted.loaded.calls.collect()["chrom"])
    # A second run reuses the cached generation and reads the same map back.
    _, again = lift_cases(tmp_path, store, source, forward)
    assert again.loaded.facets.collect().row(0, named=True)["contig_map"] == facet["contig_map"]
    lifted_count = sum(reason is None for reason in EXPECTED_REASON.values())
    assert any(warning.startswith(f"{lifted_count} lifted row(s) were breakends") for warning in lifted.loaded.warnings)
    # mate_gap's mate in a gap, mate_no_chain's on a contig the chain does not cover.
    assert "breakend mate(s) not lifted, by the mate position's reason: 1 contig_not_in_chain, 1 no_chain_block" in lifted.loaded.warnings


def test_a_lift_and_its_inverse_give_back_every_breakend_both_lifted(tmp_path: Path, setup) -> None:
    store, source, _, forward, back = setup
    resolved, there = lift_cases(tmp_path, store, source, forward)
    returned = liftover(there.loaded, back, store, target=source)
    original = by_case(resolved.calls.collect())
    home = by_case(returned.loaded.calls.collect())
    assert returned.reasons == {}
    assert set(home) == {case for case, reason in EXPECTED_REASON.items() if reason is None}
    columns = ("chrom", "pos", "end_pos", "ref", "alts", "gt", "gt_alleles", "cipos", "ciend", "sv_type")
    for case, r in home.items():
        assert {c: r[c] for c in columns} == {c: original[case][c] for c in columns}, case


def test_a_lifted_breakend_writes_and_reloads(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = by_case(lifted.loaded.calls.collect())
    # Under `source` names the lifted facet's `contig_map` is the target's own spelling, so every ALT
    # comes back verbatim; under `stored` the mate is named as CHROM is, canonically (F30).
    for style in (ContigNames.SOURCE, ContigNames.STORED):
        path = tmp_path / f"lifted-{style.value}.vcf.gz"
        write_vcf(lifted.loaded, path, contig_names=style)
        again = by_case(load(path).calls.collect())
        assert set(again) == set(calls) and calls
        for case, r in calls.items():
            alts = r["alts"] if style is ContigNames.SOURCE else [canonical_alt(alt) for alt in r["alts"]]
            assert (again[case]["chrom"], again[case]["pos"], again[case]["source_ref"], again[case]["alts"]) == (r["chrom"], r["pos"], r["ref"], alts)


def canonical_alt(alt: str) -> str:
    """A lifted ALT with its mate contig canonical (`chr2` is `2`), as `stored` names write it."""
    _, contig, _, _ = parse(alt)
    return alt if contig is None else alt.replace(f"{contig}:", f"{canonical_chrom(contig)}:")


def test_a_lifted_facet_written_with_chr_names_declares_every_mate(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward)
    calls = lifted.loaded.calls.collect()
    path = write_vcf(lifted.loaded, tmp_path / "lifted.vcf.gz", contig_names=ContigNames.CHR).path
    text = gzip.decompress(path.read_bytes()).decode()
    declared = {line.split("ID=")[1].split(",")[0].rstrip(">") for line in text.splitlines() if line.startswith("##contig=")}
    rows = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]

    # The chr style names 1, 2, 3 as chr1, chr2, chr3: the lifted mates, read canonically, then named so.
    lifted_mates = {canonical_chrom(parse(alt)[1]) for alts in calls["alts"] for alt in alts if parse(alt)[1] is not None}
    written_mates = {parse(alt)[1] for row in rows for alt in row[4].split(",") if parse(alt)[1] is not None}
    assert lifted_mates == {"1", "2", "3"}, "the fixture's mates land on all three target contigs"
    assert written_mates == {f"chr{chrom}" for chrom in lifted_mates}
    assert written_mates <= declared, "every contig a written ALT names has a ##contig line"
    # chr3 is declared although no row stands on it: only the pm_* mates point there.
    assert "3" not in set(calls["chrom"]) and "chr3" in declared
    assert {row[0] for row in rows} == {f"chr{chrom}" for chrom in set(calls["chrom"])}


def test_a_generation_without_the_mate_contigs_group_is_written_again(tmp_path: Path, setup) -> None:
    """45126a2 committed lift generations without `mate_contigs` at the `LIFT_VERSION` 0d3599a kept, and
    0d3599a's lift found one and failed with a `KeyError`. A generation holding other groups than the
    lift's is never read: the key covers the groups, and one found under it anyway (here its
    `result.json` edited, the group dropped) is written again."""
    store, source, _, forward, _ = setup
    resolved, lifted = lift_cases(tmp_path, store, source, forward)
    first = lifted.loaded.facets.collect().row(0, named=True)["contig_map"]
    (result,) = [path for path in cache_dir().glob("frames/*/*/result.json") if "mate_contigs" in json.loads(path.read_text())["files"]]
    record = json.loads(result.read_text())
    assert record["files"]["mate_contigs"]  # pm_* name chr3 as a mate, so the group is not empty
    del record["files"]["mate_contigs"]
    result.write_text(json.dumps(record))
    again = liftover(resolved, forward, store)
    assert again.loaded.facets.collect().row(0, named=True)["contig_map"] == first
    assert "mate_contigs" in json.loads(result.read_text())["files"]


def test_each_alt_of_a_record_mirrors_its_own_interval_pair(tmp_path: Path, setup) -> None:
    """Two bracket ALTs in one record with a CIPOS and a CIEND pair each (VCF 4.5's shape): CIPOS pair
    `k` mirrors with POS's strand, CIEND pair `k` with ALT `k`'s mate's, `[lo, hi]` to `[-hi, -lo]`.
    The mates of each record land on two strands, so a pair read for the wrong ALT shows."""
    store, source, _, forward, _ = setup
    ci = "SVTYPE=BND;CIPOS=-1,2,-3,4;CIEND=-5,6,-7,8"
    cases = {
        # POS on P1 (plus); mate 0 on P1 (plus), mate 1 on M2 (minus).
        "two_ci_p": bnd("two_ci_p", "1", 270, "{r}[1:300[,{r}]2:400]", gt="1/2", info=ci),
        # POS on M1 (minus); mate 0 on P2 (plus), mate 1 on M1 (minus).
        "two_ci_m": bnd("two_ci_m", "1", 600, "{r}[2:100[,{r}]1:700]", gt="1/2", info=ci),
    }
    mates = {"two_ci_p": [("1", 300), ("2", 400)], "two_ci_m": [("2", 100), ("1", 700)]}
    _, lifted = lift_cases(tmp_path, store, source, forward, cases)
    calls = by_case(lifted.loaded.calls.collect())
    assert set(calls) == set(cases)

    def mirrored(pair: list[int], minus: bool) -> list[int]:
        return [-pair[1], -pair[0]] if minus else pair

    cipos, ciend = [[-1, 2], [-3, 4]], [[-5, 6], [-7, 8]]
    for case, (chrom, pos) in (("two_ci_p", ("1", 270)), ("two_ci_m", ("1", 600))):
        own_minus = image(chrom, pos)[2] == "-"
        mate_minus = [image(*mate)[2] == "-" for mate in mates[case]]
        assert mate_minus[0] != mate_minus[1], case
        assert calls[case]["cipos"] == [v for k in range(2) for v in mirrored(cipos[k], own_minus)], case
        assert calls[case]["ciend"] == [v for k in range(2) for v in mirrored(ciend[k], mate_minus[k])], case


# VCF 4.3 §5.4.2's large insertion: a mate in angle brackets names a contig of the `##assembly` file
# (§1.6.1), not of the reference, so no chain moves it. The record's own end lifts; the mate's contig,
# position and side stay as written, and only POS's strand turns the form (`FORM_AFTER` with the mate
# on a plus strand). `CTG1` is that contig's sequence for the oracle alone (`junction`): it is in no
# FASTA, and the lift never reads it.
CTG1 = random_bases(400)
ASSEMBLY_CASES = {
    # POS on P1 (plus): the ALT comes back as written.
    "asm_plus_lr": bnd("asm_plus_lr", "1", 240, "{r}[<ctg1>:7["),
    "asm_plus_rl": bnd("asm_plus_rl", "1", 241, "]<ctg1>:214]{r}"),
    # POS on M1 (minus): the record's side flips, the mate's does not.
    **four("asm_m", "1", 640, "<ctg1>:7", info=CI),
    # Inserted bases on a minus own end: reverse-complemented to the other side of REF, as for any mate.
    "asm_ins_minus": bnd("asm_ins_minus", "1", 650, "{r}GTCA[<ctg1>:20["),
    # Two ALTs, an assembly mate beside a reference mate on a minus block (M2), which still lifts and turns.
    "asm_two_alts": bnd("asm_two_alts", "1", 280, "{r}[<ctg1>:7[,{r}]2:400]", gt="1/2"),
}


def test_an_assembly_mate_stays_as_written_and_pos_alone_turns_the_form(tmp_path: Path, setup) -> None:
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward, ASSEMBLY_CASES)
    report = by_case(lifted.report.collect())
    assert {case: r["lift_reason"] for case, r in report.items()} == {case: None for case in ASSEMBLY_CASES}
    calls = by_case(lifted.loaded.calls.collect())
    for case, r in calls.items():
        source_row = report[case]
        own = image("1", source_row["pos"])
        assert (r["chrom"], r["pos"], source_row["lift_strand"]) == (own[0].removeprefix("chr"), own[1], own[2]), case
        assert r["ref"] == (comp(source_row["ref"]) if own[2] == "-" else source_row["ref"])
        for source_alt, lifted_alt in zip(source_row["alts"], r["alts"], strict=True):
            form, mate_chrom, mate_pos, _ = parse(source_alt)
            if mate_chrom.startswith("<"):
                mate_strand, mate_at = "+", (mate_chrom, mate_pos)
            else:
                mate = image(mate_chrom, mate_pos)
                mate_strand, mate_at = mate[2], (mate[0], mate[1])
            lifted_form, lifted_chrom, lifted_pos, _ = parse(lifted_alt)
            assert lifted_form == FORM_AFTER[(form, own[2], mate_strand)], (case, source_alt, lifted_alt)
            assert (lifted_chrom, lifted_pos) == mate_at, (case, source_alt, lifted_alt)
    # On a plus own end the ALT is unchanged; on a minus one it turns, and the mate's CIEND never moves.
    assert calls["asm_plus_lr"]["alts"] == report["asm_plus_lr"]["alts"]
    assert calls["asm_m_lr"]["alts"] == [f"[<ctg1>:7[{calls['asm_m_lr']['ref']}"]
    assert (calls["asm_m_lr"]["cipos"], calls["asm_m_lr"]["ciend"]) == ([-5, 3], [-2, 7])


def test_a_lifted_assembly_mate_asserts_the_same_novel_sequence(tmp_path: Path, setup) -> None:
    """`test_a_lifted_breakend_asserts_the_same_novel_sequence`'s check, with the assembly contig's
    sequence beside both assemblies': no lift moves it, so it is the same on both sides."""
    store, source, _, forward, _ = setup
    _, lifted = lift_cases(tmp_path, store, source, forward, ASSEMBLY_CASES)
    calls = by_case(lifted.loaded.calls.collect())
    report = by_case(lifted.report.collect())
    names = {**TARGET, **{name.removeprefix("chr"): seq for name, seq in TARGET.items()}, "<ctg1>": CTG1}
    checked = 0
    for case, r in calls.items():
        source_row = report[case]
        for source_alt, lifted_alt in zip(source_row["alts"], r["alts"], strict=True):
            before = junction({**SOURCE, "<ctg1>": CTG1}, source_row["chrom"], source_row["pos"], source_alt)
            after = junction(names, r["chrom"], r["pos"], lifted_alt)
            assert after in (before, revcomp(before)), (case, source_alt, lifted_alt)
            checked += 1
    assert checked == sum(len(r["alts"]) for r in calls.values()) == len(ASSEMBLY_CASES) + 1


def test_an_assembly_mate_names_no_contig_and_lifts_back(tmp_path: Path, setup) -> None:
    store, source, _, forward, back = setup
    resolved, there = lift_cases(tmp_path, store, source, forward, ASSEMBLY_CASES)
    facet = there.loaded.facets.collect().row(0, named=True)
    # The target collection has no `<ctg1>`: the map lists the contigs rows stand on and the reference
    # mates name (asm_two_alts' 2:400 lands on chr3), never the assembly contig.
    assert {entry["source"] for entry in facet["contig_map"]} == {"chr1", "chr2", "chr3"}
    assert any("mate(s) on a contig of the ##assembly file kept as written" in warning for warning in there.loaded.warnings)
    returned = liftover(there.loaded, back, store, target=source)
    assert returned.reasons == {}
    original, home = by_case(resolved.calls.collect()), by_case(returned.loaded.calls.collect())
    columns = ("chrom", "pos", "end_pos", "ref", "alts", "gt", "gt_alleles", "cipos", "ciend", "sv_type")
    assert set(home) == set(ASSEMBLY_CASES)
    for case, r in home.items():
        assert {c: r[c] for c in columns} == {c: original[case][c] for c in columns}, case
