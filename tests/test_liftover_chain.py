"""RM6: `liftover` of blocks and symbolic SVs through a real chain. Opt-in: a reference store holding
GRCh37 and GRCh38 (`JUST_VCF_TEST_REFERENCE_DIR`) and UCSC's `hg19ToHg38.over.chain.gz`
(`JUST_VCF_TEST_CHAIN_DIR`).

The blocks and SVs are synthetic, drawn at seeded positions on the real GRCh37 autosomes; no genome is
read. Expected positions come from pyliftover's reading of the same chain, base by base, never from
what the code under test wrote.
"""

import random
from collections import Counter
from pathlib import Path

import pytest

from just_vcf import ReferenceStore, liftover, load, resolve_ref
from just_vcf.breakends import parse_breakend
from just_vcf.reference import Contig
from just_vcf.schema import canonical_chrom

CHAIN = "hg19ToHg38.over.chain.gz"
AUTOSOMES = [str(n) for n in range(1, 23)]
PRIMARY = {*AUTOSOMES, "X", "Y"}


@pytest.fixture(scope="module")
def store(test_reference_dir: Path) -> ReferenceStore:
    found = ReferenceStore.open(test_reference_dir)
    labels = {found.assembly_label(digest) for digest in found.collections()}
    if not {"GRCh37", "GRCh38"} <= labels:
        pytest.skip(f"{test_reference_dir} holds no GRCh37 and GRCh38 collections")
    return found


@pytest.fixture(scope="module")
def chain(test_chain_dir: Path) -> Path:
    path = test_chain_dir / CHAIN
    if not path.is_file():
        pytest.skip(f"{path} is absent")
    return path


def collection(store: ReferenceStore, label: str) -> str:
    (digest,) = [digest for digest in store.collections() if store.assembly_label(digest) == label]
    return digest


def by_chrom(store: ReferenceStore, digest: str) -> dict[str, Contig]:
    return {chrom: contigs[0] for chrom, contigs in store.contigs_by_chrom(digest).items() if len(contigs) == 1}


def base_fate(lifter, chrom: str, pos: int) -> tuple[str, int, str] | None:
    """pyliftover's one target of a 1-based GRCh37 base, as (chrom, 1-based pos, strand); None for
    none, several, or a target off the primary contigs (the comparison is kept to those)."""
    hits = lifter.convert_coordinate(f"chr{chrom}", pos - 1)
    if not hits or len(hits) != 1:
        return None
    name, zero, strand, _ = hits[0]
    target = name.removeprefix("chr")
    return (target, zero + 1, strand) if target in PRIMARY else None


def expected_pieces(lifter, chrom: str, start: int, end: int) -> list[tuple[str, int, int]]:
    pieces: list[list] = []
    previous = None
    for pos in range(start, end + 1):
        fate = base_fate(lifter, chrom, pos)
        if fate is None:
            previous = None
            continue
        step = 1 if fate[2] == "+" else -1
        if previous is not None and (previous[0], previous[1] + step, previous[2]) == fate:
            pieces[-1][2 if fate[2] == "+" else 1] = fate[1]
        else:
            pieces.append([fate[0], fate[1], fate[1]])
        previous = fate
    return [tuple(piece) for piece in pieces]


def test_blocks_lift_as_pyliftover_maps_their_bases(store: ReferenceStore, chain: Path, tmp_path: Path) -> None:
    pyliftover = pytest.importorskip("pyliftover")
    source = collection(store, "GRCh37")
    contigs = by_chrom(store, source)
    rng = random.Random(11)
    blocks = []
    for chrom in AUTOSOMES:
        for _ in range(6):
            start = rng.randrange(1_000_000, contigs[chrom].length - 1_000_000)
            blocks.append((chrom, start, start + rng.randrange(0, 30_000)))
    lines = [
        "##fileformat=VCFv4.2",
        *[f"##contig=<ID={chrom},length={contigs[chrom].length}>" for chrom in AUTOSOMES],
        '##ALT=<ID=NON_REF,Description="Any other allele">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
    ]
    for chrom, start, end in sorted(blocks, key=lambda block: (int(block[0]), block[1])):
        base = store.fetch(contigs[chrom], [(start - 1, start)])[0].upper()
        lines.append(f"{chrom}\t{start}\t.\t{base}\t<NON_REF>\t.\t.\tEND={end}\tGT:GQ\t0/0:30")
    vcf = tmp_path / "blocks.g.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    resolved = resolve_ref(load(vcf), store, collection=source)
    lifted = liftover(resolved, chain, store, target=collection(store, "GRCh38"))

    lifter = pyliftover.LiftOver(str(chain))
    expected: Counter = Counter()
    for chrom, start, end in blocks:
        expected.update(expected_pieces(lifter, chrom, start, end))
    got = Counter((r["chrom"], r["pos"], r["end_pos"]) for r in lifted.loaded.blocks.collect().iter_rows(named=True) if r["chrom"] in PRIMARY)
    assert got == expected
    total = sum(end - start + 1 for _, start, end in blocks)
    lifted_bases = sum(r["end_pos"] - r["pos"] + 1 for r in lifted.loaded.blocks.collect().iter_rows(named=True))
    assert lifted_bases + sum(lifted.block_reasons.values()) == total


def point_fate(lifter, chrom: str, pos: int) -> tuple[str, int, str] | str | None:
    """pyliftover's one target of a 1-based GRCh37 base; None for no target, `skip` for several or one
    off the primary contigs (those rows are left out of the comparison)."""
    hits = lifter.convert_coordinate(f"chr{chrom}", pos - 1)
    if not hits:
        return None
    if len(hits) != 1 or hits[0][0].removeprefix("chr") not in PRIMARY:
        return "skip"
    name, zero, strand, _ = hits[0]
    return name.removeprefix("chr"), zero + 1, strand


def test_deletions_lift_by_both_ends_as_pyliftover_places_them(store: ReferenceStore, chain: Path, tmp_path: Path) -> None:
    """`<DEL>`s of 50 bases to 200 kb: each end placed as pyliftover places it, and the SV kept or
    left behind by the rule (one contig and strand, in order, length within 5%)."""
    pyliftover = pytest.importorskip("pyliftover")
    source = collection(store, "GRCh37")
    target = collection(store, "GRCh38")
    contigs, target_contigs = by_chrom(store, source), by_chrom(store, target)
    rng = random.Random(13)
    svs = []
    for chrom in AUTOSOMES:
        for index in range(25):
            pos = rng.randrange(1_000_000, contigs[chrom].length - 1_000_000)
            svs.append((f"{chrom}_{index}", chrom, pos, pos + int(10 ** rng.uniform(1.7, 5.3))))
    lines = [
        "##fileformat=VCFv4.2",
        *[f"##contig=<ID={chrom},length={contigs[chrom].length}>" for chrom in AUTOSOMES],
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
    ]
    for name, chrom, pos, end in sorted(svs, key=lambda sv: (int(sv[1]), sv[2])):
        base = store.fetch(contigs[chrom], [(pos - 1, pos)])[0].upper()
        lines.append(f"{chrom}\t{pos}\t{name}\t{base}\t<DEL>\t.\tPASS\tEND={end}\tGT\t0/1")
    vcf = tmp_path / "svs.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    resolved = resolve_ref(load(vcf), store, collection=source)
    lifted = liftover(resolved, chain, store, target=target)
    report = {r["ids"][0]: r for r in lifted.report.collect().iter_rows(named=True)}
    calls = {r["ids"][0]: r for r in lifted.loaded.calls.collect().iter_rows(named=True)}

    lifter = pyliftover.LiftOver(str(chain))
    compared: Counter = Counter()
    for name, chrom, pos, end in svs:
        at_pos, at_end = point_fate(lifter, chrom, pos), point_fate(lifter, chrom, end)
        if "skip" in (at_pos, at_end):
            continue
        reason = report[name]["lift_reason"]
        if at_pos is None:
            assert reason is not None and not reason.startswith("sv_"), name
        elif at_end is None:
            assert reason == "sv_end_unlifted", name
        elif at_pos[0] != at_end[0] or at_pos[2] != at_end[2]:
            assert reason == "sv_ends_apart", name
        else:
            (target_chrom, a, strand), b = at_pos, at_end[1]
            new_pos, new_end = (a, b) if strand == "+" else (b - 1, a - 1)
            if (b < a) if strand == "+" else (b > a):
                assert reason == "sv_ends_reversed", name
            elif abs((new_end - new_pos) - (end - pos)) > 0.05 * (end - pos):
                assert reason == "sv_length_changed", name
            elif strand == "+" and store.fetch(target_contigs[target_chrom], [(new_pos - 1, new_pos)])[0].upper() != report[name]["ref"]:
                assert reason == "ref_mismatch_after_lift", name
            else:
                assert reason is None, name
                assert (calls[name]["chrom"], calls[name]["pos"], calls[name]["end_pos"]) == (target_chrom, new_pos, new_end), name
        compared[reason] += 1
    # Most of the sample is compared, and most of it lifts (the synthetic suite covers each way to stay behind).
    assert compared[None] > 100 and sum(compared.values()) > 0.9 * len(svs)


MATE_FORM = {"t[p[": "]p]t", "]p]t": "t[p[", "t]p]": "t]p]", "[p[t": "[p[t"}
"""The form a breakend's mate record writes for the same adjacency (VCF 4.3 §5.4's examples)."""

STRAND_TABLE = {
    # (form, own strand, mate strand) -> lifted form; tests/test_liftover_bnd.py's FORM_AFTER, by form.
    **{(form, "+", "+"): form for form in MATE_FORM},
    ("t[p[", "+", "-"): "t]p]", ("t[p[", "-", "+"): "[p[t", ("t[p[", "-", "-"): "]p]t",
    ("t]p]", "+", "-"): "t[p[", ("t]p]", "-", "+"): "]p]t", ("t]p]", "-", "-"): "[p[t",
    ("]p]t", "+", "-"): "[p[t", ("]p]t", "-", "+"): "t]p]", ("]p]t", "-", "-"): "t[p[",
    ("[p[t", "+", "-"): "]p]t", ("[p[t", "-", "+"): "t[p[", ("[p[t", "-", "-"): "t]p]",
}  # fmt: skip


def bnd_alt(form: str, t: str, contig: str, pos: int) -> str:
    p = f"{contig}:{pos}"
    return {"t[p[": f"{t}[{p}[", "t]p]": f"{t}]{p}]", "]p]t": f"]{p}]{t}", "[p[t": f"[{p}[{t}"}[form]


def test_translocations_lift_by_both_ends_and_back(store: ReferenceStore, chain: Path, test_chain_dir: Path, tmp_path: Path) -> None:
    """Pairs of mate records joining two autosomes, in every form: each end placed as pyliftover places
    it, the form turned by the two ends' strands (the table), the two records still each other's mates,
    and `hg38ToHg19` giving back every record both lifts kept."""
    pyliftover = pytest.importorskip("pyliftover")
    back = test_chain_dir / "hg38ToHg19.over.chain.gz"
    if not back.is_file():
        pytest.skip(f"{back} is absent")
    source, target = collection(store, "GRCh37"), collection(store, "GRCh38")
    contigs = by_chrom(store, source)
    rng = random.Random(17)
    records = []
    forms = list(MATE_FORM)
    for index in range(400):
        a, b = rng.sample(AUTOSOMES, 2)
        x = rng.randrange(1_000_000, contigs[a].length - 1_000_000)
        y = rng.randrange(1_000_000, contigs[b].length - 1_000_000)
        form = forms[index % 4]
        records.append((f"j{index}a", a, x, form, b, y, f"j{index}b"))
        records.append((f"j{index}b", b, y, MATE_FORM[form], a, x, f"j{index}a"))
    lines = [
        "##fileformat=VCFv4.3",
        *[f"##contig=<ID={chrom},length={contigs[chrom].length}>" for chrom in AUTOSOMES],
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
    ]
    for name, chrom, pos, form, mate_chrom, mate_pos, _ in sorted(records, key=lambda r: (int(r[1]), r[2])):
        base = store.fetch(contigs[chrom], [(pos - 1, pos)])[0].upper()
        lines.append(f"{chrom}\t{pos}\t{name}\t{base}\t{bnd_alt(form, base, mate_chrom, mate_pos)}\t.\tPASS\tSVTYPE=BND\tGT\t0/1")
    vcf = tmp_path / "bnd.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    resolved = resolve_ref(load(vcf), store, collection=source)
    there = liftover(resolved, chain, store, target=target)
    report = {r["ids"][0]: r for r in there.report.collect().iter_rows(named=True)}
    calls = {r["ids"][0]: r for r in there.loaded.calls.collect().iter_rows(named=True)}

    lifter = pyliftover.LiftOver(str(chain))
    compared: Counter = Counter()
    for name, chrom, pos, form, mate_chrom, mate_pos, mate_name in records:
        own, mate = point_fate(lifter, chrom, pos), point_fate(lifter, mate_chrom, mate_pos)
        if "skip" in (own, mate):
            continue
        reason = report[name]["lift_reason"]
        if own is None:
            assert reason is not None and reason != "breakend_mate_not_lifted", name
        elif mate is None:
            assert reason == "breakend_mate_not_lifted", name
        elif reason != "ref_mismatch_after_lift":
            assert reason is None, name
            lifted = parse_breakend(calls[name]["alts"][0])
            assert (calls[name]["chrom"], calls[name]["pos"]) == own[:2], name
            assert (canonical_chrom(lifted.mate_contig), lifted.mate_pos, lifted.form) == (mate[0], mate[1], STRAND_TABLE[(form, own[2], mate[2])]), name
            if mate_name in calls:
                # The mate record, lifted on its own, still stands where this one's ALT points.
                assert (calls[mate_name]["chrom"], calls[mate_name]["pos"]) == (canonical_chrom(lifted.mate_contig), lifted.mate_pos), name
                compared["pair"] += 1
            compared[own[2] + mate[2]] += 1
        compared[reason] += 1
    assert compared[None] > 500 and compared["pair"] > 400

    home = liftover(there.loaded, back, store, target=source)
    returned = {r["ids"][0]: r for r in home.loaded.calls.collect().iter_rows(named=True)}
    original = {r["ids"][0]: r for r in resolved.calls.collect().iter_rows(named=True)}
    for name, r in returned.items():
        before, after = parse_breakend(original[name]["alts"][0]), parse_breakend(r["alts"][0])
        assert (r["chrom"], r["pos"], r["ref"]) == (original[name]["chrom"], original[name]["pos"], original[name]["ref"]), name
        assert (after.form, canonical_chrom(after.mate_contig), after.mate_pos, after.bases) == (before.form, before.mate_contig, before.mate_pos, before.bases), name
    # Nearly every record both lifts keep comes back (a chain is not exactly its inverse's mirror).
    assert len(returned) > 0.95 * len(calls)
