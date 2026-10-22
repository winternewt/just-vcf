"""RM6: `normalize` splits multi-allelic sites, left-aligns and trims indels, dedups and sorts, as a
derived facet. Synthetic FASTA and VCFs under tmp_path; expectations computed from the sequence."""

import random
from collections import Counter
from math import comb
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceMismatchError, ReferenceStore, load, normalize, rejoin_split, resolve_ref, sniff, write_vcf
from just_vcf.loading import Loaded
from just_vcf.normalize import RELATION, _left_align
from just_vcf.schema import BLOCKS, CALLS

rng = random.Random(11)
_LEFT = "".join(rng.choice("ACGT") for _ in range(200))
_RIGHT = "".join(rng.choice("ACGT") for _ in range(300))
HOMOPOLYMER = 202  # 1-based start of a T run of 8, after the C at 201
TANDEM = 301  # 1-based start of six CA copies, after a G
SEQ1 = "AAAA" + _LEFT[4:200] + "C" + "TTTTTTTT" + _LEFT[: 300 - 209 - 1] + "G" + "CA" * 6 + _RIGHT
SEQ = {"1": SEQ1, "2": "".join(rng.choice("ACGT") for _ in range(600))}


def at(chrom: str, pos: int, length: int = 1) -> str:
    """The synthetic reference at 1-based `pos`, the test's own reading of it."""
    return SEQ[chrom][pos - 1 : pos - 1 + length]


def other(base: str) -> str:
    return next(b for b in "ACGT" if b != base)


def write_fasta(path: Path) -> Path:
    lines = []
    for name, seq in SEQ.items():
        lines.append(f">{name}")
        lines += [seq[i : i + 60] for i in range(0, len(seq), 60)]
    path.write_text("\n".join(lines) + "\n")
    return path


HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SEQ.items()],
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posteriors">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def record(chrom: str, pos: int, ref: str, alt: str, gt: str, vectors: bool = True) -> str:
    """A VCF line; with `vectors`, AD, PL and a GP summing to 1, all distinct, sized to the alleles."""
    alleles = len(alt.split(",")) + 1
    ploidy = len(gt.replace("|", "/").split("/"))
    genotypes = comb(alleles + ploidy - 1, ploidy)
    if not vectors:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}"
    ad = ",".join(str(10 + i) for i in range(alleles))
    pl_ = ",".join(str(7 * i) for i in range(genotypes))
    weights = [genotypes - i for i in range(genotypes)]
    gp = ",".join(f"{w / sum(weights):.6f}" for w in weights)
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:AD:PL:GP\t{gt}:{ad}:{pl_}:{gp}"


def write_vcf_text(path: Path, records: list[str]) -> Path:
    path.write_text("\n".join([*HEADER, *records]) + "\n")
    return path


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa"))
    return s


def normalized(tmp_path: Path, store: ReferenceStore, records: list[str]):
    parent = resolve_ref(load(write_vcf_text(tmp_path / "in.vcf", records)), store)
    return parent, normalize(parent, store)


def apply(chrom: str, pos: int, ref: str, alt: str) -> str:
    """The haplotype a biallelic variant makes of the reference."""
    seq = SEQ[chrom]
    assert seq[pos - 1 : pos - 1 + len(ref)].upper() == ref.upper()
    return seq[: pos - 1] + alt + seq[pos - 1 + len(ref) :]


# --- splitting ---------------------------------------------------------------------------------


def test_split_writes_dot_for_the_other_alt_never_zero(tmp_path: Path, store: ReferenceStore) -> None:
    p = 50
    ref = at("1", p)
    g, t = [b for b in "ACGT" if b != ref][:2]
    records = [
        record("1", p, ref, f"{g},{t}", "1/2"),
        record("1", p + 10, at("1", p + 10), f"{other(at('1', p + 10))},{other(other(at('1', p + 10)))}", "1/1"),
        record("1", p + 20, at("1", p + 20), ",".join(b for b in "ACGT" if b != at("1", p + 20)), "0|2"),
    ]
    _, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    assert calls.schema == CALLS
    view = {(row["pos"], row["alts"][0]): (row["gt"], row["gt_alleles"], row["phased"]) for row in calls.iter_rows(named=True)}
    assert view[(p, g)][:2] == ([1, None], [g, None])
    assert view[(p, t)][:2] == ([None, 1], [None, t])
    # 1/1 at a two-ALT site: the uncarried ALT's row holds nothing, not a hom-ref call.
    alts = [other(at("1", p + 10)), other(other(at("1", p + 10)))]
    assert view[(p + 10, alts[0])][0] == [1, 1]
    assert view[(p + 10, alts[1])][:2] == ([None, None], [None, None])
    # A phased 0|2 at a three-ALT site: REF stays REF in every row, phasing kept per slot.
    third = [b for b in "ACGT" if b != at("1", p + 20)]
    assert [view[(p + 20, b)][0] for b in third] == [[0, None], [0, 1], [0, None]]
    assert all(view[(p + 20, b)][2] == [True, True] for b in third)
    assert not any(0 in [i for i in (row["gt"] or []) if i is not None] and row["gt_alleles"][row["gt"].index(0)] != row["ref"]
                   for row in calls.iter_rows(named=True))
    assert set(calls["split_signal"].to_list()) == {"marker"}
    assert calls["source_gt"].null_count() == calls.height
    assert calls.height == 2 + 2 + 3
    assert any(w.startswith("3 multi-allelic site(s) split into 7 row(s)") for w in out.warnings)


def test_vectors_are_sliced_exactly_or_withheld(tmp_path: Path, store: ReferenceStore) -> None:
    p = 80
    ref = at("1", p)
    alts = [b for b in "ACGT" if b != ref]
    parent, out = normalized(tmp_path, store, [record("1", p, ref, ",".join(alts), "0/2")])
    source = parent.calls.collect().row(0, named=True)
    calls = {row["alts"][0]: row for row in out.calls.collect().iter_rows(named=True)}
    for k, alt in enumerate(alts, start=1):
        row = calls[alt]
        assert row["ad"] == [source["ad"][0], source["ad"][k]]
        assert row["gp"] is None and row["dosage"] is None
        if k == 2:
            # The genotype is REF/ALT2 and nothing else: the slice holds the called genotype.
            het = k * (k + 1) // 2
            assert row["gl"] == [source["gl"][0], source["gl"][het], source["gl"][het + k]]
            assert row["source_alleles"] == [ref, alt]
        else:
            # The genotype holds ALT2 in this row's `.` slot: no slice over [REF, ALT] holds it.
            assert row["gl"] is None
    assert any("GP withheld" in w for w in out.warnings)
    assert any(w.startswith("2 split row(s): GL withheld") for w in out.warnings)


def test_symbolic_sites_are_left_as_written(tmp_path: Path, store: ReferenceStore) -> None:
    p = 120
    ref = at("1", p)
    records = [record("1", p, ref, f"{other(ref)},<*>", "0/1"), record("1", p + 5, at("1", p + 5), f"{other(at('1', p + 5))},*", "1/2")]
    parent, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    assert calls.filter(pl.col("pos") == p)["alts"].to_list() == [[other(ref), "<*>"]]
    assert calls.filter(pl.col("pos") == p)["split_signal"].to_list() == [None]
    # `*` is an allele like any other here, as bcftools splits it.
    assert sorted(calls.filter(pl.col("pos") == p + 5)["alts"].to_list()) == sorted([[other(at("1", p + 5))], ["*"]])
    assert any("1 multi-allelic row(s) left as written: a symbolic allele" in w for w in out.warnings)


# --- left-alignment and trimming ----------------------------------------------------------------


def is_normalized(pos: int, ref: str, alt: str) -> bool:
    """Tan et al.'s definition: parsimonious (no common first base while both keep one) and left-aligned
    (no common last base), except at position 1, where an allele may keep its last base."""
    parsimonious = not (len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0])
    left = ref[-1] != alt[-1] or (pos == 1 and min(len(ref), len(alt)) == 1)
    return parsimonious and left


def test_indels_move_left_through_repeats(tmp_path: Path, store: ReferenceStore) -> None:
    end_t = HOMOPOLYMER + 7
    end_ca = TANDEM + 11
    cases = [
        ("1", end_t - 1, at("1", end_t - 1, 2), at("1", end_t - 1)),  # deletes the run's last T
        ("1", end_t, at("1", end_t), at("1", end_t) + "TT"),  # inserts TT at the run's end
        ("1", end_ca - 2, at("1", end_ca - 2, 3), at("1", end_ca - 2)),  # deletes the last CA
        ("1", end_ca, at("1", end_ca), at("1", end_ca) + "CACA"),  # inserts two copies after the repeat
        ("1", 160, at("1", 160, 3), other(at("1", 160)) + at("1", 161, 2)),  # MNP that trims to an SNV
    ]
    records = [record(c, p, r, a, "0/1") for c, p, r, a in cases]
    parent, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    got = [(row["pos"], row["ref"], row["alts"][0]) for row in calls.iter_rows(named=True)]
    assert len(got) == len(cases)
    haplotypes = {apply(c, p, r, a) for c, p, r, a in cases}
    assert {apply("1", p, r, a) for p, r, a in got} == haplotypes
    assert all(is_normalized(p, r, a) for p, r, a in got)
    # The run and the repeat start after a C and a G, so each indel lands on that anchor.
    assert (HOMOPOLYMER - 1, "CT", "C") in got and (HOMOPOLYMER - 1, "C", "CTT") in got
    assert (TANDEM - 1, "GCA", "G") in got and (TANDEM - 1, "G", "GCACA") in got
    assert (160, at("1", 160), other(at("1", 160))) in got
    moved = calls.filter(pl.col("pos") != 160)
    assert moved["source_ref"].null_count() == moved["ref_check"].null_count() == moved.height
    # The vectors follow the alleles: same order, new spelling.
    assert all(row["source_alleles"] == [row["ref"], *row["alts"]] for row in calls.iter_rows(named=True))
    source_gl = {apply("1", r["pos"], r["source_ref"], r["alts"][0]): r["gl"] for r in parent.calls.collect().iter_rows(named=True)}
    assert all(row["gl"] == source_gl[apply("1", row["pos"], row["ref"], row["alts"][0])] for row in calls.iter_rows(named=True))
    assert calls["end_pos"].to_list() == [p + len(r) - 1 for p, r, _ in got]
    assert any(w.startswith("4 row(s) left-aligned to another position and 1 trimmed in place") for w in out.warnings)


def test_random_indels_match_the_definition(tmp_path: Path, store: ReferenceStore) -> None:
    chooser = random.Random(5)
    seq = SEQ["1"]
    cases = {}
    for _ in range(120):
        pos = chooser.randrange(2, len(seq) - 12)
        size = chooser.randrange(1, 6)
        if chooser.random() < 0.5:
            cases[pos] = (pos, seq[pos - 1 : pos + size], seq[pos - 1])
        else:
            cases[pos] = (pos, seq[pos - 1], seq[pos - 1] + seq[pos : pos + size] if chooser.random() < 0.7 else seq[pos - 1] + "".join(chooser.choice("ACGT") for _ in range(size)))
    records = [record("1", p, r, a, "0/1", vectors=False) for p, r, a in sorted(cases.values())]
    _, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    got = [(row["pos"], row["ref"], row["alts"][0]) for row in calls.iter_rows(named=True)]
    assert all(is_normalized(p, r, a) for p, r, a in got)
    # Indels written at different places in one repeat become one row: exact duplicates are dropped.
    assert {apply("1", p, r, a) for p, r, a in got} == {apply("1", p, r, a) for p, r, a in cases.values()}
    assert len(got) == len(set(got))


def test_contig_start_keeps_the_base_after() -> None:
    # A deletion of the first A of AAAA… cannot move before position 1.
    assert _left_align(3, "AA", "A", "AA", 1) == (1, "AA", "A", True)


def test_window_too_short_asks_for_more() -> None:
    assert _left_align(10, "TT", "T", "TT", 8) is None


def test_long_repeat_reads_the_reference_again(tmp_path: Path) -> None:
    # A run longer than the first window: the deletion must travel its whole length.
    run = "T" * 500
    seq = "G" * 3 + "C" + run + "GATC" * 20
    store = ReferenceStore.open(tmp_path / "store")
    fasta = tmp_path / "long.fa"
    fasta.write_text(f">1\n{seq}\n")
    store.add_fasta(fasta)
    end = 4 + len(run)
    text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            f"##contig=<ID=1,length={len(seq)}>",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            f"1\t{end - 1}\t.\tTT\tT\t50\tPASS\t.\tGT\t0/1",
        ]
    )
    (tmp_path / "long.vcf").write_text(text + "\n")
    out = normalize(resolve_ref(load(tmp_path / "long.vcf"), store), store)
    assert out.calls.collect().select("pos", "ref", "alts").row(0) == (4, "CT", ["C"])


def test_unchecked_rows_pass_through_counted(tmp_path: Path, store: ReferenceStore) -> None:
    end_t = HOMOPOLYMER + 7
    wrong = other(at("1", end_t - 1)) + at("1", end_t)
    records = [record("1", end_t - 1, wrong, wrong[0], "0/1"), record("2", 100, at("2", 100, 2), at("2", 100), "0/1")]
    parent, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    mismatch = calls.filter(pl.col("chrom") == "1").row(0, named=True)
    assert (mismatch["pos"], mismatch["source_ref"], mismatch["ref_check"]) == (end_t - 1, wrong, "mismatch")
    assert any("1 row(s) not left-aligned or trimmed: REF differs from the reference" in w for w in out.warnings)


def test_a_dash_alt_is_counted_not_passed_over(tmp_path: Path, store: ReferenceStore) -> None:
    # Some uploads write a deletion as ALT `-` (F12). It is no VCF allele, so nothing aligns it.
    _, out = normalized(tmp_path, store, [record("1", 30, at("1", 30), "-", "0/1", vectors=False)])
    assert out.calls.collect()["alts"].to_list() == [["-"]]
    assert any("1 row(s) not left-aligned or trimmed: an ALT that is not a VCF allele" in w for w in out.warnings)


# --- dedup, determinism, the derived facet ------------------------------------------------------


def test_exact_duplicates_drop_and_conflicts_stay(tmp_path: Path, store: ReferenceStore) -> None:
    end_t = HOMOPOLYMER + 7
    same = record("1", 40, at("1", 40), other(at("1", 40)), "0/1")
    records = [
        same,
        same,
        # The same deletion written at two places in the run: one allele after left-alignment, two genotypes.
        record("1", end_t - 1, at("1", end_t - 1, 2), at("1", end_t - 1), "0/1", vectors=False),
        record("1", end_t - 3, at("1", end_t - 3, 2), at("1", end_t - 3), "1/1", vectors=False),
    ]
    parent, out = normalized(tmp_path, store, records)
    calls = out.calls.collect()
    assert parent.calls.collect().height == 4
    assert calls.filter(pl.col("pos") == 40).height == 1
    conflict = calls.filter(pl.col("pos") == HOMOPOLYMER - 1)
    assert sorted(conflict["gt"].to_list()) == [[0, 1], [1, 1]]
    assert any(w == "1 exact duplicate row(s) dropped" for w in out.warnings)
    assert any(w.startswith("2 row(s) at 1 position(s) state the same allele and differ otherwise") for w in out.warnings)


def mixed_records() -> list[str]:
    end_t = HOMOPOLYMER + 7
    ref = at("1", 90)
    return [
        record("1", 90, ref, ",".join(b for b in "ACGT" if b != ref), "1/3"),
        record("1", end_t - 1, at("1", end_t - 1, 2), at("1", end_t - 1), "0/1"),
        record("2", 300, at("2", 300, 3), at("2", 300) + "," + at("2", 300, 2), "1/2"),
        record("2", 200, at("2", 200), other(at("2", 200)), "1/1"),
    ]


def test_derived_facet_and_parent_untouched(tmp_path: Path, store: ReferenceStore) -> None:
    parent, out = normalized(tmp_path, store, mixed_records())
    before = parent.calls.collect()
    out.calls.collect()
    assert parent.calls.collect().equals(before)
    facet, source = out.facets.collect().row(0, named=True), parent.facets.collect().row(0, named=True)
    assert facet["derived_from"] == source["facet_id"] and facet["relation"] == RELATION
    assert facet["facet_id"] != source["facet_id"]
    assert set(out.calls.collect()["facet_id"]) == {facet["facet_id"]}


def test_two_runs_and_a_shuffled_input_give_identical_frames(tmp_path: Path, store: ReferenceStore) -> None:
    records = mixed_records()
    parent, first = normalized(tmp_path, store, records)
    expected = first.calls.collect()
    again = normalize(parent, store)
    assert again.calls.collect().equals(expected)
    assert again.facets.collect().equals(first.facets.collect())
    shuffled = records[:]
    random.Random(3).shuffle(shuffled)
    other_dir = tmp_path / "again"
    other_dir.mkdir()
    # Other bytes are another source, so another facet_id; every other value is the same.
    _, second = normalized(other_dir, store, shuffled)
    got = second.calls.collect()
    assert got.drop("facet_id").equals(expected.drop("facet_id"))
    keys = [(0 if row["chrom"] == "1" else 1, row["pos"]) for row in got.iter_rows(named=True)]
    assert keys == sorted(keys)


def test_unresolved_facet_is_refused(tmp_path: Path, store: ReferenceStore) -> None:
    loaded = load(write_vcf_text(tmp_path / "in.vcf", mixed_records()))
    with pytest.raises(ValueError, match="resolve_ref first"):
        normalize(loaded, store)


def test_a_store_without_the_collection_is_refused(tmp_path: Path, store: ReferenceStore) -> None:
    parent = resolve_ref(load(write_vcf_text(tmp_path / "in.vcf", mixed_records())), store)
    with pytest.raises(ReferenceMismatchError, match="does not hold collection"):
        normalize(parent, ReferenceStore.open(tmp_path / "empty"))


def test_pre_split_rows_are_left_for_rejoin_split(tmp_path: Path, store: ReferenceStore) -> None:
    end_t = HOMOPOLYMER + 7
    ref = at("1", end_t - 1, 2)
    alt = other(ref[0])
    # bcftools -m- output: one site written as two rows, the other ALT as 0; one row is an indel.
    records = [record("1", end_t - 1, ref, ref[0], "1/0", vectors=False), record("1", end_t - 1, ref[0], alt, "0/1", vectors=False)]
    parent, out = normalized(tmp_path, store, records)
    assert parent.calls.collect()["split_signal"].null_count() == 0
    assert out.calls.collect().select("pos", "source_ref").rows() == parent.calls.collect().select("pos", "source_ref").rows()
    assert any("2 row(s) flagged as one record of a pre-split site (RM270) left as written" in w for w in out.warnings)


def test_written_output_reloads_as_held(tmp_path: Path, store: ReferenceStore) -> None:
    _, out = normalized(tmp_path, store, mixed_records())
    written = write_vcf(out, tmp_path / "out.vcf.gz")
    reloaded = load(written.path).calls.collect()
    held = out.calls.collect()
    view = lambda frame: [(r["chrom"], r["pos"], r["alts"], r["gt"]) for r in frame.iter_rows(named=True)]  # noqa: E731
    assert view(reloaded) == view(held)
    assert reloaded["source_ref"].to_list() == held.select(pl.coalesce("ref", "source_ref"))["ref"].to_list()


# --- against bcftools norm (htslib's, bundled in pysam) -----------------------------------------


def test_agrees_with_bcftools_norm(tmp_path: Path, store: ReferenceStore) -> None:
    pysam = pytest.importorskip("pysam")
    bcftools = pytest.importorskip("pysam.bcftools")
    chooser = random.Random(17)
    seq = SEQ["1"]
    records = mixed_records()
    for pos in sorted(chooser.sample(range(2, len(seq) - 12), 60)):
        size = chooser.randrange(1, 5)
        deletion = seq[pos - 1 : pos + size]
        insertion = seq[pos - 1] + seq[pos : pos + size]
        snv = other(seq[pos - 1])
        kind = chooser.randrange(4)
        if kind == 0:
            records.append(record("1", pos, deletion, deletion[0], "0/1"))
        elif kind == 1:
            records.append(record("1", pos, seq[pos - 1], insertion, "1/1"))
        elif kind == 2:
            records.append(record("1", pos, deletion, f"{deletion[0]},{snv + deletion[1:]}", chooser.choice(["1/2", "0/2", "2/2", "0/1"])))
        else:
            records.append(record("1", pos, seq[pos - 1], f"{snv},{insertion}", chooser.choice(["1/2", "0|1", "2|0"])))
    lines = sorted(set(records), key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    parent = resolve_ref(load(write_vcf_text(tmp_path / "in.vcf", lines)), store)
    fasta = write_fasta(tmp_path / "plain.fa")
    pysam.faidx(str(fasta))
    exported = write_vcf(parent, tmp_path / "parent.vcf", compress=False)
    theirs_path = tmp_path / "bcftools.vcf"
    bcftools.norm("-m-any", "--multi-overlaps", ".", "-f", str(fasta), "--no-version", "-o", str(theirs_path), str(exported.path), catch_stdout=False)
    result = normalize(parent, store)
    moved = next(int(w.split()[0]) for w in result.warnings if "left-aligned to another position" in w)
    assert moved >= 10  # the comparison exercises left-alignment, not only splitting
    ours = write_vcf(result, tmp_path / "ours.vcf", compress=False)

    def sites(path: Path) -> list[tuple]:
        frame = load(path).calls.collect()
        return sorted(
            (r["chrom"], r["pos"], r["source_ref"].upper(), tuple(r["alts"]), tuple(a or "." for a in r["gt_alleles"]), tuple(r["ad"] or ()))
            for r in frame.iter_rows(named=True)
        )

    theirs, mine = sites(theirs_path), sites(ours.path)
    # bcftools keeps indels that converge on one allele as two rows; normalize drops exact duplicates.
    assert set(mine) == set(theirs)
    assert len(mine) == len(set(mine))
    # GL: ours is carried only where the slice holds the called genotype; there it equals bcftools' PL slice.
    gl = lambda path: {(r["chrom"], r["pos"], tuple(r["alts"])): r["gl"] for r in load(path).calls.collect().iter_rows(named=True)}  # noqa: E731
    their_gl, my_gl = gl(theirs_path), gl(ours.path)
    carried = {key: value for key, value in my_gl.items() if value is not None}
    # Exactly: PL -> GL lands on the Float32 nearest the decimal value (F11), the one bcftools' GL reads as.
    differ = {key: (value, their_gl[key]) for key, value in carried.items() if their_gl[key] != value}
    assert carried and not differ, differ


def test_rejoin_split_takes_back_what_did_not_move(tmp_path: Path, store: ReferenceStore) -> None:
    _, out = normalized(tmp_path, store, mixed_records())
    rejoined = rejoin_split(out).calls.collect()
    # The SNV site at 90 was split in place and joins again; the indel site at chrom 2 moved, and its
    # rows hold no source REF to join on, so they stay as normalize left them.
    site = rejoined.filter((pl.col("chrom") == "1") & (pl.col("pos") == 90))
    assert site.height == 1 and len(site["alts"].item()) == 3
    assert rejoined.filter(pl.col("chrom") == "2")["split_signal"].drop_nulls().len() >= 1


def test_real_upload_agrees_with_bcftools_on_one_contig(pgp_participants: Path, test_reference_dir: Path, tmp_path: Path) -> None:
    """Opt-in: a PGP-HMS VCF upload, chromosome 22, normalized here and by bcftools. Every row where the
    two differ is one the loader flagged as part of a pre-split site, which normalize leaves for
    rejoin_split. Aggregates only."""
    pysam = pytest.importorskip("pysam")
    bcftools = pytest.importorskip("pysam.bcftools")
    store = ReferenceStore.open(test_reference_dir)
    sizes = sorted(pgp_participants.glob("*/upload-*/*.vcf*"), key=lambda p: p.stat().st_size)
    parent = None
    for path in [p for p in sizes if p.stat().st_size > 5_000_000][:20]:
        if sniff(path).format != "vcf":
            continue
        try:
            parent = resolve_ref(load(path), store)
        except ValueError:  # ReferenceMismatchError included: not an assembly the store holds
            continue
        if parent.calls.filter(pl.col("chrom") == "22").select(pl.len()).collect().item() > 1000:
            break
        parent = None
    if parent is None:
        pytest.skip("no PGP-HMS VCF upload that resolves against the store")
    # bcftools refuses a file holding an ALT that is no VCF allele (F12's `-`), so those rows stay out of both.
    allele = r"^([ACGTNacgtn]+|\*|<[^>]+>)$"
    valid = pl.col("alts").list.eval(pl.element().str.contains(allele)).list.all()
    one = Loaded(facets=parent.facets, calls=parent.calls.filter((pl.col("chrom") == "22") & valid), blocks=pl.LazyFrame(schema=BLOCKS))
    digest = one.facets.collect()["ref_seqcol_digest"].item()
    contig = next(c for c in store.contigs(digest) if c.name in ("22", "chr22"))
    sequence = store.fetch(contig, [(0, contig.length)])[0]
    fasta = tmp_path / "22.fa"
    fasta.write_text(">22\n" + "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
    pysam.faidx(str(fasta))
    exported = write_vcf(one, tmp_path / "parent.vcf", compress=False, resolved_only=True)
    theirs_path = tmp_path / "theirs.vcf"
    bcftools.norm("-m-any", "--multi-overlaps", ".", "-f", str(fasta), "--no-version", "-o", str(theirs_path), str(exported.path), catch_stdout=False)
    result = normalize(one, store)
    ours = write_vcf(result, tmp_path / "ours.vcf", compress=False, resolved_only=True)

    def key(r: dict) -> tuple:
        return (r["pos"], (r["ref"] or r["source_ref"]).upper(), tuple(r["alts"]), tuple(a or "." for a in r["gt_alleles"] or ()))

    def sites(path: Path) -> Counter:
        return Counter(key(r) for r in load(path).calls.collect().iter_rows(named=True))

    theirs, mine = sites(theirs_path), sites(ours.path)
    held = {key(r) for r in result.calls.filter(pl.col("split_signal").is_not_null()).collect().iter_rows(named=True)}
    unexplained = [k for k in (mine - theirs) if k not in held]
    agree = sum((theirs & mine).values())
    assert agree > 1000
    assert not unexplained, f"{len(unexplained)} row(s) differ from bcftools outside pre-split sites; {agree} agree"
