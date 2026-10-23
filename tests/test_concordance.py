"""Concordance between two facets (RM277): synthetic VCF pairs on a small refget store.

Each case writes two one-sample VCFs over the same sequence, reads REF (`resolve_ref`), and compares.
Expected classes come from the fixture's own genotypes; graded values are recomputed here from the
vectors written, by enumerating the genotypes in plain Python, never copied from a run.
"""

import itertools
import math
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceMismatchError, ReferenceStore, load, mint_vrs_ids, normalize, resolve_ref, sniff
from just_vcf.concordance import SITES, SUMMARY, Agreement, concordance, summarize
from just_vcf.loading import Loaded

RNG = random.Random(277)
REPEAT_AT = 501
"""1-based start of a `CACACACACA` run in contig 1, where a deletion can be written at several positions."""


def _contig_one() -> str:
    def rnd(n: int) -> str:
        return "".join(RNG.choice("ACGT") for _ in range(n))

    head = rnd(REPEAT_AT - 1)
    # The base before the repeat and the base after it are not C or A, so the run is exactly ten bases.
    return head[:-1] + "G" + "CACACACACA" + "T" + rnd(1488)


SEQ = {"1": _contig_one(), "2": "".join(RNG.choice("ACGT") for _ in range(1500)), "X": "".join(RNG.choice("ACGT") for _ in range(800))}


def base(chrom: str, pos: int, length: int = 1) -> str:
    return SEQ[chrom][pos - 1 : pos - 1 + length]


def other(b: str, skip: str = "") -> str:
    return next(x for x in "ACGT" if x != b and x not in skip)


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    lines = []
    for name, seq in seqs.items():
        lines.append(f">{name}")
        lines += [seq[i : i + 60] for i in range(0, len(seq), 60)]
    path.write_text("\n".join(lines) + "\n")
    return path


Record = tuple[str, int, str, str, str, dict[str, str]]
"""`(chrom, pos, ref, alt, gt, extra FORMAT fields)`; INFO `END` goes in extra as `END`."""


def write_vcf(path: Path, records: list[Record], contigs: dict[str, int] | None = None) -> Path:
    contigs = contigs or {name: len(seq) for name, seq in SEQ.items()}
    keys = sorted({key for *_, extra in records for key in extra if key != "END"})
    meta = ["##fileformat=VCFv4.3", *(f"##contig=<ID={c},length={n}>" for c, n in contigs.items())]
    meta += [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
        '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posteriors">',
    ]
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = []
    for chrom, pos, ref, alt, gt, extra in records:
        info = f"END={extra['END']}" if "END" in extra else "."
        fields = ["GT", *keys]
        values = [gt, *(extra.get(key, ".") for key in keys)]
        body.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\t{':'.join(fields)}\t{':'.join(values)}")
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa", SEQ))
    return s


def facet(tmp_path: Path, name: str, records: list[Record], store: ReferenceStore, mint: bool = False) -> Loaded:
    resolved = resolve_ref(load(write_vcf(tmp_path / f"{name}.vcf", records)), store)
    return mint_vrs_ids(resolved, store) if mint else resolved


def snv(chrom: str, pos: int, gt: str, alt: str | None = None, **extra: str) -> Record:
    ref = base(chrom, pos)
    return (chrom, pos, ref, alt or other(ref), gt, extra)


def test_genotypes_compare_as_copies_of_the_allele(tmp_path: Path, store: ReferenceStore) -> None:
    # (a's genotype, b's genotype, the agreement the two state), one site each.
    cases = [
        ("0/1", "1/0", Agreement.MATCH),  # unordered: phase and slot order are not compared
        ("0|1", "1|0", Agreement.MATCH),
        ("1/1", "1/1", Agreement.MATCH),
        ("0/0", "0/0", Agreement.MATCH),
        ("0/1", "1/1", Agreement.MISMATCH),
        ("0/0", "0/1", Agreement.MISMATCH),
        ("./.", "0/1", Agreement.UNCALLED),  # the source looked and could not tell
        ("1/.", "1/1", Agreement.UNCALLED),
        ("1/.", "0/0", Agreement.MISMATCH),  # the called haplotype already contradicts
        ("1", "1/1", Agreement.PLOIDY_DIFFERS),
    ]
    positions = [20 + 30 * i for i in range(len(cases))]
    a = facet(tmp_path, "a", [snv("1", p, ga) for p, (ga, _, _) in zip(positions, cases)], store)
    b = facet(tmp_path, "b", [snv("1", p, gb) for p, (_, gb, _) in zip(positions, cases)], store)
    result = concordance(a, b)
    sites = result.sites.collect()
    assert sites.schema == SITES
    assert sites["pos"].to_list() == positions
    assert sites["agreement"].to_list() == [want.value for *_, want in cases]
    assert set(sites["key"].to_list()) == {"position"}
    assert set(sites["basis"].to_list()) == {"hard"}
    assert sites["p_same"].null_count() == sites.height
    # Copies are the allele's count in each genotype, null where a haplotype is uncalled.
    expect_a = [None if "." in ga else ga.count("1") for ga, _, _ in cases]
    assert sites["a_copies"].to_list() == expect_a


def test_a_shifted_deletion_meets_its_twin_by_vrs_id(tmp_path: Path, store: ReferenceStore) -> None:
    # One CA unit deleted from the repeat, written left-aligned in `a` and at the repeat's right end in `b`.
    left = REPEAT_AT - 1
    right = REPEAT_AT + 7
    left_ref, right_ref = base("1", left, 3), base("1", right, 3)
    seq = SEQ["1"]
    # Both spellings give the same haplotype: that is what makes them one allele.
    assert seq[: left] + seq[left + 2 :] == seq[: right] + seq[right + 2 :]
    a_rec = [("1", left, left_ref, left_ref[0], "0/1", {})]
    b_rec = [("1", right, right_ref, right_ref[0], "0/1", {})]

    unminted = concordance(facet(tmp_path, "a0", a_rec, store), facet(tmp_path, "b0", b_rec, store)).sites.collect()
    assert unminted["agreement"].to_list() == ["only_a", "only_b"]

    minted = concordance(facet(tmp_path, "a", a_rec, store, mint=True), facet(tmp_path, "b", b_rec, store, mint=True)).sites.collect()
    assert minted.height == 1
    row = minted.row(0, named=True)
    assert (row["key"], row["agreement"], row["a_pos"], row["b_pos"]) == ("vrs", "match", left, right)
    assert row["variant_class"] == "indel"
    assert row["vrs_allele_id"] is not None


def test_an_allele_the_other_facet_does_not_state_is_read_against_its_row_there(tmp_path: Path, store: ReferenceStore) -> None:
    ref = base("1", 300)
    g, t = sorted({other(ref), other(ref, skip=other(ref))})
    mono = base("1", 400)
    a = facet(tmp_path, "a", [("1", 300, ref, f"{g},{t}", "1/2", {}), snv("1", 400, "0/1"), snv("1", 500, "0/1"), snv("1", 600, "0/1")], store)
    b = facet(
        tmp_path,
        "b",
        [("1", 300, ref, g, "0/1", {}), ("1", 400, mono, ".", "0/0", {}), ("1", 500, base("1", 500), f"{other(base('1', 500))}", "./.", {})],
        store,
    )
    sites = concordance(a, b).sites.collect()
    by = {(row["pos"], row["alt"]): row for row in sites.iter_rows(named=True)}
    assert (by[(300, g)]["key"], by[(300, g)]["agreement"]) == ("position", "match")
    # `b` wrote only G at 300: its genotype holds no T, which `a` called once.
    t_row = by[(300, t)]
    assert (t_row["key"], t_row["agreement"], t_row["absence"], t_row["a_copies"], t_row["b_copies"]) == ("site", "mismatch", "other_alleles", 1, 0)
    assert by[(400, other(mono))]["agreement"] == "mismatch" and by[(400, other(mono))]["b_gt_alleles"] == [mono, mono]
    # 500: the same allele, no-called in b. 600: nothing in b at all, and b's absence says nothing.
    assert by[(500, other(base("1", 500)))]["agreement"] == "uncalled"
    lone = by[(600, other(base("1", 600)))]
    assert (lone["agreement"], lone["key"], lone["absence"], lone["basis"]) == ("only_a", None, "unknown", None)


def test_a_row_over_another_span_is_not_read_for_copies(tmp_path: Path, store: ReferenceStore) -> None:
    # `a` calls an SNV at 100; `b` wrote one row there, an MNV over two bases (a haplotype caller's
    # complex allele). Its genotype holds no allele of a one-base span, so it says nothing about the
    # SNV's copies: not a mismatch.
    ref2 = base("1", 100, 2)
    a = facet(tmp_path, "a", [snv("1", 100, "0/1")], store)
    b = facet(tmp_path, "b", [("1", 100, ref2, other(ref2[0]) + other(ref2[1]), "0/1", {})], store)
    rows = {row["alt"]: row for row in concordance(a, b).sites.collect().iter_rows(named=True)}
    lone = rows[other(base("1", 100))]
    assert (lone["agreement"], lone["key"], lone["absence"], lone["b_copies"], lone["b_gt_alleles"]) == ("only_a", None, "other_alleles", None, None)


def test_a_silent_facet_says_what_its_absence_means(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [snv("1", 100, "0/1"), snv("1", 200, "0/1"), snv("1", 900, "0/1"), snv("2", 50, "0/1")], store)
    gvcf = [
        ("1", 90, base("1", 90), "<*>", "0/0", {"END": "150"}),
        ("1", 190, base("1", 190), "<*>", "./.", {"END": "210"}),
        snv("1", 1200, "0/1"),
    ]
    b_gvcf = facet(tmp_path, "b", gvcf, store)
    by_pos = {row["pos"]: row["absence"] for row in concordance(a, b_gvcf).sites.collect().iter_rows(named=True) if row["agreement"] == "only_a"}
    assert by_pos[100] == "ref_block" and by_pos[200] == "no_call_block"
    assert by_pos[900] == "unknown"  # a gVCF decides per position, and no block covers 900

    # The same silence in a whole-genome callset without blocks reads as reference for a class it called.
    plain = facet(tmp_path, "plain", [snv("1", 1200, "0/1")], store)
    wgs = Loaded(facets=plain.facets.with_columns(pl.lit("weak_reference").cast(plain.facets.collect_schema()["absence_means"]).alias("absence_means")), calls=plain.calls, blocks=plain.blocks)
    read = {row["pos"]: row["absence"] for row in concordance(a, wgs).sites.collect().iter_rows(named=True) if row["agreement"] == "only_a"}
    assert read[100] == read[900] == "reference"
    # Contig 2 is in b's header, so its silence there counts too; a contig it never declared would not.
    assert read[50] == "reference"
    undeclared = facet(tmp_path, "nohead", [snv("1", 1200, "0/1")], store)
    calls_only = Loaded(
        facets=wgs.facets.with_columns(pl.lit(None).cast(wgs.facets.collect_schema()["header_contigs"]).alias("header_contigs")),
        calls=undeclared.calls,
        blocks=undeclared.blocks,
    )
    result = concordance(a, calls_only)
    read = {row["pos"]: row["absence"] for row in result.sites.collect().iter_rows(named=True) if row["agreement"] == "only_a"}
    assert read[50] == "unknown" and read[100] == "reference"
    assert any("neither declared nor called" in line for line in result.warnings)


def _flat(pl_values: list[int]) -> list[float]:
    weights = [10 ** (-v / 10) for v in pl_values]
    total = sum(weights)
    return [w / total for w in weights]


def _diploid_copies(n: int, allele: int) -> list[int]:
    """Copies of `allele` in each diploid genotype of `n` alleles, in VCF order."""
    order = sorted(itertools.combinations_with_replacement(range(n), 2), key=lambda jk: (jk[1], jk[0]))
    return [(j == allele) + (k == allele) for j, k in order]


def test_graded_comparison_where_both_rows_carry_one_kind_of_vector(tmp_path: Path, store: ReferenceStore) -> None:
    gp_a, gp_b = [0.1, 0.8, 0.1], [0.2, 0.6, 0.2]
    pl_a, pl_b = [20, 0, 200], [0, 10, 90]
    multi_ref = base("1", 700)
    g, t = sorted({other(multi_ref), other(multi_ref, skip=other(multi_ref))})
    pl_multi = [40, 0, 60, 15, 30, 90]  # over REF, G, T
    a = facet(
        tmp_path,
        "a",
        [
            snv("1", 100, "0/1", GP=",".join(map(str, gp_a))),
            snv("1", 200, "0/1", PL=",".join(map(str, pl_a))),
            snv("1", 300, "0/1", PL="20,0,200"),  # likelihoods here, posteriors in b: never compared
            ("1", 700, multi_ref, f"{g},{t}", "0/1", {"PL": ",".join(map(str, pl_multi))}),
        ],
        store,
    )
    b = facet(
        tmp_path,
        "b",
        [
            snv("1", 100, "0/1", GP=",".join(map(str, gp_b))),
            snv("1", 200, "0/1", PL=",".join(map(str, pl_b))),
            snv("1", 300, "0/1", GP="0.1,0.8,0.1"),
            ("1", 700, multi_ref, g, "0/1", {"PL": "30,0,50"}),
        ],
        store,
    )
    sites = {row["pos"]: row for row in concordance(a, b).sites.collect().iter_rows(named=True) if row["alt"] != t}
    assert sites[100]["basis"] == "gp"
    assert sites[100]["p_same"] == pytest.approx(sum(x * y for x, y in zip(gp_a, gp_b)))
    assert sites[200]["basis"] == "gl_flat_prior"
    assert sites[200]["p_same"] == pytest.approx(sum(x * y for x, y in zip(_flat(pl_a), _flat(pl_b))), rel=1e-6)
    assert (sites[300]["basis"], sites[300]["p_same"]) == ("hard", None)
    # The multi-allelic row's likelihoods summed over the genotypes with each number of copies of G.
    copies = _diploid_copies(3, 1)
    pa = [sum(p for p, c in zip(_flat(pl_multi), copies) if c == k) for k in range(3)]
    pb = _flat([30, 0, 50])
    assert sites[700]["basis"] == "gl_flat_prior"
    assert sites[700]["p_same"] == pytest.approx(sum(x * y for x, y in zip(pa, pb)), rel=1e-6)
    assert all(0.0 <= row["p_same"] <= 1.0 for row in sites.values() if row["p_same"] is not None)
    assert math.isclose(sum(pa), 1.0)


def test_summary_counts_every_allele_and_excluded_alleles_are_counted(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(
        tmp_path,
        "a",
        [
            snv("1", 100, "0/1"),
            snv("1", 200, "1/1"),
            ("1", 300, base("1", 300), "<DEL>", "0/1", {}),
            ("1", 400, base("1", 400), f"{other(base('1', 400))},*", "1/2", {}),
            ("1", 450, base("1", 450, 3), base("1", 450), "0/1", {}),
        ],
        store,
    )
    b = facet(tmp_path, "b", [snv("1", 100, "0/1"), snv("1", 200, "0/1"), ("1", 300, base("1", 300), f"{other(base('1', 300))},<*>", "0/1", {})], store)
    result = concordance(a, b)
    sites = result.sites.collect()
    assert result.summary.schema == SUMMARY
    assert int(result.summary["alleles"].sum()) == sites.height
    assert result.summary.equals(summarize(sites))
    excluded = {(row["facet"], row["reason"]): row["alleles"] for row in result.excluded.iter_rows(named=True)}
    assert excluded == {("a", "symbolic"): 1, ("a", "spanning_deletion"): 1, ("b", "unspecified"): 1}
    by_class = result.summary.group_by("variant_class").agg(pl.col("alleles").sum())
    assert dict(by_class.iter_rows()) == dict(sites.group_by("variant_class").len().iter_rows())


def test_region_keeps_alleles_inside_its_closed_intervals(tmp_path: Path, store: ReferenceStore) -> None:
    positions = list(range(50, 1500, 50))
    a = facet(tmp_path, "a", [snv("1", p, "0/1") for p in positions], store)
    b = facet(tmp_path, "b", [snv("1", p, "0/1") for p in positions[::2]], store)
    region = pl.DataFrame({"chrom": ["chr1", "1", "1"], "pos": [100, 140, 1000], "end_pos": [150, 400, 1000]})
    inside = {p for p in positions if 100 <= p <= 400 or p == 1000}
    sites = concordance(a, b, region=region).sites.collect()
    assert set(sites["pos"].to_list()) == inside
    assert set(concordance(a, b).sites.collect()["pos"].to_list()) == set(positions)


def test_the_result_is_deterministic_and_in_canonical_order(tmp_path: Path, store: ReferenceStore) -> None:
    records = [snv(c, p, "0/1") for c in ("X", "2", "1") for p in (700, 30, 300)]
    a = facet(tmp_path, "a", records, store)
    b = facet(tmp_path, "b", records[::-2], store)
    first = concordance(a, b).sites.collect()
    assert first.equals(concordance(a, b).sites.collect())
    keys = list(zip(first["chrom"].to_list(), first["pos"].to_list()))
    rank = {"1": 0, "2": 1, "X": 22}
    assert keys == sorted(keys, key=lambda k: (rank[k[0]], k[1]))


def test_facets_on_different_assemblies_are_refused(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [snv("1", 100, "0/1")], store)
    # Another collection whose contig 1 is another sequence of the same length: another assembly.
    changed = dict(SEQ, **{"1": SEQ["1"][:99] + other(SEQ["1"][99]) + SEQ["1"][100:]})
    store.add_fasta(write_fasta(tmp_path / "other.fa", changed))
    b_loaded = load(write_vcf(tmp_path / "b.vcf", [snv("1", 200, "0/1")]))
    digest_b = next(d for d in store.collections() if d != a.facets.collect()["ref_seqcol_digest"][0])
    b = resolve_ref(b_loaded, store, collection=digest_b)
    with pytest.raises(ReferenceMismatchError, match="lift one"):
        concordance(a, b, reference=store)


def test_two_collections_sharing_sequences_compare_the_shared_contigs(tmp_path: Path, store: ReferenceStore) -> None:
    a = facet(tmp_path, "a", [snv("1", 100, "0/1"), snv("2", 100, "0/1")], store)
    # UCSC-style names over the same contig 1 and X, contig 2 left out: one sequence, another collection.
    ucsc = {"chr1": SEQ["1"], "chrX": SEQ["X"]}
    digest = store.add_fasta(write_fasta(tmp_path / "ucsc.fa", ucsc))
    b_loaded = load(write_vcf(tmp_path / "b.vcf", [("chr1", 100, base("1", 100), other(base("1", 100)), "0/1", {})], contigs={k: len(v) for k, v in ucsc.items()}))
    b = resolve_ref(b_loaded, store, collection=digest)
    result = concordance(a, b, reference=store)
    sites = result.sites.collect()
    assert sites.select("chrom", "pos", "agreement").rows() == [("1", 100, "match")]
    assert dict(((r["facet"], r["reason"]), r["alleles"]) for r in result.excluded.iter_rows(named=True)) == {("a", "contig_not_shared"): 1}


def test_unread_ref_and_several_facets_are_refused(tmp_path: Path, store: ReferenceStore) -> None:
    raw = load(write_vcf(tmp_path / "a.vcf", [snv("1", 100, "0/1")]))
    read = resolve_ref(raw, store)
    with pytest.raises(ValueError, match="resolve_ref"):
        concordance(raw, read)
    two = Loaded(facets=pl.concat([read.facets, raw.facets]), calls=read.calls, blocks=read.blocks)
    with pytest.raises(ValueError, match="one facet per side"):
        concordance(two, read)


def test_real_upload_against_its_own_normalized_facet(pgp_participants: Path, test_reference_dir: Path) -> None:
    """Opt-in: a PGP-HMS VCF upload, chromosome 22, against `normalize` of itself, both minted. They
    state the same genotypes, so nothing may disagree: no mismatch, no ploidy difference, no allele
    on one side only, and `uncalled` only where the parent row had several ALTs (the split writes `.`
    for the other one) or an uncalled haplotype. Aggregates only."""
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
    one = Loaded(facets=parent.facets, calls=parent.calls.filter(pl.col("chrom") == "22"), blocks=parent.blocks.filter(pl.col("chrom") == "22"))
    a = mint_vrs_ids(one, store)
    b = mint_vrs_ids(normalize(one, store), store)
    result = concordance(a, b)
    sites = result.sites.collect()
    counts = dict(sites.group_by("agreement").len().iter_rows())
    assert counts.get("match", 0) > 1000
    assert not {"mismatch", "ploidy_differs", "only_a", "only_b"} & set(counts), counts
    multi = (
        one.calls.filter((pl.col("alts").list.len() > 1) | pl.col("gt_alleles").list.eval(pl.element().is_null()).list.any().fill_null(True))
        .select("pos")
        .collect()["pos"]
    )
    assert set(sites.filter(pl.col("agreement") == "uncalled")["a_pos"].to_list()) <= set(multi.to_list())
    assert int(result.summary["alleles"].sum()) == sites.height
