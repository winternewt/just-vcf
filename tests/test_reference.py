"""REF read from the reference sequence (RM3): `resolve_ref` against small synthetic FASTAs."""

import gzip
import hashlib
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceMismatchError, ReferenceStore, load, reference_dir, resolve_ref
from just_vcf.loading import Loaded
from just_vcf.reference import _span
from just_vcf.schema import CALLS, FACETS

LENGTHS = {"1": 3000, "2": 2000, "X": 1500, "MT": 400}
"""Synthetic contigs. Short, so the assembly is told by stated lengths or an asserted label."""


def sequences(seed: int = 7) -> dict[str, str]:
    rng = random.Random(seed)
    return {name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in LENGTHS.items()}


SEQ = sequences()


def write_fasta(path: Path, seqs: dict[str, str], names: dict[str, str] | None = None, lower: bool = False) -> Path:
    names = names or {}
    lines = []
    for name, seq in seqs.items():
        lines.append(f">{names.get(name, name)} synthetic")
        body = seq.lower() if lower else seq
        lines += [body[i : i + 60] for i in range(0, len(body), 60)]
    text = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


UCSC = {"1": "chr1", "2": "chr2", "X": "chrX", "MT": "chrM"}


def contig_lines(lengths: dict[str, int], names: dict[str, str] | None = None, md5: dict[str, str] | None = None) -> list[str]:
    names = names or {}
    lines = []
    for name, length in lengths.items():
        extra = f",md5={md5[name]}" if md5 and name in md5 else ""
        lines.append(f"##contig=<ID={names.get(name, name)},length={length}{extra}>")
    return lines


def vcf(tmp_path: Path, records: list[tuple[str, int, str, str, str]], contigs: list[str], name: str = "calls.vcf") -> Path:
    """A one-sample VCF of `(chrom, pos, ref, alt, gt)` records with the given `##contig` lines."""
    meta = ["##fileformat=VCFv4.2", *contigs, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{g}" for c, p, r, a, g in records]
    path = tmp_path / name
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def base(chrom: str, pos: int, length: int = 1) -> str:
    """The synthetic reference's bases at 1-based `pos`, the test's own reading of the sequence."""
    return SEQ[chrom][pos - 1 : pos - 1 + length]


def other(b: str) -> str:
    return next(x for x in "ACGT" if x != b)


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(write_fasta(tmp_path / "ref.fa", SEQ))
    return s


def frames(loaded: Loaded) -> tuple[pl.DataFrame, pl.DataFrame]:
    return pl.collect_all([loaded.facets, loaded.calls])


def assert_gt_indexes_ref(calls: pl.DataFrame) -> None:
    """For array rows: `gt_alleles` is `[ref] + alts` indexed by `gt`."""
    for row in calls.filter(pl.col("source_ref").is_null() & pl.col("gt").is_not_null()).iter_rows(named=True):
        alleles = [row["ref"], *row["alts"]]
        assert [None if i is None else alleles[i] for i in row["gt"]] == row["gt_alleles"]


def test_span_converts_one_based_pos_to_half_open() -> None:
    frame = pl.DataFrame({"pos": [1, 100], "n": [1, 3]})
    start, end = _span(pl.col("pos"), pl.col("n"))
    got = frame.select(start.alias("s"), end.alias("e"))
    assert got.rows() == [(0, 1), (99, 102)]


def test_vcf_ref_is_read_and_checked(tmp_path: Path, store: ReferenceStore) -> None:
    good_snv = ("1", 10, base("1", 10), other(base("1", 10)), "0/1")
    good_del = ("1", 50, base("1", 50, 4), base("1", 50), "1/1")
    wrong = ("2", 70, other(base("2", 70)), base("2", 70), "0/1")
    records = [good_snv, good_del, wrong]
    path = vcf(tmp_path, records, contig_lines({"1": 3000, "2": 2000}))
    loaded = load(path)
    resolved = resolve_ref(loaded, store)
    facets, calls = frames(resolved)
    parent_facets, parent_calls = frames(loaded)

    assert facets["facet_id"].to_list() == parent_facets["facet_id"].to_list()
    assert facets["ref_seqcol_digest"].to_list() == store.collections()
    assert parent_calls["ref"].null_count() == parent_calls.height  # the input is not changed
    assert calls["source_ref"].to_list() == parent_calls["source_ref"].to_list()
    expected_ref = [base(c, p, len(r)) for c, p, r, _, _ in records]
    assert calls["ref"].to_list() == expected_ref
    expected_check = ["match" if r == base(c, p, len(r)) else "mismatch" for c, p, r, _, _ in records]
    assert calls["ref_check"].cast(pl.String).to_list() == expected_check
    assert calls.drop("ref", "ref_check").equals(parent_calls.drop("ref", "ref_check"))
    assert calls.schema == CALLS and facets.schema == FACETS
    assert any("2 match, 1 mismatch" in w for w in resolved.warnings)


def test_ucsc_named_fasta_serves_an_ensembl_named_vcf_and_back(tmp_path: Path) -> None:
    ucsc_store = ReferenceStore.open(tmp_path / "ucsc")
    ucsc_store.add_fasta(write_fasta(tmp_path / "ucsc.fa", SEQ, names=UCSC))
    ensembl_store = ReferenceStore.open(tmp_path / "ensembl")
    ensembl_store.add_fasta(write_fasta(tmp_path / "ensembl.fa", SEQ))
    records = [("1", 5, base("1", 5), "A" if base("1", 5) != "A" else "C", "0/1"), ("MT", 33, base("MT", 33), "*", "1")]
    plain = vcf(tmp_path, records, contig_lines({"1": 3000, "MT": 400}), name="plain.vcf")
    prefixed = vcf(tmp_path, [(UCSC[c], *rest) for c, *rest in records], contig_lines({"1": 3000, "MT": 400}, UCSC), name="prefixed.vcf")
    for path, reference in ((plain, ucsc_store), (prefixed, ensembl_store)):
        _, calls = frames(resolve_ref(load(path), reference))
        assert calls["chrom"].to_list() == ["1", "MT"]
        assert calls["ref"].to_list() == [base("1", 5), base("MT", 33)]
        assert set(calls["ref_check"].cast(pl.String)) == {"match"}


def test_a_primary_contig_of_another_length_refuses_the_facet(tmp_path: Path, store: ReferenceStore) -> None:
    path = vcf(tmp_path, [("1", 10, base("1", 10), "N", "0/1")], contig_lines({"1": 3001}))
    with pytest.raises(ReferenceMismatchError, match="another assembly"):
        resolve_ref(load(path), store)


def test_an_mt_of_another_length_is_left_unread_and_counted(tmp_path: Path, store: ReferenceStore) -> None:
    """hg19's chrM is 16,571 bases and GRCh37's MT 16,569: the autosomes are still read."""
    records = [("1", 10, base("1", 10), "N", "0/1"), ("MT", 20, base("MT", 20), "N", "1"), ("MT", 21, base("MT", 21), "N", "1")]
    path = vcf(tmp_path, records, contig_lines({"1": 3000, "MT": 402}))
    resolved = resolve_ref(load(path), store)
    _, calls = frames(resolved)
    assert calls.filter(pl.col("chrom") == "1")["ref"].to_list() == [base("1", 10)]
    mt = calls.filter(pl.col("chrom") == "MT")
    assert mt["ref"].null_count() == 2
    assert set(mt["ref_check"].cast(pl.String)) == {"not_run"}
    assert any(w.startswith("2 row(s) left with REF unread") and "MT (stated 402, reference 400)" in w for w in resolved.warnings)


def test_md5_is_compared_when_stated(tmp_path: Path, store: ReferenceStore) -> None:
    right = {"1": hashlib.md5(SEQ["1"].encode()).hexdigest()}
    records = [("1", 10, base("1", 10), "N", "0/1")]
    ok = vcf(tmp_path, records, contig_lines({"1": 3000}, md5=right), name="ok.vcf")
    assert frames(resolve_ref(load(ok), store))[1]["ref_check"].cast(pl.String).to_list() == ["match"]
    bad = vcf(tmp_path, records, contig_lines({"1": 3000}, md5={"1": "0" * 32}), name="bad.vcf")
    with pytest.raises(ReferenceMismatchError, match="md5_mismatch"):
        resolve_ref(load(bad), store)


def test_rows_past_the_contig_end_stay_unread(tmp_path: Path, store: ReferenceStore) -> None:
    records = [("2", 1999, base("2", 1999, 2), "A", "0/1"), ("2", 2000, base("2", 2000) + "A", "A", "0/1")]
    resolved = resolve_ref(load(vcf(tmp_path, records, contig_lines({"2": 2000}))), store)
    _, calls = frames(resolved)
    assert calls["ref"].to_list() == [base("2", 1999, 2), None]
    assert calls["ref_check"].cast(pl.String).to_list() == ["match", "not_run"]
    assert any(w.startswith("1 row(s) with a source REF reaching past their contig") for w in resolved.warnings)


def test_iupac_reference_base_is_a_counted_mismatch(tmp_path: Path) -> None:
    seqs = dict(SEQ)
    seqs["1"] = seqs["1"][:99] + "R" + seqs["1"][100:]
    s = ReferenceStore.open(tmp_path / "iupac")
    s.add_fasta(write_fasta(tmp_path / "iupac.fa", seqs))
    resolved = resolve_ref(load(vcf(tmp_path, [("1", 100, "A", "G", "0/1")], contig_lines({"1": 3000}))), s)
    _, calls = frames(resolved)
    assert calls["ref"].to_list() == ["R"] and calls["source_ref"].to_list() == ["A"]
    assert calls["ref_check"].cast(pl.String).to_list() == ["mismatch"]
    assert any("1 of them against an IUPAC ambiguity code" in w for w in resolved.warnings)


def test_soft_masked_and_gzipped_fasta_reads_uppercase(tmp_path: Path) -> None:
    s = ReferenceStore.open(tmp_path / "gz")
    s.add_fasta(write_fasta(tmp_path / "ref.fa.gz", SEQ, lower=True))
    records = [("1", 7, base("1", 7, 3).lower(), "A", "0/1")]
    _, calls = frames(resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000}))), s))
    assert calls["ref"].to_list() == [base("1", 7, 3)]
    assert calls["ref_check"].cast(pl.String).to_list() == ["match"]  # VCF REF bases are case-insensitive


def test_a_fasta_added_twice_is_one_collection_and_the_store_reopens(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "ref.fa", SEQ)
    first = ReferenceStore.open(tmp_path / "s")
    digest = first.add_fasta(fasta)
    assert first.add_fasta(fasta) == digest
    reopened = ReferenceStore.open(tmp_path / "s")
    assert reopened.collections() == [digest]
    records = [("1", 10, base("1", 10), "N", "0/1")]
    _, calls = frames(resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000}))), tmp_path / "s"))
    assert calls["ref"].to_list() == [base("1", 10)]


def test_a_fasta_path_goes_into_the_configured_store(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_REFERENCE_DIR", str(tmp_path / "configured"))
    fasta = write_fasta(tmp_path / "ref.fa", SEQ)
    records = [("1", 10, base("1", 10), "N", "0/1")]
    _, calls = frames(resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000}))), fasta))
    assert calls["ref"].to_list() == [base("1", 10)]
    assert reference_dir() == tmp_path / "configured"
    assert len(ReferenceStore.open().collections()) == 1
    # No argument: the configured store, which now holds the collection.
    _, again = frames(resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000}), name="b.vcf"))))
    assert again["ref"].to_list() == [base("1", 10)]


def test_an_empty_store_is_refused(tmp_path: Path) -> None:
    records = [("1", 10, base("1", 10), "N", "0/1")]
    with pytest.raises(ReferenceMismatchError, match="holds no collection"):
        resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000}))), ReferenceStore.open(tmp_path / "empty"))


def test_two_collections_that_fit_equally_need_a_choice(tmp_path: Path) -> None:
    s = ReferenceStore.open(tmp_path / "two")
    plain = s.add_fasta(write_fasta(tmp_path / "a.fa", SEQ))
    prefixed = s.add_fasta(write_fasta(tmp_path / "b.fa", SEQ, names=UCSC))
    path = vcf(tmp_path, [("1", 10, base("1", 10), "N", "0/1")], contig_lines({"1": 3000}))
    with pytest.raises(ReferenceMismatchError, match="fit equally"):
        resolve_ref(load(path), s)
    facets, calls = frames(resolve_ref(load(path), s, collection=prefixed))
    assert facets["ref_seqcol_digest"].to_list() == [prefixed] and plain != prefixed
    assert calls["ref"].to_list() == [base("1", 10)]


def test_the_collection_that_fits_is_chosen_over_one_that_does_not(tmp_path: Path) -> None:
    s = ReferenceStore.open(tmp_path / "mixed")
    shorter = dict(SEQ, MT=SEQ["MT"][:-2])
    s.add_fasta(write_fasta(tmp_path / "short_mt.fa", shorter))
    right = s.add_fasta(write_fasta(tmp_path / "right.fa", SEQ))
    path = vcf(tmp_path, [("1", 10, base("1", 10), "N", "0/1"), ("MT", 3, base("MT", 3), "N", "1")], contig_lines({"1": 3000, "MT": 400}))
    facets, calls = frames(resolve_ref(load(path), s))
    assert facets["ref_seqcol_digest"].to_list() == [right]
    assert calls["ref"].null_count() == 0


def test_a_reference_naming_a_contig_twice_leaves_it_unread(tmp_path: Path) -> None:
    s = ReferenceStore.open(tmp_path / "dup")
    s.add_fasta(write_fasta(tmp_path / "dup.fa", {**SEQ, "chr2": SEQ["2"]}))
    records = [("1", 10, base("1", 10), "N", "0/1"), ("2", 10, base("2", 10), "N", "0/1")]
    resolved = resolve_ref(load(vcf(tmp_path, records, contig_lines({"1": 3000, "2": 2000}))), s)
    _, calls = frames(resolved)
    assert calls["ref"].to_list() == [base("1", 10), None]
    assert any("names twice: 2 (2 and chr2)" in w for w in resolved.warnings)


def test_resolving_is_deterministic_and_order_independent(tmp_path: Path, store: ReferenceStore) -> None:
    rng = random.Random(3)
    records = [("1", p, base("1", p), other(base("1", p)), "0/1") for p in sorted(rng.sample(range(1, 3000), 40))]
    records += [("2", p, base("2", p, 2), base("2", p), "1/1") for p in sorted(rng.sample(range(1, 1990), 20))]
    shuffled = records[:]
    rng.shuffle(shuffled)
    contigs = contig_lines({"1": 3000, "2": 2000})
    a = frames(resolve_ref(load(vcf(tmp_path, records, contigs, "a.vcf")), store))[1]
    b = frames(resolve_ref(load(vcf(tmp_path, shuffled, contigs, "b.vcf")), store))[1]
    assert a.drop("facet_id").equals(b.drop("facet_id"))
    again = frames(resolve_ref(load(vcf(tmp_path, records, contigs, "a.vcf")), store))[1]
    assert again.equals(a)
    assert set(a["ref_check"].cast(pl.String)) == {"match"}


def test_a_resolved_facet_is_not_resolved_again(tmp_path: Path, store: ReferenceStore) -> None:
    resolved = resolve_ref(load(vcf(tmp_path, [("1", 10, base("1", 10), "N", "0/1")], contig_lines({"1": 3000}))), store)
    with pytest.raises(ValueError, match="already read"):
        resolve_ref(resolved, store)


# Arrays: no stated lengths, so the facet's assembly label must agree with the reference's.

PREAMBLE_23 = (
    "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n"
    "#\n"
    "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
    "#\n"
    "# rsid\tchromosome\tposition\tgenotype\n"
)


def array_text(rows: list[tuple[str, str, int, str]], preamble: str = PREAMBLE_23) -> str:
    return preamble + "".join(f"{rsid}\t{chrom}\t{pos}\t{gt}\n" for rsid, chrom, pos, gt in rows)


@pytest.fixture
def grch37_store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "grch37")
    s.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    return s


def test_array_ref_comes_from_the_reference_and_alts_from_the_genotype(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    r = {p: base("1", p) for p in (10, 20, 30, 40, 50)}
    alt = {p: other(b) for p, b in r.items()}
    rows = [
        ("rs1", "1", 10, r[10] * 2),  # hom-ref
        ("rs2", "1", 20, alt[20] + r[20]),  # het, alt written first
        ("rs3", "1", 30, alt[30] * 2),  # hom-alt
        ("rs4", "1", 40, r[40] + "-"),  # half-call
        ("rs5", "1", 50, "--"),  # whole no-call
        ("i6", "1", 60, "DI"),  # 23andMe indel codes
        ("rs7", "X", 5, base("X", 5)),  # hemizygous X
        ("rs8", "25", 6, base("X", 6) * 2),  # PAR, positioned on X
        ("rs9", "MT", 7, base("MT", 7)),  # GRCh37's label cannot tell rCRS from hg19's chrM
        ("rs10", "2", 2001, "AA"),  # past the contig's end
    ]
    path = tmp_path / "genome.txt"
    path.write_text(array_text(rows))
    loaded = load(path)
    resolved = resolve_ref(loaded, grch37_store)
    _, calls = frames(resolved)
    by_id = {row["ids"][0]: row for row in calls.iter_rows(named=True)}

    assert by_id["rs1"]["ref"] == r[10] and by_id["rs1"]["alts"] == [] and by_id["rs1"]["gt"] == [0, 0]
    assert by_id["rs2"]["alts"] == [alt[20]] and by_id["rs2"]["gt"] == [1, 0]
    assert by_id["rs3"]["alts"] == [alt[30]] and by_id["rs3"]["gt"] == [1, 1]
    assert by_id["rs4"]["alts"] == [] and by_id["rs4"]["gt"] == [0, None]
    assert by_id["rs7"]["ref"] == base("X", 5) and by_id["rs7"]["gt"] == [0]
    assert by_id["rs8"]["chrom"] == "XY" and by_id["rs8"]["ref"] == base("X", 6)
    for unread in ("rs5", "i6", "rs9", "rs10"):
        assert by_id[unread]["ref"] is None and by_id[unread]["alts"] is None and by_id[unread]["gt"] is None
    assert calls["ref_check"].null_count() == calls.height  # no source REF to compare
    assert calls["source_ref"].null_count() == calls.height  # REF is never written into source_ref
    assert_gt_indexes_ref(calls)
    warnings = "\n".join(resolved.warnings)
    assert "1 array row(s) with an I/D indel code" in warnings
    assert "1 array row(s) with no allele called" in warnings
    assert "1 array row(s) past the end of their contig" in warnings
    assert "1 row(s) left with REF unread, on contig(s) the GRCh37 label cannot vouch for" in warnings
    assert "1 array row(s) (16.7% of those read) carry no REF allele" in warnings


def test_two_alts_are_sorted_and_indexed(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    ref = base("1", 100)
    a, b = sorted(x for x in "ACGT" if x != ref)[:2]
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 100, b + a)]))
    _, calls = frames(resolve_ref(load(path), grch37_store))
    row = calls.row(0, named=True)
    assert row["alts"] == [a, b] and row["gt"] == [2, 1]


def test_array_without_an_assembly_label_is_refused(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 10, "AA")], preamble="# rsid\tchromosome\tposition\tgenotype\n"))
    loaded = load(path)
    assert loaded.facets.collect()["assembly_label"].to_list() == [None]
    with pytest.raises(ReferenceMismatchError, match="neither contig lengths nor an assembly"):
        resolve_ref(loaded, grch37_store)


def test_array_on_another_assembly_label_is_refused(tmp_path: Path) -> None:
    s = ReferenceStore.open(tmp_path / "grch38")
    s.add_fasta(write_fasta(tmp_path / "grch38.fa", SEQ), assembly_label="GRCh38")
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 10, "AA")]))
    with pytest.raises(ReferenceMismatchError, match="assembly is GRCh37 and the reference is GRCh38"):
        resolve_ref(load(path), s)


def test_array_against_an_unlabelled_reference_is_refused(tmp_path: Path, store: ReferenceStore) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 10, "AA")]))
    with pytest.raises(ReferenceMismatchError, match="reference's assembly is unknown"):
        resolve_ref(load(path), store)


def test_array_resolution_keeps_canonical_order(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """Two probes at one position sort by their derived alts, as `canonical_order` would sort them."""
    ref = base("1", 200)
    x, y = sorted(b for b in "ACGT" if b != ref)[:2]
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 200, y + y), ("rs2", "1", 200, x + x)]))
    _, calls = frames(resolve_ref(load(path), grch37_store))
    assert [row[0] for row in calls["alts"].to_list()] == [x, y]

