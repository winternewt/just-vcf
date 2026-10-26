"""REF read from the reference sequence (RM3): `resolve_ref` against small synthetic FASTAs."""

import gzip
import hashlib
import random
from collections import Counter
from pathlib import Path

import polars as pl
import pytest
from gtars.refget import StorageMode, digest_sequence

from just_vcf import ReferenceMismatchError, ReferenceStore, load, mint_vrs_ids, normalize, reference_dir, resolve_ref
from just_vcf.loading import Loaded
from just_vcf.reference import Contig, _span
from just_vcf.schema import CALLS, FACETS, StrandCheck
from just_vcf.strand import MAX_FORWARD_MINUS, MIN_HETS, check_strand

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
    assert any(w.startswith("2 row(s) left with REF unread") and "MT (2 rows; length stated 402, reference 400)" in w for w in resolved.warnings)


def test_md5_is_compared_when_stated(tmp_path: Path, store: ReferenceStore) -> None:
    right = {"1": hashlib.md5(SEQ["1"].encode()).hexdigest()}
    records = [("1", 10, base("1", 10), "N", "0/1")]
    ok = vcf(tmp_path, records, contig_lines({"1": 3000}, md5=right), name="ok.vcf")
    assert frames(resolve_ref(load(ok), store))[1]["ref_check"].cast(pl.String).to_list() == ["match"]
    bad = vcf(tmp_path, records, contig_lines({"1": 3000}, md5={"1": "0" * 32}), name="bad.vcf")
    with pytest.raises(ReferenceMismatchError, match="md5_mismatch"):
        resolve_ref(load(bad), store)


WITH_Y = {**SEQ, "Y": "".join(random.Random(11).choice("ACGT") for _ in range(1200))}
"""The synthetic reference plus a Y, for the sex-chromosome rule (F20). Kept apart from `SEQ` so the
other tests' collection does not change."""


def md5_of(seqs: dict[str, str]) -> dict[str, str]:
    return {name: hashlib.md5(seq.encode()).hexdigest() for name, seq in seqs.items()}


@pytest.fixture
def store_with_y(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "store-y")
    s.add_fasta(write_fasta(tmp_path / "ref-y.fa", WITH_Y))
    return s


@pytest.mark.parametrize("sex", ["X", "Y"])
@pytest.mark.parametrize("kind", ["md5", "length"])
def test_a_sex_chromosome_that_differs_is_left_unread_when_the_autosomes_match(
    tmp_path: Path, store_with_y: ReferenceStore, sex: str, kind: str
) -> None:
    """F20: CGI's vcfBeta states a chrY of GRCh37's length with another md5. With every autosome
    matching, the autosomes are read and only the sex chromosome's rows stay unread, counted."""
    lengths = {name: len(seq) for name, seq in WITH_Y.items()}
    md5 = md5_of(WITH_Y)
    if kind == "md5":
        md5[sex] = "0" * 32
    else:
        lengths[sex] += 7
    records = [
        ("1", 10, base("1", 10), "N", "0/1"),
        ("2", 70, base("2", 70), "N", "1/1"),
        (sex, 5, WITH_Y[sex][4], "N", "1"),
        (sex, 6, WITH_Y[sex][5], "N", "1"),
        (sex, 9, WITH_Y[sex][8], "N", "1"),
    ]
    resolved = resolve_ref(load(vcf(tmp_path, records, contig_lines(lengths, md5=md5))), store_with_y)
    facets, calls = frames(resolved)
    assert facets["ref_seqcol_digest"].to_list() == store_with_y.collections()
    autosomes = calls.filter(pl.col("chrom").is_in(["1", "2"]))
    assert autosomes["ref"].to_list() == [base("1", 10), base("2", 70)]
    assert set(autosomes["ref_check"].cast(pl.String)) == {"match"}
    unread = calls.filter(pl.col("chrom") == sex)
    assert unread.height == 3 and unread["ref"].null_count() == 3
    assert set(unread["ref_check"].cast(pl.String)) == {"not_run"}
    reference = store_with_y.contigs(store_with_y.collections()[0])
    contig = next(c for c in reference if c.name == sex)
    stated = f"md5 stated {'0' * 32}, reference {contig.md5}" if kind == "md5" else f"length stated {lengths[sex]}, reference {contig.length}"
    line = next(w for w in resolved.warnings if "left with REF unread" in w)
    assert line.startswith("3 row(s) left with REF unread")
    assert f"stated {kind} differs" in line and f"{sex} (3 rows; {stated})" in line


def test_an_autosome_that_differs_still_refuses_beside_a_differing_y(tmp_path: Path, store_with_y: ReferenceStore) -> None:
    lengths = {name: len(seq) for name, seq in WITH_Y.items()}
    md5 = {**md5_of(WITH_Y), "2": "0" * 32, "Y": "0" * 32}
    records = [("1", 10, base("1", 10), "N", "0/1"), ("2", 70, base("2", 70), "N", "1/1"), ("Y", 5, WITH_Y["Y"][4], "N", "1")]
    with pytest.raises(ReferenceMismatchError, match=r"contig 2: md5_mismatch"):
        resolve_ref(load(vcf(tmp_path, records, contig_lines(lengths, md5=md5))), store_with_y)


def test_a_facet_whose_only_contig_is_a_differing_y_is_refused(tmp_path: Path, store_with_y: ReferenceStore) -> None:
    """Nothing to read: refused as a facet of which nothing matches, not read on an unchecked Y."""
    records = [("Y", 5, WITH_Y["Y"][4], "N", "1")]
    path = vcf(tmp_path, records, contig_lines({"Y": len(WITH_Y["Y"])}, md5={"Y": "0" * 32}))
    with pytest.raises(ReferenceMismatchError, match="no contig of the facet's calls could be matched"):
        resolve_ref(load(path), store_with_y)


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


MASKED = {"1": SEQ["1"][:100] + SEQ["1"][100:400].lower() + "n" * 50 + SEQ["1"][450:]}
"""Contig 1 soft-masked over 100-400 and with a lowercase N run over 400-450."""


@pytest.mark.parametrize("where", ["memory", "disk", "reopened"])
def test_fetch_returns_a_soft_masked_fasta_uppercase_as_gtars_stored_it(tmp_path: Path, where: str) -> None:
    """`fetch` keeps its "Uppercase" promise over soft-masked bases without uppercasing them itself:
    gtars' FASTA import stored them uppercase, and `fetch` returns them as read. A gtars that stopped
    uppercasing fails here, since nothing else would."""
    fasta = write_fasta(tmp_path / "masked.fa", MASKED)
    s = ReferenceStore.in_memory() if where == "memory" else ReferenceStore.open(tmp_path / "store")
    digest = s.add_fasta(fasta)
    if where == "reopened":
        s = ReferenceStore.open(tmp_path / "store")
    (contig,) = s.contigs(digest)
    ranges = [(90, 110), (380, 460), (0, contig.length)]
    assert contig.digest in s._uppercase_sequences()  # the store's record, so fetch takes gtars' text as is
    assert s.store.get_substrings(contig.digest, ranges) == [MASKED["1"][a:b].upper() for a, b in ranges]
    assert s.fetch(contig, ranges) == [MASKED["1"][a:b].upper() for a, b in ranges]


def test_fetch_uppercases_a_sequence_that_did_not_come_from_a_fasta(tmp_path: Path) -> None:
    """gtars stores a `SequenceRecord` whose `sequence` was set by hand as given in a raw store, so
    lowercase does reach its `get_substrings`; the store records no case for it, and `fetch` uppercases.
    A FASTA's collection in the same store still takes the unscanned path."""
    s = ReferenceStore.in_memory()
    s.store.set_encoding_mode(StorageMode.Raw)
    (fasta_contig,) = s.contigs(s.add_fasta(write_fasta(tmp_path / "masked.fa", MASKED)))
    sequence = "ACGTacgtnnNNACGT"
    record = digest_sequence(sequence.upper().encode(), name="hand")
    record.sequence = sequence.encode()
    s.store.add_sequence(record, force=True)
    hand = Contig(name="hand", length=len(sequence), md5=record.metadata.md5, digest=record.metadata.sha512t24u)
    ranges = [(0, 16), (4, 8), (8, 10)]
    assert s.store.get_substrings(hand.digest, ranges) == [sequence[a:b] for a, b in ranges]  # gtars kept the lowercase
    assert hand.digest not in s._uppercase_sequences()
    assert s.fetch(hand, ranges) == [sequence[a:b].upper() for a, b in ranges]
    assert fasta_contig.digest in s._uppercase_sequences()
    assert s.fetch(fasta_contig, [(90, 110)]) == [MASKED["1"][90:110].upper()]


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


def test_a_fasta_path_is_the_collection_used_even_beside_others(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_REFERENCE_DIR", str(tmp_path / "configured"))
    other_digest = ReferenceStore.open().add_fasta(write_fasta(tmp_path / "plain.fa", SEQ))
    prefixed = write_fasta(tmp_path / "prefixed.fa", SEQ, names=UCSC)
    path = vcf(tmp_path, [("1", 10, base("1", 10), "N", "0/1")], contig_lines({"1": 3000}))
    facets, calls = frames(resolve_ref(load(path), prefixed))
    assert facets["ref_seqcol_digest"].to_list() != [other_digest]
    assert facets["ref_seqcol_digest"].to_list() == [ReferenceStore.open().add_fasta(prefixed)]
    assert calls["ref"].to_list() == [base("1", 10)]


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


def test_a_vcf_without_contig_lengths_is_matched_by_its_stated_build(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    records = [("1", 10, base("1", 10), "N", "0/1"), ("MT", 5, base("MT", 5), "N", "1")]
    contigs = ["##reference=file:///refs/human_g1k_v37.fasta", "##contig=<ID=1>", "##contig=<ID=MT>"]
    loaded = load(vcf(tmp_path, records, contigs))
    assert loaded.facets.collect()["assembly_label"].to_list() == ["GRCh37"]
    resolved = resolve_ref(loaded, grch37_store)
    _, calls = frames(resolved)
    assert calls["ref"].to_list() == [base("1", 10), None]
    assert calls["ref_check"].cast(pl.String).to_list() == ["match", "not_run"]
    assert any("GRCh37 label cannot vouch for without a stated length: MT" in w for w in resolved.warnings)


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


def test_array_without_an_assembly_label_and_too_few_rows_is_refused(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """One row cannot show which build it is on (F9's row inference needs `MIN_SCORED_ROWS`)."""
    path = tmp_path / "genome.txt"
    path.write_text(array_text([("rs1", "1", 10, "AA")], preamble="# rsid\tchromosome\tposition\tgenotype\n"))
    loaded = load(path)
    assert loaded.facets.collect()["assembly_label"].to_list() == [None]
    with pytest.raises(ReferenceMismatchError, match=r"neither contig lengths nor an assembly, and its rows do not settle one \(1 SNV row"):
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


# F9: an array that states no build has it inferred from its own rows against each collection.

BARE = "# rsid\tchromosome\tposition\tgenotype\n"
"""A 23andMe-shaped preamble that names no build, as an FTDNA export names none."""

OTHER_SEQ = sequences(seed=11)
"""A second synthetic assembly: same contig names and lengths, other bases."""


def drawn_rows(seq: dict[str, str], seed: int = 3, per_contig: int = 700) -> list[tuple[str, str, int, str]]:
    """Array rows whose genotypes are drawn around `seq`'s bases: hom-ref, het and hom-alt calls."""
    rng = random.Random(seed)
    rows = []
    for chrom in ("1", "2"):
        for pos in range(1, per_contig + 1):
            ref = seq[chrom][pos - 1]
            alt = rng.choice([b for b in "ACGT" if b != ref])
            kind = rng.random()
            gt = ref * 2 if kind < 0.5 else ref + alt if kind < 0.8 else alt * 2
            rows.append((f"rs{chrom}{pos}", chrom, pos, gt))
    return rows


def share_carrying(rows: list[tuple[str, str, int, str]], seq: dict[str, str]) -> float:
    """The test's own count: rows whose genotype carries `seq`'s base at their position."""
    return sum(seq[chrom][pos - 1] in gt for _, chrom, pos, gt in rows) / len(rows)


def labelled_store(tmp_path: Path, name: str, assemblies: dict[str, dict[str, str]]) -> tuple[ReferenceStore, dict[str, str]]:
    """A store with one collection per label (a label starting `unlabelled` asserts none), and label → digest."""
    s = ReferenceStore.open(tmp_path / name)
    digests = {}
    for label, seq in assemblies.items():
        fasta = write_fasta(tmp_path / f"{name}-{label}.fa", seq)
        digests[label] = s.add_fasta(fasta, assembly_label=None if label.startswith("unlabelled") else label)
    return s, digests


def bare_array(tmp_path: Path, rows: list[tuple[str, str, int, str]], preamble: str = BARE) -> Loaded:
    path = tmp_path / "genome.txt"
    path.write_text(array_text(rows, preamble=preamble))
    return load(path)


def test_an_array_without_a_build_gets_it_from_its_rows(tmp_path: Path) -> None:
    rows = drawn_rows(SEQ)
    s, digests = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    loaded = bare_array(tmp_path, rows)
    assert loaded.facets.collect()["assembly_label"].to_list() == [None]

    resolved = resolve_ref(loaded, s)
    facets, calls = frames(resolved)
    assert facets["assembly_label"].to_list() == ["GRCh37"]
    assert facets["assembly_label_basis"].cast(pl.String).to_list() == ["inferred"]
    assert facets["ref_seqcol_digest"].to_list() == [digests["GRCh37"]]
    expected = {(chrom, pos): SEQ[chrom][pos - 1] for _, chrom, pos, _ in rows}
    assert {(c, p): r for c, p, r in calls.select("chrom", "pos", "ref").iter_rows()} == expected
    assert_gt_indexes_ref(calls)
    line = next(w for w in resolved.warnings if w.startswith("assembly inferred from the rows as GRCh37"))
    right, wrong = share_carrying(rows, SEQ), share_carrying(rows, OTHER_SEQ)
    assert f"{right:.1%} on GRCh37 ({digests['GRCh37']}, {len(rows)} rows)" in line
    assert f"{wrong:.1%} on GRCh38 ({digests['GRCh38']}, {len(rows)} rows)" in line
    assert loaded.facets.collect()["assembly_label"].to_list() == [None]  # the input is not changed


def test_the_inferred_build_follows_the_rows_not_the_store_order(tmp_path: Path) -> None:
    s, digests = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    facets, _ = frames(resolve_ref(bare_array(tmp_path, drawn_rows(OTHER_SEQ)), s))
    assert (facets["assembly_label"].item(), facets["ref_seqcol_digest"].item()) == ("GRCh38", digests["GRCh38"])


def test_one_collection_is_inferred_against_chance(tmp_path: Path) -> None:
    s, _ = labelled_store(tmp_path, "one", {"GRCh37": SEQ})
    resolved = resolve_ref(bare_array(tmp_path, drawn_rows(SEQ)), s)
    assert resolved.facets.collect()["assembly_label"].to_list() == ["GRCh37"]
    assert any("; chance " in w for w in resolved.warnings if w.startswith("assembly inferred"))


def test_rows_that_fit_no_collection_are_refused_with_the_scores(tmp_path: Path) -> None:
    s, _ = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    with pytest.raises(ReferenceMismatchError, match=r"rows do not settle one \(no collection reaches 60%.*% on GRCh3[78] "):
        resolve_ref(bare_array(tmp_path, drawn_rows(sequences(seed=23))), s)


def test_two_labels_on_the_same_bases_leave_no_margin(tmp_path: Path) -> None:
    """Two collections that differ only on MT (never scored) cannot be told apart by the rows."""
    s, _ = labelled_store(tmp_path, "same", {"GRCh37": SEQ, "GRCh38": dict(SEQ, MT=SEQ["MT"][::-1])})
    with pytest.raises(ReferenceMismatchError, match="within 20%"):
        resolve_ref(bare_array(tmp_path, drawn_rows(SEQ)), s)


def test_an_unlabelled_best_collection_names_nothing(tmp_path: Path) -> None:
    s, _ = labelled_store(tmp_path, "unl", {"unlabelled": SEQ, "GRCh38": OTHER_SEQ})
    with pytest.raises(ReferenceMismatchError, match="has no assembly label"):
        resolve_ref(bare_array(tmp_path, drawn_rows(SEQ)), s)


def test_too_few_rows_infer_nothing(tmp_path: Path) -> None:
    s, _ = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    with pytest.raises(ReferenceMismatchError, match=r"800 SNV row\(s\) could be scored on a collection; 1000 are needed"):
        resolve_ref(bare_array(tmp_path, drawn_rows(SEQ, per_contig=400)), s)


GRCH36_CHR1 = 247_249_719
"""NCBI36/hg18 chr1's length (UCSC hg18.chrom.sizes), a domain constant."""


CHR1_LENGTH = {"GRCh37": 249_250_621, "GRCh38": 248_956_422}
"""GRCh37's and GRCh38's chr1 lengths (NCBI's assembly reports), domain constants."""


def contig_1_fasta(fasta: Path, ones: list[tuple[str, str]]) -> Path:
    """A gzipped FASTA of the given contig 1s, then contig 2. A pair (name, assembly) is all N at that
    assembly's chr1 length, so the file stays small; (name, "short") is `SEQ['1']`, a length that
    names no assembly."""
    with gzip.open(fasta, "wb", compresslevel=1) as handle:
        for name, assembly in ones:
            if assembly == "short":
                handle.write(f">{name}\n{SEQ['1']}\n".encode())
                continue
            full, rest = divmod(CHR1_LENGTH[assembly], 60)
            handle.write(f">{name}\n".encode() + (b"N" * 60 + b"\n") * full + b"N" * rest + b"\n")
        handle.write(f">2\n{SEQ['2']}\n".encode())
    return fasta


def test_a_collection_serving_contig_1_twice_names_no_assembly(tmp_path: Path) -> None:
    """`chr1` at GRCh37's length, then `1` at GRCh38's: which length is contig 1's is not known, so
    the lengths name no assembly, as `contigs_by_chrom` lists both and a caller reads neither. Before,
    the later sequence overwrote the earlier one and named GRCh38.

    An asserted label is refused all the same, whichever of the two it names: each contig 1's length
    names an assembly, and one of them is not the asserted one. The lengths name none only because
    they disagree, which is no ground to trust the assertion over them. A refused label leaves no
    alias behind, so the collection stays unlabelled."""
    fasta = contig_1_fasta(tmp_path / "two_ones.fa.gz", [("chr1", "GRCh37"), ("1", "GRCh38")])
    s = ReferenceStore.in_memory()
    with pytest.raises(ReferenceMismatchError, match=r"two_ones\.fa\.gz: its contig 1 is GRCh38's length, not GRCh37's"):
        s.add_fasta(fasta, assembly_label="GRCh37")
    with pytest.raises(ReferenceMismatchError, match=r"two_ones\.fa\.gz: its contig 1 is GRCh37's length, not GRCh38's"):
        s.add_fasta(fasta, assembly_label="GRCh38")
    digest = s.add_fasta(fasta)
    assert [c.name for c in s.contigs_by_chrom(digest)["1"]] == ["chr1", "1"]
    assert s.asserted_labels(digest) == []
    assert s.assembly_label(digest) is None


def test_a_label_asserted_over_two_contig_1s_holds_where_no_length_names_another(tmp_path: Path) -> None:
    """`chr1` at GRCh37's length and `1` at a length that names nothing: no length contradicts
    GRCh37, so the asserted label is kept and names the collection, which its lengths alone do not."""
    fasta = contig_1_fasta(tmp_path / "one_named.fa.gz", [("chr1", "GRCh37"), ("1", "short")])
    s = ReferenceStore.in_memory()
    digest = s.add_fasta(fasta, assembly_label="GRCh37")
    assert s._label_from_lengths(digest) is None
    assert s.asserted_labels(digest) == ["GRCh37"]
    assert s.assembly_label(digest) == "GRCh37"


def test_a_grch36_reference_is_named_by_its_contig_1_and_cannot_vouch_for_mt(tmp_path: Path) -> None:
    """Contig 1 of GRCh36's length (all N, so the file stays small) names the collection GRCh36."""
    fasta = tmp_path / "hg18_like.fa.gz"
    full, rest = divmod(GRCH36_CHR1, 60)
    with gzip.open(fasta, "wb", compresslevel=1) as handle:
        handle.write(b">chr1\n" + (b"N" * 60 + b"\n") * full + b"N" * rest + b"\n")
        handle.write(f">chr2\n{SEQ['2']}\n>chrM\n{SEQ['MT']}\n".encode())
    s = ReferenceStore.in_memory()
    digest = s.add_fasta(fasta)
    assert s.assembly_label(digest) == "GRCh36"
    preamble = PREAMBLE_23.replace("build 37 (also known as Annotation Release 104)", "build 36")
    loaded = bare_array(tmp_path, [("rs1", "2", 10, base("2", 10) * 2), ("rs2", "MT", 7, base("MT", 7))], preamble=preamble)
    assert loaded.facets.collect()["assembly_label"].to_list() == ["GRCh36"]
    resolved = resolve_ref(loaded, s)
    _, calls = frames(resolved)
    assert dict(calls.select("chrom", "ref").iter_rows()) == {"2": base("2", 10), "MT": None}
    assert any("the GRCh36 label cannot vouch for without a stated length: MT" in w for w in resolved.warnings)


def test_the_rows_name_t2t_and_its_label_cannot_vouch_for_mt(tmp_path: Path) -> None:
    """RM286: F9 names T2T-CHM13 like any label, and `T2T-CHM13` covers two mitochondria of one length
    (the CHM13 cell line's own in `chm13v2.0.fa`, rCRS in the analysis set's `_rCRS` variant)."""
    s, digests = labelled_store(tmp_path, "two", {"GRCh38": SEQ, "T2T-CHM13": OTHER_SEQ})
    rows = [*drawn_rows(OTHER_SEQ), ("rsM", "MT", 7, OTHER_SEQ["MT"][6])]
    resolved = resolve_ref(bare_array(tmp_path, rows), s)
    facets, calls = frames(resolved)
    assert (facets["assembly_label"].item(), facets["ref_seqcol_digest"].item()) == ("T2T-CHM13", digests["T2T-CHM13"])
    assert calls.filter(pl.col("chrom") == "MT")["ref"].to_list() == [None]
    assert calls.filter(pl.col("chrom") != "MT")["ref"].null_count() == 0
    assert any("the T2T-CHM13 label cannot vouch for without a stated length: MT" in w for w in resolved.warnings)


def test_a_vcf_without_lengths_or_a_build_is_not_scored_by_its_rows(tmp_path: Path) -> None:
    """Scoring a VCF by its source REF would choose the build its `ref_check` then confirms."""
    s, _ = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    records = [("1", p, base("1", p), other(base("1", p)), "0/1") for p in range(1, 1201)]
    loaded = load(vcf(tmp_path, records, contigs=[]))
    with pytest.raises(ReferenceMismatchError, match="so no reference can be shown to be its own") as refused:
        resolve_ref(loaded, s)
    assert "rows do not settle" not in str(refused.value)


def test_a_stated_build_is_not_inferred_again(tmp_path: Path) -> None:
    """A facet that states its build keeps it and its basis, even where the rows point elsewhere."""
    s, digests = labelled_store(tmp_path, "two", {"GRCh37": SEQ, "GRCh38": OTHER_SEQ})
    resolved = resolve_ref(bare_array(tmp_path, drawn_rows(OTHER_SEQ), preamble=PREAMBLE_23), s)
    facets = resolved.facets.collect()
    assert facets["assembly_label"].to_list() == ["GRCh37"]
    assert facets["assembly_label_basis"].cast(pl.String).to_list() == ["stated"]
    assert facets["ref_seqcol_digest"].to_list() == [digests["GRCh37"]]
    assert not any(w.startswith("assembly inferred") for w in resolved.warnings)


# F13: an array export with part of its probes written on the minus strand.

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def flipped(rows: list[tuple[str, str, int, str]], every: int = 8) -> tuple[list[tuple[str, str, int, str]], set[str]]:
    """`rows` with every `every`-th genotype complemented, as a minus-strand probe writes it, and those ids."""
    out, ids = [], set()
    for n, (rsid, chrom, pos, gt) in enumerate(rows):
        if n % every == 0:
            gt = "".join(COMPLEMENT[b] for b in gt)
            ids.add(rsid)
        out.append((rsid, chrom, pos, gt))
    return out, ids


def test_a_mixed_strand_array_keeps_ref_and_withholds_the_alleles_it_cannot_read(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    true_rows = drawn_rows(SEQ)
    rows, minus = flipped(true_rows)
    resolved = resolve_ref(bare_array(tmp_path, rows, preamble=PREAMBLE_23), grch37_store)
    facets, calls = frames(resolved)
    assert facets["strand_check"].cast(pl.String).to_list() == ["mixed"]

    truth = {rsid: gt for rsid, _, _, gt in true_rows}
    expected: Counter[str] = Counter()
    for row in calls.iter_rows(named=True):
        rsid, ref = row["ids"][0], row["ref"]
        assert ref == SEQ[row["chrom"]][row["pos"] - 1], "REF is read either way"
        written = set(row["gt_alleles"])
        if len(row["gt_alleles"]) == 2 and len(written) == 2 and ref in written:
            # Strand-safe: a forward het, or an A/T or C/G het that reads the same flipped.
            assert row["alts"] == sorted(set(truth[rsid]) - {ref}), rsid
            continue
        assert row["alts"] is None and row["gt"] is None, rsid
        kind = "minus" if len(written) == 2 and COMPLEMENT[ref] in written else "neither" if len(written) == 2 else "undecidable"
        expected[kind] += 1
    assert minus and expected["minus"] and expected["undecidable"]
    warnings = "\n".join(resolved.warnings)
    assert f"{expected['minus']} array row(s) with alts and gt withheld, the facet's strand being mixed: heterozygous rows carrying" in warnings
    assert f"{expected['undecidable']} array row(s) with alts and gt withheld, the facet's strand being mixed: homozygous rows" in warnings
    assert "strand checked on " in warnings and ": mixed" in warnings

    # Minting still runs: a withheld row has its REF slot only, and asserts no ALT.
    minted = mint_vrs_ids(resolved, grch37_store).calls.collect()
    withheld = minted.filter(pl.col("alts").is_null())
    assert withheld.height == sum(expected.values())
    assert withheld["vrs_allele_ids"].list.len().to_list() == [1] * withheld.height
    assert withheld["vrs_allele_ids"].list.first().null_count() == 0


def test_a_forward_array_derives_every_allele_and_says_so(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    resolved = resolve_ref(bare_array(tmp_path, drawn_rows(SEQ), preamble=PREAMBLE_23), grch37_store)
    facets, calls = frames(resolved)
    assert facets["strand_check"].cast(pl.String).to_list() == ["forward"]
    assert calls["alts"].null_count() == 0 and calls["gt"].null_count() == 0
    assert_gt_indexes_ref(calls)
    line = next(w for w in resolved.warnings if w.startswith("strand checked on "))
    assert " 0 (0.00%) carry the complement of the reference base and not the base" in line and line.endswith(": forward")


def test_too_few_hets_leave_the_strand_unchecked_not_forward(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """An X-only export of hemizygous calls has no hets to read the strand from: null, derived as before."""
    rows = [(f"rs{p}", "X", p, base("X", p)) for p in range(1, 1201)]
    resolved = resolve_ref(bare_array(tmp_path, rows, preamble=PREAMBLE_23), grch37_store)
    facets, calls = frames(resolved)
    assert facets["strand_check"].to_list() == [None]
    assert calls["alts"].null_count() == 0
    assert any(w.startswith("strand not checked: 0 sampled heterozygous row(s), 200 needed") for w in resolved.warnings)


def test_a_vcf_gets_no_strand_check(tmp_path: Path, store: ReferenceStore) -> None:
    records = [("1", 10, base("1", 10), other(base("1", 10)), "0/1")]
    resolved = resolve_ref(load(vcf(tmp_path, records, contig_lines(LENGTHS))), store)
    assert resolved.facets.collect()["strand_check"].to_list() == [None]
    assert not any(w.startswith("strand ") for w in resolved.warnings)


def test_the_strand_threshold_is_a_share_of_hets() -> None:
    """At `MAX_FORWARD_MINUS` a facet still reads forward; one het more reads mixed."""
    hets = 1_000
    at = int(hets * MAX_FORWARD_MINUS)

    def scored(minus: int) -> StrandCheck | None:
        sample = pl.DataFrame({"alleles": [["C", "T"]] * minus + [["A", "G"]] * (hets - minus)})
        return check_strand(sample, pl.Series(["A"] * hets)).check

    assert scored(at) is StrandCheck.FORWARD
    assert scored(at + 1) is StrandCheck.MIXED
    few = MIN_HETS - 1
    small = check_strand(pl.DataFrame({"alleles": [["C", "T"]] * few}), pl.Series(["A"] * few))
    assert small.check is None and small.minus == few


# F4: Living DNA writes an indel genotype as its two allele sequences back to back (`CCA`).

LIVING_DNA = (
    "# Living DNA customer genotype data download file version: 1.0.2\n"
    "# Genotype chip: Sirius\n"
    "# Human Genome Reference Build 37 (GRCh37.p13).\n"
    "# Genotypes are presented on the forward strand.\n"
    "#\n"
    "# rsid\tchromosome\tposition\tgenotype\n"
)


def concatenated_cases() -> dict[str, tuple[int, str]]:
    """One run-together genotype per outcome, planted on SEQ["1"] after the drawn rows: name -> (pos, text).
    Each text is built from the sequence at its position, so the expected split is the construction's."""
    seq = SEQ["1"]

    def at(p: int, n: int = 1) -> str:
        return seq[p - 1 : p - 1 + n]

    positions = iter(range(1000, 2900, 7))
    cases = {}
    p = next(positions)
    b = at(p)
    x = next(c for c in "ACGT" if c not in (b, at(p + 1)))
    cases["insertion"] = (p, b + b + x)  # b / bx: b is the reference, bx not (x is not the next base)
    p = next(q for q in positions if at(q) not in at(q + 1, 4))
    cases["deletion"] = (p, at(p, 5) + at(p))  # five reference bases / their first: both equal the reference over their spans
    p = next(q for q in positions if at(q) != at(q + 1))
    cases["homozygous"] = (p, at(p, 2) * 2)
    p = next(q for q in positions if at(q) == at(q + 1))
    cases["several"] = (p, at(p) * 4)  # c/ccc and cc/cc both fit
    p = next(positions)
    y = other(at(p))
    cases["no_fit"] = (p, y + y + next(c for c in "ACGT" if c != y))
    return cases


def test_living_dna_indel_genotypes_split_against_the_reference(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    cases = concatenated_cases()
    rows = drawn_rows(SEQ, per_contig=300) + [(name, "1", pos, text) for name, (pos, text) in cases.items()]
    loaded = bare_array(tmp_path, rows, preamble=LIVING_DNA)
    assert loaded.facets.collect()["format_id"].to_list() == ["living_dna"]
    before = {r["ids"][0]: r for r in loaded.calls.collect().iter_rows(named=True) if r["ids"][0] in cases}
    assert all(r["gt_alleles"] is None for r in before.values())  # the loader cannot split them

    resolved = resolve_ref(loaded, grch37_store)
    _, calls = frames(resolved)
    got = {r["ids"][0]: r for r in calls.iter_rows(named=True) if r["ids"][0] in cases}
    seq = SEQ["1"]

    pos, text = cases["insertion"]
    ins = got["insertion"]
    assert (ins["gt_alleles"], ins["ref"], ins["alts"], ins["gt"]) == ([text[0], text[1:]], seq[pos - 1], [text[1:]], [0, 1])
    assert ins["end_pos"] == pos and ins["ref_check"] is None and ins["source_gt"] == text
    _, text = cases["deletion"]
    assert got["deletion"]["gt_alleles"] == [text[:5], text[0]]
    _, text = cases["homozygous"]
    assert got["homozygous"]["gt_alleles"] == [text[:2], text[:2]]
    for name in ("deletion", "homozygous"):
        assert (got[name]["ref"], got[name]["alts"], got[name]["gt"]) == (None, None, None), name
    for name in ("several", "no_fit"):
        assert (got[name]["gt_alleles"], got[name]["ref"]) == (None, None), name
    # The drawn SNV rows resolve as before, and every row with a REF indexes it.
    assert_gt_indexes_ref(calls)
    warnings = "\n".join(resolved.warnings)
    assert (
        "3 genotype(s) written as two allele sequences back to back split where exactly one split fits the reference "
        "(both alleles start with the same base, one equals the reference over its own span): 1 with REF and ALT read; "
        "REF left unread on 1 homozygous call(s)"
    ) in warnings
    assert "on 1 where both alleles equal the reference over their own spans" in warnings
    assert "1 genotype(s) written as two allele sequences back to back fit the reference at more than one split" in warnings
    assert "1 genotype(s) written as two allele sequences back to back fit the reference at no split" in warnings
    assert "2 array row(s) whose genotype could not be split into alleles" in warnings
    # The row with REF and ALT is an indel like any other: VRS ids are minted for both alleles.
    minted = mint_vrs_ids(resolved, grch37_store).calls.filter(pl.col("ids").list.first() == "insertion").collect().row(0, named=True)
    assert len(minted["vrs_allele_ids"]) == 2 and None not in minted["vrs_allele_ids"]


def test_normalize_trims_the_living_dna_rows_resolve_ref_gave_a_ref(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """A row of a facet read from a file, whose REF `resolve_ref` read with no source REF to check (an
    array row, a split Living DNA genotype), is aligned and trimmed like a matched VCF row.
    `bxy`/`bzy` (an MNV written with its padding) trims to the SNV `x>z` one base on; `bxy`/`by` trims
    to the deletion `bx>b` in place."""
    seq = SEQ["1"]
    p = next(q for q in range(1000, 2900) if seq[q - 1] not in seq[q : q + 2] and seq[q] != seq[q + 1])
    b, x, y = seq[p - 1 : p + 2]
    z = next(c for c in "ACGT" if c not in (b, x, y))
    q = next(r for r in range(p + 10, 2900) if seq[r - 1] not in seq[r : r + 2] and seq[r] != seq[r + 1])
    d, e, f = seq[q - 1 : q + 2]
    rows = drawn_rows(SEQ, per_contig=300) + [("mnv", "1", p, b + x + y + b + z + y), ("deletion", "1", q, d + e + f + d + f)]
    resolved = resolve_ref(bare_array(tmp_path, rows, preamble=LIVING_DNA), grch37_store)
    before = {r["ids"][0]: r for r in resolved.calls.collect().iter_rows(named=True)}
    assert (before["mnv"]["ref"], before["mnv"]["alts"], before["mnv"]["gt"]) == (b + x + y, [b + z + y], [0, 1])
    assert (before["deletion"]["ref"], before["deletion"]["alts"]) == (d + e + f, [d + f])

    normalized = normalize(resolved, grch37_store)
    after = {r["ids"][0]: r for r in normalized.calls.collect().iter_rows(named=True)}
    mnv, deletion = after["mnv"], after["deletion"]
    assert (mnv["pos"], mnv["end_pos"], mnv["ref"], mnv["alts"], mnv["gt"], mnv["gt_alleles"]) == (p + 1, p + 1, x, [z], [0, 1], [x, z])
    assert (deletion["pos"], deletion["ref"], deletion["alts"], deletion["gt_alleles"]) == (q, d + e, [d], [d + e, d])
    assert mnv["source_gt"] is None, "the trimmed row is no longer what the file wrote"
    # The array SNV rows were already as short as they go: unchanged.
    columns = ["pos", "ref", "alts", "gt", "gt_alleles"]
    snvs = resolved.calls.filter(~pl.col("ids").list.first().is_in(["mnv", "deletion"])).select(columns).collect()
    assert normalized.calls.filter(~pl.col("ids").list.first().is_in(["mnv", "deletion"])).select(columns).collect().equals(snvs)
    assert "1 row(s) left-aligned to another position and 1 trimmed in place, against the reference" in normalized.warnings


def test_only_living_dna_texts_are_split(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """Another format's long genotype text states no encoding: it stays as loaded."""
    cases = concatenated_cases()
    rows = drawn_rows(SEQ, per_contig=300) + [(name, "1", pos, text) for name, (pos, text) in cases.items()]
    resolved = resolve_ref(bare_array(tmp_path, rows, preamble=PREAMBLE_23), grch37_store)
    _, calls = frames(resolved)
    assert calls.filter(pl.col("ids").list.first().is_in(list(cases)))["gt_alleles"].null_count() == len(cases)
    assert not any("back to back" in w for w in resolved.warnings)


def test_a_mixed_strand_living_dna_file_splits_nothing(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    cases = concatenated_cases()
    rows, _ = flipped(drawn_rows(SEQ))
    rows += [(name, "1", pos, text) for name, (pos, text) in cases.items()]
    resolved = resolve_ref(bare_array(tmp_path, rows, preamble=LIVING_DNA), grch37_store)
    facets, calls = frames(resolved)
    assert facets["strand_check"].cast(pl.String).to_list() == ["mixed"]
    assert calls.filter(pl.col("ids").list.first().is_in(list(cases)))["gt_alleles"].null_count() == len(cases)
    warnings = "\n".join(resolved.warnings)
    assert f"{len(cases)} genotype(s) written as two allele sequences back to back were not split, the facet's strand being mixed" in warnings

