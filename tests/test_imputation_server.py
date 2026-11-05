"""RM287: a chip export taken the way a user takes it to an imputation server, `load` → `resolve_ref`
→ `liftover` → `write_vcf`, and the written file checked against the servers' published input rules.

The rules are `imputation_rules.check_vcf`'s, each with its source in docs/SCHEMA.md, "What an
imputation server checks". A rule just-vcf's output does not meet is an `xfail(strict=True)` with its
reason: what to do about each is the lead's decision, and a fix turns its test red until the mark goes.

A synthetic GRCh37 assembly (contigs 1, 2 and an X long enough to reach the nonPAR region), a GRCh38
built from it (contig 1 and X shifted, contig 2 reverse-complemented, so a chip row lands on a minus
chain), and the chain between them, all under tmp_path. Expected counts come from the export's own
rows, read here against the test's sequences, never from the library's frames.
"""

import random
from collections import Counter
from pathlib import Path

import polars as pl
import pytest

from imputation_rules import NON_PAR, PANEL_RULES, Outcome, RefBases, check_vcf, palindromes, read_vcf_text, strand_withheld
from just_vcf import Loaded, ReferenceStore, UnwritableRowsError, liftover, load, resolve_ref, sniff, write_vcf
from just_vcf.writing import ContigNames
from test_reference import COMPLEMENT, PREAMBLE_23, array_text, drawn_rows, flipped, write_fasta

rng = random.Random(41)


def bases(n: int) -> str:
    return "".join(rng.choices("ACGT", k=n))


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT", "TGCA"))


NONPAR_37 = NON_PAR["GRCh37"][0]
SOURCE = {"1": bases(3000), "2": bases(2000), "X": bases(NONPAR_37 + 700)}
"""GRCh37 by assertion. X reaches 700 bases into the nonPAR region the servers split chrX by."""
X_SHIFT = 90_000
FILLER = {"1": bases(100), "2": bases(50), "X": bases(X_SHIFT)}
TARGET = {"1": FILLER["1"] + SOURCE["1"], "2": FILLER["2"] + revcomp(SOURCE["2"]), "X": FILLER["X"] + SOURCE["X"]}
"""GRCh38 by assertion: X's shift carries the nonPAR rows past GRCh38's nonPAR start too."""
assert NONPAR_37 + X_SHIFT > NON_PAR["GRCh38"][0]


def chain_text() -> str:
    t1, t2, tx, sx = len(TARGET["1"]), len(TARGET["2"]), len(TARGET["X"]), len(SOURCE["X"])
    return "\n".join([
        f"chain 1000 1 3000 + 0 3000 1 {t1} + 100 {t1} 1", "3000", "",
        f"chain 900 2 2000 + 0 2000 2 {t2} - 0 2000 2", "2000", "",
        f"chain 800 X {sx} + 0 {sx} X {tx} + {X_SHIFT} {tx} 3", str(sx), "",
    ])  # fmt: skip


def lifted_pos(chrom: str, pos: int) -> int:
    """Where the construction puts a 1-based source position on the target."""
    return {"1": pos + 100, "2": 50 + 2000 - pos + 1, "X": pos + X_SHIFT}[chrom]


def reference_bases(seqs: dict[str, str]) -> RefBases:
    """The test's own reading of `seqs` under a written contig name (`chr` style, `XY` on X)."""

    def read(chrom: str, spans: list[tuple[int, int]]) -> list[str | None]:
        name = chrom.removeprefix("chr")
        seq = seqs.get({"XY": "X", "M": "MT"}.get(name, name))
        if seq is None:
            return [None] * len(spans)
        return [seq[pos - 1 : pos - 1 + length] if pos + length - 1 <= len(seq) else None for pos, length in spans]

    return read


@pytest.fixture(scope="module")
def store(tmp_path_factory: pytest.TempPathFactory) -> tuple[ReferenceStore, str, str, Path]:
    root = tmp_path_factory.mktemp("imputation")
    found = ReferenceStore.open(root / "store")
    source = found.add_fasta(write_fasta(root / "grch37.fa", SOURCE), assembly_label="GRCh37")
    target = found.add_fasta(write_fasta(root / "grch38.fa", TARGET), assembly_label="GRCh38")
    chain = root / "grch37_to_grch38.chain"
    chain.write_text(chain_text())
    return found, source, target, chain


Row = tuple[str, str, int, str]


def base(chrom: str, pos: int) -> str:
    return SOURCE["X" if chrom == "25" else chrom][pos - 1]


def other(*taken: str) -> str:
    return next(b for b in "ACGT" if b not in taken)


DRAWN = drawn_rows(SOURCE)
"""700 rows on each of contigs 1 and 2: hom-ref, het and hom-alt calls around the reference's base."""


def extras() -> list[Row]:
    """What a real export also holds: a second probe at a position, a call of two non-reference
    bases, hemizygous X calls in the nonPAR region, the pseudo-autosomal code 25, a no-call and an
    indel code. rsids must read as rsids (`rs` or `i` and digits) for the sniffer."""
    second = next(r for r in DRAWN if (r[1], r[2]) == ("2", 700))
    ref = base("1", 800)
    x_rows = [(f"rs91{n:02d}", "X", NONPAR_37 + 80 + n, base("X", NONPAR_37 + 80 + n) if n % 3 else other(base("X", NONPAR_37 + 80 + n))) for n in range(20)]
    return [
        ("rs9000001", "2", 700, second[3]),
        ("rs9000002", "1", 800, other(ref) + other(ref, other(ref))),
        *x_rows,
        ("rs9000101", "25", 1000, base("X", 1000) * 2),
        ("rs9000102", "25", 1001, base("X", 1001) + other(base("X", 1001))),
        ("rs9000201", "1", 900, "--"),
        ("i9000301", "1", 950, "DI"),
    ]


EXTRAS = extras()
FORWARD = [*DRAWN, *EXTRAS]
MIXED = [*flipped(DRAWN)[0], *EXTRAS]
"""Every eighth drawn genotype complemented, as a minus-strand probe writes it (F13)."""


def read_chip(tmp_path: Path, store: tuple[ReferenceStore, str, str, Path], rows: list[Row], name: str = "genome.txt") -> Loaded:
    found, source, _, _ = store
    path = tmp_path / name
    path.write_text(array_text(rows, preamble=PREAMBLE_23))
    return resolve_ref(load(path), found, collection=source)


def written(loaded: Loaded, path: Path, **options) -> Path:
    return write_vcf(loaded, path, resolved_only=True, **options).path


def readable(row: Row) -> bool:
    """A row whose REF `resolve_ref` reads: a called base, not `--` or a 23andMe indel code."""
    return set(row[3]) <= set("ACGT") and bool(row[3])


def called(row: Row) -> set[str]:
    return set(row[3])


def strand_safe(row: Row) -> bool:
    """F13: a het carrying the reference base reads the same on either strand."""
    return len(row[3]) == 2 and len(called(row)) == 2 and base(row[1], row[2]) in called(row)


def expected_palindromes(rows: list[Row]) -> dict[str, int]:
    counts = Counter({"het": 0, "hom_alt": 0, "alt_unstated": 0})
    for row in rows:
        ref = base(row[1], row[2])
        if called(row) == {ref}:
            counts["alt_unstated"] += 1
        elif called(row) == {ref, COMPLEMENT[ref]}:
            counts["het"] += 1
        elif called(row) == {COMPLEMENT[ref]}:
            counts["hom_alt"] += 1
    return dict(counts)


@pytest.fixture(scope="module")
def forward_37(tmp_path_factory: pytest.TempPathFactory, store) -> tuple[Loaded, Path, dict[str, Outcome]]:
    tmp_path = tmp_path_factory.mktemp("forward_37")
    loaded = read_chip(tmp_path, store, FORWARD)
    path = written(loaded, tmp_path / "chip.vcf.gz")
    return loaded, path, check_vcf(path, "GRCh37", reference_bases(SOURCE))


@pytest.fixture(scope="module")
def lifted_38(tmp_path_factory: pytest.TempPathFactory, store, forward_37) -> Loaded:
    found, _, target, chain = store
    return liftover(forward_37[0], chain, found, target=target).loaded


def test_a_chip_is_written_only_with_resolved_only(tmp_path: Path, store) -> None:
    """The typical path needs `resolved_only=True`: a no-call and an indel code have no REF to write."""
    loaded = read_chip(tmp_path, store, FORWARD)
    with pytest.raises(UnwritableRowsError) as raised:
        write_vcf(loaded, tmp_path / "chip.vcf.gz")
    assert raised.value.counts == {"ref_unresolved": sum(not readable(row) for row in FORWARD)}


FILE_RULES_MET = (
    "bgzf", "sorted", "gt_in_format", "ref_differs_from_alt", "unique_sample_ids", "chr_prefix_for_build",
    "filtered", "indel", "monomorphic", "site_call_rate", "ref_matches_reference", "chrx_ploidy", "chrx_mixed_genotypes",
)  # fmt: skip


def test_a_written_chip_meets_the_rules_a_file_can_meet(forward_37) -> None:
    _, _, report = forward_37
    assert {rule: report[rule].met for rule in FILE_RULES_MET} == dict.fromkeys(FILE_RULES_MET, True)


def test_the_rules_that_need_the_panel_are_never_reported_as_passed(forward_37) -> None:
    _, _, report = forward_37
    assert {rule: report[rule].met for rule in [*PANEL_RULES, "palindromic_sites"]} == dict.fromkeys([*PANEL_RULES, "palindromic_sites"])


def test_the_report_counts_what_the_export_holds(forward_37) -> None:
    """Each count, read here from the export's rows against the test's sequences."""
    _, _, report = forward_37
    rows = [row for row in FORWARD if readable(row)]
    hom_ref = sum(called(row) == {base(row[1], row[2])} for row in rows)
    two_alts = sum(len(called(row)) == 2 and base(row[1], row[2]) not in called(row) for row in rows)
    positions = Counter(("X" if row[1] == "25" else row[1], row[2]) for row in rows)
    x_nonpar = [row for row in rows if row[1] == "X" and row[2] >= NONPAR_37]
    assert report["multiallelic"].counts == {"records_excluded": two_alts}
    assert report["invalid_alleles"].counts == {"records_excluded": hom_ref, "alt_unstated": hom_ref}
    assert report["no_duplicate_positions"].counts == {"records_at_a_previous_position": sum(n - 1 for n in positions.values()), "dup_filter": 0}
    assert report["one_chromosome"].counts == {"chromosomes": len({row[1] for row in rows})}
    assert report["valid_chromosome"].counts == {"records_on_XY": sum(row[1] == "25" for row in rows)}
    assert report["samples_at_least_20"].counts == {"samples": 1}
    assert report["chrx_ploidy"].counts == {"nonpar_records": len(x_nonpar), "haploid_calls": len(x_nonpar), "non_haploid_after_haploid": 0}
    assert report["palindromic_sites"].counts == expected_palindromes(rows)
    assert report["vcf_version_at_most_4_2"].counts == {"v4_3": 1}


NOT_MET_37 = {
    "one_chromosome": "write_vcf writes a facet's every contig into one file; the servers take one chromosome per file (VcfFileUtil.load)",
    "valid_chromosome": "the array's pseudo-autosomal code 25 is written as contig XY, which the servers' chromosome list does not hold; as X it would be read into PAR1 or PAR2",
    "samples_at_least_20": "a person's chip is one sample; MIS2's and TOPMed's default minimum is 20, so it is uploaded beside others or not at all",
    "vcf_version_at_most_4_2": "write_vcf writes VCF 4.3 (4.5 with PSL); TOPMed's docs ask for 4.2 or lower (MIS2 states no version, and the code read checks none)",
    "invalid_alleles": "a chip's hom-ref call is written ALT `.`, which the QC drops as an invalid allele (GenomicTools.isValid): an ALT needs the probe's alleles",
    "multiallelic": "a call of two non-reference bases is written as a site with two ALTs, which the QC excludes",
    "no_duplicate_positions": "two probes at one position are written as two records (charter 8 keeps both); the docs and checkVCF flag the second",
}


@pytest.mark.parametrize("rule", [pytest.param(rule, marks=pytest.mark.xfail(strict=True, reason=reason)) for rule, reason in NOT_MET_37.items()])
def test_a_written_chip_meets_the_rule(forward_37, rule: str) -> None:
    _, _, report = forward_37
    assert report[rule].met is True, report[rule].counts


@pytest.mark.xfail(strict=True, reason="contig_names defaults to stored names, which on a GRCh38 facet is `1`, not the `chr1` the servers ask of GRCh38")
def test_a_lifted_chip_written_by_default_has_its_builds_names(tmp_path: Path, lifted_38) -> None:
    report = check_vcf(written(lifted_38, tmp_path / "chip.vcf.gz"), "GRCh38")
    assert report["chr_prefix_for_build"].met is True


def test_a_lifted_chip_with_chr_names_meets_the_build_rules(tmp_path: Path, lifted_38, forward_37) -> None:
    path = written(lifted_38, tmp_path / "chip.vcf.gz", contig_names=ContigNames.CHR)
    report = check_vcf(path, "GRCh38", reference_bases(TARGET))
    rules = ("chr_prefix_for_build", "sorted", "ref_matches_reference", "chrx_ploidy", "chrx_mixed_genotypes")
    assert {rule: report[rule].met for rule in rules} == dict.fromkeys(rules, True)
    lifted = {(r.chrom, r.pos) for r in read_vcf_text(path).records}
    source = {(f"chr{r.chrom}", lifted_pos(r.chrom, r.pos)) for r in read_vcf_text(forward_37[1]).records if r.chrom != "XY"}
    assert lifted == source, "every row but the XY ones lands where the construction puts it (no chain covers XY)"
    assert report["chrx_ploidy"].counts["nonpar_records"] == sum(row[1] == "X" for row in FORWARD), "still nonPAR on GRCh38"
    on_minus = [r for r in read_vcf_text(path).records if r.chrom == "chr2"]
    assert palindromes(on_minus) == expected_palindromes([row for row in FORWARD if row[1] == "2"]), "a minus chain keeps a palindrome one"


def test_one_contig_of_the_frame_writes_one_chromosome(tmp_path: Path, lifted_38) -> None:
    """Filtering the calls frame to one contig is how one file per chromosome is written today."""
    for chrom in ("1", "X"):
        one = lifted_38.model_copy(update={"calls": lifted_38.calls.filter(pl.col("chrom") == chrom)})
        report = check_vcf(written(one, tmp_path / f"chr{chrom}.vcf.gz", contig_names=ContigNames.CHR), "GRCh38", reference_bases(TARGET))
        rules = ("one_chromosome", "valid_chromosome", "chr_prefix_for_build", "ref_matches_reference")
        assert {rule: report[rule].met for rule in rules} == dict.fromkeys(rules, True), chrom


def test_a_diploid_x_call_after_haploid_ones_is_written_as_held_and_fails_the_job(tmp_path: Path, store) -> None:
    """The export's ploidy reaches the file unchanged (Decision 2: hemizygous X is kept as written),
    so one diploid nonPAR call after haploid ones fails the server's ploidy check."""
    pos = NONPAR_37 + 300
    het = ("rs9100099", "X", pos, base("X", pos) + other(base("X", pos)))
    path = written(read_chip(tmp_path, store, [*FORWARD, het]), tmp_path / "chip.vcf.gz")
    report = check_vcf(path, "GRCh37")
    assert report["chrx_ploidy"].met is False
    assert report["chrx_ploidy"].counts["non_haploid_after_haploid"] == 1


def test_the_helper_sees_what_a_server_refuses(tmp_path: Path, forward_37) -> None:
    """The rules a written chip meets are not met by construction of the helper: a plain-text copy
    with two records swapped, one filtered and one REF changed fails each of them."""
    loaded, _, _ = forward_37
    plain = written(loaded, tmp_path / "chip.vcf", compress=False)
    lines = plain.read_text().splitlines()
    first = next(n for n, line in enumerate(lines) if not line.startswith("#"))
    lines[first], lines[first + 1] = lines[first + 1], lines[first]
    # A record with an ALT, since the QC drops an ALT `.` as an invalid allele before reading FILTER.
    with_alt = [n for n, line in enumerate(lines) if n > first + 1 and not line.startswith("#") and line.split("\t")[4] in tuple("ACGT")]
    tiles = lines[with_alt[0]].split("\t")
    tiles[6] = "LowQual"
    lines[with_alt[0]] = "\t".join(tiles)
    tiles = lines[with_alt[1]].split("\t")
    tiles[3] = other(tiles[3], tiles[4])
    lines[with_alt[1]] = "\t".join(tiles)
    plain.write_text("\n".join(lines) + "\n")
    report = check_vcf(plain, "GRCh37", reference_bases(SOURCE))
    assert {rule: report[rule].met for rule in ("bgzf", "sorted", "filtered", "ref_matches_reference")} == dict.fromkeys(
        ("bgzf", "sorted", "filtered", "ref_matches_reference"), False
    )
    assert (report["sorted"].counts, report["filtered"].counts, report["ref_matches_reference"].counts) == (
        {"records_out_of_order": 1},
        {"records_excluded": 1},
        {"records_mismatched": 1, "records_unread": 0},
    )


def test_a_build_the_job_does_not_name_is_not_checked(forward_37) -> None:
    report = check_vcf(forward_37[1])
    assert [report[rule].met for rule in ("chr_prefix_for_build", "chrx_ploidy", "chrx_mixed_genotypes", "ref_matches_reference")] == [None] * 4


# The A/T and C/G gap (RM6; whether to settle it is desk #212). A het of a palindromic probe reads the
# same on either strand; a homozygous one does not, and the server keeps it unflipped either way.


def test_a_forward_chip_carries_its_palindromes_and_withholds_none(forward_37) -> None:
    loaded, path, report = forward_37
    rows = [row for row in FORWARD if readable(row)]
    assert loaded.facets.collect()["strand_check"].cast(pl.String).to_list() == ["forward"]
    assert report["palindromic_sites"].counts == expected_palindromes(rows)
    assert report["palindromic_sites"].counts["hom_alt"] > 0, "the fixture holds the case that depends on the strand"
    assert strand_withheld(loaded.calls.collect()) == {"ref_base_only": 0, "complement_only": 0, "other": 0}


def test_a_mixed_strand_chip_carries_only_strand_safe_palindromes(tmp_path: Path, store) -> None:
    loaded = read_chip(tmp_path, store, MIXED)
    assert loaded.facets.collect()["strand_check"].cast(pl.String).to_list() == ["mixed"]
    result = write_vcf(loaded, tmp_path / "chip.vcf.gz", resolved_only=True)
    rows = [row for row in MIXED if readable(row)]
    carried = [row for row in rows if strand_safe(row)]
    withheld = [row for row in rows if not strand_safe(row)]
    report = check_vcf(result.path, "GRCh37", reference_bases(SOURCE))
    assert report["palindromic_sites"].counts == {**expected_palindromes(carried), "alt_unstated": 0}
    assert report["palindromic_sites"].counts["hom_alt"] == 0
    expected = Counter({"ref_base_only": 0, "complement_only": 0, "other": 0})
    for row in withheld:
        ref = base(row[1], row[2])
        expected["ref_base_only" if called(row) == {ref} else "complement_only" if called(row) == {COMPLEMENT[ref]} else "other"] += 1
    assert strand_withheld(loaded.calls.collect()) == dict(expected)
    assert result.rows_skipped == {"ref_unresolved": sum(not readable(row) for row in MIXED), "alts_unresolved": len(withheld)}


# Opt-in: a real chip export under JUST_VCF_SAMPLES_DIR, against a real GRCh37 store. Aggregates only.


def store_bases(found: ReferenceStore, digest: str) -> RefBases:
    by_chrom = found.contigs_by_chrom(digest)

    def read(chrom: str, spans: list[tuple[int, int]]) -> list[str | None]:
        name = chrom.removeprefix("chr")
        contigs = by_chrom.get({"M": "MT"}.get(name, name), [])
        if len(contigs) != 1:
            return [None] * len(spans)
        (contig,) = contigs
        inside = [(pos - 1, pos - 1 + length) for pos, length in spans if pos - 1 + length <= contig.length]
        fetched = iter(found.fetch(contig, inside))
        return [next(fetched) if pos - 1 + length <= contig.length else None for pos, length in spans]

    return read


def test_the_opt_in_tests_store_reader_reads_the_test_sequences(store, forward_37) -> None:
    found, source, _, _ = store
    records = read_vcf_text(forward_37[1]).records
    for chrom in ("1", "X", "XY"):
        spans = [(r.pos, len(r.ref)) for r in records if r.chrom == chrom] + [(len(SOURCE["X"]), 2)]
        assert store_bases(found, source)(chrom, spans) == reference_bases(SOURCE)(chrom, spans), chrom


def real_chip(pgp_participants: Path) -> Path:
    candidates = sorted(pgp_participants.glob("*/upload-*/genome*Full*.txt"))[:40]
    found = next((p for p in candidates if sniff(p).format == "23andme"), None)
    if found is None:
        pytest.skip("no plain 23andMe export among the first candidates")
    return found


def summary(report: dict[str, Outcome]) -> str:
    return "\n".join(f"  {rule}: {outcome.met} {outcome.counts}" for rule, outcome in report.items())


def test_a_real_chip_export_against_the_server_rules(pgp_participants: Path, test_reference_dir: Path, tmp_path: Path) -> None:
    loaded = load(real_chip(pgp_participants))
    if loaded.facets.collect()["assembly_label"][0] != "GRCh37":
        pytest.skip("the export found is not on GRCh37")
    found = ReferenceStore.open(test_reference_dir)
    digests = [d for d in found.collections() if found.assembly_label(d) == "GRCh37"]
    if len(digests) != 1:
        pytest.skip(f"{test_reference_dir} holds no single GRCh37 collection")
    resolved = resolve_ref(loaded, found, collection=digests[0])
    result = write_vcf(resolved, tmp_path / "chip.vcf.gz", resolved_only=True)
    report = check_vcf(result.path, "GRCh37", store_bases(found, digests[0]))
    strand = resolved.facets.collect()["strand_check"].cast(pl.String)[0]
    print(f"\nwritten {result.rows_written}, skipped {result.rows_skipped}, strand {strand}, withheld {strand_withheld(resolved.calls.collect())}\n{summary(report)}")
    rules = ("bgzf", "sorted", "gt_in_format", "ref_differs_from_alt", "filtered", "indel", "site_call_rate", "ref_matches_reference")
    assert {rule: report[rule].met for rule in rules} == dict.fromkeys(rules, True)
