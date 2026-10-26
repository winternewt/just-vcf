"""F19: NCBI36's, GRCh37's and GRCh38's sequence names from NCBI's assembly reports, and their use where
a contig name meets a reference store.

Expected names and lengths are read from the reports carried in `assets/`, independently of the
package table; sequences are synthetic, of the reports' lengths. NCBI36's report is RefSeq only: every
row has GenBank `na` and Relationship `<>`, so its RefSeq accession is the sequence's one accession.
"""

import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, liftover, load, mint_vrs_ids, normalize, resolve_ref
from just_vcf.assembly import infer_assembly
from just_vcf.contig_names import PRIMARY_ACCESSIONS, SEQUENCES, scaffold_names
from just_vcf.schema import CONTIG_ALIASES, canonical_chrom

REPORTS = Path(__file__).resolve().parent.parent / "assets" / "ncbi_assembly_reports"
REPORT_FILES = {
    "GRCh37": REPORTS / "GCA_000001405.14_GRCh37.p13_assembly_report.txt",
    "GRCh38": REPORTS / "GCA_000001405.29_GRCh38.p14_assembly_report.txt",
    "GRCh36": REPORTS / "GCF_000001405.12_NCBI36_assembly_report.txt",
}


def carried_refseq(row: dict[str, str]) -> str | None:
    """The RefSeq accession the table carries for a report row: where the report says it is the GenBank
    sequence (`=`), or where the report names no GenBank sequence for it to differ from (NCBI36)."""
    if row["RefSeq-Accn"] == "na":
        return None
    return row["RefSeq-Accn"] if row["Relationship"] == "=" or row["GenBank-Accn"] == "na" else None


def report(assembly: str) -> list[dict[str, str]]:
    """The report's sequence rows, keyed by its own column names."""
    lines = REPORT_FILES[assembly].read_text().splitlines()
    columns = next(line for line in lines if line.startswith("# Sequence-Name")).removeprefix("# ").split("\t")
    return [dict(zip(columns, line.split("\t"), strict=True)) for line in lines if line and not line.startswith("#")]


def scaffold(assembly: str, ucsc: str) -> dict[str, str]:
    (row,) = [row for row in report(assembly) if row["UCSC-style-name"] == ucsc]
    return row


# --- the table ------------------------------------------------------------------------------------


def test_the_table_is_the_assembly_reports() -> None:
    expected = []
    for assembly in REPORT_FILES:
        for row in report(assembly):
            expected.append(
                (
                    assembly,
                    row["Sequence-Role"],
                    None if row["Assigned-Molecule"] == "na" else row["Assigned-Molecule"],
                    None if row["GenBank-Accn"] == "na" else row["GenBank-Accn"],
                    carried_refseq(row),
                    None if row["UCSC-style-name"] == "na" else row["UCSC-style-name"],
                    int(row["Sequence-Length"]),
                )
            )
    got = [(s.assembly, s.role, s.molecule, s.genbank, s.refseq, s.ucsc, s.length) for s in SEQUENCES]
    assert got == expected
    # Where the report says the RefSeq sequence is not the GenBank one, it gives no RefSeq accession.
    assert all(
        row["RefSeq-Accn"] == "na" for a in REPORT_FILES for row in report(a) if row["Relationship"] != "=" and row["GenBank-Accn"] != "na"
    )
    # NCBI36's report names no GenBank sequence at all, and a RefSeq accession on every row.
    assert {(row["GenBank-Accn"], row["Relationship"]) for row in report("GRCh36")} == {("na", "<>")}
    assert all(row["RefSeq-Accn"] != "na" for row in report("GRCh36"))


def test_primary_accessions_are_the_assembled_molecules() -> None:
    expected: dict[str, str] = {}
    for assembly in REPORT_FILES:
        for row in report(assembly):
            if row["Sequence-Role"] == "assembled-molecule":
                for name in (row["GenBank-Accn"], carried_refseq(row)):
                    if name not in (None, "na"):
                        expected[name] = row["Assigned-Molecule"]
    assert PRIMARY_ACCESSIONS == expected
    # GRCh37 and GRCh38: 25 molecules, GenBank and RefSeq, MT (J01415.2 / NC_012920.1) in both reports.
    # NCBI36: 24 molecules (no mitochondrion), RefSeq only.
    assert len(expected) == 2 * 25 * 2 - 2 + 24
    for name, chrom in expected.items():
        assert canonical_chrom(name) == chrom
        assert CONTIG_ALIASES[name] == chrom


def test_a_primary_accession_names_one_build() -> None:
    """An accession with its version is one sequence: no assembled molecule's accession is shared by two
    builds, except the rCRS mitochondrion GRCh37 and GRCh38 both carry, at one length. Each chromosome's
    RefSeq accession keeps its stem across builds and changes its version (`NC_000001.9`, `.10`, `.11`)."""
    builds_of: dict[str, set[str]] = {}
    lengths_of: dict[str, set[int]] = {}
    refseq_by_chrom: dict[str, dict[str, str]] = {}
    for assembly in REPORT_FILES:
        for row in report(assembly):
            if row["Sequence-Role"] != "assembled-molecule":
                continue
            for name in (row["GenBank-Accn"], carried_refseq(row)):
                if name not in (None, "na"):
                    builds_of.setdefault(name, set()).add(assembly)
                    lengths_of.setdefault(name, set()).add(int(row["Sequence-Length"]))
            refseq_by_chrom.setdefault(row["Assigned-Molecule"], {})[assembly] = row["RefSeq-Accn"]
    shared = {name for name, builds in builds_of.items() if len(builds) > 1}
    assert shared == {"J01415.2", "NC_012920.1"}
    assert builds_of["NC_012920.1"] == {"GRCh37", "GRCh38"}
    assert all(len(lengths) == 1 for lengths in lengths_of.values())
    for chrom, by_build in refseq_by_chrom.items():
        if chrom == "MT":
            continue
        assert set(by_build) == set(REPORT_FILES)
        assert len({accession.split(".")[0] for accession in by_build.values()}) == 1
        assert len(set(by_build.values())) == len(REPORT_FILES)
    assert {refseq_by_chrom["1"][a] for a in REPORT_FILES} == {"NC_000001.9", "NC_000001.10", "NC_000001.11"}


def test_one_scaffold_name_is_one_sequence() -> None:
    """Across the reports a scaffold's name points at one sequence, and a sequence has one length; only
    the assembled molecules' UCSC names (`chr1`…) differ between them. A sequence is its GenBank
    accession where any report gives one: NCBI36's report names none, and its scaffolds that GRCh37 or
    GRCh38 kept (`NT_113878.1`, GRCh37's `GL000191.1`) are found by their RefSeq accession."""
    genbank_of = {
        row["RefSeq-Accn"]: row["GenBank-Accn"]
        for assembly in REPORT_FILES
        for row in report(assembly)
        if row["Relationship"] == "=" and row["GenBank-Accn"] != "na"
    }

    def sequence(row: dict[str, str]) -> str:
        return row["GenBank-Accn"] if row["GenBank-Accn"] != "na" else genbank_of.get(row["RefSeq-Accn"], row["RefSeq-Accn"])

    accession_of: dict[str, set[str]] = {}
    length_of: dict[str, set[int]] = {}
    for assembly in REPORT_FILES:
        for row in report(assembly):
            length_of.setdefault(sequence(row), set()).add(int(row["Sequence-Length"]))
            if row["Sequence-Role"] != "assembled-molecule":
                for name in (row["GenBank-Accn"], carried_refseq(row), row["UCSC-style-name"]):
                    if name not in (None, "na"):
                        accession_of.setdefault(name, set()).add(sequence(row))
    assert all(len(found) == 1 for found in accession_of.values())
    assert all(len(found) == 1 for found in length_of.values())
    # The check reaches NCBI36: some of its scaffolds are GRCh37's too, under one sequence.
    ncbi36 = {sequence(row) for row in report("GRCh36") if row["Sequence-Role"] != "assembled-molecule"}
    later = {sequence(row) for a in ("GRCh37", "GRCh38") for row in report(a) if row["Sequence-Role"] != "assembled-molecule"}
    assert ncbi36 & later and ncbi36 - later


def test_scaffold_names_need_the_reports_length() -> None:
    row = scaffold("GRCh38", "chr1_KI270706v1_random")
    length = int(row["Sequence-Length"])
    names = {"KI270706.1", "NT_187361.1", "chr1_KI270706v1_random"}
    for name in names:
        assert scaffold_names(name, length) == names - {name}
        assert scaffold_names(name, length + 1) == frozenset()
    # One accession, two builds, two UCSC names: GRCh37 and GRCh38 place GL000194.1 differently.
    old, new = scaffold("GRCh37", "chr4_gl000194_random"), scaffold("GRCh38", "chr14_GL000194v1_random")
    assert old["GenBank-Accn"] == new["GenBank-Accn"] == "GL000194.1"
    assert {"chr4_gl000194_random", "chr14_GL000194v1_random"} <= scaffold_names("GL000194.1", int(old["Sequence-Length"]))
    # Nothing by pattern: a name the reports do not list, or another case, has no alias.
    assert scaffold_names("chr1_KI270706v2_random", length) == frozenset()
    assert scaffold_names("CHR4_GL000194_RANDOM", int(old["Sequence-Length"])) == frozenset()
    # The assembled molecules are not scaffolds: the reports' MT row (`chrM`, rCRS) vouches for nothing,
    # since hg19's chrM is another sequence.
    assert scaffold_names("chrM", 16569) == scaffold_names("J01415.2", 16569) == frozenset()
    # NCBI36's scaffolds have a RefSeq name, and two a UCSC one too (below). One GRCh37 kept gains
    # GRCh37's names; one it did not, with no UCSC name, has no other name.
    kept = next(row for row in report("GRCh36") if row["RefSeq-Accn"] == scaffold("GRCh37", "chr1_gl000191_random")["RefSeq-Accn"])
    assert scaffold_names(kept["RefSeq-Accn"], int(kept["Sequence-Length"])) == {"GL000191.1", "chr1_gl000191_random"}
    later = {row["RefSeq-Accn"] for a in ("GRCh37", "GRCh38") for row in report(a)}
    alone = next(
        row
        for row in report("GRCh36")
        if row["Sequence-Role"] != "assembled-molecule" and row["RefSeq-Accn"] not in later and row["UCSC-style-name"] == "na"
    )
    assert scaffold_names(alone["RefSeq-Accn"], int(alone["Sequence-Length"])) == frozenset()


def test_ncbi36s_two_ucsc_named_scaffolds_meet_their_refseq_accessions() -> None:
    # NCBI36's report gives a UCSC name to two scaffolds, both alt haplotypes, and to no other: each
    # name and its RefSeq accession are the other's alias, at the report's length only.
    named = [row for row in report("GRCh36") if row["Sequence-Role"] != "assembled-molecule" and row["UCSC-style-name"] != "na"]
    assert {(row["UCSC-style-name"], row["RefSeq-Accn"]) for row in named} == {("chr6_cox_hap1", "NT_113891.1"), ("chr22_h2_hap1", "NT_113959.1")}
    for row in named:
        ucsc, refseq, length = row["UCSC-style-name"], row["RefSeq-Accn"], int(row["Sequence-Length"])
        assert scaffold_names(ucsc, length) == {refseq}
        assert scaffold_names(refseq, length) == {ucsc}
        assert scaffold_names(ucsc, length + 1) == frozenset()
    # GRCh37's chr6_cox_hap2 is NT_113891.2, another sequence of another length: no build-36 name meets it.
    grch37 = scaffold("GRCh37", "chr6_cox_hap2")
    assert scaffold_names("chr6_cox_hap1", int(grch37["Sequence-Length"])) == frozenset()


ACCESSION_COLUMNS = [
    (assembly, column)
    for assembly in REPORT_FILES
    for column in ("GenBank-Accn", "RefSeq-Accn")
    if any(row[column] != "na" for row in report(assembly) if row["Sequence-Role"] == "assembled-molecule")
]
"""Every report's accession columns that name its assembled molecules: NCBI36's has no GenBank one."""


def header_vcf(tmp_path: Path, header: list[dict], first: str) -> Path:
    lines = ["##fileformat=VCFv4.2", *[f"##contig=<ID={c['name']},length={c['length']}>" for c in header]]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    lines.append(f"{first}\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
    path = tmp_path / "accessions.vcf"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_every_report_has_accession_columns_to_test() -> None:
    assert ("GRCh36", "RefSeq-Accn") in ACCESSION_COLUMNS and ("GRCh36", "GenBank-Accn") not in ACCESSION_COLUMNS
    assert len(ACCESSION_COLUMNS) == 2 * len(REPORT_FILES) - 1


@pytest.mark.parametrize(("assembly", "column"), ACCESSION_COLUMNS)
def test_an_accession_named_header_names_its_assembly(tmp_path: Path, assembly: str, column: str) -> None:
    rows = [row for row in report(assembly) if row["Sequence-Role"] == "assembled-molecule"]
    header = [{"name": row[column], "length": int(row["Sequence-Length"]), "md5": None} for row in rows]
    inferred = infer_assembly(header)
    assert inferred.label == assembly and inferred.matched == {assembly: 24}
    first = header[0]
    facet = load(header_vcf(tmp_path, header, first["name"])).facets.collect().row(0, named=True)
    assert facet["assembly_label"] == assembly
    assert {entry["source"]: entry["chrom"] for entry in facet["contig_map"]} == {first["name"]: "1"}


def test_a_build_36_header_of_mixed_names_reads_canonical_chroms(tmp_path: Path) -> None:
    """A build-36 header naming some chromosomes by RefSeq accession, the others `chrN` or bare: every
    one reads as its canonical `chrom`, and the lengths name GRCh36."""
    rows = [row for row in report("GRCh36") if row["Sequence-Role"] == "assembled-molecule"]
    spellings = [row["RefSeq-Accn"] if index % 3 == 0 else f"chr{row['Assigned-Molecule']}" if index % 3 == 1 else row["Assigned-Molecule"] for index, row in enumerate(rows)]
    header = [{"name": name, "length": int(row["Sequence-Length"]), "md5": None} for name, row in zip(spellings, rows, strict=True)]
    assert {canonical_chrom(name) for name in spellings} == {row["Assigned-Molecule"] for row in rows}
    inferred = infer_assembly(header)
    assert inferred.label == "GRCh36" and inferred.matched == {"GRCh36": len(rows)}
    accession = rows[0]["RefSeq-Accn"]
    assert accession.startswith("NC_") and spellings[0] == accession
    loaded = load(header_vcf(tmp_path, header, accession))
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["assembly_label"] == "GRCh36"
    assert {entry["source"]: entry["chrom"] for entry in facet["contig_map"]} == {accession: rows[0]["Assigned-Molecule"]}
    assert loaded.calls.collect()["chrom"].to_list() == [rows[0]["Assigned-Molecule"]]


# --- matching against a reference store -----------------------------------------------------------

rng = random.Random(19)
SMALL = scaffold("GRCh38", "chrUn_KI270394v1")  # the smallest GRCh38 scaffold
OLD = scaffold("GRCh37", "chr18_gl000207_random")  # the smallest GRCh37 one
PRIMARY = "".join(rng.choice("ACGT") for _ in range(3000))
SMALL_SEQ = "".join(rng.choice("ACGT") for _ in range(int(SMALL["Sequence-Length"])))
OLD_SEQ = "".join(rng.choice("ACGT") for _ in range(int(OLD["Sequence-Length"])))


def fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()))
    return path


def store_of(tmp_path: Path, seqs: dict[str, str], name: str = "store") -> ReferenceStore:
    store = ReferenceStore.open(tmp_path / name)
    store.add_fasta(fasta(tmp_path / f"{name}.fa", seqs))
    return store


def vcf(tmp_path: Path, contigs: dict[str, int | None], records: list[tuple[str, int, str, str]], name: str = "calls.vcf") -> Path:
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={n}>" if n is not None else f"##contig=<ID={c}>" for c, n in contigs.items()]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    lines += [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1" for c, p, r, a in records]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def alt(base: str) -> str:
    return next(b for b in "ACGT" if b != base)


ON_1 = ("1", 5, PRIMARY[4], alt(PRIMARY[4]))
"""A row on a primary contig, so the facet has something readable and is not refused outright."""


def scaffold_records(name: str, seq: str) -> list[tuple[str, int, str, str]]:
    return [(name, pos, seq[pos - 1], alt(seq[pos - 1])) for pos in (1, 37, len(seq))]


@pytest.mark.parametrize(
    ("store_name", "vcf_name"),
    [
        ("KI270394.1", "chrUn_KI270394v1"),  # Ensembl's FASTA, a UCSC-named VCF (F19's PGP-HMS case)
        ("chrUn_KI270394v1", "KI270394.1"),  # UCSC's FASTA, an Ensembl-named VCF
        ("KI270394.1", "NT_187479.1"),  # a RefSeq-named VCF
    ],
)
def test_a_scaffold_is_read_under_any_of_its_report_names(tmp_path: Path, store_name: str, vcf_name: str) -> None:
    assert {store_name, vcf_name} <= {SMALL["GenBank-Accn"], SMALL["RefSeq-Accn"], SMALL["UCSC-style-name"]}
    store = store_of(tmp_path, {"1": PRIMARY, store_name: SMALL_SEQ})
    records = scaffold_records(vcf_name, SMALL_SEQ)
    path = vcf(tmp_path, {"1": len(PRIMARY), vcf_name: len(SMALL_SEQ)}, [("1", 5, PRIMARY[4], alt(PRIMARY[4])), *records])
    resolved = resolve_ref(load(path), store)
    calls = resolved.calls.collect()
    on_scaffold = calls.filter(pl.col("chrom") == vcf_name)  # the source's spelling is kept
    assert on_scaffold["ref"].to_list() == [r for _, _, r, _ in records]
    assert set(on_scaffold["ref_check"].cast(pl.String)) == {"match"}
    assert not any("left with REF unread" in line for line in resolved.warnings)
    assert any(line.startswith(f"{len(records)} row(s) read on scaffold(s)") and f"{vcf_name} as {store_name}" in line for line in resolved.warnings)


def test_a_grch37_ucsc_scaffold_meets_ensembls_accession(tmp_path: Path) -> None:
    store = store_of(tmp_path, {"1": PRIMARY, OLD["GenBank-Accn"]: OLD_SEQ})
    records = scaffold_records(OLD["UCSC-style-name"], OLD_SEQ)
    path = vcf(tmp_path, {"1": len(PRIMARY), OLD["UCSC-style-name"]: len(OLD_SEQ)}, records)
    calls = resolve_ref(load(path), store).calls.collect()
    assert calls["ref"].to_list() == [r for _, _, r, _ in records]


def test_a_stated_length_that_is_not_the_reports_leaves_the_scaffold_unread(tmp_path: Path) -> None:
    store = store_of(tmp_path, {"1": PRIMARY, "KI270394.1": SMALL_SEQ})
    records = scaffold_records("chrUn_KI270394v1", SMALL_SEQ) + [ON_1]
    path = vcf(tmp_path, {"1": len(PRIMARY), "chrUn_KI270394v1": len(SMALL_SEQ) + 1}, records)
    resolved = resolve_ref(load(path), store)
    assert resolved.calls.collect()["ref"].null_count() == len(records) - 1
    assert any(
        line.startswith(f"{len(records) - 1} row(s) left with REF unread")
        and f"chrUn_KI270394v1 ({len(records) - 1} rows; length stated {len(SMALL_SEQ) + 1}, reference {len(SMALL_SEQ)})" in line
        for line in resolved.warnings
    )


def test_a_reference_contig_of_another_length_is_not_aliased(tmp_path: Path) -> None:
    """The store's `KI270394.1` is one base short of the report's: it is not that sequence, so the
    UCSC name does not reach it, stated length or not."""
    short = SMALL_SEQ[:-1]
    store = store_of(tmp_path, {"1": PRIMARY, "KI270394.1": short})
    records = scaffold_records("chrUn_KI270394v1", short) + [ON_1]
    for stated in (len(short), None):
        path = vcf(tmp_path, {"1": len(PRIMARY), "chrUn_KI270394v1": stated}, records, name=f"calls_{stated}.vcf")
        resolved = resolve_ref(load(path), store)
        assert resolved.calls.collect()["ref"].null_count() == len(records) - 1
        assert any("on contig(s) the reference does not have: chrUn_KI270394v1" in line for line in resolved.warnings)


def test_a_store_holding_both_names_leaves_the_scaffold_unread(tmp_path: Path) -> None:
    store = store_of(tmp_path, {"1": PRIMARY, "KI270394.1": SMALL_SEQ, "chrUn_KI270394v1": SMALL_SEQ})
    records = scaffold_records("NT_187479.1", SMALL_SEQ) + [ON_1]
    path = vcf(tmp_path, {"1": len(PRIMARY), "NT_187479.1": len(SMALL_SEQ)}, records)
    resolved = resolve_ref(load(path), store)
    assert resolved.calls.collect()["ref"].null_count() == len(records) - 1
    assert any("names twice: NT_187479.1 (KI270394.1 and chrUn_KI270394v1)" in line for line in resolved.warnings)


def test_a_resolved_scaffold_is_normalized_and_minted_against_its_store_name(tmp_path: Path) -> None:
    """`normalize` and `mint_vrs_ids` find the scaffold through the same aliases as `resolve_ref`."""
    # Bases 41-48 (1-based) are `ACACACAC`, between two `G`s.
    seq = SMALL_SEQ[:39] + "GACACACACG" + SMALL_SEQ[49:]
    assert len(seq) == len(SMALL_SEQ)
    store = store_of(tmp_path, {"1": PRIMARY, "KI270394.1": seq})
    # One `AC` deleted, written at the repeat's right end (anchor 46): normalize moves it to anchor 40.
    path = vcf(tmp_path, {"1": len(PRIMARY), "chrUn_KI270394v1": len(seq)}, [("chrUn_KI270394v1", 46, seq[45:48], seq[45])])
    normalized = normalize(resolve_ref(load(path), store), store)
    calls = normalized.calls.collect()
    assert calls.select("chrom", "pos", "ref").rows() == [("chrUn_KI270394v1", 40, "GAC")]
    assert calls["alts"].to_list() == [["G"]]
    minted = mint_vrs_ids(normalized, store).calls.collect()
    assert minted["vrs_allele_ids"][0].null_count() == 0


# --- liftover onto a scaffold ----------------------------------------------------------------------

SOURCE = "".join(rng.choice("ACGT") for _ in range(2000))


def lift_case(tmp_path: Path, q_size: int) -> tuple[ReferenceStore, Path, Path]:
    """Source contig 1; the target holds contig 1 (source 0..1000) and `KI270394.1`, onto which the
    chain maps source 1500..1600 under UCSC's name."""
    target_1 = SOURCE[:1000] + "".join(rng.choice("ACGT") for _ in range(500))
    target_scaffold = SMALL_SEQ[:200] + SOURCE[1500:1600] + SMALL_SEQ[300:]
    store = ReferenceStore.open(tmp_path / "store")
    store.add_fasta(fasta(tmp_path / "source.fa", {"1": SOURCE}))
    store.add_fasta(fasta(tmp_path / "target.fa", {"1": target_1, "KI270394.1": target_scaffold}))
    chain = tmp_path / "to_target.chain"
    chain.write_text(
        "\n".join(
            [
                f"chain 1000 1 {len(SOURCE)} + 0 1000 1 {len(target_1)} + 0 1000 1",
                "1000",
                "",
                f"chain 900 1 {len(SOURCE)} + 1500 1600 chrUn_KI270394v1 {q_size} + 200 300 2",
                "100",
                "",
            ]
        )
    )
    records = [("1", 10, SOURCE[9], alt(SOURCE[9])), ("1", 1550, SOURCE[1549], alt(SOURCE[1549]))]
    return store, chain, vcf(tmp_path, {"1": len(SOURCE)}, records)


def test_a_chain_onto_a_ucsc_scaffold_lifts_onto_ensembls_accession(tmp_path: Path) -> None:
    store, chain, path = lift_case(tmp_path, q_size=len(SMALL_SEQ))
    source_digest = next(d for d in store.collections() if any(c.name == "1" and c.length == len(SOURCE) for c in store.contigs(d)))
    resolved = resolve_ref(load(path), store, collection=source_digest)
    target = next(d for d in store.collections() if d != source_digest)
    lifted = liftover(resolved, chain, store, target=target)
    assert lifted.reasons == {}
    calls = lifted.loaded.calls.collect().sort("chrom")
    # The chain's name is kept, as a loaded file's would be (Decision 7); the base is the store's.
    assert calls.select("chrom", "pos").rows() == [("1", 10), ("chrUn_KI270394v1", 200 + (1550 - 1500))]
    assert calls["ref"].to_list() == [SOURCE[9], SOURCE[1549]]


def test_a_chain_scaffold_size_that_is_not_the_reports_is_unmatched(tmp_path: Path) -> None:
    store, chain, path = lift_case(tmp_path, q_size=len(SMALL_SEQ) + 1)
    source_digest = next(d for d in store.collections() if any(c.name == "1" and c.length == len(SOURCE) for c in store.contigs(d)))
    resolved = resolve_ref(load(path), store, collection=source_digest)
    target = next(d for d in store.collections() if d != source_digest)
    lifted = liftover(resolved, chain, store, target=target)
    assert lifted.reasons == {"target_contig_unmatched": 1}
