"""GRCh37's and GRCh38's sequence names from NCBI's assembly reports (F19).

Expected names and lengths are read from the reports carried in `assets/`, independently of the
package table.
"""

from pathlib import Path

import pytest

from just_vcf import load
from just_vcf.assembly import infer_assembly
from just_vcf.contig_names import PRIMARY_ACCESSIONS, SEQUENCES, scaffold_names
from just_vcf.schema import CONTIG_ALIASES, canonical_chrom

REPORTS = Path(__file__).resolve().parent.parent / "assets" / "ncbi_assembly_reports"
REPORT_FILES = {
    "GRCh37": REPORTS / "GCA_000001405.14_GRCh37.p13_assembly_report.txt",
    "GRCh38": REPORTS / "GCA_000001405.29_GRCh38.p14_assembly_report.txt",
}


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
            refseq = row["RefSeq-Accn"] if row["Relationship"] == "=" else "na"
            expected.append(
                (
                    assembly,
                    row["Sequence-Role"],
                    None if row["Assigned-Molecule"] == "na" else row["Assigned-Molecule"],
                    row["GenBank-Accn"],
                    None if refseq == "na" else refseq,
                    None if row["UCSC-style-name"] == "na" else row["UCSC-style-name"],
                    int(row["Sequence-Length"]),
                )
            )
    got = [(s.assembly, s.role, s.molecule, s.genbank, s.refseq, s.ucsc, s.length) for s in SEQUENCES]
    assert got == expected
    # Where the report says the RefSeq sequence is not the GenBank one, it gives no RefSeq accession.
    assert all(row["RefSeq-Accn"] == "na" for a in REPORT_FILES for row in report(a) if row["Relationship"] != "=")


def test_primary_accessions_are_the_assembled_molecules() -> None:
    expected: dict[str, str] = {}
    for assembly in REPORT_FILES:
        for row in report(assembly):
            if row["Sequence-Role"] == "assembled-molecule":
                expected[row["GenBank-Accn"]] = row["Assigned-Molecule"]
                expected[row["RefSeq-Accn"]] = row["Assigned-Molecule"]
    assert PRIMARY_ACCESSIONS == expected
    assert len(expected) == 2 * 25 * 2 - 2  # MT is J01415.2 / NC_012920.1 in both reports
    for name, chrom in expected.items():
        assert canonical_chrom(name) == chrom
        assert CONTIG_ALIASES[name] == chrom


def test_one_scaffold_name_is_one_sequence() -> None:
    """Across both reports a scaffold's name points at one GenBank accession, and an accession has one
    length; only the assembled molecules' UCSC names (`chr1`…) differ between the two."""
    accession_of: dict[str, set[str]] = {}
    length_of: dict[str, set[int]] = {}
    for assembly in REPORT_FILES:
        for row in report(assembly):
            length_of.setdefault(row["GenBank-Accn"], set()).add(int(row["Sequence-Length"]))
            if row["Sequence-Role"] != "assembled-molecule":
                for name in (row["GenBank-Accn"], row["RefSeq-Accn"], row["UCSC-style-name"]):
                    if name != "na":
                        accession_of.setdefault(name, set()).add(row["GenBank-Accn"])
    assert all(len(found) == 1 for found in accession_of.values())
    assert all(len(found) == 1 for found in length_of.values())


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


@pytest.mark.parametrize("assembly", list(REPORT_FILES))
@pytest.mark.parametrize("column", ["GenBank-Accn", "RefSeq-Accn"])
def test_an_accession_named_header_names_its_assembly(tmp_path: Path, assembly: str, column: str) -> None:
    rows = [row for row in report(assembly) if row["Sequence-Role"] == "assembled-molecule"]
    header = [{"name": row[column], "length": int(row["Sequence-Length"]), "md5": None} for row in rows]
    inferred = infer_assembly(header)
    assert inferred.label == assembly and inferred.matched == {assembly: 24}
    first = header[0]
    lines = ["##fileformat=VCFv4.2", *[f"##contig=<ID={c['name']},length={c['length']}>" for c in header]]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    lines.append(f"{first['name']}\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
    path = tmp_path / "accessions.vcf"
    path.write_text("\n".join(lines) + "\n")
    facet = load(path).facets.collect().row(0, named=True)
    assert facet["assembly_label"] == assembly
    assert {entry["source"]: entry["chrom"] for entry in facet["contig_map"]} == {first["name"]: "1"}
