"""Derive the contig-name table `just_vcf.contig_names` reads, from NCBI's assembly reports (F19).

One row per sequence of each report: the assembly label, the sequence's role, the molecule it is
assigned to, its GenBank and RefSeq accessions, its UCSC-style name and its length, every value as
the report writes it. A RefSeq accession is carried where the report says it is the same sequence as
the GenBank one (`Relationship` `=`), or where the report names no GenBank sequence for it to differ
from (`GenBank-Accn` `na`: NCBI36's report is RefSeq only, every row `<>`, so its RefSeq accession is
the sequence's one accession). A RefSeq accession beside a GenBank one it is not (`<>`) is refused.
`na` becomes an empty cell. The reports are carried verbatim in `assets/ncbi_assembly_reports/`;
docs/PROVENANCE.md names them and their sha256.

    uv run python scripts/contig_names.py \\
        --report GRCh37=assets/ncbi_assembly_reports/GCA_000001405.14_GRCh37.p13_assembly_report.txt \\
        --report GRCh38=assets/ncbi_assembly_reports/GCA_000001405.29_GRCh38.p14_assembly_report.txt \\
        --report GRCh36=assets/ncbi_assembly_reports/GCF_000001405.12_NCBI36_assembly_report.txt \\
        --out src/just_vcf/contig_names.tsv

The label is `just_vcf.assembly`'s (`GRCh36` for NCBI36); the header line records the report's own name.
"""

import logging
from pathlib import Path

import typer

logger = logging.getLogger("contig_names")

COLUMNS = ("assembly", "role", "molecule", "genbank", "refseq", "ucsc", "length")
"""`just_vcf.contig_names.COLUMNS`, which checks the file it reads against its own copy."""


def report_rows(path: Path) -> tuple[str, list[dict[str, str]]]:
    """The report's `# Assembly name:` and its sequence rows, keyed by the report's own column names."""
    lines = path.read_text().splitlines()
    name = next(line for line in lines if line.startswith("# Assembly name:")).split(":", 1)[1].strip()
    columns = next(line for line in lines if line.startswith("# Sequence-Name")).removeprefix("# ").split("\t")
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in lines if line and not line.startswith("#")]
    return name, rows


def table_rows(assembly: str, rows: list[dict[str, str]]) -> list[list[str]]:
    """The table's rows for one report, in the report's order."""
    found = []
    for row in rows:
        same = row["Relationship"] == "=" or row["GenBank-Accn"] == "na"
        if not same and row["RefSeq-Accn"] != "na":
            raise ValueError(f"{row['Sequence-Name']}: a RefSeq accession that is not the GenBank sequence")
        cells = [
            assembly,
            row["Sequence-Role"],
            row["Assigned-Molecule"],
            row["GenBank-Accn"],
            row["RefSeq-Accn"] if same else "na",
            row["UCSC-style-name"],
            row["Sequence-Length"],
        ]
        found.append(["" if cell == "na" else cell for cell in cells])
    return found


def main(
    report: list[str] = typer.Option(..., help="label=path to an NCBI assembly report; repeat per assembly"),
    out: Path = typer.Option(..., help="TSV file to write"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    header = ["# NCBI assembly reports, one row per sequence (docs/PROVENANCE.md); written by scripts/contig_names.py"]
    body: list[list[str]] = []
    for item in report:
        assembly, path = item.split("=", 1)
        name, rows = report_rows(Path(path))
        header.append(f"# {assembly}: {name}, {Path(path).name}")
        body += table_rows(assembly, rows)
        logger.info("%s: %d sequences from %s", assembly, len(rows), name)
    out.write_text("\n".join([*header, "\t".join(COLUMNS), *("\t".join(cells) for cells in body)]) + "\n")
    logger.info("%d rows written to %s", len(body), out)


if __name__ == "__main__":
    typer.run(main)
