"""GRCh37's and GRCh38's sequence names, from NCBI's assembly reports (F19).

One sequence goes by several names. UCSC calls GRCh38's unlocalized, unplaced, alt and patch scaffolds
`chr1_KI270706v1_random`, `chrUn_KI270302v1`, `chr1_KI270762v1_alt`, `chr1_KN196472v1_fix`, and
GRCh37's `chr1_gl000191_random`, `chr6_apd_hap1`; Ensembl and GenBank call the same sequences by their
GenBank accession (`KI270706.1`), RefSeq by its own (`NT_187361.1`). The package's `contig_names.tsv`
holds every sequence of NCBI's GRCh37.p13 and GRCh38.p14 assembly reports, one row each, with these
names and the sequence's length, as the reports write them (`scripts/contig_names.py`,
docs/PROVENANCE.md). No name is derived from another by its pattern: a name the reports do not list
stays unmatched.

Two uses, kept apart on purpose:

- `PRIMARY_ACCESSIONS`: the assembled molecules' GenBank and RefSeq accessions (`CM000663.2`,
  `NC_000001.11`) → their canonical `chrom`. An accession with its version names one sequence, so
  `schema.CONTIG_ALIASES` maps them as it maps T2T-CHM13's.
- `scaffold_names(name, length)`: the other names of a scaffold a reference contig holds, used only
  where a contig name is matched against a reference store. A scaffold's `chrom` keeps the source's
  spelling (docs/SCHEMA.md, Decision 7); the alias is taken only where the reference contig's length is
  the report's for that name.
"""

from dataclasses import dataclass
from importlib import resources

COLUMNS = ("assembly", "role", "molecule", "genbank", "refseq", "ucsc", "length")
"""The table's columns: the assembly label, NCBI's `Sequence-Role`, `Assigned-Molecule`, `GenBank-Accn`,
`RefSeq-Accn` (only where the report says it is the same sequence), `UCSC-style-name` and
`Sequence-Length`. An empty cell is the report's `na`."""

ASSEMBLED_MOLECULE = "assembled-molecule"


@dataclass(frozen=True)
class NamedSequence:
    """One sequence of one assembly report."""

    assembly: str
    role: str
    molecule: str | None
    genbank: str
    refseq: str | None
    ucsc: str | None
    length: int

    def names(self) -> frozenset[str]:
        """Every name the report gives the sequence."""
        return frozenset(name for name in (self.genbank, self.refseq, self.ucsc) if name)


def _read() -> tuple[NamedSequence, ...]:
    lines = [line for line in resources.files("just_vcf").joinpath("contig_names.tsv").read_text().splitlines() if not line.startswith("#")]
    if tuple(lines[0].split("\t")) != COLUMNS:
        raise ValueError(f"contig_names.tsv: columns {lines[0]!r}, expected {COLUMNS}")
    found = []
    for line in lines[1:]:
        assembly, role, molecule, genbank, refseq, ucsc, length = line.split("\t")
        found.append(NamedSequence(assembly, role, molecule or None, genbank, refseq or None, ucsc or None, int(length)))
    return tuple(found)


SEQUENCES: tuple[NamedSequence, ...] = _read()
"""Every sequence of the GRCh37.p13 and GRCh38.p14 assembly reports, in the reports' order."""

PRIMARY_ACCESSIONS: dict[str, str] = {
    name: sequence.molecule
    for sequence in SEQUENCES
    if sequence.role == ASSEMBLED_MOLECULE and sequence.molecule
    for name in (sequence.genbank, sequence.refseq)
    if name
}
"""GRCh37's and GRCh38's assembled molecules: GenBank (`CM000663.1`, `CM000663.2`) and RefSeq
(`NC_000001.10`, `NC_000001.11`) accessions, versions verbatim → canonical `chrom` (`1`…`22`, `X`, `Y`,
`MT`). The UCSC-style names of these rows are not taken from here: `chr1` is canonical already, and the
reports' `chrM` for GRCh37 is not hg19's chrM (the Yoruba sequence, 16,571 bases)."""

_SCAFFOLDS_BY_NAME: dict[str, list[NamedSequence]] = {}
for _sequence in SEQUENCES:
    if _sequence.role != ASSEMBLED_MOLECULE:
        for _name in _sequence.names():
            _SCAFFOLDS_BY_NAME.setdefault(_name, []).append(_sequence)


def scaffold_names(name: str, length: int) -> frozenset[str]:
    """The other names of the scaffold a reference contig called `name`, of `length` bases, holds.

    Its GenBank and RefSeq accessions and its UCSC-style names (one per report that lists it: GRCh37's
    `chr4_gl000194_random` and GRCh38's `chr14_GL000194v1_random` are both `GL000194.1`), `name` left
    out. Empty when no report lists `name` as a scaffold, or when its length there is not `length`.
    """
    rows = [row for row in _SCAFFOLDS_BY_NAME.get(name, []) if row.length == length]
    return frozenset(other for row in rows for other in row.names()) - {name}
