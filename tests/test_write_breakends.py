"""F30: `write_vcf` names a breakend's mate inside ALT the way it names CHROM, and declares it.

Synthetic loaded facets carry breakends of each bracket form, with mates on the row's own contig, on
another contig with rows, on a contig only an ALT names, on the mitochondrion and on a scaffold, and
one source spelling its mates another way than its CHROM. Each is written under every `ContigNames`
choice and read back through `load`. The expected names are computed from the input records and
header, never read off the output: the mate contig is read canonically, then named as that style names
a CHROM; every contig a written ALT names must have a `##contig` line.
"""

import gzip
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ContigNames, load, write_vcf
from just_vcf.breakends import parse_breakend
from just_vcf.schema import canonical_chrom

GT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"

CHR_NAMES = {**{str(n): f"chr{n}" for n in range(1, 23)}, "X": "chrX", "Y": "chrY", "MT": "chrM"}
"""The `chr` style's listed names (MT is chrM); any other contig keeps its stored name."""


def record(chrom: str, pos: int, alt: str, ref: str = "A") -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1"


# A plain-named source (CHROM `1`, `2`) whose header declares 1, 2 and MT only. Its ALTs name mates on
# 1 and 2 (rows), on 3 (no row, not declared), on MT (declared, no row), on a scaffold (neither), and one
# mate spelled `chr2` beside CHROM `2`. A single breakend, a sequence ALT and a bracket ALT that is none
# of the four forms stay as written and declare nothing.
PLAIN_CONTIGS = ["1", "2", "MT"]
PLAIN = [
    record("1", 100, "A[2:500["),  # t[p[, another contig with rows
    record("1", 200, "A]1:900]"),  # t]p], the same contig
    record("1", 300, "]chr2:700]A"),  # ]p]t, mate spelled another way than CHROM
    record("1", 400, "[3:50[A"),  # [p[t, a contig only an ALT names
    record("1", 500, "A[MT:10["),  # MT: chrM under `chr`
    record("1", 600, "A[GL000192.1:5["),  # a scaffold: no `chr` name
    record("1", 700, "A."),  # single breakend
    record("1", 800, "G"),  # sequence ALT
    record("1", 900, "A[X:5]"),  # brackets that do not pair: no form, no mate
    record("2", 500, "A]1:100]"),
    record("2", 600, "A[1:10[,A]3:20]"),  # two breakend ALTs in one record
]

# A chr-named source: CHROM, header and every mate spelled `chr`.
CHR_CONTIGS = ["chr1", "chr2", "chr3", "chrM"]
CHR = [
    record("chr1", 100, "A[chr2:500["),
    record("chr1", 200, "A]chr1:900]"),
    record("chr1", 300, "]chr3:700]A"),
    record("chr1", 400, "[chrM:50[A"),
    record("chr2", 500, "A]chr1:100]"),
]

SOURCES = {"plain": (PLAIN_CONTIGS, PLAIN), "chr": (CHR_CONTIGS, CHR)}


def source_text(contigs: list[str], records: list[str]) -> str:
    return "\n".join(["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length=100000>" for c in contigs), GT, COLUMNS, *records]) + "\n"


def written_name(style: ContigNames, canonical: str, contigs: list[str], records: list[str]) -> str:
    """How `style` names a canonical contig, worked out from the source's own header and CHROMs."""
    if style is ContigNames.CHR:
        return CHR_NAMES.get(canonical, canonical)
    if style is ContigNames.SOURCE:
        spelled = {canonical_chrom(name): name for name in [*contigs, *(line.split("\t")[0] for line in records)]}
        return spelled.get(canonical, canonical)
    return canonical


def mates(alts: list[str]) -> list[tuple[str, str, str, int]]:
    """(form, bases, mate contig as written, mate position) of each of the four forms among `alts`."""
    parsed = [parse_breakend(alt) for alt in alts]
    return [(b.form, b.bases, b.mate_contig, b.mate_pos) for b in parsed if b is not None and b.mate_contig is not None]


def others(alts: list[str]) -> list[str]:
    """The ALTs that are none of the four forms: single breakends, sequence, unpaired brackets."""
    return [alt for alt in alts if not mates([alt])]


def header_contigs(path: Path) -> set[str]:
    text = gzip.decompress(path.read_bytes()).decode()
    return {line.split("ID=")[1].split(",")[0].rstrip(">") for line in text.splitlines() if line.startswith("##contig=")}


def written_rows(path: Path) -> list[list[str]]:
    text = gzip.decompress(path.read_bytes()).decode()
    return [line.split("\t") for line in text.splitlines() if not line.startswith("#")]


@pytest.mark.parametrize("source", sorted(SOURCES))
@pytest.mark.parametrize("style", list(ContigNames))
def test_a_mate_is_named_as_chrom_is_and_declared(tmp_path: Path, source: str, style: ContigNames) -> None:
    contigs, records = SOURCES[source]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(contigs, records))
    out = write_vcf(load(path), tmp_path / "out.vcf.gz", contig_names=style).path

    rows = written_rows(out)
    assert len(rows) == len(records)
    declared = header_contigs(out)
    for given, line in zip(records, rows, strict=True):
        given_fields = given.split("\t")
        assert line[0] == written_name(style, canonical_chrom(given_fields[0]), contigs, records)
        given_alts, alts = given_fields[4].split(","), line[4].split(",")
        # Only the mate's contig changes: the form, the bases and the position are as given, and an
        # ALT that is none of the four forms is written as it came.
        expected = [(form, bases, written_name(style, canonical_chrom(mate), contigs, records), pos) for form, bases, mate, pos in mates(given_alts)]
        assert mates(alts) == expected, given
        assert others(alts) == others(given_alts)
        assert {mate for _, _, mate, _ in mates(alts)} <= declared, "every contig an ALT names has a ##contig line"

    # The header is the source's contigs, the rows' and the mates', each named once by the style;
    # the unpaired `A[X:5]` names no contig.
    canonical = {canonical_chrom(c) for c in contigs}
    canonical |= {canonical_chrom(line.split("\t")[0]) for line in records}
    canonical |= {canonical_chrom(mate) for line in records for _, _, mate, _ in mates(line.split("\t")[4].split(","))}
    assert declared == {written_name(style, c, contigs, records) for c in canonical}


@pytest.mark.parametrize("source", sorted(SOURCES))
@pytest.mark.parametrize("style", list(ContigNames))
def test_load_write_load_keeps_every_breakend(tmp_path: Path, source: str, style: ContigNames) -> None:
    contigs, records = SOURCES[source]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(contigs, records))
    first = load(path)
    again = load(write_vcf(first, tmp_path / "out.vcf.gz", contig_names=style).path)

    def adjacencies(frame: pl.DataFrame) -> set[tuple]:
        """Every ALT read: a breakend as (form, bases, canonical mate, mate position), any other as written."""
        out = set()
        for chrom, pos, alts in frame.select("chrom", "pos", "alts").iter_rows():
            for alt in alts:
                b = parse_breakend(alt)
                read = (b.form, b.bases, canonical_chrom(b.mate_contig) if b.mate_contig else None, b.mate_pos) if b else alt
                out.add((chrom, pos, read))
        return out

    before, after = first.calls.collect(), again.calls.collect()
    assert after.height == before.height
    assert adjacencies(after) == adjacencies(before)
    if style is ContigNames.SOURCE and source == "chr":
        # A source that spells every contig one way gets every ALT back verbatim.
        assert after.select("chrom", "pos", "alts").rows() == before.select("chrom", "pos", "alts").rows()


def test_one_contig_under_two_canonical_names_is_declared_once(tmp_path: Path) -> None:
    # A GRCh38 lift through UCSC's chain puts a row on the chain's `chrUn_KI270742v1` and names its mate
    # by the target FASTA's `KI270742.1`; the lifted facet's `contig_map` says they are one contig and
    # its header states none (hu1C1368's Manta file, lifted). Built here from a load, the facet made as
    # the lift makes it.
    scaffold, accession = "chrUn_KI270742v1", "KI270742.1"
    records = [record(scaffold, 100, f"A[{accession}:500["), record(scaffold, 500, f"]{accession}:100]A")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text([scaffold], records))
    loaded = load(path)
    lifted_like = loaded.model_copy(
        update={
            "facets": loaded.facets.with_columns(
                pl.lit([{"source": accession, "chrom": scaffold}]).alias("contig_map"),
                pl.lit(None, dtype=loaded.facets.collect_schema()["header_contigs"]).alias("header_contigs"),
            )
        }
    )
    for style, chroms, header in (
        (ContigNames.SOURCE, {accession}, [accession]),
        (ContigNames.STORED, {scaffold}, [accession, scaffold]),
    ):
        out = write_vcf(lifted_like, tmp_path / f"out-{style.value}.vcf.gz", contig_names=style).path
        text = gzip.decompress(out.read_bytes()).decode()
        assert sorted(line.split("ID=")[1].rstrip(">") for line in text.splitlines() if line.startswith("##contig=")) == header
        rows = written_rows(out)
        assert {row[0] for row in rows} == chroms
        assert {mate for row in rows for _, _, mate, _ in mates(row[4].split(","))} == {accession}
