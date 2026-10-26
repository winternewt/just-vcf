"""F30: `write_vcf` names a breakend's mate inside ALT the way it names CHROM, and declares it.

Synthetic loaded facets carry breakends of each bracket form, with mates on the row's own contig, on
another contig with rows, on a contig only an ALT names, on the mitochondrion and on a scaffold, and
one source spelling its mates another way than its CHROM. Each is written under every `ContigNames`
choice and read back through `load`. The expected names are computed from the input records and
header, never read off the output: the mate contig is read canonically, then named as that style names
a CHROM, except that `source` keeps an ALT's own spelling of a contig the file names nowhere else;
every contig a written ALT names must have a `##contig` line, except an `##assembly` contig (`<ctg1>`).
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

# The same chr-named records with a header that declares only the contigs with rows (chr3 and chrM are
# named by ALTs alone), and with no header contigs at all (every mate contig is).
SOURCES = {
    "plain": (PLAIN_CONTIGS, PLAIN),
    "chr": (CHR_CONTIGS, CHR),
    "chr_rows_declared": (["chr1", "chr2"], CHR),
    "chr_undeclared": ([], CHR),
}


def source_text(contigs: list[str], records: list[str]) -> str:
    return "\n".join(["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length=100000>" for c in contigs), GT, COLUMNS, *records]) + "\n"


def written_name(style: ContigNames, spelled: str, contigs: list[str], records: list[str]) -> str:
    """How `style` names a contig the source spelled `spelled` (in CHROM, the header or an ALT), worked
    out from the source's own header and CHROMs: `source` takes their spelling of the contig where
    they have one, else the spelling given."""
    canonical = canonical_chrom(spelled)
    if style is ContigNames.CHR:
        return CHR_NAMES.get(canonical, canonical)
    if style is ContigNames.SOURCE:
        known = {canonical_chrom(name): name for name in [*contigs, *(line.split("\t")[0] for line in records)]}
        return known.get(canonical, spelled)
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
        assert line[0] == written_name(style, given_fields[0], contigs, records)
        given_alts, alts = given_fields[4].split(","), line[4].split(",")
        # Only the mate's contig changes: the form, the bases and the position are as given, and an
        # ALT that is none of the four forms is written as it came.
        expected = [(form, bases, written_name(style, mate, contigs, records), pos) for form, bases, mate, pos in mates(given_alts)]
        assert mates(alts) == expected, given
        assert others(alts) == others(given_alts)
        assert {mate for _, _, mate, _ in mates(alts)} <= declared, "every contig an ALT names has a ##contig line"

    # The header is the source's contigs, the rows' and the mates', each named once by the style;
    # the unpaired `A[X:5]` names no contig.
    spelled = [*contigs, *(line.split("\t")[0] for line in records)]
    spelled += [mate for line in records for _, _, mate, _ in mates(line.split("\t")[4].split(","))]
    assert declared == {written_name(style, name, contigs, records) for name in spelled}


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
    if style is ContigNames.SOURCE and source.startswith("chr"):
        # A source that spells every contig one way gets every ALT back verbatim, whether or not its
        # header declares the contigs only an ALT names.
        assert after.select("chrom", "pos", "alts").rows() == before.select("chrom", "pos", "alts").rows()


@pytest.mark.parametrize("contigs", [[], ["chr1", "chr2"]], ids=["no_header_contigs", "rows_declared"])
def test_source_writes_a_mate_the_file_alone_names_as_it_spelled_it(tmp_path: Path, contigs: list[str]) -> None:
    # A chr-named source whose mates lie on contigs no row stands on: chr2 (declared in one header, not
    # in the other) and chr3 (declared in neither). Under `source` each comes back as the ALT spelled
    # it, and is declared under that spelling.
    records = [record("chr1", 100, "A[chr2:500["), record("chr1", 200, "A[chr3:500[")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(contigs, records))
    out = write_vcf(load(path), tmp_path / "out.vcf.gz", contig_names=ContigNames.SOURCE).path

    rows = written_rows(out)
    assert [row[:5] for row in rows] == [line.split("\t")[:5] for line in records]
    given = {line.split("\t")[0] for line in records} | {mate for line in records for _, _, mate, _ in mates([line.split("\t")[4]])}
    assert header_contigs(out) == given | set(contigs)


def test_source_prefers_the_chrom_spelling_of_a_mate_with_rows(tmp_path: Path) -> None:
    # CHROM spells contig 2 `2`; one ALT spells it `chr2`. The file is kept consistent with its own CHROM
    # column: the mate is written `2`, and contig 2 gets one `##contig` line.
    records = [record("1", 100, "]chr2:700]A"), record("2", 700, "A[1:100[")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text([], records))
    out = write_vcf(load(path), tmp_path / "out.vcf.gz", contig_names=ContigNames.SOURCE).path
    assert [row[4] for row in written_rows(out)] == ["]2:700]A", "A[1:100["]
    assert header_contigs(out) == {"1", "2"}


def test_source_refuses_a_contig_only_alts_name_spelled_two_ways(tmp_path: Path) -> None:
    # Contig 3 has no row and no header line, and the ALTs spell it `chr3` and `3`: which spelling the
    # source meant cannot be restored, as for a CHROM spelled two ways. The other styles write it.
    records = [record("1", 100, "A[chr3:500["), record("1", 200, "A[3:600[")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(["1"], records))
    loaded = load(path)
    with pytest.raises(ValueError, match="spelled 3 more than one way"):
        write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=ContigNames.SOURCE)
    out = write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=ContigNames.CHR).path
    assert [row[4] for row in written_rows(out)] == ["A[chr3:500[", "A[chr3:600["]


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


# VCF 4.3 §5.4.2's large insertion: the mate is a contig of the `##assembly` file, in angle brackets
# (§1.6.1 reads an `<ID>` CHROM the same way), not a contig of the reference.
ASSEMBLY = ["##fileformat=VCFv4.2", "##contig=<ID=13,length=400000>", "##assembly=file:///assembly.fa", GT, COLUMNS]
ASSEMBLY_RECORDS = [record("13", 321682, "]<ctg1>:7]T", ref="T"), record("13", 321683, "A[<ctg1>:214[")]


@pytest.mark.parametrize("style", list(ContigNames))
def test_an_assembly_mate_is_written_verbatim_and_declared_nowhere(tmp_path: Path, style: ContigNames) -> None:
    path = tmp_path / "in.vcf"
    path.write_text("\n".join([*ASSEMBLY, *ASSEMBLY_RECORDS]) + "\n")
    out = write_vcf(load(path), tmp_path / "out.vcf.gz", contig_names=style).path

    given = [line.split("\t") for line in ASSEMBLY_RECORDS]
    assert [row[4] for row in written_rows(out)] == [fields[4] for fields in given]
    assert header_contigs(out) == {written_name(style, fields[0], ["13"], ASSEMBLY_RECORDS) for fields in given}


@pytest.mark.parametrize("style", list(ContigNames))
def test_htslib_reads_an_assembly_mates_header_without_a_warning(tmp_path: Path, style: ContigNames, capfd: pytest.CaptureFixture) -> None:
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "in.vcf"
    path.write_text("\n".join([*ASSEMBLY, *ASSEMBLY_RECORDS]) + "\n")
    out = write_vcf(load(path), tmp_path / "out.vcf.gz", contig_names=style).path
    capfd.readouterr()
    with pysam.VariantFile(str(out)) as vcf:
        alts = [row.alts for row in vcf]
    # htslib writes its warnings to the C-level stderr, which `capfd` reads and `capsys` would not.
    err = capfd.readouterr().err
    assert alts == [(line.split("\t")[4],) for line in ASSEMBLY_RECORDS]
    assert "[W::" not in err and "[E::" not in err, err


def test_a_skipped_rows_mate_is_not_declared(tmp_path: Path) -> None:
    # The second record's GT names an ALT it does not have, so `resolved_only` leaves it out; its mate's
    # contig (7: no row, not in the source header) must not reach the written header.
    records = [record("1", 100, "A[2:500["), record("1", 200, "A[7:500[").replace("0/1", "0/2")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(["1"], records))
    written = write_vcf(load(path), tmp_path / "out.vcf.gz", resolved_only=True)
    assert written.rows_skipped == {"genotype_index_past_alts": 1}

    rows = written_rows(written.path)
    assert [row[1] for row in rows] == ["100"]
    named = {mate for row in rows for _, _, mate, _ in mates(row[4].split(","))}
    assert header_contigs(written.path) == {"1"} | named
