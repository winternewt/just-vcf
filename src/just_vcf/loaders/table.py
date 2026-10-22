"""Load the vendor genotype tables that are not 23andMe-shaped into `calls` rows.

Covers AncestryDNA, FamilyTreeDNA (both comma and space-separated exports), the vendor-neutral
`RSID,CHROMOSOME,POSITION,RESULT` table (MyHeritage and Vitagene exports land there), Color
`discovery_genotypes.csv`, deCODEme, Genos, and Illumina GenomeStudio Final Reports. Every format is read
through its `FormatSpec`: where the column names are, their alias spellings, the delimiter, the
coordinate convention and the no-call spellings. Columns are found by name, never by position, except
where the format names none (Genos without a `##Columns=` line).

What every table here shares with the 23andMe layout (SCHEMA.md Decision 2):

- no REF: the source writes none, and REF only ever comes from the assembly sequence (charter 10);
- `gt_alleles` holds what the file wrote, one entry per written allele, null for a missing one, so a
  half-call keeps its called allele and its position in the pair;
- `source_gt` keeps the genotype text. Where a format writes one allele per column (AncestryDNA), the
  two fields are concatenated without a separator (`A` + `G` gives `AG`, `0` + `0` gives `00`);
- an absent site is unknown, never reference (charter 11).

Rows that do not fit are dropped and counted, never guessed at; a genotype that cannot be split is kept
with `gt_alleles` null and its text in `source_gt`. The counts come back as warnings, one per reason.
"""

import logging
import re
from dataclasses import dataclass, field

import polars as pl

from just_vcf.formats import get_spec
from just_vcf.formats.spec import Coordinates, FormatSpec, HeaderMode
from just_vcf.loaders.text import (
    ALLELE_LETTERS,
    chrom_columns,
    ids_expr,
    lines,
    one_field_alleles,
    single_position_end,
)
from just_vcf.loaders.vcf import MultiSampleError
from just_vcf.schema import EvidenceBasis

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset(
    {
        "ancestrydna",
        "ftdna",
        "ftdna_space",
        "rsid_result_csv",
        "color_genotypes",
        "decodeme",
        "genos",
        "illumina_final_report",
    }
)
"""Format ids this loader reads."""

UNPLACED_CHROMS = frozenset({"0"})
"""Chromosome codes of probes the manifest leaves unplaced (Illumina, GSA exports). They have no
position, so they cannot be `calls` rows."""

_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

_NO_CALL_WIDTH: dict[str, int] = {"ftdna": 2, "ftdna_space": 2, "rsid_result_csv": 2, "illumina_final_report": 2}
"""Formats whose no-call spelling is a token in a field that holds two letters on every called row
(`---`, a Final Report's `-`): it stands for two missing alleles. See `one_field_alleles`."""

_REQUIRED: dict[str, tuple[str, ...]] = {
    "ancestrydna": ("rsid", "chromosome", "position", "allele1", "allele2"),
    "ftdna": ("rsid", "chromosome", "position", "genotype"),
    "ftdna_space": ("rsid", "chromosome", "position", "genotype"),
    "rsid_result_csv": ("rsid", "chromosome", "position", "genotype"),
    "color_genotypes": ("rsid", "chromosome", "position", "genotype"),
    "decodeme": ("rsid", "chromosome", "position", "strand", "genotype"),
    "genos": ("chromosome", "begin", "end", "rsid", "genotype"),
    "illumina_final_report": ("snp_name", "chromosome", "position", "genotype_plus"),
}
"""The canonical columns each format needs, by `ColumnSpec` name."""

_CONSUMED: dict[str, tuple[str, ...]] = {
    "color_genotypes": ("sample_id", "reference_version"),
    "illumina_final_report": ("sample_id",),
}
"""Columns read beyond the required ones; every other column of the file is listed as dropped."""


class TableLayoutError(ValueError):
    """The file's columns do not carry what its format needs (a position, a plus-strand genotype)."""


@dataclass
class Table:
    """What `read_table` makes of one payload."""

    calls: pl.DataFrame
    """`calls` rows without `facet_id`, plus `source_chrom` for the facet's `contig_map`."""
    warnings: list[str] = field(default_factory=list)
    dropped_columns: list[str] | None = None
    """The file's column names that no frame column carries, in file order; None where the file
    names no columns."""
    sample: str | None = None
    """The sample id the file states, if any. Identifying: never logged."""
    sample_index: int | None = None
    """The chosen sample's index among the file's sample ids, for the facet id."""
    rows_state_build: bool = False
    """Whether the rows themselves state the build (Color's per-row ReferenceVersion). Then `build` is
    the answer, None included, and no header hint overrides it."""
    build: str | None = None
    """The build the rows state, when they agree on one."""


def read_table(data: bytes, spec: FormatSpec, sample: str | None = None) -> Table:
    """Parse a decompressed payload of `spec`'s format into `calls` rows.

    `sample` picks one sample of a table with a sample column; more than one sample and no choice is
    refused with `MultiSampleError`.
    """
    if spec.id not in FORMATS:
        raise ValueError(f"{spec.id} is not a vendor table this loader reads")
    names, rows, meta = _split(data, spec)
    table = Table(calls=pl.DataFrame())
    columns = _columns(names, spec)
    missing = [name for name in _REQUIRED[spec.id] if name not in columns]
    if missing:
        raise TableLayoutError(f"{spec.id}: no column for {', '.join(missing)}{_strand_hint(spec, columns)}")
    if names is not None:
        used = {*_REQUIRED[spec.id], *_CONSUMED.get(spec.id, ())}
        table.dropped_columns = [name for name, canonical in zip(names, _canonical_names(names, spec)) if canonical not in used]

    width = len(names) if names is not None else len(spec.columns)
    if names is not None:
        repeated = pl.col("fields").list.join("\x1f").str.to_lowercase() == "\x1f".join(names).lower()
        again = rows.filter(repeated).height
        if again:
            table.warnings.append(f"{again} repeated column line(s) inside the table skipped (the file joins several exports)")
        rows = rows.filter(~repeated)
    if spec.id == "color_genotypes":
        rows = _cut_trailer(rows, columns["position"])
    wrong_width = rows.filter(pl.col("fields").list.len() != width)
    if wrong_width.height:
        widths = sorted(set(wrong_width["fields"].list.len().to_list()))
        table.warnings.append(
            f"{wrong_width.height} row(s) without exactly {width} fields dropped "
            f"(widths seen: {', '.join(map(str, widths))})"
        )
    frame = rows.filter(pl.col("fields").list.len() == width).select(
        pl.col("fields").list.get(index).alias(name) for name, index in columns.items()
    )

    frame = _choose_sample(frame, spec, sample, meta, table)

    unplaced = frame.filter(pl.col("chromosome").is_in(list(UNPLACED_CHROMS)))
    if unplaced.height:
        table.warnings.append(f"{unplaced.height} row(s) on chromosome 0 (probes with no mapped position) dropped")
    frame = frame.filter(~pl.col("chromosome").is_in(list(UNPLACED_CHROMS)))

    if spec.id == "genos":
        frame = _genos(frame, table)
    else:
        frame = _single_site(frame, spec, table)

    unreadable = frame.filter(pl.col("gt_alleles").is_null())
    if "strand" in frame.columns:
        unreadable = unreadable.filter(pl.col("strand").is_in(["+", "-"]))  # counted by `_plus_strand`
    if unreadable.height:
        lengths = sorted(set(unreadable["source_gt"].str.len_chars().to_list()))
        table.warnings.append(
            f"{unreadable.height} genotype(s) not readable as one allele per haplotype kept with gt_alleles "
            f"null and the text in source_gt (lengths seen: {', '.join(map(str, lengths))})"
        )
    if spec.id == "color_genotypes":
        table.rows_state_build = True
        table.build = _stated_build(frame, table)
    table.calls = frame.select(
        "chrom",
        "pos",
        "end_pos",
        "ids",
        "gt_alleles",
        "source_gt",
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
        "source_chrom",
    )
    return table


def samples(data: bytes, spec: FormatSpec) -> list[str]:
    """The sample ids in a table's sample column, in the order they first appear; `[]` for a format
    without one. They can be real names: treat them as identifying."""
    if spec.id not in FORMATS:
        raise ValueError(f"{spec.id} is not a vendor table this loader reads")
    names, rows, _ = _split(data, spec)
    index = _columns(names, spec).get("sample_id")
    if index is None:
        return []
    values = rows.select(pl.col("fields").list.get(index, null_on_oob=True).alias("sample_id"))
    if names is not None:
        values = values.filter(pl.col("sample_id").str.to_lowercase() != names[index].lower())
    return [value for value in values["sample_id"].unique(maintain_order=True).to_list() if value]


def spec_for(format_id: str) -> FormatSpec:
    """The registered spec for a format this loader reads."""
    if format_id not in FORMATS:
        raise ValueError(f"{format_id} is not a vendor table this loader reads")
    return get_spec(format_id)


# --- layout ----------------------------------------------------------------------------------------


def _split(data: bytes, spec: FormatSpec) -> tuple[list[str] | None, pl.DataFrame, list[str]]:
    """The file's column names (None where it names none), its body rows split into `fields`, and
    its comment lines.

    Fields are stripped of whitespace and of the double quotes FTDNA and MyHeritage put around every
    value.
    """
    text = lines(data).with_row_index("n")
    blank = pl.col("line").str.strip_chars() == ""
    comment = pl.col("line").str.starts_with(spec.comment_prefix) if spec.comment_prefix else pl.lit(False)
    meta = text.filter(comment)["line"].to_list()
    if spec.header is HeaderMode.SECTION:
        marker = rf"^\[{re.escape(spec.section or '')}\][\t,]*$"
        starts = text.filter(pl.col("line").str.contains(marker))["n"]
        if starts.is_empty():
            raise TableLayoutError(f"{spec.id}: no [{spec.section}] section")
        text = text.filter(pl.col("n") > starts[0])
    text = text.filter(~blank & ~comment)
    names: list[str] | None = None
    if spec.header in (HeaderMode.LINE, HeaderMode.SECTION):
        if text.is_empty():
            raise TableLayoutError(f"{spec.id}: no column line")
        header = text["line"][0]
        text = text.slice(1)
        delimiter = spec.delimiter if spec.delimiter not in (None, " ") else ("\t" if "\t" in header else ",")
        names = [name.strip().strip('"') for name in header.split(delimiter)]
    else:
        delimiter = spec.delimiter or "\t"
        names = _meta_columns(meta, delimiter)

    line = pl.col("line")
    if spec.delimiter == " ":
        # FTDNA's space-separated export turns into quoted comma rows near its end: each row is read by
        # its own shape.
        line = pl.when(line.str.contains(",")).then(line).otherwise(line.str.replace_all(" ", ","))
        delimiter = ","
    rows = text.select(line.str.split(delimiter).list.eval(pl.element().str.strip_chars(' \t"')).alias("fields"))
    return names, rows, meta


def _meta_columns(meta: list[str], delimiter: str) -> list[str] | None:
    """Column names stated in a `##Columns=` meta line (Genos), or None."""
    for line in meta:
        found = re.match(r"^##Columns=(.*)$", line)
        if found:
            return [name.strip() for name in found.group(1).split(delimiter)]
    return None


def _canonical_names(names: list[str], spec: FormatSpec) -> list[str | None]:
    """Each file column's canonical `ColumnSpec` name, None for a column the spec does not declare."""
    spellings: dict[str, str] = {}
    for column in spec.columns:
        for spelling in (column["name"], *column.get("aliases", ())):
            spellings[spelling.lower()] = column["name"]
    return [spellings.get(name.lower()) for name in names]


def _columns(names: list[str] | None, spec: FormatSpec) -> dict[str, int]:
    """Canonical column name → field index. A headerless format's columns come in spec order."""
    if names is None:
        return {column["name"]: index for index, column in enumerate(spec.columns)}
    found: dict[str, int] = {}
    for index, canonical in enumerate(_canonical_names(names, spec)):
        if canonical is not None and canonical not in found:
            found[canonical] = index
    return found


def _strand_hint(spec: FormatSpec, columns: dict[str, int]) -> str:
    if spec.id == "illumina_final_report" and "genotype_plus" not in columns:
        named = sorted(name for name in columns if name.startswith(("allele1_", "allele2_")))
        if named:
            return (
                f"; the alleles are given only as {', '.join(named)}, which are not the genomic plus strand, "
                "and turning them into plus-strand letters needs the chip manifest"
            )
    return ""


def _cut_trailer(rows: pl.DataFrame, position: int) -> pl.DataFrame:
    """Color's table ends at the first row whose Position is empty; a disclaimer follows it."""
    empty = rows.with_row_index("n").filter(pl.col("fields").list.get(position, null_on_oob=True).fill_null("") == "")
    return rows if empty.is_empty() else rows.slice(0, empty["n"][0])


def _choose_sample(frame: pl.DataFrame, spec: FormatSpec, sample: str | None, meta: list[str], table: Table) -> pl.DataFrame:
    """Keep the chosen sample's rows of a table with a sample column (Decision 10)."""
    if "sample_id" not in frame.columns:
        if spec.id == "genos":
            stated = next((line.split("=", 1)[1].strip() for line in meta if line.startswith("##Sample=")), None)
            table.sample = stated or None
        return frame
    samples = frame["sample_id"].unique(maintain_order=True).to_list()
    if len(samples) > 1 and sample is None:
        raise MultiSampleError(
            f"{spec.id}: {len(samples)} samples in the sample column; pass sample= to load one (just_vcf.table_samples lists them)"
        )
    if sample is not None and sample not in samples:
        raise MultiSampleError(f"{spec.id}: the requested sample is not in the sample column")
    chosen = sample if sample is not None else (samples[0] if samples else None)
    table.sample = chosen or None
    table.sample_index = samples.index(chosen) if chosen is not None else None
    return frame.filter(pl.col("sample_id") == chosen).drop("sample_id")


def _stated_build(frame: pl.DataFrame, table: Table) -> str | None:
    """The build every row states in ReferenceVersion, or None (with a warning) when rows disagree."""
    stated = sorted(set(frame["reference_version"].to_list()) - {""})
    if len(stated) > 1:
        table.warnings.append(f"rows state {len(stated)} different builds ({', '.join(stated)}); assembly left unknown")
        return None
    return stated[0] if stated and re.fullmatch(r"GRCh3[678]", stated[0]) else None


# --- genotypes -------------------------------------------------------------------------------------


def _positions(frame: pl.DataFrame, column: str, offset: int, table: Table) -> pl.DataFrame:
    """`column` as a 1-based `u32` position (`offset` added); rows that are not a count are dropped."""
    frame = frame.with_columns((pl.col(column).cast(pl.Int64, strict=False) + offset).alias("_pos"))
    fits = pl.col("_pos").is_between(0, 2**32 - 1).fill_null(False)
    bad = frame.filter(~fits).height
    if bad:
        table.warnings.append(f"{bad} row(s) with a {column} that is not a 32-bit count dropped")
    return frame.filter(fits)


def _single_site(frame: pl.DataFrame, spec: FormatSpec, table: Table) -> pl.DataFrame:
    """Every format whose row is one position with a one- or two-field genotype."""
    offset = 1 if spec.coordinates is Coordinates.ZERO_BASED_HALF_OPEN else 0
    frame = _positions(frame, "position", offset, table).with_columns(pl.col("_pos").cast(pl.UInt32).alias("pos"))
    if spec.id == "ancestrydna":
        frame = frame.with_columns(_two_field_alleles("allele1", "allele2", spec.no_call).alias("gt_alleles"),
                                   (pl.col("allele1") + pl.col("allele2")).alias("source_gt"))
    else:
        genotype = "genotype_plus" if spec.id == "illumina_final_report" else "genotype"
        alleles = one_field_alleles(genotype, spec.no_call, _NO_CALL_WIDTH.get(spec.id))
        frame = frame.with_columns(alleles.alias("gt_alleles"), pl.col(genotype).alias("source_gt"))
    if spec.id == "decodeme":
        frame = _plus_strand(frame, table)
    rsid = "snp_name" if spec.id == "illumina_final_report" else "rsid"
    return frame.with_columns(
        *chrom_columns("chromosome"),
        single_position_end("gt_alleles", "pos").alias("end_pos"),
        ids_expr(rsid).alias("ids"),
    )


def _two_field_alleles(first: str, second: str, no_call: tuple[str, ...]) -> pl.Expr:
    """`gt_alleles` from one allele per column: always two entries, a no-call spelling null.

    AncestryDNA writes a haploid call (male X and Y, MT) as the same letter twice. That is the vendor's
    convention, and it stays two entries: the list length is what the file wrote, not a ploidy.
    """
    def allele(column: str) -> pl.Expr:
        return pl.when(pl.col(column).is_in(list(no_call))).then(pl.lit(None, dtype=pl.String)).otherwise(pl.col(column))

    letter = rf"^[{ALLELE_LETTERS}]$"
    readable = pl.all_horizontal(
        (pl.col(name).str.contains(letter) | pl.col(name).is_in(list(no_call))) for name in (first, second)
    )
    return pl.when(readable).then(pl.concat_list(allele(first), allele(second))).otherwise(None)


def _plus_strand(frame: pl.DataFrame, table: Table) -> pl.DataFrame:
    """deCODEme writes each genotype on the strand in its Strand column: complement the `-` rows.

    `source_gt` keeps the letters as written. I/D codes have no complement and stay as they are. A row
    whose strand is neither `+` nor `-` cannot be put on the plus strand, so its `gt_alleles` is null.
    """
    unknown = frame.filter(~pl.col("strand").is_in(["+", "-"]))
    if unknown.height:
        table.warnings.append(f"{unknown.height} row(s) with a strand other than + or - kept with gt_alleles null")
    complemented = pl.col("gt_alleles").list.eval(pl.element().replace(_COMPLEMENT))
    return frame.with_columns(
        pl.when(pl.col("strand") == "+")
        .then(pl.col("gt_alleles"))
        .when(pl.col("strand") == "-")
        .then(complemented)
        .otherwise(None)
        .alias("gt_alleles")
    )


_GENOS_ALLELE = r"^([ACGTN]+|-)$"


def _genos(frame: pl.DataFrame, table: Table) -> pl.DataFrame:
    """Genos rows: a 0-based half-open span and alleles joined by `/`.

    A span of one or more bases is 1-based `pos = begin + 1`, `end_pos = end`. An insertion site has
    `begin == end`: it sits between two bases, and VCF places it at the base before (1-based `begin`),
    which is where `pos` goes. Its REF span is that anchor base, unknown until the reference is read, so
    `end_pos` is null. Alleles stay as written, `-` the empty allele, until RM3 reads the anchor.
    """
    frame = _positions(frame, "begin", 0, table).rename({"_pos": "_begin"})
    frame = _positions(frame, "end", 0, table).rename({"_pos": "_end"})
    backwards = frame.filter(pl.col("_end") < pl.col("_begin"))
    if backwards.height:
        table.warnings.append(f"{backwards.height} row(s) whose end lies before their begin dropped")
    frame = frame.filter(pl.col("_end") >= pl.col("_begin"))
    insertion = pl.col("_end") == pl.col("_begin")
    alleles = pl.col("genotype").str.split("/")
    readable = alleles.list.eval(pl.element().str.contains(_GENOS_ALLELE)).list.all() & (pl.col("genotype") != "")
    frame = frame.with_columns(
        pl.when(insertion).then(pl.col("_begin")).otherwise(pl.col("_begin") + 1).cast(pl.UInt32).alias("pos"),
        pl.when(insertion).then(None).otherwise(pl.col("_end")).cast(pl.UInt32).alias("end_pos"),
        pl.when(readable).then(alleles).otherwise(None).alias("gt_alleles"),
        pl.col("genotype").alias("source_gt"),
        *chrom_columns("chromosome"),
        ids_expr("rsid").alias("ids"),
    )
    sites = frame.filter(insertion).height
    if sites:
        table.warnings.append(
            f"{sites} insertion site(s) (begin == end) placed at the base before them with end_pos null; their "
            "alleles stay as written (`-` is the empty allele) until the anchor base is read from the reference"
        )
    return frame
