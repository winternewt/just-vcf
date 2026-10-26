"""Load the delimited tables that state a reference base at each site into `calls` rows (CASAVA SNP and
indel tables, SOAPsnp, Big Y, Color's variant set, GET-Evidence GFFs).

These are not array exports: the file writes the site's REF beside the genotype, as a VCF does, so a
row gets `source_ref`, `alts` and a `gt` that indexes `[source_ref] + alts` at load time, and
`resolve_ref` later compares `source_ref` with the assembly (`ref_check`). REF itself still only comes
from the sequence (charter 10).

What every format here shares (SCHEMA.md Decision 2):

- `gt_alleles` holds the alleles the file names, one entry per allele written: a one-letter call is
  one entry, a two-letter call two. The list length is what the file wrote, not an asserted ploidy;
  a consumer that needs the ploidy reads the karyotype work, not this list;
- `alts` are the called alleles other than `source_ref`, sorted by allele string: none of these
  formats states an ALT order beyond the call itself;
- `source_gt` keeps the genotype text as the file wrote it;
- rows that do not fit are dropped or kept unread and counted, one warning per reason.

Formats:

- `casava_snp_table`: CASAVA-era SNP calls, one chromosome per file. `ad` comes from the four base
  counts (they sum to `Bases Used` on every one of the 3,389,516 PGP-HMS rows), `dp` is `Total Bases`.
  The `Sample Name` column is free text, in real files a person's name: it is never read into a frame,
  logged or kept, and appears only by its column name in `dropped_fields`.
"""

import logging
import re
from dataclasses import dataclass, field

import polars as pl

from just_vcf.formats import get_spec
from just_vcf.formats.spec import FormatSpec
from just_vcf.loaders.table import TableLayoutError, _canonical_names, _columns, _split
from just_vcf.loaders.text import chrom_columns, ids_expr
from just_vcf.schema import BLOCKS, CALLS, EvidenceBasis, RefCheck, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset({"casava_snp_table", "casava_indel_table", "soapsnp_cns", "ftdna_bigy", "color_variant_set", "get_evidence_gff"})
"""Format ids this loader reads."""

_BASES = ("A", "C", "G", "T")


@dataclass
class SiteTable:
    """What `read_site_table` makes of one payload."""

    calls: pl.DataFrame
    """`calls` rows without `facet_id`, plus `source_chrom` for the facet's `contig_map`."""
    blocks: pl.DataFrame | None = None
    """`blocks` rows without `facet_id`, for a format that states reference sites; None for one that
    states none."""
    warnings: list[str] = field(default_factory=list)
    dropped_columns: list[str] = field(default_factory=list)
    """The file's fields that no frame column carries, by the file's own names."""
    read_evidence: list[str] = field(default_factory=list)
    """The VCF FORMAT keys this format's rows carry the equivalent of (`AD`, `DP`, `GQ`), for RM4's
    read-evidence gate: a callset type is inferred only for calls made from reads."""
    sample: str | None = None
    """The sample id the file states, if any. Identifying: never logged."""
    build: str | None = None
    """A build the rows state (Color's per-row ReferenceVersion), when they agree on one."""


def read_site_table(data: bytes, spec: FormatSpec) -> SiteTable:
    """Parse a decompressed payload of `spec`'s format into `calls` rows."""
    if spec.id not in FORMATS:
        raise ValueError(f"{spec.id} is not a site table this loader reads")
    return _READERS[spec.id](data, spec)


def spec_for(format_id: str) -> FormatSpec:
    """The registered spec for a format this loader reads."""
    if format_id not in FORMATS:
        raise ValueError(f"{format_id} is not a site table this loader reads")
    return get_spec(format_id)


def to_calls_frame(calls: pl.DataFrame, facet_id: str) -> pl.LazyFrame:
    """`read_site_table` rows as a `calls` frame for `facet_id`."""
    return conform(calls.drop("source_chrom").lazy().with_columns(pl.lit(facet_id).alias("facet_id")), CALLS)


def to_blocks_frame(blocks: pl.DataFrame, facet_id: str) -> pl.LazyFrame:
    """`SiteTable.blocks` rows as a `blocks` frame for `facet_id`."""
    return conform(blocks.lazy().with_columns(pl.lit(facet_id).alias("facet_id")), BLOCKS)


# --- shared ----------------------------------------------------------------------------------------


def _named_rows(
    data: bytes, spec: FormatSpec, required: tuple[str, ...], table: SiteTable, ends_at_empty: str | None = None
) -> tuple[list[str], pl.DataFrame]:
    """The file's column names, and its body rows with one column per canonical name the file has.

    `ends_at_empty` names a column whose first empty cell (or a `Raw data disclaimer` cell) ends the
    table: what follows is a footer, not rows. Rows of another width are dropped and counted. A file
    whose column line lacks a required column is refused with `TableLayoutError`.
    """
    names, rows, _ = _split(data, spec)
    if names is None:
        raise TableLayoutError(f"{spec.id}: no column line")
    columns = _columns(names, spec)
    missing = [name for name in required if name not in columns]
    if missing:
        raise TableLayoutError(f"{spec.id}: no column for {', '.join(missing)}")
    if ends_at_empty is not None:
        cell = pl.col("fields").list.get(columns[ends_at_empty], null_on_oob=True).fill_null("")
        end = rows.with_row_index("_n").filter((cell == "") | cell.str.starts_with("Raw data disclaimer"))
        if not end.is_empty():
            rows = rows.slice(0, end["_n"][0])
    wrong = rows.filter(pl.col("fields").list.len() != len(names))
    if wrong.height:
        widths = sorted(set(wrong["fields"].list.len().to_list()))
        table.warnings.append(f"{wrong.height} row(s) without exactly {len(names)} fields dropped (widths seen: {', '.join(map(str, widths))})")
    return names, rows.filter(pl.col("fields").list.len() == len(names)).select(
        pl.col("fields").list.get(index).alias(name) for name, index in columns.items()
    )


def _dropped(names: list[str], spec: FormatSpec, used: set[str]) -> list[str]:
    """The file's column names whose canonical name is not among `used`, in file order."""
    return [name for name, canonical in zip(names, _canonical_names(names, spec), strict=True) if canonical not in used]


def _positions(frame: pl.DataFrame, column: str, table: SiteTable) -> pl.DataFrame:
    """`column` as the 1-based `u32` `pos`; rows whose value is not a count above 0 are dropped."""
    frame = frame.with_columns(pl.col(column).cast(pl.Int64, strict=False).alias("_pos"))
    fits = pl.col("_pos").is_between(1, 2**32 - 1).fill_null(False)
    bad = frame.filter(~fits).height
    if bad:
        table.warnings.append(f"{bad} row(s) with a {column} that is not a position (a count from 1) dropped")
    return frame.filter(fits).with_columns(pl.col("_pos").cast(pl.UInt32).alias("pos")).drop("_pos")


def _site_genotype(frame: pl.DataFrame) -> pl.DataFrame:
    """`alts`, `gt` and `ref_check` from `source_ref` and `gt_alleles` (at most two entries).

    `alts` are the called alleles other than `source_ref`, sorted by string; `gt` indexes
    `[source_ref] + alts`, null for an uncalled haplotype. A row without `source_ref` (an array-shaped
    row) gets `alts`, `gt` and `ref_check` null: they wait for `resolve_ref`. So does a row whose
    `gt_alleles` is null (its text could not be read).
    """
    has_ref = pl.col("source_ref").is_not_null() & pl.col("gt_alleles").is_not_null()
    called = pl.col("gt_alleles").list.drop_nulls().list.unique()
    alts = called.list.set_difference(pl.concat_list(pl.col("source_ref"))).list.sort()
    frame = frame.with_columns(pl.when(has_ref).then(alts).otherwise(None).alias("alts"))

    def index(slot: int) -> pl.Expr:
        allele = pl.col("gt_alleles").list.get(slot, null_on_oob=True)
        return (
            pl.when(allele.is_null())
            .then(None)
            .when(allele == pl.col("source_ref"))
            .then(0)
            .when(allele == pl.col("alts").list.get(0, null_on_oob=True))
            .then(1)
            .when(allele == pl.col("alts").list.get(1, null_on_oob=True))
            .then(2)
            .otherwise(None)
            .cast(pl.UInt8)
        )

    gt = pl.concat_list(index(0), index(1)).list.head(pl.col("gt_alleles").list.len())
    return frame.with_columns(
        pl.when(has_ref).then(gt).otherwise(None).alias("gt"),
        pl.when(pl.col("source_ref").is_not_null()).then(pl.lit(RefCheck.NOT_RUN.value)).otherwise(None).alias("ref_check"),
    )


def _letters(column: str, pattern: str) -> pl.Expr:
    """One allele per letter of `column`, null where the text does not match `pattern`."""
    text = pl.col(column)
    return pl.when(text.str.contains(pattern)).then(text.str.split("")).otherwise(None)


def _count_unreadable(frame: pl.DataFrame, table: SiteTable, what: str) -> None:
    unreadable = frame.filter(pl.col("gt_alleles").is_null()).height
    if unreadable:
        table.warnings.append(f"{unreadable} {what} not readable as allele letters kept with gt_alleles, alts and gt null and the text in source_gt")


_SITE_COLUMNS = (
    "chrom",
    "pos",
    "end_pos",
    "ids",
    "source_ref",
    "ref_check",
    "alts",
    "source_alleles",
    "gt_alleles",
    "gt",
    "source_gt",
    "gq",
    "dp",
    "ad",
    "qual",
    "evidence_basis",
    "source_chrom",
)


def _select(frame: pl.DataFrame) -> pl.DataFrame:
    """The site columns a reader made, plus null columns for the ones it did not, in one order."""
    defaults = {
        "source_alleles": pl.List(pl.String),
        "gq": pl.UInt16,
        "dp": pl.UInt32,
        "ad": pl.List(pl.UInt32),
        "qual": pl.Float32,
    }
    missing = [pl.lit(None, dtype=dtype).alias(name) for name, dtype in defaults.items() if name not in frame.columns]
    return frame.with_columns(*missing, pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis")).select(_SITE_COLUMNS)


# --- CASAVA ----------------------------------------------------------------------------------------

_CASAVA_REQUIRED = ("chromosome", "position", "ref_base", "call", "a_bases", "c_bases", "g_bases", "t_bases", "total_bases")
_CASAVA_USED = {*_CASAVA_REQUIRED, "dbsnp", "call_type"}

_CASAVA_CALL_TYPES = {"Diff", "Het1", "Het2", "Other"}


def _casava_snp_table(data: bytes, spec: FormatSpec) -> SiteTable:
    """CASAVA SNP calls: `Ref Base` is `source_ref`, `Call` the genotype letters.

    `Call Type` restates what `Call` and `Ref Base` already say (Het1: first letter is REF, Het2: the
    second, Diff: one letter other than REF, Other: two letters, neither REF). It is checked, not
    carried: rows where the two disagree are counted, and keep what `Call` says.
    """
    table = SiteTable(calls=pl.DataFrame(), read_evidence=["AD", "DP"])
    names, frame = _named_rows(data, spec, _CASAVA_REQUIRED, table)
    table.dropped_columns = _dropped(names, spec, _CASAVA_USED)
    frame = _positions(frame, "position", table)
    ref = pl.col("ref_base").str.to_uppercase()
    counts = {base: pl.col(f"{base.lower()}_bases").cast(pl.UInt32, strict=False) for base in _BASES}
    frame = frame.with_columns(
        *chrom_columns("chromosome"),
        pl.col("pos").alias("end_pos"),
        ids_expr("dbsnp").alias("ids"),
        pl.when(ref.str.contains(r"^[ACGTN]$")).then(pl.col("ref_base")).otherwise(None).alias("source_ref"),
        _letters("call", r"^[ACGT]{1,2}$").alias("gt_alleles"),
        pl.col("call").alias("source_gt"),
        pl.col("total_bases").cast(pl.UInt32, strict=False).alias("dp"),
    )
    bad_ref = frame.filter(pl.col("source_ref").is_null()).height
    if bad_ref:
        table.warnings.append(f"{bad_ref} row(s) whose Ref Base is not one base kept with source_ref null")
    _count_unreadable(frame, table, "call(s)")
    frame = _site_genotype(frame)
    allele_list = pl.concat_list(pl.col("source_ref"), pl.col("alts"))
    count_of = pl.element().replace_strict({base: i for i, base in enumerate(_BASES)}, default=None, return_dtype=pl.Int64)
    per_base = pl.concat_list(*counts.values())
    frame = frame.with_columns(
        pl.when(pl.col("alts").is_not_null() & pl.col("source_ref").is_in(list(_BASES)))
        .then(allele_list)
        .otherwise(None)
        .alias("source_alleles"),
    ).with_columns(
        pl.col("source_alleles").list.eval(count_of).alias("_base_index"),
        per_base.alias("_counts"),
    ).with_columns(
        # Guarded (CLAUDE.md §5): a gather over all-null index lists makes a list of a Null inner type.
        pl.when(pl.col("_base_index").is_not_null()).then(pl.col("_counts").list.gather(pl.col("_base_index"), null_on_oob=True)).alias("ad"),
    )
    expected = (
        pl.when(pl.col("gt").list.len() == 1)
        .then(pl.when(pl.col("gt").list.get(0) == 0).then(pl.lit(None, dtype=pl.String)).otherwise(pl.lit("Diff")))
        .when(pl.col("gt").list.get(0) == 0)
        .then(pl.lit("Het1"))
        .when(pl.col("gt").list.get(1) == 0)
        .then(pl.lit("Het2"))
        .otherwise(pl.lit("Other"))
    )
    if "call_type" in frame.columns:
        disagree = frame.filter(pl.col("gt").is_not_null() & (expected.fill_null("") != pl.col("call_type"))).height
        if disagree:
            table.warnings.append(f"{disagree} row(s) whose Call Type disagrees with what Call and Ref Base say; the Call is kept")
        unknown = frame.filter(~pl.col("call_type").is_in(list(_CASAVA_CALL_TYPES))).height
        if unknown:
            table.warnings.append(f"{unknown} row(s) with a Call Type other than {', '.join(sorted(_CASAVA_CALL_TYPES))}")
    table.calls = _select(frame)
    return table


_INDEL_REQUIRED = ("chromosome", "position", "ref_indel", "max_gtype", "q_indel", "q_max_gtype", "ref_reads", "indel_reads")
_INDEL_USED = {*_INDEL_REQUIRED}

_INDEL_GT: dict[str, list[int]] = {"het": [0, 1], "hom": [1, 1], "ref": [0, 0]}
"""What `max_gtype` states over `[REF, indel]`."""


def _casava_indel_table(data: bytes, spec: FormatSpec) -> SiteTable:
    """CASAVA indel calls: `ref/indel` is the reference span over the called allele, `-` padding.

    `pos` is the 1-based position of the span's first base, or, for a pure insertion, of the base
    after the insertion point. Measured on the one PGP-HMS table against NCBI36 (2026-10-23, 2,000
    sampled rows each): `ref_upstream` ends at the base before `pos` on every row, the deleted bases
    start at `pos` and `ref_downstream` follows them, and an insertion's `ref_downstream` starts at
    `pos`. So the dashes are dropped: `source_ref` is the span's bases at `pos`, the ALT the called
    bases, and either can be the empty allele `""` (SCHEMA Decision 2): a deletion's ALT, or an
    insertion's REF. An insertion point sits at the base before it, `pos - 1`, with `end_pos` null,
    as Genos and CGI insertion points do; the alleles stay unanchored until `resolve_ref` reads the
    anchor base from the reference. `ref_upstream` states that base too, but REF only comes from the
    sequence (charter 10), so the flanks are not read.

    `max_gtype` gives `gt` (`het` `[0, 1]`, `hom` `[1, 1]`); `Q(max_gtype)`, the genotype's quality, is
    `gq`; `Q(indel)`, the indel's, is `qual`; `ref_reads` and `indel_reads` are `ad`.
    """
    table = SiteTable(calls=pl.DataFrame(), read_evidence=["AD", "GQ"])
    names, frame = _named_rows(data, spec, _INDEL_REQUIRED, table)
    table.dropped_columns = _dropped(names, spec, _INDEL_USED)
    frame = _positions(frame, "position", table)
    parts = pl.col("ref_indel").str.split("/")
    ref = parts.list.get(0, null_on_oob=True).str.replace_all("-", "", literal=True).str.to_uppercase()
    alt = parts.list.get(1, null_on_oob=True).str.replace_all("-", "", literal=True).str.to_uppercase()
    readable = (
        (parts.list.len() == 2)
        & ref.str.contains(r"^[ACGTN]*$")
        & alt.str.contains(r"^[ACGTN]*$")
        & (ref != alt)
        & pl.col("max_gtype").is_in(list(_INDEL_GT))
    ).fill_null(False)
    insertion = ref == ""
    unread = frame.filter(~readable).height
    if unread:
        table.warnings.append(f"{unread} row(s) whose ref/indel or max_gtype cannot be read dropped")
    frame = frame.filter(readable)
    at_start = frame.filter(insertion & (pl.col("pos") == 1)).height
    if at_start:
        table.warnings.append(f"{at_start} insertion(s) before the contig's first base dropped")
        frame = frame.filter(~(insertion & (pl.col("pos") == 1)))
    gq = pl.col("q_max_gtype").cast(pl.Int64, strict=False)
    frame = frame.with_columns(
        *chrom_columns("chromosome"),
        pl.when(insertion).then(pl.col("pos") - 1).otherwise(pl.col("pos")).cast(pl.UInt32).alias("pos"),
        pl.when(insertion).then(None).otherwise(pl.col("pos") + ref.str.len_chars().cast(pl.UInt32) - 1).cast(pl.UInt32).alias("end_pos"),
        pl.lit(None, dtype=pl.List(pl.String)).alias("ids"),
        ref.alias("source_ref"),
        pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
        pl.concat_list(alt).alias("alts"),
        pl.concat_list(ref, alt).alias("source_alleles"),
        pl.col("max_gtype").replace_strict(_INDEL_GT, default=None, return_dtype=pl.List(pl.UInt8)).alias("gt"),
        pl.col("max_gtype").alias("source_gt"),
        pl.when(gq.is_between(0, 65535)).then(gq).otherwise(None).cast(pl.UInt16).alias("gq"),
        pl.col("q_indel").cast(pl.Float32, strict=False).alias("qual"),
        pl.concat_list(pl.col("ref_reads").cast(pl.UInt32, strict=False), pl.col("indel_reads").cast(pl.UInt32, strict=False)).alias("ad"),
    ).with_columns(pl.col("source_alleles").list.gather(pl.col("gt")).alias("gt_alleles"))
    sites = frame.filter(pl.col("end_pos").is_null()).height
    if sites:
        table.warnings.append(
            f'{sites} insertion point(s) placed at the base before them with end_pos null; the alleles stay unanchored (REF "") '
            "until resolve_ref reads the anchor base from the reference"
        )
    table.calls = _select(frame)
    return table


# --- SOAPsnp -----------------------------------------------------------------------------------------

_IUPAC: dict[str, list[str]] = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["C", "G"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
}
"""The IUPAC codes SOAPsnp writes as a consensus genotype: a base, or a pair of bases (sorted). The
three-base codes (B, D, H, V) and N name no genotype of at most two alleles, so they stay unread."""

_SOAP_CARRIED = {"chromosome", "position", "ref", "consensus", "consensus_quality", "depth"}


def _soapsnp_cns(data: bytes, spec: FormatSpec) -> SiteTable:
    """SOAPsnp consensus rows: `ref` is `source_ref`, the IUPAC `consensus` the genotype.

    A pair code (R, Y, …) is the two bases it names, sorted, since the code states no order; a plain
    base is that one base. `consensus_quality`, the Phred quality of the consensus genotype, is `gq`;
    `depth` is `dp`. Rows with depth 0 (prior-driven calls at dbSNP sites, no read evidence) keep
    `dp = 0`: the depth says what they rest on.
    """
    table = SiteTable(calls=pl.DataFrame(), read_evidence=["DP", "GQ"])
    names = [column["name"] for column in spec.columns]
    rows = _split(data, spec)[1]
    wrong = rows.filter(pl.col("fields").list.len() != len(names))
    if wrong.height:
        table.warnings.append(f"{wrong.height} row(s) without exactly {len(names)} fields dropped")
    frame = rows.filter(pl.col("fields").list.len() == len(names)).select(
        pl.col("fields").list.get(index).alias(name) for index, name in enumerate(names)
    )
    table.dropped_columns = [name for name in names if name not in _SOAP_CARRIED]
    frame = _positions(frame, "position", table)
    gq = pl.col("consensus_quality").cast(pl.Int64, strict=False)
    frame = frame.with_columns(
        *chrom_columns("chromosome"),
        pl.col("pos").alias("end_pos"),
        pl.lit(None, dtype=pl.List(pl.String)).alias("ids"),
        pl.when(pl.col("ref").str.contains(r"^[ACGTN]$")).then(pl.col("ref")).otherwise(None).alias("source_ref"),
        pl.col("consensus").replace_strict(_IUPAC, default=None, return_dtype=pl.List(pl.String)).alias("gt_alleles"),
        pl.col("consensus").alias("source_gt"),
        pl.when(gq.is_between(0, 65535)).then(gq).otherwise(None).cast(pl.UInt16).alias("gq"),
        pl.col("depth").cast(pl.UInt32, strict=False).alias("dp"),
    )
    _count_unreadable(frame, table, "consensus code(s)")
    out_of_range = frame.filter(gq.is_not_null() & pl.col("gq").is_null()).height
    if out_of_range:
        table.warnings.append(f"{out_of_range} consensus quality value(s) outside 0..65535 left out of gq")
    no_reads = frame.filter(pl.col("dp") == 0).height
    if no_reads:
        table.warnings.append(f"{no_reads} row(s) with depth 0: a consensus with no read evidence, kept with dp 0")
    table.calls = _select(_site_genotype(frame))
    return table


# --- FamilyTreeDNA Big Y ---------------------------------------------------------------------------

_BIGY_USED = {"position", "reference", "genotype", "snp_name"}


def _ftdna_bigy(data: bytes, spec: FormatSpec) -> SiteTable:
    """Big Y rows: every row is chrY, the product's one chromosome (the file names none).

    Only Novel Variant rows have a position; Known SNP rows name a Y-SNP and no position, so they
    cannot be `calls` rows and are dropped with a count. `Reference` is `source_ref`; `Genotype` one
    base, `?` where the site was not called (`[null]`).
    """
    table = SiteTable(calls=pl.DataFrame())
    names, frame = _named_rows(data, spec, ("position", "reference", "genotype"), table)
    table.dropped_columns = _dropped(names, spec, _BIGY_USED)
    unplaced = frame.filter(pl.col("position").fill_null("") == "").height
    if unplaced:
        table.warnings.append(f"{unplaced} row(s) without a position (Known SNP rows, named by Y-SNP only) dropped")
    frame = frame.filter(pl.col("position").fill_null("") != "")
    frame = _positions(frame, "position", table)
    genotype = pl.col("genotype")
    frame = frame.with_columns(
        pl.lit("Y").alias("chrom"),
        pl.lit(None, dtype=pl.String).alias("source_chrom"),
        pl.col("pos").alias("end_pos"),
        ids_expr("snp_name").alias("ids") if "snp_name" in frame.columns else pl.lit(None, dtype=pl.List(pl.String)).alias("ids"),
        pl.when(pl.col("reference").str.contains(r"^[ACGTN]$")).then(pl.col("reference")).otherwise(None).alias("source_ref"),
        pl.when(genotype.is_in(list(spec.no_call)))
        .then(pl.lit([None], dtype=pl.List(pl.String)))
        .when(genotype.str.contains(r"^[ACGT]$"))
        .then(pl.concat_list(genotype))
        .otherwise(None)
        .alias("gt_alleles"),
        genotype.alias("source_gt"),
    )
    bad_ref = frame.filter(pl.col("source_ref").is_null()).height
    if bad_ref:
        table.warnings.append(f"{bad_ref} row(s) whose Reference is not one base kept with source_ref null")
    _count_unreadable(frame, table, "genotype(s)")
    table.calls = _select(_site_genotype(frame))
    return table


# --- Color variant set -----------------------------------------------------------------------------

_COLOR_USED = {"sample_id", "chromosome", "position", "ref", "alt", "zygosity", "reference_version"}

_ZYGOSITY_GT: dict[str, list[int]] = {"Heterozygous": [0, 1], "Homozygous": [1, 1]}
"""What Color's `Zygosity` word states over `[Ref, Alt]`."""


def _color_variant_set(data: bytes, spec: FormatSpec) -> SiteTable:
    """Color's reportable variants: `Ref` and `Alt` as written, `Zygosity` the genotype.

    The table ends at the first row whose SampleID is empty or `Raw data disclaimer:`; the disclaimer
    after it is not data. Each row states one ALT, so `gt` is `[0, 1]` (Heterozygous) or `[1, 1]`
    (Homozygous), and `gt_alleles` the alleles those index. Only the panel's reportable variants are
    listed, so an absent site is unknown.
    """
    table = SiteTable(calls=pl.DataFrame())
    required = ("sample_id", "chromosome", "position", "ref", "alt", "zygosity")
    names, frame = _named_rows(data, spec, required, table, ends_at_empty="sample_id")
    table.dropped_columns = _dropped(names, spec, _COLOR_USED)
    samples = frame["sample_id"].unique(maintain_order=True).to_list()
    if len(samples) > 1:
        raise TableLayoutError(f"{spec.id}: {len(samples)} sample ids in one variant set")
    table.sample = samples[0] if samples else None
    if "reference_version" in frame.columns:
        stated = sorted({re.sub(r"\.p\d+$", "", value) for value in frame["reference_version"].to_list() if value})
        if len(stated) > 1:
            table.warnings.append(f"rows state {len(stated)} different builds ({', '.join(stated)}); assembly left unknown")
        elif stated and re.fullmatch(r"GRCh3[678]", stated[0]):
            table.build = stated[0]
    frame = _positions(frame, "position", table)
    sequence = r"^[ACGTN]+$"
    readable = pl.col("ref").str.contains(sequence) & pl.col("alt").str.contains(sequence)
    gt = pl.col("zygosity").replace_strict(_ZYGOSITY_GT, default=None, return_dtype=pl.List(pl.UInt8))
    frame = frame.with_columns(
        *chrom_columns("chromosome"),
        pl.when(readable).then(pl.col("ref")).otherwise(None).alias("source_ref"),
        pl.when(readable).then(pl.concat_list("alt")).otherwise(None).alias("alts"),
        pl.when(readable).then(gt).otherwise(None).alias("gt"),
        pl.col("zygosity").alias("source_gt"),
        pl.when(readable).then(pl.col("pos") + pl.col("ref").str.len_chars().cast(pl.UInt32) - 1).otherwise(None).alias("end_pos"),
        pl.lit(None, dtype=pl.List(pl.String)).alias("ids"),
    ).with_columns(
        # Guarded (CLAUDE.md §5): a gather over all-null index lists makes a list of a Null inner type.
        pl.when(pl.col("gt").is_not_null()).then(pl.concat_list("source_ref", "alts").list.gather(pl.col("gt"), null_on_oob=True)).alias("gt_alleles"),
        pl.when(pl.col("source_ref").is_not_null()).then(pl.lit(RefCheck.NOT_RUN.value)).otherwise(None).alias("ref_check"),
    )
    unread = frame.filter(~readable.fill_null(False)).height
    if unread:
        table.warnings.append(f"{unread} row(s) whose Ref or Alt is not a base sequence kept with source_ref, alts and gt null")
    unknown = frame.filter(readable.fill_null(False) & pl.col("gt").is_null()).height
    if unknown:
        table.warnings.append(f"{unknown} row(s) with a Zygosity other than Heterozygous or Homozygous kept with gt null")
    table.calls = _select(frame)
    return table


# --- GET-Evidence GFF ------------------------------------------------------------------------------

_GFF_FIELDS = 9
_GFF_CARRIED_ATTRIBUTES = {"alleles", "ref_allele", "db_xref"}


def _get_evidence_gff(data: bytes, spec: FormatSpec) -> SiteTable:
    """GET-Evidence genome GFF rows.

    A `REF` row (attributes `.`) states hom-ref at its span: a `ref` block, evidence for that span
    only. Every other row is a site whose `alleles` attribute names the alleles called (`A/G`, or one
    letter). Where `ref_allele` states REF it is `source_ref` and the row is indexed as a VCF row;
    where it does not (the build-36 layout), the row is array-shaped: `gt_alleles` only, `alts` and `gt`
    waiting for `resolve_ref`. An empty `alleles` value is a no-call, never hom-ref. Every
    `db_xref dbsnp:rs…` id goes in `ids`.
    """
    table = SiteTable(calls=pl.DataFrame())
    rows = _split(data, spec)[1]
    wrong = rows.filter(pl.col("fields").list.len() != _GFF_FIELDS)
    if wrong.height:
        table.warnings.append(f"{wrong.height} row(s) without exactly {_GFF_FIELDS} fields dropped")
    names = [column["name"] for column in spec.columns]
    frame = rows.filter(pl.col("fields").list.len() == _GFF_FIELDS).select(
        pl.col("fields").list.get(index).alias(name) for index, name in enumerate(names)
    )
    frame = _positions(frame, "start", table)
    frame = frame.with_columns(pl.col("end").cast(pl.Int64, strict=False).alias("_end"))
    backwards = frame.filter(pl.col("_end").is_null() | (pl.col("_end") < pl.col("pos").cast(pl.Int64))).height
    if backwards:
        table.warnings.append(f"{backwards} row(s) whose end is not a position at or after their start dropped")
    frame = frame.filter(pl.col("_end").is_not_null() & (pl.col("_end") >= pl.col("pos").cast(pl.Int64)))
    frame = frame.with_columns(pl.col("_end").cast(pl.UInt32).alias("end_pos"), *chrom_columns("chromosome"))

    attributes = pl.col("attributes")
    is_ref = (pl.col("type") == "REF") & (attributes == ".")
    blocks = frame.filter(is_ref).select("chrom", "pos", "end_pos", pl.lit("ref").alias("kind"))
    table.blocks = None if blocks.is_empty() else blocks
    sites = frame.filter(~is_ref)
    keys = (
        sites.filter(attributes != ".")
        .select(attributes.str.split(";").list.eval(pl.element().str.strip_chars().str.split(" ").list.first()).alias("key"))
        .explode("key", empty_as_null=False)
        .filter(pl.col("key").fill_null("") != "")
        ["key"]
        .unique()
        .sort()
        .to_list()
    )
    table.dropped_columns = ["source", "score", "strand", "phase", *[key for key in keys if key not in _GFF_CARRIED_ATTRIBUTES]]

    def attribute(key: str) -> pl.Expr:
        return attributes.str.extract(rf"(?:^|;)\s*{key}(?: ([^;]*))?(?:;|$)", 1)

    alleles = attribute("alleles")
    stated_ref = attribute("ref_allele").str.strip_chars()
    has_alleles = attributes.str.contains(r"(?:^|;)\s*alleles(?: |;|$)")
    allele_text = alleles.fill_null("").str.strip_chars()
    readable = allele_text.str.contains(r"^[ACGTN]+(/[ACGTN]+)?$")
    sites = sites.with_columns(
        pl.when(stated_ref.str.contains(r"^[ACGTN]+$")).then(stated_ref).otherwise(None).alias("source_ref"),
        pl.when(has_alleles & (allele_text == ""))
        .then(pl.lit([None], dtype=pl.List(pl.String)))
        .when(readable)
        .then(allele_text.str.split("/"))
        .otherwise(None)
        .alias("gt_alleles"),
        pl.when(has_alleles).then(allele_text).otherwise(None).alias("source_gt"),
        pl.when(attributes.str.contains(r"dbsnp:rs\d+"))
        .then(attributes.str.extract_all(r"dbsnp:rs\d+").list.eval(pl.element().str.strip_prefix("dbsnp:")))
        .otherwise(None)
        .alias("ids"),
    )
    no_calls = sites.filter(has_alleles & (allele_text == "")).height
    if no_calls:
        table.warnings.append(f"{no_calls} row(s) with an empty alleles value read as no-calls, not as hom-ref")
    _count_unreadable(sites, table, "alleles value(s)")
    unstated = sites.filter(pl.col("source_ref").is_null()).height
    if unstated:
        table.warnings.append(f"{unstated} row(s) state no ref_allele: their alts and gt wait for REF from the reference (resolve_ref)")
    if table.blocks is not None:
        table.warnings.append(f"{table.blocks.height} REF row(s) (hom-ref at a listed site) loaded as ref blocks")
    table.calls = _select(_site_genotype(sites))
    return table


_READERS = {
    "casava_snp_table": _casava_snp_table,
    "casava_indel_table": _casava_indel_table,
    "soapsnp_cns": _soapsnp_cns,
    "ftdna_bigy": _ftdna_bigy,
    "color_variant_set": _color_variant_set,
    "get_evidence_gff": _get_evidence_gff,
}
