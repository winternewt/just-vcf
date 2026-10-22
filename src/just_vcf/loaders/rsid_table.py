"""Load genotype tables in the 23andMe layout: identifier, chromosome, position, genotype.

Covers every format of `just_vcf.formats.rsid_tables` that has positions: 23andMe (tab and
spreadsheet-resaved comma), Living DNA, SelfDecode, PLINK `--recode 23`, Genes for Good (measured and
imputed) and the two preamble-less fallbacks. Each row becomes one `calls` row (SCHEMA.md): the
genotype text verbatim in `source_gt`, one entry per written allele in `gt_alleles`, and no REF,
because an array export writes none and REF only ever comes from the assembly sequence (charter 10).

Rows that do not fit the layout are dropped and counted, never guessed at: a row with the wrong
number of fields (except one whose only extra field is an empty one right after the id, F2), a position that is not a number, a genotype with letters outside the encoding. The
counts come back as warnings, one per reason.
"""

import logging
import re

import polars as pl

from just_vcf.formats import get_spec
from just_vcf.formats.spec import FormatSpec
from just_vcf.loaders.text import chrom_columns, ids_expr, lines, one_field_alleles, single_position_end
from just_vcf.schema import CALLS, EvidenceBasis, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset(
    {
        "23andme",
        "23andme_csv",
        "living_dna",
        "selfdecode",
        "plink_23",
        "genes_for_good",
        "genes_for_good_imputed",
        "rsid_table",
        "rsid_table_headerless",
    }
)
"""Format ids this loader reads."""

_FIELDS = ("rsid", "chromosome", "position", "genotype")

_COMMENT = r"^.?\"?#"
"""A preamble line: `#`, `"#` (a spreadsheet quotes lines holding a comma), or one stray byte before
the `#` (a real Genes for Good file opens `I# Genes for Good`). No data line starts this way: ids
start with `rs`, `i`, `AX-` or a chromosome number."""


def read_calls(data: bytes, spec: FormatSpec) -> tuple[pl.DataFrame, list[str]]:
    """Parse a decompressed payload of `spec`'s format into `calls` rows without `facet_id`.

    Returns the rows and the warnings, one per reason with its count.
    """
    if spec.id not in FORMATS:
        raise ValueError(f"{spec.id} is not a 23andMe-layout table")
    delimiter = spec.delimiter or "\t"
    warnings: list[str] = []

    rows = (
        lines(data).filter(pl.col("line").str.strip_chars() != "")
        .filter(~pl.col("line").str.contains(_COMMENT))
        .with_columns(pl.col("line").str.split(delimiter).alias("fields"))
    )
    # F2: a doubled delimiter right after the id (`rs15842⇥⇥1⇥948921⇥CT`). Only that place: every such
    # row in the PGP-HMS exports has it there, and its id sits at the same chromosome and position in
    # other exports. An empty field anywhere else is not read past.
    doubled = (pl.col("fields").list.len() == len(_FIELDS) + 1) & (pl.col("fields").list.get(1) == "")
    collapsed = rows.filter(doubled).height
    if collapsed:
        warnings.append(
            f"{collapsed} row(s) with an empty field between the id and the chromosome (a doubled delimiter) "
            "read without it"
        )
        rows = rows.with_columns(
            pl.when(doubled)
            .then(pl.concat_list(pl.col("fields").list.slice(0, 1), pl.col("fields").list.slice(2)))
            .otherwise(pl.col("fields"))
            .alias("fields")
        )
    wrong_width = rows.filter(pl.col("fields").list.len() != len(_FIELDS))
    if wrong_width.height:
        widths = sorted(set(wrong_width["fields"].list.len().to_list()))
        warnings.append(
            f"{wrong_width.height} row(s) without exactly {len(_FIELDS)} fields dropped "
            f"(widths seen: {', '.join(map(str, widths))})"
        )
    rows = rows.filter(pl.col("fields").list.len() == len(_FIELDS)).select(
        pl.col("fields").list.get(i).str.strip_chars().alias(name) for i, name in enumerate(_FIELDS)
    )

    rows = rows.with_columns(pl.col("position").cast(pl.UInt32, strict=False).alias("pos"))
    bad_position = rows.filter(pl.col("pos").is_null())
    if bad_position.height:
        warnings.append(f"{bad_position.height} row(s) with a position that is not a 32-bit count dropped")
    rows = rows.filter(pl.col("pos").is_not_null())

    rows = rows.with_columns(one_field_alleles("genotype", spec.no_call).alias("gt_alleles"))
    unreadable = rows.filter(pl.col("gt_alleles").is_null())
    if unreadable.height:
        spellings = sorted(set(unreadable["genotype"].str.len_chars().to_list()))
        warnings.append(
            f"{unreadable.height} genotype(s) not readable as one allele letter per haplotype kept with "
            f"gt_alleles null and the text in source_gt (lengths seen: {', '.join(map(str, spellings))})"
        )

    calls = rows.select(
        *chrom_columns("chromosome"),
        pl.col("pos"),
        single_position_end("gt_alleles", "pos").alias("end_pos"),
        ids_expr("rsid").alias("ids"),
        pl.col("gt_alleles"),
        pl.col("genotype").alias("source_gt"),
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
    )
    return calls, warnings


def to_calls_frame(calls: pl.DataFrame, facet_id: str) -> pl.LazyFrame:
    """`read_calls` rows as a `calls` frame for `facet_id`."""
    return conform(
        calls.lazy().drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")),
        CALLS,
    )


def column_line_ok(data: bytes, spec: FormatSpec) -> bool | None:
    """Whether the file's column line names the four columns in the expected order.

    None when the format has no column line (headerless) or none was found.
    """
    head = data[: 1 << 16].decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    delimiter = spec.delimiter or "\t"
    found = re.findall(rf"^#\s*(rsid{re.escape(delimiter)}[^\n]*)$", head, flags=re.MULTILINE)
    if not found:
        return None
    names = [name.strip().lower() for name in found[-1].split(delimiter)]
    return names == list(_FIELDS)


def spec_for(format_id: str) -> FormatSpec:
    """The registered spec for a format this loader reads."""
    if format_id not in FORMATS:
        raise ValueError(f"{format_id} is not a 23andMe-layout table")
    return get_spec(format_id)
