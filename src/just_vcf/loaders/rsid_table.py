"""Load genotype tables in the 23andMe layout: identifier, chromosome, position, genotype.

Covers every format of `just_vcf.formats.rsid_tables` that has positions: 23andMe (tab and
spreadsheet-resaved comma), Living DNA, SelfDecode, PLINK `--recode 23`, Genes for Good (measured and
imputed) and the two preamble-less fallbacks. Each row becomes one `calls` row (SCHEMA.md): the
genotype text verbatim in `source_gt`, one entry per written allele in `gt_alleles`, and no REF,
because an array export writes none and REF only ever comes from the assembly sequence (charter 10).

Rows that do not fit the layout are dropped and counted, never guessed at: a row with the wrong
number of fields, a position that is not a number, a genotype with letters outside the encoding. The
counts come back as warnings, one per reason.
"""

import logging
import re

import polars as pl

from just_vcf.formats import get_spec
from just_vcf.formats.spec import FormatSpec
from just_vcf.schema import CALLS, EvidenceBasis, canonical_chrom_expr, conform

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

_ALLELE_LETTERS = "ACGTDI"
"""Letters an allele may be written with: bases, and 23andMe's D/I codes for indel probes."""


def read_calls(data: bytes, spec: FormatSpec) -> tuple[pl.DataFrame, list[str]]:
    """Parse a decompressed payload of `spec`'s format into `calls` rows without `facet_id`.

    Returns the rows and the warnings, one per reason with its count.
    """
    if spec.id not in FORMATS:
        raise ValueError(f"{spec.id} is not a 23andMe-layout table")
    lines = _lines(data)
    delimiter = spec.delimiter or "\t"
    warnings: list[str] = []

    rows = (
        lines.filter(pl.col("line").str.strip_chars() != "")
        .filter(~pl.col("line").str.contains(_COMMENT))
        .with_columns(pl.col("line").str.split(delimiter).alias("fields"))
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

    rows = rows.with_columns(_gt_alleles(spec).alias("gt_alleles"))
    unreadable = rows.filter(pl.col("gt_alleles").is_null())
    if unreadable.height:
        spellings = sorted(set(unreadable["genotype"].str.len_chars().to_list()))
        warnings.append(
            f"{unreadable.height} genotype(s) not readable as one allele letter per haplotype kept with "
            f"gt_alleles null and the text in source_gt (lengths seen: {', '.join(map(str, spellings))})"
        )

    single_bases = pl.col("gt_alleles").list.eval(pl.element().is_null() | pl.element().is_in(list("ACGT"))).list.all()
    any_called = pl.col("gt_alleles").list.eval(pl.element().is_not_null()).list.any()
    calls = rows.select(
        canonical_chrom_expr("chromosome").alias("chrom"),
        pl.col("pos"),
        pl.when(single_bases & any_called).then(pl.col("pos")).otherwise(None).alias("end_pos"),
        pl.when(pl.col("rsid").is_in(["", "."])).then(None).otherwise(pl.concat_list("rsid")).alias("ids"),
        pl.col("gt_alleles"),
        pl.col("genotype").alias("source_gt"),
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
        pl.col("chromosome").alias("source_chrom"),
    )
    return calls, warnings


def to_calls_frame(calls: pl.DataFrame, facet_id: str) -> pl.LazyFrame:
    """`read_calls` rows as a `calls` frame for `facet_id`."""
    return conform(
        calls.lazy().drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")),
        CALLS,
    )


def _lines(data: bytes) -> pl.DataFrame:
    """The payload's lines as one string column. BOM-aware; CRLF and bare CR become LF."""
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        data = data.decode("utf-16").encode("utf-8")
    elif data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
    data = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    if not data.strip():
        return pl.DataFrame({"line": []}, schema={"line": pl.String})
    return pl.read_csv(
        data,
        has_header=False,
        separator="\x1f",
        quote_char=None,
        new_columns=["line"],
        schema={"line": pl.String},
        encoding="utf8-lossy",
    ).with_columns(pl.col("line").fill_null(""))


def _gt_alleles(spec: FormatSpec) -> pl.Expr:
    """One entry per written allele; null entries for missing alleles; null when unreadable.

    `--` is two missing alleles and `-` one, because the list length is what the file wrote. A `-`
    beside a letter (`T-`) is a half-call. A genotype longer than two letters (Living DNA writes an
    indel's two allele sequences back to back) cannot be split without the site's alleles, so it is
    unreadable here and stays in `source_gt`.
    """
    genotype = pl.col("genotype")
    readable = genotype.str.contains(rf"^[{_ALLELE_LETTERS}-]{{1,2}}$")
    if spec.no_call:
        readable = readable | genotype.is_in(list(spec.no_call))
    letters = genotype.str.split("").list.eval(pl.when(pl.element() == "-").then(None).otherwise(pl.element()))
    return pl.when(readable).then(letters).otherwise(None)


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
