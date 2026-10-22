"""Reading helpers shared by the delimited-text loaders: payload lines and genotype letters."""

import polars as pl

from just_vcf.schema import canonical_chrom_expr

ALLELE_LETTERS = "ACGTDI"
"""Letters an array allele may be written with: bases, and the D/I codes of indel probes."""

_BASES = list("ACGT")


def lines(data: bytes) -> pl.DataFrame:
    """The payload's lines as one string column `line`. BOM-aware; CRLF and bare CR become LF."""
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


def one_field_alleles(column: str, no_call: tuple[str, ...], no_call_width: int | None = None) -> pl.Expr:
    """`gt_alleles` from a genotype written as letters in one field (`AG`, `A`, `T-`, `--`, `DI`).

    One entry per written letter, because the list length is what the file wrote (SCHEMA.md
    Decision 2). A `-` is one missing allele, so `T-` is a half-call and `--` two missing alleles.

    A whole no-call spelling of the format (`spec.no_call`) is one missing allele per character, unless
    `no_call_width` says how many alleles the format writes on every called row. Then the spelling is a
    token for "no call" and gives that many nulls: FTDNA's undocumented `---` and an Illumina Final
    Report's single `-` stand in fields that hold two letters on every called row, so both are two
    missing alleles, not three and not one. 23andMe writes one letter on haploid sites, so its `-` stays
    one.

    Anything else (a genotype longer than two letters, a letter outside the encoding) is null: the text
    cannot be split into haplotypes here and stays in `source_gt`.
    """
    genotype = pl.col(column)
    letters = genotype.str.split("")
    readable = genotype.str.contains(rf"^[{ALLELE_LETTERS}-]{{1,2}}$")
    called = letters.list.eval(pl.when(pl.element() == "-").then(None).otherwise(pl.element()))
    if no_call_width is None:
        missing = letters.list.eval(pl.lit(None, dtype=pl.String))
    else:
        missing = pl.lit([None] * no_call_width, dtype=pl.List(pl.String))
    return (
        pl.when(genotype.is_in(list(no_call)))
        .then(missing)
        .when(readable)
        .then(called)
        .otherwise(None)
    )


def single_position_end(alleles: str, pos: str) -> pl.Expr:
    """`end_pos` of an array call: its own position when every called allele is one base.

    Null when nothing was called, or when an allele is an indel code (`I`, `D`) whose REF span is not
    known until the reference is read.
    """
    gt = pl.col(alleles)
    single_bases = gt.list.eval(pl.element().is_null() | pl.element().is_in(_BASES)).list.all()
    any_called = gt.list.eval(pl.element().is_not_null()).list.any()
    return pl.when(single_bases & any_called).then(pl.col(pos)).otherwise(None)


def ids_expr(column: str) -> pl.Expr:
    """`ids` from a one-identifier column; null where the source wrote none."""
    return pl.when(pl.col(column).is_in(["", "."])).then(None).otherwise(pl.concat_list(column))


def chrom_columns(column: str) -> list[pl.Expr]:
    """The canonical `chrom` and the source spelling (kept for the facet's `contig_map`)."""
    return [canonical_chrom_expr(column).alias("chrom"), pl.col(column).alias("source_chrom")]
