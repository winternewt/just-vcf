"""Reading helpers shared by the text loaders: payload lines, genotype letters, and genotype text of
allele indices (the key=value table, Nirvana JSON)."""

import polars as pl

from just_vcf.schema import CALLS, canonical_chrom_expr

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


def two_field_alleles(first: str, second: str, no_call: tuple[str, ...]) -> pl.Expr:
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


READABLE_GT = r"^(\d+|\.)([/|](\d+|\.))*$"
"""Genotype text readable as allele indices: each slot an index or `.`, separated by `/` or `|`."""


GT_INDEX_MAX: int = CALLS["gt"].inner.max()
"""The largest allele index `gt` can hold (255: `gt` is a `u8` list, SCHEMA.md Decision 6)."""


def _text_indices(genotype: pl.Expr) -> pl.Expr:
    """The genotype text's allele indices, read from its digits; one too long for 64 bits as the largest."""
    biggest = pl.lit(2**63 - 1, dtype=pl.Int64)
    return genotype.str.extract_all(r"\d+").list.eval(pl.element().cast(pl.Int64, strict=False).fill_null(biggest))


def past_alts(genotype: pl.Expr, alts: pl.Expr) -> pl.Expr:
    """True on a genotype readable as allele indices (`READABLE_GT`) with an index past the row's `alts`
    ALT alleles, read from the text's digits: `1/300` over fewer than 300 ALTs is past where `gt`, a
    UInt8 list, nulls the 300, and an index too long for 64 bits is past too. An unreadable genotype is
    counted as that alone, so `1:2` is never also past (SCHEMA.md, the unreadable genotype)."""
    return genotype.str.contains(READABLE_GT) & (_text_indices(genotype).list.max() > alts)


def past_gt(genotype: pl.Expr, alts: pl.Expr) -> pl.Expr:
    """True on a genotype readable as allele indices whose every index names an allele of the row, one of
    them past what `gt` can hold (`GT_INDEX_MAX`): `0/300` over 300 ALTs, which `gt` reads `[0, null]`,
    the frame of the half-call `0/.`. A genotype past the row's ALTs is counted by `past_alts` alone, so
    each is counted under one reason at most."""
    biggest = _text_indices(genotype).list.max()
    return genotype.str.contains(READABLE_GT) & (biggest <= alts) & (biggest > GT_INDEX_MAX)


def genotype_indices(text: pl.Expr) -> pl.Expr:
    """`gt` from genotype text of allele indices (`0/1`, `1|2`, `./.`): one slot per haplotype the text
    separates, null per `.` and per index past `GT_INDEX_MAX` (counted by `past_alts` where it is past
    the row's ALTs, by `past_gt` where the row has the allele). Null where the text is not readable as
    allele indices (`READABLE_GT`): the genotype is unknown, so no slot is read, not the `0` of `0/x`
    either, which would read as the half-call `0/.` the file did not write (SCHEMA.md Decision 2).
    Gather `gt_alleles` behind a `gt` null guard: a batch whose every `gt` is null gathers a list of a
    Null inner type, which the spill's parquet sink cannot write."""
    slots = text.str.replace_all(r"\|", "/").str.split("/").list.eval(
        pl.when(pl.element() == ".").then(None).otherwise(pl.element()).cast(pl.UInt8, strict=False)
    )
    return pl.when(text.str.contains(READABLE_GT)).then(slots)


def genotype_phase(text: pl.Expr) -> pl.Expr:
    """`phased` from the same text, one entry per slot: true after `|`, false after `/`, and the first
    slot's read from the first separator, as for VCF before 4.4 (phased unless a `/` follows, so a
    haploid call is phased). Null where `genotype_indices` is."""
    seps = text.str.extract_all(r"[/|]")
    first = ~text.str.contains("/", literal=True)
    return pl.when(text.str.contains(READABLE_GT)).then(pl.concat_list(first, seps.list.eval(pl.element() == "|")))


UNREADABLE_WARNING = (
    "genotype(s) not readable as allele indices kept with gt, gt_alleles and phased null, and the text in source_gt"
)
"""The warning, after its count, for genotypes `genotype_indices` does not read."""

PAST_ALTS_WARNING = "genotype(s) index an allele past the row's ALT list; gt_alleles has a null there"
"""The warning, after its count, for genotypes `past_alts` counts, in every loader that reads allele indices."""

PAST_GT_WARNING = (
    f"genotype(s) hold an allele index past what gt can hold ({GT_INDEX_MAX}), though the row has that allele; "
    "gt and gt_alleles have a null there, as an uncalled slot does"
)
"""The warning, after its count, for genotypes `past_gt` counts, in every loader that reads allele indices."""
