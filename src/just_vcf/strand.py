"""The strand an array export's genotypes are written on, read from its own rows (F13).

An array export says its genotypes are on the forward (plus) strand, and `resolve_ref` derives `alts`
and `gt` against the reference base on that reading. Some exports write part of their probes on the
minus strand instead (one PGP-HMS FamilyTreeDNA export, about one probe in eight). On such a row the
derived alleles are wrong: a minus-strand hom-ref call reads as hom-alt.

The rows show it, for a two-allele site, through the heterozygous calls alone:

- a het carrying the reference base is safe on either reading: if its other allele is the base's
  complement (an A/T or C/G probe) it reads the same on both strands, and otherwise the minus reading
  would hold three alleles;
- a het carrying the complement of the reference base and not the base is impossible on the forward
  strand, and is a minus-strand row;
- a het carrying neither fits no strand (another build, a third allele), and is counted apart.

So `check_strand` scores a facet by the share of its sampled two-base hets that carry the complement and
not the base. On forward exports that share is near zero; on F13's export it is about one in five.

What a mixed facet cannot say, without a table of which probes the vendor wrote on which strand, is
which of its *homozygous* rows (and half-calls) are flipped: `AA` at a reference `A` is a forward
hom-ref call, or a minus-strand hom-alt call of an A/T probe; `TT` there is a minus-strand hom-ref call,
or a forward hom-alt call. Those rows keep `ref`, and their `alts` and `gt` are withheld, counted per
reason. A complement is never guessed, not even on the hets that prove a minus-strand probe: complementing
exactly the minus probes needs the vendor's strand table, and FamilyTreeDNA publishes none.
"""

from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.schema import StrandCheck

MIN_HETS = 200
"""Sampled two-base heterozygous rows needed before the check answers. FamilyTreeDNA's X-only exports
(about 2,000 mostly hemizygous rows) fall short and are left unchecked."""

MAX_FORWARD_MINUS = 0.02
"""The largest share of hets that carry the reference base's complement and not the base, at which a
facet still reads as forward. Measured on PGP-HMS array exports (docs/PREVIOUS_ISSUES.md F13)."""

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class StrandResult(BaseModel):
    """What `check_strand` found. `check` is null when it could not run."""

    model_config = ConfigDict(frozen=True)

    check: StrandCheck | None
    hets: int
    """Sampled rows with two distinct called bases, each a single base, scored on a read reference base."""
    minus: int
    """Of those, the hets carrying the reference base's complement and not the base."""
    neither: int
    """Of those, the hets carrying neither the base nor its complement."""

    @property
    def minus_share(self) -> float | None:
        return self.minus / self.hets if self.hets else None

    def summary(self) -> str:
        if self.check is None:
            return f"strand not checked: {self.hets} sampled heterozygous row(s), {MIN_HETS} needed; alleles derived as forward-strand"
        line = (
            f"strand checked on {self.hets} sampled heterozygous rows: {self.minus} ({self.minus_share:.2%}) carry the "
            f"complement of the reference base and not the base, {self.neither} carry neither: {self.check.value}"
        )
        return line


def check_strand(sample: pl.DataFrame, bases: pl.Series) -> StrandResult:
    """Score `sample_snv_rows`' rows (`alleles`: distinct called bases) against the reference `bases`
    parallel to them (null or empty where none was read)."""
    rows = sample.select("alleles").with_columns(bases.alias("_base")).filter(
        pl.col("_base").is_in(list(COMPLEMENT)).fill_null(False) & (pl.col("alleles").list.len() == 2)
    )
    complement = pl.col("_base").replace_strict(COMPLEMENT, default=None)
    scored = rows.select(
        pl.col("alleles").list.contains(pl.col("_base")).alias("_ref"),
        pl.col("alleles").list.contains(complement).alias("_comp"),
    )
    hets = scored.height
    minus = scored.filter(~pl.col("_ref") & pl.col("_comp")).height
    neither = scored.filter(~pl.col("_ref") & ~pl.col("_comp")).height
    if hets < MIN_HETS:
        return StrandResult(check=None, hets=hets, minus=minus, neither=neither)
    check = StrandCheck.MIXED if minus / hets > MAX_FORWARD_MINUS else StrandCheck.FORWARD
    return StrandResult(check=check, hets=hets, minus=minus, neither=neither)


class Withheld(str, Enum):
    """Why a mixed facet's row keeps `ref` but has its `alts` and `gt` withheld."""

    MINUS = "minus"
    """A het carrying the reference base's complement and not the base: a minus-strand probe."""
    NEITHER = "neither"
    """A het carrying neither the base nor its complement."""
    UNDECIDABLE = "undecidable"
    """A homozygous call or a half-call: forward and minus readings give different alleles."""


def strand_safe(alleles: pl.Expr, ref: pl.Expr) -> pl.Expr:
    """Whether a row's genotype reads the same on both strands or only on the forward one: fully called,
    two distinct bases, one of them the reference base."""
    called = alleles.list.drop_nulls()
    return (called.list.len() == alleles.list.len()) & (called.list.unique().list.len() == 2) & called.list.contains(ref)


def withheld_reason(alleles: pl.Expr, ref: pl.Expr) -> pl.Expr:
    """`Withheld` for a row `strand_safe` rejects, null for a safe one."""
    called = alleles.list.drop_nulls()
    distinct = called.list.unique()
    het = (called.list.len() == alleles.list.len()) & (distinct.list.len() == 2)
    complement = ref.replace_strict(COMPLEMENT, default=None)
    return (
        pl.when(strand_safe(alleles, ref))
        .then(None)
        .when(het & distinct.list.contains(complement))
        .then(pl.lit(Withheld.MINUS.value))
        .when(het)
        .then(pl.lit(Withheld.NEITHER.value))
        .otherwise(pl.lit(Withheld.UNDECIDABLE.value))
    )
