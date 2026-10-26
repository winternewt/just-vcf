"""Which intervals cover a position, where intervals may overlap one another (internal).

A gVCF's blocks and CGI's tile the genome, so the block over a position is the last one starting at or
before it, which an as-of join finds. A lifted facet's blocks do not tile: two source regions a chain
sends to one place give overlapping blocks, which `liftover` keeps (charter 8; a lifted WGS gVCF carries
hundreds of thousands), so an outer block can cover a position past an inner one's end, and several
can cover one position. `covering` finds every covering interval either way, at the as-of join's cost
where nothing overlaps; `nearest_block` picks the one block that answers for a position.
"""

from collections.abc import Iterable, Sequence

import polars as pl

from just_vcf.schema import BlockKind

_REACHES = "_reaches_next"


def covering(probes: pl.LazyFrame, intervals: pl.LazyFrame, contigs: Iterable[str], *, at: str, start: str, end: str) -> pl.LazyFrame:
    """Every pair of a `probes` row and an `intervals` row on the same contig (`chrom`) whose closed span
    `start..end` holds the probe's `at` (all three Int64): the probe's columns, then the interval's
    without its `chrom`. A probe nothing covers is left out; only `contigs` are searched.

    Per contig in start order, an interval that reaches its successor's start is one a later interval
    may lie inside or cross; only those go through a range join (an IEJoin per contig). The others
    cannot overlap one another (each ends before the next one starts), so at most one of them covers a
    position, the last of them starting at or before it, which an as-of join finds. A tiled frame has
    no interval of the first kind, so it costs the as-of join alone.
    """
    contigs = sorted(set(contigs))
    ordered = (
        intervals.filter(pl.col("chrom").is_in(contigs))
        .sort("chrom", start, end)
        .with_columns((pl.col(end) >= pl.col(start).shift(-1).over("chrom")).fill_null(False).alias(_REACHES))
    )
    apart = ordered.filter(~pl.col(_REACHES)).drop(_REACHES)
    reaching = ordered.filter(pl.col(_REACHES)).drop(_REACHES)
    tiled = (
        probes.sort("chrom", at)
        .join_asof(apart, left_on=at, right_on=start, by="chrom", strategy="backward", check_sortedness=False)
        .filter(pl.col(start).is_not_null() & (pl.col(at) <= pl.col(end)))
    )
    # One range join per contig: two inequalities and no equality run as an IEJoin, where an equality on
    # `chrom` beside them would make it a filtered cross join.
    inside = [
        probes.filter(pl.col("chrom") == contig).join_where(
            reaching.filter(pl.col("chrom") == contig).drop("chrom"), pl.col(at) >= pl.col(start), pl.col(at) <= pl.col(end)
        )
        for contig in contigs
    ]
    return pl.concat([tiled, *inside], how="diagonal_relaxed")


def nearest_block(hits: pl.DataFrame, key: str | Sequence[str], *, start: str, end: str, kind: str = "kind", then: Sequence[str] = ()) -> pl.DataFrame:
    """One row per `key` of `hits` (`covering`'s blocks over each probe): the block `reference_support`
    answers with where several cover a position. The latest starting, then the shortest, then a
    `no_call` block before a `ref` one (a reference statement does not hide a no-call stated over the
    same base), then by `then` ascending, nulls last, so no tie depends on the order a join emits rows."""
    keys = [key] if isinstance(key, str) else list(key)
    is_ref = (pl.col(kind).cast(pl.String) == BlockKind.REF.value).fill_null(True)
    return (
        hits.with_columns(is_ref.alias("_is_ref"))
        .sort([*keys, start, end, "_is_ref", *then], descending=[*([False] * len(keys)), True, False, False, *([False] * len(then))], nulls_last=True, maintain_order=True)
        .unique(keys, keep="first", maintain_order=True)
        .drop("_is_ref")
    )
