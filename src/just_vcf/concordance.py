"""Concordance between two facets of one person on one assembly (RM277, charter invariants 7 and 8).

`concordance(a, b)` compares two facets allele by allele and returns a `Concordance`: one row per
allele either facet states (`sites`), a `summary` of counts, the alleles left out of the comparison
and why. Disagreement is kept as rows, never resolved to a winner. The rules are in docs/SCHEMA.md
("Concordance"); in short:

1. **One assembly.** Both facets must have had REF read (`resolve_ref`): their `ref_seqcol_digest`
   names the collection. Equal digests compare every contig; different collections compare the
   contigs whose refget sequence is the same in both (Ensembl `1` and UCSC `chr1` of GRCh38 are one
   sequence), and a primary contig whose sequence differs refuses the pair: lift one facet first.
   Without a digest on either side, an equal stated `assembly_label` is taken with a warning
   (charter 3: a label is not the assembly); anything else is refused.
2. **The unit is one ALT allele** (SCHEMA.md Decision 1's per-allele view). Alleles meet by their
   VRS id when both facets minted one (the ALT slot of `vrs_allele_ids`, never the REF slot), which is
   independent of how each file spelled an indel, so no `normalize` is needed; otherwise by
   `(chrom, pos, REF, ALT)` as written. An allele neither meets is paired with the other facet's one
   row at that position, if it wrote exactly one (an array's hom-ref row, a DeepVariant `RefCall`, the
   other ALT of a site): that row's genotype says how many copies of this allele it holds (key `site`).
3. **Genotypes are compared as copies of the allele**, from `gt_alleles`, unordered (phase is not
   comparable across facets): `match`, `mismatch`, `ploidy_differs`, or `uncalled` where a haplotype
   is null and the called ones do not contradict. A `./.` is the source saying it could not tell,
   never a mismatch.
4. **Graded where both carry vectors** (charter 7): each side's vector becomes a distribution over
   the copies of the allele, and `p_same` is the probability the two agree, `sum_k pa(k) * pb(k)`.
   Likelihoods (`gl`) are preferred, turned into probabilities under a flat prior over the genotypes
   the vector enumerates (`gl_flat_prior`, the prior named); posteriors (`gp`) are used as given,
   renormalized to sum to 1. A likelihood is never compared with a posterior (Decision 3). Anything
   else is `hard`, and `p_same` is null: a hard call has no probability, and none is made up.
5. **One-sided alleles** say what the silent facet's silence means: a `ref` or `no_call` block over
   the position, `other_alleles` (rows at the position, none holding the allele), else the facet's
   own rule (`absence_is_reference`: `reference`, or `unknown`, never a negative).

Coordinates: every position column is 1-based (charter 13); `region` takes 1-based closed intervals.
"""

import logging
from collections import Counter
from collections.abc import Mapping
from enum import Enum
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.callset import NOT_A_VARIANT, absence_is_reference, variant_class_expr
from just_vcf.loading import Loaded
from just_vcf.reference import PRIMARY, ReferenceMismatchError, ReferenceStore, _store_for
from just_vcf.schema import CONTIG_RANK, BlockKind, VariantClass, canonical_chrom, canonical_chrom_expr, enum_dtype

logger = logging.getLogger(__name__)


class Agreement(str, Enum):
    """How the two facets' genotypes for one allele compare."""

    MATCH = "match"
    """Both fully called, the same number of copies of the allele, the same ploidy."""
    MISMATCH = "mismatch"
    """Both called and the copies differ; with a null haplotype, the called ones already contradict
    (`1/.` against `0/0`)."""
    PLOIDY_DIFFERS = "ploidy_differs"
    """Both genotypes written, with different numbers of haplotypes (a haploid X call against a diploid one)."""
    UNCALLED = "uncalled"
    """A haplotype is null on either side and the called ones do not contradict: not decided."""
    ONLY_A = "only_a"
    """Only facet `a` states the allele, and `b` gave no single row to read it against (`absence` says
    what `b`'s silence means)."""
    ONLY_B = "only_b"
    """Only facet `b` states the allele (`absence` is `a`'s silence)."""


class PairKey(str, Enum):
    """How an allele of one facet met the other facet."""

    VRS = "vrs"
    """The two rows' VRS ids for the allele are equal."""
    POSITION = "position"
    """The rows state the same `(chrom, pos, REF, ALT)`, as written (no VRS id on at least one side)."""
    SITE = "site"
    """The other facet does not state the allele, and wrote exactly one row at its position: the
    genotype of that row is read for copies of the allele."""


class AbsenceReading(str, Enum):
    """What the facet that does not state an allele says about it."""

    OTHER_ALLELES = "other_alleles"
    """It wrote rows at the position, none with this allele."""
    REF_BLOCK = "ref_block"
    """A `ref` block covers the position."""
    NO_CALL_BLOCK = "no_call_block"
    """A `no_call` block covers the position."""
    REFERENCE = "reference"
    """No row and no block: `absence_is_reference` is True for the allele's class."""
    UNKNOWN = "unknown"
    """No row and no block, and the facet's absence says nothing (`absence_is_reference` is None, or
    the facet never had the contig)."""


class ComparisonBasis(str, Enum):
    """Which evidence `p_same` was computed from."""

    HARD = "hard"
    """Genotypes only: `p_same` is null."""
    GL_FLAT_PRIOR = "gl_flat_prior"
    """Both rows' likelihoods, each turned into probabilities under a flat prior over its genotypes."""
    GP = "gp"
    """Both rows' posteriors, renormalized to sum to 1."""


SITES = pl.Schema(
    {
        "chrom": pl.String,
        "pos": pl.UInt32,
        "ref": pl.String,
        "alt": pl.String,
        "variant_class": enum_dtype(VariantClass),
        "vrs_allele_id": pl.String,
        "key": enum_dtype(PairKey),
        "agreement": enum_dtype(Agreement),
        "absence": enum_dtype(AbsenceReading),
        "basis": enum_dtype(ComparisonBasis),
        "p_same": pl.Float64,
        "a_pos": pl.UInt32,
        "a_gt_alleles": pl.List(pl.String),
        "a_copies": pl.UInt8,
        "a_filters": pl.List(pl.String),
        "a_gq": pl.UInt16,
        "b_pos": pl.UInt32,
        "b_gt_alleles": pl.List(pl.String),
        "b_copies": pl.UInt8,
        "b_filters": pl.List(pl.String),
        "b_gq": pl.UInt16,
    }
)
"""One row per allele either facet states. `pos`, `ref` and `alt` are the spelling of the facet that
states it (`a` first); `a_pos` and `b_pos` are where each facet wrote it, which differ for an indel
the two spelled at different positions and met by VRS id. `a_copies` and `b_copies` are the copies of
the allele in each genotype, null where a haplotype is uncalled."""

SUMMARY = pl.Schema(
    {
        "agreement": enum_dtype(Agreement),
        "key": enum_dtype(PairKey),
        "absence": enum_dtype(AbsenceReading),
        "variant_class": enum_dtype(VariantClass),
        "basis": enum_dtype(ComparisonBasis),
        "alleles": pl.UInt32,
        "p_same_mean": pl.Float64,
    }
)

EXCLUDED_REASONS: frozenset[str] = frozenset(
    {"unspecified", "spanning_deletion", "not_a_variant", "symbolic", "not_sequence", "ref_unread", "alts_unread", "contig_not_shared"}
)
"""Why an allele was left out of the comparison. An open vocabulary."""

_MAX_GENOTYPE_ALLELES = 64
"""The largest allele list whose diploid genotype order the graded comparison enumerates."""


class Concordance(BaseModel):
    """What `concordance` returns. A measurement over two facets, not a facet."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    a_facet_id: str
    b_facet_id: str
    sites: pl.LazyFrame
    """`SITES`' columns, one row per allele, in canonical order."""
    summary: pl.DataFrame
    """`summarize(sites)`: alleles per agreement, key, absence, variant class and basis."""
    excluded: pl.DataFrame
    """Alleles not compared, per facet (`a`, `b`) and reason (`EXCLUDED_REASONS`)."""
    warnings: tuple[str, ...] = ()


def summarize(sites: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
    """Alleles per `(agreement, key, absence, variant_class, basis)`, with the mean `p_same` of the
    graded ones; every count sums to the height of `sites`. Run it again on a filtered `sites` (PASS
    only, alleles either facet carries) to summarize that subset."""
    frame = sites.lazy() if isinstance(sites, pl.DataFrame) else sites
    groups = ["agreement", "key", "absence", "variant_class", "basis"]
    return (
        frame.group_by(groups)
        .agg(pl.len().cast(pl.UInt32).alias("alleles"), pl.col("p_same").mean().alias("p_same_mean"))
        .sort(groups, nulls_last=True)
        .collect()
        .select(SUMMARY.names())
        .cast(dict(SUMMARY))
    )


def concordance(
    a: Loaded,
    b: Loaded,
    *,
    region: pl.LazyFrame | pl.DataFrame | None = None,
    reference: "ReferenceStore | Path | str | None" = None,
) -> Concordance:
    """Compare facet `a` with facet `b`, allele by allele (module docstring, docs/SCHEMA.md).

    Each `Loaded` holds one facet. `region` keeps only the alleles whose position lies in one of its
    intervals: a frame with `chrom`, `pos` and `end_pos`, 1-based and closed, as `blocks` writes them
    (a BED's 0-based start is `pos - 1`). `reference` is the store holding both facets' collections,
    read only when the collections differ (default: `reference_dir()`).

    Raises `ReferenceMismatchError` when the two are not on one assembly (lift one first) and
    `ValueError` when a `Loaded` holds several facets or the assembly cannot be told.
    """
    fa = _one_facet(a, "a")
    fb = _one_facet(b, "b")
    warnings: list[str] = []
    shared = _shared_contigs(fa, fb, reference, warnings)
    excluded: Counter[tuple[str, str]] = Counter()
    side_a = _alleles(a.calls, shared, region, excluded, "a")
    side_b = _alleles(b.calls, shared, region, excluded, "b")

    pairs = _pair(side_a, side_b, warnings)
    paired_a = pairs.select("_a_row", "_a_k").unique()
    paired_b = pairs.select("_b_row", "_b_k").unique()
    rest_a = side_a.join(paired_a, left_on=["_row", "_k"], right_on=["_a_row", "_a_k"], how="anti")
    rest_b = side_b.join(paired_b, left_on=["_row", "_k"], right_on=["_b_row", "_b_k"], how="anti")
    alone_a = _one_sided(rest_a, b, fb, "a", "b", warnings) if rest_a.height else None
    alone_b = _one_sided(rest_b, a, fa, "b", "a", warnings) if rest_b.height else None

    pieces = [_two_sided(pairs) if pairs.height else None, alone_a, alone_b]
    present = [_conform(piece) for piece in pieces if piece is not None]
    sites = pl.concat(present, how="vertical") if present else pl.DataFrame(schema={**SITES, **_ROWS})
    sites = _ordered(sites).select(SITES.names())
    excluded_frame = pl.DataFrame(
        [{"facet": side, "reason": reason, "alleles": count} for (side, reason), count in sorted(excluded.items())],
        schema={"facet": pl.String, "reason": pl.String, "alleles": pl.UInt64},
    )
    for line in warnings:
        logger.warning("concordance %s / %s: %s", fa["facet_id"], fb["facet_id"], line)
    return Concordance(
        a_facet_id=fa["facet_id"],
        b_facet_id=fb["facet_id"],
        sites=sites.lazy(),
        summary=summarize(sites),
        excluded=excluded_frame,
        warnings=tuple(warnings),
    )


def _one_facet(loaded: Loaded, name: str) -> dict:
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"concordance takes one facet per side; `{name}` holds {facets.height}")
    return facets.row(0, named=True)


def _shared_contigs(fa: Mapping, fb: Mapping, reference: "ReferenceStore | Path | str | None", warnings: list[str]) -> set[str] | None:
    """The canonical contigs on which the two facets can be compared; None for all of them."""
    da, db = fa["ref_seqcol_digest"], fb["ref_seqcol_digest"]
    ids = f"facets {fa['facet_id']} and {fb['facet_id']}"
    if da is not None and db is not None:
        if da == db:
            return None
        store, _ = _store_for(reference)
        held = set(store.collections())
        missing = sorted({da, db} - held)
        if missing:
            raise ReferenceMismatchError(f"{ids}: the reference store does not hold collection(s) {', '.join(missing)}, which REF was read from")
        seq_a = {chrom: found[0].digest for chrom, found in store.contigs_by_chrom(da).items() if len(found) == 1}
        seq_b = {chrom: found[0].digest for chrom, found in store.contigs_by_chrom(db).items() if len(found) == 1}
        both = set(seq_a) & set(seq_b)
        differ = sorted(chrom for chrom in both & PRIMARY if seq_a[chrom] != seq_b[chrom])
        shared = {chrom for chrom in both if seq_a[chrom] == seq_b[chrom]}
        if differ or not shared & PRIMARY:
            reason = f"primary contig(s) {', '.join(differ)} differ in sequence" if differ else "no primary contig has one sequence in both"
            raise ReferenceMismatchError(f"{ids}: not on one assembly ({reason}); lift one onto the other's assembly first")
        warnings.append(f"collections {da} and {db} differ; {len(shared)} contigs share one sequence and are compared, the rest are not")
        return shared
    la, lb = fa["assembly_label"], fb["assembly_label"]
    if la is not None and la == lb:
        warnings.append(f"assembly not checked by content: REF not read on both facets, both labelled {la}")
        return None
    if la is not None and lb is not None:
        raise ReferenceMismatchError(f"{ids}: {la} against {lb}; lift one onto the other's assembly first")
    raise ValueError(f"{ids}: cannot tell they are on one assembly; run resolve_ref on both first")


def _alleles(calls: pl.LazyFrame, shared: set[str] | None, region: pl.LazyFrame | pl.DataFrame | None, excluded: Counter, side: str) -> pl.DataFrame:
    """One row per sequence ALT allele of the facet, with what the comparison reads."""
    rows = calls.select("chrom", "pos", "ref", "source_ref", "alts", "vrs_allele_ids", "gt_alleles", "source_alleles", "gl", "gp", "filters", "gq")
    rows = rows.with_row_index("_row")
    unread = rows.filter(pl.col("alts").is_null() & pl.col("gt_alleles").is_not_null()).select(pl.len()).collect().item()
    if unread:
        excluded[(side, "alts_unread")] += unread
    exploded = (
        rows.filter(pl.col("alts").list.len() > 0)
        .with_columns(pl.int_ranges(pl.col("alts").list.len(), dtype=pl.UInt32).alias("_k"))
        .explode(["alts", "_k"], empty_as_null=False)
        .rename({"alts": "alt"})
        .with_columns(pl.coalesce("ref", "source_ref").alias("_ref"))
        .with_columns(
            pl.col("vrs_allele_ids").list.get(pl.col("_k").cast(pl.Int64) + 1, null_on_oob=True).alias("_vid"),
            variant_class_expr(pl.col("_ref"), pl.col("alt")).alias("variant_class"),
        )
        .drop("ref", "source_ref", "vrs_allele_ids")
        .collect()
    )
    reason = (
        pl.when(pl.col("alt") == "<*>")
        .then(pl.lit("unspecified"))
        .when(pl.col("alt") == "*")
        .then(pl.lit("spanning_deletion"))
        .when(pl.col("alt").is_in(list(NOT_A_VARIANT)))
        .then(pl.lit("not_a_variant"))
        .when(pl.col("variant_class") == VariantClass.SV.value)
        .then(pl.lit("symbolic"))
        .when(pl.col("_ref").is_null())
        .then(pl.lit("ref_unread"))
        .when(pl.col("variant_class").is_null())
        .then(pl.lit("not_sequence"))
    )
    tagged = exploded.with_columns(reason.alias("_excluded"))
    if shared is not None:
        tagged = tagged.with_columns(
            pl.when(pl.col("_excluded").is_null() & ~pl.col("chrom").is_in(sorted(shared))).then(pl.lit("contig_not_shared")).otherwise(pl.col("_excluded")).alias("_excluded")
        )
    for reason_name, count in tagged.group_by("_excluded").len().drop_nulls("_excluded").iter_rows():
        excluded[(side, reason_name)] += count
    kept = tagged.filter(pl.col("_excluded").is_null()).drop("_excluded")
    if region is not None:
        kept = _in_region(kept, region)
    return kept


def _in_region(frame: pl.DataFrame, region: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
    """The rows of `frame` whose `pos` lies in one of `region`'s 1-based closed intervals."""
    probe = frame.with_columns(pl.col("pos").cast(pl.Int64).alias("_at")).sort("chrom", "_at")
    return probe.filter(_covered(probe, _merged_intervals(region))).drop("_at")


def _merged_intervals(region: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
    """`region`'s 1-based closed intervals on canonical contigs, overlapping or touching ones merged,
    as `chrom`, `_start`, `_end`, sorted: the one interval starting at or before a position is then
    the only one that can hold it."""
    intervals = (region.lazy() if isinstance(region, pl.DataFrame) else region).select(
        canonical_chrom_expr("chrom").alias("chrom"), pl.col("pos").cast(pl.Int64), pl.col("end_pos").cast(pl.Int64)
    )
    return (
        intervals.sort("chrom", "pos", "end_pos")
        .with_columns(pl.col("end_pos").cum_max().shift(1).over("chrom").alias("_reach"))
        .with_columns((pl.col("_reach").is_null() | (pl.col("pos") > pl.col("_reach") + 1)).cum_sum().over("chrom").alias("_run"))
        .group_by("chrom", "_run")
        .agg(pl.col("pos").min().alias("_start"), pl.col("end_pos").max().alias("_end"))
        .drop("_run")
        .sort("chrom", "_start")
        .collect()
    )


def _covered(probe: pl.DataFrame, merged: pl.DataFrame) -> pl.Series:
    """Per row of `probe` (sorted by `chrom`, `_at`; `_at` an Int64 1-based position), whether one of
    `merged`'s intervals (`_merged_intervals`) holds it."""
    hit = probe.select("chrom", "_at").join_asof(merged, left_on="_at", right_on="_start", by="chrom", strategy="backward", check_sortedness=False)
    return (hit["_end"].is_not_null() & (hit["_at"] <= hit["_end"])).rename("_covered")


def _pair(side_a: pl.DataFrame, side_b: pl.DataFrame, warnings: list[str]) -> pl.DataFrame:
    """Alleles stated by both facets: first by VRS id, then the rest by `(chrom, pos, REF, ALT)`."""
    a = side_a.rename({name: f"_a_{name.lstrip('_')}" if name.startswith("_") else f"a_{name}" for name in side_a.columns})
    b = side_b.rename({name: f"_b_{name.lstrip('_')}" if name.startswith("_") else f"b_{name}" for name in side_b.columns})
    by_vrs = a.filter(pl.col("_a_vid").is_not_null()).join(
        b.filter(pl.col("_b_vid").is_not_null()), left_on="_a_vid", right_on="_b_vid", how="inner", coalesce=False
    ).with_columns(pl.lit(PairKey.VRS.value).alias("key"))
    rest_a = a.join(by_vrs.select("_a_row", "_a_k"), on=["_a_row", "_a_k"], how="anti")
    rest_b = b.join(by_vrs.select("_b_row", "_b_k"), on=["_b_row", "_b_k"], how="anti")
    by_position = rest_a.join(
        rest_b,
        left_on=["a_chrom", "a_pos", "_a_ref", "a_alt"],
        right_on=["b_chrom", "b_pos", "_b_ref", "b_alt"],
        how="inner",
        coalesce=False,
    ).with_columns(pl.lit(PairKey.POSITION.value).alias("key"))
    pairs = pl.concat([by_vrs, by_position], how="vertical_relaxed")
    many = pairs.height - min(pairs.select("_a_row", "_a_k").unique().height, pairs.select("_b_row", "_b_k").unique().height)
    if many:
        warnings.append(f"{many} pairings beyond one per allele: a facet states the same allele in several rows, and every pairing is kept")
    return pairs


def _copies(gt: str, alt: str, max_ploidy: int) -> tuple[pl.Expr, pl.Expr, pl.Expr]:
    """Copies of `alt` among the called haplotypes of `gt`, the null haplotypes, and the ploidy."""
    hits = [(pl.col(gt).list.get(slot, null_on_oob=True) == pl.col(alt)).fill_null(False).cast(pl.Int32) for slot in range(max_ploidy)]
    called = pl.sum_horizontal(hits) if hits else pl.lit(0, dtype=pl.Int32)
    nulls = pl.col(gt).list.eval(pl.element().is_null()).list.sum().cast(pl.Int32)
    return called, nulls, pl.col(gt).list.len().cast(pl.Int32)


def _agreement(a: tuple[pl.Expr, pl.Expr, pl.Expr], b: tuple[pl.Expr, pl.Expr, pl.Expr], a_gt: str, b_gt: str) -> pl.Expr:
    ca, na, pa = a
    cb, nb, pb = b
    a_silent = pl.col(a_gt).is_null() | (pa == 0) | (na == pa)
    b_silent = pl.col(b_gt).is_null() | (pb == 0) | (nb == pb)
    return (
        pl.when(a_silent | b_silent)
        .then(pl.lit(Agreement.UNCALLED.value))
        .when(pa != pb)
        .then(pl.lit(Agreement.PLOIDY_DIFFERS.value))
        .when((na == 0) & (nb == 0))
        .then(pl.when(ca == cb).then(pl.lit(Agreement.MATCH.value)).otherwise(pl.lit(Agreement.MISMATCH.value)))
        .when((ca > cb + nb) | (cb > ca + na))
        .then(pl.lit(Agreement.MISMATCH.value))
        .otherwise(pl.lit(Agreement.UNCALLED.value))
    )


def _max_ploidy(frame: pl.DataFrame, *columns: str) -> int:
    lengths = [frame.select(pl.col(name).list.len().max()).item() for name in columns if name in frame.columns]
    return max([int(n) for n in lengths if n is not None], default=0)


def _two_sided(pairs: pl.DataFrame) -> pl.DataFrame:
    """Pairs with their agreement, copies, and graded `p_same`."""
    pairs = pairs.with_row_index("_pair")
    ploidy = _max_ploidy(pairs, "a_gt_alleles", "b_gt_alleles")
    a = _copies("a_gt_alleles", "a_alt", ploidy)
    b = _copies("b_gt_alleles", "b_alt", ploidy)
    graded = _graded(pairs)
    return (
        pairs.join(graded, on="_pair", how="left")
        .with_columns(
            pl.col("a_chrom").alias("chrom"),
            pl.col("a_pos").alias("pos"),
            pl.col("_a_ref").alias("ref"),
            pl.col("a_alt").alias("alt"),
            pl.col("a_variant_class").alias("variant_class"),
            pl.coalesce("_a_vid", "_b_vid").alias("vrs_allele_id"),
            _agreement(a, b, "a_gt_alleles", "b_gt_alleles").alias("agreement"),
            pl.lit(None, dtype=pl.String).alias("absence"),
            pl.coalesce("basis", pl.lit(ComparisonBasis.HARD.value)).alias("basis"),
            pl.when(a[1] == 0).then(a[0]).alias("a_copies"),
            pl.when(b[1] == 0).then(b[0]).alias("b_copies"),
        )
    )


def _distribution(pairs: pl.DataFrame, side: str, vector: str) -> pl.DataFrame:
    """Per pair, the probability of each number of copies of the allele under `side`'s `vector`.

    The vector is in VCF genotype order over `source_alleles` (Decision 3): haploid when it has one
    entry per allele, diploid when it has `n(n+1)/2`; any other length, a null entry, an allele list
    that does not hold the ALT exactly once, or a total of zero gives no distribution.
    """
    alleles, alt, values = f"{side}_source_alleles", f"{side}_alt", f"{side}_{vector}"
    usable = pairs.filter(pl.col(values).is_not_null() & pl.col(alleles).is_not_null())
    if usable.is_empty():
        return pl.DataFrame(schema={"_pair": pl.UInt32, "_copies": pl.Int32, "_p": pl.Float64})
    width = min(_max_ploidy(usable, alleles), _MAX_GENOTYPE_ALLELES)
    found = [pl.when(pl.col(alleles).list.get(j, null_on_oob=True) == pl.col(alt)).then(pl.lit(j, dtype=pl.Int32)) for j in range(width)]
    hits = pl.sum_horizontal([(pl.col(alleles).list.get(j, null_on_oob=True) == pl.col(alt)).fill_null(False).cast(pl.Int32) for j in range(width)])
    n = pl.col(alleles).list.len().cast(pl.Int64)
    length = pl.col(values).list.len().cast(pl.Int64)
    work = (
        usable.select(
            "_pair",
            pl.col(values).alias("_v"),
            pl.coalesce(found).alias("_i"),
            hits.alias("_hits"),
            pl.when(length == n).then(1).when(length == n * (n + 1) // 2).then(2).alias("_ploidy"),
            pl.col(values).list.eval(pl.element().is_null()).list.any().alias("_holes"),
            (n <= _MAX_GENOTYPE_ALLELES).alias("_small"),
        )
        .filter((pl.col("_hits") == 1) & pl.col("_ploidy").is_not_null() & ~pl.col("_holes") & pl.col("_small"))
        .with_columns(pl.int_ranges(pl.col("_v").list.len(), dtype=pl.Int64).alias("_g"))
        .explode(["_v", "_g"], empty_as_null=False)
    )
    order = _diploid_order()
    work = work.join(order, on="_g", how="left").with_columns(
        pl.when(pl.col("_ploidy") == 1)
        .then((pl.col("_g") == pl.col("_i")).cast(pl.Int32))
        .otherwise((pl.col("_j") == pl.col("_i")).cast(pl.Int32) + (pl.col("_k") == pl.col("_i")).cast(pl.Int32))
        .alias("_copies"),
        pl.col("_v").cast(pl.Float64),
    )
    if vector == "gl":
        # Log10 likelihoods to probabilities under a flat prior over the genotypes the vector lists.
        work = work.with_columns((pl.lit(10.0) ** (pl.col("_v") - pl.col("_v").max().over("_pair"))).alias("_w"))
    else:
        work = work.with_columns(pl.col("_v").alias("_w"))
    total = pl.col("_w").sum().over("_pair")
    work = work.filter(total > 0).with_columns((pl.col("_w") / total).alias("_w"))
    return work.group_by("_pair", "_copies").agg(pl.col("_w").sum().alias("_p"))


def _diploid_order() -> pl.DataFrame:
    """VCF's diploid genotype order: genotype `j/k` (j <= k) at index `k(k+1)/2 + j`."""
    rows = [(k * (k + 1) // 2 + j, j, k) for k in range(_MAX_GENOTYPE_ALLELES) for j in range(k + 1)]
    return pl.DataFrame(rows, schema={"_g": pl.Int64, "_j": pl.Int32, "_k": pl.Int32}, orient="row")


def _graded(pairs: pl.DataFrame) -> pl.DataFrame:
    """`basis` and `p_same` for the pairs both of whose rows carry a usable vector of one kind."""
    found = []
    taken = pl.DataFrame(schema={"_pair": pl.UInt32})
    for vector, basis in (("gl", ComparisonBasis.GL_FLAT_PRIOR), ("gp", ComparisonBasis.GP)):
        da = _distribution(pairs, "a", vector)
        db = _distribution(pairs, "b", vector)
        both = da.select("_pair").unique().join(db.select("_pair").unique(), on="_pair").join(taken, on="_pair", how="anti")
        if both.is_empty():
            continue
        products = da.join(db, on=["_pair", "_copies"], how="inner").group_by("_pair").agg((pl.col("_p") * pl.col("_p_right")).sum().alias("p_same"))
        found.append(
            both.join(products, on="_pair", how="left").with_columns(pl.col("p_same").fill_null(0.0), pl.lit(basis.value).alias("basis"))
        )
        taken = pl.concat([taken, both.select("_pair")])
    if not found:
        return pl.DataFrame(schema={"_pair": pl.UInt32, "p_same": pl.Float64, "basis": pl.String})
    return pl.concat(found).select("_pair", "p_same", "basis")


def _one_sided(rest: pl.DataFrame, other: Loaded, other_facet: Mapping, side: str, other_side: str, warnings: list[str]) -> pl.DataFrame:
    """Alleles of `side` that met no allele of `other`: read against `other`'s one row at the position,
    else against its blocks and its absence rule."""
    mine = rest.rename({name: f"{side}_{name}" for name in ("gt_alleles", "filters", "gq")}).with_columns(pl.col("pos").alias(f"{side}_pos"))
    wanted = mine.select("chrom", "pos").unique().lazy()
    rows_there = (
        other.calls.select("chrom", "pos", pl.coalesce("ref", "source_ref").str.to_uppercase().alias("_other_ref"), "gt_alleles", "filters", "gq")
        .with_row_index("_other_row")
        .join(wanted, on=["chrom", "pos"], how="semi")
        .group_by("chrom", "pos")
        .agg(
            pl.len().alias("_rows"),
            pl.col("gt_alleles").first().alias(f"{other_side}_gt_alleles"),
            pl.col("filters").first().alias(f"{other_side}_filters"),
            pl.col("gq").first().alias(f"{other_side}_gq"),
            pl.col("_other_row").first().alias(f"_{other_side}_row"),
            pl.col("_other_ref").first().alias("_other_ref"),
        )
        .collect()
    )
    mine = mine.join(rows_there, on=["chrom", "pos"], how="left")
    # Read only a row over the same REF span: a longer allele there (a haplotype caller's complex
    # allele) holds no allele of this span, so its genotype says nothing about this allele's copies.
    single = ((pl.col("_rows") == 1) & (pl.col("_other_ref") == pl.col("_ref").str.to_uppercase())).fill_null(False)
    several = (pl.col("_rows") > 1).fill_null(False)
    n_several = mine.filter(several).height
    if n_several:
        warnings.append(f"{n_several} alleles of facet {side} meet several rows of facet {other_side} at their position and none with the allele; not read against them")
    mine = mine.with_columns(
        *[pl.when(single).then(pl.col(name)).alias(name) for name in (f"{other_side}_gt_alleles", f"{other_side}_filters", f"{other_side}_gq", f"_{other_side}_row")],
        pl.when(single).then(pl.col("pos")).alias(f"{other_side}_pos"),
    )
    mine = _block_reading(mine, other)
    classes = {cls.value: absence_is_reference(other_facet, cls) for cls in VariantClass}
    reference_classes = [name for name, verdict in classes.items() if verdict is True]
    declared = _declared_contigs(other, other_facet)
    undeclared = mine.filter(pl.col("_rows").is_null() & pl.col("_block").is_null() & ~pl.col("chrom").is_in(sorted(declared))).height
    if undeclared:
        warnings.append(f"{undeclared} alleles of facet {side} lie on contigs facet {other_side} neither declared nor called; its absence there is unknown")
    absence = (
        pl.when(pl.col("_rows").is_not_null())
        .then(pl.lit(AbsenceReading.OTHER_ALLELES.value))
        .when(pl.col("_block") == BlockKind.REF.value)
        .then(pl.lit(AbsenceReading.REF_BLOCK.value))
        .when(pl.col("_block").is_not_null())
        .then(pl.lit(AbsenceReading.NO_CALL_BLOCK.value))
        .when(pl.col("chrom").is_in(sorted(declared)) & pl.col("variant_class").is_in(reference_classes))
        .then(pl.lit(AbsenceReading.REFERENCE.value))
        .otherwise(pl.lit(AbsenceReading.UNKNOWN.value))
    )
    ploidy = _max_ploidy(mine, f"{side}_gt_alleles", f"{other_side}_gt_alleles")
    own = _copies(f"{side}_gt_alleles", "alt", ploidy)
    theirs = _copies(f"{other_side}_gt_alleles", "alt", ploidy)
    pair = (own, theirs) if side == "a" else (theirs, own)
    lone = Agreement.ONLY_A.value if side == "a" else Agreement.ONLY_B.value
    return mine.with_columns(
        pl.col("_ref").alias("ref"),
        pl.col("_vid").alias("vrs_allele_id"),
        pl.when(single).then(pl.lit(PairKey.SITE.value)).alias("key"),
        pl.when(single).then(_agreement(*pair, "a_gt_alleles", "b_gt_alleles")).otherwise(pl.lit(lone)).alias("agreement"),
        absence.alias("absence"),
        pl.when(single).then(pl.lit(ComparisonBasis.HARD.value)).alias("basis"),
        pl.lit(None, dtype=pl.Float64).alias("p_same"),
        pl.when(own[1] == 0).then(own[0]).alias(f"{side}_copies"),
        pl.when(single & (theirs[1] == 0)).then(theirs[0]).alias(f"{other_side}_copies"),
        pl.col("_row").alias(f"_{side}_row"),
    )


def _block_reading(mine: pl.DataFrame, other: Loaded) -> pl.DataFrame:
    """`_block`: the kind of the other facet's block over the position, where no row of it is there."""
    blocks = other.blocks.select("chrom", pl.col("pos").cast(pl.Int64).alias("_start"), pl.col("end_pos").cast(pl.Int64).alias("_end"), pl.col("kind").cast(pl.String))
    probe = mine.with_columns(pl.col("pos").cast(pl.Int64).alias("_at"))
    needed = probe.filter(pl.col("_rows").is_null()).select("chrom").unique()
    blocks = blocks.join(needed.lazy(), on="chrom", how="semi").sort("chrom", "_start").collect()
    if blocks.is_empty():
        return mine.with_columns(pl.lit(None, dtype=pl.String).alias("_block"))
    hit = probe.sort("chrom", "_at").join_asof(blocks, left_on="_at", right_on="_start", by="chrom", strategy="backward", check_sortedness=False)
    return hit.with_columns(
        pl.when(pl.col("_rows").is_null() & pl.col("_end").is_not_null() & (pl.col("_at") <= pl.col("_end"))).then(pl.col("kind")).alias("_block")
    ).drop("_at", "_start", "_end", "kind")


def _declared_contigs(loaded: Loaded, facet: Mapping) -> set[str]:
    """Contigs the facet stated (its header) or wrote rows on: where its silence can mean anything."""
    declared = {canonical_chrom(entry["name"]) for entry in facet.get("header_contigs") or []}
    for frame in (loaded.calls, loaded.blocks):
        declared |= set(frame.select(pl.col("chrom").unique()).collect()["chrom"].to_list())
    return declared


_ROWS = pl.Schema({"_a_row": pl.UInt32, "_b_row": pl.UInt32})
"""The `calls` rows of `a` and `b` a site came from: the last tie-break of the order, then dropped."""


def _conform(piece: pl.DataFrame) -> pl.DataFrame:
    """`piece` with exactly `SITES`' columns and `_ROWS`, typed; a column it lacks is null."""
    wanted = {**SITES, **_ROWS}
    return piece.select(
        [pl.col(name).cast(dtype) if name in piece.columns else pl.lit(None, dtype=dtype).alias(name) for name, dtype in wanted.items()]
    )


def _ordered(sites: pl.DataFrame) -> pl.DataFrame:
    """Canonical order: contig rank, position, REF, ALT; then the rows each side came from."""
    rank = pl.col("chrom").replace_strict(CONTIG_RANK, default=len(CONTIG_RANK), return_dtype=pl.UInt16)
    return sites.sort([rank, pl.col("chrom"), pl.col("pos"), pl.col("ref"), pl.col("alt"), pl.col("_a_row"), pl.col("_b_row")], nulls_last=True)
