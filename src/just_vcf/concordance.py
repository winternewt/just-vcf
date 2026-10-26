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
   comparable across facets): `match`, `mismatch`, or `uncalled` where a haplotype is null and the
   called ones do not contradict. A haploid call against the same allele's diploid homozygote is a
   `match` (both hold it on every copy), with `ploidy_differs` set apart. A `./.` is the source
   saying it could not tell, never a mismatch. Each side reads its own row: another ALT of a complex
   row can hold the allele's change (`_holds`), and one that cannot be read is open, never zero
   copies (`unreadable`).
4. **Graded where both carry vectors** (charter 7): each side's vector becomes a distribution over
   the copies of the allele, and `p_same` is the probability the two agree, `sum_k pa(k) * pb(k)`.
   Likelihoods (`gl`) are preferred, turned into probabilities under a flat prior over the genotypes
   the vector enumerates (`gl_flat_prior`, the prior named); posteriors (`gp`) are used as given,
   renormalized to sum to 1. A likelihood is never compared with a posterior (Decision 3). Anything
   else is `hard`, and `p_same` is null: a hard call has no probability, and none is made up.
5. **One-sided alleles** say what the silent facet's silence means: a `ref` or `no_call` block over
   the position, `other_alleles` (rows at the position, none holding the allele), else the facet's
   own rule (`absence_is_reference`: `reference`, or `unknown`, never a negative).
6. **Structural variants meet by type and extent** (RM282), as truvari bench matches them
   (`SvMatching`, its defaults): a symbolic allele (`<DEL>`, `<DUP:TANDEM>`, `<CN0>`) pairs with at
   most one of the other facet's, of the same type (the ALT's name before any `:`), whose extent lies
   within `refdist` of its own (each widened by its CIPOS/CIEND), of similar size (`pctsize`) and
   overlap (`pctovl`); the best-scoring pair is taken first (key `extent`). Genotypes compare as in 3
   and 4. A one-sided SV reads the silent facet's SVs within `refdist` (`other_alleles`), else its
   absence rule for the `sv` class, so a facet that called no SVs leaves it `unknown`.
7. **Breakends meet by both ends and the join** (RM282), as truvari bench's `bnd_match` matches them
   (`SvMatching.bnddist`): one contig and positions within reach, the mate's contig (parsed from ALT)
   and position within reach, and the same one of the four bracket forms (key `breakend`). The unit
   is the adjacency: a pair of mates written as two records is one row. A single breakend (`.A`)
   and an ALT that does not parse are counted apart (`single_breakend`, `breakend_unparsed`).

Coordinates: every position column is 1-based (charter 13); `region` takes 1-based closed intervals.
"""

import logging
from collections import Counter
from collections.abc import Callable, Mapping
from enum import Enum
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict, Field

from just_vcf.callset import NOT_A_VARIANT, absence_is_reference, variant_class_expr
from just_vcf.intervals import covering, nearest_block
from just_vcf.loading import Loaded
from just_vcf.reference import PRIMARY, ReferenceMismatchError, ReferenceStore, _store_for
from just_vcf.schema import CONTIG_ALIASES, CONTIG_RANK, BlockKind, VariantClass, canonical_chrom, canonical_chrom_expr, enum_dtype

logger = logging.getLogger(__name__)


class Agreement(str, Enum):
    """How the two facets' genotypes for one allele compare."""

    MATCH = "match"
    """Both fully called with the same copies of the allele; where the ploidies differ, each genotype
    holds it on every haplotype, or on none (a haploid `1` against `1/1`: `ploidy_differs` says so)."""
    MISMATCH = "mismatch"
    """Both called and the copies differ; with a null haplotype, the called ones already contradict
    (`1/.` against `0/0`)."""
    UNCALLED = "uncalled"
    """A haplotype is null on either side and the called ones do not contradict: not decided."""
    UNREADABLE = "unreadable"
    """A genotype holds another allele of its row whose change, of another length, covers this allele's
    span, so whether it holds this allele cannot be read, and the rest do not contradict: not decided,
    and not the source's no-call (F21)."""
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
    EXTENT = "extent"
    """Two structural-variant alleles of one type whose extents meet under `SvMatching`'s thresholds
    (truvari's matching), each taken once."""
    BREAKEND = "breakend"
    """Two breakend adjacencies that meet under truvari's `bnd_match`: both ends within
    `SvMatching.bnddist` of each other's, on the same contigs, joined the same way; each taken once."""


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
        "ploidy_differs": pl.Boolean,
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
        "sv_kind": pl.String,
        "b_alt": pl.String,
        "a_end_pos": pl.UInt32,
        "b_end_pos": pl.UInt32,
        "size_similarity": pl.Float64,
        "reciprocal_overlap": pl.Float64,
        "ci_needed": pl.Boolean,
    }
)
"""One row per allele either facet states. `pos`, `ref` and `alt` are the spelling of the facet that
states it (`a` first); `a_pos` and `b_pos` are where each facet wrote it, which differ for an indel
the two spelled at different positions and met by VRS id. `a_copies` and `b_copies` are the copies of
the allele in each genotype, null where a haplotype is uncalled. `ploidy_differs` is True where the two
genotypes were written with different numbers of haplotypes, False where with the same, null where
either side has no genotype.

The last seven columns are a structural-variant allele's, null on any other: `sv_kind` is the type it
was matched by (`BND` for a breakend); on a pair met by extent or as breakends, `b_alt` is `b`'s ALT
as written (`<DUP>` beside `a`'s `<DUP:TANDEM>`), `size_similarity` and `reciprocal_overlap` are
truvari's measures of the two extents (null for breakends, which have none), and `ci_needed` is True
where the two met only once their confidence intervals widened them; `a_end_pos` and `b_end_pos` are
each stating side's extent (a breakend's is its own position). A breakend row is one adjacency: where
a facet wrote it as two mate records, the row shows the record of the pair that met (both on
`chrom`), else the one at the adjacency's first end in canonical order."""

SUMMARY = pl.Schema(
    {
        "agreement": enum_dtype(Agreement),
        "ploidy_differs": pl.Boolean,
        "key": enum_dtype(PairKey),
        "absence": enum_dtype(AbsenceReading),
        "variant_class": enum_dtype(VariantClass),
        "sv_kind": pl.String,
        "basis": enum_dtype(ComparisonBasis),
        "alleles": pl.UInt32,
        "p_same_mean": pl.Float64,
    }
)

EXCLUDED_REASONS: frozenset[str] = frozenset(
    {
        "unspecified",
        "spanning_deletion",
        "not_a_variant",
        "single_breakend",
        "breakend_unparsed",
        "not_sequence",
        "ref_unread",
        "alts_unread",
        "contig_not_shared",
    }
)
"""Why an allele was left out of the comparison. An open vocabulary. `single_breakend` is a breakend
with no mate (`.A`, `A.`), which truvari leaves out by default (`no_single_bnd`); `breakend_unparsed`
an ALT with brackets that is none of VCF's four breakend forms (no `contig:position`, brackets that
do not pair). `contig_not_shared` covers a breakend whose mate lies on such a contig too."""

_MAX_GENOTYPE_ALLELES = 64
"""The largest allele list whose diploid genotype order the graded comparison enumerates."""


class SvMatching(BaseModel):
    """How `concordance` matches structural-variant alleles by type and extent (docs/SCHEMA.md,
    "Concordance"). The names and defaults are truvari bench's (truvari 5.5, `VariantParams.DEFAULTS`;
    its `VariantRecord.var_match` is the rule): `refdist` bases between the two extents, `pctsize` the
    least size similarity, `pctovl` the least reciprocal overlap, `typeignore` to pair across types,
    `bnddist` the reach of a breakend's ends (its `VariantRecord.bnd_match`).
    truvari's `pctseq` has no counterpart: no symbolic allele carries its sequence in the frame, as with
    truvari's `--pctseq 0`."""

    model_config = ConfigDict(frozen=True)

    refdist: int = Field(500, ge=0)
    pctsize: float = Field(0.7, ge=0.0, le=1.0)
    pctovl: float = Field(0.0, ge=0.0, le=1.0)
    typeignore: bool = False
    bnddist: int = Field(100, ge=0)
    """truvari's `--bnddist`, for breakends: each one's position, and its mate's, is widened by this
    many bases on either side and the two windows must overlap strictly, so two ends meet when they
    lie less than `2 * bnddist` apart (plus each record's interval). At 0 only an interval can make
    two meet, as in truvari."""


class Concordance(BaseModel):
    """What `concordance` returns. A measurement over two facets, not a facet."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    a_facet_id: str
    b_facet_id: str
    sites: pl.LazyFrame
    """`SITES`' columns, one row per allele, in canonical order."""
    summary: pl.DataFrame
    """`summarize(sites)`: alleles per agreement, ploidy flag, key, absence, variant class, SV type and basis."""
    excluded: pl.DataFrame
    """Alleles not compared, per facet (`a`, `b`) and reason (`EXCLUDED_REASONS`)."""
    warnings: tuple[str, ...] = ()


def summarize(sites: pl.LazyFrame | pl.DataFrame) -> pl.DataFrame:
    """Alleles per `(agreement, ploidy_differs, key, absence, variant_class, sv_kind, basis)`, with the mean
    `p_same` of the graded ones; every count sums to the height of `sites`. Run it again on a filtered `sites` (PASS
    only, alleles either facet carries) to summarize that subset."""
    frame = sites.lazy() if isinstance(sites, pl.DataFrame) else sites
    groups = ["agreement", "ploidy_differs", "key", "absence", "variant_class", "sv_kind", "basis"]
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
    sv: SvMatching | None = None,
) -> Concordance:
    """Compare facet `a` with facet `b`, allele by allele (module docstring, docs/SCHEMA.md).

    Each `Loaded` holds one facet. `region` keeps only the alleles whose position lies in one of its
    intervals: a frame with `chrom`, `pos` and `end_pos`, 1-based and closed, as `blocks` writes them
    (a BED's 0-based start is `pos - 1`). `reference` is the store holding both facets' collections,
    read only when the collections differ (default: `reference_dir()`). `sv` sets how structural-variant
    alleles are matched (default: `SvMatching()`, truvari's defaults).

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
    is_sv = pl.col("variant_class") == VariantClass.SV.value
    is_bnd = is_sv & pl.col("_form").is_not_null()
    bnd_a, bnd_b = side_a.filter(is_bnd), side_b.filter(is_bnd)
    sv_a, sv_b = (side.filter(is_sv & ~is_bnd).drop(_BND_INPUTS) for side in (side_a, side_b))
    side_a, side_b = (side.filter(~is_sv).drop(*_SV_INPUTS, *_BND_INPUTS) for side in (side_a, side_b))
    params = sv or SvMatching()

    pairs = _pair(side_a, side_b, warnings)
    paired_a = pairs.select("_a_row", "_a_k").unique()
    paired_b = pairs.select("_b_row", "_b_k").unique()
    rest_a = side_a.join(paired_a, left_on=["_row", "_k"], right_on=["_a_row", "_a_k"], how="anti")
    rest_b = side_b.join(paired_b, left_on=["_row", "_k"], right_on=["_b_row", "_b_k"], how="anti")
    alone_a = _one_sided(rest_a, b, fb, "a", "b", warnings) if rest_a.height else None
    alone_b = _one_sided(rest_b, a, fa, "b", "a", warnings) if rest_b.height else None

    pieces = [_two_sided(pairs, warnings) if pairs.height else None, alone_a, alone_b]
    if sv_a.height or sv_b.height:
        pieces += _sv_sites(sv_a, sv_b, (a, fa), (b, fb), params, warnings)
    if bnd_a.height or bnd_b.height:
        pieces += _bnd_sites(bnd_a, bnd_b, (a, fa), (b, fb), params, warnings)
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
    rows = calls.select("chrom", "pos", "ref", "source_ref", "alts", "vrs_allele_ids", "gt_alleles", "source_alleles", "gl", "gp", "filters", "gq", *_SV_INPUTS)
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
        .with_columns(_breakend_parts(pl.col("alt")))
        .collect()
    )
    reason = (
        pl.when(pl.col("alt") == "<*>")
        .then(pl.lit("unspecified"))
        .when(pl.col("alt") == "*")
        .then(pl.lit("spanning_deletion"))
        .when(pl.col("alt").is_in(list(NOT_A_VARIANT)))
        .then(pl.lit("not_a_variant"))
        .when((pl.col("variant_class") == VariantClass.SV.value) & pl.col("alt").str.contains(_SINGLE_BREAKEND))
        .then(pl.lit("single_breakend"))
        .when((pl.col("variant_class") == VariantClass.SV.value) & ~pl.col("alt").str.starts_with("<") & pl.col("_form").is_null())
        .then(pl.lit("breakend_unparsed"))
        .when(pl.col("variant_class") == VariantClass.SV.value)
        .then(None)
        .when(pl.col("_ref").is_null())
        .then(pl.lit("ref_unread"))
        .when(pl.col("variant_class").is_null())
        .then(pl.lit("not_sequence"))
    )
    tagged = exploded.with_columns(reason.alias("_excluded"))
    if shared is not None:
        elsewhere = ~pl.col("chrom").is_in(sorted(shared)) | (pl.col("_mate_chrom").is_not_null() & ~pl.col("_mate_chrom").is_in(sorted(shared)))
        tagged = tagged.with_columns(
            pl.when(pl.col("_excluded").is_null() & elsewhere).then(pl.lit("contig_not_shared")).otherwise(pl.col("_excluded")).alias("_excluded")
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


HOLDS, HOLDS_NOT, UNREADABLE = 1, 0, -1
"""What `_holds` says of one haplotype: it holds the allele's change, it does not, or that cannot be read."""

_BASES = frozenset("ACGT")


def _change(ref: str, allele: str) -> tuple[int, int, str]:
    """`allele`'s change against `ref` once their common prefix, then suffix, is trimmed: the span
    `[start, end)` of `ref` it replaces and the bases it puts there (both empty: no change)."""
    start = 0
    limit = min(len(ref), len(allele))
    while start < limit and ref[start] == allele[start]:
        start += 1
    tail = 0
    while tail < limit - start and ref[len(ref) - 1 - tail] == allele[len(allele) - 1 - tail]:
        tail += 1
    return start, len(ref) - tail, allele[start : len(allele) - tail]


def _overlaps(first: tuple[int, int], second: tuple[int, int]) -> bool:
    """Whether two changes' spans of REF meet; an insertion (an empty span) meets a span strictly inside it."""
    (s1, e1), (s2, e2) = first, second
    if s1 < e1 and s2 < e2:
        return s1 < e2 and s2 < e1
    if s1 == e1 and s2 == e2:
        return s1 == s2
    point, (start, end) = (s1, (s2, e2)) if s1 == e1 else (s2, (s1, e1))
    return start < point < end


def _holds(ref: str, alt: str, held: str) -> int:
    """Whether the haplotype `held`, an allele of a row over `ref`, holds the change `alt` makes there.

    A haplotype caller writes a haplotype as one row (REF `CGCCC`, ALT `TGCCT,TGCCC`), so a genotype
    can hold an allele's change inside another ALT's spelling. With both trimmed against `ref`: a
    change away from the allele's span holds REF there (`HOLDS_NOT`); of one length with REF and the
    allele, it holds the change when its bases over the allele's span are the allele's; a length change
    against another length change, or against a substitution, is another allele of the site
    (`HOLDS_NOT`). The rest (a change of another length that also replaces bases over the span, or a
    base that is not ACGT there) cannot be read: `UNREADABLE`, never `HOLDS_NOT`.
    """
    ref, alt, held = ref.upper(), alt.upper(), held.upper()
    if held == alt:
        return HOLDS
    if held == ref:
        return HOLDS_NOT
    a_start, a_end, a_bases = _change(ref, alt)
    h_start, h_end, h_bases = _change(ref, held)
    if not _overlaps((a_start, a_end), (h_start, h_end)):
        return HOLDS_NOT
    if len(held) == len(ref) == len(alt):
        window = held[a_start:a_end]
        if not set(window) <= _BASES:
            return UNREADABLE
        return HOLDS if window == alt[a_start:a_end] else HOLDS_NOT

    def pure(start: int, end: int, bases: str) -> bool:
        return start == end or not bases or len(bases) == end - start

    return HOLDS_NOT if pure(a_start, a_end, a_bases) and pure(h_start, h_end, h_bases) else UNREADABLE


def _another_allele(ref: str, alt: str, held: str) -> int:
    """A structural-variant allele's reading of a haplotype that is neither it nor REF: another allele
    of its record (Canvas's `<CN2>` beside `<CN0>`), so not a copy. A symbolic allele has no bases to
    trim, so `_holds` cannot read it."""
    return HOLDS_NOT


Reader = Callable[[str, str, str], int]
"""Reads a haplotype `held` for the allele `alt` over `ref`: `HOLDS`, `HOLDS_NOT` or `UNREADABLE`."""


def _codes(frame: pl.DataFrame, alleles: str, ref: str, alt: str, reader: Reader = _holds) -> pl.Series:
    """Per row, what each entry of the list column `alleles` says of `alt`'s change over `ref`
    (`reader`: `_holds` for a sequence allele, `_another_allele` for a structural one), null where the
    entry is null; a null or empty list gives null."""
    work = frame.select(
        pl.col(alleles).alias("_h"), pl.col(ref).str.to_uppercase().alias("_r"), pl.col(alt).str.to_uppercase().alias("_x")
    ).with_row_index("_i")
    long = (
        work.filter(pl.col("_h").list.len() > 0)
        .with_columns(pl.int_ranges(pl.col("_h").list.len(), dtype=pl.UInt32).alias("_s"))
        .explode(["_h", "_s"], empty_as_null=True)
        .with_columns(pl.col("_h").str.to_uppercase())
        .with_columns(
            pl.when(pl.col("_h") == pl.col("_x"))
            .then(pl.lit(HOLDS, dtype=pl.Int8))
            .when(pl.col("_h") == pl.col("_r"))
            .then(pl.lit(HOLDS_NOT, dtype=pl.Int8))
            .alias("_code")
        )
    )
    other = long.filter(pl.col("_h").is_not_null() & pl.col("_code").is_null())
    if other.height:
        read = other.select("_i", "_s").with_columns(
            pl.Series("_read", [reader(r, x, h) for r, x, h in other.select("_r", "_x", "_h").iter_rows()], dtype=pl.Int8)
        )
        long = long.join(read, on=["_i", "_s"], how="left").with_columns(pl.coalesce("_code", "_read").alias("_code"))
    lists = long.group_by("_i").agg(pl.col("_code").sort_by("_s"))
    return work.select("_i").join(lists, on="_i", how="left").sort("_i")["_code"]


def _read(frame: pl.DataFrame, gt: str, ref: str, alt: str, name: str, reader: Reader = _holds) -> tuple[pl.DataFrame, tuple[str, str, str, str]]:
    """`frame` with the genotype `gt` read for `alt`'s change over `ref` (`reader`): the haplotypes
    holding it, the null ones, the unreadable ones and the ploidy, as columns `_{name}_c`, `_n`, `_u`, `_p`."""
    codes = _codes(frame, gt, ref, alt, reader)
    columns = (f"_{name}_c", f"_{name}_n", f"_{name}_u", f"_{name}_p")
    read = frame.with_columns(
        codes.list.eval(pl.element() == HOLDS).list.sum().fill_null(0).cast(pl.Int32).alias(columns[0]),
        codes.list.eval(pl.element().is_null()).list.sum().fill_null(0).cast(pl.Int32).alias(columns[1]),
        codes.list.eval(pl.element() == UNREADABLE).list.sum().fill_null(0).cast(pl.Int32).alias(columns[2]),
        pl.col(gt).list.len().fill_null(0).cast(pl.Int32).alias(columns[3]),
    )
    return read, columns


def _ploidy_differs(a: tuple[str, str, str, str], b: tuple[str, str, str, str]) -> pl.Expr:
    """True where both genotypes were written with different numbers of haplotypes, False where with
    the same; null where either has none."""
    pa, pb = pl.col(a[3]), pl.col(b[3])
    return pl.when((pa > 0) & (pb > 0)).then(pa != pb)


def _agreement(a: tuple[str, str, str, str], b: tuple[str, str, str, str]) -> pl.Expr:
    """`Agreement` from two genotypes read for one allele (`_read`'s columns). A null or unreadable
    haplotype may hold the allele or not; a pair is decided only where every completion agrees.

    Of one ploidy, the genotypes agree when they hold the same copies. Of two (a haploid X call against
    a diploid one), when each holds the allele on every haplotype or each on none (F22): `1` against
    `1/1` agrees, `1` against `0/1` does not. `ploidy_differs` says the ploidies differed."""
    ca, na, ua, pa = (pl.col(name) for name in a)
    cb, nb, ub, pb = (pl.col(name) for name in b)
    open_a, open_b = na + ua, nb + ub
    a_silent = (pa == 0) | (na == pa)
    b_silent = (pb == 0) | (nb == pb)
    undecided = pl.when((ua > 0) | (ub > 0)).then(pl.lit(Agreement.UNREADABLE.value)).otherwise(pl.lit(Agreement.UNCALLED.value))
    # Some completion of the open haplotypes holds the allele on every haplotype, or on none.
    every = (ca + open_a == pa) & (cb + open_b == pb)
    none = (ca == 0) & (cb == 0)
    return (
        pl.when(a_silent | b_silent)
        .then(pl.lit(Agreement.UNCALLED.value))
        .when((pa != pb) & ~(every | none))
        .then(pl.lit(Agreement.MISMATCH.value))
        .when((pa != pb) & (open_a == 0) & (open_b == 0))
        .then(pl.lit(Agreement.MATCH.value))
        .when(pa != pb)
        .then(undecided)
        .when((open_a == 0) & (open_b == 0))
        .then(pl.when(ca == cb).then(pl.lit(Agreement.MATCH.value)).otherwise(pl.lit(Agreement.MISMATCH.value)))
        .when((ca > cb + open_b) | (cb > ca + open_a))
        .then(pl.lit(Agreement.MISMATCH.value))
        .otherwise(undecided)
    )


def _max_ploidy(frame: pl.DataFrame, *columns: str) -> int:
    lengths = [frame.select(pl.col(name).list.len().max()).item() for name in columns if name in frame.columns]
    return max([int(n) for n in lengths if n is not None], default=0)


def _two_sided(pairs: pl.DataFrame, warnings: list[str]) -> pl.DataFrame:
    """Pairs with their agreement, copies, and graded `p_same`. Each side's genotype is read for its own
    row's ALT over its own row's REF: a VRS id names the trimmed allele, which another ALT of the same
    row can hold too (`_holds`)."""
    pairs = pairs.with_row_index("_pair")
    pairs, a = _read(pairs, "a_gt_alleles", "_a_ref", "a_alt", "a")
    pairs, b = _read(pairs, "b_gt_alleles", "_b_ref", "b_alt", "b")
    graded = _graded(pairs)
    sites = pairs.join(graded, on="_pair", how="left").with_columns(
        pl.col("a_chrom").alias("chrom"),
        pl.col("a_pos").alias("pos"),
        pl.col("_a_ref").alias("ref"),
        pl.col("a_alt").alias("alt"),
        pl.col("a_variant_class").alias("variant_class"),
        pl.coalesce("_a_vid", "_b_vid").alias("vrs_allele_id"),
        _agreement(a, b).alias("agreement"),
        _ploidy_differs(a, b).alias("ploidy_differs"),
        pl.lit(None, dtype=pl.String).alias("absence"),
        pl.coalesce("basis", pl.lit(ComparisonBasis.HARD.value)).alias("basis"),
        pl.when((pl.col(a[1]) == 0) & (pl.col(a[2]) == 0)).then(pl.col(a[0])).alias("a_copies"),
        pl.when((pl.col(b[1]) == 0) & (pl.col(b[2]) == 0)).then(pl.col(b[0])).alias("b_copies"),
    )
    unreadable = sites.filter(pl.col("agreement") == Agreement.UNREADABLE.value).height
    if unreadable:
        warnings.append(
            f"{unreadable} alleles met where a genotype holds another allele of its row, of another length over the allele's span: not read for copies (agreement unreadable)"
        )
    return sites


def _distribution(pairs: pl.DataFrame, side: str, vector: str, reader: Reader = _holds) -> pl.DataFrame:
    """Per pair, the probability of each share of copies of the allele under `side`'s `vector`, as
    `_dose`: the copies on a diploid scale (a haploid copy counts 2), so two ploidies compare.

    The vector is in VCF genotype order over `source_alleles` (Decision 3): haploid when it has one
    entry per allele, diploid when it has `n(n+1)/2`; any other length, a null entry, an allele list
    that does not hold the ALT exactly once, or a total of zero gives no distribution.
    """
    alleles, alt, values = f"{side}_source_alleles", f"{side}_alt", f"{side}_{vector}"
    usable = pairs.filter(pl.col(values).is_not_null() & pl.col(alleles).is_not_null())
    if usable.is_empty():
        return pl.DataFrame(schema={"_pair": pl.UInt32, "_dose": pl.Int32, "_p": pl.Float64})
    width = min(_max_ploidy(usable, alleles), _MAX_GENOTYPE_ALLELES)
    usable = usable.with_columns(_codes(usable, alleles, f"_{side}_ref", alt, reader).alias("_carry"))
    hits = pl.sum_horizontal([(pl.col(alleles).list.get(j, null_on_oob=True) == pl.col(alt)).fill_null(False).cast(pl.Int32) for j in range(width)])
    n = pl.col(alleles).list.len().cast(pl.Int64)
    length = pl.col(values).list.len().cast(pl.Int64)
    work = (
        usable.select(
            "_pair",
            pl.col(values).alias("_v"),
            "_carry",
            hits.alias("_hits"),
            pl.when(length == n).then(1).when(length == n * (n + 1) // 2).then(2).alias("_ploidy"),
            pl.col(values).list.eval(pl.element().is_null()).list.any().alias("_holes"),
            (n <= _MAX_GENOTYPE_ALLELES).alias("_small"),
            pl.col("_carry").list.eval(pl.element().is_null() | (pl.element() == UNREADABLE)).list.any().alias("_unread"),
        )
        .filter((pl.col("_hits") == 1) & pl.col("_ploidy").is_not_null() & ~pl.col("_holes") & pl.col("_small") & ~pl.col("_unread"))
        .with_columns(pl.int_ranges(pl.col("_v").list.len(), dtype=pl.Int64).alias("_g"))
        .explode(["_v", "_g"], empty_as_null=False)
    )
    order = _diploid_order()
    work = work.join(order, on="_g", how="left").with_columns(
        pl.when(pl.col("_ploidy") == 1)
        .then(pl.col("_carry").list.get(pl.col("_g"), null_on_oob=True).cast(pl.Int32))
        .otherwise(
            pl.col("_carry").list.get(pl.col("_j"), null_on_oob=True).cast(pl.Int32) + pl.col("_carry").list.get(pl.col("_k"), null_on_oob=True).cast(pl.Int32)
        )
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
    # The share of haplotypes holding the allele, on a diploid scale: haploid `k` is diploid `2k` (F22).
    work = work.with_columns((pl.col("_copies") * (2 // pl.col("_ploidy"))).cast(pl.Int32).alias("_dose"))
    return work.group_by("_pair", "_dose").agg(pl.col("_w").sum().alias("_p"))


def _diploid_order() -> pl.DataFrame:
    """VCF's diploid genotype order: genotype `j/k` (j <= k) at index `k(k+1)/2 + j`."""
    rows = [(k * (k + 1) // 2 + j, j, k) for k in range(_MAX_GENOTYPE_ALLELES) for j in range(k + 1)]
    return pl.DataFrame(rows, schema={"_g": pl.Int64, "_j": pl.Int32, "_k": pl.Int32}, orient="row")


def _graded(pairs: pl.DataFrame, reader: Reader = _holds) -> pl.DataFrame:
    """`basis` and `p_same` for the pairs both of whose rows carry a usable vector of one kind."""
    found = []
    taken = pl.DataFrame(schema={"_pair": pl.UInt32})
    for vector, basis in (("gl", ComparisonBasis.GL_FLAT_PRIOR), ("gp", ComparisonBasis.GP)):
        da = _distribution(pairs, "a", vector, reader)
        db = _distribution(pairs, "b", vector, reader)
        both = da.select("_pair").unique().join(db.select("_pair").unique(), on="_pair").join(taken, on="_pair", how="anti")
        if both.is_empty():
            continue
        products = da.join(db, on=["_pair", "_dose"], how="inner").group_by("_pair").agg((pl.col("_p") * pl.col("_p_right")).sum().alias("p_same"))
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
    mine, own = _read(mine, f"{side}_gt_alleles", "_ref", "alt", side)
    mine, theirs = _read(mine, f"{other_side}_gt_alleles", "_ref", "alt", other_side)
    pair = (own, theirs) if side == "a" else (theirs, own)
    lone = Agreement.ONLY_A.value if side == "a" else Agreement.ONLY_B.value
    return mine.with_columns(
        pl.col("_ref").alias("ref"),
        pl.col("_vid").alias("vrs_allele_id"),
        pl.when(single).then(pl.lit(PairKey.SITE.value)).alias("key"),
        pl.when(single).then(_agreement(*pair)).otherwise(pl.lit(lone)).alias("agreement"),
        pl.when(single).then(_ploidy_differs(*pair)).alias("ploidy_differs"),
        absence.alias("absence"),
        pl.when(single).then(pl.lit(ComparisonBasis.HARD.value)).alias("basis"),
        pl.lit(None, dtype=pl.Float64).alias("p_same"),
        pl.when((pl.col(own[1]) == 0) & (pl.col(own[2]) == 0)).then(pl.col(own[0])).alias(f"{side}_copies"),
        pl.when(single & (pl.col(theirs[1]) == 0) & (pl.col(theirs[2]) == 0)).then(pl.col(theirs[0])).alias(f"{other_side}_copies"),
        pl.col("_row").alias(f"_{side}_row"),
    )


def _block_reading(mine: pl.DataFrame, other: Loaded) -> pl.DataFrame:
    """`_block`: the kind of the other facet's block over the position, where no row of it is there.
    Where several cover it (a lifted facet keeps overlapping blocks), the one `reference_support`
    answers with (`intervals.nearest_block`)."""
    probe = mine.with_row_index("_probe")
    silent = probe.filter(pl.col("_rows").is_null()).select("_probe", "chrom", pl.col("pos").cast(pl.Int64).alias("_at"))
    blocks = other.blocks.select("chrom", pl.col("pos").cast(pl.Int64).alias("_start"), pl.col("end_pos").cast(pl.Int64).alias("_end"), pl.col("kind").cast(pl.String))
    hits = covering(silent.lazy(), blocks, silent["chrom"].unique().to_list(), at="_at", start="_start", end="_end").collect()
    chosen = nearest_block(hits, "_probe", start="_start", end="_end").select("_probe", pl.col("kind").alias("_block"))
    return probe.join(chosen, on="_probe", how="left", maintain_order="left").drop("_probe")


def _declared_contigs(loaded: Loaded, facet: Mapping) -> set[str]:
    """Contigs the facet stated (its header) or wrote rows on: where its silence can mean anything."""
    declared = {canonical_chrom(entry["name"]) for entry in facet.get("header_contigs") or []}
    for frame in (loaded.calls, loaded.blocks):
        declared |= set(frame.select(pl.col("chrom").unique()).collect()["chrom"].to_list())
    return declared


_SV_INPUTS = ("end_pos", "cipos", "ciend")
"""The `calls` columns only a structural-variant allele's comparison reads."""

LENGTH_KINDS: frozenset[str] = frozenset({"DEL", "DUP", "INV", "CNV"})
"""SV types whose `end_pos` is the end of the event (docs/SCHEMA.md, `end_pos`), so that `end_pos -
pos` is its size, as `write_vcf` writes SVLEN; the copy-number alleles `<CN0>`, `<CN2>`… are too. Any
other type (`<INS>`, whose length the frame does not keep; Delly's `<TRA>`) has no size."""

_COPY_NUMBER_KIND = r"^CN\d+$"


def _sv_side(side: pl.DataFrame, name: str) -> pl.DataFrame:
    """One facet's SV alleles with what matching reads, named for side `name` as `_pair` names them:
    the type, the size, the extent as truvari's half-open `[pos - 1, end_pos)` (the anchor base and the
    event; the one place 1-based closed positions turn into it), and that extent widened by the
    allele's own CIPOS and CIEND (VCF 4.5: a missing CIEND is CIPOS)."""
    k2 = pl.col("_k").cast(pl.Int64) * 2
    kind = pl.col("alt").str.extract(r"^<([^:>]+)", 1)
    lo = pl.col("cipos").list.get(k2, null_on_oob=True)
    hi = pl.coalesce(pl.col("ciend").list.get(k2 + 1, null_on_oob=True), pl.col("cipos").list.get(k2 + 1, null_on_oob=True))
    frame = side.with_columns(
        kind.alias("_kind"),
        pl.when(kind.is_in(sorted(LENGTH_KINDS)) | kind.str.contains(_COPY_NUMBER_KIND))
        .then(pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64))
        .alias("_size"),
        (pl.col("pos").cast(pl.Int64) - 1).alias("_s"),
        pl.col("end_pos").cast(pl.Int64).alias("_e"),
    ).with_columns(
        (pl.col("_s") - lo.abs().cast(pl.Int64).fill_null(0)).alias("_ws"),
        (pl.col("_e") + hi.abs().cast(pl.Int64).fill_null(0)).alias("_we"),
    )
    return frame.rename({column: f"_{name}_{column.lstrip('_')}" if column.startswith("_") else f"{name}_{column}" for column in frame.columns})


def _sv_candidates(a: pl.DataFrame, b: pl.DataFrame, params: SvMatching) -> pl.DataFrame:
    """Every pair of `a`'s and `b`'s SV alleles on one contig whose widened extents lie within
    `refdist` of each other (truvari's `overlaps(start - refdist, end + refdist, ...)`, strict), with
    truvari's measures and whether the pair passes the rest of the thresholds (`_passes`)."""
    r = params.refdist
    pairs = a.join_where(
        b, pl.col("a_chrom") == pl.col("b_chrom"), pl.col("_a_ws") - r < pl.col("_b_we"), pl.col("_b_ws") < pl.col("_a_we") + r
    )
    sa, sb = pl.col("_a_size"), pl.col("_b_size")
    # truvari's `sizesim`: a zero size counts as one base, two zeros are the same size.
    za, zb = pl.max_horizontal(sa, pl.lit(1)), pl.max_horizontal(sb, pl.lit(1))
    similarity = pl.when((sa == 0) & (sb == 0)).then(pl.lit(1.0)).otherwise(pl.min_horizontal(za, zb) / pl.max_horizontal(za, zb))
    shared = pl.min_horizontal("_a_e", "_b_e") - pl.max_horizontal("_a_s", "_b_s")
    longest = pl.max_horizontal(pl.col("_a_e") - pl.col("_a_s"), pl.col("_b_e") - pl.col("_b_s"))
    stated = (pl.col("_a_s") - r < pl.col("_b_e")) & (pl.col("_b_s") < pl.col("_a_e") + r)
    pairs = pairs.with_columns(
        pl.when(sa.is_not_null() & sb.is_not_null()).then(similarity).alias("size_similarity"),
        pl.when(shared > 0).then(shared / longest).otherwise(pl.lit(0.0)).alias("reciprocal_overlap"),
        (~stated).alias("ci_needed"),
    )
    same_type = pl.lit(True) if params.typeignore else pl.col("_a_kind") == pl.col("_b_kind")
    # A size withheld on either side (an `<INS>`) leaves the size threshold unchecked, never passed.
    sized = pl.col("size_similarity").is_null() | (pl.col("size_similarity") >= params.pctsize)
    return pairs.with_columns((same_type & sized & (pl.col("reciprocal_overlap") >= params.pctovl)).alias("_passes"))


def _pick_single(candidates: pl.DataFrame) -> pl.DataFrame:
    """truvari's `--pick single`: the passing pairs best first (score, then the nearer start and end),
    each allele taken at most once. truvari's score is `(pctseq + pctsize + pctovl) / 3` with `pctseq`
    0 for a symbolic allele; a size withheld counts 0. Ties, which truvari leaves to its sort, fall to
    each side's row order, so the pick is deterministic."""
    passing = candidates.filter(pl.col("_passes")).with_columns(
        (pl.col("size_similarity").fill_null(0.0) + pl.col("reciprocal_overlap")).alias("_score"),
        (pl.col("a_pos").cast(pl.Int64) - pl.col("b_pos").cast(pl.Int64)).abs().alias("_start"),
        (pl.col("_a_e") - pl.col("_b_e")).abs().alias("_end"),
    )
    order = ["_score", "_start", "_end", "_a_row", "_a_k", "_b_row", "_b_k"]
    passing = passing.sort(order, descending=[True, False, False, False, False, False, False])
    return _take_once(passing, ("_a_row", "_a_k"), ("_b_row", "_b_k")).drop("_score", "_start", "_end")


def _take_once(ordered: pl.DataFrame, a_key: tuple[str, ...], b_key: tuple[str, ...]) -> pl.DataFrame:
    """The rows of `ordered` (best first) taken greedily, each `a_key` and each `b_key` at most once:
    truvari's `pick_single_matches`."""
    used_a: set[tuple] = set()
    used_b: set[tuple] = set()
    keep = []
    width = len(a_key)
    for index, row in enumerate(ordered.select(*a_key, *b_key).iter_rows()):
        ka, kb = row[:width], row[width:]
        if ka in used_a or kb in used_b:
            continue
        used_a.add(ka)
        used_b.add(kb)
        keep.append(index)
    return ordered[keep]


def _sv_sites(
    sv_a: pl.DataFrame, sv_b: pl.DataFrame, a: tuple[Loaded, Mapping], b: tuple[Loaded, Mapping], params: SvMatching, warnings: list[str]
) -> list[pl.DataFrame]:
    """The SV alleles' rows of `sites`: pairs met by extent, then each side's alleles left alone."""
    side_a, side_b = _sv_side(sv_a, "a"), _sv_side(sv_b, "b")
    candidates = _sv_candidates(side_a, side_b, params)
    pairs = _pick_single(candidates)
    pieces = []
    if pairs.height:
        unsized = pairs.filter(pl.col("size_similarity").is_null()).height
        if unsized:
            warnings.append(
                f"{unsized} SV pairs met without a size: the frame keeps no length for their type (an <INS> has no SVLEN), so size similarity was not checked"
            )
        pieces.append(_sv_pair_rows(pairs, PairKey.EXTENT))
    for side, other, own in (("a", b, side_a), ("b", a, side_b)):
        taken = pairs.select(f"_{side}_row", f"_{side}_k")
        rest = own.join(taken, on=[f"_{side}_row", f"_{side}_k"], how="anti")
        if rest.height:
            near = candidates.select(f"_{side}_row", f"_{side}_k").unique()
            pieces.append(_sv_one_sided(rest, near, other, side, warnings))
    return pieces


def _sv_pair_rows(pairs: pl.DataFrame, key: PairKey) -> pl.DataFrame:
    """Two-sided structural-variant pairs (side columns as `_sv_side` names them) with their
    agreement, copies and graded `p_same`, each side's genotype read for its own ALT."""
    pairs = pairs.with_row_index("_pair")
    pairs, ra = _read(pairs, "a_gt_alleles", "_a_ref", "a_alt", "a", _another_allele)
    pairs, rb = _read(pairs, "b_gt_alleles", "_b_ref", "b_alt", "b", _another_allele)
    graded = _graded(pairs, _another_allele)
    return pairs.join(graded, on="_pair", how="left").with_columns(
        pl.col("a_chrom").alias("chrom"),
        pl.col("a_pos").alias("pos"),
        pl.col("_a_ref").alias("ref"),
        pl.col("a_alt").alias("alt"),
        pl.col("a_variant_class").alias("variant_class"),
        pl.lit(key.value).alias("key"),
        _agreement(ra, rb).alias("agreement"),
        _ploidy_differs(ra, rb).alias("ploidy_differs"),
        pl.coalesce("basis", pl.lit(ComparisonBasis.HARD.value)).alias("basis"),
        pl.when((pl.col(ra[1]) == 0) & (pl.col(ra[2]) == 0)).then(pl.col(ra[0])).alias("a_copies"),
        pl.when((pl.col(rb[1]) == 0) & (pl.col(rb[2]) == 0)).then(pl.col(rb[0])).alias("b_copies"),
        pl.col("_a_kind").alias("sv_kind"),
    )


def _sv_one_sided(
    rest: pl.DataFrame, near: pl.DataFrame, other: tuple[Loaded, Mapping], side: str, warnings: list[str], *, what: str = "SV alleles", mate: str | None = None
) -> pl.DataFrame:
    """SV alleles of `side` the other facet matched to none of its own, with what its silence means:
    SVs of its own near it (`near`'s rows: within `refdist` for a symbolic allele, within a breakend's
    position window for a breakend; `other_alleles`), else its absence rule for the `sv` class
    (`absence_is_reference`), on a contig it declared (and, for a breakend, the `mate` column's contig
    too). Its blocks are not read: a block states the bases a small-variant caller saw, and an SV's
    anchor base lies outside the event."""
    other_loaded, other_facet = other
    other_side = "b" if side == "a" else "a"
    declared = sorted(_declared_contigs(other_loaded, other_facet))
    reference = absence_is_reference(other_facet, VariantClass.SV) is True
    marked = rest.join(near.with_columns(pl.lit(True).alias("_near")), on=[f"_{side}_row", f"_{side}_k"], how="left")
    seen = pl.col(f"{side}_chrom").is_in(declared) & (pl.col(mate).is_in(declared) if mate else pl.lit(True))
    absence = (
        pl.when(pl.col("_near").fill_null(False))
        .then(pl.lit(AbsenceReading.OTHER_ALLELES.value))
        .when(pl.lit(reference) & seen)
        .then(pl.lit(AbsenceReading.REFERENCE.value))
        .otherwise(pl.lit(AbsenceReading.UNKNOWN.value))
    )
    marked = marked.with_columns(absence.alias("absence"))
    unknown = marked.filter(pl.col("absence") == AbsenceReading.UNKNOWN.value).height
    if unknown:
        why = "its absence is read as reference for SVs only on a whole-genome callset that called SVs and was not filtered"
        warnings.append(f"{unknown} {what} of facet {side} meet no SV of facet {other_side} and its silence there is unknown: {why}")
    marked, own = _read(marked, f"{side}_gt_alleles", f"_{side}_ref", f"{side}_alt", side, _another_allele)
    return marked.with_columns(
        pl.col(f"{side}_chrom").alias("chrom"),
        pl.col(f"{side}_pos").alias("pos"),
        pl.col(f"_{side}_ref").alias("ref"),
        pl.col(f"{side}_alt").alias("alt"),
        pl.col(f"{side}_variant_class").alias("variant_class"),
        pl.lit(Agreement.ONLY_A.value if side == "a" else Agreement.ONLY_B.value).alias("agreement"),
        pl.when((pl.col(own[1]) == 0) & (pl.col(own[2]) == 0)).then(pl.col(own[0])).alias(f"{side}_copies"),
        pl.col(f"_{side}_kind").alias("sv_kind"),
        pl.lit(None, dtype=pl.String).alias("b_alt"),
    )


_BND_INPUTS = ("_form", "_mate_chrom", "_mate_pos")
"""What `_alleles` parses from a breakend's ALT (`_breakend_parts`), null on any other allele."""

BREAKEND_FORMS: tuple[str, ...] = ("t[p[", "t]p]", "]p]t", "[p[t")
"""VCF 4.2-4.5's four breakend ALT forms (`t` the bases at POS, `p` the mate's `contig:position`), the
values of `_form`. They are truvari's four `(direction, strand)` pairs (`bnd_direction_strand`):
`t[p[` left/direct, `[p[t` left/complement, `t]p]` right/direct, `]p]t` right/complement."""

_SINGLE_BREAKEND = r"^\.[A-Za-z]+$|^[A-Za-z]+\.$"
"""A single breakend (`.A`, `A.`): a join to sequence the file does not place."""

_BASES_FIRST = r"^[A-Za-z]+([\[\]])([^\[\]]+):(\d+)([\[\]])$"
_BASES_LAST = r"^([\[\]])([^\[\]]+):(\d+)([\[\]])[A-Za-z]+$"
"""A breakend ALT with its bases before (`t[p[`, `t]p]`) or after (`]p]t`, `[p[t`) the mate. The
contig runs to the last `:` before the digits, so a name holding `:` (`HLA-A*01:01:01:01`) reads;
truvari's `[^\\[\\]:]+` would refuse it."""


def _breakend_parts(alt: pl.Expr) -> list[pl.Expr]:
    """`_form`, `_mate_chrom` (canonical) and `_mate_pos` of a breakend ALT, all null where the ALT is
    none of the four forms: brackets that do not pair, or no `contig:position` between them."""
    first, last = alt.str.extract_groups(_BASES_FIRST), alt.str.extract_groups(_BASES_LAST)
    f_open, f_close = first.struct.field("1"), first.struct.field("4")
    l_open, l_close = last.struct.field("1"), last.struct.field("4")
    position = pl.coalesce(first.struct.field("3"), last.struct.field("3")).cast(pl.Int64, strict=False)
    form = (
        pl.when(position.is_null())
        .then(None)
        .when(f_open.is_not_null() & (f_open == f_close))
        .then(pl.when(f_open == "[").then(pl.lit("t[p[")).otherwise(pl.lit("t]p]")))
        .when(l_open.is_not_null() & (l_open == l_close))
        .then(pl.when(l_open == "[").then(pl.lit("[p[t")).otherwise(pl.lit("]p]t")))
    )
    contig = pl.coalesce(first.struct.field("2"), last.struct.field("2")).replace(CONTIG_ALIASES)
    return [form.alias("_form"), pl.when(form.is_not_null()).then(contig).alias("_mate_chrom"), pl.when(form.is_not_null()).then(position).alias("_mate_pos")]


def _not_after(left: list[pl.Expr], right: list[pl.Expr]) -> pl.Expr:
    """Whether the tuple `left` sorts before `right` or equals it, element by element."""
    result = pl.lit(True)
    for lhs, rhs in reversed(list(zip(left, right, strict=True))):
        result = (lhs < rhs) | ((lhs == rhs) & result)
    return result


def _bnd_side(side: pl.DataFrame, name: str, bnddist: int, warnings: list[str]) -> pl.DataFrame:
    """One facet's breakend records with the adjacency each states, named for side `name` as
    `_sv_side` names its columns.

    An adjacency is two ends, each a `(contig, position, side)`: the side of the position whose
    sequence takes part in the join. A record's own end keeps the sequence left of POS in `t[p[` and
    `t]p]`, right of it in `]p]t` and `[p[t`; its mate's end keeps the sequence right of `p` in `t[p[`
    and `[p[t`, left of it in the other two. A mate record states the same two ends from the other
    side (VCF's `t[p[` beside `]p]t`, `t]p]` beside `t]p]`, `[p[t` beside `[p[t`), so records stating
    one pair of ends are one adjacency (`_adj`), whether or not the file names its mate: the frame
    keeps no MATEID or EVENT (`dropped_fields`), and neither is needed to find it. Mates whose ends
    differ by no more than their position intervals are one adjacency too (`_fold_mates`). `_first`
    marks a record standing at its adjacency's first end in canonical order (contig rank, contig,
    position, side). An EVENT of several adjacencies (an inversion's two) is several rows.

    The windows are truvari's `bnd_match` `bounds`: `pos +/- bnddist`, widened by `|CIPOS|` for the
    record's own position and `|CIEND|` for its mate's, both read at index 1 (the upper offset):
    `bounds` sets `key = 'CI' + key` before `idx = 0 if key == 'POS' else 1`, so the test never
    holds. No CIEND widens the mate by nothing; truvari does not fall back to CIPOS there."""
    hi = pl.col("_k").cast(pl.Int64) * 2 + 1
    ci = pl.col("cipos").list.get(hi, null_on_oob=True).cast(pl.Int64).abs().fill_null(0)
    ci_mate = pl.col("ciend").list.get(hi, null_on_oob=True).cast(pl.Int64).abs().fill_null(0)
    rank = {column: pl.col(column).replace_strict(CONTIG_RANK, default=len(CONTIG_RANK), return_dtype=pl.UInt16) for column in ("chrom", "_mate_chrom")}
    frame = side.with_columns(
        pl.col("pos").cast(pl.Int64).alias("_p"),
        pl.col("_mate_pos").alias("_mp"),
        pl.when(pl.col("_form").is_in(["t[p[", "t]p]"])).then(pl.lit("L")).otherwise(pl.lit("R")).alias("_side"),
        pl.when(pl.col("_form").is_in(["t[p[", "[p[t"])).then(pl.lit("R")).otherwise(pl.lit("L")).alias("_mate_side"),
        ci.alias("_ci"),
        ci_mate.alias("_ci_mate"),
        pl.lit("BND").alias("_kind"),
    )
    own = [rank["chrom"], pl.col("chrom"), pl.col("_p"), pl.col("_side")]
    mate = [rank["_mate_chrom"], pl.col("_mate_chrom"), pl.col("_mp"), pl.col("_mate_side")]
    frame = frame.with_columns(
        _not_after(own, mate).alias("_first"),
        ((pl.col("chrom") == pl.col("_mate_chrom")) & (pl.col("_p") == pl.col("_mp")) & (pl.col("_side") == pl.col("_mate_side"))).alias("_self"),
    )
    ends = {"_e1c": ("chrom", "_mate_chrom"), "_e1p": ("_p", "_mp"), "_e1s": ("_side", "_mate_side")}
    frame = frame.with_columns(
        *[pl.when(pl.col("_first")).then(pl.col(mine)).otherwise(pl.col(theirs)).alias(e) for e, (mine, theirs) in ends.items()],
        *[pl.when(pl.col("_first")).then(pl.col(theirs)).otherwise(pl.col(mine)).alias(e.replace("1", "2")) for e, (mine, theirs) in ends.items()],
    )
    keys = ["_e1c", "_e1p", "_e1s", "_e2c", "_e2p", "_e2s"]
    ids = frame.select(keys).unique().sort(keys).with_row_index("_adj")
    frame = _fold_mates(frame.join(ids, on=keys, how="left"), name, warnings)
    frame = frame.sort(["_adj", "_first", "_row", "_k"], descending=[False, True, False, False])
    frame, read = _read(frame, "gt_alleles", "_ref", "alt", "m", _another_allele)
    frame = frame.with_columns(
        (pl.int_range(pl.len()).over("_adj") == 0).alias("_rep"),
        # Records beyond one per end; an adjacency joining an end to itself (a foldback, `t]p]` at its
        # own position) has two mates standing at one end. Two foldbacks folded within their
        # intervals stand at two ends already, hence the floor at 0.
        pl.max_horizontal(
            pl.len().over("_adj").cast(pl.Int64)
            - pl.struct("chrom", "_p", "_side").n_unique().over("_adj").cast(pl.Int64)
            - (pl.col("_self").all().over("_adj") & (pl.len().over("_adj") > 1)).cast(pl.Int64),
            pl.lit(0, dtype=pl.Int64),
        ).alias("_repeats"),
        (pl.struct(*read).n_unique().over("_adj") > 1).alias("_split"),
    ).drop(*read)
    repeats = int(frame.filter(pl.col("_rep"))["_repeats"].sum())
    if repeats:
        warnings.append(f"{repeats} breakend records of facet {name} repeat an end of an adjacency another of its records states; each adjacency is counted once")
    split = frame.filter(pl.col("_rep") & pl.col("_split")).height
    if split:
        warnings.append(
            f"{split} breakend adjacencies of facet {name} whose two mate records state different genotypes; the record that met, else the one at the first end, is read"
        )
    frame = frame.drop("_repeats", "_split").with_columns(
        (pl.col("_p") - bnddist - pl.col("_ci")).alias("_lo"),
        (pl.col("_p") + bnddist + pl.col("_ci")).alias("_hi"),
        (pl.col("_mp") - bnddist - pl.col("_ci_mate")).alias("_mlo"),
        (pl.col("_mp") + bnddist + pl.col("_ci_mate")).alias("_mhi"),
        pl.concat_str("_mate_chrom", "_form", separator="\t").alias("_join"),
    )
    return frame.rename({column: f"_{name}_{column.lstrip('_')}" if column.startswith("_") else f"{name}_{column}" for column in frame.columns})


def _fold_mates(frame: pl.DataFrame, name: str, warnings: list[str]) -> pl.DataFrame:
    """`frame` (`_bnd_side`'s, with `_adj` from the exact ends) with two adjacencies made one where
    their records are mates stated within each other's position interval: each one's mate contig and
    side are the other's own, and its mate position lies in the other's `[POS + CIPOS[0], POS +
    CIPOS[1]]` (widened to hold POS). Manta writes a mate's POS at one end of the junction's
    microhomology, which CIPOS spans, and the ALT pointing at it from the other record at the other
    end, so such mates state their ends a few bases apart. Only a one-to-one pairing folds; an
    adjacency that could be the mate of several stays alone and is counted in a warning."""
    k2 = pl.col("_k").cast(pl.Int64) * 2
    lo = pl.col("cipos").list.get(k2, null_on_oob=True).cast(pl.Int64).fill_null(0)
    hi = pl.col("cipos").list.get(k2 + 1, null_on_oob=True).cast(pl.Int64).fill_null(0)
    records = frame.select(
        "_adj",
        "chrom",
        "_side",
        "_mate_chrom",
        "_mate_side",
        "_mp",
        (pl.col("_p") + pl.min_horizontal(lo, pl.lit(0))).alias("_from"),
        (pl.col("_p") + pl.max_horizontal(hi, pl.lit(0))).alias("_to"),
    )
    mates = records.rename({column: f"{column}_m" for column in records.columns})
    links = (
        records.join(mates, left_on=["_mate_chrom", "_mate_side", "chrom", "_side"], right_on=["chrom_m", "_side_m", "_mate_chrom_m", "_mate_side_m"])
        .filter(
            (pl.col("_adj") < pl.col("_adj_m"))
            & pl.col("_mp").is_between(pl.col("_from_m"), pl.col("_to_m"))
            & pl.col("_mp_m").is_between(pl.col("_from"), pl.col("_to"))
        )
        .select("_adj", "_adj_m")
        .unique()
    )
    if links.is_empty():
        return frame
    degree = pl.concat([links.select(pl.col("_adj").alias("_x")), links.select(pl.col("_adj_m").alias("_x"))]).group_by("_x").len()
    single = degree.filter(pl.col("len") == 1)["_x"].to_list()
    several = degree.filter(pl.col("len") > 1).height
    if several:
        warnings.append(
            f"{several} breakend adjacencies of facet {name} could each be the mate of several others within their position intervals; not folded, each counted alone"
        )
    folded = links.filter(pl.col("_adj").is_in(single) & pl.col("_adj_m").is_in(single))
    into = folded.select(pl.col("_adj_m").alias("_adj"), pl.col("_adj").alias("_into"))
    return frame.join(into, on="_adj", how="left").with_columns(pl.coalesce("_into", "_adj").alias("_adj")).drop("_into")


def _bnd_sites(
    bnd_a: pl.DataFrame, bnd_b: pl.DataFrame, a: tuple[Loaded, Mapping], b: tuple[Loaded, Mapping], params: SvMatching, warnings: list[str]
) -> list[pl.DataFrame]:
    """The breakend adjacencies' rows of `sites`: pairs met under truvari's `bnd_match` (truvari 5.5,
    `truvari/variant_record.py`, `VariantRecord.bnd_match`), then each side's adjacencies left alone.

    Two records meet when they lie on one contig with their position windows overlapping strictly
    (truvari's `overlaps`: `max(starts) < min(ends)`), their mates on one contig (parsed from ALT,
    canonical names here where truvari compares them as written) with those windows overlapping too,
    and their ALTs in the same form (truvari's `bnd_direction_strand`). Two adjacencies meet when any
    of their records do (truvari compares records; a pair of mates is two of them); the record pair
    shown is the best by truvari's score, `max(0, 1 - mean(|start|, |end| distance) / (2 * bnddist))`
    times 100 (0 at `bnddist` 0), then the nearer start and end, as its `MatchResult` orders them.
    Each adjacency is then taken at most once (`--pick single`), best first, ties to each side's
    record order."""
    d = params.bnddist
    side_a, side_b = _bnd_side(bnd_a, "a", d, warnings), _bnd_side(bnd_b, "b", d, warnings)
    near = side_a.join_where(side_b, pl.col("a_chrom") == pl.col("b_chrom"), pl.col("_a_lo") < pl.col("_b_hi"), pl.col("_b_lo") < pl.col("_a_hi"))
    near = near.filter(pl.max_horizontal("_a_lo", "_b_lo") < pl.min_horizontal("_a_hi", "_b_hi"))
    met = near.filter((pl.col("_a_join") == pl.col("_b_join")) & (pl.max_horizontal("_a_mlo", "_b_mlo") < pl.min_horizontal("_a_mhi", "_b_mhi")))
    start, end = (pl.col("_a_p") - pl.col("_b_p")).abs(), (pl.col("_a_mp") - pl.col("_b_mp")).abs()
    score = pl.max_horizontal(pl.lit(0.0), (1 - ((start + end) / 2) / (2 * d)) * 100) if d > 0 else pl.lit(0.0)
    met = met.with_columns(
        score.alias("_score"),
        start.alias("_start"),
        end.alias("_end"),
        ((start < 2 * d) & (end < 2 * d)).alias("_stated"),
    )
    order = ["_score", "_start", "_end", "_a_first", "_b_first", "_a_row", "_a_k", "_b_row", "_b_k"]
    met = met.sort(order, descending=[True, False, False, True, True, False, False, False, False])
    best = met.with_columns(pl.col("_stated").any().over("_a_adj", "_b_adj").alias("_any_stated")).unique(["_a_adj", "_b_adj"], keep="first", maintain_order=True)
    pairs = _take_once(best, ("_a_adj",), ("_b_adj",)).with_columns((~pl.col("_any_stated")).alias("ci_needed"))
    pieces = [_sv_pair_rows(pairs, PairKey.BREAKEND)] if pairs.height else []
    for side, other, own in (("a", b, side_a), ("b", a, side_b)):
        adj = f"_{side}_adj"
        rest = own.filter(pl.col(f"_{side}_rep")).join(pairs.select(adj), on=adj, how="anti")
        if rest.height:
            close = own.filter(pl.col(f"_{side}_rep")).join(near.select(adj).unique(), on=adj, how="semi").select(f"_{side}_row", f"_{side}_k")
            pieces.append(_sv_one_sided(rest, close, other, side, warnings, what="breakend adjacencies", mate=f"_{side}_mate_chrom"))
    return pieces


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
