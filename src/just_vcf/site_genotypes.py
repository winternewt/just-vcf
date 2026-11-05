"""Reading a facet's genotype at each site of a panel: an alternate-allele count of 0, 1 or 2, or unknown.

Shared by population inference (`just_vcf.ancestry`, RM5) and relatedness (`just_vcf.relatedness`,
RM274). Polars only, so it sits in the core; `just_vcf.ancestry.sites` turns its answer into the numpy
arrays the projections need.

just-prs read `chrom, pos, ref, alt, GT` strings and took an absent site as homozygous reference
unless told otherwise. Here the facet itself says what an absence means (`facets.absence_means`,
charter invariant 11), and a genotype is read from its literal haplotype alleles (`calls.gt_alleles`),
so an array export, which writes no REF, is read the same way as a VCF.

Per panel site, in this order:

1. **Rows at the site.** Every `calls` row of the facet at the site's `(chrom, pos)` is read against
   the panel's two alleles. A row is usable when it has no `split_signal` (in a pre-split file a `0`
   can mean another ALT; `rejoin_split` re-joins such rows first), exactly two haplotypes, and both
   of them one of the panel's two alleles (a null haplotype is neither), unless its source REF states
   no base and an ALT of it is a sequence allele (`ref_check = unstated`, `schema.unstated_sequence_alt`,
   F23): that ALT was spelled against the `N`, so `N` → `G`, `1/1`, names no base beside `G` and is
   not read. Its count is the number of haplotypes carrying the panel's effect allele. Exactly one
   usable row gives the site's count. Rows that are there but none usable, or two usable rows, make
   the site **unusable**: unknown, never 0. A no-call is the source saying it looked and could not
   tell. Where every row at the site is one spelled against an unstated REF, the site is
   **ref_unstated** instead: unknown too, and still rows at the position, so step 2 does not apply.
   Before `resolve_ref` a row's `ref_check` is `not_run` and it is read as written.
2. **No row.** What that means is the facet's `absence_means`: `weak_reference` (a whole-genome callset
   without blocks) reads it as homozygous for the panel's reference allele, count 0; `blocks` (a
   gVCF) does so only inside a `ref` block; `reference` (CGI's VCF, whose producer called what it
   omits reference) does so outside every block, a no-call block being no omission; `unknown`, or
   null, leaves it unknown. Absence is read as reference only on the autosomes, where homozygous
   means two copies. `weak_reference` and `reference` hold only where the facet's silence speaks for
   SNVs, a panel's sites (`callset.absence_is_reference(facet, "snv", force=True)`): a facet that wrote
   no SNV, or a genome cut to its indels (`content_filter = indel_only`), leaves it unknown. `force`
   is passed because a panel read has always read an SNV-only half's silence through its filter's
   hold; a class the filter cut is not lifted by it.

The panel's first allele is taken as the reference allele; a panel whose first allele is not the
assembly's REF would turn every absence into the wrong genotype. Panel loaders state where their
first allele comes from.
"""

from dataclasses import dataclass
from enum import Enum

import polars as pl

from just_vcf.callset import absence_is_reference
from just_vcf.loading import Loaded
from just_vcf.schema import CONTIG_ALIASES, AbsenceMeans, BlockKind, VariantClass, unstated_sequence_alt

AUTOSOMES: frozenset[str] = frozenset(str(n) for n in range(1, 23))
"""Canonical `chrom` names on which an absence can be read as homozygous reference."""

PANEL_SITE_COLUMNS = ("chrom", "pos", "ref_allele", "effect_allele")
"""What a panel hands to `read_site_genotypes`, one row per site in the panel's own order."""


class SiteBasis(str, Enum):
    """Where a site's genotype came from, or why there is none."""

    OBSERVED = "observed"
    """One usable row of the facet at the site."""
    ABSENT_AS_REFERENCE = "absent_as_reference"
    """No row, and the facet's `absence_means` reads the silence as homozygous reference."""
    UNUSABLE = "unusable"
    """Rows at the site, none usable or more than one usable: unknown."""
    REF_UNSTATED = "ref_unstated"
    """Rows at the site, every one with a source REF that states no base and a sequence ALT spelled
    against it (`ref_check = unstated`, F23): unknown. Still rows at the position, so the facet's
    silence is not read there."""
    ABSENT = "absent"
    """No row, and the facet's silence says nothing: unknown."""


@dataclass(frozen=True)
class Facet:
    """The facet facts reading a panel needs."""

    facet_id: str
    assembly_label: str | None
    absence_means: str | None
    silent_snv_is_reference: bool
    """`absence_is_reference(facet, "snv", force=True) is True`: the facet's silence speaks for SNVs, so
    `weak_reference` and `reference` may read an absent panel site as reference (module docstring)."""


_ABSENCE_COLUMNS = ("facet_id", "assembly_label", "absence_means", "callset_type", "variant_classes", "content_filter")


def facets_of(loaded: Loaded) -> list[Facet]:
    """Every facet of `loaded`, in its `facets` order."""
    rows = loaded.facets.select(_ABSENCE_COLUMNS).collect().iter_rows(named=True)
    return [
        Facet(
            row["facet_id"],
            row["assembly_label"],
            None if row["absence_means"] is None else str(row["absence_means"]),
            absence_is_reference(row, VariantClass.SNV, force=True) is True,
        )
        for row in rows
    ]


def facet_of(loaded: Loaded, facet_id: str | None = None) -> Facet:
    """The one facet of `loaded` to read, or the one named `facet_id`."""
    facets = facets_of(loaded)
    if facet_id is not None:
        facets = [facet for facet in facets if facet.facet_id == facet_id]
        if not facets:
            raise ValueError(f"no facet {facet_id!r} in this Loaded")
    if len(facets) != 1:
        raise ValueError(f"{len(facets)} facets in this Loaded; name the one to read with facet_id")
    return facets[0]


def panel_sites(frame: pl.DataFrame) -> pl.DataFrame:
    """A panel's site table in the shape `read_site_genotypes` reads: canonical `chrom`, `pos` as `UInt32`.

    `frame` has the columns of `PANEL_SITE_COLUMNS`. Contig names, numbers included (Privé writes
    `chr` as an integer), are mapped by the schema's alias table (Decision 7), never by stripping a
    `chr` prefix.
    """
    return frame.select(
        pl.col("chrom").cast(pl.String).replace(CONTIG_ALIASES),
        pl.col("pos").cast(pl.UInt32),
        pl.col("ref_allele").cast(pl.String),
        pl.col("effect_allele").cast(pl.String),
    )


def read_site_genotypes(loaded: Loaded, facet: Facet, sites: pl.DataFrame) -> pl.DataFrame:
    """The facet's effect-allele count at each panel site (module docstring has the rules).

    One row per site of `sites`, in its order: `site` (the row index in `sites`), `count` (`UInt8`,
    null when unknown) and `basis` (a `SiteBasis` value).
    """
    sites = sites.with_row_index("_site")
    keys = sites.lazy().select("_site", "chrom", "pos", "ref_allele", "effect_allele")
    rows = (
        loaded.calls.filter(pl.col("facet_id") == facet.facet_id)
        .select("chrom", "pos", "gt_alleles", "split_signal", _unstated=unstated_sequence_alt())
        .join(keys.select("chrom", "pos").unique(), on=["chrom", "pos"], how="semi")
        .with_row_index("_row")
    )
    pairs = keys.join(rows, on=["chrom", "pos"], how="inner").with_columns(
        _ploidy=pl.col("gt_alleles").list.len(),
        _split=pl.col("split_signal").is_not_null(),
    )
    allele = pl.col("_allele")
    per_row = (
        pairs.select("_site", "_row", "_ploidy", "_split", "_unstated", "ref_allele", "effect_allele", pl.col("gt_alleles").alias("_allele"))
        .explode("_allele", empty_as_null=True)
        .group_by("_site", "_row")
        .agg(
            pl.col("_ploidy").first(),
            pl.col("_split").first(),
            pl.col("_unstated").first(),
            _in_pair=((allele == pl.col("ref_allele")) | (allele == pl.col("effect_allele"))).fill_null(False).sum(),
            _effect=(allele == pl.col("effect_allele")).fill_null(False).sum(),
        )
        .with_columns(_usable=(~pl.col("_split")) & (~pl.col("_unstated")) & (pl.col("_ploidy") == 2) & (pl.col("_in_pair") == 2))
    )
    per_site = per_row.group_by("_site").agg(
        _n_usable=pl.col("_usable").sum(),
        _all_unstated=pl.col("_unstated").all(),
        _count=pl.col("_effect").filter(pl.col("_usable")).first(),
    )
    reference = _absent_as_reference(loaded, facet, sites)
    return (
        sites.lazy()
        .select("_site")
        .join(per_site, on="_site", how="left")
        .join(reference, on="_site", how="left")
        .select(
            pl.col("_site").alias("site"),
            pl.when(pl.col("_n_usable") == 1)
            .then(pl.col("_count").cast(pl.UInt8))
            .when(pl.col("_n_usable").is_null() & pl.col("_reference").fill_null(False))
            .then(pl.lit(0, dtype=pl.UInt8))
            .alias("count"),
            pl.when(pl.col("_n_usable") == 1)
            .then(pl.lit(SiteBasis.OBSERVED.value))
            .when((pl.col("_n_usable") == 0) & pl.col("_all_unstated"))
            .then(pl.lit(SiteBasis.REF_UNSTATED.value))
            .when(pl.col("_n_usable").is_not_null())
            .then(pl.lit(SiteBasis.UNUSABLE.value))
            .when(pl.col("_reference").fill_null(False))
            .then(pl.lit(SiteBasis.ABSENT_AS_REFERENCE.value))
            .otherwise(pl.lit(SiteBasis.ABSENT.value))
            .alias("basis"),
        )
        .sort("site")
        .collect()
    )


def _absent_as_reference(loaded: Loaded, facet: Facet, sites: pl.DataFrame) -> pl.LazyFrame:
    """`_site, _reference`: whether the facet's `absence_means` reads a silence at the site as
    reference. Only sites with no row take it; the caller applies that."""
    autosomal = pl.col("chrom").is_in(list(AUTOSOMES))
    frame = sites.lazy().select("_site", "chrom", "pos")
    silence = (AbsenceMeans.WEAK_REFERENCE.value, AbsenceMeans.REFERENCE.value)
    if facet.absence_means in silence and not facet.silent_snv_is_reference:
        # No SNV written, or the SNVs cut away (an indel-only genome half): the silence says nothing here.
        return frame.select("_site", _reference=pl.lit(False))
    if facet.absence_means == AbsenceMeans.WEAK_REFERENCE.value:
        return frame.select("_site", _reference=autosomal)
    if facet.absence_means == AbsenceMeans.BLOCKS.value:
        covered = _inside_block(loaded, facet, sites, BlockKind.REF)
        return frame.join(covered, on="_site", how="left").select("_site", _reference=autosomal & pl.col("_covered").fill_null(False))
    if facet.absence_means == AbsenceMeans.REFERENCE.value:
        # What the file omits was called reference; a site a block covers was not omitted (a no-call).
        covered = _inside_block(loaded, facet, sites, None)
        return frame.join(covered, on="_site", how="left").select("_site", _reference=autosomal & ~pl.col("_covered").fill_null(False))
    return frame.select("_site", _reference=pl.lit(False))


def _inside_block(loaded: Loaded, facet: Facet, sites: pl.DataFrame, kind: BlockKind | None) -> pl.LazyFrame:
    """`_site, _covered` for the sites a block of the facet covers, of `kind` or of any kind (blocks are
    1-based and inclusive). Each block start carries the furthest end of the blocks up to it, so a site
    inside a long block that a shorter one starts within is still covered."""
    chroms = sites["chrom"].unique().to_list()
    blocks = (
        loaded.blocks.filter(
            (pl.col("facet_id") == facet.facet_id)
            & (pl.lit(True) if kind is None else pl.col("kind") == kind.value)
            & pl.col("chrom").is_in(chroms)
        )
        .select("chrom", pl.col("pos").alias("_block_start"), pl.col("end_pos").alias("_block_end"))
        .sort("chrom", "_block_start")
        .with_columns(pl.col("_block_end").cum_max().over("chrom"))
    )
    probe = sites.lazy().select("_site", "chrom", "pos").sort("chrom", "pos")
    return (
        probe.join_asof(blocks, left_on="pos", right_on="_block_start", by="chrom", strategy="backward", check_sortedness=False)
        .filter(pl.col("_block_end") >= pl.col("pos"))
        .select("_site", _covered=pl.lit(True))
    )
