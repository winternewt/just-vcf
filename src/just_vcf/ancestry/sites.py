"""Reading a facet's genotypes at a panel's sites: one effect-allele dosage per site, or missing.

This is where the port from just-prs changes most. just-prs read `chrom, pos, ref, alt, GT` strings
and took an absent site as homozygous reference unless told otherwise. Here the facet itself says what
an absence means (`facets.absence_means`, charter invariant 11), and a genotype is read from its
literal haplotype alleles (`calls.gt_alleles`), so an array export, which writes no REF, is read the
same way as a VCF.

Per panel site, in this order:

1. **Rows at the site.** Every `calls` row of the facet at the site's `(chrom, pos)` is read against
   the panel's two alleles. A row is usable when it has no `split_signal` (in a pre-split file a `0`
   can mean another ALT; `rejoin_split` re-joins such rows first), exactly two haplotypes, and both
   of them one of the panel's two alleles (a null haplotype is neither). Its dosage is the number of
   haplotypes carrying the panel's effect allele. Exactly one usable row gives the site's dosage.
   Rows that are there but none usable, or two usable rows, make the site **unusable**: missing, never
   0. A no-call is the source saying it looked and could not tell.
2. **No row.** What that means is the facet's `absence_means`: `weak_reference` (a whole-genome callset
   without blocks) reads it as homozygous for the panel's reference allele, dosage 0; `blocks` (a
   gVCF) does so only inside a `ref` block; `reference` (CGI's VCF, whose producer called what it
   omits reference) does so outside every block, a no-call block being no omission; `unknown`, or
   null, leaves it missing. Absence is read as
   reference only on the autosomes, where homozygous means two copies.

The panel's first allele is taken as the reference allele; a panel whose first allele is not the
assembly's REF would turn every absence into the wrong genotype. Panel loaders state where their
first allele comes from.
"""

from dataclasses import dataclass

import numpy as np
import polars as pl

from just_vcf.ancestry._projection import DEFAULT_MISSING
from just_vcf.loading import Loaded
from just_vcf.schema import CONTIG_ALIASES, AbsenceMeans, BlockKind

AUTOSOMES: frozenset[str] = frozenset(str(n) for n in range(1, 23))
"""Canonical `chrom` names on which an absence can be read as homozygous reference."""

PANEL_SITE_COLUMNS = ("chrom", "pos", "ref_allele", "effect_allele")
"""What a panel hands to `site_dosages`, one row per site in the panel's own order."""


@dataclass(frozen=True)
class Facet:
    """The facet facts an ancestry estimate needs."""

    facet_id: str
    assembly_label: str | None
    absence_means: str | None


@dataclass(frozen=True)
class SiteDosages:
    """Effect-allele dosage per panel site, `DEFAULT_MISSING` where unknown, in the panel's order."""

    dosage: np.ndarray
    n_observed: int
    n_absent_as_reference: int
    n_unusable: int

    @property
    def n_used(self) -> int:
        return self.n_observed + self.n_absent_as_reference


def facet_of(loaded: Loaded, facet_id: str | None = None) -> Facet:
    """The one facet of `loaded` to read, or the one named `facet_id`."""
    facets = loaded.facets.select("facet_id", "assembly_label", "absence_means").collect()
    if facet_id is not None:
        facets = facets.filter(pl.col("facet_id") == facet_id)
        if facets.height == 0:
            raise ValueError(f"no facet {facet_id!r} in this Loaded")
    if facets.height != 1:
        raise ValueError(f"{facets.height} facets in this Loaded; name the one to read with facet_id")
    row = facets.row(0, named=True)
    absence = row["absence_means"]
    return Facet(row["facet_id"], row["assembly_label"], None if absence is None else str(absence))


def panel_sites(frame: pl.DataFrame) -> pl.DataFrame:
    """A panel's site table in the shape `site_dosages` reads: canonical `chrom`, `pos` as `UInt32`.

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


def site_dosages(loaded: Loaded, facet: Facet, sites: pl.DataFrame) -> SiteDosages:
    """The facet's effect-allele dosage at each panel site (module docstring has the rules)."""
    sites = sites.with_row_index("_site")
    keys = sites.lazy().select("_site", "chrom", "pos", "ref_allele", "effect_allele")
    rows = (
        loaded.calls.filter(pl.col("facet_id") == facet.facet_id)
        .select("chrom", "pos", "gt_alleles", "split_signal")
        .join(keys.select("chrom", "pos").unique(), on=["chrom", "pos"], how="semi")
        .with_row_index("_row")
    )
    pairs = keys.join(rows, on=["chrom", "pos"], how="inner").with_columns(
        _ploidy=pl.col("gt_alleles").list.len(),
        _split=pl.col("split_signal").is_not_null(),
    )
    allele = pl.col("_allele")
    per_row = (
        pairs.select("_site", "_row", "_ploidy", "_split", "ref_allele", "effect_allele", pl.col("gt_alleles").alias("_allele"))
        .explode("_allele", empty_as_null=True)
        .group_by("_site", "_row")
        .agg(
            pl.col("_ploidy").first(),
            pl.col("_split").first(),
            _in_pair=((allele == pl.col("ref_allele")) | (allele == pl.col("effect_allele"))).fill_null(False).sum(),
            _effect=(allele == pl.col("effect_allele")).fill_null(False).sum(),
        )
        .with_columns(
            _usable=(~pl.col("_split"))
            & (pl.col("_ploidy") == 2)
            & (pl.col("_in_pair") == 2)
        )
    )
    per_site = (
        per_row.group_by("_site")
        .agg(
            _n_usable=pl.col("_usable").sum(),
            _dosage=pl.col("_effect").filter(pl.col("_usable")).first(),
        )
        .collect()
    )

    dosage = np.full(sites.height, DEFAULT_MISSING, dtype=np.float64)
    observed = per_site.filter(pl.col("_n_usable") == 1)
    dosage[observed["_site"].to_numpy()] = observed["_dosage"].cast(pl.Float64).to_numpy()
    present = np.zeros(sites.height, dtype=bool)
    present[per_site["_site"].to_numpy()] = True

    reference = _absent_as_reference(loaded, facet, sites, present)
    dosage[reference] = 0.0
    return SiteDosages(
        dosage=dosage,
        n_observed=observed.height,
        n_absent_as_reference=int(reference.sum()),
        n_unusable=per_site.height - observed.height,
    )


def _absent_as_reference(loaded: Loaded, facet: Facet, sites: pl.DataFrame, present: np.ndarray) -> np.ndarray:
    """A mask of the sites with no row that the facet's `absence_means` reads as reference."""
    autosomal = sites["chrom"].is_in(list(AUTOSOMES)).to_numpy()
    candidates = ~present & autosomal
    if facet.absence_means == AbsenceMeans.WEAK_REFERENCE.value:
        return candidates
    if facet.absence_means == AbsenceMeans.BLOCKS.value:
        return candidates & _inside_block(loaded, facet, sites, BlockKind.REF)
    if facet.absence_means == AbsenceMeans.REFERENCE.value:
        # What the file omits was called reference; a site a block covers was not omitted (a no-call).
        return candidates & ~_inside_block(loaded, facet, sites, None)
    return np.zeros(sites.height, dtype=bool)


def _inside_block(loaded: Loaded, facet: Facet, sites: pl.DataFrame, kind: BlockKind | None) -> np.ndarray:
    """Which sites a block of the facet covers, of `kind` or of any kind (blocks are 1-based and
    inclusive). Each block start carries the furthest end of the blocks up to it, so a site inside a
    long block that a shorter one starts within is still covered."""
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
    covered = (
        probe.join_asof(blocks, left_on="pos", right_on="_block_start", by="chrom", strategy="backward", check_sortedness=False)
        .filter(pl.col("_block_end") >= pl.col("pos"))
        .select("_site")
        .collect()
    )
    mask = np.zeros(sites.height, dtype=bool)
    mask[covered["_site"].to_numpy()] = True
    return mask
