"""A facet's genotypes at an ancestry panel's sites, as the numpy arrays the projections read.

The reading itself, with its rules for rows, no-calls and absences, is the core's
`just_vcf.site_genotypes`, shared with relatedness (RM274). This module turns its per-site counts into
one effect-allele dosage per site, `DEFAULT_MISSING` where unknown.
"""

from dataclasses import dataclass

import numpy as np
import polars as pl

from just_vcf.ancestry._projection import DEFAULT_MISSING
from just_vcf.loading import Loaded
from just_vcf.site_genotypes import AUTOSOMES, PANEL_SITE_COLUMNS, Facet, SiteBasis, facet_of, panel_sites, read_site_genotypes

__all__ = ["AUTOSOMES", "PANEL_SITE_COLUMNS", "Facet", "SiteDosages", "facet_of", "panel_sites", "site_dosages"]


@dataclass(frozen=True)
class SiteDosages:
    """Effect-allele dosage per panel site, `DEFAULT_MISSING` where unknown, in the panel's order."""

    dosage: np.ndarray
    n_observed: int
    n_absent_as_reference: int
    n_unusable: int
    n_ref_unstated: int
    """Sites whose every row spells a sequence ALT against a source REF that states no base (F23)."""

    @property
    def n_used(self) -> int:
        return self.n_observed + self.n_absent_as_reference


def site_dosages(loaded: Loaded, facet: Facet, sites: pl.DataFrame) -> SiteDosages:
    """The facet's effect-allele dosage at each panel site (`just_vcf.site_genotypes` has the rules)."""
    read = read_site_genotypes(loaded, facet, sites)
    dosage = read["count"].cast(pl.Float64).fill_null(DEFAULT_MISSING).to_numpy().copy()
    basis = read["basis"]
    return SiteDosages(
        dosage=dosage,
        n_observed=int((basis == SiteBasis.OBSERVED.value).sum()),
        n_absent_as_reference=int((basis == SiteBasis.ABSENT_AS_REFERENCE.value).sum()),
        n_unusable=int((basis == SiteBasis.UNUSABLE.value).sum()),
        n_ref_unstated=int((basis == SiteBasis.REF_UNSTATED.value).sum()),
    )
