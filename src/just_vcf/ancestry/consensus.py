"""Combining estimates from several panels or methods over five continental groups.

Ported from just-prs `ancestry/__init__.py` (`to_canonical_superpops`, `bayesian_consensus`) and
`ancestry/prive.py` (`prive_group_to_continental`); docs/PROVENANCE.md.
"""

import math
from collections.abc import Sequence

from just_vcf.ancestry.result import AncestryEstimate

CANONICAL_SUPERPOPS: tuple[str, ...] = ("AFR", "AMR", "EAS", "EUR", "SAS")
"""The 1000 Genomes super-populations: the shared vocabulary a consensus is computed over."""

_LABEL_TO_CANONICAL = {"AFR": "AFR", "AMR": "AMR", "EAS": "EAS", "EUR": "EUR", "SAS": "SAS", "CSA": "SAS"}
"""Broad panel labels. HGDP+1kGP's gnomAD-style `CSA` (Central/South Asian) is `SAS`; `MID` (Middle
Eastern) has no 1000 Genomes counterpart and maps to nothing."""


def continental_of(group: str) -> str | None:
    """The continental group a panel group folds onto, or None where none fits.

    Broad 1000 Genomes and HGDP labels map directly. Privé's worldwide groups (`Africa (West)`,
    `Europe (North East)`, `Japan`, …) are rolled up by keyword, as just-prs does. Middle Eastern
    groups, and anything unrecognised, give None: their weight is dropped and the rest renormalised,
    rather than forced into a group that does not fit.
    """
    if group in _LABEL_TO_CANONICAL:
        return _LABEL_TO_CANONICAL[group]
    name = group.lower()
    if "africa" in name:
        return "AFR"
    if any(k in name for k in ("america", "latin", "hispanic")):
        return "AMR"
    if any(k in name for k in ("east asia", "asia (east", "japan", "china", "korea", "philippin")):
        return "EAS"
    if any(k in name for k in ("south asia", "asia (south", "pakistan", "bangladesh", "sri lanka", "india")):
        return "SAS"
    if any(k in name for k in ("middle east", "near east")):
        return None
    if any(
        k in name
        for k in (
            "europe", "scandinav", "united kingdom", "ireland", "italy", "finland",
            "ashkenazi", "spain", "portugal", "france", "germany",
        )
    ):
        return "EUR"
    return None


def to_continental(distribution: dict[str, float], group_to_continental: dict[str, str | None]) -> dict[str, float]:
    """Fold a distribution over panel groups onto `CANONICAL_SUPERPOPS`, renormalised.

    A distribution whose every group maps to nothing comes back all zero.
    """
    folded = dict.fromkeys(CANONICAL_SUPERPOPS, 0.0)
    for group, value in distribution.items():
        target = group_to_continental.get(group, continental_of(group))
        if target is not None:
            folded[target] += float(value)
    total = sum(folded.values())
    return {group: (value / total if total > 0 else 0.0) for group, value in folded.items()}


def consensus(
    estimates: Sequence[AncestryEstimate],
    prior: dict[str, float] | None = None,
    smoothing: float = 1e-3,
) -> AncestryEstimate:
    """One posterior over `CANONICAL_SUPERPOPS` from several estimates (a product of experts).

    Each answered estimate contributes its neighbour vote (`probabilities`) and its `proportions`,
    folded onto the five continental groups and smoothed by `smoothing`, so no single method can veto
    a group. Unknown estimates contribute nothing (an estimate that could not run is not evidence).
    The product treats every contribution as independent, which two methods on one panel are not, so
    agreement sharpens the posterior more than it should: read `probabilities` as relative support,
    not as calibrated probabilities. With no answered estimate the result is unknown (`no_estimates`).
    """
    answered = [estimate for estimate in estimates if estimate]
    panels = "+".join(sorted({estimate.panel for estimate in estimates})) or "none"
    facet_ids = tuple(sorted({facet for estimate in estimates for facet in estimate.facet_ids}))
    if not answered:
        reasons = sorted({estimate.reason for estimate in estimates if estimate.reason})
        return AncestryEstimate(
            method="consensus",
            panel=panels,
            facet_ids=facet_ids,
            reason="no_estimates",
            detail=f"none of {len(estimates)} estimates answered ({', '.join(reasons) or 'none given'})",
        )
    log_posterior = {
        group: math.log(prior[group]) if prior and prior.get(group, 0) > 0 else 0.0 for group in CANONICAL_SUPERPOPS
    }
    for estimate in answered:
        for distribution in (estimate.probabilities, estimate.proportions):
            if not distribution:
                continue
            folded = to_continental(distribution, estimate.group_to_continental)
            for group in CANONICAL_SUPERPOPS:
                log_posterior[group] += math.log(folded[group] + smoothing)
    top = max(log_posterior.values())
    scaled = {group: math.exp(value - top) for group, value in log_posterior.items()}
    total = sum(scaled.values())
    posterior = {group: round(value / total, 4) for group, value in scaled.items()}
    leader = max(CANONICAL_SUPERPOPS, key=lambda group: (posterior[group], -CANONICAL_SUPERPOPS.index(group)))
    ranked = sorted(posterior.values(), reverse=True)
    return AncestryEstimate(
        method="consensus",
        panel=panels,
        facet_ids=facet_ids,
        estimate=leader,
        probabilities=posterior,
        margin=ranked[0] - ranked[1],
        group_to_continental={group: group for group in CANONICAL_SUPERPOPS},
    )
