"""Combining estimates from several panels or methods over five continental groups.

Ported from just-prs `ancestry/__init__.py` (`to_canonical_superpops`, `bayesian_consensus`) and
`ancestry/prive.py` (`prive_group_to_continental`); docs/PROVENANCE.md. just-prs multiplies every
distribution it is given as independent evidence; here the methods of panels that share reference
samples are one witness and are averaged, and only witnesses the caller declared disjoint are
multiplied (docs/SCHEMA.md, "Combining ancestry estimates").
"""

import math
from collections.abc import Iterable, Sequence

from just_vcf.ancestry.result import AncestryEstimate, ConsensusWitness, WitnessReading

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
    *,
    disjoint_panels: Iterable[tuple[str, str]] = (),
) -> AncestryEstimate:
    """One answer over `CANONICAL_SUPERPOPS` from several estimates, one panel's methods as one witness.

    Each answered estimate gives its neighbour vote (`probabilities`) and its `proportions` as readings,
    folded onto the five continental groups. Readings from panels that share reference samples see the
    same people at the same sites, so their agreement is not independent evidence: they are one witness,
    and a witness says the mean of its readings, which gives its leader no more than the most any one
    reading gave it. Two panels are taken as having disjoint reference samples only where
    `disjoint_panels` names the pair; nothing is inferred from panel names, so by default every panel is
    one witness with every other. A witness is a connected group: panels a, b and c with only (a, b)
    declared are one witness, since c may share samples with both.

    One witness and no `prior`: the answer is its mean (`combination = mean`). Otherwise the witnesses'
    means, each smoothed by `smoothing` so that no witness can veto a group, are multiplied with the
    prior as independent evidence (`product_of_means`). A `prior` names all five groups, each above 0.
    `witnesses` lists what was pooled, each with its spread, so disagreement inside one panel stays in
    view. Read `probabilities` as relative support, not as calibrated probabilities.

    Unknown estimates contribute nothing (an estimate that could not run is not evidence). With no
    answered estimate the result is unknown (`no_estimates`); with no reading that folds onto a
    continental group, `no_continental_groups`. A consensus is not taken as an estimate to combine: pass
    the estimates it was made of.
    """
    if any(estimate.method == "consensus" for estimate in estimates):
        raise ValueError("a consensus is not combined again: pass the estimates it was made of")
    if prior is not None and (set(prior) != set(CANONICAL_SUPERPOPS) or not all(prior[g] > 0 for g in CANONICAL_SUPERPOPS)):
        raise ValueError(f"a prior gives every one of {list(CANONICAL_SUPERPOPS)} a weight above 0; got {prior}")
    known = {estimate.panel for estimate in estimates}
    declared: set[frozenset[str]] = set()
    for pair in disjoint_panels:
        first, second = pair
        if first == second:
            raise ValueError(f"panel {first!r} cannot be disjoint from itself")
        unnamed = sorted({first, second} - known)
        if unnamed:
            raise ValueError(f"disjoint_panels names {unnamed}, which no estimate is on; panels: {sorted(known)}")
        declared.add(frozenset(pair))

    panels = "+".join(sorted(known)) or "none"
    facet_ids = tuple(sorted({facet for estimate in estimates for facet in estimate.facet_ids}))
    answered = [estimate for estimate in estimates if estimate]
    if not answered:
        reasons = sorted({estimate.reason for estimate in estimates if estimate.reason})
        return AncestryEstimate(
            method="consensus",
            panel=panels,
            facet_ids=facet_ids,
            reason="no_estimates",
            detail=f"none of {len(estimates)} estimates answered ({', '.join(reasons) or 'none given'})",
        )

    witnesses = tuple(_witness(group, answered) for group in _witness_panels(sorted({e.panel for e in answered}), declared))
    speaking = [witness for witness in witnesses if witness.readings]
    if not speaking:
        return AncestryEstimate(
            method="consensus",
            panel=panels,
            facet_ids=facet_ids,
            reason="no_continental_groups",
            detail=f"none of the groups of {len(answered)} answered estimates folds onto {list(CANONICAL_SUPERPOPS)}",
            witnesses=witnesses,
        )
    if len(speaking) == 1 and prior is None:
        combination = "mean"
        posterior = dict(speaking[0].mean)
    else:
        combination = "product_of_means"
        log_posterior = {group: math.log(prior[group]) if prior else 0.0 for group in CANONICAL_SUPERPOPS}
        for witness in speaking:
            for group in CANONICAL_SUPERPOPS:
                log_posterior[group] += math.log(witness.mean[group] + smoothing)
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
        combination=combination,
        witnesses=witnesses,
    )


def _witness_panels(panels: list[str], declared: set[frozenset[str]]) -> list[tuple[str, ...]]:
    """`panels` cut into witnesses: two panels are in one unless their pair is `declared` disjoint.

    The connected groups of the graph joining every undeclared pair, each sorted, in sorted order.
    """
    parent = {panel: panel for panel in panels}

    def root(panel: str) -> str:
        while parent[panel] != panel:
            panel = parent[panel]
        return panel

    for i, first in enumerate(panels):
        for second in panels[i + 1 :]:
            if frozenset((first, second)) not in declared:
                parent[root(second)] = root(first)
    groups: dict[str, list[str]] = {}
    for panel in panels:
        groups.setdefault(root(panel), []).append(panel)
    return sorted(tuple(sorted(group)) for group in groups.values())


def _witness(panels: tuple[str, ...], answered: Sequence[AncestryEstimate]) -> ConsensusWitness:
    """The readings of every answered estimate on `panels`, their mean and their spread."""
    readings = []
    left_out = []
    for estimate in answered:
        if estimate.panel not in panels:
            continue
        for reading, distribution in (("probabilities", estimate.probabilities), ("proportions", estimate.proportions)):
            if not distribution:
                continue
            folded = to_continental(distribution, estimate.group_to_continental)
            spoke = sum(folded.values()) > 0
            entry = WitnessReading(
                panel=estimate.panel,
                method=estimate.method,
                reading=reading,
                facet_ids=estimate.facet_ids,
                distribution={group: round(value, 4) for group, value in folded.items()} if spoke else {},
            )
            (readings if spoke else left_out).append((entry, folded))
    readings.sort(key=_reading_order)
    left_out.sort(key=_reading_order)
    mean: dict[str, float] = {}
    spread: dict[str, float] = {}
    if readings:
        for group in CANONICAL_SUPERPOPS:
            values = [folded[group] for _, folded in readings]
            mean[group] = round(sum(values) / len(values), 4)
            spread[group] = round(max(values) - min(values), 4)
    return ConsensusWitness(
        panels=panels,
        readings=tuple(entry for entry, _ in readings),
        mean=mean,
        spread=spread,
        left_out=tuple(entry for entry, _ in left_out),
    )


def _reading_order(item: tuple[WitnessReading, dict[str, float]]) -> tuple:
    """A stable order for readings whatever order the estimates came in (charter invariant 12)."""
    entry, folded = item
    return (entry.panel, entry.method, entry.reading, entry.facet_ids, tuple(folded[g] for g in CANONICAL_SUPERPOPS))
