"""The facet graph's edges: the `links` frame (RM275, docs/SCHEMA.md Decision 12).

A link says how two facets relate: `from_facet_id` → `to_facet_id`, a `relation`, and the
`provenance` of the step that made it. Lineage links point from the facet a step read to the facet it
made, as Cat-VRS reads `liftover_to` (from the defining allele to the related one); a facet can have
several lineage parents (a merge, a consensus), which `facets.derived_from` cannot hold. Kin links join
facets read from one source, neither made from the other (a caller split, an appended export).

`facets.derived_from` and `facets.relation` stay for the one-parent case; `link_disagreements` states
how the two agree and lists every place they do not. This module does no I/O.
"""

from itertools import combinations

import polars as pl

from just_vcf.schema import LINKS, empty

LIFTOVER_TO = "liftover_to"
"""The Cat-VRS relation, verbatim: the `to` facet is the `from` facet lifted to another assembly."""
NORMALIZED = "normalized"
"""just-vcf's: the `to` facet is the `from` facet split, left-aligned, trimmed and deduplicated."""
SPLIT_REJOINED = "split_rejoined"
"""just-vcf's: the `to` facet is the `from` facet with its pre-split sites joined again (RM270)."""
MERGE = "merge"
"""just-vcf's: several parents, one facet whose rows are the parents' rows put together (a merge).
Documented; no step emits it yet."""
CONSENSUS = "consensus"
"""just-vcf's: several parents, one facet whose rows are computed from all of them: a posterior per
allele that weighs each parent's evidence (`consensus`, RM280). The parents keep their own calls."""
CALLER_SPLIT = "caller_split"
"""just-vcf's: facets holding one caller's rows each of a file merging several callers (Decision 8)."""
APPENDED_EXPORT = "appended_export"
"""just-vcf's: an export appended to a file after another one, read as its own facet (F3). GA4GH has
no term for two exports in one file."""

LINEAGE_RELATIONS: frozenset[str] = frozenset({LIFTOVER_TO, NORMALIZED, SPLIT_REJOINED, MERGE, CONSENSUS})
"""Relations where the `to` facet was made from the `from` facet."""
MULTI_PARENT_RELATIONS: frozenset[str] = frozenset({MERGE, CONSENSUS})
"""Lineage relations a facet holds with two or more parents, so `facets.derived_from` stays null."""
SINGLE_PARENT_RELATIONS: frozenset[str] = LINEAGE_RELATIONS - MULTI_PARENT_RELATIONS
"""Lineage relations with one parent, the ones `facets.derived_from` holds: a facet linked so is its
parent's callset, made again (lifted, normalized, rejoined)."""
KIN_RELATIONS: frozenset[str] = frozenset({CALLER_SPLIT, APPENDED_EXPORT})
"""Relations between facets of one source, neither made from the other."""
RELATIONS: frozenset[str] = LINEAGE_RELATIONS | KIN_RELATIONS
"""Every documented relation. It grows additively: `relation` is `pl.String`, not a `pl.Enum`."""

LINK_KEY = ("from_facet_id", "to_facet_id", "relation")
"""The canonical order of `links`, and what makes a link one link."""


def is_relation(value: str) -> bool:
    """Whether `value` is a documented relation."""
    return value in RELATIONS


def canonical_links(links: pl.LazyFrame) -> pl.LazyFrame:
    """`links` sorted by `LINK_KEY`, each key once (the first in that order kept)."""
    return links.sort(list(LINK_KEY), nulls_last=True, maintain_order=True).unique(subset=list(LINK_KEY), keep="first", maintain_order=True)


def _provenance(step: str, step_version: str | None, derivation_inputs: list[dict] | None) -> dict:
    return {"step": step, "step_version": step_version, "derivation_inputs": derivation_inputs or []}


def links_frame(rows: list[dict]) -> pl.LazyFrame:
    """A `links` frame from rows of its four columns, in canonical order."""
    return canonical_links(pl.LazyFrame(rows, schema=LINKS)) if rows else empty(LINKS)


def derived_links(
    parent: pl.LazyFrame,
    from_facet_id: str,
    to_facet_id: str,
    relation: str,
    *,
    step: str,
    step_version: str,
    derivation_inputs: list[dict] | None = None,
) -> pl.LazyFrame:
    """The parent's links with one lineage link added: `from_facet_id` → `to_facet_id`.

    The parent's links are carried along, so a facet three steps from its file still says how each step
    got there when only the last `Loaded` is kept. `derivation_inputs` are this step's own files (a
    liftover's chain), `[]` when it read none.
    """
    if relation not in LINEAGE_RELATIONS:
        raise ValueError(f"{relation!r} is not a lineage relation")
    added = links_frame(
        [
            {
                "from_facet_id": from_facet_id,
                "to_facet_id": to_facet_id,
                "relation": relation,
                "provenance": _provenance(step, step_version, derivation_inputs),
            }
        ]
    )
    return canonical_links(pl.concat([parent.select(LINKS.names()), added], how="vertical"))


def kin_links(facet_ids: list[str], relation: str, *, step: str = "load") -> pl.LazyFrame:
    """Kin links between facets one load read from one source, `relation` on each.

    With no order the source gives (a caller split), every pair is linked, from the lower `facet_id` to
    the higher. `facet_ids` in source order, with an order that means something (an appended export:
    the first export, then the one appended after it), is linked in that order, each to the next.
    """
    if relation not in KIN_RELATIONS:
        raise ValueError(f"{relation!r} is not a kin relation")
    if relation == APPENDED_EXPORT:
        pairs = list(zip(facet_ids, facet_ids[1:]))
    else:
        pairs = list(combinations(sorted(facet_ids), 2))
    return links_frame(
        [
            {"from_facet_id": a, "to_facet_id": b, "relation": relation, "provenance": _provenance(step, None, [])}
            for a, b in pairs
        ]
    )


def link_disagreements(facets: pl.DataFrame, links: pl.DataFrame) -> list[str]:
    """Every place `facets.derived_from` / `relation` and `links` do not agree; `[]` when they agree.

    The rules, per facet of `facets`:

    - `derived_from` set: exactly one lineage link ends at the facet, from `derived_from`, with the
      facet's `relation`, which is a lineage relation and not a multi-parent one;
    - `derived_from` null, `relation` a multi-parent relation (`merge`, `consensus`): two or more lineage
      links end at the facet, all with that relation;
    - `derived_from` null, `relation` a kin relation: no lineage link ends at the facet, and at least one
      link of that relation names it at either end;
    - both null: no lineage link ends at the facet.

    And per link: its relation is documented, and a kin link whose two facets are both in `facets` joins
    facets of one source (same `source_sha256` and `member`). A link may name a facet `facets` does not
    hold: a derived `Loaded` carries its ancestors' links.
    """
    found: list[str] = []
    unknown = sorted(set(links["relation"].to_list()) - RELATIONS)
    found += [f"relation {value!r} is not documented" for value in unknown]
    lineage = links.filter(pl.col("relation").is_in(sorted(LINEAGE_RELATIONS)))
    for facet in facets.iter_rows(named=True):
        facet_id, parent, relation = facet["facet_id"], facet["derived_from"], facet["relation"]
        into = lineage.filter(pl.col("to_facet_id") == facet_id)
        if parent is not None:
            if relation not in SINGLE_PARENT_RELATIONS:
                found.append(f"facet {facet_id}: derived_from is set but relation {relation!r} is not a one-parent lineage relation")
            if into.height != 1:
                found.append(f"facet {facet_id}: derived_from is set but {into.height} lineage link(s) end at it, not 1")
            elif into.row(0, named=True)["from_facet_id"] != parent or into.row(0, named=True)["relation"] != relation:
                link = into.row(0, named=True)
                found.append(
                    f"facet {facet_id}: derived_from {parent} with {relation!r}, but its link comes from "
                    f"{link['from_facet_id']} with {link['relation']!r}"
                )
        elif relation in MULTI_PARENT_RELATIONS:
            if into.height < 2 or set(into["relation"].to_list()) != {relation}:
                found.append(f"facet {facet_id}: relation {relation!r} needs two or more {relation!r} links into it; {into.height} lineage link(s) end at it")
        elif relation in LINEAGE_RELATIONS:
            found.append(f"facet {facet_id}: relation {relation!r} is lineage but derived_from is null")
        else:
            if into.height:
                found.append(f"facet {facet_id}: derived_from is null but {into.height} lineage link(s) end at it")
            if relation is not None and relation not in KIN_RELATIONS:
                found.append(f"facet {facet_id}: relation {relation!r} is not documented")
            elif relation is not None:
                touching = links.filter(
                    (pl.col("relation") == relation) & ((pl.col("from_facet_id") == facet_id) | (pl.col("to_facet_id") == facet_id))
                )
                if touching.height == 0:
                    found.append(f"facet {facet_id}: relation {relation!r} but no {relation!r} link names it")
    sources = {row["facet_id"]: (row["source_sha256"], row["member"]) for row in facets.select("facet_id", "source_sha256", "member").iter_rows(named=True)}
    for link in links.filter(pl.col("relation").is_in(sorted(KIN_RELATIONS))).iter_rows(named=True):
        a, b = link["from_facet_id"], link["to_facet_id"]
        if a in sources and b in sources and sources[a] != sources[b]:
            found.append(f"{link['relation']} link {a} -> {b} joins facets of two sources")
    return found

