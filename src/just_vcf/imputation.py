"""Imputed facets read back (RM287): the lineage a caller states between an imputed facet and its chip.

An imputation server makes the imputed file from a chip export (the charter: servers make the facets,
just-vcf reads what they produced). The file does not say which export went in, and just-vcf never
guesses it, so the caller states it: `state_imputed_from` links the chip's facet to the imputed one with
`imputed_from` (`links.IMPUTED_FROM`), a one-parent lineage relation, the tool and the reference panel
in the link's provenance. Stated, the two are one callset to `person_karyotype`, one source to a
person's relatedness checks, and a parent and its child to `consensus`. Unstated, they stay two
independent facets. This module does no I/O.
"""

import polars as pl

from just_vcf.links import IMPUTED_FROM, LINEAGE_RELATIONS, canonical_links, links_frame
from just_vcf.loading import Loaded
from just_vcf.person import SubjectError
from just_vcf.schema import FACETS, LINKS

PANEL_ROLE = "reference_panel"
"""The `role` of the reference panel in `derivation_inputs`, on the link's provenance and the facet."""


def state_imputed_from(
    imputed: Loaded,
    chip: Loaded,
    *,
    tool: str,
    tool_version: str | None,
    panel: str | None,
    panel_sha256: str | None,
) -> Loaded:
    """`imputed` with the caller's statement that it was imputed from `chip`'s facet.

    Each `Loaded` holds one facet. The imputed facet keeps its `facet_id`, calls and blocks: its rows
    are what the server wrote, and only its lineage is new. Its `derived_from` becomes the chip's
    `facet_id` and its `relation` `imputed_from`, and the link `chip → imputed` is added, with `step`
    the tool, `step_version` its version and `derivation_inputs` the panel (`role` `reference_panel`,
    `sha256` the panel's digest where the caller has one); the facet's own `derivation_inputs` gain the
    panel too. A `None` the caller passes stays unknown: no panel is listed, no version written. The
    chip's links come along, so a chip lifted before imputation still names its first facet.

    Refused (`ValueError`): a `Loaded` holding other than one facet, one facet stated as its own parent,
    an imputed facet that already has lineage (a facet made by a step, or one stated already: state it
    on the facet read from the server's file), a chip made from the imputed facet, an empty tool, and
    two schema versions. Two subjects (both set and different) are refused with `SubjectError`. Nothing
    is inferred: neither facet's callset type is checked, since a server's output does not always say
    it is imputed, and a low-pass genome is imputed as well as a chip.
    """
    if not tool.strip():
        raise ValueError("state the imputation tool (`tool`); pass the name the server reports")
    child, parent = _one(imputed, "imputed"), _one(chip, "chip")
    child_id, parent_id = child["facet_id"], parent["facet_id"]
    if child_id == parent_id:
        raise ValueError(f"facet {child_id} cannot be imputed from itself")
    if child["derived_from"] is not None or child["relation"] is not None:
        raise ValueError(
            f"facet {child_id} has lineage already (relation {child['relation']!r}); state imputed_from on the facet read "
            "from the server's file"
        )
    if imputed.schema_version != chip.schema_version:
        raise ValueError(f"facets of different schema versions: {imputed.schema_version}, {chip.schema_version}")
    if child["subject"] is not None and parent["subject"] is not None and child["subject"] != parent["subject"]:
        raise SubjectError(f"facets {child_id} and {parent_id} carry different subjects; an imputed facet is its chip's person")
    chip_links = chip.links.select(LINKS.names()).collect()
    if child_id in _ancestors(parent, chip_links):
        raise ValueError(f"the chip {parent_id} was made from facet {child_id}; it cannot be that facet's parent")

    inputs = [] if panel is None else [{"role": PANEL_ROLE, "name": panel, "sha256": panel_sha256}]
    link = links_frame(
        [
            {
                "from_facet_id": parent_id,
                "to_facet_id": child_id,
                "relation": IMPUTED_FROM,
                "provenance": {"step": tool, "step_version": tool_version, "derivation_inputs": inputs},
            }
        ]
    )
    facets = imputed.facets.with_columns(
        pl.lit(parent_id, dtype=FACETS["derived_from"]).alias("derived_from"),
        pl.lit(IMPUTED_FROM, dtype=FACETS["relation"]).alias("relation"),
        pl.concat_list(
            pl.col("derivation_inputs").fill_null(pl.lit([], dtype=FACETS["derivation_inputs"])),
            pl.lit(inputs, dtype=FACETS["derivation_inputs"]),
        ).alias("derivation_inputs"),
    )
    links = canonical_links(pl.concat([imputed.links.select(LINKS.names()), chip_links.lazy(), link], how="vertical"))
    return imputed.model_copy(update={"facets": facets, "links": links})


def _one(loaded: Loaded, what: str) -> dict:
    rows = loaded.facets.collect()
    if rows.height != 1:
        raise ValueError(f"the {what} Loaded holds {rows.height} facets; pass one facet")
    return rows.row(0, named=True)


def _ancestors(facet: dict, links: pl.DataFrame) -> set[str]:
    """Every facet `facet` was made from, through lineage links and its `derived_from`."""
    parents: dict[str, set[str]] = {}
    for row in links.filter(pl.col("relation").is_in(sorted(LINEAGE_RELATIONS))).iter_rows(named=True):
        parents.setdefault(row["to_facet_id"], set()).add(row["from_facet_id"])
    if facet["derived_from"] is not None:
        parents.setdefault(facet["facet_id"], set()).add(facet["derived_from"])
    seen: set[str] = set()
    todo = list(parents.get(facet["facet_id"], ()))
    while todo:
        node = todo.pop()
        if node not in seen:
            seen.add(node)
            todo += list(parents.get(node, ()))
    return seen
