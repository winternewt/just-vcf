"""GA4GH Pedigree, read as data the caller supplies (RM274). Nothing here infers a relationship.

The GA4GH Pedigree standard (https://pedigree.readthedocs.io) has three parts: a `Pedigree` with an
`id`, its `individuals` and its `relationships`; an `Individual` with an `id` and a `sex`; and a
`Relationship` saying what its `individual` is to its `relative`, the `relation` a Kinship Ontology
(KIN) term such as `KIN:027` (`isBiologicalMotherOf`). Its documented serialisation is the Phenopacket
`GA4GHPedigree` message, where an individual's `id` and `sex` sit under `subject` and a relationship
names `individual_id` and `relative_id`. `read_pedigree` reads that message as JSON, with protobuf's two JSON spellings of a field (`individual_id`, `individualId`) both accepted, into
two frames. Terms are kept verbatim (charter invariant 2): the KIN id, its label as written, the sex
value as written.

A pedigree states relationships; `just_vcf.relatedness` measures how alike genotypes are. Comparing
the two is the caller's: a stated parent whose facets read unrelated is a finding, not something this
module resolves.
"""

import json
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

PEDIGREE_INDIVIDUALS = pl.Schema({"pedigree_id": pl.String, "individual_id": pl.String, "sex": pl.String})
"""One row per individual: the `subject.id` and `subject.sex` as written."""

PEDIGREE_RELATIONSHIPS = pl.Schema(
    {
        "pedigree_id": pl.String,
        "individual_id": pl.String,
        "relation_id": pl.String,
        "relation_label": pl.String,
        "relative_id": pl.String,
    }
)
"""One row per relationship: `individual_id` is `relation` of `relative_id` (KIN's direction)."""


class PedigreeError(ValueError):
    """A pedigree document that does not have the standard's required parts."""


class Pedigree(BaseModel):
    """A supplied GA4GH Pedigree as two frames, with its own fields kept as written."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    pedigree_id: str
    individuals: pl.DataFrame
    relationships: pl.DataFrame
    index_patients: tuple[str, ...] = ()
    narrative: str | None = None
    date: str | None = None


def _field(item: dict, snake: str) -> object:
    head, *rest = snake.split("_")
    camel = head + "".join(part.title() for part in rest)
    return item.get(snake, item.get(camel))


def pedigree_from_dict(document: dict) -> Pedigree:
    """A `GA4GHPedigree` message, as parsed JSON, read into a `Pedigree`.

    Refused with `PedigreeError`: no `id`; an individual with no `subject.id` (or `id` where there is no
    `subject`); a relationship missing `individual_id`, `relation.id` or `relative_id`; a relationship
    naming an individual the pedigree does not list; two individuals with one id.
    """
    pedigree_id = document.get("id")
    if not pedigree_id:
        raise PedigreeError("a pedigree needs an id")
    people = []
    for item in document.get("individuals") or []:
        subject = item.get("subject", item)
        if not subject.get("id"):
            raise PedigreeError(f"pedigree {pedigree_id}: an individual without an id")
        sex = subject.get("sex")
        people.append({"pedigree_id": str(pedigree_id), "individual_id": str(subject["id"]), "sex": None if sex is None else str(sex)})
    individuals = pl.DataFrame(people, schema=PEDIGREE_INDIVIDUALS)
    if individuals["individual_id"].is_duplicated().any():
        raise PedigreeError(f"pedigree {pedigree_id}: two individuals share an id")
    known = set(individuals["individual_id"].to_list())
    links = []
    for item in document.get("relationships") or []:
        individual, relative, relation = _field(item, "individual_id"), _field(item, "relative_id"), item.get("relation") or {}
        if not individual or not relative or not relation.get("id"):
            raise PedigreeError(f"pedigree {pedigree_id}: a relationship needs individual_id, relation.id and relative_id")
        unknown = {str(individual), str(relative)} - known
        if unknown:
            raise PedigreeError(f"pedigree {pedigree_id}: a relationship names {', '.join(sorted(unknown))}, not among its individuals")
        label = relation.get("label")
        links.append(
            {
                "pedigree_id": str(pedigree_id),
                "individual_id": str(individual),
                "relation_id": str(relation["id"]),
                "relation_label": None if label is None else str(label),
                "relative_id": str(relative),
            }
        )
    relationships = pl.DataFrame(links, schema=PEDIGREE_RELATIONSHIPS)
    index = tuple(str(value) for value in (_field(document, "index_patients") or []))
    date = document.get("date")
    return Pedigree(
        pedigree_id=str(pedigree_id),
        individuals=individuals,
        relationships=relationships,
        index_patients=index,
        narrative=document.get("narrative"),
        date=None if date is None else str(date),
    )


def read_pedigree(path: Path | str) -> Pedigree:
    """A `GA4GHPedigree` message in a JSON file, read into a `Pedigree` (`pedigree_from_dict`)."""
    return pedigree_from_dict(json.loads(Path(path).read_text()))
