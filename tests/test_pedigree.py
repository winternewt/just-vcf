"""GA4GH Pedigree read as supplied data (RM274), from the standard's own trio and twins examples."""

import json
from pathlib import Path

import pytest

from just_vcf.pedigree import PEDIGREE_INDIVIDUALS, PEDIGREE_RELATIONSHIPS, PedigreeError, pedigree_from_dict, read_pedigree

# The standard's "Basic Trio" example (pedigree.readthedocs.io, examples), as the message's JSON.
TRIO = {
    "id": "FAM1",
    "narrative": "A Phenopacket GA4GHPedigree of a trio with an affected child",
    "date": "2022-06-23",
    "individuals": [
        {"id": 1, "subject": {"id": "MOTHER", "sex": "FEMALE"}},
        {"id": 2, "subject": {"id": "FATHER", "sex": "MALE"}},
        {"id": 3, "subject": {"id": "CHILD", "sex": "UNKNOWN"}},
    ],
    "relationships": [
        {"individual_id": "MOTHER", "relation": {"id": "KIN:027", "label": "isBiologicalMotherOf"}, "relative_id": "CHILD"},
        {"individual_id": "FATHER", "relation": {"id": "KIN:028", "label": "isBiologicalFatherOf"}, "relative_id": "CHILD"},
    ],
    "index_patients": ["CHILD"],
}


def test_the_trio_example_reads_verbatim(tmp_path: Path) -> None:
    path = tmp_path / "trio.json"
    path.write_text(json.dumps(TRIO))
    pedigree = read_pedigree(path)
    assert pedigree.individuals.schema == PEDIGREE_INDIVIDUALS and pedigree.relationships.schema == PEDIGREE_RELATIONSHIPS
    assert pedigree.individuals.select("individual_id", "sex").rows() == [
        (item["subject"]["id"], item["subject"]["sex"]) for item in TRIO["individuals"]
    ]
    assert pedigree.relationships.drop("pedigree_id").rows() == [
        (item["individual_id"], item["relation"]["id"], item["relation"]["label"], item["relative_id"]) for item in TRIO["relationships"]
    ]
    assert (pedigree.pedigree_id, pedigree.index_patients, pedigree.date) == ("FAM1", ("CHILD",), "2022-06-23")


def test_protobuf_lower_camel_names_read_alike() -> None:
    camel = {
        **{key: value for key, value in TRIO.items() if key not in ("relationships", "index_patients")},
        "relationships": [
            {"individualId": item["individual_id"], "relation": item["relation"], "relativeId": item["relative_id"]} for item in TRIO["relationships"]
        ],
        "indexPatients": ["CHILD"],
    }
    read, expected = pedigree_from_dict(camel), pedigree_from_dict(TRIO)
    assert read.relationships.equals(expected.relationships) and read.individuals.equals(expected.individuals)
    assert read.index_patients == expected.index_patients == ("CHILD",)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda d: d.pop("id"), "needs an id"),
        (lambda d: d["relationships"][0]["relation"].pop("id"), "needs individual_id, relation.id and relative_id"),
        (lambda d: d["relationships"][0].update(relative_id="UNCLE"), "names UNCLE"),
        (lambda d: d["individuals"][1]["subject"].update(id="MOTHER"), "share an id"),
    ],
)
def test_a_pedigree_without_its_required_parts_is_refused(change, message: str) -> None:
    document = json.loads(json.dumps(TRIO))
    change(document)
    with pytest.raises(PedigreeError, match=message):
        pedigree_from_dict(document)
