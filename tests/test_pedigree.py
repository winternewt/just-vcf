"""GA4GH Pedigree read as supplied data (RM274), from the standard's own trio and twins examples."""

import json
from pathlib import Path

import pytest

from just_vcf.pedigree import (
    EXPECTED_KINSHIP,
    MORE_THAN_TWO_PARENTS,
    PEDIGREE_INDIVIDUALS,
    PEDIGREE_RELATIONSHIPS,
    STATEMENTS_CONFLICT,
    TERM_WITHOUT_DEGREE,
    PedigreeError,
    expected_kinship,
    pedigree_from_dict,
    read_pedigree,
)

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


def _pedigree(people: str, links: list[tuple[str, str, str]]):
    """A pedigree of single-letter individuals and `(individual, KIN id, relative)` relationships."""
    return pedigree_from_dict(
        {
            "id": "F",
            "individuals": [{"subject": {"id": person}} for person in people],
            "relationships": [{"individual_id": a, "relation": {"id": term}, "relative_id": b} for a, term, b in links],
        }
    )


def _expected(pedigree) -> dict[frozenset[str], dict]:
    frame = expected_kinship(pedigree)
    assert frame.schema == EXPECTED_KINSHIP
    return {frozenset((row["individual_id"], row["relative_id"])): row for row in frame.iter_rows(named=True)}


def test_the_trio_expects_two_parent_child_pairs_and_unrelated_founders() -> None:
    rows = _expected(pedigree_from_dict(TRIO))
    for parent in ("MOTHER", "FATHER"):
        row = rows[frozenset((parent, "CHILD"))]
        assert (row["expected_kinship"], row["expected_bands"], row["basis"]) == (0.25, ["parent_child"], "derived")
        assert row["relation_id"] in ("KIN:027", "KIN:028")
    founders = rows[frozenset(("MOTHER", "FATHER"))]
    assert (founders["expected_kinship"], founders["expected_bands"], founders["basis"], founders["relation_id"]) == (
        0.0,
        ["unrelated"],
        "assumed",
        None,
    )


# Three generations from parent links alone: a+b have c and s; c+e have g and h; s+f have k;
# a+x have j (a half sibling of c); u is a founder linked to nobody.
FAMILY = _pedigree(
    "abcsefghjkxu",
    [
        ("a", "KIN:027", "c"), ("b", "KIN:028", "c"), ("a", "KIN:027", "s"), ("b", "KIN:028", "s"),
        ("c", "KIN:003", "g"), ("e", "KIN:003", "g"), ("h", "KIN:032", "c"), ("h", "KIN:032", "e"),
        ("s", "KIN:003", "k"), ("f", "KIN:003", "k"), ("a", "KIN:027", "j"), ("x", "KIN:028", "j"),
    ],
)  # fmt: skip


@pytest.mark.parametrize(
    ("pair", "kinship", "band"),
    [
        ("ac", 0.25, "parent_child"),  # parent and child
        ("hc", 0.25, "parent_child"),  # KIN:032 read the other way
        ("cs", 0.25, "full_siblings"),  # full siblings: the same two parents listed for both
        ("gh", 0.25, "full_siblings"),
        ("cj", 0.125, "second_degree"),  # half siblings
        ("ag", 0.125, "second_degree"),  # grandparent
        ("sg", 0.125, "second_degree"),  # aunt or uncle
        ("gk", 0.0625, "third_degree"),  # first cousins
        ("jg", 0.0625, "third_degree"),  # half aunt or uncle
        ("ab", 0.0, "unrelated"),  # founders
        ("au", 0.0, "unrelated"),
        ("ce", 0.0, "unrelated"),
    ],
)
def test_parent_links_give_every_degree(pair: str, kinship: float, band: str) -> None:
    row = _expected(FAMILY)[frozenset(pair)]
    assert (row["expected_kinship"], row["expected_bands"], row["expected_reason"]) == (kinship, [band], None)
    assert row["basis"] == ("assumed" if kinship == 0 else "derived")


def test_stated_degrees_hold_and_carry_to_relatives() -> None:
    """Siblings, cousins and twins stated without their parents: a stated sibling's child is a nephew."""
    pedigree = _pedigree(
        "abcdmnpq",
        [
            ("a", "KIN:008", "b"),  # full siblings, parents not listed
            ("a", "KIN:003", "c"),  # c is a's child, so b's nephew or niece
            ("d", "KIN:014", "a"),  # cousins
            ("m", "KIN:010", "n"),  # monozygotic twins
            ("p", "KIN:007", "q"),  # siblings, full or half
        ],
    )
    rows = _expected(pedigree)
    siblings = rows[frozenset("ab")]
    assert (siblings["expected_kinship"], siblings["expected_bands"], siblings["basis"], siblings["relation_id"]) == (
        0.25,
        ["full_siblings"],
        "stated",
        "KIN:008",
    )
    assert (rows[frozenset("bc")]["expected_kinship"], rows[frozenset("bc")]["basis"]) == (0.125, "derived")
    cousins = rows[frozenset("ad")]
    assert (cousins["expected_kinship"], cousins["expected_bands"]) == (None, ["unrelated", "third_degree"])
    assert rows[frozenset("cd")]["expected_bands"] == ["unrelated"]  # 1/32 or 1/64: below the third-degree floor
    assert (rows[frozenset("mn")]["expected_kinship"], rows[frozenset("mn")]["expected_bands"]) == (0.5, ["same_person"])
    assert rows[frozenset("pq")]["expected_bands"] == ["second_degree", "full_siblings"]


def test_a_first_degree_kinship_by_another_path_names_no_kind() -> None:
    """A twin's child is a quarter related to the other twin, by neither a parent link nor two shared
    parents: the pedigree says first degree and not which kind, so either measured kind agrees."""
    rows = _expected(_pedigree("abc", [("a", "KIN:010", "b"), ("b", "KIN:003", "c")]))
    assert (rows[frozenset("bc")]["expected_bands"], rows[frozenset("ac")]["expected_kinship"]) == (["parent_child"], 0.25)
    assert (rows[frozenset("ac")]["expected_bands"], rows[frozenset("ac")]["basis"]) == (["first_degree"], "derived")


def test_a_half_sibling_is_not_a_full_sibling() -> None:
    """One shared parent is half siblings (second degree); two are full siblings."""
    links = [("a", "KIN:027", "c"), ("b", "KIN:028", "c"), ("a", "KIN:027", "d"), ("b", "KIN:028", "d"), ("a", "KIN:027", "e")]
    rows = _expected(_pedigree("abcde", links))
    assert rows[frozenset("cd")]["expected_bands"] == ["full_siblings"]
    assert rows[frozenset("ce")]["expected_bands"] == ["second_degree"]


def test_a_child_of_full_siblings_is_closer_to_each_parent() -> None:
    pedigree = _pedigree("abc", [("a", "KIN:008", "b"), ("a", "KIN:003", "c"), ("b", "KIN:003", "c")])
    # (1/2 + 1/4) / 2: the parent's own kinship and the parents' kinship, halved.
    assert _expected(pedigree)[frozenset("ac")]["expected_kinship"] == 0.375


@pytest.mark.parametrize(
    ("links", "pair", "reason"),
    [
        ([("a", "KIN:002", "b"), ("b", "KIN:003", "c")], "ac", TERM_WITHOUT_DEGREE),  # carried to the child
        ([("a", "KIN:030", "b")], "ab", TERM_WITHOUT_DEGREE),
        ([("a", "KIN:008", "b"), ("a", "KIN:012", "b")], "ab", STATEMENTS_CONFLICT),
        ([("a", "KIN:003", "d"), ("b", "KIN:003", "d"), ("c", "KIN:038", "d")], "ad", MORE_THAN_TWO_PARENTS),
    ],
)
def test_a_pedigree_that_does_not_fix_a_kinship_says_why(links, pair: str, reason: str) -> None:
    row = _expected(_pedigree("abcd", links))[frozenset(pair)]
    assert (row["expected_kinship"], row["expected_bands"], row["expected_reason"]) == (None, None, reason)


def test_social_terms_say_nothing_about_genotypes() -> None:
    rows = _expected(_pedigree("abc", [("a", "KIN:023", "c"), ("b", "KIN:003", "c")]))  # a is c's step parent
    row = rows[frozenset("ac")]
    assert (row["expected_kinship"], row["basis"], row["relation_id"]) == (0.0, "assumed", "KIN:023")


def test_parent_links_that_loop_are_refused() -> None:
    with pytest.raises(PedigreeError, match="loop"):
        expected_kinship(_pedigree("ab", [("a", "KIN:003", "b"), ("b", "KIN:003", "a")]))
