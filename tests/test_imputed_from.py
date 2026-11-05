"""RM287: `imputed_from`, the lineage a caller states between an imputed facet and the chip it came from.

just-vcf never infers it: an imputation server makes the imputed file, and only the caller knows which
chip export went in. Stated, the imputed facet and its chip are one callset to `person_karyotype`, one
source to a person's relatedness checks, and a parent and child to `consensus`. Every fixture is
synthetic.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import Person, consensus, normalize, state_imputed_from
from just_vcf.links import IMPUTED_FROM, LIFTOVER_TO, LINEAGE_RELATIONS, MULTI_PARENT_RELATIONS, SINGLE_PARENT_RELATIONS, derived_links, link_disagreements
from just_vcf.loading import Loaded
from just_vcf.person import SubjectError, _callset_keys, _pieces
from just_vcf.schema import FACETS
from just_vcf.store import Store
from test_consensus import TRIANGLE, facet, snv
from test_consensus import store as store  # noqa: F401 - the fixture
from test_karyotype import XY, callset, callsets, made_from, weighed

TOOL = {"tool": "Minimac4", "tool_version": "1.0.2", "panel": "HRC r1.1", "panel_sha256": None}


def stated(imputed: Loaded, chip: Loaded) -> Loaded:
    return state_imputed_from(imputed, chip, **TOOL)


def facet_id(loaded: Loaded) -> str:
    return loaded.facets.select("facet_id").collect()["facet_id"].item()


def test_imputed_from_is_a_one_parent_lineage_relation() -> None:
    assert IMPUTED_FROM == "imputed_from"
    assert IMPUTED_FROM in LINEAGE_RELATIONS and IMPUTED_FROM in SINGLE_PARENT_RELATIONS
    assert IMPUTED_FROM not in MULTI_PARENT_RELATIONS


def test_the_statement_links_the_chip_to_the_imputed_facet_and_changes_no_row() -> None:
    chip, imputed = callset("c", *XY), callset("i", *XY)
    out = stated(imputed, chip)
    facets, links = out.facets.collect(), out.links.collect()
    row = facets.row(0, named=True)
    assert (row["facet_id"], row["derived_from"], row["relation"]) == ("i", "c", IMPUTED_FROM), "the id stays; the lineage is stated"
    assert row["derivation_inputs"] == [{"role": "reference_panel", "name": TOOL["panel"], "sha256": None}]
    assert facets.drop("derived_from", "relation", "derivation_inputs").equals(imputed.facets.collect().drop("derived_from", "relation", "derivation_inputs"))
    assert out.calls.collect().equals(imputed.calls.collect()) and out.blocks.collect().equals(imputed.blocks.collect())
    (link,) = links.iter_rows(named=True)
    assert (link["from_facet_id"], link["to_facet_id"], link["relation"]) == ("c", "i", IMPUTED_FROM)
    assert link["provenance"] == {
        "step": TOOL["tool"],
        "step_version": TOOL["tool_version"],
        "derivation_inputs": [{"role": "reference_panel", "name": TOOL["panel"], "sha256": None}],
    }
    assert link_disagreements(facets, links) == []


def test_an_unknown_panel_is_left_out_not_invented() -> None:
    out = state_imputed_from(callset("i", *XY), callset("c", *XY), tool="Beagle", tool_version=None, panel=None, panel_sha256=None)
    (link,) = out.links.collect().iter_rows(named=True)
    assert link["provenance"] == {"step": "Beagle", "step_version": None, "derivation_inputs": []}
    assert out.facets.collect()["derivation_inputs"].to_list() == [[]]


def test_the_chips_own_lineage_comes_along() -> None:
    """A chip lifted before imputation: the imputed facet's links reach back to the chip as loaded."""
    chip = callset("c", *XY)
    lifted = made_from(chip, "l", *XY)
    out = stated(callset("i", *XY), lifted)
    keys = {(row["from_facet_id"], row["to_facet_id"], row["relation"]) for row in out.links.collect().iter_rows(named=True)}
    assert keys == {("c", "l", "liftover_to"), ("l", "i", IMPUTED_FROM)}


def test_refusals() -> None:
    chip, imputed = callset("c", *XY), callset("i", *XY)
    with pytest.raises(ValueError, match="itself"):
        stated(chip, chip)
    with pytest.raises(ValueError, match="already"):
        stated(made_from(imputed, "l", *XY), chip)
    two = Loaded(facets=pl.concat([chip.facets, imputed.facets]), calls=pl.concat([chip.calls, imputed.calls]), blocks=chip.blocks)
    with pytest.raises(ValueError, match="one facet"):
        stated(two, chip)
    with pytest.raises(ValueError, match="one facet"):
        stated(imputed, two)
    with pytest.raises(SubjectError):
        stated(callset("i", *XY, subject="person-b"), chip)
    with pytest.raises(ValueError, match="made from"):
        stated(imputed, made_from(imputed, "x", *XY))
    with pytest.raises(ValueError, match="tool"):
        state_imputed_from(imputed, chip, tool=" ", tool_version=None, panel=None, panel_sha256=None)


def test_a_chip_and_its_stated_imputation_answer_once_in_person_karyotype() -> None:
    chip, imputed = callset("c", *XY), callset("i", *XY)
    found, _ = weighed(chip, stated(imputed, chip))
    (one,) = found.callsets
    assert (one.root_facet_id, one.members, one.standing_facet_id, one.basis) == ("c", ("c", "i"), "c", "root")
    assert found.supported_by == ("c",)
    # Unstated, the same two facets are two callsets, each a witness.
    found, _ = weighed(chip, imputed)
    assert sorted(callsets(found)) == ["c", "i"] and found.supported_by == ("c", "i")


def checked_pairs(*loaded: Loaded) -> list[frozenset[str]]:
    pairs: list[frozenset[str]] = []

    def check(a: Loaded, b: Loaded) -> bool:
        pairs.append(frozenset(facet_id(one) for one in (a, b)))
        return True

    Person.of(*loaded, check=check)
    return pairs


def test_a_person_checks_a_chip_and_its_stated_imputation_as_one_source() -> None:
    chip, imputed, genome = callset("c", *XY), callset("i", *XY), callset("g", *XY)
    assert checked_pairs(chip, stated(imputed, chip), genome) == [frozenset({"c", "g"})], "the chip stands for both"
    # Unstated, every pair of the three files is checked.
    assert sorted(checked_pairs(chip, imputed, genome), key=sorted) == [frozenset(p) for p in (("c", "g"), ("c", "i"), ("g", "i"))]
    # The imputed facet stands for its callset when the person does not hold the chip.
    assert checked_pairs(stated(imputed, chip), genome) == [frozenset({"g", "i"})]


def test_a_consensus_refuses_a_chip_beside_its_stated_imputation(tmp_path: Path, store) -> None:  # noqa: F811
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1")], store)
    with pytest.raises(ValueError, match="was made from"):
        consensus([x, stated(y, x)], atlas=TRIANGLE)
    assert consensus([x, y], atlas=TRIANGLE).facets.collect().height == 1, "unstated, the two are parents"


def test_a_consensus_refuses_two_facets_made_from_one_ancestor(tmp_path: Path, store) -> None:  # noqa: F811
    """Neither parent made from the other, but both from one chip: its reads would count twice. Two
    derivations of the chip's file (its normalized facet and another step's, here a liftover's links on
    the same rows), its normalized facet beside an imputation stated from it, and two imputations
    stated from it, are each refused naming the chip; unrelated parents are not."""
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1")], store)
    z = facet(tmp_path, "z", "Z", "1", [snv("1", 100, "0/1")], store)
    normalized = normalize(x, store)
    lifted = Loaded(
        facets=normalized.facets.with_columns(facet_id=pl.lit("lifted"), derived_from=pl.lit(facet_id(x)), relation=pl.lit(LIFTOVER_TO)),
        calls=normalized.calls.with_columns(facet_id=pl.lit("lifted")),
        blocks=normalized.blocks.with_columns(facet_id=pl.lit("lifted")),
        links=derived_links(x.links, facet_id(x), "lifted", LIFTOVER_TO, step="liftover", step_version="test"),
    )
    for parents in ([normalized, lifted], [normalized, stated(y, x)], [stated(y, x), stated(z, x)]):
        with pytest.raises(ValueError, match=rf"both made from facet\(s\) {facet_id(x)}\b"):
            consensus(parents, atlas=TRIANGLE)
    # A person reads them as one callset too, the chip's.
    pieces = _pieces((normalized, lifted, stated(y, x), stated(z, x)))
    rows = [piece.facets.collect().row(0, named=True) for piece in pieces]
    links = pl.concat([piece.links for piece in pieces]).collect()
    assert len(set(_callset_keys(rows, links))) == 1
    assert consensus([normalize(x, store), y], atlas=TRIANGLE).facets.collect().height == 1, "one parent derived, the other unrelated"


def test_a_store_keeps_the_stated_link_only_from_the_facets_first_add(tmp_path: Path) -> None:
    chip, imputed = callset("c", *XY), callset("i", *XY)
    first = Store.open(tmp_path / "first", subject="person-a")
    first.add(chip)
    first.add(stated(imputed, chip))
    relations = {(row["from_facet_id"], row["to_facet_id"], row["relation"]) for row in first.loaded().links.collect().iter_rows(named=True)}
    assert relations == {("c", "i", IMPUTED_FROM)}
    # Held already, the facet is present and never rewritten: the statement made after is not kept.
    late = Store.open(tmp_path / "late", subject="person-a")
    late.add(imputed)
    assert late.add(stated(imputed, chip)).present == ("i",)
    assert late.loaded().links.collect().height == 0


def test_the_facets_columns_are_the_schemas() -> None:
    out = stated(callset("i", *XY), callset("c", *XY))
    assert out.facets.collect_schema() == FACETS
