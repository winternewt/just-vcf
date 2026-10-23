"""RM275: the facet graph as a `links` frame, and how it agrees with `facets.derived_from` / `relation`.

Every step that sets a relation is run on a synthetic input, and its `links` are checked against its
`facets` with `link_disagreements`, the rule SCHEMA.md Decision 12 states. That function is then
tested case by case on frames built to break each clause.
"""

import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import LINKS, Loaded, link_disagreements, liftover, load, mint_vrs_ids, normalize, rejoin_split, resolve_ref
from just_vcf.links import (
    APPENDED_EXPORT,
    CALLER_SPLIT,
    KIN_RELATIONS,
    LINEAGE_RELATIONS,
    LINK_KEY,
    LIFTOVER_TO,
    MERGE,
    NORMALIZED,
    RELATIONS,
    SPLIT_REJOINED,
    derived_links,
    kin_links,
    links_frame,
)
from just_vcf.liftover import LIFT_VERSION
from just_vcf.normalize import NORMALIZE_VERSION
from just_vcf.rejoin import REJOIN_VERSION
from just_vcf.schema import FACETS
from test_liftover import lifted_from_cases, setup  # noqa: F401 (a fixture)
from test_load_fbvs import COLUMNS, FB_ONLY, FB_SHARED, META, VS_ONLY
from test_load_keyvalue import ROWS as KEYVALUE_ROWS
from test_load_keyvalue import write as write_keyvalue
from test_load_rsid import _load_text, _two_exports
from test_load_vcf import write as write_vcf_text
from test_rejoin import SPLIT
from test_rejoin import write as write_split


def frames(loaded: Loaded) -> tuple[pl.DataFrame, pl.DataFrame]:
    facets, links = pl.collect_all([loaded.facets, loaded.links])
    assert links.schema == LINKS
    return facets, links


def key(links: pl.DataFrame) -> list[tuple]:
    return list(links.select(LINK_KEY).iter_rows())


def test_the_vocabulary_is_partitioned() -> None:
    assert LINEAGE_RELATIONS.isdisjoint(KIN_RELATIONS) and LINEAGE_RELATIONS | KIN_RELATIONS == RELATIONS
    assert LIFTOVER_TO == "liftover_to", "the Cat-VRS term, verbatim"


def test_a_plain_load_has_no_links(tmp_path: Path) -> None:
    facets, links = frames(load(write_split(tmp_path, SPLIT, "split.vcf")))
    assert links.height == 0
    assert link_disagreements(facets, links) == []


def test_a_caller_split_vcf_links_its_two_facets(tmp_path: Path) -> None:
    text = "\n".join([*META, f"{COLUMNS}\tSample1", FB_SHARED, FB_ONLY, VS_ONLY]) + "\n"
    facets, links = frames(load(write_vcf_text(tmp_path, text, "merge.vcf")))
    ids = sorted(facets["facet_id"].to_list())
    assert len(ids) == 2 and set(facets["relation"]) == {CALLER_SPLIT}
    assert links.to_dicts() == [
        {
            "from_facet_id": ids[0],
            "to_facet_id": ids[1],
            "relation": CALLER_SPLIT,
            "provenance": {"step": "load", "step_version": None, "derivation_inputs": []},
        }
    ]
    assert link_disagreements(facets, links) == []


def test_a_caller_split_table_links_its_two_facets(tmp_path: Path) -> None:
    shuffled = KEYVALUE_ROWS[:]
    random.Random(3).shuffle(shuffled)
    facets, links = frames(load(write_keyvalue(tmp_path, shuffled)))
    assert key(links) == [(*sorted(facets["facet_id"].to_list()), CALLER_SPLIT)]
    assert link_disagreements(facets, links) == []


def test_an_appended_export_links_from_the_first_export(tmp_path: Path) -> None:
    text, _, _ = _two_exports()
    facets, links = frames(_load_text(tmp_path, text, "two.txt"))
    by_relation = {row["relation"]: row["facet_id"] for row in facets.iter_rows(named=True)}
    assert key(links) == [(by_relation[None], by_relation[APPENDED_EXPORT], APPENDED_EXPORT)], "file order: first, then appended"
    assert link_disagreements(facets, links) == []


def test_rejoin_links_the_parent_to_the_derived_facet(tmp_path: Path) -> None:
    split = load(write_split(tmp_path, SPLIT, "split.vcf"))
    rejoined = rejoin_split(split)
    facets, links = frames(rejoined)
    parent = split.facets.collect()["facet_id"].item()
    assert links.to_dicts() == [
        {
            "from_facet_id": parent,
            "to_facet_id": facets["facet_id"].item(),
            "relation": SPLIT_REJOINED,
            "provenance": {"step": "rejoin_split", "step_version": str(REJOIN_VERSION), "derivation_inputs": []},
        }
    ]
    assert link_disagreements(facets, links) == []
    assert rejoin_split(split).links.collect().equals(links), "deterministic"


def test_a_chain_of_steps_carries_every_link(tmp_path: Path, setup) -> None:  # noqa: F811
    """load → resolve_ref → liftover → normalize → mint_vrs_ids: each derived `Loaded` holds every link of
    its ancestry, and resolve_ref and mint_vrs_ids, which keep the facet, pass them through."""
    store, _, _, chain = setup
    resolved, lifted = lifted_from_cases(tmp_path, store, chain)
    assert resolved.links.collect().height == 0
    normalized = normalize(lifted.loaded, store)
    minted = mint_vrs_ids(normalized, store)
    root = resolved.facets.collect()["facet_id"].item()
    lift_id = lifted.loaded.facets.collect()["facet_id"].item()
    facets, links = frames(minted)
    norm_id = facets["facet_id"].item()
    assert key(links) == sorted([(root, lift_id, LIFTOVER_TO), (lift_id, norm_id, NORMALIZED)])
    steps = {row["relation"]: row["provenance"] for row in links.iter_rows(named=True)}
    chain_input = lifted.loaded.facets.collect()["derivation_inputs"].item().to_list()
    assert steps[LIFTOVER_TO] == {"step": "liftover", "step_version": str(LIFT_VERSION), "derivation_inputs": chain_input}
    assert steps[NORMALIZED] == {"step": "normalize", "step_version": str(NORMALIZE_VERSION), "derivation_inputs": []}
    # The facet's own derivation_inputs are carried from its parent; the link's are its step's own.
    assert facets["derivation_inputs"].item().to_list() == chain_input
    assert minted.links.collect().equals(normalized.links.collect())
    # The ancestors' links name facets this Loaded does not hold; the rule checks the facets it holds.
    assert link_disagreements(facets, links) == []
    lifted_facets, lifted_links = frames(lifted.loaded)
    assert key(lifted_links) == [(root, lift_id, LIFTOVER_TO)] and link_disagreements(lifted_facets, lifted_links) == []


def test_resolve_ref_passes_links_through(tmp_path: Path, setup) -> None:  # noqa: F811
    store, source, _, _ = setup
    path = tmp_path / "split.vcf"
    path.write_text("\n".join(["##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1", "1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]) + "\n")
    rejoined = rejoin_split(load(path))
    resolved = resolve_ref(rejoined, store, collection=source)
    assert resolved.links.collect().equals(rejoined.links.collect()) and resolved.links.collect().height == 1


# link_disagreements, clause by clause.

def _facet(facet_id: str, derived_from: str | None = None, relation: str | None = None, source: str = "s") -> dict:
    return {"facet_id": facet_id, "derived_from": derived_from, "relation": relation, "source_sha256": source, "member": None}


def _facets(*rows: dict) -> pl.DataFrame:
    return pl.DataFrame([{name: row.get(name) for name in FACETS.names()} for row in rows], schema=FACETS)


def _link(a: str, b: str, relation: str) -> dict:
    return {"from_facet_id": a, "to_facet_id": b, "relation": relation, "provenance": {"step": "test", "step_version": "1", "derivation_inputs": []}}


def _links(*rows: dict) -> pl.DataFrame:
    return links_frame(list(rows)).collect()


@pytest.mark.parametrize(
    ("facets", "links", "agree"),
    [
        # One parent: the link must exist, come from derived_from and carry the facet's relation.
        ([_facet("c", "p", NORMALIZED)], [_link("p", "c", NORMALIZED)], True),
        ([_facet("c", "p", NORMALIZED)], [], False),
        ([_facet("c", "p", NORMALIZED)], [_link("q", "c", NORMALIZED)], False),
        ([_facet("c", "p", NORMALIZED)], [_link("p", "c", LIFTOVER_TO)], False),
        ([_facet("c", "p", NORMALIZED)], [_link("p", "c", NORMALIZED), _link("q", "c", NORMALIZED)], False),
        ([_facet("c", "p", CALLER_SPLIT)], [_link("p", "c", NORMALIZED)], False),
        ([_facet("c", "p", MERGE)], [_link("p", "c", MERGE)], False),
        # Several parents: derived_from null, two or more merge links.
        ([_facet("m", None, MERGE)], [_link("a", "m", MERGE), _link("b", "m", MERGE)], True),
        ([_facet("m", None, MERGE)], [_link("a", "m", MERGE)], False),
        ([_facet("m", None, MERGE)], [_link("a", "m", MERGE), _link("b", "m", NORMALIZED)], False),
        ([_facet("c", None, NORMALIZED)], [_link("p", "c", NORMALIZED)], False),
        # Kin: a link of the relation names the facet; no lineage link ends at it.
        ([_facet("a", None, CALLER_SPLIT), _facet("b", None, CALLER_SPLIT)], [_link("a", "b", CALLER_SPLIT)], True),
        ([_facet("a", None, CALLER_SPLIT)], [], False),
        ([_facet("a", None, CALLER_SPLIT)], [_link("a", "b", APPENDED_EXPORT)], False),
        ([_facet("a", None, CALLER_SPLIT), _facet("b", None, CALLER_SPLIT, source="t")], [_link("a", "b", CALLER_SPLIT)], False),
        ([_facet("a", None, CALLER_SPLIT)], [_link("a", "b", CALLER_SPLIT), _link("p", "a", NORMALIZED)], False),
        # Neither: no lineage link ends at the facet; kin links from it are fine.
        ([_facet("r")], [], True),
        ([_facet("r")], [_link("r", "x", APPENDED_EXPORT)], True),
        ([_facet("r")], [_link("p", "r", LIFTOVER_TO)], False),
        # Relations outside the documented set, on a facet or a link.
        ([_facet("c", None, "consensus")], [], False),
        ([_facet("r")], [_link("r", "x", "consensus")], False),
    ],
)
def test_link_disagreements(facets: list[dict], links: list[dict], agree: bool) -> None:
    found = link_disagreements(_facets(*facets), _links(*links))
    assert (found == []) is agree, found


def test_derived_links_refuse_a_kin_relation_and_kin_links_a_lineage_one() -> None:
    with pytest.raises(ValueError, match="not a lineage relation"):
        derived_links(links_frame([]), "p", "c", CALLER_SPLIT, step="x", step_version="1")
    with pytest.raises(ValueError, match="not a kin relation"):
        kin_links(["a", "b"], NORMALIZED)


def test_link_order_does_not_depend_on_input_order() -> None:
    rows = [_link(a, b, rel) for a, b, rel in [("p", "c", NORMALIZED), ("a", "m", MERGE), ("b", "m", MERGE), ("a", "b", CALLER_SPLIT), ("c", "d", LIFTOVER_TO)]]
    expected = links_frame(rows).collect()
    assert key(expected) == sorted(key(expected))
    for seed in range(5):
        shuffled = rows[:]
        random.Random(seed).shuffle(shuffled)
        assert links_frame(shuffled).collect().equals(expected)
        parent = links_frame(shuffled)
        added = derived_links(parent, "d", "e", NORMALIZED, step="normalize", step_version="5").collect()
        assert added.equals(derived_links(links_frame(rows), "d", "e", NORMALIZED, step="normalize", step_version="5").collect())
    ids = ["c", "a", "b"]
    assert key(kin_links(ids, CALLER_SPLIT).collect()) == [("a", "b", CALLER_SPLIT), ("a", "c", CALLER_SPLIT), ("b", "c", CALLER_SPLIT)]
    assert key(kin_links(ids, APPENDED_EXPORT).collect()) == [("a", "b", APPENDED_EXPORT), ("c", "a", APPENDED_EXPORT)], "source order, sorted"


def test_a_link_already_carried_is_kept_once() -> None:
    once = derived_links(links_frame([]), "p", "c", NORMALIZED, step="normalize", step_version="5")
    twice = derived_links(once, "p", "c", NORMALIZED, step="normalize", step_version="5").collect()
    assert twice.equals(once.collect())
