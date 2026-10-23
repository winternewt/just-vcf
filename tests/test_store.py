"""RM276: the durable per-person store. Load → add → reopen gives the same frames; facets accumulate.

Every fixture is synthetic; the store is written under `tmp_path`.
"""

import hashlib
import json
import logging
import threading
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, rejoin_split
from just_vcf.links import CALLER_SPLIT, SPLIT_REJOINED, link_disagreements
from just_vcf.loading import Loaded
from just_vcf.person import Person, SubjectError
from just_vcf.schema import BLOCKS, CALLS, FACETS, LINKS, SCHEMA_VERSION
from just_vcf.store import FRAMES, Store, StoreError, StoreVersionError
from test_load_fbvs import COLUMNS, FB_ONLY, FB_SHARED, META, VS_ONLY
from test_load_vcf import vcf_text
from test_load_vcf import write as write_vcf_text
from test_person import _array_file, _vcf_file
from test_rejoin import SPLIT
from test_rejoin import write as write_split

SUBJECT = "person-a"

GVCF = [
    "1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ:MIN_DP:PL\t0/0:60:23:0,60,900",
    "1\t114\t.\tT\tTC,<NON_REF>\t213\t.\t.\tGT:GQ:MIN_DP:PL\t0/1:99:.:51,0,36,93,92,86",
    "1\t116\t.\tC\t<*>\t.\t.\tEND=120\tGT:GQ:MIN_DP:PL\t./.:.:.:.",
]


def _gvcf(tmp_path: Path, subject: str | None = SUBJECT) -> Loaded:
    text = vcf_text(GVCF, fmt=("GT", "GQ", "MIN_DP", "PL"), info=("END",))
    return load(write_vcf_text(tmp_path, text, "sample.g.vcf"), subject=subject)


def _caller_split(tmp_path: Path) -> Loaded:
    text = "\n".join([*META, f"{COLUMNS}\tSample1", FB_SHARED, FB_ONLY, VS_ONLY]) + "\n"
    return load(write_vcf_text(tmp_path, text, "merge.vcf"), subject=SUBJECT)


def _all(loaded: Loaded) -> dict[str, pl.DataFrame]:
    return dict(zip(FRAMES, pl.collect_all([loaded.facets, loaded.calls, loaded.blocks, loaded.links])))


def _files(root: Path) -> dict[str, tuple[str, int]]:
    """Every file under `facets/`: its bytes' sha256 and its inode, so a rewrite in place shows too."""
    return {
        path.relative_to(root).as_posix(): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_ino)
        for path in sorted((root / "facets").rglob("*"))
        if path.is_file()
    }


def test_load_add_reopen_gives_the_same_frames(tmp_path: Path) -> None:
    sources = [_gvcf(tmp_path), load(_array_file(tmp_path), subject=SUBJECT), _caller_split(tmp_path)]
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    for loaded in sources:
        store.add(loaded)

    reopened = Store.open(tmp_path / "store")
    back = reopened.loaded()
    assert isinstance(back.calls, pl.LazyFrame), "reads are lazy"
    got = _all(back)
    assert got["facets"].schema == FACETS and got["calls"].schema == CALLS
    assert got["blocks"].schema == BLOCKS and got["links"].schema == LINKS
    expected = Person.of(*sources)
    want_facets, want_calls, want_blocks = pl.collect_all([expected.facets, expected.calls, expected.blocks])
    assert got["facets"].equals(want_facets)
    assert got["calls"].equals(want_calls)
    assert got["blocks"].equals(want_blocks) and got["blocks"].height > 0, "the gVCF's blocks come back"
    assert got["links"].equals(_all(sources[2])["links"]), "the caller split's kin link, once"
    assert reopened.facet_ids == tuple(sorted(got["facets"]["facet_id"].to_list()))

    person = reopened.person()
    assert person.subject == SUBJECT
    assert set(person.checks["outcome"]) == {"not_run"}, "no check given: not a pass"
    assert person.calls.collect().equals(want_calls)


def test_the_manifest_records_versions_subject_and_utc_times(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    store.add(loaded)
    raw = json.loads((tmp_path / "store" / "manifest.json").read_text())
    assert raw["schema_version"] == SCHEMA_VERSION and raw["subject"] == SUBJECT and raw["store_layout"] == 1
    assert [entry["facet_id"] for entry in raw["facets"]] == loaded.facets.collect()["facet_id"].to_list()
    for stamp in (raw["created"], raw["facets"][0]["added"]):
        assert stamp.endswith("Z") or stamp.endswith("+00:00"), "ISO-8601 UTC"
    assert tuple(raw["facets"][0]["warnings"]) == loaded.warnings


def test_adding_a_facet_already_held_is_a_no_op(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    first = store.add(loaded)
    before, manifest = _files(tmp_path / "store"), store.manifest
    again = store.add(_gvcf(tmp_path))
    assert first.added == again.present and again.added == ()
    assert _files(tmp_path / "store") == before
    assert store.manifest == manifest, "the added time is the first add's"


def test_adding_a_facet_never_rewrites_another(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    store.add(_gvcf(tmp_path))
    before = _files(tmp_path / "store")
    held = store.loaded()
    store.add(load(_array_file(tmp_path), subject=SUBJECT))
    after = _files(tmp_path / "store")
    assert {name: after[name] for name in before} == before
    assert len(after) == len(before) + len(FRAMES)
    assert held.calls.collect().height > 0, "a frame read before the add still reads"


def test_a_derived_facet_keeps_its_lineage_with_or_without_its_parent(tmp_path: Path) -> None:
    split = load(write_split(tmp_path, SPLIT, "split.vcf"), subject=SUBJECT)
    rejoined = rejoin_split(split)
    parent = split.facets.collect()["facet_id"].item()

    alone = Store.open(tmp_path / "alone", subject=SUBJECT)
    alone.add(rejoined)
    links = alone.loaded().links.collect()
    assert links.select("from_facet_id", "relation").rows() == [(parent, SPLIT_REJOINED)], "the parent is named though not held"

    both = Store.open(tmp_path / "both", subject=SUBJECT)
    both.add(rejoined)
    both.add(split)
    facets, links = pl.collect_all([both.loaded().facets, both.loaded().links])
    assert links.equals(rejoined.links.collect()), "each link once"
    assert link_disagreements(facets, links) == []


def test_a_multi_facet_load_adds_each_facet_with_its_links(tmp_path: Path) -> None:
    loaded = _caller_split(tmp_path)
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    result = store.add(loaded)
    assert len(result.added) == 2
    facets, links = pl.collect_all([store.loaded().facets, store.loaded().links])
    assert set(facets["relation"]) == {CALLER_SPLIT} and links.height == 1
    assert link_disagreements(facets, links) == []


def test_the_same_input_writes_the_same_bytes(tmp_path: Path) -> None:
    loaded = _caller_split(tmp_path)
    for name in ("one", "two"):
        Store.open(tmp_path / name, subject=SUBJECT).add(loaded)
    one, two = _files(tmp_path / "one"), _files(tmp_path / "two")
    assert one.keys() == two.keys()
    assert {name: digest for name, (digest, _) in one.items()} == {name: digest for name, (digest, _) in two.items()}


def test_another_schema_version_is_refused_never_mixed(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    with pytest.raises(StoreVersionError):
        store.add(loaded.model_copy(update={"schema_version": "0.0.1"}))
    assert store.facet_ids == () and not (tmp_path / "store" / "facets").exists()

    store.add(loaded)
    manifest = tmp_path / "store" / "manifest.json"
    raw = json.loads(manifest.read_text())
    manifest.write_text(json.dumps({**raw, "schema_version": "0.0.1"}))
    reopened = Store.open(tmp_path / "store")
    assert reopened.schema_version == "0.0.1", "an old store opens for inspection"
    with pytest.raises(StoreVersionError):
        reopened.loaded()
    with pytest.raises(StoreVersionError):
        reopened.add(_gvcf(tmp_path))


def test_a_later_store_layout_is_refused(tmp_path: Path) -> None:
    Store.open(tmp_path / "store", subject=SUBJECT)
    manifest = tmp_path / "store" / "manifest.json"
    manifest.write_text(json.dumps({**json.loads(manifest.read_text()), "store_layout": 2}))
    with pytest.raises(StoreVersionError):
        Store.open(tmp_path / "store")


def test_the_store_holds_one_subject(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    for loaded in (_gvcf(tmp_path, subject=None), _gvcf(tmp_path, subject="person-b")):
        with pytest.raises(SubjectError) as refused:
            store.add(loaded)
        assert "person-b" not in str(refused.value) and SUBJECT not in str(refused.value), "facet ids, never the subject"
    assert store.facet_ids == ()
    with pytest.raises(SubjectError) as refused:
        Store.open(tmp_path / "store", subject="person-b")
    assert "person-b" not in str(refused.value)
    assert Store.open(tmp_path / "store", subject=SUBJECT).subject == SUBJECT


def test_a_store_is_created_only_where_asked_and_only_in_an_empty_directory(tmp_path: Path) -> None:
    with pytest.raises(StoreError):
        Store.open(tmp_path / "nothing")
    assert not (tmp_path / "nothing").exists()
    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "notes.txt").write_text("x")
    with pytest.raises(StoreError):
        Store.open(busy, subject=SUBJECT)
    assert sorted(child.name for child in busy.iterdir()) == ["notes.txt"]
    (tmp_path / "empty").mkdir()
    assert Store.open(tmp_path / "empty", subject=SUBJECT).facet_ids == ()


def test_an_empty_store_reads_empty_frames_and_has_no_person(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    got = _all(store.loaded())
    assert all(frame.height == 0 for frame in got.values())
    assert got["calls"].schema == CALLS
    with pytest.raises(StoreError):
        store.person()


def test_concurrent_adds_all_land_once(tmp_path: Path) -> None:
    sources = [_gvcf(tmp_path), load(_array_file(tmp_path), subject=SUBJECT), _caller_split(tmp_path)]
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    errors: list[BaseException] = []

    def add(loaded: Loaded) -> None:
        try:
            Store.open(tmp_path / "store").add(loaded)
        except BaseException as err:  # noqa: BLE001 (handed to the main thread)
            errors.append(err)

    threads = [threading.Thread(target=add, args=(loaded,)) for loaded in sources * 3]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert errors == []
    expected = sorted(facet_id for loaded in sources for facet_id in loaded.facets.collect()["facet_id"].to_list())
    assert list(store.facet_ids) == expected
    assert sorted(child.name for child in (tmp_path / "store" / "facets").iterdir()) == expected
    assert not list((tmp_path / "store").glob(".staging-*"))


def test_what_a_dead_add_left_is_reported_then_removed(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    root = tmp_path / "store"
    Store.open(root, subject=SUBJECT)
    (root / ".staging-1-dead").mkdir()
    (root / "facets" / "0").mkdir(parents=True)
    with caplog.at_level(logging.WARNING, logger="just_vcf.store"):
        store = Store.open(root)
    assert "2 facet or staging" in caplog.text
    store.add(_gvcf(tmp_path))
    assert not (root / ".staging-1-dead").exists() and not (root / "facets" / "0").exists()
    assert len(store.facet_ids) == 1
