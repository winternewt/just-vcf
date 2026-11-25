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

from just_vcf import ReferenceStore, load, mint_vrs_ids, normalize, rejoin_split, resolve_ref
from just_vcf import store as store_module
from just_vcf.links import CALLER_SPLIT, SPLIT_REJOINED, link_disagreements
from just_vcf.loading import Loaded
from just_vcf.person import Person, SubjectError
from just_vcf.schema import BLOCKS, CALLS, FACETS, LINKS, SCHEMA_VERSION
from just_vcf.reference import RESOLVE_VERSION
from just_vcf.store import FRAMES, Store, StoreConflictError, StoreError, StoreVersionError
from test_insertion_points import point_loci
from test_load_cgi import var_text
from test_load_cgi import write as write_cgi
from test_load_fbvs import COLUMNS, FB_ONLY, FB_SHARED, META, VS_ONLY
from test_load_vcf import vcf_text
from test_load_vcf import write as write_vcf_text
from test_person import _array_file, _vcf_file
from test_rejoin import SPLIT
from test_reference import SEQ, write_fasta
from test_rejoin import write as write_split
from test_spill_streams import fallback_nodes

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


def test_the_manifest_keeps_the_warnings_redacted_as_query_returns_them(tmp_path: Path) -> None:
    """The manifest is plain text beside the parquet: a warning quoting the sample, the file or a path
    is written as `query` would return it."""
    folder = tmp_path / "PlantedOwner"
    folder.mkdir()
    text = vcf_text(["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"], samples=("PlantedSample",))
    loaded = load(write_vcf_text(folder, text, "PlantedFile.vcf"), subject=SUBJECT)
    quoting = loaded.model_copy(update={"warnings": (
        *loaded.warnings,
        "2 row(s) of PlantedSample left out",
        "PlantedFile.vcf: read through /home/PlantedOwner/cache/copy.vcf",
    )})  # fmt: skip
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    store.add(quoting)
    raw = (tmp_path / "store" / "manifest.json").read_text()
    assert "Planted" not in raw
    stored = json.loads(raw)["facets"][0]["warnings"]
    assert stored[-2:] == ["2 row(s) of <sample> left out", "<file name>: read through <path>"]
    assert store.loaded().warnings[-2:] == tuple(stored[-2:])


def test_adding_a_facet_already_held_is_a_no_op(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    first = store.add(loaded)
    before, manifest = _files(tmp_path / "store"), store.manifest
    again = store.add(_gvcf(tmp_path))
    assert first.added == again.present and again.added == ()
    assert _files(tmp_path / "store") == before
    assert store.manifest == manifest, "the added time is the first add's"


def _resolved_cgi(tmp_path: Path) -> tuple[Loaded, Loaded, ReferenceStore]:
    reference = ReferenceStore.open(tmp_path / "grch37")
    reference.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    loaded = load(write_cgi(tmp_path, var_text(point_loci()), "var-GS000000001-ASM.tsv.bz2"), subject=SUBJECT)
    return loaded, resolve_ref(loaded, reference), reference


def test_a_facet_held_in_another_ref_state_is_refused_not_reported_present(tmp_path: Path) -> None:
    """`resolve_ref` keeps a facet's `facet_id` (SCHEMA.md, "How RM3 fills these slots"), so the id
    alone does not say whether the store holds the facet as loaded, or with REF read by an older
    `RESOLVE_VERSION` than this one: an add reported `present` and kept the rows it held. The facet now
    records the rule that read its REF (`resolve_version`, inherited by what is derived from it), and an
    add that meets a held facet in another REF state is refused, nothing written, naming the facet. The
    same state again is a no-op, as before."""
    loaded, resolved, reference = _resolved_cgi(tmp_path)
    assert resolved.facets.collect()["resolve_version"].to_list() == [RESOLVE_VERSION]
    assert loaded.facets.collect()["resolve_version"].to_list() == [None]
    normalized = normalize(resolved, reference)
    assert normalized.facets.collect()["resolve_version"].to_list() == [RESOLVE_VERSION]
    facet_id = resolved.facets.collect()["facet_id"].item()

    store = Store.open(tmp_path / "store", subject=SUBJECT)
    older = resolved.model_copy(update={"facets": resolved.facets.with_columns(pl.lit(RESOLVE_VERSION - 1, dtype=pl.UInt32).alias("resolve_version"))})
    assert store.add(older).added == (facet_id,)
    before, manifest = _files(tmp_path / "store"), store.manifest
    for other in (resolved, loaded):
        # REF is checked before the rows, so the refusal names the cause, not only that rows differ.
        with pytest.raises(StoreConflictError, match=f"held with REF in another state: {facet_id}"):
            store.add(other)
        assert _files(tmp_path / "store") == before and store.manifest == manifest, "nothing written"
    assert store.add(older).present == (facet_id,)

    fresh = Store.open(tmp_path / "fresh", subject=SUBJECT)
    fresh.add(loaded)
    with pytest.raises(StoreConflictError, match="REF not read"):
        fresh.add(resolved)
    # A Loaded holding a held facet in another state and a new one is refused whole.
    both = Person.of(resolved, normalized, check=None)
    with pytest.raises(StoreConflictError):
        fresh.add(Loaded(facets=both.facets, calls=both.calls, blocks=both.blocks, links=both.links, schema_version=both.schema_version))
    assert fresh.facet_ids == (facet_id,)


def _one_row_changed(frame: pl.LazyFrame, column: str) -> pl.LazyFrame:
    """`frame` with `column` of its first row changed (a null read as 0, plus one); every other row kept."""
    first = pl.int_range(pl.len()) == 0
    return frame.with_columns(pl.when(first).then(pl.col(column).fill_null(0) + 1).otherwise(pl.col(column)).alias(column))


def _refused_unwritten(store: Store, loaded: Loaded, facet_id: str, frames: set[str]) -> None:
    """`loaded` re-adds `facet_id` with other rows: refused, naming the facet and exactly `frames`, and
    nothing in the store changed."""
    before, manifest = _files(store.path), store.manifest
    with pytest.raises(StoreConflictError, match=f"held with other rows: {facet_id}") as refused:
        store.add(loaded)
    named = {name for name in FRAMES if f"{name}: " in str(refused.value)}
    assert named == frames
    assert "/" not in str(refused.value).split(". A held facet")[0], "no path in the message"
    assert _files(store.path) == before and store.manifest == manifest, "nothing written"


def test_a_re_add_whose_rows_are_the_held_ones_is_a_no_op(tmp_path: Path) -> None:
    """RM305: the rows are compared, so a re-add of the same content stays a no-op. The same file read
    from another directory differs only in `source_path`, which the comparison leaves out; the store's
    own frames read back are the held rows, links included."""
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    store.add(_gvcf(tmp_path))
    store.add(rejoin_split(load(write_split(tmp_path, SPLIT, "split.vcf"), subject=SUBJECT)))
    before, manifest = _files(tmp_path / "store"), store.manifest
    (tmp_path / "elsewhere").mkdir()
    moved = _gvcf(tmp_path / "elsewhere")
    facet_id = moved.facets.collect()["facet_id"].item()
    held = store.loaded().facets.filter(pl.col("facet_id") == facet_id).collect()
    assert moved.facets.collect()["source_path"].item() != held["source_path"].item()
    assert store.add(moved).present == (facet_id,)
    assert store.add(store.loaded()).present == store.facet_ids
    assert _files(tmp_path / "store") == before and store.manifest == manifest


def test_a_re_add_with_one_call_row_changed_is_refused(tmp_path: Path) -> None:
    """RM305: before, a facet re-added under the same `facet_id` with other rows (a loader fixed, any
    step that keeps the id) was reported `present` and the store kept its old rows without a word."""
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    facet_id = loaded.facets.collect()["facet_id"].item()
    store.add(loaded)
    assert loaded.calls.collect().height > 0
    changed = loaded.model_copy(update={"calls": _one_row_changed(loaded.calls, "qual")})
    _refused_unwritten(store, changed, facet_id, {"calls"})
    # A Loaded holding the changed facet and a new one is refused whole.
    both = Person.of(changed, load(_array_file(tmp_path), subject=SUBJECT), check=None)
    combined = Loaded(facets=both.facets, calls=both.calls, blocks=both.blocks, links=both.links, schema_version=both.schema_version)
    _refused_unwritten(store, combined, facet_id, {"calls"})
    assert store.facet_ids == (facet_id,)
    dropped = loaded.model_copy(update={"calls": loaded.calls.filter(pl.int_range(pl.len()) > 0)})
    _refused_unwritten(store, dropped, facet_id, {"calls"})


def test_a_re_add_differing_only_in_blocks_is_refused(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    gvcf = _gvcf(tmp_path)
    store.add(gvcf)
    assert gvcf.blocks.collect().height > 0
    _refused_unwritten(store, gvcf.model_copy(update={"blocks": _one_row_changed(gvcf.blocks, "gq")}), gvcf.facets.collect()["facet_id"].item(), {"blocks"})


def test_a_re_add_differing_only_in_links_is_refused(tmp_path: Path) -> None:
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    rejoined = rejoin_split(load(write_split(tmp_path, SPLIT, "split.vcf"), subject=SUBJECT))
    store.add(rejoined)
    assert rejoined.links.collect().height == 1
    relinked = rejoined.links.with_columns(pl.col("provenance").struct.with_fields(pl.lit("0").alias("step_version")))
    _refused_unwritten(store, rejoined.model_copy(update={"links": relinked}), rejoined.facets.collect()["facet_id"].item(), {"links"})


def test_a_re_add_after_minting_vrs_ids_is_refused(tmp_path: Path) -> None:
    """`mint_vrs_ids` keeps the facet (its `facets` frame as given) and fills `calls.vrs_allele_ids`: the
    store holding the facet resolved but unminted refuses the minted one, which goes into a new store."""
    _, resolved, reference = _resolved_cgi(tmp_path)
    minted = mint_vrs_ids(resolved, reference)
    facet_id = resolved.facets.collect()["facet_id"].item()
    assert minted.facets.collect().equals(resolved.facets.collect()), "minting keeps the facet and its id"
    assert minted.calls.collect()["vrs_allele_ids"].list.drop_nulls().list.len().sum() > 0, "ids were minted"
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    store.add(resolved)
    _refused_unwritten(store, minted, facet_id, {"calls"})
    assert Store.open(tmp_path / "fresh", subject=SUBJECT).add(minted).added == (facet_id,)


def test_the_re_add_comparison_streams(tmp_path: Path) -> None:
    """Both sides of the comparison are one lazy aggregate per frame that polars' streaming engine runs
    with no node falling back to the in-memory engine, so a gVCF's hundreds of millions of rows are never
    held whole."""
    store = Store.open(tmp_path / "store", subject=SUBJECT)
    loaded = _gvcf(tmp_path)
    store.add(loaded)
    facet_id = loaded.facets.collect()["facet_id"].item()
    given = store_module._facet_frames(loaded, facet_id, store_module._links_of(loaded.links.collect(), facet_id))
    for name in FRAMES:
        columns = store_module._compared(name)
        held = pl.scan_parquet(tmp_path / "store" / "facets" / facet_id / f"{name}.parquet")
        assert fallback_nodes(store_module._digest(held, columns)) == []
        assert fallback_nodes(store_module._digest(given[name], columns)) == []


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
    assert "2 facet or staging directory(ies) not in the manifest" in caplog.text
    store.add(_gvcf(tmp_path))
    assert not (root / ".staging-1-dead").exists() and not (root / "facets" / "0").exists()
    assert len(store.facet_ids) == 1
