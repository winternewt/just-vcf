"""RM276: migrating a durable store written under an older schema version to this one.

Every old store here is built by hand from docs/SCHEMA.md's history ("Changes since sign-off"), not by
running old code: a store is written now, then each of its files is rewritten as the older version's
loaders wrote it (the columns a later bump added dropped, `phase_set_kinds` cast back to its two
members, the values a later loader fills set as the older one left them), and the manifest names the
older version. The expected migrated frames are the same hand-built frames with the later columns
added null. Every fixture is synthetic; every store is under `tmp_path`.
"""

import json
import logging
import os
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load
from just_vcf.loaders.cgi import FORMATS as CGI_FORMATS
from just_vcf.loading import Loaded
from just_vcf.schema import BLOCKS, CALLS, FACETS, LINKS, SCHEMA_VERSION
from just_vcf.store import FRAMES, Store, StoreError, StoreVersionError, migrate_store
from just_vcf.store_migration import MIGRATIONS, STORE_FIRST_VERSION, migration_chain, version_key
from test_load_cgi import LOCI, var_text
from test_load_cgi import write as write_cgi
from test_load_vcf import vcf_text
from test_load_vcf import write as write_vcf_text
from test_person import _array_file
from test_store import GVCF, SUBJECT, _files

VERSIONS = ("0.0.14", "0.0.15", "0.0.16", "0.0.17", "0.0.18", "0.0.19")
"""Every schema version a store can have been written under: the store shipped under 0.0.14 (merge
6830fb7, 2026-10-23, between 0.0.14's merge 04c5db3 and 0.0.15's 8b72a46)."""

ADDED: dict[str, dict[str, tuple[str, ...]]] = {
    "0.0.15": {"calls": ("phase_set_ordinals", "phase_set_qualities")},
    "0.0.16": {"facets": ("reference", "aligner", "aligner_version", "method_steps")},
    "0.0.18": {"facets": ("gl_source",)},
    "0.0.19": {"calls": ("hq",)},
}
"""The columns each bump since the store shipped added (SCHEMA.md "Changes since sign-off")."""

KINDS_BEFORE_0017 = pl.List(pl.Enum(["ps", "psl"]))
"""`calls.phase_set_kinds` before 0.0.17 added `cgi_haplink`."""

SCHEMAS = dict(zip(FRAMES, (FACETS, CALLS, BLOCKS, LINKS)))

PHASED = [
    "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:100",
    "1\t200\t.\tC\tT\t50\tPASS\t.\tGT:PS\t1|0:100",
    "1\t300\t.\tG\tA\t50\tPASS\t.\tGT:PS\t0/1:.",
]


def _key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def _sources(tmp_path: Path) -> list[Loaded]:
    """A gVCF, a PS-phased VCF, an array export and a CGI `var`: every frame and both caveat kinds."""
    gvcf = vcf_text(GVCF, fmt=("GT", "GQ", "MIN_DP", "PL"), info=("END",))
    phased = vcf_text(PHASED, fmt=("GT", "PS"))
    return [
        load(write_vcf_text(tmp_path, gvcf, "sample.g.vcf"), subject=SUBJECT),
        load(write_vcf_text(tmp_path, phased, "phased.vcf"), subject=SUBJECT),
        load(_array_file(tmp_path), subject=SUBJECT),
        load(write_cgi(tmp_path, var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"), subject=SUBJECT),
    ]


def _facet_ids(sources: list[Loaded]) -> dict[str, str]:
    """facet_id → format_id."""
    rows = [loaded.facets.select("facet_id", "format_id").collect() for loaded in sources]
    return {facet_id: format_id for frame in rows for facet_id, format_id in frame.iter_rows()}


def _phased_id(sources: list[Loaded]) -> str:
    return sources[1].facets.collect()["facet_id"].item()


def _as_left(name: str, frame: pl.DataFrame, version: str, cgi: list[str], psl_dropped: list[str]) -> pl.DataFrame:
    """`frame` (this version's shape) with the values an older loader left: CGI rows unphased before
    0.0.17 read hapLink, and the facets in `psl_dropped` listing FORMAT/PSL among their dropped fields
    (a VCF 4.5 file whose PSL the loaders before 0.0.15 left out)."""
    if name == "calls" and _key(version) < _key("0.0.17"):
        frame = frame.with_columns(
            pl.when(pl.col("facet_id").is_in(cgi)).then(None).otherwise(pl.col(column)).alias(column)
            for column in ("phased", "phase_sets", "phase_set_kinds")
        )
    if name == "facets" and psl_dropped:
        frame = frame.with_columns(
            pl.when(pl.col("facet_id").is_in(psl_dropped))
            .then(pl.concat_list(pl.col("dropped_fields"), pl.lit(["FORMAT/PSL"])))
            .otherwise(pl.col("dropped_fields"))
            .alias("dropped_fields")
        )
    return frame


def _shaped(name: str, frame: pl.DataFrame, version: str) -> pl.DataFrame:
    """`frame` in the shape `version`'s loaders wrote: later columns dropped, the enum narrowed."""
    for bump, frames in ADDED.items():
        if _key(version) < _key(bump):
            frame = frame.drop(frames.get(name, ()))
    if name == "calls" and _key(version) < _key("0.0.17"):
        frame = frame.with_columns(pl.col("phase_set_kinds").cast(KINDS_BEFORE_0017))
    return frame


def _nulled(name: str, frame: pl.DataFrame, version: str) -> pl.DataFrame:
    """`frame` with every column a bump after `version` added set null: what a migration can give."""
    later = [column for bump, frames in ADDED.items() if _key(bump) > _key(version) for column in frames.get(name, ())]
    return frame.with_columns(pl.lit(None, dtype=frame.schema[column]).alias(column) for column in later)


def _old_store(root: Path, version: str, sources: list[Loaded], psl_dropped: list[str] = ()) -> dict[str, pl.DataFrame]:
    """A store at `version`, written as that version would have written it. Returns the frames a
    migration to this version should read back."""
    store = Store.open(root, subject=SUBJECT)
    for loaded in sources:
        store.add(loaded)
    cgi = [facet_id for facet_id, format_id in _facet_ids(sources).items() if format_id in CGI_FORMATS]
    expected: dict[str, list[pl.DataFrame]] = {name: [] for name in FRAMES}
    for facet_id in store.facet_ids:
        for name in FRAMES:
            path = root / "facets" / facet_id / f"{name}.parquet"
            left = _as_left(name, pl.read_parquet(path), version, cgi, list(psl_dropped))
            old = _shaped(name, left, version)
            path.unlink()
            old.write_parquet(path)
            expected[name].append(_nulled(name, left, version))
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({**json.loads(manifest.read_text()), "schema_version": version}))
    return {name: pl.concat(parts) if parts else pl.DataFrame(schema=SCHEMAS[name]) for name, parts in expected.items()}


def _read_back(store: Store) -> dict[str, pl.DataFrame]:
    loaded = store.loaded()
    return dict(zip(FRAMES, pl.collect_all([loaded.facets, loaded.calls, loaded.blocks, loaded.links])))


def _snapshot(root: Path) -> tuple[dict[str, tuple[str, int]], bytes]:
    return _files(root), (root / "manifest.json").read_bytes()


def test_the_chain_runs_from_the_stores_first_version_to_this_one() -> None:
    assert STORE_FIRST_VERSION == VERSIONS[0] and VERSIONS[-1] == SCHEMA_VERSION, "a new bump needs its step here"
    assert [(step.from_version, step.to_version) for step in MIGRATIONS] == list(zip(VERSIONS, VERSIONS[1:]))
    assert migration_chain(STORE_FIRST_VERSION) == MIGRATIONS
    assert migration_chain(SCHEMA_VERSION) == ()
    for name, schema in SCHEMAS.items():
        frame = _shaped(name, pl.DataFrame(schema=schema), STORE_FIRST_VERSION).lazy()
        for step in MIGRATIONS:
            frame = step.apply(name, frame)
        assert frame.collect_schema() == schema, name


def test_each_steps_added_and_retyped_columns_are_what_it_does() -> None:
    for step in MIGRATIONS:
        before = {name: _shaped(name, pl.DataFrame(schema=schema), step.from_version).schema for name, schema in SCHEMAS.items()}
        after = {name: step.apply(name, pl.LazyFrame(schema=before[name])).collect_schema() for name in FRAMES}
        added = [f"{name}.{column}" for name in FRAMES for column in after[name] if column not in before[name]]
        retyped = [
            f"{name}.{column}" for name in FRAMES for column in before[name] if after[name][column] != before[name][column]
        ]
        assert tuple(added) == step.added and tuple(retyped) == step.retyped, step.to_version


@pytest.mark.parametrize("step", MIGRATIONS, ids=lambda step: f"{step.from_version}-{step.to_version}")
def test_each_step_turns_one_versions_frames_into_the_next(tmp_path: Path, step) -> None:
    sources = _sources(tmp_path)
    cgi = [facet_id for facet_id, format_id in _facet_ids(sources).items() if format_id in CGI_FORMATS]
    current = {name: pl.concat([getattr(loaded, name).collect() for loaded in sources]) for name in FRAMES}
    for name, frame in current.items():
        left = _as_left(name, frame, step.from_version, cgi, [])
        old = _shaped(name, left, step.from_version)
        want = _shaped(name, _nulled(name, left, step.from_version), step.to_version)
        got = step.apply(name, old.lazy()).collect()
        assert got.schema == want.schema, name
        assert got.equals(want), name
    assert current["calls"]["phase_set_kinds"].explode(empty_as_null=False).drop_nulls().len() > 0, "the enum cast meets values"


@pytest.mark.parametrize("version", VERSIONS[:-1])
def test_a_store_of_each_older_version_migrates_to_this_one(tmp_path: Path, version: str) -> None:
    root = tmp_path / "store"
    sources = _sources(tmp_path)
    expected = _old_store(root, version, sources)
    store = Store.open(root)
    assert store.schema_version == version
    with pytest.raises(StoreVersionError) as refused:
        store.loaded()
    message = str(refused.value)
    assert "migrate_store" in message and version in message and SCHEMA_VERSION in message
    with pytest.raises(StoreVersionError) as refused:
        store.add(sources[0])
    assert "migrate_store" in str(refused.value)

    report = migrate_store(root)
    assert (report.from_version, report.to_version) == (version, SCHEMA_VERSION)
    chain = migration_chain(version)
    assert [(step.from_version, step.to_version) for step in report.steps] == [(s.from_version, s.to_version) for s in chain]
    assert report.facets == store.facet_ids and len(report.facets) == len(sources)
    assert store.schema_version == SCHEMA_VERSION
    got = _read_back(store)
    for name in FRAMES:
        assert got[name].schema == SCHEMAS[name], name
        assert got[name].equals(expected[name]), name
    assert not list(root.glob(".migrating-*")) and not list(root.glob(".previous-*")), "the old files go once verified"
    assert store.manifest.migrations == (report,)
    held = sources[0].facets.collect()["facet_id"].item()
    assert store.add(sources[0]).present == (held,), "the migrated store takes adds again"

    again = migrate_store(root)
    assert again.steps == () and again.from_version == SCHEMA_VERSION, "a current store is left as it is"


def test_the_report_names_what_was_filled_unknown_and_loaded_repeats_it(tmp_path: Path) -> None:
    root = tmp_path / "store"
    sources = _sources(tmp_path)
    formats = _facet_ids(sources)
    phased = _phased_id(sources)
    _old_store(root, STORE_FIRST_VERSION, sources, psl_dropped=[phased])
    report = migrate_store(root)

    vcf = tuple(sorted(facet_id for facet_id, format_id in formats.items() if format_id == "vcf"))
    cgi = tuple(sorted(facet_id for facet_id, format_id in formats.items() if format_id in CGI_FORMATS))
    assert len(vcf) == 2 and len(cgi) == 1
    filled = [(entry.version, entry.frame, entry.facet_ids) for entry in report.filled_unknown]
    assert filled == [
        ("0.0.15", "calls", (phased,)),
        ("0.0.16", "facets", vcf),
        ("0.0.17", "calls", cgi),
        ("0.0.18", "facets", vcf),
        ("0.0.19", "calls", cgi),
    ], "the array facet loses nothing; no masterVarBeta, so no read-count entry"
    assert all("load the source again" in entry.reason for entry in report.filled_unknown)

    warnings = Store.open(root).loaded().warnings
    for entry in report.filled_unknown:
        line = next(line for line in warnings if f"step to {entry.version}, {entry.frame}" in line)
        assert all(facet_id in line for facet_id in entry.facet_ids) and entry.reason in line
    assert SUBJECT not in json.dumps([*warnings, report.model_dump_json()]), "facet ids, never the subject"


def test_a_newer_store_is_refused_never_downgraded(tmp_path: Path) -> None:
    root = tmp_path / "store"
    Store.open(root, subject=SUBJECT).add(_sources(tmp_path)[0])
    newer = ".".join([*SCHEMA_VERSION.split(".")[:-1], str(int(SCHEMA_VERSION.split(".")[-1]) + 1)])
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({**json.loads(manifest.read_text()), "schema_version": newer}))
    before = _snapshot(root)
    with pytest.raises(StoreVersionError) as refused:
        Store.open(root).loaded()
    assert "upgrade just-vcf" in str(refused.value) and newer in str(refused.value)
    with pytest.raises(StoreVersionError) as refused:
        migrate_store(root)
    assert "upgrade just-vcf" in str(refused.value)
    assert _snapshot(root) == before


@pytest.mark.parametrize("version", ["0.0.1", "0.0.9", "0.0.13", "x.y"])
def test_a_version_no_store_was_written_under_has_no_migration(tmp_path: Path, version: str) -> None:
    root = tmp_path / "store"
    Store.open(root, subject=SUBJECT).add(_sources(tmp_path)[0])
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({**json.loads(manifest.read_text()), "schema_version": version}))
    assert migration_chain(version) is None
    for attempt in (lambda: Store.open(root).loaded(), lambda: migrate_store(root)):
        with pytest.raises(StoreVersionError) as refused:
            attempt()
        message = str(refused.value)
        assert "no migration" in message and STORE_FIRST_VERSION in message
        assert "upgrade" not in message, "0.0.9 is older than 0.0.19, not newer"
    assert version_key("0.0.9") < version_key("0.0.19")


def _failing(monkeypatch: pytest.MonkeyPatch, frame: str, broken: Callable[[pl.LazyFrame], pl.LazyFrame]) -> None:
    step = MIGRATIONS[-1]
    original = step.frames.get(frame, lambda lazy: lazy)
    monkeypatch.setitem(step.frames, frame, lambda lazy: broken(original(lazy)))


def test_a_migration_that_fails_its_check_or_raises_leaves_the_store_as_it_was(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "store"
    expected = _old_store(root, STORE_FIRST_VERSION, _sources(tmp_path))
    before = _snapshot(root)

    with monkeypatch.context() as patch:
        _failing(patch, "calls", lambda lazy: lazy.slice(1))
        with pytest.raises(StoreError) as refused:
            migrate_store(root)
        assert "rows" in str(refused.value)
    assert _snapshot(root) == before and not list(root.glob(".migrating-*"))

    with monkeypatch.context() as patch:
        _failing(patch, "calls", lambda lazy: lazy.with_columns(pl.col("pos") + 1))
        with pytest.raises(StoreError) as refused:
            migrate_store(root)
        assert "values of its calls frame changed" in str(refused.value)
    assert _snapshot(root) == before and not list(root.glob(".migrating-*"))

    seen: list[int] = []

    def raise_on_the_second(lazy: pl.LazyFrame) -> pl.LazyFrame:
        seen.append(1)
        if len(seen) == 2:
            raise RuntimeError("a step that fails midway")
        return lazy

    with monkeypatch.context() as patch:
        _failing(patch, "facets", raise_on_the_second)
        with pytest.raises(RuntimeError):
            migrate_store(root)
    assert _snapshot(root) == before and not list(root.glob(".migrating-*"))

    migrate_store(root)
    got = _read_back(Store.open(root))
    assert all(got[name].equals(expected[name]) for name in FRAMES)


@pytest.mark.parametrize("call", ["rename", "replace"])
def test_a_migration_cut_off_after_its_check_finishes_on_the_next_call(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, call: str
) -> None:
    """Cut off between parking the old facets and moving the new ones in (`rename`), or before the
    manifest is committed (`replace`): the store refuses reads at the old version until the next
    `migrate_store` finishes the swap."""
    root = tmp_path / "store"
    expected = _old_store(root, STORE_FIRST_VERSION, _sources(tmp_path))
    real = getattr(os, call)

    def cut(source: object, target: object, *args: object, **kwargs: object) -> None:
        moving_in = call == "rename" and Path(str(target)) == root / "facets"
        committing = call == "replace" and Path(str(target)) == root / "manifest.json"
        if moving_in or committing:
            raise OSError("cut off")
        real(source, target, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(os, call, cut)
        with pytest.raises(OSError):
            migrate_store(root)
    with caplog.at_level(logging.WARNING, logger="just_vcf.store"):
        store = Store.open(root)
    assert "migration" in caplog.text
    assert store.schema_version == STORE_FIRST_VERSION
    with pytest.raises(StoreVersionError):
        store.loaded()

    report = migrate_store(root)
    assert report.from_version == STORE_FIRST_VERSION and report.to_version == SCHEMA_VERSION
    got = _read_back(store)
    assert all(got[name].equals(expected[name]) for name in FRAMES)
    assert not list(root.glob(".migrating-*")) and not list(root.glob(".previous-*"))


def test_keep_previous_keeps_the_old_versions_files(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    root = tmp_path / "store"
    _old_store(root, STORE_FIRST_VERSION, _sources(tmp_path))
    old_files, old_manifest = _snapshot(root)
    migrate_store(root, keep_previous=True)
    previous = root / f".previous-{STORE_FIRST_VERSION}"
    assert {name: digest for name, (digest, _) in _files(previous).items()} == {
        name: digest for name, (digest, _) in old_files.items()
    }
    assert (previous / "manifest.json").read_bytes() == old_manifest
    with caplog.at_level(logging.WARNING, logger="just_vcf.store"):
        Store.open(root)
    assert "kept" in caplog.text


def test_an_empty_old_store_migrates(tmp_path: Path) -> None:
    root = tmp_path / "store"
    _old_store(root, STORE_FIRST_VERSION, [])
    report = migrate_store(root)
    assert report.facets == () and report.filled_unknown == () and len(report.steps) == len(MIGRATIONS)
    got = _read_back(Store.open(root))
    assert all(got[name].height == 0 and got[name].schema == SCHEMAS[name] for name in FRAMES)
