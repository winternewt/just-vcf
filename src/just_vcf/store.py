"""A durable store of one person's facets on disk: parquet per facet, the links, a manifest (RM276).

`Store.open(path, subject=...)` creates the store in a directory the caller names, or opens the one
there; `store.add(loaded)` adds every facet of a `Loaded`; `store.loaded()` and `store.person()` read
everything in it back, lazily. The layout (docs/SCHEMA.md, "The durable store"):

    <path>/manifest.json                 store layout, schema version, subject, created; per facet its
                                         added time and the warnings its load gave, redacted
    <path>/facets/<facet_id>/facets.parquet, calls.parquet, blocks.parquet, links.parquet
    <path>/.lock                         the OS file lock every add holds

The rules, each from the charter or the item:

- **Facets accumulate (charter 5).** A facet's directory is written once, never rewritten or removed.
  Adding a facet the store already holds, by `facet_id`, is a no-op. Adding one never touches another
  facet's files: the links live per facet (the links ending at the facet and at its ancestors, from
  the `Loaded` it came with), so no file is shared between facets. The manifest is the store's index,
  rewritten whole (temporary name, then `os.replace`) after each facet.
- **One person.** The store holds one subject, named when it is created. A facet with no subject or
  another one is refused with `SubjectError`, naming facet ids, never the subject.
- **One schema version.** The store records the `SCHEMA_VERSION` it was created under. A `Loaded` of
  another version is refused, and so is reading a store of another version, with `StoreVersionError`:
  frames of two versions are never mixed, and there is no migration yet (reload the sources into a
  new store).
- **Atomic per facet, safe under concurrent adds.** An add holds the store's lock from its check of
  what is present to its last manifest write. Each facet is sunk into a staging directory, synced, and
  renamed to `facets/<facet_id>/`; the manifest written after the rename is the commit. A staging
  directory, or a facet directory the manifest does not name, found under the lock can only be a
  writer's that died (the lock is held for the whole add), so the next add removes it. `open` reports
  one without the lock, so it may be an add still running.
- **Deterministic (charter 12).** The same `Loaded` sinks to the same parquet bytes under one polars
  release. The manifest carries times and differs.

The store is written only where the caller says; nothing here writes into the cache or the repo, and
nothing logs the subject or a path. The warnings the manifest keeps are redacted as `query` returns
them (`query.FacetRedactor`: the facets' sample, source path and member, and anything shaped like a
path), so the one plain-text file of a store quotes no name.
"""

import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict, ValidationError

from just_vcf.cache import _acquire, _release
from just_vcf.links import LINEAGE_RELATIONS, canonical_links
from just_vcf.loading import SPILL_ROW_GROUP, Loaded
from just_vcf.person import Person, SameSubjectCheck, SubjectError
from just_vcf.query import FacetRedactor
from just_vcf.schema import BLOCKS, CALLS, FACETS, LINKS, SCHEMA_VERSION, empty

logger = logging.getLogger(__name__)

STORE_LAYOUT = 1
"""The store's own layout, recorded in its manifest. A store of a later layout is refused."""

FRAMES = ("facets", "calls", "blocks", "links")
"""The parquet files each facet's directory holds, one per frame."""

_MANIFEST = "manifest.json"
_LOCK = ".lock"
_FACETS_DIR = "facets"
_STAGING = ".staging-"


class StoreError(ValueError):
    """The path is not a store and cannot become one, or the store cannot answer what was asked."""


class StoreVersionError(StoreError):
    """A schema version or store layout other than the store's: refused, never mixed."""


class StoredFacet(BaseModel):
    """One facet in the manifest."""

    model_config = ConfigDict(frozen=True)

    facet_id: str
    added: datetime
    """When the facet was committed, ISO-8601 UTC."""
    warnings: tuple[str, ...] = ()
    """The warnings of the `Loaded` the facet was added from, redacted as `query` returns them."""


class StoreManifest(BaseModel):
    """`manifest.json`: what the store holds. Times are ISO-8601 UTC."""

    model_config = ConfigDict(frozen=True)

    store_layout: int = STORE_LAYOUT
    schema_version: str
    subject: str
    """The caller's label for the person, as `facets.subject` holds it."""
    created: datetime
    facets: tuple[StoredFacet, ...] = ()
    """In `facet_id` order."""


class StoreAdded(BaseModel):
    """What `Store.add` did with each facet of the `Loaded` it was given."""

    model_config = ConfigDict(frozen=True)

    added: tuple[str, ...] = ()
    """Facet ids written by this add."""
    present: tuple[str, ...] = ()
    """Facet ids the store already held: nothing written."""


class Store:
    """One person's facets in a directory: `open`, `add`, `loaded`, `person`.

    Every read goes to the manifest on disk, so a store another process added to answers with what it
    holds now. A frame `loaded` returned keeps working: no file it reads is ever rewritten or removed.
    """

    def __init__(self, path: Path) -> None:
        self.path = path

    @classmethod
    def open(cls, path: Path | str, subject: str | None = None) -> "Store":
        """The store at `path`, created there for `subject` when there is none yet.

        Creating needs `subject` and a directory that is absent or empty: a store is never written
        into a directory holding other files. Opening an existing store with `subject` given checks
        it is the store's (`SubjectError` otherwise, naming no subject). A store of a later layout is
        refused (`StoreVersionError`); one of another schema version opens, so it can be inspected,
        but its frames cannot be read or added to.
        """
        where = Path(path)
        manifest_path = where / _MANIFEST
        if not manifest_path.is_file():
            if subject is None:
                raise StoreError("no store here yet; pass subject=... to create one")
            if where.exists() and (not where.is_dir() or any(child.name != _LOCK for child in where.iterdir())):
                raise StoreError("a new store needs an absent or empty directory; this path holds other files")
            where.mkdir(parents=True, exist_ok=True)
            with _locked(where):
                if not manifest_path.is_file():
                    _write_manifest(
                        where, StoreManifest(schema_version=SCHEMA_VERSION, subject=subject, created=_now())
                    )
        store = cls(where)
        manifest = store.manifest
        if subject is not None and manifest.subject != subject:
            raise SubjectError("the store at this path holds another subject than the one given")
        orphans = store._orphans(manifest)
        if orphans:
            logger.warning(
                "store: %d facet or staging directory(ies) not in the manifest: an add in progress, or one "
                "that did not finish, whose leftovers the next add removes: %s",
                len(orphans),
                ", ".join(orphans),
            )
        return store

    @property
    def manifest(self) -> StoreManifest:
        """The manifest as it is on disk now."""
        return _read_manifest(self.path)

    @property
    def subject(self) -> str:
        return self.manifest.subject

    @property
    def schema_version(self) -> str:
        return self.manifest.schema_version

    @property
    def facet_ids(self) -> tuple[str, ...]:
        """The facets the store holds, in `facet_id` order."""
        return tuple(entry.facet_id for entry in self.manifest.facets)

    def add(self, loaded: Loaded) -> StoreAdded:
        """Add every facet of `loaded` the store does not hold yet; a facet it holds is left as it is.

        Refused, with nothing written: a `Loaded` of another schema version (`StoreVersionError`); a
        facet with no subject or another subject than the store's (`SubjectError`); one `facet_id`
        twice in `loaded`.
        """
        manifest = self.manifest
        _check_version(manifest.schema_version, loaded.schema_version, "the Loaded")
        facets, links = pl.collect_all([loaded.facets, loaded.links])
        ids = facets["facet_id"].to_list()
        twice = sorted({facet_id for facet_id in ids if ids.count(facet_id) > 1})
        if twice:
            raise ValueError(f"facet(s) given more than once: {', '.join(twice)}")
        missing = sorted(row["facet_id"] for row in facets.iter_rows(named=True) if row["subject"] is None)
        if missing:
            raise SubjectError(f"facet(s) without a subject: {', '.join(missing)}; load them with subject=...")
        other = sorted(row["facet_id"] for row in facets.iter_rows(named=True) if row["subject"] != manifest.subject)
        if other:
            raise SubjectError(f"facet(s) {', '.join(other)} carry another subject than the store's")

        redact = FacetRedactor(facets.iter_rows(named=True))
        warnings = tuple(redact(line) for line in loaded.warnings)
        added: list[str] = []
        present: list[str] = []
        with _locked(self.path):
            manifest = self.manifest
            _check_version(manifest.schema_version, loaded.schema_version, "the Loaded")
            for orphan in self._orphans(manifest):
                # The lock is held for a whole add, so what an add left here is a dead writer's.
                logger.info("store: removing %s, left by an add that did not finish", orphan)
                shutil.rmtree(self.path / orphan, ignore_errors=True)
            held = {entry.facet_id for entry in manifest.facets}
            for facet_id in sorted(ids):
                if facet_id in held:
                    present.append(facet_id)
                    continue
                self._write_facet(loaded, facet_id, _links_of(links, facet_id))
                entry = StoredFacet(facet_id=facet_id, added=_now(), warnings=warnings)
                manifest = manifest.model_copy(
                    update={"facets": tuple(sorted((*manifest.facets, entry), key=lambda stored: stored.facet_id))}
                )
                _write_manifest(self.path, manifest)
                added.append(facet_id)
        return StoreAdded(added=tuple(added), present=tuple(present))

    def loaded(self) -> Loaded:
        """Every facet in the store as one `Loaded`: lazy scans of its files, in `facet_id` order, and
        the links of all of them, each link once (Decision 12's order). `warnings` are the facets'
        own, in facet order, each line once. An empty store gives empty frames."""
        manifest = self.manifest
        _check_version(manifest.schema_version, SCHEMA_VERSION, "this just-vcf")
        directories = [self.path / _FACETS_DIR / entry.facet_id for entry in manifest.facets]
        frames = {
            name: pl.concat([pl.scan_parquet(directory / f"{name}.parquet") for directory in directories], how="vertical")
            if directories
            else empty(schema)
            for name, schema in zip(FRAMES, (FACETS, CALLS, BLOCKS, LINKS))
        }
        warnings = list(dict.fromkeys(line for entry in manifest.facets for line in entry.warnings))
        return Loaded(
            facets=frames["facets"],
            calls=frames["calls"],
            blocks=frames["blocks"],
            links=canonical_links(frames["links"]),
            warnings=tuple(warnings),
            schema_version=manifest.schema_version,
        )

    def person(self, check: SameSubjectCheck | None = None) -> Person:
        """`Person.of` over every facet in the store, with `check` as the relatedness test. Its frames
        are the store's lazy scans, and its `links` those `loaded()` reads back. An empty store has no
        person yet: `StoreError`."""
        if not self.facet_ids:
            raise StoreError("the store holds no facets yet")
        return Person.of(self.loaded(), check=check)

    def _write_facet(self, loaded: Loaded, facet_id: str, links: pl.DataFrame) -> None:
        """Sink one facet's four frames into a staging directory, sync them, and rename it into place."""
        stage = self.path / f"{_STAGING}{os.getpid()}-{os.urandom(4).hex()}"
        stage.mkdir()
        try:
            mine = pl.col("facet_id") == facet_id
            pl.collect_all(
                [
                    loaded.facets.filter(mine).sink_parquet(stage / "facets.parquet", lazy=True),
                    loaded.calls.filter(mine).sink_parquet(stage / "calls.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
                    loaded.blocks.filter(mine).sink_parquet(stage / "blocks.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
                    links.lazy().sink_parquet(stage / "links.parquet", lazy=True),
                ]
            )
            for name in FRAMES:
                _fsync(stage / f"{name}.parquet")
            _fsync_dir(stage)
            target = self.path / _FACETS_DIR / facet_id
            target.parent.mkdir(exist_ok=True)
            os.rename(stage, target)
            _fsync_dir(target.parent)
        except BaseException:
            shutil.rmtree(stage, ignore_errors=True)
            raise

    def _orphans(self, manifest: StoreManifest) -> list[str]:
        """Staging directories, and facet directories the manifest does not name, relative to the store."""
        held = {entry.facet_id for entry in manifest.facets}
        found = sorted(child.name for child in self.path.glob(f"{_STAGING}*"))
        facets_dir = self.path / _FACETS_DIR
        if facets_dir.is_dir():
            found += sorted(f"{_FACETS_DIR}/{child.name}" for child in facets_dir.iterdir() if child.name not in held)
        return found


def _links_of(links: pl.DataFrame, facet_id: str) -> pl.DataFrame:
    """The links a facet's directory keeps: those ending at the facet or at one of its lineage
    ancestors, so the facet says how each step made it even when the store holds no ancestor."""
    lineage = links.filter(pl.col("relation").is_in(sorted(LINEAGE_RELATIONS)))
    reached = {facet_id}
    frontier = {facet_id}
    while frontier:
        parents = set(lineage.filter(pl.col("to_facet_id").is_in(sorted(frontier)))["from_facet_id"].to_list())
        frontier = parents - reached
        reached |= parents
    return canonical_links(links.lazy().filter(pl.col("to_facet_id").is_in(sorted(reached)))).collect()


def _check_version(store_version: str, other: str, what: str) -> None:
    if store_version != other:
        raise StoreVersionError(
            f"the store holds schema version {store_version} and {what} is {other}; versions are never mixed. "
            "No migration exists yet: load the sources again into a new store"
        )


class _locked:
    """The store's lock, held for the block (`just_vcf.cache`'s OS file lock)."""

    def __init__(self, path: Path) -> None:
        self.lock_path = path / _LOCK

    def __enter__(self) -> None:
        self.handle = _acquire(self.lock_path, blocking=True)

    def __exit__(self, *_: object) -> None:
        _release(self.handle)


def _read_manifest(path: Path) -> StoreManifest:
    manifest_path = path / _MANIFEST
    if not manifest_path.is_file():
        raise StoreError("no store here: there is no manifest")
    raw = manifest_path.read_bytes()
    try:
        manifest = StoreManifest.model_validate_json(raw)
    except ValidationError as err:
        raise StoreVersionError(
            f"the store's manifest is not one this just-vcf reads ({err.error_count()} error(s)); a later layout?"
        ) from err
    if manifest.store_layout != STORE_LAYOUT:
        raise StoreVersionError(f"the store has layout {manifest.store_layout}; this just-vcf reads layout {STORE_LAYOUT}")
    return manifest


def _write_manifest(path: Path, manifest: StoreManifest) -> None:
    """Written to a temporary name, synced and moved over, so a reader never sees half a manifest."""
    temporary = path / f"{_MANIFEST}.tmp-{os.getpid()}"
    temporary.write_text(manifest.model_dump_json(indent=1) + "\n")
    _fsync(temporary)
    os.replace(temporary, path / _MANIFEST)
    _fsync_dir(path)


def _fsync(path: Path) -> None:
    handle = os.open(path, os.O_RDONLY)
    try:
        os.fsync(handle)
    finally:
        os.close(handle)


def _fsync_dir(path: Path) -> None:
    """Make a rename in `path` durable. Windows cannot open a directory to sync it; there the rename
    is as durable as the filesystem makes it."""
    if os.name == "nt":
        return
    _fsync(path)


def _now() -> datetime:
    return datetime.now(timezone.utc)
