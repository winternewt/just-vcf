"""A durable store of one person's facets on disk: parquet per facet, the links, a manifest (RM276).

`Store.open(path, subject=...)` creates the store in a directory the caller names, or opens the one
there; `store.add(loaded)` adds every facet of a `Loaded`; `store.loaded()` and `store.person()` read
everything in it back, lazily. The layout (docs/SCHEMA.md, "The durable store"):

    <path>/manifest.json                 store layout, schema version, subject, created; per facet its
                                         added time and the warnings its load gave, redacted
    <path>/facets/<facet_id>/facets.parquet, calls.parquet, blocks.parquet, links.parquet
    <path>/.lock                         the OS file lock every add and migration holds

The rules, each from the charter or the item:

- **Facets accumulate (charter 5).** A facet's directory is written once, never rewritten or removed.
  Adding a facet the store already holds, by `facet_id`, is a no-op when its REF state (`REF_STATE`)
  is the held one's, and refused (`StoreConflictError`) when not: `resolve_ref` keeps the id, so the
  id alone does not say which rows the store holds. Adding one never touches another
  facet's files: the links live per facet (the links ending at the facet and at its ancestors, from
  the `Loaded` it came with), so no file is shared between facets. The manifest is the store's index,
  rewritten whole (temporary name, then `os.replace`) after each facet.
- **One person.** The store holds one subject, named when it is created. A facet with no subject or
  another one is refused with `SubjectError`, naming facet ids, never the subject.
- **One schema version.** The store records the `SCHEMA_VERSION` it was created under. A `Loaded` of
  another version is refused, and so is reading a store of another version, with `StoreVersionError`:
  frames of two versions are never mixed, and nothing is migrated on read. `migrate_store(path)`
  rewrites an older store to this version, by the steps in `just_vcf.store_migration`: every facet's
  files are written anew beside the old ones, checked against them, and swapped in under the lock,
  the manifest committing the swap; the old version's files go once it is committed. A store of a
  later version is refused (no downgrade), and so is one of a version no store was written under.
  This is the one rewrite of a facet's files: the facets and their rows stay what they were.
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

import json
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
from just_vcf.store_migration import (
    STORE_FIRST_VERSION,
    FilledUnknown,
    MigrationStep,
    StoreMigration,
    is_newer,
    migration_chain,
)

logger = logging.getLogger(__name__)

STORE_LAYOUT = 1
"""The store's own layout, recorded in its manifest. A store of a later layout is refused."""

FRAMES = ("facets", "calls", "blocks", "links")
"""The parquet files each facet's directory holds, one per frame."""

_MANIFEST = "manifest.json"
_LOCK = ".lock"
_FACETS_DIR = "facets"
_STAGING = ".staging-"
_MIGRATING = ".migrating-"
_PREVIOUS = ".previous-"
_READY = "READY"
_SCHEMAS = dict(zip(FRAMES, (FACETS, CALLS, BLOCKS, LINKS)))


class StoreError(ValueError):
    """The path is not a store and cannot become one, or the store cannot answer what was asked."""


class StoreVersionError(StoreError):
    """A schema version or store layout other than the store's: refused, never mixed."""


class StoreConflictError(StoreError):
    """A facet the store holds, given again with REF in another state: not read where the held one was,
    or read from another collection or by another `RESOLVE_VERSION`. Its files are never rewritten."""


REF_STATE = ("ref_seqcol_digest", "resolve_version")
"""The facet columns that say how a facet's REF was read. `resolve_ref` keeps the `facet_id`, so a facet
id alone does not tell a held facet's rows from those of the same facet resolved otherwise."""


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
    migrations: tuple[StoreMigration, ...] = ()
    """Every migration the store went through, oldest first, with what each filled unknown."""


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
        but its frames cannot be read or added to until `migrate_store` brings it to this version.
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
        migrating = sorted(child.name for child in where.glob(f"{_MIGRATING}*"))
        if migrating:
            logger.warning(
                "store: a migration in progress, or one that did not finish (%s): migrate_store finishes or "
                "removes it, and the store reads as the version its manifest names until then",
                ", ".join(migrating),
            )
        previous = sorted(child.name for child in where.glob(f"{_PREVIOUS}*"))
        if previous:
            logger.warning(
                "store: the files of a version a migration replaced are kept (%s); nothing reads them, and "
                "removing them is the caller's choice",
                ", ".join(previous),
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
        twice in `loaded`; a facet the store holds whose REF state (`REF_STATE`: the collection REF was
        read from and the `RESOLVE_VERSION` that read it, both null before `resolve_ref`) differs from
        the held one's (`StoreConflictError`). A held facet is never rewritten, so the rows of the same
        facet resolved otherwise go into another store.
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
            self._check_ref_state(facets.filter(pl.col("facet_id").is_in(sorted(held & set(ids)))))
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
        own, in facet order, each line once, then one line for each set of values a migration of the
        store filled unknown (`StoreMigration.filled_unknown`), so a null a migration left is not read
        as one a load wrote. An empty store gives empty frames."""
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
        warnings += [_filled_warning(migration, filled) for migration in manifest.migrations for filled in migration.filled_unknown]
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

    def _check_ref_state(self, given: pl.DataFrame) -> None:
        """Refuse `given`, facet rows the store holds, where one's REF state is not the held one's."""
        if given.is_empty():
            return
        files = [self.path / _FACETS_DIR / facet_id / "facets.parquet" for facet_id in sorted(given["facet_id"].to_list())]
        held = pl.read_parquet(files, columns=["facet_id", *REF_STATE])
        differ = []
        for row in given.select("facet_id", *REF_STATE).sort("facet_id").iter_rows(named=True):
            stored = held.filter(pl.col("facet_id") == row["facet_id"]).row(0, named=True)
            if any(row[name] != stored[name] for name in REF_STATE):
                differ.append(f"{row['facet_id']} (held: {_ref_state(stored)}; given: {_ref_state(row)})")
        if differ:
            raise StoreConflictError(
                f"facet(s) held with REF in another state: {'; '.join(differ)}. A held facet is never rewritten "
                "(resolve_ref keeps its facet_id); add these to a new store, or keep the held rows"
            )

    def _write_facet(self, loaded: Loaded, facet_id: str, links: pl.DataFrame) -> None:
        """Sink one facet's four frames into a staging directory, sync them, and rename it into place."""
        stage = self.path / f"{_STAGING}{os.getpid()}-{os.urandom(4).hex()}"
        stage.mkdir()
        try:
            mine = pl.col("facet_id") == facet_id
            _sink_facet(
                stage,
                {
                    "facets": loaded.facets.filter(mine),
                    "calls": loaded.calls.filter(mine),
                    "blocks": loaded.blocks.filter(mine),
                    "links": links.lazy(),
                },
            )
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


def _ref_state(row: dict) -> str:
    if row["ref_seqcol_digest"] is None:
        return "REF not read"
    version = "an unrecorded RESOLVE_VERSION" if row["resolve_version"] is None else f"RESOLVE_VERSION {row['resolve_version']}"
    return f"REF read from {row['ref_seqcol_digest']} by {version}"


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


def _sink_facet(directory: Path, frames: dict[str, pl.LazyFrame]) -> None:
    """Sink one facet's four frames into `directory` and sync them: an add and a migration write a
    facet's files the same way."""
    pl.collect_all(
        [
            frames["facets"].sink_parquet(directory / "facets.parquet", lazy=True),
            frames["calls"].sink_parquet(directory / "calls.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
            frames["blocks"].sink_parquet(directory / "blocks.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
            frames["links"].sink_parquet(directory / "links.parquet", lazy=True),
        ]
    )
    for name in FRAMES:
        _fsync(directory / f"{name}.parquet")
    _fsync_dir(directory)


def _check_version(store_version: str, other: str, what: str) -> None:
    if store_version == other:
        return
    if other == SCHEMA_VERSION:
        raise StoreVersionError(_version_refusal(store_version, what))
    raise StoreVersionError(f"the store holds schema version {store_version} and {what} is {other}; versions are never mixed")


def _version_refusal(store_version: str, what: str) -> str:
    """Why a store of `store_version` is neither read nor added to by this just-vcf, and what to do."""
    head = (
        f"the store holds schema version {store_version} and {what} is {SCHEMA_VERSION}; versions are never "
        "mixed, and nothing is migrated on read"
    )
    if is_newer(store_version):
        return f"{head}. A later just-vcf wrote the store, and this one does not downgrade it: upgrade just-vcf"
    if migration_chain(store_version) is None:
        return (
            f"{head}. There is no migration from {store_version}: no store was written under it (the store "
            f"shipped under {STORE_FIRST_VERSION}). Load the sources again into a new store"
        )
    return f"{head}. Migrate the store to {SCHEMA_VERSION} first, with just_vcf.migrate_store(path)"


def _filled_warning(migration: StoreMigration, filled: FilledUnknown) -> str:
    return (
        f"store migrated from schema {migration.from_version} to {migration.to_version}: at the step to "
        f"{filled.version}, {filled.frame} {', '.join(filled.columns)} on {len(filled.facet_ids)} facet(s) "
        f"({', '.join(filled.facet_ids)}): {filled.reason}"
    )


def migrate_store(path: Path | str, keep_previous: bool = False) -> StoreMigration:
    """Rewrite the store at `path`, written under an older schema version, to this one.

    The steps of `just_vcf.store_migration` are applied in order from the store's version to
    `SCHEMA_VERSION`, to every facet's four frames. Under the store's lock, each facet is written anew
    into a `.migrating-*` directory beside the old files and checked against them (this version's
    shape, the same rows, the old columns' values unchanged); only when every facet passes is the
    migration marked ready and swapped in: the old `facets/` and manifest move to `.previous-<version>`,
    the new ones move in, and the new manifest, written last, commits it. The old files are then
    removed, or kept with `keep_previous=True`. A migration cut off before it was ready leaves the store
    as it was, and the next call removes what it wrote; one cut off after finishes on the next call. A
    reader meanwhile is refused, by the version its manifest names, as before.

    Returns what was done, also recorded in the manifest's `migrations`: the steps, the facets, and
    every value filled unknown (`FilledUnknown`), which `Store.loaded()` repeats as warnings. A call
    that finishes a migration an earlier one committed but did not clean up after returns that
    migration, as the manifest records it. A store already at this version, with nothing left to
    finish, is left as it is (nothing applied, `migrated` null). Refused with `StoreVersionError`: a
    store of a later version (no downgrade), or of one no store was written under.
    """
    where = Path(path)
    _read_manifest(where)
    with _locked(where):
        finished = _finish_migration(where)
        manifest = _read_manifest(where)
        if manifest.schema_version == SCHEMA_VERSION:
            return finished or StoreMigration(from_version=SCHEMA_VERSION, to_version=SCHEMA_VERSION)
        chain = migration_chain(manifest.schema_version)
        if chain is None:
            raise StoreVersionError(_version_refusal(manifest.schema_version, "this just-vcf"))
        stage = where / f"{_MIGRATING}{os.getpid()}-{os.urandom(4).hex()}"
        (stage / _FACETS_DIR).mkdir(parents=True)
        try:
            filled = _migrate_facets(where, stage, manifest, chain)
            record = StoreMigration(
                from_version=manifest.schema_version,
                to_version=SCHEMA_VERSION,
                migrated=_now(),
                steps=tuple(step.applied() for step in chain),
                facets=tuple(entry.facet_id for entry in manifest.facets),
                filled_unknown=filled,
            )
            _fsync_dir(stage / _FACETS_DIR)
            _write_manifest(
                stage,
                manifest.model_copy(update={"schema_version": SCHEMA_VERSION, "migrations": (*manifest.migrations, record)}),
            )
            marker = stage / _READY
            marker.write_text(json.dumps({"from_version": manifest.schema_version, "keep_previous": keep_previous}) + "\n")
            _fsync(marker)
            _fsync_dir(stage)
        except BaseException:
            shutil.rmtree(stage, ignore_errors=True)
            raise
        _finish_migration(where)
    logger.info("store: migrated %d facet(s) from schema %s to %s", len(record.facets), record.from_version, record.to_version)
    return record


def _migrate_facets(where: Path, stage: Path, manifest: StoreManifest, chain: tuple[MigrationStep, ...]) -> tuple[FilledUnknown, ...]:
    """Write every facet of the store, migrated by `chain`, under `stage`, each checked against its old
    files. Returns what the steps filled unknown, step by step, on the facets each applies to."""
    applies: dict[tuple[int, int], list[str]] = {}
    for entry in manifest.facets:
        old = {name: pl.scan_parquet(where / _FACETS_DIR / entry.facet_id / f"{name}.parquet") for name in FRAMES}
        frames = dict(old)
        for at, step in enumerate(chain):
            # A caveat reads the facet row as the step finds it.
            for which, caveat in enumerate(step.caveats):
                if frames["facets"].filter(caveat.applies).select(pl.len()).collect().item():
                    applies.setdefault((at, which), []).append(entry.facet_id)
            frames = {name: step.apply(name, frame) for name, frame in frames.items()}
        target = stage / _FACETS_DIR / entry.facet_id
        target.mkdir()
        _sink_facet(target, frames)
        _verify_migrated(entry.facet_id, old, target)
    return tuple(
        FilledUnknown(
            version=step.to_version,
            frame=caveat.frame,
            columns=caveat.columns,
            facet_ids=tuple(sorted(applies[(at, which)])),
            reason=caveat.reason,
        )
        for at, step in enumerate(chain)
        for which, caveat in enumerate(step.caveats)
        if (at, which) in applies
    )


def _verify_migrated(facet_id: str, old: dict[str, pl.LazyFrame], target: Path) -> None:
    """Each migrated file has this version's shape, the old file's row count, and the old columns'
    values unchanged (a digest of every row, the new columns cast back to the old dtypes)."""
    for name in FRAMES:
        new = pl.scan_parquet(target / f"{name}.parquet")
        if new.collect_schema() != _SCHEMAS[name]:
            raise StoreError(f"migrating facet {facet_id}: its {name} frame does not have schema {SCHEMA_VERSION}'s shape")
        before = old[name].collect_schema()
        kept = before.names()
        was, now = pl.collect_all([_digest(old[name], kept), _digest(new.select(pl.col(c).cast(before[c]) for c in kept), kept)])
        if was["rows"].item() != now["rows"].item():
            raise StoreError(
                f"migrating facet {facet_id}: its {name} frame has {now['rows'].item()} rows after, {was['rows'].item()} before"
            )
        if was["digest"].item() != now["digest"].item():
            raise StoreError(f"migrating facet {facet_id}: values of its {name} frame changed")


def _digest(frame: pl.LazyFrame, columns: list[str]) -> pl.LazyFrame:
    """Row count and an order-free digest of `columns` (each row's hash, its top 32 bits summed, so the
    sum cannot overflow below four billion rows)."""
    return frame.select(pl.len().alias("rows"), (pl.struct(columns).hash(seed=0) // (1 << 32)).sum().alias("digest"))


def _finish_migration(where: Path) -> StoreMigration | None:
    """Under the store's lock: finish a migration whose files were all written and checked (its `READY`
    marker), and remove one that was cut off before. Returns the migration this call committed, or the
    one an earlier call committed (the live manifest's last) when only its cleanup was left; else None.

    Each move is skipped when it is already done, so a call cut off anywhere finishes on the next. The
    order: the old manifest copied to `.previous-<version>`, the old `facets/` moved there, the new
    `facets/` moved in, the new manifest moved over the old (the commit), then the old files removed
    unless kept, and the migration's directory last, so its marker outlives every step.
    """
    finished: StoreMigration | None = None
    for stage in sorted(where.glob(f"{_MIGRATING}*")):
        marker = stage / _READY
        if not marker.is_file():
            logger.info("store: removing %s, a migration cut off before it was ready", stage.name)
            shutil.rmtree(stage, ignore_errors=True)
            continue
        ready = json.loads(marker.read_text())
        previous = where / f"{_PREVIOUS}{ready['from_version']}"
        staged_manifest = stage / _MANIFEST
        if staged_manifest.is_file():
            previous.mkdir(exist_ok=True)
            if not (previous / _MANIFEST).is_file():
                shutil.copy2(where / _MANIFEST, previous / _MANIFEST)
                _fsync(previous / _MANIFEST)
            staged_facets = stage / _FACETS_DIR
            if staged_facets.is_dir():
                live = where / _FACETS_DIR
                if live.is_dir():
                    os.rename(live, previous / _FACETS_DIR)
                    _fsync_dir(previous)
                os.rename(staged_facets, live)
                _fsync_dir(where)
            finished = _read_manifest(stage).migrations[-1]
            os.replace(staged_manifest, where / _MANIFEST)
            _fsync_dir(where)
        else:
            # Committed by an earlier call cut off before its cleanup: the live manifest records it.
            live = _read_manifest(where)
            if live.migrations and live.migrations[-1].from_version == ready["from_version"]:
                finished = live.migrations[-1]
        if not ready["keep_previous"]:
            shutil.rmtree(previous, ignore_errors=True)
        shutil.rmtree(stage, ignore_errors=True)
    return finished


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
