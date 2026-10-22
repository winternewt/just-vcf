"""The loader's cache: what each entry is, sweeping stale versions, and `clear_cache` (RM268).

The cache (`cache_dir()`) holds four directories of entries:

- `vcf/<key>/`: a rewritten copy of an input polars-bio cannot read as it is (`readable_path`), and a
  `structural/` copy inside it for a file polars-bio refused. Kind `rewrite_copy`.
- `keyvalue/<key>/` and `cgi/<key>/`: the decompressed payload a key=value or Complete Genomics
  loader reads while it spills the rows, removed when the load ends (`work_copy`). Indexed so a copy
  a dead loader left is still traced to its source.
- `frames/<name>/`: the parquet spills the lazy frames read. `<facet_id>-<SCHEMA_VERSION>` is a load's
  spill (`frames_spill`), or a derived facet's under its own id (`rejoin_spill`, `normalize_spill`,
  `lift_spill`);
  `<facet_id>-<SCHEMA_VERSION>-ref-<digest>-<RESOLVE_VERSION>` is `resolve_ref`'s (`resolve_spill`) and
  `…-vrs-<digest>-<MINT_VERSION>` is `mint_vrs_ids`'s (`vrs_spill`). Every spill entry, `load`'s
  included, holds one directory per generation, `<key>/`, named by the content it holds or was
  computed from (`Spill`).

Each entry has two files beside it: `<name>.json`, its manifest (`CacheEntry`: kind, version keys,
the root source's sha256, the facet it holds and the facet it came from, size, created and last
used), and `<name>.lock`, an OS file lock its writer holds for as long as it writes. Nothing else
under the cache dir is ever touched here: the reference store under `reference/` is not a cache of
a genome, and only a directory whose name has an entry's shape is an entry.

An entry without a manifest was written before the index. Its origin is unknown, so `clear_cache`
never removes it; only the stale-version sweep does, and only where its name shows the version: a
`vcf/` entry without a manifest predates `REWRITE_VERSION` 7, the version that brought the index,
so it is stale by construction, and a `keyvalue/` or `cgi/` one predates their index.

Removing an entry a live `Loaded` reads breaks that frame's next `collect`. The stale-version sweep
cannot hit one: no current frame reads an old version, and it removes only entries older than the
running version, never a newer install's that shares the cache. `clear_cache` can, and that is the caller's
call. `load` or a derived step run again never can: it reuses the generation an earlier run
committed for the same content, or commits another beside it, and never rewrites or removes a
complete one.
"""

import hashlib
import json
import logging
import os
import re
import shutil
import time
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from pathlib import Path

from platformdirs import user_cache_dir
from pydantic import BaseModel, ConfigDict, ValidationError

if os.name == "nt":
    import msvcrt
else:
    import fcntl

logger = logging.getLogger(__name__)

INDEX_VERSION = 1
"""The manifest's own format."""

REWRITE_COPY = "rewrite_copy"
FRAMES_SPILL = "frames_spill"
REJOIN_SPILL = "rejoin_spill"
RESOLVE_SPILL = "resolve_spill"
VRS_SPILL = "vrs_spill"
NORMALIZE_SPILL = "normalize_spill"
LIFT_SPILL = "lift_spill"
WORK_COPY = "work_copy"
CACHE_KINDS: frozenset[str] = frozenset(
    {REWRITE_COPY, FRAMES_SPILL, REJOIN_SPILL, RESOLVE_SPILL, VRS_SPILL, NORMALIZE_SPILL, LIFT_SPILL, WORK_COPY}
)
"""What a cache entry can be. Grows as derived-facet steps that spill land."""

_COPY_DIRS = ("vcf", "keyvalue", "cgi")
"""Directories of copies of a source, each entry named by a 32-hex key that hashes the source."""
_DIRS = (*_COPY_DIRS, "frames")
_VCF_NAME = re.compile(r"[0-9a-f]{32}")
_FRAMES_NAME = re.compile(
    r"(?P<facet>[0-9a-f]{32})-(?P<schema>[0-9][0-9.]*)(?:-(?P<op>ref|vrs)-(?P<digest>[A-Za-z0-9_-]+)-(?P<step>[0-9]+)|-(?P<lift>lift))?"
)
"""A `frames/` entry's name. `-lift` is only ever a legacy name: `liftover` wrote `<facet_id>-<SCHEMA_VERSION>-lift`
before the index, and now spills under `<facet_id>-<SCHEMA_VERSION>` like the other derived facets."""
_LEGACY_STEP_KEY = {"ref": ("RESOLVE_VERSION", RESOLVE_SPILL), "vrs": ("MINT_VERSION", VRS_SPILL)}


def cache_dir() -> Path:
    """Where rewritten inputs and spills are kept: `JUST_VCF_CACHE_DIR`, else the platform's user cache."""
    configured = os.environ.get("JUST_VCF_CACHE_DIR")
    return Path(configured) if configured else Path(user_cache_dir("just-vcf"))


class CacheEntry(BaseModel):
    """The manifest beside one cache entry. Times are ISO-8601 UTC."""

    model_config = ConfigDict(frozen=True)

    index_version: int = INDEX_VERSION
    kind: str
    versions: dict[str, str]
    """The version keys the entry was written under, e.g. `{"SCHEMA_VERSION": "0.0.7"}`."""
    source_sha256: str | None
    """The root source file's sha256: the input itself for a copy, the input a derived facet goes back to."""
    source_path: str | None
    """The root source's resolved path when written, so a source deleted since can still be named."""
    member: str | None = None
    facet_id: str | None = None
    """The facet whose rows the spill holds."""
    parent_facet_id: str | None = None
    """The facet a derived spill was made from."""
    size_bytes: int | None = None
    """Bytes on disk when the last write finished; null while it is written, or when it never finished."""
    created: datetime
    last_used: datetime
    """When the entry was last written or reused: for a spill, the last run of its step, whether it
    wrote a generation or reused one; for a rewrite copy, also a later load that reused it."""
    # `kind` is checked against CACHE_KINDS where it is written (`writing`), not where it is read: a
    # manifest a newer version wrote, with a kind this one does not know, is read as data and owned by
    # no sweep here.


class Origin(BaseModel):
    """Where an entry's rows come from: what a writer knows when it writes."""

    model_config = ConfigDict(frozen=True)

    source_sha256: str | None
    source_path: str | None
    member: str | None = None
    facet_id: str | None = None
    parent_facet_id: str | None = None


def derived_origin(facet: dict, facet_id: str) -> Origin:
    """The origin of a spill a derived step writes from `facet` (a `facets` row): its root source,
    carried down every derived facet, and `facet` itself as the parent."""
    return Origin(
        source_sha256=facet["source_sha256"],
        source_path=facet["source_path"],
        member=facet["member"],
        facet_id=facet_id,
        parent_facet_id=facet["facet_id"],
    )


class CacheEntryInfo(BaseModel):
    """One entry as `clear_cache` reports it. It carries no source path: a report may be shown."""

    model_config = ConfigDict(frozen=True)

    path: str
    """Relative to the cache dir, e.g. `frames/<name>`."""
    kind: str | None
    """Null for an entry written before the index (unknown origin)."""
    versions: dict[str, str]
    source_sha256: str | None
    facet_id: str | None
    parent_facet_id: str | None
    bytes: int
    """Measured on disk now."""
    created: datetime | None
    last_used: datetime | None


class CacheClearReport(BaseModel):
    """What `clear_cache` removed, or would remove with `dry_run`, and what it left and why."""

    model_config = ConfigDict(frozen=True)

    dry_run: bool
    removed: tuple[CacheEntryInfo, ...]
    """Removed; with `dry_run`, what would be."""
    being_written: tuple[CacheEntryInfo, ...]
    """Matched, but a writer holds the entry's lock: left alone."""
    unknown_origin: tuple[CacheEntryInfo, ...]
    """Written before the index, so nothing says what they came from: never removed here. The
    stale-version sweep in `load` and the derived steps removes those whose name shows an old version."""

    @property
    def bytes_removed(self) -> int:
        return sum(entry.bytes for entry in self.removed)


class _Entry:
    """An entry on disk: its directory, manifest and lock paths, and the manifest if there is one."""

    def __init__(self, root: Path, directory: str, name: str) -> None:
        self.relative = f"{directory}/{name}"
        self.path = root / directory / name
        self.manifest_path = root / directory / f"{name}.json"
        self.lock_path = root / directory / f"{name}.lock"
        self.manifest = _read_manifest(self.manifest_path)
        self.unreadable = self.manifest is None and self.manifest_path.is_file()
        """A manifest this version cannot read (a newer one wrote it): kept, never judged stale here."""

    def size(self) -> int:
        return _size(self.path)

    def info(self) -> CacheEntryInfo:
        m = self.manifest
        return CacheEntryInfo(
            path=self.relative,
            kind=m.kind if m else None,
            versions=m.versions if m else _legacy_versions(self.relative),
            source_sha256=m.source_sha256 if m else None,
            facet_id=m.facet_id if m else None,
            parent_facet_id=m.parent_facet_id if m else None,
            bytes=self.size(),
            created=m.created if m else None,
            last_used=m.last_used if m else None,
        )


@contextmanager
def writing(entry: Path, kind: str, versions: dict[str, object], origin: Origin) -> Iterator[Path]:
    """Hold `entry`'s lock while the caller writes into it, and record its manifest.

    Waits for another writer of the same entry to finish. The manifest is written when the lock is
    taken (so an entry whose writer died can still be traced to its source, with no size) and again
    when the block finishes, with the size and the time. `created` is kept from an earlier manifest of
    the same kind and versions.
    """
    if kind not in CACHE_KINDS:
        raise ValueError(f"unknown cache entry kind {kind!r}; known: {sorted(CACHE_KINDS)}")
    entry.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = entry.parent / f"{entry.name}.json"
    handle = _acquire(entry.parent / f"{entry.name}.lock", blocking=True)
    try:
        stamped = {key: str(value) for key, value in versions.items()}
        previous = _read_manifest(manifest_path)
        now = _now()
        kept = previous is not None and previous.kind == kind and previous.versions == stamped and entry.exists()
        manifest = CacheEntry(
            kind=kind,
            versions=stamped,
            created=previous.created if kept else now,
            last_used=now,
            **origin.model_dump(),
        )
        _write_manifest(manifest_path, manifest)
        yield entry
        _write_manifest(manifest_path, manifest.model_copy(update={"size_bytes": _size(entry), "last_used": _now()}))
    finally:
        _release(handle)


class Generation(BaseModel):
    """`result.json` in one generation of a derived spill: the files its frames read, and the step's
    own warnings, so a reused generation gives back what the run that wrote it gave."""

    model_config = ConfigDict(frozen=True)

    files: dict[str, tuple[str, ...]]
    """Named groups of file names in the generation, each in the order a frame reads them."""
    warnings: tuple[str, ...] = ()


class Committed:
    """A complete generation on disk: where it is and what `result.json` says it holds."""

    def __init__(self, path: Path, record: Generation) -> None:
        self.path = path
        self.record = record

    def files(self, group: str) -> list[Path]:
        return [self.path / name for name in self.record.files[group]]

    @property
    def warnings(self) -> list[str]:
        return list(self.record.warnings)


class Spill:
    """One run of `load` or a derived step writing into its entry, under the entry's lock (`spilling`).

    The run writes into `stage`, a directory of its own inside the entry, and `commit`s it by renaming
    it to `<entry>/<key>`, the key being the content the result was computed from. A generation is
    never written into or removed once committed, so a frame an earlier run returned keeps reading its
    files whatever later runs do; a later run with the same content `find`s it and reuses it instead.
    """

    def __init__(self, entry: Path) -> None:
        self.entry = entry
        entry.mkdir(parents=True, exist_ok=True)
        for leftover in entry.glob(f"{_STAGING}*"):
            # The entry's lock is held, so a staging directory here is a writer's that died: no frame reads it.
            logger.info("cache: removing %s/%s, left by a writer that did not finish", entry.name, leftover.name)
            shutil.rmtree(leftover, ignore_errors=True)
        self.stage = entry / f"{_STAGING}{os.getpid()}-{os.urandom(4).hex()}"
        self.stage.mkdir()

    def find(self, key: str) -> Committed | None:
        """The complete generation under `key`, or None. One that cannot be read whole is removed: a
        frame over it fails already, and the run then writes it again."""
        path = self.entry / key
        if not path.is_dir():
            return None
        record = _read_generation(path / _RESULT)
        if record is not None and all((path / name).is_file() for names in record.files.values() for name in names):
            return Committed(path, record)
        logger.warning("cache: %s/%s is incomplete; written again", self.entry.name, key)
        shutil.rmtree(path)
        return None

    def commit(self, key: str, files: dict[str, list[Path]], warnings: list[str]) -> Committed:
        """Record what `stage` holds and move it into place as the generation `key`."""
        record = Generation(
            files={group: tuple(path.relative_to(self.stage).as_posix() for path in paths) for group, paths in files.items()},
            warnings=tuple(warnings),
        )
        (self.stage / _RESULT).write_text(record.model_dump_json(indent=1) + "\n")
        target = self.entry / key
        os.rename(self.stage, target)
        return Committed(target, record)

    def discard(self) -> None:
        if self.stage.exists():
            shutil.rmtree(self.stage)


SPILL_LAYOUT = 2
"""How a spill is laid out, a key in its manifest's versions: 2 is one directory per generation
(`Spill`); 1 wrote the files into the entry itself, so the sweep retires those entries. `load`'s
spill and `liftover`'s took generations after the other derived steps, at the same layout number:
their manifests before then carried no `SPILL_LAYOUT` key at all."""

_STAGING = ".staging-"
_RESULT = "result.json"


@contextmanager
def spilling(entry: Path, kind: str, versions: dict[str, object], origin: Origin) -> Iterator[Spill]:
    """`writing` for a spill: hold the entry's lock, and give the run a `Spill` to stage into.

    The manifest is updated whether the run commits a new generation or reuses one (`last_used`, size).
    The staging directory goes when the block ends, committed or not.
    """
    with writing(entry, kind, versions, origin):
        spill = Spill(entry)
        try:
            yield spill
        finally:
            spill.discard()


def content_key(files: list[Path], facet: dict) -> str:
    """What a derived result is computed from: the bytes of the staged files and the facet row.

    Under one polars release the same rows sink to the same bytes (measured 2026-10-22: 300,000 rows,
    two processes, two caches), so a run on the same input finds the generation an earlier one wrote.
    After a polars upgrade a run may write a second one beside it, which costs disk, never correctness.
    A run on other rows under the same facet id (a filtered `Loaded`) gets a key of its own: a reuse
    never crosses inputs.
    """
    digest = hashlib.sha256()
    for path in files:
        _digest_file(digest, path)
    digest.update(json.dumps(facet, sort_keys=True, default=str).encode())
    return digest.hexdigest()[:32]


def _read_generation(path: Path) -> Generation | None:
    if not path.is_file():
        return None
    try:
        return Generation.model_validate_json(path.read_bytes())
    except ValidationError:
        return None


def sweep_stale(kind: str, current: dict[str, object]) -> CacheClearReport:
    """Remove every entry of `kind` written under an older version than `current`, logging the bytes freed.

    Older means no version key ahead of `current` and at least one behind, compared numerically
    (`_older`). An entry a newer version wrote is kept: two installs of different versions can share
    the platform user cache, and the older one must not remove the newer one's live spills.

    Run by the step that owns the kind, before it writes: `load` for `rewrite_copy`, `work_copy` and `frames_spill`,
    each derived step for its own. An entry with no manifest is judged by its name where the name shows
    its version, and a `vcf/` one is stale by construction (module docstring). Always safe: no current
    lazy frame reads an entry of an old version. An entry being written is left alone.
    """
    if kind not in CACHE_KINDS:
        raise ValueError(f"unknown cache entry kind {kind!r}")
    stamped = {key: str(value) for key, value in current.items()}
    chosen = [entry for entry in _entries(cache_dir()) if _stale(entry, kind, stamped)]
    report = _remove(chosen, dry_run=False, unknown=[])
    if report.removed:
        by_kind = Counter(info.kind or "unindexed" for info in report.removed)
        logger.info(
            "cache: removed %d stale entr%s of %s (%s), %d bytes freed",
            len(report.removed),
            "y" if len(report.removed) == 1 else "ies",
            kind,
            ", ".join(f"{count} {name}" for name, count in sorted(by_kind.items())),
            report.bytes_removed,
        )
    if report.being_written:
        logger.info("cache: %d stale %s entr(ies) left alone, being written", len(report.being_written), kind)
    return report


def clear_cache(
    source: Path | str | None = None,
    older_than: timedelta | None = None,
    dry_run: bool = False,
) -> CacheClearReport:
    """Remove cache entries: those derived from `source`, those not used for `older_than`, or both.

    `source` names a source file: an existing file is hashed, and matches every entry whose root source
    has that sha256 or was read from that path; a 64-character hex string that is not a file is taken
    as the sha256 itself; any other path (a file deleted since) matches the entries read from that
    path. Derived spills (`rejoin_split`, `resolve_ref`, `normalize`, `mint_vrs_ids`) carry their root
    source, so removing a genome's copies removes them too. `older_than` matches entries whose
    `last_used` is longer ago than that. Given both, an entry must match both; given neither, every
    indexed entry matches.

    `dry_run` lists what would go and removes nothing. An entry whose writer holds its lock is left
    alone and reported in `being_written`. An entry written before the index is never removed here
    (`unknown_origin`). Only `vcf/`, `keyvalue/`, `cgi/` and `frames/` entries are considered: the
    reference store is never touched.

    Removing an entry that a live `Loaded` reads breaks that frame's next `collect`: that is the
    caller's call. Reload afterwards.
    """
    if older_than is not None and older_than < timedelta(0):
        raise ValueError("older_than must not be negative")
    matches_source = _source_matcher(source) if source is not None else None
    cutoff = _now() - older_than if older_than is not None else None
    chosen: list[_Entry] = []
    unknown: list[_Entry] = []
    for entry in _entries(cache_dir()):
        manifest = entry.manifest
        if manifest is None:
            unknown.append(entry)
            continue
        if matches_source is not None and not matches_source(manifest):
            continue
        if cutoff is not None and manifest.last_used >= cutoff:
            continue
        chosen.append(entry)
    report = _remove(chosen, dry_run=dry_run, unknown=unknown)
    logger.info(
        "cache: %s %d entr(ies), %d bytes; %d being written; %d of unknown origin kept",
        "would remove" if dry_run else "removed",
        len(report.removed),
        report.bytes_removed,
        len(report.being_written),
        len(report.unknown_origin),
    )
    return report


def _source_matcher(source: Path | str):
    path = Path(source)
    if path.is_file():
        digest = file_sha256(path)
        resolved = str(path.resolve())
        return lambda m: m.source_sha256 == digest or m.source_path == resolved
    if isinstance(source, str) and re.fullmatch(r"[0-9a-f]{64}", source):
        return lambda m: m.source_sha256 == source
    resolved = str(path.resolve())
    return lambda m: m.source_path == resolved


def _version(value: str) -> tuple[int, ...] | None:
    """A version key's value as numbers: `7` is `(7,)`, `0.0.11` is `(0, 0, 11)`; None when it is not one."""
    parts = value.split(".")
    return tuple(int(part) for part in parts) if all(part.isdigit() for part in parts) else None


def _older(versions: dict[str, str], current: dict[str, str]) -> bool:
    """Whether an entry's version keys put it behind `current`: none ahead, and at least one behind.

    Compared numerically, so `0.0.9` is behind `0.0.11`. A key `current` has and the entry lacks is
    behind: the entry was written before the key existed (`load`'s spills before `SPILL_LAYOUT`). A
    key the entry has and `current` lacks, or a value that is not a number, cannot be ordered, so the
    entry may be a newer version's and is not behind. Two installs of different versions share the
    platform user cache, and an older one must never remove a newer one's live spills.
    """
    if set(versions) - set(current):
        return False
    behind = False
    for key, now in current.items():
        if key not in versions:
            behind = True
            continue
        mine, theirs = _version(versions[key]), _version(now)
        if mine is None or theirs is None:
            return False
        width = max(len(mine), len(theirs))
        mine, theirs = mine + (0,) * (width - len(mine)), theirs + (0,) * (width - len(theirs))
        if mine > theirs:
            return False
        behind = behind or mine < theirs
    return behind


def _stale(entry: _Entry, kind: str, current: dict[str, str]) -> bool:
    if entry.manifest is not None:
        return entry.manifest.kind == kind and _older(entry.manifest.versions, current)
    if entry.unreadable or not entry.path.is_dir():
        return False
    if entry.relative.startswith("vcf/"):
        return kind == REWRITE_COPY
    if not entry.relative.startswith("frames/"):
        # `keyvalue/` and `cgi/`: every working copy written since the index has a manifest.
        return kind == WORK_COPY
    match = _FRAMES_NAME.fullmatch(entry.path.name)
    if match["lift"]:
        # Written before the index and before generations: no run of this version reads it.
        return kind == LIFT_SPILL
    if match["op"] is None:
        # Shared by load's, rejoin's and normalize's spills: only the schema version is in the name.
        return kind == FRAMES_SPILL and match["schema"] != current.get("SCHEMA_VERSION")
    key, owner = _LEGACY_STEP_KEY[match["op"]]
    return kind == owner and (match["schema"], match["step"]) != (current.get("SCHEMA_VERSION"), current.get(key))


def _legacy_versions(relative: str) -> dict[str, str]:
    """The versions an unindexed entry's name shows (none for a `vcf/` key, which hashes its version)."""
    match = _FRAMES_NAME.fullmatch(relative.split("/", 1)[1]) if relative.startswith("frames/") else None
    if match is None:
        return {}
    versions = {"SCHEMA_VERSION": match["schema"]}
    if match["op"]:
        versions[_LEGACY_STEP_KEY[match["op"]][0]] = match["step"]
    return versions


def _entries(root: Path) -> list[_Entry]:
    """Every entry under `_DIRS`: a directory or a manifest whose name has an entry's shape."""
    entries = []
    for directory in _DIRS:
        pattern = _FRAMES_NAME if directory == "frames" else _VCF_NAME
        base = root / directory
        if not base.is_dir():
            continue
        names = set()
        for child in base.iterdir():
            name = child.name[: -len(".json")] if child.suffix == ".json" and child.is_file() else child.name
            if (child.is_dir() or child.suffix == ".json") and pattern.fullmatch(name):
                names.add(name)
        entries += [_Entry(root, directory, name) for name in sorted(names)]
    return entries


def _remove(chosen: list[_Entry], *, dry_run: bool, unknown: list[_Entry]) -> CacheClearReport:
    removed: list[CacheEntryInfo] = []
    busy: list[CacheEntryInfo] = []
    for entry in chosen:
        handle = _acquire(entry.lock_path, blocking=False)
        if handle is None:
            busy.append(entry.info())
            continue
        try:
            info = entry.info()
            if not dry_run:
                if entry.path.is_dir():
                    shutil.rmtree(entry.path)
                entry.manifest_path.unlink(missing_ok=True)
                if os.name != "nt":
                    # A writer waiting on this lock finds its file gone and takes a fresh one (`_acquire`).
                    entry.lock_path.unlink(missing_ok=True)
            removed.append(info)
        finally:
            _release(handle)
        if not dry_run and os.name == "nt":
            _unlink_if_free(entry.lock_path)
    return CacheClearReport(
        dry_run=dry_run,
        removed=tuple(removed),
        being_written=tuple(busy),
        unknown_origin=tuple(entry.info() for entry in unknown),
    )


def _acquire(lock_path: Path, *, blocking: bool) -> int | None:
    """An open, locked descriptor on `lock_path`, or None when not blocking and another holds it.

    After taking the lock it checks the file is still the one at the path: a remover unlinks the lock
    file with the entry, and a writer that was waiting on the old file must not write unguarded.
    """
    while True:
        handle = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
        if not _lock(handle, blocking):
            os.close(handle)
            return None
        held = os.fstat(handle)
        current = _stat(lock_path)
        if current is not None and (current.st_dev, current.st_ino) == (held.st_dev, held.st_ino):
            return handle
        _release(handle)


def _stat(path: Path) -> os.stat_result | None:
    try:
        return os.stat(path)
    except FileNotFoundError:
        return None


def _lock(handle: int, blocking: bool) -> bool:
    if os.name == "nt":
        # Untested here: msvcrt locks a byte range, and LK_NBLCK fails at once when it is held.
        while True:
            if _try_msvcrt(handle):
                return True
            if not blocking:
                return False
            time.sleep(0.1)
    try:
        fcntl.flock(handle, fcntl.LOCK_EX if blocking else fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return False
    return True


def _try_msvcrt(handle: int) -> bool:
    try:
        msvcrt.locking(handle, msvcrt.LK_NBLCK, 1)
    except OSError:
        return False
    return True


def _release(handle: int) -> None:
    if os.name == "nt":
        os.lseek(handle, 0, os.SEEK_SET)
        msvcrt.locking(handle, msvcrt.LK_UNLCK, 1)
    else:
        fcntl.flock(handle, fcntl.LOCK_UN)
    os.close(handle)


def _unlink_if_free(lock_path: Path) -> None:
    """Windows cannot unlink an open file: the lock file goes only when no writer has it open."""
    try:
        lock_path.unlink(missing_ok=True)
    except PermissionError:
        logger.debug("cache: %s is open by a writer; left in place", lock_path.name)


def _read_manifest(path: Path) -> CacheEntry | None:
    """The manifest at `path`; None when there is none, or when this version cannot read it."""
    if not path.is_file():
        return None
    try:
        return CacheEntry.model_validate_json(path.read_bytes())
    except ValidationError as err:
        logger.info("cache: %s is not a manifest this version reads (%d error(s)); entry kept", path.name, err.error_count())
        return None


def _write_manifest(path: Path, manifest: CacheEntry) -> None:
    """Written to a temporary name and moved over, so a reader never sees half a manifest."""
    temporary = path.with_name(f"{path.name}.tmp-{os.getpid()}")
    temporary.write_text(manifest.model_dump_json(indent=1) + "\n")
    os.replace(temporary, path)


def _size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    if not path.is_dir():
        return 0
    return sum(child.stat().st_size for child in path.rglob("*") if child.is_file())


def _now() -> datetime:
    return datetime.now(timezone.utc)


def file_sha256(path: Path) -> str:
    """The sha256 of a file's bytes, hex."""
    digest = hashlib.sha256()
    _digest_file(digest, path)
    return digest.hexdigest()


def _digest_file(digest: "hashlib._Hash", path: Path) -> None:
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
