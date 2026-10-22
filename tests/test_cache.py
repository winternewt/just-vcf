"""RM268: the cache index, the stale-version sweep, and `clear_cache`.

Every test runs on the per-test cache `isolated_cache` sets up (conftest), never a real one. Sizes
are measured from the files on disk at runtime.
"""

import hashlib
import json
import logging
import os
import random
import threading
import zipfile
from datetime import timedelta
from pathlib import Path

import pytest

from just_vcf import ReferenceStore, clear_cache, load, normalize, rejoin_split, resolve_ref
from just_vcf import cache
from just_vcf.cache import CACHE_KINDS, Origin
from just_vcf.loaders.vcf import REWRITE_VERSION
from just_vcf.normalize import NORMALIZE_VERSION
from just_vcf.reference import RESOLVE_VERSION
from just_vcf.rejoin import REJOIN_VERSION
from just_vcf.schema import SCHEMA_VERSION

SEQ = "".join(random.Random(268).choice("ACGT") for _ in range(1000))


def vcf_zip(tmp_path: Path, name: str, positions: tuple[int, ...] = (100, 200, 300)) -> Path:
    """A zipped one-sample VCF on contig 1 of `SEQ`: zipped, so `load` has to copy it into the cache.

    The record at 200 is written split, one row per ALT, so `rejoin_split` has a group to join.
    """
    meta = ["##fileformat=VCFv4.2", f"##contig=<ID=1,length={len(SEQ)}>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    rows = []
    for pos in positions:
        ref = SEQ[pos - 1]
        alts = [b for b in "ACGT" if b != ref]
        if pos == 200:
            rows += [f"1\t{pos}\t.\t{ref}\t{alts[0]}\t50\tPASS\t.\tGT\t1/0", f"1\t{pos}\t.\t{ref}\t{alts[1]}\t50\tPASS\t.\tGT\t0/1"]
        else:
            rows.append(f"1\t{pos}\t.\t{ref}\t{alts[0]}\t50\tPASS\t.\tGT\t0/1")
    path = tmp_path / f"{name}.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{name}.vcf", "\n".join([*meta, header, *rows]) + "\n")
    return path


@pytest.fixture
def default_store(tmp_path: Path) -> ReferenceStore:
    """The reference store at its default place, `reference/` inside the cache dir."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">1\n" + "\n".join(SEQ[i : i + 60] for i in range(0, len(SEQ), 60)) + "\n")
    store = ReferenceStore.open()
    store.add_fasta(fasta)
    assert Path(store.path) == cache.cache_dir() / "reference"
    return store


def manifests(kind: str | None = None) -> dict[str, cache.CacheEntry]:
    root = cache.cache_dir()
    found = {}
    for directory in ("vcf", "keyvalue", "cgi", "frames"):
        for path in sorted((root / directory).glob("*.json")) if (root / directory).is_dir() else []:
            manifest = cache.CacheEntry.model_validate_json(path.read_bytes())
            if kind is None or manifest.kind == kind:
                found[f"{directory}/{path.stem}"] = manifest
    return found


def size_on_disk(path: Path) -> int:
    return sum(child.stat().st_size for child in path.rglob("*") if child.is_file())


def tree(path: Path) -> set[str]:
    return {str(child.relative_to(path)) for child in path.rglob("*")}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def full_chain(path: Path, store: ReferenceStore):
    """load → rejoin_split → resolve_ref → normalize: one entry of each kind but `vrs_spill`."""
    return normalize(resolve_ref(rejoin_split(load(path)), store), store)


def test_every_entry_records_what_it_is_and_its_root_source(tmp_path: Path, default_store: ReferenceStore) -> None:
    path = vcf_zip(tmp_path, "one")
    loaded = load(path)
    rejoined = rejoin_split(loaded)
    resolved = resolve_ref(rejoined, default_store)
    normalized = normalize(resolved, default_store)

    found = manifests()
    assert {m.kind for m in found.values()} == CACHE_KINDS - {cache.VRS_SPILL, cache.LIFT_SPILL, cache.WORK_COPY}
    assert {m.source_sha256 for m in found.values()} == {sha256(path)}
    assert {m.source_path for m in found.values()} == {str(path.resolve())}
    load_id = loaded.facets.collect()["facet_id"][0]
    rejoin_id = rejoined.facets.collect()["facet_id"][0]
    by_kind = {m.kind: (name, m) for name, m in found.items()}
    assert by_kind[cache.FRAMES_SPILL][1].facet_id == load_id
    assert (by_kind[cache.REJOIN_SPILL][1].facet_id, by_kind[cache.REJOIN_SPILL][1].parent_facet_id) == (rejoin_id, load_id)
    assert by_kind[cache.RESOLVE_SPILL][1].parent_facet_id == rejoin_id
    assert by_kind[cache.NORMALIZE_SPILL][1].facet_id == normalized.facets.collect()["facet_id"][0]
    assert by_kind[cache.NORMALIZE_SPILL][1].parent_facet_id == rejoin_id
    assert by_kind[cache.REWRITE_COPY][1].versions == {"REWRITE_VERSION": str(REWRITE_VERSION)}
    assert by_kind[cache.REJOIN_SPILL][1].versions == {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": str(REJOIN_VERSION), "SPILL_LAYOUT": str(cache.SPILL_LAYOUT)}
    assert by_kind[cache.RESOLVE_SPILL][1].versions == {"SCHEMA_VERSION": SCHEMA_VERSION, "RESOLVE_VERSION": str(RESOLVE_VERSION), "SPILL_LAYOUT": str(cache.SPILL_LAYOUT)}
    assert by_kind[cache.NORMALIZE_SPILL][1].versions == {"SCHEMA_VERSION": SCHEMA_VERSION, "NORMALIZE_VERSION": str(NORMALIZE_VERSION), "SPILL_LAYOUT": str(cache.SPILL_LAYOUT)}
    for name, manifest in found.items():
        assert manifest.size_bytes == size_on_disk(cache.cache_dir() / name)
        assert manifest.created.utcoffset() == timedelta(0) and manifest.created <= manifest.last_used
        # The lock beside it is free once the step returned.
        handle = cache._acquire(cache.cache_dir() / f"{name}.lock", blocking=False)
        assert handle is not None
        cache._release(handle)


def test_reusing_a_copy_keeps_created_and_moves_last_used(tmp_path: Path) -> None:
    path = vcf_zip(tmp_path, "one")
    load(path)
    ((name, first),) = manifests(cache.REWRITE_COPY).items()
    load(path)
    ((again, second),) = manifests(cache.REWRITE_COPY).items()
    assert again == name
    assert second.created == first.created
    assert second.last_used > first.last_used


def write_entry(directory: str, name: str, kind: str, versions: dict, origin: Origin | None = None) -> Path:
    """An entry written the way the steps write theirs, under versions the test chooses."""
    entry = cache.cache_dir() / directory / name
    with cache.writing(entry, kind, versions, origin or Origin(source_sha256="0" * 64, source_path=None)):
        entry.mkdir(parents=True, exist_ok=True)
        (entry / "calls.parquet").write_bytes(b"x" * 1000)
    return entry


def legacy_entry(directory: str, name: str, size: int = 500) -> Path:
    """An entry as code before the index wrote it: the directory alone, no manifest."""
    entry = cache.cache_dir() / directory / name
    entry.mkdir(parents=True)
    (entry / "input.vcf.gz").write_bytes(b"y" * size)
    return entry


def test_load_sweeps_the_versions_it_owns_and_nothing_else(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    root = cache.cache_dir()
    stale_copy = write_entry("vcf", "a" * 32, cache.REWRITE_COPY, {"REWRITE_VERSION": REWRITE_VERSION - 1})
    unindexed_copy = legacy_entry("vcf", "b" * 32)
    stale_spill = legacy_entry("frames", f"{'c' * 32}-0.0.1")
    stale_rejoin = write_entry("frames", f"{'d' * 32}-{SCHEMA_VERSION}", cache.REJOIN_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": 0})
    stale_resolve = legacy_entry("frames", f"{'e' * 32}-{SCHEMA_VERSION}-ref-SQ_ab-c-{RESOLVE_VERSION - 1}")
    current_unindexed = legacy_entry("frames", f"{'f' * 32}-{SCHEMA_VERSION}")
    not_an_entry = legacy_entry("frames", "notes")
    reference = legacy_entry("reference", "collections")
    expected_freed = sum(size_on_disk(p) for p in (stale_copy, unindexed_copy, stale_spill))

    with caplog.at_level(logging.INFO, logger="just_vcf.cache"):
        load(vcf_zip(tmp_path, "one"))

    assert not stale_copy.exists() and not (root / "vcf" / f"{'a' * 32}.json").exists()
    assert not unindexed_copy.exists() and not stale_spill.exists()
    # Owned by rejoin_split and resolve_ref: load leaves them for those steps.
    assert stale_rejoin.exists() and stale_resolve.exists()
    assert current_unindexed.exists() and not_an_entry.exists() and reference.exists()
    freed = [r for r in caplog.records if "stale" in r.getMessage()]
    assert sum(int(r.getMessage().rsplit(", ", 1)[1].split()[0]) for r in freed) == expected_freed


def test_each_derived_step_sweeps_its_own_kind(tmp_path: Path, default_store: ReferenceStore) -> None:
    stale_rejoin = write_entry("frames", f"{'d' * 32}-{SCHEMA_VERSION}", cache.REJOIN_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": 0})
    stale_resolve = legacy_entry("frames", f"{'e' * 32}-{SCHEMA_VERSION}-ref-SQ_ab-c-{RESOLVE_VERSION - 1}")
    current_resolve = legacy_entry("frames", f"{'e' * 32}-{SCHEMA_VERSION}-ref-SQ_ab-c-{RESOLVE_VERSION}")
    stale_normalize = write_entry(
        "frames", f"{'9' * 32}-{SCHEMA_VERSION}", cache.NORMALIZE_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION, "NORMALIZE_VERSION": 0}
    )
    loaded = rejoin_split(load(vcf_zip(tmp_path, "one")))
    assert not stale_rejoin.exists()
    assert stale_resolve.exists() and stale_normalize.exists()
    resolved = resolve_ref(loaded, default_store)
    assert not stale_resolve.exists() and current_resolve.exists() and stale_normalize.exists()
    normalize(resolved, default_store)
    assert not stale_normalize.exists()


def test_clear_by_source_removes_its_derived_spills_and_keeps_the_rest(tmp_path: Path, default_store: ReferenceStore) -> None:
    root = cache.cache_dir()
    mine, theirs = vcf_zip(tmp_path, "mine"), vcf_zip(tmp_path, "theirs", positions=(100, 200, 400))
    full_chain(mine, default_store)
    full_chain(theirs, default_store)
    unindexed = legacy_entry("frames", f"{'f' * 32}-{SCHEMA_VERSION}")
    reference_before = tree(root / "reference")
    mine_entries = {name for name, m in manifests().items() if m.source_sha256 == sha256(mine)}
    theirs_entries = {name for name, m in manifests().items() if m.source_sha256 == sha256(theirs)}
    assert len(mine_entries) == len(theirs_entries) == 5

    planned = clear_cache(mine, dry_run=True)
    assert {e.path for e in planned.removed} == mine_entries
    assert planned.bytes_removed == sum(size_on_disk(root / name) for name in mine_entries)
    assert all((root / name).exists() for name in mine_entries)

    report = clear_cache(mine)
    assert {e.path for e in report.removed} == mine_entries
    assert report.bytes_removed == planned.bytes_removed
    assert [e.path for e in report.unknown_origin] == [f"frames/{unindexed.name}"]
    assert not any((root / name).exists() or (root / f"{name}.json").exists() for name in mine_entries)
    assert set(manifests()) == theirs_entries
    assert unindexed.exists()
    assert tree(root / "reference") == reference_before
    # The other source's frames still collect, and the reference store still resolves.
    assert resolve_ref(load(theirs), default_store).calls.collect().height > 0


def test_clear_names_a_deleted_source_by_digest_or_by_path(tmp_path: Path) -> None:
    first, second = vcf_zip(tmp_path, "first"), vcf_zip(tmp_path, "second", positions=(150, 200))
    load(first)
    load(second)
    digest = sha256(first)
    first.unlink()
    second_digest = sha256(second)
    second.unlink()

    by_digest = clear_cache(digest)
    assert {e.source_sha256 for e in by_digest.removed} == {digest} and len(by_digest.removed) == 2
    by_path = clear_cache(second)
    assert {e.source_sha256 for e in by_path.removed} == {second_digest} and len(by_path.removed) == 2
    assert manifests() == {}


def test_clear_by_age_uses_last_use_and_never_touches_unknown_origin(tmp_path: Path) -> None:
    load(vcf_zip(tmp_path, "one"))
    unindexed = legacy_entry("vcf", "b" * 32)
    indexed = set(manifests())

    assert clear_cache(older_than=timedelta(days=1)).removed == ()
    assert set(manifests()) == indexed
    report = clear_cache(older_than=timedelta(0))
    assert {e.path for e in report.removed} == indexed
    assert [(e.path, e.kind, e.bytes) for e in report.unknown_origin] == [(f"vcf/{unindexed.name}", None, size_on_disk(unindexed))]
    assert unindexed.exists() and manifests() == {}
    with pytest.raises(ValueError):
        clear_cache(older_than=timedelta(seconds=-1))


def test_an_entry_being_written_is_refused(tmp_path: Path) -> None:
    root = cache.cache_dir()
    entry = root / "frames" / f"{'7' * 32}-{SCHEMA_VERSION}"
    with cache.writing(entry, cache.FRAMES_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION}, Origin(source_sha256="1" * 64, source_path=None)):
        entry.mkdir()
        (entry / "calls.parquet").write_bytes(b"z" * 10)
        report = clear_cache()
        assert [e.path for e in report.being_written] == [f"frames/{entry.name}"] and report.removed == ()
        assert entry.exists()
    assert [e.path for e in clear_cache().removed] == [f"frames/{entry.name}"]
    assert not entry.exists()


def test_a_writer_waiting_on_a_removed_entry_takes_a_fresh_lock(tmp_path: Path) -> None:
    """A remover unlinks the lock file it holds; a writer that was waiting on it must not proceed on
    the unlinked file while a second writer takes the new one."""
    lock = cache.cache_dir() / "frames" / "x.lock"
    lock.parent.mkdir(parents=True)
    held = cache._acquire(lock, blocking=True)
    got: list[int] = []
    waiter = threading.Thread(target=lambda: got.append(cache._acquire(lock, blocking=True)))
    waiter.start()
    waiter.join(0.3)
    assert waiter.is_alive()
    lock.unlink()
    cache._release(held)
    waiter.join(5)
    assert not waiter.is_alive()
    assert os.fstat(got[0]).st_ino == lock.stat().st_ino
    assert cache._acquire(lock, blocking=False) is None
    cache._release(got[0])


def test_a_manifest_a_newer_version_wrote_neither_breaks_load_nor_is_removed(tmp_path: Path) -> None:
    root = cache.cache_dir()
    newer = {"vcf": "c" * 32, "frames": f"{'c' * 32}-{SCHEMA_VERSION}"}
    for directory, name in newer.items():
        entry = legacy_entry(directory, name)
        manifest = {"index_version": 2, "kind": "liftover_spill", "versions": {"LIFTOVER_VERSION": "1"}, "created": "2026-10-22T00:00:00Z"}
        (entry.parent / f"{name}.json").write_text(json.dumps(manifest))
    load(vcf_zip(tmp_path, "one"))
    assert all((root / directory / name).exists() for directory, name in newer.items())
    report = clear_cache()
    assert {e.path for e in report.unknown_origin} == {f"{d}/{n}" for d, n in newer.items()}
    assert all((root / directory / name).exists() for directory, name in newer.items())
    with pytest.raises(ValueError):
        write_entry("frames", f"{'8' * 32}-{SCHEMA_VERSION}", "liftover_spill", {"SCHEMA_VERSION": SCHEMA_VERSION})


def test_manifest_json_is_iso_utc(tmp_path: Path) -> None:
    load(vcf_zip(tmp_path, "one"))
    for name in manifests():
        raw = json.loads((cache.cache_dir() / f"{name}.json").read_text())
        assert raw["created"].endswith("Z") or raw["created"].endswith("+00:00")
        assert "source_path" in raw
