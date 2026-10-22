"""A derived step run again must not break the frame an earlier run returned.

`rejoin_split`, `resolve_ref`, `normalize` and `mint_vrs_ids` each spill under `frames/<entry>/`. A run
stages its files and commits them as a generation named by the content it was computed from; a later
run with the same content reuses that generation, a run with other rows under the same facet id adds
one beside it, and no run rewrites or removes a committed generation. Every test runs on the per-test
cache `isolated_cache` sets up (conftest).
"""

import multiprocessing
import queue
import time
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pytest

from just_vcf import Loaded, ReferenceStore, load, mint_vrs_ids, normalize, rejoin_split, resolve_ref
from just_vcf import cache
from just_vcf.reference import RESOLVE_VERSION
from just_vcf.schema import SCHEMA_VERSION
from test_cache import default_store, manifests, size_on_disk, vcf_zip, write_entry  # noqa: F401  (default_store is a fixture)

Step = Callable[[Loaded, ReferenceStore], Loaded]

STEPS: dict[str, tuple[str, bool, Step]] = {
    # name: (kind it spills, whether its input is the resolved facet, the call)
    "rejoin_split": (cache.REJOIN_SPILL, False, lambda loaded, store: rejoin_split(loaded)),
    "resolve_ref": (cache.RESOLVE_SPILL, False, lambda loaded, store: resolve_ref(loaded, store)),
    "normalize": (cache.NORMALIZE_SPILL, True, lambda loaded, store: normalize(loaded, store)),
    "mint_vrs_ids": (cache.VRS_SPILL, True, lambda loaded, store: mint_vrs_ids(loaded, store)),
}


def step_input(path: Path, store: ReferenceStore, resolved: bool) -> Loaded:
    loaded = load(path)
    return resolve_ref(loaded, store) if resolved else loaded


def entry_of(kind: str) -> Path:
    ((name, _),) = manifests(kind).items()
    return cache.cache_dir() / name


def generations(entry: Path) -> list[Path]:
    return sorted(child for child in entry.iterdir() if child.is_dir())


def file_stamps(entry: Path) -> dict[Path, tuple[int, int]]:
    """Every file under the entry, with what a rewrite would change: its inode and its mtime."""
    return {path: (path.stat().st_ino, path.stat().st_mtime_ns) for path in entry.rglob("*") if path.is_file()}


@pytest.mark.parametrize("name", STEPS)
def test_the_same_input_again_reuses_the_spill_and_touches_no_file(name: str, tmp_path: Path, default_store: ReferenceStore) -> None:
    kind, resolved, step = STEPS[name]
    given = step_input(vcf_zip(tmp_path, "one"), default_store, resolved)
    first = step(given, default_store)
    entry = entry_of(kind)
    before = file_stamps(entry)
    rows = first.calls.collect()
    manifest = manifests(kind)[f"frames/{entry.name}"]

    second = step(given, default_store)

    assert file_stamps(entry) == before
    assert first.calls.collect().equals(rows)
    assert second.calls.collect().equals(rows)
    assert second.warnings == first.warnings
    assert len(generations(entry)) == 1
    again = manifests(kind)[f"frames/{entry.name}"]
    assert again.created == manifest.created and again.last_used > manifest.last_used
    assert again.size_bytes == size_on_disk(entry)


@pytest.mark.parametrize("name", STEPS)
def test_other_rows_under_the_same_facet_id_leave_the_first_frame_as_it_was(
    name: str, tmp_path: Path, default_store: ReferenceStore
) -> None:
    kind, resolved, step = STEPS[name]
    given = step_input(vcf_zip(tmp_path, "one"), default_store, resolved)
    first = step(given, default_store)
    rows = first.calls.collect()
    # The same facet, one unsplit site fewer: rejoin_split still has its split group to join.
    fewer = Loaded(facets=given.facets, calls=given.calls.filter(pl.col("pos") != 100), blocks=given.blocks)

    second = step(fewer, default_store)

    assert first.calls.collect().equals(rows)
    assert set(second.calls.collect()["pos"]) == set(rows["pos"]) - {100}
    assert len(generations(entry_of(kind))) == 2
    # And each input finds its own generation again.
    assert step(given, default_store).calls.collect().equals(rows)
    assert len(generations(entry_of(kind))) == 2


def test_a_staging_directory_a_dead_writer_left_is_removed(tmp_path: Path, default_store: ReferenceStore) -> None:
    given = load(vcf_zip(tmp_path, "one"))
    first = resolve_ref(given, default_store)
    rows = first.calls.collect()
    entry = entry_of(cache.RESOLVE_SPILL)
    dead = entry / ".staging-1-deadbeef"
    dead.mkdir()
    (dead / "calls-00000.parquet").write_bytes(b"half")

    resolve_ref(given, default_store)

    assert not dead.exists()
    assert [child.name for child in entry.iterdir()] == [generations(entry)[0].name]
    assert first.calls.collect().equals(rows)


def test_a_generation_missing_a_file_is_written_again(tmp_path: Path, default_store: ReferenceStore) -> None:
    given = load(vcf_zip(tmp_path, "one"))
    rows = resolve_ref(given, default_store).calls.collect()
    (generation,) = generations(entry_of(cache.RESOLVE_SPILL))
    next(generation.glob("calls-*.parquet")).unlink()

    again = resolve_ref(given, default_store)

    assert again.calls.collect().equals(rows)
    assert [g.name for g in generations(entry_of(cache.RESOLVE_SPILL))] == [generation.name]


def test_a_spill_in_the_layout_before_generations_is_swept(tmp_path: Path, default_store: ReferenceStore) -> None:
    """Same data versions, no `SPILL_LAYOUT`: files written straight into the entry, which no run of
    this version reads, so the sweep removes the entry."""
    old = write_entry(
        "frames", f"{'e' * 32}-{SCHEMA_VERSION}-ref-SQ_ab-c-{RESOLVE_VERSION}", cache.RESOLVE_SPILL,
        {"SCHEMA_VERSION": SCHEMA_VERSION, "RESOLVE_VERSION": RESOLVE_VERSION},
    )
    resolve_ref(load(vcf_zip(tmp_path, "one")), default_store)
    assert not old.exists()


def _run_chain(frames: tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame], results: "multiprocessing.Queue") -> None:
    """In a child process: every derived step on the facet the parent loaded (passed as its plans, so
    the child does not rewrite the parent's load spill), reporting each step's rows."""
    loaded = Loaded(facets=frames[0], calls=frames[1], blocks=frames[2])
    store = ReferenceStore.open()
    resolved = resolve_ref(loaded, store)
    out = {
        "rejoin_split": rejoin_split(loaded),
        "resolve_ref": resolved,
        "normalize": normalize(resolved, store),
        "mint_vrs_ids": mint_vrs_ids(resolved, store),
    }
    results.put({name: result.calls.collect().write_ipc(None).getvalue() for name, result in out.items()})


def test_two_processes_running_every_step_leave_a_live_frame_working(tmp_path: Path, default_store: ReferenceStore) -> None:
    loaded = load(vcf_zip(tmp_path, "one"))
    resolved = resolve_ref(loaded, default_store)
    live = {
        "rejoin_split": rejoin_split(loaded),
        "resolve_ref": resolved,
        "normalize": normalize(resolved, default_store),
        "mint_vrs_ids": mint_vrs_ids(resolved, default_store),
    }
    expected = {name: frame.calls.collect() for name, frame in live.items()}
    stamps = {kind: file_stamps(entry_of(kind)) for kind, _, _ in STEPS.values()}

    context = multiprocessing.get_context("spawn")
    results = context.Queue()
    children = [context.Process(target=_run_chain, args=((loaded.facets, loaded.calls, loaded.blocks), results)) for _ in range(2)]
    for child in children:
        child.start()
    reported: list[dict] = []
    reads = 0
    deadline = time.monotonic() + 120
    # Results are taken while the children run: a child that has put its result exits only once it is read.
    while len(reported) < len(children):
        assert time.monotonic() < deadline, "children did not report in time"
        assert all(child.exitcode in (None, 0) for child in children), [child.exitcode for child in children]
        for name, frame in live.items():
            assert frame.calls.collect().equals(expected[name]), name
        reads += 1
        try:
            reported.append(results.get(timeout=0.01))
        except queue.Empty:
            pass
    for child in children:
        child.join(60)
        assert child.exitcode == 0

    assert reads > 0
    for report in reported:
        for name, data in report.items():
            assert pl.read_ipc(data).equals(expected[name]), name
    for kind, before in stamps.items():
        entry = entry_of(kind)
        assert file_stamps(entry) == before, kind
        assert len(generations(entry)) == 1, kind
