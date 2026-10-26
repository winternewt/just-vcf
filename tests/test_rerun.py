"""`load` or a derived step run again must not break the frame an earlier run returned.

`load`, `rejoin_split`, `resolve_ref`, `normalize`, `mint_vrs_ids` and `liftover` each spill under
`frames/<entry>/`. A run
stages its files and commits them as a generation named by the content it was computed from; a later
run with the same content reuses that generation, a run with other rows under the same facet id adds
one beside it, and no run rewrites or removes a committed generation. Every test runs on the per-test
cache `isolated_cache` sets up (conftest).
"""

import faulthandler
import multiprocessing
import multiprocessing.connection
import os
import pickle
import shutil
import signal
import time
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pytest
from test_cache import SEQ, _schema_ahead, aged, legacy_entry, manifests, sha256, size_on_disk, vcf_zip, write_entry
from test_load_cgi import LOCI, var_text
from test_load_cgi import write as write_cgi
from test_load_keyvalue import ROWS as KEYVALUE_ROWS
from test_load_keyvalue import write as write_keyvalue

from just_vcf import (
    Loaded,
    ReferenceStore,
    cache,
    clear_cache,
    liftover,
    load,
    mint_vrs_ids,
    normalize,
    rejoin_split,
    resolve_ref,
)
from just_vcf import loading as loading_module
from just_vcf.liftover import LIFT_VERSION
from just_vcf.loaders import vcf as vcf_loader
from just_vcf.reference import RESOLVE_VERSION
from just_vcf.schema import SCHEMA_VERSION


@pytest.fixture
def default_store(tmp_path: Path) -> ReferenceStore:
    """The reference store at its default place inside the cache dir, holding `SEQ` as contig 1: a
    child process opens it with `ReferenceStore.open()`."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">1\n" + "\n".join(SEQ[i : i + 60] for i in range(0, len(SEQ), 60)) + "\n")
    store = ReferenceStore.open()
    store.add_fasta(fasta)
    return store


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
    path = vcf_zip(tmp_path, "one")
    given = step_input(path, default_store, resolved)
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
    # A fresh load of the same file, as another session would make, finds the same generation.
    fresh = step(step_input(path, default_store, resolved), default_store)
    assert fresh.calls.collect().equals(rows) and fresh.warnings == first.warnings
    assert len(generations(entry)) == 1 and file_stamps(entry) == before


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
    old = aged(write_entry(
        "frames", f"{'e' * 32}-{SCHEMA_VERSION}-ref-SQ_ab-c-{RESOLVE_VERSION}", cache.RESOLVE_SPILL,
        {"SCHEMA_VERSION": SCHEMA_VERSION, "RESOLVE_VERSION": RESOLVE_VERSION},
    ))
    resolve_ref(load(vcf_zip(tmp_path, "one")), default_store)
    assert not old.exists()


CHILD_DEADLINE = 600.0
"""Seconds the children of one two-process test may take. Alone they take a few. In the full suite at
a load average near 100, one small VCF `load` in `tests/test_callset.py` took 488 s (2026-10-23), so
the 120 s this used to allow can pass while the children still work. The deadline is what turns a
hang into a failure, not a timing claim."""


def _child(target: Callable, args: tuple, result: Path, stack: Path) -> None:
    """In a spawned child: run `target(*args)` and write what it returns to `result`.

    The result goes to a file, never a pipe: a child that has put more onto a `multiprocessing.Queue`
    than the pipe holds (64 KiB; one report here is about 57 KiB) cannot exit until the parent reads
    it, and a parent whose loop failed never does. Such a child sits in its queue's feeder thread
    forever, and pytest's own exit joins it forever (multiprocessing joins every non-daemon child).
    SIGUSR1 dumps every thread's stack to `stack`: the parent sends it when the deadline passes.
    """
    with stack.open("w") as handle:
        faulthandler.register(signal.SIGUSR1, file=handle, all_threads=True)
        data = pickle.dumps(target(*args))
        written = result.with_name(f"{result.name}.part")
        written.write_bytes(data)
        os.replace(written, result)


def _stacks(children: list, stacks: list[Path]) -> str:
    """Every live child's stack, from the dump SIGUSR1 asks of it (`_child`). A child still starting,
    before `_child` registered the dump, has no stack file, and the signal ends it."""
    alive = [(child, stack) for child, stack in zip(children, stacks, strict=True) if child.is_alive()]
    for child, _ in alive:
        os.kill(child.pid, signal.SIGUSR1)
    until = time.monotonic() + 5
    while time.monotonic() < until and not all(stack.is_file() and stack.stat().st_size for _, stack in alive):
        time.sleep(0.1)
    time.sleep(0.2)
    dumps = [(child, stack.read_text() if stack.is_file() else "") for child, stack in alive]
    return "\n".join(f"child {child.pid}:\n{dump or '(no stack written: still starting)'}" for child, dump in dumps)


def run_children(target: Callable, args: tuple, check: Callable[[], None], workdir: Path, count: int = 2) -> list:
    """Start `count` spawned children on `target(*args)`, running `check` until each has exited, and
    return what each one's `target` returned.

    Fails within `CHILD_DEADLINE` instead of hanging: past it, the message carries every live child's
    stack. A child left alive by any failure is killed before the test ends, so none outlives it.
    """
    context = multiprocessing.get_context("spawn")
    workdir.mkdir(exist_ok=True)
    results = [workdir / f"child-{n}.result" for n in range(count)]
    stacks = [workdir / f"child-{n}.stack" for n in range(count)]
    children = [
        context.Process(target=_child, args=(target, args, results[n], stacks[n]), daemon=True) for n in range(count)
    ]
    try:
        for child in children:
            child.start()
        reads = 0
        deadline = time.monotonic() + CHILD_DEADLINE
        while any(child.exitcode is None for child in children):
            if time.monotonic() > deadline:
                pytest.fail(f"children still running after {CHILD_DEADLINE:.0f} s:\n{_stacks(children, stacks)}")
            assert all(child.exitcode in (None, 0) for child in children), [child.exitcode for child in children]
            check()
            reads += 1
            multiprocessing.connection.wait([child.sentinel for child in children if child.exitcode is None], timeout=0.01)
        assert [child.exitcode for child in children] == [0] * count
        assert reads > 0
        return [pickle.loads(result.read_bytes()) for result in results]
    finally:
        for child in children:
            if child.is_alive():
                child.kill()
                child.join(10)


def _run_chain(frames: tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]) -> dict[str, bytes]:
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
    return {name: result.calls.collect().write_ipc(None).getvalue() for name, result in out.items()}


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

    def check() -> None:
        for name, frame in live.items():
            assert frame.calls.collect().equals(expected[name]), name

    reported = run_children(_run_chain, ((loaded.facets, loaded.calls, loaded.blocks),), check, tmp_path / "children")

    for report in reported:
        for name, data in report.items():
            assert pl.read_ipc(data).equals(expected[name]), name
    for kind, before in stamps.items():
        entry = entry_of(kind)
        assert file_stamps(entry) == before, kind
        assert len(generations(entry)) == 1, kind


# liftover: the same rules, over two file groups (the lifted calls and the report). The target is SEQ
# with two bases in front, a second collection added once REF is read, and the chain maps contig 1 onto
# it two bases on.
TARGET_SEQ = "GG" + SEQ


@pytest.fixture
def lift_chain(tmp_path: Path) -> Path:
    chain = tmp_path / "seq_to_target.chain"
    chain.write_text(f"chain 1000 1 {len(SEQ)} + 0 {len(SEQ)} 1 {len(TARGET_SEQ)} + 2 {len(TARGET_SEQ)} 1\n{len(SEQ)}\n\n")
    return chain


def add_target(store: ReferenceStore, tmp_path: Path) -> None:
    fasta = tmp_path / "target.fa"
    fasta.write_text(">1\n" + "\n".join(TARGET_SEQ[i : i + 60] for i in range(0, len(TARGET_SEQ), 60)) + "\n")
    store.add_fasta(fasta)


def lifted_frames(lifted) -> tuple[pl.DataFrame, pl.DataFrame]:
    return lifted.loaded.calls.collect(), lifted.report.collect()


def test_liftover_again_reuses_its_spill_and_touches_no_file(tmp_path: Path, default_store: ReferenceStore, lift_chain: Path) -> None:
    path = vcf_zip(tmp_path, "one")
    resolved = resolve_ref(load(path), default_store)
    add_target(default_store, tmp_path)
    first = liftover(resolved, lift_chain, default_store)
    entry = entry_of(cache.LIFT_SPILL)
    before = file_stamps(entry)
    calls, report = lifted_frames(first)
    manifest = manifests(cache.LIFT_SPILL)[f"frames/{entry.name}"]
    assert calls.height > 0 and report["lift_reason"].null_count() == calls.height

    second = liftover(resolved, lift_chain, default_store)

    assert file_stamps(entry) == before
    assert lifted_frames(first)[0].equals(calls) and lifted_frames(first)[1].equals(report)
    assert lifted_frames(second)[0].equals(calls) and lifted_frames(second)[1].equals(report)
    assert (second.reasons, second.loaded.warnings) == (first.reasons, first.loaded.warnings)
    assert len(generations(entry)) == 1
    again = manifests(cache.LIFT_SPILL)[f"frames/{entry.name}"]
    assert again.created == manifest.created and again.last_used > manifest.last_used
    assert again.size_bytes == size_on_disk(entry)
    assert again.versions == {"SCHEMA_VERSION": SCHEMA_VERSION, "LIFT_VERSION": str(LIFT_VERSION), "SPILL_LAYOUT": str(cache.SPILL_LAYOUT)}
    assert (again.facet_id, again.parent_facet_id) == (first.loaded.facets.collect()["facet_id"][0], resolved.facets.collect()["facet_id"][0])


def test_liftover_on_other_rows_leaves_the_first_frame_as_it_was(tmp_path: Path, default_store: ReferenceStore, lift_chain: Path) -> None:
    resolved = resolve_ref(load(vcf_zip(tmp_path, "one")), default_store)
    add_target(default_store, tmp_path)
    first = liftover(resolved, lift_chain, default_store)
    calls, report = lifted_frames(first)
    fewer = Loaded(facets=resolved.facets, calls=resolved.calls.filter(pl.col("pos") != 100), blocks=resolved.blocks)

    second = liftover(fewer, lift_chain, default_store)

    assert lifted_frames(first)[0].equals(calls) and lifted_frames(first)[1].equals(report)
    assert set(second.report.collect()["pos"]) == set(report["pos"]) - {100}
    assert len(generations(entry_of(cache.LIFT_SPILL))) == 2
    assert lifted_frames(liftover(resolved, lift_chain, default_store))[0].equals(calls)
    assert len(generations(entry_of(cache.LIFT_SPILL))) == 2


def test_a_lift_spill_from_before_the_index_is_swept(tmp_path: Path, default_store: ReferenceStore, lift_chain: Path) -> None:
    """`liftover` wrote `<facet_id>-<SCHEMA_VERSION>-lift` with no manifest; no run of this version reads it."""
    old = cache.cache_dir() / "frames" / f"{'e' * 32}-{SCHEMA_VERSION}-lift"
    old.mkdir(parents=True)
    (old / "calls-00000.parquet").write_bytes(b"x" * 100)
    aged(old)
    resolved = resolve_ref(load(vcf_zip(tmp_path, "one")), default_store)
    assert old.exists()
    add_target(default_store, tmp_path)
    liftover(resolved, lift_chain, default_store)
    assert not old.exists()


def _run_lift(frames: tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame], chain: Path) -> list[bytes]:
    """In a child process: `liftover` on the facet the parent resolved, reporting its rows and report."""
    lifted = liftover(Loaded(facets=frames[0], calls=frames[1], blocks=frames[2]), chain, ReferenceStore.open())
    return [frame.write_ipc(None).getvalue() for frame in lifted_frames(lifted)]


def test_two_processes_lifting_leave_a_live_lift_working(tmp_path: Path, default_store: ReferenceStore, lift_chain: Path) -> None:
    resolved = resolve_ref(load(vcf_zip(tmp_path, "one")), default_store)
    add_target(default_store, tmp_path)
    live = liftover(resolved, lift_chain, default_store)
    calls, report = lifted_frames(live)
    entry = entry_of(cache.LIFT_SPILL)
    before = file_stamps(entry)

    def check() -> None:
        assert lifted_frames(live)[0].equals(calls) and lifted_frames(live)[1].equals(report)

    reported = run_children(_run_lift, ((resolved.facets, resolved.calls, resolved.blocks), lift_chain), check, tmp_path / "children")

    for got_calls, got_report in reported:
        assert pl.read_ipc(got_calls).equals(calls) and pl.read_ipc(got_report).equals(report)
    assert file_stamps(entry) == before and len(generations(entry)) == 1


# load: the spill a load writes for each facet, of every loader that spills (VCF, Complete Genomics,
# key=value), follows the same rules.
SOURCES: dict[str, Callable[[Path], Path]] = {
    "vcf": lambda tmp_path: vcf_zip(tmp_path, "one"),
    "cgi": lambda tmp_path: write_cgi(tmp_path, var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"),
    "keyvalue": lambda tmp_path: write_keyvalue(tmp_path, KEYVALUE_ROWS),
}


def all_frames(loaded: Loaded) -> list[pl.DataFrame]:
    return pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])


def same(got: list[pl.DataFrame], expected: list[pl.DataFrame]) -> bool:
    return all(a.equals(b) for a, b in zip(got, expected, strict=True))


def load_entries() -> list[Path]:
    return [cache.cache_dir() / name for name in sorted(manifests(cache.FRAMES_SPILL))]


@pytest.mark.parametrize("name", SOURCES)
def test_loading_again_reuses_the_spill_and_touches_no_file(name: str, tmp_path: Path) -> None:
    path = SOURCES[name](tmp_path)
    first = load(path)
    rows = all_frames(first)
    entries = load_entries()
    assert len(entries) == rows[0].height  # one spill per facet
    before = {entry: file_stamps(entry) for entry in entries}
    made = {entry.name: manifests(cache.FRAMES_SPILL)[f"frames/{entry.name}"] for entry in entries}

    second = load(path)

    assert {entry: file_stamps(entry) for entry in entries} == before
    assert same(all_frames(first), rows) and same(all_frames(second), rows)
    assert second.warnings == first.warnings
    for entry in entries:
        assert len(generations(entry)) == 1
        again = manifests(cache.FRAMES_SPILL)[f"frames/{entry.name}"]
        assert again.created == made[entry.name].created and again.last_used > made[entry.name].last_used
        assert again.size_bytes == size_on_disk(entry)
        assert again.versions == {"SCHEMA_VERSION": SCHEMA_VERSION, "SPILL_LAYOUT": str(cache.SPILL_LAYOUT), "FRAMES_LAYOUT": str(loading_module.FRAMES_LAYOUT)}


def test_a_load_spill_in_the_layout_before_generations_is_swept(tmp_path: Path) -> None:
    """A load spill indexed before generations: its manifest has no `SPILL_LAYOUT`, so `load` sweeps it."""
    old = aged(write_entry("frames", f"{'e' * 32}-{SCHEMA_VERSION}", cache.FRAMES_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION}))
    load(vcf_zip(tmp_path, "one"))
    assert not old.exists()


def _run_load(path: Path) -> list[bytes]:
    """In a child process: load the file, reporting its frames."""
    return [frame.write_ipc(None).getvalue() for frame in all_frames(load(path))]


@pytest.mark.parametrize("name", SOURCES)
def test_two_processes_loading_leave_a_live_load_working(name: str, tmp_path: Path) -> None:
    path = SOURCES[name](tmp_path)
    live = load(path)
    rows = all_frames(live)
    entries = load_entries()
    assert len(entries) == rows[0].height  # one spill per facet
    before = {entry: file_stamps(entry) for entry in entries}

    def check() -> None:
        assert same(all_frames(live), rows), "the live frame changed"

    reported = run_children(_run_load, (path,), check, tmp_path / "children")

    for frames in reported:
        assert same([pl.read_ipc(data) for data in frames], rows)
    assert {entry: file_stamps(entry) for entry in entries} == before
    assert all(len(generations(entry)) == 1 for entry in entries)


def test_a_live_load_reads_the_spill_never_the_rewrite_copy(tmp_path: Path) -> None:
    """The zipped VCF is read through a copy under `vcf/`; once `load` returns, no frame reads it."""
    loaded = load(vcf_zip(tmp_path, "one"))
    rows = all_frames(loaded)
    ((copy, _),) = manifests(cache.REWRITE_COPY).items()
    shutil.rmtree(cache.cache_dir() / copy)
    assert same(all_frames(loaded), rows)


@pytest.mark.parametrize("name", ["cgi", "keyvalue"])
def test_clearing_a_cgi_or_keyvalue_source_removes_every_entry_from_it(name: str, tmp_path: Path) -> None:
    mine = SOURCES[name](tmp_path)
    theirs = vcf_zip(tmp_path, "theirs")
    load(mine)
    load(theirs)
    from_mine = {entry for entry, m in manifests().items() if m.source_sha256 == sha256(mine)}
    work = {entry for entry, m in manifests(cache.WORK_COPY).items()}
    assert work and work <= from_mine and {m.kind for e, m in manifests().items() if e in from_mine} == {cache.WORK_COPY, cache.FRAMES_SPILL}
    # The working copy is gone once the load returned; its manifest stays, at size 0.
    assert all(size_on_disk(cache.cache_dir() / entry) == 0 == manifests()[entry].size_bytes for entry in work)
    # A copy a loader that died left behind is found by its source.
    (cache.cache_dir() / next(iter(work)) / "payload.tsv").write_bytes(b"left by a dead loader")

    report = clear_cache(mine)

    assert {e.path for e in report.removed} == from_mine
    assert not any((cache.cache_dir() / entry).exists() for entry in from_mine)
    assert {m.source_sha256 for m in manifests().values()} == {sha256(theirs)}


def test_a_working_copy_from_before_the_index_is_swept(tmp_path: Path) -> None:
    stale = {directory: aged(legacy_entry(directory, "a" * 32)) for directory in ("keyvalue", "cgi")}
    load(vcf_zip(tmp_path, "one"))
    assert not any(entry.exists() for entry in stale.values())


# Two installs at two versions sharing one cache dir (the 2026-10-23 incident): the newer one's sweep
# removed a generation the older one's live frame read, and its next `collect` failed with
# FileNotFoundError on calls-00000.parquet. Here the parent is the older install; the child is a newer
# one, its version keys moved ahead before it loads.


def _load_as_newer(path: Path) -> list[bytes]:
    """In a child process: `load` as an install whose `SCHEMA_VERSION` and `REWRITE_VERSION` are ahead
    of the parent's, so every frames spill and copy the parent wrote is an older version's to it."""
    ahead = _schema_ahead()
    loading_module.SCHEMA_VERSION = ahead
    loading_module._FRAMES_VERSIONS = {**loading_module._FRAMES_VERSIONS, "SCHEMA_VERSION": ahead}
    vcf_loader.REWRITE_VERSION += 1
    return [frame.write_ipc(None).getvalue() for frame in all_frames(load(path))]


def test_a_newer_install_sweeping_leaves_an_older_installs_live_frame_working(tmp_path: Path) -> None:
    path = vcf_zip(tmp_path, "one")
    live = load(path)
    rows = all_frames(live)
    (live_entry,) = load_entries()
    before = file_stamps(live_entry)
    # Another file this install loaded, and nothing has used since: past the grace period, the newer
    # install retires it.
    load(vcf_zip(tmp_path, "two", positions=(150, 200)))
    (untouched,) = set(load_entries()) - {live_entry}
    aged(untouched)

    def check() -> None:
        assert same(all_frames(live), rows), "the live frame changed"

    (newer,) = run_children(_load_as_newer, (path,), check, tmp_path / "children", count=1)

    assert same(all_frames(live), rows)
    assert file_stamps(live_entry) == before
    assert not untouched.exists() and not untouched.with_name(f"{untouched.name}.json").exists()
    # The newer install spilled the same facet under its own version, beside the older one's.
    assert set(load_entries()) == {live_entry, live_entry.with_name(live_entry.name.replace(SCHEMA_VERSION, _schema_ahead()))}
    assert pl.read_ipc(newer[1]).equals(rows[1])
