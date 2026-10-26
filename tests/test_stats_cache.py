"""The cached `scan_stats` pass (F26): a load run again reads the warning counts its first run stored,
keyed by everything the pass depends on (`loaders.vcf.stats_key`), and never a result another input
or another version computed."""

import json
import math
import threading
import time
import zipfile
from pathlib import Path

import polars as pl
import pytest

from just_vcf import cache, load
from just_vcf.bgzf import BgzfWriter
from just_vcf.loaders import vcf
from test_load_vcf import vcf_text, write

RECORDS = [
    "1\t1\t.\tA\t<*>\t.\t.\tEND=50\tGT:GQ:GP\t0/0:40:0,30,40",
    "1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ:GP\t0/1:70000:0,20,30",
    "2\t5\t.\tC\tT\t.\t.\t.\tGT:GQ:GP\t1/1:30:30,20,0",
]
"""A block, an out-of-range GQ (a warning), GP (a float statistic) and two contigs."""


class Spy:
    """Counts the polars-bio stats passes `load` runs, and keeps what each returned and what `to_rows`
    was handed."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.passes: list[dict] = []
        self.handed: list[dict] = []
        scan_stats, to_rows = vcf.scan_stats, vcf.to_rows

        def counted(*args, **kwargs):
            result = scan_stats(*args, **kwargs)
            self.passes.append(result)
            return result

        def seen(records, header, stats, *args, **kwargs):
            self.handed.append(stats)
            return to_rows(records, header, stats, *args, **kwargs)

        monkeypatch.setattr(vcf, "scan_stats", counted)
        monkeypatch.setattr(vcf, "to_rows", seen)

    def runs(self, path: Path, **kwargs) -> tuple[int, tuple]:
        """Load `path`; how many passes that load ran, and what it gave (frames, facets, warnings)."""
        before = len(self.passes)
        loaded = load(path, **kwargs)
        calls, blocks, facets = pl.collect_all([loaded.calls, loaded.blocks, loaded.facets])
        return len(self.passes) - before, (calls, blocks, facets, loaded.warnings)


def same(first: tuple, second: tuple) -> bool:
    return all(a.equals(b) for a, b in zip(first[:3], second[:3], strict=True)) and first[3] == second[3]


def plain(tmp_path: Path) -> Path:
    """Read as it is: no copy, so the key takes the source's digest."""
    return write(tmp_path, vcf_text(RECORDS, fmt=("GT", "GQ", "GP"), info=("END",)))


def zipped(tmp_path: Path, records: list[str] = RECORDS) -> Path:
    """Read through a copy in the cache, so the key takes the copy's digest."""
    path = tmp_path / "calls.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("calls.vcf", vcf_text(records, fmt=("GT", "GQ", "GP"), info=("END",)))
    return path


def stats_files() -> list[Path]:
    return sorted((cache.cache_dir() / "stats").glob(f"*/{vcf.STATS_FILE}"))


@pytest.mark.parametrize("make", [plain, zipped], ids=["as_is", "copy"])
def test_a_load_run_again_reads_the_stats_its_first_run_counted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, make) -> None:
    spy = Spy(monkeypatch)
    path = make(tmp_path)
    passes, first = spy.runs(path)
    assert passes == 1
    assert "1 GQ value(s) outside 0..65535 set to null" in first[3]
    passes, second = spy.runs(path)
    assert passes == 0
    assert same(first, second)
    # The loader is handed the same values, of the same types, as the pass returned.
    counted, (fresh, cached) = spy.passes[0], spy.handed
    assert fresh == cached and counted.keys() == cached.keys()
    for key, value in counted.items():
        assert type(cached[key]) is type(value), key
        assert cached[key] == value or (isinstance(value, float) and math.isnan(value) and math.isnan(cached[key])), key
    manifest = cache.CacheEntry.model_validate_json(stats_files()[0].parent.with_suffix(".json").read_bytes())
    assert manifest.kind == cache.SCAN_STATS
    assert manifest.facet_id == first[2]["facet_id"][0]
    assert manifest.versions == {name: str(value) for name, value in vcf.stats_versions().items()}


@pytest.mark.parametrize(
    "change",
    [
        pytest.param(lambda m: m.setattr(vcf, "STATS_VERSION", vcf.STATS_VERSION + 1), id="STATS_VERSION"),
        pytest.param(lambda m: m.setattr(vcf, "POLARS_BIO_VERSION", vcf.POLARS_BIO_VERSION + ".post1"), id="polars-bio"),
        pytest.param(lambda m: m.setattr(vcf, "_STATS_CODE", "0" * 64), id="code"),
        pytest.param(lambda m: m.setattr(vcf, "SCAN_OPTIONS", {**vcf.SCAN_OPTIONS, "chunk_size": 16}), id="scan options"),
        pytest.param(
            lambda m: m.setattr(vcf, "READER_DEFINITIONS", {**vcf.READER_DEFINITIONS, (4, 3): {**vcf.READER_DEFINITIONS[(4, 3)], "INFO/XX": "Integer"}}),
            id="reader definitions",
        ),
    ],
)
def test_a_changed_key_input_runs_the_pass_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change) -> None:
    spy = Spy(monkeypatch)
    path = plain(tmp_path)
    assert spy.runs(path)[0] == 1
    with monkeypatch.context() as patched:
        change(patched)
        passes, changed = spy.runs(path)
        assert passes == 1
        assert spy.runs(path)[0] == 0
    passes, again = spy.runs(path)
    assert passes == 0 and same(changed, again)
    assert len(stats_files()) == 2


def test_another_sample_is_another_key(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy(monkeypatch)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:70000\t1/1:30", "1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ\t0/0:20\t0/1:80000"]
    path = write(tmp_path, vcf_text(records, samples=("S1", "S2"), fmt=("GT", "GQ")))
    passes, first = spy.runs(path, sample="S1")
    assert passes == 1
    passes, second = spy.runs(path, sample="S2")
    assert passes == 1
    assert first[0]["gq"].to_list() == [None, 20] and second[0]["gq"].to_list() == [30, None]
    assert spy.runs(path, sample="S1")[0] == 0 and spy.runs(path, sample="S2")[0] == 0


def test_a_copy_with_other_bytes_is_counted_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A copy written by code that changed without a `REWRITE_VERSION` bump is still the copy polars-bio
    reads: its bytes, not its name, key the stats."""
    spy = Spy(monkeypatch)
    path = zipped(tmp_path)
    passes, first = spy.runs(path)
    assert passes == 1 and first[2]["contig_map"][0].len() == 2
    (copy,) = (cache.cache_dir() / "vcf").glob("*/input.vcf.gz")
    extra = [*RECORDS, "3\t7\t.\tG\tA\t.\t.\t.\tGT:GQ:GP\t0/1:20:20,0,30"]
    with copy.open("wb") as raw, BgzfWriter(raw) as out:
        out.write(vcf_text(extra, fmt=("GT", "GQ", "GP"), info=("END",)).encode())
    passes, second = spy.runs(path)
    assert passes == 1
    assert second[2]["contig_map"][0].len() == 3 and second[0].height == first[0].height + 1


@pytest.mark.parametrize("damage", [lambda text: text[: len(text) // 2], lambda _: "\udcff not json"], ids=["truncated", "garbage"])
def test_a_stats_file_that_does_not_read_whole_is_counted_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, damage) -> None:
    spy = Spy(monkeypatch)
    path = plain(tmp_path)
    passes, first = spy.runs(path)
    (stored,) = stats_files()
    whole = stored.read_text()
    stored.write_bytes(damage(whole).encode("utf-8", errors="surrogateescape"))
    passes, second = spy.runs(path)
    assert passes == 1 and same(first, second)
    assert stored.read_text() == whole
    assert spy.runs(path)[0] == 0


def test_a_stats_file_for_another_key_is_counted_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy(monkeypatch)
    path = plain(tmp_path)
    spy.runs(path)
    (stored,) = stats_files()
    record = json.loads(stored.read_text())
    stored.write_text(json.dumps({"key": "0" * 64, "stats": {**record["stats"], "out_GQ": 99}}))
    passes, second = spy.runs(path)
    assert passes == 1 and "1 GQ value(s) outside 0..65535 set to null" in second[3]


def test_a_refused_read_is_not_stored_and_the_retry_is_reported_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy(monkeypatch)
    extra = ('##GATKCommandLine.X=<ID=X,Version=1,CommandLineOptions="reference=C:\\Users\\x\\ref.fa">',)
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"], extra_meta=extra))
    passes, first = spy.runs(path)
    assert passes == 1 and any("refused the file" in w for w in first[3])
    passes, second = spy.runs(path)
    assert passes == 0 and same(first, second)
    assert len(stats_files()) == 1


def test_a_bcf_load_run_again_reads_its_stats(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pysam = pytest.importorskip("pysam", reason="pysam (the [reference] extra) writes the BCF fixture")
    text = write(tmp_path, vcf_text(RECORDS, fmt=("GT", "GQ", "GP"), info=("END",)))
    bcf = tmp_path / "calls.bcf"
    with pysam.VariantFile(str(text)) as source, pysam.VariantFile(str(bcf), "wb", header=source.header) as out:
        for record in source:
            out.write(record)
    spy = Spy(monkeypatch)
    passes, first = spy.runs(bcf)
    assert passes == 1
    passes, second = spy.runs(bcf)
    assert passes == 0 and same(first, second)


def test_a_second_session_waits_for_the_first_pass_and_reads_its_result() -> None:
    """The entry's lock is held through the pass: a second loader of the same key waits, then reads
    the first one's result instead of counting (and never sees half a file)."""
    origin = cache.Origin(source_sha256="0" * 64, source_path=None)
    key = "a" * 64
    started, release = threading.Event(), threading.Event()
    computed: list[str] = []
    results: dict[str, dict] = {}

    def first() -> dict:
        computed.append("first")
        started.set()
        release.wait(10)
        return {"n": 1}

    def second() -> dict:
        computed.append("second")
        return {"n": 2}

    one = threading.Thread(target=lambda: results.__setitem__("first", vcf.cached_stats(key, origin, first)))
    one.start()
    assert started.wait(10)
    two = threading.Thread(target=lambda: results.__setitem__("second", vcf.cached_stats(key, origin, second)))
    two.start()
    time.sleep(0.3)
    assert two.is_alive(), "the second loader waits on the entry's lock"
    release.set()
    one.join(10)
    two.join(10)
    assert computed == ["first"]
    assert results == {"first": {"n": 1}, "second": {"n": 1}}


def test_an_older_stats_version_is_swept(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_CACHE_GRACE_DAYS", "0")
    path = plain(tmp_path)
    with monkeypatch.context() as older:
        older.setattr(vcf, "STATS_VERSION", vcf.STATS_VERSION - 1)
        load(path)
    (old,) = stats_files()
    load(path)
    assert not old.exists()
    assert len(stats_files()) == 1
