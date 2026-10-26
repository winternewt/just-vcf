"""F27 for `rejoin_split`: its spill is written in canonical order once, as `load`'s is, so the calls it
returns read their part files as they are, with no sort when they are read.

Before, the rejoined calls were one unsorted file read back one contig at a time (the contigs of the
parent's `contig_map`), a filtered scan and a sort per contig on every read. The expected calls here are
that old read, run at test time over the very bytes `rejoin_split` staged before sorting: its answer
must not change, row for row, only its plan.
"""

import random
import shutil
import zipfile
from pathlib import Path

import polars as pl
import pytest
from test_cache import aged, write_entry
from test_rerun import generations
from test_spill_order import CONTIGS, DECOYS, canonical_order_per_contig

from just_vcf import cache, load, rejoin_split
from just_vcf import loading as loading_module
from just_vcf import rejoin as rejoin_module
from just_vcf.schema import CALLS, SCHEMA_VERSION, canonical_order


def split_vcf(tmp_path: Path, seed: int = 270) -> Path:
    """A zipped gVCF over 68 contigs, most of them decoys with a row or two, its records out of order:
    sites written one row per ALT that re-join, split groups that cannot (REFs that are not prefixes of
    one another), plain calls, records sharing a site told apart only by GQ (a tie the sort must keep
    in file order), and reference blocks."""
    rng = random.Random(seed)
    meta = [
        "##fileformat=VCFv4.2",
        *[f"##contig=<ID={name},length=100000>" for name in CONTIGS],
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    ]
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    rows = []
    taken: set[tuple[str, int]] = set()
    for _ in range(600):
        chrom = rng.choice(CONTIGS[:8]) if rng.random() < 0.8 else rng.choice(DECOYS)
        pos = rng.randint(1, 2000)
        if (chrom, pos) in taken:
            continue
        taken.add((chrom, pos))
        ref = rng.choice("ACGT")
        alts = [base for base in "ACGT" if base != ref]
        gq = rng.randint(1, 99)
        kind = rng.random()
        if kind < 0.2:
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t<*>\t.\t.\tEND={pos}\tGT:GQ\t0/0:{gq}")
        elif kind < 0.45:
            # One het site of two ALTs, split as `bcftools norm -m-` writes it: re-joins.
            rows += [f"{chrom}\t{pos}\t.\t{ref}\t{alts[0]}\t50\tPASS\t.\tGT:GQ\t1/0:{gq}", f"{chrom}\t{pos}\t.\t{ref}\t{alts[1]}\t50\tPASS\t.\tGT:GQ\t0/1:{gq}"]
        elif kind < 0.55:
            # Flagged (an ALT beside an all-reference row), but the REFs are not prefixes: kept as written.
            other = rng.choice(alts)
            rows += [f"{chrom}\t{pos}\t.\t{ref}\t{alts[0]}\t50\tPASS\t.\tGT:GQ\t0/1:{gq}", f"{chrom}\t{pos}\t.\t{other}\t{ref}\t50\tPASS\t.\tGT:GQ\t0/0:{gq}"]
        elif kind < 0.65:
            # The same site twice, told apart only by GQ: a tie in every sort key, and not a split pattern.
            rows += [f"{chrom}\t{pos}\t.\t{ref}\t{alts[2]}\t50\tPASS\t.\tGT:GQ\t0/1:{value}" for value in rng.sample(range(1, 99), 2)]
        else:
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{rng.choice(alts)}\t50\tPASS\t.\tGT:GQ\t{rng.choice(['0/1', '1/1'])}:{gq}")
    rng.shuffle(rows)
    path = tmp_path / "split.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("split.g.vcf", "\n".join([*meta, columns, *rows]) + "\n")
    return path


@pytest.fixture
def staged(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[Path]:
    """Copies of every unsorted file `rejoin_split` stages, taken when it keys the spill by them, before
    anything sorts them. Spied on the key, which every layout computes from the staged bytes, so the
    oracle reads the same bytes whichever layout writes the generation."""
    copies: list[Path] = []
    real = cache.content_key

    def spying(files: list[Path], facet: dict) -> str:
        if [path.name for path in files] == ["calls.parquet"] and facet.get("relation") is None:
            copy = tmp_path / "staged" / f"{len(copies):03d}-calls.parquet"
            copy.parent.mkdir(exist_ok=True)
            shutil.copyfile(files[0], copy)
            copies.append(copy)
        return real(files, facet)

    monkeypatch.setattr(cache, "content_key", spying)
    return copies


def read_the_old_way(copy: Path, parent: pl.LazyFrame) -> pl.DataFrame:
    """What `rejoin_split` handed back before: its unsorted file sorted one contig at a time, the
    contigs being the canonical names in the parent facet's `contig_map`."""
    contigs = [item["chrom"] for item in parent.collect().row(0, named=True)["contig_map"] or []]
    return canonical_order_per_contig(pl.scan_parquet(copy), contigs).collect()


@pytest.mark.parametrize("part_rows", [None, 7, 1])
def test_rejoined_calls_are_the_old_per_contig_read_row_for_row(
    part_rows: int | None, staged: list[Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Content, row order, ties and schema, whatever the part bound: one part, parts of a few contigs,
    a part per contig (None keeps the default bound)."""
    if part_rows is not None:
        monkeypatch.setattr(loading_module, "SPILL_PART_ROWS", part_rows)
    loaded = load(split_vcf(tmp_path))
    calls = rejoin_split(loaded).calls.collect()

    (copy,) = staged
    expected = read_the_old_way(copy, loaded.facets)
    assert calls.schema == expected.schema == CALLS
    assert calls.equals(expected)
    assert calls.equals(canonical_order(calls.lazy()).collect())
    # The staged file is the parent's unflagged rows, then the re-joined ones: not in canonical order,
    # so the comparison can tell a sorted frame from an unsorted one.
    assert not pl.read_parquet(copy).equals(calls)


def test_the_split_vcf_joins_keeps_and_ties_on_every_kind_of_contig(tmp_path: Path) -> None:
    """The fixture above exercises what it claims: groups re-joined and groups kept, ties, the decoys
    as their own contigs, and rows of the parent's that `rejoin_split` passes on unflagged."""
    rejoined = rejoin_split(load(split_vcf(tmp_path)))
    calls = rejoined.calls.collect()
    assert any(" re-joined from " in warning for warning in rejoined.warnings)
    assert any("left as written: REFs that are not prefixes" in warning for warning in rejoined.warnings)
    keys = ["chrom", "pos", "end_pos", "source_ref", pl.col("alts").list.join(",")]
    assert calls.select(keys).is_duplicated().any()
    assert calls["split_signal"].is_not_null().any() and calls["split_signal"].is_null().any()
    assert len(set(calls["chrom"])) > 50


def test_rejoined_calls_read_their_spill_in_one_scan_with_no_sort(tmp_path: Path) -> None:
    """The calls' plan is a scan of the part files: no sort node, one scan whatever the number of
    contigs (68 here). The decoys with a row or two share a part with the rest (F24)."""
    rejoined = rejoin_split(load(split_vcf(tmp_path)))
    plan = rejoined.calls.explain()
    assert "SORT" not in plan
    assert plan.count("SCAN") == 1
    (facet_id,) = rejoined.facets.collect()["facet_id"]
    (generation,) = generations(cache.cache_dir() / f"frames/{facet_id}-{SCHEMA_VERSION}")
    assert [path.name for path in generation.glob("*.parquet")] == ["calls-00000.parquet"]


def test_a_rejoin_spill_in_the_unsorted_layout_is_not_read(staged: list[Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A generation an older version committed (one unsorted file, keyed by its bytes and the parent's
    facet row alone) sits in the same entry, used within the grace period so the sweep keeps it. A
    rejoin of the same rows writes a sorted generation beside it and never reads the old one."""
    path = split_vcf(tmp_path)
    loaded = load(path)
    first = rejoin_split(loaded)
    (facet_id,) = first.facets.collect()["facet_id"]
    expected = first.calls.collect()
    (copy,) = staged
    parent = loaded.facets.collect().row(0, named=True)

    monkeypatch.setenv("JUST_VCF_CACHE_DIR", str(tmp_path / "older-cache"))
    old_versions = {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": rejoin_module.REJOIN_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    entry = write_entry("frames", f"{facet_id}-{SCHEMA_VERSION}", cache.REJOIN_SPILL, old_versions)
    (entry / "calls.parquet").unlink()
    # The older version's generation: the unsorted bytes, under the key that version gave them.
    old = entry / cache.content_key([copy], parent)
    old.mkdir()
    shutil.copyfile(copy, old / "calls.parquet")
    (old / "result.json").write_text(cache.Generation(files={"calls": ("calls.parquet",)}).model_dump_json())

    again = rejoin_split(load(path))

    assert again.calls.collect().equals(expected)
    assert old in generations(entry) and len(generations(entry)) == 2
    assert "calls.parquet" not in again.calls.explain()


def test_rejoin_sweeps_a_spill_in_the_unsorted_layout(tmp_path: Path) -> None:
    """A rejoin spill whose manifest lacks `REJOIN_LAYOUT`, unused past the grace period, is removed;
    one a newer layout wrote stays."""
    current = {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": rejoin_module.REJOIN_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    before = aged(write_entry("frames", f"{'c' * 32}-{SCHEMA_VERSION}", cache.REJOIN_SPILL, current))
    newer = aged(write_entry("frames", f"{'d' * 32}-{SCHEMA_VERSION}", cache.REJOIN_SPILL, {**current, "REJOIN_LAYOUT": rejoin_module.REJOIN_LAYOUT + 1}))
    rejoin_split(load(split_vcf(tmp_path)))
    assert not before.exists() and newer.exists()
