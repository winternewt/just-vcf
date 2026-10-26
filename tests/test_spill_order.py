"""F27: `load` writes its spill in canonical order once, so the frames it returns read their part files
as they are, with no sort when they are read.

Before, each frame was one unsorted file read back one contig at a time, a filtered scan and a sort per
contig on every read (2,179 of them for a decoy-heavy GRCh38 gVCF). The expected frames here are that
old read, run at test time over the very bytes the loader staged before sorting: the answer `load`
gives must not change, row for row, only its plan.
"""

import random
import shutil
import zipfile
from pathlib import Path

import polars as pl
import pytest
from test_cache import aged, write_entry
from test_load_cgi import LOCI, var_text
from test_load_cgi import write as write_cgi
from test_load_keyvalue import ROWS as KEYVALUE_ROWS
from test_load_keyvalue import write as write_keyvalue
from test_rerun import generations

from just_vcf import cache, load
from just_vcf import loading as loading_module
from just_vcf.schema import SCHEMA_VERSION, canonical_order, canonical_order_per_contig

DECOYS = [f"chrUn_JTFH0100{index:04d}v1_decoy" for index in range(60)]
CONTIGS = ["1", "2", "10", "X", "Y", "MT", "chrUn_KI270302v1", "GL000192.1", *DECOYS]


def decoy_vcf(tmp_path: Path, seed: int = 27) -> Path:
    """A zipped gVCF over 68 contigs, most of them decoys with a row or two, its records out of order:
    calls, multi-allelic calls, records sharing a position and a REF/ALT (a tie the sort must keep in
    file order), and reference blocks."""
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
    for _ in range(600):
        chrom = rng.choice(CONTIGS[:8]) if rng.random() < 0.8 else rng.choice(DECOYS)
        pos = rng.randint(1, 500)
        ref = rng.choice("ACGT")
        alts = [base for base in "ACGT" if base != ref]
        kind = rng.random()
        if kind < 0.25:
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t<*>\t.\t.\tEND={pos + rng.randint(0, 40)}\tGT:GQ\t0/0:{rng.randint(1, 99)}")
        elif kind < 0.35:
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{alts[0]},{alts[1]}\t50\tPASS\t.\tGT:GQ\t1/2:{rng.randint(1, 99)}")
        elif kind < 0.45:
            # The same site twice, told apart only by GQ: a tie in every sort key.
            rows += [f"{chrom}\t{pos}\t.\t{ref}\t{alts[2]}\t50\tPASS\t.\tGT:GQ\t0/1:{gq}" for gq in rng.sample(range(1, 99), 2)]
        else:
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t{rng.choice(alts)}\t50\tPASS\t.\tGT:GQ\t{rng.choice(['0/1', '1/1'])}:{rng.randint(1, 99)}")
    path = tmp_path / "decoys.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("decoys.g.vcf", "\n".join([*meta, columns, *rows]) + "\n")
    return path


SOURCES = {
    "vcf": decoy_vcf,
    "cgi": lambda tmp_path: write_cgi(tmp_path, var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"),
    "keyvalue": lambda tmp_path: write_keyvalue(tmp_path, KEYVALUE_ROWS),
}


@pytest.fixture
def staged(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[Path]:
    """Copies of every unsorted file the loader stages, taken just before it sorts them into parts."""
    copies: list[Path] = []
    real = loading_module._canonical_parts

    def spying(path: Path, name: str) -> list[Path]:
        copy = tmp_path / "staged" / f"{len(copies):03d}-{name}.parquet"
        copy.parent.mkdir(exist_ok=True)
        shutil.copyfile(path, copy)
        copies.append(copy)
        return real(path, name)

    monkeypatch.setattr(loading_module, "_canonical_parts", spying)
    return copies


def read_the_old_way(copies: list[Path], name: str) -> pl.DataFrame:
    """What `load` handed back before F27: each facet's unsorted file sorted one contig at a time, the
    facets in `facet_id` order."""
    frames = []
    for copy in copies:
        if not copy.name.endswith(f"-{name}.parquet"):
            continue
        scan = pl.scan_parquet(copy)
        contigs = scan.select(pl.col("chrom").unique()).collect()["chrom"].drop_nulls().to_list()
        frames.append(canonical_order_per_contig(scan, contigs).collect())
    return pl.concat(sorted(frames, key=lambda frame: frame["facet_id"][0] if frame.height else ""), how="vertical")


@pytest.mark.parametrize("part_rows", [None, 7, 1])
@pytest.mark.parametrize("name", SOURCES)
def test_the_frames_are_the_old_per_contig_read_row_for_row(
    name: str, part_rows: int, staged: list[Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Content, row order, ties and schema, whatever the part bound: one part, parts of a few contigs,
    a part per contig (None keeps the default bound)."""
    if part_rows is not None:
        monkeypatch.setattr(loading_module, "SPILL_PART_ROWS", part_rows)
    loaded = load(SOURCES[name](tmp_path))
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])

    for got, frame in ((calls, "calls"), (blocks, "blocks")):
        expected = read_the_old_way(staged, frame)
        assert got.schema == expected.schema
        assert got.equals(expected), frame
    assert calls.equals(canonical_order(calls.lazy()).collect())
    if name == "vcf":
        # Its records are out of order, so the comparison can tell a sorted frame from an unsorted one.
        (unsorted,) = [pl.read_parquet(copy) for copy in staged if copy.name.endswith("-calls.parquet")]
        assert not unsorted.equals(calls)


def test_the_decoy_vcf_has_ties_and_blocks_on_every_kind_of_contig(staged: list[Path], tmp_path: Path) -> None:
    """The fixture above exercises what it claims: ties, blocks, and the decoys as their own contigs."""
    calls, blocks = pl.collect_all([(loaded := load(decoy_vcf(tmp_path))).calls, loaded.blocks])
    keys = ["chrom", "pos", "end_pos", "source_ref", pl.col("alts").list.join(",")]
    assert calls.select(keys).is_duplicated().any()
    assert blocks.height and calls.height
    assert len(set(calls["chrom"]) | set(blocks["chrom"])) > 50


@pytest.mark.parametrize("name", SOURCES)
def test_a_loaded_facet_reads_its_spill_in_one_scan_with_no_sort(name: str, tmp_path: Path) -> None:
    """The plan of each frame is a scan of the part files: no sort node, one scan whatever the number
    of contigs (68 in the VCF). The decoys with a row or two share a part with the rest (F24)."""
    loaded = load(SOURCES[name](tmp_path))
    facets = loaded.facets.collect().height
    # A key=value table holds no blocks: its `blocks` is an empty frame, not a scan of the spill.
    for frame, scans in ((loaded.calls, facets), (loaded.blocks, 0 if name == "keyvalue" else facets)):
        plan = frame.explain()
        assert "SORT" not in plan
        assert plan.count("SCAN") == scans
    for entry in [cache.cache_dir() / f"frames/{facet_id}-{SCHEMA_VERSION}" for facet_id in loaded.facets.collect()["facet_id"]]:
        (generation,) = generations(entry)
        assert len(list(generation.glob("calls-*.parquet"))) == 1
        assert len(list(generation.glob("blocks-*.parquet"))) == 1


def test_a_contig_over_the_part_bound_is_a_part_of_its_own(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Parts hold whole contigs: none is cut across two files, and a part over the bound holds one."""
    monkeypatch.setattr(loading_module, "SPILL_PART_ROWS", 7)
    loaded = load(decoy_vcf(tmp_path))
    (facet_id,) = loaded.facets.collect()["facet_id"]
    (generation,) = generations(cache.cache_dir() / f"frames/{facet_id}-{SCHEMA_VERSION}")
    parts = sorted(generation.glob("calls-*.parquet"))
    assert len(parts) > 1
    seen: set[str] = set()
    for part in parts:
        rows = pl.read_parquet(part)
        chroms = set(rows["chrom"])
        assert not chroms & seen
        seen |= chroms
        assert rows.height <= 7 or len(chroms) == 1


def test_a_spill_in_the_unsorted_layout_is_not_read(staged: list[Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A generation an older version committed (one unsorted file per frame, keyed by its bytes and
    the facet id alone) sits in the same entry, used within the grace period so the sweep keeps it.
    A load of the same rows writes a sorted generation beside it and never reads the old one."""
    path = decoy_vcf(tmp_path)
    first = load(path)
    (facet_id,) = first.facets.collect()["facet_id"]
    expected = pl.collect_all([first.calls, first.blocks])
    calls_copy, blocks_copy = staged

    monkeypatch.setenv("JUST_VCF_CACHE_DIR", str(tmp_path / "older-cache"))
    old_versions = {"SCHEMA_VERSION": SCHEMA_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    entry = write_entry("frames", f"{facet_id}-{SCHEMA_VERSION}", cache.FRAMES_SPILL, old_versions)
    (entry / "calls.parquet").unlink()
    # The older version's generation: the unsorted bytes, under the key that version gave them.
    old_key = cache.content_key([calls_copy, blocks_copy], {"facet_id": facet_id})
    old = entry / old_key
    old.mkdir()
    shutil.copyfile(calls_copy, old / "calls.parquet")
    shutil.copyfile(blocks_copy, old / "blocks.parquet")
    (old / "result.json").write_text(cache.Generation(files={"calls": ("calls.parquet",), "blocks": ("blocks.parquet",)}).model_dump_json())

    again = load(path)

    got = pl.collect_all([again.calls, again.blocks])
    assert all(a.equals(b) for a, b in zip(got, expected, strict=True))
    assert old in generations(entry) and len(generations(entry)) == 2
    assert "calls.parquet" not in again.calls.explain()


def test_load_sweeps_a_spill_in_the_unsorted_layout(tmp_path: Path) -> None:
    """A load spill whose manifest lacks `FRAMES_LAYOUT`, unused past the grace period, is removed."""
    before = aged(write_entry("frames", f"{'c' * 32}-{SCHEMA_VERSION}", cache.FRAMES_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}))
    newer = aged(write_entry(
        "frames", f"{'d' * 32}-{SCHEMA_VERSION}", cache.FRAMES_SPILL,
        {"SCHEMA_VERSION": SCHEMA_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT, "FRAMES_LAYOUT": loading_module.FRAMES_LAYOUT + 1},
    ))
    load(decoy_vcf(tmp_path))
    assert not before.exists() and newer.exists()
