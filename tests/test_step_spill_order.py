"""F24: `resolve_ref`, `normalize` and `liftover` read their input and write their spill a group of whole
contigs at a time (`loading.canonical_batches`), so a header full of decoys costs a few parts, where it
cost a filtered scan, a pass and a file per contig with rows.

The oracle is each step run with a part bound of one row: every contig is then a group of its own,
read in its own scan, passed through its own pass and written to its own file, which is the shape every
step had before. The same rows, row order, ties, dtypes, reasons, warnings and `contig_map` must come
back whatever the bound. (On hu416394's Manta file and bcftools gVCF the outputs of the code before
this change and after it hash the same, frame for frame: docs/PREVIOUS_ISSUES.md, F24.)

The fixture is a synthetic source assembly and a target built from it, a chain between them, and a
gVCF over 68 contigs, most of them decoys the store does not hold: calls, multi-allelic calls, indels
beside repeats, ties, reference blocks, breakends whose mates lie on any contig, contigs the chain does
not cover and contigs it covers at another size, and two source regions the chain sends to one place.
"""

import gzip
import random
import shutil
import zipfile
from pathlib import Path

import polars as pl
import pytest
from test_cache import aged, write_entry
from test_rerun import generations, manifests
from test_spill_order import DECOYS

from just_vcf import ReferenceStore, cache, liftover, load, normalize, resolve_ref
from just_vcf import loading as loading_module
from just_vcf.liftover import LIFT_LAYOUT, LIFT_VERSION
from just_vcf.normalize import NORMALIZE_LAYOUT, NORMALIZE_VERSION
from just_vcf.reference import RESOLVE_LAYOUT, RESOLVE_VERSION
from just_vcf.schema import SCHEMA_VERSION, canonical_order


def _sequences(seed: int) -> dict[str, str]:
    """The source assembly: random bases with a repeat every 200 bases, so indels have room to slide."""
    rng = random.Random(seed)
    sequences = {}
    for name, length in {"1": 3000, "2": 2000, "10": 1500, "X": 1500, "MT": 600, "GL000192.1": 800}.items():
        chars = [rng.choice("ACGT") for _ in range(length)]
        for start in range(100, length - 20, 200):
            chars[start : start + 8] = "TTTTTTTT" if start % 400 == 100 else "CACACACA"
        sequences[name] = "".join(chars)
    return sequences


def revcomp(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGT", "TGCA"))


SOURCE = _sequences(24)
UNSTORED = [*DECOYS, "chrUn_KI270302v1"]
"""Contigs the store does not hold: their rows stay unread and cannot lift."""
_DECOY_RNG = random.Random(2424)
DECOY_SEQ = {name: "".join(_DECOY_RNG.choice("ACGT") for _ in range(1000)) for name in UNSTORED}
CHAINED_DECOYS = DECOYS[:5]
"""Decoys the chain names, though the source assembly lacks them: the chain covers them at a size the
source does not have (`chain_contig_mismatch`), where the other decoys are `contig_not_in_chain`."""

# Target 1 is source 1 with source bases 1400..1420 deleted; target 2 is source 2 reverse-complemented.
TARGET = {
    "1": SOURCE["1"][:1400] + SOURCE["1"][1420:],
    "2": revcomp(SOURCE["2"]),
    "10": SOURCE["10"],
    "X": SOURCE["X"],
    "MT": SOURCE["MT"],
    "GL000192.1": SOURCE["GL000192.1"],
}


def chain_text() -> str:
    """Source 1 to target 1 with a gap; source 2 on the minus strand of target 2 from 400; source 2
    100..300 onto target 1 500..700 as well, which source 1 also lands on (overlapping lifted blocks,
    and a mismatching REF); MT as it is; GL000192.1 at a size the source does not have; five decoys to a
    contig the target lacks. 10 and X are not in the chain."""
    t1 = len(TARGET["1"])
    lines = [
        f"chain 1000 1 3000 + 0 3000 1 {t1} + 0 {t1} 1", "1400 20 0", "1580", "",
        "chain 900 2 2000 + 400 2000 2 2000 - 400 2000 2", "1600", "",
        f"chain 800 2 2000 + 100 300 1 {t1} + 500 700 3", "200", "",
        "chain 700 MT 600 + 0 600 MT 600 + 0 600 4", "600", "",
        "chain 600 GL000192.1 799 + 0 700 GL000192.1 800 + 0 700 5", "700", "",
    ]  # fmt: skip
    for index, name in enumerate(CHAINED_DECOYS):
        lines += [f"chain 500 {name} 1000 + 0 500 chrUn_q 600 + 0 500 {6 + index}", "500", ""]
    return "\n".join(lines)


def _records(seed: int) -> list[str]:
    rng = random.Random(seed)
    stored = list(SOURCE)
    rows: list[str] = []
    for _ in range(700):
        chrom = rng.choice(stored) if rng.random() < 0.7 else rng.choice(UNSTORED)
        seq = SOURCE.get(chrom) or DECOY_SEQ[chrom]
        pos = rng.randint(2, len(seq) - 10)
        kind = rng.random()
        if 0.4 <= kind < 0.58 and chrom in SOURCE and rng.random() < 0.5:
            # An indel at a repeat, where left-alignment moves it.
            pos = rng.choice(range(101, len(seq) - 20, 200)) + rng.randint(0, 7)
        ref = seq[pos - 1]
        others = [base for base in "ACGT" if base != ref]
        gq = rng.randint(1, 99)
        haploid = chrom == "MT" or (chrom in UNSTORED and rng.random() < 0.3)
        gt = "1" if haploid else rng.choice(["0/1", "1/1"])

        def line(ref: str, alt: str, gt: str = gt, info: str = ".", gq: int = gq) -> str:
            return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\tGT:GQ\t{gt}:{gq}"

        if kind < 0.2:
            end = min(len(seq), pos + rng.randint(0, 60))
            rows.append(f"{chrom}\t{pos}\t.\t{ref}\t<*>\t.\t.\tEND={end}\tGT:GQ\t{'0' if haploid else '0/0'}:{gq}")
        elif kind < 0.3:
            rows.append(line(ref, f"{others[0]},{others[1]}", "2" if haploid else "1/2"))
        elif kind < 0.4:
            # The same site twice, told apart only by GQ: a tie in every sort key.
            rows += [line(ref, others[2], gq=value) for value in rng.sample(range(1, 99), 2)]
        elif kind < 0.5:
            rows.append(line(seq[pos - 1 : pos + rng.randint(1, 3)], ref))
        elif kind < 0.58:
            rows.append(line(ref, ref + seq[pos : pos + 2]))
        elif kind < 0.64:
            mate = rng.choice(stored + UNSTORED[:3] + CHAINED_DECOYS[:1])
            mate_pos = rng.randint(1, len(SOURCE.get(mate) or DECOY_SEQ[mate]))
            rows.append(line(ref, f"{ref}[{mate}:{mate_pos}[", info="SVTYPE=BND"))
        elif kind < 0.67:
            rows.append(line(others[0], others[1]))
        else:
            rows.append(line(ref, others[0]))
    rng.shuffle(rows)
    return rows


HEADER = [
    "##fileformat=VCFv4.3",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SOURCE.items()],
    *[f"##contig=<ID={name},length={len(DECOY_SEQ[name])}>" for name in UNSTORED],
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]


def _fasta(path: Path, sequences: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    return path


@pytest.fixture(scope="module")
def inputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Path]:
    """(the zipped gVCF, the store directory, the chain)."""
    base = tmp_path_factory.mktemp("f24")
    store = ReferenceStore.open(base / "store")
    store.add_fasta(_fasta(base / "source.fa", SOURCE), assembly_label="GRCh38")
    store.add_fasta(_fasta(base / "target.fa", TARGET), assembly_label="GRCh37")
    chain = base / "source_to_target.chain.gz"
    chain.write_bytes(gzip.compress(chain_text().encode()))
    vcf = base / "decoys.zip"
    with zipfile.ZipFile(vcf, "w") as archive:
        archive.writestr("decoys.g.vcf", "\n".join([*HEADER, *_records(24)]) + "\n")
    return vcf, base / "store", chain


def run(inputs: tuple[Path, Path, Path], upto: str = "all") -> dict[str, object]:
    """Every step over the fixture (or up to `upto`), every output collected."""
    vcf, store_dir, chain = inputs
    store = ReferenceStore.open(store_dir)
    resolved = resolve_ref(load(vcf), store)
    out: dict[str, object] = {
        "resolved": resolved.calls.collect(),
        "resolved_warnings": resolved.warnings,
        "resolved_facets": resolved.facets.collect(),
    }
    if upto in ("all", "normalize"):
        normalized = normalize(resolved, store)
        out |= {"normalized": normalized.calls.collect(), "normalized_warnings": normalized.warnings}
    if upto in ("all", "lift"):
        lifted = liftover(resolved, chain, store)
        out |= {
            "lifted": lifted.loaded.calls.collect(),
            "lifted_blocks": lifted.loaded.blocks.collect(),
            "report": lifted.report.collect(),
            "reasons": lifted.reasons,
            "block_reasons": lifted.block_reasons,
            "lifted_warnings": lifted.loaded.warnings,
            "lifted_facets": lifted.loaded.facets.collect(),
        }
    if upto == "all":
        again = normalize(lifted.loaded, store)
        out |= {"normalized_lift": again.calls.collect(), "normalized_lift_warnings": again.warnings}
    return out


def same(got: dict[str, object], expected: dict[str, object]) -> None:
    assert got.keys() == expected.keys()
    for name, value in expected.items():
        if isinstance(value, pl.DataFrame):
            assert got[name].schema == value.schema, name
            assert got[name].equals(value), name
        else:
            assert got[name] == value, name


def per_contig(inputs: tuple[Path, Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    """The oracle: every step with every contig a group of its own, in a cache of its own."""
    with monkeypatch.context() as patch:
        patch.setenv("JUST_VCF_CACHE_DIR", str(tmp_path / "per-contig-cache"))
        patch.setattr(loading_module, "SPILL_PART_ROWS", 1)
        return run(inputs)


@pytest.mark.parametrize("part_rows", [None, 7, 40])
def test_every_output_is_the_per_contig_runs_row_for_row(
    part_rows: int | None, inputs: tuple[Path, Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Rows, row order, ties, dtypes, reasons, warnings and the facets (`contig_map`), whatever the part
    bound: one group, groups of a few contigs, groups of a few dozen rows (None keeps the default)."""
    expected = per_contig(inputs, tmp_path, monkeypatch)
    if part_rows is not None:
        monkeypatch.setattr(loading_module, "SPILL_PART_ROWS", part_rows)
    same(run(inputs), expected)


def test_the_fixture_exercises_what_it_claims(inputs: tuple[Path, Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Decoys left unread, contigs outside the chain and at another size in it, lifted breakends and
    mates that are not, ties, overlapping lifted blocks, conflicts at one position on several contigs,
    and indels left-aligned: the cases the groups must not change."""
    out = per_contig(inputs, tmp_path, monkeypatch)
    resolved, report = out["resolved"], out["report"]
    assert len(set(resolved["chrom"])) > 50
    assert resolved.filter(pl.col("chrom").is_in(UNSTORED))["ref"].null_count() == resolved.filter(pl.col("chrom").is_in(UNSTORED)).height > 0
    assert {"contig_not_in_chain", "chain_contig_mismatch", "ref_unread", "no_chain_block", "breakend_mate_not_lifted"} <= set(out["reasons"])
    assert report.filter(pl.col("lift_reason").is_null() & pl.col("alts").list.join(",").str.contains(r"\[")).height > 0
    keys = ["chrom", "pos", "end_pos", "source_ref", pl.col("alts").list.join(",")]
    assert resolved.select(keys).is_duplicated().any()
    assert any("lifted block(s) overlap another" in line for line in out["lifted_warnings"])
    assert any("left-aligned to another position" in line for line in out["normalized_warnings"])
    assert any("state the same allele and differ otherwise" in line for line in out["normalized_warnings"])
    # The same position, REF and ALT on two contigs: a conflict counted per contig must not join them.
    calls = out["normalized"]
    assert calls.group_by("pos", "ref", "alts").agg(pl.col("chrom").n_unique().alias("contigs"))["contigs"].max() > 1


def _parts(kind: str) -> dict[str, list[Path]]:
    (entry,) = [cache.cache_dir() / name for name in manifests(kind)]
    (generation,) = generations(entry)
    record = cache.Generation.model_validate_json((generation / "result.json").read_bytes())
    return {group: [generation / name for name in names] for group, names in record.files.items()}


def test_each_step_writes_a_few_parts_its_frames_read_in_one_scan(inputs: tuple[Path, Path, Path]) -> None:
    """At the default bound the 68 contigs share one part per frame, and each frame is one scan of its
    parts with no sort. Before, each step wrote a file per contig with rows."""
    vcf, store_dir, chain = inputs
    store = ReferenceStore.open(store_dir)
    resolved = resolve_ref(load(vcf), store)
    normalized = normalize(resolved, store)
    lifted = liftover(resolved, chain, store)
    for frame in (resolved.calls, normalized.calls, lifted.loaded.calls, lifted.loaded.blocks, lifted.report):
        plan = frame.explain()
        assert "SORT" not in plan
        assert plan.count("SCAN") == 1
    assert {group: len(paths) for group, paths in _parts(cache.RESOLVE_SPILL).items()} == {"calls": 1}
    assert {group: len(paths) for group, paths in _parts(cache.NORMALIZE_SPILL).items()} == {"calls": 1}
    lift = {group: len(paths) for group, paths in _parts(cache.LIFT_SPILL).items()}
    assert lift == {"calls": 1, "report": 1, "blocks": 1, "block_report": 1, "mate_contigs": 1}


def test_parts_hold_whole_contigs_up_to_the_bound(inputs: tuple[Path, Path, Path], monkeypatch: pytest.MonkeyPatch) -> None:
    """At a bound of 7 rows: no contig in two parts of a frame, and a part over the bound holds one
    contig. The rows `normalize` writes may outgrow their group's (a split row is several), so only
    whole contigs are asked of it."""
    monkeypatch.setattr(loading_module, "SPILL_PART_ROWS", 7)
    vcf, store_dir, chain = inputs
    store = ReferenceStore.open(store_dir)
    resolved = resolve_ref(load(vcf), store)
    normalize(resolved, store)
    liftover(resolved, chain, store)
    lift = _parts(cache.LIFT_SPILL)
    for name, paths, bounded in (
        ("resolve", _parts(cache.RESOLVE_SPILL)["calls"], True),
        ("normalize", _parts(cache.NORMALIZE_SPILL)["calls"], False),
        ("report", lift["report"], True),
        ("lifted", lift["calls"], True),
        ("lifted blocks", lift["blocks"], True),
    ):
        assert len(paths) > 1, name
        seen: set[str] = set()
        for path in paths:
            rows = pl.read_parquet(path)
            chroms = set(rows["chrom"])
            assert not chroms & seen, name
            seen |= chroms
            assert not bounded or rows.height <= 7 or len(chroms) == 1, name
            if name != "report":
                assert rows.equals(canonical_order(rows.lazy()).collect()), name


STEPS = {
    # name: (the kind it spills, its version key, its layout key, the layout key's value)
    "resolve": (cache.RESOLVE_SPILL, "RESOLVE_VERSION", "RESOLVE_LAYOUT", RESOLVE_LAYOUT),
    "normalize": (cache.NORMALIZE_SPILL, "NORMALIZE_VERSION", "NORMALIZE_LAYOUT", NORMALIZE_LAYOUT),
    "lift": (cache.LIFT_SPILL, "LIFT_VERSION", "LIFT_LAYOUT", LIFT_LAYOUT),
}
OUTPUTS = {"resolve": ("resolved",), "normalize": ("normalized",), "lift": ("lifted", "lifted_blocks", "report")}


@pytest.fixture
def keyed(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[tuple[dict, list[Path], dict]]:
    """Every generation key a step computes: its versions, copies of the staged files, the facet."""
    seen: list[tuple[dict, list[Path], dict]] = []
    real = cache.content_key

    def spying(files: list[Path], facet: dict) -> str:
        generation = facet.get("_generation")
        if generation is not None:
            copies = []
            for path in files:
                copy = tmp_path / "staged" / f"{len(seen):03d}-{path.name}"
                copy.parent.mkdir(exist_ok=True)
                shutil.copyfile(path, copy)
                copies.append(copy)
            seen.append((generation["versions"], copies, facet))
        return real(files, facet)

    monkeypatch.setattr(cache, "content_key", spying)
    return seen


@pytest.mark.parametrize("step", STEPS)
def test_a_spill_in_the_per_contig_layout_is_not_read(
    step: str, keyed: list[tuple[dict, list[Path], dict]], inputs: tuple[Path, Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A generation an older version committed, keyed as that version keyed it (its versions without
    the layout key), sits in the same entry, used within the grace period so the sweep keeps it. Its
    files are cut to one row each, so reading it would show. A run of the same input writes a generation
    of its own beside it and never reads the old one."""
    kind, version_key, layout_key, _ = STEPS[step]
    expected = run(inputs, upto=step)
    ((versions, copies, facet),) = [item for item in keyed if version_key in item[0]]
    (entry,) = [cache.cache_dir() / name for name in manifests(kind)]
    (written,) = generations(entry)
    record = cache.Generation.model_validate_json((written / "result.json").read_bytes())

    monkeypatch.setenv("JUST_VCF_CACHE_DIR", str(tmp_path / "older-cache"))
    old_versions = {name: value for name, value in versions.items() if name != layout_key}
    older = write_entry("frames", entry.name, kind, old_versions)
    (older / "calls.parquet").unlink()
    old = older / cache.content_key(copies, {**facet, "_generation": {**facet["_generation"], "versions": old_versions}})
    old.mkdir()
    for names in record.files.values():
        for name in names:
            pl.read_parquet(written / name).head(1).write_parquet(old / name)
    (old / "result.json").write_text(record.model_dump_json())
    keyed.clear()

    got = run(inputs, upto=step)

    ((_, again, _),) = [item for item in keyed if version_key in item[0]]
    # The input staged the same bytes, so the old key is the one this run would compute without the layout.
    assert [path.read_bytes() for path in again] == [path.read_bytes() for path in copies]
    for name in OUTPUTS[step]:
        assert got[name].equals(expected[name]), name
    assert old in generations(older) and len(generations(older)) == 2


@pytest.mark.parametrize("step", STEPS)
def test_a_step_sweeps_a_spill_in_the_per_contig_layout(step: str, inputs: tuple[Path, Path, Path]) -> None:
    """A spill whose manifest lacks the step's layout key, unused past the grace period, is removed; one
    a newer layout wrote stays."""
    kind, version_key, layout_key, layout = STEPS[step]
    version = {"resolve": RESOLVE_VERSION, "normalize": NORMALIZE_VERSION, "lift": LIFT_VERSION}[step]
    current = {"SCHEMA_VERSION": SCHEMA_VERSION, version_key: version, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    suffix = f"-ref-SQ_ab-c-{RESOLVE_VERSION}" if step == "resolve" else ""
    before = aged(write_entry("frames", f"{'c' * 32}-{SCHEMA_VERSION}{suffix}", kind, current))
    newer = aged(write_entry("frames", f"{'d' * 32}-{SCHEMA_VERSION}{suffix}", kind, {**current, layout_key: layout + 1}))
    run(inputs, upto="lift" if step == "lift" else "normalize")
    assert not before.exists() and newer.exists()
