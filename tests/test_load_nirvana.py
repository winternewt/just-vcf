"""Illumina Annotation Engine (Nirvana) JSON (`nirvana_json`), read through `load`.

The synthetic positions copy the shapes counted in the one PGP-HMS file (2026-10-23, counts only):
4,776,534 lines, 4,763,188 position objects, then an `omim` section of 13,343 lines; one sample;
genotypes `0/1`, `1/1`, `1/2`, haploid `1` and `0` on chrM, `0/0` reference records without
`altAlleles`, and 228 positions in the first 200,000 whose sample has no genotype.
"""

import json
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load
from just_vcf.formats import get_spec
from just_vcf.loaders.vcf import MultiSampleError
from just_vcf.schema import CALLS, canonical_chrom
from just_vcf.sniff import Sniffed, sniff
from test_load_vcf import assert_gt_consistent

HEADER = {"annotator": "Illumina Annotation Engine 1.5.1.0", "creationTime": "2016-12-02 12:36:47", "schemaVersion": 4,
          "dataVersion": "84.22.36", "samples": ["S1"]}


def _position(chrom: str, pos: int, ref: str, alts: list[str] | None, samples: list[dict], filters: list[str] | None = None) -> dict:
    found = {"chromosome": chrom, "refAllele": ref, "position": pos, "quality": 40.0, "filters": filters or ["PASS"], "samples": samples,
             "cytogeneticBand": "1p36.33", "variants": [{"vid": f"{chrom}:{pos}"}]}
    if alts is not None:
        found["altAlleles"] = alts
    return found


POSITIONS = [
    _position("chr1", 752566, "G", ["A"], [{"variantFreq": 0.5, "totalDepth": 30, "genotypeQuality": 40, "alleleDepths": [15, 15], "genotype": "0/1"}]),
    _position("chr1", 752721, "A", ["G", "T"], [{"totalDepth": 28, "genotypeQuality": 33, "alleleDepths": [0, 14, 14], "genotype": "1|2"}],
              filters=["LowGQX"]),
    _position("chr2", 1000, "ACG", ["A"], [{"totalDepth": 9, "genotypeQuality": 12, "alleleDepths": [5], "genotype": "1/1"}]),
    _position("chr2", 2000, "C", None, [{"totalDepth": 20, "genotypeQuality": 60, "genotype": "0/0"}]),
    _position("chrM", 150, "C", ["T"], [{"totalDepth": 900, "genotypeQuality": 99, "alleleDepths": [0, 900], "genotype": "1"}]),
    _position("chr3", 500, "T", ["C"], [{"totalDepth": 3}]),
]


def _write(path: Path, positions: list[dict], header: dict = HEADER, tail: str = '],"omim":[\n{"mimNumber":"100650","description":"x"}\n]}\n') -> Path:
    head = json.dumps({"header": header}, separators=(",", ":"))[:-1] + ',"positions":[\n'
    body = ",\n".join(json.dumps(p, separators=(",", ":")) for p in positions) + "\n"
    path.write_text(head + body + tail)
    return path


def test_positions_are_vcf_rows(tmp_path: Path) -> None:
    loaded = load(_write(tmp_path / "s1.json", POSITIONS))
    calls = loaded.calls.collect()
    assert calls.schema == CALLS
    assert_gt_consistent(calls)
    assert calls.height == len(POSITIONS), "the omim section holds no calls"
    rows = {(row["chrom"], row["pos"]): row for row in calls.iter_rows(named=True)}
    for position in POSITIONS:
        row = rows[(canonical_chrom(position["chromosome"]), position["position"])]
        sample = position["samples"][0]
        alts = position.get("altAlleles", [])
        assert row["source_ref"] == position["refAllele"] and row["ref"] is None and row["ref_check"] == "not_run"
        assert row["alts"] == alts
        assert row["end_pos"] == position["position"] + len(position["refAllele"]) - 1
        assert row["source_gt"] == sample.get("genotype")
        assert row["gq"] == sample.get("genotypeQuality") and row["dp"] == sample.get("totalDepth")
        depths = sample.get("alleleDepths")
        assert row["ad"] == (depths if depths is not None and len(depths) == len(alts) + 1 else None), "ad only with one value per allele"
        assert row["filters"] == position["filters"] and row["qual"] == position["quality"]
    assert rows[("1", 752721)]["phased"] == [True, True]
    assert rows[("MT", 150)]["gt"] == [1] and rows[("MT", 150)]["gt_alleles"] == ["T"]
    assert rows[("2", 2000)]["gt"] == [0, 0] and rows[("2", 2000)]["alts"] == []
    assert rows[("3", 500)]["gt"] is None
    assert any("1 position(s) whose sample has no genotype" in w for w in loaded.warnings)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["sample"] == "S1" and facet["assembly_label"] is None, "schema v4 states no assembly"
    assert "cytogeneticBand" in facet["dropped_fields"] and "variants" in facet["dropped_fields"]
    assert "samples.variantFreq" in facet["dropped_fields"]
    assert facet["contig_map"] == [{"source": "chr1", "chrom": "1"}, {"source": "chr2", "chrom": "2"},
                                   {"source": "chr3", "chrom": "3"}, {"source": "chrM", "chrom": "MT"}]


def test_the_example_loads_and_a_rerun_is_identical(tmp_path: Path) -> None:
    path = tmp_path / "example.json"
    path.write_text(get_spec("nirvana_json").example)
    first = load(path).calls.collect()
    assert first.height == 2
    assert load(path).calls.collect().equals(first)


def test_two_samples_need_a_choice(tmp_path: Path) -> None:
    two = [
        _position("chr1", 100, "A", ["G"], [{"genotype": "0/1", "totalDepth": 5}, {"genotype": "1/1", "totalDepth": 7}]),
    ]
    path = _write(tmp_path / "pair.json", two, header={**HEADER, "samples": ["S1", "S2"]})
    with pytest.raises(MultiSampleError):
        load(path)
    second = load(path, sample="S2")
    row = second.calls.collect().row(0, named=True)
    assert row["gt"] == [1, 1] and row["dp"] == 7
    assert second.facets.collect()["sample"].to_list() == ["S2"]
    first = load(path, sample="S1").facets.collect()["facet_id"].to_list()
    assert first != second.facets.collect()["facet_id"].to_list(), "each sample is its own facet"


def test_positions_without_a_genotype_load_with_the_genotype_unknown(tmp_path: Path) -> None:
    """A file whose every sample entry has no genotype loads with `gt` and `gt_alleles` null on every
    row. polars' `list.gather` over a batch whose every index list is null makes a list of a Null inner
    type, which the spill's parquet sink cannot write: the load once failed there."""
    positions = [
        _position("chr1", 100 + i, "A", ["G", "T"] if i % 2 else ["G"], [{"totalDepth": 10 + i, "genotypeQuality": 30}]) for i in range(4)
    ]
    loaded = load(_write(tmp_path / "s1.json", positions))
    calls = loaded.calls.collect()
    assert calls.schema == CALLS
    assert calls["pos"].to_list() == [p["position"] for p in positions]
    assert calls["dp"].to_list() == [p["samples"][0]["totalDepth"] for p in positions]
    assert calls["gt"].null_count() == calls["gt_alleles"].null_count() == len(positions)
    assert any(f"{len(positions)} position(s) whose sample has no genotype" in w for w in loaded.warnings)


UNREADABLE = "not readable as allele indices kept with gt, gt_alleles and phased null, and the text in source_gt"


def test_an_unreadable_genotype_is_withheld_whole(tmp_path: Path) -> None:
    """A genotype that is not allele indices keeps its text in `source_gt`, and `gt`, `gt_alleles` and
    `phased` are null: the genotype is unknown, so none of it is read, not the `0` of `0/x` either. A
    `[0, null]` would be the frame of `0/.`, a half-call the file did not write (SCHEMA.md Decision 2).
    The warning says what the row holds."""
    genotypes = ["N/A", "0/x", "x", "-", "1:2", "0/1"]
    positions = [_position("chr1", 100 + i, "A", ["G"], [{"totalDepth": 10, "genotype": gt}]) for i, gt in enumerate(genotypes)]
    loaded = load(_write(tmp_path / "s1.json", positions))
    calls = loaded.calls.collect().sort("pos")
    assert calls["source_gt"].to_list() == genotypes
    withheld = [None] * (len(genotypes) - 1)
    assert calls["gt"].to_list() == [*withheld, [0, 1]]
    assert calls["gt_alleles"].to_list() == [*withheld, ["A", "G"]]
    assert calls["phased"].to_list() == [*withheld, [False, False]]
    assert f"{len(withheld)} genotype(s) {UNREADABLE}" in loaded.warnings, loaded.warnings


def test_a_sample_of_only_unreadable_genotypes_loads(tmp_path: Path) -> None:
    """Every genotype unreadable: `gt` is null on every row, and the guarded gather keeps the spill's
    `gt_alleles` a typed list column."""
    genotypes = ["N/A", "x", "0/x"]
    positions = [_position("chr1", 100 + i, "A", ["G"], [{"totalDepth": 10, "genotype": gt}]) for i, gt in enumerate(genotypes)]
    calls = load(_write(tmp_path / "s1.json", positions)).calls.collect()
    assert calls.schema == CALLS
    assert calls["gt"].to_list() == calls["gt_alleles"].to_list() == [None] * len(genotypes)


def test_an_index_past_the_alts_is_counted_from_the_text(tmp_path: Path) -> None:
    """An index past the ALT list is counted from the genotype text's digits, as the key=value table
    counts it: `gt` is a UInt8 list, so `1/300` reads `[1, null]`, and the count taken from `gt` missed
    it. Each genotype is counted under one reason: `1:2` is unreadable, never also past the ALTs."""
    genotypes = {"1/300": "past", "0/2": "past", "0/" + "9" * 20: "past", "1:2": "unreadable", "0/1": None}
    positions = [_position("chr1", 100 + i, "A", ["G"], [{"totalDepth": 10, "genotype": gt}]) for i, gt in enumerate(genotypes)]
    loaded = load(_write(tmp_path / "s1.json", positions))
    calls = loaded.calls.collect().sort("pos")
    assert calls.filter(pl.col("source_gt") == "1/300")["gt"].to_list() == [[1, None]]
    counts = {reason: sum(1 for why in genotypes.values() if why == reason) for reason in ("unreadable", "past")}
    found = [w for w in loaded.warnings if "past the row's ALT list" in w or "not readable as allele indices" in w]
    assert sorted(found) == sorted([
        f"{counts['past']} genotype(s) index an allele past the row's ALT list; gt_alleles has a null there",
        f"{counts['unreadable']} genotype(s) {UNREADABLE}",
    ]), found


def test_real_nirvana_json(pgp_participants: Path) -> None:
    found: Sniffed | None = None
    for path in sorted(pgp_participants.glob("*/upload-*/*")):
        if path.is_file() and path.stat().st_size > 100_000_000 and (result := sniff(path)).format == "nirvana_json":
            found = result
            break
    if found is None:
        pytest.skip("no nirvana_json upload found")
    loaded = load(found)
    calls = loaded.calls.collect()
    assert calls.select(pl.col("source_ref").null_count()).item() == 0
    assert_gt_consistent(calls.head(200_000))
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["format_id"] == "nirvana_json" and facet["variant_classes"]
