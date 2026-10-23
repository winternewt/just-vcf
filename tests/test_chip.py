"""Chip detection for array facets (RM4): stated by a header, or inferred from a 23andMe probe set."""

import hashlib
import importlib.util
import json
from importlib import resources
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load
from just_vcf.chip import (
    CHIPS,
    LOW_ID_FLOOR,
    MARKER_MODULUS,
    MIN_COVERAGE,
    MIN_MARGIN,
    MIN_RECALL,
    MIN_SAMPLED_MARKERS,
    OVERLAP_FORMATS,
    ROW_BANDS_23ANDME,
    ChipOutcome,
    chip_markers,
    detect_chip,
    id_profile,
    inferred_chip,
    is_chip,
    marker_overlap,
    overlap_chip,
    stated_chip,
)
from just_vcf.formats import get_spec


def _facet(path: Path) -> dict:
    return load(path).facets.collect().row(0, named=True)


def _write(tmp_path: Path, text: str, name: str = "genome.txt") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.mark.parametrize(("version", "chip"), [("V1.0", "ancestrydna_v1"), ("V2.0", "ancestrydna_v2")])
def test_ancestrydna_states_its_array_version(tmp_path: Path, version: str, chip: str) -> None:
    example = get_spec("ancestrydna").example.replace("array version: V1.0", f"array version: {version}")
    loaded = load(_write(tmp_path, example))
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["chip"], facet["chip_basis"]) == (chip, "stated")
    assert not any(w.startswith("chip ") for w in loaded.warnings)


def test_a_stated_version_with_no_entry_leaves_the_chip_unknown_and_quotes_it(tmp_path: Path) -> None:
    example = get_spec("ancestrydna").example.replace("array version: V1.0", "array version: V3.0")
    loaded = load(_write(tmp_path, example))
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["chip"], facet["chip_basis"]) == (None, None)
    assert any(w.startswith("chip unknown: the header's AncestryDNA array version 'V3.0' names no known chip") for w in loaded.warnings)


def test_a_final_report_names_its_manifest(tmp_path: Path) -> None:
    facet = _facet(_write(tmp_path, get_spec("illumina_final_report").example))
    assert (facet["chip"], facet["chip_basis"]) == ("illumina:GSA-24v2-0_A1", "stated")


def test_a_final_report_content_line_after_data_is_not_a_statement() -> None:
    assert stated_chip("illumina_final_report", ["[Header]", "[Data]", "Content\tGSA-24v2-0_A1.bpm"]) is None


def test_living_dna_states_its_chip(tmp_path: Path) -> None:
    facet = _facet(_write(tmp_path, get_spec("living_dna").example))
    assert (facet["chip"], facet["chip_basis"]) == ("living_dna_sirius", "stated")


def test_a_statement_counts_only_for_its_own_format() -> None:
    line = "#Data was collected using AncestryDNA array version: V2.0"
    assert stated_chip("23andme", [line]) is None
    assert stated_chip("ancestrydna", [line]).chip == "ancestrydna_v2"


def test_the_vocabulary_validator() -> None:
    assert all(is_chip(chip) for chip in CHIPS)
    assert is_chip("illumina:HumanOmni2.5-8v1_A")
    assert not is_chip("illumina:")
    assert not is_chip("illumina:GSA 24")
    assert not is_chip("23andme_v6")


# 23andMe: the probe set. Synthetic frames shaped like each generation's export.


def _rows(total: int, internal: int, low: int) -> pl.DataFrame:
    """`total` rows whose first ids are `internal` internal ids (`low` of them numbered below
    `LOW_ID_FLOOR`) and rsids for the rest."""
    ids = [f"i{700_001 + n}" for n in range(low)]
    ids += [f"i{LOW_ID_FLOOR + 1 + n}" for n in range(internal - low)]
    ids += [f"rs{n + 1}" for n in range(total - internal)]
    return pl.DataFrame({"ids": [[i] for i in ids]})


GENERATIONS = {
    # chip: (rows, internal ids, of them below the floor), each inside the bands measured on PGP-HMS.
    "23andme_v1_v2": (576_516, 2_521, 0),
    "23andme_v3": (960_613, 10_709, 0),
    "23andme_v4": (602_352, 49_289, 0),
    "23andme_v5": (638_468, 17_829, 17_000),
}


@pytest.mark.parametrize("chip", sorted(GENERATIONS))
def test_a_23andme_generation_is_named_when_rows_and_ids_agree(chip: str) -> None:
    rows, internal, low = GENERATIONS[chip]
    found = inferred_chip("23andme", _rows(rows, internal, low))
    assert (found.chip, found.outcome, found.basis.value) == (chip, ChipOutcome.INFERRED, "inferred")
    assert f"{rows} rows, {internal} internal ids" in found.reason


def test_the_id_profile_is_counted_from_the_rows() -> None:
    profile = id_profile(_rows(1_000, 40, 30))
    assert profile == {"rows": 1_000, "internal": 40, "internal_share": 0.04, "low_share": 0.75}


def test_rows_of_one_generation_with_another_generations_ids_stay_unknown() -> None:
    """A v5-sized export whose internal ids look like v4's: the two features disagree."""
    found = inferred_chip("23andme", _rows(638_468, 52_000, 0))
    assert not found and found.outcome is ChipOutcome.DISAGREE
    assert "the row count points at 23andme_v5, the internal ids at 23andme_v4" in found.reason


def test_a_count_no_generation_has_stays_unknown() -> None:
    """Shaped like PGP-HMS's 1.4-million-row 23andMe exports, whose ids look like v5's: no band holds them."""
    found = inferred_chip("23andme", _rows(1_415_194, 27_582, 25_000))
    assert not found and found.outcome is ChipOutcome.NO_BAND


def test_the_bands_leave_no_overlap() -> None:
    bands = sorted(ROW_BANDS_23ANDME.values())
    assert all(high < next_low for (_, high), (next_low, _) in zip(bands, bands[1:]))


def test_other_formats_are_not_inferred() -> None:
    rows, internal, low = GENERATIONS["23andme_v5"]
    found = inferred_chip("decodeme", _rows(rows, internal, low))
    assert not found and found.outcome is ChipOutcome.NO_RULE


def test_a_resaved_23andme_table_is_read_by_its_internal_ids() -> None:
    """A preamble-less rsid table with 23andMe's `i…` ids is a re-saved 23andMe export; without them it
    could be any vendor's (an AncestryDNA export re-saved has none)."""
    rows, internal, low = GENERATIONS["23andme_v3"]
    assert inferred_chip("rsid_table_headerless", _rows(rows, internal, low)).chip == "23andme_v3"
    assert inferred_chip("rsid_table", _rows(701_478, 0, 0)).outcome is ChipOutcome.NO_RULE


def test_a_v1_v2_export_near_the_id_boundary_is_still_v1_v2() -> None:
    """The highest v1/v2 internal-id share on PGP-HMS (3,342 of 594,930 rows, 0.56%)."""
    assert inferred_chip("23andme", _rows(594_930, 3_342, 0)).chip == "23andme_v1_v2"


def test_a_loaded_23andme_export_gets_its_inferred_chip(tmp_path: Path) -> None:
    rows, internal, low = GENERATIONS["23andme_v5"]
    frame = _rows(rows, internal, low).with_row_index("n")
    body = "".join(f"{row[1][0]}\t1\t{row[0] + 1}\tAA\n" for row in frame.iter_rows())
    preamble = "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018\n# rsid\tchromosome\tposition\tgenotype\n"
    loaded = load(_write(tmp_path, preamble + body))
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["chip"], facet["chip_basis"]) == ("23andme_v5", "inferred")
    assert any(w.startswith("chip inferred from the rows as 23andme_v5: 638468 rows") for w in loaded.warnings)


def test_a_small_23andme_export_says_its_chip_is_unknown(tmp_path: Path) -> None:
    loaded = load(_write(tmp_path, get_spec("23andme").example))
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["chip"], facet["chip_basis"]) == (None, None)
    assert [w for w in loaded.warnings if w.startswith("chip ")] == [
        "chip unknown: 4 rows, 1 internal ids (25.00% of rows, 100% of them below i1000000): the row count is in no generation's band"
    ]


def test_a_vcf_has_no_chip(tmp_path: Path) -> None:
    text = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    facet = _facet(_write(tmp_path, text, "calls.vcf"))
    assert (facet["chip"], facet["chip_basis"], facet["strand_check"]) == (None, None, None)


# Marker overlap with the Illumina manifests (FTDNA, MyHeritage, Final Reports that state no chip).
# Exports are built from the packaged marker sample itself, so the expected shares follow from it.


def _markers(chip: str) -> list[int]:
    return chip_markers()[chip].to_list()


def _only(chip: str) -> list[int]:
    """Sampled rsids in this chip's manifest and in no other's."""
    others = set().union(*(set(_markers(c)) for c in chip_markers() if c != chip))
    return [n for n in _markers(chip) if n not in others]


def _off_manifest(count: int) -> list[int]:
    """Sampled rsid numbers that no manifest holds: multiples of the modulus missing from every sample."""
    held = set().union(*(set(_markers(c)) for c in chip_markers()))
    found, n = [], 0
    while len(found) < count:
        n += MARKER_MODULUS
        if n not in held:
            found.append(n)
    return found


def _export(sampled: list[int], filler: int = 0, prefix: str = "rs") -> pl.DataFrame:
    """An export's `ids`: the given sampled rsids, plus `filler` rsids outside the sample (never compared)."""
    unsampled = [n * MARKER_MODULUS + 1 for n in range(filler)]
    return pl.DataFrame({"ids": [[f"{prefix}{n}"] for n in sampled + unsampled]})


def test_the_packaged_sample_holds_only_sampled_rsids() -> None:
    markers = chip_markers()
    assert set(markers) == {"illumina_gsa", "illumina_omniexpress"} and set(markers) <= CHIPS
    for sample in markers.values():
        assert sample.len() > 5 * MIN_SAMPLED_MARKERS and (sample % MARKER_MODULUS == 0).all()
        assert sample.is_sorted() and sample.n_unique() == sample.len()


@pytest.mark.parametrize("chip", ["illumina_gsa", "illumina_omniexpress"])
def test_an_export_of_one_manifests_markers_is_that_family(chip: str) -> None:
    export = _export(_markers(chip), filler=20_000)
    found = overlap_chip(export)
    assert (found.chip, found.basis, found.outcome) == (chip, "inferred", ChipOutcome.INFERRED)
    coverage, recall = marker_overlap(export)[chip]
    assert (coverage, recall) == (1.0, 1.0), "filler outside the sample is never compared"
    for fmt in sorted(OVERLAP_FORMATS):
        assert detect_chip(fmt, [], export).chip == chip


def test_gsa_names_in_a_final_report_count_as_their_rsids() -> None:
    assert overlap_chip(_export(_markers("illumina_gsa"), prefix="GSA-rs")).chip == "illumina_gsa"


def test_a_final_report_that_states_its_manifest_is_not_inferred() -> None:
    header = ["[Header]", "Content\tGSA-24v2-0_A1.bpm", "[Data]"]
    found = detect_chip("illumina_final_report", header, _export(_markers("illumina_omniexpress")))
    assert (found.chip, found.basis) == ("illumina:GSA-24v2-0_A1", "stated")


def test_the_coverage_floor() -> None:
    """Coverage is the winning manifest's share of the export's sampled rsids; at the floor it names a chip."""
    held = _only("illumina_omniexpress")[:9_000]
    outside = round(len(held) * (1 - MIN_COVERAGE) / MIN_COVERAGE)
    at_floor = overlap_chip(_export(held + _off_manifest(outside)))
    assert marker_overlap(_export(held + _off_manifest(outside)))["illumina_omniexpress"][0] >= MIN_COVERAGE
    assert at_floor.chip == "illumina_omniexpress"
    below = overlap_chip(_export(held + _off_manifest(outside + 50)))
    assert (below.chip, below.outcome) == (None, ChipOutcome.NO_CLEAR_OVERLAP)


def test_two_families_mixed_name_neither() -> None:
    """Half of each family's own markers: each manifest holds about half, and neither leads by the margin."""
    gsa, omni = _only("illumina_gsa"), _only("illumina_omniexpress")
    export = _export(gsa[: len(gsa) // 2] + omni[: len(gsa) // 2])
    shares = marker_overlap(export)
    assert abs(shares["illumina_gsa"][0] - shares["illumina_omniexpress"][0]) < MIN_MARGIN
    assert overlap_chip(export).outcome is ChipOutcome.NO_CLEAR_OVERLAP


def test_an_extract_of_a_chip_is_not_the_chip() -> None:
    markers = _markers("illumina_omniexpress")
    extract = markers[: int(len(markers) * (MIN_RECALL - 0.05))]
    assert len(extract) >= MIN_SAMPLED_MARKERS
    coverage, recall = marker_overlap(_export(extract))["illumina_omniexpress"]
    assert coverage == 1.0 and recall < MIN_RECALL
    assert overlap_chip(_export(extract)).outcome is ChipOutcome.NO_CLEAR_OVERLAP


def test_too_few_sampled_rsids_are_not_compared() -> None:
    found = overlap_chip(_export(_markers("illumina_omniexpress")[: MIN_SAMPLED_MARKERS - 1], filler=50_000))
    assert (found.chip, found.outcome) == (None, ChipOutcome.TOO_FEW_MARKERS)
    assert overlap_chip(_export(_markers("illumina_omniexpress")[:MIN_SAMPLED_MARKERS])).outcome is not ChipOutcome.TOO_FEW_MARKERS


def test_23andme_keeps_its_own_rule_over_the_overlap() -> None:
    """A 23andMe v5 export is a custom GSA; its chip is the vendor's generation, never the overlap's family."""
    assert "23andme" not in OVERLAP_FORMATS
    found = inferred_chip("23andme", _export(_markers("illumina_gsa")))
    assert found.chip != "illumina_gsa"


def test_a_loaded_ftdna_export_gets_its_family_by_overlap(tmp_path: Path) -> None:
    body = "".join(f'"rs{n}","1","{1000 + i}","AG"\n' for i, n in enumerate(_markers("illumina_omniexpress")))
    loaded = load(_write(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\n" + body, "ftdna.csv"))
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["format_id"], facet["chip"], facet["chip_basis"]) == ("ftdna", "illumina_omniexpress", "inferred")
    assert any(w.startswith("chip inferred from the rows as illumina_omniexpress: ") for w in loaded.warnings)


def test_a_small_ftdna_export_says_why_its_chip_is_unknown(tmp_path: Path) -> None:
    loaded = load(_write(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\n" + '"rs3094315","1","752566","AG"\n', "ftdna.csv"))
    assert loaded.facets.collect()["chip"].to_list() == [None]
    assert any(w.startswith("chip unknown: 0 of the export's rsids are in the marker sample") for w in loaded.warnings)


def test_the_packaged_sample_is_derived_from_the_manifests(manifests_dir: Path) -> None:
    """Opt-in: `JUST_VCF_MANIFESTS_DIR` holds Illumina's manifest zips (docs/PROVENANCE.md names them)."""
    root = manifests_dir
    script = Path(__file__).resolve().parent.parent / "scripts" / "chip_markers.py"
    spec = importlib.util.spec_from_file_location("chip_markers_script", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    packaged = json.loads(resources.files("just_vcf").joinpath("chip_markers.json").read_text())["chips"]
    for chip, entry in packaged.items():
        found = sorted(Path(root).rglob(entry["file"]))
        if not found:
            pytest.skip(f"{entry['file']} is not under JUST_VCF_MANIFESTS_DIR")
        assert hashlib.sha256(found[0].read_bytes()).hexdigest() == entry["sha256"]
        descriptor, numbers = module.manifest_rsids(found[0])
        assert descriptor == entry["manifest"] and len(numbers) == entry["rsids"]
        assert [n for n in numbers if n % MARKER_MODULUS == 0] == entry["sample"], chip
