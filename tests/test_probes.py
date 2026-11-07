"""Chip manifests (RM281): reading a manifest, locating its indel probes, resolving an array's I/D codes."""

import io
import json
import random
import zipfile
from collections import Counter
from importlib import resources
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, normalize, probes, resolve_ref
from just_vcf.cache import cache_dir
from just_vcf.formats import get_spec
from just_vcf.reference import _orient_code_rows, resolve_inputs
from just_vcf.store import Store, StoreConflictError, migrate_store
from just_vcf.probes import (
    ID_CALIBRATION,
    MIN_ORIENTATION_PEOPLE,
    ORIENTATION_CHIPS,
    IdGate,
    IdOrientation,
    IndelSite,
    Orientation,
    ProbeSet,
    ProbeSetKind,
    available_manifests,
    catalog_digest,
    chip_absence,
    find_manifest,
    id_calibration,
    id_orientation,
    learn_probe_sets,
    locate_indel,
    orientation_digest,
    orientation_json,
    orientation_payload,
    orientation_rule,
    orientation_tables,
    place_indel,
    probe_set,
    probe_table,
    reverse_complement,
)
from test_load_rsid import _two_exports
from test_reference import PREAMBLE_23, SEQ, array_text, assert_gt_indexes_ref, write_fasta
from test_store import SUBJECT
from test_store_migration import _old_store

CONTIG = SEQ["1"]
FLANK = 60

HEADER = [
    "IlmnID", "Name", "IlmnStrand", "SNP", "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq",
    "GenomeBuild", "Chr", "MapInfo", "Ploidy", "Species", "Source", "SourceVersion", "SourceStrand", "SourceSeq",
    "TopGenomicSeq", "BeadSetID", "Exp_Clusters", "Intensity_Only", "RefStrand",
]  # fmt: skip


def probe_row(
    name: str, snp: str, chrom: str, map_info: int, source_strand: str, source_seq: str, ref_strand: str = "+", build: str = "37"
) -> list[str]:
    values = dict.fromkeys(HEADER, "")
    values.update(
        IlmnID=f"{name}-0_P_F_1", Name=name, IlmnStrand="PLUS", SNP=snp, GenomeBuild=build, Chr=chrom, MapInfo=str(map_info),
        SourceStrand=source_strand, SourceSeq=source_seq, RefStrand=ref_strand,
    )  # fmt: skip
    return [values[column] for column in HEADER]


def write_manifest(directory: Path, rows: list[list[str]], descriptor: str = "TEST-1-0_A1") -> Path:
    """An Illumina manifest zip in the layout of the published ones: `[Heading]`, `[Assay]`, `[Controls]`."""
    lines = ["Illumina, Inc.", "[Heading]", f"Descriptor File Name,{descriptor}.bpm", "Assay Format,Infinium HTS", "[Assay]", ",".join(HEADER)]
    lines += [",".join(row) for row in rows]
    lines += ["[Controls]", "0027630314,Staining,Red,DNP (High)"]
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{descriptor.lower()}-manifest-file-csv.zip"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(f"{descriptor}.csv", "\r\n".join(lines) + "\r\n")
    path.write_bytes(buffer.getvalue())
    return path


def deletion(gap: int, length: int) -> tuple[str, str]:
    """A probe whose inserted bases are the reference's own at `gap` (0-based): the reference carries them."""
    return CONTIG[gap : gap + length], f"{CONTIG[gap - FLANK : gap]}[-/{CONTIG[gap : gap + length]}]{CONTIG[gap + length : gap + length + FLANK]}"


def insertion(gap: int, inserted: str) -> str:
    """A probe whose bases go in before the reference's base at `gap`: the reference lacks them."""
    return f"{CONTIG[gap - FLANK : gap]}[-/{inserted}]{CONTIG[gap : gap + FLANK]}"


def haplotype(chrom_seq: str, pos: int, ref: str, alt: str) -> str:
    """The contig with REF at 1-based `pos` replaced by ALT, the test's own reading of a VCF allele."""
    assert chrom_seq[pos - 1 : pos - 1 + len(ref)] == ref
    return chrom_seq[: pos - 1] + alt + chrom_seq[pos - 1 + len(ref) :]


# --- locating a probe by its flanks ---------------------------------------------------------------


def test_a_deletion_site_is_found_where_the_reference_carries_the_bases() -> None:
    seq, source = deletion(500, 3)
    left, right = source.split("[")[0], source.split("]")[1]
    site = locate_indel(CONTIG[300:800], 300, left, seq, right, "PLUS")
    assert site == IndelSite(gap=500, seq=seq, reference_carries=True)


def test_an_insertion_site_is_found_where_the_flanks_meet() -> None:
    left, right = CONTIG[940:1000], CONTIG[1000:1060]
    site = locate_indel(CONTIG[700:1300], 700, left, "TTGCAG", right, "PLUS")
    assert site == IndelSite(gap=1000, seq="TTGCAG", reference_carries=False)


def test_a_minus_strand_source_sequence_is_located_on_the_plus_strand() -> None:
    seq = CONTIG[500:503]
    left, right = reverse_complement(CONTIG[503:563]), reverse_complement(CONTIG[440:500])
    site = locate_indel(CONTIG[300:800], 300, left, reverse_complement(seq), right, "MINUS")
    assert site == IndelSite(gap=500, seq=seq, reference_carries=True)


def test_a_probe_found_only_on_the_strand_its_manifest_does_not_give_is_not_located() -> None:
    seq = CONTIG[500:503]
    left, right = reverse_complement(CONTIG[503:563]), reverse_complement(CONTIG[440:500])
    assert locate_indel(CONTIG[300:800], 300, left, reverse_complement(seq), right, "PLUS") == "strand"


def test_flanks_found_nowhere_or_twice_locate_nothing() -> None:
    left, right = CONTIG[440:500], CONTIG[503:563]
    assert locate_indel(CONTIG[1500:2000], 1500, left, CONTIG[500:503], right, "PLUS") == "not_found"
    doubled = CONTIG[300:800] + CONTIG[300:800]
    assert locate_indel(doubled, 300, left, CONTIG[500:503], right, "PLUS") == "several"


# --- spelling a site at the export's position -----------------------------------------------------


def test_a_site_at_the_exports_position_is_padded_with_the_base_after_it() -> None:
    window, start = CONTIG[400:700], 400
    site = IndelSite(gap=500, seq=CONTIG[500:503], reference_carries=True)
    ref, alt = place_indel(window, start, site, pos=501)
    assert (ref, alt) == (CONTIG[500:504], CONTIG[503])
    assert haplotype(CONTIG, 501, ref, alt) == CONTIG[:500] + CONTIG[503:]


def test_a_site_after_the_exports_position_keeps_the_bases_between() -> None:
    window, start = CONTIG[400:700], 400
    site = IndelSite(gap=500, seq="TTGCAG", reference_carries=False)
    ref, alt = place_indel(window, start, site, pos=498)
    assert ref == CONTIG[497:500] and alt == CONTIG[497:500] + "TTGCAG"
    assert haplotype(CONTIG, 498, ref, alt) == CONTIG[:500] + "TTGCAG" + CONTIG[500:]


def test_a_site_in_a_repeat_is_shifted_to_the_exports_position() -> None:
    genome = "GATTACA" + "CTG" * 4 + "TTAGGC"
    # Deleting one CTG unit from the start of the run is the same haplotype as deleting the last.
    site = IndelSite(gap=7, seq="CTG", reference_carries=True)
    ref, alt = place_indel(genome, 0, site, pos=12)
    assert haplotype(genome, 12, ref, alt) == genome[:7] + genome[10:]
    assert ref[0] == genome[11]


def test_a_site_before_the_exports_position_outside_a_repeat_is_not_spelled() -> None:
    window, start = CONTIG[400:700], 400
    site = IndelSite(gap=500, seq="TTGCAG", reference_carries=False)
    if CONTIG[500] == "T":
        pytest.skip("the synthetic contig allows a shift here")
    assert place_indel(window, start, site, pos=520) == "before_pos"


# --- the manifest ---------------------------------------------------------------------------------


def test_a_manifest_is_found_by_its_descriptor_and_read_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seq, source = deletion(500, 3)
    rows = [
        probe_row("rs100", "[I/D]", "1", 501, "PLUS", source),
        probe_row("GSA-rs200", "[A/G]", "2", 77, "TOP", "ACGT[A/G]TTGA", ref_strand="-"),
        probe_row("1:900-CT", "[T/C]", "1", 900, "BOT", "ACGT[T/C]TTGA"),
        probe_row("unplaced", "[A/C]", "0", 0, "TOP", "ACGT[A/C]TTGA"),
    ]
    write_manifest(tmp_path / "manifests" / "test", rows)
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(tmp_path / "manifests"))
    [manifest] = available_manifests()
    assert manifest.descriptor == "TEST-1-0_A1" and find_manifest("TEST-1-0_A1") == manifest
    table = probe_table(manifest)
    assert table["name"].to_list() == ["rs100", "GSA-rs200", "1:900-CT", "unplaced"]
    assert table["rsid"].to_list() == ["rs100", "rs200", None, None]
    assert table["kind"].to_list() == ["indel", "snv", "snv", "snv"]
    assert table["chrom"].to_list() == ["1", "2", "1", None]
    assert table["assembly_label"].unique().to_list() == ["GRCh37"]
    indel = table.row(0, named=True)
    assert (indel["left_flank"], indel["indel_seq"], indel["right_flank"]) == (CONTIG[440:500], seq, CONTIG[503:563])
    assert table.row(1, named=True)["ref_strand"] == "-" and table.row(1, named=True)["alleles"] == ["A", "G"]
    cached = list((cache_dir() / "probes").glob(f"{manifest.sha256}-*.parquet"))
    assert len(cached) == 1
    assert probe_table(manifest).equals(table)


def test_no_configured_directory_means_no_manifest(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")
    assert available_manifests() == [] and find_manifest("GSA-24v3-0_A1") is None


# --- resolving an export's I/D codes --------------------------------------------------------------


@pytest.fixture
def grch37_store(tmp_path: Path) -> ReferenceStore:
    s = ReferenceStore.open(tmp_path / "grch37")
    s.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    return s


DELETED = CONTIG[500:503]
INSERTED = "TTGCAG"


@pytest.fixture
def manifests(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    rows = [
        probe_row("rs100", "[I/D]", "1", 501, "PLUS", deletion(500, 3)[1]),
        probe_row("GSA-rs200", "[D/I]", "1", 1001, "PLUS", insertion(1000, INSERTED)),
        probe_row("rs300", "[I/D]", "1", 1501, "PLUS", f"{'A' * 60}[-/C]{'G' * 60}"),  # flanks the reference lacks
    ]
    directory = tmp_path / "manifests"
    write_manifest(directory, rows)
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(directory))
    return directory


def export(tmp_path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    """A 23andMe export of a few rows: its chip is unknown (too few rows to name the generation)."""
    path = tmp_path / "genome.txt"
    path.write_text(array_text(rows, PREAMBLE_23))
    return path


def coded_export(tmp_path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    """An AncestryDNA export, whose header states its chip (`V2.0`: `ancestrydna_v2`, calibrated for
    I/D codes). A one-letter genotype is a half-call, the second allele written `0`."""
    body = "".join(f"{rsid}\t{chrom}\t{pos}\t{gt[0]}\t{gt[1] if len(gt) > 1 else '0'}\n" for rsid, chrom, pos, gt in rows)
    path = tmp_path / "ancestry.txt"
    path.write_text(ANCESTRY_PREAMBLE + body)
    return path


def test_id_codes_are_resolved_against_the_manifests_indel_probes(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    rows = [
        ("rs100", "1", 501, "II"),  # deletion site: I is REF
        ("rs100", "1", 501, "DI"),
        ("rs100", "1", 501, "DD"),
        ("rs200", "1", 1001, "DD"),  # insertion site: D is REF
        ("rs200", "1", 1001, "ID"),
        ("rs200", "1", 1001, "I"),  # half-call: one allele written, the other `0`
        ("rs300", "1", 1501, "II"),  # probe not located in the reference
        ("rs400", "1", 1801, "DD"),  # no probe carries the rsid
        ("i5", "1", 2001, "DI"),  # an internal id: nothing to match by
    ]
    resolved = resolve_ref(load(coded_export(tmp_path, rows)), grch37_store)
    calls = resolved.calls.collect()
    by = {(row["ids"][0], row["source_gt"].rstrip("0")): row for row in calls.iter_rows(named=True)}

    long_del, short_del = CONTIG[500:504], CONTIG[503]
    assert by[("rs100", "II")]["ref"] == long_del and by[("rs100", "II")]["alts"] == [] and by[("rs100", "II")]["gt"] == [0, 0]
    assert by[("rs100", "DI")]["alts"] == [short_del] and by[("rs100", "DI")]["gt"] == [1, 0]
    assert by[("rs100", "DI")]["gt_alleles"] == [short_del, long_del]
    assert by[("rs100", "DD")]["gt"] == [1, 1] and by[("rs100", "DD")]["end_pos"] == 501 + len(long_del) - 1
    assert haplotype(CONTIG, 501, long_del, short_del) == CONTIG[:500] + CONTIG[503:]

    short_ins, long_ins = CONTIG[1000], INSERTED + CONTIG[1000]
    assert by[("rs200", "DD")]["ref"] == short_ins and by[("rs200", "DD")]["gt"] == [0, 0]
    assert by[("rs200", "ID")]["alts"] == [long_ins] and by[("rs200", "ID")]["gt"] == [1, 0]
    assert by[("rs200", "I")]["gt"] == [1, None] and by[("rs200", "I")]["gt_alleles"] == [long_ins, None]
    assert haplotype(CONTIG, 1001, short_ins, long_ins) == CONTIG[:1000] + INSERTED + CONTIG[1000:]

    for unread in (("rs300", "II"), ("rs400", "DD"), ("i5", "DI")):
        row = by[unread]
        assert row["ref"] is None and row["gt"] is None and row["gt_alleles"] == list(unread[1]), unread
    assert all(row["source_gt"].startswith(key[1]) for key, row in by.items()), "the codes stay in source_gt"
    assert_gt_indexes_ref(calls)
    warnings = "\n".join(resolved.warnings)
    assert "6 of 9 array row(s) with an I/D indel code resolved against the indel probes of TEST-1-0_A1" in warnings
    assert "2 array row(s) with an I/D indel code left unread: no indel probe of the manifests carries one of their rsids" in warnings
    assert "1 array row(s) with an I/D indel code left unread: their probe's flanks are not found" in warnings


def test_resolved_codes_normalize_to_the_site_they_name(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    rows = [("rs100", "1", 501, "DI"), ("rs200", "1", 1001, "ID")]
    normalized = normalize(resolve_ref(load(coded_export(tmp_path, rows)), grch37_store), grch37_store).calls.collect().sort("pos")
    for row, expected in zip(normalized.iter_rows(named=True), (CONTIG[:500] + CONTIG[503:], CONTIG[:1000] + INSERTED + CONTIG[1000:])):
        assert len(row["alts"]) == 1 and row["ref"][0] == row["alts"][0][0], "anchored on its first base"
        assert haplotype(CONTIG, row["pos"], row["ref"], row["alts"][0]) == expected


def test_codes_stay_unread_without_a_manifest(tmp_path: Path, grch37_store: ReferenceStore, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")
    resolved = resolve_ref(load(coded_export(tmp_path, [("rs100", "1", 501, "DI"), ("rs1", "1", 10, CONTIG[9] * 2)])), grch37_store)
    calls = resolved.calls.collect()
    code = calls.filter(pl.col("source_gt") == "DI").row(0, named=True)
    assert code["ref"] is None and code["gt_alleles"] == ["D", "I"]
    assert any(
        w == "1 array row(s) with an I/D indel code left unread: no manifest with indel probes on GRCh37 is configured (JUST_VCF_MANIFESTS_DIR)"
        for w in resolved.warnings
    )


def test_a_manifest_on_another_build_reads_no_codes(tmp_path: Path, grch37_store: ReferenceStore, monkeypatch: pytest.MonkeyPatch) -> None:
    directory = tmp_path / "manifests36"
    write_manifest(directory, [probe_row("rs100", "[I/D]", "1", 501, "PLUS", deletion(500, 3)[1], build="36")])
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(directory))
    resolved = resolve_ref(load(coded_export(tmp_path, [("rs100", "1", 501, "DI")])), grch37_store)
    assert resolved.calls.collect()["ref"].to_list() == [None]
    assert any("no manifest with indel probes on GRCh37 is configured" in w for w in resolved.warnings)


def test_a_store_tells_codes_read_with_manifests_from_codes_left_unread(tmp_path: Path, grch37_store: ReferenceStore, monkeypatch: pytest.MonkeyPatch) -> None:
    """`resolve_ref` keeps the `facet_id` and its `RESOLVE_VERSION` is the same either way, but the rows
    differ: resolved without manifests an array's I/D rows stay unread, with them they are read. A store
    holding the first was told the second was `present`, and kept the unread rows. The facet now records
    the manifests that read it (`resolve_inputs`), so the store refuses the other state; a facet that had
    nothing for a manifest to read records none, and stays one state with or without them."""
    path = coded_export(tmp_path, [("rs100", "1", 501, "DI"), ("rs1", "1", 10, CONTIG[9] * 2)])
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")
    unread = resolve_ref(load(path, subject=SUBJECT), grch37_store)
    directory = tmp_path / "manifests"
    write_manifest(directory, [probe_row("rs100", "[I/D]", "1", 501, "PLUS", deletion(500, 3)[1])])
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(directory))
    read = resolve_ref(load(path, subject=SUBJECT), grch37_store)
    facet_id = read.facets.collect()["facet_id"].item()
    assert unread.facets.collect()["facet_id"].item() == facet_id
    code = pl.col("source_gt") == "DI"
    assert unread.calls.filter(code).collect()["ref"].to_list() == [None] and read.calls.filter(code).collect()["ref"].to_list() != [None]

    store = Store.open(tmp_path / "store", subject=SUBJECT)
    assert store.add(unread).added == (facet_id,)
    with pytest.raises(StoreConflictError, match=facet_id):
        store.add(read)
    assert store.add(unread).present == (facet_id,)

    assert unread.facets.collect()["resolve_inputs"].to_list() == [None]
    assert read.facets.collect()["resolve_inputs"].to_list() == [resolve_inputs(catalog_digest(available_manifests()), None)]
    (tmp_path / "plain").mkdir()
    plain = coded_export(tmp_path / "plain", [("rs1", "1", 10, CONTIG[9] * 2)])
    assert resolve_ref(load(plain), grch37_store).facets.collect()["resolve_inputs"].to_list() == [None], "no code row: no manifest read it"


def test_an_array_resolved_before_0021_has_its_manifests_named_unknown(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    """A 0.0.20 store recorded no manifests for the rows it holds: the step to 0.0.21 fills
    `resolve_inputs` null on a resolved array facet and names it, so the null is not read as "no
    manifest". The same facet resolved today with manifests is then another state of it."""
    resolved = resolve_ref(load(coded_export(tmp_path, [("rs100", "1", 501, "DI")]), subject=SUBJECT), grch37_store)
    root = tmp_path / "store"
    _old_store(root, "0.0.20", [resolved])
    facet_id = resolved.facets.collect()["facet_id"].item()
    report = migrate_store(root)
    assert [(entry.version, entry.frame, entry.columns, entry.facet_ids) for entry in report.filled_unknown] == [
        ("0.0.21", "facets", ("resolve_inputs",), (facet_id,))
    ]
    assert Store.open(root).loaded().facets.collect()["resolve_inputs"].to_list() == [None]
    with pytest.raises(StoreConflictError, match="no resolve_inputs recorded"):
        Store.open(root).add(resolved)


def test_codes_of_a_chip_whose_orientation_is_unverified_are_withheld(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    rows = [("rs100", "1", 501, "DI"), ("rs200", "1", 1001, "DD")]
    loaded = load(coded_export(tmp_path, rows))
    resolved = resolve_ref(loaded, grch37_store)
    assert resolved.calls.collect()["ref"].null_count() == 0, "ancestrydna_v2 is calibrated: its codes resolve"
    for chip in ("23andme_v3", "23andme_v4", "23andme_v1_v2", "illumina:TEST-1-0_A1", None):
        relabelled = loaded.model_copy(update={"facets": loaded.facets.with_columns(pl.lit(chip, dtype=pl.String).alias("chip"))})
        withheld = resolve_ref(relabelled, grch37_store)
        calls = withheld.calls.collect().sort("pos")
        assert calls["ref"].to_list() == [None, None] and calls["gt"].to_list() == [None, None], chip
        assert calls["gt_alleles"].to_list() == [["D", "I"], ["D", "D"]], chip
        line = next(w for w in withheld.warnings if "id_orientation_unverified" in w)
        assert line.startswith("2 array row(s) with an I/D indel code left unread (id_orientation_unverified)"), chip
        assert not any("resolved against the indel probes" in w for w in withheld.warnings), chip
    unknown = loaded.model_copy(update={"facets": loaded.facets.with_columns(pl.lit(None, dtype=pl.String).alias("chip"))})
    assert any("the facet's chip is unknown" in w for w in resolve_ref(unknown, grch37_store).warnings)


def test_every_calibration_names_its_evidence() -> None:
    by_mode = {mode: {chip for chip, entry in ID_CALIBRATION.items() if entry.mode is mode} for mode in IdGate}
    assert by_mode[IdGate.CHIP_WIDE] == {"23andme_v5", "ancestrydna_v2"}
    assert by_mode[IdGate.PER_PROBE] == {"23andme_v3", "23andme_v4"} == ORIENTATION_CHIPS
    assert by_mode[IdGate.WITHHELD] == {"23andme_v1_v2"}
    assert all("2026-10-23" in entry.evidence for entry in ID_CALIBRATION.values())
    assert all("2026-11-07" in ID_CALIBRATION[chip].evidence for chip in ORIENTATION_CHIPS)
    assert not id_calibration("living_dna_sirius").resolve and not id_calibration(None).resolve
    assert {chip for chip, entry in ID_CALIBRATION.items() if entry.resolve} == by_mode[IdGate.CHIP_WIDE] | by_mode[IdGate.PER_PROBE]


# --- per-probe orientation (23andMe v3, v4) -------------------------------------------------------


@pytest.mark.parametrize(
    ("hom_ref", "hom_alt", "expected"),
    [
        (3, 0, Orientation.SAME),
        (13, 0, Orientation.SAME),
        (0, 3, Orientation.FLIPPED),
        (0, 30, Orientation.FLIPPED),
        (2, 0, None),  # too few people
        (0, 2, None),
        (3, 1, None),  # both kinds of evidence
        (1, 3, None),
        (11, 2, None),
        (0, 0, None),
    ],
)
def test_the_orientation_rule_case_by_case(hom_ref: int, hom_alt: int, expected: Orientation | None) -> None:
    assert orientation_rule(hom_ref, hom_alt) is expected


def test_the_rule_is_min_orientation_people_one_way_and_none_the_other() -> None:
    """The rule against its own definition over a grid, so a change to the threshold moves both."""
    for hom_ref in range(MIN_ORIENTATION_PEOPLE + 3):
        for hom_alt in range(MIN_ORIENTATION_PEOPLE + 3):
            same = hom_ref >= MIN_ORIENTATION_PEOPLE and hom_alt == 0
            flipped = hom_alt >= MIN_ORIENTATION_PEOPLE and hom_ref == 0
            assert orientation_rule(hom_ref, hom_alt) is (Orientation.SAME if same else Orientation.FLIPPED if flipped else None)


EVIDENCE = {
    "23andme_v4": {"rs9": (4, 0), "rs10": (0, 3), "rs2": (2, 0), "rs7": (3, 1), "rs1": (0, 0)},
    "23andme_v3": {"rs5": (0, 5)},
}


def test_the_payload_lists_only_what_the_rule_gives_and_writes_deterministically() -> None:
    payload = orientation_payload(EVIDENCE, {"23andme_v4": 13, "23andme_v3": 23}, "synthetic evidence", "2026-11-07")
    tables = orientation_tables(json.loads(orientation_json(payload)))
    expected = {
        chip: {rsid: orientation_rule(*counts) for rsid, counts in by_rsid.items() if orientation_rule(*counts) is not None}
        for chip, by_rsid in EVIDENCE.items()
    }
    assert {chip: {rsid: entry.orientation for rsid, entry in table.items()} for chip, table in tables.items()} == expected
    assert tables["23andme_v4"]["rs10"] == IdOrientation(Orientation.FLIPPED, 0, 3)
    assert payload["chips"]["23andme_v4"]["pairs"] == 13 and payload["provenance"]["min_orientation_people"] == MIN_ORIENTATION_PEOPLE

    reordered = {chip: dict(reversed(list(by_rsid.items()))) for chip, by_rsid in reversed(list(EVIDENCE.items()))}
    again = orientation_payload(reordered, {"23andme_v3": 23, "23andme_v4": 13}, "synthetic evidence", "2026-11-07")
    text = orientation_json(payload)
    assert orientation_json(again) == text and text.endswith("}\n") and "\n" not in text[:-1]
    assert list(tables["23andme_v4"]) == sorted(tables["23andme_v4"])

    empty = orientation_payload({}, {}, "nothing", "2026-11-07")
    assert orientation_tables(empty) == {chip: {} for chip in ORIENTATION_CHIPS}


def test_a_table_the_rule_does_not_give_is_refused() -> None:
    def payload() -> dict:
        return json.loads(orientation_json(orientation_payload(EVIDENCE, {}, "synthetic evidence", "2026-11-07")))

    edited = payload()
    edited["chips"]["23andme_v4"]["probes"]["rs10"]["orientation"] = "same"  # 0 hom-ref, 3 hom-alt: flipped
    with pytest.raises(ValueError, match="rs10 is listed 'same'"):
        orientation_tables(edited)
    edited = payload()
    added = {**edited["chips"]["23andme_v4"]["probes"], "rs7": {"orientation": "same", "hom_ref_people": 3, "hom_alt_people": 1}}
    edited["chips"]["23andme_v4"]["probes"] = dict(sorted(added.items()))
    with pytest.raises(ValueError, match="rs7"):
        orientation_tables(edited)
    edited = payload()
    probes = edited["chips"]["23andme_v4"]["probes"]
    edited["chips"]["23andme_v4"]["probes"] = dict(reversed(list(probes.items())))
    with pytest.raises(ValueError, match="not sorted"):
        orientation_tables(edited)
    edited = payload()
    del edited["chips"]["23andme_v3"]
    with pytest.raises(ValueError, match="covers"):
        orientation_tables(edited)
    edited = payload()
    edited["provenance"]["min_orientation_people"] = MIN_ORIENTATION_PEOPLE + 1
    with pytest.raises(ValueError, match="min_orientation_people"):
        orientation_tables(edited)
    with pytest.raises(ValueError, match="23andme_v5"):
        orientation_payload({"23andme_v5": {"rs1": (3, 0)}}, {}, "pooled", "2026-11-07")


def test_the_packaged_table_validates() -> None:
    text = resources.files("just_vcf").joinpath("id_orientation.json").read_text()
    payload = json.loads(text)
    assert orientation_json(payload) == text, "written as orientation_json writes it"
    tables = orientation_tables(payload)
    assert set(tables) == ORIENTATION_CHIPS
    assert {"learned_from", "date", "rule", "min_orientation_people", "variant_window"} <= set(payload["provenance"])
    for chip in ORIENTATION_CHIPS:
        assert id_orientation(chip) == tables[chip]
        assert payload["chips"][chip]["pairs"] >= 0
    assert id_orientation("23andme_v5") == {}, "a chip resolved chip-wide has no table"


def with_table(monkeypatch: pytest.MonkeyPatch, evidence: dict[str, dict[str, tuple[int, int]]]) -> None:
    """The packaged table replaced by one learned from `evidence`, read through the same validation."""
    payload = json.loads(orientation_json(orientation_payload(evidence, {}, "synthetic evidence", "2026-11-07")))
    monkeypatch.setattr(probes, "_packaged_tables", lambda: orientation_tables(payload))


@pytest.fixture
def oriented_manifests(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Three locatable probes: rs100 a deletion site, rs200 an insertion site, rs500 another deletion."""
    rows = [
        probe_row("rs100", "[I/D]", "1", 501, "PLUS", deletion(500, 3)[1]),
        probe_row("GSA-rs200", "[D/I]", "1", 1001, "PLUS", insertion(1000, INSERTED)),
        probe_row("rs500", "[I/D]", "1", 2001, "PLUS", deletion(2000, 2)[1]),
    ]
    directory = tmp_path / "manifests"
    write_manifest(directory, rows)
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(directory))
    return directory


ORIENTED_ROWS = [
    ("rs100", "1", 501, "II"),
    ("rs100", "1", 501, "DI"),
    ("rs200", "1", 1001, "DD"),
    ("rs200", "1", 1001, "ID"),
    ("rs500", "1", 2001, "DI"),
    ("rs500", "1", 2001, "DD"),
]


def as_chip(loaded, chip: str | None):
    return loaded.model_copy(update={"facets": loaded.facets.with_columns(pl.lit(chip, dtype=pl.String).alias("chip"))})


def by_codes(resolved) -> dict[tuple[str, str], dict]:
    return {(row["ids"][0], row["source_gt"]): row for row in resolved.calls.collect().iter_rows(named=True)}


@pytest.mark.parametrize("chip", sorted(ORIENTATION_CHIPS))
def test_a_per_probe_chip_reads_each_listed_probe_its_way_and_withholds_the_rest(
    tmp_path: Path, grch37_store: ReferenceStore, oriented_manifests: Path, monkeypatch: pytest.MonkeyPatch, chip: str
) -> None:
    with_table(monkeypatch, {chip: {"rs100": (3, 0), "rs200": (0, 4)}})
    loaded = load(coded_export(tmp_path, ORIENTED_ROWS))
    as_written = by_codes(resolve_ref(as_chip(loaded, "23andme_v5"), grch37_store))
    oriented = resolve_ref(as_chip(loaded, chip), grch37_store)
    rows = by_codes(oriented)
    slots = ("ref", "alts", "gt", "end_pos", "gt_alleles")

    # rs100 `same`: read as written, as a chip-wide chip reads it.
    long_del, short_del = CONTIG[500:504], CONTIG[503]
    for codes in ("II", "DI"):
        assert [rows[("rs100", codes)][s] for s in slots] == [as_written[("rs100", codes)][s] for s in slots]
    assert rows[("rs100", "II")]["ref"] == long_del and rows[("rs100", "II")]["gt"] == [0, 0]
    assert haplotype(CONTIG, 501, long_del, short_del) == CONTIG[:500] + CONTIG[503:]

    # rs200 `flipped`: the export's D is Illumina's I. Over an insertion site (the reference lacks the
    # bases) Illumina's D is REF, so the export's DD is hom-alt with the bases inserted, its ID het.
    short_ins, long_ins = CONTIG[1000], INSERTED + CONTIG[1000]
    flipped_dd, flipped_id = rows[("rs200", "DD")], rows[("rs200", "ID")]
    assert (flipped_dd["ref"], flipped_dd["alts"], flipped_dd["gt"], flipped_dd["gt_alleles"]) == (short_ins, [long_ins], [1, 1], [long_ins, long_ins])
    assert (flipped_id["gt"], flipped_id["gt_alleles"]) == ([0, 1], [short_ins, long_ins])
    assert haplotype(CONTIG, 1001, short_ins, long_ins) == CONTIG[:1000] + INSERTED + CONTIG[1000:]
    # The same rows as a chip-wide chip reads the codes swapped by hand.
    (tmp_path / "swapped").mkdir()
    by_hand = by_codes(resolve_ref(as_chip(load(coded_export(tmp_path / "swapped", [("rs200", "1", 1001, "II"), ("rs200", "1", 1001, "DI")])), "23andme_v5"), grch37_store))
    for codes, written in {"DD": "II", "ID": "DI"}.items():
        assert [rows[("rs200", codes)][s] for s in slots] == [by_hand[("rs200", written)][s] for s in slots]
    assert all(row["source_gt"] == codes for (_, codes), row in rows.items()), "the codes stay in source_gt as the export wrote them"

    # rs500 is not listed: left as loaded, though its probe would resolve it.
    for codes in ("DI", "DD"):
        assert as_written[("rs500", codes)]["ref"] == CONTIG[2000:2003]
        row = rows[("rs500", codes)]
        assert row["ref"] is None and row["gt"] is None and row["gt_alleles"] == list(codes)
    assert_gt_indexes_ref(oriented.calls.collect())

    warnings = "\n".join(oriented.warnings)
    assert "4 of 6 array row(s) with an I/D indel code resolved against the indel probes of TEST-1-0_A1" in warnings
    assert "2 of them with I and D swapped first" in warnings
    assert (
        f"2 array row(s) with an I/D indel code left unread (id_orientation_unverified): {chip} is read per probe, and no WGS "
        f"evidence of this probe's orientation is in its table (id_orientation.json, 2 probe(s) of {chip} listed)"
    ) in warnings


def test_a_chip_resolved_chip_wide_ignores_the_orientation_table(
    tmp_path: Path, grch37_store: ReferenceStore, oriented_manifests: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    loaded = load(coded_export(tmp_path, ORIENTED_ROWS))
    before = {chip: resolve_ref(as_chip(loaded, chip), grch37_store) for chip in ("ancestrydna_v2", "23andme_v5")}
    with_table(monkeypatch, {"23andme_v4": {"rs100": (0, 3), "rs200": (0, 3)}})
    for chip, resolved in before.items():
        after = resolve_ref(as_chip(loaded, chip), grch37_store)
        assert after.calls.collect().equals(resolved.calls.collect()), chip
        assert after.warnings == resolved.warnings, chip
        assert after.facets.collect()["resolve_inputs"].to_list() == [resolve_inputs(catalog_digest(available_manifests()), None)], chip
        assert not any("swapped" in w for w in after.warnings), chip


def test_the_orientation_table_is_part_of_a_per_probe_facets_resolve_inputs(
    tmp_path: Path, grch37_store: ReferenceStore, oriented_manifests: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Two tables that read a probe two ways give two REF states of one facet; a table whose counts
    alone differ reads the same rows, and is the same state."""
    loaded = as_chip(load(coded_export(tmp_path, ORIENTED_ROWS)), "23andme_v4")

    def inputs(evidence: dict[str, tuple[int, int]]) -> str | None:
        with_table(monkeypatch, {"23andme_v4": evidence})
        return resolve_ref(loaded, grch37_store).facets.collect()["resolve_inputs"].item()

    same, flipped, more_people = inputs({"rs200": (3, 0)}), inputs({"rs200": (0, 3)}), inputs({"rs200": (9, 0)})
    chip_wide = resolve_inputs(catalog_digest(available_manifests()), None)
    assert len({same, flipped, chip_wide}) == 3 and None not in (same, flipped)
    assert more_people == same
    with_table(monkeypatch, {"23andme_v4": {"rs200": (3, 0)}})
    assert same == resolve_inputs(catalog_digest(available_manifests(), orientation_digest(id_orientation("23andme_v4"))), None)


@pytest.mark.parametrize("chip", sorted(ORIENTATION_CHIPS))
def test_an_empty_table_withholds_a_per_probe_chip_as_the_chip_wide_gate_did(
    tmp_path: Path, grch37_store: ReferenceStore, oriented_manifests: Path, monkeypatch: pytest.MonkeyPatch, chip: str
) -> None:
    """Until the table is learned, a v3 or v4 export resolves as it did when the chip was withheld:
    every code row as loaded and counted, no manifest read, no manifest named in its REF state."""
    with_table(monkeypatch, {})
    loaded = load(coded_export(tmp_path, ORIENTED_ROWS))
    withheld = resolve_ref(as_chip(loaded, chip), grch37_store)
    calls = withheld.calls.collect()
    assert calls["ref"].null_count() == calls.height and calls["gt"].null_count() == calls.height
    assert calls.sort("pos", "source_gt")["gt_alleles"].to_list() == [list(gt) for _, _, _, gt in sorted(ORIENTED_ROWS, key=lambda r: (r[2], r[3]))]
    code_lines = [w for w in withheld.warnings if "I/D indel code" in w]
    assert len(code_lines) == 1 and code_lines[0].startswith("6 array row(s) with an I/D indel code left unread (id_orientation_unverified)")
    assert f"which lists no probe of {chip} yet" in code_lines[0]
    assert withheld.facets.collect()["resolve_inputs"].to_list() == [None]
    assert list((cache_dir() / "probes").glob("*.parquet")) == [], "no manifest was read"


def test_rows_are_oriented_by_the_table_and_a_row_listed_both_ways_is_withheld() -> None:
    table = {"rs1": IdOrientation(Orientation.SAME, 3, 0), "rs2": IdOrientation(Orientation.FLIPPED, 0, 3)}
    rows = pl.DataFrame(
        {
            "_row": [0, 1, 2, 3, 4],
            "ids": [["rs1"], ["rs2"], ["rs1", "rs2"], ["rs9"], ["rs2"]],
            "gt_alleles": [["I", "D"], ["I", "D"], ["I", "I"], ["D", "D"], ["I", None]],
        },
        schema={"_row": pl.UInt32, "ids": pl.List(pl.String), "gt_alleles": pl.List(pl.String)},
    )
    counts: Counter[str] = Counter()
    oriented, flipped = _orient_code_rows(rows, table, counts)
    assert dict(zip(oriented["_row"].to_list(), oriented["gt_alleles"].to_list())) == {0: ["I", "D"], 1: ["D", "I"], 4: ["D", None]}
    assert flipped == {1, 4}
    assert counts == Counter({"code_id_orientation_unverified": 1, "code_orientation_ambiguous": 1})


# --- probe sets and the cause of an absence -------------------------------------------------------

ANCESTRY_PREAMBLE = (
    "#AncestryDNA raw data download\n"
    "#This file was generated by AncestryDNA at: 01/01/2019 00:00:00 UTC\n"
    "#Data was collected using AncestryDNA array version: V2.0\n"
    "#Data is formatted using AncestryDNA converter version: V1.0\n"
    "#of the SNP using human reference build 37.1 coordinates.  Columns four and five \n"
    "rsid\tchromosome\tposition\tallele1\tallele2\n"
)


def ancestry_export(path: Path, positions: list[int]) -> Path:
    path.write_text(ANCESTRY_PREAMBLE + "".join(f"rs{p}\t1\t{p}\tA\tG\n" for p in positions))
    return path


def test_a_probe_set_is_learned_from_a_chips_exports_and_read_back(tmp_path: Path) -> None:
    exports = [ancestry_export(tmp_path / f"e{i}.txt", [100, 200 + (i == 2), 300 + i]) for i in range(3)]
    [learned] = learn_probe_sets(exports)
    assert (learned.chip, learned.kind, learned.assembly_label) == ("ancestrydna_v2", ProbeSetKind.LEARNED, "GRCh37")
    assert sorted(learned.sites.rows()) == [("1", 100), ("1", 200)], "a site one export of three carries is not on the chip"
    read = probe_set("ancestrydna_v2", "GRCh37")
    assert isinstance(read, ProbeSet) and read.sites.equals(learned.sites) and read.source == "learned from 3 export(s)"


def test_too_few_exports_learn_nothing(tmp_path: Path) -> None:
    exports = [ancestry_export(tmp_path / f"e{i}.txt", [100, 200]) for i in range(2)]
    assert learn_probe_sets(exports) == []
    assert probe_set("ancestrydna_v2", "GRCh37") == "no probe set for ancestrydna_v2 on GRCh37 has been learned (learn_probe_sets)"


def test_an_absence_on_the_chip_is_a_no_call_and_off_it_not_on_the_chip(tmp_path: Path) -> None:
    learn_probe_sets([ancestry_export(tmp_path / f"e{i}.txt", [100, 200, 300]) for i in range(3)])
    loaded = load(ancestry_export(tmp_path / "person.txt", [100, 300]))
    asked = pl.DataFrame({"chrom": ["1", "1", "1", "2"], "pos": [100, 200, 250, 100]})
    answer = chip_absence(loaded, asked)
    assert answer.sites["absence_cause"].to_list() == ["in_export", "no_call", "not_on_chip", "not_on_chip"]
    assert answer.probe_set is not None and answer.reason.startswith("read against the learned probe set of ancestrydna_v2 on GRCh37")


def test_a_load_of_several_facets_is_refused(tmp_path: Path) -> None:
    """An appended export is two facets: one facet's chip read against both facets' rows would mix them."""
    text, _, _ = _two_exports()
    path = tmp_path / "appended.txt"
    path.write_text(text)
    appended = load(path)
    assert appended.facets.collect().height == 2
    with pytest.raises(ValueError, match="chip_absence takes one facet; this Loaded has 2"):
        chip_absence(appended, pl.DataFrame({"chrom": ["1"], "pos": [1001]}))
    exports = [ancestry_export(tmp_path / f"e{i}.txt", [100, 200]) for i in range(3)]
    with pytest.raises(ValueError, match="learn_probe_sets takes one facet; this Loaded has 2"):
        learn_probe_sets([*exports, appended])
    assert not (cache_dir() / "probes" / "learned").exists(), "nothing is written before the refusal"


def test_a_written_no_call_is_a_no_call_not_a_site_in_the_export(tmp_path: Path) -> None:
    learn_probe_sets([ancestry_export(tmp_path / f"e{i}.txt", [100, 200, 300]) for i in range(3)])
    person = tmp_path / "person.txt"
    person.write_text(ANCESTRY_PREAMBLE + "rs100\t1\t100\tA\tG\nrs200\t1\t200\t0\t0\nrs300\t1\t300\tA\t0\n")
    asked = pl.DataFrame({"chrom": ["1", "1", "1"], "pos": [100, 200, 300]})
    assert chip_absence(load(person), asked).sites["absence_cause"].to_list() == ["in_export", "no_call", "in_export"]
    # 23andMe's spelling, and no probe set at all: the written row itself shows the chip probed there.
    twentythree = export(tmp_path, [("rs100", "1", 100, "AG"), ("rs200", "1", 200, "--"), ("rs201", "1", 200, "--")])
    answer = chip_absence(load(twentythree), pl.DataFrame({"chrom": ["1", "1", "1"], "pos": [100, 200, 300]}))
    assert answer.probe_set is None
    assert answer.sites["absence_cause"].to_list() == ["in_export", "no_call", None]


@pytest.mark.parametrize(
    "format_id",
    ["vcf_no_meta", "casava_snp_table", "casava_indel_table", "soapsnp_cns", "ftdna_bigy",
     "nirvana_json", "color_variant_set", "get_evidence_gff", "cgi_var", "cgi_master_var"],
)  # fmt: skip
def test_a_sequencing_facet_that_states_no_vcf_version_is_not_an_array(tmp_path: Path, format_id: str) -> None:
    """These formats state no VCF version and load with no callset type, as an array export can; only
    a genotype table is one. A site such a file has no row at is not a chip's no-call."""
    path = tmp_path / f"{format_id}.txt"
    path.write_text(get_spec(format_id).example)
    loaded = load(path)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["format_id"] == format_id and facet["vcf_version"] is None and facet["callset_type"] is None
    first = loaded.calls.select("chrom", pl.col("pos").cast(pl.Int64)).head(1).collect()
    answer = chip_absence(loaded, pl.concat([first, pl.DataFrame({"chrom": ["1"], "pos": [1]})]))
    assert answer.reason == "the facet is not an array export" and answer.probe_set is None
    assert answer.sites["absence_cause"].null_count() == answer.sites.height == 2


def test_a_site_with_a_genotype_in_any_of_its_rows_is_in_the_export(tmp_path: Path) -> None:
    twentythree = export(tmp_path, [("rs200", "1", 200, "--"), ("rs201", "1", 200, "CT")])
    answer = chip_absence(load(twentythree), pl.DataFrame({"chrom": ["1"], "pos": [200]}))
    assert answer.sites["absence_cause"].to_list() == ["in_export"]


def test_without_a_probe_set_an_absence_has_no_cause(tmp_path: Path) -> None:
    loaded = load(ancestry_export(tmp_path / "person.txt", [100, 300]))
    answer = chip_absence(loaded, pl.DataFrame({"chrom": ["1", "1"], "pos": [100, 200]}))
    assert answer.sites["absence_cause"].to_list() == ["in_export", None]
    assert answer.probe_set is None and "has been learned" in answer.reason


def test_a_manifest_named_by_the_chip_is_its_probe_set(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rows = [
        probe_row("rs1", "[A/G]", "1", 752566, "TOP", "ACGT[A/G]TTGA"),
        probe_row("rs2", "[A/G]", "1", 752721, "TOP", "ACGT[A/G]TTGA"),
        probe_row("rs3", "[A/G]", "1", 900000, "TOP", "ACGT[A/G]TTGA", build="36"),
        probe_row("rs4", "[A/G]", "0", 0, "TOP", "ACGT[A/G]TTGA"),
    ]
    write_manifest(tmp_path / "manifests", rows, descriptor="GSA-24v2-0_A1")
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(tmp_path / "manifests"))
    found = probe_set("illumina:GSA-24v2-0_A1", "GRCh37")
    assert isinstance(found, ProbeSet) and found.kind is ProbeSetKind.MANIFEST
    assert found.sites.rows() == [("1", 752566), ("1", 752721)], "only the probes the manifest maps to that build"
    assert probe_set("illumina:GSA-24v3-0_A1", "GRCh37") == "the manifest GSA-24v3-0_A1 is not configured (JUST_VCF_MANIFESTS_DIR)"
    assert probe_set("illumina:GSA-24v2-0_A1", "GRCh38") == "the manifest GSA-24v2-0_A1 maps no probe to GRCh38"


def test_a_platform_family_has_no_probe_set() -> None:
    reason = probe_set("illumina_omniexpress", "GRCh37")
    assert isinstance(reason, str) and "platform family" in reason
    assert probe_set(None, "GRCh37") == "the facet's chip is unknown"


# --- a mixed-strand export read by its probes' allele pairs ---------------------------------------

STRAND_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def probed_rows(seed: int = 5, per_contig: int = 700) -> tuple[list[tuple[str, str, int, str]], dict[str, list[str]]]:
    """Array rows drawn around the synthetic reference, forward strand, and each probe's plus-strand pair.

    One probe in ten is an A/T or C/G probe, which reads alike on both strands."""
    rng = random.Random(seed)
    rows, pairs = [], {}
    for chrom in ("1", "2"):
        for pos in range(1, per_contig + 1):
            ref = SEQ[chrom][pos - 1]
            palindromic = rng.random() < 0.1
            alt = STRAND_COMPLEMENT[ref] if palindromic else rng.choice([b for b in "ACGT" if b not in (ref, STRAND_COMPLEMENT[ref])])
            kind = rng.random()
            gt = ref * 2 if kind < 0.5 else ref + alt if kind < 0.8 else alt * 2
            rsid = f"rs{chrom}{pos}"
            rows.append((rsid, chrom, pos, gt))
            pairs[rsid] = sorted([ref, alt])
    return rows, pairs


def complemented_every(rows: list[tuple[str, str, int, str]], every: int = 6) -> list[tuple[str, str, int, str]]:
    """`rows` with every `every`-th genotype written as a minus-strand probe writes it: complemented."""
    return [(r, c, p, "".join(STRAND_COMPLEMENT[b] for b in gt) if n % every == 0 else gt) for n, (r, c, p, gt) in enumerate(rows)]


def test_a_mixed_strand_export_is_read_by_its_probes_allele_pairs(tmp_path: Path, grch37_store: ReferenceStore, monkeypatch: pytest.MonkeyPatch) -> None:
    truth, pairs = probed_rows()
    rows = complemented_every(truth)
    manifest_rows = []
    for n, (rsid, chrom, pos, _) in enumerate(truth):
        # Half the probes are designed on the minus strand: the manifest writes their pair complemented.
        design = pairs[rsid] if n % 2 else [STRAND_COMPLEMENT[b] for b in pairs[rsid]]
        manifest_rows.append(probe_row(rsid, f"[{design[0]}/{design[1]}]", chrom, pos, "TOP", "ACGT[A/G]TTGA", ref_strand="+" if n % 2 else "-"))
    write_manifest(tmp_path / "manifests", manifest_rows)
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(tmp_path / "manifests"))

    resolved = resolve_ref(load(export(tmp_path, rows)), grch37_store)
    facets, calls = resolved.facets.collect(), resolved.calls.collect()
    assert facets["strand_check"].cast(pl.String).to_list() == ["mixed"]
    by_id = {row["ids"][0]: row for row in calls.iter_rows(named=True)}
    written = {rsid: gt for rsid, _, _, gt in rows}
    read = flipped_back = 0
    for rsid, _, _, gt in truth:
        row, ref = by_id[rsid], by_id[rsid]["ref"]
        if pairs[rsid] in (["A", "T"], ["C", "G"]):
            safe = len(set(written[rsid])) == 2 and ref in written[rsid]
            assert (row["alts"] is not None) == safe, rsid  # only strand-safe hets are read on such a probe
            continue
        assert row["alts"] == sorted(set(gt) - {ref}), rsid
        assert row["gt_alleles"] == list(gt), "a minus-strand row is complemented back to the forward strand"
        assert row["source_gt"] == written[rsid], "the export's text is kept"
        read += 1
        flipped_back += written[rsid] != gt
    assert read and flipped_back
    assert_gt_indexes_ref(calls)
    line = next(w for w in resolved.warnings if "read by their probe's two alleles on the plus strand" in w)
    assert "still withheld: " in line and " on A/T or C/G probes" in line


def test_a_mixed_strand_export_without_manifests_withholds_as_before(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    resolved = resolve_ref(load(export(tmp_path, complemented_every(probed_rows()[0]))), grch37_store)
    calls = resolved.calls.collect()
    homozygous = calls.filter(pl.col("gt_alleles").list.n_unique() == 1)
    assert homozygous["alts"].null_count() == homozygous.height, "a homozygous row cannot be read without the probe's pair"
    assert "strand of withheld rows not read by probe: no manifest is configured (JUST_VCF_MANIFESTS_DIR)" in resolved.warnings
