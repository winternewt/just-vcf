"""Chip manifests (RM281): reading a manifest, locating its indel probes, resolving an array's I/D codes."""

import io
import random
import zipfile
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, normalize, resolve_ref
from just_vcf.cache import cache_dir
from just_vcf.probes import (
    IndelSite,
    ProbeSet,
    ProbeSetKind,
    available_manifests,
    chip_absence,
    find_manifest,
    learn_probe_sets,
    locate_indel,
    place_indel,
    probe_set,
    probe_table,
    reverse_complement,
)
from test_reference import PREAMBLE_23, SEQ, array_text, assert_gt_indexes_ref, write_fasta

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
    path = tmp_path / "genome.txt"
    path.write_text(array_text(rows, PREAMBLE_23))
    return path


def test_id_codes_are_resolved_against_the_manifests_indel_probes(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    rows = [
        ("rs100", "1", 501, "II"),  # deletion site: I is REF
        ("rs100", "1", 501, "DI"),
        ("rs100", "1", 501, "DD"),
        ("rs200", "1", 1001, "DD"),  # insertion site: D is REF
        ("rs200", "1", 1001, "ID"),
        ("rs200", "1", 1001, "I"),  # half-call or hemizygous: one letter
        ("rs300", "1", 1501, "II"),  # probe not located in the reference
        ("rs400", "1", 1801, "DD"),  # no probe carries the rsid
        ("i5", "1", 2001, "DI"),  # an internal id: nothing to match by
    ]
    resolved = resolve_ref(load(export(tmp_path, rows)), grch37_store)
    calls = resolved.calls.collect()
    by = {(row["ids"][0], row["source_gt"]): row for row in calls.iter_rows(named=True)}

    long_del, short_del = CONTIG[500:504], CONTIG[503]
    assert by[("rs100", "II")]["ref"] == long_del and by[("rs100", "II")]["alts"] == [] and by[("rs100", "II")]["gt"] == [0, 0]
    assert by[("rs100", "DI")]["alts"] == [short_del] and by[("rs100", "DI")]["gt"] == [1, 0]
    assert by[("rs100", "DI")]["gt_alleles"] == [short_del, long_del]
    assert by[("rs100", "DD")]["gt"] == [1, 1] and by[("rs100", "DD")]["end_pos"] == 501 + len(long_del) - 1
    assert haplotype(CONTIG, 501, long_del, short_del) == CONTIG[:500] + CONTIG[503:]

    short_ins, long_ins = CONTIG[1000], INSERTED + CONTIG[1000]
    assert by[("rs200", "DD")]["ref"] == short_ins and by[("rs200", "DD")]["gt"] == [0, 0]
    assert by[("rs200", "ID")]["alts"] == [long_ins] and by[("rs200", "ID")]["gt"] == [1, 0]
    assert by[("rs200", "I")]["gt"] == [1] and by[("rs200", "I")]["gt_alleles"] == [long_ins]
    assert haplotype(CONTIG, 1001, short_ins, long_ins) == CONTIG[:1000] + INSERTED + CONTIG[1000:]

    for unread in (("rs300", "II"), ("rs400", "DD"), ("i5", "DI")):
        row = by[unread]
        assert row["ref"] is None and row["gt"] is None and row["gt_alleles"] == list(unread[1]), unread
    assert all(row["source_gt"] == key[1] for key, row in by.items()), "the codes stay in source_gt"
    assert_gt_indexes_ref(calls)
    warnings = "\n".join(resolved.warnings)
    assert "6 of 9 array row(s) with an I/D indel code resolved against the indel probes of TEST-1-0_A1" in warnings
    assert "2 array row(s) with an I/D indel code left unread: no indel probe of the manifests carries one of their rsids" in warnings
    assert "1 array row(s) with an I/D indel code left unread: their probe's flanks are not found" in warnings


def test_resolved_codes_normalize_to_the_site_they_name(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    rows = [("rs100", "1", 501, "DI"), ("rs200", "1", 1001, "ID")]
    normalized = normalize(resolve_ref(load(export(tmp_path, rows)), grch37_store), grch37_store).calls.collect().sort("pos")
    for row, expected in zip(normalized.iter_rows(named=True), (CONTIG[:500] + CONTIG[503:], CONTIG[:1000] + INSERTED + CONTIG[1000:])):
        assert len(row["alts"]) == 1 and row["ref"][0] == row["alts"][0][0], "anchored on its first base"
        assert haplotype(CONTIG, row["pos"], row["ref"], row["alts"][0]) == expected


def test_codes_stay_unread_without_a_manifest(tmp_path: Path, grch37_store: ReferenceStore, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")
    resolved = resolve_ref(load(export(tmp_path, [("rs100", "1", 501, "DI"), ("rs1", "1", 10, CONTIG[9] * 2)])), grch37_store)
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
    resolved = resolve_ref(load(export(tmp_path, [("rs100", "1", 501, "DI")])), grch37_store)
    assert resolved.calls.collect()["ref"].to_list() == [None]
    assert any("no manifest with indel probes on GRCh37 is configured" in w for w in resolved.warnings)


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
    exports = [ancestry_export(tmp_path / f"e{i}.txt", [100, 200, 300 + i]) for i in range(3)]
    [learned] = learn_probe_sets(exports)
    assert (learned.chip, learned.kind, learned.assembly_label) == ("ancestrydna_v2", ProbeSetKind.LEARNED, "GRCh37")
    assert sorted(learned.sites.rows()) == [("1", 100), ("1", 200), ("1", 300), ("1", 301), ("1", 302)]
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
