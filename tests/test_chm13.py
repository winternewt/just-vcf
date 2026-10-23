"""T2T-CHM13 v2.0 against the real assembly (RM286). Opt-in: a reference store holding GRCh38 and
T2T-CHM13 (`JUST_VCF_TEST_REFERENCE_DIR`) and UCSC's hg38↔hs1 chains (`JUST_VCF_TEST_CHAIN_DIR`).

Expected values come from the store's sequence, the carried NCBI assembly report and pyliftover's
reading of the same chain, never from what the code under test wrote.
"""

import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, liftover, load, resolve_ref
from just_vcf.assembly import ASSEMBLY_LENGTHS
from just_vcf.reference import Contig
from just_vcf.schema import canonical_chrom

REPORT = Path(__file__).resolve().parent.parent / "assets" / "ncbi_assembly_reports" / "GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt"
TO_CHM13 = "hg38ToHs1.over.chain.gz"
TO_GRCH38 = "hs1ToHg38.over.chain.gz"
AUTOSOMES = [str(n) for n in range(1, 23)]


def report_lengths() -> dict[str, int]:
    """UCSC-style name → length, from the carried assembly report."""
    rows = [line.split("\t") for line in REPORT.read_text().splitlines() if line and not line.startswith("#")]
    return {row[9]: int(row[8]) for row in rows}


@pytest.fixture(scope="module")
def store(test_reference_dir: Path) -> ReferenceStore:
    found = ReferenceStore.open(test_reference_dir)
    labels = {found.assembly_label(digest) for digest in found.collections()}
    if not {"GRCh38", "T2T-CHM13"} <= labels:
        pytest.skip(f"{test_reference_dir} holds no GRCh38 and T2T-CHM13 collections")
    return found


def collection(store: ReferenceStore, label: str) -> str:
    (digest,) = [digest for digest in store.collections() if store.assembly_label(digest) == label]
    return digest


def by_chrom(store: ReferenceStore, digest: str) -> dict[str, Contig]:
    return {chrom: contigs[0] for chrom, contigs in store.contigs_by_chrom(digest).items() if len(contigs) == 1}


def test_the_t2t_collection_is_chm13v2_and_its_own(store: ReferenceStore) -> None:
    digest = collection(store, "T2T-CHM13")
    contigs = store.contigs(digest)
    assert {contig.name: contig.length for contig in contigs} == report_lengths()
    for contig in contigs:
        chrom = canonical_chrom(contig.name)
        assert {contig.length} <= ASSEMBLY_LENGTHS["T2T-CHM13"][chrom]
    # Charter 3: another assembly is another collection, and no T2T sequence is a GRCh38 one.
    grch38 = collection(store, "GRCh38")
    assert digest != grch38
    assert not {c.digest for c in contigs} & {c.digest for c in store.contigs(grch38)}


def sampled_positions(contigs: dict[str, Contig], per_chrom: int, seed: int) -> list[tuple[str, int]]:
    rng = random.Random(seed)
    return [(chrom, rng.randrange(1_000_000, contigs[chrom].length - 1_000_000)) for chrom in AUTOSOMES for _ in range(per_chrom)]


def bases_at(store: ReferenceStore, contigs: dict[str, Contig], positions: list[tuple[str, int]]) -> list[str]:
    out = []
    for chrom, pos in positions:
        out.append(store.fetch(contigs[chrom], [(pos - 1, pos)])[0].upper())
    return out


def test_an_array_drawn_on_t2t_is_named_t2t_by_its_rows(store: ReferenceStore, tmp_path: Path) -> None:
    """F9 against the real collections: rows whose genotypes are drawn around CHM13's bases."""
    digest = collection(store, "T2T-CHM13")
    contigs = by_chrom(store, digest)
    positions = sampled_positions(contigs, 80, seed=5)
    rng = random.Random(9)
    lines = ["# rsid\tchromosome\tposition\tgenotype"]
    for index, ((chrom, pos), ref) in enumerate(zip(positions, bases_at(store, contigs, positions), strict=True)):
        if ref not in "ACGT":
            continue
        alt = rng.choice([b for b in "ACGT" if b != ref])
        kind = rng.random()
        gt = ref * 2 if kind < 0.5 else ref + alt if kind < 0.8 else alt * 2
        lines.append(f"rs{index}\t{chrom}\t{pos}\t{gt}")
    path = tmp_path / "genome.txt"
    path.write_text("\n".join(lines) + "\n")
    loaded = load(path)
    assert loaded.facets.collect()["assembly_label"].to_list() == [None]
    resolved = resolve_ref(loaded, store)
    facets = resolved.facets.collect()
    assert (facets["assembly_label"].item(), facets["ref_seqcol_digest"].item()) == ("T2T-CHM13", digest)
    line = next(w for w in resolved.warnings if w.startswith("assembly inferred from the rows as T2T-CHM13"))
    assert f"on GRCh38 ({collection(store, 'GRCh38')}" in line


def grch38_vcf(store: ReferenceStore, tmp_path: Path, per_chrom: int = 40) -> Path:
    """A GRCh38 VCF of SNVs at sampled positions, REF read from the store, contigs as the store names them."""
    digest = collection(store, "GRCh38")
    contigs = by_chrom(store, digest)
    positions = sampled_positions(contigs, per_chrom, seed=11)
    meta = ["##fileformat=VCFv4.2", *[f"##contig=<ID={contigs[c].name},length={contigs[c].length}>" for c in AUTOSOMES]]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = []
    for (chrom, pos), ref in sorted(zip(positions, bases_at(store, contigs, positions), strict=True), key=lambda item: (AUTOSOMES.index(item[0][0]), item[0][1])):
        if ref in "ACGT":
            alt = next(b for b in "ACGT" if b != ref)
            body.append(f"{contigs[chrom].name}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1")
    path = tmp_path / "grch38.vcf"
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def test_grch38_to_t2t_and_back_agree_with_pyliftover(store: ReferenceStore, test_chain_dir: Path, tmp_path: Path) -> None:
    pyliftover = pytest.importorskip("pyliftover")
    to_chm13, to_grch38 = test_chain_dir / TO_CHM13, test_chain_dir / TO_GRCH38
    if not (to_chm13.is_file() and to_grch38.is_file()):
        pytest.skip(f"{test_chain_dir} lacks {TO_CHM13} or {TO_GRCH38}")
    parent = resolve_ref(load(grch38_vcf(store, tmp_path)), store)
    assert parent.facets.collect()["assembly_label"].item() == "GRCh38"
    lifted = liftover(parent, to_chm13, store)
    facet = lifted.loaded.facets.collect().row(0, named=True)
    assert (facet["assembly_label"], facet["ref_seqcol_digest"], facet["relation"]) == ("T2T-CHM13", collection(store, "T2T-CHM13"), "liftover_to")

    check_against_pyliftover(lifted.report.collect(), pyliftover.LiftOver(str(to_chm13)))
    # REF on the target was read from CHM13: where it is neither allele the row stays behind, counted.
    assert set(lifted.reasons) <= {"ref_mismatch_after_lift", "no_chain_block", "multiple_targets"}

    back = liftover(lifted.loaded, to_grch38, store)
    assert back.loaded.facets.collect()["assembly_label"].item() == "GRCh38"
    check_against_pyliftover(back.report.collect(), pyliftover.LiftOver(str(to_grch38)))
    original = {(c, p) for c, p in parent.calls.select("chrom", "pos").collect().iter_rows()}
    round_trip = {(c, p) for c, p in back.loaded.calls.select("chrom", "pos").collect().iter_rows()}
    assert round_trip and len(round_trip & original) > 0


def check_against_pyliftover(report: pl.DataFrame, lifter) -> None:
    """Every lifted row lands where pyliftover puts its position (0-based in, 0-based out), on the same
    strand; a row with no block has no pyliftover hit; one with several targets has several distinct hits."""
    lifted = 0
    for chrom, pos, lift_chrom, lift_pos, strand, reason in report.select("chrom", "pos", "lift_chrom", "lift_pos", "lift_strand", "lift_reason").iter_rows():
        hits = {(canonical_chrom(name), zero + 1, hit_strand) for name, zero, hit_strand, _ in lifter.convert_coordinate(f"chr{chrom}", pos - 1)}
        if reason is None:
            lifted += 1
            assert hits == {(lift_chrom, lift_pos, strand)}, (chrom, pos)
        elif reason == "no_chain_block":
            assert hits == set(), (chrom, pos)
        elif reason == "multiple_targets":
            assert len(hits) > 1, (chrom, pos)
        elif reason == "ref_mismatch_after_lift":
            assert len(hits) == 1, (chrom, pos)
    assert lifted == report.filter(pl.col("lift_reason").is_null()).height > 0
