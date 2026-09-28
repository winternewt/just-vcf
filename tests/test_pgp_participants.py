"""Tests for scripts/pgp_participants.py: manifest parsing, collection triage and per-person dedup."""

from pathlib import Path

import pytest

import pgp_leech as pl
import pgp_participants as pp


def test_manifest_paths_join_stream_and_unescape_spaces() -> None:
    # Shape of a real CGI delivery manifest: a root stream and a nested one, with block locators first.
    text = (
        ". a642c3ecd26956b88e6feb08d12a5c65+67108864 0:20:version 20:1024:README.2.0.1.txt\n"
        "./GS02995-DNA_C02/ASM 3e527e7a22595cd8bcbc80399cf3afc9+67108864 0:3000:var-GS000036655-ASM.tsv.bz2"
        " 3000:10:odd\\040name.txt\n"
    )
    assert pp.parse_manifest(text) == [
        ("version", 20),
        ("README.2.0.1.txt", 1024),
        ("GS02995-DNA_C02/ASM/var-GS000036655-ASM.tsv.bz2", 3000),
        ("GS02995-DNA_C02/ASM/odd name.txt", 10),
    ]


@pytest.mark.parametrize(
    ("name", "kind"),
    [("AD9JZ4FY9DC.vcf", "vcf"), ("WGC069886D vcf", "vcf"), ("hu577152.GS000043235-DID", "cgi"),
     ("copy of Raw sequencing data GS03052-DNA_H01", "cgi"), ("GS000039609-ASM", "cgi"),
     ("AD9JZ4FY9DC.bam", None), ("WGC069886D bam", None),
     ("Output 54e7cb61 of GETEvidenceReportHTML of su92l-d1hrv-zdxey2zsibvspyt", None)],
)
def test_collection_names_seen_on_the_cluster(name: str, kind: str | None) -> None:
    assert pp.classify_collection(name) == kind


@pytest.mark.parametrize(
    ("base", "keep"),
    [("var-GS000036655-ASM.tsv.bz2", True), ("masterVarBeta-GS000036655-ASM.tsv.bz2", True),
     ("vcfBeta-GS000036655-ASM.vcf.bz2", True), ("summary-GS000036655-ASM.tsv", True),
     ("var-GS000036655-ASM.tsv_with_wellcount.bz2", False), ("dbSNPAnnotated-GS000036655-ASM.tsv.bz2", False),
     ("coverageRefScore-chr1-GS000037344-ASM.tsv.bz2", False), ("evidenceDnbs-chr2-GS000037344-ASM.tsv.bz2", False)],
)
def test_cgi_keep_list(base: str, keep: bool) -> None:
    assert bool(pp.CGI_KEEP.match(base)) is keep


@pytest.mark.parametrize(
    ("name", "reads"),
    [("AE2CH6SK4DG-EXT.chr1.bam - BAM", True), ("WGC069886D.bam.tgz) - BAM", True),
     ("AE2CH6SK4DG-EXT - VCF", False), ("huA2692E-veritas-gVCF-2019.vcf.bz2", False)],
)
def test_reads_are_recognised_by_upload_name(name: str, reads: bool) -> None:
    assert bool(pp.READS_BY_NAME.search(name)) is reads


def _arv(coll: str, name: str, file: str, size: int, hu: str = "hu3F7E26") -> pp.Item:
    return pp.Item(hu, "arvados", coll, name, "vcf", file, size, f"https://x/{coll}/{file}")


def test_plan_keeps_the_original_over_a_copy_and_folds_same_size_upload() -> None:
    copy = _arv("ac2it-4zz18-aaaa", "copy of AD9JZ4FY9DC.vcf", "AD9JZ4FY9DC.vcf.gz", 700)
    orig = _arv("ac2it-4zz18-zzzz", "AD9JZ4FY9DC.vcf", "AD9JZ4FY9DC.vcf.gz", 700)
    upload_same = pp.Item("hu3F7E26", "upload", "4001", "AD9JZ4FY9DC - VCF", "Veritas Genetics", None, 700, None)
    upload_other = pp.Item("hu3F7E26", "upload", "4002", "23andMe", "23andMe", None, 5_000, None)
    other_person = _arv("ac2it-4zz18-bbbb", "AD9JZ4FY9DC.vcf", "AD9JZ4FY9DC.vcf.gz", 700, hu="hu715C3F")
    kept, dup = pp.plan([copy, upload_same, orig, upload_other, other_person])
    assert set(kept) == {orig, upload_other, other_person}
    assert dup == {copy: orig, upload_same: orig}


def test_migration_moves_by_type_files_into_participant_dirs(tmp_path: Path) -> None:
    old = tmp_path / "23andMe" / "hu23502F" / "2" / "hu23502F_20110106182158.zip"
    old.parent.mkdir(parents=True)
    old.write_bytes(b"PK")
    (tmp_path / "manifest.jsonl").write_text("{}\n")
    root = tmp_path / "participants"
    assert pp.migrate_by_type_layout(tmp_path, root) == 1
    assert (root / "hu23502F" / "upload-2" / old.name).read_bytes() == b"PK"
    assert not (tmp_path / "23andMe").exists()
    assert (tmp_path / "manifest.legacy-by-type.jsonl").exists()
    # A migrated upload directory satisfies the leecher's offline resume check for that row.
    row = pl.ListingRow(2, "hu23502F", None, 1215, "2011-01-06", "23andMe", "Participant", "SNP Data from 23andMe",
                        "e9e7b1d959510ac01dd0107df2e5cc54+69", 2)
    assert pl._already_present(row, root / "hu23502F" / "upload-2") is not None


def test_plan_folds_an_upload_whose_listed_size_is_wrong_by_its_sample_id() -> None:
    # Real case, hu6942AE: upload 2021 "WGC069895D - VCF" is listed at 1 byte; the bytes equal the dataset file.
    arv = pp.Item("hu6942AE", "arvados", "ac2it-4zz18-hphf6foa8ax1aru", "WGC069895D vcf", "vcf", "WGC069895D.vcf.gz",
                  344669369, "https://x")
    upload = pp.Item("hu6942AE", "upload", "2021", "WGC069895D - VCF", "Veritas Genetics", None, 1, None)
    kept, dup = pp.plan([upload, arv])
    assert kept == [arv] and dup == {upload: arv}
