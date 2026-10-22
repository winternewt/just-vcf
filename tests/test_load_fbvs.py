"""PGP-HMS freeBayes+VarScan2 merges (`vcf_freebayes_varscan2`): one facet per caller, split on INFO/SF.

The layouts are the three docs/PGP_HMS.md counted over the arvados merges: two populated sample
columns (13 of 110), one declared sample (48), and two declared samples with one carried (50). A
PGP-HMS merge read 2026-10-22 had SF `0,1` on 84% of its first 300,000 rows.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff
from test_load_vcf import assert_gt_consistent, write

META = [
    "##fileformat=VCFv4.1",
    "##fileDate=20160919",
    "##source=freeBayes v0.9.21-24-g840b412-dirty",
    "##reference=/vbin/ref/hg19.fa",
    "##phasing=none",
    "##source=VarScan2",
    "##contig=<ID=chr1,length=249250621>",
    '##INFO=<ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">',
    '##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled">',
    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">',
    '##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">',
    "##source_20160919.1=vcf-merge(r840) Sample_X-EXT.chr1.FB.vcf.gz Sample_X-EXT.chr1.VS.vcf.gz",
    "##sourceFiles_20160919.1=0:Sample_X-EXT.chr1.FB.vcf.gz,1:Sample_X-EXT.chr1.VS.vcf.gz",
]
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
FB_SHARED = "chr1\t752566\t.\tG\tA\t425.3\tPASS\tAC=2;SF=0,1\tGT:GL\t0/1:-30.1,0,-25.2"
FB_ONLY = "chr1\t752721\t.\tA\tG\t96.5\t.\tAC=1;SF=0\tGT:GL\t1/1:-50,-9,0"
VS_ONLY = "chr1\t776546\t.\tA\tG\t.\tPASS\tAC=1;SF=1\tGT:GQ:AD:FREQ\t0/1:41:8:28.57%"
VS_FILTERED = "chr1\t800000\t.\tC\tT\t.\tstr10\tAC=1;SF=1f\tGT:GQ:AD:FREQ\t0/1:12:3:10%"


def frames(path: Path):
    loaded = load(path)
    facets, calls, blocks = pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])
    by_caller = {row["caller"]: row for row in facets.iter_rows(named=True)}
    rows = {caller: calls.filter(pl.col("facet_id") == facet["facet_id"]) for caller, facet in by_caller.items()}
    return loaded, facets, by_caller, rows, blocks


def test_one_data_column(tmp_path: Path) -> None:
    text = "\n".join([*META, f"{COLUMNS}\tSample1", FB_SHARED, FB_ONLY, VS_ONLY, VS_FILTERED]) + "\n"
    path = write(tmp_path, text, "merge.vcf")
    assert sniff(path).format == "vcf_freebayes_varscan2"
    loaded, facets, by_caller, rows, blocks = frames(path)
    assert facets.height == 2 and blocks.height == 0
    assert facets["facet_id"].to_list() == sorted(facets["facet_id"].to_list())
    fb, vs = by_caller["freeBayes"], by_caller["VarScan2"]
    assert (fb["caller_version"], vs["caller_version"]) == ("0.9.21-24-g840b412-dirty", None)
    assert {fb["relation"], vs["relation"]} == {"caller_split"} and fb["derived_from"] is None
    assert fb["source_sha256"] == vs["source_sha256"] and fb["vendor"] == "PGP-HMS"
    assert (fb["absence_means"], vs["absence_means"]) == ("unknown", "unknown")
    assert rows["freeBayes"]["pos"].to_list() == [752566, 752721]
    assert rows["VarScan2"]["pos"].to_list() == [776546, 800000], "the shared site carries freeBayes's genotype only"
    assert rows["freeBayes"]["gl"].to_list()[0] == pytest.approx([-30.1, 0.0, -25.2], rel=1e-5)
    assert rows["VarScan2"]["gq"].to_list() == [41, 12] and rows["VarScan2"]["ad"].to_list() == [None, None]
    assert rows["VarScan2"]["filters"].to_list() == [["PASS"], ["str10"]]
    assert any(w.startswith("VarScan2 facet: one data column: at the 1 site(s)") for w in loaded.warnings)
    assert "FORMAT/AD" in fb["dropped_fields"], "VarScan2's AD is the ALT depth (Number=1)"
    calls = loaded.calls.collect()
    assert calls.height == 4 and calls["facet_id"].to_list() == sorted(calls["facet_id"].to_list())
    assert_gt_consistent(calls)


def test_two_populated_columns(tmp_path: Path) -> None:
    records = [
        "chr1\t752566\t.\tG\tA\t425.3\tPASS\tAC=2;SF=0,1\tGT:GL:GQ\t0/1:-30.1,0,-25.2:.\t1/1:.:35",
        "chr1\t752721\t.\tA\tG\t96.5\t.\tAC=1;SF=0\tGT:GL:GQ\t1/1:-50,-9,0:.\t.",
        "chr1\t776546\t.\tA\tG\t.\tPASS\tAC=1;SF=1\tGT:GL:GQ\t.\t0/1:.:41",
    ]
    text = "\n".join([*META, f"{COLUMNS}\tunknown\tSample1", *records]) + "\n"
    loaded, facets, by_caller, rows, _ = frames(write(tmp_path, text, "two.vcf"))
    assert (by_caller["freeBayes"]["sample"], by_caller["VarScan2"]["sample"]) == ("unknown", "Sample1")
    assert rows["freeBayes"]["gt"].to_list() == [[0, 1], [1, 1]]
    assert rows["VarScan2"]["pos"].to_list() == [752566, 776546], "with a column per caller, VarScan2 keeps shared sites"
    assert rows["VarScan2"]["gt"].to_list() == [[1, 1], [0, 1]], "where the callers disagree, both are kept"
    assert not any("one data column" in w for w in loaded.warnings)


def test_two_declared_columns_one_carried(tmp_path: Path) -> None:
    text = "\n".join([*META, f"{COLUMNS}\tunknown\tSample1", FB_SHARED, FB_ONLY, VS_ONLY]) + "\n"
    loaded, _, by_caller, rows, _ = frames(write(tmp_path, text, "declared_two.vcf"))
    assert (by_caller["freeBayes"]["sample"], by_caller["VarScan2"]["sample"]) == ("unknown", "unknown")
    assert rows["VarScan2"]["pos"].to_list() == [776546]
    assert rows["VarScan2"]["gt"].to_list() == [[0, 1]], "the one carried column holds VarScan2's lone sites"
    assert any("one data column" in w for w in loaded.warnings)


def test_a_source_file_name_without_a_caller_token_leaves_the_caller_null(tmp_path: Path) -> None:
    meta = [line.replace(".VS.", ".x.") for line in META]
    text = "\n".join([*meta, f"{COLUMNS}\tSample1", FB_SHARED, VS_ONLY]) + "\n"
    loaded = load(write(tmp_path, text, "unnamed.vcf"))
    facets = loaded.facets.collect()
    assert set(facets["caller"].to_list()) == {"freeBayes", None}
    assert any(w.startswith("unnamed caller facet:") for w in loaded.warnings)


def test_real_merge(pgp_participants: Path) -> None:
    """The smallest PGP-HMS merge on disk, counts only: the facets cover every row once, except the
    shared sites a one-column merge gives to freeBayes alone."""
    candidates = sorted(pgp_participants.glob("*/arvados-*/*.vcf.gz"), key=lambda p: p.stat().st_size)
    found = next((p for p in candidates[:20] if sniff(p).format == "vcf_freebayes_varscan2"), None)
    if found is None:
        pytest.skip("no PGP-HMS freeBayes+VarScan2 merge on disk")
    loaded = load(found)
    facets = loaded.facets.collect()
    counts = loaded.calls.group_by("facet_id").len().collect()
    by_id = dict(zip(counts["facet_id"], counts["len"], strict=True))
    callers = dict(zip(facets["caller"], facets["facet_id"], strict=True))
    assert set(callers) == {"freeBayes", "VarScan2"}
    assert by_id[callers["freeBayes"]] > by_id[callers["VarScan2"]] > 0
    assert loaded.calls.filter(pl.col("end_pos") < pl.col("pos")).select(pl.len()).collect().item() == 0
