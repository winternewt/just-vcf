"""RM270 D: `rejoin_split` puts a pre-split VCF's rows back into one site each."""

import random
from pathlib import Path

import polars as pl

from just_vcf import load, rejoin_split
from just_vcf.rejoin import RELATION
from just_vcf.schema import CALLS

HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=249250621>",
    "##contig=<ID=2,length=243199373>",
    '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]

# One person's multi-allelic sites, as a caller writes them, and as `bcftools norm -m-` splits them:
# one row per ALT, GT slots kept, the other ALT written as 0.
UNSPLIT = [
    "1\t100\trs1\tA\tG,T\t50\tPASS\t.\tGT\t1/2",
    "1\t200\t.\tC\tG,T\t50\tPASS\t.\tGT\t1/1",
    "1\t300\t.\tA\tG,C,T\t50\tPASS\t.\tGT\t0/2",
    "1\t400\t.\tAT\tA,GT\t50\tPASS\t.\tGT\t1/2",
    "1\t500\t.\tA\tG,T\t50\tPASS\t.\tGT\t0|1",
    "2\t600\t.\tG\tA\t50\tPASS\t.\tGT\t0/1",
]
SPLIT = [
    "1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t1/0",
    "1\t100\trs1\tA\tT\t50\tPASS\t.\tGT\t0/1",
    "1\t200\t.\tC\tG\t50\tPASS\t.\tGT\t1/1",
    "1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/0",
    "1\t300\t.\tA\tG\t50\tPASS\t.\tGT\t0/0",
    "1\t300\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
    "1\t300\t.\tA\tT\t50\tPASS\t.\tGT\t0/0",
    "1\t400\t.\tAT\tA\t50\tPASS\t.\tGT\t1/0",
    "1\t400\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",  # trimmed: AT>GT written as A>G
    "1\t500\t.\tA\tG\t50\tPASS\t.\tGT\t0|1",
    "1\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0|0",
    "2\t600\t.\tG\tA\t50\tPASS\t.\tGT\t0/1",
]


def write(tmp_path: Path, records: list[str], name: str) -> Path:
    path = tmp_path / name
    path.write_text("\n".join([*HEADER, *records]) + "\n")
    return path


def site_view(calls: pl.DataFrame) -> list[tuple]:
    return [
        (row["chrom"], row["pos"], row["source_ref"], frozenset(row["alts"]), tuple(row["gt_alleles"]))
        for row in calls.iter_rows(named=True)
    ]


def test_rejoined_split_file_reads_as_the_unsplit_one(tmp_path: Path) -> None:
    unsplit = load(write(tmp_path, UNSPLIT, "unsplit.vcf")).calls.collect()
    split = load(write(tmp_path, SPLIT, "split.vcf"))
    rejoined = rejoin_split(split)
    calls = rejoined.calls.collect()
    assert calls.schema == CALLS
    # ALT order is not recoverable from split rows, so alts compare as sets; gt_alleles exactly.
    assert site_view(calls) == site_view(unsplit)
    assert calls.filter(pl.col("pos") == 100)["ids"].to_list() == [["rs1"]]
    # The joined ALTs are sorted by allele, whatever order the rows came in (here the A>G row sorts first).
    assert calls.filter(pl.col("pos") == 400)["alts"].to_list() == [["A", "GT"]]
    assert calls["split_signal"].null_count() == calls.height
    joined = calls.filter(pl.col("pos").is_in([100, 200, 300, 400, 500]))
    assert joined["gl"].null_count() == joined["ad"].null_count() == joined["source_gt"].null_count() == joined.height
    assert calls.filter(pl.col("pos") == 500)["phased"].to_list() == [[True, True]]
    assert any(w.startswith("5 split site(s) re-joined from 11 row(s)") for w in rejoined.warnings)


def test_the_parent_is_untouched_and_the_facet_is_derived(tmp_path: Path) -> None:
    split = load(write(tmp_path, SPLIT, "split.vcf"))
    before = split.calls.collect()
    rejoined = rejoin_split(split)
    parent, derived = split.facets.collect().row(0, named=True), rejoined.facets.collect().row(0, named=True)
    assert derived["derived_from"] == parent["facet_id"] and derived["relation"] == RELATION
    assert derived["facet_id"] != parent["facet_id"]
    assert set(rejoined.calls.collect()["facet_id"]) == {derived["facet_id"]}
    assert split.calls.collect().equals(before)
    assert rejoin_split(split).facets.collect().equals(rejoined.facets.collect()), "deterministic id"


def test_a_file_with_nothing_flagged_passes_through(tmp_path: Path) -> None:
    loaded = load(write(tmp_path, UNSPLIT, "unsplit.vcf"))
    rejoined = rejoin_split(loaded)
    assert rejoined.calls.collect().drop("facet_id").equals(loaded.calls.collect().drop("facet_id"))
    assert rejoined.warnings == loaded.warnings


def test_row_order_does_not_change_the_result(tmp_path: Path) -> None:
    expected = rejoin_split(load(write(tmp_path, SPLIT, "split.vcf"))).calls.collect().drop("facet_id")
    for seed in range(5):
        shuffled = SPLIT[:]
        random.Random(seed).shuffle(shuffled)
        got = rejoin_split(load(write(tmp_path, shuffled, f"shuffled{seed}.vcf"))).calls.collect().drop("facet_id")
        assert got.equals(expected), seed


def test_groups_that_cannot_be_read_as_one_site_stay_as_written(tmp_path: Path) -> None:
    tag = "OLD_MULTIALLELIC=1:{}:A/G/T"
    records = [
        f"1\t100\t.\tA\tG\t50\tPASS\t{tag.format(100)}\tGT\t0/1",  # two ALTs claim slot 1
        f"1\t100\t.\tA\tT\t50\tPASS\t{tag.format(100)}\tGT\t0/1",
        f"1\t200\t.\tA\tG\t50\tPASS\t{tag.format(200)}\tGT\t1",  # haploid beside diploid
        f"1\t200\t.\tA\tT\t50\tPASS\t{tag.format(200)}\tGT\t0/1",
        f"1\t300\t.\tA\tG\t50\tPASS\t{tag.format(300)}\tGT\t1/.",  # vt writes the other ALT as .
        f"1\t300\t.\tA\tT\t50\tPASS\t{tag.format(300)}\tGT\t./1",
        f"1\t400\t.\tA\tG\t50\tPASS\t{tag.format(400)}\tGT\t0/1",  # one allele reached from two rows
        f"1\t400\t.\tAT\tGT\t50\tPASS\t{tag.format(400)}\tGT\t0/1",
        f"1\t500\t.\tA\tG\t40\tq10\t{tag.format(500)}\tGT\t1/0",  # filters disagree: union; QUAL: null
        f"1\t500\t.\tA\tT\t50\tPASS\t{tag.format(500)}\tGT\t0/1",
        f"1\t600\t.\tA\tG\t50\tPASS\t{tag.format(600)}\tGT\t1/.",  # a slot no row called stays null
        f"1\t600\t.\tA\tT\t50\tPASS\t{tag.format(600)}\tGT\t./.",
        f"1\t700\t.\tA\tG\t50\tPASS\t{tag.format(700)}\tGT\t0/2",  # an index past the row's one ALT
        f"1\t700\t.\tA\tT\t50\tPASS\t{tag.format(700)}\tGT\t0/0",
        f"1\t800\t.\tAT\tA\t50\tPASS\t{tag.format(800)}\tGT\t1/0",  # REFs that are not prefixes
        f"1\t800\t.\tAG\tA\t50\tPASS\t{tag.format(800)}\tGT\t0/1",
        f"1\t900\t.\tA\t<DEL>\t50\tPASS\t{tag.format(900)}\tGT\t1/0",  # symbolic
        f"1\t900\t.\tA\tT\t50\tPASS\t{tag.format(900)}\tGT\t0/1",
    ]
    rejoined = rejoin_split(load(write(tmp_path, records, "edge.vcf")))
    calls = rejoined.calls.collect()
    by_pos = {pos: group for (pos,), group in calls.group_by(["pos"])}
    assert by_pos[100].height == 2 and set(by_pos[100]["split_signal"]) == {"marker"}
    assert by_pos[200].height == 2 and set(by_pos[200]["split_signal"]) == {"marker"}
    assert by_pos[300]["gt_alleles"].to_list() == [["G", "T"]]
    assert by_pos[400]["alts"].to_list() == [["GT"]] and by_pos[400]["gt_alleles"].to_list() == [["AT", "GT"]]
    assert by_pos[500]["filters"].to_list() == [["q10"]] and by_pos[500]["qual"].to_list() == [None]
    assert by_pos[600]["gt_alleles"].to_list() == [["G", None]]
    for pos in (100, 200, 700, 800, 900):
        assert by_pos[pos].height == 2 and set(by_pos[pos]["split_signal"]) == {"marker"}, pos
    assert set(rejoined.warnings) >= {
        "1 split site(s) left as written: two ALTs in one genotype slot",
        "1 split site(s) left as written: rows of different ploidy",
        "1 split site(s) left as written: a genotype index past the row's ALTs",
        "1 split site(s) left as written: REFs that are not prefixes of one another",
        "1 split site(s) left as written: a symbolic allele",
    }


def test_a_row_without_a_genotype_stays_as_written(tmp_path: Path) -> None:
    header = HEADER[:-1] + ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">', HEADER[-1]]
    records = [
        "1\t100\t.\tA\tG\t50\tPASS\tOLD_MULTIALLELIC=1:100:A/G/T\tGT\t1/0",
        "1\t100\t.\tA\tT\t50\tPASS\tOLD_MULTIALLELIC=1:100:A/G/T\tDP\t7",
    ]
    path = tmp_path / "nogt.vcf"
    path.write_text("\n".join([*header, *records]) + "\n")
    rejoined = rejoin_split(load(path))
    assert rejoined.calls.collect().height == 2
    assert "1 split site(s) left as written: a row has no genotype" in rejoined.warnings
