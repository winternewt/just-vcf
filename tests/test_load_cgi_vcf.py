"""Complete Genomics `vcfBeta` (`vcf_cgi`): its record classes, half-calls and GL scale.

The synthetic head copies the shapes counted in a PGP-HMS vcfBeta (2026-10-22, counts only): per
400,000 rows, 61,660 `<CGA_CNVWIN>` windows, 1,576 `<CGA_NOCALL>` spans with END, about 164,000 ALT `.`
records with GT `./.` and AN=0 (no-called spans over their REF), a few thousand half-calls, and no
ALT `.` record called reference. In 126,604 of 126,604 rows with GQ, GQ was the gap between the two
best GL values, which is why GL is read in phred units.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff
from test_load_vcf import assert_gt_consistent, write

HEADER = [
    "##fileformat=VCFv4.1",
    "##center=Complete Genomics",
    "##source=CGAPipeline_2.5.0.37",
    "##source_GENOME_REFERENCE=NCBI build 37",
    "##source_MAX_PLOIDY=10",
    "##phasing=partial",
    "##reference=ftp://ftp.completegenomics.com/ReferenceFiles/build37.fa.bz2",
    "##contig=<ID=1,length=249250621,assembly=B37>",
    "##contig=<ID=X,length=155270560,assembly=B37>",
    '##ALT=<ID=CGA_CNVWIN,Description="Copy number analysis window">',
    '##ALT=<ID=CGA_NOCALL,Description="No-called record">',
    '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CGA_WINEND,Number=1,Type=Integer,Description="End of coverage window">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">',
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filters">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    '##FORMAT=<ID=GL,Number=.,Type=Integer,Description="Genotype Likelihood">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">',
    '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)">',
    '##FORMAT=<ID=CGA_GP,Number=1,Type=Float,Description="Depth of coverage for 2k window GC normalized to mean">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGS000040194-ASM",
]
FMT = "GT:PS:FT:GQ:GL:DP:AD"
RECORDS = [
    "1\t1\t.\tN\t<CGA_NOCALL>\t.\t.\tEND=10000;NS=1;AN=0\tGT:PS\t./.:.",
    "1\t10001\t.\tT\t<CGA_CNVWIN>\t.\t.\tNS=1;CGA_WINEND=12000\tGT:CGA_GP\t.:0.97",
    "1\t10001\t.\tTAACCC\t.\t.\t.\tNS=1;AN=0\tGT:PS\t./.:.",
    f"1\t752566\trs3094315\tG\tA\t.\t.\tNS=1;AN=2;AC=1\t{FMT}\t1|0:752566:PASS:253:-253,0,-253:70:25,45",
    f"1\t752721\trs3131972\tA\tG\t.\t.\tNS=1;AN=1;AC=1\t{FMT}\t1|.:752566:PASS:.:-102,0,0:7:6,.",
    f"1\t760000\t.\tC\tT\t.\t.\tNS=1;AN=2;AC=2\t{FMT}\t1/1:.:VQLOW:73:-170,-73,0:5:3,3",
    "1\t770000\t.\tGA\t.\t.\t.\tNS=1;AN=1\tGT:PS\t.|0:752566",
    f"1\t780000\t.\tC\tA,T\t.\t.\tNS=1;AN=2;AC=1,1\t{FMT}\t1|2:752566:PASS:40:-200,-90,-120,-60,0,-40:30:14,16",
    "1\t790000\t.\tC\t.\t.\t.\tNS=1;AN=0\tGT:PS\t.|.:.",
    "X\t2699509\t.\tAAAAGAAGTTAA\t.\t.\t.\tNS=1;AN=0\tGT:PS\t.:.",  # haploid no-call on male X
]


@pytest.fixture
def cgi_vcf(tmp_path: Path) -> Path:
    return write(tmp_path, "\n".join([*HEADER, *RECORDS]) + "\n", "vcfBeta-GS000040194-ASM.vcf")


def test_record_classes(cgi_vcf: Path) -> None:
    assert sniff(cgi_vcf).format == "vcf_cgi"
    loaded = load(cgi_vcf)
    facets, calls, blocks = pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])
    assert blocks.select("chrom", "pos", "end_pos", "kind").rows() == [
        ("1", 1, 10000, "no_call"),  # <CGA_NOCALL> to END
        ("1", 10001, 10006, "no_call"),  # ALT . with ./.: the REF's span
        ("1", 790000, 790000, "no_call"),
        ("X", 2699509, 2699520, "no_call"),  # a haploid `.`, which polars-bio returns as a null GT
    ]
    assert calls["pos"].to_list() == [752566, 752721, 760000, 770000, 780000], "the CNV window is left out"
    assert "1 <CGA_CNVWIN> record(s) left out: copy-number windows, not genotypes" in loaded.warnings
    facet = facets.row(0, named=True)
    assert (facet["format_id"], facet["vendor"], facet["caller"], facet["caller_version"]) == (
        "vcf_cgi", "Complete Genomics", "CGAPipeline", "2.5.0.37",
    )
    assert facet["absence_means"] == "blocks" and facet["ploidy_statements"] == ["##source_MAX_PLOIDY=10"]
    assert "FORMAT/AD" in facet["dropped_fields"]
    assert_gt_consistent(calls)


def test_half_calls_phase_and_filters(cgi_vcf: Path) -> None:
    rows = {row["pos"]: row for row in load(cgi_vcf).calls.collect().iter_rows(named=True)}
    assert rows[752721]["gt_alleles"] == ["G", None] and rows[752721]["phased"] == [True, True]
    assert rows[752721]["phase_sets"] == ["752566", "752566"]
    assert rows[770000]["gt_alleles"] == [None, "GA"] and rows[770000]["alts"] == []
    assert rows[770000]["end_pos"] == 770001, "a half-called reference span stays a call over its REF"
    assert rows[760000]["sample_filters"] == ["VQLOW"] and rows[760000]["filters"] is None
    assert all(row["ad"] is None for row in rows.values()), "CGI's AD counts reads per GT slot"


def test_gl_is_read_in_phred_units(cgi_vcf: Path) -> None:
    loaded = load(cgi_vcf)
    rows = {row["pos"]: row for row in loaded.calls.collect().iter_rows(named=True)}
    assert rows[752566]["gl"] == pytest.approx([-25.3, 0.0, -25.3])
    assert rows[760000]["gl"] == pytest.approx([-17.0, -7.3, 0.0])
    assert rows[780000]["gl"] == pytest.approx([-20.0, -9.0, -12.0, -6.0, 0.0, -4.0])
    for row in rows.values():
        if row["gq"] is not None and row["gl"] is not None:
            best, second = sorted(row["gl"], reverse=True)[:2]
            assert row["gq"] == pytest.approx(10 * (best - second)), "GQ is the gap between the two best, in phred"
    assert any(w.startswith("CGAPipeline writes GL in phred units") for w in loaded.warnings)


def test_a_cgi_vcf_without_the_centre_line_gets_the_same_rules(tmp_path: Path) -> None:
    """Uploads can lose `##center`: the file then sniffs as `vcf`, and the rules still follow the
    declared alleles and the CGAPipeline source line."""
    header = [line for line in HEADER if not line.startswith("##center")]
    path = write(tmp_path, "\n".join([*header, *RECORDS]) + "\n", "upload.vcf")
    assert sniff(path).format == "vcf"
    loaded = load(path)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    assert blocks.height == 4 and calls.height == 5
    assert calls.filter(pl.col("pos") == 760000)["gl"].to_list()[0] == pytest.approx([-17.0, -7.3, 0.0])
    assert loaded.facets.collect()["vendor"][0] is None, "no line names the vendor"


def test_real_cgi_vcf(pgp_participants: Path) -> None:
    """The smallest PGP-HMS vcfBeta on disk, by counts only: no CNV window in calls, every block
    no-called, GQ always the gap between the two best `gl` values."""
    candidates = sorted(pgp_participants.glob("*/arvados-*/vcfBeta-*.vcf.bz2"), key=lambda p: p.stat().st_size)
    if not candidates:
        pytest.skip("no PGP-HMS vcfBeta on disk")
    loaded = load(candidates[0])
    summary = loaded.calls.select(
        pl.len().alias("calls"),
        pl.col("alts").list.contains("<CGA_CNVWIN>").sum().alias("windows"),
        pl.col("alts").list.contains("<CGA_NOCALL>").sum().alias("nocall_calls"),
        ((pl.col("gl").list.max() - pl.col("gl").list.eval(pl.element().sort(descending=True)).list.get(1, null_on_oob=True)) * 10 - pl.col("gq").cast(pl.Float64))
        .abs()
        .gt(0.01)
        .sum()
        .alias("gq_not_gap"),
        (pl.col("end_pos") < pl.col("pos")).sum().alias("end_before_pos"),
    ).collect().row(0, named=True)
    kinds = loaded.blocks.group_by("kind").len().collect()
    assert summary["calls"] > 0 and summary["windows"] == 0 and summary["nocall_calls"] == 0
    assert summary["gq_not_gap"] == 0 and summary["end_before_pos"] == 0
    assert set(kinds["kind"].to_list()) == {"no_call"}
    no_genotype = loaded.calls.filter(pl.col("gt").is_null() & (pl.col("alts").list.len() == 0)).select(pl.len()).collect().item()
    assert no_genotype == 0, "every ALT . record without a called allele is a no-call block"
