"""Complete Genomics native callsets (`cgi_var`, `cgi_master_var`): loci into calls and blocks.

The synthetic loci copy shapes counted in a PGP-HMS `var` and `masterVarBeta` of one assembly
(2026-10-22, counts only): a het SNP written as one `snp` and one `ref` haplotype row; a locus whose
haplotypes are split differently (`TAG>CAA` against `TA` + `G>A`); a deletion with an empty allele; a
zero-length insertion; a no-called stretch followed by an insertion on one haplotype (masterVarBeta
writes it `?A`); half-called reference; zero-length no-calls (346,474 in that masterVarBeta); haploid
male X; `PAR-called-in-X` on chrY; format 2.0's `varQuality`.
"""

import bz2
import random
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff
from just_vcf.schema import BLOCKS, CALLS, canonical_order
from test_load_vcf import assert_gt_consistent

VAR_COLUMNS = (
    "locus", "ploidy", "allele", "chromosome", "begin", "end", "varType", "reference", "alleleSeq",
    "varScoreVAF", "varScoreEAF", "varFilter", "hapLink", "xRef", "alleleFreq", "alternativeCalls",
)
MASTER_COLUMNS = (
    "locus", "ploidy", "chromosome", "begin", "end", "zygosity", "varType", "reference", "allele1Seq", "allele2Seq",
    "allele1VarScoreVAF", "allele2VarScoreVAF", "allele1VarScoreEAF", "allele2VarScoreEAF", "allele1VarFilter",
    "allele2VarFilter", "allele1HapLink", "allele2HapLink", "allele1XRef", "allele2XRef", "allele1Freq", "allele2Freq",
    "allele1AlternativeCalls", "allele2AlternativeCalls", "evidenceIntervalId", "allele1ReadCount", "allele2ReadCount",
    "referenceAlleleReadCount", "totalReadCount",
)


def meta(kind: str, version: str = "2.5", title: bool = True) -> str:
    lines = [
        *(["#TITLE\tWhole Human Genome Sequencing"] if title else []),
        "#ASSEMBLY_ID\tGS000000001-ASM",
        "#GENERATED_BY\tcgatools",
        "#GENOME_REFERENCE\tNCBI build 37",
        "#SOFTWARE_VERSION\t2.5.0.48",
        f"#FORMAT_VERSION\t{version}",
        "#GENERATED_BY\tdbsnptool",
        f"#TYPE\t{kind}",
    ]
    return "\n".join(lines) + "\n\n"


@dataclass(frozen=True)
class Seg:
    """One var row of one haplotype: 0-based half-open span, varType, reference, sequence."""

    begin: int
    end: int
    var_type: str
    reference: str
    seq: str
    filter: str = "PASS"
    xref: str = ""


@dataclass(frozen=True)
class Locus:
    chrom: str
    begin: int
    end: int
    ploidy: int = 2
    uniform: str | None = None
    """varType of a one-row (`allele all`) locus."""
    haps: tuple[tuple[Seg, ...], ...] = ()
    total_reads: str = ""


LOCI = [
    Locus("chr1", 0, 10000, uniform="no-ref"),
    Locus("chr1", 10000, 10468, uniform="no-call"),
    Locus("chr1", 10468, 752565, uniform="ref"),
    Locus("chr1", 752565, 752566, haps=((Seg(752565, 752566, "snp", "G", "A", xref="dbsnp.86:rs3094315"),), (Seg(752565, 752566, "ref", "G", "G"),)), total_reads="26"),
    Locus("chr1", 752566, 760000, uniform="ref"),
    Locus(
        "chr1", 760000, 760003,
        haps=(
            (Seg(760000, 760003, "sub", "TAG", "CAA", xref="dbsnp.89:rs1638317;dbsnp.100:rs2854683"),),
            (Seg(760000, 760002, "ref", "TA", "TA"), Seg(760002, 760003, "snp", "G", "A", xref="dbsnp.100:rs2854683")),
        ),
        total_reads="208",
    ),
    Locus("chr1", 760003, 770000, uniform="ref"),
    Locus("chr1", 770000, 770001, haps=((Seg(770000, 770001, "ref", "T", "T"),), (Seg(770000, 770001, "del", "T", "", filter="VQLOW"),))),
    Locus("chr1", 770001, 780000, uniform="ref"),
    Locus("chr1", 780000, 780000, haps=((Seg(780000, 780000, "ref", "", ""),), (Seg(780000, 780000, "ins", "", "GG"),))),
    Locus("chr1", 780000, 790000, uniform="ref"),
    Locus(
        "chr1", 790000, 790001,
        haps=(
            (Seg(790000, 790001, "snp", "T", "A", filter="VQLOW"),),
            (Seg(790000, 790001, "no-call", "T", "?", filter=""), Seg(790001, 790001, "ins", "", "A")),
        ),
    ),
    Locus("chr1", 790001, 800000, uniform="ref"),
    Locus("chr1", 800000, 800001, haps=((Seg(800000, 800001, "ref", "C", "C"),), (Seg(800000, 800001, "no-call", "C", "?", filter=""),))),
    Locus("chr1", 800001, 810000, uniform="ref"),
    Locus(
        "chr1", 810000, 810002,
        haps=((Seg(810000, 810002, "no-call", "AC", "?", filter=""),), (Seg(810000, 810002, "no-call-rc", "AC", "?", filter=""),)),
    ),
    Locus("chr1", 810002, 820000, uniform="ref"),
    Locus("chr1", 820000, 820000, uniform="no-call"),
    Locus("chr1", 820000, 830000, uniform="ref"),
    Locus("chrX", 2699520, 2699521, ploidy=1, haps=((Seg(2699520, 2699521, "snp", "A", "G"),),)),
    Locus("chrX", 2699521, 2700000, ploidy=1, uniform="ref"),
    Locus("chrY", 10000, 2649520, uniform="PAR-called-in-X"),
    Locus("chrM", 0, 16569, ploidy=1, uniform="ref"),
]


def var_text(loci: list[Locus], version: str = "2.5", title: bool = True) -> str:
    columns = list(VAR_COLUMNS)
    if version == "2.0":
        columns = [*columns[:11], "varQuality", "hapLink", "xRef"]
    lines = [">" + "\t".join(columns)]
    for number, locus in enumerate(loci, 1):
        if locus.uniform:
            seq = "=" if locus.uniform == "ref" else "?"
            row = [str(number), str(locus.ploidy), "all", locus.chrom, str(locus.begin), str(locus.end), locus.uniform, "=", seq]
            lines.append("\t".join(row + [""] * (len(columns) - len(row))))
            continue
        for allele, segments in enumerate(locus.haps, 1):
            for seg in segments:
                score = "" if seg.seq == "?" else "99"
                filt = seg.filter if version == "2.5" else ("VQHIGH" if seg.filter == "PASS" else seg.filter)
                row = [str(number), str(locus.ploidy), str(allele), locus.chrom, str(seg.begin), str(seg.end), seg.var_type,
                       seg.reference, seg.seq, score, score, filt, "", seg.xref]
                lines.append("\t".join(row + [""] * (len(columns) - len(row))))
    return meta("VAR-ANNOTATION", version, title) + "\n".join(lines) + "\n"


def joined(segments: tuple[Seg, ...]) -> str:
    text = "".join(seg.seq for seg in segments)
    while "??" in text:
        text = text.replace("??", "?")
    return text


def master_var_text(loci: list[Locus]) -> str:
    lines = [">" + "\t".join(MASTER_COLUMNS)]
    for number, locus in enumerate(loci, 1):
        if locus.uniform:
            seq = "=" if locus.uniform == "ref" else "?"
            zygosity = ("hom" if locus.ploidy == 2 else "hap") if locus.uniform == "ref" else "no-call"
            var_type = locus.uniform if locus.uniform in ("ref", "no-ref", "PAR-called-in-X") else "complex"
            seqs = [seq, seq if locus.ploidy == 2 else ""]
            row = [str(number), str(locus.ploidy), locus.chrom, str(locus.begin), str(locus.end), zygosity, var_type, "=", *seqs]
        else:
            reference = "".join(seg.reference for seg in locus.haps[0])
            # masterVarBeta tends to write the variant haplotype first: reverse the slots of a het-ref.
            haps = list(locus.haps)[::-1] if len(locus.haps) == 2 and joined(locus.haps[0]) != reference else list(locus.haps)
            seqs = [joined(h) for h in haps] + ([""] if locus.ploidy == 1 else [])
            filters = [";".join(sorted({s.filter for s in h if s.filter})) for h in haps] + ([""] if locus.ploidy == 1 else [])
            xrefs = [";".join(s.xref for s in h if s.xref) for h in haps] + ([""] if locus.ploidy == 1 else [])
            row = [str(number), str(locus.ploidy), locus.chrom, str(locus.begin), str(locus.end), "het-ref", "complex", reference,
                   *seqs, "99", "99", "99", "99", *filters, "", "", *xrefs, "", "", "", "", "1", "1", "1", "1", locus.total_reads]
        lines.append("\t".join(row + [""] * (len(MASTER_COLUMNS) - len(row))))
    return meta("VAR-OLPL") + "\n".join(lines) + "\n"


def write(tmp_path: Path, text: str, name: str) -> Path:
    path = tmp_path / name
    path.write_bytes(bz2.compress(text.encode()))
    return path


def frames(path: Path):
    loaded = load(path)
    facets, calls, blocks = pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])
    return loaded, facets.row(0, named=True), calls, blocks


@pytest.fixture
def var(tmp_path: Path):
    return frames(write(tmp_path, var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"))


def by_pos(calls: pl.DataFrame) -> dict[tuple[str, int], dict]:
    return {(r["chrom"], r["pos"]): r for r in calls.iter_rows(named=True)}


def test_facet(var) -> None:
    loaded, facet, calls, blocks = var
    assert calls.schema == CALLS and blocks.schema == BLOCKS
    assert (facet["format_id"], facet["vendor"], facet["caller"], facet["caller_version"]) == (
        "cgi_var", "Complete Genomics", "CGAPipeline", "2.5.0.48",
    )
    assert facet["sample"] == "GS000000001-ASM", "the assembly id, as CGI's own VCF names its sample column"
    assert (facet["assembly_label"], facet["assembly_label_basis"]) == ("GRCh37", "stated")
    assert (facet["callset_type"], facet["absence_means"]) == ("wgs", "blocks")
    assert facet["dropped_fields"] == ["varScoreVAF", "varScoreEAF", "hapLink", "alleleFreq", "alternativeCalls"]
    assert {(m["source"], m["chrom"]) for m in facet["contig_map"]} == {("chr1", "1"), ("chrX", "X"), ("chrY", "Y"), ("chrM", "MT")}
    assert facet["vcf_version"] is None and facet["ploidy_statements"] is None
    assert "1 chrY stretch(es) called on chrX (PAR-called-in-X) kept as no-call blocks" in loaded.warnings


def test_blocks_are_the_uniform_loci_of_positive_length(var) -> None:
    _, _, _, blocks = var
    expected = [
        (locus.chrom.removeprefix("chr").replace("M", "MT"), locus.begin + 1, locus.end, "ref" if locus.uniform == "ref" else "no_call")
        for locus in LOCI
        if locus.uniform and locus.end > locus.begin
    ] + [("1", 810001, 810002, "no_call")]  # the locus whose haplotypes are all no-called
    rank = {"1": 0, "X": 1, "Y": 2, "MT": 3}
    assert blocks.select("chrom", "pos", "end_pos", "kind").rows() == sorted(expected, key=lambda r: (rank[r[0]], r[1]))


def test_variant_loci(var) -> None:
    _, _, calls, _ = var
    rows = by_pos(calls)
    snp = rows[("1", 752566)]
    assert (snp["end_pos"], snp["source_ref"], snp["alts"], snp["gt"], snp["gt_alleles"]) == (752566, "G", ["A"], [1, 0], ["A", "G"])
    assert snp["ids"] == ["rs3094315"] and snp["sample_filters"] == ["PASS"] and snp["ref"] is None
    assert snp["ref_check"] == "not_run" and snp["phased"] is None and snp["gl"] is None and snp["source_gt"] is None
    sub = rows[("1", 760001)]
    assert (sub["end_pos"], sub["source_ref"], sub["alts"], sub["gt_alleles"]) == (760003, "TAG", ["CAA", "TAA"], ["CAA", "TAA"])
    assert sub["ids"] == ["rs1638317", "rs2854683"], "rs ids of both haplotypes, once each, sorted"
    deletion = rows[("1", 770001)]
    assert (deletion["end_pos"], deletion["source_ref"], deletion["alts"], deletion["gt_alleles"]) == (770001, "T", [""], ["T", ""])
    assert deletion["sample_filters"] == ["VQLOW"], "PASS is left out beside another filter"
    insertion = rows[("1", 780000)]
    assert (insertion["end_pos"], insertion["source_ref"], insertion["alts"], insertion["gt"]) == (None, "", ["GG"], [0, 1])
    partial = rows[("1", 790001)]
    assert (partial["source_ref"], partial["alts"], partial["gt_alleles"]) == ("T", ["?A", "A"], ["A", "?A"])
    half = rows[("1", 800001)]
    assert (half["alts"], half["gt"], half["gt_alleles"]) == ([], [0, None], ["C", None])
    unknown_insertion = rows[("1", 820000)]
    assert (unknown_insertion["source_ref"], unknown_insertion["alts"], unknown_insertion["gt"]) == ("", [], [None, None])
    assert unknown_insertion["end_pos"] is None, "an insertion point CGI could not call stays a row"
    haploid = rows[("X", 2699521)]
    assert (haploid["gt"], haploid["gt_alleles"]) == ([1], ["G"])
    assert ("1", 810001) not in rows, "a locus with every haplotype no-called is a block"
    assert calls.height == 8
    assert_gt_consistent(calls)


def test_canonical_order_and_determinism(tmp_path: Path) -> None:
    first = load(write(tmp_path, var_text(LOCI), "a.tsv.bz2"))
    again = load(write(tmp_path, var_text(LOCI), "b.tsv.bz2"))
    calls = first.calls.collect()
    assert calls.equals(canonical_order(first.calls).collect())
    assert first.blocks.collect().equals(canonical_order(first.blocks).collect())
    assert calls.equals(again.calls.collect()), "the same bytes give the same rows and facet id"


def test_shuffled_loci_give_the_same_rows(tmp_path: Path) -> None:
    text = var_text(LOCI)
    head, body = text.split(">", 1)
    column_line, *rows = body.splitlines()
    random.Random(3).shuffle(rows)
    shuffled = head + ">" + column_line + "\n" + "\n".join(rows) + "\n"
    a = load(write(tmp_path, text, "a.tsv.bz2"))
    b = load(write(tmp_path, shuffled, "b.tsv.bz2"))
    assert a.calls.drop("facet_id").collect().equals(b.calls.drop("facet_id").collect())
    assert a.blocks.drop("facet_id").collect().equals(b.blocks.drop("facet_id").collect())


def test_master_var_gives_the_same_sites_blocks_and_allele_sets(tmp_path: Path, var) -> None:
    _, _, var_calls, var_blocks = var
    path = write(tmp_path, master_var_text(LOCI), "masterVarBeta-GS000000001-ASM.tsv.bz2")
    assert sniff(path).format == "cgi_master_var"
    loaded, facet, calls, blocks = frames(path)
    assert facet["caller_version"] == "2.5.0.48" and facet["sample"] == "GS000000001-ASM"
    assert "totalReadCount" not in facet["dropped_fields"] and "allele1ReadCount" in facet["dropped_fields"]
    assert blocks.drop("facet_id").equals(var_blocks.drop("facet_id"))

    def sites(frame: pl.DataFrame) -> pl.DataFrame:
        return frame.select(
            "chrom", "pos", "end_pos", "source_ref", "alts", "ids", "sample_filters",
            pl.col("gt_alleles").list.sort(nulls_last=True).alias("allele_set"),
        )

    assert sites(calls).equals(sites(var_calls))
    slots = calls.join(var_calls, on=["chrom", "pos"], suffix="_var")
    assert (slots["gt_alleles"] != slots["gt_alleles_var"]).any(), "the slots do come in another order"
    rows = by_pos(calls)
    assert rows[("1", 752566)]["dp"] == 26 and rows[("1", 760001)]["dp"] == 208, "totalReadCount"
    assert all(r["dp"] is None for r in var_calls.iter_rows(named=True)), "var has no read counts"
    assert_gt_consistent(calls)


def test_format_2_0_and_no_title(tmp_path: Path) -> None:
    loaded, facet, calls, _ = frames(write(tmp_path, var_text(LOCI, version="2.0", title=False), "old.tsv.bz2"))
    assert facet["callset_type"] is None, "no #TITLE says what the file is"
    assert facet["absence_means"] == "blocks"
    assert by_pos(calls)[("1", 752566)]["sample_filters"] == ["VQHIGH"], "varQuality's vocabulary is kept, not mapped"
    assert facet["dropped_fields"] == ["varScoreVAF", "varScoreEAF", "hapLink"]


def test_bad_rows_and_loci_are_counted(tmp_path: Path) -> None:
    gap = Locus("chr1", 900000, 900003, haps=((Seg(900000, 900001, "ref", "A", "A"), Seg(900002, 900003, "snp", "C", "T")), (Seg(900000, 900003, "ref", "AGC", "AGC"),)))
    expanded = Locus("chr1", 910000, 910001, uniform="snp")
    text = var_text([*LOCI, gap]).rstrip("\n") + "\n9x\t2\tall\tchr1\tbegin\t5\tref\t=\t=\n"
    loaded, _, calls, _ = frames(write(tmp_path, text, "bad.tsv.bz2"))
    assert "1 row(s) whose locus, begin or end is not a number, or whose end lies before its begin, dropped" in loaded.warnings
    assert "1 locus(es) with haplotype rows that do not tile the locus dropped" in loaded.warnings
    assert ("1", 900001) not in by_pos(calls)
    # An `all` row with a variant varType stands for every haplotype.
    row = var_text([expanded]).splitlines()[-1].split("\t")
    row[7], row[8] = "C", "T"
    _, _, calls, _ = frames(write(tmp_path, var_text([]) + "\t".join(row) + "\n", "all.tsv.bz2"))
    assert by_pos(calls)[("1", 910001)]["gt_alleles"] == ["T", "T"]

