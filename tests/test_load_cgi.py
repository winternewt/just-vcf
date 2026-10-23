"""Complete Genomics native callsets (`cgi_var`, `cgi_master_var`): loci into calls and blocks.

The synthetic loci copy shapes counted in a PGP-HMS `var` and `masterVarBeta` of one assembly
(2026-10-22, counts only): a het SNP written as one `snp` and one `ref` haplotype row; a locus whose
haplotypes are split differently (`TAG>CAA` against `TA` + `G>A`); a deletion with an empty allele; a
zero-length insertion; a no-called stretch followed by an insertion on one haplotype (masterVarBeta
writes it `?A`); half-called reference; zero-length no-calls (346,474 in that masterVarBeta); haploid
male X; `PAR-called-in-X` on chrY; format 2.0's `varQuality`.
"""

import bz2
import os
import random
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

from just_vcf import Loaded, callset, load, sniff
from just_vcf.loaders import cgi
from just_vcf.schema import BLOCKS, CALLS, CallsetType, canonical_order
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


def meta(kind: str, version: str = "2.5", title: bool = True, software: str = "2.5.0.48") -> str:
    lines = [
        *(["#TITLE\tWhole Human Genome Sequencing"] if title else []),
        "#ASSEMBLY_ID\tGS000000001-ASM",
        "#GENERATED_BY\tcgatools",
        "#GENOME_REFERENCE\tNCBI build 37",
        f"#SOFTWARE_VERSION\t{software}",
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
    link: str = ""
    """hapLink: the id of the haplotype this row is on; empty when it is linked to nothing."""


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
    Locus("chr1", 810002, 815000, uniform="ref"),
    Locus(
        "chr1", 815000, 815004,
        haps=(
            (Seg(815000, 815004, "no-call", "ACGT", "?", filter=""),),
            (Seg(815000, 815002, "no-call", "AC", "?", filter=""), Seg(815002, 815004, "no-call-ri", "GT", "?", filter="")),
        ),
    ),
    Locus("chr1", 815004, 820000, uniform="ref"),
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
    elif version.startswith("1."):
        columns = [*columns[:9], "totalScore", "hapLink", "xRef"]
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
                head = [str(number), str(locus.ploidy), str(allele), locus.chrom, str(seg.begin), str(seg.end), seg.var_type, seg.reference, seg.seq]
                # 1.x: one totalScore, no filter column.
                row = [*head, score, seg.link, seg.xref] if version.startswith("1.") else [*head, score, score, filt, seg.link, seg.xref]
                lines.append("\t".join(row + [""] * (len(columns) - len(row))))
    software = "1.11.0.16" if version.startswith("1.") else "2.5.0.48"
    return meta("VAR-ANNOTATION", version, title, software) + "\n".join(lines) + "\n"


def joined(segments: tuple[Seg, ...]) -> str:
    """A haplotype as masterVarBeta writes it: the sub-rows' sequences run together, `??` included."""
    return "".join(seg.seq for seg in segments)


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
            links = [next((s.link for s in h if s.link), "") for h in haps] + ([""] if locus.ploidy == 1 else [])
            xrefs = [";".join(s.xref for s in h if s.xref) for h in haps] + ([""] if locus.ploidy == 1 else [])
            row = [str(number), str(locus.ploidy), locus.chrom, str(locus.begin), str(locus.end), "het-ref", "complex", reference,
                   *seqs, "99", "99", "99", "99", *filters, *links, *xrefs, "", "", "", "", "1", "1", "1", "1", locus.total_reads]
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
    assert facet["dropped_fields"] == ["varScoreVAF", "varScoreEAF", "alleleFreq", "alternativeCalls"]
    assert {(m["source"], m["chrom"]) for m in facet["contig_map"]} == {("chr1", "1"), ("chrX", "X"), ("chrY", "Y"), ("chrM", "MT")}
    assert facet["vcf_version"] is None and facet["ploidy_statements"] is None
    assert "1 chrY stretch(es) called on chrX (PAR-called-in-X) kept as no-call blocks" in loaded.warnings


def test_blocks_are_the_uniform_loci_of_positive_length(var) -> None:
    _, _, _, blocks = var
    expected = [
        (locus.chrom.removeprefix("chr").replace("M", "MT"), locus.begin + 1, locus.end, "ref" if locus.uniform == "ref" else "no_call")
        for locus in LOCI
        if locus.uniform and locus.end > locus.begin
    ] + [("1", 810001, 810002, "no_call"), ("1", 815001, 815004, "no_call")]  # loci with every haplotype no-called
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
    assert not {"allele1ReadCount", "allele2ReadCount", "referenceAlleleReadCount", "totalReadCount"} & set(facet["dropped_fields"])
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


def _genome_loci(snvs: int, deletions: int, insertions: int) -> list[Locus]:
    """Het loci 100 bp apart on chr1: SNVs, one-base deletions written `T` to `""` and two-base
    insertions written `""` to `GG`, as CGI writes them, with no anchor base."""
    loci = []
    for i in range(snvs + deletions + insertions):
        b = 1_000_000 + 100 * i
        if i < snvs:
            loci.append(Locus("chr1", b, b + 1, haps=((Seg(b, b + 1, "snp", "A", "G"),), (Seg(b, b + 1, "ref", "A", "A"),))))
        elif i < snvs + deletions:
            loci.append(Locus("chr1", b, b + 1, haps=((Seg(b, b + 1, "ref", "T", "T"),), (Seg(b, b + 1, "del", "T", ""),))))
        else:
            loci.append(Locus("chr1", b, b, haps=((Seg(b, b, "ref", "", ""),), (Seg(b, b, "ins", "", "GG"),))))
    return loci


def test_unanchored_indels_count_so_a_genome_does_not_read_snv_only(tmp_path: Path) -> None:
    """Decision 4: CGI's deletions (`T` to `""`) and insertions (`""` to `GG`) are indels. Counted, the
    rows show more than one class; if the empty allele were no class, 10,000 SNVs alone would read
    `snv_only` and hold back every absence."""
    snvs, deletions, insertions = 10_000, 300, 300
    for text, name in (
        (var_text(_genome_loci(snvs, deletions, insertions)), "var-GS000000002-ASM.tsv.bz2"),
        (master_var_text(_genome_loci(snvs, deletions, insertions)), "masterVarBeta-GS000000002-ASM.tsv.bz2"),
    ):
        _, facet, calls, _ = frames(write(tmp_path, text, name))
        assert calls.filter(pl.col("alts") == [""]).height == deletions and calls.filter(pl.col("source_ref") == "").height == insertions
        assert facet["variant_classes"] == [{"variant_class": "snv", "alleles": snvs}, {"variant_class": "indel", "alleles": deletions + insertions}], name
        assert (facet["content_filter"], facet["content_filter_basis"]) == (None, None), name


def test_the_fixture_facet_counts_its_classes(var) -> None:
    """`LOCI`, by the classes of its called ALTs: partial calls (`?` inside a sequence) are no class."""
    _, facet, calls, _ = var
    alleles = [(r["source_ref"], alt) for r in calls.iter_rows(named=True) for alt in r["alts"] if "?" not in alt]
    found = {"snv": 0, "mnv": 0, "indel": 0}
    for ref, alt in alleles:
        found["indel" if len(ref) != len(alt) else "snv" if len(alt) == 1 else "mnv"] += 1
    assert facet["variant_classes"] == [{"variant_class": cls, "alleles": n} for cls, n in found.items() if n]
    assert facet["content_filter"] is None, "too few alleles to show a filter, and nothing else says one"


def test_without_a_title_the_rows_give_the_callset_type(tmp_path: Path) -> None:
    """RM4 on a native file whose header is silent: the spread of its calls names the class, as for a
    VCF. Absence stays the blocks' to say, whatever the class."""
    loci = _genome_loci(callset.MIN_PANEL_SITES * 2, 0, 0)
    untitled = frames(write(tmp_path, var_text(loci, title=False), "untitled.tsv.bz2"))[1]
    expected = callset.classify(pl.LazyFrame({"chrom": ["1"] * len(loci), "pos": [lo.begin + 1 for lo in loci], "alts": [["G"]] * len(loci)}), ["DP"], None)
    assert expected.callset_type is CallsetType.PANEL, expected.reason
    assert (untitled["callset_type"], untitled["absence_means"]) == ("panel", "blocks")
    titled = frames(write(tmp_path, var_text(loci), "titled.tsv.bz2"))[1]
    assert (titled["callset_type"], titled["absence_means"]) == ("wgs", "blocks"), "a title that says so is kept"


def test_format_2_0_and_no_title(tmp_path: Path) -> None:
    loaded, facet, calls, _ = frames(write(tmp_path, var_text(LOCI, version="2.0", title=False), "old.tsv.bz2"))
    assert facet["callset_type"] is None, "no #TITLE says what the file is, and too few sites for the rows to"
    assert facet["absence_means"] == "blocks"
    assert by_pos(calls)[("1", 752566)]["sample_filters"] == ["VQHIGH"], "varQuality's vocabulary is kept, not mapped"
    assert facet["dropped_fields"] == ["varScoreVAF", "varScoreEAF"]


def test_format_1_x_reads_as_2_5_does_less_the_filters(tmp_path: Path, var) -> None:
    """Format 1.6/1.7 (CGA 1.11, 1.12): twelve columns, one `totalScore` where 2.x has two scores and a
    filter, and xRef may name COSMIC beside dbSNP. The same loci give the same sites, blocks, alleles,
    rs ids and phase sets; `sample_filters` is null, since the file has no filter to say anything."""
    _, _, v25_calls, v25_blocks = var
    cosmic = [
        Locus(lo.chrom, lo.begin, lo.end, lo.ploidy, lo.uniform, tuple(tuple(Seg(**{**s.__dict__, "xref": f"COSMIC:ncv_id:1;{s.xref}" if s.xref else ""}) for s in h) for h in lo.haps), lo.total_reads)
        for lo in LOCI
    ]
    path = write(tmp_path, var_text(cosmic, version="1.6", title=False), "var-GS000000006-ASM.tsv.bz2")
    assert sniff(path).format == "cgi_var"
    loaded, facet, calls, blocks = frames(path)
    assert (facet["caller_version"], facet["dropped_fields"]) == ("1.11.0.16", ["totalScore"])
    assert facet["callset_type"] is None, "no #TITLE, and too few sites for the rows to say"
    assert blocks.drop("facet_id").equals(v25_blocks.drop("facet_id")), "no-call-rc and no-call-ri are no calls, as in 2.5"
    columns = ["chrom", "pos", "end_pos", "source_ref", "alts", "gt", "gt_alleles", "ids", "phased"]
    assert calls.select(columns).equals(v25_calls.select(columns))
    assert calls["sample_filters"].null_count() == calls.height
    assert not [w for w in loaded.warnings if "dropped" in w]
    linked = frames(write(tmp_path, var_text(LINKED, version="1.7"), "var-GS000000007-ASM.tsv.bz2"))[2]
    assert {(r["chrom"], r["pos"], a, ps) for r in linked.iter_rows(named=True) if r["phase_sets"] for a, ps in zip(r["gt_alleles"], r["phase_sets"]) if ps} == linked_slots(LINKED)


# --- masterVarBeta's read counts as `ad` -------------------------------------------------------------


@dataclass(frozen=True)
class Counted:
    """One masterVarBeta variant locus with its read counts, and the `ad` it should give (None:
    withheld, for `reason`). Counts are text as the file writes them; `""` is an empty field."""

    begin: int
    end: int
    zygosity: str
    reference: str
    seqs: tuple[str, str]
    counts: tuple[str, str, str]
    """allele1ReadCount, allele2ReadCount, referenceAlleleReadCount."""
    ad: list[int] | None
    reason: str | None = None
    ploidy: int = 2


COUNTED = [
    Counted(100, 101, "het-ref", "G", ("A", "G"), ("7", "9", "9"), [9, 7]),
    Counted(200, 201, "hom", "G", ("T", "T"), ("12", "12", "0"), [0, 12]),
    # The slots name the ALTs in the other order than `alts` sorts them: the counts follow the alleles.
    Counted(300, 301, "het-alt", "G", ("T", "C"), ("5", "6", "1"), [1, 6, 5]),
    Counted(400, 401, "hap", "G", ("C", ""), ("8", "", "0"), [0, 8], ploidy=1),
    Counted(500, 501, "half", "G", ("A", "?"), ("4", "", "2"), [2, 4]),
    Counted(600, 600, "het-ref", "", ("GG", ""), ("3", "10", "10"), [10, 3]),  # an insertion point
    Counted(700, 701, "het-ref", "T", ("", "T"), ("5", "6", "6"), [6, 5]),  # a deletion
    Counted(800, 801, "half", "G", ("G", "?"), ("", "", ""), None, "a read count the file leaves empty"),
    Counted(850, 850, "no-call", "", ("?", "?"), ("", "", ""), None, "no haplotype called"),  # an uncalled insertion point
    Counted(900, 902, "half", "GA", ("?A", "GA"), ("2", "9", "9"), None, "a partly called allele"),
    Counted(1000, 1001, "hom", "G", ("T", "T"), ("12", "11", "0"), None, "two haplotypes of one allele with different counts"),
    Counted(1100, 1101, "het-ref", "G", ("A", "G"), ("7", "9", "8"), None, "a reference haplotype whose count is not referenceAlleleReadCount"),
]


def counted_master_var(loci: list[Counted]) -> str:
    lines = [">" + "\t".join(MASTER_COLUMNS)]
    for number, locus in enumerate(loci, 1):
        row = [str(number), str(locus.ploidy), "chr1", str(locus.begin), str(locus.end), locus.zygosity, "complex", locus.reference,
               *locus.seqs, *([""] * 4), "PASS", "" if locus.ploidy == 1 else "PASS", *([""] * 9), *locus.counts, "30"]
        lines.append("\t".join(row))
    return meta("VAR-OLPL") + "\n".join(lines) + "\n"


def test_master_var_read_counts_are_ad_over_ref_and_alts(tmp_path: Path) -> None:
    """`ad` is indexed over `source_alleles = [source_ref] + alts`: REF gets referenceAlleleReadCount,
    each ALT the count of the slot that carries it. Where a count is missing, a slot's count
    contradicts another, or an allele is only partly called, `ad` is withheld and the reason counted."""
    loaded, facet, calls, _ = frames(write(tmp_path, counted_master_var(COUNTED), "masterVarBeta-GS000000003-ASM.tsv.bz2"))
    assert not {"allele1ReadCount", "allele2ReadCount", "referenceAlleleReadCount", "totalReadCount"} & set(facet["dropped_fields"])
    rows = by_pos(calls)
    for locus in COUNTED:
        row = rows[("1", locus.begin if locus.begin == locus.end else locus.begin + 1)]
        assert row["ad"] == locus.ad, locus
        if locus.ad is None:
            assert row["source_alleles"] is None, locus
        else:
            assert row["source_alleles"] == [row["source_ref"], *row["alts"]], locus
            assert len(row["ad"]) == len(row["source_alleles"])
        assert row["dp"] == 30
    reasons = {locus.reason for locus in COUNTED if locus.reason}
    assert reasons <= cgi.AD_WITHHELD_REASONS
    assert {w for w in loaded.warnings if "without `ad`" in w} == {f"1 row(s) without `ad`: {reason}" for reason in reasons}
    # var has no read counts: nothing to index, and nothing withheld to report.
    var_loaded, _, var_calls, _ = frames(write(tmp_path, var_text(LOCI), "var-GS000000003-ASM.tsv.bz2"))
    assert var_calls["ad"].null_count() == var_calls.height and var_calls["source_alleles"].null_count() == var_calls.height
    assert not [w for w in var_loaded.warnings if "without `ad`" in w]


# --- hapLink as phase sets ---------------------------------------------------------------------------


def snp_locus(chrom: str, begin: int, ref: str, alt: str, links: tuple[str, str] = ("", ""), alt_slot: int = 0) -> Locus:
    """A het SNV locus whose ALT sits on haplotype `alt_slot + 1`, each haplotype with its hapLink."""
    haps = [(Seg(begin, begin + 1, "ref", ref, ref, link=links[0]),), (Seg(begin, begin + 1, "ref", ref, ref, link=links[1]),)]
    haps[alt_slot] = (Seg(begin, begin + 1, "snp", ref, alt, link=links[alt_slot]),)
    return Locus(chrom, begin, begin + 1, haps=tuple(haps))


# Two molecules of one 2.5 phase block (`Phased_1_1_1` and `_2`), the ALT on either slot; a 2.0-style
# integer id over a haplotype of two sub-rows, and on a haplotype CGI could not call; the same integer
# on chr2, which is another molecule; and a locus linked to nothing.
LINKED = [
    snp_locus("chr1", 9, "A", "G", ("Phased_1_1_1", "Phased_1_1_2")),
    snp_locus("chr1", 19, "C", "T", ("Phased_1_1_2", "Phased_1_1_1"), alt_slot=1),
    Locus("chr1", 29, 32, haps=((Seg(29, 31, "ref", "GA", "GA", link="7"), Seg(31, 32, "snp", "T", "C", link="7")), (Seg(29, 32, "ref", "GAT", "GAT"),))),
    snp_locus("chr1", 39, "A", "G"),
    Locus("chr1", 49, 50, haps=((Seg(49, 50, "no-call", "A", "?", filter="", link="7"),), (Seg(49, 50, "snp", "A", "G", link="8"),))),
    snp_locus("chr2", 99, "A", "G", ("7", "")),
]


def linked_slots(loci: list[Locus]) -> set[tuple]:
    """`(chrom, pos, allele, hapLink)` for every linked haplotype, from the loci as written; a
    haplotype written `?` is the allele None."""
    out = set()
    for locus in loci:
        for segments in locus.haps:
            link = next((s.link for s in segments if s.link), "")
            if link:
                seq = joined(segments)
                out.add((locus.chrom.removeprefix("chr"), locus.begin + 1, None if seq == "?" else seq, link))
    return out


def test_hap_links_are_phase_sets_of_their_own_kind(tmp_path: Path) -> None:
    """A hapLink names a molecule, whatever slot it sits in: both files give the same sets of alleles,
    though masterVarBeta writes the ALT haplotype first. A linked slot is phased with the id as its
    set, kind `cgi_haplink`; an unlinked one says nothing (null), and a row with none has no phase."""
    expected = linked_slots(LINKED)
    for name, text in (("var-GS000000004-ASM.tsv.bz2", var_text(LINKED)), ("masterVarBeta-GS000000004-ASM.tsv.bz2", master_var_text(LINKED))):
        _, facet, calls, _ = frames(write(tmp_path, text, name))
        assert not {"hapLink", "allele1HapLink", "allele2HapLink"} & set(facet["dropped_fields"]), name
        got = {
            (row["chrom"], row["pos"], allele, ps)
            for row in calls.iter_rows(named=True) if row["phase_sets"] is not None
            for allele, ps in zip(row["gt_alleles"], row["phase_sets"]) if ps is not None
        }
        assert got == expected, name
        for row in calls.iter_rows(named=True):
            if row["phase_sets"] is None:
                assert row["phased"] is None and row["phase_set_kinds"] is None, (name, row["pos"])
                continue
            assert row["phased"] == [True if ps is not None else None for ps in row["phase_sets"]], (name, row["pos"])
            assert row["phase_set_kinds"] == ["cgi_haplink" if ps is not None else None for ps in row["phase_sets"]], (name, row["pos"])
        assert by_pos(calls)[("1", 40)]["phased"] is None, "linked to nothing"


def test_a_haplotype_naming_two_hap_links_gets_no_phase(tmp_path: Path) -> None:
    torn = Locus("chr1", 29, 32, haps=((Seg(29, 31, "ref", "GA", "GA", link="7"), Seg(31, 32, "snp", "T", "C", link="9")), (Seg(29, 32, "ref", "GAT", "GAT", link="8"),)))
    loaded, _, calls, _ = frames(write(tmp_path, var_text([torn]), "torn.tsv.bz2"))
    row = by_pos(calls)[("1", 30)]
    assert (row["gt_alleles"], row["phase_sets"], row["phased"]) == (["GAC", "GAT"], [None, "8"], [None, True]), "the locus stays; one slot loses its phase"
    assert "1 haplotype(s) whose rows name more than one hapLink left without a phase set" in loaded.warnings


# --- var against CGI's own VCF of the same assembly --------------------------------------------------


_KEY = pl.Struct({"key_pos": pl.Int64, "key_ref": pl.String, "set": pl.List(pl.String), "set_called": pl.List(pl.String)})


def _key(row: dict) -> dict:
    """A site as both files can be compared: the bases every called allele and REF share at the start
    are cut and the position moves past them, so CGI's VCF anchor base (`CT>C` against var's `T>""`)
    and a base var's locus already shares (`T>TAGAG`) both fall away. `set` is the sorted genotype, a
    no-called haplotype `.`; `set_called` also reads a partly called one (`?A`) as `.`."""
    alleles = row["gt_alleles"] or []
    basis = [row["source_ref"], *(a for a in alleles if a is not None and "?" not in a)]
    cut = len(os.path.commonprefix(basis))
    return {
        "key_pos": row["start"] + cut,
        "key_ref": row["source_ref"][cut:],
        "set": sorted("." if a is None else a[cut:] for a in alleles),
        "set_called": sorted("." if a is None or "?" in a else a[cut:] for a in alleles),
    }


def _sites(calls: pl.LazyFrame, prefix: str) -> pl.LazyFrame:
    """`_key` over every row. The first base of a row's span is `pos`, except at a CGI insertion point,
    which sits at the base before it (`source_ref` empty)."""
    start = pl.col("pos").cast(pl.Int64) + (pl.col("source_ref") == "").cast(pl.Int64)
    keyed = pl.struct(start.alias("start"), "source_ref", "gt_alleles").map_elements(_key, return_dtype=_KEY)
    return calls.select(
        "chrom",
        keyed.alias("_k"),
        pl.col("gt_alleles").list.eval(pl.element().str.contains("?", literal=True)).list.any().fill_null(False).alias(f"{prefix}_partial"),
        ((pl.col("source_ref") == "") & pl.col("gt").list.eval(pl.element().is_null()).list.all()).alias(f"{prefix}_unknown_insertion"),
    ).unnest("_k").rename({"set": f"{prefix}_set", "set_called": f"{prefix}_set_called"})


def _bases(intervals: dict[str, pl.LazyFrame]) -> pl.DataFrame:
    """Bases per combination of interval sets (each `chrom, pos, end_pos`, 1-based inclusive, disjoint
    within a set), by a sweep over their breakpoints."""
    events = pl.concat(
        [
            frame.filter(pl.col("end_pos") >= pl.col("pos")).select(
                "chrom",
                pl.concat_list(pl.col("pos").cast(pl.Int64), pl.col("end_pos").cast(pl.Int64) + 1).alias("at"),
                pl.lit([1, -1]).alias("delta"),
            )
            .explode("at", "delta", empty_as_null=False)  # two entries each, never empty
            .with_columns(pl.lit(name).alias("set"))
            for name, frame in intervals.items()
        ]
    )
    depth = (
        events.group_by("chrom", "at")
        .agg(*[pl.col("delta").filter(pl.col("set") == name).sum().alias(name) for name in intervals])
        .sort("chrom", "at")
        .with_columns(*[pl.col(name).cum_sum().over("chrom") for name in intervals])
        .with_columns((pl.col("at").shift(-1).over("chrom") - pl.col("at")).alias("bases"))
    )
    return (
        depth.filter(pl.col("bases").is_not_null())
        .group_by(*[(pl.col(name) > 0).alias(name) for name in intervals])
        .agg(pl.col("bases").sum())
        .collect()
    )


def compare(var: Loaded, vcf: Loaded) -> dict[str, int]:
    """How a `var` facet and the vcfBeta facet of the same assembly agree, as counts per class."""
    sites = _sites(var.calls, "var").join(_sites(vcf.calls, "vcf"), on=["chrom", "key_pos", "key_ref"], how="full", coalesce=True)
    in_var, in_vcf = pl.col("var_set").is_not_null(), pl.col("vcf_set").is_not_null()
    same = pl.col("var_set") == pl.col("vcf_set")
    partial = ~same & (pl.col("var_set_called") == pl.col("vcf_set_called"))
    var_none_called = pl.col("var_set_called").list.eval(pl.element() == ".").list.all()
    classes = sites.select(
        (in_var & in_vcf & same).sum().alias("same"),
        (in_var & in_vcf & partial).sum().alias("partial_as_no_call"),
        (in_var & in_vcf & ~same & ~partial).sum().alias("genotype_differs"),
        (in_var & ~in_vcf & pl.col("var_unknown_insertion")).sum().alias("unknown_insertion_var_only"),
        (in_var & ~in_vcf & ~pl.col("var_unknown_insertion") & var_none_called).sum().alias("partial_only_var_only"),
        (in_var & ~in_vcf & ~pl.col("var_unknown_insertion") & ~var_none_called).sum().alias("var_only"),
        (~in_var & in_vcf).sum().alias("vcf_only"),
    ).collect().row(0, named=True)
    covered = _bases(
        {
            "var_ref": var.blocks.filter(pl.col("kind") == "ref").select("chrom", "pos", "end_pos"),
            "var_other": pl.concat([var.blocks.filter(pl.col("kind") == "no_call").select("chrom", "pos", "end_pos"), var.calls.select("chrom", "pos", "end_pos")]),
            "vcf": pl.concat([vcf.blocks.select("chrom", "pos", "end_pos"), vcf.calls.select("chrom", "pos", "end_pos")]),
        }
    )

    def bases(**where: bool) -> int:
        rows = covered.filter(pl.all_horizontal(pl.col(k) == v for k, v in where.items()))
        return int(rows["bases"].sum())

    return {
        **classes,
        "bases_var_ref": bases(var_ref=True),
        "bases_var_ref_vcf_covered": bases(var_ref=True, vcf=True),
        "bases_var_other": bases(var_other=True),
        "bases_var_other_vcf_silent": bases(var_other=True, vcf=False),
        "bases_vcf_outside_var": bases(var_ref=False, var_other=False, vcf=True),
    }


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


# --- Real files (opt-in) -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real_assembly(pgp_participants: Path) -> dict[str, Loaded]:
    """The smallest PGP-HMS `var` on disk with the `masterVarBeta` and `vcfBeta` of the same assembly,
    loaded once for the tests below."""
    for var_path in sorted(pgp_participants.glob("*/arvados-*/var-*-ASM.tsv.bz2"), key=lambda p: p.stat().st_size):
        assembly = var_path.name.removeprefix("var-").removesuffix(".tsv.bz2")
        siblings = {
            "master_var": var_path.with_name(f"masterVarBeta-{assembly}.tsv.bz2"),
            "vcf": var_path.with_name(f"vcfBeta-{assembly}.vcf.bz2"),
        }
        if all(path.is_file() for path in siblings.values()):
            return {"var": load(var_path), **{name: load(path) for name, path in siblings.items()}}
    pytest.skip("no PGP-HMS CGI assembly with var, masterVarBeta and vcfBeta on disk")


def test_real_var_and_master_var_agree(real_assembly: dict[str, Loaded]) -> None:
    """Counts only: the two native files of one assembly give the same blocks and the same sites with
    the same allele sets."""
    var, master = real_assembly["var"], real_assembly["master_var"]
    for loaded in (var, master):
        facet = loaded.facets.collect().row(0, named=True)
        assert (facet["caller"], facet["absence_means"]) == ("CGAPipeline", "blocks")
        # A whole genome, unanchored indels counted: it shows no content filter.
        classes = {entry["variant_class"]: entry["alleles"] for entry in facet["variant_classes"]}
        assert classes["indel"] >= 0.05 * sum(classes.values()) and facet["content_filter"] is None
    columns = ["chrom", "pos", "end_pos", "kind"]
    assert var.blocks.select(columns).collect().equals(master.blocks.select(columns).collect())

    def sites(loaded: Loaded) -> pl.LazyFrame:
        return loaded.calls.select(
            "chrom", "pos", "source_ref", "alts",
            pl.col("gt_alleles").list.eval(pl.element().fill_null(".")).list.sort().alias("allele_set"),
        )

    keys = ["chrom", "pos", "source_ref", "alts", "allele_set"]
    only_var = sites(var).join(sites(master), on=keys, how="anti", nulls_equal=True).select(pl.len()).collect().item()
    only_master = sites(master).join(sites(var), on=keys, how="anti", nulls_equal=True).select(pl.len()).collect().item()
    assert (only_var, only_master) == (0, 0)


def test_real_var_converges_with_its_vcf_beta(real_assembly: dict[str, Loaded]) -> None:
    """Counts only: where `var` and CGI's own VCF of one assembly hold a site, the allele sets agree,
    and the disagreements fall in classes the two representations explain.

    On the smallest assembly (2026-10-22): 4,138,184 sites the same; 27,505 where vcfBeta writes a
    partly called haplotype `.`; 10,039 var loci called only partly, which vcfBeta writes as no-calls
    (blocks); 346,474 uncalled insertion points, which vcfBeta writes as a no-call on the anchor base.
    Left: 101 genotypes that differ, 1,119 var-only and 445 vcfBeta-only sites (0.04%), where vcfBeta
    joins neighbouring loci into one record or places a repeat's indel elsewhere. Coverage: vcfBeta
    covers every base var does not call reference, and 904,739 of var's 2.78G reference bases, the
    anchor bases of its indels and no-calls.
    """
    assert real_assembly["vcf"].facets.collect()["absence_means"].to_list() == ["reference"]
    counts = compare(real_assembly["var"], real_assembly["vcf"])
    matched = counts["same"] + counts["partial_as_no_call"] + counts["genotype_differs"]
    assert counts["same"] > 0.99 * matched, counts
    assert counts["partial_as_no_call"] > 0, "vcfBeta writes a partly called haplotype as `.`"
    assert counts["partial_only_var_only"] > 0, "and a locus called only partly as a no-call, which is a block"
    assert counts["unknown_insertion_var_only"] > 0, "vcfBeta writes an uncalled insertion point as a block on its anchor"
    unexplained = counts["genotype_differs"] + counts["var_only"] + counts["vcf_only"]
    assert unexplained < 0.001 * matched, counts
    # vcfBeta writes no reference-called record: what it leaves out is what var calls reference.
    assert counts["bases_var_ref_vcf_covered"] < 0.001 * counts["bases_var_ref"], counts
    assert counts["bases_var_other_vcf_silent"] < 0.001 * counts["bases_var_other"], counts



def _snv_slots(calls: pl.LazyFrame, molecule: list[pl.Expr]) -> pl.LazyFrame:
    """`(chrom, pos, allele)` of every phased slot of a single-base row whose allele sits in one slot
    only (a hom row says nothing about which molecule is which), with the molecule it is on."""
    slots = (
        calls.filter((pl.col("source_ref").str.len_chars() == 1) & pl.col("phased").is_not_null())
        .with_columns(pl.int_ranges(0, pl.col("gt_alleles").list.len(), dtype=pl.UInt32).alias("_h"))
        .explode("_h", empty_as_null=False)
        .with_columns(
            pl.col("gt_alleles").list.get(pl.col("_h")).alias("allele"),
            pl.col("phased").list.get(pl.col("_h")).alias("_phased"),
            pl.col("phase_sets").list.get(pl.col("_h"), null_on_oob=True).alias("_set"),
            pl.col("phase_sets").list.drop_nulls().list.first().alias("_row_set"),
        )
        .filter(pl.col("_phased").fill_null(False) & (pl.col("allele").str.len_chars() == 1))
    )
    once = pl.len().over("chrom", "pos", "allele") == 1
    return slots.filter(once).select("chrom", "pos", "allele", *molecule)


def test_real_hap_links_are_the_molecules_of_vcf_beta_phase_sets(real_assembly: dict[str, Loaded]) -> None:
    """Counts only: CGI's own VCF of the assembly writes the same phasing as PS sets, slot `h` of every
    record of a set one molecule. Matched at single-base sites, each hapLink is one such molecule and
    each molecule one hapLink, in var and in masterVarBeta alike.

    On the smallest assembly (2026-10-23, the raw files): 4,190,488 linked single-base slots match a
    vcfBeta slot, and their 25,010 hapLinks and 25,010 molecules pair one to one. About 2% do not
    match: 90,181 sit where vcfBeta writes no record (it joins the locus with a neighbour), 52,763
    where it writes the other haplotype `.`."""
    vcf = _snv_slots(real_assembly["vcf"].calls, [pl.struct("_row_set", "_h").alias("vcf_molecule")])
    for name in ("var", "master_var"):
        links = _snv_slots(real_assembly[name].calls, [pl.col("_set").alias("link")])
        joined = links.join(vcf, on=["chrom", "pos", "allele"], how="inner")
        counts = pl.collect_all(
            [
                joined.select(pl.len()),
                links.select(pl.len()),
                joined.group_by("chrom", "link").agg(pl.col("vcf_molecule").n_unique().alias("n")).select((pl.col("n") > 1).sum(), pl.len()),
                joined.group_by("chrom", "vcf_molecule").agg(pl.col("link").n_unique().alias("n")).select((pl.col("n") > 1).sum(), pl.len()),
            ]
        )
        matched, total = counts[0].item(), counts[1].item()
        (split_links, n_links), (merged_molecules, n_molecules) = counts[2].row(0), counts[3].row(0)
        assert matched > 0.95 * total, (name, matched, total)
        assert n_links > 0 and split_links == 0, (name, split_links, n_links)
        assert merged_molecules == 0 and n_molecules == n_links, (name, merged_molecules, n_molecules, n_links)


def _format_version(path: Path) -> str | None:
    """A CGI file's `#FORMAT_VERSION`, read from its head; None where the head names none."""
    with bz2.open(path, "rt", errors="replace") as text:
        for _, line in zip(range(40), text):
            if line.startswith("#FORMAT_VERSION\t"):
                return line.split("\t", 1)[1].strip()
            if line.startswith(">"):
                break
    return None


def test_a_real_format_1_x_var_loads(pgp_participants: Path) -> None:
    """Counts only. The PGP-HMS census (2026-10-23) has two 1.x `var` files, formats 1.6 and 1.7 (CGA
    1.11, 1.12), both in GET-Evidence collections; the first in path order is loaded whole. It reads
    as a genome: whole-genome classes, blocks of both kinds, rs ids and hapLink sets, nothing dropped."""
    found = next(
        (path for path in sorted(pgp_participants.glob("*/evidence-*/*.tsv.bz2")) if (_format_version(path) or "").startswith("1.")),
        None,
    )
    if found is None:
        pytest.skip("no CGI format 1.x var file on disk")
    loaded = load(found)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["format_id"] == "cgi_var" and facet["caller_version"].startswith("1.")
    assert facet["dropped_fields"] == ["totalScore"]
    assert facet["absence_means"] == "blocks" and facet["callset_type"] == "wgs"
    classes = {entry["variant_class"]: entry["alleles"] for entry in facet["variant_classes"]}
    assert classes["snv"] > 1_000_000 and classes["indel"] >= 0.05 * sum(classes.values())
    counts = loaded.calls.select(
        pl.len().alias("rows"),
        pl.col("sample_filters").is_not_null().sum().alias("filtered"),
        pl.col("ids").is_not_null().sum().alias("with_ids"),
        pl.col("phase_sets").is_not_null().sum().alias("linked"),
    ).collect().row(0, named=True)
    assert counts["filtered"] == 0, "1.x has no filter column"
    assert counts["with_ids"] > 0.3 * counts["rows"] and counts["linked"] > 0, counts
    kinds = dict(loaded.blocks.group_by("kind").len().collect().iter_rows())
    assert kinds.get("ref", 0) > 0 and kinds.get("no_call", 0) > 0
    assert not [w for w in loaded.warnings if "dropped" in w], loaded.warnings
