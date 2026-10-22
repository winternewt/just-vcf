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

from just_vcf import Loaded, load, sniff
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


def test_format_2_0_and_no_title(tmp_path: Path) -> None:
    loaded, facet, calls, _ = frames(write(tmp_path, var_text(LOCI, version="2.0", title=False), "old.tsv.bz2"))
    assert facet["callset_type"] is None, "no #TITLE says what the file is"
    assert facet["absence_means"] == "blocks"
    assert by_pos(calls)[("1", 752566)]["sample_filters"] == ["VQHIGH"], "varQuality's vocabulary is kept, not mapped"
    assert facet["dropped_fields"] == ["varScoreVAF", "varScoreEAF", "hapLink"]


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

