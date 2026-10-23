"""Tables that state REF at each site (`just_vcf.loaders.site_tables`), read through `load`.

The synthetic rows copy the shapes counted in the PGP-HMS files (2026-10-23, counts only). CASAVA:
24 tables of one participant, one chromosome each, 3,389,516 rows; Call Type Diff (one letter)
1,372,307, Het1 1,163,955, Het2 852,040, Other 1,214; A+C+G+T Bases equal Bases Used on every row.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load
from just_vcf.formats import get_spec
from just_vcf.loaders.site_tables import FORMATS as SITE_FORMATS
from just_vcf.loading import LOADABLE_FORMATS
from just_vcf.schema import CALLS, canonical_chrom
from just_vcf.sniff import Sniffed, open_payload, sniff
from test_load_vcf import assert_gt_consistent


def test_every_site_table_format_is_loadable() -> None:
    assert SITE_FORMATS <= LOADABLE_FORMATS


@pytest.mark.parametrize("format_id", sorted(SITE_FORMATS))
def test_every_example_loads_with_the_calls_schema(tmp_path: Path, format_id: str) -> None:
    path = tmp_path / "example.txt"
    path.write_text(get_spec(format_id).example)
    loaded = load(path)
    calls = loaded.calls.collect()
    assert calls.schema == CALLS
    assert calls.height > 0
    assert calls["ref"].null_count() == calls.height, "REF is never read from the file (charter 10)"
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["format_id"] == format_id
    again = load(path)
    assert again.calls.collect().equals(calls), "the same file loads to the same rows"


# --- CASAVA ------------------------------------------------------------------------------------------

CASAVA_HEADER = (
    "Index\tSample Name\tChr Id\tPosition\tA Bases\tC Bases\tG Bases\tT Bases\tCall\tTotal Bases\tBases Used\t"
    "Call 1 Score\tCall 2 Score\tRef Base\tCall Type\tdbSNP SNP Name\tGenes\r\n"
)
CASAVA_ROWS = [
    # index, chrom, pos, A, C, G, T, call, total, used, ref, type, dbsnp
    ("1", "1", "752566", 14, 0, 12, 0, "GA", 28, 26, "G", "Het1", "rs3094315"),
    ("2", "1", "752721", 0, 0, 31, 0, "G", 33, 31, "A", "Diff", "rs3131972"),
    ("3", "1", "798959", 9, 0, 0, 11, "TA", 20, 20, "A", "Het2", ""),
    ("4", "1", "800007", 0, 7, 0, 8, "CT", 15, 15, "G", "Other", "rs1"),
    ("5", "M", "150", 0, 0, 0, 900, "T", 905, 900, "C", "Diff", "rs62581312"),
]


def _casava(rows: list[tuple]) -> str:
    lines = [
        "\t".join([index, "Jane Doe", chrom, pos, *map(str, (a, c, g, t)), call, str(total), str(used), "120.5", "98.2", ref, kind, dbsnp, ""])
        for index, chrom, pos, a, c, g, t, call, total, used, ref, kind, dbsnp in rows
    ]
    return CASAVA_HEADER + "\r\n".join(lines) + "\r\n"


def test_casava_calls_index_the_stated_ref_and_carry_base_counts(tmp_path: Path) -> None:
    path = tmp_path / "s_1_snps.txt"
    path.write_bytes(_casava(CASAVA_ROWS).encode())
    loaded = load(path)
    calls = loaded.calls.collect()
    assert_gt_consistent(calls)
    by_pos = {row["pos"]: row for row in calls.iter_rows(named=True)}
    for _, chrom, pos, a, c, g, t, call, total, used, ref, _, dbsnp in CASAVA_ROWS:
        row = by_pos[int(pos)]
        letters = list(call)
        alts = sorted(set(letters) - {ref})
        assert row["chrom"] == canonical_chrom(chrom)
        assert row["source_ref"] == ref and row["ref"] is None and row["ref_check"] == "not_run"
        assert row["gt_alleles"] == letters, "one entry per letter written"
        assert row["alts"] == alts
        assert row["gt"] == [([ref] + alts).index(x) for x in letters]
        counts = {"A": a, "C": c, "G": g, "T": t}
        assert row["source_alleles"] == [ref] + alts
        assert row["ad"] == [counts[x] for x in [ref] + alts]
        assert row["dp"] == total
        assert row["ids"] == ([dbsnp] if dbsnp else None)
        assert row["end_pos"] == row["pos"]
    facet = loaded.facets.collect().row(0, named=True)
    assert "Sample Name" in facet["dropped_fields"] and "Bases Used" in facet["dropped_fields"]
    assert facet["sample"] is None, "the free-text sample name is never read into a frame"
    assert not any("Jane" in w for w in loaded.warnings)
    assert facet["absence_means"] == "unknown"
    assert facet["assembly_label"] is None, "the file states no build"
    assert {"MT"} <= set(calls["chrom"].to_list()), "M is mapped to MT explicitly"


def test_casava_call_type_that_disagrees_is_counted(tmp_path: Path) -> None:
    rows = [*CASAVA_ROWS[:2], ("9", "1", "900000", 5, 0, 5, 0, "GA", 10, 10, "G", "Het2", "")]
    path = tmp_path / "s_1_snps.txt"
    path.write_bytes(_casava(rows).encode())
    loaded = load(path)
    assert any("1 row(s) whose Call Type disagrees" in w for w in loaded.warnings)
    row = loaded.calls.filter(pl.col("pos") == 900000).collect().row(0, named=True)
    assert row["gt"] == [0, 1], "the Call is kept"


def test_casava_unreadable_call_is_kept_unread(tmp_path: Path) -> None:
    rows = [*CASAVA_ROWS[:1], ("9", "1", "900000", 5, 0, 5, 0, "GAT", 10, 10, "G", "Other", "")]
    path = tmp_path / "s_1_snps.txt"
    path.write_bytes(_casava(rows).encode())
    loaded = load(path)
    row = loaded.calls.filter(pl.col("pos") == 900000).collect().row(0, named=True)
    assert row["gt_alleles"] is None and row["gt"] is None and row["alts"] is None
    assert row["source_gt"] == "GAT"
    assert any("1 call(s) not readable" in w for w in loaded.warnings)


INDEL_HEADER = (
    "Index\tSample\tChr\tpos\tCIGAR\tref_upstream\tref/indel\tref_downstream\tQ(indel)\tmax_gtype\tQ(max_gtype)\t"
    "max2_gtype\tbp1_reads\tref_reads\tindel_reads\tother_reads\trepeat_unit\tref_repeat_count\tindel_repeat_count\r\n"
)


def _indel_row(index: int, chrom: str, pos: int, cigar: str, ref_indel: str, q: int, gtype: str, gq: int, ref_reads: int, indel_reads: int) -> str:
    return "\t".join(map(str, [index, "Jane Doe", chrom, pos, cigar, "GTTTACTACA", ref_indel, "CTTTGTAGTA", q, gtype, gq, "het", 33,
                               ref_reads, indel_reads, 1, "G", 1, 0]))


def test_casava_indels_drop_the_padding_and_place_insertion_points(tmp_path: Path) -> None:
    rows = [
        (0, "20", 10000012, "1D", "G/-", 3070, "hom", 76, 0, 23),
        (1, "20", 10000058, "1I", "-/T", 2887, "het", 33, 9, 11),
        (2, "20", 10000100, "2D", "AG/--", 900, "het", 40, 7, 6),
        (3, "M", 300, "1D1M1I", "C-A/-TA", 120, "hom", 20, 0, 9),
    ]
    path = tmp_path / "indels.txt"
    path.write_bytes((INDEL_HEADER + "".join(_indel_row(*row) + "\r\n" for row in rows)).encode())
    loaded = load(path)
    calls = loaded.calls.collect()
    assert_gt_consistent(calls)
    got = {(row["chrom"], row["pos"]): row for row in calls.iter_rows(named=True)}
    for _, chrom, pos, _, ref_indel, q, gtype, gq, ref_reads, indel_reads in rows:
        ref, alt = (side.replace("-", "") for side in ref_indel.split("/"))
        at = pos - 1 if ref == "" else pos
        row = got[(canonical_chrom(chrom), at)]
        assert (row["source_ref"], row["alts"]) == (ref, [alt])
        assert row["end_pos"] == (None if ref == "" else pos + len(ref) - 1)
        assert row["gt"] == {"het": [0, 1], "hom": [1, 1]}[gtype] and row["source_gt"] == gtype
        assert (row["gq"], row["qual"], row["ad"]) == (gq, q, [ref_reads, indel_reads])
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["content_filter"] == "indel_only" and "Sample" in facet["dropped_fields"] and facet["sample"] is None
    assert any("1 insertion point(s)" in w for w in loaded.warnings)


# --- SOAPsnp -----------------------------------------------------------------------------------------

IUPAC = {"A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"], "R": ["A", "G"], "Y": ["C", "T"], "K": ["G", "T"],
         "M": ["A", "C"], "S": ["C", "G"], "W": ["A", "T"]}


def _soap_row(chrom: str, pos: int, ref: str, consensus: str, quality: int, depth: int) -> str:
    return "\t".join(map(str, [chrom, pos, ref, consensus, quality, "G", 36, 12, 12, "A", 30, 10, 10, depth, 0.412, "1.00000", 1, 15882]))


def test_soapsnp_iupac_codes_quality_and_depth(tmp_path: Path) -> None:
    rows = [("chr1", 742429, "G", "R", 34, 22), ("chr1", 758311, "G", "A", 87, 19), ("chrX", 9000, "C", "Y", 41, 0),
            ("chrM", 150, "C", "T", 99, 900), ("chr2", 500, "A", "B", 20, 9)]
    path = tmp_path / "cns.snp"
    path.write_text("".join(_soap_row(*row) + "\n" for row in rows))
    loaded = load(path)
    calls = loaded.calls.collect()
    assert_gt_consistent(calls)
    by_pos = {row["pos"]: row for row in calls.iter_rows(named=True)}
    for chrom, pos, ref, consensus, quality, depth in rows:
        row = by_pos[pos]
        assert row["chrom"] == canonical_chrom(chrom) and row["source_ref"] == ref and row["source_gt"] == consensus
        assert row["gq"] == quality and row["dp"] == depth
        if consensus in IUPAC:
            assert row["gt_alleles"] == IUPAC[consensus]
            assert row["alts"] == sorted(set(IUPAC[consensus]) - {ref})
        else:
            assert row["gt_alleles"] is None and row["gt"] is None, "a three-base code names no genotype"
    assert any("1 row(s) with depth 0" in w for w in loaded.warnings)
    assert any("1 consensus code(s) not readable" in w for w in loaded.warnings)


# --- FamilyTreeDNA Big Y -----------------------------------------------------------------------------


def test_bigy_known_snps_without_positions_are_dropped_and_rows_are_chr_y(tmp_path: Path) -> None:
    text = (
        '"Type","Position","SNPName","Derived","OnTree","Reference","Genotype"\n'
        '"Novel Variant","10061836","","","","A","A"\n'
        '"Known SNP","","A14825","?","Yes","A","?"\n'
        '"Novel Variant","10061838","","","","A","?"\n'
        "Novel Variant,10061900,,,,C,T\n"
    )
    path = tmp_path / "bigy.csv"
    path.write_text(text)
    loaded = load(path)
    calls = loaded.calls.collect()
    assert_gt_consistent(calls)
    assert calls["chrom"].unique().to_list() == ["Y"]
    rows = {row["pos"]: row for row in calls.iter_rows(named=True)}
    assert set(rows) == {10061836, 10061838, 10061900}
    assert rows[10061836]["gt"] == [0] and rows[10061836]["alts"] == []
    assert rows[10061838]["gt_alleles"] == [None] and rows[10061838]["gt"] == [None]
    assert rows[10061900]["gt"] == [1] and rows[10061900]["alts"] == ["T"]
    assert any("1 row(s) without a position" in w for w in loaded.warnings)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["assembly_label"] is None, "Big Y states no build"
    assert facet["contig_map"] == [], "the file names no contig"


# --- Color variant set -------------------------------------------------------------------------------

COLOR_HEAD_NEW = (
    "SampleID,Chromosome,Position,Ref,Alt,Gene,Type,cHGVS,pHGVS,Transcript,BIC,Zygosity,StructuralVariant,"
    "Classification,Comments,TestType,ReferenceVersion,AnalysisCompleted\r\n"
)
COLOR_HEAD_OLD = (
    "SampleID,Chromosome,Pos,Ref,Alt,Gene,Type,cHGVS,pHGVS,Exon,Transcript,Zygosity,isStructuralVariant,"
    "Classification,TestOrdered,ReferenceVersion\r\n"
)


def test_color_variant_set_both_layouts_and_the_footer(tmp_path: Path) -> None:
    new = COLOR_HEAD_NEW + (
        "K1,1,55523033,A,G,PCSK9,SNV,c.1,p.1,ENST1,,Homozygous,,Benign,,T,GRCh37.p12,2025-06-23\r\n"
        "K1,13,32906729,CA,C,BRCA2,indel,c.2,p.2,ENST2,,Heterozygous,,Benign,,T,GRCh37.p12,2025-06-23\r\n"
        "\r\nRaw data disclaimer:\r\n\"This file, which is not a report,\r\nspans lines, with, commas\"\r\n"
    )
    old = COLOR_HEAD_OLD + (
        "K2,1,55523033,A,G,PCSK9,SNV,c.1,p.1,4,ENST1,Homozygous,,Benign,T,GRCh37\r\n"
        "K2,13,32906729,CA,C,BRCA2,indel,c.2,p.2,5,ENST2,Heterozygous,,Benign,T,GRCh37\r\n"
    )
    frames = []
    for name, text in (("new.csv", new), ("old.csv", old)):
        path = tmp_path / name
        path.write_text(text)
        loaded = load(path)
        assert loaded.warnings == (), "the footer is not counted as rows"
        calls = loaded.calls.collect()
        assert_gt_consistent(calls)
        facet = loaded.facets.collect().row(0, named=True)
        assert facet["assembly_label"] == "GRCh37" and facet["absence_means"] == "unknown"
        frames.append(calls.select("chrom", "pos", "end_pos", "source_ref", "alts", "gt", "gt_alleles"))
        assert calls.filter(pl.col("pos") == 32906729)["end_pos"].to_list() == [32906730], "the REF span"
    assert frames[0].equals(frames[1]), "both layouts give the same rows"
    assert frames[0].filter(pl.col("pos") == 55523033)["gt"].to_list() == [[1, 1]]


def test_color_variant_set_unknown_zygosity_keeps_gt_null(tmp_path: Path) -> None:
    path = tmp_path / "v.csv"
    path.write_text(COLOR_HEAD_NEW + "K1,1,55523033,A,G,PCSK9,SNV,c.1,p.1,ENST1,,Hemizygous,,Benign,,T,GRCh37,2025-06-23\r\n")
    loaded = load(path)
    row = loaded.calls.collect().row(0, named=True)
    assert row["gt"] is None and row["gt_alleles"] is None and row["source_gt"] == "Hemizygous"
    assert any("Zygosity other than" in w for w in loaded.warnings)


# --- GET-Evidence GFF --------------------------------------------------------------------------------


def test_get_evidence_gff_ref_rows_are_blocks_and_ref_allele_indexes(tmp_path: Path) -> None:
    text = (
        "##genome-build b37\n# File creation date: 2019-01-15\n"
        "chr1\tCGI\tREF\t752566\t752566\t.\t+\t.\t.\n"
        "chr1\tCGI\tSNP\t752721\t752721\t.\t+\t.\talleles A/G;db_xref dbsnp:rs3131972,dbsnp:rs61770172;ref_allele A\n"
        "chr1\tCGI\tSNP\t776546\t776546\t.\t+\t.\talleles G;db_xref dbsnp:rs12124819;ref_allele A;amino_acid X\n"
        "chr10\tCGI\tSNP\t100\t100\t.\t+\t.\talleles ;ref_allele C\n"
    )
    path = tmp_path / "genome.gff"
    path.write_text(text)
    loaded = load(path)
    calls, blocks = loaded.calls.collect(), loaded.blocks.collect()
    assert_gt_consistent(calls)
    assert blocks.select("chrom", "pos", "end_pos", "kind").rows() == [("1", 752566, 752566, "ref")]
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["absence_means"] == "blocks", "only a listed REF site is reference"
    assert facet["assembly_label"] == "GRCh37"
    assert "amino_acid" in facet["dropped_fields"]
    rows = {row["pos"]: row for row in calls.iter_rows(named=True)}
    assert rows[752721]["ids"] == ["rs3131972", "rs61770172"] and rows[752721]["gt"] == [0, 1]
    assert rows[776546]["gt_alleles"] == ["G"] and rows[776546]["alts"] == ["G"] and rows[776546]["gt"] == [1]
    assert rows[100]["gt_alleles"] == [None], "an empty alleles value is a no-call, not hom-ref"
    assert [row["chrom"] for row in calls.iter_rows(named=True)] == ["1", "1", "10"], "karyotypic order"


def test_get_evidence_gff_without_ref_allele_is_array_shaped(tmp_path: Path) -> None:
    path = tmp_path / "b36.gff"
    path.write_text("##genome-build 36\nchr1\t-\t-\t1000\t1000\t.\t+\t.\talleles A/G\nchrX\t-\t-\t2000\t2000\t.\t+\t.\talleles T\n")
    loaded = load(path)
    calls = loaded.calls.collect()
    assert calls["source_ref"].null_count() == 2 and calls["alts"].null_count() == 2 and calls["gt"].null_count() == 2
    assert calls["gt_alleles"].to_list() == [["A", "G"], ["T"]]
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["absence_means"] == "unknown" and facet["variant_classes"] is None
    assert facet["assembly_label"] == "GRCh36"


# --- Real files (opt-in) -----------------------------------------------------------------------------


@pytest.fixture(scope="session")
def site_uploads(pgp_participants: Path) -> dict[str, Sniffed]:
    """The first participant file of each site-table format, by content (sniff), never by name."""
    found: dict[str, Sniffed] = {}
    for path in sorted(pgp_participants.glob("*/upload-*/*")):
        if len(found) == len(SITE_FORMATS):
            break
        if not path.is_file() or path.stat().st_size > 60_000_000:
            continue
        result = sniff(path)
        if result.format in SITE_FORMATS and result.format not in found:
            found[result.format] = result
    return found


def test_real_casava_table_loads_every_row(site_uploads: dict[str, Sniffed]) -> None:
    sniffed = site_uploads.get("casava_snp_table")
    if sniffed is None:
        pytest.skip("no casava_snp_table upload found")
    with open_payload(sniffed) as handle:
        text = handle.read().decode()
    rows = [line.split("\t") for line in text.replace("\r", "").split("\n")[1:] if line]
    calls = load(sniffed).calls.collect()
    assert calls.height == len(rows)
    assert_gt_consistent(calls)
    assert sorted(calls["source_ref"].to_list()) == sorted(r[13] for r in rows)
    used = {int(r[3]): int(r[10]) for r in rows}
    assert all(row["ad"] is not None and sum(row["ad"]) <= used[row["pos"]] for row in calls.select("pos", "ad").iter_rows(named=True)), (
        "the called alleles' base counts are among the bases used"
    )


def _data_lines(fmt: str, text: str) -> list[str]:
    """The test's own reading of which lines of a real file are data rows."""
    lines = [line for line in text.replace("\r", "").split("\n") if line.strip()]
    if fmt == "get_evidence_gff":
        return [line for line in lines if not line.startswith("#")]
    if fmt == "soapsnp_cns":
        return lines
    body = lines[1:]
    if fmt == "ftdna_bigy":
        return [line for line in body if line.split(",")[1].strip('"')]
    if fmt == "color_variant_set":
        end = next((i for i, line in enumerate(body) if not line.split(",")[0] or line.startswith(("Raw data", '"'))), len(body))
        return body[:end]
    return body


@pytest.mark.parametrize("format_id", sorted(SITE_FORMATS))
def test_real_site_table_loads_every_data_row(site_uploads: dict[str, Sniffed], format_id: str) -> None:
    sniffed = site_uploads.get(format_id)
    if sniffed is None:
        pytest.skip(f"no {format_id} upload found")
    with open_payload(sniffed) as handle:
        text = handle.read().decode("utf-8", errors="replace")
    loaded = load(sniffed)
    calls, blocks = loaded.calls.collect(), loaded.blocks.collect()
    assert calls.schema == CALLS
    assert calls.height + blocks.height == len(_data_lines(format_id, text)), "every data row is one calls or blocks row"
    assert_gt_consistent(calls)
    assert calls["ref"].null_count() == calls.height
