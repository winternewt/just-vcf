"""The vendor-table loader (`just_vcf.loaders.table`): AncestryDNA, FTDNA, the RSID/RESULT table, Color,
deCODEme, Genos and Illumina Final Reports, read through `load`."""

import random
import zipfile
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

from just_vcf import MultiSampleError, TableLayoutError, UnsupportedFormatError, load, table_samples
from just_vcf.formats import get_spec
from just_vcf.loaders.table import FORMATS as TABLE_FORMATS
from just_vcf.loading import LOADABLE_FORMATS, Loaded
from just_vcf.schema import CALLS, FACETS, canonical_chrom
from just_vcf.sniff import DATA_SHAPE_LINES, Sniffed, open_payload, sniff


@dataclass(frozen=True)
class Layout:
    """The test's own reading of a format's data rows: field indices, the position offset, the split."""

    rsid: int
    chrom: int
    pos: int
    genotype: tuple[int, ...]
    offset: int = 0
    delimiter: str = "\t"
    skip_to: str | None = None
    """Data rows start after the line that starts with this text, and after one column line."""
    header: bool = True


LAYOUTS: dict[str, Layout] = {
    "ancestrydna": Layout(rsid=0, chrom=1, pos=2, genotype=(3, 4)),
    "ftdna": Layout(rsid=0, chrom=1, pos=2, genotype=(3,), delimiter=","),
    "ftdna_space": Layout(rsid=0, chrom=1, pos=2, genotype=(3,), delimiter=" "),
    "rsid_result_csv": Layout(rsid=0, chrom=1, pos=2, genotype=(3,), delimiter=","),
    "color_genotypes": Layout(rsid=3, chrom=1, pos=2, genotype=(4,), offset=1, delimiter=","),
    "decodeme": Layout(rsid=0, chrom=2, pos=3, genotype=(5,), delimiter=","),
    "genos": Layout(rsid=3, chrom=0, pos=1, genotype=(4,), offset=1, header=False),
    "illumina_final_report": Layout(rsid=0, chrom=1, pos=2, genotype=(3,), skip_to="[Data]"),
}


def _expected_rows(text: str, layout: Layout) -> list[dict]:
    """Data rows of `text` as the test reads them: id, canonical chrom, 1-based pos, genotype text."""
    lines = [line for line in text.replace("\r\n", "\n").split("\n") if line.strip()]
    if layout.skip_to:
        lines = lines[next(i for i, line in enumerate(lines) if line.startswith(layout.skip_to)) + 1 :]
    lines = [line for line in lines if not line.startswith("#")]
    if layout.header:
        lines = lines[1:]  # the column line
    rows = []
    for line in lines:
        fields = [field.strip('"') for field in line.split(layout.delimiter)]
        if fields[layout.chrom] == "0" or not fields[layout.pos].isdigit():
            continue  # unplaced probes, and a column line repeated where two exports were joined
        begin = int(fields[layout.pos])
        offset = layout.offset
        if layout.header is False and int(fields[2]) == begin:
            offset = 0  # a Genos insertion site sits at the base before it
        rows.append(
            {
                "id": fields[layout.rsid],
                "chrom": canonical_chrom(fields[layout.chrom]),
                "pos": begin + offset,
                "source_gt": "".join(fields[i] for i in layout.genotype),
            }
        )
    return rows


def _load_text(tmp_path: Path, text: str | bytes, name: str = "genome.txt", sample: str | None = None) -> Loaded:
    path = tmp_path / name
    path.write_bytes(text if isinstance(text, bytes) else text.encode())
    return load(path, sample=sample)


def _by_id(loaded: Loaded) -> dict[str, dict]:
    return {row["ids"][0]: row for row in loaded.calls.collect().iter_rows(named=True)}


def test_every_table_format_is_loadable() -> None:
    assert TABLE_FORMATS <= LOADABLE_FORMATS
    assert set(LAYOUTS) == set(TABLE_FORMATS)


@pytest.mark.parametrize("format_id", sorted(TABLE_FORMATS))
def test_every_example_round_trips(tmp_path: Path, format_id: str) -> None:
    spec = get_spec(format_id)
    expected = _expected_rows(spec.example, LAYOUTS[format_id])
    loaded = _load_text(tmp_path, spec.example)
    calls = loaded.calls.collect()
    facet = loaded.facets.collect()

    assert calls.schema == CALLS
    assert facet.schema == FACETS
    assert facet["format_id"].to_list() == [format_id]
    assert facet["vendor"].to_list() == [None if spec.vendor == "unknown" else spec.vendor]
    assert facet["absence_means"].to_list() == ["unknown"], "an absent site in a genotype table is unknown"
    got = sorted(zip(calls["ids"].list.first(), calls["chrom"], calls["pos"], calls["source_gt"]))
    want = sorted((r["id"], r["chrom"], r["pos"], r["source_gt"]) for r in expected)
    assert got == want
    no_call_rows = {r["id"] for r in expected if r["source_gt"] in {*spec.no_call, "00", "--", "-"}}
    all_null = calls.filter(pl.col("gt_alleles").list.eval(pl.element().is_null()).list.all())
    assert set(all_null["ids"].list.first().to_list()) == no_call_rows
    assert calls["ref"].null_count() == calls.height, "a genotype table writes no REF; REF comes from the sequence"
    assert calls["source_ref"].null_count() == calls.height
    assert calls["phased"].null_count() == calls.height, "a genotype table says nothing about phase"
    assert calls["evidence_basis"].to_list() == ["assay"] * calls.height


def test_load_is_deterministic_and_order_free(tmp_path: Path) -> None:
    header = "#AncestryDNA raw data download\nrsid\tchromosome\tposition\tallele1\tallele2\n"
    rows = [f"rs{i}\t{random.Random(i).choice(['1', '2', '10', '23', '26'])}\t{1000 + i}\tA\tG\n" for i in range(200)]
    first = _load_text(tmp_path, header + "".join(rows), "one.txt")
    again = load(tmp_path / "one.txt")
    assert first.calls.collect().equals(again.calls.collect())
    assert first.facets.collect().equals(again.facets.collect())
    shuffled = rows[:]
    random.Random(7).shuffle(shuffled)
    other = _load_text(tmp_path, header + "".join(shuffled), "two.txt")
    assert other.calls.collect().drop("facet_id").equals(first.calls.collect().drop("facet_id"))


# --- AncestryDNA -------------------------------------------------------------------------------------

ANCESTRY_HEAD = (
    "#AncestryDNA raw data download\n"
    "#Data was collected using AncestryDNA array version: V2.0\n"
    "#of the SNP using human reference build 37.1 coordinates.\n"
    "rsid\tchromosome\tposition\tallele1\tallele2\n"
)


def test_ancestrydna_haploid_letters_stay_doubled_and_codes_map(tmp_path: Path) -> None:
    body = (
        "rs1\t23\t100\tA\tA\n"
        "rs2\t24\t200\tC\tC\n"
        "rs3\t25\t300\tA\tG\n"
        "rs4\t26\t400\tT\tT\n"
        "rs5\t1\t500\tA\t0\n"
        "rs6\t1\t600\tD\tI\n"
        "rs7\t1\t700\t0\t0\n"
    )
    loaded = _load_text(tmp_path, ANCESTRY_HEAD + body)
    by_id = _by_id(loaded)
    assert {k: (v["chrom"], v["gt_alleles"]) for k, v in by_id.items()} == {
        "rs1": ("X", ["A", "A"]),
        "rs2": ("Y", ["C", "C"]),
        "rs3": ("XY", ["A", "G"]),
        "rs4": ("MT", ["T", "T"]),
        "rs5": ("1", ["A", None]),
        "rs6": ("1", ["D", "I"]),
        "rs7": ("1", [None, None]),
    }
    assert [by_id[k]["source_gt"] for k in ("rs1", "rs5", "rs7")] == ["AA", "A0", "00"]
    assert by_id["rs6"]["end_pos"] is None, "an I/D code's REF span is unknown until the reference is read"
    assert by_id["rs5"]["end_pos"] == 500
    mapping = loaded.facets.collect()["contig_map"][0].to_list()
    assert {m["source"]: m["chrom"] for m in mapping} == {"1": "1", "23": "X", "24": "Y", "25": "XY", "26": "MT"}
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["assembly_label"], facet["assembly_label_basis"], facet["callset_type"]) == ("GRCh37", "stated", "array")
    assert loaded.warnings == ()


def test_ancestrydna_mixed_line_endings_and_bad_rows(tmp_path: Path) -> None:
    body = "rs1\t1\t100\tA\tG\r\nrs2\t1\t200\tA\r\nrs3\t1\tx\tA\tG\r\nrs4\t1\t400\tN\tA\r\n"
    loaded = _load_text(tmp_path, ANCESTRY_HEAD + body)
    by_id = _by_id(loaded)
    assert set(by_id) == {"rs1", "rs4"}
    assert by_id["rs1"]["gt_alleles"] == ["A", "G"]
    assert by_id["rs4"]["gt_alleles"] is None and by_id["rs4"]["source_gt"] == "NA"
    assert len(loaded.warnings) == 3
    assert any(w.startswith("1 row(s) without exactly 5 fields") for w in loaded.warnings)
    assert any(w.startswith("1 row(s) with a position") for w in loaded.warnings)
    assert any(w.startswith("1 genotype(s) not readable") for w in loaded.warnings)


# --- FTDNA and the RSID/RESULT table -----------------------------------------------------------------


def test_ftdna_half_calls_no_calls_and_hemizygous_letters(tmp_path: Path) -> None:
    body = (
        '"rs1","1","100","-A"\n'
        '"rs2","1","200","A-"\n'
        '"rs3","1","300","---"\n'
        '"rs4","1","400","00"\n'
        '"rs5","X","500","G"\n'
        "rs6,1,600,CT\n"
    )
    loaded = _load_text(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\n" + body, "ftdna.csv")
    assert sniff(tmp_path / "ftdna.csv").format == "ftdna"
    by_id = _by_id(loaded)
    assert {k: v["gt_alleles"] for k, v in by_id.items()} == {
        "rs1": [None, "A"],
        "rs2": ["A", None],
        "rs3": [None, None],
        "rs4": [None, None],
        "rs5": ["G"],
        "rs6": ["C", "T"],
    }
    assert by_id["rs3"]["source_gt"] == "---"
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["assembly_label"] is None, "FTDNA states no build; it is never read off positions"
    assert loaded.warnings == ()


def test_ftdna_files_joined_with_their_column_line_repeated(tmp_path: Path) -> None:
    """Seven PGP-HMS FTDNA uploads are an autosomal and an X export joined, the column line repeated."""
    text = 'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\nRSID,CHROMOSOME,POSITION,RESULT\n"rs2","X","200","AA"\n'
    loaded = _load_text(tmp_path, text, "joined.csv")
    assert set(_by_id(loaded)) == {"rs1", "rs2"}
    assert loaded.warnings == ("1 repeated column line(s) inside the table skipped (the file joins several exports)",)


def test_ftdna_space_export_with_a_quoted_comma_tail(tmp_path: Path) -> None:
    filler = "".join(f"rs{1000 + i} 1 {1000 + i} AG\r\n" for i in range(DATA_SHAPE_LINES))
    text = (
        "RSID,CHROMOSOME,POSITION,RESULT\r\n"
        "rs1 1 100 AG\r\n"
        + filler
        + "rs2 22 200 -C\r\n"
        "rs3 22 300 --\r\n"
        '"rs4","22","400","TT"\r\n'
        '"rs5","22","500","CG"\r\n'
    )
    loaded = _load_text(tmp_path, text, "space.csv")
    assert sniff(tmp_path / "space.csv").format == "ftdna_space"
    by_id = {k: v for k, v in _by_id(loaded).items() if int(k[2:]) < 1000}
    assert len(_by_id(loaded)) == 5 + DATA_SHAPE_LINES
    assert {k: (v["chrom"], v["pos"], v["gt_alleles"]) for k, v in by_id.items()} == {
        "rs1": ("1", 100, ["A", "G"]),
        "rs2": ("22", 200, [None, "C"]),
        "rs3": ("22", 300, [None, None]),
        "rs4": ("22", 400, ["T", "T"]),
        "rs5": ("22", 500, ["C", "G"]),
    }
    assert loaded.warnings == ()


def test_rsid_result_table_with_a_preamble_and_unplaced_probes(tmp_path: Path) -> None:
    text = (
        "# A raw data export\n"
        "# genotypes on the forward strand\n"
        "RSID,CHROMOSOME,POSITION,RESULT\n"
        '"rs1","1","100","AA"\n'
        '"GSA-rs2","2","200","DI"\n'
        '"kgp3","0","0","--"\n'
        '"kgp4","0","0","AG"\n'
        '"rs5","XY","500","CT"\n'
    )
    loaded = _load_text(tmp_path, text, "export.csv")
    assert sniff(tmp_path / "export.csv").format == "rsid_result_csv"
    by_id = _by_id(loaded)
    assert set(by_id) == {"rs1", "GSA-rs2", "rs5"}
    assert by_id["GSA-rs2"]["gt_alleles"] == ["D", "I"] and by_id["GSA-rs2"]["end_pos"] is None
    assert by_id["rs5"]["chrom"] == "XY"
    assert loaded.warnings == ("2 row(s) on chromosome 0 (probes with no mapped position) dropped",)
    assert loaded.facets.collect()["vendor"][0] is None, "the vendor is not claimed for this layout"


def test_a_table_whose_every_row_is_dropped_loads_empty(tmp_path: Path) -> None:
    loaded = _load_text(tmp_path, 'RSID,CHROMOSOME,POSITION,RESULT\n"kgp1","0","0","--"\n', "unplaced.csv")
    calls = loaded.calls.collect()
    assert calls.height == 0 and calls.schema == CALLS
    facet = loaded.facets.collect()
    assert facet.height == 1 and facet["contig_map"][0].to_list() == []
    assert loaded.warnings == ("1 row(s) on chromosome 0 (probes with no mapped position) dropped",)


# --- Color -------------------------------------------------------------------------------------------

COLOR_HEADER = "SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion\n"
COLOR_TRAILER = 'Raw data disclaimer:,,,,,\n"This file, which has commas, is not data, and it ends here.",,,,,\n'


def test_color_positions_are_shifted_to_one_based_and_the_trailer_is_cut(tmp_path: Path) -> None:
    body = "K1,2,136608645,rs4988235,AG,GRCh37\nK1,12,112241765,rs671,GG,GRCh37\n"
    loaded = _load_text(tmp_path, COLOR_HEADER + body + COLOR_TRAILER, "discovery_genotypes.csv")
    by_id = _by_id(loaded)
    assert {k: v["pos"] for k, v in by_id.items()} == {"rs4988235": 136608646, "rs671": 112241766}
    assert all(v["end_pos"] == v["pos"] for v in by_id.values())
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["assembly_label"], facet["assembly_label_basis"]) == ("GRCh37", "stated")
    assert facet["sample"] == "K1"
    assert facet["callset_type"] is None, "Color does not say which assay made these calls"
    assert loaded.warnings == (), "the disclaimer after the table is not a malformed row"


def test_color_build_is_read_from_the_rows(tmp_path: Path) -> None:
    same = _load_text(tmp_path, COLOR_HEADER + "K1,2,100,rs1,AG,GRCh38\n", "a.csv")
    assert same.facets.collect()["assembly_label"][0] == "GRCh38"
    mixed = _load_text(tmp_path, COLOR_HEADER + "K1,2,100,rs1,AG,GRCh37\nK1,2,200,rs2,AG,GRCh38\n", "b.csv")
    assert mixed.facets.collect()["assembly_label"][0] is None
    assert mixed.warnings == ("rows state 2 different builds (GRCh37, GRCh38); assembly left unknown",)


def test_a_table_with_two_sample_ids_needs_a_choice(tmp_path: Path) -> None:
    text = COLOR_HEADER + "K1,2,100,rs1,AG,GRCh37\nK2,2,100,rs1,GG,GRCh37\n"
    with pytest.raises(MultiSampleError) as refused:
        _load_text(tmp_path, text, "two.csv")
    assert "K1" not in str(refused.value) and "K2" not in str(refused.value), "sample ids can identify a person"
    chosen = load(tmp_path / "two.csv", sample="K2")
    assert chosen.calls.collect()["gt_alleles"].to_list() == [["G", "G"]]
    assert chosen.facets.collect()["sample"][0] == "K2"
    other = load(tmp_path / "two.csv", sample="K1")
    assert other.facets.collect()["facet_id"][0] != chosen.facets.collect()["facet_id"][0]
    with pytest.raises(MultiSampleError):
        load(tmp_path / "two.csv", sample="K3")


# --- deCODEme ----------------------------------------------------------------------------------------

DECODEME_HEADER = "Name,Variation,Chromosome,Position,Strand,YourCode\n"


def test_decodeme_minus_strand_rows_are_complemented(tmp_path: Path) -> None:
    body = (
        "rs1,A/G,1,100,+,AG\n"
        "rs2,C/T,1,200,-,CC\n"
        "rs3,C/T,1,300,-,CT\n"
        "rs4,C/T,1,400,-,--\n"
        "rs5,D/I,1,500,-,DI\n"
        "rs6,A/G,1,600,?,AG\n"
    )
    loaded = _load_text(tmp_path, DECODEME_HEADER + body, "genotypes.csv")
    by_id = _by_id(loaded)
    assert {k: v["gt_alleles"] for k, v in by_id.items()} == {
        "rs1": ["A", "G"],
        "rs2": ["G", "G"],
        "rs3": ["G", "A"],
        "rs4": [None, None],
        "rs5": ["D", "I"],
        "rs6": None,
    }
    assert by_id["rs2"]["source_gt"] == "CC", "source_gt keeps the letters as written"
    assert by_id["rs6"]["source_gt"] == "AG"
    assert loaded.warnings == ("1 row(s) with a strand other than + or - kept with gt_alleles null",)
    assert loaded.facets.collect()["dropped_fields"][0].to_list() == ["Variation"]


def test_decodeme_build_comes_from_the_info_member_of_its_zip(tmp_path: Path) -> None:
    path = tmp_path / "decodeme.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("genotypes.csv", get_spec("decodeme").example)
        zf.writestr("info.txt", get_spec("decodeme_info").example)
    sniffed = sniff(path)
    assert sniffed.format == "decodeme"
    facet = load(sniffed).facets.collect().row(0, named=True)
    assert (facet["assembly_label"], facet["assembly_label_basis"]) == ("GRCh36", "stated")
    assert facet["member"] == "genotypes.csv"

    alone = tmp_path / "genotypes.csv"
    alone.write_text(get_spec("decodeme").example)
    assert load(alone).facets.collect()["assembly_label"][0] is None, "no info member, no build"


# --- Genos -------------------------------------------------------------------------------------------

GENOS_HEAD = (
    "##FileFormat=Genos\n"
    "##ReferenceGenome=UCSC.hg19\n"
    "##Sample=KIT1\n"
    "##Columns=CHROM\tBEGIN\tEND\tID\tGENOTYPE\n"
)


def test_genos_spans_insertions_and_hemizygous_calls(tmp_path: Path) -> None:
    body = (
        "chr1\t99\t100\trs1\tC/T\n"
        "chr1\t199\t201\trs2\tAT/GC\n"
        "chr1\t299\t300\trs3\tT/-\n"
        "chr1\t400\t400\trs4\t-/A\n"
        "chrX\t499\t500\trs5\tG\n"
        "chrM\t599\t600\trs6\tA/A\n"
        "chr1\t699\t700\trs7\tA/?\n"
    )
    loaded = _load_text(tmp_path, GENOS_HEAD + body, "genos.tsv")
    by_id = _by_id(loaded)
    assert {k: (v["chrom"], v["pos"], v["end_pos"], v["gt_alleles"]) for k, v in by_id.items()} == {
        "rs1": ("1", 100, 100, ["C", "T"]),
        "rs2": ("1", 200, 201, ["AT", "GC"]),
        "rs3": ("1", 300, 300, ["T", "-"]),
        "rs4": ("1", 400, None, ["-", "A"]),
        "rs5": ("X", 500, 500, ["G"]),
        "rs6": ("MT", 600, 600, ["A", "A"]),
        "rs7": ("1", 700, 700, None),
    }
    assert by_id["rs4"]["source_gt"] == "-/A"
    assert len(loaded.warnings) == 2
    assert any(w.startswith("1 insertion site(s) (begin == end)") for w in loaded.warnings)
    assert any(w.startswith("1 genotype(s) not readable") for w in loaded.warnings)
    facet = loaded.facets.collect().row(0, named=True)
    assert (facet["sample"], facet["assembly_label"], facet["callset_type"]) == ("KIT1", "GRCh37", None)


def test_genos_with_other_columns_is_refused(tmp_path: Path) -> None:
    head = GENOS_HEAD.replace("CHROM\tBEGIN\tEND\tID\tGENOTYPE", "CHROM\tSTART\tSTOP\tID\tGENOTYPE")
    with pytest.raises(TableLayoutError, match="begin, end"):
        _load_text(tmp_path, head + "chr1\t99\t100\trs1\tC/T\n", "genos.tsv")


# --- Illumina Final Report ---------------------------------------------------------------------------

REPORT_HEAD = "[Header]\nGSGT Version\t2.0.4\nContent\tGSA-24v2-0_A1.bpm\nNum Samples\t2\n[Data]\n"


def test_final_report_tab_and_comma_exports_load_alike(tmp_path: Path) -> None:
    spec = get_spec("illumina_final_report")
    tab = _load_text(tmp_path, spec.example, "tab.txt")
    comma_text = spec.example.replace("\t", ",").replace("[Header]\n", "[Header],,,\n").replace("[Data]\n", "[Data],,,\n")
    comma = _load_text(tmp_path, comma_text, "comma.csv")
    assert sniff(tmp_path / "comma.csv").format == "illumina_final_report"
    assert comma.calls.collect().drop("facet_id").equals(tab.calls.collect().drop("facet_id"))
    dropped = tab.facets.collect()["dropped_fields"][0].to_list()
    assert dropped == ["GenCall Score", "LogR", "BFreq", "X", "Y", "Xraw", "Yraw"]


def test_final_report_with_two_samples_and_unplaced_probes(tmp_path: Path) -> None:
    text = REPORT_HEAD + (
        "SNP Name\tSample ID\tChr\tMapInfo\tAlleles - Plus\tGC Score\n"
        "rs1\tS1\t1\t100\tAG\t0.9\n"
        "rs1\tS2\t1\t100\tGG\t0.9\n"
        "rs2\tS1\t0\t0\tAA\t0.8\n"
        "rs3\tS1\tXY\t300\t-\t0.0\n"
        "rs4\tS1\tMT\t400\tII\t0.7\n"
    )
    with pytest.raises(MultiSampleError):
        _load_text(tmp_path, text, "report.txt")
    loaded = load(tmp_path / "report.txt", sample="S1")
    by_id = _by_id(loaded)
    assert {k: (v["chrom"], v["gt_alleles"]) for k, v in by_id.items()} == {
        "rs1": ("1", ["A", "G"]),
        "rs3": ("XY", [None, None]),
        "rs4": ("MT", ["I", "I"]),
    }
    assert loaded.warnings == ("1 row(s) on chromosome 0 (probes with no mapped position) dropped",)
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["sample"] == "S1" and facet["assembly_label"] is None, "a manifest name is not a build statement"
    assert facet["dropped_fields"] == ["GC Score"]


def test_table_samples_lists_the_sample_column(tmp_path: Path) -> None:
    rows = [("rs1", "S2"), ("rs1", "S1"), ("rs2", "S2"), ("rs3", "S3"), ("rs3", "S1")]
    body = "".join(f"{rsid}\t{sample}\t1\t100\tAG\t0.9\n" for rsid, sample in rows)
    text = REPORT_HEAD + "SNP Name\tSample ID\tChr\tMapInfo\tAlleles - Plus\tGC Score\n" + body
    path = tmp_path / "report.txt"
    path.write_text(text)
    first_seen = list(dict.fromkeys(sample for _, sample in rows))
    assert table_samples(path) == first_seen
    with pytest.raises(MultiSampleError, match="table_samples"):
        load(path)
    for sample in first_seen:
        loaded = load(path, sample=sample)
        assert loaded.calls.select(pl.len()).collect().item() == sum(s == sample for _, s in rows)
    color = tmp_path / "color.csv"
    color.write_text(COLOR_HEADER + "K1,2,100,rs1,AG,GRCh37\n")
    assert table_samples(color) == ["K1"]
    ancestry = tmp_path / "ancestry.txt"
    ancestry.write_text(get_spec("ancestrydna").example)
    assert table_samples(ancestry) == [], "no sample column, nothing to choose"
    vcf = tmp_path / "calls.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")
    with pytest.raises(UnsupportedFormatError):
        table_samples(vcf)


def test_final_report_without_plus_strand_alleles_is_refused(tmp_path: Path) -> None:
    text = REPORT_HEAD + (
        "SNP Name\tSample ID\tAllele1 - Forward\tAllele2 - Forward\tGC Score\n"
        "rs1\tS1\tA\tG\t0.9\n"
    )
    with pytest.raises(TableLayoutError, match="not the genomic plus strand") as refused:
        _load_text(tmp_path, text, "forward.txt")
    message = str(refused.value)
    assert "no column for chromosome, position, genotype_plus" in message
    assert "the positions and plus-strand alleles would come from the chip's manifest, GSA-24v2-0_A1.bpm" in message
    unnamed = text.replace("Content\tGSA-24v2-0_A1.bpm\n", "")
    with pytest.raises(TableLayoutError, match="which the header does not name"):
        _load_text(tmp_path, unnamed, "unnamed.txt")


# --- Real files (opt-in) -----------------------------------------------------------------------------

REAL_LAYOUTS = {
    "ancestrydna": LAYOUTS["ancestrydna"],
    "ftdna": LAYOUTS["ftdna"],
    "decodeme": LAYOUTS["decodeme"],
    "genos": LAYOUTS["genos"],
    "illumina_final_report": LAYOUTS["illumina_final_report"],
}


@pytest.fixture(scope="session")
def upload_formats(pgp_participants: Path) -> dict[str, Sniffed]:
    """The first participant upload of each table format, by content (sniff), never by name."""
    found: dict[str, Sniffed] = {}
    for path in sorted(pgp_participants.glob("*/upload-*/*")):
        if len(found) == len(REAL_LAYOUTS):
            break
        if not path.is_file() or path.stat().st_size > 60_000_000:
            continue
        result = sniff(path)
        if result.format in REAL_LAYOUTS and result.format not in found:
            with open_payload(result) as handle:
                if result.format == "illumina_final_report" and b"Alleles - Plus" not in handle.read(1 << 16):
                    continue
            found[result.format] = result
    return found


@pytest.mark.parametrize("format_id", sorted(REAL_LAYOUTS))
def test_real_upload_loads_every_data_row(upload_formats: dict[str, Sniffed], format_id: str) -> None:
    sniffed = upload_formats.get(format_id)
    if sniffed is None:
        pytest.skip(f"no {format_id} upload found")
    with open_payload(sniffed) as handle:
        text = handle.read().decode("utf-8", errors="replace")
    expected = _expected_rows(text, REAL_LAYOUTS[format_id])
    loaded = load(sniffed)
    calls = loaded.calls.collect()
    assert calls.height == len(expected), "every placed data row becomes one calls row"
    assert sorted(calls["source_gt"].to_list()) == sorted(r["source_gt"] for r in expected)
    assert calls["ref"].null_count() == calls.height
    if format_id == "decodeme":
        facet = loaded.facets.collect().row(0, named=True)
        assert facet["assembly_label"] == "GRCh36", "deCODEme states its build in the info member of its zip"
