"""BCF loads into the same frames as the text VCF it was written from (pysam writes the fixtures)."""

import gzip
import struct
import zipfile
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff, vcf_samples
from just_vcf.bgzf import BgzfWriter
from just_vcf.loaders.vcf import HeaderRefused
from test_load_vcf import FORMAT_LINES, INFO_LINES, assert_gt_consistent, vcf_text, write

pysam = pytest.importorskip("pysam", reason="pysam (the [reference] extra) writes the BCF fixtures")

RECORDS = [
    "1\t100\trs1;rs2\tA\tG,T\t50\tq10\t.\tGT:GQ:DP:AD:PL\t1/2:30:7:0,3,4:90,40,0,50,20,80",
    "1\t200\t.\tACG\tA\t.\tPASS\t.\tGT:GQ:DP:AD:PL\t0|1:20:9:5,4:30,0,60",
    "1\t300\t.\tG\t<*>\t.\t.\tEND=313\tGT:GQ:DP:AD:PL\t0/0:60:.:.:0,60,900",
    "2\t5\t.\tT\t.\t.\t.\t.\tGT:GQ:DP:AD:PL\t./.:.:.:.:.",
    "2\t10\t.\tC\tT\t.\t.\t.\tGT:GQ:DP:AD:PL\t1:.:.:.:.",
]
FMT = ("GT", "GQ", "DP", "AD", "PL")


def to_bcf(text_vcf: Path, target: Path, mode: str = "wb") -> Path:
    with pysam.VariantFile(str(text_vcf)) as source, pysam.VariantFile(str(target), mode, header=source.header) as out:
        for record in source:
            out.write(record)
    return target


def collected(path: Path) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    loaded = load(path)
    facets, calls, blocks = pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])
    return facets, calls, blocks


@pytest.fixture
def text_vcf(tmp_path: Path) -> Path:
    q10 = ('##FILTER=<ID=q10,Description="Quality below 10">',)
    return write(tmp_path, vcf_text(RECORDS, fmt=FMT, info=("END",), extra_meta=q10), "calls.vcf")


def test_bcf_loads_the_same_rows_as_its_text(tmp_path: Path, text_vcf: Path) -> None:
    bcf = to_bcf(text_vcf, tmp_path / "calls.bcf")
    assert sniff(bcf).format == "bcf"
    text_facets, text_calls, text_blocks = collected(text_vcf)
    facets, calls, blocks = collected(bcf)
    assert facets["format_id"].to_list() == ["bcf"] and facets["facet_id"][0] != text_facets["facet_id"][0]
    assert calls.height == len(RECORDS) - 1 and blocks.height == 1
    assert calls.drop("facet_id", "source_gt").equals(text_calls.drop("facet_id", "source_gt"))
    assert blocks.drop("facet_id").equals(text_blocks.drop("facet_id"))
    assert calls["source_gt"].null_count() == calls.height, "a BCF wrote no genotype text"
    assert text_calls["source_gt"].null_count() == 0
    assert_gt_consistent(calls)
    keep = ["header_contigs", "vcf_version", "sample", "dropped_fields", "absence_means", "ploidy_statements"]
    assert facets.select(keep).equals(text_facets.select(keep))
    assert vcf_samples(bcf) == ["S1"]


@pytest.mark.parametrize("how", ["renamed", "uncompressed", "zipped"])
def test_a_bcf_under_another_name_or_layer_is_copied_and_read(tmp_path: Path, text_vcf: Path, how: str) -> None:
    _, reference, _ = collected(to_bcf(text_vcf, tmp_path / "reference.bcf"))
    if how == "renamed":
        path = to_bcf(text_vcf, tmp_path / "upload.dat")
    elif how == "uncompressed":
        path = to_bcf(text_vcf, tmp_path / "raw.bcf", mode="wb0")
    else:
        inner = to_bcf(text_vcf, tmp_path / "inner.bcf")
        path = tmp_path / "upload.zip"
        with zipfile.ZipFile(path, "w") as zf:
            zf.write(inner, "genome.bcf")
    sniffed = sniff(path)
    assert sniffed.format == "bcf", sniffed
    _, calls, _ = collected(path)
    assert calls.drop("facet_id").equals(reference.drop("facet_id"))


def test_multi_sample_bcf_needs_a_choice(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:30\t1/1:40"]
    text = write(tmp_path, vcf_text(records, samples=("mother", "child"), fmt=("GT", "GQ")), "trio.vcf")
    bcf = to_bcf(text, tmp_path / "trio.bcf")
    assert vcf_samples(bcf) == ["mother", "child"]
    child = load(bcf, sample="child").calls.collect()
    assert (child["gt"][0].to_list(), child["gq"][0]) == ([1, 1], 40)


def test_contig_dictionary_follows_idx(tmp_path: Path) -> None:
    """A header whose `##contig` lines carry IDX out of line order: the first record's contig is read
    by its dictionary index, so the coordinate check still compares the right names."""
    meta = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=2,length=243199373,IDX=1>",
        "##contig=<ID=1,length=249250621,IDX=0>",
        INFO_LINES["END"],
        FORMAT_LINES["GT"],
    ]
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    text = write(tmp_path, "\n".join([*meta, columns, "1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1", "2\t7\t.\tC\tT\t.\t.\t.\tGT\t1/1"]) + "\n")
    bcf = to_bcf(text, tmp_path / "idx.bcf")
    calls = load(bcf).calls.collect()
    assert list(zip(calls["chrom"], calls["pos"], strict=True)) == [("1", 100), ("2", 7)]


def test_a_bcf_header_states_its_build_as_a_text_header_does(tmp_path: Path) -> None:
    """`sniff` stops at the BCF magic, so the loader reads the header's build lines itself (found by
    the facet gate, tests/test_facet_gate.py)."""
    text = write(tmp_path, vcf_text(RECORDS[1:2], fmt=FMT, extra_meta=("##reference=file:///ref/hs37d5.fa",)), "built.vcf")
    bcf = to_bcf(text, tmp_path / "built.bcf")
    assert sniff(bcf).build is None, "the binary head carries no build hint"
    facets = [load(path).facets.collect().row(0, named=True) for path in (text, bcf)]
    assert [(f["assembly_label"], f["assembly_label_basis"]) for f in facets] == [("GRCh37", "stated")] * 2


def _retype_bcf_gt(bcf: Path, kind: str) -> None:
    """Rewrite a BCF's header text to declare GT `Type=kind`: htslib will not write such a BCF (it
    cannot parse a GT `0/1` under Integer), but a BCF's records encode GT the same way whatever the line."""
    payload = gzip.decompress(bcf.read_bytes())
    (length,) = struct.unpack("<I", payload[5:9])
    text = payload[9 : 9 + length].replace(b"ID=GT,Number=1,Type=String", f"ID=GT,Number=1,Type={kind}".encode())
    with bcf.open("wb") as raw, BgzfWriter(raw) as out:
        out.write(payload[:5] + struct.pack("<I", len(text)) + text + payload[9 + length :])


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3"])
def test_a_bcf_declaring_gt_another_type_is_refused_naming_it(tmp_path: Path, version: str) -> None:
    """A text VCF whose GT is declared `Integer` is read through a copy declaring it String, and its
    warning says so. A BCF is read as stored, with no copy, and polars-bio refuses its header (under
    4.3 with "invalid record", naming nothing), so the load refuses it naming GT and its Type."""
    record = ["1\t100\t.\tA\tG\t50\t.\t.\tGT\t0/1"]
    bcf = to_bcf(write(tmp_path, vcf_text(record, version=version), "gt.vcf"), tmp_path / "gt.bcf")
    _retype_bcf_gt(bcf, "Integer")
    with pytest.raises(HeaderRefused, match=r"BCF declaring FORMAT/GT Type=Integer"):
        load(bcf)
    meta = ('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="x">',)
    text = load(write(tmp_path, vcf_text(record, fmt=(), extra_meta=meta, version=version), "gt_integer.vcf"))
    assert text.calls.collect()["gt"].to_list() == [[0, 1]]
    assert [w for w in text.warnings if w.startswith("FORMAT/GT declared Type=Integer")], text.warnings


def test_a_bcf_float_integer_key_is_read_as_stored_with_no_retype_warning(tmp_path: Path) -> None:
    """A BCF stores a Float as a 32-bit float and its header cannot be rewritten, so no copy declares
    the key String: the text VCF's retype warning would name a read that never happened."""
    meta = ('##FORMAT=<ID=DP,Number=1,Type=Float,Description="x">',)
    text = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t50\t.\t.\tGT:DP\t0/1:30"], extra_meta=meta), "float.vcf")
    bcf = to_bcf(text, tmp_path / "float.bcf")
    retyped = [[w for w in loaded.warnings if w.startswith("integer key(s) declared Type=Float")] for loaded in (load(text), load(bcf))]
    assert [len(found) for found in retyped] == [1, 0], retyped
    assert load(bcf).calls.collect()["dp"].to_list() == [30]
