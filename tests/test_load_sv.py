"""Structural-variant and copy-number facets (RM282): the SV columns, the extent of copy-number
alleles, copy-number segments kept out of `blocks`, the VRS reason, and the round trip through
`write_vcf`. Synthetic records in the shapes of Manta, DRAGEN's CNV caller, Canvas, DELLY and the
VCF 4.5 examples; the PGP-HMS files these shapes come from are counted in docs/SCHEMA.md."""

import gzip
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, normalize, resolve_ref, vrs_withheld, write_vcf
from just_vcf.schema import CALLS
from just_vcf.vrs import WITHHELD_REASONS
from test_load_vcf import frames, vcf_text, write
from test_normalize import HEADER, HOMOPOLYMER, at, store  # noqa: F401  (store is a fixture)

SV_META = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">',
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype">',
)
"""A VCF 4.2 SV caller's declarations, as Manta, DELLY and Canvas write them (CIPOS Number=2, CN Integer)."""

SV_COLUMNS = ["pos", "end_pos", "alts", "sv_type", "cipos", "ciend", "copy_number", "alt_copy_numbers", "gt"]

RECORDS = [
    # Manta: a deletion with both END and SVLEN, an imprecise interval on each side.
    "1\t1000\tdel1\tT\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1300;SVLEN=-300;CIPOS=-5,5;CIEND=-10,12;IMPRECISE\tGT\t0/1",
    "1\t2000\tdup1\tA\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP;END=2100;SVLEN=100\tGT\t1/1",
    "1\t3000\tinv1\tC\t<INV>\t.\tPASS\tSVTYPE=INV;END=3500\tGT\t0/1",
    "1\t4000\tins1\tG\t<INS>\t.\tPASS\tSVTYPE=INS;END=4000;SVLEN=250\tGT\t0/1",
    # A breakend pair: each end anchors at its own position, whatever its mate says.
    "1\t5000\tbnd1\tA\tA]2:7000]\t.\tPASS\tSVTYPE=BND;MATEID=bnd2;CIPOS=0,6\tGT\t0/1",
    "2\t7000\tbnd2\tC\tC]1:5000]\t.\tPASS\tSVTYPE=BND;MATEID=bnd1;CIPOS=0,6\tGT\t0/1",
    # DRAGEN's CNV caller: a <DEL>/<DUP> with SVTYPE CNV and the sample's copy number in FORMAT.
    "1\t6000\tcnv1\tN\t<DEL>\t.\tPASS\tSVTYPE=CNV;END=6999\tGT:CN\t0/1:1",
    # Canvas: copy-number alleles <CNn> to END, one interval for the record over two ALTs, and a
    # copy-neutral segment written ALT `.` to END, which is a copy-number call, not a reference block.
    "1\t8000\tcanvas1\tN\t<CN0>,<CN2>\t.\tPASS\tSVTYPE=LOH;END=8999;CIPOS=-20,30\tGT:CN\t1/2:2",
    "1\t9000\tcanvas2\tN\t.\t.\tPASS\tEND=9999\tGT:CN\t0/0:2",
    "1\t10000\tcanvas3\tN\t<CN3>\t.\tPASS\tSVTYPE=CNV;END=10999\tGT:CN\t0/1:3",
    # An interval that is not one pair per ALT is withheld and counted.
    "1\t11000\tbad\tT\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=11050;CIPOS=-5,5,7\tGT\t0/1",
    # An Illumina gVCF reference block (ALT `.` with END and no copy number) stays a block.
    "1\t12000\t.\tA\t.\t.\tPASS\tEND=12099\tGT\t0/0",
    # A sequence allele carries none of the SV columns.
    "1\t13000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
]


def sv_vcf(tmp_path: Path, records: list[str] = RECORDS, name: str = "sample.sv.vcf") -> Path:
    return write(tmp_path, vcf_text(records, info=("END", "SVLEN"), extra_meta=SV_META), name)


def by_id(calls: pl.DataFrame) -> dict[str, dict]:
    return {row["ids"][0]: row for row in calls.iter_rows(named=True) if row["ids"]}


def test_sv_columns_are_read_as_written(tmp_path: Path) -> None:
    calls, blocks, warnings = frames(sv_vcf(tmp_path))
    assert calls.schema == CALLS
    rows = by_id(calls)

    # SVTYPE as the caller wrote it: DRAGEN's CNV beside its <DEL>, Canvas's LOH beside <CN0>,<CN2>.
    assert {key: rows[key]["sv_type"] for key in rows} == {
        "del1": "DEL", "dup1": "DUP", "inv1": "INV", "ins1": "INS", "bnd1": "BND", "bnd2": "BND",
        "cnv1": "CNV", "canvas1": "LOH", "canvas2": None, "canvas3": "CNV", "bad": "DEL",
    }  # fmt: skip
    # One pair per ALT; a record's one pair is that of each of its ALTs; CIEND never borrowed from CIPOS.
    assert rows["del1"]["cipos"] == [-5, 5] and rows["del1"]["ciend"] == [-10, 12]
    assert rows["bnd1"]["cipos"] == [0, 6] and rows["bnd1"]["ciend"] is None
    assert rows["canvas1"]["cipos"] == [-20, 30, -20, 30]
    assert rows["dup1"]["cipos"] is None and rows["bad"]["cipos"] is None
    assert any(w.startswith("1 CIPOS value(s) not two per ALT") for w in warnings), warnings
    # The sample's copy number, from FORMAT CN.
    assert {key: rows[key]["copy_number"] for key in ("cnv1", "canvas1", "canvas2", "canvas3", "del1")} == {
        "cnv1": 1.0, "canvas1": 2.0, "canvas2": 2.0, "canvas3": 3.0, "del1": None,
    }  # fmt: skip

    # Copy-number alleles end at END; a copy-neutral segment is a call over its extent.
    ends = {key: rows[key]["end_pos"] for key in rows}
    assert ends == {
        "del1": 1300, "dup1": 2100, "inv1": 3500, "ins1": 4000, "bnd1": 5000, "bnd2": 7000,
        "cnv1": 6999, "canvas1": 8999, "canvas2": 9999, "canvas3": 10999, "bad": 11050,
    }  # fmt: skip
    assert rows["canvas2"]["alts"] == [] and rows["canvas2"]["gt"] == [0, 0]
    # Only the gVCF block is a block: a segment with a copy number says nothing about small variants.
    assert blocks.select("pos", "end_pos", "kind").rows() == [(12000, 12099, "ref")]
    sequence = calls.filter(pl.col("pos") == 13000).row(0, named=True)
    assert [sequence[c] for c in ("sv_type", "cipos", "ciend", "copy_number", "alt_copy_numbers")] == [None] * 5


def test_sv_facet_counts_its_classes_and_its_absence_is_unknown(tmp_path: Path) -> None:
    loaded = load(sv_vcf(tmp_path))
    facet = loaded.facets.collect().row(0, named=True)
    calls = loaded.calls.collect()
    sv_alleles = int(calls.filter(pl.col("pos") != 13000)["alts"].list.len().sum())
    assert facet["variant_classes"] == [{"variant_class": "snv", "alleles": 1}, {"variant_class": "sv", "alleles": sv_alleles}]
    # The gVCF block in this synthetic file still makes absence a matter of blocks; a pure SV file
    # with copy-neutral segments has none (below).
    no_block = [r for r in RECORDS if "12000" not in r]
    pure = load(sv_vcf(tmp_path, no_block, "pure.sv.vcf"))
    assert pure.blocks.collect().height == 0
    assert pure.facets.collect()["absence_means"].to_list() == ["unknown"]


def test_vcf_45_per_allele_copy_numbers(tmp_path: Path) -> None:
    """VCF 4.5's own example: INFO CN per ALT allele beside FORMAT CN, the sample's total."""
    meta = (
        '##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number of allele">',
        '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS">',
        '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">',
    )
    records = [
        "1\t100\ta\tT\t<CNV>,<CNV>\t.\t.\tSVLEN=30,30;CN=1,2;CIPOS=-5,5,0,0\tGT:CN\t1/2:3",
        "1\t200\tb\tT\t<CNV>\t.\t.\tSVLEN=30;CN=1.5\tGT:CN\t0/1:2.5",
        "1\t300\tc\tT\t<CNV>,<DEL>\t.\t.\tSVLEN=30,30;CN=1\tGT\t1/2",  # not one per ALT: withheld
    ]
    path = write(tmp_path, vcf_text(records, version="VCFv4.5", info=("SVLEN",), extra_meta=meta))
    calls, _, warnings = frames(path)
    rows = by_id(calls)
    assert rows["a"]["alt_copy_numbers"] == [1.0, 2.0] and rows["a"]["copy_number"] == 3.0
    assert rows["a"]["cipos"] == [-5, 5, 0, 0]
    assert rows["b"]["alt_copy_numbers"] == [1.5] and rows["b"]["copy_number"] == 2.5
    assert rows["c"]["alt_copy_numbers"] is None
    assert rows["a"]["end_pos"] == 130
    assert any(w.startswith("1 INFO CN value(s) not one per ALT") for w in warnings), warnings


def test_sv_alleles_are_withheld_from_vrs_with_their_reason(tmp_path: Path) -> None:
    """An SV allele or breakend has no VRS Allele id: `sv_no_vrs_model`, never approximated (charter 4)."""
    calls = load(sv_vcf(tmp_path)).calls
    reasons = vrs_withheld(calls).collect()
    alts = {row["ids"][0]: row["vrs_withheld"][1:] for row in reasons.iter_rows(named=True) if row["ids"]}
    assert alts["del1"] == alts["bnd1"] == alts["cnv1"] == ["sv_no_vrs_model"]
    assert alts["canvas1"] == ["sv_no_vrs_model", "sv_no_vrs_model"]
    assert alts["canvas2"] == []  # a copy-neutral segment has no ALT to mint
    assert "sv_no_vrs_model" in WITHHELD_REASONS
    assert reasons["vrs_allele_ids"].null_count() == reasons.height


def test_a_complete_genomics_no_call_is_not_an_sv(tmp_path: Path) -> None:
    """`<CGA_NOCALL>` without END stays in `calls` and names no variant: `symbolic`, not an SV reason."""
    meta = ('##ALT=<ID=CGA_NOCALL,Description="No-called record">',)
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\t<CGA_NOCALL>\t.\t.\t.\tGT\t./."], extra_meta=meta))
    reasons = vrs_withheld(load(path).calls).collect()
    assert reasons["vrs_withheld"].to_list() == [["ref_unchecked", "symbolic"]]


@pytest.mark.parametrize("compress", [True, False])
def test_sv_rows_round_trip_through_write_vcf(tmp_path: Path, compress: bool) -> None:
    no_block = [r for r in RECORDS if "12000" not in r]
    first = load(sv_vcf(tmp_path, no_block))
    written = write_vcf(first, tmp_path / ("out.vcf.gz" if compress else "out.vcf"), compress=compress)
    back = load(written.path)
    before = first.calls.collect().select(SV_COLUMNS)
    after = back.calls.collect().select(SV_COLUMNS)
    assert after.equals(before), (before, after)
    assert back.blocks.collect().height == 0
    if not compress:
        text = written.path.read_text()
        assert "##INFO=<ID=SVTYPE,Number=1,Type=String" in text
        assert "##FORMAT=<ID=CN,Number=1,Type=Integer" in text
        assert "SVTYPE=LOH;CIPOS=-20,30,-20,30" in text


def info_of(text: str) -> dict[str, dict[str, str]]:
    """Each written record's INFO as key → value, by its ID."""
    out = {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        info = {} if fields[7] == "." else dict(item.split("=", 1) if "=" in item else (item, "") for item in fields[7].split(";"))
        out[fields[2]] = info
    return out


# Records whose ALTs mix length-type SV alleles, beside RECORDS: one length-type ALT beside an <INS>,
# two whose own lengths the frame cannot tell apart, and two that VCF 4.5 says share one (FORMAT CN).
MIXED = [
    "1\t14000\tdelins\tT\t<DEL>,<INS>\t.\tPASS\tSVTYPE=DEL;END=14100;SVLEN=-100,40\tGT\t1/2",
    "1\t15000\tdeldup\tT\t<DEL>,<DUP>\t.\tPASS\tEND=15200;SVLEN=-150,200\tGT\t1/2",
    "1\t16000\tdeldupcn\tT\t<DEL>,<DUP>\t.\tPASS\tEND=16300\tGT:CN\t1/2:2",
]


def test_svlen_is_written_as_vcf_43_states_it(tmp_path: Path) -> None:
    """SVLEN on each length-type SV allele, `end_pos - pos`, negative for a deletion (VCF 4.3's
    examples), `.` on any other ALT. Where the source wrote SVLEN for a one-SV-allele record, the
    value written is the source's; an <INS>'s inserted length is not in the frame and is not made up."""
    no_block = [r for r in RECORDS if "12000" not in r]
    first = load(sv_vcf(tmp_path, [*no_block, *MIXED]))
    written = write_vcf(first, tmp_path / "out.vcf", compress=False)
    text = written.path.read_text()
    assert '##INFO=<ID=SVLEN,Number=.,Type=Integer,' in text
    svlen = {key: info.get("SVLEN") for key, info in info_of(text).items() if key != "."}
    assert svlen == {
        # As the source wrote it (Manta's del1 and dup1), else END - POS from the source's END.
        "del1": "-300", "dup1": "100", "inv1": "500", "cnv1": "-999", "bad": "-50",
        # Copy-number alleles share the segment's length (VCF 4.5); a copy-neutral segment has no ALT.
        "canvas1": "999,999", "canvas3": "999", "canvas2": None,
        # An <INS> and a breakend have no length in the frame.
        "ins1": None, "bnd1": None, "bnd2": None,
        "delins": "-100,.",
        # <DEL>,<DUP> lengths are not kept apart in the frame: left out, counted.
        "deldup": None,
        # ... unless FORMAT CN says they share one (VCF 4.5).
        "deldupcn": "-300,300",
    }  # fmt: skip
    assert any(w.startswith("1 row(s) with SV alleles whose own lengths the frame does not keep apart") for w in written.warnings), written.warnings

    back = load(written.path)
    assert back.calls.collect().select(SV_COLUMNS).equals(first.calls.collect().select(SV_COLUMNS))
    # END and the SVLEN written beside it agree, so the reload counts no disagreement.
    assert not any("disagrees with POS + |SVLEN|" in w for w in back.warnings), back.warnings


def test_htslib_reads_the_written_svlen(tmp_path: Path) -> None:
    """Opt-in external validation: htslib (pysam, the `[reference]` extra) parses every SVLEN as written."""
    pysam = pytest.importorskip("pysam")
    no_block = [r for r in RECORDS if "12000" not in r]
    written = write_vcf(load(sv_vcf(tmp_path, [*no_block, *MIXED])), tmp_path / "out.vcf.gz")
    text = gzip.decompress(written.path.read_bytes()).decode()
    expected = {
        key: tuple(None if value == "." else int(value) for value in info["SVLEN"].split(",")) if "SVLEN" in info else None
        for key, info in info_of(text).items()
    }
    with pysam.VariantFile(str(written.path)) as reader:
        got = {record.id or ".": record.info.get("SVLEN") for record in reader}
    assert got == expected and any(value is not None for value in got.values())


def test_svlen_is_written_unsigned_per_alt_under_vcf_45(tmp_path: Path) -> None:
    """A facet holding a PSL set is written as VCF 4.5, where SVLEN is `Number=A` and never negative."""
    meta = (
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">',
        '##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number of allele">',
        '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">',
        '##FORMAT=<ID=PSL,Number=P,Type=String,Description="Phase set list">',
    )
    records = [
        "1\t5\tsnv\tT\tG\t.\tPASS\t.\tGT:PSL\t0|1:.,1*5*1",
        "1\t100\tcnv\tT\t<CNV>,<CNV>\t.\t.\tSVLEN=30,30;CN=1,2\tGT:CN\t1/2:3",
        "1\t200\tdel\tT\t<DEL>\t.\t.\tSVLEN=300\tGT\t0/1",
        "1\t600\tins\tT\t<INS>\t.\t.\tSVLEN=40\tGT\t0/1",
    ]
    first = load(write(tmp_path, vcf_text(records, version="VCFv4.5", extra_meta=meta), "v45.vcf"))
    written = write_vcf(first, tmp_path / "out.vcf", compress=False)
    text = written.path.read_text()
    assert "##fileformat=VCFv4.5\n" in text
    assert '##INFO=<ID=SVLEN,Number=A,Type=Integer,' in text
    svlen = {key: info.get("SVLEN") for key, info in info_of(text).items()}
    assert svlen == {"snv": None, "cnv": "30,30", "del": "300", "ins": None}
    back = load(written.path)
    assert back.calls.collect().select(SV_COLUMNS).equals(first.calls.collect().select(SV_COLUMNS))


def test_a_vcf_43_copy_number_declared_float_is_not_carried(tmp_path: Path) -> None:
    """VCF 4.3 reserves FORMAT CN as Integer, and polars-bio refuses a 4.3 header declaring it Float.
    The file still loads, without CN, and says so; under 4.5, where CN is Float, it is carried."""
    meta = ('##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">',)
    records = ["1\t100\ta\tT\t<DEL>\t.\t.\tEND=200\tGT:CN\t0/1:1"]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records, version="VCFv4.3", info=("END",), extra_meta=meta), "v43.vcf"))
    assert calls["copy_number"].to_list() == [None] and calls["end_pos"].to_list() == [200]
    assert any("not carried: FORMAT/CN (Type=Float, 4.3 reserves Integer)" in w for w in warnings), warnings
    assert "FORMAT/CN" in load(tmp_path / "v43.vcf").facets.collect()["dropped_fields"][0]
    calls, _, _ = frames(write(tmp_path, vcf_text(records, version="VCFv4.5", info=("END",), extra_meta=meta), "v45.vcf"))
    assert calls["copy_number"].to_list() == [1.0]


def test_an_undeclared_copy_number_is_declared_as_its_version_reserves(tmp_path: Path) -> None:
    records = ["1\t100\ta\tT\t<DEL>\t.\t.\tEND=200;CN=1\tGT:CN\t0/1:1"]
    for version in ("VCFv4.3", "VCFv4.2"):
        calls, _, warnings = frames(write(tmp_path, vcf_text(records, version=version, info=("END",)), f"{version}.vcf"))
        assert calls.select("copy_number", "alt_copy_numbers").rows() == [(1.0, [1.0])], version
        assert any("declared by the loader: INFO/CN, FORMAT/CN" in w for w in warnings), warnings


def test_a_fractional_copy_number_is_left_out_of_the_export_and_counted(tmp_path: Path) -> None:
    meta = ('##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">',)
    records = ["1\t100\ta\tT\t<CNV>\t.\t.\tSVLEN=100\tGT:CN\t0/1:2.5", "1\t300\tb\tT\t<CNV>\t.\t.\tSVLEN=100\tGT:CN\t0/1:3"]
    first = load(write(tmp_path, vcf_text(records, version="VCFv4.5", info=("SVLEN",), extra_meta=meta)))
    written = write_vcf(first, tmp_path / "out.vcf", compress=False)
    assert any(w.startswith("1 row(s) had a copy number that is not a whole number") for w in written.warnings), written.warnings
    back = load(written.path).calls.collect()
    assert back.select("pos", "copy_number").rows() == [(100, None), (300, 3.0)]


def test_normalize_slices_the_per_alt_columns_and_clears_moved_intervals(tmp_path: Path, store: ReferenceStore) -> None:
    """A split row keeps its own ALT's interval and copy number; a row left-alignment moved loses its
    interval, an offset from the POS the source wrote."""
    sv_info = [
        '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS">',
        '##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number of allele">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    ]
    p = 50
    ref = at("1", p)
    g, t = [b for b in "ACGT" if b != ref][:2]
    end_t = HOMOPOLYMER + 7
    records = [
        f"1\t{p}\t.\t{ref}\t{g},{t}\t50\tPASS\tCIPOS=-1,1,-2,2;CN=1,2;SVTYPE=SNV\tGT\t1/2",
        f"1\t{end_t - 1}\t.\t{at('1', end_t - 1, 2)}\t{at('1', end_t - 1)}\t50\tPASS\tCIPOS=0,3;SVTYPE=DEL\tGT\t0/1",
    ]
    path = tmp_path / "in.vcf"
    path.write_text("\n".join([*HEADER[:-1], *sv_info, HEADER[-1], *records]) + "\n")
    out = normalize(resolve_ref(load(path), store), store).calls.collect()
    rows = {(row["pos"], row["alts"][0]): row for row in out.iter_rows(named=True)}
    assert (rows[(p, g)]["cipos"], rows[(p, g)]["alt_copy_numbers"]) == ([-1, 1], [1.0])
    assert (rows[(p, t)]["cipos"], rows[(p, t)]["alt_copy_numbers"]) == ([-2, 2], [2.0])
    assert rows[(p, g)]["sv_type"] == rows[(p, t)]["sv_type"] == "SNV"
    moved = rows[(HOMOPOLYMER - 1, "C")]
    assert moved["cipos"] is None and moved["sv_type"] == "DEL"
