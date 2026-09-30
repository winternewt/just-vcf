import gzip
import zipfile
from pathlib import Path

import polars as pl
import pytest

from just_vcf import MultiSampleError, load, vcf_samples
from just_vcf.schema import BLOCKS, CALLS

FORMAT_LINES = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    "PL": '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "GP": '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posteriors">',
    "PP": '##FORMAT=<ID=PP,Number=G,Type=Integer,Description="Phred posteriors">',
    "PS": '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "FT": '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">',
    "LEN": '##FORMAT=<ID=LEN,Number=1,Type=Integer,Description="Block length">',
    "MIN_DP": '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Block minimum depth">',
}
INFO_LINES = {
    "END": '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    "SVLEN": '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">',
    "DP": '##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth">',
}


def vcf_text(records: list[str], samples: tuple[str, ...] = ("S1",), fmt: tuple[str, ...] = ("GT",),
             info: tuple[str, ...] = (), version: str = "VCFv4.2", extra_meta: tuple[str, ...] = ()) -> str:
    meta = [f"##fileformat={version}", "##contig=<ID=1,length=249250621>", "##contig=<ID=2,length=243199373>"]
    meta += [INFO_LINES[k] for k in info] + [FORMAT_LINES[k] for k in fmt] + list(extra_meta)
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] + (["FORMAT", *samples] if samples else [])
    return "\n".join([*meta, "\t".join(columns), *records]) + "\n"


def write(tmp_path: Path, text: str, name: str = "calls.vcf") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def frames(path: Path, **kwargs) -> tuple[pl.DataFrame, pl.DataFrame, tuple[str, ...]]:
    loaded = load(path, **kwargs)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    return calls, blocks, loaded.warnings


def assert_gt_consistent(calls: pl.DataFrame) -> None:
    """`gt_alleles` is `[source_ref] + alts` indexed by `gt`, row by row."""
    for row in calls.iter_rows(named=True):
        if row["gt"] is None:
            continue
        alleles = [row["source_ref"], *row["alts"]]
        assert row["gt_alleles"] == [None if i is None else alleles[i] for i in row["gt"]], row


def test_single_sample_fields(tmp_path: Path) -> None:
    records = [
        "1\t100\trs1;rs2\tA\tG,T\t50\tq10;s50\tDP=99\tGT:GQ:DP:AD:PL:FT\t1/2:30:7:0,3,4:90,40,0,50,20,80:LowGQ",
        "1\t200\t.\tACG\tA\t.\tPASS\t.\tGT:GQ:DP:AD:PL:FT\t0/1:20:9:5,4:30,0,60:PASS",
        "2\t300\t.\tT\t.\t.\t.\t.\tGT:GQ:DP:AD:PL:FT\t0/0:.:.:.:.:.",
    ]
    path = write(tmp_path, vcf_text(records, fmt=("GT", "GQ", "DP", "AD", "PL", "FT"), info=("DP",)))
    calls, blocks, warnings = frames(path)
    assert calls.schema == CALLS and blocks.schema == BLOCKS and blocks.height == 0
    assert calls.height == len(records) and warnings == ()
    first = calls.row(0, named=True)
    assert first["ids"] == ["rs1", "rs2"] and first["filters"] == ["q10", "s50"]
    assert first["alts"] == ["G", "T"] and first["gt"] == [1, 2] and first["gt_alleles"] == ["G", "T"]
    assert first["dp"] == 7, "FORMAT/DP, not INFO/DP"
    assert first["ad"] == [0, 3, 4] and first["sample_filters"] == ["LowGQ"]
    assert first["gl"] == pytest.approx([-9.0, -4.0, 0.0, -5.0, -2.0, -8.0])
    assert first["source_alleles"] == ["A", "G", "T"] and first["ref"] is None and first["ref_check"] == "not_run"
    second = calls.row(1, named=True)
    assert second["end_pos"] == 202 and second["ids"] is None and second["filters"] == ["PASS"]
    third = calls.row(2, named=True)
    assert third["alts"] == [] and third["filters"] is None and third["source_alleles"] is None
    assert_gt_consistent(calls)


def test_multi_sample_needs_a_choice(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:30\t1/1:40"]
    path = write(tmp_path, vcf_text(records, samples=("mother", "child"), fmt=("GT", "GQ")))
    assert vcf_samples(path) == ["mother", "child"]
    with pytest.raises(MultiSampleError):
        load(path)
    with pytest.raises(MultiSampleError):
        load(path, sample="father")
    child, _, _ = frames(path, sample="child")
    mother, _, _ = frames(path, sample="mother")
    assert (child["gt"][0].to_list(), child["gq"][0]) == ([1, 1], 40)
    assert (mother["gt"][0].to_list(), mother["gq"][0]) == ([0, 1], 30)
    assert load(path, sample="child").facets.collect()["facet_id"][0] != load(path, sample="mother").facets.collect()["facet_id"][0]


def test_blocks_are_all_star_records_asserting_reference(tmp_path: Path) -> None:
    records = [
        "1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ:MIN_DP:PL\t0/0:60:23:0,60,900",
        "1\t114\t.\tT\tTC,<NON_REF>\t213\t.\t.\tGT:GQ:MIN_DP:PL\t0/1:99:.:51,0,36,93,92,86",
        "1\t116\t.\tC\t<*>\t.\t.\tEND=120\tGT:GQ:MIN_DP:PL\t./.:.:.:.",
        "1\t121\t.\tA\t<*>\t.\t.\tEND=125\tGT:GQ:MIN_DP:PL\t0/1:5:3:40,0,10",
    ]
    path = write(tmp_path, vcf_text(records, fmt=("GT", "GQ", "MIN_DP", "PL"), info=("END",)))
    calls, blocks, _ = frames(path)
    assert blocks["pos"].to_list() == [100, 116]
    assert blocks["end_pos"].to_list() == [113, 120]
    assert blocks["kind"].to_list() == ["ref", "no_call"]
    assert blocks["dp_min"].to_list() == [23, None]
    assert blocks["gl"][0].to_list() == pytest.approx([0.0, -6.0, -90.0])
    assert calls["pos"].to_list() == [114, 121], "a variant with <*> and a GT calling <*> stay calls"
    assert calls["alts"].to_list() == [["TC", "<*>"], ["<*>"]]
    assert load(path).facets.collect()["absence_means"][0] == "blocks"


def test_end_pos_follows_the_vcf_45_rule(tmp_path: Path) -> None:
    records = [
        "1\t200\t.\tT\t<DEL>\t.\t.\tSVLEN=30\tGT:LEN\t0/1:.",
        "1\t300\t.\tT\t<DUP:TANDEM>\t.\t.\tSVLEN=-12\tGT:LEN\t0/1:.",
        "1\t400\t.\tG\t<INS>\t.\t.\tSVLEN=500\tGT:LEN\t0/1:.",
        "1\t500\t.\tA\tA[2:321682[\t.\t.\t.\tGT:LEN\t0/1:.",
        "1\t600\t.\tC\t<*>\t.\t.\t.\tGT:LEN\t0/0:14",
        "1\t700\t.\tCTTT\tC\t.\t.\t.\tGT:LEN\t0/1:.",
    ]
    path = write(tmp_path, vcf_text(records, version="VCFv4.5", fmt=("GT", "LEN"), info=("SVLEN", "END")))
    calls, blocks, warnings = frames(path)
    ends = dict(zip(calls["pos"].to_list(), calls["end_pos"].to_list(), strict=True))
    assert ends == {200: 230, 300: 312, 400: 400, 500: 500, 700: 703}
    assert blocks["end_pos"].to_list() == [613]
    assert any("P4" in w for w in warnings)
    assert calls["source_gt"].null_count() == calls.height, "4.5: polars-bio rewrites GT, so it is withheld"
    assert {row[0] for row in calls["phased"].to_list()} == {None}


def test_phasing_and_phase_sets(tmp_path: Path) -> None:
    records = [
        "1\t100\t.\tA\tG\t.\t.\t.\tGT:PS\t0|1:100",
        "1\t200\t.\tC\tT\t.\t.\t.\tGT:PS\t1|0:.",
        "1\t300\t.\tG\tA\t.\t.\t.\tGT:PS\t0/1:.",
        "2\t100\t.\tA\tG\t.\t.\t.\tGT:PS\t1|0:100",
    ]
    calls, _, _ = frames(write(tmp_path, vcf_text(records, fmt=("GT", "PS"))))
    rows = {(r["chrom"], r["pos"]): r for r in calls.iter_rows(named=True)}
    assert rows[("1", 100)]["phased"] == [True, True] and rows[("1", 100)]["phase_sets"] == ["100", "100"]
    assert rows[("1", 200)]["phase_sets"] == [None, None] and rows[("1", 200)]["phase_set_kinds"] == ["ps", "ps"]
    assert rows[("1", 300)]["phased"] == [False, False] and rows[("1", 300)]["phase_sets"] is None
    assert rows[("2", 100)]["phase_sets"] == ["100", "100"], "same PS on another chromosome; keyed with chrom downstream"
    assert [r["source_gt"] for r in rows.values()] == ["0|1", "1|0", "0/1", "1|0"]


def test_header_repairs(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"]
    sorted_meta = vcf_text(records).split("\n")
    sorted_meta = sorted_meta[1:3] + [sorted_meta[0]] + sorted_meta[3:]
    calls, _, warnings = frames(write(tmp_path, "\n".join(sorted_meta), "late.vcf"))
    assert calls["gt"].to_list() == [[0, 1]] and any("reordered" in w for w in warnings)

    undeclared = vcf_text(records, fmt=())
    calls, _, warnings = frames(write(tmp_path, undeclared, "nogt.vcf"))
    assert calls["gt"].to_list() == [[0, 1]] and any("GT" in w and "not declared" in w for w in warnings)


def test_out_of_range_values_are_counted(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:300", "1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ\t0/1:40"]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"))))
    assert calls["gq"].to_list() == [None, 40]
    assert warnings == ("1 GQ value(s) outside 0..255 set to null",)


@pytest.mark.parametrize(
    ("version", "values", "expected", "warned"),
    [
        ("VCFv4.2", "0,20,30", [1.0, 0.01, 0.001], False),
        ("VCFv4.3", "0.9,0.09,0.01", [0.9, 0.09, 0.01], False),
        ("VCFv4.3", "0,20,30", None, True),
        ("VCFv4.2", "0.9,0.09,0.01", None, True),
    ],
)
def test_gp_scale_follows_the_header(tmp_path: Path, version: str, values: str, expected, warned: bool) -> None:
    records = [f"1\t100\t.\tA\tG\t.\t.\t.\tGT:GP\t0/0:{values}"]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records, version=version, fmt=("GT", "GP"))))
    got = calls["gp"][0]
    assert (got is None) if expected is None else got.to_list() == pytest.approx(expected, rel=1e-5)
    assert any("GP scale" in w for w in warnings) is warned


def test_pp_converts_to_gp(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:PP\t0/0:0,10,20"]
    calls, _, _ = frames(write(tmp_path, vcf_text(records, fmt=("GT", "PP"))))
    assert calls["gp"][0].to_list() == pytest.approx([1.0, 0.1, 0.01], rel=1e-5)


def test_frames_can_be_collected_together_and_again(tmp_path: Path) -> None:
    records = [f"{c}\t{p}\t.\tA\tG\t.\t.\t.\tGT\t0/1" for c, p in [("2", 5), ("1", 9), ("1", 3)]]
    records.insert(0, "1\t1\t.\tA\t<*>\t.\t.\tEND=2\tGT\t0/0")
    path = write(tmp_path, vcf_text(records, info=("END",)))
    loaded = load(path)
    first = pl.collect_all([loaded.calls, loaded.blocks])
    second = pl.collect_all([loaded.calls, loaded.blocks])
    assert all(a.equals(b) for a, b in zip(first, second, strict=True))
    assert list(zip(first[0]["chrom"], first[0]["pos"])) == [("1", 3), ("1", 9), ("2", 5)]
    again = load(path)
    assert again.calls.collect().equals(first[0]) and again.facets.collect().equals(loaded.facets.collect())


def test_archived_vcf_is_rewritten_and_loaded(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1", "1\t200\t.\tC\tT\t.\t.\t.\tGT\t1/1"]
    path = tmp_path / "upload.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("calls.vcf.gz", gzip.compress(vcf_text(records).encode()))
    loaded = load(path)
    assert loaded.calls.collect()["gt"].to_list() == [[0, 1], [1, 1]]
    assert loaded.facets.collect()["member"][0] == "calls.vcf.gz"


def test_real_vcf_upload(pgp_participants: Path) -> None:
    from just_vcf.sniff import sniff

    candidates = sorted(pgp_participants.glob("*/upload-*/*.vcf*"), key=lambda p: p.stat().st_size)
    found = next((p for p in candidates[:200] if p.stat().st_size > 10_000 and sniff(p).format == "vcf"), None)
    if found is None:
        pytest.skip("no plain VCF upload found")
    opener = gzip.open if found.read_bytes()[:2] == b"\x1f\x8b" else open
    with opener(found, "rt", errors="replace") as handle:
        data_lines = sum(1 for line in handle if line.strip() and not line.startswith("#"))
    loaded = load(found, sample=vcf_samples(found)[0] if vcf_samples(found) else None)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    assert calls.height + blocks.height == data_lines
    assert_gt_consistent(calls)
    assert calls.filter(pl.col("end_pos") < pl.col("pos")).height == 0


def test_malformed_field_the_frames_do_not_carry_does_not_fail_the_load(tmp_path: Path) -> None:
    extra = ('##FORMAT=<ID=XAF,Number=1,Type=Float,Description="A vendor annotation">',)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:XAF\t0/1:not_a_float"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["gt"].to_list() == [[0, 1]]
    assert loaded.facets.collect()["dropped_fields"][0].to_list() == ["FORMAT/XAF"]


def test_header_polars_bio_refuses_is_retried_with_structural_lines(tmp_path: Path) -> None:
    extra = ('##GATKCommandLine.X=<ID=X,Version=1,CommandLineOptions="reference=C:\\Users\\x\\ref.fa">',)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["gt"].to_list() == [[0, 1]]
    assert any("refused the file" in w for w in loaded.warnings)


def test_a_scalar_key_declared_as_a_list_keeps_its_value(tmp_path: Path) -> None:
    extra = ('##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Platypus-style GQ">',)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:45", "1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ\t0/1:12,30"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["gq"].to_list() == [45, 12]
    assert loaded.warnings == ("1 GQ value(s) held more than one number; the first was kept",)


def test_duplicate_declarations_are_dropped_in_the_retry_copy(tmp_path: Path) -> None:
    extra = ('##FILTER=<ID=LowQual,Description="one">', '##FILTER=<ID=LowQual,Description="two">')
    records = ["1\t100\t.\tA\tG\t.\tLowQual\t.\tGT\t0/1"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["filters"].to_list() == [["LowQual"]]
    assert any("refused the file" in w for w in loaded.warnings)
