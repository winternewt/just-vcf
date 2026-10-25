import csv
import gzip
import io
import re
import struct
import zipfile
from pathlib import Path

import polars as pl
import pytest

from just_vcf import MultiSampleError, NoSampleError, load, rejoin_split, vcf_samples, write_vcf
from just_vcf.loaders.vcf import FORMAT_KEYS, HeaderRefused, undo_resave
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


def test_sites_only_vcf_is_refused(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t50\tPASS\t.", "2\t200\trs1\tC\tT,G\t.\t.\t."]
    path = write(tmp_path, vcf_text(records, samples=(), fmt=()))
    assert vcf_samples(path) == []
    with pytest.raises(NoSampleError):
        load(path)
    with pytest.raises(MultiSampleError):
        load(path, sample="S1")


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


PSL_LINES = (
    '##FORMAT=<ID=PSL,Number=P,Type=String,Description="Phase set list">',
    '##FORMAT=<ID=PSO,Number=P,Type=Integer,Description="Phase set list ordinal">',
    '##FORMAT=<ID=PSQ,Number=P,Type=Integer,Description="Phase set list quality">',
)


def test_psl_phase_sets_are_read_per_allele(tmp_path: Path) -> None:
    """VCF 4.5's PSL, PSO and PSQ (`Number=P`, which polars-bio refuses, P10): one entry per GT allele."""
    records = [
        "1\t5\t.\tT\tG\t.\tPASS\t.\tGT:PSL\t|0/1:1*5*1,.",  # slot 0 named: phased though P4 drops the `|`
        "1\t10\t.\tA\tT,G\t.\tPASS\t.\tGT:PSL:PSO:PSQ\t|1/2|1:1*10*1,.,1*5*1:1,.,2:30,.,40",
        "2\t15\t.\tG\tC\t.\tPASS\t.\tGT:PSL\t1|0:.,1*10*1",  # slot 0's `.` says nothing: unknown, in no set
        "2\t20\t.\tG\tC\t.\tPASS\t.\tGT:PSL\t1|1:1*5*1,.",  # slot 1 phased with `.`: the implicit set
        "2\t25\t.\tG\tC\t.\tPASS\t.\tGT:PSL\t1/1:.,x",  # a name on an unphased allele is not read
        "2\t30\t.\tG\tC\t.\tPASS\t.\tGT:PSL\t1|1:x",  # one entry for two alleles: phase withheld
        "2\t35\t.\tG\tC\t.\tPASS\t.\tGT:PS:PSL\t1|1:35:y,y",  # PS and PSL both: PSL read
    ]
    path = write(tmp_path, vcf_text(records, fmt=("GT", "PS"), version="VCFv4.5", extra_meta=PSL_LINES))
    calls, _, warnings = frames(path)
    rows = {(r["chrom"], r["pos"]): r for r in calls.iter_rows(named=True)}
    columns = ("phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities")

    def phase(key: tuple) -> tuple:
        return tuple(rows[key][name] for name in columns)

    assert phase(("1", 5)) == ([True, False], ["1*5*1", None], ["psl", None], None, None)
    assert phase(("1", 10)) == ([True, False, True], ["1*10*1", None, "1*5*1"], ["psl", None, "psl"], [1, None, 2], [30, None, 40])
    assert phase(("2", 15)) == ([None, True], [None, "1*10*1"], [None, "psl"], None, None)
    assert phase(("2", 20)) == ([True, True], ["1*5*1", None], ["psl", "ps"], None, None)
    assert phase(("2", 25)) == ([None, False], None, None, None, None)
    assert phase(("2", 30)) == ([None, None], None, None, None, None)
    assert phase(("2", 35)) == ([True, True], ["y", "y"], ["psl", "psl"], None, None)
    assert calls["gt"].to_list() == [[0, 1], [1, 2, 1], [1, 0], [1, 1], [1, 1], [1, 1], [1, 1]]
    assert any("Number=P" in w and "PSL, PSO, PSQ" in w for w in warnings)
    assert any(w.startswith("1 record(s) name a PSL phase set on an unphased allele") for w in warnings)
    assert any(w.startswith("1 record(s) whose PSL is not one value per allele") for w in warnings)
    assert any(w.startswith("1 record(s) carry both PS and PSL") for w in warnings)


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.4", "VCFv4.5"])
@pytest.mark.parametrize(("psl", "pso", "psq"), [("P", "1", "1"), ("1", "P", "P"), ("P", "2", "2"), ("1", "1", "1")])
def test_phase_set_lists_declared_with_a_fixed_count_are_read_per_allele(tmp_path: Path, version: str, psl: str, pso: str, psq: str) -> None:
    """VCF defines PSL, PSO and PSQ as one value per GT allele (`Number=P`). A file declaring one of them
    `Number=1` writes two values on a diploid row, which polars-bio cannot parse under that declaration;
    the loader reads them as it reads `P`, each row's count checked against GT and PSL, so a haploid row
    and a diploid one both load. A `Number=2` declaration, which loads as it is, reads the same."""
    meta = tuple(
        f'##FORMAT=<ID={key},Number={number},Type={kind},Description="x">'
        for key, number, kind in (("PSL", psl, "String"), ("PSO", pso, "Integer"), ("PSQ", psq, "Integer"))
    )
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tGT:PSL:PSO:PSQ\t1|1:a,a:1,2:30,40", "1\t200\t.\tA\tG\t50\t.\t.\tGT:PSL:PSO:PSQ\t1:a:3:50"]
    loaded = load(write(tmp_path, vcf_text(records, version=version, extra_meta=meta)))
    calls, facet = loaded.calls.collect().sort("pos"), loaded.facets.collect().row(0, named=True)
    assert calls["phase_sets"].to_list() == [["a", "a"], ["a"]]
    assert calls["phase_set_ordinals"].to_list() == [[1, 2], [3]]
    assert calls["phase_set_qualities"].to_list() == [[30, 40], [50]]
    assert not {"FORMAT/PSL", "FORMAT/PSO", "FORMAT/PSQ"} & set(facet["dropped_fields"] or [])
    one = [key for key, number in (("PSL", psl), ("PSO", pso), ("PSQ", psq)) if number == "1"]
    stated = [w for w in loaded.warnings if w.startswith("phase-set-list key(s) declared Number=1")]
    assert stated == ([] if not one else [stated[0]]) and (not one or stated[0].endswith(", ".join(one))), loaded.warnings


def test_a_file_without_psl_keeps_its_phase_set_lists_null(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:PS\t0|1:100"]
    calls, _, _ = frames(write(tmp_path, vcf_text(records, fmt=("GT", "PS"))))
    assert calls["phase_set_ordinals"].to_list() == calls["phase_set_qualities"].to_list() == [None]


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
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:70000", "1\t200\t.\tC\tT\t.\t.\t.\tGT:GQ\t0/1:300"]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"))))
    assert calls["gq"].to_list() == [None, 300]
    assert warnings == ("1 GQ value(s) outside 0..65535 set to null",)


@pytest.mark.parametrize(
    ("version", "values", "expected", "warning"),
    [
        ("VCFv4.2", "0,20,30", [1.0, 0.01, 0.001], None),
        ("VCFv4.3", "0.9,0.09,0.01", [0.9, 0.09, 0.01], None),
        ("VCFv4.3", "0,20,30", [1.0, 0.01, 0.001], "GP written as phred under a VCFv4.3 header"),
        ("VCFv4.2", "0.9,0.09,0.01", [0.9, 0.09, 0.01], "GP written as probabilities under a VCFv4.2 header"),
        ("VCFv4.2", "0.5,0.5,0.5", None, "GP is ambiguous"),
        # F10: values the numbers alone cannot settle are what a 4.3+ header defines GP as.
        ("VCFv4.3", "0.5,0.5,0.5", [0.5, 0.5, 0.5], "GP read as probabilities, as the VCFv4.3 header defines it: 1 of 1"),
        ("VCFv4.5", "1,0.1,0.01", [1.0, 0.1, 0.01], "GP read as probabilities, as the VCFv4.5 header defines it"),
        (None, "0.5,0.5,0.5", None, "GP is ambiguous: every value is within [0, 1] but 1 of 1 rows do not sum to 1, and the header states no VCF version"),
        (None, "0,20,30", [1.0, 0.01, 0.001], None),
    ],
)
def test_gp_scale_is_read_from_decisive_values(tmp_path: Path, version: str | None, values: str, expected, warning) -> None:
    records = [f"1\t100\t.\tA\tG\t.\t.\t.\tGT:GP\t0/0:{values}"]
    text = vcf_text(records, version=version or "VCFv4.2", fmt=("GT", "GP"))
    if version is None:
        text = "\n".join(line for line in text.split("\n") if not line.startswith("##"))
    calls, _, warnings = frames(write(tmp_path, text))
    got = calls["gp"][0]
    assert (got is None) if expected is None else got.to_list() == pytest.approx(expected, rel=1e-5)
    assert [w for w in warnings if w.startswith("GP")] == ([] if warning is None else [next(w for w in warnings if w.startswith(warning))])


def test_gp_from_phred_survives_its_own_export(tmp_path: Path) -> None:
    """F10: phred GP converted on load does not sum to 1; `write_vcf` writes it as 4.3 probabilities,
    and the reload keeps it. Tolerance: `rel=1e-6`, since the writer's shortest float text reads back
    as the same float32."""
    records = [
        "1\t100\t.\tA\tG\t.\t.\t.\tGT:GP\t0/0:0,10,20",
        "1\t200\t.\tC\tT\t.\t.\t.\tGT:GP\t0/1:13,0,7",
        "1\t300\t.\tG\tA\t.\t.\t.\tGT:GP\t1/1:30,20,0",
    ]
    first = load(write(tmp_path, vcf_text(records, fmt=("GT", "GP"))))
    original = first.calls.collect()
    assert original["gp"].null_count() == 0
    assert (original["gp"].list.sum() - 1).abs().max() > 0.02, "the case F10 is about: posteriors that do not sum to 1"
    again = load(write_vcf(first, tmp_path / "out.vcf.gz").path)
    back = again.calls.collect()
    assert back["gp"].null_count() == 0 and not any(w.startswith("GP is ambiguous") for w in again.warnings)
    for before, after in zip(original["gp"].to_list(), back["gp"].to_list(), strict=True):
        assert after == pytest.approx(before, rel=1e-6)
    assert back["source_alleles"].to_list() == original["source_alleles"].to_list()


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


# A carried key declared as VCF reserves it, with a value that does not parse under that declaration.
# polars-bio names the key for a list-valued field: (meta line, the bad sample or INFO field, the key).
_UNPARSED = [
    ('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="h">', "GT:HQ\t0/1:{}", "30.5,40", "FORMAT/HQ"),
    ('##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="h">', "GT:HQ\t0/1:{}", "x,40", "FORMAT/HQ"),
    ('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="h">', "GT:AD\t0/1:{}", "x,4", "FORMAT/AD"),
    ('##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="h">', "SVLEN={}\tGT\t0/1", "x", "INFO/SVLEN"),
]


def _unparsed_record(pos: int, field: str, value: str) -> str:
    """A record whose INFO or sample column carries `field` with `value`."""
    if field.startswith("GT"):
        return f"1\t{pos}\t.\tA\tG\t50\t.\t.\t{field.format(value)}"
    return f"1\t{pos}\t.\tA\t<DEL>\t50\t.\t{field.format(value)}"


@pytest.mark.parametrize(("meta", "field", "value", "name"), _UNPARSED)
def test_a_carried_key_polars_bio_names_as_unparsable_is_not_carried(tmp_path: Path, meta: str, field: str, value: str, name: str) -> None:
    """A VCF 4.2 file declaring a carried key as VCF reserves it, with a value that does not parse under
    that declaration, loads without the key where polars-bio's error names it: listed in
    `dropped_fields`, one warning naming it, and every other key still read."""
    records = [_unparsed_record(100, field, value), _unparsed_record(200, field, "5"), "1\t300\t.\tC\tT\t50\t.\t.\tGT:GQ\t1/1:40"]
    loaded = load(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"), extra_meta=(meta,))))
    calls, facet = loaded.calls.collect().sort("pos"), loaded.facets.collect().row(0, named=True)
    assert calls["gt"].to_list() == [[0, 1], [0, 1], [1, 1]] and calls["gq"].to_list() == [None, None, 40]
    assert name in facet["dropped_fields"]
    named = [w for w in loaded.warnings if name in w and "cannot parse" in w]
    assert len(named) == 1 and "not carried" in named[0], loaded.warnings


def test_each_carried_key_polars_bio_cannot_parse_is_left_out_in_turn(tmp_path: Path) -> None:
    """Two carried keys with values that do not parse: polars-bio names one per read, so each is left
    out in its own retry, and the file loads without both."""
    meta = ('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="h">', '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="h">')
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tGT:AD:HQ\t0/1:x,4:30,40", "1\t200\t.\tA\tG\t50\t.\t.\tGT:AD:HQ\t0/1:3,4:30,x"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=meta)))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls["gt"].to_list() == [[0, 1], [0, 1]] and calls["ad"].to_list() == [None, None] and calls["hq"].to_list() == [None, None]
    assert {"FORMAT/AD", "FORMAT/HQ"} <= set(facet["dropped_fields"])
    assert sorted(name for name in ("FORMAT/AD", "FORMAT/HQ") for w in loaded.warnings if name in w and "cannot parse" in w) == ["FORMAT/AD", "FORMAT/HQ"]


def test_an_unparsable_carried_key_polars_bio_does_not_name_still_refuses(tmp_path: Path) -> None:
    """polars-bio names no key when a single-value field fails to parse (GQ `30.5` under Integer), and
    the loader never guesses which key it was: the load is refused, as before."""
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ\t0/1:30.5"]
    with pytest.raises(HeaderRefused, match="Error reading FORMAT field"):
        load(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"))))


def test_header_polars_bio_refuses_is_retried_with_structural_lines(tmp_path: Path) -> None:
    extra = ('##GATKCommandLine.X=<ID=X,Version=1,CommandLineOptions="reference=C:\\Users\\x\\ref.fa">',)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["gt"].to_list() == [[0, 1]]
    assert any("refused the file" in w for w in loaded.warnings)


def test_split_marker_survives_the_structural_retry(tmp_path: Path) -> None:
    extra = (
        '##GATKCommandLine.X=<ID=X,Version=1,CommandLineOptions="reference=C:\\Users\\x\\ref.fa">',
        '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">',
    )
    records = ["1\t100\t.\tA\tG\t.\t.\tOLD_MULTIALLELIC=1:100:A/G/T\tGT\t1/."]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert any("refused the file" in w for w in loaded.warnings)
    assert loaded.calls.collect()["split_signal"].to_list() == ["marker"]


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


def test_records_at_pos_zero_are_left_out_and_counted(tmp_path: Path) -> None:
    records = [
        "NA\t0\trs2828663\tN\tA,C\t.\t.\t.\tGT\t1/1",
        "NA\t0\tkgp9608035\tN\tA,G\t.\t.\t.\tGT\t1/2",
        "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1",
    ]
    loaded = load(write(tmp_path, vcf_text(records)))
    assert loaded.calls.collect()["ids"].to_list() == [["rs1"]]
    assert "2 record(s) at POS 0 left out: polars-bio cannot read them" in loaded.warnings


def test_illumina_gvcf_blocks(tmp_path: Path) -> None:
    records = [
        "1\t100\t.\tA\t.\t.\tPASS\tEND=150\tGT:GQ\t0/0:1200",
        "1\t151\t.\tC\t.\t.\tPASS\t.\tGT:GQ\t0/0:30",
        "1\t152\t.\tG\tT\t50\tPASS\t.\tGT:GQ\t0/1:99",
    ]
    loaded = load(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"), info=("END",))))
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    assert blocks.select("pos", "end_pos", "kind", "gq").rows() == [(100, 150, "ref", 1200)]
    assert calls["pos"].to_list() == [151, 152], "a single monomorphic record without END stays a call"
    assert loaded.facets.collect()["absence_means"][0] == "blocks"


def test_split_rows_are_flagged_and_left_as_written(tmp_path: Path) -> None:
    """RM270: rows that look like one multi-allelic site written one row per ALT get `split_signal`;
    `gt` and `gt_alleles` stay as the file wrote them (option B)."""
    records = [
        "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/0",  # a split 1/2: each 0 is the other ALT
        "1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1",
        "1\t200\t.\tC\tG\t50\tPASS\t.\tGT\t1/1",  # a split 1/1 at C>G,T: the T row is all zero
        "1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/0",
        "1\t300\t.\tC\tG,T\t50\tPASS\t.\tGT\t1/2",  # unsplit
        "1\t400\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
        "1\t400\t.\tA\t<*>\t.\t.\tEND=410\tGT\t0/0",  # a block at a call's position: not a call row
        "1\t500\t.\tA\tG\t50\tPASS\t.\tGT\t1",  # haploid: the all-zero rule is any ploidy
        "1\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0",
        "1\t600\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",  # two het rows sharing one slot: not the pattern
        "1\t600\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1",
        "1\t700\t.\tA\tG\t50\tPASS\t.\tGT\t0|1",  # a real SNV and indel on two haplotypes look split
        "1\t700\t.\tAT\tA\t50\tPASS\t.\tGT\t1|0",
        "1\t800\t.\tA\tG\t50\tPASS\t.\tGT\t0/0",  # hom-ref at a split site: every 0 there is REF
        "1\t800\t.\tA\tT\t50\tPASS\t.\tGT\t0/0",
    ]
    path = write(tmp_path, vcf_text(records, info=("END",)))
    calls, blocks, warnings = frames(path)
    assert calls.schema == CALLS and calls.height == len(records) - 1 and blocks.height == 1
    assert_gt_consistent(calls)
    signals = {(row["pos"], row["alts"][0]): row["split_signal"] for row in calls.iter_rows(named=True)}
    assert signals == {
        (100, "G"): "slot_pattern", (100, "T"): "slot_pattern",
        (200, "G"): "all_zero_beside_alt", (200, "T"): "all_zero_beside_alt",
        (300, "G"): None, (400, "G"): None,
        (500, "G"): "all_zero_beside_alt", (500, "T"): "all_zero_beside_alt",
        (600, "G"): None, (600, "AT"): None,
        (700, "G"): "slot_pattern", (700, "A"): "slot_pattern",
        (800, "G"): None, (800, "T"): None,
    }
    counted = [re.match(r"(\d+) row\(s\) at (\d+) position\(s\) .* \((\w+)\):", w) for w in warnings if "RM270" in w]
    assert {m.group(3): (int(m.group(1)), int(m.group(2))) for m in counted} == {
        "slot_pattern": (4, 2), "all_zero_beside_alt": (4, 2),
    }
    assert load(path).facets.collect()["split_marker"][0] is None


def test_split_markers_flag_their_rows(tmp_path: Path) -> None:
    """A bcftools `--old-rec-tag` key, whatever its name, marks a row only with its USED_ALT_IDX field;
    vt's OLD_MULTIALLELIC marks every row it is on. Both are read, so neither is a dropped field."""
    bcftools = '##INFO=<ID=ORIGINAL,Number=1,Type=String,Description="Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX">'
    vt = '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">'
    records = [
        "1\t100\t.\tA\tG\t50\tPASS\tORIGINAL=1|100|A|G,T|1\tGT\t0/1",  # a split 0/1, which no pattern fits
        "1\t100\t.\tA\tT\t50\tPASS\tORIGINAL=1|100|A|G,T|2\tGT\t0/0",
        "1\t200\t.\tC\tG\t50\tPASS\tORIGINAL=1|199|TC|TG\tGT\t0/1",  # realigned only: four fields
        "1\t300\t.\tG\tA\t50\tPASS\tOLD_MULTIALLELIC=1:300:G/A/C\tGT\t1/.",
        "1\t400\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
    ]
    path = write(tmp_path, vcf_text(records, extra_meta=(bcftools, vt)))
    loaded = load(path)
    calls = loaded.calls.collect()
    signals = {row["pos"]: row["split_signal"] for row in calls.iter_rows(named=True) if row["alts"] != ["T"]}
    assert signals == {100: "marker", 200: None, 300: "marker", 400: None}
    assert calls.filter(pl.col("alts") == ["T"])["split_signal"].to_list() == ["marker"]
    facet = loaded.facets.collect().row(0, named=True)
    assert facet["split_marker"] == "OLD_MULTIALLELIC,ORIGINAL"
    assert not {"INFO/ORIGINAL", "INFO/OLD_MULTIALLELIC"} & set(facet["dropped_fields"] or [])
    assert any(w.startswith("3 row(s) tagged as split") for w in loaded.warnings)


def test_a_declared_marker_that_tags_nothing_leaves_the_facet_null(tmp_path: Path) -> None:
    bcftools = '##INFO=<ID=ORIGINAL,Number=1,Type=String,Description="Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX">'
    records = ["1\t200\t.\tC\tG\t50\tPASS\tORIGINAL=1|199|TC|TG\tGT\t0/1"]  # realigned, not split
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=(bcftools,))))
    assert loaded.facets.collect()["split_marker"][0] is None
    assert loaded.calls.collect()["split_signal"].to_list() == [None] and loaded.warnings == ()


def test_vcf_no_meta_loads_under_a_synthesised_header(tmp_path: Path) -> None:
    """A `#CHROM` line with no meta lines (`vcf_no_meta`): the carried keys the rows use are declared
    by the loader, the rest are listed as dropped, and the rows equal the same body under a header."""
    body = [
        "chr1\t752566\t.\tG\tA\t50.2\t.\tAB=0.5;SF=0\tGT:RO:AO:GQ\t0/1:7:5:40",
        "chr1\t776546\t.\tA\tG\t.\tPASS\tADP=28;SF=1\tGT:FREQ\t0/1:28.57%",
        "chr2\t100\t.\tC\t<*>\t.\t.\tEND=180;SF=0\tGT:GQ\t0/0:30",
    ]
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1"
    bare = write(tmp_path, "\n".join([columns, *body]) + "\n", "bare.vcf")
    from just_vcf import sniff

    assert sniff(bare).format == "vcf_no_meta"
    loaded = load(bare)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    facet = loaded.facets.collect().row(0, named=True)
    assert calls["gt"].to_list() == [[0, 1], [0, 1]] and calls["gq"].to_list() == [40, None]
    assert blocks.select("chrom", "pos", "end_pos", "kind").rows() == [("2", 100, 180, "ref")]
    assert (facet["format_id"], facet["vcf_version"], facet["header_contigs"]) == ("vcf_no_meta", None, [])
    assert facet["dropped_fields"] == ["INFO/AB", "INFO/ADP", "INFO/SF", "FORMAT/AO", "FORMAT/FREQ", "FORMAT/RO"]
    assert any(w.startswith("no meta lines") for w in loaded.warnings)
    assert any("INFO/END" in w and "FORMAT/GT" in w and "not declared" in w for w in loaded.warnings)

    headed = write(tmp_path, vcf_text([], fmt=("GT", "GQ"), info=("END",)).split("#CHROM")[0] + columns + "\n" + "\n".join(body) + "\n", "headed.vcf")
    reference = load(headed)
    ref_calls, ref_blocks = pl.collect_all([reference.calls, reference.blocks])
    assert calls.drop("facet_id").equals(ref_calls.drop("facet_id"))
    assert blocks.drop("facet_id").equals(ref_blocks.drop("facet_id"))


def test_an_undeclared_end_is_declared_so_a_block_keeps_its_extent(tmp_path: Path) -> None:
    records = ["1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ\t0/0:60", "1\t200\t.\tA\t.\t.\tPASS\tEND=250\tGT\t0/0"]
    loaded = load(write(tmp_path, vcf_text(records, fmt=("GT", "GQ"))))
    blocks = loaded.blocks.collect()
    assert blocks.select("pos", "end_pos").rows() == [(100, 113), (200, 250)]
    assert any("INFO/END" in w and "not declared" in w for w in loaded.warnings)


def test_a_dash_alt_is_read_as_the_empty_allele_and_counted(tmp_path: Path) -> None:
    # F12: ALT `-` is an unanchored deletion of the REF bases; the loader has no reference to anchor it.
    # A zero-length allele has one spelling, "" (as CGI's native deletions and VRS's empty literal).
    records = [
        "1\t100\t.\tCT\t-\t.\t.\t.\tGT\t0/1",
        "1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1",
        "1\t300\t.\tG\t-\t.\t.\t.\tGT\t0/0",
        "1\t400\t.\tTA\tT,-\t.\t.\t.\tGT\t1/2",
    ]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records)))
    rows = {row["pos"]: row for row in calls.iter_rows(named=True)}
    assert (rows[100]["source_ref"], rows[100]["alts"], rows[100]["end_pos"]) == ("CT", [""], 101)
    assert (rows[100]["gt_alleles"], rows[100]["source_gt"]) == (["CT", ""], "0/1"), "source_gt keeps the file's text"
    assert (rows[400]["alts"], rows[400]["gt_alleles"]) == (["T", ""], ["T", ""])
    assert not any("-" in (row["alts"] or []) for row in rows.values())
    assert_gt_consistent(calls)
    assert any(w.startswith("3 record(s) write a deletion as ALT `-`") for w in warnings), warnings


def test_vectors_not_parallel_to_their_alleles_are_withheld(tmp_path: Path) -> None:
    """`ad` is one value per allele and `gl`/`gp` one per genotype. A row whose vector has another
    length is withheld and counted, never truncated; a haploid row needs one value per allele."""
    records = [
        "1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t0/1:5,4:30,0,60",  # fits
        "1\t200\t.\tC\tT\t.\t.\t.\tGT:AD:PL\t0/1:4:30,0",  # VarScan2-style ALT depth; two PL for three genotypes
        "1\t300\t.\tG\tA,C\t.\t.\t.\tGT:AD:PL\t1/2:0,3,4:90,40,0,50,20,80",  # fits, multi-allelic
        "1\t400\t.\tT\tC\t.\t.\t.\tGT:AD:PL\t1:0,9:90,0",  # haploid: two genotypes
        "1\t500\t.\tT\tC\t.\t.\t.\tGT:AD:GP\t0/1:3,3:0.1,0.8",  # GP short by one
    ]
    loaded = load(write(tmp_path, vcf_text(records, fmt=("GT", "AD", "PL", "GP"))))
    calls = loaded.calls.collect()
    rows = {row["pos"]: row for row in calls.iter_rows(named=True)}
    assert rows[100]["ad"] == [5, 4] and rows[100]["gl"] == pytest.approx([-3.0, 0.0, -6.0])
    assert rows[200]["ad"] is None and rows[200]["gl"] is None and rows[200]["source_alleles"] is None
    assert rows[300]["ad"] == [0, 3, 4] and len(rows[300]["gl"]) == 6
    assert rows[400]["ad"] == [0, 9] and rows[400]["gl"] == pytest.approx([-9.0, 0.0])
    assert rows[500]["gp"] is None and rows[500]["ad"] == [3, 3]
    assert set(loaded.warnings) >= {
        "1 AD value(s) not one per allele set to null",
        "1 PL value(s) not one per genotype set to null",
        "1 GP value(s) not one per genotype set to null",
    }


def test_pl_and_the_same_likelihood_as_gl_load_equal(tmp_path: Path) -> None:
    """F11: PL `p` and GL `-p/10` are one likelihood, so they load as the same Float32: the one nearest
    to `-p/10`, as a GL file's own value is parsed."""
    values = list(range(0, 300))
    triples = [values[i : i + 3] for i in range(0, len(values), 3)]
    pl_records = [f"1\t{100 + i}\t.\tA\tG\t.\t.\t.\tGT:PL\t0/1:{','.join(map(str, t))}" for i, t in enumerate(triples)]
    gl_records = [f"1\t{100 + i}\t.\tA\tG\t.\t.\t.\tGT:GL\t0/1:{','.join(repr(-p / 10) for p in t)}" for i, t in enumerate(triples)]
    gl_line = '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods">'
    from_pl, _, _ = frames(write(tmp_path, vcf_text(pl_records, fmt=("GT", "PL")), "pl.vcf"))
    from_gl, _, _ = frames(write(tmp_path, vcf_text(gl_records, fmt=("GT",), extra_meta=(gl_line,)), "gl.vcf"))
    assert from_pl["gl"].to_list() == from_gl["gl"].to_list()
    nearest = [[struct.unpack("f", struct.pack("f", -p / 10))[0] for p in t] for t in triples]
    assert from_pl["gl"].to_list() == nearest


def test_a_vector_key_declared_with_a_fixed_count_is_not_carried(tmp_path: Path) -> None:
    """VarScan2 declares AD as Number=1 (the ALT depth) and CGI as Number=2 (per GT slot): neither is
    one value per allele, so `ad` stays null and the key is listed as dropped."""
    extra = ('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)">',)
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT:AD\t1/1:3,3"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=extra)))
    assert loaded.calls.collect()["ad"].to_list() == [None]
    assert "FORMAT/AD" in loaded.facets.collect()["dropped_fields"][0].to_list()
    assert any("not carried: AD (Number=2)" in w for w in loaded.warnings)


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.5"])
def test_symbolic_sv_end_pos_reads_end_then_svlen(tmp_path: Path, version: str) -> None:
    """SCHEMA.md `end_pos`: a length-type symbolic SV ends at INFO END where the record gives one, else at
    POS + |SVLEN|; `<INS>`, breakends and other symbolic alleles stay at the anchor, whatever END says."""
    records = [
        "1\t200\t.\tT\t<DEL>\t.\t.\tEND=230\tGT\t0/1",  # END only
        "1\t300\t.\tT\t<DUP:TANDEM>\t.\t.\tSVLEN=12\tGT\t0/1",  # SVLEN only
        "1\t400\t.\tG\t<DEL>\t.\t.\tEND=450;SVLEN=-50\tGT\t0/1",  # both, agreeing
        "1\t500\t.\tA\t<CNV>\t.\t.\tEND=600;SVLEN=101\tGT\t0/1",  # both, disagreeing: END wins, counted
        "1\t700\t.\tC\t<INS>\t.\t.\tEND=760;SVLEN=60\tGT\t0/1",  # an insertion spans no reference
        "1\t800\t.\tG\t<TRA>\t.\t.\tEND=5000\tGT\t0/1",  # Delly: END is the partner's position on CHR2
        "1\t900\t.\tA\t<DEL>,<INS>\t.\t.\tSVLEN=10,4\tGT\t1/2",  # the largest end over the site
        "1\t1000\t.\tC\t<*>\t.\t.\tEND=1050\tGT\t0/0",  # a gVCF block, unchanged
        "1\t1100\t.\tCTTT\tC\t.\t.\t.\tGT\t0/1",  # a sequence allele: the REF span
    ]
    path = write(tmp_path, vcf_text(records, version=version, info=("END", "SVLEN")))
    calls, blocks, warnings = frames(path)
    ends = dict(zip(calls["pos"].to_list(), calls["end_pos"].to_list(), strict=True))
    assert ends == {200: 230, 300: 312, 400: 450, 500: 600, 700: 700, 800: 800, 900: 910, 1100: 1103}
    assert blocks["end_pos"].to_list() == [1050]
    assert "1 symbolic SV record(s) whose INFO END disagrees with POS + |SVLEN|: end_pos follows END" in warnings
    assert any(w.startswith("2 record(s) carry an INFO END that is not their extent") for w in warnings), warnings


@pytest.mark.parametrize("number", ["1", "."])
def test_svlen_declared_with_another_number_still_sets_end_pos(tmp_path: Path, number: str) -> None:
    """Older SV callers declare SVLEN `Number=1` (polars-bio then returns a scalar) or `Number=.`."""
    records = ["1\t200\t.\tT\t<DEL>\t.\t.\tSVLEN=-30\tGT\t0/1", "1\t300\t.\tT\t<DUP>\t.\t.\tSVLEN=12\tGT\t0/1"]
    svlen = f'##INFO=<ID=SVLEN,Number={number},Type=Integer,Description="SV length">'
    calls, _, _ = frames(write(tmp_path, vcf_text(records, extra_meta=(svlen,))))
    assert calls["end_pos"].to_list() == [230, 312]


def test_symbolic_sv_extent_survives_its_own_export(tmp_path: Path) -> None:
    """`write_vcf` writes a symbolic SV's extent as END and no SVLEN, so the reload has to read END."""
    records = ["1\t200\t.\tT\t<DEL>\t.\t.\tSVLEN=-30\tGT\t0/1", "1\t300\t.\tA\tG\t.\t.\t.\tGT\t0/1"]
    first = load(write(tmp_path, vcf_text(records, info=("SVLEN",))))
    back = load(write_vcf(first, tmp_path / "out.vcf.gz").path).calls.collect()
    assert back["end_pos"].to_list() == first.calls.collect()["end_pos"].to_list() == [230, 300]


def test_a_column_line_written_with_two_hashes_is_read(tmp_path: Path) -> None:
    """F7: Veritas SV and CNV uploads write `##CHROM` and put a third `#` on most meta lines. The
    column line is read; a `###` line is a meta line whose key starts with `#` and declares nothing."""
    text = "\n".join(
        [
            "##fileformat=VCFv4.1",
            '###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "###contig=<ID=chr1,length=249250621,assembly=hg19>",
            "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            "chr1\t100\t.\tC\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=150\tGT\t0/1",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1/1",
        ]
    ) + "\n"
    loaded = load(write(tmp_path, text, "sv.vcf"))
    calls = loaded.calls.collect()
    facet = loaded.facets.collect().row(0, named=True)
    assert calls.select("chrom", "pos", "end_pos", "alts", "gt_alleles").rows() == [
        ("1", 100, 150, ["<DEL>"], ["C", "<DEL>"]),
        ("1", 300, 300, ["G"], ["G", "G"]),
    ]
    assert facet["header_contigs"] == [], "a ###contig line is not a contig declaration"
    assert facet["sample"] == "S1"
    assert any(w.startswith("the column line is written ##CHROM") and "2 meta line(s) starting ###" in w for w in loaded.warnings)
    assert any("FORMAT/GT" in w for w in loaded.warnings), "GT is undeclared, so the loader declares it"


def test_a_column_line_after_the_real_one_is_not_reread(tmp_path: Path) -> None:
    text = vcf_text(["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"]).replace(
        "#CHROM", "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOTHER\n#CHROM", 1
    )
    loaded = load(write(tmp_path, text))
    assert loaded.facets.collect()["sample"].to_list() == ["S1"]
    assert not any("##CHROM" in w for w in loaded.warnings)


def test_an_alt_that_is_no_vcf_allele_has_no_variant_class(tmp_path: Path) -> None:
    """F15: a Veritas CNV upload writes ALT `CNV`. It is kept verbatim and counted, never an indel."""
    records = ["1\t100\t.\tG\tCNV\t.\t.\t.\tGT\t0/1", "1\t200\t.\tA\tAT\t.\t.\t.\tGT\t0/1", "1\t300\t.\tA\tR\t.\t.\t.\tGT\t0/1"]
    loaded = load(write(tmp_path, vcf_text(records)))
    calls = loaded.calls.collect()
    assert calls["alts"].to_list() == [["CNV"], ["AT"], ["R"]]
    assert loaded.facets.collect()["variant_classes"].to_list() == [[{"variant_class": "indel", "alleles": 1}]]
    assert any(w.startswith("2 record(s) carry an ALT that is no VCF allele") for w in loaded.warnings), loaded.warnings


def _spreadsheet_resave(text: str) -> bytes:
    """`text` as a spreadsheet re-saves a tab file: CSV quoting of fields holding a comma or a quote,
    lines padded with empty fields to the widest row, bare-CR line ends. The test's own encoder."""
    rows = [line.split("\t") for line in text.rstrip("\n").split("\n")]
    width = max(len(row) for row in rows)
    out = []
    for row in rows:
        cells = ['"' + cell.replace('"', '""') + '"' if ("," in cell or '"' in cell) else cell for cell in row]
        out.append("\t".join(cells + [""] * (width - len(cells))))
    return "\r".join(out).encode() + b"\r"


def test_a_vcf_resaved_by_a_spreadsheet_loads_as_the_text_it_was_saved_from(tmp_path: Path) -> None:
    """F7: a VCF re-saved by a spreadsheet: quoted `"##…""…""` meta lines padded with tabs, quoted
    INFO and sample fields, bare-CR line ends."""
    records = [
        '1\t100\trs1\tA\tG\t50\tPASS\tEFF=A(B|C),D(E);NOTE="x"\tGT:AD:PL\t0/1:10,12:40,0,50',
        "1\t200\t.\tC\tCT,T\t60\tLowQual\tEFF=F\tGT:AD:PL\t1/2:1,5,6:90,30,80,20,0,70",
        "2\t300\t.\tG\t.\t.\t.\t.\tGT:AD:PL\t./.",
    ]
    meta = (
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##INFO=<ID=EFF,Number=.,Type=String,Description="Effects, as \'A(B|C)\'">',
        '##INFO=<ID=NOTE,Number=1,Type=String,Description="Note">',
    )
    plain = vcf_text(records, fmt=("GT", "AD", "PL"), extra_meta=meta)
    resaved = write(tmp_path, "", "resaved.vcf")
    resaved.write_bytes(_spreadsheet_resave(plain))
    assert b'"##INFO=<ID=EFF' in resaved.read_bytes() and b"\n" not in resaved.read_bytes()

    expected, _, _ = frames(write(tmp_path, plain, "plain.vcf"))
    loaded = load(resaved)
    calls = loaded.calls.collect()
    facet = loaded.facets.collect().row(0, named=True)
    columns = ["chrom", "pos", "ids", "source_ref", "alts", "gt_alleles", "ad", "gl", "qual", "filters"]
    assert calls.select(columns).rows() == expected.select(columns).rows()
    assert calls.filter(pl.col("pos") == 100)["ad"].to_list() == [[10, 12]]
    assert "INFO/EFF" in facet["dropped_fields"] and "INFO/NOTE" in facet["dropped_fields"]
    assert vcf_samples(resaved) == ["S1"]
    assert any(w.startswith("the file was re-saved by a spreadsheet") for w in loaded.warnings)


def test_a_resaved_vcf_with_a_field_past_the_csv_limit_loads(tmp_path: Path) -> None:
    """F7's re-save with one INFO field of about 200 KiB, quoted for its commas: Python's `csv` stops at
    131,072 bytes a field, a VCF field has no such limit."""
    effects = ",".join(f"EFFECT_{i:06d}(MODIFIER|gene{i})" for i in range(7000))
    assert len(effects) > 200_000 > csv.field_size_limit()
    records = [
        f"1\t100\trs1\tA\tG\t50\tPASS\tEFF={effects}\tGT:AD\t0/1:10,12",
        "1\t200\t.\tC\tT\t60\tPASS\tEFF=F\tGT:AD\t1/1:0,9",
    ]
    meta = ('##INFO=<ID=EFF,Number=.,Type=String,Description="Effects, as \'A(B|C)\'">',)
    plain = vcf_text(records, fmt=("GT", "AD"), extra_meta=meta)
    resaved = write(tmp_path, "", "resaved.vcf")
    resaved.write_bytes(_spreadsheet_resave(plain))
    assert f'"EFF={effects}"'.encode() in resaved.read_bytes()

    expected, _, _ = frames(write(tmp_path, plain, "plain.vcf"))
    loaded = load(resaved)
    columns = ["chrom", "pos", "ids", "source_ref", "alts", "gt_alleles", "ad", "qual", "filters"]
    assert loaded.calls.collect().select(columns).rows() == expected.select(columns).rows()
    assert any(w.startswith("the file was re-saved by a spreadsheet") for w in loaded.warnings)


def test_undo_resave_reads_quotes_as_csv_does(tmp_path: Path) -> None:
    """The unquoting against `csv` itself, on the shapes a spreadsheet writes (fields under its limit)."""
    lines = [
        '"##INFO=<ID=EFF,Number=.,Type=String,Description=""Effects, as \'A(B|C)\'"">"\t\t\t',
        "#CHROM\tPOS\tID",
        '1\t"NOTE=""x"";EFF=A,B"\t""\t"a""""b"\tplain',
        '2\t"""quoted"""\t.',
    ]
    text = "\r".join(lines) + "\r"
    ours = undo_resave(io.BytesIO(text.encode())).read().decode().split("\n")[:-1]
    theirs = []
    for row in csv.reader(io.StringIO(text, newline=None), delimiter="\t", quotechar='"'):
        if row and row[0].startswith("##"):
            while len(row) > 1 and row[-1] == "":
                row.pop()
        theirs.append("\t".join(row))
    assert ours == theirs
    assert ours[2] == '1\tNOTE="x";EFF=A,B\t\ta""b\tplain'


def test_crlf_and_plain_vcfs_are_not_taken_for_a_resave(tmp_path: Path) -> None:
    text = vcf_text(["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"])
    for name, body in (("lf.vcf", text.encode()), ("crlf.vcf", text.replace("\n", "\r\n").encode())):
        path = tmp_path / name
        path.write_bytes(body)
        assert not any("spreadsheet" in w for w in load(path).warnings), name


GL_FLOAT = '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods">'


@pytest.mark.parametrize(
    ("fmt", "extra", "sample", "expected"),
    [
        (("GT", "PL"), (), "0/1:40,0,300", "pl"),
        (("GT",), (GL_FLOAT,), "0/1:-4,0,-30", "gl"),
        (("GT", "PL"), (GL_FLOAT,), "0/1:40,0,300:-4,0,-30", "gl_over_pl"),
        (("GT", "GQ"), (), "0/1:30", None),
        # A GL with a fixed count is another quantity under the name, not carried: PL is the only source.
        (("GT", "PL"), ('##FORMAT=<ID=GL,Number=1,Type=Float,Description="Something else">',), "0/1:40,0,300:-4", "pl"),
    ],
    ids=["pl", "gl", "both", "neither", "gl_fixed_count"],
)
def test_the_facet_names_the_key_gl_is_read_from(tmp_path: Path, fmt: tuple[str, ...], extra: tuple[str, ...], sample: str,
                                                 expected: str | None) -> None:  # fmt: skip
    keys = ":".join([*fmt, *("GL" for line in extra if "ID=GL" in line)])
    path = write(tmp_path, vcf_text([f"1\t100\t.\tA\tG\t50\t.\t.\t{keys}\t{sample}"], fmt=fmt, extra_meta=extra))
    loaded = load(path)
    facets, calls = pl.collect_all([loaded.facets, loaded.calls])
    assert facets["gl_source"].to_list() == [expected]
    # Null means the facet carries no likelihoods: no row has `gl`.
    assert (calls["gl"].is_null().all()) is (expected is None)


def test_an_undeclared_pl_the_rows_use_is_a_pl_facet(tmp_path: Path) -> None:
    """The loader declares a key the rows use without a declaration, so it is read, and named."""
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t50\t.\t.\tGT:PL\t0/1:40,0,300"]))
    loaded = load(path)
    assert loaded.facets.collect()["gl_source"].to_list() == ["pl"]
    assert loaded.calls.collect()["gl"].to_list()[0] == pytest.approx([-4.0, 0.0, -30.0])


def test_a_derived_facet_keeps_its_parents_gl_source(tmp_path: Path) -> None:
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t50\t.\t.\tGT:PL\t0/1:40,0,300"], fmt=("GT", "PL")))
    loaded = load(path)
    derived = rejoin_split(loaded).facets.collect()
    assert derived["relation"].to_list() == ["split_rejoined"]
    assert derived["gl_source"].to_list() == loaded.facets.collect()["gl_source"].to_list() == ["pl"]


# --- FORMAT HQ, haplotype qualities, as `hq` -------------------------------------------------------

HQ_LINE = '##FORMAT=<ID=HQ,Number={number},Type=Integer,Description="Haplotype quality">'

# (GT:HQ as the file writes it, the `hq` the frame holds, what the row adds to the counts). One row
# per rule in SCHEMA.md's "A VCF's haplotype qualities": as written where one value per GT allele; a
# Number=2 field on a haploid row with its second value missing (CGI writes `51,.`) read for the one
# allele; every other count withheld; a value outside u16 or on an uncalled allele left out; null
# where no value is left.
HQ_ROWS = [
    ("0/1:51,40", [51, 40], None),
    ("0|0:48,47", [48, 47], None),
    ("1:51,.", [51], None),
    ("1:51", [51], None),
    ("1/1:.", None, None),
    ("1/1:.,.", None, None),
    ("1/1:.,33", [None, 33], None),
    ("1/1:70000,-3", None, "out"),
    ("0/1:70000,30", [None, 30], "out"),
    ("./1:20,30", [None, 30], "uncalled"),
    ("1:51,40", None, "unfit"),
    ("0/1:30", None, "unfit"),
    ("1/1/1:5,6,7", [5, 6, 7], None),
]


def _hq_records(rows: list[tuple]) -> list[str]:
    return [f"1\t{100 + i}\t.\tA\tG\t50\t.\t.\tGT:HQ\t{row[0]}" for i, row in enumerate(rows)]


def test_hq_is_read_per_gt_allele(tmp_path: Path) -> None:
    extra = (HQ_LINE.format(number=2),)
    path = write(tmp_path, vcf_text([*_hq_records(HQ_ROWS), "1\t200\t.\tA\tG\t50\t.\t.\tGT\t1/1"], extra_meta=extra))
    loaded = load(path)
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls.schema["hq"] == pl.List(pl.UInt16)
    by_pos = dict(zip(calls["pos"].to_list(), calls["hq"].to_list()))
    assert [by_pos[100 + i] for i in range(len(HQ_ROWS))] == [hq for _, hq, _ in HQ_ROWS]
    assert by_pos[200] is None, "a row without HQ has no score"
    # Every element left is on a called allele: the frame never scores a `.`.
    for row in calls.filter(pl.col("hq").is_not_null()).iter_rows(named=True):
        assert len(row["hq"]) == len(row["gt"])
        assert all(score is None for score, allele in zip(row["hq"], row["gt"]) if allele is None)
    values = [v for sample, _, kind in HQ_ROWS if kind == "out" for v in sample.split(":")[1].split(",")]
    out = sum(not 0 <= int(v) <= 2**16 - 1 for v in values if v != ".")
    unfit = sum(kind == "unfit" for *_, kind in HQ_ROWS)
    assert set(loaded.warnings) == {
        f"{unfit} HQ value(s) not one per GT allele set to null",
        f"{out} HQ value(s) outside 0..65535 set to null",
        "1 HQ value(s) on an uncalled GT allele not read",
    }
    assert "FORMAT/HQ" not in (facet["dropped_fields"] or [])


@pytest.mark.parametrize(("number", "version"), [(".", "VCFv4.2"), ("P", "VCFv4.4")])
def test_hq_declared_with_any_number_or_one_per_allele_is_read(tmp_path: Path, number: str, version: str) -> None:
    rows = [("0/1:51,40", [51, 40]), ("1:51", [51]), ("1/1/1:5,6,7", [5, 6, 7]), ("1:51,40", None)]
    path = write(tmp_path, vcf_text(_hq_records(rows), version=version, extra_meta=(HQ_LINE.format(number=number),)))
    loaded = load(path)
    assert loaded.calls.collect().sort("pos")["hq"].to_list() == [hq for _, hq in rows]
    assert "1 HQ value(s) not one per GT allele set to null" in loaded.warnings


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3", "VCFv4.5"])
@pytest.mark.parametrize(("number", "value"), [("1", "51"), ("1", "51,40"), ("G", "51,40,7"), ("R", "51,40,7"), ("3", "51,40,7")])
def test_hq_declared_with_another_number_is_not_carried(tmp_path: Path, number: str, value: str, version: str) -> None:
    """A fixed count other than VCF's 2, or a per-allele or per-genotype one, is another quantity under
    the name: not read into `hq`, listed in `dropped_fields`, and the load does not fail on it, even
    where the row holds more values than the Number allows (`Number=1` with `51,40`), which polars-bio
    cannot parse under that declaration."""
    record = f"1\t100\t.\tA\tG\t50\t.\t.\tGT:HQ\t0/1:{value}"
    path = write(tmp_path, vcf_text([record], version=version, extra_meta=(HQ_LINE.format(number=number),)))
    loaded = load(path)
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls["hq"].to_list() == [None] and calls["gt"].to_list() == [[0, 1]]
    assert "FORMAT/HQ" in facet["dropped_fields"]
    assert any(f"HQ (Number={number})" in w and "not carried" in w for w in loaded.warnings)


@pytest.mark.parametrize(
    ("kind", "rows"),
    [
        # (GT:HQ, the `hq` read, why a value is left out). A whole number is read as it is; one with a
        # fraction is never rounded or cut to fit u16, and text that is no number is no score.
        ("Float", [("0/1:30.5,40", [None, 40], "fraction"), ("0/1:30.0,40.0", [30, 40], None), ("0/1:.,40.9", None, "fraction"),
                   ("0/1:nan,40", [None, 40], "text"), ("./1:30.5,40", [None, 40], "uncalled"), ("0/1:1e30,40", [None, 40], "out")]),
        ("String", [("0/1:30,x", [30, None], "text"), ("0/1:30.5,40", [None, 40], "fraction"), ("0/1:30,40", [30, 40], None),
                    ("0/1:.,x", None, "text"), ("./1:x,40", [None, 40], "uncalled"), ("0/1:100000000000000000000,40", [None, 40], "out")]),
    ],
)
def test_hq_that_is_not_a_whole_number_is_left_out_and_counted(tmp_path: Path, kind: str, rows: list[tuple]) -> None:
    """VCF 4.2 reserves nothing, so polars-bio reads HQ under whatever Type the header declares. A
    score is a whole phred number: a value with a fraction is left out, not truncated, and text that is
    no number is left out too, each reason counted in its own warning."""
    meta = (f'##FORMAT=<ID=HQ,Number=2,Type={kind},Description="Haplotype quality">',)
    loaded = load(write(tmp_path, vcf_text(_hq_records(rows), extra_meta=meta)))
    assert loaded.calls.collect().sort("pos")["hq"].to_list() == [hq for _, hq, _ in rows]
    counts = {reason: sum(1 for *_, why in rows if why == reason) for reason in ("fraction", "text", "uncalled", "out")}
    hq_warnings = {w for w in loaded.warnings if "HQ" in w}
    assert hq_warnings == {
        f"{counts['out']} HQ value(s) outside 0..65535 set to null",
        f"{counts['fraction']} HQ value(s) not a whole number left out",
        f"{counts['text']} HQ value(s) not a number left out",
        f"{counts['uncalled']} HQ value(s) on an uncalled GT allele not read",
    }


_PSL_LINE = '##FORMAT=<ID=PSL,Number=P,Type=String,Description="Phase set list">'
_END_LINE = INFO_LINES["END"]

# Each carried key VCF types Integer that the frames read as whole numbers, other than HQ: (key, its
# namespace and Number, other meta lines the record needs, the record with value `v` at `pos`, the frame
# and column it lands in, the value read for a whole `n`, and the value there when `v` is left out).
_WHOLE_KEYS = [
    ("GQ", "FORMAT", "1", (), lambda p, v: f"1\t{p}\t.\tA\tG\t50\t.\t.\tGT:GQ\t0/1:{v}", "calls", "gq", lambda p, n: n, lambda p: None),
    ("DP", "FORMAT", "1", (), lambda p, v: f"1\t{p}\t.\tA\tG\t50\t.\t.\tGT:DP\t0/1:{v}", "calls", "dp", lambda p, n: n, lambda p: None),
    ("AD", "FORMAT", "R", (), lambda p, v: f"1\t{p}\t.\tA\tG\t50\t.\t.\tGT:AD\t0/1:{v},4", "calls", "ad", lambda p, n: [n, 4], lambda p: [None, 4]),
    ("PSO", "FORMAT", "P", (_PSL_LINE,), lambda p, v: f"1\t{p}\t.\tA\tG\t50\t.\t.\tGT:PSL:PSO\t0|1:a,a:{v},4", "calls", "phase_set_ordinals",
     lambda p, n: [n, 4], lambda p: [None, 4]),
    ("MIN_DP", "FORMAT", "1", (_END_LINE,), lambda p, v: f"1\t{p}\t.\tA\t<*>\t.\t.\tEND={p + 5}\tGT:MIN_DP\t0/0:{v}", "blocks", "dp_min",
     lambda p, n: n, lambda p: None),
    ("LEN", "FORMAT", "1", (), lambda p, v: f"1\t{p}\t.\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:{v}", "blocks", "end_pos", lambda p, n: p + n - 1, lambda p: p),
    ("END", "INFO", "1", (), lambda p, v: f"1\t{p}\t.\tA\t<DEL>\t50\t.\tEND={p}{v}\tGT\t0/1", "calls", "end_pos", lambda p, n: int(f"{p}{n}"), lambda p: p),
    ("SVLEN", "INFO", "A", (), lambda p, v: f"1\t{p}\t.\tA\t<DEL>\t50\t.\tSVLEN=-{v}\tGT\t0/1", "calls", "end_pos", lambda p, n: p + n, lambda p: p),
    ("CIPOS", "INFO", ".", (), lambda p, v: f"1\t{p}\t.\tA\t<DEL>\t50\t.\tCIPOS={v},5\tGT\t0/1", "calls", "cipos", lambda p, n: [n, 5], lambda p: [None, 5]),
]

# (Type, the values written, each one's whole number or the reason it is left out). polars-bio parses a
# Float, so text there fails the read; under String, text reaches the loader.
_WHOLE_VALUES = {
    "Float": [("30.5", "not a whole number"), ("30.0", 30), ("31", 31)],
    "String": [("30.5", "not a whole number"), ("x", "not a number"), ("31", 31)],
}


@pytest.mark.parametrize("kind", ["Float", "String"])
@pytest.mark.parametrize(("key", "namespace", "number", "meta", "record", "frame", "column", "read", "unread"), _WHOLE_KEYS)
def test_an_integer_key_declared_float_or_string_keeps_only_whole_numbers(
    tmp_path: Path, kind: str, key: str, namespace: str, number: str, meta: tuple, record, frame: str, column: str, read, unread
) -> None:
    """VCF 4.2 reserves no Type, so a file can declare GQ, DP, AD and the other integer keys Float or
    String. A whole number is read as it is (`30.0` is 30); a number with a fraction is left out, never
    truncated, and text that is no number is left out too, each reason one warning with its count."""
    declared = f'##{namespace}=<ID={key},Number={number},Type={kind},Description="x">'
    values = _WHOLE_VALUES[kind]
    positions = [1000 + 100 * i for i in range(len(values))]
    records = [record(p, v) for p, (v, _) in zip(positions, values)]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=(*meta, declared))))
    rows = (loaded.calls if frame == "calls" else loaded.blocks).collect().sort("pos")
    assert rows["pos"].to_list() == positions
    expected = [read(p, n) if isinstance(n, int) else unread(p) for p, (_, n) in zip(positions, values)]
    assert rows[column].to_list() == expected
    label = key if namespace == "FORMAT" else f"INFO {key}"
    reasons = [n for _, n in values if isinstance(n, str)]
    stated = {w for w in loaded.warnings if w.startswith(tuple(f"{i} {label} value(s) not a" for i in range(1, 9)))}
    assert stated == {f"{reasons.count(r)} {label} value(s) {r} left out" for r in set(reasons)}, loaded.warnings


@pytest.mark.parametrize(("namespace", "key", "record", "column"), [
    ("FORMAT", "PS", "1\t100\t.\tA\tG\t50\t.\t.\tGT:PS\t0|1:123456789", "phase_sets"),
    ("INFO", "END", "1\t100\t.\tA\t<DEL>\t50\t.\tEND=123456789\tGT\t0/1", "end_pos"),
])
def test_a_whole_number_a_float_cannot_hold_exactly_is_left_out(tmp_path: Path, namespace: str, key: str, record: str, column: str) -> None:
    """polars-bio reads a Float as a 32-bit float, which holds every whole number only up to 2^24:
    123456789 comes back 123456792. A position or a phase set read so would be one the file did not
    write, so such a value is left out and counted."""
    declared = f'##{namespace}=<ID={key},Number=1,Type=Float,Description="x">'
    loaded = load(write(tmp_path, vcf_text([record], extra_meta=(declared,))))
    (value,) = loaded.calls.collect()[column].to_list()
    assert value != 123456792 and value != ["123456792", "123456792"]
    label = key if namespace == "FORMAT" else f"INFO {key}"
    assert any(w.startswith(f"1 {label} value(s) too large to read exactly as a whole number left out") for w in loaded.warnings), loaded.warnings


@pytest.mark.parametrize(("kind", "value", "reason"), [("String", "abc", "not a number"), ("Float", "7.5", "not a whole number"), ("String", "7", None)])
def test_a_phase_set_that_is_not_read_withholds_the_phase(tmp_path: Path, kind: str, value: str, reason: str | None) -> None:
    """A PS left out does not put the record in its chromosome's implicit phase set, which a phased
    record without PS is in: the set it named is unknown, so its phase is withheld. An unphased record
    keeps its `phased` as written."""
    meta = (f'##FORMAT=<ID=PS,Number=1,Type={kind},Description="Phase set">',)
    records = [f"1\t100\t.\tA\tG\t50\t.\t.\tGT:PS\t0|1:{value}", f"1\t200\t.\tA\tG\t50\t.\t.\tGT:PS\t0/1:{value}"]
    loaded = load(write(tmp_path, vcf_text(records, extra_meta=meta)))
    calls = loaded.calls.collect().sort("pos")
    if reason is None:
        assert calls["phased"].to_list() == [[True, True], [False, False]] and calls["phase_sets"].to_list()[0] == ["7", "7"]
        return
    assert calls["phased"].to_list() == [[None, None], [False, False]]
    assert calls["phase_sets"].to_list() == [None, None]
    assert any(w.startswith(f"2 PS value(s) {reason} left out") and "phase" in w for w in loaded.warnings), loaded.warnings


@pytest.mark.parametrize("version", ["VCFv4.3", "VCFv4.4", "VCFv4.5"])
def test_hq_declared_with_a_type_its_version_reserves_otherwise_is_not_carried(tmp_path: Path, version: str) -> None:
    """VCF 4.3 on reserves HQ as Integer, and polars-bio refuses a header of those versions declaring it
    Float: the file loads without HQ, says so, and lists it in `dropped_fields`. Under 4.2, where
    nothing is reserved, the same declaration is read (a whole number as it is)."""
    meta = ('##FORMAT=<ID=HQ,Number=2,Type=Float,Description="Haplotype quality">',)
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tGT:HQ\t0/1:30.5,40.5", "1\t200\t.\tA\tG\t50\t.\t.\tGT:HQ\t0/1:30,40"]
    loaded = load(write(tmp_path, vcf_text(records, version=version, extra_meta=meta)))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls["hq"].to_list() == [None, None] and calls["gt"].to_list() == [[0, 1], [0, 1]]
    assert "FORMAT/HQ" in facet["dropped_fields"]
    reserved = version.removeprefix("VCFv")
    assert any(f"not carried: FORMAT/HQ (Type=Float, {reserved} reserves Integer)" in w for w in loaded.warnings), loaded.warnings
    older = load(write(tmp_path, vcf_text(records[1:], version="VCFv4.2", extra_meta=meta), "v42.vcf"))
    assert older.calls.collect()["hq"].to_list() == [[30, 40]]


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3", "VCFv4.5"])
@pytest.mark.parametrize(("key", "value", "column"), [("AD", "5,5", "ad"), ("PL", "0,10,20", "gl"), ("GL", "0,-1,-2", "gl"), ("GP", "0,0,1", "gp"), ("PP", "0,10,20", "gp")])
def test_a_vector_declared_one_value_with_several_is_not_carried(tmp_path: Path, key: str, value: str, column: str, version: str) -> None:
    """A carried vector declared `Number=1` names another quantity (`VECTOR_NUMBERS`), and a row writing
    several values there cannot be parsed under that declaration: the file still loads, without the key."""
    kind = "Float" if key in ("GL", "GP") else "Integer"
    meta = (f'##FORMAT=<ID={key},Number=1,Type={kind},Description="x">',)
    loaded = load(write(tmp_path, vcf_text([f"1\t100\t.\tA\tG\t50\t.\t.\tGT:{key}\t0/1:{value}"], version=version, extra_meta=meta)))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls[column].to_list() == [None] and calls["gt"].to_list() == [[0, 1]]
    assert f"FORMAT/{key}" in facet["dropped_fields"]
    assert any(f"{key} (Number=1)" in w and "not carried" in w for w in loaded.warnings), loaded.warnings


# Each carried numeric FORMAT key, a value it reads, and the `calls` column that holds it. Declared
# with the other numeric Type below, the value still parses, so only the header can refuse it.
_TYPED_KEYS = [
    ("GQ", "1", "30", "gq"),
    ("HQ", "2", "30,40", "hq"),
    ("DP", "1", "10", "dp"),
    ("AD", "R", "5,5", "ad"),
    ("PL", "G", "0,10,20", "gl"),
    ("GL", "G", "0,-1,-2", "gl"),
    ("GP", "G", "0,0,1", "gp"),
    ("PP", "G", "0,10,20", "gp"),
    ("DS", "A", "1", "dosage"),
    ("PS", "1", "7", "phase_sets"),
    ("CN", "1", "2", "copy_number"),
]


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3", "VCFv4.4", "VCFv4.5"])
@pytest.mark.parametrize(("key", "number", "value", "column"), _TYPED_KEYS)
def test_a_carried_key_declared_with_another_type_loads(tmp_path: Path, key: str, number: str, value: str, column: str, version: str) -> None:
    """polars-bio refuses a VCF 4.3+ header declaring a key its version reserves with another Type. The
    file loads all the same: the key is either read or not carried and listed in `dropped_fields`, never
    both, never neither. VCF 4.2 reserves nothing, so there every key is read."""
    other = "Integer" if FORMAT_KEYS[key][1] == "Float" else "Float"
    meta = (f'##FORMAT=<ID={key},Number={number},Type={other},Description="x">',)
    record = f"1\t100\t.\tA\tG\t50\t.\t.\tGT:{key}\t0|1:{value}"
    loaded = load(write(tmp_path, vcf_text([record], version=version, extra_meta=meta)))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    dropped = f"FORMAT/{key}" in (facet["dropped_fields"] or [])
    (value,) = calls[column].to_list()
    read = value is not None and (not isinstance(value, list) or any(v is not None for v in value))
    assert dropped != read, (value, facet["dropped_fields"])
    assert dropped == any(f"FORMAT/{key} (Type={other}" in w for w in loaded.warnings), loaded.warnings
    if version == "VCFv4.2":
        assert not dropped


@pytest.mark.parametrize("version", ["VCFv4.3", "VCFv4.4", "VCFv4.5"])
@pytest.mark.parametrize(("key", "number", "value"), [("END", "1", "200"), ("SVLEN", "A", "-100"), ("CIPOS", ".", "-1,1"), ("CIEND", ".", "-1,1")])
def test_a_carried_info_key_declared_with_another_type_loads(tmp_path: Path, key: str, number: str, value: str, version: str) -> None:
    """As for FORMAT: an INFO key the version reserves as Integer, declared Float, is not carried, and
    the file loads."""
    meta = (f'##INFO=<ID={key},Number={number},Type=Float,Description="x">',)
    record = f"1\t100\t.\tA\t<DEL>\t50\t.\t{key}={value}\tGT\t0/1"
    loaded = load(write(tmp_path, vcf_text([record], version=version, extra_meta=meta)))
    assert loaded.calls.collect()["gt"].to_list() == [[0, 1]]
    assert f"INFO/{key}" in loaded.facets.collect()["dropped_fields"][0]
    assert any(f"INFO/{key} (Type=Float" in w and "not carried" in w for w in loaded.warnings), loaded.warnings


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3", "VCFv4.5"])
@pytest.mark.parametrize("kind", ["Integer", "Float", "Character"])
@pytest.mark.parametrize(("namespace", "key", "record", "column"), [
    ("FORMAT", "FT", "1\t100\t.\tA\tG\t50\t.\t.\tGT:FT\t0/1:PASS", "sample_filters"),
    ("INFO", "SVTYPE", "1\t100\t.\tA\t<DEL>\t50\t.\tSVTYPE=DEL\tGT\t0/1", "sv_type"),
])
def test_a_text_key_declared_as_a_number_is_not_carried(
    tmp_path: Path, namespace: str, key: str, record: str, column: str, kind: str, version: str
) -> None:
    """FT and SVTYPE hold text (`PASS`, `DEL`). Declared Integer, Float or Character, the value does not parse, and
    under VCF 4.2, which reserves no Type, polars-bio fails the whole read on the first one: the file
    loads without the key, which is in `dropped_fields`, with one warning that names it."""
    declared = f'##{namespace}=<ID={key},Number=1,Type={kind},Description="x">'
    loaded = load(write(tmp_path, vcf_text([record], version=version, extra_meta=(declared,))))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    assert calls["gt"].to_list() == [[0, 1]] and calls[column].to_list() == [None]
    assert f"{namespace}/{key}" in facet["dropped_fields"]
    assert sum(f"{namespace}/{key} (Type={kind}" in w and "not carried" in w for w in loaded.warnings) == 1, loaded.warnings


@pytest.mark.parametrize("version", ["VCFv4.2", "VCFv4.3", "VCFv4.4", "VCFv4.5"])
@pytest.mark.parametrize("kind", ["Integer", "Float", "Character"])
def test_gt_declared_with_another_type_is_read_as_the_genotype(tmp_path: Path, kind: str, version: str) -> None:
    """GT is the one key a facet cannot do without. Declared Integer, Float or Character, polars-bio
    refuses a VCF 4.3+ header, and a retry leaving GT out had every row without a genotype (and failed
    the spill). The copy declares GT String and the rows are read as written, in every version, with a
    warning naming the Type; GT is never in `dropped_fields`. The same file declaring String is the
    expected answer."""
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ\t0/1:30", "1\t200\t.\tC\tT,G\t50\t.\t.\tGT:GQ\t1|2:40", "2\t300\t.\tT\tC\t50\t.\t.\tGT:GQ\t./.:."]
    written = vcf_text(records, fmt=("GT", "GQ"), version=version)
    expected = load(write(tmp_path, written, "string.vcf")).calls.collect()
    declared = FORMAT_LINES["GT"].replace("Type=String", f"Type={kind}")
    loaded = load(write(tmp_path, written.replace(FORMAT_LINES["GT"], declared)))
    calls, facet = loaded.calls.collect(), loaded.facets.collect().row(0, named=True)
    for column in ("gt", "gt_alleles", "phased", "gq", "source_gt"):
        assert calls[column].to_list() == expected[column].to_list(), column
    assert calls["gt"].null_count() == 0
    assert "FORMAT/GT" not in (facet["dropped_fields"] or [])
    assert sum(f"FORMAT/GT declared Type={kind}" in w for w in loaded.warnings) == 1, loaded.warnings
    assert not any("refused" in w or "not carried" in w for w in loaded.warnings), loaded.warnings


@pytest.mark.parametrize("fmt", [("DP",), ("GT", "DP")])
def test_records_without_gt_load_with_the_genotype_unknown(tmp_path: Path, fmt: tuple[str, ...]) -> None:
    """A sample column that carries no GT (none declared, or declared and unused) leaves every record's
    genotype unknown: the rows load with `gt` and `gt_alleles` null, never a guessed call. A run of such
    records once failed the spill (polars' `list.gather` over all-null indices). Without a GT declared,
    a warning says no genotype is read."""
    records = ["1\t100\t.\tA\tG\t50\t.\t.\tDP\t12", "1\t200\t.\tC\tT\t50\t.\t.\tDP\t13"]
    calls, _, warnings = frames(write(tmp_path, vcf_text(records, fmt=fmt)))
    assert calls["pos"].to_list() == [100, 200] and calls["dp"].to_list() == [12, 13]
    assert calls["gt"].to_list() == [None, None] and calls["gt_alleles"].to_list() == [None, None]
    assert any(w.startswith("no GT:") for w in warnings) == ("GT" not in fmt), warnings


def test_an_undeclared_hq_the_rows_use_is_read(tmp_path: Path) -> None:
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t50\t.\t.\tGT:HQ\t0|1:51,40"]))
    loaded = load(path)
    assert loaded.calls.collect()["hq"].to_list() == [[51, 40]]
    assert any("FORMAT/HQ" in w and "not declared" in w for w in loaded.warnings)


def test_hq_is_never_filled_from_gq(tmp_path: Path) -> None:
    path = write(tmp_path, vcf_text(["1\t100\t.\tA\tG\t50\t.\t.\tGT:GQ\t0/1:40"], fmt=("GT", "GQ"), extra_meta=(HQ_LINE.format(number=2),)))
    calls = load(path).calls.collect()
    assert calls["gq"].to_list() == [40] and calls["hq"].to_list() == [None]
