"""A CGI file is read one locus range at a time (`cgi.locus_ranges`, `cgi.range_frames`), so a whole
`var` no longer holds every haplotype in memory at once (17.5 GB for a whole PGP-HMS `var`, 2026-10-26).

The expected frames and warnings here are the read `load` made before: the whole file in one lazy plan,
kept in this file as the oracle and run at test time over the same file. Ranges must change the plan
only: the same rows, order, dtypes and warnings, whatever the range size.
"""

import bz2
import random
import threading
import time
from pathlib import Path

import polars as pl
import pytest
from test_load_cgi import (
    COUNTED,
    LINKED,
    LOCI,
    NOT_WHOLE,
    SCORED,
    VAF_DIFFERS,
    Locus,
    Seg,
    counted_master_var,
    master_var_text,
    var_text,
    write,
)

from just_vcf import load, sniff
from just_vcf.loaders import cgi
from just_vcf.schema import canonical_order
from just_vcf.sniff import open_payload


def whole_file(path: Path, body: Path) -> tuple[pl.DataFrame, pl.DataFrame, list[str]]:
    """The calls and blocks (without `facet_id`, in canonical order) and the loader's own warnings, as
    `_load_cgi_copy` made them before ranges: the whole file in one plan."""
    sniffed = sniff(path)
    with open_payload(sniffed) as handle:
        header = cgi.copy_body(handle, body)
    rows = cgi.prepare(cgi.scan(body, header), sniffed.format)
    if sniffed.format == cgi.VAR:
        uniform, segments = cgi.var_parts(rows)
        haplotypes = cgi.var_haplotypes(segments)
    else:
        uniform, haplotypes = cgi.master_var_parts(rows)
    checks = cgi.locus_checks(haplotypes)
    good = haplotypes.join(checks.filter(pl.col("problem").is_null()).select("locus"), on="locus", how="semi")
    calls, none_called = cgi.call_rows(good)
    stretches = pl.concat([uniform.select("chromosome", "begin", "end", "varType"), none_called], how="vertical")
    blocks = cgi.block_rows(stretches)
    counts = pl.collect_all(
        [
            rows.select(cgi.malformed(rows).sum().alias("malformed")),
            checks.filter(pl.col("problem").is_not_null()).group_by("problem").len().sort("problem"),
            uniform.group_by("varType").len().sort("varType"),
            calls.filter(pl.col("ad_withheld").is_not_null()).group_by("ad_withheld").len().sort("ad_withheld"),
            checks.filter(pl.col("problem").is_null()).select(pl.col("link_conflicts").sum()),
            good.select(
                (pl.col("score").is_not_null() & ~pl.col("score").is_between(0, cgi.HQ_MAX)).sum().alias("hq_out"),
                (pl.col("vaf_differs") & pl.col("score").is_between(0, cgi.HQ_MAX)).fill_null(False).sum().alias("vaf_differs"),
                pl.col("score_not_whole").sum().alias("not_whole"),
                pl.col("score_not_number").sum().alias("not_number"),
            ),
        ]
    )
    warnings = []
    if counts[0]["malformed"][0]:
        warnings.append(f"{counts[0]['malformed'][0]} row(s) whose locus, begin or end is not a number, or whose end lies before its begin, dropped")
    warnings += [f"{n} locus(es) with {problem} dropped" for problem, n in counts[1].iter_rows()]
    warnings += [f"{n} row(s) without `ad`: {reason}" for reason, n in counts[3].iter_rows()]
    if counts[4].item():
        warnings.append(f"{counts[4].item()} haplotype(s) whose rows name more than one hapLink left without a phase set")
    scores = counts[5].row(0, named=True)
    if scores["hq_out"]:
        warnings.append(f"{scores['hq_out']} haplotype score(s) outside 0..{cgi.HQ_MAX} left out of `hq`")
    if scores["not_whole"]:
        warnings.append(f"{scores['not_whole']} haplotype score value(s) not a whole number left out of `hq`")
    if scores["not_number"]:
        warnings.append(f"{scores['not_number']} haplotype score value(s) not a number left out of `hq`")
    if scores["vaf_differs"]:
        warnings.append(f"{scores['vaf_differs']} haplotype(s) whose varScoreVAF differs from the varScoreEAF read into `hq`")
    by_type = dict(counts[2].iter_rows())
    if by_type.get("PAR-called-in-X"):
        warnings.append(f"{by_type['PAR-called-in-X']} chrY stretch(es) called on chrX (PAR-called-in-X) kept as no-call blocks")
    calls_frame, blocks_frame = cgi.to_frames(calls, blocks, "oracle")
    expected_calls, expected_blocks = pl.collect_all([canonical_order(calls_frame), canonical_order(blocks_frame)])
    body.unlink()
    return expected_calls.drop("facet_id"), expected_blocks.drop("facet_id"), warnings


def shuffled(text: str, seed: int = 5) -> str:
    """The same file with its data rows out of locus order."""
    head, body = text.split(">", 1)
    column_line, *rows = body.splitlines()
    random.Random(seed).shuffle(rows)
    return head + ">" + column_line + "\n" + "\n".join(rows) + "\n"


def one_locus_on_two_contigs(text: str) -> str:
    """The chrX locus renumbered to share the first chr1 variant locus's number: no CGI file does
    this, but a range is cut by locus number, so it must read as the whole file did."""
    first = next(line.split("\t")[0] for line in text.splitlines() if "\tchr1\t" in line and "\tall\t" not in line)
    return "\n".join(
        "\t".join([first, *line.split("\t")[1:]]) if "\tchrX\t" in line and "\tall\t" not in line else line for line in text.splitlines()
    ) + "\n"


BAD = [
    *LOCI,
    Locus("chr1", 900000, 900003, haps=((Seg(900000, 900001, "ref", "A", "A"), Seg(900002, 900003, "snp", "C", "T")), (Seg(900000, 900003, "ref", "AGC", "AGC"),))),
    Locus("chr1", 920029, 920032, haps=((Seg(920029, 920031, "ref", "GA", "GA", link="7"), Seg(920031, 920032, "snp", "T", "C", link="9")), (Seg(920029, 920032, "ref", "GAT", "GAT", link="8"),))),
]
"""LOCI with a locus whose sub-rows leave a gap and a haplotype naming two hapLinks."""

FILES = {
    "var": lambda: (var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"),
    "var shuffled": lambda: (shuffled(var_text(LOCI)), "var-GS000000001-ASM.tsv.bz2"),
    "var format 2.0": lambda: (var_text(LOCI, version="2.0"), "var-GS000000001-ASM.tsv.bz2"),
    "var bad rows and loci": lambda: (var_text(BAD).rstrip("\n") + "\n9x\t2\tall\tchr1\tbegin\t5\tref\t=\t=\n", "bad.tsv.bz2"),
    "var one locus on two contigs": lambda: (one_locus_on_two_contigs(var_text(LOCI)), "var-GS000000001-ASM.tsv.bz2"),
    "var scores": lambda: (var_text([*SCORED, *VAF_DIFFERS[:1], *NOT_WHOLE[2:4]]), "var-GS000000009-ASM.tsv.bz2"),
    "var links": lambda: (var_text(LINKED), "var-GS000000004-ASM.tsv.bz2"),
    "masterVarBeta": lambda: (master_var_text(LOCI), "masterVarBeta-GS000000001-ASM.tsv.bz2"),
    "masterVarBeta shuffled": lambda: (shuffled(master_var_text(LOCI)), "masterVarBeta-GS000000001-ASM.tsv.bz2"),
    "masterVarBeta read counts": lambda: (counted_master_var(COUNTED), "masterVarBeta-GS000000003-ASM.tsv.bz2"),
    "masterVarBeta scores": lambda: (master_var_text([*SCORED, *NOT_WHOLE]), "masterVarBeta-GS000000009-ASM.tsv.bz2"),
}


@pytest.mark.parametrize("size", [1, 3, 7, None], ids=["1", "3", "7", "default"])
@pytest.mark.parametrize("name", list(FILES))
def test_ranges_give_the_rows_and_warnings_the_whole_file_gave(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str, size: int | None) -> None:
    text, file_name = FILES[name]()
    path = write(tmp_path, text, file_name)
    expected_calls, expected_blocks, expected_warnings = whole_file(path, tmp_path / "oracle-body.tsv")
    if size is not None:
        monkeypatch.setattr(cgi, "CHUNK_ROWS", size)
    loaded = load(path)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    assert calls.drop("facet_id").equals(expected_calls), "same rows, order and values"
    assert blocks.drop("facet_id").equals(expected_blocks)
    assert calls.schema == loaded.calls.collect_schema() and calls.drop("facet_id").schema == expected_calls.schema, "same dtypes"
    assert blocks.drop("facet_id").schema == expected_blocks.schema
    # The loader's own warnings come first, as before; the variant-class ones follow.
    assert list(loaded.warnings[: len(expected_warnings)]) == expected_warnings
    assert not [w for w in loaded.warnings[len(expected_warnings):] if w in expected_warnings]


def test_the_ranges_hold_whole_loci_and_every_row_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With ranges of 3 rows, each range is read on its own: fewer than 3 rows before its last locus,
    no locus in two ranges, every well-formed row in one. Fails where the file is read in one plan,
    since there is no range to read."""
    text = shuffled(var_text(LOCI))
    seen: list[pl.DataFrame] = []
    original = cgi.range_frames

    def spy(rows: pl.DataFrame, format_id: str):
        seen.append(rows)
        return original(rows, format_id)

    monkeypatch.setattr(cgi, "range_frames", spy)
    monkeypatch.setattr(cgi, "CHUNK_ROWS", 3)
    load(write(tmp_path, text, "var-GS000000001-ASM.tsv.bz2"))
    rows = text.split(">", 1)[1].splitlines()[1:]
    loci_rows: dict[int, int] = {}
    for line in rows:
        loci_rows[int(line.split("\t")[0])] = loci_rows.get(int(line.split("\t")[0]), 0) + 1
    assert len(seen) > 1
    per_range = [set(frame["locus"].to_list()) for frame in seen]
    for frame, loci in zip(seen, per_range):
        assert frame.height - loci_rows[max(loci)] < 3, (frame.height, loci)
        assert frame.height == sum(loci_rows[locus] for locus in loci), "a range holds its loci whole"
    assert sum(len(loci) for loci in per_range) == len(set().union(*per_range)) == len(loci_rows), "each locus in one range"


def test_a_range_is_cut_by_locus_number_not_file_order() -> None:
    rows = pl.LazyFrame({"locus": [5, 1, 1, 9, 2, None, 2, 2, 7]})
    ranges = cgi.locus_ranges(rows, 3)
    # Rows per locus 1:2, 2:3, 5:1, 7:1, 9:1 start 0, 2, 5, 6 and 7 rows in: ranges start at 0, 3 and 6.
    assert ranges == [(1, 5), (5, 7), (7, None)]
    held = [rows.filter(cgi.in_range(low, high)).collect()["locus"].to_list() for low, high in ranges]
    assert sorted(x for part in held for x in part) == [1, 1, 2, 2, 2, 5, 7, 9], "every numbered row once; a null locus is malformed"
    assert cgi.locus_ranges(pl.LazyFrame({"locus": pl.Series([], dtype=pl.Int64)}), 3) == [(None, None)]


# --- Real files (opt-in) -----------------------------------------------------------------------------


def _smallest_var(pgp_participants: Path) -> Path:
    found = sorted(pgp_participants.glob("*/arvados-*/var-*-ASM.tsv.bz2"), key=lambda p: p.stat().st_size)
    if not found:
        pytest.skip("no PGP-HMS CGI var on disk")
    return found[0]


def test_a_real_var_head_reads_the_same_in_ranges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pgp_participants: Path) -> None:
    """Counts only. The first 300,000 data rows of the smallest PGP-HMS `var` (consented public), read
    in ranges of 20,000 rows, give the frames and warnings the whole-file plan gives."""
    head = []
    with bz2.open(_smallest_var(pgp_participants), "rt") as text:
        for line in text:
            head.append(line)
            if line.startswith(">"):
                break
        head += [line for _, line in zip(range(300_000), text)]
    path = write(tmp_path, "".join(head), "var-head-ASM.tsv.bz2")
    expected_calls, expected_blocks, expected_warnings = whole_file(path, tmp_path / "oracle-body.tsv")
    monkeypatch.setattr(cgi, "CHUNK_ROWS", 20_000)
    loaded = load(path)
    calls, blocks = pl.collect_all([loaded.calls, loaded.blocks])
    assert expected_calls.height > 10_000
    assert calls.drop("facet_id").equals(expected_calls) and blocks.drop("facet_id").equals(expected_blocks)
    assert list(loaded.warnings[: len(expected_warnings)]) == expected_warnings


def _rss() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def test_a_real_whole_var_loads_within_a_bounded_peak(pgp_participants: Path) -> None:
    """Counts only. The smallest PGP-HMS `var` loaded whole stays within 6 GB above the process's
    memory before the load. Measured 2026-10-26 on a 233 MB `var` (18.9M rows): 17.5 GB peak read in
    one plan, 3.3 GB in ranges (process total)."""
    peak = [0]
    done = threading.Event()

    def sample() -> None:
        while not done.is_set():
            peak[0] = max(peak[0], _rss())
            time.sleep(0.05)

    before = _rss()
    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    loaded = load(_smallest_var(pgp_participants))
    done.set()
    sampler.join()
    assert loaded.calls.select(pl.len()).collect().item() > 1_000_000
    assert peak[0] - before < 6 * 2**30, f"peak {(peak[0] - before) / 2**30:.1f} GB above the start"
