import random
import re
from pathlib import Path

import polars as pl
import pytest

from just_vcf.schema import (
    BLOCKS,
    CALLS,
    CONTIG_ALIASES,
    CONTIG_RANK,
    FACETS,
    canonical_chrom,
    canonical_chrom_expr,
    canonical_order,
    conform,
    contig_map,
)

SCHEMA_DOC = Path(__file__).resolve().parent.parent / "docs" / "SCHEMA.md"


def _section(heading: str) -> str:
    text = SCHEMA_DOC.read_text()
    start = text.index(heading)
    end = text.find("\n## ", start + len(heading))
    return text[start:end]


def _table_names(section: str) -> list[str]:
    """First-column names of the markdown table rows in `section`, split on `, `, in order."""
    names: list[str] = []
    for line in section.splitlines():
        match = re.match(r"^\| (`[^|]+) \|", line)
        if match:
            names += re.findall(r"`([a-z_0-9]+)`", match.group(1))
    return names


def test_calls_columns_match_the_doc() -> None:
    assert _table_names(_section("## Decision 6")) == CALLS.names()


def test_facets_columns_match_the_doc() -> None:
    assert _table_names(_section("## Decision 5")) == FACETS.names()


def test_blocks_columns_match_the_doc() -> None:
    listed = re.search(r"`blocks` is small in columns: `([^`]+)`", _section("## Decision 4")).group(1)
    names = [re.match(r"[a-z_]+", part.strip()).group(0) for part in listed.split(",")]
    assert names == BLOCKS.names()


@pytest.mark.parametrize(
    ("source", "canonical"),
    [
        ("chrM", "MT"),
        ("M", "MT"),
        ("MT", "MT"),
        ("chrMT", "MT"),
        ("chr1", "1"),
        ("22", "22"),
        ("chrX", "X"),
        ("23", "X"),
        ("24", "Y"),
        ("25", "XY"),
        ("26", "MT"),
        ("chrUn_KI270302v1", "chrUn_KI270302v1"),
        ("GL000192.1", "GL000192.1"),
        ("chr1_KI270706v1_random", "chr1_KI270706v1_random"),
    ],
)
def test_canonical_chrom(source: str, canonical: str) -> None:
    assert canonical_chrom(source) == canonical


def test_expression_agrees_with_function() -> None:
    names = [*CONTIG_ALIASES, "chrUn_KI270302v1", "HLA-A*01:01:01:01", "chrEBV"]
    frame = pl.DataFrame({"name": names}).with_columns(canonical_chrom_expr("name").alias("chrom"))
    assert frame["chrom"].to_list() == [canonical_chrom(name) for name in names]


def test_contig_map_is_sorted_and_complete() -> None:
    seen = ["chr2", "chrM", "chr2", "chr10", "chrUn_KI270302v1"]
    mapping = contig_map(seen)
    assert [entry["source"] for entry in mapping] == sorted(set(seen))
    assert {entry["source"]: entry["chrom"] for entry in mapping} == {name: canonical_chrom(name) for name in seen}


def _calls(rows: list[tuple]) -> pl.LazyFrame:
    frame = pl.LazyFrame(
        rows, schema={"facet_id": pl.String, "chrom": pl.String, "pos": pl.UInt32, "end_pos": pl.UInt32,
                      "source_ref": pl.String, "alts": pl.List(pl.String)}, orient="row",
    )
    return conform(frame, CALLS)


ROWS = [
    ("f1", "10", 5, 5, "A", ["G"]),
    ("f1", "2", 7, 7, "C", ["T"]),
    ("f1", "2", 7, 7, "C", ["A"]),
    ("f1", "2", 7, None, None, None),
    ("f1", "2", 7, 9, "CTT", ["C"]),
    ("f1", "MT", 1, 1, "G", ["A"]),
    ("f1", "XY", 3, 3, "A", ["C"]),
    ("f1", "X", 3, 3, "A", ["C"]),
    ("f1", "chrUn_KI270302v1", 1, 1, "A", ["T"]),
    ("f1", "GL000192.1", 1, 1, "A", ["T"]),
    ("f0", "22", 1, 1, "A", ["T"]),
]


def test_canonical_order_ranks_contigs_numerically_then_by_name() -> None:
    ordered = canonical_order(_calls(ROWS)).collect()
    f1 = ordered.filter(pl.col("facet_id") == "f1")
    chroms = list(dict.fromkeys(f1["chrom"].to_list()))
    assert chroms == ["2", "10", "X", "MT", "XY", "GL000192.1", "chrUn_KI270302v1"]
    assert ordered["facet_id"][0] == "f0"
    at_2_7 = f1.filter((pl.col("chrom") == "2") & (pl.col("pos") == 7))
    assert at_2_7["end_pos"].to_list() == [7, 7, 9, None]
    assert at_2_7["alts"].to_list()[:2] == [["A"], ["T"]]


def test_canonical_order_is_independent_of_input_order() -> None:
    expected = canonical_order(_calls(ROWS)).collect()
    for seed in range(20):
        shuffled = ROWS[:]
        random.Random(seed).shuffle(shuffled)
        assert canonical_order(_calls(shuffled)).collect().equals(expected)


def test_contig_rank_covers_primary_and_par() -> None:
    assert list(CONTIG_RANK) == [*[str(n) for n in range(1, 23)], "X", "Y", "MT", "XY"]


def test_conform_adds_missing_columns_as_nulls() -> None:
    frame = conform(pl.LazyFrame({"facet_id": ["f"], "chrom": ["1"], "pos": [1]}), CALLS).collect()
    assert frame.schema == CALLS
    assert set(frame.select(pl.exclude("facet_id", "chrom", "pos")).null_count().row(0)) == {1}


def test_conform_refuses_a_position_past_u32() -> None:
    with pytest.raises(pl.exceptions.InvalidOperationError):
        conform(pl.LazyFrame({"facet_id": ["f"], "chrom": ["1"], "pos": [2**32]}), CALLS).collect()


def test_conform_refuses_unknown_columns() -> None:
    with pytest.raises(ValueError, match="not in the schema"):
        conform(pl.LazyFrame({"facet_id": ["f"], "stray": [1]}), CALLS)


def test_per_contig_order_equals_the_global_order(tmp_path: Path) -> None:
    from just_vcf.schema import canonical_order_per_contig

    rng = random.Random(11)
    contigs = ["1", "2", "10", "X", "MT", "XY", "chrUn_KI270302v1", "GL000192.1"]
    rows = [
        ("f", rng.choice(contigs), rng.randrange(1, 50), rng.choice([None, 60, 61]), rng.choice([None, "A", "C"]),
         rng.choice([None, ["G"], ["T"], ["G", "T"]]))
        for _ in range(3000)
    ]
    frame = _calls([(f, c, p, None if e is None else max(e, p), r, a) for f, c, p, e, r, a in rows])
    spill = tmp_path / "calls.parquet"
    frame.sink_parquet(spill, row_group_size=100)
    expected = canonical_order(frame).collect()
    got = canonical_order_per_contig(pl.scan_parquet(spill), sorted({r[1] for r in rows})).collect()
    assert got.equals(expected)
