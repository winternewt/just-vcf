"""`intervals.covering` against a brute-force scan: every interval over every probe, on frames that tile,
nest, cross, share a start, touch by one base, repeat exactly, or lie apart, across several contigs.
`intervals.contig_join` against the one `join_where` with the contig equality beside the inequalities,
which it replaces."""

import random

import polars as pl

from just_vcf.intervals import contig_join, covering

SHAPES = {
    # (start, end) per contig, in no particular order: each shape the as-of join and the range join split.
    "1": [(1, 10), (11, 20), (21, 30)],  # tiled: no interval reaches the next
    "2": [(1, 100), (20, 30), (25, 40), (50, 50), (101, 120)],  # nested, crossing, a one-base interval, then apart
    "3": [(5, 10), (5, 15), (10, 12), (12, 12), (12, 20)],  # one start shared, an end touching the next start
    "X": [(1, 10), (1, 10), (30, 40), (41, 41)],  # an exact repeat, then adjacent intervals apart
}


def brute_force(probes: pl.DataFrame, intervals: pl.DataFrame) -> list[tuple]:
    return sorted(
        (probe["_probe"], interval["id"])
        for probe in probes.iter_rows(named=True)
        for interval in intervals.iter_rows(named=True)
        if probe["chrom"] == interval["chrom"] and interval["start"] <= probe["at"] <= interval["end"]
    )


def frames(shapes: dict[str, list[tuple[int, int]]], probe_contigs: list[str]) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows = [(chrom, start, end) for chrom, spans in shapes.items() for start, end in spans]
    intervals = pl.DataFrame(rows, schema={"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}, orient="row").with_row_index("id")
    top = max(end for _, _, end in rows) + 2
    probes = pl.DataFrame(
        [(chrom, at) for chrom in probe_contigs for at in range(0, top)], schema={"chrom": pl.String, "at": pl.Int64}, orient="row"
    ).with_row_index("_probe")
    return probes, intervals


def found(probes: pl.DataFrame, intervals: pl.DataFrame, contigs: list[str]) -> list[tuple]:
    hits = covering(probes.lazy(), intervals.lazy(), contigs, at="at", start="start", end="end")
    assert set(hits.columns) == {"_probe", "chrom", "at", "id", "start", "end"}
    return sorted(hits.select("_probe", "id").rows())


def test_every_covering_interval_is_found_for_every_shape() -> None:
    """Each probe position on each contig, against the scan: nothing missed, nothing found twice."""
    probes, intervals = frames(SHAPES, list(SHAPES))
    expected = brute_force(probes, intervals)
    assert found(probes, intervals, list(SHAPES)) == expected
    # Shapes the as-of join alone gets wrong are in the set: a nested interval hides its outer one past
    # its end, a shared start hides the shorter one, and an exact repeat is found twice.
    where = {probe: (chrom, at) for probe, chrom, at in probes.rows()}
    per_position = {where[probe]: sum(1 for p, _ in expected if p == probe) for probe in {p for p, _ in expected}}
    assert (per_position[("2", 35)], per_position[("2", 45)], per_position[("3", 12)], per_position[("X", 5)]) == (2, 1, 4, 2)


def test_only_the_asked_contigs_are_searched() -> None:
    probes, intervals = frames(SHAPES, list(SHAPES))
    asked = ["2", "X"]
    expected = brute_force(probes.filter(pl.col("chrom").is_in(asked)), intervals)
    assert found(probes, intervals, asked) == expected


def test_seeded_random_intervals_match_the_scan() -> None:
    """Random spans, many of them overlapping, on three contigs: the same pairs as the scan."""
    rng = random.Random(20261026)
    shapes = {}
    for chrom in ("1", "7", "MT"):
        spans = []
        for _ in range(60):
            start = rng.randint(1, 400)
            spans.append((start, start + rng.choice([0, 0, 1, 3, 10, 40, 150])))
        shapes[chrom] = spans
    probes, intervals = frames(shapes, list(shapes))
    assert found(probes, intervals, list(shapes)) == brute_force(probes, intervals)


def windows(rng: random.Random, name: str, contigs: list[str | None], rows: int) -> pl.DataFrame:
    """`rows` windows on `contigs` (None a row with no contig), with a row number and a list column, as
    the concordance sides carry."""
    starts = [rng.randint(1, 600) for _ in range(rows)]
    return pl.DataFrame(
        {
            f"{name}_chrom": [rng.choice(contigs) for _ in range(rows)],
            f"_{name}_lo": starts,
            f"_{name}_hi": [start + rng.choice([0, 1, 5, 40, 200]) for start in starts],
            f"_{name}_row": list(range(rows)),
            f"{name}_tags": [[rng.randint(0, 3)] * rng.randint(0, 2) for _ in range(rows)],
        },
        schema={
            f"{name}_chrom": pl.String,
            f"_{name}_lo": pl.Int64,
            f"_{name}_hi": pl.Int64,
            f"_{name}_row": pl.UInt32,
            f"{name}_tags": pl.List(pl.Int64),
        },
    )


def test_contig_join_gives_what_the_equality_join_gave() -> None:
    """Seeded random windows on shared, unshared and missing contigs: the same rows, columns and dtypes
    as `join_where` with the contig equality, for strict and widened inequalities alike (the breakend
    windows' and the SV extents' shapes), and an empty frame of the same schema when no contig is shared."""
    rng = random.Random(20261027)
    cases = [
        (pl.col("_a_lo") < pl.col("_b_hi"), pl.col("_b_lo") < pl.col("_a_hi")),
        (pl.col("_a_lo") - 25 < pl.col("_b_hi"), pl.col("_b_lo") < pl.col("_a_hi") + 25),
    ]
    pairs = 0
    for trial in range(20):
        a = windows(rng, "a", ["1", "2", "X", None], rng.randint(0, 80))
        b = windows(rng, "b", ["1", "2", "MT", None], rng.randint(0, 80))
        for inequalities in cases:
            expected = a.join_where(b, pl.col("a_chrom") == pl.col("b_chrom"), *inequalities)
            got = contig_join(a, b, "a_chrom", "b_chrom", *inequalities)
            assert got.schema == expected.schema, trial
            assert got.sort(got.columns).equals(expected.sort(expected.columns)), trial
            pairs += got.height
    assert pairs > 500  # the trials pair plenty of rows, so the comparison is not of empty frames
    a, b = windows(rng, "a", ["1"], 10), windows(rng, "b", ["2"], 10)
    expected = a.join_where(b, pl.col("a_chrom") == pl.col("b_chrom"), *cases[0])
    got = contig_join(a, b, "a_chrom", "b_chrom", *cases[0])
    assert got.is_empty() and got.schema == expected.schema
