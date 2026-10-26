"""F26: every plan `load` spills runs in polars' streaming engine, with no node falling back to the
in-memory engine.

polars 1.44 runs `replace` with a mapping (the canonical `chrom`) as a columnar function the streaming
engine cannot split into batches, an in-memory fallback, and holds every other column of the plan
beside it until that function has the whole column. On PGP-HMS hu416394's gVCF (9,651,143 calls) the
`calls` sink alone held 8.0 GB for an output of 2.6 GB. The fallback is read off the plan's physical
graph, as polars draws it, so the test does not depend on a measurement of memory.
"""

import re
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pytest
from test_load_cgi import LOCI, var_text
from test_load_cgi import write as write_cgi
from test_load_keyvalue import ROWS as KEYVALUE_ROWS
from test_load_keyvalue import write as write_keyvalue
from test_load_nirvana import POSITIONS
from test_load_nirvana import _write as write_nirvana
from test_spill_order import decoy_vcf

from just_vcf import load
from just_vcf import loading as loading_module
from just_vcf.schema import CONTIG_ALIASES, canonical_chrom, canonical_chrom_expr

FALLBACK_COLOUR = '"0.0 0.3 1.0"'
"""The fill colour polars' physical graph gives a node its legend calls "in-memory engine fallback"."""


def fallback_nodes(frame: pl.LazyFrame) -> list[str]:
    """The kinds of the nodes in `frame`'s streaming plan that fall back to the in-memory engine."""
    graph = frame.show_graph(engine="streaming", plan_stage="physical", raw_output=True, show=False)
    assert "in-memory engine fallback" in graph, "polars' graph no longer names its fallback colour: read the legend again"
    labels = re.findall(r'label="([^"]*)",style=filled,fillcolor=' + re.escape(FALLBACK_COLOUR), graph)
    return [label.split("\\n", 1)[0] for label in labels]


SOURCES: dict[str, Callable[[Path], Path]] = {
    "gvcf": decoy_vcf,
    "cgi": lambda tmp_path: write_cgi(tmp_path, var_text(LOCI), "var-GS000000001-ASM.tsv.bz2"),
    "keyvalue": lambda tmp_path: write_keyvalue(tmp_path, KEYVALUE_ROWS),
    "nirvana": lambda tmp_path: write_nirvana(tmp_path / "s1.json", POSITIONS),
}


@pytest.mark.parametrize("name", SOURCES)
def test_every_spilled_plan_streams_without_an_in_memory_fallback(name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spill_files = loading_module._spill_files
    plans: list[tuple[str, pl.LazyFrame]] = []

    def watched(facet_id: str, calls: pl.LazyFrame, blocks: pl.LazyFrame, origin) -> tuple[list[Path], list[Path]]:
        plans.extend([("calls", calls), ("blocks", blocks)])
        return spill_files(facet_id, calls, blocks, origin)

    monkeypatch.setattr(loading_module, "_spill_files", watched)
    load(SOURCES[name](tmp_path))
    assert plans, "the loader spilled nothing"
    assert {frame: nodes for frame, plan in plans if (nodes := fallback_nodes(plan))} == {}


NAMES = [*CONTIG_ALIASES, *CONTIG_ALIASES.values(), "chrUn_KI270302v1", "GL000192.1", "chr1_KI270706v1_random", "", "chrm", "M "]


@pytest.mark.parametrize("engine", ["streaming", "in-memory"])
def test_the_canonical_chrom_is_the_mapping_row_for_row(engine: str) -> None:
    """Every listed spelling, every canonical name, names the list leaves alone, and nulls: the same
    value `canonical_chrom` gives, and the value `replace` over the whole column gave."""
    names = pl.LazyFrame({"chrom": [*NAMES, None, *reversed(NAMES)]}, schema={"chrom": pl.String})
    got = names.select(canonical_chrom_expr("chrom")).collect(engine=engine)["chrom"]
    expected = [None if name is None else canonical_chrom(name) for name in names.collect()["chrom"]]
    assert got.dtype == pl.String
    assert got.to_list() == expected
    assert got.equals(names.select(pl.col("chrom").replace(CONTIG_ALIASES)).collect()["chrom"])


def test_the_canonical_chrom_streams_beside_other_columns() -> None:
    frame = pl.LazyFrame({"chrom": ["chr1", "2"], "pos": [1, 2]}).select(canonical_chrom_expr("chrom"), "pos")
    assert fallback_nodes(frame) == []
    assert fallback_nodes(pl.LazyFrame({"chrom": ["chr1"], "pos": [1]}).select(pl.col("chrom").replace(CONTIG_ALIASES), "pos")) == [
        "columnar-function"
    ], "polars streams `replace` now: canonical_chrom_expr can be the plain expression again"
