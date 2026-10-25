"""What sits in each cell of two facets' confidence bands, by aggregate only (docs/CONSENSUS.md).

    uv run python scripts/band_cells.py A B --reference <refget store> --giab <GIAB v3.6 directory> \
        [--b-caller freeBayes] [--array <array export>] [--out data/output/band_cells]

`A` and `B` are two facets of one person (a file that loads as several facets names the one to take
with `--a-caller` / `--b-caller`). Both have their REF read and VRS ids minted, as
`scripts/method_atlas.py` does, and are compared contig by contig with `concordance`. Kept: the SNV
alleles either side carries, both called, every stating row passing its filters (the atlas's
`carried_pass`). Each is banded by each side's own confidence (`atlas.confidence_band`), as the atlas's
cells are, and counted by what the additive band fit assumes does not matter: the copies of the
allele each side called, whether A's row carries the sample filter `VQLOW` (CGI's), GIAB `easy` or
`difficult`, B's depth, and, with `--array`, the copies of the allele in that facet's genotype at the
same position (an independent third reading where it has a row).

Written to `<out>.tsv`: one row per `(a_band, b_band, a_copies, b_copies, a_vqlow, region, b_dp,
array_copies)` with `alleles` and `mismatch`. No path, sample name or position is written.
"""

import logging
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from just_vcf import Loaded, ReferenceStore, concordance, load, mint_vrs_ids, resolve_ref
from just_vcf.atlas import _confidence_at, _passes, confidence_band, giab_stratifications
from just_vcf.concordance import _covered

logger = logging.getLogger("band_cells")
app = typer.Typer(add_completion=False)

DP_BINS = [10, 20, 30, 40, 60]
"""B's read depth bins (`calls.dp`)."""


def _one(loaded: Loaded, caller: str | None) -> Loaded:
    if caller is None:
        return loaded
    facets = loaded.facets.filter(pl.col("caller") == caller)
    ids = facets.select("facet_id").collect()["facet_id"].to_list()
    if len(ids) != 1:
        raise ValueError(f"{len(ids)} facets with caller {caller}")
    return Loaded(facets=facets, calls=loaded.calls.filter(pl.col("facet_id") == ids[0]), blocks=loaded.blocks.filter(pl.col("facet_id") == ids[0]))


def _on(loaded: Loaded, chrom: str) -> Loaded:
    return Loaded(facets=loaded.facets, calls=loaded.calls.filter(pl.col("chrom") == chrom), blocks=loaded.blocks.filter(pl.col("chrom") == chrom))


def _contigs(loaded: Loaded) -> set[str]:
    return set(loaded.calls.select(pl.col("chrom").unique()).collect()["chrom"].to_list())


@app.command()
def main(
    a: Annotated[Path, typer.Argument(help="Facet A (banded rows: CGI var, say)")],
    b: Annotated[Path, typer.Argument(help="Facet B")],
    reference: Annotated[Path, typer.Option(help="The refget store")],
    giab: Annotated[Path, typer.Option(help="The GIAB v3.6 stratifications directory")],
    assembly: Annotated[str, typer.Option(help="The assembly of both facets")] = "GRCh37",
    a_caller: Annotated[str | None, typer.Option(help="Facet A's caller, where its file loads as several")] = None,
    b_caller: Annotated[str | None, typer.Option(help="Facet B's caller, where its file loads as several")] = None,
    array: Annotated[Path | None, typer.Option(help="A third facet read at the same positions (an array export)")] = None,
    out: Annotated[Path, typer.Option(help="Output stem")] = Path("data/output/band_cells"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s", force=True)
    for name in ("just_vcf.concordance", "just_vcf.reference", "just_vcf.vrs"):
        logging.getLogger(name).setLevel(logging.ERROR)
    store = ReferenceStore.open(reference)
    regions = giab_stratifications(giab, assembly, ["easy", "difficult"])
    fa = mint_vrs_ids(resolve_ref(_one(load(a), a_caller), store), store)
    fb = mint_vrs_ids(resolve_ref(_one(load(b), b_caller), store), store)
    third = None
    if array is not None:
        third = load(array).calls.select("chrom", "pos", pl.col("gt_alleles").alias("_x")).collect().unique(["chrom", "pos"], keep="none")
    parts = []
    for chrom in sorted(_contigs(fa) & _contigs(fb)):
        pa, pb = _on(fa, chrom), _on(fb, chrom)
        sites = concordance(pa, pb, reference=store).sites.collect()
        carries = (
            pl.col("a_gt_alleles").list.contains(pl.col("alt")).fill_null(False)
            | pl.col("b_gt_alleles").list.contains(pl.col("alt")).fill_null(False)
            | (pl.col("a_copies") > 0).fill_null(False)
            | (pl.col("b_copies") > 0).fill_null(False)
        )
        sites = sites.filter(
            carries & _passes("a") & _passes("b") & (pl.col("variant_class").cast(pl.String) == "snv")
            & pl.col("agreement").cast(pl.String).is_in(["match", "mismatch"])
        )  # fmt: skip
        for side, piece in (("a", pa), ("b", pb)):
            at = _confidence_at(piece).select("chrom", pl.col("pos").cast(pl.UInt32).alias(f"{side}_pos"), pl.col("confidence").alias(f"_{side}_conf"))
            sites = sites.join(at, on=["chrom", f"{side}_pos"], how="left").with_columns(confidence_band(pl.col(f"_{side}_conf")).alias(f"{side}_band"))
        rows_a = pa.calls.select("chrom", pl.col("pos").alias("a_pos"), pl.col("sample_filters").alias("_a_sf")).collect().unique(["chrom", "a_pos"], keep="first", maintain_order=True)
        rows_b = pb.calls.select("chrom", pl.col("pos").alias("b_pos"), pl.col("dp").alias("_b_dp")).collect().unique(["chrom", "b_pos"], keep="first", maintain_order=True)
        sites = sites.join(rows_a, on=["chrom", "a_pos"], how="left").join(rows_b, on=["chrom", "b_pos"], how="left")
        if third is not None:
            sites = sites.join(third, on=["chrom", "pos"], how="left")
            called = ~pl.col("_x").list.eval(pl.element().is_null()).list.any()
            sites = sites.with_columns(
                pl.when(pl.col("_x").is_null()).then(pl.lit("absent")).when(~called).then(pl.lit("no_call"))
                .otherwise(_alt_copies().cast(pl.String)).alias("array_copies")
            )  # fmt: skip
        else:
            sites = sites.with_columns(pl.lit(None, dtype=pl.String).alias("array_copies"))
        probe = sites.with_columns(pl.col("pos").cast(pl.Int64).alias("_at")).sort("chrom", "_at")
        probe = probe.with_columns(
            pl.when(_covered(probe, regions.intervals["easy"])).then(pl.lit("easy"))
            .when(_covered(probe, regions.intervals["difficult"])).then(pl.lit("difficult"))
            .otherwise(pl.lit("neither")).alias("region"),
            pl.col("_a_sf").list.contains("VQLOW").fill_null(False).alias("a_vqlow"),
            pl.col("_b_dp").cut(DP_BINS).cast(pl.String).alias("b_dp"),
        )  # fmt: skip
        keys = ["a_band", "b_band", "a_copies", "b_copies", "a_vqlow", "region", "b_dp", "array_copies"]
        parts.append(probe.group_by(keys).agg(pl.len().alias("alleles"), (pl.col("agreement").cast(pl.String) == "mismatch").sum().alias("mismatch")))
        logger.info("contig %s: %d alleles", chrom, probe.height)
    keys = ["a_band", "b_band", "a_copies", "b_copies", "a_vqlow", "region", "b_dp", "array_copies"]
    table = pl.concat(parts).group_by(keys).agg(pl.col("alleles", "mismatch").sum()).sort(keys, nulls_last=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.write_csv(out.with_suffix(".tsv"), separator="\t")
    logger.info("%d cells, %d alleles", table.height, table["alleles"].sum())


def _alt_copies() -> pl.Expr:
    """The copies of `alt` among the third facet's called alleles (`_x`), SNVs being one base."""
    first = (pl.col("_x").list.get(0, null_on_oob=True) == pl.col("alt")).cast(pl.Int32).fill_null(0)
    second = (pl.col("_x").list.get(1, null_on_oob=True) == pl.col("alt")).cast(pl.Int32).fill_null(0)
    return first + second


if __name__ == "__main__":
    app()
