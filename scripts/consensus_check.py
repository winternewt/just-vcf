"""Check a consensus posterior (RM280) on one person's facets, by aggregate only.

    uv run python scripts/consensus_check.py PLAN.json [--contigs 22] [--out data/output/consensus_check]

The plan is a JSON file kept out of the repo (it names files), for example under `data/interim/`:

    {"reference": "<refget store>", "giab": "<GIAB v3.6 directory>", "atlas": "<method_atlas.parquet>",
     "assembly": "GRCh37",
     "hold_out": [{"a_caller": "...", "a_caller_version": "...", "b_caller": "...", "b_caller_version": "..."}],
     "facets": [{"method": "CGI vcfBeta", "path": "<file>", "caller": null}, ...],
     "leave_out": "<the method of the facet held out for calibration>"}

`hold_out` drops the atlas rows of these method pairs before the rates are taken: the rows the person
under test contributed, so the rates come from other people. A file that loads as several facets (a
PGP freeBayes + VarScan2 merge) names the one to take with `caller`, or none to take them all. Each
facet's REF is read where `resolve_ref` can; the rest are compared by assembly label.

Two runs per contig: the consensus over every facet, and one without `leave_out`. Written to
`<out>.agreement.tsv` (per source method and variant class: the alleles either side carries, how
often the consensus's most probable genotype matches the source's), `<out>.calibration.tsv` (the
held-out facet against the consensus made without it, per confidence bin: the predicted rate of
disagreement against the observed one), `<out>.support.tsv` (rows per supported / contradicted /
silent count) and `<out>.rates.tsv` (the rates each method was weighted by). Aggregates only: no
path, sample name or position is written. Autosomes by default.
"""

import json
import logging
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from just_vcf import Loaded, ReferenceMismatchError, ReferenceStore, concordance, consensus, load, resolve_ref
from just_vcf.atlas import giab_stratifications, wilson_interval
from just_vcf.consensus import method_error_rates

logger = logging.getLogger("consensus_check")
app = typer.Typer(add_completion=False)

BINS = [0.0, 10.0, 20.0, 30.0, 40.0, 60.0, float("inf")]
"""Confidence bins (phred) for the calibration table."""


def _only(loaded: Loaded, caller: str | None) -> Loaded:
    if caller is None:
        return loaded
    facets = loaded.facets.filter(pl.col("caller") == caller)
    ids = facets.select("facet_id").collect()["facet_id"].to_list()
    if len(ids) != 1:
        raise ValueError(f"{len(ids)} facets with caller {caller}")
    return Loaded(facets=facets, calls=loaded.calls.filter(pl.col("facet_id") == ids[0]), blocks=loaded.blocks.filter(pl.col("facet_id") == ids[0]))


def _on(loaded: Loaded, chrom: str) -> Loaded:
    return Loaded(
        facets=loaded.facets, calls=loaded.calls.filter(pl.col("chrom") == chrom), blocks=loaded.blocks.filter(pl.col("chrom") == chrom), links=loaded.links
    )


def _one(loaded: Loaded, facet_id: str) -> Loaded:
    return Loaded(
        facets=loaded.facets.filter(pl.col("facet_id") == facet_id),
        calls=loaded.calls.filter(pl.col("facet_id") == facet_id),
        blocks=loaded.blocks.filter(pl.col("facet_id") == facet_id),
    )


def _carried(sites: pl.DataFrame) -> pl.DataFrame:
    carries = pl.col("a_gt_alleles").list.contains(pl.col("alt")).fill_null(False) | pl.col("b_gt_alleles").list.contains(pl.col("alt")).fill_null(False)
    return sites.filter(carries)


@app.command()
def main(
    plan: Annotated[Path, typer.Argument(help="The plan JSON")],
    out: Annotated[Path, typer.Option(help="Output stem")] = Path("data/output/consensus_check"),
    contigs: Annotated[list[str] | None, typer.Option(help="Contigs (canonical names); autosomes by default")] = None,
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s", force=True)
    for name in ("just_vcf.consensus", "just_vcf.concordance", "just_vcf.reference"):
        logging.getLogger(name).setLevel(logging.ERROR)
    spec = json.loads(plan.read_text())
    store = ReferenceStore.open(Path(spec["reference"]))
    atlas = pl.read_parquet(spec["atlas"])
    for pair in spec.get("hold_out", []):
        same = [pl.col(k).eq_missing(pl.lit(v, dtype=pl.String)) for k, v in pair.items()]
        atlas = atlas.filter(~pl.all_horizontal(same))
    regions = giab_stratifications(spec["giab"], spec["assembly"])
    chroms = contigs or [str(n) for n in range(1, 23)]

    methods: dict[str, Loaded] = {}
    for entry in spec["facets"]:
        whole = _only(load(Path(entry["path"])), entry.get("caller"))
        try:
            whole = resolve_ref(whole, store)
        except ReferenceMismatchError:
            logger.info("%s: REF not read; compared by assembly label", entry["method"])
        methods[entry["method"]] = whole
    held = spec["leave_out"]
    ids = {method: loaded.facets.select("facet_id").collect()["facet_id"].to_list() for method, loaded in methods.items()}

    agreement, calibration, support, rates = [], [], [], None
    for chrom in chroms:
        logger.info("contig %s", chrom)
        parts = {method: _on(loaded, chrom) for method, loaded in methods.items()}
        full = consensus(list(parts.values()), atlas=atlas, regions=regions)
        rates = full.rates if rates is None else rates
        support.append(full.support.collect().group_by("variant_class", "facets", "supported", "contradicted", "silent", "defaulted").len())
        for method, part in parts.items():
            for facet_id in ids[method]:
                sites = _carried(concordance(full, _one(part, facet_id)).sites.collect())
                agreement.append(
                    sites.group_by("variant_class").agg(
                        pl.len().alias("carried"),
                        (pl.col("agreement") == "match").sum().alias("match"),
                        (pl.col("agreement") == "mismatch").sum().alias("mismatch"),
                        (pl.col("agreement") == "ploidy_differs").sum().alias("ploidy_differs"),
                        (pl.col("agreement") == "only_a").sum().alias("consensus_only"),
                        (pl.col("agreement") == "only_b").sum().alias("source_only"),
                    ).with_columns(pl.lit(method).alias("method"))
                )
        others = [part for method, part in parts.items() if method != held]
        loo = consensus(others, atlas=atlas, regions=regions)
        held_rates = method_error_rates(atlas, by="caller")
        held_caller = methods[held].facets.select("caller").collect()["caller"][0]
        e = dict(
            zip(
                *held_rates.filter((pl.col("caller") == held_caller) & (pl.col("region_class") == "all")).select(pl.col("variant_class").cast(pl.String), "genotype_error")
                .to_dict(as_series=False).values()
            )
        )  # fmt: skip
        sites = concordance(loo, _one(parts[held], ids[held][0])).sites.collect()
        both = _carried(sites).filter(pl.col("agreement").cast(pl.String).is_in(["match", "mismatch"]) & (pl.col("variant_class").cast(pl.String) == "snv"))
        sup = loo.support.collect().select("chrom", "pos", "ref", "alt", "copies", "confidence", "ploidy")
        posterior = loo.calls.collect().select("chrom", "pos", pl.col("source_ref").alias("ref"), pl.col("alts").list.first().alias("alt"), pl.col("gp").cast(pl.List(pl.Float64)))
        joined = both.join(sup, on=["chrom", "pos", "ref", "alt"], how="inner").join(posterior, on=["chrom", "pos", "ref", "alt"], how="inner")
        e_snv = e.get("snv")
        p_map = pl.col("gp").list.get(pl.col("copies").cast(pl.Int64))
        predicted = 1 - (p_map * (1 - e_snv) + (1 - p_map) * e_snv / pl.col("ploidy").cast(pl.Float64))
        calibration.append(
            joined.with_columns(pl.col("confidence").cut(BINS[1:-1], left_closed=True).alias("bin"), predicted.alias("_predicted"))
            .group_by("bin")
            .agg(pl.len().alias("alleles"), (pl.col("agreement") == "mismatch").sum().alias("mismatch"), pl.col("_predicted").sum().alias("predicted_sum"))
        )

    ag = (
        pl.concat(agreement).group_by("method", "variant_class").agg(pl.exclude("method", "variant_class").sum())
        .with_columns((pl.col("match") / (pl.col("match") + pl.col("mismatch") + pl.col("ploidy_differs"))).alias("agreement_share"))
        .sort("method", "variant_class")
    )  # fmt: skip
    cal = pl.concat(calibration).with_columns(pl.col("bin").cast(pl.String)).group_by("bin").agg(pl.exclude("bin").sum()).sort("bin")
    cal = cal.with_columns((pl.col("mismatch") / pl.col("alleles")).alias("observed"), (pl.col("predicted_sum") / pl.col("alleles")).alias("predicted"))
    intervals = [wilson_interval(int(m), int(n)) for m, n in zip(cal["mismatch"], cal["alleles"])]
    cal = cal.with_columns(pl.Series("observed_ci_low", [i[0] for i in intervals]), pl.Series("observed_ci_high", [i[1] for i in intervals]))
    sp = pl.concat(support).group_by(pl.exclude("len")).agg(pl.col("len").sum()).sort(pl.exclude("len"))
    out.parent.mkdir(parents=True, exist_ok=True)
    ag.write_csv(out.with_name(out.name + ".agreement.tsv"), separator="\t")
    cal.write_csv(out.with_name(out.name + ".calibration.tsv"), separator="\t")
    sp.with_columns(pl.col("variant_class").cast(pl.String)).write_csv(out.with_name(out.name + ".support.tsv"), separator="\t")
    if rates is not None:
        names = {fid: method for method, fids in ids.items() for fid in fids}
        rates.with_columns(pl.col("facet_id").replace_strict(names).alias("method")).drop("facet_id").with_columns(
            pl.col("variant_class").cast(pl.String), pl.col("rate_basis").cast(pl.String)
        ).write_csv(out.with_name(out.name + ".rates.tsv"), separator="\t")
    with pl.Config(tbl_rows=60, tbl_cols=20, tbl_width_chars=200):
        print(ag)
        print(cal)


if __name__ == "__main__":
    app()
