"""Build a method-error atlas (RM279) from a plan of same-person facet pairs.

    uv run python scripts/method_atlas.py PLAN.json [--contigs 22 X] [--out data/output/method_atlas]

The plan is a JSON file kept out of the repo (it names private files), for example under
`data/interim/`:

    {"reference": "<refget store>", "giab": "<GIAB v3.6 directory>",
     "pairs": [{"group": "same reads, two callers", "a": "<facet file>", "b": "<facet file>"}]}

A file that loads into several facets (a merge split per caller) names the one to take with
`"a_caller"` or `"b_caller"`. Where `resolve_ref` refuses either facet (no collection in the store fits
it), the pair is compared as loaded, by its equal assembly label (`concordance` warns);
`<out>.pairs.tsv` says which pairs were (`assembly_check`: `content` or `label`).

Each facet is loaded and its REF read once, however many pairs name it. A pair is compared only when
`same_person` answers True: a pair it calls unrelated, or cannot tell, is skipped and the verdict
written down, so the atlas never pools two people's differences as a method's error. Each compared
pair's counts are written to `<out>.pairs/<n>.parquet` as soon as they are done, and a rerun reuses
them. The output holds aggregates only: `<out>.parquet` and `<out>.tsv` (the pooled atlas) and
`<out>.pairs.tsv` (per pair: group, the two methods, the verdict, alleles compared). No path, sample
name or position is written.
"""

import json
import logging
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from just_vcf import Loaded, ReferenceMismatchError, ReferenceStore, load, resolve_ref, same_person
from just_vcf.atlas import ALL, giab_stratifications, pair_counts, pool_counts

logger = logging.getLogger("method_atlas")
app = typer.Typer(add_completion=False)


def _on(loaded: Loaded, contigs: list[str] | None) -> Loaded:
    if not contigs:
        return loaded
    return Loaded(facets=loaded.facets, calls=loaded.calls.filter(pl.col("chrom").is_in(contigs)), blocks=loaded.blocks.filter(pl.col("chrom").is_in(contigs)))


@app.command()
def main(
    plan: Annotated[Path, typer.Argument(help="The plan JSON (reference, giab, pairs)")],
    out: Annotated[Path, typer.Option(help="Output stem")] = Path("data/output/method_atlas"),
    contigs: Annotated[list[str] | None, typer.Option(help="Compare only these contigs (canonical names)")] = None,
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s", force=True)
    spec = json.loads(plan.read_text())
    store = ReferenceStore.open(Path(spec["reference"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    per_pair = out.with_name(out.name + ".pairs")
    per_pair.mkdir(exist_ok=True)
    loaded: dict[tuple[str, str | None], Loaded] = {}
    regions = {}

    def facet(path: str, caller: str | None) -> Loaded:
        """The file's one facet, or the one by `caller` where the file loads into several."""
        if (path, caller) not in loaded:
            whole = _on(load(Path(path)), contigs)
            loaded[(path, caller)] = whole if caller is None else _by_caller(whole, caller)
        return loaded[(path, caller)]

    def pair(entry: dict) -> tuple[Loaded, Loaded, str]:
        """Both facets with REF read, or, where resolve_ref refuses either (the store holds no
        collection that fits it), both as loaded: concordance then compares them by their equal
        assembly label, with a warning, and the ledger says so."""
        a, b = facet(entry["a"], entry.get("a_caller")), facet(entry["b"], entry.get("b_caller"))
        try:
            return resolve_ref(a, store), resolve_ref(b, store), "content"
        except ReferenceMismatchError as refused:
            logger.warning("pair %s: REF not read (%s); compared by assembly label", entry["group"], str(refused).split(";")[0][:200])
            return a, b, "label"

    ledger = []
    for index, entry in enumerate(spec["pairs"]):
        done = per_pair / f"{index:03d}.parquet"
        row = {"pair": index, "group": entry["group"]}
        if done.exists():
            counts = pl.read_parquet(done)
            sidecar = done.with_suffix(".json")
            check = json.loads(sidecar.read_text())["assembly_check"] if sidecar.exists() else None
            ledger.append(row | {"verdict": "same_person", "assembly_check": check, "alleles": _alleles(counts)} | _methods(counts))
            continue
        a, b, check = pair(entry)
        row["assembly_check"] = check
        verdict = same_person(a, b)
        if verdict is not True:
            logger.warning("pair %d (%s): same_person answered %s; skipped", index, entry["group"], verdict)
            ledger.append(row | {"verdict": "not_same_person" if verdict is False else "unknown", "alleles": None})
            continue
        label = a.facets.collect()["assembly_label"].item()
        if label not in regions:
            regions[label] = giab_stratifications(Path(spec["giab"]), label)
        counts = pair_counts(a, b, regions[label], reference=store, by_contig=True)
        counts.write_parquet(done)
        done.with_suffix(".json").write_text(json.dumps({"assembly_check": check}))
        logger.info("pair %d (%s): %d alleles", index, entry["group"], _alleles(counts))
        ledger.append(row | {"verdict": "same_person", "alleles": _alleles(counts)} | _methods(counts))

    frames = [pl.read_parquet(path) for path in sorted(per_pair.glob("*.parquet"))]
    atlas = pool_counts(frames)
    atlas.write_parquet(out.with_suffix(".parquet"))
    atlas.with_columns(pl.col("variant_class").cast(pl.String)).write_csv(out.with_suffix(".tsv"), separator="\t")
    pl.DataFrame(ledger, infer_schema_length=None).write_csv(out.with_name(out.name + ".pairs.tsv"), separator="\t")
    logger.info("atlas: %d rows from %d pairs", atlas.height, len(frames))


def _by_caller(loaded: Loaded, caller: str) -> Loaded:
    """The facet of `loaded` whose caller is `caller` (a file split per caller loads one facet each)."""
    facets = loaded.facets.filter(pl.col("caller") == caller)
    ids = facets.select("facet_id").collect()["facet_id"].to_list()
    if len(ids) != 1:
        raise ValueError(f"{len(ids)} facets with caller {caller}")
    return Loaded(facets=facets, calls=loaded.calls.filter(pl.col("facet_id") == ids[0]), blocks=loaded.blocks.filter(pl.col("facet_id") == ids[0]))


def _alleles(counts: pl.DataFrame) -> int:
    return int(counts.filter((pl.col("region_class") == ALL) & (pl.col("scope") == "carried"))["alleles"].sum())


def _methods(counts: pl.DataFrame) -> dict[str, str | None]:
    return counts.select("a_caller", "a_caller_version", "b_caller", "b_caller_version").row(0, named=True)


if __name__ == "__main__":
    app()
