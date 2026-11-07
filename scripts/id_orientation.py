"""Learn which way 23andMe v3 and v4 write each indel probe's `I`/`D` codes, from whole genomes (RM281).

The package file `src/just_vcf/id_orientation.json` lists, per chip and rsid, whether the chip's `I` is
Illumina's `I` (`same`) or Illumina's `D` (`flipped`); `resolve_ref` reads a v3 or v4 code row only when
its probe is listed. This script writes it from pairs of a person's 23andMe export and the same
person's Complete Genomics `var` genome (any whole-genome callset the library loads with reference
blocks will do).

Per pair: the export's code rows are resolved by rsid against the manifests' indel probes with the
codes read as written, whatever the chip's gate says (`just_vcf.reference._ungated_code_rows`, the
resolution `resolve_ref` runs behind the gate). Of the rows that resolve homozygous, a row is reference
evidence where a `ref` block of the genome covers its whole span and no call of the genome carrying an
ALT lies within `ORIENTATION_VARIANT_WINDOW` bases of it: the person is hom-ref there, so a row read
hom-ref says the chip's codes mean what Illumina's do, a row read hom-alt that they mean the other way
round. Each person counts once per chip, probe and kind, however many exports they have. The table is
then the rule's (`just_vcf.probes.MIN_ORIENTATION_PEOPLE`): at least that many people one way and none
the other. A hom-alt row over a variant of the genome is counted in the log and never as evidence.

`--pairs` is a TSV with a header and the columns `participant`, `chip` (`23andme_v3` or `23andme_v4`),
`export` and `genome` (paths). The chip is the column's: an export too small to name its chip is read as
the column says; one whose facet names another chip is skipped. Logs counts only, never a row, a
position, an rsid or a path, so it is safe on a private genome.

    uv run python scripts/id_orientation.py --pairs pairs.tsv --manifests $JUST_VCF_MANIFESTS_DIR \\
        --reference $JUST_VCF_REFERENCE_DIR --out src/just_vcf/id_orientation.json
"""

import csv
import logging
import os
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
import typer

from just_vcf import load
from just_vcf.intervals import contig_join, covering
from just_vcf.loading import Loaded
from just_vcf.probes import (
    ORIENTATION_CHIPS,
    ORIENTATION_VARIANT_WINDOW,
    orientation_json,
    orientation_payload,
    orientation_tables,
)
from just_vcf.reference import _ungated_code_rows
from just_vcf.schema import BlockKind

logger = logging.getLogger("id_orientation")

PAIR_COLUMNS = ("participant", "chip", "export", "genome")

DEFAULT_SOURCE = "PGP-HMS participants' public 23andMe exports and Complete Genomics var genomes (CC0)"

HOM_REF, HOM_ALT = "hom_ref", "hom_alt"


def read_pairs(path: Path) -> list[dict[str, str]]:
    """The TSV's rows, each with the four `PAIR_COLUMNS`. Raises `ValueError` when a column is missing."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(PAIR_COLUMNS) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path.name}: missing column(s) {missing}; the columns are {list(PAIR_COLUMNS)}")
        return [{column: row[column].strip() for column in PAIR_COLUMNS} for row in reader]


def homozygous_sites(export: Loaded, reference: Path | None) -> tuple[pl.DataFrame, Counter[str]]:
    """The export's code rows that resolve homozygous with the codes as written (`chrom`, `pos`,
    `end_pos`, `rsid`, `kind` hom-ref or hom-alt), and the resolution's counts. A half-call, a haploid
    row and a het are not homozygous."""
    rows, counts = _ungated_code_rows(export, reference)
    gt = pl.col("gt")
    homozygous = rows.filter((gt.list.len() == 2) & (gt.list.drop_nulls().list.len() == 2) & (gt.list.n_unique() == 1))
    kind = pl.when(gt.list.first() == 0).then(pl.lit(HOM_REF)).otherwise(pl.lit(HOM_ALT))
    sites = homozygous.select(
        "chrom", pl.col("pos").cast(pl.Int64), pl.col("end_pos").cast(pl.Int64), "rsid", kind.alias("kind")
    ).with_row_index("_site")
    counts["homozygous"] += sites.height
    return sites, counts


def reference_evidence(sites: pl.DataFrame, genome: Loaded) -> pl.DataFrame:
    """`sites` with `covered` (a `ref` block of the genome covers the whole span) and `near_variant` (a
    call of the genome carrying an ALT within `ORIENTATION_VARIANT_WINDOW` bases of it). Only the
    sites' contigs are read from the genome."""
    contigs = sorted(set(sites["chrom"].drop_nulls().to_list()))
    blocks = genome.blocks.filter((pl.col("kind") == BlockKind.REF.value) & pl.col("chrom").is_in(contigs)).select(
        "chrom", pl.col("pos").cast(pl.Int64).alias("_block_pos"), pl.col("end_pos").cast(pl.Int64).alias("_block_end")
    )
    over = covering(sites.select("_site", "chrom", "pos", "end_pos"), blocks, contigs, at="pos", start="_block_pos", end="_block_end")
    covered = set(over.filter(pl.col("_block_end") >= pl.col("end_pos"))["_site"].to_list())
    variants = (
        genome.calls.filter(pl.col("chrom").is_in(contigs) & (pl.col("alts").list.len() > 0))
        .select(
            pl.col("chrom").alias("_variant_chrom"),
            pl.col("pos").cast(pl.Int64).alias("_variant_pos"),
            pl.coalesce("end_pos", "pos").cast(pl.Int64).alias("_variant_end"),
        )
        .collect()
    )
    window = ORIENTATION_VARIANT_WINDOW
    near = contig_join(
        sites.select("_site", "chrom", "pos", "end_pos"), variants, "chrom", "_variant_chrom",
        pl.col("_variant_pos") <= pl.col("end_pos") + window, pl.col("_variant_end") >= pl.col("pos") - window,
    )  # fmt: skip
    near_sites = set(near["_site"].to_list())
    return sites.with_columns(
        pl.col("_site").is_in(sorted(covered)).alias("covered"), pl.col("_site").is_in(sorted(near_sites)).alias("near_variant")
    )


def tally_pair(
    participant: str, chip: str, export: Loaded, genome: Loaded, reference: Path | None, evidence: dict[str, dict[tuple[str, str], set[str]]]
) -> Counter[str]:
    """One pair's reference evidence added to `evidence` (kind → (chip, rsid) → participants); returns
    the pair's counts (what resolved, what was evidence, what was not and why)."""
    sites, counts = homozygous_sites(export, reference)
    read = reference_evidence(sites, genome)
    usable = pl.col("covered") & ~pl.col("near_variant")
    for kind, rsid in read.filter(usable).select("kind", "rsid").iter_rows():
        evidence[kind].setdefault((chip, rsid), set()).add(participant)
    for kind in (HOM_REF, HOM_ALT):
        of_kind = read.filter(pl.col("kind") == kind)
        counts[f"{kind} over a reference block"] += of_kind.filter(usable).height
        counts[f"{kind} near a variant"] += of_kind.filter(pl.col("near_variant")).height
        counts[f"{kind} not covered by a reference block"] += of_kind.filter(~pl.col("covered") & ~pl.col("near_variant")).height
    return counts


def table_evidence(evidence: dict[str, dict[tuple[str, str], set[str]]]) -> dict[str, dict[str, tuple[int, int]]]:
    """Per chip and rsid, the people hom-ref over a reference block and the people hom-alt over one."""
    keys = sorted(set(evidence[HOM_REF]) | set(evidence[HOM_ALT]))
    table: dict[str, dict[str, tuple[int, int]]] = {}
    for chip, rsid in keys:
        table.setdefault(chip, {})[rsid] = (len(evidence[HOM_REF].get((chip, rsid), ())), len(evidence[HOM_ALT].get((chip, rsid), ())))
    return table


def learn(pairs: list[dict[str, str]], reference: Path | None, source: str) -> tuple[dict, Counter[str]]:
    """The orientation payload learned from `pairs`, one pair loaded at a time, and the run's counts."""
    evidence: dict[str, dict[tuple[str, str], set[str]]] = {HOM_REF: {}, HOM_ALT: {}}
    read_per_chip: Counter[str] = Counter()
    people: set[str] = set()
    tally: Counter[str] = Counter()
    for number, pair in enumerate(pairs, 1):
        chip = pair["chip"]
        if chip not in ORIENTATION_CHIPS:
            tally["pairs skipped: chip not one the table covers"] += 1
            continue
        export, genome = load(pair["export"]), load(pair["genome"])
        export_facets, genome_facets = export.facets.collect(), genome.facets.collect()
        if export_facets.height != 1 or genome_facets.height != 1:
            tally["pairs skipped: a file loads as several facets"] += 1
            continue
        export_facet, genome_facet = export_facets.row(0, named=True), genome_facets.row(0, named=True)
        if export_facet["chip"] is not None and export_facet["chip"] != chip:
            tally["pairs skipped: the export's facet names another chip"] += 1
            continue
        if export_facet["assembly_label"] is None or export_facet["assembly_label"] != genome_facet["assembly_label"]:
            tally["pairs skipped: export and genome not on one stated assembly"] += 1
            continue
        counts = tally_pair(pair["participant"], chip, export, genome, reference, evidence)
        logger.info("pair %d of %d (%s): %s", number, len(pairs), chip, ", ".join(f"{key} {value}" for key, value in sorted(counts.items())))
        tally.update({f"{chip}: {key}": value for key, value in counts.items()})
        read_per_chip[chip] += 1
        people.add(pair["participant"])
    learned_from = f"{source}: {sum(read_per_chip.values())} pair(s) of an export and the same person's whole genome, {len(people)} people"
    payload = orientation_payload(table_evidence(evidence), read_per_chip, learned_from, datetime.now(UTC).date().isoformat())
    return payload, tally


def main(
    pairs: Path = typer.Option(..., exists=True, dir_okay=False, help="TSV: participant, chip, export, genome"),
    manifests: Path = typer.Option(..., exists=True, file_okay=False, help="directory of Illumina's manifest zips"),
    reference: Path = typer.Option(None, exists=True, file_okay=False, help="refget store directory (default: the library's)"),
    out: Path = typer.Option(..., dir_okay=False, help="JSON file to write"),
    source: str = typer.Option(DEFAULT_SOURCE, help="what the pairs are, written into the table's provenance"),
) -> None:
    # `force`: importing polars-bio (through just_vcf) installs a root handler, which a plain
    # basicConfig leaves in place at WARNING, and this script's counts would never print.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    logging.getLogger("just_vcf").setLevel(logging.ERROR)
    # The library reads the manifests from this variable and nowhere else.
    os.environ["JUST_VCF_MANIFESTS_DIR"] = str(manifests)
    payload, tally = learn(read_pairs(pairs), reference, source)
    orientation_tables(payload)  # the table the library will read, checked before it is written
    for key, value in sorted(tally.items()):
        logger.info("run: %s: %d", key, value)
    for chip, entry in sorted(payload["chips"].items()):
        kinds = Counter(probe["orientation"] for probe in entry["probes"].values())
        logger.info("table: %s from %d pair(s): %d same, %d flipped", chip, entry["pairs"], kinds["same"], kinds["flipped"])
    out.write_text(orientation_json(payload))


if __name__ == "__main__":
    typer.run(main)
