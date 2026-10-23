"""Learn the probe sets of chips whose vendor publishes no manifest, from exports under a directory (RM281).

Walks the directory, sniffs every file, loads those in an array format whose exports write a row per
probe (23andMe, AncestryDNA, Living DNA), and hands them to `just_vcf.probes.learn_probe_sets`, which
writes one set per chip and assembly into the cache (`probes/learned/`), read from then on by
`probe_set` and `chip_absence`. Logs counts only: participant uploads can carry real names in file names.

    uv run python scripts/learn_probe_sets.py $JUST_VCF_SAMPLES_DIR/pgp-hms/participants
"""

import logging
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import typer

from just_vcf import load
from just_vcf.loading import Loaded
from just_vcf.probes import learn_probe_sets
from just_vcf.sniff import sniff

logger = logging.getLogger("learn_probe_sets")

FORMATS = frozenset({"23andme", "23andme_csv", "ancestrydna", "living_dna"})
"""Array formats that write a row for every probe of the chip, no-calls included."""

SKIP_SUFFIXES = (".part", ".jsonl")


def exports(root: Path, tally: Counter[str]) -> Iterator[Loaded]:
    for path in sorted(p for p in root.rglob("*") if p.is_file() and not p.name.endswith(SKIP_SUFFIXES)):
        sniffed = sniff(path)
        if sniffed.format not in FORMATS:
            continue
        try:
            loaded = load(sniffed)
        except Exception as error:  # noqa: BLE001 - one unreadable upload must not end the walk; counted by type
            tally[f"not loaded: {type(error).__name__}"] += 1
            continue
        facets = loaded.facets.collect()
        if facets.height != 1:
            # An appended export (F3): which facet's rows are the chip's export is not known.
            tally["skipped: several facets in one file"] += 1
            continue
        facet = facets.row(0, named=True)
        tally[f"{facet['chip'] or 'chip unknown'} on {facet['assembly_label'] or 'assembly unknown'}"] += 1
        yield loaded


def main(
    root: Path = typer.Argument(..., exists=True, file_okay=False, help="directory to walk"),
    min_exports: int = typer.Option(3, help="exports one chip needs on one assembly before its set is written"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("just_vcf").setLevel(logging.ERROR)
    tally: Counter[str] = Counter()
    written = learn_probe_sets(exports(root, tally), min_exports=min_exports)
    for label, count in sorted(tally.items()):
        logger.info("exports: %s: %d", label, count)
    for found in written:
        logger.info("learned: %s on %s, %d sites, %s", found.chip, found.assembly_label, found.sites.height, found.source)


if __name__ == "__main__":
    typer.run(main)
