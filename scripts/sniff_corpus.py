"""Sniff every file under a directory and report what the format registry makes of it.

Writes one TSV row per file (path, outcome, format or kind, layers, member, build, reason) and logs an
aggregate by outcome and format. Partial downloads (`*.part`) and manifests are skipped.

    uv run python scripts/sniff_corpus.py $JUST_VCF_SAMPLES_DIR/pgp-hms/participants --out census.tsv
"""

import csv
import logging
from collections import Counter
from pathlib import Path

import typer

from just_vcf.sniff import sniff

logger = logging.getLogger("sniff_corpus")

SKIP_SUFFIXES = (".part", ".jsonl")


def main(
    root: Path = typer.Argument(..., exists=True, file_okay=False, help="directory to walk"),
    out: Path = typer.Option(..., help="TSV to write"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    paths = sorted(p for p in root.rglob("*") if p.is_file() and not p.name.endswith(SKIP_SUFFIXES))
    tally: Counter[tuple[str, str]] = Counter()
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "outcome", "format_or_kind", "layers", "member", "build", "members", "reason"])
        for path in paths:
            result = sniff(path)
            label = result.format or result.kind or ""
            tally[(result.outcome.value, label)] += 1
            writer.writerow([
                path.relative_to(root),
                result.outcome.value,
                label,
                "+".join(result.layers),
                result.member or "",
                result.build or "",
                len(result.members),
                result.reason or "; ".join(result.warnings),
            ])
    for (outcome, label), count in sorted(tally.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1])):
        logger.info("%6d  %-17s %s", count, outcome, label)
    logger.info("%d files, report in %s", len(paths), out)


if __name__ == "__main__":
    typer.run(main)
