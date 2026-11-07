"""Derive the marker sample `just_vcf.chip` compares an export with, from Illumina's manifest files.

Each manifest's `[Assay]` section names its probes. A probe named by an rsid (`rs123`, or GSA's
`GSA-rs123` and `exm-rs123`) gives that rsid; every other name (positions, `rs123_MNV`, `rs123.1`)
gives none. The sample keeps the rsids whose number is divisible by `MODULUS`, so an export's ids are
sampled by the same predicate and the two samples compare like the full sets. Each sampled rsid's
chromosome comes from the same rows' `Chr` column, written as `sample_chrom`, a list parallel to
`sample` (RM4: recall is measured over the chromosomes an export carries). The result is the
package file `src/just_vcf/chip_markers.json`; docs/PROVENANCE.md names the manifests and their sha256.

    uv run python scripts/chip_markers.py \\
        --manifest illumina_gsa=GSA-24v3-0-A1-manifest-file-csv.zip \\
        --manifest illumina_omniexpress=infinium-omniexpress-24-v1-2-a1-manifest-file-csv.zip \\
        --out src/just_vcf/chip_markers.json
"""

import csv
import hashlib
import io
import json
import logging
import re
import zipfile
from collections import Counter
from pathlib import Path

import typer

from just_vcf.chip import MARKER_CHROMS, MARKER_MODULUS

logger = logging.getLogger("chip_markers")

_RSID_NAME = re.compile(r"^(?:GSA-|exm-)?rs(\d+)$")

CHROMS: dict[str, str] = {**{name: name for name in MARKER_CHROMS}, "XY": "X"}
"""A manifest `Chr` value → the chromosome written for it. Illumina's `XY` (the pseudo-autosomal region)
is written `X`; `0` (unplaced) and any value not listed are written null."""


def manifest_rsids(path: Path) -> tuple[str, list[int], dict[int, str | None]]:
    """The manifest's descriptor name (`GSA-24v3-0_A1`), the distinct rsid numbers its probes name, and
    each rsid's chromosome (`CHROMS`; null when unplaced).

    An rsid named by several probes takes their chromosome when every placed probe agrees, and null when
    placed probes disagree: a marker of unknown chromosome counts in no recall. The disagreements are
    logged with a count, as are the `Chr` values written null.
    """
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
        if len(members) != 1:
            raise ValueError(f"{path.name}: expected one CSV member, found {len(members)}")
        with archive.open(members[0]) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            descriptor, placed, section, name_index, chr_index = None, {}, None, None, None
            unplaced: Counter[str] = Counter()
            for row in csv.reader(text):
                if not row:
                    continue
                if row[0].startswith("[") and row[0].rstrip().endswith("]"):
                    section, name_index = row[0].strip(), None
                    continue
                if section == "[Heading]" and row[0] == "Descriptor File Name":
                    descriptor = row[1].removesuffix(".bpm")
                elif section == "[Assay]":
                    if name_index is None:
                        name_index, chr_index = row.index("Name"), row.index("Chr")
                    elif found := _RSID_NAME.match(row[name_index]):
                        chroms = placed.setdefault(int(found.group(1)), set())
                        if (chrom := CHROMS.get(row[chr_index])) is None:
                            unplaced[row[chr_index]] += 1
                        else:
                            chroms.add(chrom)
    if descriptor is None:
        raise ValueError(f"{path.name}: no Descriptor File Name in [Heading]")
    disagree = sum(len(chroms) > 1 for chroms in placed.values())
    logger.info("%s: Chr values written null on rsid probes: %s", path.name, dict(sorted(unplaced.items())) or "none")
    logger.info("%s: %d rsids on probes whose chromosomes disagree, written null", path.name, disagree)
    numbers = sorted(placed)
    return descriptor, numbers, {n: next(iter(placed[n])) if len(placed[n]) == 1 else None for n in numbers}


def main(
    manifest: list[str] = typer.Option(..., help="chip=path to an Illumina manifest CSV zip; repeat per chip"),
    out: Path = typer.Option(..., help="JSON file to write"),
) -> None:
    # `force`: importing polars-bio (through just_vcf) installs a root handler, which a plain
    # basicConfig leaves in place at WARNING, and this script's counts would never print.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True)
    chips = {}
    for item in manifest:
        chip, _, path_text = item.partition("=")
        path = Path(path_text)
        descriptor, numbers, chroms = manifest_rsids(path)
        sample = [n for n in numbers if n % MARKER_MODULUS == 0]
        sample_chrom = [chroms[n] for n in sample]
        chips[chip] = {
            "manifest": descriptor,
            "file": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rsids": len(numbers),
            "sample": sample,
            "sample_chrom": sample_chrom,
        }
        classes = Counter("null" if c is None else c if c in ("X", "Y", "MT") else "autosome" for c in sample_chrom)
        logger.info("%s: %s, %d rsids, %d sampled, by chromosome %s", chip, descriptor, len(numbers), len(sample), dict(sorted(classes.items())))
    payload = {"modulus": MARKER_MODULUS, "chips": dict(sorted(chips.items()))}
    out.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    typer.run(main)
