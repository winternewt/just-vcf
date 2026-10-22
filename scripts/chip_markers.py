"""Derive the marker sample `just_vcf.chip` compares an export with, from Illumina's manifest files.

Each manifest's `[Assay]` section names its probes. A probe named by an rsid (`rs123`, or GSA's
`GSA-rs123` and `exm-rs123`) gives that rsid; every other name (positions, `rs123_MNV`, `rs123.1`)
gives none. The sample keeps the rsids whose number is divisible by `MODULUS`, so an export's ids are
sampled by the same predicate and the two samples compare like the full sets. The result is the
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
from pathlib import Path

import typer

from just_vcf.chip import MARKER_MODULUS

logger = logging.getLogger("chip_markers")

_RSID_NAME = re.compile(r"^(?:GSA-|exm-)?rs(\d+)$")


def manifest_rsids(path: Path) -> tuple[str, list[int]]:
    """The manifest's descriptor name (`GSA-24v3-0_A1`) and the distinct rsid numbers its probes name."""
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
        if len(members) != 1:
            raise ValueError(f"{path.name}: expected one CSV member, found {len(members)}")
        with archive.open(members[0]) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            descriptor, numbers, section, name_index = None, set(), None, None
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
                        name_index = row.index("Name")
                    elif found := _RSID_NAME.match(row[name_index]):
                        numbers.add(int(found.group(1)))
    if descriptor is None:
        raise ValueError(f"{path.name}: no Descriptor File Name in [Heading]")
    return descriptor, sorted(numbers)


def main(
    manifest: list[str] = typer.Option(..., help="chip=path to an Illumina manifest CSV zip; repeat per chip"),
    out: Path = typer.Option(..., help="JSON file to write"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    chips = {}
    for item in manifest:
        chip, _, path_text = item.partition("=")
        path = Path(path_text)
        descriptor, numbers = manifest_rsids(path)
        sample = [n for n in numbers if n % MARKER_MODULUS == 0]
        chips[chip] = {
            "manifest": descriptor,
            "file": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rsids": len(numbers),
            "sample": sample,
        }
        logger.info("%s: %s, %d rsids, %d sampled", chip, descriptor, len(numbers), len(sample))
    payload = {"modulus": MARKER_MODULUS, "chips": dict(sorted(chips.items()))}
    out.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    typer.run(main)
