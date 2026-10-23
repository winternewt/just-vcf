"""Derive the site panel `just_vcf.relatedness` reads, from somalier's published site VCFs.

somalier (Brent Pedersen, MIT) publishes one sites VCF per build: about 17,766 common SNVs from gnomAD,
chosen for relatedness and sample checks. This keeps the sites that a genotype can be read at the same
way on an array and on a sequenced genome:

- on an autosome (an absence is read as homozygous reference only there, and X and Y are haploid in
  half of the people);
- one base REF, one base ALT, both A, C, G or T;
- not strand-ambiguous (A/T, C/G): an array that writes some probes on the minus strand (F13) gives the
  same two letters, and the pair cannot tell a flipped probe from a genotype;
- REF equal to the base the reference store holds there, for the assembly the VCF is on: an absence is
  read as homozygous for this REF, so a wrong one turns every silence into the wrong genotype.

The result is the package file `src/just_vcf/relatedness_sites.json`; docs/PROVENANCE.md names the
VCFs, their sha256 and the counts each rule dropped.

    uv run python scripts/relatedness_sites.py \\
        --vcf GRCh37=sites.GRCh37.vcf.gz --vcf GRCh38=sites.hg38.vcf.gz \\
        --reference /path/to/refget-store --out src/just_vcf/relatedness_sites.json
"""

import gzip
import hashlib
import json
import logging
from collections import Counter
from pathlib import Path

import polars as pl
import typer

from just_vcf.reference import ReferenceStore
from just_vcf.relatedness import PANEL_VERSION
from just_vcf.schema import CONTIG_RANK, canonical_chrom
from just_vcf.site_genotypes import AUTOSOMES

logger = logging.getLogger("relatedness_sites")

_BASES = frozenset("ACGT")
_PALINDROMES = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


def read_sites(path: Path) -> list[tuple[str, int, str, str]]:
    """`(source chrom, pos, REF, ALT)` for every record of a sites VCF, in file order."""
    records = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt = line.split("\t", 5)[:5]
            records.append((chrom, int(pos), ref, alt))
    return records


def collection_for(store: ReferenceStore, assembly: str) -> str:
    """The one collection in the store whose assembly is `assembly`."""
    found = [digest for digest in store.collections() if store.assembly_label(digest) == assembly]
    if len(found) != 1:
        raise ValueError(f"the store holds {len(found)} collections of {assembly}; need exactly one")
    return found[0]


def derive(records: list[tuple[str, int, str, str]], store: ReferenceStore, assembly: str) -> tuple[pl.DataFrame, Counter[str]]:
    """The kept sites as `chrom, pos, ref, alt` in canonical order, and a count per dropped reason."""
    dropped: Counter[str] = Counter()
    kept: list[tuple[str, int, str, str]] = []
    for source_chrom, pos, ref, alt in records:
        chrom = canonical_chrom(source_chrom)
        if chrom not in AUTOSOMES:
            dropped["not on an autosome"] += 1
        elif len(ref) != 1 or len(alt) != 1 or ref not in _BASES or alt not in _BASES:
            dropped["not a one-base ACGT SNV"] += 1
        elif (ref, alt) in _PALINDROMES:
            dropped["strand-ambiguous (A/T or C/G)"] += 1
        else:
            kept.append((chrom, pos, ref, alt))
    digest = collection_for(store, assembly)
    contigs = store.contigs_by_chrom(digest)
    checked: list[tuple[str, int, str, str]] = []
    by_chrom: dict[str, list[tuple[str, int, str, str]]] = {}
    for site in kept:
        by_chrom.setdefault(site[0], []).append(site)
    for chrom, sites in by_chrom.items():
        serving = contigs.get(chrom, [])
        if len(serving) != 1:
            dropped["contig not in the reference"] += len(sites)
            continue
        bases = store.fetch(serving[0], [(pos - 1, pos) for _, pos, _, _ in sites])
        for site, base in zip(sites, bases):
            if base == site[2]:
                checked.append(site)
            else:
                dropped["REF is not the reference base"] += 1
    frame = pl.DataFrame(checked, schema={"chrom": pl.String, "pos": pl.UInt32, "ref": pl.String, "alt": pl.String}, orient="row")
    frame = frame.unique(["chrom", "pos"], keep="none", maintain_order=True)
    if (duplicates := len(checked) - frame.height):
        dropped["two records at one position"] += duplicates
    order = pl.col("chrom").replace_strict(CONTIG_RANK, return_dtype=pl.UInt16)
    return frame.sort(order, "pos"), dropped


def main(
    vcf: list[str] = typer.Option(..., help="ASSEMBLY=PATH of a somalier sites VCF, once per assembly"),
    reference: Path = typer.Option(..., help="a refget store holding each assembly"),
    out: Path = typer.Option(..., help="the JSON to write"),
) -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    store = ReferenceStore.open(reference)
    assemblies: dict[str, dict] = {}
    for item in vcf:
        assembly, path_text = item.split("=", 1)
        path = Path(path_text)
        records = read_sites(path)
        frame, dropped = derive(records, store, assembly)
        logger.info("%s: %d of %d sites kept; dropped %s", assembly, frame.height, len(records), dict(sorted(dropped.items())))
        assemblies[assembly] = {
            "file": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "records": len(records),
            "dropped": dict(sorted(dropped.items())),
            "chrom": frame["chrom"].to_list(),
            "pos": frame["pos"].to_list(),
            "ref": "".join(frame["ref"].to_list()),
            "alt": "".join(frame["alt"].to_list()),
        }
    payload = {"version": PANEL_VERSION, "source": "somalier sites VCFs (MIT), https://github.com/brentp/somalier/releases", "assemblies": assemblies}
    out.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    typer.run(main)
