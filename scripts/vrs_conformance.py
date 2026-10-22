"""Dogfood RM271 on a real VCF: load, resolve REF, mint VRS ids, and check a sample against vrs-python.

    uv run python scripts/vrs_conformance.py VCF STORE [--sample N] [--seed S]

Prints aggregates only (counts, timings, peak memory, the warnings' counted lines), never a row, a
sample name or the input's file name, so it can run on any genome. The sample is drawn from the
minted alleles with a fixed seed and spread over the kinds that take different paths through the
normalization (substitution, insertion, deletion, the REF slot), each checked against vrs-python's
id for the same position, REF and allele on the same store.
"""

import argparse
import random
import resource
import time
from pathlib import Path

import polars as pl
from vrs_oracle import StoreProxy, oracle_id

from just_vcf import ReferenceStore, load, mint_vrs_ids, resolve_ref


def _kind(ref: str, allele: str, slot: int) -> str:
    if slot == 0:
        return "ref"
    if len(ref) == len(allele):
        return "substitution"
    return "insertion" if len(allele) > len(ref) else "deletion"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=Path)
    parser.add_argument("store", type=Path)
    parser.add_argument("--sample", type=int, default=1000, help="alleles checked per kind")
    parser.add_argument("--seed", type=int, default=271)
    parser.add_argument("--sample-name", default=None, help="for a multi-sample file")
    args = parser.parse_args()

    store = ReferenceStore.open(args.store)
    clock = time.perf_counter()
    loaded = load(args.vcf, sample=args.sample_name)
    rows = loaded.calls.select(pl.len()).collect().item()
    loaded_at = time.perf_counter()
    resolved = resolve_ref(loaded, store)
    resolved.calls.select(pl.len()).collect()
    resolved_at = time.perf_counter()
    minted = mint_vrs_ids(resolved, store)
    called_at = time.perf_counter()
    calls = minted.calls.select("chrom", "pos", "ref", "alts", "vrs_allele_ids").collect()
    minted_at = time.perf_counter()
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

    print(f"rows {rows}")
    print(
        f"load {loaded_at - clock:.1f} s, resolve_ref {resolved_at - loaded_at:.1f} s, "
        f"mint_vrs_ids {called_at - resolved_at:.1f} s, reading the minted calls back {minted_at - called_at:.1f} s"
    )
    print(f"peak RSS {peak_mb:.0f} MB")
    for warning in minted.warnings[len(resolved.warnings) :]:
        print("  ", warning)

    slots = (
        calls.filter(pl.col("vrs_allele_ids").is_not_null())
        .with_columns(pl.concat_list("ref", "alts").alias("allele"))
        .with_row_index("_row")
        .explode("allele", "vrs_allele_ids", empty_as_null=False)
        .with_columns(pl.int_range(pl.len()).over("_row").alias("slot"))
        .filter(pl.col("vrs_allele_ids").is_not_null())
    )
    (digest,) = minted.facets.select("ref_seqcol_digest").collect().row(0)
    accession = {chrom: f"SQ.{contigs[0].digest}" for chrom, contigs in store.contigs_by_chrom(digest).items() if len(contigs) == 1}
    proxy = StoreProxy(store, [digest])
    rng = random.Random(args.seed)
    by_kind: dict[str, list[tuple]] = {}
    for row in slots.select("chrom", "pos", "ref", "allele", "slot", "vrs_allele_ids").iter_rows():
        by_kind.setdefault(_kind(row[2], row[3], row[4]), []).append(row)
    checked = mismatched = 0
    for kind in sorted(by_kind):
        population = by_kind[kind]
        sample = rng.sample(population, min(args.sample, len(population)))
        bad = sum(oracle_id(proxy, accession[chrom], pos, ref, allele) != got for chrom, pos, ref, allele, _, got in sample)
        print(f"oracle {kind}: {len(sample)} checked of {len(population)} minted, {bad} mismatch")
        checked += len(sample)
        mismatched += bad
    print(f"oracle total: {checked} checked, {mismatched} mismatch")


if __name__ == "__main__":
    main()
