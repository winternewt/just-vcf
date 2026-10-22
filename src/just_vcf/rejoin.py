"""Re-join the rows a pre-split VCF wrote one per ALT into one site each (RM270, docs/SCHEMA.md).

The result is a derived facet (charter 5): the parent is never changed. Only rows the loader flagged
in `split_signal` take part, grouped by position. A group is joined from its rows alone, as
`bcftools norm -m+` does: no reference sequence and no reads. What split rows cannot give back (the
likelihood vectors over the whole allele list) is withheld, since rebuilding it would mean calling
genotypes again.
"""

import hashlib
import json
from collections import Counter
from itertools import groupby

import polars as pl

from just_vcf import cache
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.schema import CALLS, SCHEMA_VERSION, EvidenceBasis, RefCheck, canonical_order_per_contig

REJOIN_VERSION = 1
"""Part of the derived `facet_id`: bumped whenever the join rule changes what it emits."""

RELATION = "split_rejoined"
"""The derived facet's `relation`."""

_KEPT_WHEN_EQUAL = ("phased", "phase_sets", "phase_set_kinds", "dp", "qual")
"""Columns a joined row keeps only when every row of its group agrees; otherwise null (not known)."""


def rejoin_split(loaded: Loaded) -> Loaded:
    """A derived facet in which each position's split-flagged rows are one site again.

    Per position, the rows with a `split_signal` are joined when they can be read as one site:
    REF is the longest `source_ref` (the others must be prefixes of it, and each ALT is extended by the
    REF bases it lacks); `alts` are the distinct extended alleles, sorted; each GT slot takes the one
    allele the rows place there, `0` when none places an ALT and a row says `0`, and null when every row
    is null there. A group that cannot be read so (two different ALTs in one slot, rows of different
    ploidy, a symbolic allele, REFs that disagree) keeps its rows as written, counted in the warnings.
    """
    cache.sweep_stale(cache.REJOIN_SPILL, _VERSIONS)
    parent = loaded.facets.collect()
    if parent.height != 1:
        raise ValueError(f"rejoin_split takes one facet; this Loaded has {parent.height}")
    facet = parent.row(0, named=True)
    facet_id = _derived_id(facet["facet_id"])

    flagged = loaded.calls.filter(pl.col("split_signal").is_not_null()).collect()
    joined: list[dict] = []
    kept: list[dict] = []
    unjoined: Counter[str] = Counter()
    # Plain rows grouped in Python: a DataFrame per group costs far more than the join itself.
    for _, group in groupby(flagged.sort("chrom", "pos").iter_rows(named=True), key=lambda row: (row["chrom"], row["pos"])):
        rows = list(group)
        result = _join(rows)
        if isinstance(result, str):
            unjoined[result] += 1
            kept += rows
        else:
            joined.append(result)

    warnings = []
    if joined:
        warnings.append(f"{len(joined)} split site(s) re-joined from {flagged.height - len(kept)} row(s)")
    warnings += [f"{count} split site(s) left as written: {reason}" for reason, count in sorted(unjoined.items())]

    if flagged.height:
        # Re-sorting the parent's lazy frame per contig would replay every contig's sort for each one.
        # The rows are spilled once instead, as `load` does, and read back one contig at a time.
        rest = loaded.calls.filter(pl.col("split_signal").is_null())
        rebuilt = pl.DataFrame(kept + joined, schema=CALLS).lazy()
        spill = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
        with cache.writing(spill, cache.REJOIN_SPILL, _VERSIONS, cache.derived_origin(facet, facet_id)):
            spill.mkdir(parents=True, exist_ok=True)
            pl.concat([rest, rebuilt], how="vertical").with_columns(pl.lit(facet_id).alias("facet_id")).sink_parquet(
                spill / "calls.parquet", row_group_size=SPILL_ROW_GROUP
            )
        contigs = [entry["chrom"] for entry in facet["contig_map"] or []]
        calls = canonical_order_per_contig(pl.scan_parquet(spill / "calls.parquet"), contigs)
    else:
        calls = loaded.calls.with_columns(pl.lit(facet_id).alias("facet_id"))
    facets = parent.lazy().with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.lit(facet["facet_id"]).alias("derived_from"),
        pl.lit(RELATION).alias("relation"),
    )
    return Loaded(
        facets=facets,
        calls=calls,
        blocks=loaded.blocks.with_columns(pl.lit(facet_id).alias("facet_id")),
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )


_VERSIONS = {"SCHEMA_VERSION": SCHEMA_VERSION, "REJOIN_VERSION": REJOIN_VERSION}
"""The version keys of a rejoin spill, in its cache manifest."""


def _derived_id(parent_id: str) -> str:
    key = json.dumps({"derived_from": parent_id, "operation": RELATION, "version": REJOIN_VERSION}, sort_keys=True)
    return hashlib.sha256(key.encode()).hexdigest()[:32]


def _join(rows: list[dict]) -> dict | str:
    """One joined row, or the reason the group stays as written."""
    if any(row["gt"] is None for row in rows):
        return "a row has no genotype"
    if any(row["source_ref"] is None for row in rows):
        return "a row has no source REF (normalize moved it)"
    if any(allele.startswith("<") or allele == "*" for row in rows for allele in row["alts"] or []):
        return "a symbolic allele"
    ref = max((row["source_ref"] for row in rows), key=len)
    if any(not ref.startswith(row["source_ref"]) for row in rows):
        return "REFs that are not prefixes of one another"
    if len({len(row["gt"]) for row in rows}) > 1:
        return "rows of different ploidy"
    extended = [[allele + ref[len(row["source_ref"]):] for allele in row["alts"] or []] for row in rows]
    alts = sorted({allele for row_alts in extended for allele in row_alts})

    gt: list[int | None] = []
    for slot in range(len(rows[0]["gt"])):
        placed = set()
        for row, row_alts in zip(rows, extended):
            index = row["gt"][slot]
            if index is not None and index > 0:
                if index > len(row_alts):
                    return "a genotype index past the row's ALTs"
                placed.add(row_alts[index - 1])
        if len(placed) > 1:
            return "two ALTs in one genotype slot"
        if placed:
            gt.append(alts.index(placed.pop()) + 1)
        elif any(row["gt"][slot] == 0 for row in rows):
            gt.append(0)
        else:
            gt.append(None)

    alleles = [ref, *alts]
    first = rows[0]  # every row of a group shares chrom and pos
    joined = {name: None for name in CALLS.names()}
    joined.update(
        chrom=first["chrom"],
        pos=first["pos"],
        end_pos=first["pos"] + len(ref) - 1,
        ids=sorted({value for row in rows for value in row["ids"] or []}) or None,
        source_ref=ref,
        ref_check=RefCheck.NOT_RUN.value,
        alts=alts,
        gt=gt,
        gt_alleles=[None if index is None else alleles[index] for index in gt],
        evidence_basis=EvidenceBasis.ASSAY.value,
        filters=_filters(row["filters"] for row in rows),
        sample_filters=_filters(row["sample_filters"] for row in rows),
    )
    for name in _KEPT_WHEN_EQUAL:
        values = [row[name] for row in rows]
        joined[name] = values[0] if all(value == values[0] for value in values) else None
    return joined


def _filters(lists) -> list[str] | None:
    """The union of the rows' filters, sorted; `PASS` only when no row failed one. Null when no row said."""
    found = {value for values in lists for value in values or []}
    if len(found) > 1:
        found.discard("PASS")
    return sorted(found) or None
