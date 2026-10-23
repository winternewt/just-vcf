"""Lift a facet to another assembly through a chain file, as a derived facet (RM6, charter invariant 5).

`liftover(loaded, chain, reference)` returns a `Lifted`: the derived facet (`relation = liftover_to`,
the Cat-VRS term, `derived_from` the parent, which is never changed) and a report with one row per
parent row saying where it went or why it did not. The rules are in docs/SCHEMA.md ("Lifting"); in
short:

1. **The chain is checked against both assemblies first.** A chain header states each contig's size
   on both sides (`tSize`, `qSize`). The source sizes must equal the collection the parent's REF was
   read from, the target sizes the target collection, contig by contig through `canonical_chrom`. A
   primary contig (1-22, X, Y) that differs refuses the lift: the chain is for another assembly. Any
   other contig that differs leaves its rows unlifted, counted (UCSC hg19's chrM, 16,571 bases,
   against GRCh37's rCRS MT).
2. **Positions.** A row's REF span must lie inside one ungapped block, and inside no block with
   another target (chains that overlap on the source side with different targets make a position's
   target ambiguous; blocks repeating one mapping count once, F14).
   Otherwise it is unlifted: no block at all, a span that crosses a block edge, or more than one
   target.
3. **Alleles.** On a minus-strand block every allele is reverse-complemented. An indel whose alleles
   all start with one base (VCF's padding base) then has that base at their end instead, so it is
   re-anchored: the base is dropped from every allele and the target's base before the lifted span
   is put in front, which moves the row one base left. The indel is not left-aligned on the target
   here; `normalize` on the lifted facet does that.
   REF is then read from the target sequence over the lifted span (charter 10: never carried over).
   Where it equals the parent's REF the row keeps its alleles; where it equals one of the ALTs, REF
   and that ALT swap and the genotype indices follow (`bcftools +liftover`'s and Picard's swap);
   where it equals none, the row is unlifted (`ref_mismatch_after_lift`).
4. **Withheld in this slice**, each counted: a minus-strand indel whose lifted span starts at the
   target contig's first base (no base before it to anchor on), symbolic alleles and `*`, rows
   flagged as one record of a pre-split site, and the `blocks` frame.

A lifted row's `source_ref`, `ref_check`, `source_gt` and VRS ids are null: they described the
source's representation on the source assembly. `gl`, `gp` and `ad` stay indexed over
`source_alleles`, never reordered (SCHEMA.md Decision 2); a swapped row's `dosage` is withheld, since it
counted the other allele. Re-mint ids on the target with `mint_vrs_ids`.

Coordinates: `pos` is 1-based; chains and the store are 0-based half-open. `reference._span` converts
a row to the store's range and `_to_pos` converts a lifted start back, the only two places (charter 13).
"""

import gzip
import hashlib
import json
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf import cache
from just_vcf.links import LIFTOVER_TO, derived_links
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.reference import PRIMARY, Contig, ReferenceMismatchError, ReferenceStore, _derive_array_alleles, _span, _store_for
from just_vcf.schema import CALLS, CONTIG_RANK, FACETS, SCHEMA_VERSION, AbsenceMeans, AssemblyLabelBasis, RefCheck, canonical_chrom, canonical_order, empty

logger = logging.getLogger(__name__)

LIFT_VERSION = 2
"""Part of the derived `facet_id`: bumped whenever the rule changes what `liftover` emits."""

RELATION = LIFTOVER_TO
"""The derived facet's `relation`: the Cat-VRS term, verbatim (charter, GA4GH from the start)."""

UNLIFTED_REASONS: frozenset[str] = frozenset(
    {
        # What the parent row is.
        "ref_unread",
        "source_ref_mismatch",
        "allele_not_sequence",
        "split_flagged",
        # What the chain says about its position.
        "contig_not_in_chain",
        "chain_contig_mismatch",
        "no_chain_block",
        "split_across_blocks",
        "multiple_targets",
        "target_contig_unmatched",
        # What the target says about its alleles.
        "indel_minus_strand",
        "ref_mismatch_after_lift",
    }
)
"""Why a parent row is not in the lifted facet; the report's `lift_reason`. An open vocabulary (str)."""

REASON_TEXT: dict[str, str] = {
    "ref_unread": "REF was not read from the source assembly (resolve_ref left it unread)",
    "source_ref_mismatch": "the source REF disagrees with the source assembly (ref_check mismatch)",
    "allele_not_sequence": "an allele that is not a base sequence (symbolic, `*`, `<*>`, a breakend); SV extents are not lifted",
    "split_flagged": "flagged as one record of a pre-split site (RM270); run rejoin_split first",
    "contig_not_in_chain": "a contig the chain does not cover",
    "chain_contig_mismatch": "a contig whose size in the chain differs from the source assembly's",
    "no_chain_block": "no chain block covers it (a gap in the alignment)",
    "split_across_blocks": "its REF span crosses the edge of a chain block",
    "multiple_targets": "it lies in more than one chain block (overlapping chains), so its target is ambiguous",
    "target_contig_unmatched": "the chain's target contig is not in the target collection, or its size differs",
    "indel_minus_strand": "an indel on a minus-strand block lifted to the target contig's first base, with no base before it to anchor on",
    "ref_mismatch_after_lift": "the target's REF equals none of the row's alleles",
}

PARENT_COLUMNS_KEPT = ("ids", "phased", "phase_sets", "phase_set_kinds", "gl", "gp", "dosage", "evidence_basis", "gq", "dp", "ad", "qual", "filters", "sample_filters", "split_signal")
"""`calls` columns a lifted row carries from its parent as they were (phase sets and dosage aside, below)."""

REPORT = pl.Schema(
    {
        **CALLS,
        "lift_chrom": pl.String,
        "lift_pos": pl.UInt32,
        "lift_strand": pl.String,
        "lift_swapped": pl.Boolean,
        "lift_reason": pl.String,
    }
)
"""`Lifted.report`: the parent's `calls` row, then where it went (`lift_chrom`, 1-based `lift_pos`,
the block's strand `+`/`-`, whether REF and an ALT swapped) or why not (`lift_reason`, null when lifted)."""

_SEQUENCE = r"^[ACGTN]+$"
_COMPLEMENT = (["A", "C", "G", "T"], ["T", "G", "C", "A"])
_FAR = 1 << 40
"""An end past every contig, closing the last gap segment of a chain's source contig."""


class Lifted(BaseModel):
    """What `liftover` returns: the derived facet, and the fate of every parent row."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    loaded: Loaded
    """The lifted facet (`relation = liftover_to`)."""
    report: pl.LazyFrame
    """One row per parent `calls` row, in the parent's order, with `REPORT`'s columns: CrossMap's
    `.unmap` (the rows with a `lift_reason`) and bcftools' `--write-src` (the rows without) in one frame."""
    reasons: dict[str, int]
    """Unlifted rows per reason, sorted by reason."""


@dataclass(frozen=True)
class Chain:
    """A chain file read into blocks, with each side's contig sizes."""

    name: str
    """The file's name (not its path: a path is a label, and can carry a person's name)."""
    sha256: str
    blocks: pl.DataFrame
    """One row per ungapped block: `t_chrom` (canonical), `t_start`, `t_end` (0-based half-open on the
    source), `q_chrom` (the target's name as the chain writes it), `q_start` (0-based on the chain's
    target strand), `q_size`, `q_minus`, `chain` (the chain's index in the file)."""
    source_sizes: dict[str, set[int]]
    """Source contig sizes by canonical name; more than one size means the chain contradicts itself."""
    target_sizes: dict[str, set[int]]
    """Target contig sizes by the chain's own names."""


def read_chain(path: Path | str) -> Chain:
    """Read a UCSC chain file (plain or gzipped, by magic bytes) into blocks.

    The format: a header `chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd
    id`, then lines `size dt dq` and a last `size`. `t` is the source, `q` the target; `tStrand` is
    always `+`, and on a `-` chain the q coordinates count on the reverse strand.
    """
    path = Path(path)
    data = path.read_bytes()
    text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
    columns: dict[str, list] = {name: [] for name in ("t_chrom", "t_start", "t_end", "q_chrom", "q_start", "q_size", "q_minus", "chain")}
    source_sizes: dict[str, set[int]] = {}
    target_sizes: dict[str, set[int]] = {}
    chain = -1
    t = q = 0
    t_chrom = q_chrom = ""
    q_size = 0
    q_minus = False
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "chain":
            if len(parts) < 12 or parts[4] != "+":
                raise ValueError(f"{path.name}: not a chain header (or a source strand other than +): {line[:80]}")
            chain += 1
            t_chrom, q_chrom = canonical_chrom(parts[2]), parts[7]
            source_sizes.setdefault(t_chrom, set()).add(int(parts[3]))
            target_sizes.setdefault(q_chrom, set()).add(int(parts[8]))
            q_size, q_minus = int(parts[8]), parts[9] == "-"
            t, q = int(parts[5]), int(parts[10])
            continue
        size = int(parts[0])
        for name, value in zip(columns, (t_chrom, t, t + size, q_chrom, q, q_size, q_minus, chain)):
            columns[name].append(value)
        if len(parts) == 3:
            t += size + int(parts[1])
            q += size + int(parts[2])
    blocks = pl.DataFrame(
        columns,
        schema={"t_chrom": pl.String, "t_start": pl.Int64, "t_end": pl.Int64, "q_chrom": pl.String, "q_start": pl.Int64, "q_size": pl.Int64, "q_minus": pl.Boolean, "chain": pl.Int64},
    ).with_row_index("block")
    return Chain(name=path.name, sha256=hashlib.sha256(data).hexdigest(), blocks=blocks, source_sizes=source_sizes, target_sizes=target_sizes)


def _segments(blocks: pl.DataFrame) -> pl.DataFrame:
    """One source contig's axis cut where the set of targets changes: `seg_start`, how many distinct
    targets cover the segment (`depth`), a block giving it when exactly one does, and running counts of
    covered and multiply-covered segments, so a span over several segments is judged by two lookups.

    A target is a mapping, not a block: blocks with the same target contig, strand and offset
    (`q_start - t_start`) send every shared base to the same place, so they count once. Ensembl's
    GRCh38_to_GRCh37 chain lists every block twice, in two chains (F14); counting blocks would make
    the whole genome ambiguous.
    """
    keys = list(zip(blocks["q_chrom"], blocks["q_minus"], (blocks["q_start"] - blocks["t_start"])))
    ids = blocks["block"].to_list()
    events = sorted([(s, 1, i) for i, s in enumerate(blocks["t_start"])] + [(e, -1, i) for i, e in enumerate(blocks["t_end"])])
    starts, depths, singles = [0], [0], [None]
    active: dict[int, tuple] = {}
    index = 0
    while index < len(events):
        at = events[index][0]
        while index < len(events) and events[index][0] == at:
            _, step, i = events[index]
            if step == 1:
                active[i] = keys[i]
            else:
                active.pop(i, None)
            index += 1
        targets = {key: i for i, key in sorted(active.items())}
        depth = len(targets)
        single = next(iter(targets.values())) if depth == 1 else None
        if at == starts[-1]:
            depths[-1], singles[-1] = depth, single
            continue
        # A cut where the one target stays the same (overlapping blocks of one mapping) is no edge.
        previous = singles[-1]
        if depth == depths[-1] and (depth == 0 or (depth == 1 and keys[previous] == keys[single])):
            continue
        starts.append(at)
        depths.append(depth)
        singles.append(single)
    blocks_of = [None if i is None else ids[i] for i in singles]
    frame = pl.DataFrame({"seg_start": starts, "depth": depths, "seg_block": blocks_of}, schema={"seg_start": pl.Int64, "depth": pl.Int64, "seg_block": pl.UInt32})
    return frame.with_row_index("seg").with_columns(
        (pl.col("depth") >= 1).cast(pl.Int64).cum_sum().alias("covered_cum"),
        (pl.col("depth") >= 2).cast(pl.Int64).cum_sum().alias("multi_cum"),
    )


def _to_pos(start: pl.Expr) -> pl.Expr:
    """A lifted 0-based start as a 1-based `pos` (charter 13: the way back from `reference._span`)."""
    return (start + 1).cast(pl.UInt32)


def _complement(allele: pl.Expr, minus: pl.Expr) -> pl.Expr:
    """An allele on the target strand: uppercase, and reverse-complemented on a minus-strand block."""
    upper = allele.str.to_uppercase()
    return pl.when(minus).then(upper.str.reverse().str.replace_many(*_COMPLEMENT)).otherwise(upper)


def _complement_list(alleles: pl.Expr, minus: pl.Expr) -> pl.Expr:
    upper = alleles.list.eval(pl.element().str.to_uppercase())
    flipped = alleles.list.eval(pl.element().str.to_uppercase().str.reverse().str.replace_many(*_COMPLEMENT))
    return pl.when(minus).then(flipped).otherwise(upper)


def liftover(
    loaded: Loaded,
    chain: Path | str | Chain,
    reference: "ReferenceStore | Path | str | None" = None,
    *,
    target: str | None = None,
) -> Lifted:
    """The facet lifted to another assembly through `chain`, as a derived facet, and a report of every row.

    `loaded` holds one facet whose REF `resolve_ref` read; `reference` is a store (or its directory;
    default the store at `reference_dir()`) holding that collection and the target. `target` names the
    target collection's digest; without it, the one collection whose contig sizes match the chain's
    target side on every primary contig is used, and anything else is refused. The input is not changed.
    """
    cache.sweep_stale(cache.LIFT_SPILL, _VERSIONS)
    parent = loaded.facets.collect()
    if parent.height != 1:
        raise ValueError(f"liftover takes one facet; this Loaded has {parent.height}")
    facet = parent.row(0, named=True)
    source_digest = facet["ref_seqcol_digest"]
    if source_digest is None:
        raise ValueError(f"facet {facet['facet_id']}: REF has not been read; run resolve_ref first, then liftover")
    store, _ = _store_for(reference)
    if source_digest not in store.collections():
        raise ReferenceMismatchError(f"facet {facet['facet_id']}: the store does not hold collection {source_digest}, which REF was read from")
    chain = chain if isinstance(chain, Chain) else read_chain(chain)
    source_ok = _check_source(chain, store, source_digest, facet["facet_id"])
    target_digest = _choose_target(chain, store, target, facet["facet_id"])
    target_contigs = _target_contigs(chain, store, target_digest)

    facet_id = _derived_id(facet["facet_id"], chain.sha256, target_digest)
    blocks = loaded.blocks.select(pl.len()).collect().item()
    entry = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
    with cache.spilling(entry, cache.LIFT_SPILL, _VERSIONS, cache.derived_origin(facet, facet_id)) as spill:
        # Streamed to parquet once, so each contig below is a filtered scan; its bytes key the result,
        # with what else the warnings name, so a run on the same input reuses what an earlier one wrote.
        staged = spill.stage / "input.parquet"
        loaded.calls.sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
        named = {"chain_name": chain.name, "blocks": blocks, "labels": [store.assembly_label(source_digest), store.assembly_label(target_digest)]}
        key = cache.content_key([staged], {**facet, "_lift": named})
        done = spill.find(key)
        if done is None:
            done = _lift_into(spill, key, staged, chain, source_ok, store, target_contigs, facet_id, blocks, source_digest, target_digest)
    calls_parts = done.files("calls")
    report_parts = done.files("report")
    calls = pl.scan_parquet(calls_parts) if calls_parts else empty(CALLS)
    report = pl.scan_parquet(report_parts) if report_parts else empty(REPORT)
    # Read back from the report, so a reused generation gives the same counts as the run that wrote it.
    found = dict(report.filter(pl.col("lift_reason").is_not_null()).group_by("lift_reason").len().collect().iter_rows())
    reasons = {reason: found[reason] for reason in sorted(found)}
    warnings = done.warnings
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)

    present = calls.select(pl.col("chrom").unique()).collect()["chrom"].to_list()
    names = {chrom: contig.name for chrom, contig in target_contigs.items()}
    label = store.assembly_label(target_digest)
    chain_input = [{"role": "chain", "name": chain.name, "sha256": chain.sha256}]
    facets = parent.lazy().with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.lit(facet["facet_id"]).alias("derived_from"),
        pl.lit(RELATION).alias("relation"),
        pl.lit(label, dtype=pl.String).alias("assembly_label"),
        pl.lit(AssemblyLabelBasis.INFERRED.value if label else None, dtype=FACETS["assembly_label_basis"]).alias("assembly_label_basis"),
        pl.lit(None, dtype=pl.String).alias("seqcol_digest"),
        pl.lit(None, dtype=pl.String).alias("sorted_name_length_pairs_digest"),
        pl.lit(target_digest).alias("ref_seqcol_digest"),
        pl.lit(AbsenceMeans.UNKNOWN.value, dtype=FACETS["absence_means"]).alias("absence_means"),
        pl.lit(None, dtype=FACETS["header_contigs"]).alias("header_contigs"),
        pl.lit(
            sorted(({"source": names.get(chrom, chrom), "chrom": chrom} for chrom in present), key=lambda entry: entry["source"]),
            dtype=FACETS["contig_map"],
        ).alias("contig_map"),
        pl.lit(chain_input, dtype=FACETS["derivation_inputs"]).alias("derivation_inputs"),
    )
    derived = Loaded(
        facets=facets.select(FACETS.names()),
        calls=calls,
        blocks=loaded.blocks.clear().with_columns(pl.lit(facet_id).alias("facet_id")),
        links=derived_links(
            loaded.links, facet["facet_id"], facet_id, RELATION, step="liftover", step_version=str(LIFT_VERSION), derivation_inputs=chain_input
        ),
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )
    return Lifted(loaded=derived, report=report, reasons=reasons)


_VERSIONS = {"SCHEMA_VERSION": SCHEMA_VERSION, "LIFT_VERSION": LIFT_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
"""The version keys of a lift spill, in its cache manifest."""


def _lift_into(
    spill: cache.Spill,
    key: str,
    staged: Path,
    chain: Chain,
    source_ok: dict[str, bool],
    store: ReferenceStore,
    target_contigs: dict[str, Contig],
    facet_id: str,
    blocks: int,
    source_digest: str,
    target_digest: str,
) -> cache.Committed:
    """Lift the staged input contig by contig into the spill's stage, and commit it as the generation `key`."""
    source = pl.scan_parquet(staged)
    chroms = source.select(pl.col("chrom").unique()).collect()["chrom"].to_list()
    counts: Counter[str] = Counter()
    lifted_parts: list[Path] = []
    report_parts: list[Path] = []
    for index, chrom in enumerate(sorted(chroms, key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))):
        rows = source.filter(pl.col("chrom") == chrom).collect()
        report, lifted = _lift_contig(rows, chrom, chain, source_ok.get(chrom), store, target_contigs, facet_id, counts)
        report_parts.append(spill.stage / f"report-{index:05d}.parquet")
        report.write_parquet(report_parts[-1], row_group_size=SPILL_ROW_GROUP)
        if lifted.height:
            lifted_parts.append(spill.stage / f"lifted-{index:05d}.parquet")
            lifted.write_parquet(lifted_parts[-1], row_group_size=SPILL_ROW_GROUP)
    staged.unlink()

    calls_parts: list[Path] = []
    if lifted_parts:
        scan = pl.scan_parquet(lifted_parts)
        targets = scan.select(pl.col("chrom").unique()).collect()["chrom"].to_list()
        for index, chrom in enumerate(sorted(targets, key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))):
            part = spill.stage / f"calls-{index:05d}.parquet"
            canonical_order(scan.filter(pl.col("chrom") == chrom)).collect().write_parquet(part, row_group_size=SPILL_ROW_GROUP)
            calls_parts.append(part)
        for part in lifted_parts:
            part.unlink()
    reasons = {reason: counts[f"reason:{reason}"] for reason in sorted(UNLIFTED_REASONS) if counts[f"reason:{reason}"]}
    warnings = _warnings(counts, reasons, blocks, chain, store, source_digest, target_digest)
    return spill.commit(key, {"calls": calls_parts, "report": report_parts}, warnings)


def _derived_id(parent_id: str, chain_sha256: str, target: str) -> str:
    key = json.dumps(
        {"chain": chain_sha256, "derived_from": parent_id, "operation": RELATION, "reference": target, "version": LIFT_VERSION}, sort_keys=True
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]


def _single(store: ReferenceStore, digest: str) -> dict[str, Contig]:
    """The collection's contigs by canonical name; a name it holds twice is left out."""
    return {chrom: contigs[0] for chrom, contigs in store.contigs_by_chrom(digest).items() if len(contigs) == 1}


def _check_source(chain: Chain, store: ReferenceStore, digest: str, facet_id: str) -> dict[str, bool]:
    """Per source contig the chain covers: whether its size there is the source collection's.

    Raises `ReferenceMismatchError` when a primary contig differs: the chain starts from another assembly.
    """
    contigs = _single(store, digest)
    result: dict[str, bool] = {}
    for chrom, sizes in chain.source_sizes.items():
        ok = chrom in contigs and sizes == {contigs[chrom].length}
        if not ok and chrom in PRIMARY and chrom in contigs:
            raise ReferenceMismatchError(
                f"facet {facet_id}: chain {chain.name} gives contig {chrom} {', '.join(map(str, sorted(sizes)))} bases on its source side, "
                f"the facet's assembly {contigs[chrom].length}: the chain starts from another assembly"
            )
        result[chrom] = ok
    return result


def _primary_fit(chain: Chain, store: ReferenceStore, digest: str) -> tuple[int, int]:
    """(primary contigs whose target size the collection matches, primary contigs where it differs)."""
    contigs = _single(store, digest)
    matched = differing = 0
    for name, sizes in chain.target_sizes.items():
        chrom = canonical_chrom(name)
        if chrom not in PRIMARY or chrom not in contigs:
            continue
        if sizes == {contigs[chrom].length}:
            matched += 1
        else:
            differing += 1
    return matched, differing


def _choose_target(chain: Chain, store: ReferenceStore, target: str | None, facet_id: str) -> str:
    if target is not None:
        if target not in store.collections():
            raise ReferenceMismatchError(f"facet {facet_id}: the store does not hold target collection {target}")
        matched, differing = _primary_fit(chain, store, target)
        if differing or not matched:
            raise ReferenceMismatchError(
                f"facet {facet_id}: chain {chain.name}'s target sizes do not fit collection {target} "
                f"({matched} primary contig(s) match, {differing} differ)"
            )
        return target
    fits = []
    for digest in store.collections():
        matched, differing = _primary_fit(chain, store, digest)
        if matched and not differing:
            fits.append(digest)
    if len(fits) != 1:
        found = ", ".join(fits) if fits else "none"
        raise ReferenceMismatchError(
            f"facet {facet_id}: {len(fits)} collection(s) in the store fit chain {chain.name}'s target side ({found}); pass target= to choose"
        )
    return fits[0]


def _target_contigs(chain: Chain, store: ReferenceStore, digest: str) -> dict[str, Contig]:
    """The target collection's contigs by canonical name, kept only where the chain's size for it agrees."""
    contigs = _single(store, digest)
    sizes: dict[str, set[int]] = {}
    for name, found in chain.target_sizes.items():
        sizes.setdefault(canonical_chrom(name), set()).update(found)
    return {chrom: contig for chrom, contig in contigs.items() if sizes.get(chrom) == {contig.length}}


def _lift_contig(
    rows: pl.DataFrame,
    chrom: str,
    chain: Chain,
    contig_ok: bool | None,
    store: ReferenceStore,
    targets: dict[str, Contig],
    facet_id: str,
    counts: Counter[str],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """One source contig's rows: (the report rows, the lifted `calls` rows)."""
    counts["rows"] += rows.height
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    not_sequence = alts.list.eval(~pl.element().str.to_uppercase().str.contains(_SEQUENCE)).list.any()
    reason = (
        pl.when(pl.col("ref").is_null())
        .then(pl.lit("ref_unread"))
        .when((pl.col("ref_check") == RefCheck.MISMATCH.value).fill_null(False))
        .then(pl.lit("source_ref_mismatch"))
        .when(not_sequence)
        .then(pl.lit("allele_not_sequence"))
        .when(pl.col("split_signal").is_not_null())
        .then(pl.lit("split_flagged"))
    )
    if contig_ok is None:
        reason = reason.otherwise(pl.lit("contig_not_in_chain"))
    elif not contig_ok:
        reason = reason.otherwise(pl.lit("chain_contig_mismatch"))
    start, end = _span(pl.col("pos"), pl.col("ref").str.len_chars())
    work = rows.with_row_index("_row").with_columns(reason.alias("lift_reason"), start.alias("_s0"), end.alias("_e0"))

    blocks = chain.blocks.filter(pl.col("t_chrom") == chrom)
    candidates = work.filter(pl.col("lift_reason").is_null())
    if candidates.height and blocks.height:
        candidates = _place(candidates, blocks, targets)
        candidates = _read_target(candidates, store, targets)
        work = work.join(candidates.select("_row", "_reason2", "_q_chrom", "_q0", "_minus", "_reanchor", "_t_ref", "_swap"), on="_row", how="left", maintain_order="left")
        work = work.with_columns(pl.coalesce("lift_reason", "_reason2").alias("lift_reason"))
    else:
        work = work.with_columns(
            pl.lit(None, dtype=pl.String).alias("_q_chrom"),
            pl.lit(None, dtype=pl.Int64).alias("_q0"),
            pl.lit(None, dtype=pl.Boolean).alias("_minus"),
            pl.lit(None, dtype=pl.Boolean).alias("_reanchor"),
            pl.lit(None, dtype=pl.String).alias("_t_ref"),
            pl.lit(None, dtype=pl.Boolean).alias("_swap"),
        )
        if candidates.height:
            work = work.with_columns(pl.col("lift_reason").fill_null(pl.lit("no_chain_block")))

    for name, count in work.filter(pl.col("lift_reason").is_not_null()).group_by("lift_reason").len().iter_rows():
        counts[f"reason:{name}"] += count
    ok = pl.col("lift_reason").is_null()
    report = work.select(
        *CALLS.names(),
        pl.when(ok).then(pl.col("_q_chrom")).alias("lift_chrom"),
        pl.when(ok).then(_to_pos(pl.col("_q0"))).alias("lift_pos"),
        pl.when(pl.col("_minus").is_not_null()).then(pl.when(pl.col("_minus")).then(pl.lit("-")).otherwise(pl.lit("+"))).alias("lift_strand"),
        pl.when(ok).then(pl.col("_swap")).alias("lift_swapped"),
        pl.col("lift_reason"),
    )
    lifted = _lifted_rows(work.filter(ok), facet_id, chrom, counts)
    return report, lifted


def _place(rows: pl.DataFrame, blocks: pl.DataFrame, targets: dict[str, Contig]) -> pl.DataFrame:
    """Each row's block, target contig and 0-based target start, or why it has none (`_reason2`)."""
    segments = _segments(blocks)
    by_start = rows.select("_row", "_s0", "_e0").sort("_s0")
    first = by_start.join_asof(
        segments.rename({"seg": "seg_s", "depth": "depth_s", "seg_block": "block_s", "covered_cum": "cov_s", "multi_cum": "multi_s"}),
        left_on="_s0",
        right_on="seg_start",
        strategy="backward",
    )
    last = (
        rows.select("_row", (pl.col("_e0") - 1).alias("_last"))
        .sort("_last")
        .join_asof(segments.select(pl.col("seg_start"), pl.col("seg").alias("seg_e"), pl.col("covered_cum").alias("cov_e"), pl.col("multi_cum").alias("multi_e")), left_on="_last", right_on="seg_start", strategy="backward")
        .select("_row", "seg_e", "cov_e", "multi_e")
    )
    placed = first.join(last, on="_row", how="left").sort("_row")
    # Running counts include the segment itself, so a range [s, e] holds cum[e] - cum[s] + (s's own).
    own_cov = (pl.col("depth_s") >= 1).cast(pl.Int64)
    own_multi = (pl.col("depth_s") >= 2).cast(pl.Int64)
    multi = pl.col("multi_e") - pl.col("multi_s") + own_multi
    covered = pl.col("cov_e") - pl.col("cov_s") + own_cov
    one_segment = pl.col("seg_s") == pl.col("seg_e")
    reason = (
        pl.when(multi > 0)
        .then(pl.lit("multiple_targets"))
        .when(one_segment & (pl.col("depth_s") == 1))
        .then(None)
        .when(covered > 0)
        .then(pl.lit("split_across_blocks"))
        .otherwise(pl.lit("no_chain_block"))
    )
    placed = placed.with_columns(reason.alias("_reason2")).join(
        blocks.select(pl.col("block").alias("block_s"), "t_start", "q_chrom", "q_start", "q_size", "q_minus"), on="block_s", how="left"
    )
    length = pl.col("_e0") - pl.col("_s0")
    offset = pl.col("_s0") - pl.col("t_start")
    q0 = pl.when(pl.col("q_minus")).then(pl.col("q_size") - (pl.col("q_start") + offset + length)).otherwise(pl.col("q_start") + offset)
    target = pl.col("q_chrom").replace_strict({name: canonical_chrom(name) for name in placed["q_chrom"].drop_nulls().unique().to_list()}, default=None, return_dtype=pl.String)
    placed = placed.with_columns(
        pl.when(pl.col("_reason2").is_null()).then(target).alias("_q_chrom"),
        pl.when(pl.col("_reason2").is_null()).then(q0).alias("_q0"),
        pl.when(pl.col("_reason2").is_null()).then(pl.col("q_minus")).alias("_minus"),
    )
    placed = placed.with_columns(
        pl.when(pl.col("_reason2").is_null() & ~pl.col("_q_chrom").is_in(list(targets))).then(pl.lit("target_contig_unmatched")).otherwise(pl.col("_reason2")).alias("_reason2")
    )
    return rows.join(placed.select("_row", "_reason2", "_q_chrom", "_q0", "_minus"), on="_row", how="left", maintain_order="left")


def _read_target(rows: pl.DataFrame, store: ReferenceStore, targets: dict[str, Contig]) -> pl.DataFrame:
    """The target's REF over each placed row's lifted span, and whether the row keeps or swaps its alleles.

    A minus-strand indel whose alleles share their first base is re-anchored (`_reanchor`): its span
    moves one base left on the target, so the base read first is the new padding base.
    """
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    lengths = alts.list.eval(pl.element().str.len_chars())
    ref_length = pl.col("ref").str.len_chars()
    indel = ((lengths.list.min() != ref_length) | (lengths.list.max() != ref_length)).fill_null(False)
    first_bases = pl.concat_list(pl.col("ref"), alts).list.eval(pl.element().str.to_uppercase().str.slice(0, 1))
    padded = first_bases.list.n_unique() == 1
    reanchor = pl.col("_reason2").is_null() & pl.col("_minus") & indel & padded
    rows = rows.with_columns(
        pl.when(reanchor & (pl.col("_q0") == 0)).then(pl.lit("indel_minus_strand")).otherwise(pl.col("_reason2")).alias("_reason2"),
        (reanchor & (pl.col("_q0") > 0)).fill_null(False).alias("_reanchor"),
        pl.lit(None, dtype=pl.String).alias("_t_ref"),
    )
    rows = rows.with_columns(pl.when(pl.col("_reanchor")).then(pl.col("_q0") - 1).otherwise(pl.col("_q0")).alias("_q0"))
    readable = rows.filter(pl.col("_reason2").is_null())
    pieces = []
    for chrom in sorted(readable["_q_chrom"].unique().to_list()):
        part = readable.filter(pl.col("_q_chrom") == chrom).select("_row", "_q0", (pl.col("_q0") + pl.col("_e0") - pl.col("_s0")).alias("_q1"))
        spans = part.select("_q0", "_q1").unique().sort("_q0", "_q1")
        bases = store.fetch(targets[chrom], list(zip(spans["_q0"].to_list(), spans["_q1"].to_list())))
        table = spans.with_columns(pl.Series("_bases", [text.upper() for text in bases], dtype=pl.String))
        pieces.append(part.join(table, on=["_q0", "_q1"], how="left").select("_row", "_bases"))
    if pieces:
        read = pl.concat(pieces)
        rows = rows.join(read, on="_row", how="left", maintain_order="left").with_columns(pl.col("_bases").alias("_t_ref")).drop("_bases")
    lifted_ref = _complement(pl.col("ref"), pl.col("_minus"))
    lifted_alts = _complement_list(alts, pl.col("_minus"))
    # A re-anchored row's alleles are the target's first base, then each complemented allele less its
    # last base (the old padding base). The first base read is that new padding base by construction,
    # so the alleles are compared without it.
    reanchored = pl.col("_reanchor")
    target = pl.when(reanchored).then(pl.col("_t_ref").str.slice(1)).otherwise(pl.col("_t_ref"))
    ref_body = pl.when(reanchored).then(_drop_last(lifted_ref)).otherwise(lifted_ref)
    alts_body = pl.when(reanchored).then(lifted_alts.list.eval(_drop_last(pl.element()))).otherwise(lifted_alts)
    keep = target == ref_body
    swap = alts_body.list.contains(target)
    open_ = pl.col("_reason2").is_null()
    reason = (
        pl.when(open_ & pl.col("_t_ref").is_null())
        .then(pl.lit("target_contig_unmatched"))
        .when(open_ & ~keep & ~swap)
        .then(pl.lit("ref_mismatch_after_lift"))
        .otherwise(pl.col("_reason2"))
    )
    return rows.with_columns(reason.alias("_reason2"), pl.when(open_ & pl.col("_t_ref").is_not_null()).then(~keep).alias("_swap"))


def _drop_last(allele: pl.Expr) -> pl.Expr:
    return allele.str.slice(0, allele.str.len_chars() - 1)


def _reanchored(rows: pl.DataFrame) -> pl.DataFrame:
    """Re-anchored rows' alleles on the target: the new padding base (REF's first, read from the
    target), then each complemented allele less its old padding base. Few rows, so plain Python."""

    def moved(anchor: str, allele: str | None) -> str | None:
        return None if allele is None else anchor + allele[:-1]

    columns: dict[str, list] = {name: [] for name in ("_lifted_ref", "alts", "gt_alleles", "source_alleles")}
    for ref, lifted_ref, alts, gt_alleles, source_alleles in rows.select("ref", "_lifted_ref", "alts", "gt_alleles", "source_alleles").iter_rows():
        anchor = ref[0]
        columns["_lifted_ref"].append(moved(anchor, lifted_ref))
        columns["alts"].append([moved(anchor, allele) for allele in alts])
        for name, alleles in (("gt_alleles", gt_alleles), ("source_alleles", source_alleles)):
            columns[name].append(None if alleles is None else [moved(anchor, allele) for allele in alleles])
    return rows.with_columns(
        pl.Series("_lifted_ref", columns["_lifted_ref"], dtype=pl.String),
        *[pl.Series(name, columns[name], dtype=CALLS[name]) for name in ("alts", "gt_alleles", "source_alleles")],
    )


def _lifted_rows(work: pl.DataFrame, facet_id: str, source_chrom: str, counts: Counter[str]) -> pl.DataFrame:
    """The lifted `calls` rows from placed, read, kept-or-swapped parent rows."""
    if work.is_empty():
        return pl.DataFrame(schema=CALLS)
    minus = pl.col("_minus")
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    lifted_ref = _complement(pl.col("ref"), minus)
    moved_contig = pl.col("_q_chrom") != source_chrom
    array = pl.col("source_ref").is_null()
    counts["lifted"] += work.height
    counts["minus"] += int(work["_minus"].sum())
    counts["reanchored"] += int(work["_reanchor"].sum())
    counts["swapped"] += int(work["_swap"].sum())
    counts["other_contig"] += int(work.select(moved_contig.sum()).item())
    counts["dosage_withheld"] += int(work.select((pl.col("_swap") & pl.col("dosage").is_not_null()).sum()).item())
    counts["phase_cleared"] += int(work.select((moved_contig & pl.col("phased").is_not_null()).sum()).item())
    out = work.with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.col("_q_chrom").alias("chrom"),
        _to_pos(pl.col("_q0")).alias("pos"),
        (pl.col("_q0") + pl.col("ref").str.len_chars()).cast(pl.UInt32).alias("end_pos"),
        pl.col("_t_ref").alias("ref"),
        lifted_ref.alias("_lifted_ref"),
        _complement_list(alts, minus).alias("alts"),
        pl.when(pl.col("gt_alleles").is_not_null()).then(_complement_list(pl.col("gt_alleles"), minus)).alias("gt_alleles"),
        pl.when(pl.col("source_alleles").is_not_null()).then(_complement_list(pl.col("source_alleles"), minus)).alias("source_alleles"),
        pl.lit(None, dtype=CALLS["source_ref"]).alias("source_ref"),
        pl.lit(None, dtype=CALLS["ref_check"]).alias("ref_check"),
        pl.lit(None, dtype=CALLS["source_gt"]).alias("source_gt"),
        pl.lit(None, dtype=CALLS["vrs_allele_ids"]).alias("vrs_allele_ids"),
        pl.lit(None, dtype=CALLS["vrs_errors"]).alias("vrs_errors"),
        pl.when(pl.col("_swap")).then(None).otherwise(pl.col("dosage")).alias("dosage"),
        *[pl.when(moved_contig).then(None).otherwise(pl.col(name)).alias(name) for name in ("phased", "phase_sets", "phase_set_kinds")],
        array.alias("_array"),
    )
    if out["_reanchor"].any():
        out = pl.concat([out.filter(~pl.col("_reanchor")), _reanchored(out.filter(pl.col("_reanchor")))], how="vertical")
    swapped = out.filter(pl.col("_swap") & ~pl.col("_array"))
    if swapped.height:
        out = pl.concat([out.filter(~(pl.col("_swap") & ~pl.col("_array"))), _swap(swapped)], how="vertical")
    # Array rows (no source REF) take their ALTs from the genotype against the new REF, as resolve_ref
    # derives them: the observed alleles other than REF, sorted, with `gt` indexing `[ref] + alts`.
    redo = pl.col("_array") & (pl.col("_swap") | pl.col("_minus"))
    if out.select(redo.any()).item():
        out = _derive_array_alleles(out.with_columns(redo.alias("_array")), Counter())
    return out.select(CALLS.names())


def _swap(rows: pl.DataFrame) -> pl.DataFrame:
    """VCF rows whose target REF is one of their ALTs: that ALT becomes REF and the old REF takes its
    place among the ALTs; the genotype's indices follow. Few rows, so plain Python."""
    alts_out: list[list[str]] = []
    gts_out: list[list[int | None] | None] = []
    for old_ref, alts, gt, new_ref in rows.select("_lifted_ref", "alts", "gt", "ref").iter_rows():
        k = alts.index(new_ref) + 1
        swapped = list(alts)
        swapped[k - 1] = old_ref
        alts_out.append(swapped)
        gts_out.append(None if gt is None else [None if g is None else k if g == 0 else 0 if g == k else g for g in gt])
    return rows.with_columns(pl.Series("alts", alts_out, dtype=CALLS["alts"]), pl.Series("gt", gts_out, dtype=CALLS["gt"]))


def _warnings(counts: Counter[str], reasons: dict[str, int], blocks: int, chain: Chain, store: ReferenceStore, source: str, target: str) -> list[str]:
    source_label = store.assembly_label(source) or source
    target_label = store.assembly_label(target) or target
    lines = [
        f"{counts['lifted']} of {counts['rows']} row(s) lifted from {source_label} to {target_label} through chain {chain.name} "
        f"(sha256 {chain.sha256}); REF read from collection {target}"
    ]
    if counts["swapped"]:
        lines.append(f"{counts['swapped']} lifted row(s) had the target's REF as an ALT: REF and that ALT swapped, the genotype re-indexed")
    if counts["minus"]:
        lines.append(f"{counts['minus']} lifted row(s) on minus-strand chain blocks: alleles reverse-complemented")
    if counts["reanchored"]:
        lines.append(
            f"{counts['reanchored']} of them indels re-anchored on the target's base before the lifted span (one base left); "
            "run normalize on the lifted facet to left-align them on the target"
        )
    if counts["dosage_withheld"]:
        lines.append(f"{counts['dosage_withheld']} swapped row(s): dosage withheld (it counted the other allele)")
    if counts["other_contig"]:
        lines.append(f"{counts['other_contig']} row(s) lifted to another contig than their source's; phase sets cleared on {counts['phase_cleared']}")
    for reason, count in reasons.items():
        lines.append(f"{count} row(s) not lifted: {REASON_TEXT[reason]} [{reason}]")
    if blocks:
        lines.append(f"{blocks} block(s) not lifted (lifting blocks is not built); the lifted facet's absence_means is unknown")
    return lines
