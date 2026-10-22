"""Normalize a facet: split multi-allelic sites, left-align and trim indels, dedup, sort (RM6).

`normalize(loaded, reference)` returns a derived facet (charter 5: `relation = normalized`,
`derived_from` the parent, which is never changed). It needs a facet whose REF was read
(`resolve_ref`), and reads further bases from the same collection of the same `ReferenceStore`.
The rules are in docs/SCHEMA.md ("Normalizing"); in short:

1. **Split.** A row with two or more ALTs becomes one row per ALT, as `bcftools norm -m-any
   --multi-overlaps .`, vt decompose and Glow do: a GT slot holding another ALT is written null
   (`.`), never `0`, so a split row never claims a REF haplotype the person does not have. Each split
   row carries `split_signal = marker`. `ad` is sliced to `[REF, ALT]`; `gl` is sliced to the
   genotypes over `[REF, ALT]` (haploid or diploid) only where no GT slot holds another ALT, because
   otherwise the slice leaves out the genotype the caller chose; `gp` and `dosage` are withheld (a
   posterior over part of the genotypes no longer sums to 1, and making it do so changes the
   values). A site holding a symbolic allele (`<DEL>`, `<*>`, a breakend) is left as written.
2. **Left-align and trim** each biallelic sequence row whose REF matched the reference
   (`ref_check = match`), with the standard algorithm (Tan et al. 2015, vt, `bcftools norm -f`):
   drop a common last base, and when an allele runs empty, extend both to the left with the
   reference base before them; repeat; then drop common first bases while both alleles keep one. A
   row that moved or changed has `source_ref`, `ref_check`, `source_gt` and its VRS ids null: they
   described the representation the source wrote, at the position it wrote it, and keeping
   `ref_check = match` on bases read from the very sequence it checks would be a check agreeing with
   itself. `ref` holds the new REF, read from the reference. Rows whose REF was not shown to match
   pass through unchanged, counted.
3. **Dedup and sort.** Rows identical in every column are kept once, counted. Rows for the same
   allele at the same position that differ in anything else are all kept, with a warning
   (charter 8: disagreement is data, never resolved by picking a winner). The result is in canonical
   order (SCHEMA.md Decision 11).

Rows the loader flagged as one record of a pre-split site (RM270) are left as written: a `0` there
may mean another ALT, and moving one row of such a site would keep `rejoin_split` from finding it.
Run `rejoin_split` before `resolve_ref` and `normalize`.

Coordinates: `pos` is 1-based; the store takes 0-based half-open ranges, converted in `_window`.
"""

import hashlib
import json
import logging
import shutil
from collections import Counter
from pathlib import Path

import polars as pl

from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.reference import Contig, ReferenceMismatchError, ReferenceStore
from just_vcf.schema import CALLS, CONTIG_RANK, SCHEMA_VERSION, RefCheck, SplitSignal, canonical_chrom, canonical_order, empty

logger = logging.getLogger(__name__)

NORMALIZE_VERSION = 1
"""Part of the derived `facet_id`: bumped whenever the rule changes what `normalize` emits."""

RELATION = "normalized"
"""The derived facet's `relation`."""

_FIRST_WINDOW = 100
"""Reference bases read to the left of an indel on the first pass; a row that needs more is read again
with a window `_WINDOW_GROWTH` times wider, until it reaches the contig's start."""
_WINDOW_GROWTH = 16

_POSITIONAL = ("source_ref", "ref_check", "source_gt", "vrs_allele_ids", "vrs_errors")
"""Columns that describe the source's representation at the source's position: null on a row that
left-alignment or trimming changed."""


def normalize(loaded: Loaded, reference: "ReferenceStore | Path | str | None" = None) -> Loaded:
    """A derived facet with multi-allelic sites split, indels left-aligned and trimmed, exact
    duplicates dropped, and rows in canonical order.

    `loaded` holds one facet whose REF was read by `resolve_ref`. `reference` is the store it was read
    from (a `ReferenceStore` or its directory; default: the store at `reference_dir()`); it must hold
    the collection the facet names in `ref_seqcol_digest`. The input is not changed. What could not
    be normalized is left as written, each reason counted once in the warnings.
    """
    parent = loaded.facets.collect()
    if parent.height != 1:
        raise ValueError(f"normalize takes one facet; this Loaded has {parent.height}")
    facet = parent.row(0, named=True)
    digest = facet["ref_seqcol_digest"]
    if digest is None:
        raise ValueError(f"facet {facet['facet_id']}: REF has not been read; run resolve_ref first, then normalize")
    store = reference if isinstance(reference, ReferenceStore) else ReferenceStore.open(reference)
    if digest not in store.collections():
        raise ReferenceMismatchError(f"facet {facet['facet_id']}: the store does not hold collection {digest}, which REF was read from")
    contigs = _contigs_by_chrom(store, digest)

    facet_id = _derived_id(facet["facet_id"], digest)
    spill = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
    if spill.exists():
        shutil.rmtree(spill)
    spill.mkdir(parents=True)
    # Streamed to parquet once, so each contig below is a filtered scan, not a replay of the input's plan.
    loaded.calls.sink_parquet(spill / "input.parquet", row_group_size=SPILL_ROW_GROUP)
    source = pl.scan_parquet(spill / "input.parquet")
    chroms = source.select(pl.col("chrom").unique()).collect()["chrom"].to_list()

    counts: Counter[str] = Counter()
    parts: list[Path] = []
    for index, chrom in enumerate(sorted(chroms, key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))):
        rows = source.filter(pl.col("chrom") == chrom).collect()
        normalized = _normalize_contig(rows.with_columns(pl.lit(facet_id).alias("facet_id")), store, contigs.get(chrom), counts)
        part = spill / f"calls-{index:05d}.parquet"
        normalized.write_parquet(part, row_group_size=SPILL_ROW_GROUP)
        parts.append(part)
    (spill / "input.parquet").unlink()

    warnings = _warnings(counts)
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    facets = parent.lazy().with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.lit(facet["facet_id"]).alias("derived_from"),
        pl.lit(RELATION).alias("relation"),
    )
    return Loaded(
        facets=facets,
        calls=pl.scan_parquet(parts) if parts else empty(CALLS),
        blocks=loaded.blocks.with_columns(pl.lit(facet_id).alias("facet_id")),
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )


def _derived_id(parent_id: str, digest: str) -> str:
    key = json.dumps(
        {"derived_from": parent_id, "operation": RELATION, "reference": digest, "version": NORMALIZE_VERSION}, sort_keys=True
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]


def _contigs_by_chrom(store: ReferenceStore, digest: str) -> dict[str, Contig]:
    """The collection's contigs by canonical name; a name the collection holds twice is left out."""
    found: dict[str, list[Contig]] = {}
    for contig in store.contigs(digest):
        found.setdefault(canonical_chrom(contig.name), []).append(contig)
    return {chrom: matches[0] for chrom, matches in found.items() if len(matches) == 1}


def _normalize_contig(rows: pl.DataFrame, store: ReferenceStore, contig: Contig | None, counts: Counter[str]) -> pl.DataFrame:
    """One contig's rows split, aligned, deduplicated and sorted."""
    if rows.is_empty():
        return rows
    flagged = pl.col("split_signal").is_not_null()
    held = rows.filter(flagged)
    counts["flagged"] += held.height
    work = _align(_split(rows.filter(~flagged), counts), store, contig, counts)
    ordered = canonical_order(pl.concat([work, held], how="vertical").lazy()).collect()
    unique = ordered.unique(maintain_order=True)
    counts["duplicates"] += ordered.height - unique.height
    written_ref = pl.coalesce(pl.col("ref"), pl.col("source_ref"))
    conflicts = (
        unique.group_by(pl.col("pos"), written_ref.alias("_ref"), pl.col("alts"))
        .agg(pl.len().alias("rows"))
        .filter(pl.col("rows") > 1)
    )
    counts["conflict_sites"] += conflicts.height
    counts["conflict_rows"] += int(conflicts["rows"].sum()) if conflicts.height else 0
    return unique


_SYMBOLIC = r"^<|\[|\]"
"""A symbolic allele (`<DEL>`, `<*>`) or a breakend."""


def _slots(ploidy: int) -> list[pl.Expr]:
    return [pl.col("gt").list.get(slot, null_on_oob=True) for slot in range(ploidy)]


def _per_slot(values: list[pl.Expr]) -> pl.Expr:
    """A list of one value per GT slot, as long as the row's own `gt`; null where `gt` is null."""
    return pl.when(pl.col("gt").is_not_null()).then(pl.concat_list(values).list.head(pl.col("gt").list.len()))


def _split(rows: pl.DataFrame, counts: Counter[str]) -> pl.DataFrame:
    """Every multi-allelic row as one row per ALT; rows that are biallelic, or cannot be split, as they are."""
    multi = (pl.col("alts").list.len() >= 2).fill_null(False)
    symbolic = pl.col("alts").list.eval(pl.element().str.contains(_SYMBOLIC)).list.any()
    unindexed = pl.col("gt").is_null() & pl.col("gt_alleles").is_not_null()
    n_symbolic, n_unindexed = rows.select((multi & symbolic).sum().alias("a"), (multi & ~symbolic & unindexed).sum().alias("b")).row(0)
    counts["unsplit_symbolic"] += n_symbolic
    counts["unsplit_unindexed"] += n_unindexed
    target = multi & ~symbolic & ~unindexed
    work = rows.filter(target)
    if work.is_empty():
        return rows
    counts["split_sites"] += work.height
    ploidy = work.select(pl.col("gt").list.len().max()).item() or 0
    n = pl.col("alts").list.len().cast(pl.Int64)
    k = pl.col("_k")
    slots = _slots(ploidy)
    pieces = work.with_columns(pl.int_ranges(1, n + 1, dtype=pl.Int64).alias("_k")).explode("_k", empty_as_null=False)
    counts["split_rows"] += pieces.height
    alt = pl.col("alts").list.get(k - 1)
    other_alt = pl.any_horizontal([(g.is_not_null() & (g != 0) & (g != k)).fill_null(False) for g in slots]) if slots else pl.lit(False)
    indexable = (pl.col("source_alleles") == pl.concat_list(pl.col("source_ref"), pl.col("alts"))).fill_null(False)
    pieces = pieces.with_columns(other_alt.alias("_other"), indexable.alias("_indexable"))

    def pair(name: str) -> pl.Expr:
        column = pl.col(name)
        return pl.when(column.list.len() == n + 1).then(pl.concat_list(column.list.get(0), column.list.get(k)))

    het = k * (k + 1) // 2
    gl = pl.col("gl")
    usable = pl.col("_indexable") & ~pl.col("_other")
    sliced_gl = (
        pl.when(usable & (pl.col("gt").list.len() == 1) & (gl.list.len() == n + 1))
        .then(pl.concat_list(gl.list.get(0), gl.list.get(k)))
        .when(usable & (pl.col("gt").list.len() == 2) & (gl.list.len() == (n + 1) * (n + 2) // 2))
        .then(pl.concat_list(gl.list.get(0), gl.list.get(het), gl.list.get(het + k)))
    )
    ref = pl.coalesce(pl.col("ref"), pl.col("source_ref"))
    gt = _per_slot([pl.when(g == 0).then(0).when(g == k).then(1).cast(pl.UInt8) for g in slots]) if slots else pl.col("gt")
    gt_alleles = (
        _per_slot([pl.when((g == 0) | (g == k)).then(pl.col("gt_alleles").list.get(s, null_on_oob=True)) for s, g in enumerate(slots)])
        if slots
        else pl.col("gt_alleles")
    )
    pieces = pieces.with_columns(
        pl.concat_list(alt).alias("alts"),
        (pl.col("pos").cast(pl.Int64) + ref.str.len_chars().cast(pl.Int64) - 1).cast(pl.UInt32).alias("end_pos"),
        gt.alias("gt"),
        gt_alleles.alias("gt_alleles"),
        pair("vrs_allele_ids").alias("vrs_allele_ids"),
        pair("vrs_errors").alias("vrs_errors"),
        pl.when(pl.col("_indexable")).then(pair("ad")).alias("ad"),
        sliced_gl.alias("gl"),
        pl.lit(None, dtype=CALLS["gp"]).alias("gp"),
        pl.lit(None, dtype=CALLS["dosage"]).alias("dosage"),
        pl.lit(None, dtype=CALLS["source_gt"]).alias("source_gt"),
        pl.lit(SplitSignal.MARKER.value, dtype=CALLS["split_signal"]).alias("split_signal"),
        pl.col("gl").is_not_null().alias("_had_gl"),
        (pl.col("gl").is_not_null() | pl.col("ad").is_not_null()).alias("_had_vector"),
        pl.col("gp").is_not_null().alias("_had_gp"),
    )
    pieces = pieces.with_columns(
        pl.when(pl.col("ad").is_not_null() | pl.col("gl").is_not_null())
        .then(pl.concat_list(pl.col("source_alleles").list.get(0), pl.col("alts").list.first()))
        .alias("source_alleles")
    )
    tallies = pieces.select(
        (pl.col("_indexable") & pl.col("_had_gl") & pl.col("_other")).sum().alias("a"),
        (pl.col("_indexable") & pl.col("_had_gl") & ~pl.col("_other") & pl.col("gl").is_null()).sum().alias("b"),
        (~pl.col("_indexable") & pl.col("_had_vector")).sum().alias("c"),
        pl.col("_had_gp").sum().alias("d"),
    ).row(0)
    for key, value in zip(("gl_other_alt", "gl_unsliceable", "vectors_unindexed", "gp_withheld"), tallies):
        counts[key] += value
    return pl.concat([rows.filter(~target), pieces.select(CALLS.names())], how="vertical")


def _align(rows: pl.DataFrame, store: ReferenceStore, contig: Contig | None, counts: Counter[str]) -> pl.DataFrame:
    """Every biallelic sequence row whose REF matched, left-aligned and trimmed; the rest as they are,
    counted by why they were not aligned."""
    alt = pl.col("alts").list.first()
    written = pl.coalesce(pl.col("ref"), pl.col("source_ref"))
    shape = (
        (pl.col("alts").list.len() == 1)
        & alt.str.to_uppercase().str.contains(r"^[ACGTN]+$")
        & ((written.str.len_chars() > 1) | (alt.str.len_chars() > 1))
    ).fill_null(False)
    reason = (
        pl.when((pl.col("ref_check") != RefCheck.MATCH.value).fill_null(True) | pl.col("ref").is_null())
        .then(pl.col("ref_check").cast(pl.String).fill_null("unchecked"))
        .when(pl.col("gt").is_null() & pl.col("gt_alleles").is_not_null())
        .then(pl.lit("unindexed"))
        .when((pl.col("gt").list.max() > 1).fill_null(False))
        .then(pl.lit("index_past_alts"))
    )
    # An ALT that is neither bases, `*`, symbolic nor a breakend (a `-` written for a deletion) is not a
    # VCF allele; it cannot be aligned, and is counted rather than passed over in silence.
    not_an_allele = ((pl.col("alts").list.len() == 1) & ~alt.str.to_uppercase().str.contains(r"^([ACGTN]+|\*)$") & ~alt.str.contains(_SYMBOLIC)).fill_null(False)
    rows = rows.with_row_index("_row").with_columns(
        pl.when(shape).then(reason).when(not_an_allele).then(pl.lit("not_an_allele")).alias("_reason"), shape.alias("_shape")
    )
    for name, count in rows.filter(pl.col("_reason").is_not_null()).group_by("_reason").len().iter_rows():
        counts[f"unaligned_{name}"] += count
    candidates = rows.filter(pl.col("_shape") & pl.col("_reason").is_null()).select("_row", "pos", "ref", alt.str.to_uppercase().alias("_alt"))
    rows = rows.drop("_reason", "_shape")
    if candidates.is_empty():
        return rows.drop("_row")
    if contig is None:
        # REF matched, so the contig was read once; a collection naming it twice is not read again.
        counts["unaligned_contig"] += candidates.height
        return rows.drop("_row")

    changed: dict[str, list] = {"_row": [], "_pos": [], "_ref": [], "_alt": []}
    pending = list(candidates.iter_rows())
    width = _FIRST_WINDOW
    while pending:
        windows = _windows(store, contig, [(pos, width) for _, pos, _, _ in pending])
        retry = []
        for (index, pos, ref, alt_text), (window, start) in zip(pending, windows):
            result = _left_align(pos, ref, alt_text, window, start)
            if result is None:
                retry.append((index, pos, ref, alt_text))
                continue
            new_pos, new_ref, new_alt, at_start = result
            counts["aligned_at_start"] += at_start
            if (new_pos, new_ref, new_alt) != (pos, ref, alt_text):
                for key, value in zip(changed, (index, new_pos, new_ref, new_alt)):
                    changed[key].append(value)
        pending = retry
        width *= _WINDOW_GROWTH
    if not changed["_row"]:
        return rows.drop("_row")
    moved = pl.DataFrame(changed, schema={"_row": pl.UInt32, "_pos": pl.UInt32, "_ref": pl.String, "_alt": pl.String})
    return _moved(rows.join(moved, on="_row", how="left", maintain_order="left"), counts).drop("_row")


def _windows(store: ReferenceStore, contig: Contig, requests: list[tuple[int, int]]) -> list[tuple[str, int]]:
    """For each `(pos, width)`: up to `width` reference bases ending just before 1-based `pos`, and the
    1-based position of the window's first base. A window at the contig's start is empty."""
    starts = [max(1, pos - width) for pos, width in requests]
    # 1-based [start, pos - 1] is 0-based half-open [start - 1, pos - 1).
    ranges = [(start - 1, pos - 1) for start, (pos, _) in zip(starts, requests)]
    reads = [index for index, (begin, end) in enumerate(ranges) if end > begin]
    bases = store.fetch(contig, [ranges[index] for index in reads])
    windows = [""] * len(requests)
    for index, text in zip(reads, bases):
        windows[index] = text.upper()
    return list(zip(windows, starts))


def _left_align(pos: int, ref: str, alt: str, window: str, start: int) -> tuple[int, str, str, bool] | None:
    """`(pos, REF, ALT, reached the contig start)` left-aligned and trimmed, or None when the reference
    `window` (bases at 1-based `start` .. `pos - 1`) ran out before the allele stopped moving.

    At the contig's start an allele that would run empty keeps its last base instead, so a deletion
    at position 1 stays anchored on the base after it, as VCF allows.
    """
    r, a, p = ref, alt, pos
    at_start = False
    while True:
        changed = False
        if r and a and r[-1] == a[-1]:
            if (len(r) == 1 or len(a) == 1) and p == 1:
                at_start = True
            else:
                r, a = r[:-1], a[:-1]
                changed = True
        if not r or not a:
            if p - 1 < start:
                return None
            base = window[p - 1 - start]
            r, a, p = base + r, base + a, p - 1
            changed = True
        if not changed:
            break
    while len(r) >= 2 and len(a) >= 2 and r[0] == a[0]:
        r, a, p = r[1:], a[1:], p + 1
    return p, r, a, at_start


def _moved(rows: pl.DataFrame, counts: Counter[str]) -> pl.DataFrame:
    """Rows with a new position or spelling (`_pos`, `_ref`, `_alt`) rewritten to it, positional columns
    cleared; the vectors follow when they were indexed over the row's own two alleles."""
    hit = pl.col("_pos").is_not_null()
    follow = (
        pl.col("source_alleles").list.eval(pl.element().str.to_uppercase())
        == pl.concat_list(pl.col("ref"), pl.col("alts").list.first().str.to_uppercase())
    ).fill_null(False)
    rows = rows.with_columns(follow.alias("_follow"))
    moved, trimmed, unindexed = rows.select(
        (hit & (pl.col("_pos") != pl.col("pos"))).sum().alias("a"),
        (hit & (pl.col("_pos") == pl.col("pos"))).sum().alias("b"),
        (hit & pl.col("source_alleles").is_not_null() & ~pl.col("_follow")).sum().alias("c"),
    ).row(0)
    counts["aligned_moved"] += moved
    counts["aligned_trimmed"] += trimmed
    counts["vectors_unindexed"] += unindexed
    ploidy = rows.filter(hit).select(pl.col("gt").list.len().max()).item() or 0
    slots = _slots(ploidy)
    gt_alleles = _per_slot([pl.when(g == 0).then(pl.col("_ref")).when(g == 1).then(pl.col("_alt")) for g in slots]) if slots else pl.col("gt_alleles")
    drop_vectors = hit & ~pl.col("_follow")
    return rows.with_columns(
        pl.when(hit).then(pl.col("_pos")).otherwise(pl.col("pos")).alias("pos"),
        pl.when(hit)
        .then((pl.col("_pos").cast(pl.Int64) + pl.col("_ref").str.len_chars().cast(pl.Int64) - 1).cast(pl.UInt32))
        .otherwise(pl.col("end_pos"))
        .alias("end_pos"),
        pl.when(hit).then(pl.col("_ref")).otherwise(pl.col("ref")).alias("ref"),
        pl.when(hit).then(pl.concat_list(pl.col("_alt"))).otherwise(pl.col("alts")).alias("alts"),
        pl.when(hit).then(gt_alleles).otherwise(pl.col("gt_alleles")).alias("gt_alleles"),
        pl.when(hit & pl.col("_follow"))
        .then(pl.concat_list(pl.col("_ref"), pl.col("_alt")))
        .when(hit)
        .then(None)
        .otherwise(pl.col("source_alleles"))
        .alias("source_alleles"),
        *[pl.when(drop_vectors).then(None).otherwise(pl.col(name)).alias(name) for name in ("gl", "gp", "ad")],
        *[pl.when(hit).then(None).otherwise(pl.col(name)).alias(name) for name in _POSITIONAL],
    ).drop("_pos", "_ref", "_alt", "_follow")


def _warnings(counts: Counter[str]) -> list[str]:
    lines = []
    if counts["split_sites"]:
        lines.append(
            f"{counts['split_sites']} multi-allelic site(s) split into {counts['split_rows']} row(s), one per ALT; "
            "a GT slot holding another ALT is written `.`, and each row carries split_signal = marker"
        )
    if counts["gp_withheld"]:
        lines.append(f"{counts['gp_withheld']} split row(s): GP withheld (a posterior over part of the genotypes is not one)")
    if counts["gl_other_alt"]:
        lines.append(
            f"{counts['gl_other_alt']} split row(s): GL withheld, since the genotype holds another ALT and the slice "
            "would leave out the genotype the caller chose"
        )
    if counts["gl_unsliceable"]:
        lines.append(f"{counts['gl_unsliceable']} split row(s): GL withheld (not haploid or diploid, or not one value per genotype)")
    if counts["vectors_unindexed"]:
        lines.append(f"{counts['vectors_unindexed']} row(s): AD/GL/GP withheld, indexed over another allele list than the row's")
    moved, trimmed = counts["aligned_moved"], counts["aligned_trimmed"]
    if moved or trimmed:
        lines.append(f"{moved} row(s) left-aligned to another position and {trimmed} trimmed in place, against the reference")
    if counts["aligned_at_start"]:
        lines.append(f"{counts['aligned_at_start']} row(s) reached their contig's start and stay anchored on the base after")
    reasons = {
        "unaligned_not_run": "REF not compared with the reference (ref_check not_run)",
        "unaligned_mismatch": "REF differs from the reference (ref_check mismatch)",
        "unaligned_unchecked": "no source REF was checked",
        "unaligned_no_ref": "no REF at all",
        "unaligned_unindexed": "a genotype that is not indexed",
        "unaligned_index_past_alts": "a genotype index past the row's ALTs",
        "unaligned_contig": "a contig the collection names twice",
        "unaligned_not_an_allele": "an ALT that is not a VCF allele (such as `-` for a deletion, F12)",
    }
    for key, reason in reasons.items():
        if counts[key]:
            lines.append(f"{counts[key]} row(s) not left-aligned or trimmed: {reason}")
    if counts["unsplit_symbolic"]:
        lines.append(f"{counts['unsplit_symbolic']} multi-allelic row(s) left as written: a symbolic allele (`<*>`, `<DEL>`, a breakend)")
    if counts["unsplit_unindexed"]:
        lines.append(f"{counts['unsplit_unindexed']} multi-allelic row(s) left as written: a genotype that is not indexed")
    if counts["flagged"]:
        lines.append(
            f"{counts['flagged']} row(s) flagged as one record of a pre-split site (RM270) left as written; "
            "run rejoin_split before resolve_ref and normalize"
        )
    if counts["duplicates"]:
        lines.append(f"{counts['duplicates']} exact duplicate row(s) dropped")
    if counts["conflict_sites"]:
        lines.append(
            f"{counts['conflict_rows']} row(s) at {counts['conflict_sites']} position(s) state the same allele and differ "
            "otherwise; all kept, none chosen"
        )
    return lines
