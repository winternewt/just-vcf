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

_SEQUENCE = frozenset("ACGTN")

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
    counts["flagged"] += rows.select(flagged.sum()).item()
    shaped = (pl.col("alts").list.len() >= 2) | _may_change()
    work_mask = shaped.fill_null(False) & ~flagged
    work = rows.filter(work_mask)
    rest = rows.filter(~work_mask)

    pieces: list[dict] = []
    for row in work.iter_rows(named=True):
        pieces += _split(row, counts)
    aligned = _align_all(pieces, store, contig, counts)
    combined = pl.concat([rest, pl.DataFrame(aligned, schema=CALLS)], how="vertical") if aligned else rest
    ordered = canonical_order(combined.lazy()).collect()
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


def _may_change() -> pl.Expr:
    """True on a biallelic row that left-alignment or trimming could change: an allele longer than one base."""
    ref = pl.coalesce(pl.col("ref"), pl.col("source_ref"))
    alt = pl.col("alts").list.first()
    return (pl.col("alts").list.len() == 1) & ((ref.str.len_chars() > 1) | (alt.str.len_chars() > 1))


def _symbolic(allele: str) -> bool:
    return allele.startswith("<") or "[" in allele or "]" in allele


def _split(row: dict, counts: Counter[str]) -> list[dict]:
    """The row as one row per ALT; a row that is already biallelic, or cannot be split, as itself."""
    alts = row["alts"]
    if len(alts) < 2:
        return [row]
    if any(_symbolic(allele) for allele in alts):
        counts["unsplit_symbolic"] += 1
        return [row]
    if row["gt"] is None and row["gt_alleles"] is not None:
        counts["unsplit_unindexed"] += 1
        return [row]
    counts["split_sites"] += 1
    gt = row["gt"]
    ploidy = len(gt) if gt is not None else None
    ref = row["ref"] if row["ref"] is not None else row["source_ref"]
    alleles = row["source_alleles"]
    pieces = []
    for k, alt in enumerate(alts, start=1):
        piece = dict(row)
        other_alt = gt is not None and any(index is not None and index not in (0, k) for index in gt)
        piece.update(
            alts=[alt],
            end_pos=row["pos"] + len(ref) - 1,
            gt=None if gt is None else [None if index is None else 0 if index == 0 else 1 if index == k else None for index in gt],
            gt_alleles=None
            if row["gt_alleles"] is None or gt is None
            else [allele if index in (0, k) else None for allele, index in zip(row["gt_alleles"], gt)],
            source_gt=None,
            vrs_allele_ids=_pair(row["vrs_allele_ids"], k, len(alts)),
            vrs_errors=_pair(row["vrs_errors"], k, len(alts)),
            gp=None,
            dosage=None,
            split_signal=SplitSignal.MARKER.value,
        )
        slot = alleles.index(alt) if alleles is not None and alleles.count(alt) == 1 and alleles.index(alt) > 0 else None
        ad = _pair(row["ad"], slot, len(alleles) - 1) if slot is not None else None
        gl = None
        if slot is not None and row["gl"] is not None:
            if other_alt:
                counts["gl_other_alt"] += 1
            else:
                gl = _genotype_slice(row["gl"], slot, len(alleles), ploidy)
                if gl is None:
                    counts["gl_unsliceable"] += 1
        elif row["gl"] is not None or row["ad"] is not None:
            counts["vectors_unindexed"] += 1
        if row["gp"] is not None:
            counts["gp_withheld"] += 1
        piece.update(ad=ad, gl=gl, source_alleles=[alleles[0], alt] if ad is not None or gl is not None else None)
        pieces.append(piece)
    counts["split_rows"] += len(pieces)
    return pieces


def _pair(values: list | None, k: int, alt_count: int) -> list | None:
    """`[values[0], values[k]]` from a list parallel to REF plus `alt_count` ALTs; null when it is not one."""
    if values is None or len(values) != alt_count + 1:
        return None
    return [values[0], values[k]]


def _genotype_slice(vector: list, slot: int, allele_count: int, ploidy: int | None) -> list | None:
    """The entries of a genotype-ordered vector (VCF Number=G) for the genotypes over alleles 0 and `slot`.

    Haploid: genotypes `0`, `slot`. Diploid: `0/0`, `0/slot`, `slot/slot`, at VCF's index
    `j + k(k+1)/2` for `j <= k`. Null for any other ploidy, or a vector of another length.
    """
    if ploidy == 1 and len(vector) == allele_count:
        return [vector[0], vector[slot]]
    if ploidy == 2 and len(vector) == allele_count * (allele_count + 1) // 2:
        het = slot * (slot + 1) // 2
        return [vector[0], vector[het], vector[het + slot]]
    return None


def _align_all(rows: list[dict], store: ReferenceStore, contig: Contig | None, counts: Counter[str]) -> list[dict]:
    """Every row that can be left-aligned and trimmed, done; the rest as they are, counted by reason."""
    pending: list[tuple[int, int, str, str]] = []  # (row index, pos, REF, ALT), uppercase
    for index, row in enumerate(rows):
        reason = _unalignable(row)
        if reason is None:
            pending.append((index, row["pos"], row["ref"], row["alts"][0].upper()))
        elif reason != "fixed":
            counts[f"unaligned_{reason}"] += 1
    if pending and contig is None:
        # REF matched, so the contig was read once; a collection naming it twice cannot be read again.
        counts["unaligned_contig"] += len(pending)
        return rows
    width = _FIRST_WINDOW
    while pending:
        windows = _windows(store, contig, [(pos, width) for _, pos, _, _ in pending])
        retry = []
        for (index, pos, ref, alt), (window, start) in zip(pending, windows):
            result = _left_align(pos, ref, alt, window, start)
            if result is None:
                retry.append((index, pos, ref, alt))
                continue
            new_pos, new_ref, new_alt, at_start = result
            counts["aligned_at_start"] += at_start
            if (new_pos, new_ref, new_alt) != (pos, ref, alt):
                rows[index] = _moved(rows[index], new_pos, new_ref, new_alt, counts)
        pending = retry
        width *= _WINDOW_GROWTH
    return rows


def _unalignable(row: dict) -> str | None:
    """Why a row is not left-aligned; `fixed` when there is nothing to align; None when it can be."""
    alts = row["alts"]
    if len(alts) != 1:
        return "multiallelic"
    ref = row["ref"] if row["ref"] is not None else row["source_ref"]
    alt = alts[0]
    if ref is None:
        return "no_ref"
    if len(ref) == 1 and len(alt) == 1:
        return "fixed"
    if alt == "*" or _symbolic(alt) or not set(alt.upper()) <= _SEQUENCE:
        return "fixed"
    if row["ref_check"] != RefCheck.MATCH.value or row["ref"] is None:
        return str(row["ref_check"]) if row["ref_check"] is not None else "unchecked"
    if row["gt"] is None and row["gt_alleles"] is not None:
        return "unindexed"
    if row["gt"] is not None and any(index is not None and index > 1 for index in row["gt"]):
        return "index_past_alts"
    return None


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


def _moved(row: dict, pos: int, ref: str, alt: str, counts: Counter[str]) -> dict:
    """The row at its new position and representation, positional columns cleared."""
    counts["aligned_moved" if pos != row["pos"] else "aligned_trimmed"] += 1
    moved = dict(row)
    moved.update({name: None for name in _POSITIONAL})
    gt = row["gt"]
    old = [row["ref"], row["alts"][0]]
    alleles = row["source_alleles"]
    vectors_follow = alleles is not None and [allele.upper() for allele in alleles] == [allele.upper() for allele in old]
    if alleles is not None and not vectors_follow:
        counts["vectors_unindexed"] += 1
    moved.update(
        pos=pos,
        end_pos=pos + len(ref) - 1,
        ref=ref,
        alts=[alt],
        gt_alleles=None if gt is None else [None if index is None else (ref, alt)[index] for index in gt],
        source_alleles=[ref, alt] if vectors_follow else None,
    )
    if not vectors_follow:
        moved.update(gl=None, gp=None, ad=None)
    return moved


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
