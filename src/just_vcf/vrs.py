"""VRS Allele ids (RM271, charter invariant 4).

`mint_vrs_ids(loaded, reference)` fills `calls.vrs_allele_ids` for a facet whose REF `resolve_ref`
has read: one id per allele, parallel to `[ref] + alts` (gnomAD's convention), each the VRS 2.x id of
the fully justified allele on the refget sequence the facet's REF was read from. The collection is
`facets.ref_seqcol_digest`; each contig's accession (`SQ.…`) is the store's.

An id that cannot be minted conformantly is withheld, never approximated. Why one is missing is
computed by `vrs_withheld`, over the frame; `vrs_errors` holds only what minting itself reported.

How an allele becomes an id:

- `_justify` applies the VRS 2.x normalization, step by step as vrs-python numbers it, and returns
  the allele's location and state. gtars rolls an indel over its region of ambiguity
  (`normalize_allele`); the choice of state (a `LiteralSequenceExpression`, or a
  `ReferenceLengthExpression` for a deletion, a repeat insertion and the reference allele) is made
  here, because gtars 0.10 always returns a literal sequence. The reference allele and single-base
  substitutions, most of a genome, need no bases and take the same rule as Polars expressions.
- `_location_json` and `_allele_json` serialize the objects as the spec's canonical JSON (keys
  sorted, no whitespace, inherent properties only: a `ReferenceLengthExpression`'s `sequence` is not
  digested), and gtars' sha512t24u digests them. That is the only place ids are serialized.

The test suite compares every kind of allele with vrs-python's ids on the same sequence.

Coordinates: `pos` is 1-based; VRS locations are 0-based interbase. `reference._span` is the one place
that converts between them (charter 13), for the store and for VRS alike.
"""

import logging
import re
from collections import Counter
from pathlib import Path

import polars as pl
from gtars.refget import sha512t24u_digest
from gtars.vrs import normalize_allele

from just_vcf import cache
from just_vcf.callset import NOT_A_VARIANT
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.reference import Contig, ReferenceMismatchError, ReferenceStore, _span, _store_for
from just_vcf.schema import CALLS, SCHEMA_VERSION, RefCheck, empty

logger = logging.getLogger(__name__)

MINT_VERSION = 3
"""Part of the spill key: bumped whenever the rule changes what `mint_vrs_ids` writes."""

MINT_ERRORS: frozenset[str] = frozenset({"no_contig", "refused"})
"""What minting itself reports in `vrs_errors`: the collection has no single sequence for the
contig (`no_contig`), or gtars refused the allele (`refused`, its message counted in the warnings)."""

WITHHELD_REASONS: frozenset[str] = frozenset(
    {
        # What the allele is, read from its string.
        "sv_no_vrs_model",
        "symbolic",
        "unspecified",
        "spanning_deletion",
        "array_code",
        "partial_call",
        "not_sequence",
        # Pipeline state, read from the row.
        "telomere",
        "ref_unresolved",
        "ref_unchecked",
        "ref_mismatch",
        # Nothing forbids minting, and it has not run.
        "not_run",
    }
) | MINT_ERRORS
"""Every reason `vrs_withheld` can give. An open vocabulary: a new kind of facet may add one."""

_SEQUENCE = r"^[ACGTN]*$"
"""A sequence allele. The empty allele (`schema.EMPTY_ALLELE`) is one: VRS's zero-length literal, the
ALT of a deletion written without its anchor base, minted like the anchored spelling of the same event."""
_ID_LIST = pl.List(pl.String)
_READ = ["chrom", "pos", "ref", "source_ref", "ref_check", "alts", "gt_alleles"]
"""The `calls` columns minting reads."""


def is_withheld_reason(value: str) -> bool:
    """Whether `value` is one of `WITHHELD_REASONS`."""
    return value in WITHHELD_REASONS


def _alleles() -> pl.Expr:
    """`[ref] + alts` per row; a null `alts` (an array row whose REF was not read) counts as none."""
    return pl.concat_list(pl.col("ref"), pl.col("alts").fill_null(pl.lit([], dtype=_ID_LIST)))


def _allele_kind(allele: pl.Expr) -> pl.Expr:
    """Why an allele cannot have an id, read from its string alone; null for a sequence allele."""
    upper = allele.str.to_uppercase()
    symbolic = upper.str.starts_with("<") | upper.str.contains(r"[\[\]]") | upper.str.starts_with(".") | upper.str.ends_with(".")
    return (
        pl.when(upper.is_null())
        .then(None)
        .when(upper == "<*>")
        .then(pl.lit("unspecified"))
        .when(upper == "*")
        .then(pl.lit("spanning_deletion"))
        .when(symbolic & upper.is_in(list(NOT_A_VARIANT)))
        .then(pl.lit("symbolic"))
        .when(symbolic)
        .then(pl.lit("sv_no_vrs_model"))
        .when(upper.str.contains("?", literal=True))
        .then(pl.lit("partial_call"))
        .when(upper.is_in(["I", "D"]) & pl.col("source_ref").is_null())
        .then(pl.lit("array_code"))
        .when(~upper.str.contains(_SEQUENCE))
        .then(pl.lit("not_sequence"))
        .otherwise(None)
    )


def _row_reason() -> pl.Expr:
    """Why no allele of the row can have an id, read from the row's other columns; null when none."""
    array_code = pl.col("source_ref").is_null() & pl.col("gt_alleles").list.eval(pl.element().is_in(["I", "D"])).list.any().fill_null(False)
    return (
        pl.when(array_code & pl.col("ref").is_null())
        .then(pl.lit("array_code"))
        .when(pl.col("pos") == 0)
        .then(pl.lit("telomere"))
        .when(pl.col("ref").is_null() & pl.col("source_ref").is_null())
        .then(pl.lit("ref_unresolved"))
        .when(pl.col("ref").is_null())
        .then(pl.lit("ref_unchecked"))
        .when(pl.col("ref_check") == RefCheck.MISMATCH.value)
        .then(pl.lit("ref_mismatch"))
        .otherwise(None)
    )


def _slots(calls: pl.LazyFrame) -> pl.LazyFrame:
    """One row per allele slot of `[ref] + alts`, from a frame with a `_row` index: `_slot`,
    `_allele`, and why nothing could mint it before minting runs (`_static`, null when it may try)."""
    return (
        calls.with_columns(_alleles().alias("_allele"), _row_reason().alias("_row_reason"))
        .with_columns(pl.int_ranges(0, pl.col("_allele").list.len(), dtype=pl.UInt32).alias("_slot"))
        .explode("_allele", "_slot", empty_as_null=False)  # never empty: the REF slot is always there
        .with_columns(pl.coalesce(_allele_kind(pl.col("_allele")), pl.col("_row_reason")).alias("_static"))
    )


def vrs_withheld(calls: pl.LazyFrame) -> pl.LazyFrame:
    """`calls` with a `vrs_withheld` column: why each allele's id is missing, parallel to `[ref] + alts`.

    Null where the id exists. The first reason that applies wins, in this order: what the allele is
    (`sv_no_vrs_model`, `symbolic`, `unspecified`, `spanning_deletion`, `array_code`, `partial_call`,
    `not_sequence`),
    the row's state (`telomere`, `ref_unresolved`, `ref_unchecked`, `ref_mismatch`), what minting
    reported (the `vrs_errors` element), else `not_run`. A computed column, never stored: RM3 fills
    `ref`, and a stored reason would still say "unchecked".
    """
    indexed = calls.with_row_index("_row")
    slots = _slots(indexed).with_columns(
        pl.col("vrs_allele_ids").list.get(pl.col("_slot"), null_on_oob=True).alias("_id"),
        pl.col("vrs_errors").list.get(pl.col("_slot"), null_on_oob=True).alias("_error"),
    )
    reason = (
        pl.when(pl.col("_id").is_not_null())
        .then(None)
        .otherwise(pl.coalesce(pl.col("_static"), pl.col("_error"), pl.lit("not_run")))
    )
    per_row = slots.group_by("_row", maintain_order=True).agg(reason.alias("vrs_withheld"))
    return indexed.join(per_row, on="_row", how="left", maintain_order="left").drop("_row")


def mint_vrs_ids(loaded: Loaded, reference: "ReferenceStore | Path | str | None" = None) -> Loaded:
    """The facet with `vrs_allele_ids` minted, and `vrs_errors` where minting refused an allele.

    `loaded` must have been through `resolve_ref`: its `facets.ref_seqcol_digest` names the
    collection REF was read from, and `reference` (a `ReferenceStore`, a store directory, or None for
    `reference_dir()`) must hold it. Only an allele whose row has `ref` read and not contradicted by
    the source (`ref_check` is not `mismatch`) is minted. The result is the same facet, its rows in
    the same order; the input is not changed. Withheld alleles are counted per reason in the
    warnings (`vrs_withheld`).

    The input is spilled to parquet once, then read one contig at a time, only the columns minting
    reads; the ids are written to their own files and put beside the spilled rows by position. Memory
    holds one contig's alleles and, while its indels are justified, its sequence.

    Raises `ValueError` when REF was never read, and `ReferenceMismatchError` when the store does not
    hold the facet's collection.
    """
    versions = {"SCHEMA_VERSION": SCHEMA_VERSION, "MINT_VERSION": MINT_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    cache.sweep_stale(cache.VRS_SPILL, versions)
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"mint_vrs_ids takes one facet; this Loaded has {facets.height}")
    facet = facets.row(0, named=True)
    digest = facet["ref_seqcol_digest"]
    if digest is None:
        raise ValueError(f"facet {facet['facet_id']}: REF has not been read; run resolve_ref first")
    if isinstance(reference, (Path, str)) and Path(reference).is_file():
        raise ValueError(f"{reference}: mint_vrs_ids reads a reference store, not a FASTA; pass the store resolve_ref used")
    store, _ = _store_for(reference)
    if digest not in store.collections():
        raise ReferenceMismatchError(f"facet {facet['facet_id']}: the reference store does not hold collection {digest}, which REF was read from")
    by_chrom = store.contigs_by_chrom(digest)

    entry = cache_dir() / "frames" / f"{facet['facet_id']}-{SCHEMA_VERSION}-vrs-{digest}-{MINT_VERSION}"
    with cache.spilling(entry, cache.VRS_SPILL, versions, cache.derived_origin(facet, facet["facet_id"])) as spill:
        # The input is spilled once as it is: the result reads its rows from this file, so they are the
        # rows the ids were minted from, in the same order, whatever plan built the input. Its bytes
        # key the result: a run on the same input reuses what an earlier one wrote.
        staged = spill.stage / "calls.parquet"
        loaded.calls.drop("vrs_allele_ids", "vrs_errors").sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
        key = cache.content_key([staged], facet)
        done = spill.find(key)
        if done is None:
            tally = _Tally()
            parts: list[Path] = []
            calls_spill = pl.scan_parquet(staged)
            # Walked run by run of one contig, in row order: a loader groups rows by contig (Decision 11),
            # so a run is a contig, and the ids come out in the input's order with nothing to sort.
            runs = calls_spill.select(pl.col("chrom").rle()).collect().unnest("chrom")
            batch: list[pl.DataFrame] = []
            offset = 0
            for length, chrom in runs.select("len", "value").iter_rows():
                rows = calls_spill.slice(offset, length).select(_READ).collect().with_row_index("_row", offset=offset)
                matches = by_chrom.get(chrom, [])
                batch.append(_mint_rows(rows, store, matches[0] if len(matches) == 1 else None, tally))
                offset += length
                if sum(frame.height for frame in batch) >= _IDS_PER_PART:
                    parts.append(_write_part(batch, spill.stage, len(parts)))
                    batch = []
            if batch:
                parts.append(_write_part(batch, spill.stage, len(parts)))
            done = spill.commit(key, {"calls": [staged], "ids": parts}, tally.warnings(digest))
    # Read from where the generation was committed, never from the staging directory.
    calls_spill = pl.scan_parquet(done.files("calls")[0])
    parts = done.files("ids")

    warnings = done.warnings
    for warning in warnings:
        logger.warning("%s: %s", facet["facet_id"], warning)
    calls = (
        pl.concat([calls_spill, pl.scan_parquet(parts).drop("_row")], how="horizontal", strict=True).select(CALLS.names())
        if parts
        else empty(CALLS)
    )
    return Loaded(
        facets=loaded.facets,
        calls=calls,
        blocks=loaded.blocks,
        links=loaded.links,
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )


_IDS_PER_PART = 1_000_000
"""Rows of ids gathered before a part file is written: few files, and a bounded batch in memory."""


def _write_part(batch: list[pl.DataFrame], spill: Path, index: int) -> Path:
    part = spill / f"ids-{index:05d}.parquet"
    pl.concat(batch).write_parquet(part, row_group_size=SPILL_ROW_GROUP)
    return part


class _Tally:
    """Counts over every allele slot, for the warnings."""

    def __init__(self) -> None:
        self.minted: Counter[str] = Counter()
        self.withheld: Counter[str] = Counter()
        self.refusals: Counter[str] = Counter()

    def warnings(self, digest: str) -> list[str]:
        slots = sum(self.minted.values()) + sum(self.withheld.values())
        lines = [
            f"VRS ids minted on collection {digest} for {sum(self.minted.values())} of {slots} allele(s): "
            f"{self.minted['ref']} REF, {self.minted['alt']} ALT"
        ]
        for reason, count in sorted(self.withheld.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"{count} allele(s) without a VRS id: {_EXPLAIN.get(reason, reason)} ({reason})")
        for message, count in sorted(self.refusals.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"{count} allele(s) refused by gtars: {message}")
        return lines


_EXPLAIN = {
    "sv_no_vrs_model": (
        "a structural variant (a symbolic allele or a breakend): no VRS Allele holds one, and VRS 2's "
        "copy-number and adjacency objects are not minted"
    ),
    "symbolic": "a symbolic allele that names no variant (Complete Genomics' no-call)",
    "unspecified": "the unspecified allele <*>",
    "spanning_deletion": "the spanning deletion *, called at another record",
    "array_code": "a 23andMe I/D code (resolving it needs the chip manifest)",
    "partial_call": "an allele with bases the caller could not call (?)",
    "not_sequence": "an allele with characters other than A, C, G, T and N",
    "telomere": "a telomere record at POS 0",
    "ref_unresolved": "no REF: the source wrote none and the reference was not read",
    "ref_unchecked": "REF not read from the reference (a contig it could not match, or a span past the contig's end)",
    "ref_mismatch": "the source REF disagrees with the reference",
    "no_contig": "the collection has no single sequence for the contig",
    "refused": "gtars refused the allele",
    "not_run": "not minted",
}

_JUSTIFIED = {"_start": pl.Int64, "_end": pl.Int64, "_length": pl.Int64, "_cycle": pl.Int64, "_sequence": pl.String, "_error": pl.String}
"""What `_justify` gives per allele: the location, and the state (a literal `_sequence`, or a
reference `_length` with its repeat `_cycle`), or `_error` when gtars refused it."""


def _mint_rows(rows: pl.DataFrame, store: ReferenceStore, contig: Contig | None, tally: _Tally) -> pl.DataFrame:
    """One contig's `_row`, `vrs_allele_ids` and `vrs_errors`, in `_row` order."""
    ref_length = pl.col("ref").str.len_chars().cast(pl.Int64)
    start, _ = _span(pl.col("pos"), ref_length)
    slots = (
        _slots(rows.lazy())
        .with_columns(start.alias("_start"), pl.col("_allele").str.to_uppercase())  # VCF bases are case-insensitive
        .collect()
    )
    todo = slots.filter(pl.col("_static").is_null()).select("_row", "_slot", "_start", "ref", "_allele")
    # The reference allele (2.a) and single-base substitutions (2.b) need no bases and no trimming.
    identity = pl.col("_allele") == pl.col("ref")
    snv = (ref_length == 1) & (pl.col("_allele").str.len_chars() == 1)
    simple = todo.filter(identity | snv).select(
        "_row",
        "_slot",
        "_start",
        (pl.col("_start") + ref_length).alias("_end"),
        pl.when(identity).then(ref_length).alias("_length"),
        pl.when(identity).then(ref_length).alias("_cycle"),
        pl.when(identity).then(None).otherwise(pl.col("_allele")).alias("_sequence"),
        pl.lit(None, dtype=pl.String).alias("_error"),
    )
    rest = todo.filter(~(identity | snv))
    if contig is None:
        simple = simple.with_columns(pl.lit("no_contig").alias("_error"))
        justified = rest.select("_row", "_slot").with_columns(
            [pl.lit(None, dtype=dtype).alias(name) for name, dtype in _JUSTIFIED.items() if name != "_error"] + [pl.lit("no_contig").alias("_error")]
        )
    else:
        sequence = _ContigSequence(store, contig)
        results = [_justify(*args, sequence, tally) for args in zip(rest["_start"].to_list(), rest["ref"].to_list(), rest["_allele"].to_list())]
        justified = rest.select("_row", "_slot").hstack(pl.DataFrame(results, schema=_JUSTIFIED, orient="row"))
    states = pl.concat([simple, justified.select(simple.columns)])
    ids = _identify(states, f"SQ.{contig.digest}") if contig is not None else states.with_columns(pl.lit(None, dtype=pl.String).alias("_id"))

    slots = slots.join(ids.select("_row", "_slot", "_id", "_error"), on=["_row", "_slot"], how="left").sort("_row", "_slot")
    tally.minted["ref"] += slots.filter(pl.col("_id").is_not_null() & (pl.col("_slot") == 0)).height
    tally.minted["alt"] += slots.filter(pl.col("_id").is_not_null() & (pl.col("_slot") > 0)).height
    reasons = slots.filter(pl.col("_id").is_null()).select(pl.coalesce("_static", "_error", pl.lit("not_run")).alias("r"))["r"]
    tally.withheld.update(dict(reasons.value_counts().iter_rows()))
    # A list with no element set is null as a whole: nothing was minted, or nothing was refused.
    return (
        slots.group_by("_row", maintain_order=True)
        .agg(pl.col("_id").alias("vrs_allele_ids"), pl.col("_error").alias("vrs_errors"))
        .with_columns(
            pl.when(pl.col(name).list.drop_nulls().list.len() > 0).then(pl.col(name)).otherwise(None).cast(_ID_LIST).alias(name)
            for name in ("vrs_allele_ids", "vrs_errors")
        )
    )


def _location_json(accession: str) -> pl.Expr:
    """A SequenceLocation over `_start`/`_end` on `accession`, as canonical JSON."""
    reference = '"sequenceReference":{"refgetAccession":"' + accession + '","type":"SequenceReference"}'
    return pl.concat_str(
        pl.lit('{"end":'), pl.col("_end").cast(pl.String), pl.lit("," + reference + ',"start":'), pl.col("_start").cast(pl.String), pl.lit(',"type":"SequenceLocation"}')
    )


def _allele_json() -> pl.Expr:
    """An Allele at the location digest `_location` with the state in `_sequence` or `_length`/`_cycle`.

    The accession is base64url and a sequence is `[ACGTN]`, so nothing needs JSON escaping.
    """
    state = (
        pl.when(pl.col("_sequence").is_not_null())
        .then(pl.concat_str(pl.lit('{"sequence":"'), pl.col("_sequence"), pl.lit('","type":"LiteralSequenceExpression"}')))
        .otherwise(
            pl.concat_str(
                pl.lit('{"length":'),
                pl.col("_length").cast(pl.String),
                pl.lit(',"repeatSubunitLength":'),
                pl.col("_cycle").cast(pl.String),
                pl.lit(',"type":"ReferenceLengthExpression"}'),
            )
        )
    )
    return pl.concat_str(pl.lit('{"location":"'), pl.col("_location"), pl.lit('","state":'), state, pl.lit(',"type":"Allele"}'))


def _digest(texts: pl.Series) -> pl.Series:
    """gtars' sha512t24u of each string."""
    return pl.Series([sha512t24u_digest(text) for text in texts.to_list()], dtype=pl.String)


def _identify(states: pl.DataFrame, accession: str) -> pl.DataFrame:
    """`states` with `_id`, the VRS Allele id, where the state is known."""
    ok = states.filter(pl.col("_error").is_null())
    spans = ok.select("_start", "_end").unique(maintain_order=True)
    spans = spans.with_columns(_digest(spans.select(_location_json(accession)).to_series()).alias("_location"))
    ok = ok.join(spans, on=["_start", "_end"], how="left", maintain_order="left")
    ok = ok.with_columns(("ga4gh:VA." + _digest(ok.select(_allele_json()).to_series())).alias("_id")).drop("_location")
    refused = states.filter(pl.col("_error").is_not_null()).with_columns(pl.lit(None, dtype=pl.String).alias("_id"))
    return pl.concat([ok, refused])


class _ContigSequence:
    """A contig's bases, read from the store once, on the first allele that needs them (an indel)."""

    def __init__(self, store: ReferenceStore, contig: Contig) -> None:
        self.store = store
        self.contig = contig
        self._bases: str | None = None

    def bases(self) -> str:
        if self._bases is None:
            self._bases = self.store.fetch(self.contig, [(0, self.contig.length)])[0]
        return self._bases


def _justify(start: int, ref: str, alt: str, sequence: _ContigSequence, tally: _Tally) -> tuple:
    """`ref` → `alt` at 0-based `start`, fully justified: a row of `_JUSTIFIED`.

    The VRS 2.x normalization, numbered as vrs-python numbers its steps. `ref` is the reference's own
    bases over the span (`resolve_ref` read them). Steps 2.a and 2.b for single bases are taken in
    `_mint_rows` before this is called; they are repeated here for any length.
    """
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    trimmed_ref, trimmed_alt = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(trimmed_ref), len(trimmed_alt)) and trimmed_ref[-1 - suffix] == trimmed_alt[-1 - suffix]:
        suffix += 1
    trimmed_ref = trimmed_ref[: len(trimmed_ref) - suffix]
    trimmed_alt = trimmed_alt[: len(trimmed_alt) - suffix]

    # 2.a: the reference allele, on its own span, as a ReferenceLengthExpression.
    if not trimmed_ref and not trimmed_alt:
        return start, start + len(ref), len(ref), len(ref), None, None
    # 2.b: a substitution, on the trimmed span; no bases needed.
    if trimmed_ref and trimmed_alt:
        return start + prefix, start + prefix + len(trimmed_ref), None, None, trimmed_alt, None

    # 3-4: an insertion or deletion, rolled over its whole region of ambiguity by gtars.
    seed_length = len(trimmed_ref) or len(trimmed_alt)
    bases = sequence.bases()
    try:
        rolled = normalize_allele(bases, start, ref, alt)
    except ValueError as refusal:
        tally.refusals[re.sub(r"[0-9]+", "#", str(refusal).split(":")[0])] += 1  # one line per kind, not per position
        return None, None, None, None, None, "refused"
    new_start, new_end, extended_alt = int(rolled["start"]), int(rolled["end"]), str(rolled["allele"])
    extended_ref = bases[new_start:new_end]
    # 5.a: an insertion the reference does not repeat.
    if not extended_ref:
        return new_start, new_end, None, None, extended_alt, None
    # 5.b: a deletion.
    if len(extended_alt) < len(extended_ref):
        return new_start, new_end, len(extended_alt), seed_length, None, None
    # 5.c: an insertion inside a repeat: a ReferenceLengthExpression when the inserted bases continue
    # a cycle of the reference, the longest cycle that divides the inserted length first.
    if len(extended_alt) > len(extended_ref):
        for cycle in _factors_descending(seed_length):
            if cycle > len(extended_ref):
                continue
            template = extended_ref[len(extended_ref) - cycle :]
            tail = extended_alt[len(extended_ref) :]
            if all(base == template[i % cycle] for i, base in enumerate(tail)):
                return new_start, new_end, len(extended_alt), cycle, None, None
    # 5.c.3 and 5.e: a literal sequence over the expanded span.
    return new_start, new_end, None, None, extended_alt, None


def _factors_descending(n: int) -> list[int]:
    """Every factor of `n`, largest first."""
    return sorted((f for f in range(1, n + 1) if n % f == 0), reverse=True)
