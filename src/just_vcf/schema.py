"""The frames every loader emits: `facets`, `calls` and `blocks` (docs/SCHEMA.md).

This module is the contract in code: column names, dtypes, the closed vocabularies, contig naming
(SCHEMA.md Decision 7) and the canonical row order (Decision 11). It does no I/O. A test compares the
column lists here with the tables in docs/SCHEMA.md, so the two cannot drift apart silently.

Closed vocabularies are `str` enums with a matching `pl.Enum` dtype. Vocabularies that grow
additively (`format_id`, `relation`, `vrs_errors`) stay `pl.String`, because adding a category to a
`pl.Enum` changes the frame's dtype.
"""

from enum import Enum

import polars as pl

from just_vcf.formats.spec import Family

SCHEMA_VERSION = "0.0.4"
"""The frame schema version. Written into parquet metadata as `just_vcf.schema_version` on export."""

SCHEMA_VERSION_KEY = "just_vcf.schema_version"


class RefCheck(str, Enum):
    """The outcome of comparing `source_ref` with `ref` read from the assembly sequence."""

    NOT_RUN = "not_run"
    MATCH = "match"
    MISMATCH = "mismatch"


class EvidenceBasis(str, Enum):
    """Where a row's evidence came from. Which quantities it has is said by `gl`, `gp` and `dosage`."""

    ASSAY = "assay"
    """The source wrote a record at this position."""
    REF_BLOCK = "ref_block"
    """Materialised from a reference block."""
    ABSENCE = "absence"
    """Materialised from a whole-genome callset's silence (charter 11)."""


class BlockKind(str, Enum):
    """What a row of the `blocks` frame covers."""

    REF = "ref"
    NO_CALL = "no_call"


class PhaseSetKind(str, Enum):
    """Which VCF field a phase set came from. PS and PSL sets never merge, even when spelled alike."""

    PS = "ps"
    PSL = "psl"


class AbsenceMeans(str, Enum):
    """What a position absent from a facet means (SCHEMA.md Decision 4)."""

    UNKNOWN = "unknown"
    WEAK_REFERENCE = "weak_reference"
    BLOCKS = "blocks"


class CallsetType(str, Enum):
    """What kind of callset a facet is (RM4). Null when not known."""

    WGS = "wgs"
    WES = "wes"
    PANEL = "panel"
    ARRAY = "array"
    IMPUTED = "imputed"


class SplitSignal(str, Enum):
    """Why a `calls` row looks like one record of a multi-allelic site written one row per ALT (RM270).

    In such a row a `0` can mean another ALT, not REF. Where more than one fits, the first listed wins.
    """

    MARKER = "marker"
    """The row carries a split marker: INFO `OLD_MULTIALLELIC` (vt, Glow), or a bcftools
    `--old-rec-tag` value with its `USED_ALT_IDX` field."""
    SLOT_PATTERN = "slot_pattern"
    """Two or more non-reference rows at one position, each GT slot non-zero in at most one of them."""
    ALL_ZERO_BESIDE_ALT = "all_zero_beside_alt"
    """A fully called all-`0` row beside a non-reference row at the same position."""


class AssemblyLabelBasis(str, Enum):
    """How a facet's `assembly_label` was obtained."""

    STATED = "stated"
    INFERRED = "inferred"


def enum_dtype(vocabulary: type[Enum]) -> pl.Enum:
    """The `pl.Enum` dtype for a closed vocabulary, categories in declaration order."""
    return pl.Enum([member.value for member in vocabulary])


_CONTIG_MAP = pl.List(pl.Struct({"source": pl.String, "chrom": pl.String}))

FACETS = pl.Schema(
    {
        "facet_id": pl.String,
        "source_sha256": pl.String,
        "source_path": pl.String,
        "member": pl.String,
        "format_id": pl.String,
        "family": enum_dtype(Family),
        "vendor": pl.String,
        "caller": pl.String,
        "caller_version": pl.String,
        "sample": pl.String,
        "assembly_label": pl.String,
        "assembly_label_basis": enum_dtype(AssemblyLabelBasis),
        "seqcol_digest": pl.String,
        "ref_seqcol_digest": pl.String,
        "callset_type": enum_dtype(CallsetType),
        "absence_means": enum_dtype(AbsenceMeans),
        "phasing_method": pl.String,
        "derived_from": pl.String,
        "relation": pl.String,
        "contig_map": _CONTIG_MAP,
        "header_contigs": pl.List(pl.Struct({"name": pl.String, "length": pl.UInt64, "md5": pl.String})),
        "vcf_version": pl.String,
        "sample_meta": pl.List(pl.String),
        "dropped_fields": pl.List(pl.String),
        "split_marker": pl.String,
    }
)
"""One row per facet (SCHEMA.md Decision 5)."""

CALLS = pl.Schema(
    {
        "facet_id": pl.String,
        "chrom": pl.String,
        "pos": pl.UInt32,
        "end_pos": pl.UInt32,
        "ids": pl.List(pl.String),
        "ref": pl.String,
        "source_ref": pl.String,
        "ref_check": enum_dtype(RefCheck),
        "alts": pl.List(pl.String),
        "source_alleles": pl.List(pl.String),
        "vrs_allele_ids": pl.List(pl.String),
        "vrs_errors": pl.List(pl.String),
        "gt_alleles": pl.List(pl.String),
        "gt": pl.List(pl.UInt8),
        "source_gt": pl.String,
        "phased": pl.List(pl.Boolean),
        "phase_sets": pl.List(pl.String),
        "phase_set_kinds": pl.List(enum_dtype(PhaseSetKind)),
        "gl": pl.List(pl.Float32),
        "gp": pl.List(pl.Float32),
        "dosage": pl.Float32,
        "evidence_basis": enum_dtype(EvidenceBasis),
        "gq": pl.UInt16,
        "dp": pl.UInt32,
        "ad": pl.List(pl.UInt32),
        "qual": pl.Float32,
        "filters": pl.List(pl.String),
        "sample_filters": pl.List(pl.String),
        "split_signal": enum_dtype(SplitSignal),
    }
)
"""One row per site per facet (SCHEMA.md Decision 6)."""

BLOCKS = pl.Schema(
    {
        "facet_id": pl.String,
        "chrom": pl.String,
        "pos": pl.UInt32,
        "end_pos": pl.UInt32,
        "kind": enum_dtype(BlockKind),
        "gq": pl.UInt16,
        "dp_min": pl.UInt32,
        "gl": pl.List(pl.Float32),
    }
)
"""Reference and no-call intervals, for formats that have them (SCHEMA.md Decision 4)."""


# Contig names (Decision 7). Every accepted spelling is listed; nothing is stripped by pattern, so
# `chrM` cannot turn into `M` and fall out of a join against `MT`.
_PRIMARY = [str(n) for n in range(1, 23)] + ["X", "Y"]
CONTIG_ALIASES: dict[str, str] = {
    **{name: name for name in _PRIMARY},
    **{f"{prefix}{name}": name for prefix in ("chr", "Chr", "CHR") for name in _PRIMARY},
    **{spelling: "MT" for spelling in ("MT", "M", "chrM", "chrMT", "ChrM", "CHRM")},
    # PLINK array codes: 23 X, 24 Y, 25 the X pseudo-autosomal region, 26 mitochondrion.
    "23": "X",
    "24": "Y",
    "25": "XY",
    "26": "MT",
    "XY": "XY",
}
"""Source spelling → canonical `chrom`. A name not listed keeps its source spelling."""

CONTIG_RANK: dict[str, int] = {name: rank for rank, name in enumerate([*_PRIMARY, "MT", "XY"])}
"""Sort rank of the canonical names; any other contig sorts after them, by name."""


def canonical_chrom(name: str) -> str:
    """The canonical `chrom` for a source contig name."""
    return CONTIG_ALIASES.get(name, name)


def canonical_chrom_expr(column: str) -> pl.Expr:
    """`canonical_chrom` as a Polars expression over a string column."""
    return pl.col(column).replace(CONTIG_ALIASES)


def contig_map(source_names: list[str]) -> list[dict[str, str]]:
    """The facet's `contig_map`: every source name seen, with its canonical name, sorted by source."""
    return [{"source": name, "chrom": canonical_chrom(name)} for name in sorted(set(source_names))]


def canonical_order(frame: pl.LazyFrame) -> pl.LazyFrame:
    """Sort a `calls` or `blocks` frame into the canonical order (Decision 11).

    `facet_id`, then contig rank (the canonical names first, then other names as strings), then
    `pos`, `end_pos`, and for `calls` `source_ref` and `alts`. Nulls sort last. Until a facet's
    sequence collection is known, the rank is the fixed one in `CONTIG_RANK`.
    """
    names = frame.collect_schema().names()
    rank = pl.col("chrom").replace_strict(CONTIG_RANK, default=len(CONTIG_RANK), return_dtype=pl.UInt16)
    keys: list[pl.Expr] = [pl.col("facet_id"), rank, pl.col("chrom"), pl.col("pos"), pl.col("end_pos")]
    if "source_ref" in names:
        keys += [pl.col("source_ref"), pl.col("alts").list.join(",")]
    return frame.sort(keys, nulls_last=True, maintain_order=True)


def canonical_order_per_contig(frame: pl.LazyFrame, contigs: list[str]) -> pl.LazyFrame:
    """`canonical_order` for one facet, one contig at a time, so memory stays within one contig.

    `frame` should be a scan that can skip to a contig cheaply: a parquet file whose rows are grouped
    by contig, where row-group statistics let the filter skip the rest. `contigs` are the canonical
    `chrom` values in the frame. The result equals `canonical_order(frame)` row for row (tested).
    """
    ordered = sorted(set(contigs), key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))
    names = frame.collect_schema().names()
    keys: list[pl.Expr] = [pl.col("pos"), pl.col("end_pos")]
    if "source_ref" in names:
        keys += [pl.col("source_ref"), pl.col("alts").list.join(",")]
    parts = [frame.filter(pl.col("chrom") == name).sort(keys, nulls_last=True, maintain_order=True) for name in ordered]
    return pl.concat(parts, how="vertical", parallel=False) if parts else frame.clear()


def empty(schema: pl.Schema) -> pl.LazyFrame:
    """An empty frame of `schema`."""
    return pl.LazyFrame(schema=schema)


def conform(frame: pl.LazyFrame, schema: pl.Schema) -> pl.LazyFrame:
    """`frame` with exactly `schema`'s columns, in order: missing columns are added as nulls.

    A column present in `frame` with another dtype is cast strictly, so a value that does not fit
    (a position past `u32`) raises instead of turning into a null.
    """
    present = frame.collect_schema()
    extra = [name for name in present.names() if name not in schema]
    if extra:
        raise ValueError(f"columns not in the schema: {', '.join(extra)}")
    return frame.select(
        [
            pl.col(name).cast(dtype, strict=True) if name in present else pl.lit(None, dtype=dtype).alias(name)
            for name, dtype in schema.items()
        ]
    )
