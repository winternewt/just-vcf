"""The frames every loader emits: `facets`, `calls` and `blocks`, and the `links` between facets (docs/SCHEMA.md).

This module is the contract in code: column names, dtypes, the closed vocabularies, contig naming
(SCHEMA.md Decision 7) and the canonical row order (Decision 11). It does no I/O of its own (GRCh37's and
GRCh38's accessions come from the package table `just_vcf.contig_names` reads). A test compares the
column lists here with the tables in docs/SCHEMA.md, so the two cannot drift apart silently.

Closed vocabularies are `str` enums with a matching `pl.Enum` dtype. Vocabularies that grow
additively (`format_id`, `relation`, `vrs_errors`) stay `pl.String`, because adding a category to a
`pl.Enum` changes the frame's dtype.
"""

from enum import Enum

import polars as pl

from just_vcf.contig_names import PRIMARY_ACCESSIONS
from just_vcf.formats.spec import Family

SCHEMA_VERSION = "0.0.21"
"""The frame schema version. Written into parquet metadata as `just_vcf.schema_version` on export."""

SCHEMA_VERSION_KEY = "just_vcf.schema_version"

EMPTY_ALLELE = ""
"""The one spelling of a zero-length allele in `ref`, `source_ref`, `alts`, `source_alleles` and
`gt_alleles`, as VRS writes a zero-length literal: an insertion point's REF (CGI's native files), a
deletion written without its anchor base (CGI's native files, an upload's ALT `-`, F12), and Genos's
`-`. `source_gt` keeps what the file wrote. Such a row is unanchored: `resolve_ref` reads an insertion
point's anchor base and writes it into every allele (RM3), `normalize` anchors a deletion once
`resolve_ref` has read its REF, and `write_vcf` refuses the row until then (VCF has no empty allele).
An anchored insertion point keeps `source_ref = ""`, what the file wrote."""


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
    """Which field a phase set came from. Sets of different kinds never merge, even when spelled alike."""

    PS = "ps"
    PSL = "psl"
    CGI_HAPLINK = "cgi_haplink"
    """Complete Genomics' `hapLink` (native `var` and `masterVarBeta`): an id on one haplotype of a
    locus, shared by the haplotypes of other loci on the same molecule. Like a PSL name it is one
    haplotype, not a GT slot; unlike PSL it is keyed per chromosome. just-vcf's own term; VCF has none."""


class AbsenceMeans(str, Enum):
    """What a position absent from a facet means (SCHEMA.md Decision 4)."""

    UNKNOWN = "unknown"
    WEAK_REFERENCE = "weak_reference"
    BLOCKS = "blocks"
    REFERENCE = "reference"
    """Every position the file omits was called reference by its producer: Complete Genomics' VCF,
    which writes its variants and its no-called spans and leaves out what it called reference. Stronger
    than `weak_reference`, whose silence is only the caller not reporting a variant. A position a block
    covers is not omitted: the block says what it is. just-vcf's own term; GA4GH has none."""


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


class VariantClass(str, Enum):
    """The kind of change one ALT allele makes against its REF."""

    SNV = "snv"
    MNV = "mnv"
    INDEL = "indel"
    SV = "sv"
    """A symbolic ALT (`<DEL>`, `<DUP:TANDEM>`, `<CNV>`) or a breakend. A deletion written out base
    by base is an `indel` whatever its length: the class follows how the file writes the allele."""


class ContentFilter(str, Enum):
    """The one variant class a facet's rows were cut down to (SCHEMA.md Decision 4). Null: none found."""

    SNV_ONLY = "snv_only"
    INDEL_ONLY = "indel_only"
    SV_ONLY = "sv_only"


CONTENT_FILTER_KEEPS: dict[ContentFilter, VariantClass] = {
    ContentFilter.SNV_ONLY: VariantClass.SNV,
    ContentFilter.INDEL_ONLY: VariantClass.INDEL,
    ContentFilter.SV_ONLY: VariantClass.SV,
}
"""The class each content filter keeps; every other class was cut."""


class ContentFilterBasis(str, Enum):
    """What `content_filter` rests on, strongest first."""

    ROWS = "rows"
    """The rows themselves: one class, the others a residue."""
    STATED = "stated"
    """A header line records the filtering command (`bcftools view -v snps`, GATK SelectVariants)."""
    FILE_NAME = "file_name"
    """Only the file name (`.snp.`, `.indel.`): indirect evidence."""


class AssemblyLabelBasis(str, Enum):
    """How a facet's `assembly_label` was obtained."""

    STATED = "stated"
    INFERRED = "inferred"


class ChipBasis(str, Enum):
    """How a facet's `chip` was obtained: written in the file's header, or read off its probe set."""

    STATED = "stated"
    INFERRED = "inferred"


class StrandCheck(str, Enum):
    """What `resolve_ref` found about the strand an array facet's genotypes are written on (F13).

    Null when the check did not run: not an array, REF not read yet, or too few heterozygous rows.
    """

    FORWARD = "forward"
    """Heterozygous rows carry the reference base, as forward-strand genotypes must."""
    MIXED = "mixed"
    """A share of heterozygous rows carry the complement of the reference base and not the base itself:
    part of the probes is written on the minus strand."""


class GlSource(str, Enum):
    """The FORMAT key(s) a VCF facet's `gl` is read from, by the header's declarations (Decision 3).

    GL and PL state one quantity, likelihoods, on two scales, and load equal into `gl`; this says which
    was written. Null: the facet carries no likelihoods (not a VCF, or a header declaring neither key
    with a genotype count). just-vcf's own vocabulary."""

    GL = "gl"
    """GL alone: log10 likelihoods as written (CGAPipeline's phred-unit GL as GL / 10)."""
    PL = "pl"
    """PL alone: phred-scaled likelihoods, read as `-PL/10`."""
    GL_OVER_PL = "gl_over_pl"
    """Both declared: a row's `gl` is its GL, and its PL only where the GL is missing or does not fit
    the genotypes. Which key a given row took is not kept."""


def enum_dtype(vocabulary: type[Enum]) -> pl.Enum:
    """The `pl.Enum` dtype for a closed vocabulary, categories in declaration order."""
    return pl.Enum([member.value for member in vocabulary])


STEP_ROLES: frozenset[str] = frozenset({"align", "call", "filter", "normalize", "merge", "lift", "annotate", "phase", "other"})
"""What a step of `facets.method_steps` did (SCHEMA.md Decision 5). An open vocabulary that grows
additively, so the column is `pl.String`: `other` holds a step the header records whose role is none
of the rest (a format conversion, a gvcftools pass)."""

METHOD_STEP = pl.Struct(
    {"tool": pl.String, "version": pl.String, "role": pl.String, "options": pl.List(pl.String), "line": pl.String}
)
"""One step of the method chain a header records, in the order the steps ran (`just_vcf.source`)."""

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
        "ploidy_statements": pl.List(pl.String),
        "sample": pl.String,
        "assembly_label": pl.String,
        "assembly_label_basis": enum_dtype(AssemblyLabelBasis),
        "seqcol_digest": pl.String,
        "sorted_name_length_pairs_digest": pl.String,
        "ref_seqcol_digest": pl.String,
        "resolve_version": pl.UInt32,
        "resolve_inputs": pl.String,
        "callset_type": enum_dtype(CallsetType),
        "absence_means": enum_dtype(AbsenceMeans),
        "variant_classes": pl.List(pl.Struct({"variant_class": enum_dtype(VariantClass), "alleles": pl.UInt64})),
        "content_filter": enum_dtype(ContentFilter),
        "content_filter_basis": enum_dtype(ContentFilterBasis),
        "phasing_method": pl.String,
        "reference": pl.String,
        "aligner": pl.String,
        "aligner_version": pl.String,
        "method_steps": pl.List(METHOD_STEP),
        "derived_from": pl.String,
        "relation": pl.String,
        "contig_map": _CONTIG_MAP,
        "header_contigs": pl.List(pl.Struct({"name": pl.String, "length": pl.UInt64, "md5": pl.String})),
        "vcf_version": pl.String,
        "sample_meta": pl.List(pl.String),
        "dropped_fields": pl.List(pl.String),
        "split_marker": pl.String,
        "gl_source": enum_dtype(GlSource),
        "derivation_inputs": pl.List(pl.Struct({"role": pl.String, "name": pl.String, "sha256": pl.String})),
        "chip": pl.String,
        "chip_basis": enum_dtype(ChipBasis),
        "strand_check": enum_dtype(StrandCheck),
        "subject": pl.String,
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
        "phase_set_ordinals": pl.List(pl.UInt32),
        "phase_set_qualities": pl.List(pl.UInt16),
        "gl": pl.List(pl.Float32),
        "gp": pl.List(pl.Float32),
        "dosage": pl.Float32,
        "evidence_basis": enum_dtype(EvidenceBasis),
        "gq": pl.UInt16,
        "hq": pl.List(pl.UInt16),
        "dp": pl.UInt32,
        "ad": pl.List(pl.UInt32),
        "qual": pl.Float32,
        "filters": pl.List(pl.String),
        "sample_filters": pl.List(pl.String),
        "split_signal": enum_dtype(SplitSignal),
        "sv_type": pl.String,
        "cipos": pl.List(pl.Int32),
        "ciend": pl.List(pl.Int32),
        "copy_number": pl.Float32,
        "alt_copy_numbers": pl.List(pl.Float32),
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

_DERIVATION_INPUTS = pl.List(pl.Struct({"role": pl.String, "name": pl.String, "sha256": pl.String}))

LINKS = pl.Schema(
    {
        "from_facet_id": pl.String,
        "to_facet_id": pl.String,
        "relation": pl.String,
        "provenance": pl.Struct({"step": pl.String, "step_version": pl.String, "derivation_inputs": _DERIVATION_INPUTS}),
    }
)
"""One row per link between two facets: the facet graph's edges (SCHEMA.md Decision 12, `just_vcf.links`)."""


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

CHM13_ACCESSIONS: dict[str, str] = {
    # GenBank, GCA_009914755.4: chr1 is CP068277.2 counting down to chr22 CP068256.2 and chrX CP068255.2.
    **{f"CP0682{77 - index:02d}.2": name for index, name in enumerate([*_PRIMARY[:22], "X"])},
    "CP086569.2": "Y",
    "CP068254.1": "MT",
    # RefSeq, GCF_009914755.1: NC_060925.1 (chr1) to NC_060948.1 (chrY). RefSeq has no mitochondrion.
    **{f"NC_0609{25 + index}.1": name for index, name in enumerate(_PRIMARY)},
}
"""T2T-CHM13v2.0's sequence accessions → canonical `chrom`, versions verbatim, from NCBI's assembly
report for GCA_009914755.4 (docs/PROVENANCE.md). An accession names one sequence of one assembly, so
mapping it is safe in a way a bare `chr` strip is not. GRCh37's and GRCh38's assembled molecules map
the same way (`just_vcf.contig_names.PRIMARY_ACCESSIONS`, F19); their scaffolds keep their source
name and are matched against a reference store by `contig_names.scaffold_names`."""

CONTIG_ALIASES.update(CHM13_ACCESSIONS)
CONTIG_ALIASES.update(PRIMARY_ACCESSIONS)

CONTIG_RANK: dict[str, int] = {name: rank for rank, name in enumerate([*_PRIMARY, "MT", "XY"])}
"""Sort rank of the canonical names; any other contig sorts after them, by name."""


def canonical_chrom(name: str) -> str:
    """The canonical `chrom` for a source contig name."""
    return CONTIG_ALIASES.get(name, name)


def canonical_chrom_expr(column: str) -> pl.Expr:
    """`canonical_chrom` as a Polars expression over a string column, applied one batch at a time.

    polars 1.44's streaming engine runs `replace` with a mapping as a columnar function it cannot
    split into batches (an in-memory fallback), and holds every other column of the plan beside it
    until that function has the whole column: `load`'s `calls` sink over a 9.65M-call gVCF held 8 GB
    for 2.6 GB of rows (F26). The mapping is elementwise, so each batch is mapped on its own and the
    plan streams; the values are the ones `replace` over the whole column gives."""
    return pl.col(column).map_batches(_canonical_batch, return_dtype=pl.String, is_elementwise=True)


def _canonical_batch(names: pl.Series) -> pl.Series:
    return names.cast(pl.String).replace(CONTIG_ALIASES)


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
