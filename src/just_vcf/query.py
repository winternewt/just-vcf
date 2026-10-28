"""Questions about a genome answered without its rows (RM283; CLAUDE.md §2, docs/TRUST.md).

`load` gives the frames, and code that may read rows uses them. This module is the helper beside it
for agents and apps: `query(source)` answers counts, summary statistics, the facet's provenance, the
load's warnings and whether the facet carries one given variant, and nothing it returns can hold a row.
Every answer is a frozen pydantic model built from counts, numbers, the schema's enums and bounded
strings; `extra="forbid"` keeps a field from being added on the way out. What never appears:

- the facet's `sample`, `source_path`, `member` (an archive member's file name: participant uploads
  carry real names there), `sample_meta`, `ploidy_statements`, `header_contigs` and
  `derivation_inputs`: header and command lines can carry names and paths. `source_sha256` is left
  out too; `facet_id` already names the load. Any of those strings found inside a text that is
  returned (a warning, a caller line) is replaced by a placeholder, and so is anything shaped like a
  file path;
- a position, an allele or a genotype the caller did not ask about. Counts group by contig, variant
  class, filter set, genotype class, REF check, evidence basis and block kind, never by position, and
  there is no region filter: counts narrowed region by region would find a row by bisection.

A count below some k is not suppressed. This is one person's genome, so a small count identifies no
one the genome does not already name; with no position dimension a count cannot be narrowed onto a
row; and a suppressed cell next to its total is recovered by subtraction, so suppression would only
look like protection while breaking the totals apps add up.

This is a convention for a session, not a sandbox: the `Query` holds the `Loaded` in a private slot,
and Python code can still reach it. What it guarantees is that an agent working through these answers
never reads a row into its context, which is the exit docs/TRUST.md models (the model provider).
"""

import re
from collections.abc import Iterable, Mapping, Sequence
from enum import Enum
from pathlib import Path, PurePath
from typing import Annotated

import polars as pl
from pydantic import BaseModel, ConfigDict, Field, StringConstraints

from just_vcf.callset import absence_is_reference, variant_class_expr
from just_vcf.loading import Loaded, load
from just_vcf.schema import (
    CONTIG_RANK,
    AbsenceMeans,
    AssemblyLabelBasis,
    CallsetType,
    ChipBasis,
    ContentFilter,
    ContentFilterBasis,
    EvidenceBasis,
    RefCheck,
    SplitSignal,
    StrandCheck,
    VariantClass,
    canonical_chrom,
)
from just_vcf.formats.spec import Family
from just_vcf.sniff import Sniffed

MAX_LABEL = 128
"""The longest string a key or provenance field may carry. Contig, filter and caller names are short;
a longer one is refused, not cut."""

MAX_REASON = 400
"""The longest warning reason returned; a longer one is shortened, ending in `…`."""

Label = Annotated[str, StringConstraints(max_length=MAX_LABEL)]
Reason = Annotated[str, StringConstraints(max_length=MAX_REASON)]
Digest = Annotated[str, StringConstraints(max_length=MAX_LABEL, pattern=r"^[A-Za-z0-9._:-]*$")]

VRS_ALLELE_ID = re.compile(r"^ga4gh:VA\.[A-Za-z0-9_-]{32}$")
SEQUENCE_ALLELE = re.compile(r"^[ACGTN]*$")
_PATH_LIKE = re.compile(r"(?<![\w/\\])(?:(?:~|[A-Za-z]:)(?:[\\/]+[\w.@+-]+)+|(?:[\\/]+[\w.@+-]+){2,})[\\/]?")
"""A path: two components or more (`/home/x/y`, a URL's `///home/x`, `\\\\host\\share`), or one after
`~` or a drive letter (`~/x`, `C:\\x`). A file name alone is not one: the facet's own names cover it."""

PROVENANCE_TEXT = ("format_id", "vendor", "caller", "caller_version", "phasing_method", "assembly_label", "relation", "vcf_version", "chip", "split_marker", "reference", "aligner", "aligner_version")
"""Facet columns holding free text from the file's header or its caller (a `##DeepVariant_version`
can be a path), passed through `FacetRedactor` wherever they are shown: here and in the CLI."""
_COUNTED = re.compile(r"^(?:(?P<scope>[^:]{1,60} facet): )?(?P<count>\d+) (?P<rest>.+)$", re.DOTALL)

SAMPLE_PLACEHOLDER = "<sample>"
PATH_PLACEHOLDER = "<path>"
NAME_PLACEHOLDER = "<file name>"


class _Answer(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class Frame(str, Enum):
    """Which frame a count or a statistic reads."""

    CALLS = "calls"
    BLOCKS = "blocks"


class Dimension(str, Enum):
    """What a count may be grouped by. None of them is a position."""

    CHROM = "chrom"
    """The canonical contig (Decision 7), both frames."""
    VARIANT_CLASS = "variant_class"
    """The class of each ALT allele (`calls` only); grouping by it counts ALT alleles, not rows."""
    FILTER = "filter"
    """The row's site filters, sorted and joined with `;` (`calls` only); null = the source said `.`."""
    GENOTYPE_CLASS = "genotype_class"
    """`GenotypeClass` of the row (`calls` only)."""
    REF_CHECK = "ref_check"
    """`calls` only."""
    EVIDENCE_BASIS = "evidence_basis"
    """`calls` only."""
    BLOCK_KIND = "block_kind"
    """`ref` or `no_call` (`blocks` only)."""


class CountUnit(str, Enum):
    """What one unit of a count is."""

    ROWS = "rows"
    ALT_ALLELES = "alt_alleles"
    BLOCKS = "blocks"


class StatField(str, Enum):
    """The numeric columns `stats` summarises."""

    GQ = "gq"
    DP = "dp"
    """`calls.dp`, this sample's depth."""
    QUAL = "qual"
    DP_MIN = "dp_min"
    """`blocks.dp_min`, a block's minimum depth."""


_FRAME_DIMENSIONS = {
    Frame.CALLS: {Dimension.CHROM, Dimension.VARIANT_CLASS, Dimension.FILTER, Dimension.GENOTYPE_CLASS, Dimension.REF_CHECK, Dimension.EVIDENCE_BASIS},
    Frame.BLOCKS: {Dimension.CHROM, Dimension.BLOCK_KIND},
}
_FRAME_FIELDS = {Frame.CALLS: {StatField.GQ, StatField.DP, StatField.QUAL}, Frame.BLOCKS: {StatField.GQ, StatField.DP_MIN}}


class GenotypeClass(str, Enum):
    """A `calls` row's genotype read from `gt` (indices, `0` = REF). just-vcf's own vocabulary.

    In a row with a `split_signal` a `0` can mean another ALT (RM270), so `hom_ref` and `het` there may
    be `het_alt` and `hom_alt`; run `rejoin_split` first where that matters.
    """

    HOM_REF = "hom_ref"
    HET = "het"
    """REF and one ALT."""
    HOM_ALT = "hom_alt"
    HET_ALT = "het_alt"
    """Two or more different ALTs."""
    HAPLOID_REF = "haploid_ref"
    HAPLOID_ALT = "haploid_alt"
    PARTIAL = "partial"
    """Some haplotypes called, some not (`0/.`, `1/.`)."""
    NO_CALL = "no_call"
    UNINDEXED = "unindexed"
    """The alleles are stated (`gt_alleles`) but not indexed against REF: an array row before
    `resolve_ref`, or one it withheld (a mixed-strand facet)."""
    NOT_STATED = "not_stated"
    """The row states no genotype."""


class CountGroup(_Answer):
    """One group of a count: a value per dimension, in the order asked. Null is the frame's null."""

    values: tuple[Label | None, ...]
    count: int = Field(ge=0)
    bases: int | None = Field(default=None, ge=0)
    """For blocks: the bases they cover, `end_pos - pos + 1` summed. Null for calls."""


class CountTable(_Answer):
    """A facet's rows (or ALT alleles) counted by the dimensions asked, groups in a fixed order."""

    facet_id: Digest
    frame: Frame
    unit: CountUnit
    dimensions: tuple[Dimension, ...]
    total: int = Field(ge=0)
    groups: tuple[CountGroup, ...]
    rows_without_alts: int | None = Field(default=None, ge=0)
    """When counting ALT alleles: rows whose `alts` are not stated (an array row before `resolve_ref`),
    which hold no allele to count. Null otherwise."""


class NumericSummary(_Answer):
    """Summary of one numeric column over one group. Every statistic is null when nothing is present."""

    values: tuple[Label | None, ...]
    present: int = Field(ge=0)
    missing: int = Field(ge=0)
    minimum: float | None
    p25: float | None
    median: float | None
    p75: float | None
    maximum: float | None
    mean: float | None


class StatTable(_Answer):
    facet_id: Digest
    frame: Frame
    field: StatField
    dimensions: tuple[Dimension, ...]
    groups: tuple[NumericSummary, ...]


class VariantClassCount(_Answer):
    variant_class: VariantClass
    alleles: int = Field(ge=0)


class FacetProvenance(_Answer):
    """A facet's provenance (docs/SCHEMA.md Decision 5), the fields that name no person and no path."""

    facet_id: Digest
    format_id: Label | None
    family: Family | None
    vendor: Label | None
    caller: Label | None
    caller_version: Label | None
    phasing_method: Label | None
    reference: Label | None
    aligner: Label | None
    aligner_version: Label | None
    assembly_label: Label | None
    assembly_label_basis: AssemblyLabelBasis | None
    seqcol_digest: Digest | None
    sorted_name_length_pairs_digest: Digest | None
    ref_seqcol_digest: Digest | None
    callset_type: CallsetType | None
    absence_means: AbsenceMeans | None
    variant_classes: tuple[VariantClassCount, ...] | None
    content_filter: ContentFilter | None
    content_filter_basis: ContentFilterBasis | None
    derived_from: Digest | None
    relation: Label | None
    vcf_version: Label | None
    chip: Label | None
    chip_basis: ChipBasis | None
    strand_check: StrandCheck | None
    split_marker: Label | None
    contigs: tuple[Label, ...] | None
    """The canonical contig names of `contig_map`; the source spellings are not returned."""
    dropped_fields: tuple[Label, ...] | None


class WarningGroup(_Answer):
    """The load's warnings that share a reason. `count` sums their leading counts where every line
    had one ("12 row(s) … dropped"), and is null where the reason states no count."""

    reason: Reason
    count: int | None = Field(default=None, ge=0)
    lines: int = Field(ge=1)


class PresenceBasis(str, Enum):
    """What a `Presence` answer rests on."""

    STATED = "stated"
    """A `calls` row states the asked allele; its genotype decides. Never a row whose source REF states
    no base (`ref_check = unstated`): an ALT spelled against an `N` is not stated against the site."""
    ROWS_DISAGREE = "rows_disagree"
    """Several rows state the allele and their genotypes disagree (charter 8)."""
    SITE_GENOTYPED = "site_genotyped"
    """A row at the same position and REF states other alleles only. A row whose source REF states no
    base (`ref_check = unstated`) never answers so: the site then reads `overlapping_call`."""
    REF_MISMATCH = "ref_mismatch"
    """A row at the position, but the file's REF there contradicts the reference (`ref_check`)."""
    ALLELES_UNRESOLVED = "alleles_unresolved"
    """A row at the position whose REF is not read (an array before `resolve_ref`): its alleles cannot be
    set against the asked ones without strand and REF."""
    OVERLAPPING_CALL = "overlapping_call"
    """Rows overlap the asked span without starting at it with the same REF (an indel nearby, an MNV), or
    the row at it has a source REF that states no base (`ref_check = unstated`): its ALTs were spelled
    against an `N`, so it states neither the asked allele nor the site's base."""
    REF_BLOCK = "ref_block"
    """A reference block covers the asked span."""
    NO_CALL_BLOCK = "no_call_block"
    """A no-call block reaches into the asked span."""
    ABSENCE = "absence"
    """Nothing covers the position; `absence_means` says what that is worth."""
    NOT_FOUND_BY_ID = "not_found_by_id"
    """Asked by VRS id only, and no row carries the id: without a position there is nothing else to read."""


class Presence(_Answer):
    """Whether one facet carries the asked allele, graded, never the row.

    `carried` is Kleene: `True` (a stated genotype holds it), `False` (evidence of its absence: the
    site genotyped without it, a reference block, or a silence `absence_is_reference` accepts) or
    `None` (unknown). `basis` says which. `False` is only ever answered for an SNV: an indel or an MNV
    the facet does not state may be written there in another representation (shifted, split, merged),
    so its absence is unknown. Rejected: answering an indel's absence when no call lies within some
    window of it, since a repeat can shift the same indel further than any fixed window. No position,
    allele or genotype other than the asked one is returned.
    """

    facet_id: Digest
    carried: bool | None
    basis: PresenceBasis
    copies: int | None = Field(default=None, ge=0)
    """Haplotypes carrying the asked allele, where the genotype states it."""
    ploidy: int | None = Field(default=None, ge=0)
    genotype_class: GenotypeClass | None = None
    """The stating row's genotype class; null with a `split_signal`, where a `0` may be another ALT."""
    evidence_basis: EvidenceBasis | None = None
    gq: int | None = Field(default=None, ge=0)
    dp: int | None = Field(default=None, ge=0)
    filters_passed: bool | None = None
    """`True` for `["PASS"]`, `False` for any other filter, null when the source said `.`."""
    split_signal: SplitSignal | None = None
    absence_means: AbsenceMeans | None = None
    rows: int = Field(default=0, ge=0)
    """How many rows the answer read (rows stating the allele, or rows at the site)."""

    def __bool__(self) -> bool:
        return self.carried is True


def genotype_class_expr(gt: pl.Expr = pl.col("gt"), gt_alleles: pl.Expr = pl.col("gt_alleles")) -> pl.Expr:
    """`GenotypeClass` of each row, as a string, from `gt` (and `gt_alleles` where `gt` is null)."""
    length = gt.list.len()
    missing = gt.list.eval(pl.element().is_null()).list.sum()
    zeros = gt.list.eval(pl.element() == 0).list.sum()
    distinct_alts = gt.list.eval(pl.element().drop_nulls().filter(pl.element() != 0)).list.n_unique()
    return (
        pl.when(gt.is_null() & gt_alleles.is_null())
        .then(pl.lit(GenotypeClass.NOT_STATED.value))
        .when(gt.is_null())
        .then(pl.lit(GenotypeClass.UNINDEXED.value))
        .when(length == 0)
        .then(pl.lit(GenotypeClass.NOT_STATED.value))
        .when(missing == length)
        .then(pl.lit(GenotypeClass.NO_CALL.value))
        .when(missing > 0)
        .then(pl.lit(GenotypeClass.PARTIAL.value))
        .when(length == 1)
        .then(pl.when(zeros == 1).then(pl.lit(GenotypeClass.HAPLOID_REF.value)).otherwise(pl.lit(GenotypeClass.HAPLOID_ALT.value)))
        .when(zeros == length)
        .then(pl.lit(GenotypeClass.HOM_REF.value))
        .when(distinct_alts >= 2)
        .then(pl.lit(GenotypeClass.HET_ALT.value))
        .when(zeros > 0)
        .then(pl.lit(GenotypeClass.HET.value))
        .otherwise(pl.lit(GenotypeClass.HOM_ALT.value))
    )


def _genotype_class(gt: list[int | None] | None, gt_alleles: list[str | None] | None) -> GenotypeClass:
    """`genotype_class_expr` for one row, in Python (the presence answer reads a handful of rows)."""
    value = pl.DataFrame({"gt": [gt], "gt_alleles": [gt_alleles]}, schema={"gt": pl.List(pl.UInt8), "gt_alleles": pl.List(pl.String)})
    return GenotypeClass(value.select(genotype_class_expr()).item())


def _dimension_expr(dimension: Dimension) -> pl.Expr:
    if dimension is Dimension.CHROM:
        return pl.col("chrom")
    if dimension is Dimension.VARIANT_CLASS:
        return variant_class_expr(pl.coalesce("ref", "source_ref"), pl.col("alts"))
    if dimension is Dimension.FILTER:
        return pl.col("filters").list.sort().list.join(";")
    if dimension is Dimension.GENOTYPE_CLASS:
        return genotype_class_expr()
    if dimension is Dimension.BLOCK_KIND:
        return pl.col("kind").cast(pl.String)
    return pl.col(dimension.value).cast(pl.String)


def _sort_key(values: tuple[str | None, ...], dimensions: tuple[Dimension, ...]) -> tuple:
    """Groups in a fixed order: contigs by rank, enum members in declaration order, others by name."""
    order = {
        Dimension.VARIANT_CLASS: [m.value for m in VariantClass],
        Dimension.GENOTYPE_CLASS: [m.value for m in GenotypeClass],
        Dimension.REF_CHECK: [m.value for m in RefCheck],
        Dimension.EVIDENCE_BASIS: [m.value for m in EvidenceBasis],
        Dimension.BLOCK_KIND: ["ref", "no_call"],
    }
    key = []
    for value, dimension in zip(values, dimensions, strict=True):
        if value is None:
            key.append((2, 0, ""))
        elif dimension is Dimension.CHROM:
            key.append((0, CONTIG_RANK.get(value, len(CONTIG_RANK)), value))
        elif dimension in order and value in order[dimension]:
            key.append((0, order[dimension].index(value), value))
        else:
            key.append((1, 0, value))
    return tuple(key)


def redact_paths(text: str) -> str:
    """`text` with anything shaped like a path replaced by `PATH_PLACEHOLDER`."""
    return _PATH_LIKE.sub(PATH_PLACEHOLDER, text)


class FacetRedactor:
    """Replaces the facets' identifying strings, and anything shaped like a path, in returned text.

    The one redaction of header and library text: `Query` uses it, and so does the CLI, before it
    replaces the labels its own run knows (the input's path, `--sample`)."""

    def __init__(self, facets: Iterable[Mapping[str, object]]) -> None:
        tokens: dict[str, str] = {}
        for facet in facets:
            sample = facet.get("sample")
            if isinstance(sample, str):
                tokens[sample] = SAMPLE_PLACEHOLDER
            for column in ("source_path", "member"):
                value = facet.get(column)
                if isinstance(value, str) and value:
                    tokens[value] = PATH_PLACEHOLDER
                    tokens.setdefault(PurePath(value).name, NAME_PLACEHOLDER)
        # Longest first, so a path is replaced whole before its file name. A token shorter than three
        # characters ("S1") names no one and would eat ordinary words. The file name's stem is not a
        # token: stems are words ("calls", "genome") the warnings use, and the loaders quote no names.
        ordered = sorted(((t, r) for t, r in tokens.items() if len(t) >= 3), key=lambda item: -len(item[0]))
        self._tokens = [(re.compile(rf"(?<!\w){re.escape(t)}(?!\w)"), r) for t, r in ordered]

    def __call__(self, text: str) -> str:
        for pattern, replacement in self._tokens:
            text = pattern.sub(replacement, text)
        return redact_paths(text)

    def optional(self, text: object) -> str | None:
        return None if text is None else self(str(text))


class Query:
    """Aggregate answers about one load. Make it with `query(...)`; every method returns models only."""

    __slots__ = ("_loaded", "_facets", "_redact")

    def __init__(self, loaded: Loaded) -> None:
        self._loaded = loaded
        self._facets: list[dict] = loaded.facets.collect().to_dicts()
        self._redact = FacetRedactor(self._facets)

    def __repr__(self) -> str:
        return f"Query(facets={tuple(f['facet_id'] for f in self._facets)!r})"

    @property
    def facet_ids(self) -> tuple[str, ...]:
        return tuple(facet["facet_id"] for facet in self._facets)

    def facets(self) -> tuple[FacetProvenance, ...]:
        """Each facet's provenance, without the sample, the paths and the header's own lines."""
        out = []
        for facet in self._facets:
            text = {
                name: self._redact.optional(facet[name])
                for name in PROVENANCE_TEXT
            }
            classes = facet["variant_classes"]
            contig_map = facet["contig_map"]
            dropped = facet["dropped_fields"]
            out.append(
                FacetProvenance(
                    facet_id=facet["facet_id"],
                    family=facet["family"],
                    assembly_label_basis=facet["assembly_label_basis"],
                    seqcol_digest=facet["seqcol_digest"],
                    sorted_name_length_pairs_digest=facet["sorted_name_length_pairs_digest"],
                    ref_seqcol_digest=facet["ref_seqcol_digest"],
                    callset_type=facet["callset_type"],
                    absence_means=facet["absence_means"],
                    variant_classes=None if classes is None else tuple(VariantClassCount(**entry) for entry in classes),
                    content_filter=facet["content_filter"],
                    content_filter_basis=facet["content_filter_basis"],
                    derived_from=facet["derived_from"],
                    chip_basis=facet["chip_basis"],
                    strand_check=facet["strand_check"],
                    contigs=None if contig_map is None else tuple(sorted({e["chrom"] for e in contig_map}, key=lambda c: (CONTIG_RANK.get(c, len(CONTIG_RANK)), c))),
                    dropped_fields=None if dropped is None else tuple(self._redact(name) for name in dropped),
                    **text,
                )
            )
        return tuple(out)

    def counts(self, by: Sequence[Dimension | str] = (Dimension.CHROM,), frame: Frame | str = Frame.CALLS) -> tuple[CountTable, ...]:
        """Rows of `frame` per facet, grouped by `by` (any number of dimensions, none a position).

        Grouping by `variant_class` counts ALT alleles, as `facets.variant_classes` does; the null class
        holds ALTs that are no variant (`<*>`, `*`, letters outside `ACGTN`). Blocks also report the
        bases they cover.
        """
        frame = Frame(frame)
        dimensions = self._dimensions(by, frame)
        source = self._frame(frame)
        if Dimension.VARIANT_CLASS in dimensions:
            unit = CountUnit.ALT_ALLELES
            unstated = source.filter(pl.col("alts").is_null()).group_by("facet_id").agg(pl.len().alias("n"))
            source = source.filter(pl.col("alts").is_not_null()).explode("alts", empty_as_null=False)
        else:
            unit = CountUnit.BLOCKS if frame is Frame.BLOCKS else CountUnit.ROWS
            unstated = None
        keyed = source.select(pl.col("facet_id"), *[_dimension_expr(d).alias(f"k{i}") for i, d in enumerate(dimensions)], *([(pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64) + 1).alias("bases")] if frame is Frame.BLOCKS else []))
        aggregates = [pl.len().alias("count")] + ([pl.col("bases").sum().alias("bases")] if frame is Frame.BLOCKS else [])
        grouped = keyed.group_by("facet_id", *[f"k{i}" for i in range(len(dimensions))]).agg(aggregates)
        collected = pl.collect_all([grouped] + ([unstated] if unstated is not None else []))
        rows = collected[0].to_dicts()
        no_alts = {r["facet_id"]: r["n"] for r in collected[1].to_dicts()} if unstated is not None else {}
        tables = []
        for facet_id in self.facet_ids:
            groups = []
            for row in rows:
                if row["facet_id"] != facet_id:
                    continue
                values = tuple(self._redact.optional(row[f"k{i}"]) for i in range(len(dimensions)))
                groups.append(CountGroup(values=values, count=row["count"], bases=row.get("bases")))
            groups.sort(key=lambda g: _sort_key(g.values, dimensions))
            tables.append(
                CountTable(
                    facet_id=facet_id,
                    frame=frame,
                    unit=unit,
                    dimensions=dimensions,
                    total=sum(g.count for g in groups),
                    groups=tuple(groups),
                    rows_without_alts=no_alts.get(facet_id, 0) if unit is CountUnit.ALT_ALLELES else None,
                )
            )
        return tuple(tables)

    def stats(self, field: StatField | str, frame: Frame | str | None = None, by: Sequence[Dimension | str] = ()) -> tuple[StatTable, ...]:
        """Present/missing counts, quartiles, extremes and mean of `field` per facet and group.

        `frame` defaults to the one frame that has the field (`dp_min` is a block's, `dp` and `qual` a
        call's); `gq` is in both, and defaults to `calls`.
        """
        field = StatField(field)
        frame = Frame(frame) if frame is not None else (Frame.BLOCKS if field is StatField.DP_MIN else Frame.CALLS)
        if field not in _FRAME_FIELDS[frame]:
            raise ValueError(f"{frame.value} has no {field.value}; it has {', '.join(sorted(f.value for f in _FRAME_FIELDS[frame]))}")
        dimensions = self._dimensions(by, frame)
        if Dimension.VARIANT_CLASS in dimensions:
            raise ValueError("variant_class counts ALT alleles, and a row's statistic belongs to the row: group stats by another dimension")
        value = pl.col(field.value).cast(pl.Float64)
        keys = [f"k{i}" for i in range(len(dimensions))]
        grouped = (
            self._frame(frame)
            .select(pl.col("facet_id"), *[_dimension_expr(d).alias(f"k{i}") for i, d in enumerate(dimensions)], value.alias("v"))
            .group_by("facet_id", *keys)
            .agg(
                pl.col("v").count().alias("present"),
                pl.col("v").null_count().alias("missing"),
                pl.col("v").min().alias("minimum"),
                pl.col("v").quantile(0.25, "linear").alias("p25"),
                pl.col("v").median().alias("median"),
                pl.col("v").quantile(0.75, "linear").alias("p75"),
                pl.col("v").max().alias("maximum"),
                pl.col("v").mean().alias("mean"),
            )
            .collect()
        )
        tables = []
        for facet_id in self.facet_ids:
            groups = []
            for row in grouped.filter(pl.col("facet_id") == facet_id).to_dicts():
                values = tuple(self._redact.optional(row[k]) for k in keys)
                numbers = {name: row[name] for name in ("minimum", "p25", "median", "p75", "maximum", "mean")}
                groups.append(NumericSummary(values=values, present=row["present"], missing=row["missing"], **numbers))
            groups.sort(key=lambda g: _sort_key(g.values, dimensions))
            tables.append(StatTable(facet_id=facet_id, frame=frame, field=field, dimensions=dimensions, groups=tuple(groups)))
        return tuple(tables)

    def warnings(self) -> tuple[WarningGroup, ...]:
        """The load's warnings, one group per reason, in the order the reasons first appear."""
        groups: dict[str, list[int | None]] = {}
        for line in self._loaded.warnings:
            text = self._redact(line)
            match = _COUNTED.match(text)
            if match:
                reason = (f"{match['scope']}: " if match["scope"] else "") + match["rest"]
                count: int | None = int(match["count"])
            else:
                reason, count = text, None
            if len(reason) > MAX_REASON:
                reason = reason[: MAX_REASON - 1] + "…"
            groups.setdefault(reason, []).append(count)
        out = []
        for reason, counts in groups.items():
            total = None if any(c is None for c in counts) else sum(c for c in counts if c is not None)
            out.append(WarningGroup(reason=reason, count=total, lines=len(counts)))
        return tuple(out)

    def carries(
        self,
        *,
        vrs_id: str | None = None,
        chrom: str | None = None,
        pos: int | None = None,
        ref: str | None = None,
        alt: str | None = None,
    ) -> tuple[Presence, ...]:
        """Whether each facet carries one ALT allele, asked by VRS Allele id, by `chrom`/`pos`/`ref`/`alt`
        (1-based VCF coordinates, sequence alleles), or both.

        Alleles are matched as the facet spells them: ask a facet `normalize` produced, with the allele
        left-aligned, or ask by VRS id after `mint_vrs_ids`. An array facet answers once `resolve_ref`
        has read its REF; before that its rows are `alleles_unresolved`.
        """
        coordinates = (chrom, pos, ref, alt)
        if vrs_id is None and all(v is None for v in coordinates):
            raise ValueError("ask by vrs_id, by chrom/pos/ref/alt, or both")
        if vrs_id is not None and not VRS_ALLELE_ID.match(vrs_id):
            raise ValueError("vrs_id is not a VRS Allele id (ga4gh:VA.<32 characters>)")
        by_position = all(v is not None for v in coordinates)
        if not by_position and any(v is not None for v in coordinates):
            raise ValueError("chrom, pos, ref and alt go together")
        if by_position:
            ref, alt = str(ref).upper(), str(alt).upper()
            if not (SEQUENCE_ALLELE.match(ref) and SEQUENCE_ALLELE.match(alt)) or ref == alt or int(pos) < 1:  # type: ignore[arg-type]
                raise ValueError("ref and alt must be different sequence alleles (ACGTN) and pos a 1-based position")
            chrom, pos = canonical_chrom(str(chrom)), int(pos)  # type: ignore[arg-type]
        span_end = (int(pos) + max(len(ref), 1) - 1) if by_position else None  # type: ignore[arg-type]
        calls_filter = pl.lit(False)
        if vrs_id is not None:
            calls_filter = calls_filter | pl.col("vrs_allele_ids").list.slice(1).list.contains(vrs_id)
        if by_position:
            calls_filter = calls_filter | ((pl.col("chrom") == chrom) & (pl.col("pos") <= span_end) & (pl.coalesce("end_pos", "pos") >= pos))
        columns = ["facet_id", "chrom", "pos", "ref", "source_ref", "ref_check", "alts", "vrs_allele_ids", "gt", "gt_alleles", "gq", "dp", "filters", "split_signal", "evidence_basis"]
        frames = [self._loaded.calls.filter(calls_filter).select(columns)]
        if by_position:
            frames.append(self._loaded.blocks.filter((pl.col("chrom") == chrom) & (pl.col("pos") <= span_end) & (pl.col("end_pos") >= pos)).select("facet_id", "pos", "end_pos", "kind", "gq"))
        collected = pl.collect_all(frames)
        calls = collected[0].to_dicts()
        blocks = collected[1].to_dicts() if by_position else []
        asked = _Asked(vrs_id, chrom, pos, ref, alt, span_end) if by_position else _Asked(vrs_id, None, None, None, None, None)
        return tuple(
            _presence(facet, [r for r in calls if r["facet_id"] == facet["facet_id"]], [b for b in blocks if b["facet_id"] == facet["facet_id"]], asked)
            for facet in self._facets
        )

    def _frame(self, frame: Frame) -> pl.LazyFrame:
        return self._loaded.calls if frame is Frame.CALLS else self._loaded.blocks

    @staticmethod
    def _dimensions(by: Sequence[Dimension | str], frame: Frame) -> tuple[Dimension, ...]:
        if isinstance(by, (str, Dimension)):
            by = (by,)
        dimensions = tuple(Dimension(d) for d in by)
        wrong = [d.value for d in dimensions if d not in _FRAME_DIMENSIONS[frame]]
        if wrong:
            raise ValueError(f"{frame.value} cannot be grouped by {', '.join(wrong)}")
        if len(set(dimensions)) != len(dimensions):
            raise ValueError("a dimension is named twice")
        return dimensions


class _Asked:
    __slots__ = ("vrs_id", "chrom", "pos", "ref", "alt", "span_end")

    def __init__(self, vrs_id: str | None, chrom: str | None, pos: int | None, ref: str | None, alt: str | None, span_end: int | None) -> None:
        self.vrs_id, self.chrom, self.pos, self.ref, self.alt, self.span_end = vrs_id, chrom, pos, ref, alt, span_end

    @property
    def by_position(self) -> bool:
        return self.pos is not None

    @property
    def snv(self) -> bool:
        return self.ref is not None and self.alt is not None and len(self.ref) == 1 and len(self.alt) == 1


def _row_ref(row: Mapping) -> str | None:
    value = row["ref"] if row["ref"] is not None else row["source_ref"]
    return None if value is None else value.upper()


def _allele_index(row: Mapping, asked: _Asked) -> int | None:
    """The asked allele's index in the row's `[REF] + alts`, when the row states it.

    A row whose source REF states no base (`ref_check = unstated`, F23) states none of its sequence ALTs:
    each was spelled against the `N`, not against `ref` (`schema.unstated_sequence_alt`). Read beside
    the reference's base, `N` → `G` with GT `0/0` would answer `hom_ref` for a base the source never
    named. The asked allele is always a sequence, so such a row states nothing asked here.
    """
    unstated = row["ref_check"] == RefCheck.UNSTATED.value
    if asked.vrs_id is not None and row["vrs_allele_ids"]:
        for index, value in enumerate(row["vrs_allele_ids"][1:], start=1):
            if value == asked.vrs_id and not (unstated and SEQUENCE_ALLELE.match(row["alts"][index - 1].upper())):
                return index
    stated_ref = row["ref_check"] not in (RefCheck.MISMATCH.value, RefCheck.UNSTATED.value)
    if asked.by_position and row["pos"] == asked.pos and row["alts"] and _row_ref(row) == asked.ref and stated_ref:
        for index, value in enumerate(row["alts"], start=1):
            if value.upper() == asked.alt:
                return index
    return None


def _filters_passed(filters: list[str] | None) -> bool | None:
    if filters is None:
        return None
    return filters == ["PASS"]


def _presence(facet: Mapping, calls: list[dict], blocks: list[dict], asked: _Asked) -> Presence:
    facet_id = facet["facet_id"]
    absence_means = facet["absence_means"]
    stating = [(row, index) for row in calls if (index := _allele_index(row, asked)) is not None]
    if stating:
        answers = []
        for row, index in stating:
            gt = row["gt"]
            called = [g for g in gt or [] if g is not None]
            copies = sum(1 for g in called if g == index) if gt else None
            if not gt:
                carried = None
            elif copies:
                carried = True
            else:
                carried = False if len(called) == len(gt) else None
            answers.append((row, carried, copies))
        if len({a[1] for a in answers}) > 1 or len({a[2] for a in answers}) > 1:
            return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.ROWS_DISAGREE, rows=len(answers), absence_means=absence_means)
        row, carried, copies = answers[0]
        split = row["split_signal"]
        return Presence(
            facet_id=facet_id,
            carried=carried,
            basis=PresenceBasis.STATED,
            copies=copies,
            ploidy=len(row["gt"]) if row["gt"] else None,
            genotype_class=None if split is not None else _genotype_class(row["gt"], row["gt_alleles"]),
            evidence_basis=row["evidence_basis"],
            gq=row["gq"] if len(answers) == 1 else None,
            dp=row["dp"] if len(answers) == 1 else None,
            filters_passed=_filters_passed(row["filters"]) if len(answers) == 1 else None,
            split_signal=split,
            absence_means=absence_means,
            rows=len(answers),
        )
    if not asked.by_position:
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.NOT_FOUND_BY_ID, absence_means=absence_means)
    at_site = [row for row in calls if row["chrom"] == asked.chrom and row["pos"] == asked.pos]
    if any(_row_ref(row) is None for row in at_site):
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.ALLELES_UNRESOLVED, rows=len(at_site), absence_means=absence_means)
    if any(row["ref_check"] == RefCheck.MISMATCH.value for row in at_site):
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.REF_MISMATCH, rows=len(at_site), absence_means=absence_means)
    # A source REF of only N (`unstated`, F23) says nothing about which base the site holds, and its
    # sequence ALTs were spelled against that N (above): its row answers neither for an ALT nor for the
    # site, so the asked allele reads `overlapping_call`, carried unknown.
    same_ref = [row for row in at_site if _row_ref(row) == asked.ref and row["ref_check"] != RefCheck.UNSTATED.value]
    if same_ref and len(same_ref) == len(calls):
        # Every row the span reaches is at this site with this REF, and none states the allele.
        fully_called = all(row["gt"] and all(g is not None for g in row["gt"]) and row["split_signal"] is None for row in same_ref)
        return Presence(
            facet_id=facet_id,
            carried=False if (fully_called and asked.snv) else None,
            basis=PresenceBasis.SITE_GENOTYPED,
            evidence_basis=same_ref[0]["evidence_basis"] if len(same_ref) == 1 else None,
            rows=len(same_ref),
            absence_means=absence_means,
        )
    if calls:
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.OVERLAPPING_CALL, rows=len(calls), absence_means=absence_means)
    if any(block["kind"] == "no_call" for block in blocks):
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.NO_CALL_BLOCK, absence_means=absence_means)
    covering = [b for b in blocks if b["kind"] == "ref" and b["pos"] <= asked.pos and b["end_pos"] >= asked.span_end]
    if covering:
        return Presence(
            facet_id=facet_id,
            carried=False if asked.snv else None,
            basis=PresenceBasis.REF_BLOCK,
            evidence_basis=EvidenceBasis.REF_BLOCK,
            gq=covering[0]["gq"] if len(covering) == 1 else None,
            absence_means=absence_means,
        )
    if blocks:
        # Reference blocks reach into the span without covering all of it; the rest is silence.
        return Presence(facet_id=facet_id, carried=None, basis=PresenceBasis.REF_BLOCK, evidence_basis=EvidenceBasis.REF_BLOCK, absence_means=absence_means)
    reference = absence_is_reference(facet, VariantClass.SNV) if asked.snv else None
    return Presence(
        facet_id=facet_id,
        carried=False if reference is True else None,
        basis=PresenceBasis.ABSENCE,
        evidence_basis=EvidenceBasis.ABSENCE if reference is True else None,
        absence_means=absence_means,
    )


def query(source: Loaded | Path | str | Sniffed, sample: str | None = None) -> Query:
    """Aggregate-only answers about a genome: `source` is a `Loaded`, or a file `load` reads here.

    Passing the file keeps the `Loaded` out of the caller's hands altogether. `sample` picks one sample
    of a multi-sample file, as in `load`.
    """
    return Query(source if isinstance(source, Loaded) else load(source, sample))
