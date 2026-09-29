"""Declarative description of an input file format.

A `FormatSpec` says how to recognise a format from the head of a decompressed file and what a loader
needs to read it: delimiter, comment prefix, where the column names are, the columns with their dtypes
and alternative spellings, the coordinate convention and the genotype encoding. Recognition and
loading both read the spec; neither hardcodes a vendor.

The shape follows `cgm_format`'s schema modules, with one difference required by just-vcf's rules:
the set of format ids grows additively, so it is a `frozenset[str]` assembled from the registry
(`just_vcf.formats.FORMATS`), not an enum.
"""

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import NotRequired, TypedDict

import polars as pl

_FORMAT_ID = re.compile(r"^[a-z0-9][a-z0-9_]*$")


class Family(str, Enum):
    """What kind of content a format carries. Bound to the schema, so a closed enum."""

    GENOTYPE_TABLE = "genotype_table"
    """Genotypes at predefined sites (arrays, imputation output). An absent site is unknown."""
    VARIANT_CALLS = "variant_calls"
    """Discovery calls in VCF or a VCF-like layout."""
    CALLSET_NATIVE = "callset_native"
    """A vendor's native whole-genome callset that states reference and no-call regions (CGI var)."""
    ANNOTATED_VARIANTS = "annotated_variants"
    """A variant list whose purpose is annotation (clinical reports, GFF, annotator JSON)."""
    METADATA = "metadata"
    """Describes a callset without holding its calls (CGI summary)."""


class Coordinates(str, Enum):
    """The position convention of a format's coordinate columns."""

    ONE_BASED = "1-based"
    ZERO_BASED_HALF_OPEN = "0-based-half-open"


class HeaderMode(str, Enum):
    """Where a delimited format keeps its column names."""

    NONE = "none"
    """No column names; `FormatSpec.columns` gives them in file order."""
    LINE = "line"
    """The first line that is not a comment."""
    COMMENTED = "commented"
    """The last comment line before the data (23andMe `# rsid\\tchromosome…`)."""
    VCF = "vcf"
    """VCF meta lines and a `#CHROM` line; read through polars-bio."""


class ColumnSpec(TypedDict):
    """One column of a delimited format."""

    name: str
    """The canonical spelling, used in the loaded frame."""
    dtype: type[pl.DataType] | pl.DataType
    description: str
    aliases: NotRequired[tuple[str, ...]]
    """Other spellings the column has had across export versions; normalised to `name` on load."""


@dataclass(frozen=True)
class Detection:
    """How to recognise a format from a decoded head.

    A spec matches when every `patterns` regex is found in the head (`re.search`, multiline), no
    `excludes` regex is found, and, when `data_shape` is set, one of the first data lines (lines not
    starting with the spec's comment prefix) matches it in full. Patterns are compiled once.
    """

    patterns: tuple[str, ...] = ()
    excludes: tuple[str, ...] = ()
    data_shape: str | None = None

    def __post_init__(self) -> None:
        if not self.patterns and self.data_shape is None:
            raise ValueError("a Detection needs at least one pattern or a data_shape")
        for regex in (*self.patterns, *self.excludes, *filter(None, [self.data_shape])):
            re.compile(regex)


@dataclass(frozen=True)
class BuildHint:
    """A header statement that names the genome build. `build` is what a match means."""

    pattern: str
    build: str
    """`GRCh36`, `GRCh37` or `GRCh38`."""


@dataclass(frozen=True)
class FormatSpec:
    """One recognisable input format and what a loader needs to read it."""

    id: str
    family: Family
    vendor: str
    description: str
    detection: Detection
    comment_prefix: str | None = "#"
    delimiter: str | None = "\t"
    """None for formats that are not delimited text (JSON) or that polars-bio reads (VCF)."""
    header: HeaderMode = HeaderMode.LINE
    columns: tuple[ColumnSpec, ...] = ()
    coordinates: Coordinates | None = Coordinates.ONE_BASED
    genotype_encoding: str | None = None
    """How a genotype is written, in prose: allele letters, I/D codes, hemizygous single letters."""
    no_call: tuple[str, ...] = ()
    """Genotype spellings that mean no call (`--`, `00`, `./.`)."""
    build_hints: tuple[BuildHint, ...] = ()
    """Header statements naming the build. Never a vendor default: no hint means unknown."""
    notes: tuple[str, ...] = field(default=())
    """Gotchas a loader must handle, learned from real files."""
    example: str | None = None
    """A synthetic head of this format that must detect as it and as no earlier registry entry.
    Real rsids and positions, never a person's genotypes. None only for binary formats."""

    def __post_init__(self) -> None:
        if not _FORMAT_ID.match(self.id):
            raise ValueError(f"format id {self.id!r} must be lowercase snake_case")
        if self.header is HeaderMode.NONE and not self.columns:
            raise ValueError(f"{self.id}: a headerless format must list its columns")
