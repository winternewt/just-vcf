"""Load Complete Genomics' native callsets, `var` and `masterVarBeta`, into `calls` and `blocks` rows.

Both files tile the genome in loci with 0-based half-open `begin`/`end` (docs/FORMATS.md). A uniform
locus is a block: `ref` states reference, `no-call` and `no-ref` (an assembly gap) state that nothing
was called. A variant locus becomes one `calls` row. `var` writes it as one or more rows per haplotype
(`allele` 1, 2), whose sequences are joined in `begin` order (`TA` + `A` is the haplotype `TAA`; a
no-called stretch `?` then an inserted `A` is `?A`); `masterVarBeta` writes the joined sequences
itself, one row per locus. Both are read into the same per-haplotype shape and from there the same
code makes the rows, so the two files of one assembly give the same sites, blocks and allele sets;
the haplotype slots can come in another order (masterVarBeta tends to put the variant first).

What a `calls` row of a CGI locus holds (SCHEMA.md Decision 2):

- `source_ref` is the locus reference as written, `alts` the distinct non-reference haplotype
  sequences sorted by string (CGI writes haplotype slots, not an ALT order), `gt` their indices and
  `gt_alleles` the sequences, null for a haplotype written `?` (a half-call keeps its called slot). A
  sequence with `?` inside (`GC?`, `?A`) is a partial call, kept verbatim.
- No anchor base: CGI writes no base before an insertion or deletion, and REF only ever comes from the
  assembly sequence (charter 10). So a deletion's allele is `""`, and a locus of zero length (an
  insertion point) has `source_ref = ""`, sits at the base before it (`pos = begin`, where VCF puts
  it) and has `end_pos` null until RM3 reads that base, as Genos insertion sites do. An insertion
  point CGI could not call (a zero-length no-call) is such a row too, with every haplotype null.
- `sample_filters` are the haplotypes' `varFilter` tokens (`varQuality` before format 2.5, kept
  verbatim), with PASS left out beside another token, as CGI's own VCF writes FT. `ids` are the rs
  ids of `xRef`, sorted. `dp` is masterVarBeta's `totalReadCount`. The allele scores are no likelihood, and
  `hapLink` does not say which haplotype slot is which across loci, so `gl`, `gq` and `phased` stay
  null; those columns are listed in `facets.dropped_fields`.
"""

import logging
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO

import polars as pl

from just_vcf.schema import BLOCKS, CALLS, BlockKind, EvidenceBasis, RefCheck, canonical_chrom_expr, conform

logger = logging.getLogger(__name__)

VAR = "cgi_var"
MASTER_VAR = "cgi_master_var"
FORMATS: frozenset[str] = frozenset({VAR, MASTER_VAR})
"""Format ids this loader reads."""

CALLER = "CGAPipeline"
"""The caller CGI's own VCF names in `##source=CGAPipeline_<version>`; its version is the native
files' `#SOFTWARE_VERSION`."""

READ_EVIDENCE = ("DP",)
"""What `callset.classify`'s read-evidence gate is given for a native file, which has no FORMAT keys.
The gate keeps array and imputed VCFs, whose genotypes come from no reads, from reading as a genome;
a `var` or `masterVarBeta` file is CGAPipeline's own output from reads by construction (allele
scores, masterVarBeta's read counts). The key=value loader passes DP for the same reason. Rejected:
mapping the file's score and count columns onto FORMAT keys, which would name a likelihood the
scores are not."""

NO_CALL = "?"
"""A haplotype CGI could not call, or a stretch of one inside a partly called sequence."""

SAME = "="
"""A sequence CGI did not spell out because it equals the reference (blocks only)."""

_NO_CALL_TYPES = frozenset({"no-call", "no-call-rc", "no-call-ri", "no-ref", "PAR-called-in-X"})
"""varType values of a stretch where nothing was called. `no-ref` is an assembly gap (reference N);
`PAR-called-in-X` is a chrY pseudo-autosomal stretch whose calls are on chrX."""

_FILTER = {VAR: ("varFilter",), MASTER_VAR: ("allele1VarFilter", "allele2VarFilter")}
_FILTER_ALIASES = {"varQuality": "varFilter", "allele1VarQuality": "allele1VarFilter", "allele2VarQuality": "allele2VarFilter"}
"""Format 2.0 and 2.2 name the filter column `…VarQuality` (VQHIGH/VQLOW). The values are kept."""

_CARRIED = {
    VAR: frozenset({"locus", "ploidy", "allele", "chromosome", "begin", "end", "varType", "reference", "alleleSeq", "varFilter", "xRef"}),
    MASTER_VAR: frozenset(
        {
            "locus", "ploidy", "chromosome", "begin", "end", "zygosity", "varType", "reference", "allele1Seq",
            "allele2Seq", "allele1VarFilter", "allele2VarFilter", "allele1XRef", "allele2XRef", "totalReadCount",
        }
    ),
}
"""The columns a frame carries; every other column of the file is listed in `dropped_fields`."""


@dataclass
class Header:
    """The `#KEY\\tvalue` metadata and the `>` column line of a CGI file."""

    meta: dict[str, list[str]] = field(default_factory=dict)
    """Values per key, in file order (`#GENERATED_BY` appears twice). `#CUSTOMER_SAMPLE_ID` and
    `#SAMPLE` can identify the person: never logged."""
    columns: list[str] = field(default_factory=list)
    """Column names as the file writes them (`varQuality` in format 2.0)."""

    def first(self, key: str) -> str | None:
        values = self.meta.get(key)
        return values[0].strip() if values else None


def states_whole_genome(header: Header) -> bool:
    """Whether the header says the file is a whole-genome callset (`#TITLE Whole Human Genome
    Sequencing`; 28 of the 250 PGP-HMS files have no `#TITLE`, and their rows are classified)."""
    return re.search(r"\bWhole (Human )?Genome\b", header.first("TITLE") or "", flags=re.IGNORECASE) is not None


def copy_body(handle: BinaryIO, target: Path) -> Header:
    """Read the header off `handle` and write the data rows after the `>` line to `target`."""
    header = Header()
    for raw in handle:
        line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
        if line.startswith(">"):
            header.columns = line[1:].split("\t")
            break
        if line.startswith("#") and "\t" in line:
            key, value = line[1:].split("\t", 1)
            header.meta.setdefault(key, []).append(value)
    with target.open("wb") as out:
        shutil.copyfileobj(handle, out, 16 << 20)
    return header


def scan(path: Path, header: Header) -> pl.LazyFrame:
    """The data rows as text, columns named as in format 2.5 (`varQuality` read as `varFilter`)."""
    names = [_FILTER_ALIASES.get(name, name) for name in header.columns]
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=names,
        schema={name: pl.String for name in names},
        quote_char=None,
        truncate_ragged_lines=True,
        empty_string_is_null=False,
    )


def dropped_fields(format_id: str, header: Header) -> list[str]:
    """The file's column names that no frame column carries, in file order."""
    return [name for name in header.columns if _FILTER_ALIASES.get(name, name) not in _CARRIED[format_id]]


# --- per-haplotype rows ----------------------------------------------------------------------------


def _with_optional(rows: pl.LazyFrame, names: list[str]) -> pl.LazyFrame:
    present = set(rows.collect_schema().names())
    return rows.with_columns(pl.lit(None, dtype=pl.String).alias(name) for name in names if name not in present)


def _typed(rows: pl.LazyFrame) -> pl.LazyFrame:
    return rows.with_columns(
        pl.col("locus").cast(pl.Int64, strict=False),
        pl.col("begin").cast(pl.Int64, strict=False),
        pl.col("end").cast(pl.Int64, strict=False),
    )


def malformed(rows: pl.LazyFrame) -> pl.Expr:
    """A row whose locus or coordinates are not numbers, or whose end lies before its begin."""
    return (
        pl.col("locus").is_null() | pl.col("begin").is_null() | pl.col("end").is_null() | (pl.col("end") < pl.col("begin"))
    ).fill_null(True)


def prepare(rows: pl.LazyFrame, format_id: str) -> pl.LazyFrame:
    """`scan` output with the optional columns added (null) and the locus and coordinates typed.
    Rows that are then `malformed` are the caller's to count and leave out."""
    optional = ["varFilter", "xRef"] if format_id == VAR else [*_FILTER[MASTER_VAR], "allele1XRef", "allele2XRef", "totalReadCount"]
    return _typed(_with_optional(rows, optional))


def var_parts(rows: pl.LazyFrame) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """A prepared `var` file's uniform rows and its haplotype rows.

    A row with `allele` `all` is uniform when its varType is ref or a no-call kind; one with another
    varType says every haplotype carries its sequence, so it stands for one row per haplotype.
    """
    rows = rows.filter(~malformed(rows))
    every = rows.filter(pl.col("allele") == "all")
    uniform = pl.col("varType").is_in(["ref", *_NO_CALL_TYPES]) & (pl.col("end") > pl.col("begin"))
    expanded = (
        _spelled(every.filter(~uniform), ["alleleSeq"])
        .with_columns(pl.int_ranges(1, pl.col("ploidy").cast(pl.Int64, strict=False).fill_null(0) + 1).alias("allele"))
        .explode("allele", empty_as_null=False)
        .filter(pl.col("allele").is_not_null())
        .with_columns(pl.col("allele").cast(pl.String))
    )
    return every.filter(uniform), pl.concat([rows.filter(pl.col("allele") != "all"), expanded], how="vertical")


def _spelled(rows: pl.LazyFrame, sequences: list[str]) -> pl.LazyFrame:
    """Rows that become haplotypes, with `=` spelled out: a zero-length span's reference is empty, and
    a sequence written `=` is the reference."""
    reference = pl.when(pl.col("begin") == pl.col("end")).then(pl.lit("")).otherwise(pl.col("reference"))
    rows = rows.with_columns(reference.alias("reference"))
    return rows.with_columns(
        pl.when(pl.col(name) == SAME).then(pl.col("reference")).otherwise(pl.col(name)).alias(name) for name in sequences
    )


def var_haplotypes(segments: pl.LazyFrame) -> pl.LazyFrame:
    """One row per haplotype of each `var` locus: its sub-rows joined in `begin` order.

    Carries `tiles`, False where the haplotype's sub-rows leave a gap or overlap; the caller drops
    such loci. A no-called stretch is one `?`, however many sub-rows wrote it.
    """
    ordered = segments.sort("locus", "allele", "begin", "end")
    joined = ordered.group_by("locus", "allele", maintain_order=True).agg(
        pl.col("chromosome").first(),
        pl.col("ploidy").first(),
        pl.col("begin").min(),
        pl.col("end").max(),
        pl.col("alleleSeq").str.join("").str.replace_all(r"\?+", NO_CALL).alias("seq"),
        pl.col("reference").str.join("").alias("reference"),
        (pl.col("begin").slice(1) == pl.col("end").slice(0, pl.len() - 1)).all().alias("tiles"),
        pl.col("varFilter").filter(pl.col("alleleSeq") != NO_CALL).str.split(";").list.explode(keep_nulls=False, empty_as_null=False).alias("filters"),
        pl.col("xRef").str.join(";").alias("xref"),
    )
    return joined.with_columns(
        pl.col("allele").cast(pl.Int64, strict=False).alias("allele"),
        pl.lit(None, dtype=pl.UInt32).alias("dp"),
    )


def master_var_parts(rows: pl.LazyFrame) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """A prepared `masterVarBeta` file's uniform loci and its haplotype rows (allele 1 and, on a diploid locus,
    allele 2), in the shape `var_haplotypes` returns."""
    rows = rows.filter(~malformed(rows))
    one = pl.col("allele1Seq")
    two = pl.col("allele2Seq")
    haploid = pl.col("ploidy") == "1"
    uniform = (
        pl.col("varType").is_in(list(_NO_CALL_TYPES))
        | (pl.col("varType") == "ref") & (one == SAME)
        | (one == NO_CALL) & (haploid | (two == NO_CALL))
    ) & (pl.col("end") > pl.col("begin"))
    blocks = rows.filter(uniform)
    loci =_spelled(rows.filter(~uniform), ["allele1Seq", "allele2Seq"])
    alleles = []
    for number in (1, 2):
        alleles.append(
            loci.filter(pl.lit(number == 1) | ~haploid).select(
                "locus",
                pl.lit(number, dtype=pl.Int64).alias("allele"),
                "chromosome",
                "ploidy",
                "begin",
                "end",
                # A no-called stretch is one `?`, as var_haplotypes joins it (masterVarBeta can write `??`).
                pl.col(f"allele{number}Seq").str.replace_all(r"\?+", NO_CALL).alias("seq"),
                "reference",
                pl.lit(True).alias("tiles"),
                pl.when(pl.col(f"allele{number}Seq").str.contains(r"^\?+$"))
                .then(pl.lit([], dtype=pl.List(pl.String)))
                .otherwise(pl.col(f"allele{number}VarFilter").str.split(";"))
                .alias("filters"),
                pl.col(f"allele{number}XRef").alias("xref"),
                pl.col("totalReadCount").cast(pl.UInt32, strict=False).alias("dp"),
            )
        )
    return blocks, pl.concat(alleles, how="vertical")


# --- rows ------------------------------------------------------------------------------------------


def block_rows(uniform: pl.LazyFrame) -> pl.LazyFrame:
    """`blocks` rows (without `facet_id`) from uniform rows of a positive length, plus `source_chrom`."""
    return uniform.filter(pl.col("end") > pl.col("begin")).select(
        canonical_chrom_expr("chromosome").alias("chrom"),
        pl.col("chromosome").alias("source_chrom"),
        (pl.col("begin") + 1).cast(pl.UInt32).alias("pos"),
        pl.col("end").cast(pl.UInt32).alias("end_pos"),
        pl.when(pl.col("varType") == "ref")
        .then(pl.lit(BlockKind.REF.value))
        .otherwise(pl.lit(BlockKind.NO_CALL.value))
        .alias("kind"),
    )


def locus_checks(haplotypes: pl.LazyFrame) -> pl.LazyFrame:
    """Per locus: whether it can become a row, and why not. Columns `locus`, `problem` (null = fine)."""
    expected = pl.col("ploidy").first().cast(pl.Int64, strict=False)
    return haplotypes.group_by("locus").agg(
        pl.when(expected.is_null() | (expected < 1))
        .then(pl.lit("a ploidy that is not a positive number"))
        .when((pl.col("allele").sort() != pl.int_range(1, pl.len() + 1)).any() | (pl.len() != expected))
        .then(pl.lit("haplotype numbers that are not 1 to the ploidy"))
        .when(~pl.col("tiles").all() | (pl.col("begin").n_unique() > 1) | (pl.col("end").n_unique() > 1))
        .then(pl.lit("haplotype rows that do not tile the locus"))
        .when(pl.col("reference").n_unique() > 1)
        .then(pl.lit("haplotypes that disagree on the reference"))
        .when(pl.col("reference").first().str.contains(SAME, literal=True) | pl.col("seq").str.contains(SAME, literal=True).any())
        .then(pl.lit("a sequence written `=` inside a variant locus"))
        .otherwise(None)
        .alias("problem")
    )


def call_rows(haplotypes: pl.LazyFrame) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """`calls` rows (without `facet_id`) from the haplotype rows of loci that passed `locus_checks`,
    and the loci of a positive length that called no haplotype (no-call blocks).

    A locus of zero length is an insertion point: it covers no base, so it cannot be a block, and
    stays a row whatever it called. A whole no-call there (`gt` all null) says CGI looked for an
    insertion and could not tell.
    """
    called = pl.col("seq") != NO_CALL
    is_alt = called & (pl.col("seq") != pl.col("reference"))
    ranked = haplotypes.with_columns(
        pl.when(is_alt).then(pl.col("seq")).rank("dense").over("locus").alias("_rank"),
        (called.any().over("locus") | (pl.col("end") == pl.col("begin"))).alias("_row"),
    )
    none_called = ranked.filter(~pl.col("_row") & (pl.col("allele") == 1)).select(
        "chromosome", "begin", "end", pl.lit("no-call").alias("varType")
    )
    gt = pl.when(~called).then(None).when(~is_alt).then(0).otherwise(pl.col("_rank")).cast(pl.UInt8)
    loci = (
        ranked.filter(pl.col("_row"))
        .sort("locus", "allele")
        .group_by("locus", maintain_order=True)
        .agg(
            pl.col("chromosome").first(),
            pl.col("begin").first(),
            pl.col("end").first(),
            pl.col("reference").first().alias("source_ref"),
            pl.col("seq").filter(is_alt).unique().sort().alias("alts"),
            gt.alias("gt"),
            pl.when(called).then(pl.col("seq")).otherwise(None).alias("gt_alleles"),
            pl.col("filters").filter(called).list.explode(keep_nulls=False, empty_as_null=False).alias("_filters"),
            pl.col("xref").str.join(";").alias("_xref"),
            pl.col("dp").first(),
        )
    )
    tokens = pl.col("_filters").list.eval(pl.element().filter(pl.element().is_not_null() & (pl.element() != ""))).list.unique().list.sort()
    not_pass = tokens.list.eval(pl.element().filter(pl.element() != "PASS"))
    # Sorted: each haplotype writes its own xRef, so their order says nothing.
    ids = pl.col("_xref").str.extract_all(r"rs\d+").list.unique().list.sort()
    insertion = pl.col("end") == pl.col("begin")
    rows = loci.select(
        canonical_chrom_expr("chromosome").alias("chrom"),
        pl.col("chromosome").alias("source_chrom"),
        pl.when(insertion).then(pl.col("begin")).otherwise(pl.col("begin") + 1).cast(pl.UInt32).alias("pos"),
        pl.when(insertion).then(None).otherwise(pl.col("end")).cast(pl.UInt32).alias("end_pos"),
        pl.when(ids.list.len() > 0).then(ids).otherwise(None).alias("ids"),
        "source_ref",
        pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
        "alts",
        "gt_alleles",
        "gt",
        pl.when(not_pass.list.len() > 0).then(not_pass).when(tokens.list.len() > 0).then(tokens).otherwise(None).alias("sample_filters"),
        "dp",
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
    )
    return rows, none_called


def to_frames(calls: pl.LazyFrame, blocks: pl.LazyFrame, facet_id: str) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    return (
        conform(calls.drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")), CALLS),
        conform(blocks.drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")), BLOCKS),
    )
