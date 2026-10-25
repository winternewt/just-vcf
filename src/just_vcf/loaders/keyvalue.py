"""Load PGP-HMS's freeBayes+VarScan2 key=value table (`freebayes_varscan2_keyvalue`), one facet per caller.

The table sits beside PGP's merged VCF in the same arvados collection (`<id>-orig.vcf.gz`, which is not
a VCF). Each row is one site: CHROM, POS, REF, ALT, then four `KEY=value` fields per caller in a fixed
order (genotype, reference allele count, reference allele quality, total depth; freeBayes first, then
VarScan2), then an annotation tail of any length. A caller that did not report the site writes `NULL`
in all four of its fields. Unlike the one-column merges (SCHEMA.md Decision 8), it keeps VarScan2's
genotype where both callers found a site.

Each caller's rows become one facet (`relation = caller_split`): `gt` indexes `[REF] + ALT` as the
genotype text does, `source_gt` keeps that text, and `dp` is the caller's total depth. The reference
allele counts and qualities have no column (a count for REF alone is not an `ad` vector), and the
annotation tail is not read; both are listed in `facets.dropped_fields`.
"""

import logging
import shutil
from pathlib import Path
from typing import BinaryIO

import polars as pl

from just_vcf.schema import CALLS, EvidenceBasis, RefCheck, canonical_chrom_expr, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset({"freebayes_varscan2_keyvalue"})
"""Format ids this loader reads."""

CALLERS: tuple[tuple[int, str, str], ...] = ((0, "FB", "freeBayes"), (1, "VS", "VarScan2"))
"""Each caller's source index (as INFO/SF numbers it in the sibling merge), key prefix and name."""

_SITE = ("chromosome", "position", "ref", "alt")
_KEYS = ("GenoType", "RefAlleleCount", "ReferenceAlleleQ", "TotalDepth")
_FIXED = (*_SITE, *(f"{prefix}{key}" for _, prefix, _ in CALLERS for key in _KEYS))
"""The twelve fixed fields: the site, then each caller's four keyed fields."""

NULL = "NULL"
"""The table's spelling of a missing value."""

_TAIL_ROWS = 1000
"""Rows whose annotation tail is read to name its keys in `dropped_fields`."""


def copy_payload(handle: BinaryIO, target: Path) -> list[bytes]:
    """Write the payload to `target` as it is and return its first rows, for the tail's key names."""
    head = [handle.readline() for _ in range(_TAIL_ROWS)]
    with target.open("wb") as out:
        out.writelines(head)
        shutil.copyfileobj(handle, out, 16 << 20)
    return [line for line in head if line]


def tail_keys(head: list[bytes]) -> list[str]:
    """The annotation tail's key names in the first rows (`KG`, `EFF_HGVS`, bare tokens such as
    `EFF_RESULTS`), sorted."""
    keys: set[str] = set()
    for raw in head:
        for field in raw.rstrip(b"\r\n").split(b"\t")[len(_FIXED) :]:
            key = field.split(b"=", 1)[0].decode("utf-8", errors="replace").strip()
            if key:
                keys.add(key)
    return sorted(keys)


def dropped_fields(head: list[bytes]) -> list[str]:
    """Every field the frames do not carry: the reference counts and qualities, and the tail's keys."""
    carried = {"GenoType", "TotalDepth"}
    fixed = [f"{prefix}{key}" for _, prefix, _ in CALLERS for key in _KEYS if key not in carried]
    return sorted(fixed) + [key for key in tail_keys(head) if key not in fixed]


def scan(path: Path) -> pl.LazyFrame:
    """The twelve fixed fields of every row, as text with each `KEY=` prefix still on."""
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=list(_FIXED),
        schema={name: pl.String for name in _FIXED},
        quote_char=None,
        truncate_ragged_lines=True,
        empty_string_is_null=False,
    )


def malformed(rows: pl.LazyFrame) -> pl.Expr:
    """True on a row whose keyed fields do not carry their keys in the fixed order, or whose position
    is not a 32-bit count."""
    keyed = [pl.col(name).str.starts_with(f"{name}=") for name in _FIXED[len(_SITE) :]]
    position = pl.col("position").cast(pl.UInt32, strict=False).is_not_null()
    return ~(pl.all_horizontal(keyed) & position).fill_null(False)


def _value(name: str) -> pl.Expr:
    """A keyed field's value, null where the table writes NULL (or freeBayes's `.`)."""
    value = pl.col(name).str.strip_prefix(f"{name}=")
    return pl.when(value.is_in([NULL, "."])).then(None).otherwise(value)


def caller_rows(rows: pl.LazyFrame, prefix: str) -> pl.LazyFrame:
    """One caller's `calls` rows (without `facet_id`): the well-formed rows where its genotype is
    not NULL, plus `source_chrom` for the facet's `contig_map`."""
    genotype = _value(f"{prefix}GenoType")
    text = pl.col("source_gt")
    seps = text.str.extract_all(r"[/|]")
    gt = text.str.replace_all(r"\|", "/").str.split("/").list.eval(
        pl.when(pl.element() == ".").then(None).otherwise(pl.element()).cast(pl.UInt8, strict=False)
    )
    alleles = pl.concat_list(pl.col("source_ref"), pl.col("alts"))
    return (
        rows.filter(~malformed(rows))
        .with_columns(genotype.alias("source_gt"))
        .filter(pl.col("source_gt").is_not_null())
        .select(
            canonical_chrom_expr("chromosome").alias("chrom"),
            pl.col("chromosome").alias("source_chrom"),
            pl.col("position").cast(pl.UInt32).alias("pos"),
            (pl.col("position").cast(pl.UInt32) + pl.col("ref").str.len_chars().cast(pl.UInt32) - 1).alias("end_pos"),
            pl.col("ref").alias("source_ref"),
            pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
            pl.col("alt").str.split(",").alias("alts"),
            pl.col("source_gt"),
            gt.alias("gt"),
            # The first haplotype's phase is read from the first separator, as for VCF before 4.4.
            pl.concat_list(~text.str.contains("/", literal=True), seps.list.eval(pl.element() == "|")).alias("phased"),
            _value(f"{prefix}TotalDepth").cast(pl.UInt32, strict=False).alias("dp"),
            pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
        )
        .with_columns(alleles.list.gather(pl.col("gt"), null_on_oob=True).alias("gt_alleles"))
    )


def stats(rows: pl.LazyFrame) -> dict:
    """Counts the warnings need, over every row: malformed rows, genotypes per caller, and genotype
    indices past a row's ALT list."""
    alts = pl.col("alt").str.split(",").list.len()
    exprs = [malformed(rows).sum().alias("malformed")]
    for _, prefix, _ in CALLERS:
        genotype = _value(f"{prefix}GenoType")
        indices = genotype.str.extract_all(r"\d+").list.eval(pl.element().cast(pl.Int64)).list.max()
        exprs += [
            (genotype.is_not_null() & ~malformed(rows)).sum().alias(f"{prefix}_rows"),
            ((indices > alts) & ~malformed(rows)).fill_null(False).sum().alias(f"{prefix}_past_alts"),
            (genotype.is_not_null() & ~genotype.str.contains(r"^(\d+|\.)([/|](\d+|\.))*$") & ~malformed(rows))
            .fill_null(False)
            .sum()
            .alias(f"{prefix}_unreadable"),
        ]
    exprs.append(pl.col("chromosome").filter(~malformed(rows)).unique().sort().alias("chroms").implode())
    return rows.select(exprs).collect().row(0, named=True)


def warnings_for(counts: dict, prefix: str) -> list[str]:
    """One caller's warnings, one per reason with its count."""
    found = []
    if counts[f"{prefix}_past_alts"]:
        found.append(f"{counts[f'{prefix}_past_alts']} genotype(s) index an allele past the row's ALT list; gt_alleles has a null there")
    if counts[f"{prefix}_unreadable"]:
        found.append(f"{counts[f'{prefix}_unreadable']} genotype(s) not readable as allele indices kept with a null in gt and gt_alleles for each slot that is not an index, and the text in source_gt")
    return found


def to_calls_frame(calls: pl.LazyFrame, facet_id: str) -> pl.LazyFrame:
    """`caller_rows` output as a `calls` frame for `facet_id`."""
    return conform(calls.drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")), CALLS)
