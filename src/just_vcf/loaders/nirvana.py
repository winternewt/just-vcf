"""Load Illumina Annotation Engine (Nirvana) JSON (`nirvana_json`) into `calls` rows.

The file is one JSON object written a line at a time: line 1 is the header object and the opening of
`"positions":[`, each following line one position object ending in `,`, until a line starting with `]`
closes the array (gene or OMIM sections may follow; they hold no calls). Each position carries the VCF
record it was annotated from: `chromosome`, `position` (1-based VCF POS), `refAllele` and `altAlleles`
(VCF REF and ALT, with the anchor base), `quality`, `filters`, and one `samples[]` entry per sample of
the header's `samples`, with the VCF GT text in `genotype`, `genotypeQuality`, `totalDepth` and
`alleleDepths`. The `variants[]` list repeats the alleles trimmed, for annotation; it is not read.

So a row is a VCF row (SCHEMA.md Decision 2): `source_ref` the REF, `alts` the ALTs (`[]` for a
position without `altAlleles`, a reference record), `gt` the GT's indices, `gt_alleles` the alleles they
name, `phased` read from the separators as for a VCF before 4.4, `gq`, `dp` and `ad` (only where it has
one value per allele). The header names no assembly in schema v4, and no contig lengths, so the build
stays unknown unless a header states one.
"""

import json
import logging
from pathlib import Path
from typing import BinaryIO

import polars as pl

from just_vcf.loaders.vcf import MultiSampleError, NoSampleError
from just_vcf.schema import CALLS, EvidenceBasis, RefCheck, canonical_chrom_expr, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset({"nirvana_json"})
"""Format ids this loader reads."""

READ_EVIDENCE = ("AD", "DP", "GQ")
"""The VCF FORMAT keys a position's sample fields carry the equivalent of, for RM4's read-evidence gate."""

_SAMPLE = pl.Struct(
    {
        "genotype": pl.String,
        "genotypeQuality": pl.Int64,
        "totalDepth": pl.Int64,
        "alleleDepths": pl.List(pl.Int64),
    }
)
_POSITION = {
    "chromosome": pl.String,
    "position": pl.Int64,
    "refAllele": pl.String,
    "altAlleles": pl.List(pl.String),
    "quality": pl.Float64,
    "filters": pl.List(pl.String),
    "samples": pl.List(_SAMPLE),
}
"""The position fields read; every other key is listed in `dropped_fields`."""

_CARRIED_SAMPLE_KEYS = {"genotype", "genotypeQuality", "totalDepth", "alleleDepths"}

_HEAD_ROWS = 1000
"""Position lines whose keys are read to name the dropped fields."""


class NirvanaLayoutError(ValueError):
    """The first line is not a Nirvana header object opening a `positions` array."""


def copy_positions(handle: BinaryIO, target: Path) -> tuple[dict, list[bytes], int]:
    """Write the position objects to `target` as NDJSON, one per line, and return the header object,
    the first position lines (for the dropped fields) and the count of other lines inside the array."""
    first = handle.readline().rstrip(b"\r\n")
    marker = b',"positions":['
    if not first.endswith(marker):
        raise NirvanaLayoutError('nirvana_json: line 1 does not end by opening the "positions" array')
    header = json.loads(first[: -len(marker)] + b"}").get("header", {})
    head: list[bytes] = []
    other = 0
    with target.open("wb") as out:
        for raw in handle:
            if raw.startswith(b"]"):
                break
            if not raw.startswith(b"{"):
                other += 1
                continue
            line = raw.rstrip(b"\r\n").rstrip(b",")
            if len(head) < _HEAD_ROWS:
                head.append(line)
            out.write(line + b"\n")
    return header, head, other


def samples(header: dict) -> list[str]:
    """The header's sample names, in order. They can be real names: treat them as identifying."""
    return [str(name) for name in header.get("samples", [])]


def sample_index(header: dict, sample: str | None) -> int:
    """The index of the sample to load: the only one, or the one named."""
    names = samples(header)
    if not names:
        raise NoSampleError("nirvana_json: the header lists no sample")
    if sample is None:
        if len(names) > 1:
            raise MultiSampleError(f"nirvana_json: {len(names)} samples; pass sample= to load one")
        return 0
    if sample not in names:
        raise MultiSampleError("nirvana_json: the requested sample is not in the header")
    return names.index(sample)


def dropped_fields(head: list[bytes]) -> list[str]:
    """The position and sample keys of the first rows that no frame column carries, sorted."""
    keys: set[str] = set()
    for line in head:
        position = json.loads(line)
        keys |= {key for key in position if key not in _POSITION}
        for entry in position.get("samples") or []:
            keys |= {f"samples.{key}" for key in entry if key not in _CARRIED_SAMPLE_KEYS}
    return sorted(keys)


def scan(path: Path) -> pl.LazyFrame:
    """The read position fields of every line of the NDJSON copy."""
    return pl.scan_ndjson(path, schema=_POSITION)


def calls_rows(positions: pl.LazyFrame, index: int) -> pl.LazyFrame:
    """`calls` rows (without `facet_id`) for sample `index`, plus `source_chrom`."""
    entry = pl.col("samples").list.get(index, null_on_oob=True)
    text = pl.col("source_gt")
    alts = pl.col("altAlleles").fill_null(pl.lit([], dtype=pl.List(pl.String))).list.eval(
        pl.when(pl.element() == "<NON_REF>").then(pl.lit("<*>")).otherwise(pl.element())
    )
    gt = text.str.replace_all(r"\|", "/").str.split("/").list.eval(
        pl.when(pl.element() == ".").then(None).otherwise(pl.element()).cast(pl.UInt8, strict=False)
    )
    alleles = pl.concat_list(pl.col("source_ref"), pl.col("alts"))
    depths = entry.struct.field("alleleDepths")
    return (
        positions.select(
            canonical_chrom_expr("chromosome").alias("chrom"),
            pl.col("chromosome").alias("source_chrom"),
            pl.col("position").cast(pl.UInt32).alias("pos"),
            (pl.col("position") + pl.col("refAllele").str.len_chars() - 1).cast(pl.UInt32).alias("end_pos"),
            pl.col("refAllele").alias("source_ref"),
            pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
            alts.alias("alts"),
            entry.struct.field("genotype").alias("source_gt"),
            entry.struct.field("genotypeQuality").alias("_gq"),
            entry.struct.field("totalDepth").cast(pl.UInt32, strict=False).alias("dp"),
            depths.alias("_ad"),
            pl.col("quality").cast(pl.Float32).alias("qual"),
            pl.col("filters").alias("filters"),
            pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
        )
        .with_columns(
            gt.alias("gt"),
            pl.concat_list(~text.str.contains("/", literal=True), text.str.extract_all(r"[/|]").list.eval(pl.element() == "|")).alias("phased"),
            pl.when(pl.col("_gq").is_between(0, 65535)).then(pl.col("_gq")).otherwise(None).cast(pl.UInt16).alias("gq"),
        )
        .with_columns(
            # Guarded: `list.gather` over a batch whose every index list is null makes a list of a Null inner
            # type, which the spill's parquet sink cannot write (every position with no genotype failed the load).
            pl.when(pl.col("gt").is_not_null()).then(alleles.list.gather(pl.col("gt"), null_on_oob=True)).alias("gt_alleles"),
            pl.when(pl.col("_ad").list.len() == pl.col("alts").list.len() + 1)
            .then(pl.col("_ad").cast(pl.List(pl.UInt32), strict=False))
            .otherwise(None)
            .alias("ad"),
        )
        .with_columns(pl.when(pl.col("ad").is_not_null()).then(alleles).otherwise(None).alias("source_alleles"))
        .drop("_gq", "_ad")
    )


def stats(rows: pl.LazyFrame) -> dict:
    """Counts the warnings need."""
    return (
        rows.select(
            pl.len().alias("rows"),
            pl.col("source_gt").is_null().sum().alias("no_genotype"),
            (pl.col("gt").list.max() > pl.col("alts").list.len()).fill_null(False).sum().alias("past_alts"),
            (pl.col("source_gt").is_not_null() & ~pl.col("source_gt").str.contains(r"^(\d+|\.)([/|](\d+|\.))*$")).sum().alias("unreadable"),
            pl.col("chrom").is_null().sum().alias("no_chrom"),
        )
        .collect()
        .row(0, named=True)
    )


def warnings_for(counts: dict, other_lines: int) -> list[str]:
    found = []
    if other_lines:
        found.append(f"{other_lines} line(s) inside the positions array that are not a position object skipped")
    if counts["no_genotype"]:
        found.append(f"{counts['no_genotype']} position(s) whose sample has no genotype kept with gt null")
    if counts["unreadable"]:
        found.append(f"{counts['unreadable']} genotype(s) not readable as allele indices kept with a null in gt and gt_alleles for each slot that is not an index, and the text in source_gt")
    if counts["past_alts"]:
        found.append(f"{counts['past_alts']} genotype(s) index an allele past the row's ALT list; gt_alleles has a null there")
    return found


def to_calls_frame(calls: pl.LazyFrame, facet_id: str) -> pl.LazyFrame:
    """`calls_rows` output as a `calls` frame for `facet_id`."""
    return conform(calls.drop("source_chrom").with_columns(pl.lit(facet_id).alias("facet_id")), CALLS)
