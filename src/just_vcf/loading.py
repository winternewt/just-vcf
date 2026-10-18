"""Load a genotype file into the `facets`, `calls` and `blocks` frames (docs/SCHEMA.md).

`load(path)` sniffs the file, reopens the payload `sniff` recognised and hands it to the loader for
that format. A format without a loader is refused with `UnsupportedFormatError`, never read by a
loader for a look-alike format.
"""

import hashlib
import json
import logging
import os
from pathlib import Path

import polars as pl
from platformdirs import user_cache_dir
from pydantic import BaseModel, ConfigDict

from just_vcf.loaders import rsid_table, vcf
from just_vcf.loaders.vcf import MultiSampleError
from just_vcf.schema import (
    BLOCKS,
    FACETS,
    SCHEMA_VERSION,
    AbsenceMeans,
    AssemblyLabelBasis,
    CallsetType,
    canonical_chrom,
    canonical_order,
    canonical_order_per_contig,
    contig_map,
    empty,
)
from just_vcf.sniff import Sniffed, open_payload, sniff

logger = logging.getLogger(__name__)

LOADABLE_FORMATS: frozenset[str] = rsid_table.FORMATS | vcf.FORMATS
"""Format ids `load` can read. Grows as loaders land."""

_SPILL_ROW_GROUP = 200_000
"""Rows per parquet row group in a spill. Small enough that a contig filter skips most of the file."""

_IMPUTED_FORMATS = frozenset({"genes_for_good_imputed"})
_UNKNOWN_CALLSET_FORMATS = frozenset({"plink_23"})
"""PLINK's 23andMe layout is written from any source (the one real file looks imputed or sequenced)."""


class UnsupportedFormatError(ValueError):
    """The file was recognised, but no loader reads its format yet, or it was not recognised at all."""

    def __init__(self, sniffed: Sniffed) -> None:
        self.sniffed = sniffed
        what = sniffed.format or sniffed.kind or sniffed.reason or sniffed.outcome.value
        super().__init__(f"{sniffed.path}: no loader for {what} ({sniffed.outcome.value})")


class Loaded(BaseModel):
    """One loaded input: its facet and its rows. The frames are lazy."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    facets: pl.LazyFrame
    calls: pl.LazyFrame
    blocks: pl.LazyFrame
    warnings: tuple[str, ...] = ()
    """Rows the loader could not take as written, one line per reason with a count."""
    schema_version: str = SCHEMA_VERSION


def load(source: Path | str | Sniffed, sample: str | None = None) -> Loaded:
    """Load `source` (a path, or a `sniff` result for it) into the three frames.

    `sample` picks one sample of a multi-sample VCF; a VCF with several samples and no choice is
    refused with `MultiSampleError`. It is ignored for formats that hold one person by construction.
    """
    sniffed = source if isinstance(source, Sniffed) else sniff(Path(source))
    if not sniffed or sniffed.format not in LOADABLE_FORMATS:
        raise UnsupportedFormatError(sniffed)
    if sniffed.format in vcf.FORMATS:
        return _load_vcf(sniffed, sample)
    spec = rsid_table.spec_for(sniffed.format)
    with open_payload(sniffed) as handle:
        data = handle.read()

    calls, warnings = rsid_table.read_calls(data, spec)
    if spec.columns and rsid_table.column_line_ok(data, spec) is False:
        warnings.append("the column line does not read rsid, chromosome, position, genotype; columns taken by position")
    source_sha256 = _sha256(sniffed.path)
    facet_id = _facet_id(source_sha256, sniffed.member, spec.id, options={})

    callset = (
        CallsetType.IMPUTED
        if spec.id in _IMPUTED_FORMATS
        else None
        if spec.id in _UNKNOWN_CALLSET_FORMATS
        else CallsetType.ARRAY
    )
    facet = {
        "facet_id": facet_id,
        "source_sha256": source_sha256,
        "source_path": str(sniffed.path),
        "member": sniffed.member,
        "format_id": spec.id,
        "family": spec.family.value,
        "vendor": None if spec.vendor == "unknown" else spec.vendor,
        "assembly_label": sniffed.build,
        "assembly_label_basis": AssemblyLabelBasis.STATED.value if sniffed.build else None,
        "callset_type": callset.value if callset else None,
        "absence_means": AbsenceMeans.UNKNOWN.value,
        "contig_map": contig_map(calls["source_chrom"].unique().to_list()),
    }
    facets = pl.LazyFrame([facet], schema={name: FACETS[name] for name in facet}).select(
        [pl.col(name) if name in facet else pl.lit(None, dtype=dtype).alias(name) for name, dtype in FACETS.items()]
    )
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return Loaded(
        facets=facets,
        calls=canonical_order(rsid_table.to_calls_frame(calls, facet_id)),
        blocks=empty(BLOCKS),
        warnings=tuple(warnings),
    )


def cache_dir() -> Path:
    """Where rewritten inputs are kept: `JUST_VCF_CACHE_DIR`, else the platform's user cache."""
    configured = os.environ.get("JUST_VCF_CACHE_DIR")
    return Path(configured) if configured else Path(user_cache_dir("just-vcf"))


def vcf_samples(source: Path | str | Sniffed) -> list[str]:
    """The sample names of a VCF, in file order. They can be real names: treat them as identifying."""
    sniffed = source if isinstance(source, Sniffed) else sniff(Path(source))
    if not sniffed or sniffed.format not in vcf.FORMATS:
        raise UnsupportedFormatError(sniffed)
    with open_payload(sniffed) as handle:
        header, _ = vcf.read_header(handle)
    return header.samples


def _load_vcf(sniffed: Sniffed, sample: str | None) -> Loaded:
    with open_payload(sniffed) as handle:
        header, data = vcf.read_header(handle)
    if len(header.samples) > 1 and sample is None:
        raise vcf.MultiSampleError(
            f"{sniffed.path}: {len(header.samples)} samples; pass sample= to load one "
            "(just_vcf.vcf_samples lists them)"
        )
    if sample is not None and sample not in header.samples:
        raise vcf.MultiSampleError(f"{sniffed.path}: the requested sample is not in the file")
    chosen = sample if sample is not None else (header.samples[0] if header.samples else None)

    warnings: list[str] = []
    missing = vcf.undeclared_format_keys(header, data)
    if missing:
        warnings.append(f"FORMAT key(s) used but not declared, declared by the loader: {', '.join(missing)}")
        header.format.update({key: vcf.FORMAT_KEYS[key][0] for key in missing})
    if not header.format_first:
        warnings.append("##fileformat was not the first line; the header was reordered for reading")
    source_sha256 = _sha256(sniffed.path)
    cache = cache_dir() / "vcf" / hashlib.sha256(f"{source_sha256}!{sniffed.member}!{vcf.REWRITE_VERSION}".encode()).hexdigest()[:32]
    path = vcf.readable_path(sniffed.path, sniffed.layers, lambda: open_payload(sniffed), header, missing, cache)
    try:
        records = vcf.read_records(path, header, chosen, data[0] if data else None)
        stats = vcf.scan_stats(records)
    except vcf.HeaderRefused as err:
        warnings.append(
            f"polars-bio refused the file ({err}); read again with a header of only the carried INFO and "
            "FORMAT keys and the FILTER, ALT and contig lines"
        )
        path = vcf.readable_path(sniffed.path, sniffed.layers, lambda: open_payload(sniffed), header, missing, cache, structural_only=True)
        first = next((line for line in data if line.split(b"\t", 2)[1:2] != [b"0"]), None)
        records = vcf.read_records(path, header, chosen, first)
        stats = vcf.scan_stats(records)
        if left_out := vcf.pos_zero_left_out(cache):
            warnings.append(f"{left_out} record(s) at POS 0 left out: polars-bio cannot read them")
    calls, blocks, row_warnings = vcf.to_rows(records, header, stats)
    warnings += row_warnings
    if (header.version_tuple or (0, 0)) >= (4, 4):
        warnings.append("VCF 4.4+: polars-bio drops a leading GT phase indicator (P4), so phased[0] and source_gt are withheld")

    facet_id = _facet_id(source_sha256, sniffed.member, sniffed.format, options={"sample_index": header.samples.index(chosen) if chosen else None})
    calls_frame, blocks_frame = vcf.to_frames(calls, blocks, facet_id)
    # A global sort holds every column in memory (11 GB for a 6.4M-call VCF). The rows are spilled
    # unsorted in one streaming pass instead, and read back one contig at a time in canonical order.
    spill = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
    spill.mkdir(parents=True, exist_ok=True)
    pl.collect_all(
        [
            calls_frame.sink_parquet(spill / "calls.parquet", row_group_size=_SPILL_ROW_GROUP, lazy=True),
            blocks_frame.sink_parquet(spill / "blocks.parquet", row_group_size=_SPILL_ROW_GROUP, lazy=True),
        ]
    )
    contigs = [canonical_chrom(name) for name in stats["chroms"]]
    spilled = pl.scan_parquet(spill / "calls.parquet")
    flagged, marker_rows = vcf.split_groups(spilled, stats["ploidy_max"] or 1)
    warnings += vcf.split_warnings(flagged, marker_rows, header)
    facet = {
        "facet_id": facet_id,
        "source_sha256": source_sha256,
        "source_path": str(sniffed.path),
        "member": sniffed.member,
        "format_id": sniffed.format,
        "family": sniffed.family.value if sniffed.family else None,
        "sample": chosen,
        "assembly_label": sniffed.build,
        "assembly_label_basis": AssemblyLabelBasis.STATED.value if sniffed.build else None,
        "absence_means": (AbsenceMeans.BLOCKS if stats["any_block"] else AbsenceMeans.UNKNOWN).value,
        "contig_map": contig_map(stats["chroms"]),
        "header_contigs": header.contigs,
        "vcf_version": header.version,
        "sample_meta": header.sample_meta,
        "dropped_fields": vcf.dropped_fields(header),
        "split_marker": ",".join(sorted(header.split_markers)) if marker_rows else None,
    }
    facets = pl.LazyFrame([facet], schema={name: FACETS[name] for name in facet}).select(
        [pl.col(name) if name in facet else pl.lit(None, dtype=dtype).alias(name) for name, dtype in FACETS.items()]
    )
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return Loaded(
        facets=facets,
        calls=canonical_order_per_contig(vcf.with_split_groups(spilled, flagged), contigs),
        blocks=canonical_order_per_contig(pl.scan_parquet(spill / "blocks.parquet"), contigs),
        warnings=tuple(warnings),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _facet_id(source_sha256: str, member: str | None, format_id: str, options: dict) -> str:
    """Deterministic facet id: the source bytes, the member, the format and the load options."""
    key = json.dumps(
        {"sha256": source_sha256, "member": member, "format": format_id, "options": options},
        sort_keys=True,
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]
