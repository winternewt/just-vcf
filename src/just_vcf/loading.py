"""Load a genotype file into the `facets`, `calls` and `blocks` frames (docs/SCHEMA.md).

`load(path)` sniffs the file, reopens the payload `sniff` recognised and hands it to the loader for
that format. A format without a loader is refused with `UnsupportedFormatError`, never read by a
loader for a look-alike format.
"""

import hashlib
import json
import logging
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf import cache, callset
from just_vcf.assembly import infer_assembly, sorted_name_length_pairs_digest
from just_vcf.cache import cache_dir, file_sha256
from just_vcf.formats import get_spec
from just_vcf.formats.spec import Family
from just_vcf.loaders import rsid_table, table, vcf
from just_vcf.loaders.vcf import MultiSampleError, NoSampleError
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
from just_vcf.source import detect_source

logger = logging.getLogger(__name__)

LOADABLE_FORMATS: frozenset[str] = rsid_table.FORMATS | table.FORMATS | vcf.FORMATS
"""Format ids `load` can read. Grows as loaders land."""

SPILL_ROW_GROUP = 200_000
"""Rows per parquet row group in a spill. Small enough that a contig filter skips most of the file."""

_IMPUTED_FORMATS = frozenset({"genes_for_good_imputed"})
_UNKNOWN_CALLSET_FORMATS = frozenset({"plink_23", "genos", "color_genotypes"})
"""Formats whose callset type the file does not say. PLINK's 23andMe layout is written from any source
(the one real file looks imputed or sequenced). Genos reports exome sequencing at a request list of
sites, and Color does not say which assay its trait SNPs come from; RM4 infers either."""


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

    `sample` picks one sample of a multi-sample VCF or of a table with a sample column (an Illumina
    Final Report); several samples and no choice is refused with `MultiSampleError`, and a sites-only
    VCF, which holds no sample, with `NoSampleError`. A table whose columns lack what its format needs (a
    Final Report without positions or plus-strand alleles) is refused with `TableLayoutError`. `sample`
    is ignored for formats that hold one person by construction.

    A VCF's rows are spilled to parquet under `cache_dir()`, and an input polars-bio cannot read as
    it is is first copied there; the lazy frames read those files. Each run first removes the cache
    entries an older `REWRITE_VERSION` or `SCHEMA_VERSION` wrote (`just_vcf.cache.sweep_stale`).
    """
    cache.sweep_stale(cache.REWRITE_COPY, {"REWRITE_VERSION": vcf.REWRITE_VERSION})
    cache.sweep_stale(cache.FRAMES_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION})
    sniffed = source if isinstance(source, Sniffed) else sniff(Path(source))
    if not sniffed or sniffed.format not in LOADABLE_FORMATS:
        raise UnsupportedFormatError(sniffed)
    if sniffed.format in vcf.FORMATS:
        return _load_vcf(sniffed, sample)
    if sniffed.format in table.FORMATS:
        return _load_table(sniffed, sample)
    return _load_rsid_table(sniffed)


def _load_rsid_table(sniffed: Sniffed) -> Loaded:
    spec = rsid_table.spec_for(sniffed.format)
    with open_payload(sniffed) as handle:
        data = handle.read()
    calls, warnings = rsid_table.read_calls(data, spec)
    if spec.columns and rsid_table.column_line_ok(data, spec) is False:
        warnings.append("the column line does not read rsid, chromosome, position, genotype; columns taken by position")
    return _array_loaded(sniffed, calls, warnings, build=sniffed.build, options={}, extra={})


def _load_table(sniffed: Sniffed, sample: str | None) -> Loaded:
    spec = table.spec_for(sniffed.format)
    with open_payload(sniffed) as handle:
        data = handle.read()
    read = table.read_table(data, spec, sample)
    build = read.build if read.rows_state_build else sniffed.build or _sibling_build(sniffed)
    options = {} if read.sample_index is None else {"sample_index": read.sample_index}
    extra = {"sample": read.sample, "dropped_fields": read.dropped_columns}
    return _array_loaded(sniffed, read.calls, read.warnings, build=build, options=options, extra=extra)


def _sibling_build(sniffed: Sniffed) -> str | None:
    """The build a metadata member beside this one in the same archive states (deCODEme's info file)."""
    stated = sorted({m.build for m in sniffed.members if m.family is Family.METADATA and m.build})
    return stated[0] if len(stated) == 1 else None


def _array_loaded(
    sniffed: Sniffed,
    calls: pl.DataFrame,
    warnings: list[str],
    *,
    build: str | None,
    options: dict,
    extra: dict,
) -> Loaded:
    """The three frames of a genotype table: an array-like facet, the calls, no blocks."""
    source_sha256 = file_sha256(sniffed.path)
    facet_id = _facet_id(source_sha256, sniffed.member, sniffed.format, options=options)
    callset = (
        CallsetType.IMPUTED
        if sniffed.format in _IMPUTED_FORMATS
        else None
        if sniffed.format in _UNKNOWN_CALLSET_FORMATS
        else CallsetType.ARRAY
    )
    spec_vendor = get_spec(sniffed.format).vendor
    facet = {
        "facet_id": facet_id,
        "source_sha256": source_sha256,
        "source_path": str(sniffed.path),
        "member": sniffed.member,
        "format_id": sniffed.format,
        "family": sniffed.family.value if sniffed.family else None,
        "vendor": None if spec_vendor == "unknown" else spec_vendor,
        "assembly_label": build,
        "assembly_label_basis": AssemblyLabelBasis.STATED.value if build else None,
        "callset_type": callset.value if callset else None,
        "absence_means": AbsenceMeans.UNKNOWN.value,
        "contig_map": contig_map(calls["source_chrom"].unique().to_list()),
        **{name: value for name, value in extra.items() if value is not None},
    }
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return Loaded(
        facets=_facet_frame(facet),
        calls=canonical_order(rsid_table.to_calls_frame(calls, facet_id)),
        blocks=empty(BLOCKS),
        warnings=tuple(warnings),
    )


def _facet_frame(facet: dict) -> pl.LazyFrame:
    """One `facets` row from the columns a loader knows; every other column null."""
    return pl.LazyFrame([facet], schema={name: FACETS[name] for name in facet}).select(
        [pl.col(name) if name in facet else pl.lit(None, dtype=dtype).alias(name) for name, dtype in FACETS.items()]
    )


def vcf_samples(source: Path | str | Sniffed) -> list[str]:
    """The sample names of a VCF, in file order. They can be real names: treat them as identifying."""
    sniffed = source if isinstance(source, Sniffed) else sniff(Path(source))
    if not sniffed or sniffed.format not in vcf.FORMATS:
        raise UnsupportedFormatError(sniffed)
    with open_payload(sniffed) as handle:
        header = vcf.read_bcf_header(handle)[0] if sniffed.format == vcf.BCF else vcf.read_header(handle)[0]
    return header.samples


def _load_vcf(sniffed: Sniffed, sample: str | None) -> Loaded:
    binary = sniffed.format == vcf.BCF
    with open_payload(sniffed) as handle:
        if binary:
            header, first_record = vcf.read_bcf_header(handle)
            data: list[bytes] = []
        else:
            header, data = vcf.read_header(handle)
            first_record = data[0] if data else None
    warnings: list[str] = []
    if sniffed.format == vcf.FREEBAYES_VARSCAN2:
        # One person, two callers (SCHEMA Decision 8): one facet per caller, whatever the columns.
        parts = vcf.caller_parts(header, data)
    else:
        if len(header.samples) > 1 and sample is None:
            raise vcf.MultiSampleError(
                f"{sniffed.path}: {len(header.samples)} samples; pass sample= to load one "
                "(just_vcf.vcf_samples lists them)"
            )
        if sample is not None and sample not in header.samples:
            raise vcf.MultiSampleError(f"{sniffed.path}: the requested sample is not in the file")
        if not header.samples:
            raise vcf.NoSampleError(f"{sniffed.path}: a sites-only VCF, with no sample column to load")
        chosen = sample if sample is not None else header.samples[0]
        parts = [vcf.Part(sample=chosen, options={"sample_index": header.samples.index(chosen)})]

    missing = vcf.undeclared_keys(header, data)
    if missing:
        warnings.append(f"key(s) used but not declared, declared by the loader: {', '.join(missing)}")
        vcf.declare(header, missing)
    if unfit := vcf.unfit_keys(header):
        warnings.append(
            "FORMAT key(s) declared with a Number that is not one value per allele or genotype, not carried: "
            + ", ".join(f"{key} (Number={number})" for key, number in sorted(unfit.items()))
        )
    if not header.meta:
        warnings.append(
            "no meta lines: read under a synthesised header (##fileformat=VCFv4.2 and the carried keys the "
            "rows use); the file states no VCF version, no contigs and no build"
        )
    elif not header.format_first:
        warnings.append("##fileformat was not the first line; the header was reordered for reading")
    source_sha256 = file_sha256(sniffed.path)
    copy = cache_dir() / "vcf" / hashlib.sha256(f"{source_sha256}!{sniffed.member}!{vcf.REWRITE_VERSION}".encode()).hexdigest()[:32]
    origin = cache.Origin(source_sha256=source_sha256, source_path=str(sniffed.path.resolve()), member=sniffed.member)
    source = detect_source(header.meta, get_spec(sniffed.format).vendor)
    warnings += source.warnings
    gl_phred = source.caller == "CGAPipeline"
    if gl_phred:
        warnings.append("CGAPipeline writes GL in phred units (GQ is the gap between the two best): read as GL / 10")
    if (header.version_tuple or (0, 0)) >= (4, 4):
        warnings.append("VCF 4.4+: polars-bio drops a leading GT phase indicator (P4), so phased[0] and source_gt are withheld")
    versions = {caller.name: caller.version for caller in source.callers}
    label, basis, label_warnings = _assembly_label(sniffed.build, header.contigs)
    warnings += label_warnings

    facets: list[dict] = []
    facet_stats: list[dict] = []
    calls_parts: list[pl.LazyFrame] = []
    blocks_parts: list[pl.LazyFrame] = []
    for part in parts:
        part_warnings: list[str] = []
        if binary:
            # A BCF declares every key it uses, and its header cannot be rewritten without re-encoding
            # the records, so a refusal from polars-bio is final here.
            path = vcf.readable_bcf_path(sniffed.path, sniffed.layers, lambda: open_payload(sniffed), copy, origin)
            records = part.select(vcf.read_records(path, header, part.sample, first_record, extra_info=part.info))
            stats = vcf.scan_stats(records)
        else:
            records, stats = _read_text_records(sniffed, header, missing, data, copy, origin, part, part_warnings)
        calls, blocks, row_warnings = vcf.to_rows(records, header, stats, text_gt=not binary, gl_phred=gl_phred)
        part_warnings += row_warnings

        facet_id = _facet_id(source_sha256, sniffed.member, sniffed.format, options=part.options)
        calls_frame, blocks_frame = vcf.to_frames(calls, blocks, facet_id)
        # A global sort holds every column in memory (11 GB for a 6.4M-call VCF). The rows are spilled
        # unsorted in one streaming pass instead, and read back one contig at a time in canonical order.
        spill = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
        with cache.writing(spill, cache.FRAMES_SPILL, {"SCHEMA_VERSION": SCHEMA_VERSION}, origin.model_copy(update={"facet_id": facet_id})):
            spill.mkdir(parents=True, exist_ok=True)
            pl.collect_all(
                [
                    calls_frame.sink_parquet(spill / "calls.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
                    blocks_frame.sink_parquet(spill / "blocks.parquet", row_group_size=SPILL_ROW_GROUP, lazy=True),
                ]
            )
        contigs = [canonical_chrom(name) for name in stats["chroms"]]
        spilled = pl.scan_parquet(spill / "calls.parquet")
        flagged, marker_rows = vcf.split_groups(spilled, stats["ploidy_max"] or 1)
        part_warnings += vcf.split_warnings(flagged, marker_rows, header)
        if source.callset_type:
            kind = source.callset_type
        else:
            inference = callset.classify(spilled, list(header.format), callset.imputation_marker(header.meta, list(header.info)))
            logger.info("%s: callset type %s (%s)", facet_id, inference.callset_type.value if inference else "unknown", inference.reason)
            kind = inference.callset_type if inference else None
        caller = part.caller if part.split else source.caller
        facets.append(
            {
                "facet_id": facet_id,
                "source_sha256": source_sha256,
                "source_path": str(sniffed.path),
                "member": sniffed.member,
                "format_id": sniffed.format,
                "family": sniffed.family.value if sniffed.family else None,
                "vendor": source.vendor,
                "caller": caller,
                "caller_version": versions.get(caller) if part.split else source.caller_version,
                "ploidy_statements": source.ploidy_statements,
                # What the header states (array or imputed content) wins; RM4 infers the rest from the rows.
                "callset_type": kind.value if kind else None,
                "phasing_method": source.phasing_method,
                "relation": vcf.CALLER_SPLIT if part.split else None,
                "sample": part.sample,
                "assembly_label": label,
                "assembly_label_basis": basis.value if basis else None,
                "sorted_name_length_pairs_digest": sorted_name_length_pairs_digest(header.contigs),
                # A one-column merge facet holds only the sites its caller alone found: absent there is unknown.
                "absence_means": (AbsenceMeans.UNKNOWN if part.sf_alone else callset.absence_means(kind, stats["any_block"])).value,
                "contig_map": contig_map(stats["chroms"]),
                "header_contigs": header.contigs,
                "vcf_version": header.version,
                "sample_meta": header.sample_meta,
                "dropped_fields": vcf.dropped_fields(header, data),
                "split_marker": ",".join(sorted(header.split_markers)) if marker_rows else None,
            }
        )
        # The flagged positions join above the per-contig sort, so a filter still reaches the parquet scan.
        calls_parts.append(vcf.with_split_groups(canonical_order_per_contig(spilled, contigs), flagged))
        blocks_parts.append(canonical_order_per_contig(pl.scan_parquet(spill / "blocks.parquet"), contigs))
        if part.sf_alone:
            shared = facet_stats[0].get("shared_sources", 0)
            part_warnings.append(
                f"one data column: at the {shared} site(s) another caller also found, the merge kept that caller's "
                "genotype, so this facet holds only the sites it alone found; absent here is unknown"
            )
        facet_stats.append(stats)
        warnings += [f"{part.caller or 'unnamed caller'} facet: {w}" if part.split else w for w in part_warnings]

    order = sorted(range(len(facets)), key=lambda i: facets[i]["facet_id"])
    frame = pl.LazyFrame([facets[i] for i in order], schema={name: FACETS[name] for name in facets[0]}).select(
        [pl.col(name) if name in facets[0] else pl.lit(None, dtype=dtype).alias(name) for name, dtype in FACETS.items()]
    )
    for warning in warnings:
        logger.warning("%s: %s", facets[order[0]]["facet_id"], warning)
    return Loaded(
        facets=frame,
        calls=_concat([calls_parts[i] for i in order]),
        blocks=_concat([blocks_parts[i] for i in order]),
        warnings=tuple(warnings),
    )


def _concat(frames: list[pl.LazyFrame]) -> pl.LazyFrame:
    """One facet's frame as it is; several, one after another (they are already in facet order)."""
    return frames[0] if len(frames) == 1 else pl.concat(frames, how="vertical", parallel=False)


def _read_text_records(
    sniffed: Sniffed,
    header: vcf.Header,
    missing: list[str],
    data: list[bytes],
    copy: Path,
    origin: cache.Origin,
    part: "vcf.Part",
    warnings: list[str],
) -> tuple[pl.LazyFrame, dict]:
    """A text VCF's records for one part and their stats, read once more from a structural copy if
    polars-bio refuses the first read (P7, P8). The retry is reported in `warnings`."""
    path = vcf.readable_path(sniffed.path, sniffed.layers, lambda: open_payload(sniffed), header, missing, copy, origin)
    try:
        records = part.select(vcf.read_records(path, header, part.sample, data[0] if data else None, extra_info=part.info))
        return records, vcf.scan_stats(records)
    except vcf.HeaderRefused as err:
        warnings.append(
            f"polars-bio refused the file ({err}); read again with a header of only the carried INFO and "
            "FORMAT keys and the FILTER, ALT and contig lines"
        )
    path = vcf.readable_path(sniffed.path, sniffed.layers, lambda: open_payload(sniffed), header, missing, copy, origin, structural_only=True)
    first = next((line for line in data if line.split(b"\t", 2)[1:2] != [b"0"]), None)
    records = part.select(vcf.read_records(path, header, part.sample, first, extra_info=part.info))
    stats = vcf.scan_stats(records)
    if left_out := vcf.pos_zero_left_out(copy):
        warnings.append(f"{left_out} record(s) at POS 0 left out: polars-bio cannot read them")
    return records, stats


def _assembly_label(stated: str | None, contigs: list[dict]) -> tuple[str | None, AssemblyLabelBasis | None, list[str]]:
    """The facet's assembly label and its basis, plus a warning for each disagreement.

    A label the header states is kept as stated, even when the contig lengths point elsewhere: the
    disagreement is reported, never resolved here. Without a stated label, the contig lengths give
    one only when they settle it (`infer_assembly`); otherwise the label stays unknown.
    """
    inferred = infer_assembly(contigs)
    warnings = []
    if len(inferred.matched) > 1:
        warnings.append(f"##contig lengths match more than one assembly ({inferred.reason})")
    if stated and inferred and inferred.label != stated:
        warnings.append(
            f"the header states {stated}, but {inferred.matched[inferred.label]} primary chromosome length(s) "
            f"match {inferred.label}; the stated label is kept"
        )
    if stated:
        return stated, AssemblyLabelBasis.STATED, warnings
    if inferred:
        return inferred.label, AssemblyLabelBasis.INFERRED, warnings
    return None, None, warnings


def _facet_id(source_sha256: str, member: str | None, format_id: str, options: dict) -> str:
    """Deterministic facet id: the source bytes, the member, the format and the load options."""
    key = json.dumps(
        {"sha256": source_sha256, "member": member, "format": format_id, "options": options},
        sort_keys=True,
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]
