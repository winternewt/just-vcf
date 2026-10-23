"""Export one facet's calls as a VCF Zarr store (RM276; docs/SCHEMA.md, "VCF Zarr export").

`write_vcf_zarr(loaded, path)` writes a VCF Zarr 0.4 store in Zarr format 2, which polars-bio reads
back with `scan_vcf_zarr`. It is an export at the edge, like `write_vcf`, and follows its rules for
which rows can be written: REF is `ref`, else `source_ref`, never the genotype (charter 10), and a row
VCF cannot hold is refused, counted per reason, unless `resolved_only=True`.

Written with the standard library alone: the metadata is JSON, strings are `vlen-utf8`, chunks are
`gzip` without a timestamp, so two writes of one frame give the same bytes (charter 12) and no extra
is needed. polars-bio's reader takes `gzip` (and no compressor); it refuses `zlib`.

**Blocks are not representable.** The VCF Zarr spec has no reference or no-call blocks (a gVCF's
`<*>` records): a facet's `blocks` stay in the parquet store (`just_vcf.store`) and are counted in
`blocks_not_written` and a warning, never written as calls.

What is written:

- `contig_id` (and `contig_length` where the facet's header gave every length), `variant_contig`,
  `variant_position` (1-based), `variant_length` (`end_pos - pos + 1`, so a symbolic ALT keeps its
  span), `variant_allele` (REF then ALTs, padded with `""`), `variant_id` with `variant_id_mask`,
  `variant_quality`, `filter_id` (`PASS` first) and `variant_filter`;
- `sample_id`, `call_genotype` (`-1` an uncalled haplotype, `-2` past the row's ploidy) and
  `call_genotype_phased` (true where every separator after the first slot is `|`, as `write_vcf`
  writes GT);
- the FORMAT fields any row has: `call_GQ`, `call_DP`, `call_AD`, `call_GL`, `call_GP`, `call_DS`,
  `call_PS` (PS-kind phase sets, as integers), `call_FT`, `call_CN`; and the INFO fields any row has:
  `variant_SVTYPE`, `variant_JUST_VCF_SPLIT`. `call_AD`, `call_GL` and `call_GP` are written only
  where they are indexed over the written REF and ALTs, as in `write_vcf`.

Missing values are the spec's: `-1` for integers, the NaN `0x7F800001` for floats, `"."` for
strings; padding is `-2`, the NaN `0x7F800002` and `""`. The root's attributes name the facet
(`facet_id`, `source_sha256`, `format_id`, `derived_from`, the schema version) and never its
`source_path`. Fields not written (`cipos`, `ciend`, `alt_copy_numbers`, `dropped_fields`) are counted
in a warning where a row has them.
"""

import gzip
import json
import logging
import os
import shutil
import struct
import sys
from array import array
from importlib.metadata import version
from math import prod
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.loaders.vcf import SPLIT_INFO_KEY
from just_vcf.loading import Loaded
from just_vcf.schema import SCHEMA_VERSION
from just_vcf.writing import (
    ContigNames,
    UnsortedRowsError,
    UnwritableRowsError,
    _contig_names,
    _Names,
    _ps_text,
    _unwritable,
    _UNWRITABLE_REASONS,
    _vectors_indexable,
    _written_ref,
)

logger = logging.getLogger(__name__)

VCF_ZARR_VERSION = "0.4"
"""The VCF Zarr spec version written, the one polars-bio reads."""

VARIANTS_CHUNK_SIZE = 10_000
"""Variants per chunk."""

_INT_MISSING, _INT_FILL = -1, -2
_FLOAT_MISSING, _FLOAT_FILL = 0x7F800001, 0x7F800002
"""The spec's float sentinels, as the bits of a float32 NaN."""
_STR_MISSING, _STR_FILL = ".", ""
_COMPRESSOR = {"id": "gzip", "level": 6}
_VLEN = [{"id": "vlen-utf8"}]
_TYPECODES = {"|i1": "b", "<i2": "h", "<i4": "i", "<i8": "q"}

_DESCRIPTIONS = {
    "GQ": "Conditional genotype quality",
    "DP": "Read depth of this sample",
    "AD": "Read depth for each allele",
    "GL": "Genotype likelihoods, log10",
    "GP": "Genotype posterior probabilities",
    "DS": "Expected count of the alternate allele",
    "PS": "Phase set",
    "FT": "Sample filter",
    "CN": "Copy number of this sample",
}


class VcfZarrWritten(BaseModel):
    """What `write_vcf_zarr` wrote."""

    model_config = ConfigDict(frozen=True)

    path: Path
    rows_written: int
    rows_skipped: dict[str, int]
    """Rows left out under `resolved_only=True`, per reason. Empty when nothing was left out."""
    blocks_not_written: int
    """Every block of the facet: VCF Zarr has no reference or no-call blocks."""
    arrays: tuple[str, ...]
    """The arrays written, sorted."""
    warnings: tuple[str, ...] = ()


def write_vcf_zarr(
    loaded: Loaded,
    path: Path | str,
    *,
    contig_names: ContigNames = ContigNames.STORED,
    sample: str | None = None,
    resolved_only: bool = False,
    variants_chunk_size: int = VARIANTS_CHUNK_SIZE,
) -> VcfZarrWritten:
    """Write `loaded`'s one facet's calls to `path` as a VCF Zarr 0.4 store (a directory).

    `path` must not exist: a store is written whole beside it and renamed into place, never over
    another. `sample` and `contig_names` are `write_vcf`'s. The rows must be in the order `load`
    emits them; rows VCF cannot hold raise `UnwritableRowsError` unless `resolved_only=True`. The
    selected columns of the rows written are held in memory while the chunks are encoded.
    """
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"{path.name} exists; a VCF Zarr store is never written over another")
    if variants_chunk_size < 1:
        raise ValueError("variants_chunk_size must be at least 1")
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"write_vcf_zarr takes one facet; this Loaded has {facets.height}")
    facet = facets.row(0, named=True)
    facet_id = facet["facet_id"]
    sample_name = sample if sample is not None else facet["sample"] if facet["sample"] is not None else facet_id

    calls = loaded.calls.with_columns(_unwritable().alias("_unwritable"))
    counts = calls.select(
        *[(pl.col("_unwritable") == reason).sum().alias(reason) for reason in _UNWRITABLE_REASONS],
        *[pl.col(name).list.len().fill_null(0).gt(0).sum().alias(f"n_{name}") for name in ("cipos", "ciend", "alt_copy_numbers")],
    ).collect().row(0, named=True)
    skipped = {reason: counts[reason] for reason in _UNWRITABLE_REASONS if counts[reason]}
    if skipped and not resolved_only:
        raise UnwritableRowsError(facet_id, skipped)
    rows = (calls.filter(pl.col("_unwritable").is_null()) if skipped else calls).select(_columns()).collect()
    _check_order(rows)
    block_count = loaded.blocks.select(pl.len()).collect().item()

    names = _contig_names(contig_names, facet, rows["chrom"].unique().to_list())
    contig_index = {chrom: index for index, chrom in enumerate(names.mapping)}
    warnings: list[str] = []
    if skipped:
        warnings.append(f"{sum(skipped.values())} row(s) not written: " + ", ".join(f"{n} {r}" for r, n in skipped.items()))
    if block_count:
        warnings.append(
            f"{block_count} reference/no-call block(s) not written: VCF Zarr has no blocks; they stay in the parquet store"
        )
    not_written = {name: counts[f"n_{name}"] for name in ("cipos", "ciend", "alt_copy_numbers") if counts[f"n_{name}"]}
    if not_written:
        warnings.append("row(s) with fields VCF Zarr export does not write: " + ", ".join(f"{n} {k}" for k, n in not_written.items()))
    vectors_dropped = rows.filter(pl.col("_has_vector") & ~pl.col("_indexable")).height
    if vectors_dropped:
        warnings.append(
            f"{vectors_dropped} row(s) had GL/GP/AD indexed over another allele list than the written REF and ALTs; "
            "those vectors were left out"
        )
    mixed = rows.filter(pl.col("_mixed_phase")).height
    if mixed:
        warnings.append(f"{mixed} row(s) had phased and unphased slots in one genotype; written unphased")

    partial = path.with_name(path.name + ".partial")
    shutil.rmtree(partial, ignore_errors=True)
    try:
        written = _write_store(partial, rows, facet, sample_name, names, contig_index, block_count, variants_chunk_size)
    except BaseException:
        shutil.rmtree(partial, ignore_errors=True)
        raise
    os.rename(partial, path)
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return VcfZarrWritten(
        path=path,
        rows_written=rows.height,
        rows_skipped=skipped,
        blocks_not_written=block_count,
        arrays=tuple(sorted(written)),
        warnings=tuple(warnings),
    )


def _columns() -> list[pl.Expr]:
    """What the export reads from each row."""
    indexable = _vectors_indexable()
    phased = pl.col("phased").list.slice(1)
    return [
        pl.col("chrom"),
        pl.col("pos"),
        (pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64) + 1).alias("_length"),
        pl.concat_list(_written_ref(), pl.col("alts")).alias("_alleles"),
        pl.col("ids").list.join(";").alias("_id"),
        pl.col("qual"),
        pl.col("filters"),
        pl.col("gt"),
        (phased.list.len().fill_null(0).gt(0) & phased.list.eval(pl.element().fill_null(False)).list.all()).alias("_phased"),
        (phased.list.eval(pl.element().fill_null(False)).list.any() & ~phased.list.eval(pl.element().fill_null(False)).list.all())
        .fill_null(False)
        .alias("_mixed_phase"),
        pl.col("gq"),
        pl.col("dp"),
        pl.when(indexable).then(pl.col("ad")).alias("_AD"),
        pl.when(indexable).then(pl.col("gl")).alias("_GL"),
        pl.when(indexable).then(pl.col("gp")).alias("_GP"),
        indexable.alias("_indexable"),
        pl.any_horizontal(pl.col(name).list.len().fill_null(0) > 0 for name in ("gl", "gp", "ad")).alias("_has_vector"),
        pl.col("dosage"),
        _ps_text().cast(pl.Int64, strict=False).alias("_PS"),
        pl.when(pl.col("sample_filters").list.len() > 0).then(pl.col("sample_filters").list.join(";")).alias("_FT"),
        pl.col("copy_number"),
        pl.col("sv_type"),
        pl.col("split_signal").is_not_null().alias("_split"),
    ]


def _check_order(rows: pl.DataFrame) -> None:
    """Rows grouped by contig, positions ascending, as `load` emits them; `UnsortedRowsError` otherwise."""
    same = pl.col("chrom") == pl.col("chrom").shift(1)
    seen = rows.select(
        (same & (pl.col("pos") < pl.col("pos").shift(1))).fill_null(False).any().alias("back"),
        (~same).fill_null(True).sum().alias("runs"),
        pl.col("chrom").n_unique().alias("contigs"),
    ).row(0, named=True)
    if rows.height and (seen["back"] or seen["runs"] != seen["contigs"]):
        raise UnsortedRowsError("the rows are not grouped by contig with positions ascending")


class _Pad:
    """Padding past a row's own length: the spec's fill value, distinct from missing."""


_PAD = _Pad()


def _write_store(
    root: Path,
    rows: pl.DataFrame,
    facet: dict,
    sample_name: str,
    names: _Names,
    contig_index: dict[str, int],
    block_count: int,
    chunk: int,
) -> list[str]:
    """Every array of the store under `root`; returns their names."""
    root.mkdir(parents=True)
    _write_json(root / ".zgroup", {"zarr_format": 2})
    _write_json(
        root / ".zattrs",
        {
            "vcf_zarr_version": VCF_ZARR_VERSION,
            "source": f"just-vcf-{version('just-vcf')}",
            "vcf_meta_information": [["fileformat", "VCFv4.3"]],
            "just_vcf": {
                "schema_version": SCHEMA_VERSION,
                "facet_id": facet["facet_id"],
                "source_sha256": facet["source_sha256"],
                "format_id": facet["format_id"],
                "derived_from": facet["derived_from"],
                "blocks_not_written": block_count,
            },
        },
    )
    variants = rows.height
    alleles = max(rows["_alleles"].list.len().max() or 1, 1)
    ploidy = max(rows["gt"].list.len().max() or 1, 1)
    gt_max = rows["gt"].list.max().max()
    filters = sorted({name for value in rows["filters"].drop_nulls().to_list() for name in value} - {"PASS"})
    filter_ids = ["PASS", *filters]
    contigs = names.contigs
    written: list[str] = []

    def array_of(name: str, dims: list[str], dtype: str, width: tuple[int, ...], encode, fill, description: str | None = None) -> None:
        """One array over the variants: `encode(part)` gives a chunk's flat values, row-major."""
        _write_meta(root / name, [variants, *width], [chunk, *width], dtype, fill, dims, description)
        per_row = prod(width)
        for index, start in enumerate(range(0, variants, chunk)):
            key = ".".join([str(index), *("0" for _ in width)])
            values = encode(rows.slice(start, chunk))
            # Zarr v2 stores every chunk whole: the last one is padded past the variants.
            values.extend([_PAD] * (chunk * per_row - len(values)))
            _write_chunk(root / name / key, _encode(dtype, values))
        written.append(name)

    def fixed(name: str, dims: list[str], dtype: str, values: list, fill) -> None:
        """A small array of its own dimension, in one chunk."""
        _write_meta(root / name, [len(values)], [max(len(values), 1)], dtype, fill, dims, None)
        if values:
            _write_chunk(root / name / "0", _encode(dtype, values))
        written.append(name)

    fixed("contig_id", ["contigs"], "|O", [contig["name"] for contig in contigs], _STR_FILL)
    lengths = [contig["length"] for contig in contigs]
    if any(length is not None for length in lengths):
        fixed("contig_length", ["contigs"], "<i8", lengths, _INT_FILL)
    fixed("filter_id", ["filters"], "|O", filter_ids, _STR_FILL)
    fixed("sample_id", ["samples"], "|O", [sample_name], _STR_FILL)

    contig_dtype = "|i1" if len(contigs) < 2**7 else "<i2" if len(contigs) < 2**15 else "<i4"
    array_of("variant_contig", ["variants"], contig_dtype, (), lambda part: [contig_index[c] for c in part["chrom"].to_list()], None)
    array_of("variant_position", ["variants"], "<i4", (), lambda part: part["pos"].to_list(), None, "The reference position")
    array_of(
        "variant_length", ["variants"], "<i4", (), lambda part: part["_length"].to_list(), None, "The length of the variant measured in bases"
    )
    array_of(
        "variant_allele",
        ["variants", "alleles"],
        "|O",
        (alleles,),
        lambda part: _padded(part["_alleles"], alleles),
        _STR_FILL,
        "List of the reference and alternate alleles",
    )
    array_of("variant_id", ["variants"], "|O", (), lambda part: [value or None for value in part["_id"].to_list()], None)
    array_of("variant_id_mask", ["variants"], "|b1", (), lambda part: [not value for value in part["_id"].to_list()], None)
    array_of("variant_quality", ["variants"], "<f4", (), lambda part: part["qual"].to_list(), "NaN", "Phred-scaled quality score")
    array_of(
        "variant_filter",
        ["variants", "filters"],
        "|b1",
        (len(filter_ids),),
        lambda part: [name in (value or ()) for value in part["filters"].to_list() for name in filter_ids],
        None,
        "Filter status of the variant",
    )
    if rows["sv_type"].is_not_null().any():
        array_of("variant_SVTYPE", ["variants"], "|O", (), lambda part: part["sv_type"].to_list(), None)
    if rows["_split"].any():
        array_of(f"variant_{SPLIT_INFO_KEY}", ["variants"], "|b1", (), lambda part: part["_split"].to_list(), None)

    gt_dtype = "|i1" if (gt_max or 0) < 2**7 else "<i2"
    array_of("call_genotype", ["variants", "samples", "ploidy"], gt_dtype, (1, ploidy), lambda part: _padded(part["gt"], ploidy), _INT_FILL)
    array_of("call_genotype_phased", ["variants", "samples"], "|b1", (1,), lambda part: part["_phased"].to_list(), None)
    ps_max = rows["_PS"].abs().max() or 0
    ps_dtype = "<i4" if ps_max < 2**31 else "<i8"
    scalars = (("GQ", "gq", "<i4"), ("DP", "dp", "<i4"), ("DS", "dosage", "<f4"), ("PS", "_PS", ps_dtype), ("FT", "_FT", "|O"), ("CN", "copy_number", "<f4"))
    for key, column, dtype in scalars:
        if rows[column].is_not_null().any():
            array_of(
                f"call_{key}",
                ["variants", "samples"],
                dtype,
                (1,),
                lambda part, column=column: part[column].to_list(),
                _fill(dtype),
                _DESCRIPTIONS[key],
            )
    for key, dtype, dim in (("AD", "<i4", "alleles"), ("GL", "<f4", "genotypes"), ("GP", "<f4", "genotypes")):
        width = rows[f"_{key}"].list.len().max() or 0
        if width:
            array_of(
                f"call_{key}",
                ["variants", "samples", dim],
                dtype,
                (1, width),
                lambda part, key=key, width=width: _padded(part[f"_{key}"], width),
                _fill(dtype),
                _DESCRIPTIONS[key],
            )
    return written


def _fill(dtype: str) -> object:
    """The `.zarray` fill value of a dtype: the spec's padding."""
    return _STR_FILL if dtype == "|O" else "NaN" if dtype == "<f4" else _INT_FILL


def _padded(column: pl.Series, width: int) -> list:
    """A list column as flat values, `width` per row: `None` for a missing element (a null list is
    missing throughout), `_PAD` past the row's own length."""
    flat: list = []
    for value in column.to_list():
        if value is None:
            flat.extend([None] * width)
        else:
            flat.extend(value)
            flat.extend([_PAD] * (width - len(value)))
    return flat


def _encode(dtype: str, values: list) -> bytes:
    """Flat values as a chunk's raw bytes, little-endian, with the spec's sentinels in place of `None`
    (missing) and `_PAD` (fill)."""
    if dtype == "|O":
        parts = [struct.pack("<I", len(values))]
        for value in values:
            text = (_STR_FILL if value is _PAD else _STR_MISSING if value is None else value).encode()
            parts.append(struct.pack("<I", len(text)))
            parts.append(text)
        return b"".join(parts)
    if dtype == "|b1":
        return bytes(1 if value is True else 0 for value in values)
    if dtype == "<f4":
        floats = array("f", (0.0 if value is None or value is _PAD else value for value in values))
        encoded = array("I")
        encoded.frombytes(floats.tobytes())
        for index, value in enumerate(values):
            if value is None:
                encoded[index] = _FLOAT_MISSING
            elif value is _PAD:
                encoded[index] = _FLOAT_FILL
    else:
        encoded = array(_TYPECODES[dtype], (_INT_MISSING if value is None else _INT_FILL if value is _PAD else value for value in values))
    if sys.byteorder == "big":
        encoded.byteswap()
    return encoded.tobytes()


def _write_meta(directory: Path, shape: list[int], chunks: list[int], dtype: str, fill: object, dims: list[str], description: str | None) -> None:
    directory.mkdir()
    _write_json(
        directory / ".zarray",
        {
            "zarr_format": 2,
            "shape": shape,
            "chunks": chunks,
            "dtype": dtype,
            "compressor": _COMPRESSOR,
            "fill_value": fill,
            "order": "C",
            "filters": _VLEN if dtype == "|O" else None,
            "dimension_separator": ".",
        },
    )
    _write_json(directory / ".zattrs", {"_ARRAY_DIMENSIONS": dims, **({"description": description} if description else {})})


def _write_chunk(path: Path, raw: bytes) -> None:
    path.write_bytes(gzip.compress(raw, compresslevel=_COMPRESSOR["level"], mtime=0))


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=1, sort_keys=True) + "\n")
