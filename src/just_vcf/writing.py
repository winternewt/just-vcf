"""Write one facet back out as a VCF: the export at the edge (charter, RM6).

`write_vcf(loaded, path)` writes a `Loaded`'s `calls` rows as a single-sample VCF 4.3, BGZF-compressed
by default. The frames stay the product; the file is one way out of them, written by this module's own
text writer rather than polars-bio's `sink_vcf`, because the writer has to own every byte: the GT text
with its per-haplotype phasing, the float formatting that makes two writes of one frame identical
(charter 12), and the declaration of every key it writes (polars-bio drops an undeclared column
silently on write, P13).

What is written, and the rules it follows:

- **REF** is `ref` (read from the assembly sequence) when there is one, else `source_ref` (what the
  file wrote). It is never taken from the genotype (charter 10). A row with neither, a row without
  `alts` (an array row before RM3, or one RM3 withheld in a mixed-strand facet, F13), a row whose
  `ref_check` is `mismatch`, a row with an ALT that is no VCF allele (a `-` written for a deletion,
  F12, until `normalize` anchors it), a row with a genotype that
  is not indexed (`gt` null beside `gt_alleles`) and a row whose `gt` points past its ALTs cannot be
  written as valid VCF. The writer refuses them with `UnwritableRowsError`, counted per reason, unless
  the caller passes `resolved_only=True`; then the rest is written and the counts are reported.
- **FORMAT** is per row: a key is listed where the row has a value for it, so a field the source
  left out stays out. GT is rebuilt from `gt` and `phased` (`source_gt` is null on derived facets and
  VCF 4.4+ files). `gl` is written as GL: the frame holds log10 likelihoods and does not remember
  whether the source wrote PL. `gl`, `gp` and `ad` are written only when `source_alleles`, the order
  they are indexed over, is the written `[REF] + ALTs`; otherwise they are left out and counted.
- **VCF 4.3**, not 4.4+: GP is a probability from 4.3 on, which is what `gp` holds, and the writer
  emits no leading phase indicator, so the first haplotype's phase follows 4.3's implicit rule.
- **INFO** carries `END` where `end_pos` is not the REF span (symbolic alleles), and
  `JUST_VCF_SPLIT` on rows flagged as one record of a split multi-allelic site (`split_signal`:
  `normalize`'s split rows, and rows RM270 flagged). Those rows are written as the frame holds them,
  with a warning and the tag, so a reader is told that a `0` there may mean another ALT. `load` reads
  the tag back as a marker (`split_signal = marker`), so a written normalized facet keeps its flag.
  `rejoin_split` joins them first, if that is what the caller wants.
- **Blocks** are not written. A gVCF block record needs the reference base at its start, which the
  `blocks` frame does not hold until RM3 reads the sequence; the count is reported.
- **Contig names** are written as stored (the canonical `chrom`), with the `chr` prefix (`MT` becomes
  `chrM`), or as the source spelled them (`contig_map`). The header's `##contig` lines follow the same
  names and the rows' order.

- **A tabix index** (`<path>.tbi`) is written beside BGZF output by default, from the byte offsets
  the writer itself produced (`just_vcf.tabix`), so the file is region-queryable by htslib without a
  second pass over it. Each record's interval is the one htslib's VCF preset reads back: from `POS`,
  to INFO `END` where the line writes one, else over REF. A plain-text write has no index, and a
  write without one removes a stale `<path>.tbi`.

The header names the facet (`facet_id`, `source_sha256`, `format_id`, the parent of a derived facet)
and this package's version. It never names `source_path` or the source's `##SAMPLE` lines, which can
carry a real name; the sample column is the facet's `sample`, as given.
"""

import hashlib
import io
import logging
import os
from dataclasses import dataclass
from enum import Enum
from importlib.metadata import version
from pathlib import Path
from typing import BinaryIO

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.bgzf import BgzfWriter
from just_vcf.loaders.vcf import FORMAT_KEYS, SPLIT_INFO_KEY
from just_vcf.loading import Loaded
from just_vcf.schema import CONTIG_RANK, SCHEMA_VERSION, RefCheck, canonical_chrom
from just_vcf.tabix import TabixIndexer

logger = logging.getLogger(__name__)

VCF_VERSION = "VCFv4.3"

_FORMAT_DESCRIPTIONS = {
    "GT": "Genotype",
    "GQ": "Conditional genotype quality",
    "DP": "Read depth of this sample",
    "AD": "Read depth for each allele",
    "GL": "Genotype likelihoods, log10",
    "GP": "Genotype posterior probabilities",
    "DS": "Expected count of the alternate allele",
    "FT": "Sample filter",
    "PS": "Phase set",
}
"""The FORMAT keys the writer can emit, in the order they are written. Number and Type come from the
loader's `FORMAT_KEYS`, so a written file declares each key the way the loader reads it back."""

_NO_STRINGS = pl.lit([], dtype=pl.List(pl.String))

_UNWRITABLE_REASONS = (
    "ref_unresolved",
    "ref_mismatch",
    "alts_unresolved",
    "alt_not_an_allele",
    "genotype_unindexed",
    "genotype_index_past_alts",
)

_VCF_ALLELE = r"^([ACGTNacgtn]+|\*|<[^<>]+>|\.[ACGTNacgtn]+|[ACGTNacgtn]+\.)$|\[|\]"
"""An ALT VCF can write: bases, `*`, a symbolic allele, a single breakend (`.A`, `A.`) or a breakend
(`A[2:500[`). Anything else (F12's `-` for a deletion) is not one, and htslib refuses the file."""


class ContigNames(str, Enum):
    """How contig names are spelled in a written VCF."""

    STORED = "stored"
    """The canonical `chrom`: `1`–`22`, `X`, `Y`, `MT`, `XY`; other contigs as the source spelled them."""
    CHR = "chr"
    """`chr1`–`chr22`, `chrX`, `chrY`, `chrM`; every other contig as stored."""
    SOURCE = "source"
    """As the source file spelled them, from the facet's `contig_map`."""


_CHR_NAMES = {**{str(n): f"chr{n}" for n in range(1, 23)}, "X": "chrX", "Y": "chrY", "MT": "chrM"}
"""Canonical → `chr`-prefixed, listed name by name: `MT` is `chrM`, and nothing is prefixed by pattern."""


class UnwritableRowsError(ValueError):
    """Rows that cannot be written as valid VCF, counted per reason. `resolved_only=True` skips them."""

    def __init__(self, facet_id: str, counts: dict[str, int]) -> None:
        self.counts = counts
        detail = ", ".join(f"{count} {reason}" for reason, count in counts.items())
        super().__init__(
            f"{facet_id}: {sum(counts.values())} row(s) cannot be written as VCF ({detail}); "
            "pass resolved_only=True to write the rest"
        )


class UnsortedRowsError(ValueError):
    """The rows are not grouped by contig with positions ascending, which a VCF requires."""


class VcfWritten(BaseModel):
    """What `write_vcf` wrote."""

    model_config = ConfigDict(frozen=True)

    path: Path
    sha256: str
    """Digest of the bytes written, compressed or not."""
    rows_written: int
    rows_skipped: dict[str, int]
    """Rows left out under `resolved_only=True`, per reason. Empty when nothing was left out."""
    blocks_not_written: int
    warnings: tuple[str, ...] = ()
    index_path: Path | None = None
    """The tabix index written beside a BGZF file (`<path>.tbi`); None when none was written."""


class _Hashing:
    """A binary sink that digests what passes through it."""

    def __init__(self, raw: BinaryIO) -> None:
        self._raw = raw
        self.digest = hashlib.sha256()

    def write(self, data: bytes) -> int:
        self.digest.update(data)
        return self._raw.write(data)


def write_vcf(
    loaded: Loaded,
    path: Path | str,
    *,
    contig_names: ContigNames = ContigNames.STORED,
    sample: str | None = None,
    resolved_only: bool = False,
    compress: bool = True,
    index: bool | None = None,
) -> VcfWritten:
    """Write `loaded`'s one facet to `path` as a VCF 4.3 (BGZF unless `compress=False`).

    `sample` names the sample column; the default is the facet's `sample`, then its `facet_id` for a
    source that names none (an array export). Rows that cannot be valid VCF raise
    `UnwritableRowsError` unless `resolved_only=True`. The rows must be in the order `load` emits
    them (grouped by contig, positions ascending); otherwise `UnsortedRowsError`. The same frames and
    options give the same bytes.

    `index` writes a tabix index, `<path>.tbi`, beside the file; the default is to write one whenever
    the output is BGZF. A plain-text file cannot be indexed (`index=True` with `compress=False` is
    refused), and a write without an index removes a `<path>.tbi` left from an earlier write, which
    would no longer describe the file.
    """
    path = Path(path)
    if index and not compress:
        raise ValueError("a tabix index needs BGZF output; pass compress=True or leave index unset")
    indexed = compress if index is None else index
    index_path = path.with_name(path.name + ".tbi")
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"write_vcf takes one facet; this Loaded has {facets.height}")
    facet = facets.row(0, named=True)
    facet_id = facet["facet_id"]
    sample_name = sample if sample is not None else facet["sample"] if facet["sample"] is not None else facet_id
    if not sample_name or any(char in sample_name for char in "\t\n\r"):
        raise ValueError("the sample name must be non-empty and hold no tab or line break")

    calls = loaded.calls.with_columns(_unwritable().alias("_unwritable"))
    stats = _scan(calls)
    skipped = {reason: stats[f"n_{reason}"] for reason in _UNWRITABLE_REASONS if stats[f"n_{reason}"]}
    if skipped and not resolved_only:
        raise UnwritableRowsError(facet_id, skipped)
    blocks = loaded.blocks.select(pl.len()).collect().item()

    names = _contig_names(contig_names, facet, stats["chroms"])
    warnings: list[str] = []
    if skipped:
        warnings.append(f"{sum(skipped.values())} row(s) not written: " + ", ".join(f"{n} {r}" for r, n in skipped.items()))
    if stats["n_vectors_dropped"]:
        warnings.append(
            f"{stats['n_vectors_dropped']} row(s) had GL/GP/AD indexed over another allele list than the written "
            "REF and ALTs; those vectors were left out"
        )
    if stats["n_ps_dropped"]:
        warnings.append(
            f"{stats['n_ps_dropped']} row(s) had a phase set VCF 4.3 cannot write as PS (a PSL name, or not an "
            "integer); left out"
        )
    if stats["n_split"]:
        warnings.append(
            f"{stats['n_split']} row(s) look split one row per ALT (RM270) and were written as held, tagged "
            f"INFO/{SPLIT_INFO_KEY}: a `0` there may mean another ALT, not REF; rejoin_split joins them first"
        )
    if blocks:
        warnings.append(
            f"{blocks} reference/no-call block(s) not written: a gVCF block needs the REF base at its start, "
            "which the blocks frame does not hold before RM3"
        )

    declared = [key for key in _FORMAT_DESCRIPTIONS if stats[f"has_{key}"]]
    header = _header(facet, sample_name, names, declared, stats)
    rows = calls.filter(pl.col("_unwritable").is_null()) if skipped else calls
    lines = rows.select(
        _line(names.mapping, stats["ploidy"]).alias("line"),
        "chrom",
        "pos",
        *_index_columns(names.mapping),
    )

    # Written beside the target and moved over it at the end, so a refused or failed write leaves no
    # half-written VCF (or index) behind under the requested name.
    partial = path.with_name(path.name + ".partial")
    index_partial = index_path.with_name(index_path.name + ".partial")
    indexer = TabixIndexer() if indexed else None
    try:
        with partial.open("wb") as raw:
            sink = _Hashing(raw)
            if compress:
                with BgzfWriter(sink, level=6) as out:
                    written = _write_body(out, header, lines, indexer)
                if indexer is not None:
                    with index_partial.open("wb") as raw_index, BgzfWriter(raw_index, level=6) as index_out:
                        index_out.write(indexer.build(out.virtual_offset))
            else:
                written = _write_body(sink, header, lines, None)
    except BaseException:
        partial.unlink(missing_ok=True)
        index_partial.unlink(missing_ok=True)
        raise
    os.replace(partial, path)
    if indexer is not None:
        os.replace(index_partial, index_path)
    else:
        index_path.unlink(missing_ok=True)

    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return VcfWritten(
        path=path,
        sha256=sink.digest.hexdigest(),
        rows_written=written,
        rows_skipped=skipped,
        blocks_not_written=blocks,
        warnings=tuple(warnings),
        index_path=index_path if indexer is not None else None,
    )


def _write_body(out: BinaryIO | BgzfWriter | _Hashing, header: str, lines: pl.LazyFrame, indexer: TabixIndexer | None) -> int:
    """The header, then every line batch by batch, checking the order as it goes and handing each
    record's interval and byte offsets to `indexer`. Returns the row count."""
    head = header.encode()
    out.write(head)
    offset = len(head)
    order = _OrderCheck()
    written = 0
    for batch in lines.collect_batches():
        order.check(batch)
        buffer = io.BytesIO()
        batch.select("line").write_csv(buffer, include_header=False, quote_style="never", line_terminator="\n")
        data = buffer.getvalue()
        out.write(data)
        if indexer is not None:
            size = pl.col("line").str.len_bytes().cast(pl.Int64) + 1
            located = batch.select(
                pl.col("_ix_name").alias("name"),
                pl.col("_ix_beg").alias("beg"),
                pl.col("_ix_end").alias("end"),
                (offset + size.cum_sum()).alias("u_end"),
            ).with_columns((pl.col("u_end") - batch.select(size).to_series()).alias("u_start"))
            if batch.height and located["u_end"][-1] != offset + len(data):
                raise RuntimeError("the written lines' sizes do not add up to the bytes written; the index would be wrong")
            indexer.add(located)
        offset += len(data)
        written += batch.height
    return written


def _index_columns(mapping: dict[str, str]) -> list[pl.Expr]:
    """Each record's contig as written and its interval as htslib's VCF preset reads it back: 0-based
    from `POS - 1`, to INFO `END` where the line writes one, else to `POS - 1 + len(REF)`."""
    ref_span_end = pl.col("pos").cast(pl.Int64) + _written_ref().str.len_chars().cast(pl.Int64) - 1
    end = pl.col("end_pos").cast(pl.Int64)
    return [
        pl.col("chrom").replace(mapping).alias("_ix_name"),
        (pl.col("pos").cast(pl.Int64) - 1).alias("_ix_beg"),
        pl.when((end != ref_span_end).fill_null(False)).then(end).otherwise(ref_span_end).alias("_ix_end"),
    ]


def _written_ref() -> pl.Expr:
    """REF as written: `ref` from the sequence when known, else what the file wrote."""
    return pl.coalesce(pl.col("ref"), pl.col("source_ref"))


def _unwritable() -> pl.Expr:
    """Why a row cannot be written as valid VCF, first reason first; null when it can."""
    return (
        pl.when(_written_ref().is_null())
        .then(pl.lit("ref_unresolved"))
        .when(pl.col("ref_check") == RefCheck.MISMATCH.value)
        .then(pl.lit("ref_mismatch"))
        .when(pl.col("alts").is_null())
        .then(pl.lit("alts_unresolved"))
        .when(~pl.col("alts").list.eval(pl.element().str.contains(_VCF_ALLELE)).list.all())
        .then(pl.lit("alt_not_an_allele"))
        .when(pl.col("gt").is_null() & pl.col("gt_alleles").is_not_null())
        .then(pl.lit("genotype_unindexed"))
        .when(pl.col("gt").list.max().cast(pl.Int64) > pl.col("alts").list.len().cast(pl.Int64))
        .then(pl.lit("genotype_index_past_alts"))
        .otherwise(None)
    )


def _vectors_indexable() -> pl.Expr:
    """True where `gl`, `gp` and `ad` are indexed over the written `[REF] + ALTs`."""
    return (pl.col("source_alleles") == pl.concat_list(_written_ref(), pl.col("alts"))).fill_null(False)


def _has_vector() -> pl.Expr:
    return pl.any_horizontal(pl.col(name).list.len().fill_null(0) > 0 for name in ("gl", "gp", "ad"))


def _ps_text() -> pl.Expr:
    """PS: the row's PS-kind phase set, as an integer. Null where there is none or it is a PSL name."""
    value = pl.col("phase_sets").list.drop_nulls().list.first()
    psl = pl.col("phase_set_kinds").list.contains("psl").fill_null(False)
    return pl.when(psl).then(None).otherwise(value.cast(pl.Int64, strict=False)).cast(pl.String)


def _number_text(value: pl.Expr) -> pl.Expr:
    """A number as text: the shortest form that reads back to the same value, and `0.0` for `-0.0`
    (a PL of 0 becomes a GL of -0.0, which is the same likelihood)."""
    return pl.when(value == 0).then(value.abs()).otherwise(value).cast(pl.String)


def _list_text(name: str) -> pl.Expr:
    """A list column as VCF's comma-separated text, `.` per missing element; null when empty or null."""
    column = pl.col(name)
    text = column.list.eval(_number_text(pl.element()).fill_null(".")).list.join(",")
    return pl.when(column.list.len() > 0).then(text).otherwise(None)


def _gt_text(ploidy: int) -> pl.Expr:
    """GT from `gt` and `phased`: `.` per uncalled haplotype, `|` before a phased one, else `/`."""
    gt = pl.col("gt")
    parts = [gt.list.get(0, null_on_oob=True).cast(pl.String).fill_null(".")]
    for slot in range(1, ploidy):
        separator = pl.when(pl.col("phased").list.get(slot, null_on_oob=True).fill_null(False)).then(pl.lit("|")).otherwise(pl.lit("/"))
        allele = gt.list.get(slot, null_on_oob=True).cast(pl.String).fill_null(".")
        parts.append(pl.when(gt.list.len() > slot).then(pl.concat_str(separator, allele)))
    return pl.when(gt.list.len() > 0).then(pl.concat_str(parts, ignore_nulls=True))


def _format_values(ploidy: int) -> dict[str, pl.Expr]:
    """Each FORMAT key's text per row, null where the row has no value for it."""
    vectors = _vectors_indexable()
    return {
        "GT": _gt_text(ploidy),
        "GQ": pl.col("gq").cast(pl.String),
        "DP": pl.col("dp").cast(pl.String),
        "AD": pl.when(vectors).then(_list_text("ad")),
        "GL": pl.when(vectors).then(_list_text("gl")),
        "GP": pl.when(vectors).then(_list_text("gp")),
        "DS": _number_text(pl.col("dosage")),
        "FT": pl.when(pl.col("sample_filters").list.len() > 0).then(pl.col("sample_filters").list.join(";")),
        "PS": _ps_text(),
    }


def _scan(calls: pl.LazyFrame) -> dict:
    """One pass for the header and the warnings: which keys occur, the filters, the contigs, the counts."""
    ploidy = calls.select(pl.col("gt").list.len().max()).collect().item() or 1
    values = _format_values(ploidy)
    ref_span_end = pl.col("pos").cast(pl.Int64) + _written_ref().str.len_chars().cast(pl.Int64) - 1
    row = (
        calls.select(
            *[(pl.col("_unwritable") == reason).sum().alias(f"n_{reason}") for reason in _UNWRITABLE_REASONS],
            *[text.is_not_null().any().alias(f"has_{key}") for key, text in values.items()],
            (_has_vector() & ~_vectors_indexable() & pl.col("_unwritable").is_null()).sum().alias("n_vectors_dropped"),
            (pl.col("phase_sets").list.drop_nulls().list.len().fill_null(0).gt(0) & _ps_text().is_null()).sum().alias("n_ps_dropped"),
            pl.col("split_signal").is_not_null().sum().alias("n_split"),
            (pl.col("end_pos").cast(pl.Int64) != ref_span_end).fill_null(False).any().alias("has_END"),
            pl.concat_list(pl.col("filters").fill_null(_NO_STRINGS), pl.col("sample_filters").fill_null(_NO_STRINGS))
            .explode(empty_as_null=False)
            .drop_nulls()
            .unique()
            .sort()
            .implode()
            .alias("filters"),
            pl.col("chrom").unique().sort().implode().alias("chroms"),
        )
        .collect()
        .row(0, named=True)
    )
    row["ploidy"] = ploidy
    return row


@dataclass(frozen=True)
class _Names:
    """Canonical `chrom` → written name, and the `##contig` lines in the rows' order."""

    mapping: dict[str, str]
    contigs: list[dict]


def _contig_names(style: ContigNames, facet: dict, chroms: list[str]) -> _Names:
    """The written name of every contig in the rows and in the facet's `##contig` lines."""
    if style is ContigNames.SOURCE:
        # `contig_map` holds the contigs that have rows; a `##contig` with none is spelled by the header.
        by_chrom: dict[str, set[str]] = {}
        for entry in facet["contig_map"] or []:
            by_chrom.setdefault(entry["chrom"], set()).add(entry["source"])
        for contig in facet["header_contigs"] or []:
            by_chrom.setdefault(canonical_chrom(contig["name"]), set()).add(contig["name"])
        clashes = sorted(chrom for chrom, sources in by_chrom.items() if len(sources) > 1)
        if clashes:
            raise ValueError(f"the source spelled {', '.join(clashes)} more than one way; its spelling cannot be restored")
        source = {chrom: sources.pop() for chrom, sources in by_chrom.items()}

    def written(chrom: str) -> str:
        if style is ContigNames.CHR:
            return _CHR_NAMES.get(chrom, chrom)
        if style is ContigNames.SOURCE:
            return source.get(chrom, chrom)
        return chrom

    declared: dict[str, dict] = {}
    for contig in facet["header_contigs"] or []:
        chrom = canonical_chrom(contig["name"])
        declared.setdefault(chrom, {"length": contig["length"], "md5": contig["md5"]})
    for chrom in chroms:
        declared.setdefault(chrom, {"length": None, "md5": None})
    ordered = sorted(declared, key=lambda chrom: (CONTIG_RANK.get(chrom, len(CONTIG_RANK)), chrom))
    return _Names(
        mapping={chrom: written(chrom) for chrom in ordered},
        contigs=[{"name": written(chrom), **declared[chrom]} for chrom in ordered],
    )


def _header(facet: dict, sample: str, names: _Names, declared: list[str], stats: dict) -> str:
    lines = [
        f"##fileformat={VCF_VERSION}",
        f"##source=just-vcf {version('just-vcf')}",
        f"##just_vcf_schema_version={SCHEMA_VERSION}",
        f"##just_vcf_facet_id={facet['facet_id']}",
    ]
    for key in ("source_sha256", "format_id", "derived_from", "relation", "assembly_label"):
        if facet[key] is not None:
            lines.append(f"##just_vcf_{key}={facet[key]}")
    filters = stats["filters"] or []
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
    lines += [
        f'##FILTER=<ID={name},Description="Carried from the source; its description was not kept">'
        for name in filters
        if name != "PASS"
    ]
    if stats["has_END"]:
        lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the record">')
    if stats["n_split"]:
        lines.append(
            f'##INFO=<ID={SPLIT_INFO_KEY},Number=1,Type=String,Description="This row looks like one record of a '
            'multi-allelic site written one row per ALT (just-vcf RM270): a 0 in GT may mean another ALT, not REF">'
        )
    for key in declared:
        number, kind = FORMAT_KEYS[key]
        lines.append(f'##FORMAT=<ID={key},Number={number},Type={kind},Description="{_FORMAT_DESCRIPTIONS[key]}">')
    for contig in names.contigs:
        attrs = f"ID={contig['name']}"
        if contig["length"] is not None:
            attrs += f",length={contig['length']}"
        if contig["md5"]:
            attrs += f",md5={contig['md5']}"
        lines.append(f"##contig=<{attrs}>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]))
    return "\n".join(lines) + "\n"


def _line(mapping: dict[str, str], ploidy: int) -> pl.Expr:
    """One VCF data line per row."""
    ref = _written_ref()
    alts = pl.col("alts")
    ref_span_end = pl.col("pos").cast(pl.Int64) + ref.str.len_chars().cast(pl.Int64) - 1
    end = pl.when(pl.col("end_pos").cast(pl.Int64) != ref_span_end).then(pl.concat_str(pl.lit("END="), pl.col("end_pos").cast(pl.String)))
    split = pl.when(pl.col("split_signal").is_not_null()).then(
        pl.concat_str(pl.lit(f"{SPLIT_INFO_KEY}="), pl.col("split_signal").cast(pl.String))
    )
    info = pl.concat_str([end, split], separator=";", ignore_nulls=True)
    values = _format_values(ploidy)
    keys = pl.concat_str([pl.when(text.is_not_null()).then(pl.lit(key)) for key, text in values.items()], separator=":", ignore_nulls=True)
    sample = pl.concat_str(list(values.values()), separator=":", ignore_nulls=True)
    # A row with no FORMAT value at all still needs both columns: GT, uncalled.
    nothing = keys == ""
    return pl.concat_str(
        [
            pl.col("chrom").replace(mapping),
            pl.col("pos").cast(pl.String),
            pl.when(pl.col("ids").list.len() > 0).then(pl.col("ids").list.join(";")).otherwise(pl.lit(".")),
            ref,
            pl.when(alts.list.len() > 0).then(alts.list.join(",")).otherwise(pl.lit(".")),
            _number_text(pl.col("qual")).fill_null("."),
            pl.when(pl.col("filters").list.len() > 0).then(pl.col("filters").list.join(";")).otherwise(pl.lit(".")),
            pl.when(info == "").then(pl.lit(".")).otherwise(info),
            pl.when(nothing).then(pl.lit("GT")).otherwise(keys),
            pl.when(nothing).then(pl.lit(".")).otherwise(sample),
        ],
        separator="\t",
    )


class _OrderCheck:
    """Raises `UnsortedRowsError` when a contig comes back after another, or a position goes down."""

    def __init__(self) -> None:
        self._seen: set[str] = set()
        self._last: tuple[str, int] | None = None

    def check(self, batch: pl.DataFrame) -> None:
        if batch.is_empty():
            return
        chroms = batch["chrom"]
        positions = batch["pos"].cast(pl.Int64)
        starts = (chroms != chroms.shift(1)).fill_null(True)
        if ((positions.diff().fill_null(0) < 0) & ~starts).any():
            raise UnsortedRowsError("positions go down within a contig; write the frames in the order load emits them")
        for index in starts.arg_true().to_list():
            chrom, pos = chroms[index], positions[index]
            if index == 0 and self._last is not None and chrom == self._last[0]:
                if pos < self._last[1]:
                    raise UnsortedRowsError(f"positions go down within contig {chrom}")
                continue
            if chrom in self._seen:
                raise UnsortedRowsError(f"contig {chrom} comes back after another; the rows are not grouped by contig")
            self._seen.add(chrom)
        self._last = (chroms[-1], positions[-1])
