"""Coverage as a facet: what was callable, read from depth summaries or computed over alignments (RM278).

A coverage facet says how deep the reads were, interval by interval. It is evidence for absence: a
whole-genome callset that wrote no variant at a position the reads covered well says more than one
that wrote none where no read landed. It is read two ways:

- **`read_coverage(path)`** reads what a depth tool already wrote: mosdepth's `.per-base.bed.gz`,
  `.regions.bed.gz`, `.quantized.bed.gz` and `.mosdepth.summary.txt`, and `samtools depth` and
  `samtools coverage` output. The kind is told from the content, never from the file name, and
  compression from the magic bytes.
- **`alignment_coverage(path)`** runs polars-bio's `depth()` over a BAM or CRAM (charter: alignment
  processing is the substrate's, used as is; just-vcf walks no CIGAR of its own) and folds its
  per-base depth into fixed windows while it streams, so the result is a few million rows for a
  human genome whatever the depth.

Both return a `Coverage`: a `CoverageFacet` (what was read and how) and `intervals`, a lazy frame in
the shape `COVERAGE_INTERVALS`. A position no interval covers is unknown, never zero: a file can be a
region cut, and `samtools depth` leaves zero-depth bases out unless asked not to. The one source whose
silence is measured zero is polars-bio's depth over a BAM, whose header names every contig; its
windows tile those contigs with the zeros filled in.

`callable_sites` and `callable_at` turn a coverage facet into evidence per position: `True` when the
interval's depth bounds show at least `min_depth` reads at every base of it, `False` when they show
fewer at every base, `None` with a reason otherwise. A mean over an interval bounds no single base, so
a regions file or a contig summary answers `None` for every position: what was callable cannot be
read from an average.

polars-bio 0.36's `depth()` writes its depth as Int16 by a wrapping cast from a 32-bit count, so a
base deeper than 32,767 reads comes out negative (and one deeper than 65,535 as a small positive
number, which cannot be told apart). Negative depths are left out of every window, counted in
`unreadable_bases` and in a warning. Its accumulator holds a 32-bit array the length of every contig
it touched for the whole file before it emits anything, about 4 bytes per base of the genome, so the
peak is bounded by the genome, not by a window: 13.3 GB on a 27 GB GRCh38 BAM (docs/ROADMAP.md
RM278). `depth()` takes no region and drops pushed-down filters, so a smaller bound cannot be set
from here. Its run-length block mode is not used: it builds every contig's batches before emitting
the first, so its peak also grows with the reads, and it drops the depth of a contig whose reads
reach its last base.
"""

import gzip
import hashlib
import io
import json
import logging
import re
from collections.abc import Iterator
from enum import Enum
from importlib.metadata import version
from pathlib import Path

import polars as pl
import polars_bio as pb
from pydantic import BaseModel, ConfigDict

from just_vcf import cache
from just_vcf.assembly import infer_assembly
from just_vcf.cache import cache_dir
from just_vcf.schema import CONTIG_RANK, AssemblyLabelBasis, canonical_chrom, canonical_chrom_expr

logger = logging.getLogger(__name__)

COVERAGE_VERSION = 1
"""Bumped when the same input would give different intervals; part of every `facet_id` and of the
cache key of `alignment_coverage`'s spill."""

POLARS_BIO_VERSION = version("polars-bio")
"""The polars-bio whose `depth()` computed a spill: part of its cache key and its manifest, so a
spill one release wrote is never served under another (a fix to the Int16 wrap changes the windows
without a `COVERAGE_VERSION` bump here)."""

DEFAULT_WINDOW = 1_000
"""Window size for per-base sources folded into windows: 3.1 million rows for a human genome."""

DEFAULT_FILTER_FLAG = 1796
"""Reads `alignment_coverage` leaves out: unmapped, secondary, QC-failed, duplicate. mosdepth's and
polars-bio's default, so a depth computed here and a mosdepth file read here count the same reads."""

DEFAULT_MIN_MAPPING_QUALITY = 0
"""mosdepth's and polars-bio's default: every mapping quality counts."""

COVERAGE_INTERVALS: dict[str, pl.DataType] = {
    "chrom": pl.String(),
    "pos": pl.UInt32(),
    "end_pos": pl.UInt32(),
    "depth_mean": pl.Float32(),
    "depth_min": pl.UInt32(),
    "depth_max": pl.UInt32(),
    "bases": pl.UInt32(),
    "label": pl.String(),
}
"""The `intervals` frame of a `Coverage`. `pos` and `end_pos` are 1-based and inclusive, as in
`calls` and `blocks`; `chrom` is canonical (Decision 7).

- `depth_mean`: mean depth over the interval's described bases; null where the source gives none (a
  quantized bin).
- `depth_min`, `depth_max`: the lowest and highest depth at any described base. Null where the source
  does not bound them (a regions mean); `depth_max` is also null for an open-ended bin (`4:inf`).
- `bases`: how many of the interval's bases the source described. Fewer than the interval's length
  means a window only partly described (a `samtools depth` file without zero rows, a wrapped depth
  left out), and its bounds then speak only for the described bases.
- `label`: a quantized bin's or a named region's label, verbatim; null otherwise.
"""


class CoverageFormat(str, Enum):
    """What a coverage facet was read from. Told by shape, so a name covers every tool writing it."""

    PER_BASE_BED = "per_base_bed"
    """BED4 runs of one depth: mosdepth `.per-base.bed.gz`, bedtools `genomecov -bg`."""
    REGIONS_BED = "regions_bed"
    """BED4 or BED5 with a mean depth last: mosdepth `.regions.bed.gz` (`--by`)."""
    QUANTIZED_BED = "quantized_bed"
    """BED4 runs of a depth bin: mosdepth `.quantized.bed.gz` (`--quantize`)."""
    MOSDEPTH_SUMMARY = "mosdepth_summary"
    """mosdepth `.mosdepth.summary.txt`: length, bases, mean, min and max per contig."""
    SAMTOOLS_DEPTH = "samtools_depth"
    """`samtools depth` of one alignment file: contig, 1-based position, depth."""
    SAMTOOLS_COVERAGE = "samtools_coverage"
    """`samtools coverage`: one row per contig or region, with its mean depth."""
    ALIGNMENT_DEPTH = "alignment_depth"
    """polars-bio `depth()` over a BAM or CRAM, folded into windows here."""


class Resolution(str, Enum):
    """What one interval states about the bases inside it."""

    BASE = "base"
    """Every interval has one depth, exact at each of its bases."""
    WINDOW = "window"
    """Fixed windows with the mean, min and max over each, folded here from per-base depth."""
    BIN = "bin"
    """Runs whose depth falls in one bin: bounds, no mean."""
    REGION = "region"
    """A mean over each interval, no bounds; mosdepth's summary adds the bounds per contig."""


class CoverageFacet(BaseModel):
    """What a coverage facet was read from and how. just-vcf's own shape (charter 2)."""

    model_config = ConfigDict(frozen=True)

    facet_id: str
    """Deterministic: the source's digest (for an alignment file, its header, index and size), the
    format, the options and `COVERAGE_VERSION`."""
    format_id: CoverageFormat
    source_path: str
    """A label, not identity."""
    source_sha256: str | None
    """The source bytes' sha256. Null for an alignment file: hashing tens of gigabytes to name what
    `depth()` read would double the cost; `facet_id` uses its header, index and size instead."""
    assembly_label: str | None
    assembly_label_basis: AssemblyLabelBasis | None
    """`stated` when the caller passed it, `inferred` from contig lengths (an alignment header,
    mosdepth's summary). Null when nothing says."""
    contig_lengths: list[dict] | None
    """`{source, chrom, length}` per contig, where the source states lengths; else null."""
    resolution: Resolution
    window: int | None
    """The window size where intervals were folded into windows here."""
    options: dict[str, int]
    """The depth options `alignment_coverage` ran with (`filter_flag`, `min_mapping_quality`); empty
    for a file, whose options the file does not record."""
    zero_filled: bool
    """True when a base inside a listed contig and outside every source interval was measured as
    zero depth and is written so (alignment depth). False: such a base is unknown."""
    dropped_fields: list[str]
    """Source columns no interval column carries (`samtools coverage`'s `numreads`, …)."""
    warnings: list[str]
    subject: str | None = None
    """The caller's label for the person whose reads these are (`subject=`), as `load` writes
    `facets.subject` (RM273): supplied, never inferred. Not identity: `facet_id` does not digest it, and
    `alignment_coverage`'s cache spill does not key on it. `Person.of(..., coverage=...)` groups by it."""


class Coverage(BaseModel):
    """A coverage facet and its intervals (`COVERAGE_INTERVALS`), sorted by contig rank and `pos`."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    facet: CoverageFacet
    intervals: pl.LazyFrame


class UnsupportedCoverageError(ValueError):
    """The file is not a coverage output this module reads, or reading it would mean guessing.

    The message names the kind of input, never its path: a path is a label and can carry a person's
    name (as `DamagedBgzfError`)."""


# --- Reading depth tools' outputs --------------------------------------------------------------

_CHUNK_BYTES = 64 << 20

_MOSDEPTH_SUMMARY_HEADER = ["chrom", "length", "bases", "mean", "min", "max"]
_SAMTOOLS_COVERAGE_HEADER = ["#rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
_BIN_LABEL = re.compile(r"^(\d+):(\d+|inf)$")


def read_coverage(
    path: Path | str,
    *,
    window: int | None = DEFAULT_WINDOW,
    assembly_label: str | None = None,
    kind: CoverageFormat | None = None,
    subject: str | None = None,
) -> Coverage:
    """Read a depth tool's output into a coverage facet (module docstring for the formats).

    `window` folds per-base sources (`per_base_bed`, `samtools_depth`) into fixed windows while they
    stream, so a whole genome's runs never sit in memory; None keeps each source row as an interval,
    read whole, for small or targeted files. The other formats keep their own intervals.
    `assembly_label` states the build where the file cannot (no format but mosdepth's summary gives
    contig lengths); a stated label disagreeing with the lengths is kept, with a warning.
    `kind` states the format where the content cannot tell it (a headerless three-column file whose
    rows read as intervals is refused without it).
    `subject` is the caller's label for the person the reads are of (`CoverageFacet.subject`).
    """
    path = Path(path)
    _check_subject(subject)
    if window is not None and window < 1:
        raise ValueError("window must be a positive number of bases, or None")
    if kind is CoverageFormat.ALIGNMENT_DEPTH:
        raise ValueError("alignment depth is computed by alignment_coverage, not read from a file")
    lines = _head_lines(path)
    kind = kind or _coverage_kind(lines)
    sha256 = _sha256(path)
    warnings: list[str] = []
    dropped: list[str] = []
    lengths: list[dict] | None = None
    folded = window if kind in (CoverageFormat.PER_BASE_BED, CoverageFormat.SAMTOOLS_DEPTH) else None

    if kind is CoverageFormat.MOSDEPTH_SUMMARY:
        intervals, lengths, skipped = _read_summary(path)
        if skipped:
            warnings.append(f"{skipped:,} summary row(s) for the total or for --by regions carry no interval and were left out")
        resolution = Resolution.REGION
    elif kind is CoverageFormat.SAMTOOLS_COVERAGE:
        intervals = _read_samtools_coverage(path)
        dropped = ["numreads", "covbases", "coverage", "meanbaseq", "meanmapq"]
        resolution = Resolution.REGION
    elif kind is CoverageFormat.REGIONS_BED:
        intervals = _read_regions(path, named=len(lines[0].split("\t")) >= 5)
        resolution = Resolution.REGION
    elif kind is CoverageFormat.QUANTIZED_BED:
        intervals, unbounded = _read_quantized(path)
        if unbounded:
            warnings.append(f"{unbounded:,} quantized interval(s) carry a bin label that names no depth range (a renamed bin); their depth is unknown")
        resolution = Resolution.BIN
    elif folded is not None:
        partial = _fold_chunks(_depth_chunks(path, kind), folded)
        intervals = _finish_windows(partial, folded, lengths=None, zero_filled=False)
        resolution = Resolution.WINDOW
    else:
        intervals = _read_runs(path, kind)
        resolution = Resolution.BASE

    label, basis, label_warnings = _assembly(assembly_label, lengths)
    warnings += label_warnings
    for line in warnings:
        logger.warning("%s: %s", kind.value, line)
    options: dict[str, int] = {}
    facet = CoverageFacet(
        facet_id=_facet_id({"sha256": sha256, "format": kind.value, "window": folded, "assembly_label": assembly_label}),
        format_id=kind,
        source_path=str(path),
        source_sha256=sha256,
        assembly_label=label,
        assembly_label_basis=basis,
        contig_lengths=lengths,
        resolution=resolution,
        window=folded,
        options=options,
        zero_filled=False,
        dropped_fields=dropped,
        warnings=warnings,
        subject=subject,
    )
    return Coverage(facet=facet, intervals=_ordered(intervals).lazy())


def _check_subject(subject: str | None) -> None:
    """`load`'s rule for the label: a non-empty one, or None."""
    if subject is not None and not subject.strip():
        raise ValueError("subject must be a non-empty label, or None for no subject")


def _open_text(path: Path) -> io.BufferedIOBase:
    """The file's bytes, decompressed when the magic bytes say gzip (BGZF included)."""
    with path.open("rb") as handle:
        magic = handle.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else path.open("rb")  # type: ignore[return-value]


def _head_lines(path: Path, count: int = 50) -> list[str]:
    lines: list[str] = []
    with _open_text(path) as handle:
        for raw in handle:
            text = raw.decode("utf-8", errors="replace").rstrip("\r\n")
            if text:
                lines.append(text)
            if len(lines) >= count:
                break
    if not lines:
        raise UnsupportedCoverageError("the coverage input is an empty file")
    return lines


def _is_int(value: str) -> bool:
    return value.isdigit()


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _coverage_kind(lines: list[str]) -> CoverageFormat:
    """Which coverage output `lines` (the file's first non-empty lines) come from."""
    first = lines[0].split("\t")
    if first == _MOSDEPTH_SUMMARY_HEADER:
        return CoverageFormat.MOSDEPTH_SUMMARY
    if first == _SAMTOOLS_COVERAGE_HEADER:
        return CoverageFormat.SAMTOOLS_COVERAGE
    body = [line.split("\t") for line in lines if not line.startswith("#")]
    if not body:
        raise UnsupportedCoverageError("no data rows in the first lines")
    widths = {len(fields) for fields in body}
    if len(widths) != 1:
        raise UnsupportedCoverageError(f"rows of different widths ({sorted(widths)}) in the first lines")
    width = widths.pop()
    if lines[0].startswith("#CHROM\tPOS") or (width == 3 and all(_is_int(f[1]) and _is_int(f[2]) for f in body)):
        if width != 3:
            raise UnsupportedCoverageError(
                f"samtools depth output with {width - 2} depth columns: one alignment file per facet; "
                "run samtools depth on one file, or read each column as its own file"
            )
        if not lines[0].startswith("#CHROM") and all(int(f[2]) > int(f[1]) for f in body):
            raise UnsupportedCoverageError(
                "three columns without a header whose third always exceeds the second: a BED3 of intervals "
                "or a samtools depth output, the rows cannot tell; pass kind=CoverageFormat.SAMTOOLS_DEPTH if it is one"
            )
        return CoverageFormat.SAMTOOLS_DEPTH
    if lines[0].startswith("#chrom\tstart\tend\tregion"):
        raise UnsupportedCoverageError("a mosdepth thresholds file: counts per threshold, not a depth (not read yet)")
    if not all(_is_int(f[1]) and _is_int(f[2]) for f in body):
        raise UnsupportedCoverageError("not a BED of depths: the second and third columns are not all positions")
    last = [f[-1] for f in body]
    if width == 4 and all(_is_int(v) for v in last):
        return CoverageFormat.PER_BASE_BED
    if width in (4, 5) and all(_is_number(v) for v in last):
        return CoverageFormat.REGIONS_BED
    if width == 4:
        return CoverageFormat.QUANTIZED_BED
    raise UnsupportedCoverageError(f"a BED of {width} columns whose last is not a depth")


def _read_tsv(path: Path, names: list[str], dtypes: list[pl.DataType], *, skip_header: bool) -> pl.DataFrame:
    with _open_text(path) as handle:
        data = handle.read()
    return pl.read_csv(
        data,
        separator="\t",
        has_header=False,
        skip_rows=1 if skip_header else 0,
        comment_prefix=None if skip_header else "#",
        new_columns=names,
        schema_overrides=dict(zip(names, dtypes)),
    )


def _read_summary(path: Path) -> tuple[pl.DataFrame, list[dict], int]:
    frame = _read_tsv(path, _MOSDEPTH_SUMMARY_HEADER, [pl.String(), pl.UInt64(), pl.UInt64(), pl.Float64(), pl.UInt32(), pl.UInt32()], skip_header=True)
    contigs = frame.filter((pl.col("chrom") != "total") & ~pl.col("chrom").str.ends_with("_region"))
    skipped = frame.height - contigs.height
    lengths = [{"source": row["chrom"], "chrom": canonical_chrom(row["chrom"]), "length": int(row["length"])} for row in contigs.iter_rows(named=True)]
    intervals = contigs.select(
        chrom=pl.col("chrom"),
        pos=pl.lit(1, dtype=pl.UInt32),
        end_pos=pl.col("length").cast(pl.UInt32),
        depth_mean=pl.col("mean").cast(pl.Float32),
        depth_min=pl.col("min"),
        depth_max=pl.col("max"),
        bases=pl.col("length").cast(pl.UInt32),
        label=pl.lit(None, dtype=pl.String),
    )
    return intervals, lengths, skipped


def _read_samtools_coverage(path: Path) -> pl.DataFrame:
    names = ["chrom", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
    frame = _read_tsv(path, names, [pl.String(), pl.UInt32(), pl.UInt32(), pl.UInt64(), pl.UInt64(), pl.Float64(), pl.Float64(), pl.Float64(), pl.Float64()], skip_header=True)
    return frame.select(
        chrom=pl.col("chrom"),
        pos=pl.col("startpos"),
        end_pos=pl.col("endpos"),
        depth_mean=pl.col("meandepth").cast(pl.Float32),
        depth_min=pl.lit(None, dtype=pl.UInt32),
        depth_max=pl.lit(None, dtype=pl.UInt32),
        bases=(pl.col("endpos") - pl.col("startpos") + 1).cast(pl.UInt32),
        label=pl.lit(None, dtype=pl.String),
    )


def _read_regions(path: Path, *, named: bool) -> pl.DataFrame:
    names = ["chrom", "start", "end", "name", "mean"] if named else ["chrom", "start", "end", "mean"]
    dtypes: list[pl.DataType] = [pl.String(), pl.UInt32(), pl.UInt32(), pl.String(), pl.Float64()] if named else [pl.String(), pl.UInt32(), pl.UInt32(), pl.Float64()]
    frame = _read_tsv(path, names, dtypes, skip_header=False)
    return frame.select(
        chrom=pl.col("chrom"),
        pos=(pl.col("start") + 1).cast(pl.UInt32),
        end_pos=pl.col("end"),
        depth_mean=pl.col("mean").cast(pl.Float32),
        depth_min=pl.lit(None, dtype=pl.UInt32),
        depth_max=pl.lit(None, dtype=pl.UInt32),
        bases=(pl.col("end") - pl.col("start")).cast(pl.UInt32),
        label=pl.col("name") if named else pl.lit(None, dtype=pl.String),
    )


def _read_quantized(path: Path) -> tuple[pl.DataFrame, int]:
    """mosdepth's bins are `lower:upper`, half-open (`0:1` is depth 0, `4:inf` four or more)."""
    frame = _read_tsv(path, ["chrom", "start", "end", "label"], [pl.String(), pl.UInt32(), pl.UInt32(), pl.String()], skip_header=False)
    parts = pl.col("label").str.extract_groups(_BIN_LABEL.pattern)
    out = frame.with_columns(parts.alias("bin")).select(
        chrom=pl.col("chrom"),
        pos=(pl.col("start") + 1).cast(pl.UInt32),
        end_pos=pl.col("end"),
        depth_mean=pl.lit(None, dtype=pl.Float32),
        depth_min=pl.col("bin").struct.field("1").cast(pl.UInt32),
        depth_max=pl.when(pl.col("bin").struct.field("2") != "inf").then(pl.col("bin").struct.field("2").cast(pl.UInt32, strict=False) - 1),
        bases=(pl.col("end") - pl.col("start")).cast(pl.UInt32),
        label=pl.col("label"),
    )
    return out, int(out["depth_min"].null_count())


def _read_runs(path: Path, kind: CoverageFormat) -> pl.DataFrame:
    """Every source row as one interval of one depth (`window=None`)."""
    chunks = list(_depth_chunks(path, kind))
    runs = pl.concat(chunks) if chunks else pl.DataFrame(schema={"chrom": pl.String, "pos": pl.Int64, "end_pos": pl.Int64, "depth": pl.Int64})
    return runs.select(
        chrom=pl.col("chrom"),
        pos=pl.col("pos").cast(pl.UInt32),
        end_pos=pl.col("end_pos").cast(pl.UInt32),
        depth_mean=pl.col("depth").cast(pl.Float32),
        depth_min=pl.col("depth").cast(pl.UInt32),
        depth_max=pl.col("depth").cast(pl.UInt32),
        bases=(pl.col("end_pos") - pl.col("pos") + 1).cast(pl.UInt32),
        label=pl.lit(None, dtype=pl.String),
    )


def _depth_chunks(path: Path, kind: CoverageFormat) -> Iterator[pl.DataFrame]:
    """A per-base source as runs `chrom, pos, end_pos, depth` (1-based, inclusive), chunk by chunk."""
    if kind is CoverageFormat.PER_BASE_BED:
        names, dtypes = ["chrom", "start", "end", "depth"], [pl.String(), pl.Int64(), pl.Int64(), pl.Int64()]
    else:
        names, dtypes = ["chrom", "pos", "depth"], [pl.String(), pl.Int64(), pl.Int64()]
    with _open_text(path) as handle:
        tail = b""
        while True:
            block = handle.read(_CHUNK_BYTES)
            data = tail + block
            if not block:
                cut = len(data)
            else:
                cut = data.rfind(b"\n") + 1
                if cut == 0:
                    tail = data
                    continue
            tail = data[cut:]
            if data[:cut].strip():
                frame = pl.read_csv(
                    data[:cut], separator="\t", has_header=False, comment_prefix="#", new_columns=names, schema_overrides=dict(zip(names, dtypes))
                )
                if kind is CoverageFormat.PER_BASE_BED:
                    yield frame.select("chrom", pos=pl.col("start") + 1, end_pos=pl.col("end"), depth=pl.col("depth"))
                else:
                    yield frame.select("chrom", "pos", end_pos=pl.col("pos"), depth=pl.col("depth"))
            if not block:
                return


# --- Windows ------------------------------------------------------------------------------------

_PARTIAL_AGG = ("depth_sum", "bases", "depth_min", "depth_max", "unreadable_bases")


def _window_partials(runs: pl.LazyFrame | pl.DataFrame, window: int) -> pl.LazyFrame:
    """Runs (`chrom, pos, end_pos, depth`, 1-based inclusive) folded into partial window sums.

    A run crossing window edges is cut at them. A negative depth (polars-bio's Int16 wrap) is left out
    of the sums and bounds and counted in `unreadable_bases`."""
    lf = runs.lazy()
    first = (pl.col("pos") - 1) // window
    last = (pl.col("end_pos") - 1) // window
    pieces = (
        lf.with_columns(window=pl.int_ranges(first, last + 1, dtype=pl.Int64))
        .explode("window", empty_as_null=False)
        .with_columns(
            n=pl.min_horizontal(pl.col("end_pos"), (pl.col("window") + 1) * window)
            - pl.max_horizontal(pl.col("pos"), pl.col("window") * window + 1)
            + 1,
            readable=pl.col("depth") >= 0,
        )
    )
    return pieces.group_by("chrom", "window").agg(
        depth_sum=(pl.col("depth") * pl.col("n")).filter(pl.col("readable")).sum().cast(pl.Int64),
        bases=pl.col("n").filter(pl.col("readable")).sum().cast(pl.Int64),
        depth_min=pl.col("depth").filter(pl.col("readable")).min().cast(pl.Int64),
        depth_max=pl.col("depth").filter(pl.col("readable")).max().cast(pl.Int64),
        unreadable_bases=pl.col("n").filter(~pl.col("readable")).sum().cast(pl.Int64),
    )


def _combine(partials: pl.DataFrame) -> pl.DataFrame:
    return partials.group_by("chrom", "window").agg(
        pl.col("depth_sum").sum(),
        pl.col("bases").sum(),
        pl.col("depth_min").min(),
        pl.col("depth_max").max(),
        pl.col("unreadable_bases").sum(),
    )


def _fold_chunks(chunks: Iterator[pl.DataFrame], window: int) -> pl.DataFrame:
    """Fold chunk after chunk into window sums, re-combining as they come so memory stays bounded."""
    folded: pl.DataFrame | None = None
    for chunk in chunks:
        part = _window_partials(chunk, window).collect()
        folded = part if folded is None else _combine(pl.concat([folded, part]))
    if folded is None:
        return pl.DataFrame(schema={"chrom": pl.String, "window": pl.Int64, **dict.fromkeys(_PARTIAL_AGG, pl.Int64)})
    return folded


def _finish_windows(partial: pl.DataFrame, window: int, *, lengths: list[dict] | None, zero_filled: bool) -> pl.DataFrame:
    """Window sums as `COVERAGE_INTERVALS` rows, source contig names made canonical.

    With `lengths`, every window of every listed contig is written and the last one ends at the
    contig's end; with `zero_filled`, a base no run described is depth zero. Without lengths only the
    windows a run reached are written, each `window` bases long."""
    if lengths is not None:
        tiles = pl.DataFrame(
            {"chrom": [c["source"] for c in lengths], "length": [c["length"] for c in lengths]},
            schema={"chrom": pl.String, "length": pl.Int64},
        ).with_columns(window=pl.int_ranges(0, (pl.col("length") + window - 1) // window, dtype=pl.Int64)).explode("window", empty_as_null=False)
        frame = tiles.join(partial, on=["chrom", "window"], how="left").with_columns(
            pl.col("bases", "depth_sum", "unreadable_bases").fill_null(0)
        )
        end = pl.min_horizontal((pl.col("window") + 1) * window, pl.col("length"))
    else:
        frame = partial
        end = (pl.col("window") + 1) * window
    frame = frame.with_columns(pos=pl.col("window") * window + 1, end_pos=end)
    if zero_filled:
        missing = pl.col("end_pos") - pl.col("pos") + 1 - pl.col("bases") - pl.col("unreadable_bases")
        frame = frame.with_columns(
            depth_min=pl.when(missing > 0).then(pl.lit(0, dtype=pl.Int64)).otherwise(pl.col("depth_min")),
            depth_max=pl.col("depth_max").fill_null(pl.when(missing > 0).then(pl.lit(0, dtype=pl.Int64))),
            bases=pl.col("bases") + missing,
        )
    return frame.select(
        chrom=canonical_chrom_expr("chrom"),
        pos=pl.col("pos").cast(pl.UInt32),
        end_pos=pl.col("end_pos").cast(pl.UInt32),
        depth_mean=pl.when(pl.col("bases") > 0).then(pl.col("depth_sum") / pl.col("bases")).cast(pl.Float32),
        depth_min=pl.col("depth_min").cast(pl.UInt32),
        depth_max=pl.col("depth_max").cast(pl.UInt32),
        bases=pl.col("bases").cast(pl.UInt32),
        label=pl.lit(None, dtype=pl.String),
        unreadable_bases=pl.col("unreadable_bases"),
    )


def _ordered(intervals: pl.DataFrame) -> pl.DataFrame:
    """Canonical chrom, then contig rank, `pos`, `end_pos` (charter 12), in the `COVERAGE_INTERVALS` shape."""
    frame = intervals.with_columns(chrom=canonical_chrom_expr("chrom"))
    rank = pl.col("chrom").replace_strict(CONTIG_RANK, default=len(CONTIG_RANK), return_dtype=pl.UInt16)
    return frame.sort([rank, pl.col("chrom"), pl.col("pos"), pl.col("end_pos")], maintain_order=True).select(
        [pl.col(name).cast(dtype) for name, dtype in COVERAGE_INTERVALS.items()]
    )


# --- Depth over alignments (polars-bio) ------------------------------------------------------------


def alignment_coverage(
    path: Path | str,
    *,
    window: int = DEFAULT_WINDOW,
    min_mapping_quality: int = DEFAULT_MIN_MAPPING_QUALITY,
    filter_flag: int = DEFAULT_FILTER_FLAG,
    assembly_label: str | None = None,
    subject: str | None = None,
) -> Coverage:
    """Depth over a BAM or CRAM by polars-bio's `depth()`, in fixed windows (module docstring).

    The per-base rows `depth(per_base=True)` streams are folded into `window`-base windows by
    Polars' streaming engine and never held whole; the windows tile every contig the header names, with uncovered bases as
    zero (`zero_filled`). The read filters are `depth()`'s own and are recorded in `options`. A CRAM
    needs its reference where htslib-style lookup finds it; polars-bio is given none here.

    The windows are spilled to the cache (`coverage/` under `JUST_VCF_CACHE_DIR`, `just_vcf.cache`), so
    the file is read once: a later call with the same file, `window`, `min_mapping_quality` and
    `filter_flag`, under the same `COVERAGE_VERSION` and polars-bio, reads the spill. The file is known
    by what names it in `facet_id` (header, index, size); `assembly_label` is left out of the key, since
    no window depends on it. A spill made with other parameters or versions is never served.

    `subject` is the caller's label for the person the reads are of (`CoverageFacet.subject`); it is
    carried on the facet and keys nothing.
    """
    path = Path(path)
    _check_subject(subject)
    if window < 1:
        raise ValueError("window must be a positive number of bases")
    versions = {"COVERAGE_VERSION": COVERAGE_VERSION, "POLARS_BIO_VERSION": POLARS_BIO_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    cache.sweep_stale(cache.COVERAGE_SPILL, versions)
    kind = _alignment_kind(path)
    scan = pb.scan_cram(str(path)) if kind == "cram" else pb.scan_bam(str(path))
    header = {key: value for key, value in pb.get_metadata(scan).get("header", {}).items() if not key.startswith("_")}
    references = json.loads(header.get("reference_sequences") or "[]")
    if not references:
        raise UnsupportedCoverageError(f"the {kind.upper()} header lists no reference sequence, so no window can be placed")
    lengths = [{"source": ref["name"], "chrom": canonical_chrom(ref["name"]), "length": int(ref["length"])} for ref in references]
    index = _index_sha256(path)
    options = {"filter_flag": filter_flag, "min_mapping_quality": min_mapping_quality}
    source = {"format": CoverageFormat.ALIGNMENT_DEPTH.value, "kind": kind, "size": path.stat().st_size, "header": header, "index": index,
              "options": options, "window": window}
    facet_id = _facet_id({**source, "assembly_label": assembly_label})
    key = _facet_id({**source, "polars_bio": POLARS_BIO_VERSION})
    origin = cache.Origin(source_sha256=None, source_path=str(path.resolve()), facet_id=facet_id)
    with cache.spilling(cache_dir() / "coverage" / key, cache.COVERAGE_SPILL, versions, origin) as spill:
        done = spill.find(key)
        if done is None:
            finished, depth_warnings = _alignment_windows(path, window, filter_flag, min_mapping_quality, lengths)
            part = spill.stage / "intervals.parquet"
            finished.write_parquet(part)
            done = spill.commit(key, {"intervals": [part]}, depth_warnings)
    warnings = done.warnings
    label, basis, label_warnings = _assembly(assembly_label, lengths)
    warnings += label_warnings
    for line in warnings:
        logger.warning("alignment depth: %s", line)
    facet = CoverageFacet(
        facet_id=facet_id,
        format_id=CoverageFormat.ALIGNMENT_DEPTH,
        source_path=str(path),
        source_sha256=None,
        assembly_label=label,
        assembly_label_basis=basis,
        contig_lengths=lengths,
        resolution=Resolution.WINDOW,
        window=window,
        options=options,
        zero_filled=True,
        dropped_fields=[],
        warnings=warnings,
        subject=subject,
    )
    return Coverage(facet=facet, intervals=pl.scan_parquet(done.files("intervals")))


def _alignment_windows(
    path: Path, window: int, filter_flag: int, min_mapping_quality: int, lengths: list[dict]
) -> tuple[pl.DataFrame, list[str]]:
    """polars-bio's per-base depth over `path`, folded into `COVERAGE_INTERVALS` windows, and the
    warning for wrapped depths: the one full read of the file."""
    # Per-base rows, not run-length blocks: polars-bio 0.36 builds every contig's block batches before
    # it emits the first, so block mode's peak grows with the reads (22.6 GB on a 27 GB BAM, 33.6 GB on
    # a 48 GB one); per-base mode emits contig by contig from its dense arrays (13.3 GB on the first).
    bases = pb.depth(
        str(path), filter_flag=filter_flag, min_mapping_quality=min_mapping_quality, use_zero_based=False, per_base=True
    )
    readable = pl.col("coverage") >= 0
    depth = pl.when(readable).then(pl.col("coverage").cast(pl.Int64))
    partial = (
        bases.group_by(chrom=pl.col("contig"), window=(pl.col("pos").cast(pl.Int64) - 1) // window)
        .agg(
            depth_sum=depth.sum().cast(pl.Int64),
            bases=readable.sum().cast(pl.Int64),
            depth_min=depth.min(),
            depth_max=depth.max(),
            unreadable_bases=(~readable).sum().cast(pl.Int64),
        )
        .collect(engine="streaming")
    )
    finished = _finish_windows(partial, window, lengths=lengths, zero_filled=True)
    warnings: list[str] = []
    unreadable = int(finished["unreadable_bases"].sum())
    if unreadable:
        windows = int((finished["unreadable_bases"] > 0).sum())
        warnings.append(
            f"{unreadable:,} base(s) in {windows:,} window(s) came back with a negative depth (polars-bio writes depth as Int16, "
            "wrapping past 32,767); they are left out of those windows' depth, which then describe fewer bases"
        )
    return _ordered(finished), warnings


def _alignment_kind(path: Path) -> str:
    """`bam` or `cram`, from the magic bytes (a BAM is BGZF holding `BAM\\1`)."""
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic == b"CRAM":
        return "cram"
    if magic[:2] == b"\x1f\x8b":
        with gzip.open(path, "rb") as payload:
            if payload.read(4) == b"BAM\x01":
                return "bam"
    raise UnsupportedCoverageError("the alignment input is not a BAM or CRAM by its magic bytes")


def _index_sha256(path: Path) -> str | None:
    for candidate in (Path(f"{path}.bai"), Path(f"{path}.csi"), Path(f"{path}.crai"), path.with_suffix(".bai"), path.with_suffix(".crai")):
        if candidate.is_file():
            return _sha256(candidate)
    return None


# --- Shared ---------------------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _facet_id(key: dict) -> str:
    return hashlib.sha256(json.dumps({**key, "coverage_version": COVERAGE_VERSION}, sort_keys=True, default=str).encode()).hexdigest()[:32]


def _assembly(stated: str | None, lengths: list[dict] | None) -> tuple[str | None, AssemblyLabelBasis | None, list[str]]:
    inferred = infer_assembly([{"name": c["source"], "length": c["length"]} for c in lengths]) if lengths else None
    warnings: list[str] = []
    if stated is not None:
        if inferred and inferred.label != stated:
            warnings.append(f"the stated build {stated} disagrees with the contig lengths, which read {inferred.label}; the stated one is kept")
        return stated, AssemblyLabelBasis.STATED, warnings
    if inferred:
        return inferred.label, AssemblyLabelBasis.INFERRED, warnings
    return None, None, warnings


# --- Callable --------------------------------------------------------------------------------------

CALLABLE_REASONS: frozenset[str] = frozenset({"no_interval", "no_bounds", "straddles", "partly_described"})
"""Why a position's callability is unknown.

- `no_interval`: no interval of the coverage facet covers it.
- `no_bounds`: the interval states a mean or a renamed bin, which bounds no single base.
- `straddles`: the interval's lowest depth is under `min_depth` and its highest is not: some of its
  bases reach it and some do not, and the interval does not say which.
- `partly_described`: the bounds clear the threshold, but the source described only some of the
  interval's bases, and the position may be one of the others.
"""


class Callability(BaseModel):
    """Whether the reads covered one position deeply enough. Truthy only when shown callable."""

    model_config = ConfigDict(frozen=True)

    chrom: str
    pos: int
    min_depth: int
    callable: bool | None
    reason: str | None
    """One of `CALLABLE_REASONS` when `callable` is None."""
    depth_min: int | None
    depth_max: int | None
    depth_mean: float | None

    def __bool__(self) -> bool:
        return self.callable is True


def callable_sites(coverage: Coverage, sites: pl.LazyFrame | pl.DataFrame, *, min_depth: int) -> pl.LazyFrame:
    """`sites` (`chrom`, canonical, and 1-based `pos`) with `callable` and `callable_reason` added.

    `callable` is True where the covering interval shows at least `min_depth` reads at every base it
    describes and describes them all, False where it shows fewer at every base, null otherwise with
    `callable_reason` from `CALLABLE_REASONS`. `min_depth` has no default: what depth makes a base
    callable is the caller's choice (GATK CallableLoci, for one, uses a minimum of 4). The depth
    columns of the covering interval come along as `coverage_depth_min`, `_max` and `_mean`.
    """
    if min_depth < 0:
        raise ValueError("min_depth must not be negative")
    intervals = coverage.intervals.select(
        "chrom",
        interval_pos=pl.col("pos").cast(pl.Int64),
        interval_end=pl.col("end_pos").cast(pl.Int64),
        coverage_depth_min=pl.col("depth_min"),
        coverage_depth_max=pl.col("depth_max"),
        coverage_depth_mean=pl.col("depth_mean"),
        interval_bases=pl.col("bases"),
    ).sort("chrom", "interval_pos")
    queried = sites.lazy().with_columns(_site_row=pl.int_range(pl.len(), dtype=pl.UInt64), _pos=pl.col("pos").cast(pl.Int64)).sort("chrom", "_pos")
    joined = queried.join_asof(intervals, left_on="_pos", right_on="interval_pos", by="chrom", strategy="backward", check_sortedness=False)
    inside = pl.col("interval_end").is_not_null() & (pl.col("_pos") <= pl.col("interval_end"))
    whole = pl.col("interval_bases").cast(pl.Int64) >= pl.col("interval_end") - pl.col("interval_pos") + 1
    low, high = pl.col("coverage_depth_min"), pl.col("coverage_depth_max")
    shown_true = low.is_not_null() & (low >= min_depth)
    shown_false = high.is_not_null() & (high < min_depth)
    callable_expr = pl.when(~inside).then(None).when(shown_true & whole).then(True).when(shown_false & whole).then(False).otherwise(None)
    reason = (
        pl.when(~inside)
        .then(pl.lit("no_interval"))
        .when((shown_true | shown_false) & whole)
        .then(None)
        .when(shown_true | shown_false)
        .then(pl.lit("partly_described"))
        .when(low.is_null())
        .then(pl.lit("no_bounds"))
        .otherwise(pl.lit("straddles"))
    )
    return (
        joined.with_columns(
            callable=callable_expr.cast(pl.Boolean),
            callable_reason=reason.cast(pl.String),
            coverage_depth_min=pl.when(inside).then(low),
            coverage_depth_max=pl.when(inside).then(high),
            coverage_depth_mean=pl.when(inside).then(pl.col("coverage_depth_mean")),
        )
        .sort("_site_row")
        .drop("_site_row", "_pos", "interval_pos", "interval_end", "interval_bases")
    )


def callable_at(coverage: Coverage, chrom: str, pos: int, *, min_depth: int) -> Callability:
    """`callable_sites` for one position; `chrom` is made canonical first."""
    row = (
        callable_sites(coverage, pl.DataFrame({"chrom": [canonical_chrom(chrom)], "pos": [pos]}, schema={"chrom": pl.String, "pos": pl.Int64}), min_depth=min_depth)
        .collect()
        .row(0, named=True)
    )
    return Callability(
        chrom=row["chrom"],
        pos=pos,
        min_depth=min_depth,
        callable=row["callable"],
        reason=row["callable_reason"],
        depth_min=row["coverage_depth_min"],
        depth_max=row["coverage_depth_max"],
        depth_mean=row["coverage_depth_mean"],
    )

