"""A tabix index (`.tbi`) for a BGZF VCF this package writes, built while the lines are written.

The format is htslib's (the tabix paper and htslib's `tbx.c`/`hts.c`): a magic `TBI\\1`, the VCF
preset's column layout, the contig names in the order the file first uses them, then per contig the
binning index (UCSC's five-level bins over 0-based half-open intervals, each bin a list of chunks of
virtual offsets) and the linear index (per 16 kb window, the smallest virtual offset of a record that
overlaps it, an empty window taking the next one's). The whole file is BGZF-compressed.

A record's interval is the one htslib's VCF preset gives it: from `POS - 1`, to INFO `END` where the
record writes one, else to `POS - 1 + len(REF)`. The writer hands that interval over with each line,
so the index never re-parses its own output.

What is left out, both optional in the format and both only size: htslib's `compress_binning`, which
moves the chunks of a small bin into its parent. The index stays valid without it. The pseudo-bin
37450 with each contig's offset range and record count is written, as htslib does.
"""

import struct
from collections.abc import Callable

import polars as pl

MIN_SHIFT = 14
"""Linear-index window: 16 kb."""
MAX_POSITION = 1 << 29
"""The largest coordinate a `.tbi` can bin (five levels over 2^29); a longer contig needs CSI."""
META_BIN = 37450
"""htslib's pseudo-bin holding a contig's offset range and record count."""

_PRESET_VCF = 2
_META_CHAR = ord("#")


class TabixLimitError(ValueError):
    """A record past 2^29, which a `.tbi` cannot index."""


def reg2bin(beg: pl.Expr, end: pl.Expr) -> pl.Expr:
    """The smallest bin holding `[beg, end)`, 0-based half-open (SAM spec §5.3, `reg2bin`)."""
    last = end - 1
    expr = pl.lit(0, dtype=pl.Int64)
    # Built from the widest level in, so the innermost `when` that holds wins.
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        size = 1 << shift
        expr = pl.when(beg // size == last // size).then(offset + beg // size).otherwise(expr)
    return expr


class TabixIndexer:
    """Collects, batch by batch, what a `.tbi` needs from the lines as they are written.

    Each batch gives every record's contig name as written, its 0-based `beg` and exclusive `end`,
    and its uncompressed start and end offsets in the file. Offsets stay uncompressed until `build`,
    which takes the function turning them into virtual offsets (`BgzfWriter.virtual_offset`).
    """

    def __init__(self) -> None:
        self._names: list[str] = []
        self._tid: dict[str, int] = {}
        self._chunks: list[dict[int, list[tuple[int, int]]]] = []
        self._windows: list[dict[int, int]] = []
        self._span: list[list[int]] = []
        """Per contig: first record's start offset, last record's end offset, records."""

    def add(self, batch: pl.DataFrame) -> None:
        """`batch` holds `name`, `beg`, `end`, `u_start`, `u_end` (Int64, offsets uncompressed)."""
        if batch.is_empty():
            return
        if batch["end"].max() > MAX_POSITION:
            raise TabixLimitError(f"a record ends past {MAX_POSITION}, which a .tbi cannot index; pass index=False")
        for name in batch["name"].unique(maintain_order=True).to_list():
            if name not in self._tid:
                self._tid[name] = len(self._names)
                self._names.append(name)
                self._chunks.append({})
                self._windows.append({})
                self._span.append([-1, 0, 0])
        work = batch.with_columns(
            pl.col("name").replace_strict(self._tid, return_dtype=pl.Int64).alias("tid"),
            pl.max_horizontal(pl.col("end"), pl.col("beg") + 1).alias("end"),
        ).with_columns(reg2bin(pl.col("beg"), pl.col("end")).alias("bin"))

        runs = (
            work.with_columns(pl.struct("tid", "bin").rle_id().alias("run"))
            .group_by("run", maintain_order=True)
            .agg(pl.col("tid").first(), pl.col("bin").first(), pl.col("u_start").first(), pl.col("u_end").last())
        )
        for tid, bin_, start, end in runs.select("tid", "bin", "u_start", "u_end").iter_rows():
            self._chunks[tid].setdefault(bin_, []).append((start, end))

        windows = (
            work.select(
                "tid",
                "u_start",
                pl.int_ranges(pl.col("beg") // (1 << MIN_SHIFT), (pl.col("end") - 1) // (1 << MIN_SHIFT) + 1).alias("window"),
            )
            .explode("window", empty_as_null=False)
            .group_by("tid", "window")
            .agg(pl.col("u_start").min())
        )
        # Batches come in file order, so a window set by an earlier batch already holds the smaller offset.
        for tid, window, start in windows.iter_rows():
            self._windows[tid].setdefault(window, start)

        spans = work.group_by("tid", maintain_order=True).agg(pl.col("u_start").first(), pl.col("u_end").last(), pl.len())
        for tid, start, end, count in spans.iter_rows():
            span = self._span[tid]
            if span[0] < 0:
                span[0] = start
            span[1] = end
            span[2] += count

    def build(self, virtual: Callable[[int], int]) -> bytes:
        """The uncompressed `.tbi` content, offsets turned virtual by `virtual`."""
        names = b"".join(name.encode() + b"\0" for name in self._names)
        out = [b"TBI\x01", struct.pack("<7i", len(self._names), _PRESET_VCF, 1, 2, 0, _META_CHAR, 0), struct.pack("<i", len(names)), names]
        for tid in range(len(self._names)):
            bins = {bin_: _merged([(virtual(start), virtual(end)) for start, end in chunks]) for bin_, chunks in self._chunks[tid].items()}
            first, last, count = self._span[tid]
            bins[META_BIN] = [(virtual(first), virtual(last)), (count, 0)]
            out.append(struct.pack("<i", len(bins)))
            for bin_ in sorted(bins):
                out.append(struct.pack("<Ii", bin_, len(bins[bin_])))
                out.extend(struct.pack("<QQ", start, end) for start, end in bins[bin_])
            # A window no record overlaps takes the next window's offset: every record overlapping a
            # later window starts after it. That is htslib's own fill (checked against its index).
            windows = self._windows[tid]
            linear = [0] * (max(windows) + 1)
            filled = virtual(windows[len(linear) - 1])
            for window in reversed(range(len(linear))):
                if window in windows:
                    filled = virtual(windows[window])
                linear[window] = filled
            out.append(struct.pack("<i", len(linear)))
            out.append(struct.pack(f"<{len(linear)}Q", *linear))
        out.append(struct.pack("<Q", 0))
        return b"".join(out)


def _merged(chunks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """A bin's chunks, in file order, with a chunk joined to the one before when it starts in the
    BGZF block where that one ends (htslib's rule: reading on through one block costs nothing)."""
    merged: list[tuple[int, int]] = []
    for start, end in chunks:
        if merged and merged[-1][1] >> 16 == start >> 16:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged
