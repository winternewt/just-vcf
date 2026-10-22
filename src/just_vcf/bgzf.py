"""A BGZF writer on the standard library, for the loader's cached copies, and a block walker that
finds where a BGZF file is damaged.

BGZF (SAM spec §4.1) is gzip written as independent blocks of at most 64 KiB, each carrying its own
size in a `BC` extra subfield, and ended by a fixed empty block. polars-bio reads it as fast as plain
text, where a plain gzip copy reads about 2.7 times slower (measured 2026-09-30). Blocks are compressed
on a thread pool; `zlib` releases the GIL while it works.
"""

import struct
import zlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import BinaryIO

BLOCK_PAYLOAD = 65280
"""Uncompressed bytes per block: the size htslib uses, so a compressed block always fits in 64 KiB."""

EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
"""The empty block that ends every BGZF file (SAM spec §4.1.2)."""

_BATCH = 256


def compress_block(payload: bytes, level: int = 1) -> bytes:
    """One BGZF block holding `payload` (at most `BLOCK_PAYLOAD` bytes)."""
    compressor = zlib.compressobj(level, zlib.DEFLATED, -15)
    deflated = compressor.compress(payload) + compressor.flush()
    block_size = 18 + len(deflated) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", block_size - 1)
    return header + deflated + struct.pack("<II", zlib.crc32(payload), len(payload))


class BgzfWriter:
    """A write-only binary stream that writes BGZF to `raw`. Closing it writes the end block.

    Every block but the last holds exactly `BLOCK_PAYLOAD` uncompressed bytes, so the uncompressed
    offset of any byte written names its block; `block_starts` records where each block begins in
    the compressed stream, which is all a tabix index needs (`virtual_offset`).
    """

    def __init__(self, raw: BinaryIO, level: int = 1, threads: int = 8) -> None:
        self._raw = raw
        self._level = level
        self._pool = ThreadPoolExecutor(threads)
        self._buffer = bytearray()
        self._pending: list[bytes] = []
        self._compressed = 0
        self._written = 0
        self.block_starts: list[int] = []
        """Compressed offset of each block written, the end block last once closed."""

    def write(self, data: bytes) -> int:
        self._written += len(data)
        self._buffer += data
        while len(self._buffer) >= BLOCK_PAYLOAD:
            self._pending.append(bytes(self._buffer[:BLOCK_PAYLOAD]))
            del self._buffer[:BLOCK_PAYLOAD]
            if len(self._pending) >= _BATCH:
                self._flush_pending()
        return len(data)

    def close(self) -> None:
        if self._buffer:
            self._pending.append(bytes(self._buffer))
            self._buffer.clear()
        self._flush_pending()
        self.block_starts.append(self._compressed)
        self._raw.write(EOF_BLOCK)
        self._compressed += len(EOF_BLOCK)
        self._pool.shutdown()

    def virtual_offset(self, offset: int) -> int:
        """The BGZF virtual offset (SAM spec §4.1.1) of uncompressed byte `offset`, once closed:
        the compressed start of its block, shifted 16 bits, plus the offset inside the block. The end
        of a block is the start of the next, as htslib's reader reports it, so the offset just past
        the last byte written is the end block's start."""
        if offset == self._written:
            return self.block_starts[-1] << 16
        block, within = divmod(offset, BLOCK_PAYLOAD)
        return (self.block_starts[block] << 16) | within

    def _flush_pending(self) -> None:
        for block in self._pool.map(lambda chunk: compress_block(chunk, self._level), self._pending):
            self.block_starts.append(self._compressed)
            self._raw.write(block)
            self._compressed += len(block)
        self._pending.clear()

    def __enter__(self) -> "BgzfWriter":
        return self

    def __exit__(self, kind: type[BaseException] | None, value: BaseException | None, tb: TracebackType | None) -> None:
        if kind is None:
            self.close()
        else:
            self._pool.shutdown(cancel_futures=True)


@dataclass(frozen=True)
class BlockDamage:
    """The first block of a BGZF file that does not read back."""

    offset: int
    """The block's compressed byte offset in the file (where its gzip header starts)."""
    reason: str
    """What is wrong with it: `not_a_block` (no gzip header with a `BC` subfield where a block must
    start), `truncated` (the file ends inside the block), `inflate` (the deflate data is invalid),
    `crc32` or `isize` (it inflates, to other bytes or another length than the block states)."""
    detail: str
    """The specifics, for the message: zlib's text, or the stated and found values."""


_FIXED = struct.Struct("<4BI2BH")
"""A gzip member's first 12 bytes: ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN (when FLG has FEXTRA)."""


def find_damage(path: Path, chunk: int = 1 << 20) -> BlockDamage | None:
    """Walk `path` block by block and return the first block that does not read back, or None.

    Each block's header must be gzip's, with FEXTRA and a `BC` subfield giving the block's size; its
    deflate data must inflate to exactly ISIZE bytes whose CRC32 is the one stated. A file that ends
    cleanly after a block is not damaged, whether or not it carries the empty end block. Reads the
    whole file, so it is for the failure path, not for every load.
    """
    with path.open("rb", buffering=chunk) as handle:
        offset = 0
        while True:
            head = handle.read(_FIXED.size)
            if not head:
                return None
            if len(head) < _FIXED.size:
                return BlockDamage(offset, "truncated", f"{len(head)} byte(s) where a block header starts")
            id1, id2, method, flags, _, _, _, xlen = _FIXED.unpack(head)
            if (id1, id2, method) != (0x1F, 0x8B, 8) or not flags & 4:
                return BlockDamage(offset, "not_a_block", "no gzip header with an extra field")
            extra = handle.read(xlen)
            if len(extra) < xlen:
                return BlockDamage(offset, "truncated", "the file ends inside the block header")
            bsize = _block_size(extra)
            if bsize is None:
                return BlockDamage(offset, "not_a_block", "no BC subfield in the gzip extra field")
            rest = bsize + 1 - _FIXED.size - xlen
            if rest < 8:
                return BlockDamage(offset, "not_a_block", f"a block size of {bsize + 1} bytes is too short")
            body = handle.read(rest)
            if len(body) < rest:
                return BlockDamage(offset, "truncated", f"the file ends {rest - len(body)} byte(s) short of the block's end")
            crc, isize = struct.unpack("<II", body[-8:])
            inflater = zlib.decompressobj(-15)
            try:
                data = inflater.decompress(body[:-8]) + inflater.flush()
            except zlib.error as err:
                return BlockDamage(offset, "inflate", str(err))
            if not inflater.eof:
                return BlockDamage(offset, "inflate", "the deflate stream does not end inside the block")
            if len(data) != isize:
                return BlockDamage(offset, "isize", f"inflates to {len(data)} bytes, the block states {isize}")
            if zlib.crc32(data) != crc:
                return BlockDamage(offset, "crc32", f"CRC32 {zlib.crc32(data):08x}, the block states {crc:08x}")
            offset += bsize + 1


def _block_size(extra: bytes) -> int | None:
    """BSIZE (the block's total size less one) from a gzip extra field, or None without a `BC` subfield."""
    at = 0
    while at + 4 <= len(extra):
        si1, si2, length = extra[at], extra[at + 1], struct.unpack_from("<H", extra, at + 2)[0]
        if (si1, si2, length) == (66, 67, 2) and at + 6 <= len(extra):
            return struct.unpack_from("<H", extra, at + 4)[0]
        at += 4 + length
    return None
