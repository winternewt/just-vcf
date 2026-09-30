"""A BGZF writer on the standard library, for the loader's cached copies.

BGZF (SAM spec §4.1) is gzip written as independent blocks of at most 64 KiB, each carrying its own
size in a `BC` extra subfield, and ended by a fixed empty block. polars-bio reads it as fast as plain
text, where a plain gzip copy reads about 2.7 times slower (measured 2026-09-30). Blocks are compressed
on a thread pool; `zlib` releases the GIL while it works.
"""

import struct
import zlib
from concurrent.futures import ThreadPoolExecutor
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
    """A write-only binary stream that writes BGZF to `raw`. Closing it writes the end block."""

    def __init__(self, raw: BinaryIO, level: int = 1, threads: int = 8) -> None:
        self._raw = raw
        self._level = level
        self._pool = ThreadPoolExecutor(threads)
        self._buffer = bytearray()
        self._pending: list[bytes] = []

    def write(self, data: bytes) -> int:
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
        self._raw.write(EOF_BLOCK)
        self._pool.shutdown()

    def _flush_pending(self) -> None:
        for block in self._pool.map(lambda chunk: compress_block(chunk, self._level), self._pending):
            self._raw.write(block)
        self._pending.clear()

    def __enter__(self) -> "BgzfWriter":
        return self

    def __exit__(self, kind: type[BaseException] | None, value: BaseException | None, tb: TracebackType | None) -> None:
        if kind is None:
            self.close()
        else:
            self._pool.shutdown(cancel_futures=True)
