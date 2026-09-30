import gzip
import io
import random
from pathlib import Path

from just_vcf.bgzf import BLOCK_PAYLOAD, EOF_BLOCK, BgzfWriter
from just_vcf.sniff import sniff


def _write(payload: bytes, chunk: int) -> bytes:
    raw = io.BytesIO()
    with BgzfWriter(raw) as out:
        for start in range(0, len(payload), chunk):
            out.write(payload[start : start + chunk])
    return raw.getvalue()


def test_round_trip_across_block_boundaries() -> None:
    rng = random.Random(3)
    for size in (0, 1, BLOCK_PAYLOAD - 1, BLOCK_PAYLOAD, BLOCK_PAYLOAD + 1, 300 * BLOCK_PAYLOAD + 17):
        payload = bytes(rng.getrandbits(8) for _ in range(min(size, 4096))) * (size // 4096 + 1)
        payload = payload[:size]
        for chunk in (1000, 70000):
            data = _write(payload, chunk)
            assert gzip.decompress(data) == payload, (size, chunk)
            assert data.endswith(EOF_BLOCK)


def test_every_block_fits_in_64_kib_and_states_its_size() -> None:
    data = _write(bytes(random.Random(5).getrandbits(8) for _ in range(5 * BLOCK_PAYLOAD)), 1 << 20)
    offset, blocks = 0, 0
    while offset < len(data):
        assert data[offset : offset + 4] == b"\x1f\x8b\x08\x04" and data[offset + 12 : offset + 14] == b"BC"
        size = int.from_bytes(data[offset + 16 : offset + 18], "little") + 1
        assert size <= 65536
        offset += size
        blocks += 1
    assert offset == len(data) and blocks == 6


def test_sniff_sees_bgzf(tmp_path: Path) -> None:
    path = tmp_path / "calls.vcf.gz"
    path.write_bytes(_write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tG\t.\t.\t.\n", 7))
    assert sniff(path).layers == ("bgzf",)
