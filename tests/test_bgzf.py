import gzip
import io
import random
import struct
from pathlib import Path

import polars as pl
import pytest

from just_vcf import DamagedBgzfError, bgzf, load
from just_vcf.bgzf import BLOCK_PAYLOAD, EOF_BLOCK, BgzfWriter, find_damage
from just_vcf.cache import file_sha256
from just_vcf.loaders import vcf
from just_vcf.loading import _facet_id
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


# --- find_damage: the first block that does not read back ----------------------------------------


def _vcf_bgzf(path: Path, rows: int = 60_000) -> list[int]:
    """A BGZF VCF of `rows` SNVs at `path`; returns the compressed start of every block, the end block's last."""
    rng = random.Random(9)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=100000000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        *(f"1\t{1000 + 10 * i}\t.\t{rng.choice('ACGT')}\t{rng.choice('ACGT')}\t50\tPASS\t.\tGT\t0/1" for i in range(rows)),
    ]
    with path.open("wb") as raw, BgzfWriter(raw) as out:
        out.write(("\n".join(lines) + "\n").encode())
    return out.block_starts


def _damaged(path: Path, starts: list[int], block: int, mode: str) -> int:
    """Damage block `block` of the file at `path` in one of three ways; returns the block's offset."""
    data = bytearray(path.read_bytes())
    begin, end = starts[block], starts[block + 1]
    if mode == "inflate":
        for at in range(begin + 40, begin + 60):
            data[at] ^= 0xFF
    elif mode == "isize":
        data[end - 4 : end] = struct.pack("<I", 12345)
    elif mode == "crc32":
        data[end - 8] ^= 0xFF
    path.write_bytes(bytes(data))
    return begin


def test_an_intact_file_has_no_damage(tmp_path: Path) -> None:
    path = tmp_path / "ok.vcf.gz"
    _vcf_bgzf(path)
    assert find_damage(path) is None
    # Without the end block it still reads to the end: a missing EOF marker alone is not damage.
    bare = tmp_path / "no_eof.vcf.gz"
    bare.write_bytes(path.read_bytes()[: -len(EOF_BLOCK)])
    assert find_damage(bare) is None


def test_htslib_bgzip_output_has_no_damage(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    ours = tmp_path / "ours.vcf.gz"
    _vcf_bgzf(ours)
    plain = tmp_path / "plain.vcf"
    plain.write_bytes(gzip.decompress(ours.read_bytes()))
    pysam.tabix_compress(str(plain), str(tmp_path / "theirs.vcf.gz"))
    assert find_damage(tmp_path / "theirs.vcf.gz") is None


@pytest.mark.parametrize("mode", ["inflate", "isize", "crc32"])
def test_a_damaged_block_is_found_at_its_offset(tmp_path: Path, mode: str) -> None:
    path = tmp_path / "bad.vcf.gz"
    starts = _vcf_bgzf(path)
    offset = _damaged(path, starts, len(starts) // 2, mode)
    damage = find_damage(path)
    assert damage is not None and (damage.offset, damage.reason) == (offset, mode)


def test_a_cut_file_is_truncated_at_the_last_block_it_starts(tmp_path: Path) -> None:
    path = tmp_path / "cut.vcf.gz"
    starts = _vcf_bgzf(path)
    path.write_bytes(path.read_bytes()[: starts[5] + 30])
    damage = find_damage(path)
    assert damage is not None and (damage.offset, damage.reason) == (starts[5], "truncated")
    junk = tmp_path / "junk.vcf.gz"
    _vcf_bgzf(junk)
    junk.write_bytes(junk.read_bytes() + b"not a block at all")
    damage = find_damage(junk)
    assert damage is not None and (damage.offset, damage.reason) == (junk.stat().st_size - 18, "not_a_block")


# --- load names the damage ----------------------------------------------------------------------


def test_load_names_the_facet_and_the_block_never_the_path(tmp_path: Path) -> None:
    """The starling upload in docs/PGP_HMS.md, synthetic: one block in the middle does not inflate.
    polars-bio's own error says neither where nor which facet."""
    path = tmp_path / "person_name_genome.vcf.gz"
    starts = _vcf_bgzf(path)
    offset = _damaged(path, starts, len(starts) // 2, "inflate")
    with pytest.raises(DamagedBgzfError) as caught:
        load(path)
    err = caught.value
    facet_id = _facet_id(file_sha256(path), None, "vcf", options={"sample_index": 0})
    assert (err.facet_id, err.offset, err.reason) == (facet_id, offset, "inflate")
    assert f"facet {facet_id}" in str(err) and f"byte {offset}" in str(err)
    assert "person_name" not in str(err) and str(tmp_path) not in str(err)
    assert not any("person_name" in str(value) for value in vars(err).values())
    assert isinstance(err, ValueError) and err.__cause__ is not None


def test_damage_in_the_header_names_the_source(tmp_path: Path) -> None:
    path = tmp_path / "head.vcf.gz"
    starts = _vcf_bgzf(path)
    offset = _damaged(path, starts, 0, "crc32")
    with pytest.raises(DamagedBgzfError) as caught:
        load(path)
    assert (caught.value.facet_id, caught.value.offset, caught.value.source_sha256) == (None, offset, file_sha256(path))
    assert f"source sha256 {file_sha256(path)}" in str(caught.value)


def test_damage_in_plain_gzip_keeps_the_readers_error(tmp_path: Path) -> None:
    """Only BGZF has blocks to name; a plain gzip member fails as the reader says."""
    good = tmp_path / "good.vcf.gz"
    _vcf_bgzf(good)
    path = tmp_path / "plain.vcf.gz"
    data = bytearray(gzip.compress(gzip.decompress(good.read_bytes())))
    for at in range(len(data) // 2, len(data) // 2 + 20):
        data[at] ^= 0xFF
    path.write_bytes(bytes(data))
    assert sniff(path).layers == ("gzip",)
    with pytest.raises(Exception) as caught:
        load(path)
    assert not isinstance(caught.value, DamagedBgzfError)


def test_a_read_failure_that_is_not_damage_walks_the_file_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An intact BGZF input whose read fails for another reason is walked once, not once per level
    that names damage, and fails with the reader's own error."""
    path = tmp_path / "intact.vcf.gz"
    _vcf_bgzf(path, rows=2_000)
    walks: list[Path] = []

    def counted(target: Path, *args, **kwargs):
        walks.append(target)
        return find_damage(target, *args, **kwargs)

    def refuse(_records):
        raise pl.exceptions.ComputeError("a read failure that is not damage")

    monkeypatch.setattr(bgzf, "find_damage", counted)
    monkeypatch.setattr(vcf, "scan_stats", refuse)
    with pytest.raises(pl.exceptions.ComputeError, match="not damage"):
        load(path)
    assert len(walks) == 1
