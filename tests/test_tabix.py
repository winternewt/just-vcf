"""RM6: `write_vcf` writes a tabix index (`.tbi`) beside its BGZF output, built from the offsets it wrote.

The index is read back two ways: by a small reader of the format here (layout, names, the linear
index computed from the frame), and, opt-in through pysam (the `[reference]` extra), by htslib itself:
region queries against the frame, and the linear index against the one htslib builds for the same file.
"""

import gzip
import random
import shutil
import struct
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ContigNames, load, write_vcf
from just_vcf.bgzf import BLOCK_PAYLOAD, BgzfWriter
from just_vcf.tabix import META_BIN, reg2bin

rng = random.Random(11)

HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=249250621>",
    "##contig=<ID=2,length=243199373>",
    "##contig=<ID=MT,length=16569>",
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
]


def records() -> list[str]:
    """SNVs spread over five megabases of contig 1, deletions, symbolic deletions whose END reaches
    over many windows and up the bin levels, and a few rows on 2 and MT: enough for several BGZF blocks."""
    rows: list[tuple[str, int, str]] = []
    for pos in sorted(rng.sample(range(1, 5_000_000), 6000)):
        rows.append(("1", pos, f"1\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"))
    for pos in sorted(rng.sample(range(100, 5_000_000), 40)):
        rows.append(("1", pos, f"1\t{pos}\t.\tACGTACGTAC\tA\t.\tPASS\t.\tGT\t0/1"))
    for pos, length in ((40_000, 20_000), (1_000_000, 300_000), (2_500_000, 3_000_000), (4_999_000, 70)):
        rows.append(("1", pos, f"1\t{pos}\t.\tC\t<DEL>\t.\tPASS\tSVLEN=-{length};END={pos + length}\tGT\t0/1"))
    for pos in (10, 16_384, 16_385, 99_999):
        rows.append(("2", pos, f"2\t{pos}\t.\tT\tC\t.\tPASS\t.\tGT\t1/1"))
    rows.append(("MT", 150, "MT\t150\t.\tC\tT\t.\tPASS\t.\tGT\t1"))
    rank = {"1": 0, "2": 1, "MT": 2}
    return [line for _, _, line in sorted(rows, key=lambda r: (rank[r[0]], r[1]))]


@pytest.fixture
def written(tmp_path: Path):
    source = tmp_path / "in.vcf"
    source.write_text("\n".join([*HEADER, *records()]) + "\n")
    loaded = load(source)
    return loaded, write_vcf(loaded, tmp_path / "out.vcf.gz", contig_names=ContigNames.CHR)


def read_tbi(path: Path) -> dict:
    """The `.tbi` format, read field by field (htslib's tabix layout)."""
    data = gzip.decompress(path.read_bytes())
    assert data[:4] == b"TBI\x01"
    n_ref, preset, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", data, 4)
    at = 36
    names = data[at : at + l_nm].split(b"\0")[:-1]
    at += l_nm
    refs = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", data, at)
        at += 4
        bins = {}
        for _ in range(n_bin):
            bin_, n_chunk = struct.unpack_from("<Ii", data, at)
            at += 8
            bins[bin_] = [struct.unpack_from("<QQ", data, at + 16 * k) for k in range(n_chunk)]
            at += 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", data, at)
        at += 4
        linear = list(struct.unpack_from(f"<{n_intv}Q", data, at))
        at += 8 * n_intv
        refs.append({"bins": bins, "linear": linear})
    return {
        "layout": (preset, col_seq, col_beg, col_end, meta, skip),
        "names": [name.decode() for name in names],
        "refs": refs,
        "rest": data[at:],
    }


def intervals(text: str) -> list[tuple[str, int, int, int]]:
    """Each data line's contig, 0-based start and end as htslib's VCF preset reads them, and its byte offset."""
    found = []
    offset = 0
    for line in text.splitlines(keepends=True):
        if not line.startswith("#"):
            fields = line.rstrip("\n").split("\t")
            beg = int(fields[1]) - 1
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            end = int(info["END"]) if "END" in info else beg + len(fields[3])
            found.append((fields[0], beg, end, offset))
        offset += len(line.encode())
    return found


def test_bgzf_output_gets_a_tabix_index(written) -> None:
    loaded, out = written
    assert out.index_path == out.path.with_name("out.vcf.gz.tbi") and out.index_path.exists()
    index = read_tbi(out.index_path)
    # The VCF preset: sequence column 1, begin column 2, no end column, `#` meta lines.
    assert index["layout"] == (2, 1, 2, 0, ord("#"), 0)
    assert index["names"] == ["chr1", "chr2", "chrM"], "contigs as written, in the order the file first uses them"
    assert index["rest"] == struct.pack("<Q", 0)

    text = gzip.decompress(out.path.read_bytes()).decode()
    lines = intervals(text)
    assert len(lines) == loaded.calls.collect().height
    # Offsets are virtual: the compressed block's start, then the offset inside it. Recompute the block
    # starts from the file itself.
    raw = out.path.read_bytes()
    starts, at = [], 0
    while at < len(raw):
        starts.append(at)
        at += struct.unpack_from("<H", raw, at + 16)[0] + 1

    def virtual(u: int) -> int:
        return (starts[u // BLOCK_PAYLOAD] << 16) | (u % BLOCK_PAYLOAD)

    assert len(starts) > 3, "the fixture spans several BGZF blocks"
    for tid, name in enumerate(index["names"]):
        mine = [(beg, end, virtual(u)) for chrom, beg, end, u in lines if chrom == name]
        ref = index["refs"][tid]
        # Linear index: per 16 kb window, the smallest offset of a record overlapping it; an empty
        # window takes the smallest offset of any record after it.
        expected: dict[int, int] = {}
        for beg, end, v in mine:
            for window in range(beg >> 14, ((end - 1) >> 14) + 1):
                expected.setdefault(window, v)
        linear = [min(v for w, v in expected.items() if w >= window) for window in range(max(expected) + 1)]
        assert ref["linear"] == linear
        # Every record is inside a chunk of its own bin, and the pseudo-bin counts the contig's records.
        for beg, end, v in mine:
            bin_ = pl.select(reg2bin(pl.lit(beg, dtype=pl.Int64), pl.lit(end, dtype=pl.Int64))).item()
            assert any(start <= v < stop for start, stop in ref["bins"][bin_]), (name, beg, end)
        assert ref["bins"][META_BIN][1] == (len(mine), 0)


def test_reg2bin_matches_the_spec_examples() -> None:
    cases = {(0, 1): 4681, (16_383, 16_384): 4681, (16_384, 16_385): 4682, (0, 16_385): 585, (0, 1 << 17): 585,
             (0, (1 << 17) + 1): 73, (0, 1 << 26): 1, (0, (1 << 26) + 1): 0, ((1 << 26) - 1, (1 << 26) + 1): 0}  # fmt: skip
    frame = pl.DataFrame({"beg": [k[0] for k in cases], "end": [k[1] for k in cases]})
    assert frame.select(reg2bin(pl.col("beg"), pl.col("end")).alias("bin"))["bin"].to_list() == list(cases.values())


def test_no_index_for_plain_text_and_a_stale_one_is_removed(written, tmp_path: Path) -> None:
    loaded, out = written
    stale = out.index_path.read_bytes()
    bare = write_vcf(loaded, out.path, index=False, contig_names=ContigNames.CHR)
    assert bare.index_path is None and not out.index_path.exists(), "an index left from an earlier write no longer describes the file"
    assert bare.sha256 == out.sha256
    again = write_vcf(loaded, out.path, contig_names=ContigNames.CHR)
    assert again.index_path.read_bytes() == stale, "the same frames give the same index bytes"
    plain = write_vcf(loaded, tmp_path / "out.vcf", compress=False)
    assert plain.index_path is None and not (tmp_path / "out.vcf.tbi").exists()
    with pytest.raises(ValueError, match="needs BGZF"):
        write_vcf(loaded, tmp_path / "other.vcf", compress=False, index=True)
    assert not (tmp_path / "other.vcf").exists()


def test_virtual_offsets_name_the_block_and_the_byte_inside_it() -> None:
    class Sink:
        def __init__(self) -> None:
            self.data = bytearray()

        def write(self, chunk: bytes) -> int:
            self.data += chunk
            return len(chunk)

    sink = Sink()
    payload = bytes(rng.getrandbits(8) for _ in range(2 * BLOCK_PAYLOAD + 10))
    with BgzfWriter(sink) as out:
        out.write(payload)
    raw = bytes(sink.data)
    for u in (0, 5, BLOCK_PAYLOAD - 1, BLOCK_PAYLOAD, 2 * BLOCK_PAYLOAD + 9, 2 * BLOCK_PAYLOAD + 10):
        v = out.virtual_offset(u)
        block, within = v >> 16, v & 0xFFFF
        size = struct.unpack_from("<H", raw, block + 16)[0] + 1
        content = gzip.decompress(raw[block : block + size])
        assert content[within:] == payload[u : u + len(content) - within]
    # The end of the data is the end block's start, as htslib's reader reports a position at a block's end.
    assert out.virtual_offset(len(payload)) == out.block_starts[3] << 16
    assert len(out.block_starts) == 4, "three data blocks, then the end block"


def test_htslib_answers_region_queries_from_the_index(written, tmp_path: Path) -> None:
    """Opt-in: htslib (pysam) fetches regions through our index, and builds the same linear index."""
    pysam = pytest.importorskip("pysam")
    loaded, out = written
    text = gzip.decompress(out.path.read_bytes()).decode()
    lines = intervals(text)
    with pysam.TabixFile(str(out.path)) as tabix:
        assert list(tabix.contigs) == ["chr1", "chr2", "chrM"]
        windows = [("chr1", s, s + w) for s, w in ((rng.randrange(0, 5_400_000), rng.choice([1, 50, 5_000, 200_000])) for _ in range(300))]
        windows += [("chr1", 0, 6_000_000), ("chr2", 16_383, 16_384), ("chr2", 0, 10), ("chrM", 149, 150), ("chr1", 3_000_000, 3_000_001)]
        for chrom, start, stop in windows:
            got = sorted(int(line.split("\t")[1]) for line in tabix.fetch(chrom, start, stop))
            want = sorted(beg + 1 for name, beg, end, _ in lines if name == chrom and beg < stop and end > start)
            assert got == want, (chrom, start, stop)
    with pysam.VariantFile(str(out.path)) as reader:
        assert sum(1 for _ in reader.fetch("chr1", 2_600_000, 2_600_001)) >= 1, "the symbolic deletion spanning it, by its END"

    # htslib's own index of a copy: the linear index must be the same (bins may differ only where
    # htslib folds a small bin into its parent, which this writer does not).
    copy = tmp_path / "copy.vcf.gz"
    shutil.copyfile(out.path, copy)
    pysam.tabix_index(str(copy), preset="vcf", force=True, keep_original=True)
    theirs, ours = read_tbi(copy.with_name("copy.vcf.gz.tbi")), read_tbi(out.index_path)
    assert theirs["names"] == ours["names"] and theirs["layout"] == ours["layout"]
    for mine, htslib in zip(ours["refs"], theirs["refs"]):
        assert mine["linear"] == htslib["linear"]
        assert mine["bins"][META_BIN] == htslib["bins"][META_BIN]
