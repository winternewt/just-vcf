"""Container peeling and outcome kinds of `sniff`, built from synthetic files in tmp_path."""

import bz2
import gzip
import io
import lzma
import struct
import tarfile
import zipfile
import zlib
from pathlib import Path

from just_vcf.formats import get_spec
from just_vcf.sniff import HEAD_BYTE_CAP, Sniffed, SniffOutcome, sniff

TWENTYTHREE = get_spec("23andme").example
VCF = get_spec("vcf").example


def _bgzf_block(payload: bytes) -> bytes:
    """One BGZF block (SAM spec §4.1): gzip with FEXTRA and a BC subfield holding the block size."""
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    deflated = compressor.compress(payload) + compressor.flush()
    block_size = 18 + len(deflated) + 8
    header = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, block_size - 1)
    return header + deflated + struct.pack("<II", zlib.crc32(payload), len(payload))


def _write(path: Path, data: bytes) -> Path:
    path.write_bytes(data)
    return path


def test_plain_text_is_detected_with_build_from_header(tmp_path: Path) -> None:
    result = sniff(_write(tmp_path / "a.txt", TWENTYTHREE.encode()))
    assert result
    assert (result.format, result.layers, result.build) == ("23andme", (), "GRCh37")
    assert "reference human assembly build 37" in result.build_evidence


def test_build_is_unknown_when_the_header_does_not_say(tmp_path: Path) -> None:
    no_build = "\n".join(line for line in TWENTYTHREE.splitlines() if "assembly build" not in line)
    result = sniff(_write(tmp_path / "a.txt", no_build.encode()))
    assert result.format == "23andme"
    assert result.build is None and result.build_evidence is None


def test_conflicting_build_hints_withhold_the_build(tmp_path: Path) -> None:
    both = TWENTYTHREE.replace("#\n# rsid", "# reference human assembly build 36\n# rsid", 1)
    result = sniff(_write(tmp_path / "a.txt", both.encode()))
    assert result.format == "23andme"
    assert result.build is None
    assert any("conflicting builds" in w for w in result.warnings)


def test_compression_is_detected_by_magic_not_extension(tmp_path: Path) -> None:
    payload = VCF.encode()
    cases = {
        "gzip": gzip.compress(payload),
        "bgzf": _bgzf_block(payload) + _bgzf_block(b""),
        "bz2": bz2.compress(payload),
        "xz": lzma.compress(payload),
    }
    for layer, data in cases.items():
        result = sniff(_write(tmp_path / f"{layer}.vcf", data))
        assert (result.format, result.layers) == ("vcf", (layer,)), layer


def test_zip_skips_macos_forks_and_reports_the_member(tmp_path: Path) -> None:
    path = tmp_path / "export.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("__MACOSX/._genome.txt", b"\x00\x05\x16\x07" + b"\x00" * 4096)
        zf.writestr("._genome.txt", b"\x00\x05\x16\x07" + b"\x00" * 4096)
        zf.writestr("genome.txt", TWENTYTHREE)
    result = sniff(path)
    assert (result.format, result.layers, result.member, result.members) == ("23andme", ("zip",), "genome.txt", ())


def test_nested_zip_then_gzip(tmp_path: Path) -> None:
    path = tmp_path / "calls.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("calls.vcf.gz", gzip.compress(VCF.encode()))
    result = sniff(path)
    assert (result.format, result.layers, result.member) == ("vcf", ("zip", "gzip"), "calls.vcf.gz")


def test_multi_member_archive_leads_with_largest_detected_and_keeps_all(tmp_path: Path) -> None:
    path = tmp_path / "bundle.tar.gz"
    small_vcf = VCF.encode()
    large_vcf = VCF.encode() + b"".join(b"1\t%d\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n" % (800000 + i) for i in range(50))
    notes = b"just some notes\n"
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        for name, data in (("b/small.vcf", small_vcf), ("a/large.vcf", large_vcf), ("readme.txt", notes)):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    result = sniff(_write(path, gzip.compress(buf.getvalue())))
    assert result.format == "vcf"
    assert result.member == "a/large.vcf"
    assert result.layers == ("gzip", "tar")
    sizes = {"a/large.vcf": len(large_vcf), "b/small.vcf": len(small_vcf), "readme.txt": len(notes)}
    assert [m.member for m in result.members] == sorted(sizes, key=lambda n: (-sizes[n], n))
    assert {m.member: m.outcome for m in result.members}["readme.txt"] is SniffOutcome.UNKNOWN


def test_members_with_different_formats_are_warned_about(tmp_path: Path) -> None:
    path = tmp_path / "mixed.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("calls.vcf", VCF)
        zf.writestr("genome.txt", TWENTYTHREE)
    result = sniff(path)
    assert {m.format for m in result.members} == {"vcf", "23andme"}
    assert any("different formats" in w for w in result.warnings)


def test_empty_zip_is_unreadable_not_unknown(tmp_path: Path) -> None:
    path = tmp_path / "empty.zip"
    with zipfile.ZipFile(path, "w"):
        pass
    result = sniff(path)
    assert result.outcome is SniffOutcome.UNREADABLE
    assert not result
    assert result.reason == "archive has no data members"


def test_office_documents_are_not_variant_data(tmp_path: Path) -> None:
    path = tmp_path / "sheet.xlsx"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("[Content_Types].xml", "<Types/>")
        zf.writestr("xl/workbook.xml", "<workbook/>")
    result = sniff(path)
    assert (result.outcome, result.kind) == (SniffOutcome.NOT_VARIANT_DATA, "ooxml")


def test_binary_and_text_kinds(tmp_path: Path) -> None:
    cases = {
        "pdf": b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n",
        "jpeg": b"\xff\xd8\xff\xe0\x00\x10JFIF",
        "fasta": b">chrM\nGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG\n",
        "gedcom": b"0 HEAD\n1 SOUR FTW\n",
    }
    for kind, data in cases.items():
        result = sniff(_write(tmp_path / f"f.{kind}", data))
        assert (result.outcome, result.kind) == (SniffOutcome.NOT_VARIANT_DATA, kind), kind
    bam = sniff(_write(tmp_path / "reads", _bgzf_block(b"BAM\x01" + b"\x00" * 64) + _bgzf_block(b"")))
    assert (bam.outcome, bam.kind, bam.layers) == (SniffOutcome.NOT_VARIANT_DATA, "bam", ("bgzf",))


def test_bcf_is_detected_by_magic(tmp_path: Path) -> None:
    result = sniff(_write(tmp_path / "calls", _bgzf_block(b"BCF\x02\x02" + b"\x00" * 64) + _bgzf_block(b"")))
    assert (result.outcome, result.format, result.layers) == (SniffOutcome.DETECTED, "bcf", ("bgzf",))


def test_unknown_empty_binary_and_missing_are_different_negatives(tmp_path: Path) -> None:
    unknown = sniff(_write(tmp_path / "notes.txt", b"dear diary\nnothing genomic here\n"))
    empty = sniff(_write(tmp_path / "empty.txt", b""))
    binary = sniff(_write(tmp_path / "blob", b"\x01\x02\x00\x00\x07" * 100))
    missing = sniff(tmp_path / "missing")
    assert [(r.outcome, r.reason) for r in (unknown, empty, binary, missing)] == [
        (SniffOutcome.UNKNOWN, "no registered format matched"),
        (SniffOutcome.UNKNOWN, "empty"),
        (SniffOutcome.UNKNOWN, "binary content"),
        (SniffOutcome.UNREADABLE, "not a file"),
    ]
    assert not any((unknown, empty, binary, missing))


def test_truncated_gzip_is_unreadable(tmp_path: Path) -> None:
    data = gzip.compress(VCF.encode() * 50)
    result = sniff(_write(tmp_path / "cut.gz", data[:20]))
    assert result.outcome is SniffOutcome.UNREADABLE
    assert result.layers == ("gzip",)


def test_utf16_crlf_and_bom_exports_are_read(tmp_path: Path) -> None:
    crlf = TWENTYTHREE.replace("\n", "\r\n")
    for name, data in (
        ("utf16", crlf.encode("utf-16")),
        ("bom", b"\xef\xbb\xbf" + crlf.encode()),
        ("latin1", TWENTYTHREE.replace("build 37", "build 37 \xe9").encode("latin-1")),
    ):
        result = sniff(_write(tmp_path / name, data))
        assert result.format == "23andme", name


def test_head_is_bounded_and_ends_on_a_whole_line(tmp_path: Path) -> None:
    """A header longer than the cap still sniffs from what fits, without reading the rest."""
    contigs = "".join(f"##contig=<ID=HLA-{i:06d},length=3000>\n" for i in range(HEAD_BYTE_CAP // 30))
    vcf = VCF.replace("##reference", contigs + "##reference", 1)
    result = sniff(_write(tmp_path / "long.vcf", gzip.compress(vcf.encode())))
    assert result.format == "vcf"
    assert result.build is None  # the ##reference line sits beyond the head cap: not seen, so not claimed


def test_sniff_is_deterministic(tmp_path: Path) -> None:
    path = tmp_path / "mixed.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("z.vcf", VCF)
        zf.writestr("a.txt", TWENTYTHREE)
        zf.writestr("m.txt", TWENTYTHREE)
    first, second = sniff(path), sniff(path)
    assert first == second
    assert isinstance(first, Sniffed)
