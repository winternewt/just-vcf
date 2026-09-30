"""Recognise an input file from its content, never from its name or label.

`sniff(path)` peels compression and archive layers by magic bytes (gzip and BGZF, bz2, xz, zip, tar,
nested), reads a bounded head of the payload and matches it against the format registry in order.
It never decompresses a whole file. The answer is a `Sniffed` result: truthy only when a format was
recognised, and otherwise saying which kind of negative it is.
"""

import bz2
import gzip
import io
import logging
import lzma
import re
import tarfile
import zipfile
from collections.abc import Callable, Iterator
from contextlib import ExitStack, contextmanager
from enum import Enum
from pathlib import Path
from typing import BinaryIO

from pydantic import BaseModel, ConfigDict

from just_vcf.formats import BCF_FORMAT, REGISTRY
from just_vcf.formats.spec import Family, FormatSpec

logger = logging.getLogger(__name__)

HEAD_BYTE_CAP = 1 << 20
"""Most decompressed bytes read from one payload. VCF headers with many `##contig` lines need room."""
HEAD_DATA_LINES = 50
"""Reading stops early once this many lines not starting with `#` are in the head."""
MAX_DEPTH = 4
"""Most nested layers followed (zip → tar → gzip → payload)."""
MAX_MEMBERS = 64
"""Most archive members sniffed; the rest are counted in a warning."""
NESTED_ZIP_CAP = 256 << 20
"""A zip inside another layer must be buffered to be read; larger ones are reported unreadable."""
TAR_SCAN_CAP = 1 << 30
"""A tar stream is read sequentially; listing stops after this many bytes of it."""
DATA_SHAPE_LINES = 5
"""A spec's `data_shape` is tried against this many leading data lines."""

_CHUNK = 64 << 10

# Payloads that are not variant data, recognised so they are reported as what they are rather than
# as unknown. Checked on raw bytes of the innermost layer.
_BINARY_KINDS: tuple[tuple[bytes, str], ...] = (
    (b"BAM\x01", "bam"),
    (b"CRAM", "cram"),
    (b"BAI\x01", "bai"),
    (b"TBI\x01", "tbi"),
    (b"CSI\x01", "csi"),
    (b"%PDF", "pdf"),
    (b"\x89PNG", "png"),
    (b"\xff\xd8\xff", "jpeg"),
    (b"GIF8", "gif"),
    (b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1", "ole2"),
    (b"Rar!", "rar"),
    (b"7z\xbc\xaf\x27\x1c", "7z"),
    (b"bplist00", "bplist"),
    (b"L\x00\x00\x00\x01\x14\x02\x00", "lnk"),
    (b"{\\rtf", "rtf"),
)

# Text payloads that are not variant data. Tried only after every format spec failed, on the first
# non-blank line of the head.
_TEXT_KINDS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"^@(HD|SQ|RG|PG)\t"), "sam"),
    (re.compile(r"^>"), "fasta"),
    (re.compile(r"^@\S"), "fastq"),
    (re.compile(r"^<!doctype html|^<html", re.IGNORECASE), "html"),
    (re.compile(r"^<\?xml"), "xml"),
    (re.compile(r"^0 HEAD\s*$"), "gedcom"),
    (re.compile(r"^Seq-submit ::= \{"), "genbank_asn1"),
    # Genetic-genealogy results that hold no positioned calls: Y-STR repeat counts, Y-SNP panels
    # by marker name, the rCRS/RSRS reference table, prose haplogroup reports, ancestry segments.
    (re.compile(r"^DYS393,DYS390,DYS19,"), "ftdna_ystr"),
    (re.compile(r'^"?SnpName"?,"?Positive"?,"?Test"?\s*$'), "ftdna_ysnp_panel"),
    (re.compile(r"^Position,rCRS,RSRS\s*$"), "mtdna_reference_table"),
    (re.compile(r"^\s*Type : mtDNA"), "genographic_report"),
    (re.compile(r"^Ancestry,Copy,Chromosome,Start Point,End Point"), "ancestry_segments"),
)

# Zip methods the standard library cannot read, by APPNOTE number. Deflate64 is what Windows
# Explorer writes for members over 2 GiB.
_ZIP_METHODS = {9: "Deflate64 (method 9)", 98: "PPMd (method 98)", 99: "AES encryption (method 99)"}

_ZIP_JUNK = re.compile(r"(^|/)(__MACOSX/|\._|\.DS_Store$|Thumbs\.db$)")


class SniffOutcome(str, Enum):
    """The kind of answer `sniff` gave."""

    DETECTED = "detected"
    """A registered format matched."""
    NOT_VARIANT_DATA = "not_variant_data"
    """Recognised as something that holds no variant calls (a PDF, a BAM, an index); see `kind`."""
    UNKNOWN = "unknown"
    """Readable, but no registered format matched."""
    UNREADABLE = "unreadable"
    """The container or payload could not be opened; see `reason`."""


class Sniffed(BaseModel):
    """What `sniff` found. Truthy only when a format was detected."""

    model_config = ConfigDict(frozen=True)

    path: Path
    outcome: SniffOutcome
    format: str | None = None
    family: Family | None = None
    kind: str | None = None
    """For `NOT_VARIANT_DATA`: what the payload is (`bam`, `pdf`, `ooxml`…)."""
    layers: tuple[str, ...] = ()
    """Container layers from the outside in: `zip`, `tar`, `gzip`, `bgzf`, `bz2`, `xz`."""
    member: str | None = None
    """The archive member the answer describes, innermost last, joined with `!`."""
    build: str | None = None
    """The build the header names. None means the header did not say, not GRCh37."""
    build_evidence: str | None = None
    evidence: tuple[str, ...] = ()
    """The detection patterns that matched."""
    also_matched: tuple[str, ...] = ()
    """Later registry formats that matched as well; the first match wins."""
    reason: str | None = None
    warnings: tuple[str, ...] = ()
    members: tuple["Sniffed", ...] = ()
    """For an archive with several data members: one result per member, largest first."""

    def __bool__(self) -> bool:
        return self.outcome is SniffOutcome.DETECTED


Opener = Callable[[ExitStack], BinaryIO]


def sniff(path: Path | str) -> Sniffed:
    """Recognise the format of `path` from its content."""
    path = Path(path)
    if not path.is_file():
        return Sniffed(path=path, outcome=SniffOutcome.UNREADABLE, reason="not a file")

    def open_raw(stack: ExitStack) -> BinaryIO:
        return stack.enter_context(path.open("rb"))

    try:
        result = _sniff_layer(path, open_raw, (), None, 0)
    except OSError as err:
        result = Sniffed(path=path, outcome=SniffOutcome.UNREADABLE, reason=f"{type(err).__name__}: {err}")
    logger.debug("sniffed %s: %s %s", path, result.outcome.value, result.format or result.kind or result.reason)
    return result


@contextmanager
def open_payload(result: Sniffed) -> Iterator[BinaryIO]:
    """The decompressed payload `result` describes, as a stream, for a loader to read in full.

    Replays `result.layers` from the outside in and opens `result.member` inside each archive, so
    the bytes are exactly the ones `sniff` matched. Only a detected result has a payload to open.
    """
    if not result:
        raise ValueError(f"{result.path}: nothing to open, sniff outcome is {result.outcome.value}")
    members = result.member.split("!") if result.member else []
    with ExitStack() as stack:
        handle: BinaryIO = stack.enter_context(result.path.open("rb"))
        for depth, layer in enumerate(result.layers):
            if layer in ("gzip", "bgzf"):
                handle = stack.enter_context(gzip.GzipFile(fileobj=handle))
            elif layer == "bz2":
                handle = stack.enter_context(bz2.BZ2File(handle))
            elif layer == "xz":
                handle = stack.enter_context(lzma.LZMAFile(handle))
            elif layer == "zip":
                handle = _open_zip_member(stack, handle, members.pop(0), nested=depth > 0)
            elif layer == "tar":
                handle = _open_tar_member(stack, handle, members.pop(0))
            else:
                raise ValueError(f"unknown layer {layer!r}")
        yield handle


def _open_zip_member(stack: ExitStack, handle: BinaryIO, name: str, nested: bool) -> BinaryIO:
    if nested:
        buffered = handle.read(NESTED_ZIP_CAP + 1)
        if len(buffered) > NESTED_ZIP_CAP:
            raise OSError("nested zip too large to buffer")
        handle = io.BytesIO(buffered)
    zf = stack.enter_context(zipfile.ZipFile(handle))
    info = zf.getinfo(name)
    if info.compress_type in _ZIP_METHODS:
        raise OSError(f"zip compression {_ZIP_METHODS[info.compress_type]} is not supported")
    return stack.enter_context(zf.open(info))


def _open_tar_member(stack: ExitStack, handle: BinaryIO, name: str) -> BinaryIO:
    tf = stack.enter_context(tarfile.open(fileobj=handle, mode="r|"))
    for info in tf:
        if info.name == name and info.isfile():
            extracted = tf.extractfile(info)
            if extracted is None:
                break
            return extracted
    raise OSError(f"tar member {name!r} not found")


def _sniff_layer(path: Path, opener: Opener, layers: tuple[str, ...], member: str | None, depth: int) -> Sniffed:
    base = {"path": path, "layers": layers, "member": member}
    if depth > MAX_DEPTH:
        return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason=f"more than {MAX_DEPTH} nested layers")

    try:
        with ExitStack() as stack:
            magic = _read_exactly(opener(stack), 512)
    except (OSError, EOFError, lzma.LZMAError) as err:
        return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason=f"{type(err).__name__}: {err}")

    try:
        if magic[:2] == b"\x1f\x8b":
            layer = "bgzf" if _is_bgzf(magic) else "gzip"
            return _sniff_layer(path, _wrap(opener, lambda f: gzip.GzipFile(fileobj=f)), (*layers, layer), member, depth + 1)
        if magic[:3] == b"BZh":
            return _sniff_layer(path, _wrap(opener, bz2.BZ2File), (*layers, "bz2"), member, depth + 1)
        if magic[:6] == b"\xfd7zXZ\x00":
            return _sniff_layer(path, _wrap(opener, lzma.LZMAFile), (*layers, "xz"), member, depth + 1)
        if magic[:4] == b"\x28\xb5\x2f\xfd":
            return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason="zstd compression is not supported")
        if magic[:4] in (b"PK\x03\x04", b"PK\x05\x06"):
            return _sniff_zip(path, opener, layers, member, depth)
        if magic[257:262] == b"ustar":
            return _sniff_tar(path, opener, layers, member, depth)
    except (OSError, EOFError, zipfile.BadZipFile, tarfile.TarError, lzma.LZMAError) as err:
        return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason=f"{type(err).__name__}: {err}")

    if magic[:4] == b"BCF\x02":
        spec = next(s for s in REGISTRY if s.id == BCF_FORMAT)
        return Sniffed(**base, outcome=SniffOutcome.DETECTED, format=spec.id, family=spec.family, evidence=("BCF\\2",))
    for prefix, kind in _BINARY_KINDS:
        if magic.startswith(prefix):
            return Sniffed(**base, outcome=SniffOutcome.NOT_VARIANT_DATA, kind=kind)

    try:
        with ExitStack() as stack:
            raw = _read_head(opener(stack))
    except (OSError, EOFError, lzma.LZMAError) as err:
        return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason=f"{type(err).__name__}: {err}")
    if not raw.strip():
        return Sniffed(**base, outcome=SniffOutcome.UNKNOWN, reason="empty")
    text = _decode(raw)
    if text is None:
        return Sniffed(**base, outcome=SniffOutcome.UNKNOWN, reason="binary content")
    # CRLF and CR-only exports exist; patterns see `\n` only, so `$` needs no `\r?`.
    return _match(base, text.replace("\r\n", "\n").replace("\r", "\n"))


def _match(base: dict, text: str) -> Sniffed:
    lines = text.splitlines()
    matched: list[tuple[FormatSpec, tuple[str, ...]]] = []
    for spec in REGISTRY:
        evidence = _spec_matches(spec, text, lines)
        if evidence is not None:
            matched.append((spec, evidence))
    if matched:
        spec, evidence = matched[0]
        build, build_evidence, warnings = _build_from_hints(spec, text)
        return Sniffed(
            **base,
            outcome=SniffOutcome.DETECTED,
            format=spec.id,
            family=spec.family,
            evidence=evidence,
            also_matched=tuple(s.id for s, _ in matched[1:]),
            build=build,
            build_evidence=build_evidence,
            warnings=warnings,
        )
    first = next((line for line in lines if line.strip()), "")
    for regex, kind in _TEXT_KINDS:
        if regex.search(first):
            return Sniffed(**base, outcome=SniffOutcome.NOT_VARIANT_DATA, kind=kind)
    return Sniffed(**base, outcome=SniffOutcome.UNKNOWN, reason="no registered format matched")


def _spec_matches(spec: FormatSpec, text: str, lines: list[str]) -> tuple[str, ...] | None:
    """The evidence for `spec` matching, or None when it does not."""
    detection = spec.detection
    if any(re.search(p, text, re.MULTILINE) for p in detection.excludes):
        return None
    if not all(re.search(p, text, re.MULTILINE) for p in detection.patterns):
        return None
    evidence = detection.patterns
    if detection.data_shape is not None:
        prefix = spec.comment_prefix
        data = [line for line in lines if line.strip() and not (prefix and line.startswith(prefix))]
        if not any(re.fullmatch(detection.data_shape, line) for line in data[:DATA_SHAPE_LINES]):
            return None
        evidence = (*evidence, detection.data_shape)
    return evidence


def _build_from_hints(spec: FormatSpec, text: str) -> tuple[str | None, str | None, tuple[str, ...]]:
    """The build the header names. Disagreeing hints make the answer unknown, with a warning."""
    found: list[tuple[str, str]] = []
    for hint in spec.build_hints:
        m = re.search(hint.pattern, text, re.MULTILINE)
        if m:
            found.append((hint.build, m.group(0).strip()))
    builds = sorted({build for build, _ in found})
    if not builds:
        return None, None, ()
    if len(builds) > 1:
        return None, None, (f"header names conflicting builds: {', '.join(builds)}",)
    return builds[0], found[0][1], ()


def _sniff_zip(path: Path, opener: Opener, layers: tuple[str, ...], member: str | None, depth: int) -> Sniffed:
    base = {"path": path, "layers": (*layers, "zip"), "member": member}
    with ExitStack() as stack:
        handle = opener(stack)
        if not layers:
            zf = stack.enter_context(zipfile.ZipFile(handle))
        else:
            buffered = handle.read(NESTED_ZIP_CAP + 1)
            if len(buffered) > NESTED_ZIP_CAP:
                return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason="nested zip too large to buffer")
            zf = stack.enter_context(zipfile.ZipFile(io.BytesIO(buffered)))
        names = zf.namelist()
        if "[Content_Types].xml" in names:
            return Sniffed(**base, outcome=SniffOutcome.NOT_VARIANT_DATA, kind="ooxml")
        if "mimetype" in names:
            return Sniffed(**base, outcome=SniffOutcome.NOT_VARIANT_DATA, kind="odf")
        if any(n.startswith("Index/") and n.endswith(".iwa") for n in names) or "index.xml" in names:
            return Sniffed(**base, outcome=SniffOutcome.NOT_VARIANT_DATA, kind="iwork")
        infos = [i for i in zf.infolist() if not i.is_dir() and not _ZIP_JUNK.search(i.filename)]
        if not infos:
            return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason="archive has no data members")
        infos.sort(key=lambda i: (-i.file_size, i.filename))

        def member_opener(info: zipfile.ZipInfo) -> Opener:
            return lambda s: s.enter_context(zf.open(info))

        results = []
        for info in infos[:MAX_MEMBERS]:
            name = info.filename if member is None else f"{member}!{info.filename}"
            try:
                results.append(_sniff_layer(path, member_opener(info), base["layers"], name, depth + 1))
            except NotImplementedError:
                method = _ZIP_METHODS.get(info.compress_type, f"method {info.compress_type}")
                reason = f"zip compression {method} is not supported"
                results.append(Sniffed(**{**base, "member": name}, outcome=SniffOutcome.UNREADABLE, reason=reason))
        return _aggregate(base, results, len(infos))


def _sniff_tar(path: Path, opener: Opener, layers: tuple[str, ...], member: str | None, depth: int) -> Sniffed:
    base = {"path": path, "layers": (*layers, "tar"), "member": member}
    heads: list[tuple[str, int, bytes]] = []
    total = 0
    truncated = False
    with ExitStack() as stack:
        with tarfile.open(fileobj=opener(stack), mode="r|") as tf:
            for info in tf:
                if tf.offset > TAR_SCAN_CAP:
                    truncated = True
                    break
                if not info.isfile() or _ZIP_JUNK.search(info.name):
                    continue
                total += 1
                if len(heads) < MAX_MEMBERS:
                    extracted = tf.extractfile(info)
                    heads.append((info.name, info.size, extracted.read(HEAD_BYTE_CAP) if extracted else b""))
    if not heads:
        return Sniffed(**base, outcome=SniffOutcome.UNREADABLE, reason="archive has no data members")
    heads.sort(key=lambda h: (-h[1], h[0]))

    def head_opener(data: bytes) -> Opener:
        return lambda s: io.BytesIO(data)

    results = []
    for name, _, data in heads:
        full = name if member is None else f"{member}!{name}"
        results.append(_sniff_layer(path, head_opener(data), base["layers"], full, depth + 1))
    result = _aggregate(base, results, total)
    if truncated:
        note = f"tar listing stopped after {TAR_SCAN_CAP >> 20} MiB; later members not sniffed"
        result = result.model_copy(update={"warnings": (*result.warnings, note)})
    return result


def _aggregate(base: dict, results: list[Sniffed], total: int) -> Sniffed:
    """One answer for an archive: its only member's, or the largest detected member's with the rest attached."""
    if len(results) == 1 and total == 1:
        return results[0]
    warnings: list[str] = []
    if total > len(results):
        warnings.append(f"{total - len(results)} archive member(s) beyond the first {MAX_MEMBERS} not sniffed")
    detected = [r for r in results if r]
    formats = sorted({r.format for r in detected if r.format})
    if len(formats) > 1:
        warnings.append(f"members hold different formats: {', '.join(formats)}")
    lead = detected[0] if detected else results[0]
    return lead.model_copy(update={"members": tuple(results), "warnings": (*lead.warnings, *warnings)})


def _wrap(opener: Opener, decompressor: Callable[[BinaryIO], BinaryIO]) -> Opener:
    def wrapped(stack: ExitStack) -> BinaryIO:
        return stack.enter_context(decompressor(opener(stack)))

    return wrapped


def _is_bgzf(magic: bytes) -> bool:
    """BGZF is gzip with FEXTRA set and a `BC` subfield (SAM spec §4.1)."""
    return len(magic) >= 18 and magic[3] & 0x04 != 0 and magic[12:14] == b"BC"


def _read_exactly(handle: BinaryIO, size: int) -> bytes:
    """Up to `size` bytes; a decompressing stream may return fewer per read."""
    parts: list[bytes] = []
    remaining = size
    while remaining > 0:
        chunk = handle.read(remaining)
        if not chunk:
            break
        parts.append(chunk)
        remaining -= len(chunk)
    return b"".join(parts)


def _read_head(handle: BinaryIO) -> bytes:
    """Whole lines from the start, until `HEAD_DATA_LINES` non-`#` lines or `HEAD_BYTE_CAP` bytes."""
    buf = bytearray()
    at_eof = False
    while len(buf) < HEAD_BYTE_CAP:
        chunk = handle.read(_CHUNK)
        if not chunk:
            at_eof = True
            break
        buf += chunk
        data_lines = sum(1 for line in buf.split(b"\n")[:-1] if line and not line.startswith(b"#"))
        if data_lines >= HEAD_DATA_LINES:
            break
    head = bytes(buf[:HEAD_BYTE_CAP])
    if not at_eof and b"\n" in head:
        head = head[: head.rindex(b"\n") + 1]
    return head


def _decode(raw: bytes) -> str | None:
    """Text from a head, or None when it is binary. BOM-aware; latin-1 when UTF-8 fails."""
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16", errors="replace")
    if b"\x00" in raw[:8192]:
        return None
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as err:
        if err.start >= len(raw) - 4:
            return raw[: err.start].decode("utf-8-sig")
        return raw.decode("latin-1")
