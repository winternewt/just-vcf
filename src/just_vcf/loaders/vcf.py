"""Load one sample of a VCF or gVCF into `calls` and `blocks` rows (docs/SCHEMA.md).

polars-bio reads the records. This module reads the header itself, because the loader needs what
polars-bio does not hand back: which namespace a key was declared in (polars-bio names a column `DP`
or `fmt_DP` depending on what else the header declares, P5), the file's VCF version, the `##contig`
lines, and the `##SAMPLE` and `##PEDIGREE` lines.

A header polars-bio would misread is rewritten into a temporary copy first: meta lines are put
`##fileformat` first (P2), and a FORMAT key the rows use without declaring it gets the declaration
VCF 4.5 reserves for it, since polars-bio drops undeclared keys, GT included (P1). Every rewrite is
reported in the warnings.

Coordinates stay 1-based (`use_zero_based=False` is passed, never assumed). `end_pos` is computed
here from VCF 4.5's END rule, not taken from polars-bio, whose `end` is one short for symbolic SVs
(P3). polars-bio drops a GT's leading phase indicator (P4), so in a VCF 4.4+ file the first
haplotype's `phased` entry is withheld rather than reported from the rewritten text, except where a
PSL name says it is phased (`_phase_sets`).
"""

import hashlib
import io
import logging
import re
import shutil
import struct
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import BinaryIO

import polars as pl
import polars_bio as pb
from polars_bio.sql import SQL

from just_vcf import cache as cache_index
from just_vcf.bgzf import BgzfWriter
from just_vcf.links import CALLER_SPLIT  # noqa: F401 (kept importable here; defined in just_vcf.links)
from just_vcf.schema import BLOCKS, CALLS, EMPTY_ALLELE, BlockKind, EvidenceBasis, PhaseSetKind, RefCheck, SplitSignal, canonical_chrom_expr, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset({"vcf", "bcf", "vcf_no_meta", "vcf_cgi", "vcf_freebayes_varscan2"})
"""Format ids this loader reads. The dialects with their own ids join as their quirks are handled."""

BCF = "bcf"
"""Binary VCF: the header is read from the binary layout here, the records by polars-bio."""

FREEBAYES_VARSCAN2 = "vcf_freebayes_varscan2"
"""PGP-HMS's vcf-merge of freeBayes and VarScan2: one facet per caller, split on INFO/SF."""

SELECTION_INFO = ("SF",)
"""INFO keys read only to select a part's rows. The structural copy made for a refused header keeps
them, or the retry could not split the file."""

_MERGE_CALLERS = {"FB": "freeBayes", "VS": "VarScan2"}
"""The caller tokens PGP's pipeline puts in the merged inputs' names (`chr1.FB.vcf.gz`), as the
merge's `##sourceFiles` line records them."""

FORMAT_KEYS: dict[str, tuple[str, str]] = {
    "GT": ("1", "String"),
    "GQ": ("1", "Integer"),
    "DP": ("1", "Integer"),
    "AD": ("R", "Integer"),
    "PL": ("G", "Integer"),
    "GL": ("G", "Float"),
    "GP": ("G", "Float"),
    "PP": ("G", "Integer"),
    "DS": ("A", "Float"),
    "FT": ("1", "String"),
    "PS": ("1", "Integer"),
    "PSL": ("P", "String"),
    "PSO": ("P", "Integer"),
    "PSQ": ("P", "Integer"),
    "LEN": ("1", "Integer"),
    "MIN_DP": ("1", "Integer"),
    "CN": ("1", "Float"),
}
"""FORMAT keys the schema carries, with the Number and Type VCF 4.5 reserves for them. `P` (one value
per allele in GT, VCF 4.4+) is refused by polars-bio's reader (P10), so a copy declares those keys
`Number=.` (`READ_AS_ANY_NUMBER`)."""

READ_AS_ANY_NUMBER = frozenset({"P"})
"""FORMAT Numbers polars-bio refuses (P10) that the loader reads through a copy declaring `Number=.`:
the values are the same lists, and the length is checked against GT here."""

VECTOR_NUMBERS: dict[str, str] = {"AD": "R", "PL": "G", "GL": "G", "GP": "G", "PP": "G"}
"""Carried vectors that are parallel to the alleles (R) or the genotypes (G). A file declaring one of
them with a fixed count is using the name for something else: VarScan2's AD (Number=1) is the ALT
depth only, CGI's AD (Number=2) counts reads per GT slot, not per allele. Those are not carried."""

INFO_KEYS: dict[str, tuple[str, str]] = {
    "END": ("1", "Integer"),
    "SVLEN": ("A", "Integer"),
    "SVTYPE": ("1", "String"),
    "CIPOS": (".", "Integer"),
    "CIEND": (".", "Integer"),
    "CN": ("A", "Float"),
}
"""INFO keys the loader reads, with the Number and Type VCF 4.5 gives them: END and SVLEN for
`end_pos`, and the structural-variant keys the SV columns carry (RM282). VCF 4.5 deprecates SVTYPE
as redundant with ALT; it is read because callers write a type ALT does not say (DRAGEN's `CNV` on a
`<DEL>`, Canvas's `LOH`)."""

V43_RESERVED_TYPES: dict[str, str] = {"FORMAT/CN": "Integer", "INFO/CN": "Integer"}
"""The carried keys whose Type VCF 4.3 reserves differently from VCF 4.5 (CN is Float from 4.4 on).
polars-bio's reader (noodles) refuses a VCF 4.3 header declaring one of them with another Type, so
such a key is not carried (`refused_types`), and the loader declares an undeclared one so."""

_LIST_INFO = ("SVLEN", "CIPOS", "CIEND", "CN")
"""INFO keys read as lists. One declared `Number=1` (VCF 4.2's INFO CN, an older caller's SVLEN) is
read as the one-value list of a one-ALT record."""

SPLIT_MARKER_ID = "OLD_MULTIALLELIC"
"""The INFO key vt decompose and Glow write on a record split off a multi-allelic site."""

SPLIT_INFO_KEY = "JUST_VCF_SPLIT"
"""The INFO key `write_vcf` puts on a row flagged as one record of a split multi-allelic site
(`normalize`'s split rows, and rows RM270 flagged). Read back as a marker, whatever its value."""

BCFTOOLS_OLD_REC = "Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX"
"""The Description bcftools norm writes for its `--old-rec-tag` key, whatever the key is named
(vcfnorm.c). Only a split record gets the fifth, `USED_ALT_IDX`, field; a realigned one has four."""

CGI_NOCALL = "<CGA_NOCALL>"
"""Complete Genomics' no-called span: INFO END gives its extent. It goes to `blocks` as `no_call`."""

DASH_DELETION = "-"
"""An ALT some uploads write for a deletion (F12): the REF bases at POS deleted, with no anchor base.
VCF has no such allele. The loader reads it as the empty allele (`schema.EMPTY_ALLELE`), the one
spelling of a zero-length allele, since anchoring it needs the reference base before POS; `normalize`
anchors it, and `write_vcf` refuses it until then."""

_VCF_ALLELE_OR_DASH = r"^([ACGTNacgtn]+|\*|<[^<>]+>|\.[ACGTNacgtn]+|[ACGTNacgtn]+\.|-)$|\[|\]"
"""An ALT VCF 4.5 allows (bases, `*`, a symbolic allele, a breakend), or the dash deletion F12 counts
apart. Anything else is counted as no allele (F15)."""

CGI_CNV_WINDOW = "<CGA_CNVWIN>"
"""Complete Genomics' 2 kb copy-number window: coverage, not a genotype. Left out of the frames."""

_SV_LENGTH_TYPES = ("DEL", "DUP", "INV", "CNV", r"CN\d+")
"""Symbolic allele types whose extent is `SVLEN` bases after POS, or END. INS and breakends end at
POS. `<CN0>`, `<CN2>`… are copy-number alleles (Canvas, the 1000 Genomes SV release): a segment of
this contig to END, as `<CNV>` is."""

_SV_PATTERN = rf"^<({'|'.join(_SV_LENGTH_TYPES)})(:[^>]*)?>$"
"""An ALT of one of `_SV_LENGTH_TYPES`, subtypes included (`<DUP:TANDEM>`, `<DEL:ME:ALU>`)."""

_REGISTRATION_LOCK = Lock()
_REGISTERED: set[str] = set()

REWRITE_VERSION = 9
"""Bumped whenever `readable_path` changes what it writes, so a cached copy from older code is not reused.
7: every copy has a manifest beside it (RM268), so a copy without one is from an older version.
8: the SV keys (SVTYPE, CIPOS, CIEND, CN) are carried, so a copy declares them where the rows use
them undeclared, and a structural copy keeps them (RM282); an older copy lacks the declarations,
and polars-bio panics on a requested key its header does not declare.
9: a FORMAT key declared `Number=P` (PSL, PSO, PSQ) is declared `Number=.` in the copy (P10), and the
phase-set-list keys are carried, so a structural copy keeps them."""

HASHED_COLUMN_LINE = "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
"""The VCF column line written with two `#` (F7). Read as the column line only where no `#CHROM` line
came before it; the copy polars-bio reads writes it with one."""

_HEADER_READ_CAP = 64 << 20
_DATA_LINES_CHECKED = 1000


class MultiSampleError(ValueError):
    """The file holds more than one sample and none was chosen."""


class NoSampleError(ValueError):
    """The file holds no sample (a sites-only VCF): it has sites, but no person's genotype to load."""


class HeaderRefused(ValueError):
    """polars-bio refused the file: its header (P8), or a field value it could not parse (P7)."""


_STRUCTURAL_META = ("##fileformat=", "##FILTER=", "##ALT=", "##contig=")
"""Meta lines kept in a structural-only copy, besides the INFO and FORMAT keys the frames carry.
Everything else is read by `read_header`, not by polars-bio."""


def _keep_in_structural_copy(line: str, refused: dict[str, str]) -> bool:
    if line.startswith(_STRUCTURAL_META):
        return True
    declared = re.match(r"^##(INFO|FORMAT)=<(.*)>$", line)
    if not declared:
        return False
    kind, body = declared.groups()
    attrs = _attrs(body)
    if f"{kind}/{attrs.get('ID')}" in refused:
        return False
    if kind == "INFO":
        return attrs.get("ID") in INFO_KEYS or attrs.get("ID") in SELECTION_INFO or _split_marker_kind(attrs) is not None
    return attrs.get("ID") in FORMAT_KEYS


def _attrs(body: str) -> dict[str, str]:
    """The key=value pairs of a structured meta line, quotes left on."""
    return dict(re.findall(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,]*)', body))


def _split_marker_kind(attrs: dict[str, str]) -> str | None:
    """`vt`, `just_vcf` or `bcftools` when a declared INFO key marks records split off a multi-allelic site."""
    if attrs.get("ID") == SPLIT_MARKER_ID:
        return "vt"
    if attrs.get("ID") == SPLIT_INFO_KEY:
        return "just_vcf"
    if attrs.get("Description", "").strip('"') == BCFTOOLS_OLD_REC:
        return "bcftools"
    return None


@dataclass
class Header:
    """What the loader reads from a VCF header."""

    meta: list[str]
    columns: list[str]
    samples: list[str]
    version: str | None
    info: dict[str, str] = field(default_factory=dict)
    """Declared INFO key → Number."""
    format: dict[str, str] = field(default_factory=dict)
    """Declared FORMAT key → Number."""
    types: dict[str, str] = field(default_factory=dict)
    """Declared key, namespaced (`FORMAT/CN`) → Type."""
    contigs: list[dict] = field(default_factory=list)
    sample_meta: list[str] = field(default_factory=list)
    format_first: bool = True
    hashed_column_line: bool = False
    resaved: bool = False
    """The file was re-saved by a spreadsheet (F7): read through `undo_resave`, so polars-bio needs the
    rewritten copy. Set by the caller that chose to read it so."""
    """The column line was written `##CHROM` (F7: Veritas SV and CNV uploads), so polars-bio cannot
    read the file as it is."""
    split_markers: dict[str, str] = field(default_factory=dict)
    alt_ids: set[str] = field(default_factory=set)
    """Symbolic alleles the header declares with `##ALT` (`CGA_NOCALL`, `NON_REF`, `DEL`…)."""
    """Declared INFO keys that mark split records → `vt`, `just_vcf` or `bcftools` (`_split_marker_kind`)."""

    @property
    def version_tuple(self) -> tuple[int, int] | None:
        match = re.match(r"VCFv(\d+)\.(\d+)", self.version or "")
        return (int(match.group(1)), int(match.group(2))) if match else None


def read_header(handle: BinaryIO) -> tuple[Header, list[bytes]]:
    """The header, and up to `_DATA_LINES_CHECKED` data lines after it."""
    meta: list[str] = []
    columns: list[str] = []
    data: list[bytes] = []
    read = 0
    hashed: str | None = None
    for raw in handle:
        read += len(raw)
        line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
        if line.startswith("##"):
            if hashed is None and line.startswith(HASHED_COLUMN_LINE):
                hashed = line
            meta.append(line)
        elif line.startswith("#"):
            columns = line[1:].split("\t")
        elif line:
            data.append(raw)
            if len(data) >= _DATA_LINES_CHECKED:
                break
        if read > _HEADER_READ_CAP and not columns:
            raise ValueError("no #CHROM line in the first 64 MiB")
    from_hashed = not columns and hashed is not None
    if from_hashed:
        # F7: the column line with one `#` too many, and no real one. It names VCF's fixed columns
        # exactly, so it is read as the column line and the rows after it as the records.
        columns = hashed[2:].split("\t")
        meta.remove(hashed)
    if not columns:
        raise ValueError("no #CHROM line")
    header = _parse_header(meta, columns)
    header.hashed_column_line = from_hashed
    return header, data


RESAVE_HEAD_BYTES = 64 << 10
"""How much of a payload `is_resaved` looks at."""

_QUOTED_META = re.compile(rb'(\A|[\r\n])"##')


def is_resaved(head: bytes) -> bool:
    """Whether a VCF payload's head shows a spreadsheet's re-save (F7): a meta line wrapped in quotes
    (`"##FILTER=<…""Low quality"">"`), or lines ended by a bare CR with no LF at all."""
    return bool(_QUOTED_META.search(head)) or (b"\r" in head and b"\n" not in head)


class _Lines(io.RawIOBase):
    """A readable byte stream over an iterator of byte chunks."""

    def __init__(self, chunks: Iterator[bytes]) -> None:
        self._chunks = chunks
        self._left = b""

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        while not self._left:
            chunk = next(self._chunks, None)
            if chunk is None:
                return 0
            self._left = chunk
        size = min(len(buffer), len(self._left))
        buffer[:size] = self._left[:size]
        self._left = self._left[size:]
        return size


def undo_resave(handle: BinaryIO) -> BinaryIO:
    """A spreadsheet-re-saved VCF payload as the tab-separated text it was saved from (F7).

    The re-save is CSV quoting applied to a tab file: a field holding a comma or a quote is wrapped in
    quotes with its own quotes doubled (an INFO with `EFF=…,…`, a sample field whose AD is `12,0`), meta
    lines are padded with empty fields to the widest row, and lines may end in a bare CR. Universal
    newlines undo the CRs; each line is cut at its tabs and a field wrapped in quotes is unwrapped, its
    doubled quotes made single (`_unquote`); trailing empty fields are dropped from `##` lines only.
    What the spreadsheet did to a value's text (a QUAL written `6980.9` for `6980.90`) cannot be undone
    or seen.

    Not `csv`: its reader refuses a field over `csv.field_size_limit()` (128 KiB by default), and a VCF
    field (SnpEff's EFF on a busy site) has no limit. The limit is process-wide, so raising it would
    change every other reader in the process, and a set-and-restore cannot be scoped to a reader that
    is consumed lazily. A VCF field holds no tab and no line break, so the quoting needs no state
    across fields or lines.
    """
    text = io.TextIOWrapper(handle, encoding="utf-8", errors="replace", newline=None)

    def lines() -> Iterator[bytes]:
        for line in text:
            row = [_unquote(field) for field in line.rstrip("\n").split("\t")]
            if row[0].startswith("##"):
                while len(row) > 1 and row[-1] == "":
                    row.pop()
            yield ("\t".join(row) + "\n").encode()

    return io.BufferedReader(_Lines(lines()), 1 << 20)


def _unquote(field: str) -> str:
    """A spreadsheet's CSV quoting taken off one field (outer quotes dropped, doubled quotes made
    single); a field not wrapped in quotes as it is."""
    if len(field) >= 2 and field[0] == '"' and field[-1] == '"':
        return field[1:-1].replace('""', '"')
    return field


def read_bcf_header(handle: BinaryIO) -> tuple[Header, bytes | None]:
    """A BCF's header, and its first record's CHROM and 1-based POS as a `CHROM\tPOS` line.

    `handle` is the decompressed payload (`open_payload` peels the BGZF layer). The layout is BCF 2.x:
    the magic `BCF` with major version 2 and a minor version byte, a little-endian u32 header length, the VCF header text NUL-terminated, then
    records whose shared part starts with `l_shared`, `l_indiv`, the contig's dictionary index and the
    0-based POS. That position is what `read_records` checks polars-bio's coordinates against.
    """
    magic = handle.read(5)
    if magic[:4] != b"BCF\x02":
        raise ValueError("not a BCF 2 payload")
    (length,) = struct.unpack("<I", handle.read(4))
    text = handle.read(length).rstrip(b"\0").decode("utf-8", errors="replace")
    meta = [line for line in text.splitlines() if line.startswith("##")]
    columns = next((line[1:].split("\t") for line in text.splitlines() if line.startswith("#") and not line.startswith("##")), [])
    if not columns:
        raise ValueError("no #CHROM line in the BCF header")
    header = _parse_header(meta, columns)
    shared = handle.read(16)
    if len(shared) < 16:
        return header, None
    _, _, contig, pos = struct.unpack("<IIii", shared)
    names = _contig_dictionary(meta)
    if not 0 <= contig < len(names) or names[contig] is None:
        raise ValueError(f"the first BCF record names contig index {contig}, which the header does not declare")
    return header, f"{names[contig]}\t{pos + 1}".encode()


def _contig_dictionary(meta: list[str]) -> list[str | None]:
    """BCF's contig dictionary: `##contig` lines in order, or by their `IDX` where the header gives one."""
    entries = []
    for line in meta:
        declared = re.match(r"^##contig=<(.*)>$", line)
        if declared and "ID" in (attrs := _attrs(declared.group(1))):
            entries.append((attrs["ID"], int(attrs["IDX"]) if attrs.get("IDX", "").isdigit() else None))
    if not any(idx is not None for _, idx in entries):
        return [name for name, _ in entries]
    names: list[str | None] = [None] * (max(idx for _, idx in entries if idx is not None) + 1)
    for name, idx in entries:
        if idx is not None:
            names[idx] = name
    return names


def _parse_header(meta: list[str], columns: list[str]) -> Header:
    header = Header(
        meta=meta,
        columns=columns,
        samples=columns[9:],
        version=next((m.split("=", 1)[1] for m in meta if m.startswith("##fileformat=")), None),
        format_first=bool(meta) and meta[0].startswith("##fileformat="),
    )
    for line in meta:
        structured = re.match(r"^##(\w+)=<(.*)>$", line)
        if not structured:
            continue
        kind, body = structured.groups()
        attrs = _attrs(body)
        if kind in ("INFO", "FORMAT") and "ID" in attrs:
            header.types[f"{kind}/{attrs['ID']}"] = attrs.get("Type", "")
        if kind == "INFO" and "ID" in attrs:
            header.info[attrs["ID"]] = attrs.get("Number", ".")
            if marker := _split_marker_kind(attrs):
                header.split_markers[attrs["ID"]] = marker
        elif kind == "FORMAT" and "ID" in attrs:
            header.format[attrs["ID"]] = attrs.get("Number", ".")
        elif kind == "ALT" and "ID" in attrs:
            header.alt_ids.add(attrs["ID"])
        elif kind == "contig" and "ID" in attrs:
            length = attrs.get("length")
            header.contigs.append(
                {"name": attrs["ID"], "length": int(length) if length and length.isdigit() else None, "md5": attrs.get("md5")}
            )
        elif kind in ("SAMPLE", "PEDIGREE"):
            header.sample_meta.append(line)
    return header


def row_keys(data: list[bytes]) -> tuple[set[str], set[str]]:
    """The INFO and FORMAT keys the rows in `data` use."""
    info: set[str] = set()
    fmt: set[str] = set()
    for raw in data:
        fields = raw.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")
        if len(fields) > 7 and fields[7] not in ("", "."):
            info.update(item.split("=", 1)[0] for item in fields[7].split(";") if item)
        if len(fields) > 8:
            fmt.update(fields[8].split(":"))
    return info, fmt


def undeclared_keys(header: Header, data: list[bytes]) -> list[str]:
    """Carried keys the rows use without the header declaring them, namespaced (`INFO/END`, `FORMAT/GT`)."""
    info, fmt = row_keys(data)
    return sorted(f"INFO/{key}" for key in info if key in INFO_KEYS and key not in header.info) + sorted(
        f"FORMAT/{key}" for key in fmt if key in FORMAT_KEYS and key not in header.format
    )


def declare(header: Header, missing: list[str]) -> None:
    """Record the loader's declarations of `missing` keys in `header`, as `readable_path` writes them."""
    for name in missing:
        namespace, key = name.split("/", 1)
        header.types[name] = declared_type(header, name)
        if namespace == "INFO":
            header.info[key] = INFO_KEYS[key][0]
        else:
            header.format[key] = FORMAT_KEYS[key][0]


def declared_type(header: Header, name: str) -> str:
    """The Type the loader declares an undeclared carried key with (`INFO/CN`): VCF 4.5's, except where
    the file's VCF 4.3 reserves another one (`V43_RESERVED_TYPES`)."""
    namespace, key = name.split("/", 1)
    if header.version_tuple == (4, 3) and name in V43_RESERVED_TYPES:
        return V43_RESERVED_TYPES[name]
    return (INFO_KEYS if namespace == "INFO" else FORMAT_KEYS)[key][1]


def refused_types(header: Header) -> dict[str, str]:
    """Carried keys a VCF 4.3 header declares with a Type its reserved definition does not have
    (`V43_RESERVED_TYPES`) → that Type. polars-bio refuses the header, so they are not carried."""
    if header.version_tuple != (4, 3):
        return {}
    return {
        name: header.types[name]
        for name, reserved in V43_RESERVED_TYPES.items()
        if name in header.types and header.types[name] != reserved
    }


@dataclass
class Part:
    """One facet's share of a VCF: which sample column, which rows, and what makes its id."""

    sample: str
    options: dict
    caller: str | None = None
    split: bool = False
    """A caller's share of a merged file (`relation = caller_split`)."""
    sf: int | None = None
    """Keep the rows whose INFO/SF lists this source index."""
    sf_alone: bool = False
    """Keep them only where SF lists no other index (the merge kept another caller's genotype there)."""
    merge_indices: tuple[int, ...] = ()
    """Every part's source index, set on a merge's first part only: the records none of the parts
    takes are counted once per file (`unsourced`)."""

    @property
    def info(self) -> tuple[str, ...]:
        return SELECTION_INFO if self.sf is not None else ()

    def select(self, records: pl.LazyFrame) -> pl.LazyFrame:
        """`records` restricted to this part's rows. An `f` after an index marks the record as
        filtered in that source; the record is still that source's."""
        if self.sf is None:
            return records
        indices = pl.col("info_SF").list.eval(pl.element().str.strip_chars_end("f"))
        keep = indices.list.contains(str(self.sf))
        if self.sf_alone:
            keep = keep & indices.list.eval(pl.element() != str(self.sf)).list.any().not_()
        return records.filter(keep.fill_null(False))

    def unsourced(self, records: pl.LazyFrame) -> pl.LazyFrame | None:
        """The records of a merge that no part takes, which land in no facet: INFO/SF absent or empty,
        or naming no `##sourceFiles` index. None on every part but the merge's first, so the file's
        count is taken once (`scan_stats`)."""
        if not self.merge_indices:
            return None
        indices = pl.col("info_SF").list.eval(pl.element().str.strip_chars_end("f"))
        claimed = indices.list.eval(pl.element().is_in([str(index) for index in self.merge_indices])).list.any()
        return records.filter(~claimed.fill_null(False))


def caller_parts(header: Header, data: list[bytes]) -> list[Part]:
    """The parts of a freeBayes+VarScan2 merge: one per source index of its `##sourceFiles` lines.

    The caller of an index is read from the merged input's name on those lines (`.FB.`, `.VS.`);
    lines that disagree, or a name without a token, leave it null. Two layouts (FORMATS.md): with one
    data column per caller, each index reads its own column, in `##sourceFiles` order; with one data
    column, it holds the first index's genotype wherever that caller found the site, so a later index
    keeps only the rows where it was the sole caller, and says what it lost.
    """
    callers: dict[int, set[str | None]] = {}
    for line in header.meta:
        if re.match(r"^##sourceFiles_\d{8}\.\d+=", line):
            for index, name in re.findall(r"(\d+):([^,]+)", line.split("=", 1)[1]):
                token = next((t for t in _MERGE_CALLERS if re.search(rf"[._]{t}[._]", name)), None)
                callers.setdefault(int(index), set()).add(_MERGE_CALLERS.get(token) if token else None)
    indices = sorted(callers) or [0, 1]
    columns = len(header.samples)
    per_caller = columns >= len(indices) and bool(data) and all(len(raw.split(b"\t")) >= 9 + columns for raw in data)
    parts = []
    for position, index in enumerate(indices):
        named = callers.get(index, {None})
        caller = next(iter(named)) if len(named) == 1 else None
        lone = not per_caller and position > 0
        parts.append(
            Part(
                sample=header.samples[position] if per_caller else header.samples[0],
                options={"caller_index": index},
                caller=caller,
                split=True,
                sf=index,
                sf_alone=lone,
                merge_indices=tuple(indices) if position == 0 else (),
            )
        )
    return parts


def any_number_keys(header: Header) -> list[str]:
    """FORMAT keys the header declares with a Number polars-bio refuses and the loader reads through a
    copy declaring `Number=.` (`READ_AS_ANY_NUMBER`, P10). Any such key, carried or not: polars-bio
    refuses the whole header for one."""
    return sorted(key for key, number in header.format.items() if number in READ_AS_ANY_NUMBER)


def _readable_number(line: str) -> str:
    """A `##FORMAT` line with a Number polars-bio refuses (`READ_AS_ANY_NUMBER`) declared `Number=.`."""
    if not line.startswith("##FORMAT=<"):
        return line
    numbers = "|".join(re.escape(number) for number in sorted(READ_AS_ANY_NUMBER))
    return re.sub(rf"(?<=[<,])Number=(?:{numbers})(?=[,>])", "Number=.", line)


def unfit_keys(header: Header) -> dict[str, str]:
    """Carried vector keys the header declares with a Number they cannot be read under → that Number."""
    return {
        key: header.format[key]
        for key, number in VECTOR_NUMBERS.items()
        if key in header.format and header.format[key] not in (number, ".")
    }


def readable_path(
    path: Path,
    layers: tuple[str, ...],
    opener,
    header: Header,
    missing: list[str],
    cache: Path,
    origin: cache_index.Origin,
    structural_only: bool = False,
) -> Path:
    """A path polars-bio can read: the file itself when possible, else a rewritten copy under `cache`.

    The copy is kept, because the frames `load` returns are lazy and read it when collected. It is
    named after the source's digest and reused when complete, and written as BGZF: about a quarter of
    the plain size, and polars-bio reads it as fast as plain text. `structural_only` keeps only the meta
    lines polars-bio needs and the INFO and FORMAT keys the frames carry, for a file it refused (P7, P8):
    polars-bio skips a key its header does not declare (P1), so a malformed value in a key the frames
    do not carry is never parsed.

    `cache` is the entry (`vcf/<key>`); it is written, or reused, under its lock, with `origin` in its
    manifest (`just_vcf.cache.writing`).
    """
    as_is = (
        header.format_first
        and not header.hashed_column_line
        and not header.resaved
        and not missing
        and not structural_only
        and not any_number_keys(header)
    )
    if layers in ((), ("gzip",), ("bgzf",)) and as_is:
        return path
    with cache_index.writing(cache, cache_index.REWRITE_COPY, {"REWRITE_VERSION": REWRITE_VERSION}, origin):
        return _write_copy(opener, header, missing, cache / "structural" if structural_only else cache, structural_only)


def _write_copy(opener, header: Header, missing: list[str], cache: Path, structural_only: bool) -> Path:
    """The copy `readable_path` describes, in the directory `cache`; reused when it is complete."""
    target = cache / "input.vcf.gz"
    done = cache / "input.vcf.gz.complete"
    if done.exists():
        return target
    cache.mkdir(parents=True, exist_ok=True)
    declared: set[tuple[str, str]] = set()
    refused = refused_types(header)
    with opener() as source, target.open("wb") as raw, BgzfWriter(raw) as out:
        out.write(f"##fileformat={header.version or 'VCFv4.2'}\n".encode())
        for line in header.meta:
            if line.startswith("##fileformat=") or (structural_only and not _keep_in_structural_copy(line, refused)):
                continue
            key = re.match(r"^##(\w+)=<ID=([^,>]+)", line)
            if structural_only and key:
                if key.groups() in declared:
                    continue
                declared.add(key.groups())
            out.write(f"{_readable_number(line)}\n".encode())
        for name in missing:
            namespace, key = name.split("/", 1)
            number = (INFO_KEYS if namespace == "INFO" else FORMAT_KEYS)[key][0]
            number = "." if number in READ_AS_ANY_NUMBER else number
            kind = declared_type(header, name)
            out.write(
                f'##{namespace}=<ID={key},Number={number},Type={kind},'
                f'Description="Declared by just-vcf: used by the rows, missing from the header">\n'.encode()
            )
        out.write(("#" + "\t".join(header.columns) + "\n").encode())
        if structural_only:
            # polars-bio refuses POS 0 ("Missing variant start"); Illumina array exports put unmapped
            # probes there (CHROM NA). They are left out of the copy and counted.
            left_out = 0
            for raw in source:
                if raw.startswith(b"#"):
                    continue
                fields = raw.split(b"\t", 2)
                if len(fields) > 1 and fields[1] == b"0":
                    left_out += 1
                    continue
                out.write(raw)
            (cache / "pos_zero_left_out").write_text(str(left_out))
        else:
            for raw in source:
                if not raw.startswith(b"#"):
                    out.write(raw)
                    break
            shutil.copyfileobj(source, out, 16 << 20)
    done.touch()
    return target


def readable_bcf_path(path: Path, layers: tuple[str, ...], opener, cache: Path, origin: cache_index.Origin) -> Path:
    """A path polars-bio reads as BCF: it picks BCF by the `.bcf` name and refuses any other name
    ("empty input", probed with 0.36). The file itself when it is named so and not wrapped in another
    layer, else a copy under `cache` named `input.bcf`: the bytes as they are when only BGZF wraps
    them, else the decompressed payload written as BGZF. Written under the entry's lock, as `readable_path`."""
    if layers in ((), ("bgzf",)) and path.name.lower().endswith(".bcf"):
        return path
    with cache_index.writing(cache, cache_index.REWRITE_COPY, {"REWRITE_VERSION": REWRITE_VERSION}, origin):
        return _write_bcf_copy(path, layers, opener, cache)


def _write_bcf_copy(path: Path, layers: tuple[str, ...], opener, cache: Path) -> Path:
    target = cache / "input.bcf"
    done = cache / "input.bcf.complete"
    if done.exists():
        return target
    cache.mkdir(parents=True, exist_ok=True)
    if layers in ((), ("bgzf",)):
        shutil.copyfile(path, target)
    else:
        with opener() as source, target.open("wb") as raw, BgzfWriter(raw) as out:
            shutil.copyfileobj(source, out, 16 << 20)
    done.touch()
    return target


def pos_zero_left_out(cache: Path) -> int:
    """Records at POS 0 left out of the structural copy under `cache`."""
    marker = cache / "structural" / "pos_zero_left_out"
    return int(marker.read_text()) if marker.exists() else 0


def scan(path: Path, info_fields: list[str], format_fields: list[str], sample: str | None) -> pl.LazyFrame:
    """A polars-bio scan of `path` that is safe to collect more than once, and concurrently.

    `pb.scan_vcf` registers a table named after the file stem each time a plan runs, so two scans of
    one file in one plan (or `pl.collect_all` over frames from the same file) panic with "table
    already exists" (P6). Each source is registered once, under a name derived from its path, size,
    modification time and the fields and sample read, and queried by that name. Ported from just-prs
    `vcf.scan_vcf` (docs/PROVENANCE.md).

    Only the carried FORMAT keys and the chosen sample are parsed, so a malformed field the frames do
    not carry cannot fail the load. The public `register_vcf` takes neither, so polars-bio's own
    `SQL._register_variant` is called (P7). A named registration returns 0-based `start`; it is
    converted to the 1-based `POS` here, and `read_records` checks it against the first record.
    """
    stat = path.stat()
    key = f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}:{sorted(info_fields)}:{sorted(format_fields)}:{sample}"
    name = f"just_vcf_{hashlib.sha256(key.encode()).hexdigest()[:24]}"
    with _REGISTRATION_LOCK:
        if name not in _REGISTERED:
            _register(
                path=str(path),
                name=name,
                info_fields=info_fields,
                format_fields=format_fields,
                samples=[sample] if sample is not None else None,
                genotype_output="string",
                chunk_size=8,
                concurrent_fetches=1,
                allow_anonymous=True,
                max_retries=5,
                timeout=300,
                enable_request_payer=False,
                compression_type="auto",
            )
            _REGISTERED.add(name)
    return pb.sql(f'SELECT * FROM "{name}"').with_columns((pl.col("start") + 1).alias("start"))


def _collect(plan: pl.LazyFrame) -> pl.DataFrame:
    """`plan` collected, with a field value polars-bio cannot parse raised as `HeaderRefused` (P7)."""
    try:
        return plan.collect(engine="streaming")
    except pl.exceptions.ComputeError as err:
        if refused := _refused(err):
            raise refused from err
        raise


def _collect_all(plans: list[pl.LazyFrame]) -> list[pl.DataFrame]:
    """`plans` collected together, as `_collect` collects one."""
    try:
        return pl.collect_all(plans, engine="streaming")
    except pl.exceptions.ComputeError as err:
        if refused := _refused(err):
            raise refused from err
        raise


def _refused(err: pl.exceptions.ComputeError) -> HeaderRefused | None:
    """polars-bio's failure to parse a field value, as `HeaderRefused`; None for any other failure."""
    if ("Error reading" in str(err) and "field" in str(err)) or "Missing variant start" in str(err):
        return HeaderRefused(str(err).splitlines()[-1][-160:])
    return None


def _register(**kwargs) -> None:
    """`SQL._register_variant`, with a refused header raised as `HeaderRefused`."""
    try:
        SQL._register_variant(**kwargs)
    except ValueError as err:
        if "Failed to register table" in str(err):
            raise HeaderRefused(str(err)) from err
        raise


def read_records(
    path: Path, header: Header, sample: str, first_record: bytes | None, extra_info: tuple[str, ...] = ()
) -> pl.LazyFrame:
    """The records of one sample with flat, namespaced columns.

    INFO keys become `info_<KEY>` and FORMAT keys `fmt_<KEY>`, whatever polars-bio names them. Every
    carried key is present, as a typed null when the file does not declare it. `extra_info` are INFO
    keys read for selecting rows (`SF`), as lists of strings.
    """
    refused = refused_types(header)
    info = [key for key in INFO_KEYS if key in header.info and f"INFO/{key}" not in refused]
    info += sorted(header.split_markers) + list(extra_info)
    unfit = unfit_keys(header)
    fmt = [key for key in FORMAT_KEYS if key in header.format and key not in unfit and f"FORMAT/{key}" not in refused]
    multi = len(header.samples) > 1
    frame = scan(path, info, fmt, sample if multi else None)
    if first_record is not None:
        fields = first_record.decode("utf-8", errors="replace").split("\t")
        got = _collect(frame.select("chrom", "start").head(1)).row(0)
        if (got[0], got[1]) != (fields[0], int(fields[1])):
            raise RuntimeError(
                f"polars-bio returned {got[0]}:{got[1]} for the first record, which the file writes as "
                f"{fields[0]}:{fields[1]}; the coordinate convention changed upstream"
            )
    names = frame.collect_schema().names()
    renames = {key: f"info_{key}" for key in info}
    if "genotypes" in names:
        present = [f.name for f in frame.collect_schema()["genotypes"].fields]
        frame = frame.with_columns(
            pl.col("genotypes").struct.field(key).list.get(0, null_on_oob=True).alias(f"fmt_{key}")
            for key in fmt
            if key in present
        ).drop("genotypes")
    else:
        for key in fmt:
            renames[f"fmt_{key}" if key in info else key] = f"fmt_{key}"
    frame = frame.rename({old: new for old, new in renames.items() if old in names})
    schema = frame.collect_schema()
    present = set(schema.names())
    reshaped = []
    for key, (number, _) in FORMAT_KEYS.items():
        name = f"fmt_{key}"
        if name not in present:
            continue
        is_list = isinstance(schema[name], pl.List)
        if number == "1" and is_list:
            # Declared with another Number (Platypus writes GQ as Number=.): keep the value, and count
            # the rows where there was more than one (`scan_stats`).
            reshaped += [pl.col(name).list.first().alias(name), (pl.col(name).list.len() > 1).alias(f"_multi_{key}")]
        elif number != "1" and not is_list:
            reshaped.append(pl.concat_list(pl.col(name)).alias(name))
    for key in _LIST_INFO:
        name = f"info_{key}"
        if name in present and not isinstance(schema[name], pl.List):
            # Declared Number=1, as older SV callers do: one value, read as the list of a one-ALT record.
            reshaped.append(pl.concat_list(pl.col(name)).alias(name))
    if "info_SVTYPE" in present and isinstance(schema["info_SVTYPE"], pl.List):
        reshaped.append(pl.col("info_SVTYPE").list.first().alias("info_SVTYPE"))
    frame = frame.with_columns(reshaped)
    # A marker declared with a Number other than 1 comes back as a list: its first value is enough.
    frame = frame.with_columns(
        (pl.col(name).list.first() if isinstance(schema[name], pl.List) else pl.col(name)).cast(pl.String).alias(name)
        for name in (f"info_{key}" for key in header.split_markers)
        if name in present
    )
    frame = frame.with_columns(
        (pl.col(name).cast(pl.List(pl.String)) if isinstance(schema[name], pl.List) else pl.col(name).cast(pl.String).str.split(",")).alias(name)
        for name in (f"info_{key}" for key in extra_info)
        if name in present
    )
    return frame.with_columns(
        [pl.lit(None, dtype=_dtype(number, kind)).alias(f"fmt_{key}") for key, (number, kind) in FORMAT_KEYS.items() if f"fmt_{key}" not in present]
        + [
            pl.lit(None, dtype=pl.List(_dtype("1", kind)) if key in _LIST_INFO else _dtype("1", kind)).alias(f"info_{key}")
            for key, (_, kind) in INFO_KEYS.items()
            if f"info_{key}" not in present
        ]
    ).with_columns(pl.col("fmt_GT").cast(pl.String), pl.col("fmt_FT").cast(pl.String))


_BOUNDED = {"GQ": 2**16 - 1, "DP": 2**32 - 1, "MIN_DP": 2**32 - 1}


def scan_stats(records: pl.LazyFrame, unsourced: pl.LazyFrame | None = None) -> dict:
    """One streaming pass for what the warnings and the facet need: out-of-range counts, GP's
    largest value, the source contig names, and whether any record is a reference block. It reads
    GT, which makes polars-bio parse the sample column, so a field value it cannot parse fails here.

    `unsourced` (`Part.unsourced`, a merge's first part only) is counted in the same collect, as
    `unsourced`: the merge's records no part takes, which `records`, already one part's, no longer
    holds. It is a second scan of the file, read beside the first."""
    out_of_range = [
        ((pl.col(f"fmt_{key}").cast(pl.Int64, strict=False) < 0) | (pl.col(f"fmt_{key}").cast(pl.Int64, strict=False) > top))
        .sum()
        .alias(f"out_{key}")
        for key, top in _BOUNDED.items()
    ]
    alts = _alts_expr()
    plan = records.select(
        *out_of_range,
        pl.col("fmt_GP").cast(pl.List(pl.Float64), strict=False).list.max().max().alias("gp_max"),
        *_end_stats(alts, pl.col("ref")),
        _is_block_record(alts).any().alias("any_block"),
        (alts.list.len().eq(1) & alts.list.first().eq(CGI_CNV_WINDOW)).sum().alias("cnv_windows"),
        alts.list.contains(DASH_DELETION).sum().alias("dash_alts"),
        *[
            (pl.col(f"info_{key}").is_not_null() & _ci_pairs(f"info_{key}", alts).is_null()).sum().alias(f"unfit_{key}")
            for key in ("CIPOS", "CIEND")
        ],
        (pl.col("info_CN").is_not_null() & _per_alt("info_CN", alts).is_null()).sum().alias("unfit_info_CN"),
        alts.list.eval(~pl.element().str.contains(_VCF_ALLELE_OR_DASH)).list.any().sum().alias("not_allele_alts"),
        pl.col("fmt_GP").cast(pl.List(pl.Float64), strict=False).is_not_null().sum().alias("gp_rows"),
        (pl.col("fmt_GP").cast(pl.List(pl.Float64), strict=False).list.sum() - 1.0).abs().gt(0.02).sum().alias("gp_off_one"),
        pl.col("chrom").unique().implode().alias("chroms"),
        # Reading GT makes polars-bio parse every sample field here, where a value it cannot parse
        # can still be retried (P7), rather than when the caller collects the frames.
        pl.col("fmt_GT").null_count().alias("gt_missing"),
        (_gt_body().str.count_matches(r"[/|]").max() + 1).alias("ploidy_max"),
        *_psl_stats(),
        *[pl.col(name).sum().alias(name) for name in records.collect_schema().names() if name.startswith("_multi_")],
        *[
            (pl.col(f"fmt_{key}").is_not_null() & ~_fits(pl.col(f"fmt_{key}"), key, alts)).sum().alias(f"unfit_{key}")
            for key in VECTOR_NUMBERS
        ],
        *(
            [pl.col("info_SF").list.eval(pl.element().str.strip_chars_end("f")).list.n_unique().gt(1).sum().alias("shared_sources")]
            if "info_SF" in records.collect_schema().names()
            else []
        ),
    )
    if unsourced is None:
        row = _collect(plan).row(0, named=True)
    else:
        found, left_out = _collect_all([plan, unsourced.select(pl.len().alias("unsourced"))])
        row = {**found.row(0, named=True), "unsourced": left_out.item()}
    row["chroms"] = sorted(row["chroms"] or [])
    return row


def _fits(vector: pl.Expr, key: str, alts: pl.Expr) -> pl.Expr:
    """True where `vector` has the length its Number requires on this row: one value per allele
    (R), or one per genotype of the row's ploidy (G; checked for haploid and diploid GT, and diploid
    when GT is missing, as the spec says to assume)."""
    alleles = alts.list.len().cast(pl.Int64) + 1
    if VECTOR_NUMBERS[key] == "R":
        return vector.list.len() == alleles
    ploidy = (_gt_body().str.count_matches(r"[/|]") + 1).fill_null(2)
    genotypes = (
        pl.when(ploidy == 1).then(alleles).when(ploidy == 2).then(alleles * (alleles + 1) // 2).otherwise(None)
    )
    return (vector.list.len() == genotypes).fill_null(True)


def _gt_body() -> pl.Expr:
    """GT text without a VCF 4.4 leading phase indicator."""
    return pl.col("fmt_GT").str.replace_all(r"^[/|]", "")


def _gt_expr() -> pl.Expr:
    """GT's allele indices, a null per `.`; null where the record has no GT."""
    return (
        pl.when(pl.col("fmt_GT").is_null())
        .then(None)
        .otherwise(
            _gt_body().str.replace_all(r"\|", "/").str.split("/").list.eval(
                pl.when(pl.element() == ".").then(None).otherwise(pl.element()).cast(pl.UInt8, strict=False)
            )
        )
    )


def _tagged(key: str, kind: str) -> pl.Expr:
    """True on a record INFO `key` tags as split off a multi-allelic site: any vt or just-vcf value, or
    a bcftools value with its fifth, `USED_ALT_IDX`, field."""
    column = pl.col(f"info_{key}")
    return column.str.count_matches(r"\|").fill_null(0) == 4 if kind == "bcftools" else column.is_not_null()


def _marker_expr(header: Header) -> pl.Expr:
    """True on a record any declared split marker tags."""
    tagged = [_tagged(key, kind) for key, kind in sorted(header.split_markers.items())]
    return pl.any_horizontal(tagged) if tagged else pl.lit(False)


def split_groups(calls: pl.LazyFrame, ploidy: int) -> tuple[pl.DataFrame, int]:
    """Positions whose rows look like one multi-allelic site written one row per ALT (RM270).

    `calls` is a facet's `calls` frame before this pass, with `split_signal` set on marker-tagged rows
    only; blocks are never in it. Returns the flagged positions, `(chrom, pos, rows, split_group)`,
    where `split_group` is the pattern a position fits, and the number of marker-tagged rows. One
    pass over three columns; the result holds only flagged positions.
    """
    slots = range(max(ploidy, 1))
    nonzero = [pl.col("gt").list.get(i, null_on_oob=True).fill_null(0) > 0 for i in slots]
    all_zero = pl.col("gt").list.eval(pl.element().is_not_null() & (pl.element() == 0)).list.all() & (pl.col("gt").list.len() > 0)
    per_row = calls.select(
        "chrom",
        "pos",
        pl.any_horizontal(nonzero).alias("nonref"),
        all_zero.fill_null(False).alias("all_zero"),
        *[flag.alias(f"slot_{i}") for i, flag in zip(slots, nonzero)],
    )
    groups = (
        per_row.group_by("chrom", "pos")
        .agg(
            pl.len().alias("rows"),
            pl.col("nonref").sum().alias("nonref"),
            pl.col("all_zero").sum().alias("all_zero"),
            *[pl.col(f"slot_{i}").sum().alias(f"slot_{i}") for i in slots],
        )
        .filter(pl.col("rows") > 1)
        .select(
            "chrom",
            "pos",
            "rows",
            pl.when((pl.col("nonref") >= 2) & (pl.max_horizontal(f"slot_{i}" for i in slots) <= 1))
            .then(pl.lit(SplitSignal.SLOT_PATTERN.value))
            .when((pl.col("all_zero") >= 1) & (pl.col("nonref") >= 1))
            .then(pl.lit(SplitSignal.ALL_ZERO_BESIDE_ALT.value))
            .alias("split_group"),
        )
        .filter(pl.col("split_group").is_not_null())
    )
    marked = calls.select((pl.col("split_signal") == SplitSignal.MARKER.value).sum())
    flagged, marker_rows = pl.collect_all([groups, marked])
    return flagged.sort("chrom", "pos"), marker_rows.item()


def with_split_groups(calls: pl.LazyFrame, flagged: pl.DataFrame) -> pl.LazyFrame:
    """`calls` with each flagged position's pattern in `split_signal`, where no marker set it first."""
    if flagged.is_empty():
        return calls
    return (
        calls.join(flagged.lazy().select("chrom", "pos", "split_group"), on=["chrom", "pos"], how="left", maintain_order="left")
        .with_columns(pl.coalesce(pl.col("split_signal"), pl.col("split_group").cast(CALLS["split_signal"])).alias("split_signal"))
        .drop("split_group")
    )


def split_warnings(flagged: pl.DataFrame, marker_rows: int, header: Header) -> list[str]:
    """One warning per RM270 signal that fired, with its count."""
    warnings = []
    if marker_rows:
        warnings.append(
            f"{marker_rows} row(s) tagged as split off a multi-allelic site (INFO {', '.join(sorted(header.split_markers))}): "
            "a `0` there may mean another ALT, not REF (RM270)"
        )
    for signal in (SplitSignal.SLOT_PATTERN, SplitSignal.ALL_ZERO_BESIDE_ALT):
        group = flagged.filter(pl.col("split_group") == signal.value)
        if group.height:
            warnings.append(
                f"{group['rows'].sum()} row(s) at {group.height} position(s) look split one row per ALT ({signal.value}): "
                "a `0` there may mean another ALT, not REF (RM270)"
            )
    return warnings


def _alts_expr() -> pl.Expr:
    return (
        pl.when(pl.col("alt") == "")
        .then(pl.lit([], dtype=pl.List(pl.String)))
        .otherwise(pl.col("alt").str.split("|"))
        .list.eval(pl.when(pl.element() == "<NON_REF>").then(pl.lit("<*>")).otherwise(pl.element()))
    )


def to_rows(
    records: pl.LazyFrame, header: Header, stats: dict, text_gt: bool = True, gl_phred: bool = False
) -> tuple[pl.LazyFrame, pl.LazyFrame, list[str]]:
    """Lazy `calls` and `blocks` rows (without `facet_id`) from `read_records` output, plus warnings.

    `text_gt` is False for a BCF: its genotypes are binary, so the GT text polars-bio renders is
    not something the file wrote, and `source_gt` stays null. `gl_phred` is True for Complete
    Genomics' CGAPipeline, which writes GL as integers in phred units (GQ is the gap between the two
    best values), so `gl` is GL / 10.

    Complete Genomics' record classes are read by their declared alleles, whatever the format id:
    `<CGA_CNVWIN>` rows are copy-number windows and are left out; `<CGA_NOCALL>` rows are no-called
    spans to END; and where the header declares `CGA_NOCALL`, a record with ALT `.`, no END and a GT
    with no allele called is a no-called span over its REF, so it goes to `blocks` too.
    """
    warnings = [
        f"{stats[f'out_{key}']} {key} value(s) outside 0..{top} set to null"
        for key, top in _BOUNDED.items()
        if stats[f"out_{key}"]
    ] + [
        f"{stats[f'unfit_{key}']} {key} value(s) not one per {'allele' if number == 'R' else 'genotype'} set to null"
        for key, number in VECTOR_NUMBERS.items()
        if stats.get(f"unfit_{key}")
    ] + [
        f"{count} {name.removeprefix('_multi_')} value(s) held more than one number; the first was kept"
        for name, count in stats.items()
        if name.startswith("_multi_") and count
    ]
    if stats.get("sv_end_disagrees"):
        warnings.append(f"{stats['sv_end_disagrees']} symbolic SV record(s) whose INFO END disagrees with POS + |SVLEN|: end_pos follows END")
    if stats.get("end_not_read"):
        warnings.append(
            f"{stats['end_not_read']} record(s) carry an INFO END that is not their extent (an <INS>, a breakend, "
            "a sequence allele or a symbolic allele of no known span): end_pos does not read it"
        )
    for key in ("CIPOS", "CIEND"):
        if stats.get(f"unfit_{key}"):
            warnings.append(
                f"{stats[f'unfit_{key}']} {key} value(s) not two per ALT (VCF 4.5) nor one pair for the record "
                "(VCF 4.2's Number=2) set to null"
            )
    if stats.get("unfit_info_CN"):
        warnings.append(f"{stats['unfit_info_CN']} INFO CN value(s) not one per ALT set to null")
    if stats.get("dash_alts"):
        warnings.append(
            f"{stats['dash_alts']} record(s) write a deletion as ALT `-`, which VCF has no allele for (F12): read as "
            'the empty allele "", a deletion of the REF bases with no anchor base; normalize anchors it on the base '
            "before, and write_vcf refuses it until then"
        )
    if stats.get("not_allele_alts"):
        warnings.append(
            f"{stats['not_allele_alts']} record(s) carry an ALT that is no VCF allele (neither bases ACGTN, "
            "a symbolic <…> allele, a breakend nor *), F15: kept as written, with no variant class"
        )
    v44 = (header.version_tuple or (0, 0)) >= (4, 4)
    # The file's REF goes to `source_ref` only. `ref` is read from the assembly sequence or stays null
    # (charter 10), so polars-bio's `ref` column must not reach the frame under that name.
    # F12's `-` becomes the empty allele here, after the stats counted it as written, and after `.` became [].
    alts = _alts_expr().list.eval(pl.when(pl.element() == DASH_DELETION).then(pl.lit(EMPTY_ALLELE)).otherwise(pl.element()))
    rows = records.with_columns(alts.alias("alts"), pl.col("ref").alias("source_ref")).drop("ref")
    if stats.get("cnv_windows"):
        warnings.append(f"{stats['cnv_windows']} {CGI_CNV_WINDOW} record(s) left out: copy-number windows, not genotypes")
        rows = rows.filter(~(pl.col("alts").list.len().eq(1) & pl.col("alts").list.first().eq(CGI_CNV_WINDOW)))
    rows = rows.with_columns(
        _end_pos(pl.col("alts"), pl.col("source_ref")).alias("end_pos"),
        _is_block_record(pl.col("alts")).alias("_block"),
    )

    gt_text = pl.col("fmt_GT")
    seps = _gt_body().str.extract_all(r"[/|]")
    first_phased = pl.lit(None, dtype=pl.Boolean) if v44 else pl.when(seps.list.contains("/")).then(pl.lit(False)).otherwise(pl.lit(True))
    rows = rows.with_columns(
        _gt_expr().alias("gt"),
        pl.when(gt_text.is_null()).then(None).otherwise(pl.concat_list(first_phased, seps.list.eval(pl.element() == "|"))).alias("phased"),
    )
    allele_list = pl.concat_list(pl.col("source_ref"), pl.col("alts"))
    rows = rows.with_columns(allele_list.list.gather(pl.col("gt"), null_on_oob=True).alias("gt_alleles"))
    # A record is a block only when it asserts reference or nothing: a GT calling `<*>` stays a call.
    reference_or_nothing = pl.col("gt").is_null() | pl.col("gt").list.eval(pl.element().is_null() | (pl.element() == 0)).list.all()
    rows = rows.with_columns((pl.col("_block") & reference_or_nothing).alias("_block"))
    if CGI_NOCALL.strip("<>") in header.alt_ids:
        # polars-bio returns a haploid no-call `.` as a null GT (male X and Y), so a null GT on such a
        # record is a no-call too: CGI writes GT on every record.
        no_call_span = (
            pl.col("alts").list.len().eq(0)
            & pl.col("info_END").is_null()
            & (pl.col("gt").is_null() | pl.col("gt").list.eval(pl.element().is_null()).list.all())
        )
        rows = rows.with_columns((pl.col("_block") | no_call_span.fill_null(False)).alias("_block"))
    # RM270: a marker tags its own row here; the patterns need whole positions, read from the spill.
    rows = rows.with_columns(
        pl.when(_marker_expr(header)).then(pl.lit(SplitSignal.MARKER.value)).alias("split_signal")
    )
    rows = _phase_sets(rows, psl="PSL" in header.format, ploidy=stats.get("ploidy_max") or 1, v44=v44)
    warnings += psl_warnings(stats)

    gp, gp_warning = _gp(header, stats)
    if gp_warning:
        warnings.append(gp_warning)
    # A vector of the wrong length for its row cannot be indexed over the alleles: withheld (counted
    # in `scan_stats`), never truncated or padded.
    rows = rows.with_columns(
        pl.when(_fits(pl.col(f"fmt_{key}"), key, pl.col("alts"))).then(pl.col(f"fmt_{key}")).otherwise(None).alias(f"fmt_{key}")
        for key in VECTOR_NUMBERS
    )
    rows = rows.with_columns(
        # Divided in Float64 and rounded to Float32 once, so PL 21 gives the Float32 nearest -2.1, the
        # value a file writing GL -2.1 loads as (F11); a Float32 division lands one step off.
        pl.coalesce(
            (
                pl.col("fmt_GL").cast(pl.List(pl.Float64), strict=False).list.eval(pl.element() / 10.0).cast(pl.List(pl.Float32))
                if gl_phred
                else pl.col("fmt_GL").cast(pl.List(pl.Float32), strict=False)
            ),
            pl.col("fmt_PL").cast(pl.List(pl.Float64), strict=False).list.eval(pl.element() / -10.0).cast(pl.List(pl.Float32)),
        ).alias("gl"),
        gp.alias("gp"),
        pl.when(allele_list.list.len() == 2)
        .then(pl.col("fmt_DS").cast(pl.List(pl.Float32), strict=False).list.first())
        .otherwise(None)
        .alias("dosage"),
        _bounded("GQ", pl.UInt16).alias("gq"),
        _bounded("DP", pl.UInt32).alias("dp"),
        _bounded("MIN_DP", pl.UInt32).alias("dp_min"),
        pl.col("fmt_AD").cast(pl.List(pl.UInt32), strict=False).alias("ad"),
        pl.when(pl.col("id") == "").then(None).otherwise(pl.col("id").str.split(";")).alias("ids"),
        pl.when(pl.col("filter") == "").then(None).otherwise(pl.col("filter").str.split(";")).alias("filters"),
        pl.when(pl.col("fmt_FT").is_in(["", "."])).then(None).otherwise(pl.col("fmt_FT").str.split(";")).alias("sample_filters"),
        pl.col("qual").cast(pl.Float32).alias("qual"),
        (pl.lit(None, dtype=pl.String) if v44 or not text_gt else gt_text).alias("source_gt"),
        canonical_chrom_expr("chrom").alias("chrom"),
        pl.col("start").alias("pos"),
        pl.when(pl.col("info_SVTYPE").is_in(["", "."])).then(None).otherwise(pl.col("info_SVTYPE")).alias("sv_type"),
        _ci_pairs("info_CIPOS", pl.col("alts")).alias("cipos"),
        _ci_pairs("info_CIEND", pl.col("alts")).alias("ciend"),
        pl.col("fmt_CN").cast(pl.Float32, strict=False).alias("copy_number"),
        _per_alt("info_CN", pl.col("alts")).alias("alt_copy_numbers"),
    )
    full_vector = pl.any_horizontal(pl.col("gl").is_not_null(), pl.col("gp").is_not_null(), pl.col("ad").is_not_null())
    rows = rows.with_columns(
        pl.when(full_vector).then(allele_list).otherwise(None).alias("source_alleles"),
        pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
    )
    present = set(rows.collect_schema().names())
    calls = rows.filter(~pl.col("_block")).select([name for name in CALLS.names() if name != "facet_id" and name in present])
    blocks = rows.filter(pl.col("_block")).select(
        "chrom",
        "pos",
        "end_pos",
        pl.when(pl.col("gt").list.eval(pl.element().is_not_null()).list.any())
        .then(pl.lit(BlockKind.REF.value))
        .otherwise(pl.lit(BlockKind.NO_CALL.value))
        .alias("kind"),
        "gq",
        "dp_min",
        "gl",
    )
    return calls, blocks, warnings


def _is_block_record(alts: pl.Expr) -> pl.Expr:
    """A reference-block record, before its GT is looked at: every ALT is `<*>` (GATK, VCF 4.5), or
    there is no ALT and INFO/END gives an extent (Illumina's gVCF: `ALT .` with `END`), or it is
    Complete Genomics' `<CGA_NOCALL>` span with END.

    An `ALT .` record to END that carries a copy number (FORMAT or INFO CN) is a copy-number segment,
    not a block: Canvas writes its copy-neutral segments so (`GT 0/0`, `CN 2`). It states the copy
    number of the region and nothing about small variants in it, so it is a `calls` row with
    `alts = []`, and the facet's absence does not become `blocks`."""
    all_star = alts.list.len().gt(0) & alts.list.eval(pl.element() == "<*>").list.all()
    copy_number = pl.col("fmt_CN").is_not_null() | pl.col("info_CN").is_not_null()
    illumina = alts.list.len().eq(0) & pl.col("info_END").is_not_null() & ~copy_number
    cgi = alts.list.len().eq(1) & alts.list.first().eq(CGI_NOCALL) & pl.col("info_END").is_not_null()
    return all_star | illumina | cgi


def _dtype(number: str, kind: str) -> pl.DataType:
    scalar = {"Integer": pl.Int32, "Float": pl.Float32, "String": pl.String}[kind]
    return scalar if number == "1" else pl.List(scalar)


def _ci_pairs(name: str, alts: pl.Expr) -> pl.Expr:
    """A confidence interval (INFO CIPOS or CIEND) as two values per ALT, VCF 4.5's shape: as written
    when the file wrote one pair per ALT; VCF 4.2's one pair for the record (`Number=2`) repeated for
    each ALT, since every ALT of a record shares its POS; null for any other count (counted in
    `scan_stats`). A record with no ALT (a copy-neutral segment) keeps its one pair."""
    values = pl.col(name).cast(pl.List(pl.Int32), strict=False)
    width = pl.max_horizontal(alts.list.len().cast(pl.Int64), pl.lit(1, dtype=pl.Int64))
    repeated = values.list.gather(pl.int_ranges(0, 2 * width).list.eval(pl.element() % 2))
    return pl.when(values.list.len() == 2 * width).then(values).when(values.list.len() == 2).then(repeated).otherwise(None)


def _per_alt(name: str, alts: pl.Expr) -> pl.Expr:
    """An INFO list with one value per ALT (`Number=A`) as Float32, null where its count is not that."""
    values = pl.col(name).cast(pl.List(pl.Float32), strict=False)
    return pl.when(values.list.len() == alts.list.len()).then(values).otherwise(None)


def _has_sv(alts: pl.Expr) -> pl.Expr:
    """True where an ALT is a symbolic SV that spans reference bases (`_SV_LENGTH_TYPES`)."""
    return alts.list.eval(pl.element().str.contains(_SV_PATTERN)).list.any()


def _svlen_end(alts: pl.Expr) -> pl.Expr:
    """`POS + |SVLEN|`, the largest over the length-type SV alleles that have an SVLEN; null where none
    has one. An SVLEN that is not one value per ALT cannot be matched to its alleles and is not read."""
    svlen = pl.col("info_SVLEN").cast(pl.List(pl.Int64), strict=False)
    fitting = pl.when(svlen.list.len() == alts.list.len()).then(svlen).otherwise(None)
    is_sv = alts.list.eval(pl.when(pl.element().str.contains(_SV_PATTERN)).then(pl.lit(1, dtype=pl.Int64)).otherwise(None))
    return pl.col("start").cast(pl.Int64) + (is_sv * fitting.list.eval(pl.element().abs())).list.max()


def _end_pos(alts: pl.Expr, ref: pl.Expr) -> pl.Expr:
    """VCF 4.5's computed END, per allele, then the largest over the site (SCHEMA.md `end_pos`).

    A length-type SV allele ends at the record's INFO END where it gives one, else at POS + |SVLEN|.
    END is read only there, and for `<*>` blocks (after FORMAT LEN), `<CGA_NOCALL>` spans and ALT `.`
    blocks: an `<INS>` or a breakend ends at its anchor, and another symbolic allele's END is not known
    to be on this contig (Delly's `<TRA>` writes the partner's position on CHR2 there).
    """
    ref_end = pl.col("start") + ref.str.len_chars().cast(pl.UInt32) - 1
    end = pl.col("info_END").cast(pl.Int64)
    block_end = (
        pl.when(pl.col("fmt_LEN").is_not_null())
        .then(pl.col("start") + pl.col("fmt_LEN").cast(pl.UInt32) - 1)
        .otherwise(end)
    )
    return pl.max_horizontal(
        ref_end,
        pl.when(alts.list.contains("<*>")).then(block_end).otherwise(None),
        pl.when(alts.list.contains(CGI_NOCALL)).then(end).otherwise(None),
        pl.when(alts.list.len() == 0).then(end).otherwise(None),
        pl.when(_has_sv(alts)).then(pl.coalesce(end, _svlen_end(alts))).otherwise(None),
    ).cast(pl.UInt32)


def _end_stats(alts: pl.Expr, ref: pl.Expr) -> list[pl.Expr]:
    """Counts for the warnings: SV records whose END and SVLEN disagree, and records whose END is
    not read because it is not their extent."""
    end = pl.col("info_END").cast(pl.Int64)
    ref_end = pl.col("start").cast(pl.Int64) + ref.str.len_chars().cast(pl.Int64) - 1
    has_sv = _has_sv(alts)
    read = has_sv | alts.list.contains("<*>") | alts.list.contains(CGI_NOCALL) | (alts.list.len() == 0)
    return [
        (has_sv & (end != _svlen_end(alts))).fill_null(False).sum().alias("sv_end_disagrees"),
        (~read & (end != ref_end)).fill_null(False).sum().alias("end_not_read"),
    ]


def _psl_names() -> pl.Expr:
    """FORMAT PSL's names, a null per `.` (and per empty name)."""
    return pl.col("fmt_PSL").list.eval(pl.when(pl.element().is_in(["", "."])).then(None).otherwise(pl.element()))


def _psl_present() -> pl.Expr:
    """True on a record whose PSL names at least one phase set."""
    return _psl_names().list.eval(pl.element().is_not_null()).list.any().fill_null(False)


def _psl_fits() -> pl.Expr:
    """True where PSL has one value per allele in GT (`Number=P`); a record without GT has no slots."""
    slots = _gt_body().str.count_matches(r"[/|]") + 1
    return (pl.col("fmt_PSL").list.len().cast(pl.Int64) == slots.cast(pl.Int64)).fill_null(False)


def _psl_stats() -> list[pl.Expr]:
    """Counts for the warnings about PSL, PSO and PSQ (`scan_stats`)."""
    present = _psl_present()
    fits = _psl_fits()
    # Slots after the first: `/` before the allele is unphased; a PSL name there breaks VCF 4.5's rule.
    # (Null where the lengths differ: list arithmetic is evaluated on every row, guarded or not.)
    unphased = pl.when(fits).then(_gt_body().str.extract_all(r"[/|]").list.eval((pl.element() == "/").cast(pl.Int64)))
    named = pl.when(fits).then(_psl_names().list.slice(1).list.eval(pl.element().is_not_null().cast(pl.Int64)))
    psq = pl.col("fmt_PSQ").cast(pl.List(pl.Int64), strict=False)
    pso = pl.col("fmt_PSO").cast(pl.List(pl.Int64), strict=False)
    return [
        present.sum().alias("psl_rows"),
        (present & ~fits).sum().alias("psl_unfit"),
        (present & fits & pl.col("fmt_PS").is_not_null()).sum().alias("psl_and_ps"),
        (present & fits & ((unphased * named).list.sum() > 0).fill_null(False)).sum().alias("psl_on_unphased"),
        (present & fits & (psq.list.len() != pl.col("fmt_PSL").list.len()).fill_null(False)).sum().alias("psq_unfit"),
        (present & fits & (pso.list.len() != pl.col("fmt_PSL").list.len()).fill_null(False)).sum().alias("pso_unfit"),
        (present & fits).cast(pl.Int64).mul(
            psq.list.eval(((pl.element() < 0) | (pl.element() > _BOUNDED["GQ"])).cast(pl.Int64)).list.sum().fill_null(0)
            + pso.list.eval(((pl.element() < 0) | (pl.element() > 2**32 - 1)).cast(pl.Int64)).list.sum().fill_null(0)
        ).sum().alias("psl_values_out"),
    ]


def psl_warnings(stats: dict) -> list[str]:
    """One warning per PSL rule a record broke, with its count."""
    warnings = []
    if stats.get("psl_unfit"):
        warnings.append(
            f"{stats['psl_unfit']} record(s) whose PSL is not one value per allele in GT (Number=P): their phase is "
            "withheld (phased null), since the set each phased allele is in cannot be read"
        )
    if stats.get("psl_and_ps"):
        warnings.append(
            f"{stats['psl_and_ps']} record(s) carry both PS and PSL, which VCF 4.5 forbids: PSL read, PS ignored"
        )
    if stats.get("psl_on_unphased"):
        warnings.append(
            f"{stats['psl_on_unphased']} record(s) name a PSL phase set on an unphased allele (VCF 4.5 says `.`): "
            "the name is not read, the allele stays unphased"
        )
    for key in ("PSO", "PSQ"):
        if stats.get(f"{key.lower()}_unfit"):
            warnings.append(f"{stats[f'{key.lower()}_unfit']} {key} value(s) not one per PSL entry set to null")
    if stats.get("psl_values_out"):
        warnings.append(f"{stats['psl_values_out']} PSO or PSQ value(s) out of range (negative, or too large) set to null")
    return warnings


def _phase_sets(rows: pl.LazyFrame, psl: bool = False, ploidy: int = 2, v44: bool = False) -> pl.LazyFrame:
    """`phase_sets` and `phase_set_kinds`, and from PSL `phase_set_ordinals` and `phase_set_qualities`.

    PS: the PS value on each phased haplotype, kind `ps`. Written without a join, as list arithmetic:
    a phased haplotype carries `PS + 1`, an unphased one 0, which then becomes null; a null PS (the
    chromosome's implicit set) stays null.

    PSL (`psl`, the file declares it; VCF 4.4+): on a record whose PSL names a set, each slot is read
    from its own entry (docs/SCHEMA.md, Decision 2). A phased slot with a name is in that set, kind
    `psl`, with the PSO and PSQ entries beside it; a phased slot with `.` is in the chromosome's
    implicit set (kind `ps`, set null: it has no PS subfield either); an unphased slot is in none,
    and a name there is not read (counted). In a VCF 4.4+ file the first slot's phase is unknown
    (P4), and a name there is what says it is phased: VCF 4.5 gives an unphased allele `.`. A PSL
    that is not one value per GT allele withholds the record's phase (counted); PS is ignored where
    PSL names a set (counted). `ploidy` is the file's largest, the slots written out one by one.
    """
    encoded = pl.col("phased").list.eval(pl.element().cast(pl.Int64)) * (pl.col("fmt_PS").cast(pl.Int64) + 1)
    phased_any = pl.col("phased").list.eval(pl.element().fill_null(False)).list.any()
    ps_sets = (
        pl.when(phased_any)
        .then(encoded.list.eval(pl.when(pl.element() > 0).then((pl.element() - 1).cast(pl.String)).otherwise(None)))
        .otherwise(None)
    )
    ps_kinds = (
        pl.when(phased_any)
        .then(pl.col("phased").list.eval(pl.when(pl.element()).then(pl.lit(PhaseSetKind.PS.value)).otherwise(None)))
        .otherwise(None)
    )
    if not psl:
        return rows.with_columns(
            ps_sets.alias("phase_sets"),
            ps_kinds.alias("phase_set_kinds"),
            pl.lit(None, dtype=pl.List(pl.UInt32)).alias("phase_set_ordinals"),
            pl.lit(None, dtype=pl.List(pl.UInt16)).alias("phase_set_qualities"),
        )
    names = _psl_names()
    width = pl.col("phased").list.len()
    pso = pl.col("fmt_PSO").cast(pl.List(pl.Int64), strict=False)
    psq = pl.col("fmt_PSQ").cast(pl.List(pl.Int64), strict=False)
    pso = pl.when(pso.list.len() == pl.col("fmt_PSL").list.len()).then(pso)
    psq = pl.when(psq.list.len() == pl.col("fmt_PSL").list.len()).then(psq)
    phased, sets, kinds, ordinals, qualities = [], [], [], [], []
    for slot in range(max(ploidy, 1)):
        was = pl.col("phased").list.get(slot, null_on_oob=True)
        name = names.list.get(slot, null_on_oob=True)
        now = pl.when(was.is_null() & name.is_not_null() & pl.lit(slot == 0 and v44)).then(pl.lit(True)).otherwise(was)
        named = (now & name.is_not_null()).fill_null(False)
        phased.append(now)
        sets.append(pl.when(named).then(name))
        kinds.append(pl.when(named).then(pl.lit(PhaseSetKind.PSL.value)).when(now.fill_null(False)).then(pl.lit(PhaseSetKind.PS.value)))
        ordinals.append(pl.when(named).then(pso.list.get(slot, null_on_oob=True)))
        qualities.append(pl.when(named).then(psq.list.get(slot, null_on_oob=True)))

    def per_slot(values: list[pl.Expr]) -> pl.Expr:
        return pl.concat_list(values).list.head(width)

    def bounded(values: list[pl.Expr], top: int, dtype: pl.DataType) -> pl.Expr:
        listed = per_slot(values).list.eval(pl.when((pl.element() >= 0) & (pl.element() <= top)).then(pl.element()).cast(dtype))
        return pl.when(listed.list.eval(pl.element().is_not_null()).list.any()).then(listed)

    # As for PS: a record with no phased slot has no phase-set lists.
    read = _psl_present() & _psl_fits()
    in_some_set = read & per_slot(kinds).list.eval(pl.element().is_not_null()).list.any().fill_null(False)
    unfit = _psl_present() & ~_psl_fits()
    return rows.with_columns(
        pl.when(read)
        .then(per_slot(phased))
        .when(unfit)
        .then(pl.col("phased").list.eval(pl.when(pl.element().is_null()).then(None).otherwise(None).cast(pl.Boolean)))
        .otherwise(pl.col("phased"))
        .alias("_phased"),
        pl.when(in_some_set).then(per_slot(sets)).when(read | unfit).then(None).otherwise(ps_sets).alias("phase_sets"),
        pl.when(in_some_set).then(per_slot(kinds)).when(read | unfit).then(None).otherwise(ps_kinds).alias("phase_set_kinds"),
        pl.when(read).then(bounded(ordinals, 2**32 - 1, pl.UInt32)).alias("phase_set_ordinals"),
        pl.when(read).then(bounded(qualities, _BOUNDED["GQ"], pl.UInt16)).alias("phase_set_qualities"),
    ).with_columns(pl.col("_phased").alias("phased")).drop("_phased")


def _gp(header: Header, stats: dict) -> tuple[pl.Expr, str | None]:
    """`gp` from GP or PP. GP's scale is read from the values where they settle it, else from the
    header's VCF version where it defines one, else withheld.

    VCF before 4.3 wrote GP phred-scaled, and 4.3 made it a probability. Real files do not follow the
    version: imputation tools wrote probabilities under 4.1 and 4.2 headers. So values above 1 mean
    phred, and values that are all within [0, 1] and sum to 1 per row mean probabilities, whatever
    the header says; the warnings name a header that disagreed. Values within [0, 1] that do not sum
    to 1 (posteriors converted from rounded phred, as just-vcf's own 4.3 export writes them) are
    probabilities under a 4.3+ header, which defines GP so, with the rows counted in a warning. Under
    an older header, or one that states no version, they settle nothing, and GP is not used.
    """
    from_pp = pl.col("fmt_PP").cast(pl.List(pl.Float32), strict=False).list.eval(pl.lit(10.0) ** (pl.element() / -10.0))
    largest = stats["gp_max"]
    if largest is None:
        return from_pp, None
    raw = pl.col("fmt_GP").cast(pl.List(pl.Float32), strict=False)
    version = header.version_tuple
    defined = None if version is None else ("phred" if version < (4, 3) else "probabilities")
    off_one = f"{stats['gp_off_one']} of {stats['gp_rows']} rows do not sum to 1"
    if largest > 1.0 + 1e-6:
        found, gp = "phred", raw.list.eval(pl.lit(10.0) ** (pl.element() / -10.0))
    elif stats["gp_off_one"] == 0:
        found, gp = "probabilities", raw
    elif defined == "probabilities":
        return pl.coalesce(raw, from_pp), f"GP read as probabilities, as the {header.version} header defines it: {off_one}"
    else:
        stated = "states no VCF version" if version is None else f"is {header.version}, where GP was phred-scaled"
        return from_pp, f"GP is ambiguous: every value is within [0, 1] but {off_one}, and the header {stated}; GP not used"
    note = None if defined in (None, found) else f"GP written as {found} under a {header.version} header; read as {found}"
    return pl.coalesce(gp, from_pp), note


def _bounded(key: str, dtype: pl.DataType) -> pl.Expr:
    """FORMAT `key` cast to `dtype`; values outside it become null (counted by `scan_stats`)."""
    value = pl.col(f"fmt_{key}").cast(pl.Int64, strict=False)
    return pl.when((value >= 0) & (value <= _BOUNDED[key])).then(value).otherwise(None).cast(dtype)


def dropped_fields(header: Header, data: list[bytes] = ()) -> list[str]:
    """Keys present in the file that the frames do not carry, namespaced: the declared ones, and the
    ones the rows in `data` use without a declaration (polars-bio drops those unread, P1)."""
    info, fmt = row_keys(list(data))
    info |= set(header.info)
    fmt |= set(header.format)
    refused = refused_types(header)
    dropped_info = [
        f"INFO/{key}" for key in info if (key not in INFO_KEYS and key not in header.split_markers) or f"INFO/{key}" in refused
    ]
    dropped_fmt = [f"FORMAT/{key}" for key in fmt if key not in FORMAT_KEYS or key in unfit_keys(header) or f"FORMAT/{key}" in refused]
    return sorted(dropped_info) + sorted(dropped_fmt)


def to_frames(calls: pl.LazyFrame, blocks: pl.LazyFrame, facet_id: str) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    return (
        conform(calls.with_columns(pl.lit(facet_id).alias("facet_id")), CALLS),
        conform(blocks.with_columns(pl.lit(facet_id).alias("facet_id")), BLOCKS),
    )
