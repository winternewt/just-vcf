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
haplotype's `phased` entry is withheld rather than reported from the rewritten text.
"""

import hashlib
import logging
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import BinaryIO

import polars as pl
import polars_bio as pb
from polars_bio.sql import SQL

from just_vcf.bgzf import BgzfWriter
from just_vcf.schema import BLOCKS, CALLS, BlockKind, EvidenceBasis, PhaseSetKind, RefCheck, SplitSignal, canonical_chrom_expr, conform

logger = logging.getLogger(__name__)

FORMATS: frozenset[str] = frozenset({"vcf"})
"""Format ids this loader reads. The dialects with their own ids join as their quirks are handled."""

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
    "LEN": ("1", "Integer"),
    "MIN_DP": ("1", "Integer"),
}
"""FORMAT keys the schema carries, with the Number and Type VCF 4.5 reserves for them."""

INFO_KEYS = ("END", "SVLEN")
"""INFO keys the loader reads, for `end_pos` only."""

SPLIT_MARKER_ID = "OLD_MULTIALLELIC"
"""The INFO key vt decompose and Glow write on a record split off a multi-allelic site."""

BCFTOOLS_OLD_REC = "Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX"
"""The Description bcftools norm writes for its `--old-rec-tag` key, whatever the key is named
(vcfnorm.c). Only a split record gets the fifth, `USED_ALT_IDX`, field; a realigned one has four."""

_SV_LENGTH_TYPES = ("DEL", "DUP", "INV", "CNV")
"""Symbolic allele types whose extent is `SVLEN` bases after POS. INS and breakends end at POS."""

_REGISTRATION_LOCK = Lock()
_REGISTERED: set[str] = set()

REWRITE_VERSION = 5
"""Bumped whenever `readable_path` changes what it writes, so a cached copy from older code is not reused."""

_HEADER_READ_CAP = 64 << 20
_DATA_LINES_CHECKED = 1000


class MultiSampleError(ValueError):
    """The file holds more than one sample and none was chosen."""


class HeaderRefused(ValueError):
    """polars-bio refused the file: its header (P8), or a field value it could not parse (P7)."""


_STRUCTURAL_META = ("##fileformat=", "##FILTER=", "##ALT=", "##contig=")
"""Meta lines kept in a structural-only copy, besides the INFO and FORMAT keys the frames carry.
Everything else is read by `read_header`, not by polars-bio."""


def _keep_in_structural_copy(line: str) -> bool:
    if line.startswith(_STRUCTURAL_META):
        return True
    declared = re.match(r"^##(INFO|FORMAT)=<(.*)>$", line)
    if not declared:
        return False
    kind, body = declared.groups()
    attrs = _attrs(body)
    if kind == "INFO":
        return attrs.get("ID") in INFO_KEYS or _split_marker_kind(attrs) is not None
    return attrs.get("ID") in FORMAT_KEYS


def _attrs(body: str) -> dict[str, str]:
    """The key=value pairs of a structured meta line, quotes left on."""
    return dict(re.findall(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,]*)', body))


def _split_marker_kind(attrs: dict[str, str]) -> str | None:
    """`vt` or `bcftools` when a declared INFO key marks records split off a multi-allelic site."""
    if attrs.get("ID") == SPLIT_MARKER_ID:
        return "vt"
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
    contigs: list[dict] = field(default_factory=list)
    sample_meta: list[str] = field(default_factory=list)
    format_first: bool = True
    split_markers: dict[str, str] = field(default_factory=dict)
    """Declared INFO keys that mark split records → `vt` or `bcftools` (`_split_marker_kind`)."""

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
    for raw in handle:
        read += len(raw)
        line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
        if line.startswith("##"):
            meta.append(line)
        elif line.startswith("#"):
            columns = line[1:].split("\t")
        elif line:
            data.append(raw)
            if len(data) >= _DATA_LINES_CHECKED:
                break
        if read > _HEADER_READ_CAP and not columns:
            raise ValueError("no #CHROM line in the first 64 MiB")
    if not columns:
        raise ValueError("no #CHROM line")
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
        if kind == "INFO" and "ID" in attrs:
            header.info[attrs["ID"]] = attrs.get("Number", ".")
            if marker := _split_marker_kind(attrs):
                header.split_markers[attrs["ID"]] = marker
        elif kind == "FORMAT" and "ID" in attrs:
            header.format[attrs["ID"]] = attrs.get("Number", ".")
        elif kind == "contig" and "ID" in attrs:
            length = attrs.get("length")
            header.contigs.append(
                {"name": attrs["ID"], "length": int(length) if length and length.isdigit() else None, "md5": attrs.get("md5")}
            )
        elif kind in ("SAMPLE", "PEDIGREE"):
            header.sample_meta.append(line)
    return header, data


def undeclared_format_keys(header: Header, data: list[bytes]) -> list[str]:
    """Carried FORMAT keys the rows use without the header declaring them."""
    used: set[str] = set()
    for raw in data:
        fields = raw.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")
        if len(fields) > 8:
            used.update(fields[8].split(":"))
    return sorted(key for key in used if key in FORMAT_KEYS and key not in header.format)


def readable_path(
    path: Path, layers: tuple[str, ...], opener, header: Header, missing: list[str], cache: Path, structural_only: bool = False
) -> Path:
    """A path polars-bio can read: the file itself when possible, else a rewritten copy under `cache`.

    The copy is kept, because the frames `load` returns are lazy and read it when collected. It is
    named after the source's digest and reused when complete, and written as BGZF: about a quarter of
    the plain size, and polars-bio reads it as fast as plain text. `structural_only` keeps only the meta
    lines polars-bio needs and the INFO and FORMAT keys the frames carry, for a file it refused (P7, P8):
    polars-bio skips a key its header does not declare (P1), so a malformed value in a key the frames
    do not carry is never parsed.
    """
    if layers in ((), ("gzip",), ("bgzf",)) and header.format_first and not missing and not structural_only:
        return path
    if structural_only:
        cache = cache / "structural"
    target = cache / "input.vcf.gz"
    done = cache / "input.vcf.gz.complete"
    if done.exists():
        return target
    cache.mkdir(parents=True, exist_ok=True)
    declared: set[tuple[str, str]] = set()
    with opener() as source, target.open("wb") as raw, BgzfWriter(raw) as out:
        out.write(f"##fileformat={header.version or 'VCFv4.2'}\n".encode())
        for line in header.meta:
            if line.startswith("##fileformat=") or (structural_only and not _keep_in_structural_copy(line)):
                continue
            key = re.match(r"^##(\w+)=<ID=([^,>]+)", line)
            if structural_only and key:
                if key.groups() in declared:
                    continue
                declared.add(key.groups())
            out.write(f"{line}\n".encode())
        for key in missing:
            number, kind = FORMAT_KEYS[key]
            out.write(
                f'##FORMAT=<ID={key},Number={number},Type={kind},'
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
        if ("Error reading" in str(err) and "field" in str(err)) or "Missing variant start" in str(err):
            raise HeaderRefused(str(err).splitlines()[-1][-160:]) from err
        raise


def _register(**kwargs) -> None:
    """`SQL._register_variant`, with a refused header raised as `HeaderRefused`."""
    try:
        SQL._register_variant(**kwargs)
    except ValueError as err:
        if "Failed to register table" in str(err):
            raise HeaderRefused(str(err)) from err
        raise


def read_records(path: Path, header: Header, sample: str | None, first_record: bytes | None) -> pl.LazyFrame:
    """The records of one sample (or none, for a sites-only file) with flat, namespaced columns.

    INFO keys become `info_<KEY>` and FORMAT keys `fmt_<KEY>`, whatever polars-bio names them. Every
    carried key is present, as a typed null when the file does not declare it.
    """
    info = [key for key in INFO_KEYS if key in header.info] + sorted(header.split_markers)
    fmt = [key for key in FORMAT_KEYS if key in header.format] if sample is not None else []
    multi = sample is not None and len(header.samples) > 1
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
    frame = frame.with_columns(reshaped)
    # A marker declared with a Number other than 1 comes back as a list: its first value is enough.
    frame = frame.with_columns(
        (pl.col(name).list.first() if isinstance(schema[name], pl.List) else pl.col(name)).cast(pl.String).alias(name)
        for name in (f"info_{key}" for key in header.split_markers)
        if name in present
    )
    return frame.with_columns(
        [pl.lit(None, dtype=_dtype(number, kind)).alias(f"fmt_{key}") for key, (number, kind) in FORMAT_KEYS.items() if f"fmt_{key}" not in present]
        + [pl.lit(None, dtype=pl.List(pl.Int32) if key == "SVLEN" else pl.Int32).alias(f"info_{key}") for key in INFO_KEYS if f"info_{key}" not in present]
    ).with_columns(pl.col("fmt_GT").cast(pl.String), pl.col("fmt_FT").cast(pl.String))


_BOUNDED = {"GQ": 2**16 - 1, "DP": 2**32 - 1, "MIN_DP": 2**32 - 1}


def scan_stats(records: pl.LazyFrame) -> dict:
    """One streaming pass for what the warnings and the facet need: out-of-range counts, GP's
    largest value, the source contig names, and whether any record is a reference block. It reads
    GT, which makes polars-bio parse the sample column, so a field value it cannot parse fails here."""
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
        _is_block_record(alts).any().alias("any_block"),
        pl.col("fmt_GP").cast(pl.List(pl.Float64), strict=False).is_not_null().sum().alias("gp_rows"),
        (pl.col("fmt_GP").cast(pl.List(pl.Float64), strict=False).list.sum() - 1.0).abs().gt(0.02).sum().alias("gp_off_one"),
        pl.col("chrom").unique().implode().alias("chroms"),
        # Reading GT makes polars-bio parse every sample field here, where a value it cannot parse
        # can still be retried (P7), rather than when the caller collects the frames.
        pl.col("fmt_GT").null_count().alias("gt_missing"),
        (_gt_body().str.count_matches(r"[/|]").max() + 1).alias("ploidy_max"),
        *[pl.col(name).sum().alias(name) for name in records.collect_schema().names() if name.startswith("_multi_")],
    )
    row = _collect(plan).row(0, named=True)
    row["chroms"] = sorted(row["chroms"] or [])
    return row


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
    """True on a record INFO `key` tags as split off a multi-allelic site: any vt value, or a
    bcftools value with its fifth, `USED_ALT_IDX`, field."""
    column = pl.col(f"info_{key}")
    return column.is_not_null() if kind == "vt" else column.str.count_matches(r"\|").fill_null(0) == 4


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


def to_rows(records: pl.LazyFrame, header: Header, stats: dict) -> tuple[pl.LazyFrame, pl.LazyFrame, list[str]]:
    """Lazy `calls` and `blocks` rows (without `facet_id`) from `read_records` output, plus warnings."""
    warnings = [
        f"{stats[f'out_{key}']} {key} value(s) outside 0..{top} set to null"
        for key, top in _BOUNDED.items()
        if stats[f"out_{key}"]
    ] + [
        f"{count} {name.removeprefix('_multi_')} value(s) held more than one number; the first was kept"
        for name, count in stats.items()
        if name.startswith("_multi_") and count
    ]
    v44 = (header.version_tuple or (0, 0)) >= (4, 4)
    # The file's REF goes to `source_ref` only. `ref` is read from the assembly sequence or stays null
    # (charter 10), so polars-bio's `ref` column must not reach the frame under that name.
    rows = records.with_columns(_alts_expr().alias("alts"), pl.col("ref").alias("source_ref")).drop("ref")
    rows = rows.with_columns(
        _end_pos(sv_declared="SVLEN" in header.info).alias("end_pos"),
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
    # RM270: a marker tags its own row here; the patterns need whole positions, read from the spill.
    rows = rows.with_columns(
        pl.when(_marker_expr(header)).then(pl.lit(SplitSignal.MARKER.value)).alias("split_signal")
    )
    rows = _phase_sets(rows)

    gp, gp_warning = _gp(header, stats)
    if gp_warning:
        warnings.append(gp_warning)
    rows = rows.with_columns(
        pl.coalesce(
            pl.col("fmt_GL").cast(pl.List(pl.Float32), strict=False),
            pl.col("fmt_PL").cast(pl.List(pl.Float32), strict=False).list.eval(pl.element() / -10.0),
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
        (pl.lit(None, dtype=pl.String) if v44 else gt_text).alias("source_gt"),
        canonical_chrom_expr("chrom").alias("chrom"),
        pl.col("start").alias("pos"),
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
    there is no ALT and INFO/END gives an extent (Illumina's gVCF: `ALT .` with `END`)."""
    all_star = alts.list.len().gt(0) & alts.list.eval(pl.element() == "<*>").list.all()
    illumina = alts.list.len().eq(0) & pl.col("info_END").is_not_null()
    return all_star | illumina


def _dtype(number: str, kind: str) -> pl.DataType:
    scalar = {"Integer": pl.Int32, "Float": pl.Float32, "String": pl.String}[kind]
    return scalar if number == "1" else pl.List(scalar)


def _end_pos(sv_declared: bool) -> pl.Expr:
    """VCF 4.5's computed END, per allele, then the largest over the site (SCHEMA.md `end_pos`)."""
    ref_end = pl.col("start") + pl.col("source_ref").str.len_chars().cast(pl.UInt32) - 1
    block_end = (
        pl.when(pl.col("fmt_LEN").is_not_null())
        .then(pl.col("start") + pl.col("fmt_LEN").cast(pl.UInt32) - 1)
        .otherwise(pl.col("info_END").cast(pl.UInt32))
    )
    sv_pattern = rf"^<({'|'.join(_SV_LENGTH_TYPES)})(:[^>]*)?>$"
    has_sv = pl.col("alts").list.eval(pl.element().str.contains(sv_pattern)).list.any()
    sv_end = (
        pl.struct("start", "alts", "info_SVLEN").map_elements(
            lambda row: _sv_end(row["start"], row["alts"], row["info_SVLEN"], sv_pattern), return_dtype=pl.UInt32
        )
        if sv_declared
        else pl.lit(None, dtype=pl.UInt32)
    )
    return pl.max_horizontal(
        ref_end,
        pl.when(pl.col("alts").list.contains("<*>")).then(block_end).otherwise(None),
        pl.when(pl.col("alts").list.len() == 0).then(pl.col("info_END").cast(pl.UInt32)).otherwise(None),
        pl.when(has_sv).then(sv_end).otherwise(None),
    ).cast(pl.UInt32)


def _sv_end(start: int, alts: list[str] | None, svlen: list[int | None] | None, pattern: str) -> int | None:
    alts = alts or []
    ends = [
        start + abs(length)
        for allele, length in zip(alts, svlen or [None] * len(alts), strict=False)
        if length is not None and re.match(pattern, allele)
    ]
    return max(ends) if ends else None


def _phase_sets(rows: pl.LazyFrame) -> pl.LazyFrame:
    """`phase_sets` and `phase_set_kinds`: the PS value on each phased haplotype, kind `ps`.

    Written without a join, as list arithmetic: a phased haplotype carries `PS + 1`, an unphased one
    0, which then becomes null; a null PS (the chromosome's implicit set) stays null.
    """
    encoded = pl.col("phased").list.eval(pl.element().cast(pl.Int64)) * (pl.col("fmt_PS").cast(pl.Int64) + 1)
    phased_any = pl.col("phased").list.eval(pl.element().fill_null(False)).list.any()
    return rows.with_columns(
        pl.when(phased_any)
        .then(encoded.list.eval(pl.when(pl.element() > 0).then((pl.element() - 1).cast(pl.String)).otherwise(None)))
        .otherwise(None)
        .alias("phase_sets"),
        pl.when(phased_any)
        .then(pl.col("phased").list.eval(pl.when(pl.element()).then(pl.lit(PhaseSetKind.PS.value)).otherwise(None)))
        .otherwise(None)
        .alias("phase_set_kinds"),
    )


def _gp(header: Header, stats: dict) -> tuple[pl.Expr, str | None]:
    """`gp` from GP or PP. GP's scale is read from the values when they settle it, else withheld.

    VCF before 4.3 wrote GP phred-scaled, and 4.3 made it a probability. Real files do not follow the
    version: imputation tools wrote probabilities under 4.1 and 4.2 headers. So values above 1 mean
    phred, and values that are all within [0, 1] and sum to 1 per row mean probabilities, whatever
    the header says; the warnings name a header that disagreed. Anything else is ambiguous, and GP
    is not used.
    """
    from_pp = pl.col("fmt_PP").cast(pl.List(pl.Float32), strict=False).list.eval(pl.lit(10.0) ** (pl.element() / -10.0))
    largest = stats["gp_max"]
    if largest is None:
        return from_pp, None
    raw = pl.col("fmt_GP").cast(pl.List(pl.Float32), strict=False)
    header_says = "phred" if (header.version_tuple or (4, 3)) < (4, 3) else "probabilities"
    if largest > 1.0 + 1e-6:
        found, gp = "phred", raw.list.eval(pl.lit(10.0) ** (pl.element() / -10.0))
    elif stats["gp_off_one"] == 0:
        found, gp = "probabilities", raw
    else:
        return from_pp, (
            f"GP is ambiguous: every value is within [0, 1] but {stats['gp_off_one']} of {stats['gp_rows']} "
            "rows do not sum to 1; GP not used"
        )
    note = None if found == header_says else f"GP written as {found} under a {header.version} header; read as {found}"
    return pl.coalesce(gp, from_pp), note


def _bounded(key: str, dtype: pl.DataType) -> pl.Expr:
    """FORMAT `key` cast to `dtype`; values outside it become null (counted by `scan_stats`)."""
    value = pl.col(f"fmt_{key}").cast(pl.Int64, strict=False)
    return pl.when((value >= 0) & (value <= _BOUNDED[key])).then(value).otherwise(None).cast(dtype)


def dropped_fields(header: Header) -> list[str]:
    """Declared keys the frames do not carry, namespaced."""
    info = [f"INFO/{key}" for key in header.info if key not in INFO_KEYS and key not in header.split_markers]
    fmt = [f"FORMAT/{key}" for key in header.format if key not in FORMAT_KEYS]
    return sorted(info) + sorted(fmt)


def to_frames(calls: pl.LazyFrame, blocks: pl.LazyFrame, facet_id: str) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    return (
        conform(calls.with_columns(pl.lit(facet_id).alias("facet_id")), CALLS),
        conform(blocks.with_columns(pl.lit(facet_id).alias("facet_id")), BLOCKS),
    )
