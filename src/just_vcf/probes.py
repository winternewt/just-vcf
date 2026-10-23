"""Chip manifests: what each probe of a genotyping chip tests (RM281).

A vendor manifest lists a chip's probes: name, position, the two alleles on the probe's design strand,
which strand that is against the reference, and the probe's source sequence with the variant between
its flanks (`…GAGA[-/AGCG]GCGA…`). Illumina publishes its manifests (docs/PROVENANCE.md); 23andMe and
AncestryDNA publish none. The manifests are not shipped: they are 70–80 MB zips under Illumina's terms,
so a configured directory holds them (`JUST_VCF_MANIFESTS_DIR`), each read once into a probe table in
the cache (`probes/<sha256>-<PROBES_VERSION>.parquet`, beside the reference store and, like it, not a
cache entry `clear_cache` manages). Nothing configured, or a manifest absent: the answers that need it
are unknown, and say why.

**Indel probes and an array's `I`/`D` codes.** An export writes an indel probe's genotype as `I` (the
allele with the inserted bases) and `D` (the allele without), which names no sequence. An Illumina
indel probe's source sequence names the inserted bases `S` and 60 or so bases of flank on each side.
`locate_indel` finds the site in the reference by those flanks, on either strand, within a window
around the manifest's position: the flanks with `S` between them say the reference carries `S` (a
deletion site, so `I` is REF), the flanks joined say it lacks it (an insertion site, `D` is REF). Only
one hit anywhere in the window is taken, and only on the strand the manifest's `SourceStrand` gives
the source sequence on; anything else is unlocated, counted. REF is then read from the reference, so
the manifest contributes the inserted bases and the locator, never REF (charter 10). `place_indel`
spells the site with REF starting at the export's own position, shifting it along a repeat where the
export placed it elsewhere in one.

Measured on GSA-24v3-0 A1 against Ensembl GRCh37 (2026-10-23): of its 10,118 indel probes, 9,618 are
found by their flanks at exactly one place; on the strand `SourceStrand` gives for all but 4. The
manifest's `MapInfo` is not the anchor: the gap lies one base before it for 5,594 probes, at or after
it, up to several bases on in a repeat, for the rest.

Coordinates: `map_info` and the exports' `pos` are 1-based; a gap is the 0-based index of the first
base after the left flank, so an insertion goes before that base and a deletion starts at it.
"""

import csv
import hashlib
import io
import logging
import os
import re
import zipfile
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import polars as pl

from just_vcf.cache import cache_dir, file_sha256
from just_vcf.schema import canonical_chrom

logger = logging.getLogger(__name__)

PROBES_VERSION = 1
"""Part of a probe table's file name: bumped whenever `probe_table` would read a manifest differently."""

MANIFEST_BUILDS: dict[str, str] = {"36": "GRCh36", "37": "GRCh37", "38": "GRCh38"}
"""A manifest's `GenomeBuild` column as an assembly label. A manifest states its build per probe."""

PROBES_SCHEMA: dict[str, pl.DataType] = {
    "name": pl.String,
    "rsid": pl.String,
    "chrom": pl.String,
    "map_info": pl.Int64,
    "assembly_label": pl.String,
    "kind": pl.String,
    "alleles": pl.List(pl.String),
    "ilmn_strand": pl.String,
    "ref_strand": pl.String,
    "source_strand": pl.String,
    "left_flank": pl.String,
    "indel_seq": pl.String,
    "right_flank": pl.String,
}
"""A probe table, one row per probe in manifest order. `rsid` is the rsid the probe's name gives
(`rs123`, `GSA-rs123`, `exm-rs123`), else null; `chrom` canonical, null for `Chr 0` (unmapped); `kind`
`snv`, `indel` or `other`; `alleles` the `SNP` column's two letters as written (design strand for an
SNV, `I`/`D` for an indel); the flanks and `indel_seq` only for an indel, on the source strand."""

_RSID_NAME = re.compile(r"^(?:GSA-|exm-)?(rs\d+)$")
_SOURCE_INDEL = re.compile(r"^([ACGTN]*)\[-/([ACGT]+)\]([ACGTN]*)$")
_COLUMNS = ("Name", "IlmnStrand", "SNP", "GenomeBuild", "Chr", "MapInfo", "SourceStrand", "SourceSeq", "RefStrand")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

SEARCH_MARGIN = 200
"""Bases searched beyond the probe's flanks on each side of `MapInfo` when locating an indel."""

MAX_PREFIX = 64
"""The most reference bases `place_indel` puts between an export's position and the site it names."""


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Manifest:
    """One manifest file in the configured directory."""

    descriptor: str
    """The `Descriptor File Name` its `[Heading]` gives, without `.bpm` (`GSA-24v3-0_A1`)."""
    path: Path
    sha256: str


def manifests_dir() -> Path | None:
    """`JUST_VCF_MANIFESTS_DIR`, or None when it is unset or empty."""
    configured = os.environ.get("JUST_VCF_MANIFESTS_DIR")
    return Path(configured) if configured else None


def available_manifests() -> list[Manifest]:
    """The manifest zips under the configured directory, by descriptor. Empty when none is configured."""
    root = manifests_dir()
    if root is None or not root.is_dir():
        return []
    return list(_index(str(root.resolve())).values())


def find_manifest(descriptor: str) -> Manifest | None:
    """The configured manifest of that descriptor (`GSA-24v3-0_A1`), or None."""
    root = manifests_dir()
    if root is None or not root.is_dir():
        return None
    return _index(str(root.resolve())).get(descriptor)


@cache
def _index(root: str) -> dict[str, Manifest]:
    found: dict[str, Manifest] = {}
    for path in sorted(Path(root).rglob("*.zip")):
        descriptor = _descriptor(path)
        if descriptor is None:
            continue
        if descriptor in found:
            logger.warning("manifests: %s and %s both hold %s; the first is read", found[descriptor].path.name, path.name, descriptor)
            continue
        found[descriptor] = Manifest(descriptor=descriptor, path=path, sha256=file_sha256(path))
    return found


def _member(archive: zipfile.ZipFile) -> str | None:
    members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
    return members[0] if len(members) == 1 else None


def _descriptor(path: Path) -> str | None:
    """The zip's manifest descriptor, read from its `[Heading]`; None for a zip that holds no manifest."""
    try:
        with zipfile.ZipFile(path) as archive:
            member = _member(archive)
            if member is None:
                return None
            with archive.open(member) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
                for line_number, row in enumerate(csv.reader(text)):
                    if row and row[0] == "Descriptor File Name" and len(row) > 1:
                        return row[1].strip().removesuffix(".bpm")
                    if line_number > 50 or (row and row[0].startswith("[Assay]")):
                        return None
    except zipfile.BadZipFile:
        return None
    return None


def probe_table(manifest: Manifest) -> pl.DataFrame:
    """The manifest's probes (`PROBES_SCHEMA`), read once into the cache and from there afterwards."""
    target = cache_dir() / "probes" / f"{manifest.sha256}-{PROBES_VERSION}.parquet"
    if target.is_file():
        return pl.read_parquet(target)
    table = _read_manifest(manifest.path)
    target.parent.mkdir(parents=True, exist_ok=True)
    staged = target.with_name(f".{target.name}.{os.getpid()}-{os.urandom(4).hex()}")
    table.write_parquet(staged)
    os.replace(staged, target)
    logger.info("manifests: %s read into the probe cache (%d probes)", manifest.descriptor, table.height)
    return table


def _read_manifest(path: Path) -> pl.DataFrame:
    rows: dict[str, list] = {name: [] for name in PROBES_SCHEMA}
    with zipfile.ZipFile(path) as archive:
        member = _member(archive)
        if member is None:
            raise ValueError(f"{path.name}: expected one CSV member")
        with archive.open(member) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            section, index = None, None
            for row in csv.reader(text):
                if not row:
                    continue
                if row[0].startswith("[") and row[0].rstrip().endswith("]"):
                    section, index = row[0].strip(), None
                    continue
                if section != "[Assay]":
                    continue
                if index is None:
                    index = {name: row.index(name) for name in _COLUMNS}
                    continue
                _append(rows, {name: row[position].strip() for name, position in index.items()})
    return pl.DataFrame(rows, schema=PROBES_SCHEMA)


def _append(rows: dict[str, list], probe: dict[str, str]) -> None:
    snp = probe["SNP"].strip("[]").split("/")
    indel = sorted(snp) == ["D", "I"]
    found = _SOURCE_INDEL.match(probe["SourceSeq"].upper()) if indel else None
    rsid = _RSID_NAME.match(probe["Name"])
    chrom = probe["Chr"]
    values = {
        "name": probe["Name"],
        "rsid": rsid.group(1) if rsid else None,
        "chrom": None if chrom in ("0", "") else canonical_chrom(chrom),
        "map_info": int(probe["MapInfo"]) if probe["MapInfo"].isdigit() else None,
        "assembly_label": MANIFEST_BUILDS.get(probe["GenomeBuild"]),
        "kind": "indel" if indel else ("snv" if len(snp) == 2 and all(len(a) == 1 and a in "ACGT" for a in snp) else "other"),
        "alleles": snp,
        "ilmn_strand": probe["IlmnStrand"] or None,
        "ref_strand": probe["RefStrand"] or None,
        "source_strand": probe["SourceStrand"] or None,
        "left_flank": found.group(1) if found else None,
        "indel_seq": found.group(2) if found else None,
        "right_flank": found.group(3) if found else None,
    }
    for name, value in values.items():
        rows[name].append(value)


def indel_catalog(assembly_label: str) -> tuple[pl.DataFrame, list[Manifest]]:
    """The indel probes of every configured manifest on `assembly_label` that an rsid names, and the
    manifests read. Columns: `rsid`, `chrom`, `map_info`, `left_flank`, `indel_seq`, `right_flank`,
    `source_strand`, `manifest`. A probe repeated across manifests is kept once."""
    parts, used = [], []
    for manifest in available_manifests():
        probes = probe_table(manifest).filter(
            (pl.col("kind") == "indel") & (pl.col("assembly_label") == assembly_label) & pl.col("rsid").is_not_null()
            & pl.col("chrom").is_not_null() & pl.col("map_info").is_not_null() & pl.col("indel_seq").is_not_null()
        )  # fmt: skip
        if probes.is_empty():
            continue
        used.append(manifest)
        parts.append(
            probes.select("rsid", "chrom", "map_info", "left_flank", "indel_seq", "right_flank", "source_strand").with_columns(
                pl.lit(manifest.descriptor).alias("manifest")
            )
        )
    if not parts:
        return pl.DataFrame(schema={"rsid": pl.String, "chrom": pl.String, "map_info": pl.Int64, "left_flank": pl.String,
                                    "indel_seq": pl.String, "right_flank": pl.String, "source_strand": pl.String,
                                    "manifest": pl.String}), used  # fmt: skip
    catalog = pl.concat(parts).unique(subset=["rsid", "chrom", "map_info", "left_flank", "indel_seq", "right_flank", "source_strand"], keep="first", maintain_order=True)
    return catalog, used


def catalog_digest(manifests: list[Manifest]) -> str:
    """What a result computed with these manifests depends on: their bytes and this module's version."""
    text = ",".join(sorted(manifest.sha256 for manifest in manifests)) + f";{PROBES_VERSION}"
    return hashlib.sha256(text.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class IndelSite:
    """An indel probe's site in the reference: the inserted bases on the plus strand, the gap before which
    they go (0-based), and whether the reference carries them (a deletion site) or not (an insertion site)."""

    gap: int
    seq: str
    reference_carries: bool


def window_for(map_info: int, left: str, seq: str, right: str, contig_length: int) -> tuple[int, int]:
    """The 0-based half-open reference span `locate_indel` searches for a probe."""
    reach = len(left) + len(seq) + len(right) + SEARCH_MARGIN
    return max(0, map_info - 1 - reach), min(contig_length, map_info + reach)


def locate_indel(window: str, start: int, left: str, seq: str, right: str, source_strand: str | None) -> IndelSite | str:
    """The site of an indel probe in `window` (reference bases from 0-based `start`), or why there is none:
    `not_found` (no flank match), `several` (more than one), `strand` (found only on the strand the
    manifest does not give its source sequence on)."""
    window = window.upper()
    hits: list[tuple[str, IndelSite]] = []
    for strand, (l, s, r) in (("+", (left, seq, right)), ("-", (reverse_complement(right), reverse_complement(seq), reverse_complement(left)))):
        for carries, probe in ((True, l + s + r), (False, l + r)):
            at = window.find(probe)
            while at != -1:
                hits.append((strand, IndelSite(gap=start + at + len(l), seq=s, reference_carries=carries)))
                at = window.find(probe, at + 1)
    if not hits:
        return "not_found"
    if len(hits) > 1:
        return "several"
    strand, site = hits[0]
    expected = {"PLUS": "+", "MINUS": "-"}.get((source_strand or "").upper())
    return site if strand == expected else "strand"


def place_indel(window: str, start: int, site: IndelSite, pos: int) -> tuple[str, str] | str:
    """The site spelled as `(REF, ALT)` with REF starting at 1-based `pos`, or why it cannot be:
    `before_pos` (the site lies before `pos` and no shift along a repeat brings it there), `far` (more
    than `MAX_PREFIX` bases after it), `window` (the window does not reach).

    The site is first shifted along a repeat towards `pos` (an equivalent spelling of the same
    haplotype), then the reference bases from `pos` to the gap are put before both alleles. Where none
    are, both alleles end with the base after the indel, so neither is empty: VCF's padding base,
    written after the event, which `normalize` moves to the left as it aligns.
    """
    window = window.upper()
    target = pos - 1
    gap, seq, n = site.gap, site.seq, len(site.seq)

    def base(at: int) -> str | None:
        return window[at - start] if 0 <= at - start < len(window) else None

    def right_shift(g: int, s: str) -> tuple[int, str] | None:
        if site.reference_carries:
            return (g + 1, window[g + 1 - start : g + 1 + n - start]) if base(g) is not None and base(g) == base(g + n) else None
        return (g + 1, s[1:] + s[0]) if base(g) is not None and s[0] == base(g) else None

    def left_shift(g: int, s: str) -> tuple[int, str] | None:
        if site.reference_carries:
            return (g - 1, window[g - 1 - start : g - 1 + n - start]) if base(g - 1) is not None and base(g - 1) == base(g - 1 + n) else None
        return (g - 1, s[-1] + s[:-1]) if base(g - 1) is not None and s[-1] == base(g - 1) else None

    while gap < target and (shifted := right_shift(gap, seq)) is not None:
        gap, seq = shifted
    while gap > target and (shifted := left_shift(gap, seq)) is not None:
        gap, seq = shifted
    if gap < target:
        return "before_pos"
    if gap - target > MAX_PREFIX:
        return "far"
    after = gap + n if site.reference_carries else gap
    if target < start or after + 1 > start + len(window):
        return "window"
    prefix = window[target - start : gap - start]
    pad = "" if prefix else window[after - start]
    long_allele, short_allele = prefix + seq + pad, prefix + pad
    return (long_allele, short_allele) if site.reference_carries else (short_allele, long_allele)

