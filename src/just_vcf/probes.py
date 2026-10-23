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

**Probe sets and the cause of an absence.** Which sites a chip has turns an array's absence, still an
unknown genotype, into one with a cause (`chip_absence`): `no_call` where the chip has a probe and the
export reports nothing, `not_on_chip` where it has none. A chip's probe set is its own manifest where a
Final Report names one (`illumina:<manifest>`, configured), or one learned from exports of the chip
(`learn_probe_sets`, written into the cache) for the vendors that publish none. A platform family named
by marker overlap has none: which version's probes it carries is not known.

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
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from functools import cache
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.cache import cache_dir, file_sha256
from just_vcf.chip import ILLUMINA_PREFIX
from just_vcf.loading import Loaded, load
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


_BASE_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def snv_catalog() -> tuple[pl.DataFrame, list[Manifest]]:
    """The two alleles of every configured manifest's rsid-named SNV probe on the plus strand, and the
    manifests read. Columns: `rsid`, `chrom`, `plus` (the pair, sorted).

    The `SNP` column gives a probe's alleles on its design strand; `RefStrand` says whether that strand
    is the reference's plus strand (`+`, the pair as written) or not (`-`, its complement). Checked on
    20,000 sampled SNV probes of each manifest against GRCh37: the base at `MapInfo` is in the pair so
    derived for all but 2 of GSA-24v3-0 A1's and 2 of OmniExpress-24 v1.2 A1's. A pair is a property of
    the variant, not of the build, so it is kept whatever build the manifest maps to; the reader checks
    it against its own reference base. An rsid whose probes give two different pairs is left out.
    """
    parts, used = [], []
    for manifest in available_manifests():
        probes = probe_table(manifest).filter(
            (pl.col("kind") == "snv") & pl.col("rsid").is_not_null() & pl.col("chrom").is_not_null() & pl.col("ref_strand").is_in(["+", "-"])
        )
        if probes.is_empty():
            continue
        used.append(manifest)
        plus = (
            pl.when(pl.col("ref_strand") == "+")
            .then(pl.col("alleles"))
            .otherwise(pl.col("alleles").list.eval(pl.element().replace_strict(_BASE_COMPLEMENT, default=None)))
        )
        parts.append(probes.select("rsid", "chrom", plus.list.sort().alias("plus")))
    if not parts:
        return pl.DataFrame(schema={"rsid": pl.String, "chrom": pl.String, "plus": pl.List(pl.String)}), used
    pairs = pl.concat(parts).unique()
    single = pairs.group_by("rsid", "chrom").len().filter(pl.col("len") == 1).drop("len")
    return pairs.join(single, on=["rsid", "chrom"], how="semi").sort("chrom", "rsid"), used


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


# --- probe sets: which sites a chip has ------------------------------------------------------------


class ProbeSetKind(str, Enum):
    """Where a chip's probe set comes from."""

    MANIFEST = "manifest"
    """The chip's own manifest: an Illumina chip a Final Report names (`illumina:<manifest>`), on disk."""
    LEARNED = "learned"
    """The sites the exports of this chip were seen to carry (`learn_probe_sets`), for a vendor that
    publishes no manifest. As complete as the exports it was learned from."""


class AbsenceCause(str, Enum):
    """What a site's absence from an array export means, given the chip's probe set. Null = unknown."""

    IN_EXPORT = "in_export"
    """Not absent: the export has a row at the site (a called genotype or a written no-call)."""
    NO_CALL = "no_call"
    """The chip has a probe at the site and the export reports nothing there: a no-call the export
    dropped, or a probe the vendor left out of this release. Unknown genotype, never reference."""
    NOT_ON_CHIP = "not_on_chip"
    """No probe of the chip is at the site: the export could not have said anything about it."""


class ProbeSet(BaseModel):
    """A chip's sites on one assembly: 1-based `pos` on canonical `chrom`."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    chip: str
    kind: ProbeSetKind
    source: str
    """The manifest's descriptor and sha256, or how many exports a learned set was read from."""
    assembly_label: str
    sites: pl.DataFrame


class ChipAbsence(BaseModel):
    """`chip_absence`'s answer: each asked site with its `absence_cause`, the probe set it was read
    against, and why the cause is null where it is."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    sites: pl.DataFrame
    probe_set: ProbeSet | None
    reason: str


FAMILY_CHIPS: frozenset[str] = frozenset({"illumina_gsa", "illumina_omniexpress"})
"""Chips named as a platform family by marker overlap (`just_vcf.chip`), whose version is not known."""

LEARNED_VERSION = 1
"""Part of a learned probe set's file name: bumped whenever `learn_probe_sets` would learn differently."""

MIN_LEARNING_EXPORTS = 3
"""Exports of one chip on one assembly a learned probe set needs: fewer cannot tell a chip's site from
one export's oddity, and the set is not written."""

MIN_EXPORT_SHARE = 0.5
"""Share of a chip's exports that must carry a site for it to be on the chip. Not the union: on PGP-HMS,
one of 81 23andMe v1/v2 exports labelled build 36 carried 578,794 positions no other export has (an
export in another coordinate system, or another chip), which the union would have put on the chip."""


def _on_chip(table: pl.DataFrame) -> pl.DataFrame:
    return table.filter(pl.col("exports") >= pl.col("of_exports") * MIN_EXPORT_SHARE).select("chrom", "pos")


def _learned_path(chip: str, assembly_label: str) -> Path:
    return cache_dir() / "probes" / "learned" / f"{chip}-{assembly_label}-{LEARNED_VERSION}.parquet"


def probe_set(chip: str | None, assembly_label: str | None) -> ProbeSet | str:
    """The probe set of `chip` on `assembly_label`, or why there is none.

    An `illumina:<manifest>` chip takes that manifest's mapped probes, when the manifest is configured.
    Any other chip takes the set `learn_probe_sets` wrote into the cache. A platform family named by
    marker overlap (`illumina_gsa`, `illumina_omniexpress`) has none: its version is not known, and
    another version's manifest would put probes on the chip that it may not carry.
    """
    if chip is None:
        return "the facet's chip is unknown"
    if assembly_label is None:
        return "the facet's assembly is unknown"
    if chip.startswith(ILLUMINA_PREFIX):
        descriptor = chip[len(ILLUMINA_PREFIX) :]
        manifest = find_manifest(descriptor)
        if manifest is None:
            return f"the manifest {descriptor} is not configured (JUST_VCF_MANIFESTS_DIR)"
        sites = probe_table(manifest).filter(
            (pl.col("assembly_label") == assembly_label) & pl.col("chrom").is_not_null() & pl.col("map_info").is_not_null()
        ).select("chrom", pl.col("map_info").alias("pos")).unique().sort("chrom", "pos")  # fmt: skip
        if sites.is_empty():
            return f"the manifest {descriptor} maps no probe to {assembly_label}"
        return ProbeSet(chip=chip, kind=ProbeSetKind.MANIFEST, source=f"{descriptor} (sha256 {manifest.sha256[:12]}…)",
                        assembly_label=assembly_label, sites=sites)  # fmt: skip
    if chip in FAMILY_CHIPS:
        return f"{chip} is a platform family named by marker overlap; which version's probes the chip carries is not known"
    path = _learned_path(chip, assembly_label)
    if not path.is_file():
        return f"no probe set for {chip} on {assembly_label} has been learned (learn_probe_sets)"
    table = pl.read_parquet(path)
    exports = int(table["of_exports"].max() or 0)
    return ProbeSet(chip=chip, kind=ProbeSetKind.LEARNED, source=f"learned from {exports} export(s)", assembly_label=assembly_label,
                    sites=_on_chip(table))  # fmt: skip


def chip_absence(loaded: Loaded, sites: pl.DataFrame) -> ChipAbsence:
    """Why each of `sites` (`chrom` canonical, `pos` 1-based) is or is not in an array facet's export.

    `absence_cause` is `in_export` where the export has a row there, else `no_call` where the facet's
    chip has a probe there and `not_on_chip` where it has none; null for every site when the facet has
    no probe set (not an array, chip or assembly unknown, nothing configured or learned), with the
    reason. Either absence is an unknown genotype (Decision 4); the cause says why it is unknown.
    """
    facet = loaded.facets.collect().row(0, named=True)
    asked = sites.select(pl.col("chrom").cast(pl.String), pl.col("pos").cast(pl.Int64)).unique(maintain_order=True)
    if facet["vcf_version"] is not None or facet["callset_type"] not in (None, "array"):
        return ChipAbsence(sites=asked.with_columns(pl.lit(None, dtype=pl.String).alias("absence_cause")), probe_set=None,
                           reason="the facet is not an array export")  # fmt: skip
    found = probe_set(facet["chip"], facet["assembly_label"])
    present = loaded.calls.select("chrom", pl.col("pos").cast(pl.Int64)).unique().collect().with_columns(pl.lit(True).alias("_in"))
    marked = asked.join(present, on=["chrom", "pos"], how="left", maintain_order="left").with_columns(pl.col("_in").fill_null(False))
    if isinstance(found, str):
        cause = pl.when(pl.col("_in")).then(pl.lit(AbsenceCause.IN_EXPORT.value)).otherwise(None)
        return ChipAbsence(sites=marked.with_columns(cause.alias("absence_cause")).drop("_in"), probe_set=None, reason=found)
    on_chip = marked.join(found.sites.with_columns(pl.lit(True).alias("_on")), on=["chrom", "pos"], how="left", maintain_order="left")
    cause = (
        pl.when(pl.col("_in")).then(pl.lit(AbsenceCause.IN_EXPORT.value))
        .when(pl.col("_on").fill_null(False)).then(pl.lit(AbsenceCause.NO_CALL.value))
        .otherwise(pl.lit(AbsenceCause.NOT_ON_CHIP.value))
    )  # fmt: skip
    return ChipAbsence(sites=on_chip.with_columns(cause.alias("absence_cause")).drop("_in", "_on"), probe_set=found,
                       reason=f"read against the {found.kind.value} probe set of {found.chip} on {found.assembly_label}: {found.source}")  # fmt: skip


def learn_probe_sets(exports: Iterable[Path | str | Loaded], min_exports: int = MIN_LEARNING_EXPORTS) -> list[ProbeSet]:
    """Learn each chip's sites from array exports of it, and write them into the cache for `probe_set`.

    For a vendor that publishes no manifest (23andMe, AncestryDNA, Living DNA), the sites a chip has are
    the positions its exports carry: 23andMe and AncestryDNA write a row for every probe, no-calls
    included. Each export counts under its facet's `chip` and `assembly_label`; an export with either
    unknown, a VCF, a platform family or a chip with its own manifest (`illumina:`) is passed over. A
    set holds every position its exports carry, each with the number of exports that carry it
    (`exports`; `of_exports` is how many were read), and is written only from `min_exports` exports
    or more; a site is on the chip when at least `MIN_EXPORT_SHARE` of them carry it. Learning again
    replaces the set. Returns the sets written.
    """
    seen: dict[tuple[str, str], pl.DataFrame] = {}
    counted: dict[tuple[str, str], int] = {}
    for export in exports:
        loaded = export if isinstance(export, Loaded) else load(export)
        facet = loaded.facets.collect().row(0, named=True)
        chip, label = facet["chip"], facet["assembly_label"]
        if facet["vcf_version"] is not None or chip is None or label is None or chip in FAMILY_CHIPS or chip.startswith(ILLUMINA_PREFIX):
            continue
        sites = loaded.calls.select("chrom", pl.col("pos").cast(pl.Int64)).unique().collect().with_columns(pl.lit(1, dtype=pl.UInt32).alias("exports"))
        key = (chip, label)
        counted[key] = counted.get(key, 0) + 1
        before = seen.get(key)
        seen[key] = sites if before is None else pl.concat([before, sites]).group_by("chrom", "pos").agg(pl.col("exports").sum())
    written = []
    for (chip, label), sites in sorted(seen.items()):
        if counted[(chip, label)] < min_exports:
            logger.info("probe sets: %s on %s seen in %d export(s), %d needed; not written", chip, label, counted[(chip, label)], min_exports)
            continue
        table = sites.sort("chrom", "pos").with_columns(pl.lit(counted[(chip, label)], dtype=pl.UInt32).alias("of_exports"))
        target = _learned_path(chip, label)
        target.parent.mkdir(parents=True, exist_ok=True)
        staged = target.with_name(f".{target.name}.{os.getpid()}-{os.urandom(4).hex()}")
        table.write_parquet(staged)
        os.replace(staged, target)
        written.append(ProbeSet(chip=chip, kind=ProbeSetKind.LEARNED, source=f"learned from {counted[(chip, label)]} export(s)",
                                assembly_label=label, sites=_on_chip(table)))  # fmt: skip
    return written
