"""Which genotyping chip an array export came from (RM4; charter invariant 9: unknown unless shown).

Two kinds of evidence, recorded apart in `facets.chip_basis`:

- **Stated**: the export's own header names the chip. AncestryDNA writes `array version: V1.0`, an
  Illumina GenomeStudio Final Report names its manifest on the `Content` line of its `[Header]`, and
  Living DNA writes `Genotype chip:`. A statement this module cannot map leaves the chip unknown and
  is quoted in a warning, so the vocabulary can grow from it.
- **Inferred**: a 23andMe export names no chip, but each generation's probe set has its own shape.
  Two features of the rows are read independently, the row count and the internal-id (`i…`) profile,
  and a generation is named only when both point at it. Measured on 341 public PGP-HMS 23andMe exports
  (docs/ROADMAP.md RM4). The bands come from just-prs `detect_chip_generation` (row counts per
  generation), re-measured; just-prs' fallback ("GSA v3" for any count it does not know) is not taken:
  a count no band holds is unknown. A table in 23andMe's layout without its preamble (a re-saved
  export) is read the same way when it carries 23andMe's internal ids.
- **Inferred by marker overlap**: FamilyTreeDNA and MyHeritage exports (`ftdna`, `ftdna_space`,
  `rsid_result_csv`) and a Final Report whose header names no manifest state nothing, so their rsids
  are compared with Illumina's own manifests (`chip_markers.json`, a sample of each manifest's rsids
  derived by `scripts/chip_markers.py`; docs/PROVENANCE.md). A chip is named only when one manifest
  holds nearly all of the export's sampled rsids, clearly more than any other does, and the export
  holds a good part of the manifest's: an extract, an export on a chip no manifest here describes, or
  two chips' content mixed stays unknown. Thresholds measured on PGP-HMS (docs/ROADMAP.md RM4). The
  name is the platform's family (`illumina_omniexpress`, `illumina_gsa`), not the manifest compared
  with: one manifest per family is on hand, and a Final Report that states GSA v2 reaches 0.97 against
  the GSA v3 manifest, so the overlap does not tell versions apart. A row count names nothing here.

The vocabulary is open (`CHIPS`, a `frozenset[str]` with a validator, stored as `str`): every vendor
release adds a chip. A Final Report's manifest is kept as Illumina names it, `illumina:<manifest>`.
"""

import json
import re
from enum import Enum
from functools import cache
from importlib import resources

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.schema import ChipBasis

CHIPS: frozenset[str] = frozenset(
    {
        "23andme_v1_v2",
        "23andme_v3",
        "23andme_v4",
        "23andme_v5",
        "ancestrydna_v1",
        "ancestrydna_v2",
        "illumina_gsa",
        "illumina_omniexpress",
        "living_dna_sirius",
    }
)
"""Chips named by vendor and generation. `23andme_v1_v2` is either of the two HumanHap550-based
generations: their exports are not told apart by their rows. `illumina_gsa` and
`illumina_omniexpress` are Illumina platform families, named only by marker overlap: the content of
the Global Screening Array or of the OmniExpress BeadChip, whichever version."""

ILLUMINA_PREFIX = "illumina:"
"""A chip given as the Illumina manifest a Final Report names (`illumina:GSA-24v2-0_A1`), verbatim."""


def is_chip(value: str) -> bool:
    """Whether `value` belongs to the chip vocabulary: one of `CHIPS`, or an `illumina:` manifest name."""
    if value in CHIPS:
        return True
    return value.startswith(ILLUMINA_PREFIX) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._\-]*", value[len(ILLUMINA_PREFIX) :]) is not None


_ANCESTRY_VERSIONS = {"V1.0": "ancestrydna_v1", "V2.0": "ancestrydna_v2"}
_LIVING_DNA_CHIPS = {"sirius": "living_dna_sirius"}

_ANCESTRY = re.compile(r"^#\s*Data was collected using AncestryDNA array version:\s*(\S+)\s*$")
_LIVING_DNA = re.compile(r"^#\s*Genotype chip:\s*(.+?)\s*$")
_CONTENT = re.compile(r"^Content[\t,]+(?:[\t,]*)([^\t,]+?)\s*[\t,]*$")


class ChipOutcome(str, Enum):
    """Why `detect_chip` gave the answer it gave. Only `STATED` and `INFERRED` carry a chip."""

    STATED = "stated"
    INFERRED = "inferred"
    UNMAPPED = "unmapped"
    """The header names a chip this module has no entry for."""
    NO_RULE = "no_rule"
    """The format states no chip and has no probe-set rule."""
    NO_BAND = "no_band"
    """The row count lies in no generation's band."""
    DISAGREE = "disagree"
    """The row count and the internal-id profile point at different generations, or the id profile at none."""
    TOO_FEW_MARKERS = "too_few_markers"
    """Too few of the export's rsids fall in the marker sample to compare it with a manifest."""
    NO_CLEAR_OVERLAP = "no_clear_overlap"
    """No manifest holds nearly all of the export's sampled rsids, clearly ahead of the others, with a good
    part of its own in the export."""


class ChipInference(BaseModel):
    """The chip an array export came from, or why it is unknown. Truthy only with a chip."""

    model_config = ConfigDict(frozen=True)

    chip: str | None
    basis: ChipBasis | None
    outcome: ChipOutcome
    reason: str

    def __bool__(self) -> bool:
        return self.chip is not None


def final_report_content(header: list[str]) -> str | None:
    """The manifest a Final Report's `Content` line names (`HumanOmni2.5-8v1_A.bpm`), as written."""
    for line in header:
        line = line.rstrip("\r\n")
        if line.startswith("[Data]"):
            break
        if found := _CONTENT.match(line):
            return found.group(1)
    return None


def stated_chip(format_id: str, header: list[str]) -> ChipInference | None:
    """The chip the header names, or None when it names none. `header` is the file's leading lines."""
    for line in header:
        line = line.rstrip("\r\n")
        if line.startswith("[Data]"):
            break
        if format_id == "ancestrydna" and (found := _ANCESTRY.match(line)):
            return _mapped(found.group(1), _ANCESTRY_VERSIONS.get(found.group(1)), "AncestryDNA array version")
        if format_id == "living_dna" and (found := _LIVING_DNA.match(line)):
            return _mapped(found.group(1), _LIVING_DNA_CHIPS.get(found.group(1).lower()), "Living DNA genotype chip")
        if format_id == "illumina_final_report" and (found := _CONTENT.match(line)):
            manifest = found.group(1).removesuffix(".bpm").removesuffix(".BPM")
            chip = ILLUMINA_PREFIX + manifest
            return _mapped(found.group(1), chip if is_chip(chip) else None, "Final Report Content")
    return None


def _mapped(stated: str, chip: str | None, where: str) -> ChipInference:
    if chip is None:
        return ChipInference(chip=None, basis=None, outcome=ChipOutcome.UNMAPPED, reason=f"the header's {where} {stated!r} names no known chip")
    return ChipInference(chip=chip, basis=ChipBasis.STATED, outcome=ChipOutcome.STATED, reason=f"the header's {where} is {stated!r}")


# --- 23andMe: the probe set -------------------------------------------------------------------------

ROW_BANDS_23ANDME: dict[str, tuple[int, int]] = {
    "23andme_v1_v2": (560_000, 597_000),
    "23andme_v4": (597_001, 625_000),
    "23andme_v5": (625_001, 660_000),
    "23andme_v3": (940_000, 1_020_000),
}
"""Rows in the export (no-calls included) per generation, inclusive. Measured on PGP-HMS (RM4)."""

LOW_ID_FLOOR = 1_000_000
"""Internal ids below this (`i700001` and up) are the v5 numbering; older generations start at `i1000001`."""


def id_profile(calls: pl.DataFrame) -> dict[str, float | int]:
    """The internal-id features of an export's rows: their count, their share of rows, and the share of
    them numbered below `LOW_ID_FLOOR`. Reads the first id of each row."""
    first = calls.select(pl.col("ids").list.first().alias("id"))["id"]
    internal = first.filter(first.str.contains(r"^i\d+$").fill_null(False)).str.slice(1).cast(pl.Int64)
    count = internal.len()
    low = int((internal < LOW_ID_FLOOR).sum()) if count else 0
    return {"rows": calls.height, "internal": count, "internal_share": count / calls.height if calls.height else 0.0, "low_share": low / count if count else 0.0}


def _generation_by_ids(profile: dict[str, float | int]) -> str | None:
    """The 23andMe generation the internal-id profile alone points at, or None.

    v5 numbers most of its internal ids from `i700001` (2.7–2.8% of its rows on PGP-HMS); v4 carries
    many (7.7–9.8%); v3 1.0–1.4%; v1 and v2 0.3–0.6%. The bands split each measured gap near its middle.
    """
    share, low = profile["internal_share"], profile["low_share"]
    if low >= 0.5:
        return "23andme_v5" if 0.01 <= share <= 0.06 else None
    if low > 0.1:
        return None
    if 0.05 <= share <= 0.15:
        return "23andme_v4"
    if 0.008 < share <= 0.03:
        return "23andme_v3"
    if share <= 0.008:
        return "23andme_v1_v2"
    return None


_23ANDME_FORMATS = frozenset({"23andme", "23andme_csv"})
_23ANDME_LAYOUT = frozenset({"rsid_table", "rsid_table_headerless"})
"""Tables in 23andMe's layout without its preamble, which is what a re-saved 23andMe export looks like.
They are read as 23andMe only when they carry its internal ids (`i…`), which no other vendor writes."""


def _generation_by_rows(rows: int) -> str | None:
    return next((chip for chip, (low, high) in ROW_BANDS_23ANDME.items() if low <= rows <= high), None)


def inferred_chip(format_id: str, calls: pl.DataFrame) -> ChipInference:
    """The chip the export's probe set points at, from rows a loader read (`ids` per row)."""
    if format_id in OVERLAP_FORMATS:
        return overlap_chip(calls)
    profile = id_profile(calls)
    if format_id not in _23ANDME_FORMATS and not (format_id in _23ANDME_LAYOUT and profile["internal"]):
        return ChipInference(chip=None, basis=None, outcome=ChipOutcome.NO_RULE, reason=f"{format_id} states no chip and has no probe-set rule")
    features = (
        f"{profile['rows']} rows, {profile['internal']} internal ids ({profile['internal_share']:.2%} of rows, "
        f"{profile['low_share']:.0%} of them below i{LOW_ID_FLOOR})"
    )
    by_rows = _generation_by_rows(int(profile["rows"]))
    by_ids = _generation_by_ids(profile)
    if by_rows is None:
        return ChipInference(chip=None, basis=None, outcome=ChipOutcome.NO_BAND, reason=f"{features}: the row count is in no generation's band")
    if by_ids != by_rows:
        return ChipInference(
            chip=None, basis=None, outcome=ChipOutcome.DISAGREE,
            reason=f"{features}: the row count points at {by_rows}, the internal ids at {by_ids or 'none'}",
        )  # fmt: skip
    return ChipInference(chip=by_rows, basis=ChipBasis.INFERRED, outcome=ChipOutcome.INFERRED, reason=f"{features}: row count and internal ids both point at {by_rows}")


# --- marker overlap with Illumina manifests ---------------------------------------------------------

OVERLAP_FORMATS = frozenset({"ftdna", "ftdna_space", "rsid_result_csv", "illumina_final_report"})
"""Formats whose chip is inferred by marker overlap: they state none (a Final Report only when its
header names no manifest; `detect_chip` takes a statement first). 23andMe exports keep their own rule,
which names the vendor's generation: 23andMe's v5 is a custom GSA, and the overlap would call it
`illumina_gsa`. AncestryDNA states its version."""

MARKER_MODULUS = 61
"""The marker sample: rsids whose number is divisible by this. A prime, so no decimal or binary pattern
in how rsids are numbered lines up with it. Applied alike to the manifests (`chip_markers.json`) and to
an export's ids, so the two samples compare as the full sets do. Measured on 606 PGP-HMS array exports
of 61,000 rsids or more, the sampled shares differ from the full sets' by at most 0.011 (RM4)."""

MIN_SAMPLED_MARKERS = 1_000
"""Sampled rsids an export needs before it is compared (about 61,000 rsids in the export). On PGP-HMS the
FTDNA X-chromosome exports have 245-270, every whole export 8,964 or more."""

MIN_COVERAGE = 0.90
"""Share of the export's sampled rsids the winning manifest must hold. On PGP-HMS: 0.98-0.99 for the 57
FTDNA exports on OmniExpress; at most 0.33 for the 2010 Affymetrix ones; 0.84 against GSA for the one
rsid CSV (vendor unstated), which therefore reads unknown."""

MIN_MARGIN = 0.30
"""How far the winning manifest's share must stand above every other manifest's (0.73 or more on PGP-HMS
where a chip is named: the two families share about a quarter of their rsids)."""

MIN_RECALL = 0.40
"""Share of the winning manifest's sampled rsids the export must hold: an extract of a chip is not the chip.
0.94-0.99 for the FTDNA exports on OmniExpress. Set below the 0.51 of the one GSA Final Report on PGP-HMS
(it states its chip, `GSA-24v2-0_A1`, and names half its probes by position, not by rsid), so a report
like it that did not state its chip would still be read."""

_RSID_NAME = r"^(?:GSA-|exm-)?rs(\d+)$"
"""An id naming an rsid, as the manifests and the Final Reports write it (`rs123`, `GSA-rs123`)."""


@cache
def chip_markers() -> dict[str, pl.Series]:
    """Each chip family's sampled rsid numbers, from the package's `chip_markers.json`."""
    payload = json.loads(resources.files("just_vcf").joinpath("chip_markers.json").read_text())
    if payload["modulus"] != MARKER_MODULUS:
        raise ValueError(f"chip_markers.json samples 1 in {payload['modulus']}, this module 1 in {MARKER_MODULUS}")
    return {chip: pl.Series(chip, entry["sample"], dtype=pl.UInt64) for chip, entry in sorted(payload["chips"].items())}


def sampled_rsids(calls: pl.DataFrame) -> pl.Series:
    """The distinct rsid numbers among the export's ids that fall in the marker sample."""
    numbers = calls.select(pl.col("ids").explode(empty_as_null=False).str.extract(_RSID_NAME, 1).cast(pl.UInt64).alias("rs"))["rs"]
    numbers = numbers.drop_nulls().unique()
    return numbers.filter(numbers % MARKER_MODULUS == 0).sort()


def marker_overlap(calls: pl.DataFrame) -> dict[str, tuple[float, float]]:
    """Per chip family: the share of the export's sampled rsids in its manifest (coverage), and the
    share of its manifest's sampled rsids in the export (recall). Empty when the export has none."""
    sampled = sampled_rsids(calls)
    if sampled.is_empty():
        return {}
    shares = {}
    for chip, markers in chip_markers().items():
        shared = int(sampled.is_in(markers.implode()).sum())
        shares[chip] = (shared / sampled.len(), shared / markers.len())
    return shares


def overlap_chip(calls: pl.DataFrame) -> ChipInference:
    """The chip family whose manifest the export's rsids clearly belong to, or why none is named."""
    count = sampled_rsids(calls).len()
    if count < MIN_SAMPLED_MARKERS:
        return ChipInference(
            chip=None, basis=None, outcome=ChipOutcome.TOO_FEW_MARKERS,
            reason=f"{count} of the export's rsids are in the marker sample (1 in {MARKER_MODULUS}), fewer than {MIN_SAMPLED_MARKERS}",
        )  # fmt: skip
    shares = marker_overlap(calls)
    ranked = sorted(shares.items(), key=lambda item: (-item[1][0], item[0]))
    (best, (coverage, recall)), runner_up = ranked[0], max((share[0] for _, share in ranked[1:]), default=0.0)
    features = f"{count} sampled rsids; " + ", ".join(f"{chip} holds {c:.1%} of them, the export {r:.1%} of its" for chip, (c, r) in ranked)
    if coverage >= MIN_COVERAGE and coverage - runner_up >= MIN_MARGIN and recall >= MIN_RECALL:
        return ChipInference(chip=best, basis=ChipBasis.INFERRED, outcome=ChipOutcome.INFERRED, reason=f"{features}: the overlap points at {best}")
    return ChipInference(chip=None, basis=None, outcome=ChipOutcome.NO_CLEAR_OVERLAP, reason=f"{features}: no manifest clearly holds the export")


def detect_chip(format_id: str, header: list[str], calls: pl.DataFrame) -> ChipInference:
    """The chip of an array export: stated by its header where it says, else inferred from its rows."""
    stated = stated_chip(format_id, header)
    return stated if stated is not None else inferred_chip(format_id, calls)


def header_lines(data: bytes, limit: int = 200) -> list[str]:
    """The payload's first `limit` lines, decoded leniently: where a header states a chip."""
    head = data[: 1 << 16]
    if head.startswith((b"\xff\xfe", b"\xfe\xff")):
        text = data[: 1 << 17].decode("utf-16", errors="replace")
    else:
        text = head.decode("utf-8", errors="replace").removeprefix("﻿")
    return text.replace("\r\n", "\n").replace("\r", "\n").split("\n")[:limit]
