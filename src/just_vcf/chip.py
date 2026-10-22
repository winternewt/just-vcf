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

FamilyTreeDNA and MyHeritage exports state nothing and get no rule: every FamilyTreeDNA export in
PGP-HMS is on OmniExpress (by marker overlap with Illumina's manifest), so nothing here measures what a
later chip's export looks like, and naming one by its row count alone would be the guess this module
refuses (RM4).

The vocabulary is open (`CHIPS`, a `frozenset[str]` with a validator, stored as `str`): every vendor
release adds a chip. A Final Report's manifest is kept as Illumina names it, `illumina:<manifest>`.
"""

import re
from enum import Enum

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
        "living_dna_sirius",
    }
)
"""Chips named by vendor and generation. `23andme_v1_v2` is either of the two HumanHap550-based
generations: their exports are not told apart by their rows."""

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
