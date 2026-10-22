"""What kind of callset a VCF facet is, and so what an absent position means (RM4; charter 9, 11).

The scope test is ported from just-dna-lite `annotation/restoration.py` `detect_callset_scope`
(docs/PROVENANCE.md): a whole-genome callset carries millions of sites, and nearly all of its span
lies within one flank of a call. That file separated whole-genome from "targeted" only. Here the
targeted side is split into exome and panel, and every band between two classes answers unknown, so
weak evidence never lands in the nearest class.

The classes, by what the calls look like:

- `wgs`: at least `MIN_WGS_SITES` sites on the primary assembly and at least `MIN_WGS_BREADTH` of the
  span within `FLANK_BP` of a call (the ported rule, thresholds unchanged).
- `wes`: an exome's spread. Between `MIN_WES_SITES` and `MAX_WES_SITES` sites, on every autosome, and
  clustered (breadth at most `MAX_TARGETED_BREADTH`).
- `panel`: calls confined to a small region set. Between `MIN_PANEL_SITES` and `MAX_PANEL_SITES` sites,
  all within `MAX_PANEL_REACH_BP` of genome. A region extract of a larger callset (`bcftools view -r`)
  looks the same and cannot be told apart from the calls: either way absence outside the regions means
  nothing was looked at, which is what the facet then says.
- `imputed`: the header names an imputation tool or declares an imputation INFO key. Imputed sites are
  a panel's sites, whatever their spread.

Each sequencing class also needs the header to declare a read-evidence FORMAT key (AD, DP, GQ, PL, GL,
MIN_DP): a genome-wide file of genotypes with no read evidence at all is what imputation and array
conversions write, so its spread alone is not taken as a whole-genome callset.
"""

import re
from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.schema import AbsenceMeans, CallsetType

FLANK_BP = 10_000
"""just-dna-lite's `restoration_max_flank_bp` default, the flank the breadth is measured at."""

MIN_WGS_SITES = 1_000_000
"""just-dna-lite `MIN_WGS_SITES`: about 4x below the WGS floor it measured (4.3–4.7M sites)."""

MIN_WGS_BREADTH = 0.75
"""just-dna-lite `MIN_WGS_BREADTH`: WGS measured 0.942–0.950 there, a clustered callset 0.21."""

MIN_WES_SITES = 10_000
MAX_WES_SITES = 400_000
"""An exome carries roughly 20,000–150,000 sites, more with off-target calls; the bounds sit below and
above that with room on both sides. 400,000 to 1,000,000 is unknown."""

MAX_TARGETED_BREADTH = 0.6
"""An exome's breadth stays well under a genome's 0.94; 0.6 to 0.75 is unknown."""

MIN_PANEL_SITES = 50
MAX_PANEL_SITES = 5_000
"""Fewer than 50 sites (a test file, a single-gene extract) say too little to name a class; 5,000 to
10,000 is the unknown band between a large panel and a small exome."""

MAX_PANEL_REACH_BP = 50_000_000
"""Genome within `FLANK_BP` of a panel's calls. A few hundred genes, each with its flanks, stay far
below it; an exome's calls reach hundreds of megabases."""

READ_EVIDENCE_KEYS = frozenset({"AD", "DP", "GQ", "PL", "GL", "MIN_DP"})
"""FORMAT keys a caller that read the sample's reads writes; imputation and array conversions do not."""

IMPUTATION_INFO_KEYS = frozenset({"R2", "ER2", "DR2", "IMP", "IMPUTED", "TYPED"})
"""INFO keys minimac (R2, ER2, IMPUTED, TYPED) and Beagle (DR2, IMP) declare."""

_IMPUTATION_SOURCE = re.compile(r"^##source=.*(?i:minimac|beagle|impute[245]?|glimpse|loimpute)", re.MULTILINE)
"""A `##source` line naming an imputation tool (Gencove writes `loimpute`)."""

_PRIMARY = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]
_AUTOSOMES = frozenset(str(n) for n in range(1, 23))


class CallsetInference(BaseModel):
    """The callset type the evidence supports, or why there is none. Truthy only with a type."""

    model_config = ConfigDict(frozen=True)

    callset_type: CallsetType | None
    reason: str
    sites: int | None = None
    """Distinct `(chrom, pos)` on the primary assembly; None when the spread was not measured."""
    breadth: float | None = None
    reach_bp: int | None = None
    autosomes: int | None = None

    def __bool__(self) -> bool:
        return self.callset_type is not None


def imputation_marker(meta: list[str], info_keys: list[str]) -> str | None:
    """What in the header says the file is imputation output, or None."""
    declared = sorted(IMPUTATION_INFO_KEYS & set(info_keys))
    if declared:
        return f"the header declares imputation INFO key(s) {', '.join(declared)}"
    for line in meta:
        if _IMPUTATION_SOURCE.match(line):
            return f"the header names an imputation tool ({line[:60]})"
    return None


def spread(calls: pl.LazyFrame, flank_bp: int = FLANK_BP) -> dict:
    """Site count, breadth, reach and autosome count of a `calls` frame's small-variant sites.

    A site is a distinct `(chrom, pos)` on the primary assembly (canonical `chrom`) of a row that
    writes a sequence allele, or none (a monomorphic record). Rows whose ALTs are all symbolic
    (`<DEL>`, `<*>`, breakends, `*`) do not count: a genome-wide structural-variant callset says
    nothing about small variants, and must not read as a whole-genome one. Alt contigs and decoys do
    not count either; they carry a handful of calls each.

    Breadth is just-dna-lite's `_coverage_breadth`: each gap between neighbouring sites counts up to
    `2 * flank_bp`, summed (the reach), divided by the span the sites spread over per contig.
    """
    sequence_allele = pl.col("alts").list.eval(pl.element().str.contains(r"^[A-Za-z]+$")).list.any()
    primary = (
        calls.filter(pl.col("chrom").is_in(_PRIMARY) & ((pl.col("alts").list.len() == 0) | sequence_allele))
        .select("chrom", "pos")
        .unique()
    )
    per_contig = (
        primary.group_by("chrom")
        .agg(
            pl.len().alias("n"),
            (pl.col("pos").max().cast(pl.Int64) - pl.col("pos").min().cast(pl.Int64)).alias("span"),
            pl.col("pos").sort().cast(pl.Int64).diff().clip(upper_bound=2 * flank_bp).sum().alias("reach"),
        )
        .collect()
    )
    total = int(per_contig["n"].sum() or 0)
    span = int(per_contig["span"].sum() or 0)
    reach = int(per_contig["reach"].sum() or 0)
    return {
        "sites": total,
        "breadth": reach / span if span else 0.0,
        "reach_bp": reach,
        "autosomes": int(per_contig.filter(pl.col("chrom").is_in(_AUTOSOMES)).height),
    }


def classify(calls: pl.LazyFrame, format_keys: list[str], imputed_by: str | None) -> CallsetInference:
    """The callset type of a VCF facet from its `calls` (`chrom`, `pos`, `alts`) and its header."""
    if imputed_by:
        return CallsetInference(callset_type=CallsetType.IMPUTED, reason=imputed_by)
    measured = spread(calls)
    n, breadth, reach, autosomes = measured["sites"], measured["breadth"], measured["reach_bp"], measured["autosomes"]
    where = f"{n:,} site(s) on the primary assembly, {breadth:.1%} of the span within {FLANK_BP:,} bp of a call, {reach:,} bp reached, {autosomes} autosome(s)"

    def answer(callset_type: CallsetType | None, why: str) -> CallsetInference:
        return CallsetInference(callset_type=callset_type, reason=f"{why}: {where}", **measured)

    evidence = sorted(READ_EVIDENCE_KEYS & set(format_keys))
    if not evidence:
        return answer(None, "unknown, the header declares no read-evidence FORMAT key (AD, DP, GQ, PL, GL, MIN_DP), so the calls may not come from reads")
    if n >= MIN_WGS_SITES and breadth >= MIN_WGS_BREADTH:
        return answer(CallsetType.WGS, "whole genome")
    if MIN_WES_SITES <= n <= MAX_WES_SITES and breadth <= MAX_TARGETED_BREADTH and autosomes == len(_AUTOSOMES):
        return answer(CallsetType.WES, "exome")
    if MIN_PANEL_SITES <= n <= MAX_PANEL_SITES and reach <= MAX_PANEL_REACH_BP:
        return answer(CallsetType.PANEL, "panel or region extract")
    return answer(None, "unknown, the spread fits no class with room to spare")


def absence_means(callset_type: CallsetType | None, has_blocks: bool) -> AbsenceMeans:
    """What a position absent from the facet means (SCHEMA.md Decision 4).

    Reference blocks decide it wherever they exist. Without them only a whole-genome callset makes
    absence weak evidence of reference; an exome, a panel, an array or imputed content says nothing
    outside what it looked at, and neither does a callset of unknown type.
    """
    if has_blocks:
        return AbsenceMeans.BLOCKS
    return AbsenceMeans.WEAK_REFERENCE if callset_type is CallsetType.WGS else AbsenceMeans.UNKNOWN
