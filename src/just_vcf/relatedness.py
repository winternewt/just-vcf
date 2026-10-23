"""How closely two facets' genotypes say their people are related (RM274), on the frames, in Polars.

The estimator is KING-robust (Manichaikul et al. 2010, "Robust relationship inference in genome-wide
association studies", Bioinformatics 26:2867), the between-family form PLINK 2's `--make-king` uses:

    kinship = 1/2 - (H_a + H_b) / (4 H_min) + (H_ab - 2 O) / (2 H_min)

over the panel sites both facets have a genotype at, where `H_a` and `H_b` count each facet's
heterozygous sites, `H_min` is the smaller of the two, `H_ab` counts sites heterozygous in both and `O`
sites homozygous for opposite alleles (IBS0). It needs no allele frequency, so a person's ancestry does
not move it the way it moves an estimator that does. On this scale one person (or a monozygotic twin)
is 0.5, a parent, child or full sibling 0.25, a grandparent, half sibling, aunt or uncle 0.125, a first
cousin 0.0625, and unrelated people about 0.

**Sites.** A panel of autosomal, strand-unambiguous SNVs on one assembly (`SitePanel`). The package
carries one for GRCh37 and one for GRCh38, derived from somalier's published site VCFs by
`scripts/relatedness_sites.py` (docs/PROVENANCE.md): common sites, so an array, an exome and a whole
genome all have a share of them. Each facet's genotype at a site is read by
`just_vcf.site_genotypes`, the same reading population inference uses: a row's literal alleles against
the panel's two, an absence as homozygous reference only where the facet's `absence_means` says so.
Sites are matched by position and alleles on one assembly; two facets whose `assembly_label`s differ,
or that have none, are not compared (`assembly_mismatch`, `assembly_unknown`). Where both facets carry
minted VRS ids for the panel's ALT at a site (`mint_vrs_ids`), the ids are compared too, and a site
where they differ is counted (`n_vrs_differ`): on one assembly they never should.

**Bands.** KING's inference criteria, on the kinship scale: above 2^-1.5 (0.354) one person, from 2^-2.5
(0.177) first degree, from 2^-3.5 (0.0884) second degree, from 2^-4.5 (0.0442) third degree, below that
unrelated. A band is named only when the kinship's interval (two jackknife standard errors either side,
over `JACKKNIFE_BLOCKS` contiguous blocks of panel sites) lies inside it; an interval across a boundary
is `unknown` with the reason `between_bands`; a facet homozygous reference at none of the shared sites
chose them by its own genotype, and the band is `unknown`, `variant_sites_only`, its kinship kept.
The bands are just-vcf's own vocabulary: GA4GH Pedigree's KIN ontology names relationships (parent
of, sibling of), not degrees inferred from genotypes.

A facet's genotypes say nothing about who someone is: this measures only how alike two facets'
genotypes are. Same-person (`same_person`) is the check `just_vcf.person` runs on facets a caller says
are one person.
"""

import json
import math
from collections.abc import Sequence
from enum import Enum
from functools import cache
from importlib import resources

import polars as pl
from pydantic import BaseModel, ConfigDict, field_validator

from just_vcf.loading import Loaded
from just_vcf.site_genotypes import AUTOSOMES, Facet, facet_of, facets_of, read_site_genotypes

PANEL_VERSION = 1
"""The layout of `relatedness_sites.json`; `scripts/relatedness_sites.py` writes it."""

PANEL_COLUMNS = ("chrom", "pos", "ref_allele", "effect_allele")

KING_METHOD = "king_robust"
SOMALIER_METHOD = "somalier"
RELATEDNESS_METHODS: frozenset[str] = frozenset({KING_METHOD, SOMALIER_METHOD})

ASSEMBLY_UNKNOWN = "assembly_unknown"
"""A facet has no `assembly_label`: its sites cannot be placed on a panel."""
ASSEMBLY_MISMATCH = "assembly_mismatch"
"""The two facets are on different assemblies: positions do not match across them."""
NO_PANEL = "no_panel"
"""No panel for the facets' assembly: the package carries none and the caller passed none, or passed
one on another assembly."""
TOO_FEW_SITES = "too_few_sites"
"""Fewer than `min_sites` panel sites where both facets have a genotype."""
NO_HETEROZYGOUS_SITES = "no_heterozygous_sites"
"""A facet is heterozygous at none of the shared sites: KING-robust divides by that count."""
BETWEEN_BANDS = "between_bands"
"""The kinship's interval crosses a band boundary."""
VARIANT_SITES_ONLY = "variant_sites_only"
"""A facet is homozygous reference at none of the shared sites: it has a genotype only where it carries
the ALT (variant rows, and a silence that says nothing), so its person's genotype chose the sites, which
biases KING-robust in a way the jackknife SE does not capture (up, on the exome measured). The kinship and its SE are
kept; the band is not named. docs/RELATEDNESS.md has the measurement (a Genos exome gVCF whose blocks are all no-call)."""
RELATEDNESS_REASONS: frozenset[str] = frozenset(
    {ASSEMBLY_UNKNOWN, ASSEMBLY_MISMATCH, NO_PANEL, TOO_FEW_SITES, NO_HETEROZYGOUS_SITES, BETWEEN_BANDS, VARIANT_SITES_ONLY}
)
"""Why a relatedness is `unknown`. Grows additively."""

DEFAULT_MIN_SITES = 500
"""A guard against a kinship from almost nothing, not a calibrated threshold: the jackknife standard
error is what says how far the sites support the value."""
JACKKNIFE_BLOCKS = 50
INTERVAL_SES = 2.0
"""How many standard errors either side of the kinship the band's interval reaches."""

SAME_PERSON_FLOOR = 2**-1.5
FIRST_DEGREE_FLOOR = 2**-2.5
SECOND_DEGREE_FLOOR = 2**-3.5
THIRD_DEGREE_FLOOR = 2**-4.5


class RelationBand(str, Enum):
    """The degree of relationship a kinship falls in (KING's inference criteria). Just-vcf's own terms."""

    SAME_PERSON = "same_person"
    """One person, or a monozygotic twin: genotypes cannot tell the two apart."""
    FIRST_DEGREE = "first_degree"
    SECOND_DEGREE = "second_degree"
    THIRD_DEGREE = "third_degree"
    UNRELATED = "unrelated"
    UNKNOWN = "unknown"


def band_of(kinship: float) -> RelationBand:
    """The band a kinship value falls in."""
    if kinship > SAME_PERSON_FLOOR:
        return RelationBand.SAME_PERSON
    if kinship > FIRST_DEGREE_FLOOR:
        return RelationBand.FIRST_DEGREE
    if kinship > SECOND_DEGREE_FLOOR:
        return RelationBand.SECOND_DEGREE
    if kinship > THIRD_DEGREE_FLOOR:
        return RelationBand.THIRD_DEGREE
    return RelationBand.UNRELATED


class SitePanel(BaseModel):
    """Autosomal, strand-unambiguous biallelic SNVs on one assembly, `ref_allele` the assembly's REF."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    assembly_label: str
    sites: pl.DataFrame
    """`chrom, pos, ref_allele, effect_allele`, canonical `chrom`, one site per position, in panel order."""
    source: str

    @field_validator("sites")
    @classmethod
    def _valid(cls, sites: pl.DataFrame) -> pl.DataFrame:
        missing = [name for name in PANEL_COLUMNS if name not in sites.columns]
        if missing:
            raise ValueError(f"a panel needs the columns {', '.join(PANEL_COLUMNS)}; missing {', '.join(missing)}")
        sites = sites.select(PANEL_COLUMNS)
        bases = pl.lit(["A", "C", "G", "T"])
        problems = {
            "not on an autosome": ~pl.col("chrom").is_in(list(AUTOSOMES)),
            "not a one-base ACGT SNV": ~(pl.col("ref_allele").is_in(bases) & pl.col("effect_allele").is_in(bases))
            | (pl.col("ref_allele") == pl.col("effect_allele")),
            "strand-ambiguous (A/T or C/G)": pl.concat_str("ref_allele", "effect_allele").is_in(["AT", "TA", "CG", "GC"]),
            "two sites at one position": pl.struct("chrom", "pos").is_duplicated(),
        }
        counts = sites.select(**{name: expr.sum() for name, expr in problems.items()}).row(0, named=True)
        found = {name: count for name, count in counts.items() if count}
        if found:
            raise ValueError("not a relatedness panel: " + "; ".join(f"{count} sites {name}" for name, count in found.items()))
        return sites


@cache
def _packaged() -> dict:
    payload = json.loads(resources.files("just_vcf").joinpath("relatedness_sites.json").read_text())
    if payload["version"] != PANEL_VERSION:
        raise ValueError(f"relatedness_sites.json has layout {payload['version']}, this module reads {PANEL_VERSION}")
    return payload


def packaged_panel(assembly_label: str) -> SitePanel | None:
    """The panel the package carries for `assembly_label` (GRCh37, GRCh38), or None."""
    payload = _packaged()
    entry = payload["assemblies"].get(assembly_label)
    if entry is None:
        return None
    sites = pl.DataFrame(
        {
            "chrom": entry["chrom"],
            "pos": pl.Series(entry["pos"], dtype=pl.UInt32),
            "ref_allele": list(entry["ref"]),
            "effect_allele": list(entry["alt"]),
        }
    )
    return SitePanel(assembly_label=assembly_label, sites=sites, source=f"{payload['source']}: {entry['file']}")


def packaged_assemblies() -> tuple[str, ...]:
    """The assemblies the package carries a panel for."""
    return tuple(sorted(_packaged()["assemblies"]))


class Relatedness(BaseModel):
    """How alike two facets' genotypes are, and the relationship band that puts them in, or unknown."""

    model_config = ConfigDict(frozen=True)

    facet_a: str
    facet_b: str
    method: str
    band: RelationBand
    reason: str | None = None
    """Why `band` is `unknown`; one of `RELATEDNESS_REASONS`. None otherwise."""
    assembly_label: str | None = None
    panel: str | None = None
    kinship: float | None = None
    kinship_se: float | None = None
    """Delete-one-block jackknife standard error over `JACKKNIFE_BLOCKS` blocks of panel sites."""
    ibs0: float | None = None
    """Share of the shared sites where the two are homozygous for opposite alleles."""
    ibs2: float | None = None
    """Share of the shared sites where the two genotypes are the same."""
    n_sites: int | None = None
    """Panel sites where both facets have a genotype."""
    n_het_a: int | None = None
    n_het_b: int | None = None
    n_het_both: int | None = None
    n_ibs0: int | None = None
    n_ibs2: int | None = None
    n_vrs_compared: int | None = None
    """Shared sites where both facets carry a minted VRS id for the panel's ALT."""
    n_vrs_differ: int | None = None
    somalier_relatedness: float | None = None
    """somalier's own `relatedness` column, verbatim (twice the kinship scale), for `method = somalier`."""

    @field_validator("method")
    @classmethod
    def _method(cls, value: str) -> str:
        if value not in RELATEDNESS_METHODS:
            raise ValueError(f"unknown relatedness method {value!r}")
        return value

    @field_validator("reason")
    @classmethod
    def _reason(cls, value: str | None) -> str | None:
        if value is not None and value not in RELATEDNESS_REASONS:
            raise ValueError(f"unknown relatedness reason {value!r}")
        return value


RELATEDNESS_FRAME = pl.Schema(
    {
        "facet_a": pl.String,
        "facet_b": pl.String,
        "method": pl.String,
        "band": pl.Enum([band.value for band in RelationBand]),
        "reason": pl.String,
        "assembly_label": pl.String,
        "panel": pl.String,
        "kinship": pl.Float64,
        "kinship_se": pl.Float64,
        "ibs0": pl.Float64,
        "ibs2": pl.Float64,
        "n_sites": pl.UInt32,
        "n_het_a": pl.UInt32,
        "n_het_b": pl.UInt32,
        "n_het_both": pl.UInt32,
        "n_ibs0": pl.UInt32,
        "n_ibs2": pl.UInt32,
        "n_vrs_compared": pl.UInt32,
        "n_vrs_differ": pl.UInt32,
        "somalier_relatedness": pl.Float64,
    }
)
"""The frame `pairwise_relatedness` returns: one row per pair of facets, the fields of `Relatedness`."""


def king_robust(n_het_a: int, n_het_b: int, n_het_both: int, n_ibs0: int) -> float | None:
    """KING-robust's between-family kinship from the four counts, or None when a facet has no het site."""
    smaller = min(n_het_a, n_het_b)
    if smaller == 0:
        return None
    return 0.5 - (n_het_a + n_het_b) / (4 * smaller) + (n_het_both - 2 * n_ibs0) / (2 * smaller)


def interval_band(kinship: float, se: float) -> RelationBand:
    """The band the whole interval `kinship ± INTERVAL_SES·se` lies in, else `unknown`."""
    low, high = band_of(kinship - INTERVAL_SES * se), band_of(kinship + INTERVAL_SES * se)
    return low if low is high else RelationBand.UNKNOWN


def relatedness(
    a: Loaded,
    b: Loaded,
    panel: SitePanel | None = None,
    *,
    facet_a: str | None = None,
    facet_b: str | None = None,
    min_sites: int = DEFAULT_MIN_SITES,
) -> Relatedness:
    """KING-robust relatedness of one facet of `a` and one of `b` (the module docstring has the method).

    `a` and `b` may be the same `Loaded` (two facets of one file, like F3's two exports) when
    `facet_a` and `facet_b` name them. `panel` defaults to the package's panel for the facets' assembly.
    """
    first, second = facet_of(a, facet_a), facet_of(b, facet_b)
    frame = _pairs([(a, first), (b, second)], panel, min_sites, [(0, 1)])
    return Relatedness(**frame.row(0, named=True))


def pairwise_relatedness(
    loaded: Sequence[Loaded],
    panel: SitePanel | None = None,
    *,
    min_sites: int = DEFAULT_MIN_SITES,
    pairs: Sequence[tuple[int, int]] | None = None,
) -> pl.DataFrame:
    """KING-robust relatedness of every pair of facets across `loaded`, one row per pair.

    Facets are taken in order, each `Loaded`'s facets in its `facets` order; a pair is `(i, j)` with
    `i < j`. `pairs` limits the rows to those index pairs, in the order given. Each facet's calls are
    read once. Columns are `RELATEDNESS_FRAME`.
    """
    facets = [(item, facet) for item in loaded for facet in facets_of(item)]
    if pairs is None:
        pairs = [(i, j) for i in range(len(facets)) for j in range(i + 1, len(facets))]
    return _pairs(facets, panel, min_sites, list(pairs))


def same_person(a: Loaded, b: Loaded, panel: SitePanel | None = None, *, min_sites: int = DEFAULT_MIN_SITES) -> bool | None:
    """Whether the one facet of `a` and the one facet of `b` are one person's: True, False, or None.

    True when the kinship's interval lies wholly above KING's same-person floor (0.354), False when it
    lies wholly below it, None when it crosses it or relatedness is unknown for another reason (no
    shared assembly, no panel, too few shared sites). A monozygotic twin reads True.
    """
    found = relatedness(a, b, panel, min_sites=min_sites)
    if found.kinship is None or found.kinship_se is None:
        return None
    low = found.kinship - INTERVAL_SES * found.kinship_se
    high = found.kinship + INTERVAL_SES * found.kinship_se
    if low > SAME_PERSON_FLOOR:
        return True
    if high <= SAME_PERSON_FLOOR:
        return False
    return None


def _unknown(first: Facet, second: Facet, reason: str, **known: object) -> dict:
    return {"facet_a": first.facet_id, "facet_b": second.facet_id, "method": KING_METHOD, "band": RelationBand.UNKNOWN.value, "reason": reason, **known}


def _pairs(facets: list[tuple[Loaded, Facet]], panel: SitePanel | None, min_sites: int, pairs: list[tuple[int, int]]) -> pl.DataFrame:
    panels: dict[str, SitePanel | None] = {}
    reads: dict[int, pl.DataFrame] = {}

    def panel_for(assembly: str) -> SitePanel | None:
        if assembly not in panels:
            panels[assembly] = panel if panel is not None and panel.assembly_label == assembly else (None if panel is not None else packaged_panel(assembly))
        return panels[assembly]

    def read(index: int, sites: SitePanel) -> pl.DataFrame:
        if index not in reads:
            loaded, facet = facets[index]
            reads[index] = _known_genotypes(loaded, facet, sites)
        return reads[index]

    rows: list[dict] = []
    for i, j in pairs:
        first, second = facets[i][1], facets[j][1]
        if first.assembly_label is None or second.assembly_label is None:
            rows.append(_unknown(first, second, ASSEMBLY_UNKNOWN))
            continue
        if first.assembly_label != second.assembly_label:
            rows.append(_unknown(first, second, ASSEMBLY_MISMATCH))
            continue
        assembly = first.assembly_label
        sites = panel_for(assembly)
        if sites is None:
            rows.append(_unknown(first, second, NO_PANEL, assembly_label=assembly))
            continue
        rows.append(_compare(first, second, read(i, sites), read(j, sites), sites, min_sites))
    return pl.DataFrame(rows, schema=RELATEDNESS_FRAME) if rows else pl.DataFrame(schema=RELATEDNESS_FRAME)


def _known_genotypes(loaded: Loaded, facet: Facet, panel: SitePanel) -> pl.DataFrame:
    """`site, count, vrs_id` at the panel sites where the facet has a genotype."""
    read = read_site_genotypes(loaded, facet, panel.sites).filter(pl.col("count").is_not_null())
    return read.join(_effect_vrs_ids(loaded, facet, panel), on="site", how="left").select("site", "count", "vrs_id")


def _effect_vrs_ids(loaded: Loaded, facet: Facet, panel: SitePanel) -> pl.DataFrame:
    """`site, vrs_id`: the minted VRS id of the panel's ALT on the facet's row at the site, where one is."""
    keys = panel.sites.with_row_index("site").lazy()
    return (
        loaded.calls.filter((pl.col("facet_id") == facet.facet_id) & pl.col("vrs_allele_ids").is_not_null())
        .select("chrom", "pos", "alts", pl.col("vrs_allele_ids").list.slice(1))  # the REF slot comes first
        .filter(pl.col("alts").list.len() == pl.col("vrs_allele_ids").list.len())
        .join(keys, on=["chrom", "pos"], how="inner")
        .select("site", "alts", "vrs_allele_ids", "effect_allele")
        .explode("alts", "vrs_allele_ids", empty_as_null=True)
        .filter((pl.col("alts") == pl.col("effect_allele")) & pl.col("vrs_allele_ids").is_not_null())
        .group_by("site")
        .agg(pl.col("vrs_allele_ids").unique())
        .filter(pl.col("vrs_allele_ids").list.len() == 1)
        .select("site", vrs_id=pl.col("vrs_allele_ids").list.first())
        .collect()
    )


def _compare(first: Facet, second: Facet, a: pl.DataFrame, b: pl.DataFrame, panel: SitePanel, min_sites: int) -> dict:
    blocks = max(1, math.ceil(panel.sites.height / JACKKNIFE_BLOCKS))
    joined = a.join(b, on="site", how="inner", suffix="_b").select(
        (pl.col("site") // blocks).alias("block"),
        het_a=pl.col("count") == 1,
        het_b=pl.col("count_b") == 1,
        het_both=(pl.col("count") == 1) & (pl.col("count_b") == 1),
        ibs0=(pl.col("count").cast(pl.Int16) - pl.col("count_b").cast(pl.Int16)).abs() == 2,
        ibs2=pl.col("count") == pl.col("count_b"),
        hom_ref_a=pl.col("count") == 0,
        hom_ref_b=pl.col("count_b") == 0,
        vrs_compared=pl.col("vrs_id").is_not_null() & pl.col("vrs_id_b").is_not_null(),
        vrs_differ=pl.col("vrs_id") != pl.col("vrs_id_b"),
    )
    per_block = joined.group_by("block").agg(
        n=pl.len(),
        het_a=pl.col("het_a").sum(),
        het_b=pl.col("het_b").sum(),
        het_both=pl.col("het_both").sum(),
        ibs0=pl.col("ibs0").sum(),
        ibs2=pl.col("ibs2").sum(),
        hom_ref_a=pl.col("hom_ref_a").sum(),
        hom_ref_b=pl.col("hom_ref_b").sum(),
        vrs_compared=pl.col("vrs_compared").sum(),
        vrs_differ=pl.col("vrs_differ").fill_null(False).sum(),
    ).sort("block")  # group_by's order varies; the jackknife's float sums must not (charter 12)
    totals = {name: int(per_block[name].sum()) for name in per_block.columns if name != "block"}
    known = {
        "assembly_label": first.assembly_label,
        "panel": panel.source,
        "n_sites": totals["n"],
        "n_het_a": totals["het_a"],
        "n_het_b": totals["het_b"],
        "n_het_both": totals["het_both"],
        "n_ibs0": totals["ibs0"],
        "n_ibs2": totals["ibs2"],
        "n_vrs_compared": totals["vrs_compared"],
        "n_vrs_differ": totals["vrs_differ"],
    }
    if totals["n"] < min_sites:
        return _unknown(first, second, TOO_FEW_SITES, **known)
    kinship = king_robust(totals["het_a"], totals["het_b"], totals["het_both"], totals["ibs0"])
    if kinship is None:
        return _unknown(first, second, NO_HETEROZYGOUS_SITES, **known)
    known |= {"ibs0": totals["ibs0"] / totals["n"], "ibs2": totals["ibs2"] / totals["n"], "kinship": kinship}
    se = _jackknife_se(per_block, totals)
    known["kinship_se"] = se
    band = interval_band(kinship, se) if se is not None else band_of(kinship)
    reason = BETWEEN_BANDS if band is RelationBand.UNKNOWN else None
    if totals["hom_ref_a"] == 0 or totals["hom_ref_b"] == 0:
        band, reason = RelationBand.UNKNOWN, VARIANT_SITES_ONLY
    return {"facet_a": first.facet_id, "facet_b": second.facet_id, "method": KING_METHOD, "band": band.value, "reason": reason, **known}


def _jackknife_se(per_block: pl.DataFrame, totals: dict[str, int]) -> float | None:
    """Delete-one-block jackknife standard error of the kinship; None with fewer than two blocks."""
    if per_block.height < 2:
        return None
    estimates = []
    for row in per_block.iter_rows(named=True):
        left = king_robust(
            totals["het_a"] - row["het_a"], totals["het_b"] - row["het_b"], totals["het_both"] - row["het_both"], totals["ibs0"] - row["ibs0"]
        )
        if left is None:
            return None
        estimates.append(left)
    count = len(estimates)
    mean = sum(estimates) / count
    return math.sqrt((count - 1) / count * sum((value - mean) ** 2 for value in estimates))


__all__ = [
    "ASSEMBLY_MISMATCH",
    "ASSEMBLY_UNKNOWN",
    "BETWEEN_BANDS",
    "DEFAULT_MIN_SITES",
    "NO_HETEROZYGOUS_SITES",
    "NO_PANEL",
    "RELATEDNESS_FRAME",
    "RELATEDNESS_METHODS",
    "RELATEDNESS_REASONS",
    "TOO_FEW_SITES",
    "RelationBand",
    "Relatedness",
    "SitePanel",
    "interval_band",
    "band_of",
    "king_robust",
    "packaged_assemblies",
    "packaged_panel",
    "pairwise_relatedness",
    "relatedness",
    "same_person",
]
