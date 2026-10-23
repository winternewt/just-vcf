"""A method-error atlas: how often two methods disagree, by region of the genome (RM279, charter 8).

`method_atlas(pairs, regions)` compares pairs of facets of one person (`concordance`) and pools the
result per (method of each side, region class, variant class, scope): how many alleles each pair
compared, how often both sides called a genotype and it differed, how often only one side carried
the allele, and the mean `p_same` where the comparison was graded, each rate with a binomial
confidence interval. The output never holds a per-site row, only counts and rates. The rules are in
docs/SCHEMA.md ("A method-error atlas"); in short:

1. **Regions are a `RegionSet`**: named classes of 1-based closed intervals on one assembly, merged
   once when the set is built and reused for every pair. `giab_stratifications` reads the GIAB
   v3.6 stratification BEDs (easy and difficult, low mappability, segmental duplications,
   homopolymers and tandem repeats by length, GC extremes, X and Y). The class `all` is no region:
   every compared allele.
2. **One comparison per pair.** `concordance(a, b)` runs once, and each allele is placed in every
   class whose intervals hold its `pos`, so an indel the two facets spelled at two positions stays
   one paired allele (cutting each facet per class would split it at an edge).
3. **Rates are over carried alleles**: an allele counts when at least one side's genotype holds a
   copy of it. A stated but not carried allele (a DeepVariant `RefCall`, a split row's `0`) says
   nothing about a method's error. Scope `carried_pass` keeps those whose every stating row passed
   its filters (`PASS`, or no filter written).
4. **`p_same` is pooled per basis**, likelihoods and posteriors apart (Decision 3).
5. **Rows by confidence band.** Beside each pooled row (`a_band = b_band = all`), the same counts
   split by both sides' own confidence in the allele (`genotype_evidence`'s `confidence` of the row
   stating it, in `CONFIDENCE_BANDS`), one row per cell of the two bands where alleles fall: how
   often two methods disagree where each was sure, and where it was not.
"""

import logging
import math
from collections.abc import Iterable, Mapping
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.concordance import SITES, AbsenceReading, Agreement, ComparisonBasis, _covered, _merged_intervals, concordance
from just_vcf.evidence import _calls_evidence
from just_vcf.loading import Loaded
from just_vcf.reference import ReferenceStore
from just_vcf.schema import VariantClass, enum_dtype

logger = logging.getLogger(__name__)

ALL = "all"
"""The region class that is no region: every allele the pair compared."""

GIAB_V36_FILES: dict[str, str] = {
    "easy": "{asm}_notinalldifficultregions.bed.gz",
    "difficult": "{asm}_alldifficultregions.bed.gz",
    "low_mappability": "{asm}_lowmappabilityall.bed.gz",
    "segdups": "{asm}_segdups.bed.gz",
    "homopolymer_4to6": "{asm}_SimpleRepeat_homopolymer_4to6_slop5.bed.gz",
    "homopolymer_7to11": "{asm}_SimpleRepeat_homopolymer_7to11_slop5.bed.gz",
    "homopolymer_ge12": "{asm}_SimpleRepeat_homopolymer_ge12_slop5.bed.gz",
    "tandem_repeat_le50bp": "{asm}_AllTandemRepeats_le50bp_slop5.bed.gz",
    "tandem_repeat_51to200bp": "{asm}_AllTandemRepeats_51to200bp_slop5.bed.gz",
    "tandem_repeat_201to10000bp": "{asm}_AllTandemRepeats_201to10000bp_slop5.bed.gz",
    "tandem_repeat_ge10001bp": "{asm}_AllTandemRepeats_ge10001bp_slop5.bed.gz",
    "gc_lt25_or_gt65": "{asm}_gclt25orgt65_slop50.bed.gz",
    "x_par": "{asm}_chrX_PAR.bed.gz",
    "x_nonpar": "{asm}_chrX_nonPAR.bed.gz",
    "y_nonpar": "{asm}_chrY_nonPAR.bed.gz",
}
"""Each GIAB v3.6 region class and its BED's file name (`{asm}` is `GRCh37` or `GRCh38`). The names
are the release's; where the files sit under the directory does not matter."""

REGION_CLASSES: frozenset[str] = frozenset({ALL, *GIAB_V36_FILES})
"""The region classes `giab_stratifications` knows. A `RegionSet` built from frames may name others."""

SCOPES: frozenset[str] = frozenset({"carried", "carried_pass"})
"""Which alleles a row counts: those either side carries, and of those, the ones every stating row
passed its filters."""

CI_Z = 1.959963984540054
"""The normal quantile of the 95% Wilson score interval every rate carries."""

ALL_BANDS = "all"
"""The confidence band that is no band: every allele, whatever either side's confidence."""

NO_CONFIDENCE = "none"
"""The band of an allele whose row states no confidence (a genotype alone, `gt_only`), or where the
side wrote no row."""

BAND_EDGES: tuple[float, ...] = (10.0, 20.0, 30.0, 40.0, 60.0, 100.0)
"""Where the confidence bands (phred-like, `genotype_evidence`'s `confidence`) break, each band
closed below: `[0, 10)`, `[10, 20)`, ... `[100, inf)`."""

CONFIDENCE_BANDS: tuple[str, ...] = ("0-10", "10-20", "20-30", "30-40", "40-60", "60-100", "100+", NO_CONFIDENCE)
"""The bands of one side's confidence in an allele, in order: what `a_band` and `b_band` hold beside
`all`."""


def confidence_band(confidence: pl.Expr) -> pl.Expr:
    """The `CONFIDENCE_BANDS` label of a phred-like confidence; `none` where it is null."""
    labels = CONFIDENCE_BANDS[:-1]
    band = pl.lit(labels[0])
    for edge, label in zip(BAND_EDGES, labels[1:]):
        band = pl.when(confidence >= edge).then(pl.lit(label)).otherwise(band)
    return pl.when(confidence.is_null() | confidence.is_nan()).then(pl.lit(NO_CONFIDENCE)).otherwise(band)

_GIAB_ASSEMBLIES = ("GRCh37", "GRCh38")


class RegionSet(BaseModel):
    """Named classes of intervals on one assembly, merged once (`from_frames`) and reused per pair."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    """What the set is, as the atlas's `region_set` column states it (`GIAB v3.6`)."""
    assembly_label: str
    """The assembly the intervals are on; a pair on another one is refused."""
    intervals: dict[str, pl.DataFrame]
    """Per class, its merged intervals (`chrom`, `_start`, `_end`, 1-based closed), sorted."""

    @classmethod
    def from_frames(cls, name: str, assembly_label: str, frames: Mapping[str, pl.DataFrame | pl.LazyFrame]) -> "RegionSet":
        """A set from one frame per class, each with `chrom`, `pos` and `end_pos`, 1-based and closed as
        `blocks` writes them (a BED's 0-based start is `pos - 1`), any contig spelling."""
        if ALL in frames:
            raise ValueError(f"`{ALL}` is the class with no region; it cannot be given intervals")
        return cls(name=name, assembly_label=assembly_label, intervals={label: _merged_intervals(frame) for label, frame in sorted(frames.items())})

    @property
    def classes(self) -> tuple[str, ...]:
        """`all`, then the set's classes in name order."""
        return (ALL, *sorted(self.intervals))


def read_bed(path: Path | str) -> pl.LazyFrame:
    """A BED file's intervals as `chrom`, `pos`, `end_pos`, 1-based closed (`pos` is the BED start + 1).
    Header, track and browser lines are skipped; columns beyond the third are ignored."""
    # read_csv, not scan_csv: it decompresses a gzipped BED, as the GIAB releases are.
    frame = pl.read_csv(Path(path), separator="\t", has_header=False, comment_prefix="#", infer_schema=False, truncate_ragged_lines=True)
    return (
        frame.lazy()
        .select(pl.nth(0).alias("chrom"), pl.nth(1).alias("start"), pl.nth(2).alias("end"))
        .filter(~pl.col("chrom").str.starts_with("track") & ~pl.col("chrom").str.starts_with("browser"))
        .select("chrom", (pl.col("start").cast(pl.Int64) + 1).alias("pos"), pl.col("end").cast(pl.Int64).alias("end_pos"))
    )


def giab_stratifications(directory: Path | str, assembly_label: str, classes: Iterable[str] | None = None) -> RegionSet:
    """The GIAB v3.6 stratifications for `assembly_label` (`GRCh37` or `GRCh38`) read from the BEDs
    under `directory`, found by file name at any depth (`GIAB_V36_FILES`). `classes` picks some (all
    by default). A missing file is refused with `FileNotFoundError` naming it: a class is never
    silently left out. The files come from
    `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/`."""
    if assembly_label not in _GIAB_ASSEMBLIES:
        raise ValueError(f"GIAB v3.6 stratifications are read here for {', '.join(_GIAB_ASSEMBLIES)}, not {assembly_label}")
    wanted = sorted(set(classes) if classes is not None else set(GIAB_V36_FILES))
    unknown = sorted(set(wanted) - set(GIAB_V36_FILES))
    if unknown:
        raise ValueError(f"unknown GIAB region class(es): {', '.join(unknown)}; known: {', '.join(sorted(GIAB_V36_FILES))}")
    root = Path(directory)
    found: dict[str, Path] = {}
    missing = []
    for label in wanted:
        name = GIAB_V36_FILES[label].format(asm=assembly_label)
        hits = sorted(root.rglob(name))
        if not hits:
            missing.append(name)
        else:
            found[label] = hits[0]
    if missing:
        raise FileNotFoundError(f"{root}: no {', '.join(missing)}")
    return RegionSet.from_frames(f"GIAB v3.6 {assembly_label}", assembly_label, {label: read_bed(path) for label, path in found.items()})


COUNTS = pl.Schema(
    {
        "region_set": pl.String,
        "assembly_label": pl.String,
        "a_caller": pl.String,
        "a_caller_version": pl.String,
        "b_caller": pl.String,
        "b_caller_version": pl.String,
        "region_class": pl.String,
        "variant_class": enum_dtype(VariantClass),
        "scope": pl.String,
        "a_band": pl.String,
        "b_band": pl.String,
        "pairs": pl.UInt32,
        "alleles": pl.UInt64,
        "match": pl.UInt64,
        "mismatch": pl.UInt64,
        "ploidy_differs": pl.UInt64,
        "uncalled": pl.UInt64,
        "unreadable": pl.UInt64,
        "only_a": pl.UInt64,
        "only_b": pl.UInt64,
        "one_sided_reference": pl.UInt64,
        "graded_gl": pl.UInt64,
        "p_same_gl_sum": pl.Float64,
        "graded_gp": pl.UInt64,
        "p_same_gp_sum": pl.Float64,
    }
)
"""Per pair (`pair_counts`) or pooled (`pool_counts`), the counts the atlas's rates come from."""

_KEYS = ["region_set", "assembly_label", "a_caller", "a_caller_version", "b_caller", "b_caller_version", "region_class", "variant_class", "scope", "a_band", "b_band"]
_SUMMED = [name for name in COUNTS.names() if name not in _KEYS]

ATLAS = pl.Schema(
    {
        **{name: COUNTS[name] for name in _KEYS},
        "pairs": pl.UInt32,
        "alleles": pl.UInt64,
        "match": pl.UInt64,
        "mismatch": pl.UInt64,
        "ploidy_differs": pl.UInt64,
        "called_both": pl.UInt64,
        "discordant": pl.UInt64,
        "discordance_rate": pl.Float64,
        "discordance_ci_low": pl.Float64,
        "discordance_ci_high": pl.Float64,
        "only_a": pl.UInt64,
        "only_b": pl.UInt64,
        "one_sided_rate": pl.Float64,
        "one_sided_ci_low": pl.Float64,
        "one_sided_ci_high": pl.Float64,
        "one_sided_reference": pl.UInt64,
        "uncalled": pl.UInt64,
        "unreadable": pl.UInt64,
        "graded_gl": pl.UInt64,
        "p_same_gl_mean": pl.Float64,
        "graded_gp": pl.UInt64,
        "p_same_gp_mean": pl.Float64,
    }
)
"""One row per (region set, assembly, each side's caller and version, region class, variant class,
scope). docs/SCHEMA.md ("A method-error atlas") defines every column."""


def method_atlas(
    pairs: Iterable[tuple[Loaded, Loaded]],
    regions: RegionSet,
    *,
    reference: "ReferenceStore | Path | str | None" = None,
) -> pl.DataFrame:
    """The atlas over `pairs` (each two facets of one person on `regions`' assembly): `pair_counts`
    per pair, pooled by `pool_counts`. A pair `concordance` refuses raises as it does."""
    return pool_counts([pair_counts(a, b, regions, reference=reference) for a, b in pairs])


def pair_counts(
    a: Loaded,
    b: Loaded,
    regions: RegionSet,
    *,
    reference: "ReferenceStore | Path | str | None" = None,
    by_contig: bool = False,
) -> pl.DataFrame:
    """`COUNTS` for one pair: `concordance(a, b)` once, each allele placed in every class of `regions`
    holding its `pos`. The two sides are put in method order (caller, then version), `only_a` and
    `only_b` swapped with them, so a pair counts the same whichever facet is given first.

    `by_contig` runs `concordance` once per contig either facet has rows on and sums the counts: the
    same counts (concordance pairs alleles within a position), with one contig's sites in memory at a
    time instead of a whole genome's."""
    fa = _facet(a)
    fb = _facet(b)
    for facet in (fa, fb):
        if facet["assembly_label"] != regions.assembly_label:
            raise ValueError(f"facet {facet['facet_id']} is on {facet['assembly_label']}, the regions on {regions.assembly_label}; lift one first")
    if by_contig:
        contigs = sorted(set(_contigs(a)) | set(_contigs(b)))
        chunks = [_pair_chunk(_on(a, chrom), _on(b, chrom), regions, reference) for chrom in contigs]
        keys = ["region_class", "variant_class", "scope", "a_band", "b_band"]
        counts = pl.concat(chunks).group_by(keys).agg(pl.exclude(keys).sum()) if chunks else _counts(_NO_SITES, regions)
    else:
        counts = _pair_chunk(a, b, regions, reference)
    swap = _method(fb) < _method(fa)
    if swap:
        fa, fb = fb, fa
    counts = counts.with_columns(
        pl.lit(regions.name).alias("region_set"),
        pl.lit(regions.assembly_label).alias("assembly_label"),
        pl.lit(fa["caller"], dtype=pl.String).alias("a_caller"),
        pl.lit(fa["caller_version"], dtype=pl.String).alias("a_caller_version"),
        pl.lit(fb["caller"], dtype=pl.String).alias("b_caller"),
        pl.lit(fb["caller_version"], dtype=pl.String).alias("b_caller_version"),
        pl.lit(1, dtype=pl.UInt32).alias("pairs"),
    )
    if swap:
        counts = counts.rename({"only_a": "only_b", "only_b": "only_a", "a_band": "b_band", "b_band": "a_band"})
    return _sorted(counts.select(COUNTS.names()).cast(dict(COUNTS)))


def _pair_chunk(a: Loaded, b: Loaded, regions: RegionSet, reference: "ReferenceStore | Path | str | None") -> pl.DataFrame:
    """One comparison's counts, each allele banded by each side's confidence in it."""
    sites = concordance(a, b, reference=reference).sites.collect()
    return _counts(sites, regions, (_confidence_at(a), _confidence_at(b)))


def _confidence_at(loaded: Loaded) -> pl.DataFrame:
    """Per `(chrom, pos)` the facet wrote a row at, the lowest `genotype_evidence` confidence among its
    rows there (null where none states one): what an allele the facet wrote there is banded by."""
    graded = _calls_evidence(loaded.calls, False, loaded.facets).select("chrom", "pos", "confidence")
    return graded.group_by("chrom", "pos").agg(pl.col("confidence").min()).collect()


def pool_counts(counts: Iterable[pl.DataFrame]) -> pl.DataFrame:
    """The `ATLAS` from `COUNTS` frames (one per pair, or already pooled ones): counts summed per key,
    rates and their Wilson intervals computed from the sums. A rate whose denominator is zero is null,
    never zero."""
    frames = list(counts)
    if not frames:
        return pl.DataFrame(schema=ATLAS)
    # Counts written before the bands carry no band columns: they are the pooled rows.
    frames = [frame.with_columns(*[pl.lit(ALL_BANDS).alias(name) for name in ("a_band", "b_band") if name not in frame.columns]) for frame in frames]
    pooled = pl.concat([frame.select(COUNTS.names()).cast(dict(COUNTS)) for frame in frames]).group_by(_KEYS).agg(pl.col(_SUMMED).sum())
    called = pl.col("match") + pl.col("mismatch")
    discordant = pl.col("mismatch")
    one_sided = pl.col("only_a") + pl.col("only_b")
    atlas = pooled.with_columns(
        called.alias("called_both"),
        discordant.alias("discordant"),
        *_wilson(discordant, called, "discordance"),
        *_wilson(one_sided, pl.col("alleles"), "one_sided"),
        pl.when(pl.col("graded_gl") > 0).then(pl.col("p_same_gl_sum") / pl.col("graded_gl")).alias("p_same_gl_mean"),
        pl.when(pl.col("graded_gp") > 0).then(pl.col("p_same_gp_sum") / pl.col("graded_gp")).alias("p_same_gp_mean"),
    )
    return _sorted(atlas.select(ATLAS.names()).cast(dict(ATLAS)))


def _wilson(successes: pl.Expr, trials: pl.Expr, name: str) -> list[pl.Expr]:
    """`{name}_rate` and its 95% Wilson score interval; all three null when `trials` is zero."""
    n = trials.cast(pl.Float64)
    p = successes.cast(pl.Float64) / n
    z2 = CI_Z * CI_Z
    centre = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = CI_Z * ((p * (1 - p) / n + z2 / (4 * n * n)).sqrt()) / (1 + z2 / n)
    some = trials > 0
    return [
        pl.when(some).then(p).alias(f"{name}_rate"),
        pl.when(some).then((centre - half).clip(lower_bound=0.0)).alias(f"{name}_ci_low"),
        pl.when(some).then((centre + half).clip(upper_bound=1.0)).alias(f"{name}_ci_high"),
    ]


def wilson_interval(successes: int, trials: int) -> tuple[float, float] | None:
    """The 95% Wilson score interval of `successes` out of `trials`; None when `trials` is zero."""
    if trials == 0:
        return None
    p = successes / trials
    z2 = CI_Z * CI_Z
    centre = (p + z2 / (2 * trials)) / (1 + z2 / trials)
    half = CI_Z * math.sqrt(p * (1 - p) / trials + z2 / (4 * trials * trials)) / (1 + z2 / trials)
    return max(0.0, centre - half), min(1.0, centre + half)


def _facet(loaded: Loaded) -> dict:
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"an atlas pair takes one facet per side; this Loaded holds {facets.height}")
    return facets.row(0, named=True)


_NO_SITES = pl.DataFrame(schema=SITES)


def _contigs(loaded: Loaded) -> list[str]:
    return loaded.calls.select(pl.col("chrom").unique()).collect()["chrom"].to_list()


def _on(loaded: Loaded, chrom: str) -> Loaded:
    """The facet's rows and blocks on one contig; its `facets` row as it is."""
    return Loaded(facets=loaded.facets, calls=loaded.calls.filter(pl.col("chrom") == chrom), blocks=loaded.blocks.filter(pl.col("chrom") == chrom))


def _method(facet: Mapping) -> tuple[str, str]:
    """The sort key of a side's method; unknown sorts first."""
    return (facet["caller"] or "", facet["caller_version"] or "")


def _passes(side: str) -> pl.Expr:
    """Whether `side`'s row passed its filters: no row there, no filter written, or `PASS` only."""
    filters = pl.col(f"{side}_filters")
    return filters.is_null() | filters.list.eval(~pl.element().is_in(["PASS", ""])).list.any().not_()


def _counts(sites: pl.DataFrame, regions: RegionSet, confidence: tuple[pl.DataFrame, pl.DataFrame] | None = None) -> pl.DataFrame:
    """Counts per (region class, variant class, scope), pooled (`a_band = b_band = all`) and, where both
    sides' confidence at their rows is given (`_confidence_at`), per cell of the two sides' bands."""
    # A genotype carries the allele where it spells it, or where it holds it under another ALT of its
    # row (a complex row met by VRS id, whose copies `concordance` read from the change).
    spelled = pl.col("a_gt_alleles").list.contains(pl.col("alt")).fill_null(False) | pl.col("b_gt_alleles").list.contains(pl.col("alt")).fill_null(False)
    carries = spelled | (pl.col("a_copies") > 0).fill_null(False) | (pl.col("b_copies") > 0).fill_null(False)
    probe = sites.filter(carries).with_columns(pl.col("pos").cast(pl.Int64).alias("_at"), (_passes("a") & _passes("b")).alias("_pass"))
    if confidence is not None:
        for side, frame in zip("ab", confidence):
            at = frame.select("chrom", pl.col("pos").cast(pl.UInt32).alias(f"{side}_pos"), pl.col("confidence").alias(f"_{side}_confidence"))
            probe = probe.join(at, on=["chrom", f"{side}_pos"], how="left", maintain_order="left").with_columns(
                confidence_band(pl.col(f"_{side}_confidence")).alias(f"{side}_band")
            )
    probe = probe.sort("chrom", "_at")
    members = {ALL: pl.Series("_covered", [True] * probe.height, dtype=pl.Boolean)}
    for label, merged in regions.intervals.items():
        members[label] = _covered(probe, merged)
    pieces = []
    for label in regions.classes:
        inside = probe.filter(members[label])
        for scope, kept in (("carried", inside), ("carried_pass", inside.filter(pl.col("_pass")))):
            tags = [pl.lit(label).alias("region_class"), pl.lit(scope).alias("scope")]
            pieces.append(_tally(kept).with_columns(*tags, pl.lit(ALL_BANDS).alias("a_band"), pl.lit(ALL_BANDS).alias("b_band")))
            if confidence is not None:
                pieces.append(_tally(kept, by=["a_band", "b_band"]).with_columns(*tags))
    columns = pieces[0].columns
    return pl.concat([piece.select(columns) for piece in pieces], how="vertical")


def _tally(sites: pl.DataFrame, by: list[str] | None = None) -> pl.DataFrame:
    """Counts per variant class, one row for every class even where it has no allele; with `by`, per
    variant class and value of those columns, only where alleles fall."""
    agreement = pl.col("agreement").cast(pl.String)
    gl = (pl.col("basis").cast(pl.String) == ComparisonBasis.GL_FLAT_PRIOR.value) & pl.col("p_same").is_not_null()
    gp = (pl.col("basis").cast(pl.String) == ComparisonBasis.GP.value) & pl.col("p_same").is_not_null()
    one_sided_reference = pl.col("absence").cast(pl.String).is_in([AbsenceReading.REFERENCE.value, AbsenceReading.REF_BLOCK.value])
    tallied = sites.group_by(["variant_class", *(by or [])]).agg(
        pl.len().cast(pl.UInt64).alias("alleles"),
        *[_count(agreement == value.value).alias(value.value) for value in Agreement],
        _count(pl.col("ploidy_differs") & agreement.is_in([Agreement.MATCH.value, Agreement.MISMATCH.value])).alias("ploidy_differs"),
        _count(one_sided_reference & agreement.is_in([Agreement.ONLY_A.value, Agreement.ONLY_B.value])).alias("one_sided_reference"),
        _count(gl).alias("graded_gl"),
        pl.col("p_same").filter(gl).sum().cast(pl.Float64).alias("p_same_gl_sum"),
        _count(gp).alias("graded_gp"),
        pl.col("p_same").filter(gp).sum().cast(pl.Float64).alias("p_same_gp_sum"),
    )
    if by is not None:
        small = tallied.filter(pl.col("variant_class").cast(pl.String) != VariantClass.SV.value)
        return small.with_columns(pl.col("variant_class").cast(enum_dtype(VariantClass)))
    every = pl.DataFrame({"variant_class": [cls.value for cls in VariantClass if cls is not VariantClass.SV]}, schema={"variant_class": enum_dtype(VariantClass)})
    return every.join(tallied.with_columns(pl.col("variant_class").cast(enum_dtype(VariantClass))), on="variant_class", how="left").fill_null(0)


def _count(condition: pl.Expr) -> pl.Expr:
    return condition.fill_null(False).sum().cast(pl.UInt64)


def _sorted(frame: pl.DataFrame) -> pl.DataFrame:
    order = {label: rank for rank, label in enumerate([ALL, *GIAB_V36_FILES])}
    rank = pl.col("region_class").replace_strict(order, default=len(order), return_dtype=pl.UInt16)
    bands = {label: rank for rank, label in enumerate([ALL_BANDS, *CONFIDENCE_BANDS])}
    band = [pl.col(name).replace_strict(bands, default=len(bands), return_dtype=pl.UInt16) for name in ("a_band", "b_band")]
    return frame.sort(
        [pl.col(name) for name in _KEYS[:6]] + [rank, pl.col("region_class"), pl.col("variant_class"), pl.col("scope"), *band],
        nulls_last=True,
    )
