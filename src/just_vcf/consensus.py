"""A consensus posterior across one person's facets (RM280; charter invariants 5, 7, 8, 9).

`consensus(person, atlas=...)` reads every facet of one person at the union of the alleles they state
on one assembly, and gives each allele a posterior over the copies of it the person carries, combining
each facet's evidence as independent and weighting it by that facet's method error rate from a
method-error atlas (RM279). The result is one more facet, derived from all of them (`relation =
consensus`): the parents keep their own calls, and each consensus row says how many facets supported
its most probable genotype, how many contradicted it, and how many were silent (charter 8). The rules
are in docs/SCHEMA.md ("A consensus posterior"); in short:

1. **The unit is one ALT allele**, as in `concordance` and the atlas: a row is `(chrom, pos, REF,
   ALT)`, and its vector is over the copies of that ALT (0, 1, 2 diploid; 0, 1 haploid), which is VCF
   genotype order over `[REF, ALT]`. A multi-allelic site gives one row per ALT, each a marginal.
2. **What each facet says** about an allele: a row stating it (its vector where it has one, else its
   genotype); one row at the position over the same REF span without it (zero copies, `site`); a
   `ref` block over it (zero copies, graded by the block's likelihoods where it has them); its silence
   where `absence_is_reference` says it reads as reference; otherwise nothing.
3. **An error rate turns a statement into a likelihood**: `L(k) = (1 - e)·q(k) + e·(1 - q(k))/P`,
   where `q` is the facet's own distribution over copies (a genotype is a point mass), `P` the ploidy.
   A statement gets the method's genotype error, a block or a silence its silence error
   (`method_error_rates`). A method the atlas never measured gets a stated default, and the row says
   how many of its facets did.
4. **Posterior** `∝ prior(k) · ∏ L_f(k)`, computed in log space. The prior is flat over the copies.
"""

import hashlib
import json
import logging
import math
from collections import Counter
from collections.abc import Iterable
from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.atlas import ALL, CI_Z, RegionSet
from just_vcf.callset import NOT_A_VARIANT, absence_is_reference, variant_class_expr
from just_vcf.concordance import _covered, _declared_contigs, _distribution
from just_vcf.evidence import GradeBasis, _calls_evidence
from just_vcf.links import CONSENSUS, LINEAGE_RELATIONS, canonical_links, links_frame
from just_vcf.loading import Loaded
from just_vcf.person import Person, SubjectError, _pieces
from just_vcf.reference import ReferenceMismatchError
from just_vcf.schema import BLOCKS, CALLS, CONTIG_RANK, FACETS, LINKS, SCHEMA_VERSION, AbsenceMeans, EvidenceBasis, RefCheck, VariantClass, enum_dtype

logger = logging.getLogger(__name__)

__all__ = [
    "CONSENSUS_SUPPORT",
    "CONSENSUS_VERSION",
    "FACET_RATES",
    "METHOD_RATES",
    "Consensus",
    "GpReading",
    "RateBasis",
    "consensus",
    "method_error_rates",
]

CONSENSUS_VERSION = 1
"""Part of the consensus `facet_id` and its links' `step_version`: bumped whenever the rule changes
what it emits."""

CALLER = "just-vcf consensus"
"""The consensus facet's `caller`: the rule that made its rows, versioned by `CONSENSUS_VERSION`."""

_SMALL = (VariantClass.SNV.value, VariantClass.MNV.value, VariantClass.INDEL.value)
_LN10 = math.log(10.0)


class RateBasis(str, Enum):
    """Where a facet's error rate for one region class and variant class came from."""

    METHOD = "method"
    """Atlas rows of the facet's caller at its version."""
    CALLER = "caller"
    """Atlas rows of the facet's caller at other versions, pooled (the facet's version has none)."""
    DEFAULT = "default"
    """No atlas row for the caller: the stated default (`consensus(default_rate=...)`, else the worst
    rate any measured method has in that class)."""


class GpReading(str, Enum):
    """How a row whose evidence is a posterior (GP) is read (docs/SCHEMA.md Decision 3)."""

    CALL = "call"
    """As the genotype it calls, weighted by the method's error rate: the posterior's own prior is not
    multiplied in. The default."""
    LIKELIHOOD = "likelihood"
    """Its vector multiplied in as if it were a likelihood. The source's prior (a population panel)
    then counts once per such facet; the support frame counts how many facets were read so."""


METHOD_RATES = pl.Schema(
    {
        "caller": pl.String,
        "caller_version": pl.String,
        "region_class": pl.String,
        "variant_class": enum_dtype(VariantClass),
        "genotype_error": pl.Float64,
        "genotype_attribution": pl.String,
        "silence_error": pl.Float64,
        "partners": pl.UInt32,
        "called_both": pl.UInt64,
    }
)
"""What `method_error_rates` returns: one method's error rates per region class and variant class.
`caller_version` is null for rates pooled over a caller's versions (`by="caller"`)."""

FACET_RATES = pl.Schema(
    {
        "facet_id": pl.String,
        "caller": pl.String,
        "caller_version": pl.String,
        "region_class": pl.String,
        "variant_class": enum_dtype(VariantClass),
        "genotype_error": pl.Float64,
        "silence_error": pl.Float64,
        "rate_basis": enum_dtype(RateBasis),
    }
)
"""`Consensus.rates`: the rates each parent facet was weighted by, per region class and variant class."""

CONSENSUS_SUPPORT = pl.Schema(
    {
        "chrom": pl.String,
        "pos": pl.UInt32,
        "ref": pl.String,
        "alt": pl.String,
        "variant_class": enum_dtype(VariantClass),
        "ploidy": pl.UInt8,
        "copies": pl.UInt8,
        "confidence": pl.Float64,
        "facets": pl.UInt32,
        "supported": pl.UInt32,
        "contradicted": pl.UInt32,
        "silent": pl.UInt32,
        "defaulted": pl.UInt32,
        "gp_as_likelihood": pl.UInt32,
    }
)
"""`Consensus.support`: one row per consensus row, in its order. `copies` is the most probable number
of copies of `alt`, `confidence` `-10·log10(1 - p)` of it from the log-space posterior (the `gp`
column is Float32 and can round a tiny remainder to zero). `supported` + `contradicted` + `silent` =
`facets`, the parents. `defaulted` counts the facets that spoke with a default rate,
`gp_as_likelihood` those whose posterior was multiplied in (`GpReading.LIKELIHOOD`)."""


class Consensus(Loaded):
    """What `consensus` returns: the derived facet as a `Loaded`, plus how each row came about."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    support: pl.LazyFrame
    """`CONSENSUS_SUPPORT`: per row, the copies it gives, its confidence and who agreed."""
    rates: pl.DataFrame
    """`FACET_RATES`: the rates each parent was weighted by."""


# --------------------------------------------------------------------------------------------------
# Rates


def _floor(n: pl.Expr) -> pl.Expr:
    """The 95% Wilson upper bound of zero events in `n`: the least rate `n` trials can rule out."""
    z2 = CI_Z * CI_Z
    return pl.lit(z2) / (n.cast(pl.Float64) + z2)


def method_error_rates(atlas: pl.DataFrame, *, scope: str = "carried_pass", by: str = "method") -> pl.DataFrame:
    """Per method, its error rates by region class and variant class, from an atlas's pair rates.

    The atlas states how often two methods disagree, not which one erred. With independent rare
    errors a pair's discordance is about the sum of the two methods' error rates, so:

    - **`genotype_error`**, the rate at which a method's stated genotype is wrong: where the method
      forms a triangle with two partners that were also compared with each other, `(d_xy + d_xz -
      d_yz) / 2`, averaged over its triangles and kept at or below its lowest pair discordance
      (`genotype_attribution = triangle`); otherwise half of each pair's discordance, averaged over
      its partners (`half_pair`);
    - **`silence_error`**, the rate at which a method is silent where a variant is: the share of the
      alleles a partner carries that the method left one-sided, the lowest over its partners. It is an
      upper bound: a partner's false calls count against the method too.

    No rate goes below the 95% Wilson upper bound of zero events in the alleles it rests on, so a
    method that never disagreed in a million alleles gets about 4e-6, never zero. Rows of `scope` only;
    rows of an unknown caller are left out (an unknown method matches no facet). `by="caller"` pools
    each caller's versions into one method and drops the pairs of one caller with itself, whose
    disagreement is version drift. Assemblies are pooled. Returns `METHOD_RATES`.
    """
    if by not in ("method", "caller"):
        raise ValueError(f"by is 'method' or 'caller', not {by!r}")
    rows = atlas.filter((pl.col("scope") == scope) & pl.col("a_caller").is_not_null() & pl.col("b_caller").is_not_null())
    if by == "caller":
        rows = rows.filter(pl.col("a_caller") != pl.col("b_caller")).with_columns(
            pl.lit(None, dtype=pl.String).alias("a_caller_version"), pl.lit(None, dtype=pl.String).alias("b_caller_version")
        )
    if rows.is_empty():
        return pl.DataFrame(schema=METHOD_RATES)
    sep, unknown = "\x1f", "\x00"
    key_a = pl.concat_str([pl.col("a_caller"), pl.col("a_caller_version").fill_null(unknown)], separator=sep)
    key_b = pl.concat_str([pl.col("b_caller"), pl.col("b_caller_version").fill_null(unknown)], separator=sep)
    pooled = (
        rows.with_columns(key_a.alias("_a"), key_b.alias("_b"), pl.col("variant_class").cast(pl.String))
        .group_by("_a", "_b", "region_class", "variant_class")
        .agg(pl.col("called_both", "discordant", "alleles", "only_a", "only_b").sum())
        .with_columns(pl.when(pl.col("called_both") > 0).then(pl.col("discordant") / pl.col("called_both")).alias("_d"))
    )
    cls = ["region_class", "variant_class"]
    forward = pooled.select(
        pl.col("_a").alias("_x"), pl.col("_b").alias("_y"), *cls, "_d", pl.col("called_both").alias("_n"),
        pl.col("only_b").alias("_silent"), (pl.col("alleles") - pl.col("only_a")).alias("_carried"),
    )  # fmt: skip
    backward = pooled.select(
        pl.col("_b").alias("_x"), pl.col("_a").alias("_y"), *cls, "_d", pl.col("called_both").alias("_n"),
        pl.col("only_a").alias("_silent"), (pl.col("alleles") - pl.col("only_b")).alias("_carried"),
    )  # fmt: skip
    edges = pl.concat([forward, backward]).filter(pl.col("_x") != pl.col("_y"))
    measured = edges.filter(pl.col("_d").is_not_null())
    halves = measured.group_by("_x", *cls).agg(
        pl.max_horizontal(pl.col("_d") / 2, _floor(pl.col("_n"))).mean().alias("_half"),
        pl.col("_d").min().alias("_bound"),
        pl.len().cast(pl.UInt32).alias("partners"),
        pl.col("_n").min().cast(pl.UInt64).alias("called_both"),
    )
    second = measured.select(pl.col("_x"), pl.col("_y").alias("_z"), *cls, pl.col("_d").alias("_d_xz"), pl.col("_n").alias("_n_xz"))
    opposite = measured.select(pl.col("_x").alias("_y"), pl.col("_y").alias("_z"), *cls, pl.col("_d").alias("_d_yz"), pl.col("_n").alias("_n_yz"))
    triangles = (
        measured.join(second, on=["_x", *cls])
        .filter(pl.col("_y") < pl.col("_z"))
        .join(opposite, on=["_y", "_z", *cls])
        .with_columns(
            pl.max_horizontal(
                (pl.col("_d") + pl.col("_d_xz") - pl.col("_d_yz")) / 2, _floor(pl.min_horizontal("_n", "_n_xz", "_n_yz"))
            ).alias("_t")
        )
        .group_by("_x", *cls)
        .agg(pl.col("_t").mean().alias("_triangle"))
    )
    silence = (
        edges.filter(pl.col("_carried") > 0)
        .group_by("_x", *cls)
        .agg(pl.max_horizontal(pl.col("_silent") / pl.col("_carried"), _floor(pl.col("_carried"))).min().alias("silence_error"))
    )
    methods = edges.select("_x", *cls).unique()
    rates = (
        methods.join(halves, on=["_x", *cls], how="left")
        .join(triangles, on=["_x", *cls], how="left")
        .join(silence, on=["_x", *cls], how="left")
        .with_columns(
            pl.when(pl.col("_triangle").is_not_null()).then(pl.min_horizontal("_triangle", "_bound")).otherwise(pl.col("_half")).alias("genotype_error"),
            pl.when(pl.col("_triangle").is_not_null()).then(pl.lit("triangle")).when(pl.col("_half").is_not_null()).then(pl.lit("half_pair")).alias("genotype_attribution"),
            pl.col("_x").str.split_exact(sep, 1).struct.rename_fields(["caller", "_version"]).alias("_m"),
        )
        .unnest("_m")
        .with_columns(pl.when(pl.col("_version") != unknown).then(pl.col("_version")).alias("caller_version"))
        .with_columns(pl.col("partners").fill_null(0), pl.col("called_both").fill_null(0))
    )
    order = {label: rank for rank, label in enumerate([ALL, *sorted(set(rates["region_class"].to_list()) - {ALL})])}
    return (
        rates.select(METHOD_RATES.names())
        .cast(dict(METHOD_RATES))
        .sort(["caller", "caller_version", pl.col("region_class").replace_strict(order, return_dtype=pl.UInt16), "variant_class"], nulls_last=True)
    )


def _facet_rates(facets: list[dict], atlas: pl.DataFrame | None, classes: tuple[str, ...], default_rate: float | None, scope: str) -> pl.DataFrame:
    """`FACET_RATES` for every parent: its method's rates, else its caller's, else the default."""
    grid = pl.DataFrame(
        [{"region_class": c, "variant_class": v} for c in classes for v in _SMALL],
        schema={"region_class": pl.String, "variant_class": enum_dtype(VariantClass)},
    )
    empty = pl.DataFrame(schema=METHOD_RATES)
    by_method = method_error_rates(atlas, scope=scope) if atlas is not None else empty
    by_caller = method_error_rates(atlas, scope=scope, by="caller") if atlas is not None else empty
    if default_rate is not None:
        if not 0.0 < default_rate < 1.0:
            raise ValueError(f"default_rate must lie strictly between 0 and 1, not {default_rate}")
        default = grid.with_columns(pl.lit(default_rate).alias("_dg"), pl.lit(default_rate).alias("_ds"))
    else:
        worst = by_method.group_by("region_class", "variant_class").agg(
            pl.col("genotype_error").max().alias("_dg"), pl.col("silence_error").max().alias("_ds")
        )
        default = grid.join(worst, on=["region_class", "variant_class"], how="left")
    out = []
    for facet in facets:
        caller, version = facet["caller"], facet["caller_version"]
        own = by_method.filter((pl.col("caller") == caller) & pl.col("caller_version").eq_missing(pl.lit(version, dtype=pl.String))) if caller else empty
        kin = by_caller.filter(pl.col("caller") == caller) if caller else empty
        rows = (
            grid.join(own.select("region_class", "variant_class", pl.col("genotype_error").alias("_mg"), pl.col("silence_error").alias("_ms")), on=["region_class", "variant_class"], how="left")
            .join(kin.select("region_class", "variant_class", pl.col("genotype_error").alias("_cg"), pl.col("silence_error").alias("_cs")), on=["region_class", "variant_class"], how="left")
            .join(default, on=["region_class", "variant_class"], how="left")
        )
        measured_m = pl.col("_mg").is_not_null() & pl.col("_ms").is_not_null()
        measured_c = pl.col("_cg").is_not_null() & pl.col("_cs").is_not_null()
        rows = rows.with_columns(
            pl.when(measured_m).then(pl.col("_mg")).when(measured_c).then(pl.col("_cg")).otherwise(pl.col("_dg")).alias("genotype_error"),
            pl.when(measured_m).then(pl.col("_ms")).when(measured_c).then(pl.col("_cs")).otherwise(pl.col("_ds")).alias("silence_error"),
            pl.when(measured_m).then(pl.lit(RateBasis.METHOD.value)).when(measured_c).then(pl.lit(RateBasis.CALLER.value)).otherwise(pl.lit(RateBasis.DEFAULT.value)).alias("rate_basis"),
            pl.lit(facet["facet_id"]).alias("facet_id"),
            pl.lit(caller, dtype=pl.String).alias("caller"),
            pl.lit(version, dtype=pl.String).alias("caller_version"),
        )
        missing = rows.filter(pl.col("genotype_error").is_null() | pl.col("silence_error").is_null())
        if missing.height:
            where = ", ".join(sorted({f"{r['region_class']}/{r['variant_class']}" for r in missing.iter_rows(named=True)}))
            raise ValueError(
                f"facet {facet['facet_id']} has no atlas rate in {where} and no measured method gives a default there; pass default_rate=..."
            )
        out.append(rows.select(FACET_RATES.names()).cast(dict(FACET_RATES)))
    return pl.concat(out)


# --------------------------------------------------------------------------------------------------
# Inputs


def _parents(person: Person | Loaded | Iterable[Loaded]) -> tuple[list[Loaded], pl.LazyFrame, list[str]]:
    """One `Loaded` per parent facet, sorted by `facet_id`, the links they carried, and warnings."""
    warnings: list[str] = []
    if isinstance(person, Person):
        if person.one_person is False:
            raise SubjectError(f"the relatedness check read the facets of subject {person.subject!r} as more than one person; a consensus needs one")
        if person.one_person is None:
            warnings.append("the facets were not all checked as one person (a pair not run or not told); the consensus takes them as given")
        whole = Loaded(facets=person.facets, calls=person.calls, blocks=person.blocks, links=person.links, schema_version=person.schema_version)
        pieces = _pieces((whole,))
    else:
        given = (person,) if isinstance(person, Loaded) else tuple(person)
        pieces = _pieces(given)
    rows = [piece.facets.collect().row(0, named=True) for piece in pieces]
    ids = [row["facet_id"] for row in rows]
    twice = sorted({i for i in ids if ids.count(i) > 1})
    if twice:
        raise ValueError(f"facet(s) given more than once: {', '.join(twice)}")
    if len(pieces) < 2:
        raise ValueError(f"a consensus needs two facets or more; {len(pieces)} given")
    subjects = sorted({row["subject"] for row in rows if row["subject"] is not None})
    if len(subjects) > 1:
        raise SubjectError(f"facets under {len(subjects)} different subjects cannot make one consensus")
    already = sorted(row["facet_id"] for row in rows if row["relation"] == CONSENSUS)
    if already:
        raise ValueError(f"facet(s) {', '.join(already)} are consensus facets; their parents' evidence would count twice")
    order = sorted(range(len(pieces)), key=lambda i: ids[i])
    pieces = [pieces[i] for i in order]
    links = canonical_links(pl.concat([piece.links.select(LINKS.names()) for piece in pieces], how="vertical"))
    _refuse_lineage([rows[i] for i in order], links.collect())
    return pieces, links, warnings


def _refuse_lineage(rows: list[dict], links: pl.DataFrame) -> None:
    """Refuse a parent made from another parent: the same reads would count twice."""
    parents: dict[str, set[str]] = {}
    for link in links.filter(pl.col("relation").is_in(sorted(LINEAGE_RELATIONS))).iter_rows(named=True):
        parents.setdefault(link["to_facet_id"], set()).add(link["from_facet_id"])
    for row in rows:
        if row["derived_from"] is not None:
            parents.setdefault(row["facet_id"], set()).add(row["derived_from"])
    given = {row["facet_id"] for row in rows}
    for row in rows:
        seen: set[str] = set()
        todo = list(parents.get(row["facet_id"], ()))
        while todo:
            ancestor = todo.pop()
            if ancestor in seen:
                continue
            seen.add(ancestor)
            todo += list(parents.get(ancestor, ()))
        hit = sorted(seen & given)
        if hit:
            raise ValueError(f"facet {row['facet_id']} was made from facet(s) {', '.join(hit)}, also given; pass one of them, not both")


def _assembly(rows: list[dict], warnings: list[str]) -> tuple[str | None, str | None]:
    """The parents' one assembly (`assembly_label`, `ref_seqcol_digest`), or a refusal."""
    digests = {row["ref_seqcol_digest"] for row in rows}
    labels = {row["assembly_label"] for row in rows}
    if None not in digests:
        if len(digests) > 1:
            raise ReferenceMismatchError(f"the facets' REF was read from {len(digests)} collections; lift them onto one assembly first")
        label = labels.pop() if len(labels) == 1 else None
        return label, digests.pop()
    if len(labels) == 1 and None not in labels:
        label = labels.pop()
        warnings.append(f"assembly not checked by content: REF not read on every facet, all labelled {label}")
        digest = digests - {None}
        if len(digest) > 1:
            raise ReferenceMismatchError(f"the facets' REF was read from {len(digest)} collections; lift them onto one assembly first")
        return label, None
    if len(labels - {None}) > 1:
        raise ReferenceMismatchError(f"the facets are on {', '.join(sorted(labels - {None}))}; lift them onto one assembly first")
    raise ValueError("cannot tell the facets are on one assembly; run resolve_ref on them first")


# --------------------------------------------------------------------------------------------------
# Readings


_SILENT = ("filtered", "no_call", "covered", "several_rows", "other_span", "no_call_block", "unknown", "ploidy_differs")
"""Why a facet said nothing about an allele (the counts in the warnings)."""


def _passes(filters: pl.Expr) -> pl.Expr:
    return filters.is_null() | filters.list.eval(~pl.element().is_in(["PASS", ""])).list.any().not_()


def _evidence(piece: Loaded, chrom: str, gq_on_reference: bool) -> pl.DataFrame:
    """The facet's rows on `chrom`, graded as `genotype_evidence` grades them, with what matching needs."""
    calls = piece.calls.filter(pl.col("chrom") == chrom).with_row_index("_row").with_columns(
        pl.coalesce("ref", "source_ref").str.to_uppercase().alias("_ref_key"), pl.col("ref").is_not_null().alias("_ref_read")
    )
    keep = ("_row", "_ref_key", "_ref_read", "alts", "gt_alleles", "filters")
    return _calls_evidence(calls, gq_on_reference, piece.facets, keep=keep).collect()


def _stated_alleles(evidence: pl.DataFrame) -> pl.DataFrame:
    """One row per sequence ALT allele a row states, with the row's evidence."""
    exploded = (
        evidence.filter(pl.col("alts").list.len() > 0)
        .explode("alts", empty_as_null=False)
        .rename({"alts": "alt"})
        .filter(pl.col("_ref_key").is_not_null() & ~pl.col("alt").is_in(sorted(NOT_A_VARIANT)))
        .with_columns(variant_class_expr(pl.col("_ref_key"), pl.col("alt")).alias("variant_class"))
    )
    return exploded.filter(pl.col("variant_class").is_in(list(_SMALL)))


def _copies_of(gt_alleles: str, alt: str) -> pl.Expr:
    """Copies of the row's `alt` among the haplotypes of `gt_alleles` (two slots at most are counted)."""
    slots = [(pl.col(gt_alleles).list.get(i, null_on_oob=True) == pl.col(alt)).fill_null(False).cast(pl.Int8) for i in range(2)]
    return pl.sum_horizontal(slots).cast(pl.Int8)


def _fully_called(gt_alleles: str) -> pl.Expr:
    col = pl.col(gt_alleles)
    return (col.is_not_null() & col.list.len().is_between(1, 2) & ~col.list.eval(pl.element().is_null()).list.any()).fill_null(False)


def _vector_copies(frame: pl.DataFrame) -> pl.DataFrame:
    """Per `_site`, the distribution over copies of `alt` from `genotype_probabilities` over
    `genotype_alleles` (VCF genotype order, haploid or diploid), as `_q0`, `_q1`, `_q2` and `_ploidy`."""
    work = frame.select(
        pl.col("_site").alias("_pair"),
        pl.col("genotype_alleles").alias("x_source_alleles"),
        pl.col("alt").alias("x_alt"),
        pl.col("genotype_probabilities").alias("x_p"),
    )
    dist = _distribution(work, "x", "p")
    if dist.is_empty():
        return pl.DataFrame(schema={"_site": pl.UInt32, "_q0": pl.Float64, "_q1": pl.Float64, "_q2": pl.Float64, "_ploidy": pl.Int8})
    n = work.select("_pair", pl.col("x_source_alleles").list.len().alias("_n"), pl.col("x_p").list.len().alias("_len"))
    wide = dist.pivot(on="_copies", index="_pair", values="_p", aggregate_function="sum")
    for k in range(3):
        if str(k) not in wide.columns:
            wide = wide.with_columns(pl.lit(None, dtype=pl.Float64).alias(str(k)))
    return (
        wide.join(n, on="_pair")
        .select(
            pl.col("_pair").alias("_site"),
            pl.col("0").fill_null(0.0).alias("_q0"),
            pl.col("1").fill_null(0.0).alias("_q1"),
            pl.when(pl.col("_len") != pl.col("_n")).then(pl.col("2").fill_null(0.0)).alias("_q2"),
            pl.when(pl.col("_len") == pl.col("_n")).then(1).otherwise(2).cast(pl.Int8).alias("_ploidy"),
        )
    )


_READING = pl.Schema(
    {"_site": pl.UInt32, "_reading": pl.String, "_ploidy": pl.Int8, "_hard": pl.Int8, "_q0": pl.Float64, "_q1": pl.Float64, "_q2": pl.Float64, "_gp_lik": pl.Boolean}
)
"""One facet's reading of one allele: `_reading` is `stated`, `site`, `ref_block`, `absence` or one of
`_SILENT`; a hard reading has `_hard` copies (its ploidy, or null where the source did not say), a
graded one `_q0.._q2`."""


def _readings(
    sites: pl.DataFrame,
    evidence: pl.DataFrame,
    stated: pl.DataFrame,
    piece: Loaded,
    facet: dict,
    chrom: str,
    declared: set[str],
    *,
    pass_only: bool,
    gp_as: GpReading,
    force: bool,
    duplicates: Counter,
) -> pl.DataFrame:
    """What one facet says about each allele of `sites` on `chrom` (`_READING`)."""
    key = ["pos", "_ref_key", "alt"]
    mine = stated.sort("_row").unique(subset=key, keep="first", maintain_order=True)
    duplicates[facet["facet_id"]] += stated.height - mine.height
    hit = sites.select("_site", *key).join(mine, on=key, how="inner")
    vector = pl.col("grade_basis").cast(pl.String).is_in([GradeBasis.GL.value, GradeBasis.GQ.value]) | (
        (pl.col("grade_basis").cast(pl.String) == GradeBasis.GP.value) & pl.lit(gp_as is GpReading.LIKELIHOOD)
    )
    hit = hit.with_columns(
        (_passes(pl.col("filters")) | pl.lit(not pass_only)).alias("_pass"),
        (vector & pl.col("genotype_probabilities").is_not_null()).alias("_vector"),
        _fully_called("gt_alleles").alias("_called"),
        _copies_of("gt_alleles", "alt").alias("_count"),
        pl.col("gt_alleles").list.len().cast(pl.Int8).alias("_gt_ploidy"),
        (pl.col("grade_basis").cast(pl.String) == GradeBasis.GP.value).alias("_gp_lik"),
    )
    graded = _vector_copies(hit.filter(pl.col("_pass") & pl.col("_vector")))
    hit = hit.join(graded, on="_site", how="left")
    stated_reading = hit.select(
        "_site",
        pl.when(~pl.col("_pass")).then(pl.lit("filtered"))
        .when(pl.col("_q0").is_not_null() | pl.col("_called")).then(pl.lit("stated"))
        .otherwise(pl.lit("no_call")).alias("_reading"),
        pl.when(pl.col("_q0").is_not_null()).then(pl.col("_ploidy")).when(pl.col("_called")).then(pl.col("_gt_ploidy")).alias("_ploidy"),
        pl.when(pl.col("_q0").is_null() & pl.col("_called")).then(pl.col("_count")).alias("_hard"),
        "_q0", "_q1", "_q2",
        (pl.col("_q0").is_not_null() & pl.col("_gp_lik")).alias("_gp_lik"),
    )  # fmt: skip

    rest = sites.join(hit.select("_site"), on="_site", how="anti")
    at = evidence.group_by("pos").agg(
        pl.len().alias("_rows"),
        pl.col("_ref_key").first().alias("_other_ref"),
        pl.col("end_pos").first().alias("_other_end"),
        pl.col("gt_alleles").first().alias("_other_gt"),
        pl.col("filters").first().alias("_other_filters"),
    )
    rest = rest.join(at, on="pos", how="left")
    same_span = (pl.col("_other_ref") == pl.col("_ref_key")) | (
        pl.col("_other_ref").is_null() & (pl.col("variant_class").cast(pl.String) == VariantClass.SNV.value)
        & (pl.col("_other_end").is_null() | (pl.col("_other_end").cast(pl.Int64) == pl.col("pos").cast(pl.Int64)))
    )  # fmt: skip
    one = (pl.col("_rows") == 1).fill_null(False) & same_span.fill_null(False)
    rest = rest.with_columns(
        _fully_called("_other_gt").alias("_called"),
        _copies_of("_other_gt", "alt").alias("_count"),
        pl.col("_other_gt").list.len().cast(pl.Int8).alias("_gt_ploidy"),
        (_passes(pl.col("_other_filters")) | pl.lit(not pass_only)).alias("_pass"),
    )
    site_reading = rest.filter(pl.col("_rows").is_not_null()).select(
        "_site",
        pl.when(pl.col("_rows") > 1).then(pl.lit("several_rows")).when(~one).then(pl.lit("other_span"))
        .when(~pl.col("_pass")).then(pl.lit("filtered"))
        .when(pl.col("_called")).then(pl.lit("site"))
        .otherwise(pl.lit("no_call")).alias("_reading"),
        pl.when(one & pl.col("_pass") & pl.col("_called")).then(pl.col("_gt_ploidy")).alias("_ploidy"),
        pl.when(one & pl.col("_pass") & pl.col("_called")).then(pl.col("_count")).alias("_hard"),
    )  # fmt: skip

    silent = rest.filter(pl.col("_rows").is_null()).select("_site", "pos", "variant_class").with_columns(pl.col("pos").cast(pl.Int64).alias("_at"))
    # Inside a record that starts before the position and spans it (a deletion): not absent, and no
    # reading of this allele either.
    spans = (
        evidence.select(pl.col("pos").cast(pl.Int64).alias("_start"), pl.coalesce("end_pos", "pos").cast(pl.Int64).alias("_end"))
        .sort("_start")
        .with_columns(pl.col("_end").cum_max().alias("_reach"))
    )
    probe = silent.sort("_at").join_asof(spans.select("_start", "_reach"), left_on="_at", right_on="_start", strategy="backward", allow_exact_matches=False)
    probe = probe.with_columns((pl.col("_reach").is_not_null() & (pl.col("_reach") >= pl.col("_at"))).alias("_covered"))
    blocks = (
        piece.blocks.filter(pl.col("chrom") == chrom)
        .select(pl.col("pos").cast(pl.Int64).alias("_bstart"), pl.col("end_pos").cast(pl.Int64).alias("_bend"), pl.col("kind").cast(pl.String).alias("_kind"), pl.col("gl").cast(pl.List(pl.Float64)).alias("_gl"))
        .sort("_bstart")
        .collect()
    )  # fmt: skip
    if blocks.height:
        probe = probe.join_asof(blocks, left_on="_at", right_on="_bstart", strategy="backward")
        probe = probe.with_columns(pl.when(pl.col("_bend") >= pl.col("_at")).then(pl.col("_kind")).alias("_kind"))
    else:
        probe = probe.with_columns(pl.lit(None, dtype=pl.String).alias("_kind"), pl.lit(None, dtype=pl.List(pl.Float64)).alias("_gl"))
    reads_reference = [cls for cls in _SMALL if absence_is_reference(facet, cls, force=force) is True] if chrom in declared else []
    gl = pl.col("_gl")
    usable = (gl.list.len().is_between(2, 3) & gl.list.eval(pl.element().is_finite() | (pl.element() == float("-inf"))).list.all() & gl.list.max().is_finite()).fill_null(False)
    weights = gl.list.eval(((pl.element() - pl.element().max()) * _LN10).exp())
    total = weights.list.sum()
    is_ref = (pl.col("_kind") == "ref").fill_null(False)
    graded_block = is_ref & usable
    absence_reading = probe.select(
        "_site",
        pl.when(pl.col("_covered")).then(pl.lit("covered"))
        .when(is_ref).then(pl.lit("ref_block"))
        .when(pl.col("_kind").is_not_null()).then(pl.lit("no_call_block"))
        .when(pl.col("variant_class").cast(pl.String).is_in(reads_reference)).then(pl.lit("absence"))
        .otherwise(pl.lit("unknown")).alias("_reading"),
        pl.when(~pl.col("_covered") & graded_block).then(gl.list.len() - 1).cast(pl.Int8).alias("_ploidy"),
        pl.when(~pl.col("_covered") & graded_block).then(None)
        .when(~pl.col("_covered") & (is_ref | (pl.col("_kind").is_null() & pl.col("variant_class").cast(pl.String).is_in(reads_reference)))).then(0)
        .cast(pl.Int8).alias("_hard"),
        pl.when(~pl.col("_covered") & graded_block).then(weights.list.get(0) / total).alias("_q0"),
        pl.when(~pl.col("_covered") & graded_block).then(weights.list.get(1) / total).alias("_q1"),
        pl.when(~pl.col("_covered") & graded_block & (gl.list.len() == 3)).then(weights.list.get(2, null_on_oob=True) / total).alias("_q2"),
    )  # fmt: skip
    pieces = [frame for frame in (stated_reading, site_reading, absence_reading) if frame.height]
    if not pieces:
        return pl.DataFrame(schema=_READING)
    return pl.concat(
        [frame.select([pl.col(name).cast(dtype) if name in frame.columns else pl.lit(None, dtype=dtype).alias(name) for name, dtype in _READING.items()]) for frame in pieces]
    ).with_columns(pl.col("_gp_lik").fill_null(False))


# --------------------------------------------------------------------------------------------------
# The posterior


def _site_rates(sites: pl.DataFrame, rates: pl.DataFrame, regions: RegionSet | None) -> pl.DataFrame:
    """Per allele and parent, the genotype and silence error rates it is weighted by: the highest rate
    among the region classes holding the allele (`all` where none does), and whether any of those
    came from the default."""
    probe = sites.select("_site", "chrom", "variant_class", pl.col("pos").cast(pl.Int64).alias("_at")).sort("chrom", "_at")
    classes = sorted(regions.intervals) if regions is not None else []
    members = {label: _covered(probe, regions.intervals[label]) for label in classes} if regions is not None else {}
    held = probe.with_columns(*[series.alias(f"_in_{label}") for label, series in members.items()])
    none_held = ~pl.any_horizontal([pl.col(f"_in_{label}") for label in classes]) if classes else pl.lit(True)
    out = []
    for facet_id in rates["facet_id"].unique(maintain_order=True).to_list():
        own = rates.filter(pl.col("facet_id") == facet_id)

        def lookup(label: str, column: str) -> pl.Expr:
            table = own.filter(pl.col("region_class") == label)
            mapping = dict(zip(table["variant_class"].cast(pl.String).to_list(), table[column].to_list()))
            return pl.col("variant_class").cast(pl.String).replace_strict(mapping, default=None, return_dtype=pl.Float64)

        def defaulted(label: str) -> pl.Expr:
            table = own.filter(pl.col("region_class") == label)
            mapping = dict(zip(table["variant_class"].cast(pl.String).to_list(), (b == RateBasis.DEFAULT.value for b in table["rate_basis"].cast(pl.String).to_list())))
            return pl.col("variant_class").cast(pl.String).replace_strict(mapping, default=False, return_dtype=pl.Boolean)

        def pick(column: str) -> pl.Expr:
            inside = [pl.when(pl.col(f"_in_{label}")).then(lookup(label, column)) for label in classes]
            return pl.when(none_held).then(lookup(ALL, column)).otherwise(pl.max_horizontal(inside) if inside else lookup(ALL, column))

        flag = pl.when(none_held).then(defaulted(ALL)).otherwise(pl.any_horizontal([pl.col(f"_in_{label}") & defaulted(label) for label in classes]) if classes else defaulted(ALL))
        out.append(held.select("_site", pl.lit(facet_id).alias("_facet"), pick("genotype_error").alias("_ge"), pick("silence_error").alias("_se"), flag.alias("_defaulted")))
    return pl.concat(out)


def _posterior(readings: pl.DataFrame, prior: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Per allele: ploidy, the log posterior's three entries and the most probable copies; and per
    reading, whether it supported, contradicted or was silent."""
    speaking = ~pl.col("_reading").is_in(list(_SILENT))
    ploidy = (
        readings.filter(speaking & pl.col("_ploidy").is_in([1, 2]))
        .group_by("_site")
        .agg((pl.col("_ploidy") == 1).sum().alias("_p1"), (pl.col("_ploidy") == 2).sum().alias("_p2"))
        .select("_site", pl.when(pl.col("_p1") > pl.col("_p2")).then(1).otherwise(2).cast(pl.Int8).alias("_P"))
    )
    work = readings.join(ploidy, on="_site", how="left").with_columns(pl.col("_P").fill_null(2))
    P = pl.col("_P")
    work = work.with_columns(
        pl.when(speaking & pl.col("_ploidy").is_not_null() & (pl.col("_ploidy") != P)).then(pl.lit("ploidy_differs")).otherwise(pl.col("_reading")).alias("_reading")
    )
    speaking = ~pl.col("_reading").is_in(list(_SILENT))
    hard = pl.col("_hard")
    q = [pl.coalesce(pl.col(f"_q{k}"), (hard == k).cast(pl.Float64)) for k in range(3)]
    e = pl.when(pl.col("_reading").is_in(["stated", "site"])).then(pl.col("_ge")).otherwise(pl.col("_se"))
    likelihood = [((1 - e) * q[k] + e * (1 - q[k]) / P.cast(pl.Float64)).log() for k in range(3)]
    work = work.with_columns(
        *[pl.when(speaking & (pl.lit(k) <= P)).then(likelihood[k]).alias(f"_l{k}") for k in range(3)],
        # The copies this facet says: the top of its distribution, the lowest on a tie.
        pl.when(speaking).then(
            pl.when(q[0] >= pl.max_horizontal(q[1], pl.when(P == 2).then(q[2]).otherwise(-1.0))).then(0)
            .when(q[1] >= pl.when(P == 2).then(q[2]).otherwise(-1.0)).then(1)
            .otherwise(2)
        ).cast(pl.Int8).alias("_says"),
    )  # fmt: skip
    if prior != "flat":
        raise ValueError(f"prior is 'flat' (every number of copies equally likely beforehand), not {prior!r}")
    sums = work.group_by("_site").agg(
        pl.col("_P").first(),
        *[pl.col(f"_l{k}").sum().alias(f"_s{k}") for k in range(3)],
    )
    s = [pl.col(f"_s{k}") for k in range(3)]
    s2 = pl.when(pl.col("_P") == 2).then(s[2]).otherwise(float("-inf"))
    top = pl.max_horizontal(s[0], s[1], s2)
    copies = pl.when(s[0] >= pl.max_horizontal(s[1], s2)).then(0).when(s[1] >= s2).then(1).otherwise(2).cast(pl.Int8)
    sums = sums.with_columns(
        top.alias("_top"), copies.alias("_copies")
    ).with_columns(
        *[((pl.col(f"_s{k}") - pl.col("_top")) if k < 2 else (s2 - pl.col("_top"))).exp().alias(f"_w{k}") for k in range(3)]
    )  # fmt: skip
    total = pl.col("_w0") + pl.col("_w1") + pl.col("_w2")
    others = total - pl.when(pl.col("_copies") == 0).then(pl.col("_w0")).when(pl.col("_copies") == 1).then(pl.col("_w1")).otherwise(pl.col("_w2"))
    sums = sums.with_columns(
        *[(pl.col(f"_w{k}") / total).alias(f"_p{k}") for k in range(3)],
        (-10.0 * (others.log() - total.log()) / _LN10).alias("confidence"),
    )
    votes = work.join(sums.select("_site", "_copies"), on="_site").select(
        "_site", "_facet", "_defaulted", "_gp_lik", "_reading",
        pl.when(~speaking).then(pl.lit("silent")).when(pl.col("_says") == pl.col("_copies")).then(pl.lit("supported")).otherwise(pl.lit("contradicted")).alias("_vote"),
    )  # fmt: skip
    return sums, votes


# --------------------------------------------------------------------------------------------------
# The facet


def consensus(
    person: Person | Loaded | Iterable[Loaded],
    *,
    atlas: pl.DataFrame | None,
    regions: RegionSet | None = None,
    prior: str = "flat",
    default_rate: float | None = None,
    scope: str = "carried_pass",
    pass_only: bool = True,
    gp_as: GpReading | str = GpReading.CALL,
    gq_on_reference: bool = False,
    force: bool = False,
) -> Consensus:
    """A consensus facet over one person's facets: per allele, a posterior over its copies (module
    docstring, docs/SCHEMA.md "A consensus posterior").

    `person` is a `Person` (refused when its check read two people) or `Loaded`s of one subject, two
    facets or more, none made from another (the same reads would count twice), all on one assembly
    (equal `ref_seqcol_digest`, else one equal `assembly_label` with a warning; otherwise lift first).
    `atlas` is a `method_atlas` table; its pair rates become per-method rates (`method_error_rates`,
    rows of `scope`). `regions` places each allele in its region classes, and a facet is weighted by
    the highest rate among them; without it, by its `all` rate. A facet whose caller the atlas never
    measured gets `default_rate`, or, when that is None, the worst rate a measured method has in the
    class; `atlas=None` needs `default_rate`. `pass_only` reads a row that failed its filters as
    silent (and takes the atlas rows of scope `carried_pass`); `gp_as` says how a posterior is read
    (`GpReading`); `gq_on_reference` and `force` are `genotype_evidence`'s and `absence_is_reference`'s.

    The rows are every sequence ALT allele a parent states in a row it would be read by, on every
    contig, in canonical order. Each is a `calls` row of the consensus facet: `alts = [alt]`,
    `source_alleles = [REF, ALT]`, `gp` the posterior over 0..ploidy copies (VCF genotype order over
    `[REF, ALT]`), `gt` the most probable, `ref` the REF where a parent read it from the reference.
    """
    gp_as = GpReading(gp_as)
    if atlas is None and default_rate is None:
        raise ValueError("without an atlas every facet is weighted by the default rate; pass default_rate=...")
    pieces, parent_links, warnings = _parents(person)
    rows = [piece.facets.collect().row(0, named=True) for piece in pieces]
    label, digest = _assembly(rows, warnings)
    if regions is not None and regions.assembly_label != label:
        raise ValueError(f"the regions are on {regions.assembly_label}, the facets on {label}; lift one first")
    classes = regions.classes if regions is not None else (ALL,)
    scope = "carried_pass" if pass_only else ("carried" if scope == "carried_pass" else scope)
    rates = _facet_rates(rows, atlas, classes, default_rate, scope)
    defaulted = sorted({r["facet_id"] for r in rates.filter(pl.col("rate_basis") == RateBasis.DEFAULT.value).iter_rows(named=True)})
    if defaulted:
        warnings.append(f"{len(defaulted)} facet(s) weighted by a default rate in some class (no atlas row for their caller): {', '.join(defaulted)}")

    options = {"prior": prior, "default_rate": default_rate, "scope": scope, "pass_only": pass_only, "gp_as": gp_as.value,
               "gq_on_reference": gq_on_reference, "force": force, "regions": regions.name if regions is not None else None}  # fmt: skip
    rates_digest = hashlib.sha256(rates.write_csv().encode()).hexdigest()
    facet_id = hashlib.sha256(
        json.dumps({"derived_from": [r["facet_id"] for r in rows], "operation": CONSENSUS, "version": CONSENSUS_VERSION, "options": options, "rates": rates_digest}, sort_keys=True).encode()
    ).hexdigest()[:32]

    contigs = sorted(
        {chrom for piece in pieces for chrom in piece.calls.select(pl.col("chrom").unique()).collect()["chrom"].to_list()},
        key=lambda c: (CONTIG_RANK.get(c, len(CONTIG_RANK)), c),
    )
    declared = {row["facet_id"]: _declared_contigs(piece, row) for piece, row in zip(pieces, rows)}
    silent_counts: Counter[tuple[str, str]] = Counter()
    duplicates: Counter[str] = Counter()
    calls_out, support_out = [], []
    for chrom in contigs:
        evidence = [_evidence(piece, chrom, gq_on_reference) for piece in pieces]
        stated = [_stated_alleles(frame) for frame in evidence]
        reads = [
            frame.filter(
                _passes(pl.col("filters")) if pass_only else pl.lit(True),
                ~pl.col("grade_basis").cast(pl.String).is_in([GradeBasis.NO_CALL.value, GradeBasis.DOSAGE.value]),
            )
            for frame in stated
        ]
        union = pl.concat([frame.select("pos", "_ref_key", "alt", "variant_class", pl.col("_ref_read").alias("_read")) for frame in reads]) if reads else None
        if union is None or union.is_empty():
            continue
        sites = (
            union.group_by("pos", "_ref_key", "alt", "variant_class").agg(pl.col("_read").any())
            .sort("pos", "_ref_key", "alt")
            .with_columns(pl.int_range(pl.len(), dtype=pl.UInt32).alias("_site"), pl.lit(chrom).alias("chrom"))
        )  # fmt: skip
        site_rates = _site_rates(sites, rates, regions)
        per_facet = []
        for piece, row, frame, own in zip(pieces, rows, evidence, stated):
            reading = _readings(
                sites, frame, own, piece, row, chrom, declared[row["facet_id"]], pass_only=pass_only, gp_as=gp_as, force=force, duplicates=duplicates
            ).with_columns(pl.lit(row["facet_id"]).alias("_facet"))
            per_facet.append(reading)
        readings = pl.concat(per_facet).join(site_rates, on=["_site", "_facet"], how="left")
        sums, votes = _posterior(readings, prior)
        for facet, reason, count in votes.filter(pl.col("_vote") == "silent").group_by("_facet", "_reading").len().iter_rows():
            silent_counts[(facet, reason)] += count
        tally = votes.group_by("_site").agg(
            pl.len().cast(pl.UInt32).alias("facets"),
            *[(pl.col("_vote") == v).sum().cast(pl.UInt32).alias(v) for v in ("supported", "contradicted", "silent")],
            ((pl.col("_vote") != "silent") & pl.col("_defaulted").fill_null(False)).sum().cast(pl.UInt32).alias("defaulted"),
            ((pl.col("_vote") != "silent") & pl.col("_gp_lik")).sum().cast(pl.UInt32).alias("gp_as_likelihood"),
        )
        done = sites.join(sums, on="_site").join(tally, on="_site").sort("_site")
        support_out.append(
            done.select(
                "chrom", "pos", pl.col("_ref_key").alias("ref"), "alt", "variant_class", pl.col("_P").alias("ploidy"), pl.col("_copies").alias("copies"),
                "confidence", "facets", "supported", "contradicted", "silent", "defaulted", "gp_as_likelihood",
            ).cast(dict(CONSENSUS_SUPPORT))
        )  # fmt: skip
        calls_out.append(_calls(done, facet_id))

    for facet, count in sorted(duplicates.items()):
        if count:
            warnings.append(f"facet {facet}: {count} allele(s) stated by several rows; the first row in file order was read")
    by_facet: dict[str, list[str]] = {}
    for (facet, reason), count in sorted(silent_counts.items()):
        by_facet.setdefault(facet, []).append(f"{reason} {count}")
    for facet, parts in by_facet.items():
        warnings.append(f"facet {facet} silent on alleles others stated: {', '.join(parts)}")
    for line in warnings:
        logger.warning("consensus %s: %s", facet_id, line)

    calls = pl.concat(calls_out) if calls_out else pl.DataFrame(schema=CALLS)
    support = pl.concat(support_out) if support_out else pl.DataFrame(schema=CONSENSUS_SUPPORT)
    facet_row = _facet_row(facet_id, rows, label, digest, calls, rates_digest)
    added = links_frame(
        [
            {"from_facet_id": row["facet_id"], "to_facet_id": facet_id, "relation": CONSENSUS,
             "provenance": {"step": "consensus", "step_version": str(CONSENSUS_VERSION), "derivation_inputs": []}}
            for row in rows
        ]
    )  # fmt: skip
    return Consensus(
        facets=facet_row.lazy(),
        calls=calls.lazy(),
        blocks=pl.LazyFrame(schema=BLOCKS),
        links=canonical_links(pl.concat([parent_links, added], how="vertical")),
        warnings=tuple(warnings),
        schema_version=SCHEMA_VERSION,
        support=support.lazy(),
        rates=rates,
    )


def _calls(done: pl.DataFrame, facet_id: str) -> pl.DataFrame:
    """The consensus facet's `calls` rows for one contig's alleles."""
    P, copies = pl.col("_P"), pl.col("_copies")
    gt = (
        pl.when(P == 1).then(pl.concat_list(copies.cast(pl.UInt8)))
        .when(copies == 0).then(pl.lit([0, 0], dtype=pl.List(pl.UInt8)))
        .when(copies == 1).then(pl.lit([0, 1], dtype=pl.List(pl.UInt8)))
        .otherwise(pl.lit([1, 1], dtype=pl.List(pl.UInt8)))
    )  # fmt: skip
    probabilities = pl.when(P == 1).then(pl.concat_list("_p0", "_p1")).otherwise(pl.concat_list("_p0", "_p1", "_p2"))
    ref = pl.col("_ref_key")
    frame = done.select(
        pl.lit(facet_id).alias("facet_id"),
        "chrom",
        pl.col("pos").cast(pl.UInt32),
        (pl.col("pos").cast(pl.Int64) + pl.max_horizontal(ref.str.len_chars().cast(pl.Int64), pl.lit(1)) - 1).cast(pl.UInt32).alias("end_pos"),
        pl.when(pl.col("_read")).then(ref).alias("ref"),
        ref.alias("source_ref"),
        pl.lit(RefCheck.NOT_RUN.value).alias("ref_check"),
        pl.concat_list("alt").alias("alts"),
        pl.concat_list(ref, "alt").alias("source_alleles"),
        gt.alias("gt"),
        probabilities.cast(pl.List(pl.Float32)).alias("gp"),
        pl.lit(EvidenceBasis.ASSAY.value).alias("evidence_basis"),
        ref.alias("_r"),
        pl.col("alt").alias("_a"),
    )
    first = pl.when(pl.col("gt").list.get(0) == 0).then(pl.col("_r")).otherwise(pl.col("_a"))
    second = pl.when(pl.col("gt").list.get(1, null_on_oob=True) == 0).then(pl.col("_r")).otherwise(pl.col("_a"))
    gt_alleles = pl.when(pl.col("gt").list.len() == 1).then(pl.concat_list(first)).otherwise(pl.concat_list(first, second))
    frame = frame.with_columns(gt_alleles.alias("gt_alleles"))
    return frame.select([pl.col(name).cast(dtype) if name in frame.columns else pl.lit(None, dtype=dtype).alias(name) for name, dtype in CALLS.items()])


def _facet_row(facet_id: str, rows: list[dict], label: str | None, digest: str | None, calls: pl.DataFrame, rates_digest: str) -> pl.DataFrame:
    """The consensus facet's `facets` row."""

    def shared(name: str) -> object:
        values = {json.dumps(row[name], sort_keys=True, default=str) for row in rows}
        return rows[0][name] if len(values) == 1 else None

    counted = (
        calls.with_columns(variant_class_expr(pl.col("source_ref"), pl.col("alts").list.first()).alias("_v"))
        .group_by("_v").len()
    )  # fmt: skip
    classes = sorted(({"variant_class": v, "alleles": n} for v, n in counted.iter_rows() if v is not None), key=lambda item: _SMALL.index(item["variant_class"]))
    row = {name: None for name in FACETS.names()}
    row.update(
        facet_id=facet_id,
        caller=CALLER,
        caller_version=str(CONSENSUS_VERSION),
        assembly_label=label,
        assembly_label_basis=shared("assembly_label_basis") if label is not None else None,
        ref_seqcol_digest=digest,
        callset_type=shared("callset_type"),
        absence_means=AbsenceMeans.UNKNOWN.value,
        variant_classes=classes,
        method_steps=[{"tool": CALLER, "version": str(CONSENSUS_VERSION), "role": "merge", "options": [], "line": None}],
        relation=CONSENSUS,
        contig_map=[{"source": chrom, "chrom": chrom} for chrom in sorted(set(calls["chrom"].to_list()))] if calls.height else [],
        derivation_inputs=[{"role": "error_rates", "name": "method_error_rates", "sha256": rates_digest}],
        subject=shared("subject"),
    )
    return pl.DataFrame([row], schema=FACETS)
