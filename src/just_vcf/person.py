"""Facets of one person, grouped by the subject the caller supplied (RM273, docs/SCHEMA.md).

The caller says whose a file is with `load(..., subject=...)`. This module groups facets by that label
and checks the claim; it never infers a subject and never moves a facet to another one. The check is a
relatedness test over pairs of facets, passed in as a callable (`SameSubjectCheck`): `True` says one
person, `False` two, `None` that it could not tell. A `False` is a warning naming both facets, and the
facets stay where the caller put them. The default check is `just_vcf.relatedness.same_person` (RM274),
KING-robust kinship over the package's site panel; `check=None` runs none, and every pair is then
`not_run`, which is not a pass.

A person also carries the coverage facets (RM278) the caller labels with its subject
(`read_coverage(..., subject=...)`, `alignment_coverage(..., subject=...)`). They are grouped by that
label alone: a depth file holds no genotype, so no relatedness check reads one.
"""

import logging
from collections import defaultdict
from collections.abc import Callable, Iterable
from enum import Enum
from itertools import combinations

import polars as pl
from pydantic import BaseModel, ConfigDict, Field

from just_vcf.coverage import Coverage
from just_vcf.links import MULTI_PARENT_RELATIONS, SINGLE_PARENT_RELATIONS, canonical_links
from just_vcf.loading import Loaded
from just_vcf.pedigree import BAND_ORDER, EXPECTED_KINSHIP, PEDIGREE_REASONS, Pedigree, expected_kinship
from just_vcf.relatedness import (
    BETWEEN_BANDS,
    DEFAULT_MIN_SITES,
    INTERVAL_SES,
    RELATEDNESS_FRAME,
    RELATEDNESS_REASONS,
    RelationBand,
    SitePanel,
    band_of,
    finest_bands,
    first_degree_kind,
    pairwise_relatedness,
    same_person,
)
from just_vcf.schema import LINKS, SCHEMA_VERSION, empty, enum_dtype

logger = logging.getLogger(__name__)

SameSubjectCheck = Callable[[Loaded, Loaded], bool | None]
"""A relatedness check of two facets, each passed as a `Loaded` holding that one facet: `True` = one
person, `False` = two people, `None` = could not tell."""


class SubjectCheckOutcome(str, Enum):
    """What the relatedness check said about one pair of facets under one subject."""

    SAME = "same"
    DIFFERENT = "different"
    COULD_NOT_TELL = "could_not_tell"
    """The check ran and answered `None`."""
    NOT_RUN = "not_run"
    """No check was given."""


SUBJECT_CHECKS = pl.Schema(
    {"facet_id": pl.String, "other_facet_id": pl.String, "outcome": enum_dtype(SubjectCheckOutcome)}
)
"""`Person.checks`: one row per pair of sources checked, by the facet that stood for each source."""


class SubjectError(ValueError):
    """Facets that cannot be grouped as given: a facet without a subject, or two subjects in one person."""


class Person(BaseModel):
    """The facets one subject label groups, their rows, and what the relatedness check said.

    The frames are the facets' own, concatenated in `facet_id` order, so `calls` and `blocks` keep the
    canonical order (SCHEMA.md Decision 11). `links` are every link the facets' loads carried, once
    each, in `canonical_links` order (Decision 12).
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    subject: str
    """The caller's label. It identifies a person if the caller makes it so; nothing here logs it."""
    facets: pl.LazyFrame
    calls: pl.LazyFrame
    blocks: pl.LazyFrame
    links: pl.LazyFrame = Field(default_factory=lambda: empty(LINKS))
    checks: pl.DataFrame
    """`SUBJECT_CHECKS`: each pair of the person's sources and what the check said about it."""
    warnings: tuple[str, ...] = ()
    """One line per reason, naming facet ids, never the subject."""
    schema_version: str = SCHEMA_VERSION
    coverage: tuple[Coverage, ...] = ()
    """The coverage facets under this subject, in `facet_id` order (RM278). Grouped by the caller's
    label and nothing else: they are in no row of `checks`, since no relatedness test reads depth.
    `just_vcf.karyotype.person_sex_chromosome_evidence` pairs them with the callsets."""

    @property
    def one_person(self) -> bool | None:
        """Kleene AND over `checks`: `False` when a pair read as two people, `True` only when every pair
        was checked and read as one, `None` otherwise (a pair not run or not told, or one source alone,
        where there is nothing to check)."""
        outcomes = set(self.checks["outcome"].to_list())
        if SubjectCheckOutcome.DIFFERENT.value in outcomes:
            return False
        if outcomes == {SubjectCheckOutcome.SAME.value}:
            return True
        return None

    @classmethod
    def of(
        cls, *loaded: Loaded, check: SameSubjectCheck | None = same_person, coverage: Iterable[Coverage] = ()
    ) -> "Person":
        """One person from facets that all carry one subject.

        A facet without a subject, or facets under two subjects, are refused with `SubjectError`; so is
        one `facet_id` given twice, whose rows would then count twice. Facets of one source (the same
        bytes, archive member and sample: a file and the facets derived from it, or the facets one file
        loads as) are one person by construction and are not checked against each other; nor are
        facets joined by one-parent lineage across files (an imputed facet and the chip the caller
        stated it was imputed from, `state_imputed_from`), whose rows are not independent evidence.
        Every other pair of such callsets is checked once, each standing as a facet read from a file
        with no lineage (the chip before its imputation), else its lowest `facet_id`.

        `coverage` is coverage facets of the same subject, carried on `Person.coverage`. One without a
        subject, or under another subject than the facets, is refused with `SubjectError`; one
        `facet_id` given twice with `ValueError`. Their label is the only check they get.
        """
        pieces = _pieces(loaded)
        if not pieces:
            raise SubjectError("no facets given")
        facets = [piece.facets.collect().row(0, named=True) for piece in pieces]
        missing = sorted(row["facet_id"] for row in facets if row["subject"] is None)
        if missing:
            raise SubjectError(f"facet(s) without a subject: {', '.join(missing)}; load them with subject=...")
        by_subject: dict[str, list[str]] = defaultdict(list)
        for row in facets:
            by_subject[row["subject"]].append(row["facet_id"])
        if len(by_subject) > 1:
            groups = "; ".join(", ".join(sorted(ids)) for ids in sorted(by_subject.values()))
            raise SubjectError(f"facets under {len(by_subject)} different subjects cannot be one person: {groups}")
        ids = [row["facet_id"] for row in facets]
        twice = sorted({facet_id for facet_id in ids if ids.count(facet_id) > 1})
        if twice:
            raise ValueError(f"facet(s) given more than once: {', '.join(twice)}")
        versions = sorted({piece.schema_version for piece in pieces})
        if len(versions) > 1:
            raise ValueError(f"facets of different schema versions cannot be grouped: {', '.join(versions)}")
        held = _coverage_of(facets[0]["subject"], tuple(coverage))

        order = sorted(range(len(pieces)), key=lambda i: ids[i])
        pieces = [pieces[i] for i in order]
        facets = [facets[i] for i in order]
        links = canonical_links(pl.concat([piece.links for piece in pieces], how="vertical"))
        checks, warnings = _check(pieces, facets, links.collect(), check)
        for warning in warnings:
            logger.warning("%s", warning)
        return cls(
            subject=facets[0]["subject"],
            facets=pl.concat([piece.facets for piece in pieces], how="vertical"),
            calls=pl.concat([piece.calls for piece in pieces], how="vertical"),
            blocks=pl.concat([piece.blocks for piece in pieces], how="vertical"),
            links=links,
            checks=checks,
            warnings=tuple(warnings),
            schema_version=versions[0],
            coverage=held,
        )


def _unlabelled_coverage(coverage: tuple[Coverage, ...]) -> None:
    missing = sorted(item.facet.facet_id for item in coverage if item.facet.subject is None)
    if missing:
        raise SubjectError(f"coverage facet(s) without a subject: {', '.join(missing)}; read them with subject=...")


def _coverage_of(subject: str, coverage: tuple[Coverage, ...]) -> tuple[Coverage, ...]:
    """The coverage facets of `subject`, in `facet_id` order; anything else is refused."""
    _unlabelled_coverage(coverage)
    other = sorted(item.facet.facet_id for item in coverage if item.facet.subject != subject)
    if other:
        raise SubjectError(f"coverage facet(s) under another subject than the facets: {', '.join(other)}")
    ids = [item.facet.facet_id for item in coverage]
    twice = sorted({facet_id for facet_id in ids if ids.count(facet_id) > 1})
    if twice:
        raise ValueError(f"coverage facet(s) given more than once: {', '.join(twice)}")
    return tuple(sorted(coverage, key=lambda item: item.facet.facet_id))


class PedigreeAgreement(str, Enum):
    """What the measured relatedness of two facets says about the kinship a pedigree expects."""

    AGREES = "agrees"
    """Every band the measured interval reaches is one the pedigree allows."""
    DISAGREES = "disagrees"
    """No band the measured interval reaches is one the pedigree allows. A first-degree kind is a band
    of its own here: a stated parent who reads `full_siblings` disagrees."""
    CANNOT_TELL = "cannot_tell"
    """The interval reaches allowed and other bands, the first degree's kind is not told where the
    pedigree allows one, or one side is unknown; `reason` says which (`PEDIGREE_CHECK_REASONS`)."""


NO_FACET = "no_facet"
"""A pedigree individual has no facet under its id as subject."""
KIND_NOT_TOLD = "kind_not_told"
"""The kinship's interval lies wholly in the first degree, the pedigree allows one of its two kinds, and
the IBS0 count does not tell which kind the pair is (its interval is too thin, or there is none)."""
PEDIGREE_CHECK_REASONS: frozenset[str] = frozenset({NO_FACET, KIND_NOT_TOLD}) | PEDIGREE_REASONS | RELATEDNESS_REASONS
"""Why a `PEDIGREE_CHECKS` row is `cannot_tell`: its own reasons, a pedigree's (`expected_reason`), or
the measured relatedness's (`between_bands` also when the interval reaches allowed and other bands).
Grows additively; `check_pedigree` refuses to write any other."""

PEDIGREE_CHECKS = pl.Schema(
    {
        **EXPECTED_KINSHIP,
        "facet_a": pl.String,
        "facet_b": pl.String,
        "band": RELATEDNESS_FRAME["band"],
        "kinship": pl.Float64,
        "kinship_se": pl.Float64,
        "ibs0": pl.Float64,
        "n_sites": pl.UInt32,
        "relatedness_reason": pl.String,
        "agreement": enum_dtype(PedigreeAgreement),
        "reason": pl.String,
    }
)
"""`check_pedigree`: `EXPECTED_KINSHIP`'s pair, a facet of each (`facet_a` the individual's, `facet_b`
the relative's), the measured `relatedness` (`band`, `kinship`, `kinship_se`, `ibs0`, `n_sites`, and
its `reason` as `relatedness_reason`), and the `agreement`, with `reason` when it is `cannot_tell`."""


def _measured_ladder(kinship: float, se: float | None) -> list[RelationBand]:
    """The bands of `BAND_ORDER` the interval `kinship ± INTERVAL_SES·se` reaches (just the kinship's
    band without an SE), the first degree as one band."""
    if se is None:
        return [band_of(kinship)]
    low, high = BAND_ORDER.index(band_of(kinship - INTERVAL_SES * se)), BAND_ORDER.index(band_of(kinship + INTERVAL_SES * se))
    return list(BAND_ORDER[low : high + 1])


def _agreement(expected: dict, measured: dict) -> tuple[PedigreeAgreement, str | None]:
    """The agreement of a measured pair with what the pedigree expects, and the reason when it is
    `cannot_tell`. The measured bands are the finest the interval reaches, the first degree as the kind
    its IBS0 count tells, or both kinds where it does not (no count tells none). Where the interval lies
    wholly in the first degree and only the kind is not told, the reason is `kind_not_told`; where it
    reaches allowed and other bands, `between_bands`."""
    if expected["expected_reason"] is not None:
        return PedigreeAgreement.CANNOT_TELL, expected["expected_reason"]
    # A kinship kept beside a reason other than its interval's (`variant_sites_only`) names no band.
    if measured["kinship"] is None or measured["reason"] not in (None, BETWEEN_BANDS):
        return PedigreeAgreement.CANNOT_TELL, measured["reason"]
    allowed = set().union(*(finest_bands(RelationBand(value)) for value in expected["expected_bands"]))
    ladder = _measured_ladder(measured["kinship"], measured["kinship_se"])
    kind = first_degree_kind(measured["n_ibs0"], measured["n_sites"])
    reached = set().union(*(finest_bands(kind if band is RelationBand.FIRST_DEGREE else band) for band in ladder))
    if reached <= allowed:
        return PedigreeAgreement.AGREES, None
    if not reached & allowed:
        return PedigreeAgreement.DISAGREES, None
    if ladder == [RelationBand.FIRST_DEGREE] and kind is RelationBand.FIRST_DEGREE:
        return PedigreeAgreement.CANNOT_TELL, KIND_NOT_TOLD
    return PedigreeAgreement.CANNOT_TELL, BETWEEN_BANDS


def check_pedigree(
    pedigree: Pedigree,
    *sources: Person | Loaded,
    panel: SitePanel | None = None,
    min_sites: int = DEFAULT_MIN_SITES,
) -> pl.DataFrame:
    """A supplied pedigree set beside the measured relatedness of its individuals' facets, as
    `PEDIGREE_CHECKS`: one row per pair of individuals (`expected_kinship`) and pair of their facets.

    An individual's facets are those whose `subject` is its `individual_id` (RM273), each source once,
    standing as `Person.checks` has it. A pair with an individual that has none is one row,
    `cannot_tell`, `no_facet`. Relatedness is `just_vcf.relatedness`'s KING-robust estimate (`panel`
    and `min_sites` passed on). The pedigree is never changed: a disagreement is a row, not a
    correction (charter invariant 8). Facets under subjects the pedigree does not list are not read.
    """
    expected = expected_kinship(pedigree)
    pieces: list[Loaded] = []
    for source in sources:
        if isinstance(source, Person):
            source = Loaded(facets=source.facets, calls=source.calls, blocks=source.blocks, links=source.links, schema_version=source.schema_version)
        pieces.extend(_pieces((source,)))
    rows = [piece.facets.collect().row(0, named=True) for piece in pieces]
    missing = sorted(row["facet_id"] for row in rows if row["subject"] is None)
    if missing:
        raise SubjectError(f"facet(s) without a subject: {', '.join(missing)}; load them with subject=...")
    listed = set(pedigree.individuals["individual_id"].to_list())
    by_subject: dict[str, list[int]] = defaultdict(list)
    for index in sorted(range(len(rows)), key=lambda i: rows[i]["facet_id"]):
        if rows[index]["subject"] in listed:
            by_subject[rows[index]["subject"]].append(index)
    links = canonical_links(pl.concat([piece.links for piece in pieces], how="vertical")).collect() if pieces else empty(LINKS).collect()
    standing = {
        subject: [indices[i] for i in _standing([rows[index] for index in indices], links)] for subject, indices in by_subject.items()
    }

    used = sorted({index for indices in standing.values() for index in indices})
    position = {index: place for place, index in enumerate(used)}
    plan: list[tuple[dict, int | None, int | None]] = []
    for pair in expected.iter_rows(named=True):
        first, second = standing.get(pair["individual_id"], []), standing.get(pair["relative_id"], [])
        if not first or not second:
            plan.append((pair, None, None))
            continue
        plan.extend((pair, a, b) for a in first for b in second)
    measured_pairs = sorted({(min(position[a], position[b]), max(position[a], position[b])) for _, a, b in plan if a is not None})
    measured = pairwise_relatedness([pieces[index] for index in used], panel, min_sites=min_sites, pairs=measured_pairs)
    by_pair = {frozenset((row["facet_a"], row["facet_b"])): row for row in measured.iter_rows(named=True)}

    out = []
    for pair, a, b in plan:
        if a is None:
            out.append({**pair, "agreement": PedigreeAgreement.CANNOT_TELL.value, "reason": NO_FACET})
            continue
        facet_a, facet_b = rows[a]["facet_id"], rows[b]["facet_id"]
        found = by_pair[frozenset((facet_a, facet_b))]
        agreement, reason = _agreement(pair, found)
        out.append(
            {
                **pair,
                "facet_a": facet_a,
                "facet_b": facet_b,
                "band": found["band"],
                "kinship": found["kinship"],
                "kinship_se": found["kinship_se"],
                "ibs0": found["ibs0"],
                "n_sites": found["n_sites"],
                "relatedness_reason": found["reason"],
                "agreement": agreement.value,
                "reason": reason,
            }
        )
    unknown = sorted({row["reason"] for row in out if row["reason"] is not None} - PEDIGREE_CHECK_REASONS)
    if unknown:
        raise ValueError(f"pedigree check reason(s) outside PEDIGREE_CHECK_REASONS: {', '.join(unknown)}")
    return pl.DataFrame(out, schema=PEDIGREE_CHECKS)


def _pedigree_warnings(checks: pl.DataFrame, facet_ids: set[str]) -> list[str]:
    """One warning for the disagreeing pairs a person's facets are in, naming facet ids and bands."""
    pairs = [
        f"{row['facet_a']}/{row['facet_b']} (expected {' or '.join(row['expected_bands'])}, measured {row['band']})"
        for row in checks.filter(pl.col("agreement") == PedigreeAgreement.DISAGREES.value).iter_rows(named=True)
        if row["facet_a"] in facet_ids or row["facet_b"] in facet_ids
    ]
    return [f"{len(pairs)} pair(s) of facets disagree with the supplied pedigree, kept as supplied: {'; '.join(pairs)}"] if pairs else []


def people(
    *loaded: Loaded,
    check: SameSubjectCheck | None = same_person,
    pedigree: Pedigree | None = None,
    coverage: Iterable[Coverage] = (),
) -> dict[str, Person]:
    """Every facet given, grouped into one `Person` per subject, sorted by subject.

    `coverage` is coverage facets, each going to the person its `subject` names. One without a
    subject, or under a subject no facet carries, is refused with `SubjectError`.

    The grouping is the caller's labels and nothing else. A facet without a subject is refused with
    `SubjectError`. Facets of one source under different subjects (one file loaded twice under two
    labels) are allowed and warned in each person they appear in, naming the facet ids.

    With a `pedigree`, its individuals are matched to subjects by id and `check_pedigree` runs (whatever
    `check` is); each person whose facets are in a pair that disagrees with it gets a warning naming
    the facet ids and the bands. Nothing is moved or relabelled.
    """
    pieces = _pieces(loaded)
    rows = [piece.facets.collect().row(0, named=True) for piece in pieces]
    missing = sorted(row["facet_id"] for row in rows if row["subject"] is None)
    if missing:
        raise SubjectError(f"facet(s) without a subject: {', '.join(missing)}; load them with subject=...")
    groups: dict[str, list[Loaded]] = defaultdict(list)
    subjects_of_source: dict[tuple, set[str]] = defaultdict(set)
    for piece, row in zip(pieces, rows):
        groups[row["subject"]].append(piece)
        subjects_of_source[_source(row)].add(row["subject"])
    coverage = tuple(coverage)
    _unlabelled_coverage(coverage)
    stray = sorted(item.facet.facet_id for item in coverage if item.facet.subject not in groups)
    if stray:
        raise SubjectError(f"coverage facet(s) under a subject no facet carries: {', '.join(stray)}")

    result: dict[str, Person] = {}
    for subject in sorted(groups):
        person = Person.of(*groups[subject], check=check, coverage=[item for item in coverage if item.facet.subject == subject])
        shared = _shared_sources(subject, rows, subjects_of_source)
        for warning in shared:
            logger.warning("%s", warning)
        result[subject] = person.model_copy(update={"warnings": person.warnings + tuple(shared)}) if shared else person
    if pedigree is not None:
        checks = check_pedigree(pedigree, *loaded)
        for subject, person in result.items():
            ids = set(person.facets.select("facet_id").collect()["facet_id"].to_list())
            found = _pedigree_warnings(checks, ids)
            for warning in found:
                logger.warning("%s", warning)
            if found:
                result[subject] = person.model_copy(update={"warnings": person.warnings + tuple(found)})
    return result


def _pieces(loaded: tuple[Loaded, ...]) -> list[Loaded]:
    """Each facet of each `Loaded` as a `Loaded` of its own, with the load's links whole: a derived
    facet's lineage runs through links that do not name it (its parent's own steps)."""
    pieces = []
    for item in loaded:
        for facet_id in item.facets.select("facet_id").collect()["facet_id"].to_list():
            pieces.append(
                Loaded(
                    facets=item.facets.filter(pl.col("facet_id") == facet_id),
                    calls=item.calls.filter(pl.col("facet_id") == facet_id),
                    blocks=item.blocks.filter(pl.col("facet_id") == facet_id),
                    links=item.links,
                    warnings=item.warnings,
                    schema_version=item.schema_version,
                )
            )
    return pieces


def _source(facet: dict) -> tuple:
    """What a facet's data was read from: the bytes, the archive member and the sample column. A facet
    made from several (a merge, a consensus) was read from no file, so it is a source of its own."""
    if facet["relation"] in MULTI_PARENT_RELATIONS:
        return ("derived", facet["facet_id"])
    return (facet["source_sha256"], facet["member"], facet["sample"])


def _callset_keys(facets: list[dict], links: pl.DataFrame) -> list[str]:
    """Per facet, in `facets` order, a key shared by every facet of one source and every facet joined to
    it by one-parent lineage (`links.SINGLE_PARENT_RELATIONS`, from `links` and `derived_from`), over
    every id the links name, as `karyotype.person_karyotype` groups a callset. A liftover or a normalize
    shares its parent's source already; the lineage joins what a caller stated across two files (an
    imputed facet and its chip, `imputed_from`)."""
    leader: dict[str, str] = {}

    def find(node: str) -> str:
        leader.setdefault(node, node)
        while leader[node] != node:
            leader[node] = leader[leader[node]]
            node = leader[node]
        return node

    def join(a: str, b: str) -> None:
        first, second = find(a), find(b)
        if first != second:
            leader[max(first, second)] = min(first, second)

    edges = {
        (row["from_facet_id"], row["to_facet_id"])
        for row in links.filter(pl.col("relation").is_in(sorted(SINGLE_PARENT_RELATIONS))).iter_rows(named=True)
    }
    edges |= {(facet["derived_from"], facet["facet_id"]) for facet in facets if facet["derived_from"] is not None}
    for facet in facets:
        join(f"facet:{facet['facet_id']}", f"source:{_source(facet)!r}")
    for a, b in sorted(edges):
        join(f"facet:{a}", f"facet:{b}")
    return [find(f"facet:{facet['facet_id']}") for facet in facets]


def _standing(facets: list[dict], links: pl.DataFrame) -> list[int]:
    """The index of the facet standing for each callset (`_callset_keys`), in `facets` order (which is
    `facet_id` order)."""
    keys = _callset_keys(facets, links)
    standing: dict[str, int] = {}
    for index, facet in enumerate(facets):
        # A facet made from several of them is no evidence of its own: checked against its parents it
        # would agree with what it was computed from.
        if facet["relation"] in MULTI_PARENT_RELATIONS:
            continue
        key = keys[index]
        held = standing.get(key)
        # A facet read from the file stands for its source over one derived from it; facets are in
        # `facet_id` order, so ties go to the lowest id.
        if held is None or (facets[held]["derived_from"] is not None and facet["derived_from"] is None):
            standing[key] = index
    return sorted(standing.values())


def _check(
    pieces: list[Loaded], facets: list[dict], links: pl.DataFrame, check: SameSubjectCheck | None
) -> tuple[pl.DataFrame, list[str]]:
    """The pair checks of one person's callsets, and the warnings they give."""
    stand = _standing(facets, links)

    rows: list[dict] = []
    for a, b in combinations(stand, 2):
        if check is None:
            outcome = SubjectCheckOutcome.NOT_RUN
        else:
            answer = check(pieces[a], pieces[b])
            if answer is not None and not isinstance(answer, bool):
                raise TypeError(f"a subject check must answer True, False or None, not {type(answer).__name__}")
            outcome = (
                SubjectCheckOutcome.COULD_NOT_TELL
                if answer is None
                else SubjectCheckOutcome.SAME
                if answer
                else SubjectCheckOutcome.DIFFERENT
            )
        rows.append({"facet_id": facets[a]["facet_id"], "other_facet_id": facets[b]["facet_id"], "outcome": outcome.value})
    checks = pl.DataFrame(rows, schema=SUBJECT_CHECKS)

    warnings = []
    for outcome, text in (
        (SubjectCheckOutcome.DIFFERENT, "read as two people by the relatedness check, kept under one subject as supplied"),
        (SubjectCheckOutcome.COULD_NOT_TELL, "the relatedness check could not tell apart or together"),
    ):
        pairs = [f"{row['facet_id']}/{row['other_facet_id']}" for row in rows if row["outcome"] == outcome.value]
        if pairs:
            warnings.append(f"{len(pairs)} pair(s) of facets {text}: {', '.join(pairs)}")
    return checks, warnings


def _shared_sources(subject: str, rows: list[dict], subjects_of_source: dict[tuple, set[str]]) -> list[str]:
    """A warning when facets of this subject share their source with facets under another subject."""
    own = sorted({row["facet_id"] for row in rows if row["subject"] == subject and len(subjects_of_source[_source(row)]) > 1})
    if not own:
        return []
    sources = {_source(row) for row in rows if row["facet_id"] in own and row["subject"] == subject}
    others = sorted({row["facet_id"] for row in rows if row["subject"] != subject and _source(row) in sources})
    count = len(set().union(*(subjects_of_source[key] for key in sources)) - {subject})
    if others == own:
        return [f"facet(s) {', '.join(own)} given under {count + 1} subjects"]
    return [
        f"facet(s) {', '.join(own)} share their source (bytes, member, sample) with facet(s) {', '.join(others)} "
        f"under {count} other subject(s)"
    ]
