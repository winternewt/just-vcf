"""Facets of one person, grouped by the subject the caller supplied (RM273, docs/SCHEMA.md).

The caller says whose a file is with `load(..., subject=...)`. This module groups facets by that label
and checks the claim; it never infers a subject and never moves a facet to another one. The check is a
relatedness test over pairs of facets, passed in as a callable (`SameSubjectCheck`): `True` says one
person, `False` two, `None` that it could not tell. A `False` is a warning naming both facets, and the
facets stay where the caller put them. The default check is `just_vcf.relatedness.same_person` (RM274),
KING-robust kinship over the package's site panel; `check=None` runs none, and every pair is then
`not_run`, which is not a pass.
"""

import logging
from collections import defaultdict
from collections.abc import Callable
from enum import Enum
from itertools import combinations

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.loading import Loaded
from just_vcf.relatedness import same_person
from just_vcf.schema import SCHEMA_VERSION, enum_dtype

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
    canonical order (SCHEMA.md Decision 11).
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    subject: str
    """The caller's label. It identifies a person if the caller makes it so; nothing here logs it."""
    facets: pl.LazyFrame
    calls: pl.LazyFrame
    blocks: pl.LazyFrame
    checks: pl.DataFrame
    """`SUBJECT_CHECKS`: each pair of the person's sources and what the check said about it."""
    warnings: tuple[str, ...] = ()
    """One line per reason, naming facet ids, never the subject."""
    schema_version: str = SCHEMA_VERSION

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
    def of(cls, *loaded: Loaded, check: SameSubjectCheck | None = same_person) -> "Person":
        """One person from facets that all carry one subject.

        A facet without a subject, or facets under two subjects, are refused with `SubjectError`; so is
        one `facet_id` given twice, whose rows would then count twice. Facets of one source (the same
        bytes, archive member and sample: a file and the facets derived from it, or the facets one file
        loads as) are one person by construction and are not checked against each other. Every other
        pair of sources is checked once, each source standing as its facet read from the file (else
        its lowest `facet_id`).
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

        order = sorted(range(len(pieces)), key=lambda i: ids[i])
        pieces = [pieces[i] for i in order]
        facets = [facets[i] for i in order]
        checks, warnings = _check(pieces, facets, check)
        for warning in warnings:
            logger.warning("%s", warning)
        return cls(
            subject=facets[0]["subject"],
            facets=pl.concat([piece.facets for piece in pieces], how="vertical"),
            calls=pl.concat([piece.calls for piece in pieces], how="vertical"),
            blocks=pl.concat([piece.blocks for piece in pieces], how="vertical"),
            checks=checks,
            warnings=tuple(warnings),
            schema_version=versions[0],
        )


def people(*loaded: Loaded, check: SameSubjectCheck | None = same_person) -> dict[str, Person]:
    """Every facet given, grouped into one `Person` per subject, sorted by subject.

    The grouping is the caller's labels and nothing else. A facet without a subject is refused with
    `SubjectError`. Facets of one source under different subjects (one file loaded twice under two
    labels) are allowed and warned in each person they appear in, naming the facet ids.
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

    result: dict[str, Person] = {}
    for subject in sorted(groups):
        person = Person.of(*groups[subject], check=check)
        shared = _shared_sources(subject, rows, subjects_of_source)
        for warning in shared:
            logger.warning("%s", warning)
        result[subject] = person.model_copy(update={"warnings": person.warnings + tuple(shared)}) if shared else person
    return result


def _pieces(loaded: tuple[Loaded, ...]) -> list[Loaded]:
    """Each facet of each `Loaded` as a `Loaded` of its own."""
    pieces = []
    for item in loaded:
        for facet_id in item.facets.select("facet_id").collect()["facet_id"].to_list():
            pieces.append(
                Loaded(
                    facets=item.facets.filter(pl.col("facet_id") == facet_id),
                    calls=item.calls.filter(pl.col("facet_id") == facet_id),
                    blocks=item.blocks.filter(pl.col("facet_id") == facet_id),
                    warnings=item.warnings,
                    schema_version=item.schema_version,
                )
            )
    return pieces


def _source(facet: dict) -> tuple:
    """What a facet's data was read from: the bytes, the archive member and the sample column."""
    return (facet["source_sha256"], facet["member"], facet["sample"])


def _check(pieces: list[Loaded], facets: list[dict], check: SameSubjectCheck | None) -> tuple[pl.DataFrame, list[str]]:
    """The pair checks of one person's sources, and the warnings they give."""
    standing: dict[tuple, int] = {}
    for index, facet in enumerate(facets):
        key = _source(facet)
        held = standing.get(key)
        # A facet read from the file stands for its source over one derived from it; facets are in
        # `facet_id` order, so ties go to the lowest id.
        if held is None or (facets[held]["derived_from"] is not None and facet["derived_from"] is None):
            standing[key] = index
    stand = sorted(standing.values())

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
