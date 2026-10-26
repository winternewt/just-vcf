"""GA4GH Pedigree, read as data the caller supplies (RM274). Nothing here infers a relationship.

The GA4GH Pedigree standard (https://pedigree.readthedocs.io) has three parts: a `Pedigree` with an
`id`, its `individuals` and its `relationships`; an `Individual` with an `id` and a `sex`; and a
`Relationship` saying what its `individual` is to its `relative`, the `relation` a Kinship Ontology
(KIN) term such as `KIN:027` (`isBiologicalMotherOf`). Its documented serialisation is the Phenopacket
`GA4GHPedigree` message, where an individual's `id` and `sex` sit under `subject` and a relationship
names `individual_id` and `relative_id`. `read_pedigree` reads that message as JSON, with protobuf's two JSON spellings of a field (`individual_id`, `individualId`) both accepted, into
two frames. Terms are kept verbatim (charter invariant 2): the KIN id, its label as written, the sex
value as written.

A pedigree states relationships; `just_vcf.relatedness` measures how alike genotypes are.
`expected_kinship` reads what a pedigree says each pair's kinship should be, and
`just_vcf.person.check_pedigree` sets that beside the measured band. A stated parent whose facets read
unrelated is a finding, never a correction: nothing here changes a pedigree (charter invariant 8).

**Expected kinship** (just-vcf's own reading of KIN, charter invariant 2; KIN has no kinship values).
On KING's scale, one person 1/2, the kinship coefficient of the pair:

- *Parent links* are KIN:003 (`isBiologicalParentOf`) and the terms under it (KIN:004 sperm donor,
  KIN:006 surrogate ovum donor, KIN:027 mother, KIN:028 father, KIN:038 ovum donor), and KIN:032
  (`isBiologicalChildOf`) read the other way. Through them every pair's kinship is the pedigree
  coefficient, `phi(x, y) = (phi(x, p1) + phi(x, p2)) / 2` over the later one's parents, a person's own
  `(1 + phi(mother, father)) / 2`. A parent the pedigree does not list is a founder. This is what
  makes grandparents, siblings, aunts, uncles and cousins follow from parent links alone. Basis `derived`.
- *Stated degrees* hold where a pair is related by a term that fixes its kinship, and they carry into
  the recursion, so a stated sibling's child is a stated nephew: full sibling (KIN:008) and polyzygotic
  twin (KIN:011) 1/4; monozygotic twin (KIN:010) 1/2; a multiple-birth sibling of unknown zygosity
  (KIN:009) 1/4 or 1/2; a biological sibling (KIN:007) 1/4 or 1/8; half sibling (KIN:012, 054, 055)
  1/8; grandparent (KIN:017, 036, 052, 053) 1/8; parental sibling (KIN:013, 046, 058 to 061, a sibling
  of a parent, full or half) 1/8 or 1/16; cousin (KIN:014 to 016, a child of a parental sibling) 1/16
  or 1/32; great-grandparent (KIN:018, 047) 1/16. Where a value is a set, every band in it is
  expected. A stated term takes precedence over the parent links for its pair. KIN files twins under
  full siblings; a full sibling stated as KIN:008 is read as 1/4, so a twin reads as disagreeing
  unless it is stated as one. Basis `stated`.
- *Terms with no degree* (KIN:001 relative, KIN:002 biological relative, KIN:030 and 051
  consanguineous partner) make the pair's kinship unknown (`term_without_degree`), and every pair
  reached through it. Social, legal and gestational terms (KIN:005, 019 to 026, 029, 050) say nothing
  about genotypes and are not read.
- *No path*: founders of one pedigree are taken as unrelated, the convention every pedigree kinship
  computation rests on. Such a pair expects `unrelated`, basis `assumed`, so the pedigree's silence is
  visible as an assumption.

*First degree, by kind.* A pair whose allowed kinship falls in the first degree is expected as one of
`relatedness`'s two first-degree bands where the pedigree says which: `parent_child` where one is a
listed biological parent of the other (a parent and child share an allele at every site, inbred or
not); `full_siblings` where both have the same two listed parents, or a stated term allows the pair a
kinship of 1/4 (every such term is a sibling term: KIN:007, 008, 009, 011). A stated term takes
precedence here too. A first-degree kinship reached any other way (a twin's child and the other twin)
is `first_degree`, either kind. just-vcf's own reading, as the rest of this list is.

A child listed with more than two biological parents makes every pair reached through it unknown
(`more_than_two_parents`); two statements about one pair that share no value make it
`statements_conflict`. A pedigree whose parent links loop is refused (`PedigreeError`).
"""

import json
from enum import Enum
from fractions import Fraction
from functools import cache
from itertools import combinations
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.relatedness import RelationBand, band_of, finest_bands

PEDIGREE_INDIVIDUALS = pl.Schema({"pedigree_id": pl.String, "individual_id": pl.String, "sex": pl.String})
"""One row per individual: the `subject.id` and `subject.sex` as written."""

PEDIGREE_RELATIONSHIPS = pl.Schema(
    {
        "pedigree_id": pl.String,
        "individual_id": pl.String,
        "relation_id": pl.String,
        "relation_label": pl.String,
        "relative_id": pl.String,
    }
)
"""One row per relationship: `individual_id` is `relation` of `relative_id` (KIN's direction)."""


class PedigreeError(ValueError):
    """A pedigree document that does not have the standard's required parts."""


class Pedigree(BaseModel):
    """A supplied GA4GH Pedigree as two frames, with its own fields kept as written."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    pedigree_id: str
    individuals: pl.DataFrame
    relationships: pl.DataFrame
    index_patients: tuple[str, ...] = ()
    narrative: str | None = None
    date: str | None = None


def _field(item: dict, snake: str) -> object:
    head, *rest = snake.split("_")
    camel = head + "".join(part.title() for part in rest)
    return item.get(snake, item.get(camel))


def pedigree_from_dict(document: dict) -> Pedigree:
    """A `GA4GHPedigree` message, as parsed JSON, read into a `Pedigree`.

    Refused with `PedigreeError`: no `id`; an individual with no `subject.id` (or `id` where there is no
    `subject`); a relationship missing `individual_id`, `relation.id` or `relative_id`; a relationship
    naming an individual the pedigree does not list; two individuals with one id.
    """
    pedigree_id = document.get("id")
    if not pedigree_id:
        raise PedigreeError("a pedigree needs an id")
    people = []
    for item in document.get("individuals") or []:
        subject = item.get("subject", item)
        if not subject.get("id"):
            raise PedigreeError(f"pedigree {pedigree_id}: an individual without an id")
        sex = subject.get("sex")
        people.append({"pedigree_id": str(pedigree_id), "individual_id": str(subject["id"]), "sex": None if sex is None else str(sex)})
    individuals = pl.DataFrame(people, schema=PEDIGREE_INDIVIDUALS)
    if individuals["individual_id"].is_duplicated().any():
        raise PedigreeError(f"pedigree {pedigree_id}: two individuals share an id")
    known = set(individuals["individual_id"].to_list())
    links = []
    for item in document.get("relationships") or []:
        individual, relative, relation = _field(item, "individual_id"), _field(item, "relative_id"), item.get("relation") or {}
        if not individual or not relative or not relation.get("id"):
            raise PedigreeError(f"pedigree {pedigree_id}: a relationship needs individual_id, relation.id and relative_id")
        unknown = {str(individual), str(relative)} - known
        if unknown:
            raise PedigreeError(f"pedigree {pedigree_id}: a relationship names {', '.join(sorted(unknown))}, not among its individuals")
        label = relation.get("label")
        links.append(
            {
                "pedigree_id": str(pedigree_id),
                "individual_id": str(individual),
                "relation_id": str(relation["id"]),
                "relation_label": None if label is None else str(label),
                "relative_id": str(relative),
            }
        )
    relationships = pl.DataFrame(links, schema=PEDIGREE_RELATIONSHIPS)
    index = tuple(str(value) for value in (_field(document, "index_patients") or []))
    date = document.get("date")
    return Pedigree(
        pedigree_id=str(pedigree_id),
        individuals=individuals,
        relationships=relationships,
        index_patients=index,
        narrative=document.get("narrative"),
        date=None if date is None else str(date),
    )


def read_pedigree(path: Path | str) -> Pedigree:
    """A `GA4GHPedigree` message in a JSON file, read into a `Pedigree` (`pedigree_from_dict`)."""
    return pedigree_from_dict(json.loads(Path(path).read_text()))


PARENT_TERMS: frozenset[str] = frozenset({"KIN:003", "KIN:004", "KIN:006", "KIN:027", "KIN:028", "KIN:038"})
"""`individual_id` is a biological parent of `relative_id`."""
CHILD_TERMS: frozenset[str] = frozenset({"KIN:032"})
"""`individual_id` is a biological child of `relative_id`."""

_HALF, _QUARTER, _EIGHTH, _SIXTEENTH, _THIRTY_SECOND = (Fraction(1, n) for n in (2, 4, 8, 16, 32))
_ZERO = frozenset({Fraction(0)})
STATED_KINSHIP: dict[str, frozenset[Fraction]] = {
    **dict.fromkeys(("KIN:008", "KIN:011"), frozenset({_QUARTER})),
    "KIN:010": frozenset({_HALF}),
    "KIN:009": frozenset({_QUARTER, _HALF}),
    "KIN:007": frozenset({_QUARTER, _EIGHTH}),
    **dict.fromkeys(("KIN:012", "KIN:054", "KIN:055"), frozenset({_EIGHTH})),
    **dict.fromkeys(("KIN:017", "KIN:036", "KIN:052", "KIN:053"), frozenset({_EIGHTH})),
    **dict.fromkeys(("KIN:013", "KIN:046", "KIN:058", "KIN:059", "KIN:060", "KIN:061"), frozenset({_EIGHTH, _SIXTEENTH})),
    **dict.fromkeys(("KIN:014", "KIN:015", "KIN:016"), frozenset({_SIXTEENTH, _THIRTY_SECOND})),
    **dict.fromkeys(("KIN:018", "KIN:047"), frozenset({_SIXTEENTH})),
}
"""The kinships a KIN term allows for its pair (module docstring). just-vcf's own reading."""
NO_DEGREE_TERMS: frozenset[str] = frozenset({"KIN:001", "KIN:002", "KIN:030", "KIN:051"})
"""Terms that say a pair is related without saying how closely."""

TERM_WITHOUT_DEGREE = "term_without_degree"
"""A term relating the pair, or a pair on the way, says it is related but not how closely."""
MORE_THAN_TWO_PARENTS = "more_than_two_parents"
"""A child on the way has more than two biological parents listed."""
STATEMENTS_CONFLICT = "statements_conflict"
"""Two statements about one pair allow no common kinship."""
PEDIGREE_REASONS: frozenset[str] = frozenset({TERM_WITHOUT_DEGREE, MORE_THAN_TWO_PARENTS, STATEMENTS_CONFLICT})
"""Why a pedigree gives a pair no expected kinship. Grows additively."""

BAND_ORDER: tuple[RelationBand, ...] = (
    RelationBand.UNRELATED,
    RelationBand.THIRD_DEGREE,
    RelationBand.SECOND_DEGREE,
    RelationBand.FIRST_DEGREE,
    RelationBand.SAME_PERSON,
)
"""The named bands of the kinship ladder, from least to most related. `parent_child` and
`full_siblings` stand at `first_degree`'s rung (`relatedness.finest_bands`)."""


def _rung(band: RelationBand) -> int:
    """`band`'s place on `BAND_ORDER`, a first-degree kind at the first degree's."""
    return BAND_ORDER.index(RelationBand.FIRST_DEGREE if band in finest_bands(RelationBand.FIRST_DEGREE) else band)


class ExpectedBasis(str, Enum):
    """Where a pair's expected kinship comes from."""

    STATED = "stated"
    """A term relating the two fixes it."""
    DERIVED = "derived"
    """Computed through parent links (and stated degrees further up)."""
    ASSUMED = "assumed"
    """No biological path: founders of one pedigree are taken as unrelated."""


EXPECTED_KINSHIP = pl.Schema(
    {
        "pedigree_id": pl.String,
        "individual_id": pl.String,
        "relative_id": pl.String,
        "relation_id": pl.String,
        "relation_label": pl.String,
        "basis": pl.Enum([basis.value for basis in ExpectedBasis]),
        "expected_kinship": pl.Float64,
        "expected_bands": pl.List(pl.Enum([band.value for band in RelationBand])),
        "expected_reason": pl.String,
    }
)
"""`expected_kinship`: one row per pair of individuals, in the pedigree's order (`individual_id` the
earlier). `relation_id`/`relation_label` are the first statement relating the two, verbatim, null when
none does. `expected_kinship` is null unless one value is allowed; `expected_bands` is every band the
allowed values fall in (least related first), null with `expected_reason` when unknown."""

_Kin = frozenset[Fraction] | str
"""The kinships a pair may have, or why it is unknown."""


def _half_sum(first: _Kin, second: _Kin) -> _Kin:
    if isinstance(first, str):
        return first
    if isinstance(second, str):
        return second
    return frozenset((a + b) / 2 for a in first for b in second)


def _parents(pedigree: Pedigree) -> dict[str, list[str]]:
    parents: dict[str, list[str]] = {person: [] for person in pedigree.individuals["individual_id"].to_list()}
    for individual, relation, relative in pedigree.relationships.select("individual_id", "relation_id", "relative_id").iter_rows():
        if relation in PARENT_TERMS:
            parent, child = individual, relative
        elif relation in CHILD_TERMS:
            parent, child = relative, individual
        else:
            continue
        if parent not in parents[child]:
            parents[child].append(parent)
    return parents


def _order(pedigree_id: str, parents: dict[str, list[str]]) -> dict[str, int]:
    """Each individual's place, every parent before its children; parent links that loop are refused."""
    order: dict[str, int] = {}
    visiting: set[str] = set()

    def visit(person: str) -> None:
        if person in order:
            return
        if person in visiting:
            raise PedigreeError(f"pedigree {pedigree_id}: its parent links loop")
        visiting.add(person)
        for parent in parents[person]:
            visit(parent)
        visiting.discard(person)
        order[person] = len(order)

    for person in parents:
        visit(person)
    return order


def _stated(pedigree: Pedigree) -> dict[frozenset[str], _Kin]:
    """The kinships stated terms allow, per pair: the common values of every statement about it."""
    stated: dict[frozenset[str], _Kin] = {}
    for individual, relation, relative in pedigree.relationships.select("individual_id", "relation_id", "relative_id").iter_rows():
        if relation in STATED_KINSHIP:
            value: _Kin = STATED_KINSHIP[relation]
        elif relation in NO_DEGREE_TERMS:
            value = TERM_WITHOUT_DEGREE
        else:
            continue
        key = frozenset((individual, relative))
        held = stated.get(key)
        if held is None or (isinstance(value, str) and not isinstance(held, str)):
            stated[key] = value
        elif not isinstance(held, str) and not isinstance(value, str):
            stated[key] = (held & value) or STATEMENTS_CONFLICT
    return stated


def expected_kinship(pedigree: Pedigree) -> pl.DataFrame:
    """What `pedigree` says each pair of its individuals' kinship should be (module docstring), as
    `EXPECTED_KINSHIP`. A pedigree whose parent links loop is refused with `PedigreeError`."""
    parents = _parents(pedigree)
    order = _order(pedigree.pedigree_id, parents)
    stated = _stated(pedigree)

    @cache
    def phi(x: str, y: str) -> _Kin:
        if order[x] > order[y]:
            x, y = y, x
        if x == y:
            listed = parents[x]
            if len(listed) > 2:
                return MORE_THAN_TWO_PARENTS
            inbreeding = phi(*listed) if len(listed) == 2 else _ZERO
            return inbreeding if isinstance(inbreeding, str) else frozenset((1 + value) / 2 for value in inbreeding)
        if (key := frozenset((x, y))) in stated:
            return stated[key]
        # y comes later, so it is never an ancestor of x: expand y through its parents, or x when y is
        # a founder. A parent the pedigree does not list is a founder, unrelated to everyone.
        for person, other in ((y, x), (x, y)):
            listed = parents[person]
            if len(listed) > 2:
                return MORE_THAN_TWO_PARENTS
            if listed:
                halves = [phi(parent, other) for parent in listed] + [_ZERO] * (2 - len(listed))
                return _half_sum(*halves)
        return _ZERO

    first_statement: dict[frozenset[str], tuple[str, str | None]] = {}
    for individual, relation, label, relative in pedigree.relationships.select(
        "individual_id", "relation_id", "relation_label", "relative_id"
    ).iter_rows():
        first_statement.setdefault(frozenset((individual, relative)), (relation, label))

    def first_degree_kind(x: str, y: str, key: frozenset[str]) -> RelationBand:
        """Which first-degree band the pedigree expects of the pair, or `first_degree` when it does not say."""
        if key in stated:
            return RelationBand.FULL_SIBLINGS  # every term allowing 1/4 is a sibling term
        if x in parents[y] or y in parents[x]:
            return RelationBand.PARENT_CHILD
        if len(parents[x]) == 2 and set(parents[x]) == set(parents[y]):
            return RelationBand.FULL_SIBLINGS
        return RelationBand.FIRST_DEGREE

    rows = []
    for individual, relative in combinations(pedigree.individuals["individual_id"].to_list(), 2):
        key = frozenset((individual, relative))
        value = phi(individual, relative)
        relation_id, relation_label = first_statement.get(key, (None, None))
        kinship, bands, reason = None, None, None
        if isinstance(value, str):
            basis, reason = (ExpectedBasis.STATED if key in stated else ExpectedBasis.DERIVED), value
        else:
            basis = ExpectedBasis.STATED if key in stated else ExpectedBasis.ASSUMED if value == _ZERO else ExpectedBasis.DERIVED
            kinship = float(next(iter(value))) if len(value) == 1 else None
            found = {band_of(float(v)) for v in value}
            if RelationBand.FIRST_DEGREE in found:
                found = (found - {RelationBand.FIRST_DEGREE}) | {first_degree_kind(individual, relative, key)}
            bands = [band.value for band in sorted(found, key=_rung)]
        rows.append(
            {
                "pedigree_id": pedigree.pedigree_id,
                "individual_id": individual,
                "relative_id": relative,
                "relation_id": relation_id,
                "relation_label": relation_label,
                "basis": basis.value,
                "expected_kinship": kinship,
                "expected_bands": bands,
                "expected_reason": reason,
            }
        )
    return pl.DataFrame(rows, schema=EXPECTED_KINSHIP)
