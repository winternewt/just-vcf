"""RM273: a person (subject) in the schema. The caller supplies the subject; nothing is inferred."""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, rejoin_split
from just_vcf.formats import get_spec
from just_vcf.loading import Loaded
from just_vcf.person import SUBJECT_CHECKS, Person, SubjectError, people
from just_vcf.schema import FACETS

VCF_HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=249250621>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def _array_file(tmp_path: Path, name: str = "genome.txt", format_id: str = "23andme") -> Path:
    path = tmp_path / name
    path.write_text(get_spec(format_id).example)
    return path


def _vcf_file(tmp_path: Path, samples: list[str], name: str = "calls.vcf") -> Path:
    """A small synthetic VCF; each sample's genotype differs so the samples' rows can be told apart."""
    genotypes = ["0/1", "1/1", "0/0", "1|0"]
    records = [
        "\t".join(["1", str(pos), ".", "A", "G", "50", "PASS", ".", "GT", *(genotypes[(i + j) % 4] for j in range(len(samples)))])
        for i, pos in enumerate((100, 200, 300))
    ]
    path = tmp_path / name
    path.write_text("\n".join([*VCF_HEADER, "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]), *records]) + "\n")
    return path


def _subjects(loaded: Loaded) -> list[str | None]:
    return loaded.facets.collect()["subject"].to_list()


# --- the column and `load` ---------------------------------------------------------------------


def test_subject_is_a_facets_column() -> None:
    assert FACETS["subject"] == pl.String


@pytest.mark.parametrize("kind", ["array", "vcf"])
def test_load_writes_the_supplied_subject_and_null_without_one(tmp_path: Path, kind: str) -> None:
    path = _array_file(tmp_path) if kind == "array" else _vcf_file(tmp_path, ["S1"])
    labelled = load(path, subject="person-a")
    plain = load(path)
    assert _subjects(labelled) == ["person-a"]
    assert _subjects(plain) == [None], "nothing is inferred: no subject supplied means none"
    assert labelled.facets.collect().schema == FACETS


def test_the_subject_is_never_read_from_the_sample_name(tmp_path: Path) -> None:
    path = _vcf_file(tmp_path, ["S1"])
    facet = load(path).facets.collect().row(0, named=True)
    assert facet["sample"] == "S1" and facet["subject"] is None


def test_facet_id_does_not_digest_the_subject(tmp_path: Path) -> None:
    """The facet id identifies the data, not the claim about whose it is."""
    for path in (_array_file(tmp_path), _vcf_file(tmp_path, ["S1"])):
        ids = {tuple(load(path, subject=s).facets.collect()["facet_id"].to_list()) for s in (None, "person-a", "person-b")}
        assert len(ids) == 1
    calls_a = load(_array_file(tmp_path), subject="person-a").calls.collect()
    calls_b = load(_array_file(tmp_path), subject="person-b").calls.collect()
    assert calls_a.equals(calls_b), "the subject is a facet column only; the rows do not change"


def test_every_facet_of_a_multi_facet_load_gets_the_subject(tmp_path: Path) -> None:
    """A file loaded as several facets (an appended export, F3) is one claim: every facet carries it."""
    path = tmp_path / "two.txt"
    first = "".join(f"rs{i}\t{1 + i % 22}\t{1000 + i}\t{['AG', 'CC', '--', 'TT'][i % 4]}\n" for i in range(1, 400))
    alleles = ["A", "C", "G", "T", "0"]
    second = "".join(f"rs{i}\t{1 + i % 26}\t{1000 + i}\t{alleles[i % 5]}\t{alleles[(i * 3) % 5]}\n" for i in range(1, 1201))
    path.write_text(first + second)
    loaded = load(path, subject="person-a")
    facets = loaded.facets.collect()
    assert facets.height == 2, "the fixture must load as two facets for this test to mean anything"
    assert _subjects(loaded) == ["person-a", "person-a"]


@pytest.mark.parametrize("bad", ["", "   "])
def test_an_empty_subject_is_refused(tmp_path: Path, bad: str) -> None:
    with pytest.raises(ValueError, match="subject"):
        load(_array_file(tmp_path), subject=bad)


def test_a_derived_facet_inherits_the_subject(tmp_path: Path) -> None:
    parent = load(_array_file(tmp_path), subject="person-a")
    derived = rejoin_split(parent)
    child = derived.facets.collect().row(0, named=True)
    assert child["derived_from"] == parent.facets.collect()["facet_id"].item()
    assert child["subject"] == "person-a"


# --- grouping and the check ----------------------------------------------------------------------


class Recorder:
    """A relatedness check that answers from a table keyed by the pair of facet ids, and records calls."""

    def __init__(self, answers: dict[frozenset[str], bool | None] | None = None, default: bool | None = True) -> None:
        self.answers = answers or {}
        self.default = default
        self.calls: list[tuple[str, str]] = []

    def __call__(self, a: Loaded, b: Loaded) -> bool | None:
        ids = []
        for one in (a, b):
            facets = one.facets.collect()
            assert facets.height == 1, "the check is handed one facet at a time"
            assert one.calls.collect()["facet_id"].unique().to_list() == facets["facet_id"].to_list()
            ids.append(facets["facet_id"].item())
        self.calls.append((ids[0], ids[1]))
        return self.answers.get(frozenset(ids), self.default)


def _facet_id(loaded: Loaded) -> str:
    return loaded.facets.collect()["facet_id"].item()


def _two_sources(tmp_path: Path, subject: str = "person-a") -> tuple[Loaded, Loaded]:
    return load(_array_file(tmp_path), subject=subject), load(_vcf_file(tmp_path, ["S1"]), subject=subject)


def test_a_person_holds_its_facets_rows_in_facet_id_order(tmp_path: Path) -> None:
    array, vcf = _two_sources(tmp_path)
    person = Person.of(vcf, array)
    facets = person.facets.collect()
    assert person.subject == "person-a"
    assert facets["facet_id"].to_list() == sorted([_facet_id(array), _facet_id(vcf)])
    calls = person.calls.collect()
    expected = pl.concat([array.calls.collect(), vcf.calls.collect()])
    assert calls.height == expected.height
    assert calls["facet_id"].to_list() == sorted(calls["facet_id"].to_list()), "facet_id is the first canonical key"
    for one in (array, vcf):
        own = calls.filter(pl.col("facet_id") == _facet_id(one))
        assert own.equals(one.calls.collect()), "a facet's rows come through unchanged and in their order"


def test_without_a_check_every_pair_is_not_run_and_nothing_passes(tmp_path: Path) -> None:
    array, vcf = _two_sources(tmp_path)
    person = Person.of(array, vcf)
    assert person.checks.schema == SUBJECT_CHECKS
    assert person.checks["outcome"].to_list() == ["not_run"]
    assert person.one_person is None, "a check that did not run is not a check that passed"
    assert person.warnings == ()


@pytest.mark.parametrize(
    ("answer", "outcome", "one_person"),
    [(True, "same", True), (False, "different", False), (None, "could_not_tell", None)],
)
def test_the_check_answer_is_kept_kleene(tmp_path: Path, answer: bool | None, outcome: str, one_person: bool | None) -> None:
    array, vcf = _two_sources(tmp_path)
    check = Recorder(default=answer)
    person = Person.of(array, vcf, check=check)
    assert len(check.calls) == 1
    assert person.checks["outcome"].to_list() == [outcome]
    assert person.one_person is one_person
    assert person.facets.collect().height == 2, "never refused, never regrouped, whatever the check says"
    if answer is True:
        assert person.warnings == ()
    else:
        (warning,) = person.warnings
        assert _facet_id(array) in warning and _facet_id(vcf) in warning
        assert "person-a" not in warning, "warnings name facet ids, never the subject"


def test_one_different_pair_makes_the_person_false_and_one_unknown_cannot_make_it_true(tmp_path: Path) -> None:
    array, vcf = _two_sources(tmp_path)
    other = load(_array_file(tmp_path, "other.txt", "23andme_csv"), subject="person-a")
    ids = [_facet_id(x) for x in (array, vcf, other)]
    unknown = Person.of(array, vcf, other, check=Recorder({frozenset(ids[:2]): None}))
    assert sorted(unknown.checks["outcome"].to_list()) == ["could_not_tell", "same", "same"]
    assert unknown.one_person is None
    different = Person.of(array, vcf, other, check=Recorder({frozenset(ids[:2]): None, frozenset(ids[1:]): False}))
    assert different.one_person is False
    assert len(different.warnings) == 2, "one line per reason"


def test_facets_of_one_source_are_not_checked_against_each_other(tmp_path: Path) -> None:
    """A parent and the facet derived from it are one person's data by construction."""
    array, vcf = _two_sources(tmp_path)
    derived = rejoin_split(array)
    check = Recorder()
    person = Person.of(array, derived, vcf, check=check)
    assert person.facets.collect().height == 3
    assert len(check.calls) == 1
    assert {frozenset(pair) for pair in check.calls} == {frozenset({_facet_id(array), _facet_id(vcf)})}, (
        "the source stands as its facet read from the file, not the derived one"
    )
    alone = Person.of(array, derived, check=check)
    assert alone.checks.height == 0 and alone.one_person is None


def test_samples_of_one_file_are_different_sources(tmp_path: Path) -> None:
    """A trio VCF: one file, three people. Same bytes, different sample columns."""
    path = _vcf_file(tmp_path, ["S1", "S2"])
    first, second = (load(path, sample=s, subject="person-a") for s in ("S1", "S2"))
    check = Recorder(default=False)
    person = Person.of(first, second, check=check)
    assert len(check.calls) == 1 and person.one_person is False


def test_grouping_refuses_what_it_cannot_group(tmp_path: Path) -> None:
    array, vcf = _two_sources(tmp_path)
    unlabelled = load(_array_file(tmp_path, "other.txt", "23andme_csv"))
    with pytest.raises(SubjectError, match=_facet_id(unlabelled)):
        Person.of(array, unlabelled)
    with pytest.raises(SubjectError, match=_facet_id(unlabelled)):
        people(array, unlabelled)
    elsewhere = load(_vcf_file(tmp_path, ["S1"], "b.vcf"), subject="person-b")
    with pytest.raises(SubjectError, match="2 different subjects") as refused:
        Person.of(array, elsewhere)
    assert "person-a" not in str(refused.value) and "person-b" not in str(refused.value)
    with pytest.raises(ValueError, match="more than once"):
        Person.of(array, load(_array_file(tmp_path), subject="person-a"))
    with pytest.raises(SubjectError):
        Person.of()


def test_people_groups_by_the_supplied_label_only(tmp_path: Path) -> None:
    a_array, a_vcf = _two_sources(tmp_path, "person-a")
    b = load(_array_file(tmp_path, "b.txt", "23andme_csv"), subject="person-b")
    check = Recorder(default=False)
    grouped = people(b, a_vcf, a_array, check=check)
    assert list(grouped) == ["person-a", "person-b"]
    assert set(grouped["person-a"].facets.collect()["facet_id"]) == {_facet_id(a_array), _facet_id(a_vcf)}
    assert grouped["person-b"].facets.collect()["facet_id"].to_list() == [_facet_id(b)]
    assert grouped["person-a"].one_person is False, "the check said two people, and the grouping stands"
    assert {frozenset(pair) for pair in check.calls} == {frozenset({_facet_id(a_array), _facet_id(a_vcf)})}, (
        "people never checks across subjects"
    )


def test_one_file_under_two_subjects_is_allowed_and_warned(tmp_path: Path) -> None:
    path = _array_file(tmp_path)
    as_a, as_b = load(path, subject="person-a"), load(path, subject="person-b")
    unrelated = load(_vcf_file(tmp_path, ["S1"]), subject="person-c")
    grouped = people(as_a, as_b, unrelated)
    facet_id = _facet_id(as_a)
    assert _facet_id(as_b) == facet_id
    for subject in ("person-a", "person-b"):
        (warning,) = grouped[subject].warnings
        assert facet_id in warning and "1 other subject" in warning
        assert "person-a" not in warning and "person-b" not in warning
    assert grouped["person-c"].warnings == ()


def test_two_samples_of_one_file_under_two_subjects_are_not_warned(tmp_path: Path) -> None:
    path = _vcf_file(tmp_path, ["S1", "S2"])
    grouped = people(load(path, sample="S1", subject="child"), load(path, sample="S2", subject="mother"))
    assert all(person.warnings == () for person in grouped.values())


def test_a_check_must_answer_kleene(tmp_path: Path) -> None:
    array, vcf = _two_sources(tmp_path)
    with pytest.raises(TypeError, match="True, False or None"):
        Person.of(array, vcf, check=lambda a, b: 0.9)  # type: ignore[arg-type, return-value]
