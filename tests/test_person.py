"""RM273: a person (subject) in the schema. The caller supplies the subject; nothing is inferred."""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, rejoin_split
from just_vcf.formats import get_spec
from just_vcf.loading import Loaded
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
