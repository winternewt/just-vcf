"""The registry as a whole: ids, examples and detection order."""

import re

import pytest

from just_vcf.formats import BCF_FORMAT, FORMATS, REGISTRY, get_spec, validate_format
from just_vcf.formats.spec import HeaderMode
from just_vcf.reference import CHR1_LENGTHS
from just_vcf.sniff import SniffOutcome, sniff

TEXT_SPECS = [spec for spec in REGISTRY if spec.id != BCF_FORMAT]


def test_ids_are_unique_and_match_the_vocabulary() -> None:
    ids = [spec.id for spec in REGISTRY]
    assert len(ids) == len(set(ids))
    assert set(ids) == FORMATS


def test_validate_format_rejects_unregistered_ids() -> None:
    for format_id in sorted(FORMATS):
        assert validate_format(format_id) == format_id
        assert get_spec(format_id).id == format_id
    with pytest.raises(ValueError):
        validate_format("not_a_format")


@pytest.mark.parametrize("spec", TEXT_SPECS, ids=lambda s: s.id)
def test_every_text_spec_carries_an_example(spec) -> None:
    assert spec.example, f"{spec.id} needs a synthetic example head"


@pytest.mark.parametrize("spec", TEXT_SPECS, ids=lambda s: s.id)
def test_example_detects_as_its_own_format_and_no_earlier_one(spec, tmp_path) -> None:
    """The ordering test: first match wins, so no earlier entry may capture a later entry's example."""
    path = tmp_path / "example"
    path.write_text(spec.example)
    result = sniff(path)
    assert result.outcome is SniffOutcome.DETECTED, result.reason
    assert result.format == spec.id
    assert result.family is spec.family


@pytest.mark.parametrize("spec", TEXT_SPECS, ids=lambda s: s.id)
def test_example_rows_have_the_declared_column_count(spec) -> None:
    """A delimited spec's example data lines split into as many fields as it declares columns."""
    if spec.delimiter is None or not spec.columns or spec.header is HeaderMode.VCF:
        pytest.skip("not a delimited table with declared columns")
    lines = [line for line in spec.example.splitlines() if line.strip()]
    if spec.comment_prefix:
        lines = [line for line in lines if not line.startswith(spec.comment_prefix)]
    if spec.header is HeaderMode.LINE:
        lines = lines[1:]
    assert lines, "example has no data lines"
    widths = {len(line.split(spec.delimiter)) for line in lines}
    assert widths == {len(spec.columns)}


@pytest.mark.parametrize("spec", REGISTRY, ids=lambda s: s.id)
def test_build_hints_compile_and_name_a_known_build(spec) -> None:
    for hint in spec.build_hints:
        re.compile(hint.pattern)
        # A label a reference collection can carry (contig 1's length names it), so a stated build can be checked.
        assert hint.build in set(CHR1_LENGTHS.values())
