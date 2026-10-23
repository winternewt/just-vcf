"""docs/FORMAT_REGISTRY.md against the code it is generated from (scripts/format_registry.py)."""

import importlib.util
import re
import sys
from pathlib import Path

import pytest

from just_vcf import LOADABLE_FORMATS
from just_vcf.formats import FORMATS, REGISTRY
from just_vcf.formats.spec import Coordinates, Family, HeaderMode
from just_vcf.sniff import SniffOutcome, sniff

_ROOT = Path(__file__).resolve().parent.parent
_PATH = _ROOT / "scripts" / "format_registry.py"
_spec = importlib.util.spec_from_file_location("format_registry", _PATH)
registry_doc = importlib.util.module_from_spec(_spec)
sys.modules["format_registry"] = registry_doc
_spec.loader.exec_module(registry_doc)

SHAPES_WITH_EXPECTATION = [shape for shape in registry_doc.SHAPES if shape.expect is not None]


@pytest.fixture(scope="module")
def rendered() -> str:
    return registry_doc.render()


def test_the_committed_doc_is_what_the_code_renders(rendered: str) -> None:
    """Fails when a spec, a sniff rule, a loader or a shape changed and the doc was not regenerated."""
    committed = registry_doc.DOC.read_text(encoding="utf-8")
    assert committed == rendered, "docs/FORMAT_REGISTRY.md is stale: run uv run python scripts/format_registry.py"


def test_every_registered_format_has_its_section_in_order(rendered: str) -> None:
    headings = re.findall(r"^### `([a-z0-9_]+)`$", rendered, flags=re.MULTILINE)
    assert headings == [spec.id for spec in REGISTRY]
    assert set(headings) == FORMATS


def test_every_format_names_whether_it_loads(rendered: str) -> None:
    sections = re.split(r"^### `", rendered, flags=re.MULTILINE)[1:]
    for section in sections:
        format_id = section.split("`", 1)[0]
        if format_id not in FORMATS:
            continue
        loads = "- **Loads:** yes, through `just_vcf.loaders." in section
        assert loads == (format_id in LOADABLE_FORMATS), format_id
    loaders = {format_id for _, formats in registry_doc.LOADERS for format_id in formats}
    assert loaders == LOADABLE_FORMATS


def test_every_vocabulary_member_is_described() -> None:
    assert set(registry_doc.FAMILY_MEANING) == set(Family)
    assert set(registry_doc.COORDINATES_MEANING) == set(Coordinates)
    assert set(registry_doc.HEADER_MEANING) == set(HeaderMode)


def test_every_finding_f1_to_f17_is_accounted_for(rendered: str) -> None:
    findings = {f for shape in registry_doc.SHAPES for f in shape.findings} | {f for f, _ in registry_doc.NOT_SHAPES}
    assert findings == {f"F{n}" for n in range(1, 18)}
    for n in range(1, 18):
        assert re.search(rf"\bF{n}\b", rendered), f"F{n}"


def test_shapes_have_unique_keys_a_known_status_and_an_example_or_a_reason() -> None:
    keys = [shape.key for shape in registry_doc.SHAPES]
    assert len(keys) == len(set(keys))
    for shape in registry_doc.SHAPES:
        assert shape.status in registry_doc.STATUSES, shape.key
        assert shape.example is not None or shape.members or shape.example_note, shape.key


@pytest.mark.parametrize("shape", SHAPES_WITH_EXPECTATION, ids=lambda s: s.key)
def test_shape_example_sniffs_as_the_doc_says(shape, tmp_path: Path) -> None:
    expect = shape.expect
    result = sniff(registry_doc.materialize(shape, tmp_path))
    assert result.outcome is SniffOutcome(expect.outcome), result.reason
    assert result.format == expect.format
    if expect.check_build:
        assert result.build == expect.build
    if expect.layers:
        assert result.layers == expect.layers
    if expect.reason_contains:
        assert expect.reason_contains in (result.reason or "")


def test_every_cited_test_exists() -> None:
    defined = {
        name
        for path in (_ROOT / "tests").rglob("test_*.py")
        for name in re.findall(r"^def (test_\w+)\(", path.read_text(encoding="utf-8"), flags=re.MULTILINE)
    }
    cited = {test for shape in registry_doc.SHAPES for test in shape.tests}
    assert cited - defined == set()
