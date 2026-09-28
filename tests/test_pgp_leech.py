"""Tests for scripts/pgp_leech.py. Live tests read the public listing and skip when offline."""

import importlib.util
import re
import sys
from pathlib import Path

import httpx
import pytest

_PATH = Path(__file__).resolve().parent.parent / "scripts" / "pgp_leech.py"
_spec = importlib.util.spec_from_file_location("pgp_leech", _PATH)
leech = importlib.util.module_from_spec(_spec)
sys.modules["pgp_leech"] = leech
_spec.loader.exec_module(leech)

# A real row of the public 23andMe listing (download 4208), as the site served it on 2026-09-28.
ROW_4208 = leech.ListingRow(
    download_id=4208,
    participant="huC8B936",
    pgp_alias=None,
    participant_numeric_id=17756,
    published="2026-04-14",
    data_type="23andMe",
    source="Participant",
    name="James",
    locator="745d71146d59a622dc9f936edf97db77+99",
    listed_size=5906047,
)


@pytest.fixture(scope="module")
def listing_html() -> str:
    try:
        resp = httpx.get(
            f"{leech.PGP_BASE}{leech.LISTING_PATH}?data_type=23andMe",
            headers={"User-Agent": leech.USER_AGENT},
            timeout=60,
        )
    except httpx.TransportError as err:
        pytest.skip(f"PGP-HMS unreachable: {err}")
    resp.raise_for_status()
    return resp.text


def test_listing_parses_every_file_row(listing_html: str) -> None:
    rows, dropped = leech.parse_listing(listing_html)
    # Independent count of the row marker in the raw page: every marked row is parsed or accounted for.
    marked = len(re.findall(r"<tr data-file-row", listing_html))
    assert len(rows) + sum(dropped.values()) == marked
    assert not dropped
    ids = [r.download_id for r in rows]
    assert ids == sorted(ids) and len(set(ids)) == len(ids)
    assert all(re.fullmatch(r"\d{4}-\d{2}-\d{2}", r.published) for r in rows)
    assert all(re.fullmatch(r"hu[0-9A-F]+", r.participant) for r in rows)
    aliased = [r for r in rows if r.pgp_alias is not None]
    assert aliased, "the listing has stopped showing PGP enrolment aliases; re-check the participant cell"
    assert all(re.fullmatch(r"PGP\d+", r.pgp_alias) for r in aliased)
    # One participant keeps one alias across all their files.
    assert all(len({r.pgp_alias for r in aliased if r.participant == p}) == 1 for p in {r.participant for r in aliased})
    assert all(r.listed_size > 0 for r in rows)
    assert {r.data_type for r in rows} == {"23andMe"}


def test_listing_row_matches_recorded_values(listing_html: str) -> None:
    rows, _ = leech.parse_listing(listing_html)
    by_id = {r.download_id: r for r in rows}
    if 4208 not in by_id:
        pytest.skip("download 4208 is no longer listed")
    assert by_id[4208] == ROW_4208


def test_collection_listing_size_matches_row() -> None:
    fetcher = leech.Fetcher(delay_s=0)
    try:
        with fetcher.open(f"{leech.PGP_BASE}/user_file/download/{ROW_4208.download_id}", follow=False) as resp:
            location = resp.headers.get("location")
        if location is None:
            pytest.skip("download 4208 is no longer served through Arvados")
        page = fetcher.get(location)
    except httpx.TransportError as err:
        pytest.skip(f"PGP-HMS or Arvados unreachable: {err}")
    files = leech.parse_collection(page.text, str(page.url))
    assert sum(f.size for f in files) == ROW_4208.listed_size
    assert all(f.url.startswith(str(page.url).rstrip("/")) for f in files)


def test_already_present_needs_the_full_listed_size(tmp_path: Path) -> None:
    target = tmp_path / "23andMe" / ROW_4208.participant / str(ROW_4208.download_id)
    target.mkdir(parents=True)
    assert leech._already_present(ROW_4208, target) is None
    (target / "export.zip").write_bytes(b"\0" * (ROW_4208.listed_size - 1))
    assert leech._already_present(ROW_4208, target) is None
    (target / "export.zip").write_bytes(b"\0" * ROW_4208.listed_size)
    (target / "other.zip.part").write_bytes(b"\0" * 10)  # an interrupted file never counts
    present = leech._already_present(ROW_4208, target)
    assert present is not None and [r.file_name for r in present] == ["export.zip"]
    assert all(present)


def test_failed_results_are_falsy() -> None:
    for outcome in leech.Outcome:
        result = leech.FileResult(ROW_4208, outcome)
        assert bool(result) is (outcome in (leech.Outcome.DOWNLOADED, leech.Outcome.ALREADY_PRESENT))


def test_manifest_is_sorted_merged_and_order_independent(tmp_path: Path) -> None:
    older = leech.FileResult(leech.ListingRow(**{**ROW_4208.__dict__, "download_id": 2}), leech.Outcome.HTTP_ERROR, "a.zip")
    newer = leech.FileResult(ROW_4208, leech.Outcome.DOWNLOADED, "b.zip", "/x/b.zip", 5906047)
    first, second = tmp_path / "one.jsonl", tmp_path / "two.jsonl"
    leech.write_manifest(first, [newer, older])
    leech.write_manifest(second, [older, newer])
    assert first.read_bytes() == second.read_bytes()
    retried = leech.FileResult(older.row, leech.Outcome.DOWNLOADED, "a.zip", "/x/a.zip", 5052615)
    leech.write_manifest(first, [retried])
    lines = first.read_text().splitlines()
    assert len(lines) == 2
    assert [('"outcome": "downloaded"' in line) for line in lines] == [True, True]


@pytest.mark.parametrize("name", ["../../etc/passwd", "..", "a/b\\c", ".hidden"])
def test_safe_path_component_cannot_escape(tmp_path: Path, name: str) -> None:
    cleaned = leech._safe(name)
    assert "/" not in cleaned and "\\" not in cleaned and not cleaned.startswith(".")
    assert (tmp_path / cleaned).resolve().parent == tmp_path.resolve()
