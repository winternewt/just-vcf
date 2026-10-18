"""Tests for scripts/item-next.py: which ids count as taken, and where a reservation row lands.

The script reads its configuration from the environment at import time and its name has a hyphen, so it
is driven as a subprocess against a docs tree under tmp_path, the way the runbook runs it.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "item-next.py"


def run(docs: Path, *args: str) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, "INDEX": str(docs / "ROADMAP.md"), "SCAN_DIR": str(docs), "ITEM_PREFIX": "RM"}
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], env=env, capture_output=True, text=True, check=True
    )


def headed_ids(*texts: str) -> set[int]:
    return {int(n) for t in texts for n in re.findall(r"^## RM(\d+)\b", t, flags=re.M)}


def test_a_sibling_id_quoted_in_prose_is_not_taken(tmp_path: Path) -> None:
    # The incident: RM3's body quoted just-dna-format's RM267, and the allocator handed out RM268.
    roadmap = (
        "# Roadmap\n\n## Items\n\n"
        "## RM1 — Plan\nopen\n\n"
        "## RM3 — REF lookup\n- just-dna-format mints ids offline. Their RM267 shows a variant table\n"
        "  used as the REF source producing wrong events.\n"
    )
    history = "# Roadmap history\n\n## RM2 — Shipped thing\ndone\n"
    (tmp_path / "ROADMAP.md").write_text(roadmap)
    (tmp_path / "ROADMAP_HISTORY.md").write_text(history)

    out = run(tmp_path, "--dry-run").stdout.strip()

    assert out == f"RM{max(headed_ids(roadmap, history)) + 1}"


def test_reservation_and_tombstone_rows_still_count(tmp_path: Path) -> None:
    (tmp_path / "ROADMAP.md").write_text("# Roadmap\n\n## Items\n\n## RM1 — Plan\nopen\n")
    first = run(tmp_path).stdout.strip()
    run(tmp_path, "--release", first)
    second = run(tmp_path).stdout.strip()

    assert first == "RM2"
    assert second == "RM3", "a released id is spent, never handed out again"


def test_reservation_lands_in_the_anchor_section_not_in_a_later_items_list(tmp_path: Path) -> None:
    # `## Items` is empty; the old walk ran past `## RM1` and `## RM2` into RM2's bullet list.
    (tmp_path / "ROADMAP.md").write_text(
        "# Roadmap\n\n## Items\n\n## RM1 — Plan\nopen\n\n## RM2 — Loading\n- first bullet\n- second\n"
    )
    number = run(tmp_path, "--note", "a probe").stdout.strip()

    lines = (tmp_path / "ROADMAP.md").read_text().splitlines()
    row = next(i for i, line in enumerate(lines) if f"**{number}**" in line)
    items = lines.index("## Items")
    following = next(i for i in range(items + 1, len(lines)) if lines[i].startswith("## "))
    assert items < row < following, "\n".join(lines)
    assert lines[row - 1] == "", "a list row directly under prose or a heading reads as part of it"
    assert lines[row + 1] == "" and lines[row + 2].startswith("## "), "one blank line, then the next item"
