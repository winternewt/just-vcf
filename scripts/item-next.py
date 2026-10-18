#!/usr/bin/env python3
"""Allocate the next tracked-item id, and **reserve** it so two sessions cannot take the same one.

See TRIAGE_LOOP.md § 2. This is the companion to `triage-state.py --next`: that one allocates the
*inbox* ids a reporter writes, this one allocates the ids of the tracked items route **(a)** files.

**An index is not an allocator.** A complete index of your tracked items answers *where does item 47
live*; it cannot answer *which number is mine*, because reading a number claims nothing. So the usual
procedure — grep the highest id, add one, write the entry — leaves a window between the read and the
write, and that window is exactly where a second session reads. Both write the same number.

**That is a reproduced incident, not a hypothesis.** In the project this pattern came from, two agent
sessions sharing one working tree filed different work under the same id a minute apart, and the tree
carried two entries pointing at different items; one had to be renumbered afterwards, and the choice of
which was decided by how many other documents already referenced each. Grepping cannot fix it. The claim
has to be a single atomic write.

So this does both halves in one critical section:

* **Scan** every markdown file under `SCAN_DIR` for the lines that *own* an id, and compute the next
  number. The whole tree rather than the index alone, because an item outlives its place in the index:
  shipped ones move to a history file under their own heading.
* **Reserve** it by appending a placeholder row to `INDEX` while holding an advisory lock, so the next
  scan sees the number as taken — the placeholder contains it.

**The lock is on the directory, never on the index file.** Two reasons, and the second is measured
rather than assumed. A lockfile left behind by exactly the crash it guards against would block every
later run, and the staleness rule that repairs that is a clock — whereas `flock` dies with the process,
so there is nothing to go stale and nothing to clean up. And `flock` binds an **inode**: an editor or an
atomic writer that renames a new file over your index leaves the holder locking an unlinked inode while
a second process opens the new file and acquires immediately. Verified in a sandbox before this was
written; locking the file you are appending to looks correct and excludes nothing.

**An id is taken where a line owns it, never where prose mentions it.** Three forms own one: an item's
heading (`## RM47 — …`, in the index or the history file), a reservation row and a tombstone row. A
mention in a body does not, because nothing in free prose tells ours from a sibling's: just-vcf's RM3
quoted just-dna-format's "Their RM267", and the first cut of this scan, which counted every mention,
handed out RM268. Under the allocator a number is claimed only by a row this tool writes, so an id
that owns no line was never claimed here.

**A reservation is a visible row in the index, not a side-car.** A state file the index cannot see is
how a number goes missing, which is the failure a complete index exists to prevent. So a claimed and
then abandoned number is legible rather than silently burned: replace the row with the item's real entry
when you write it, or `--release` it.

**`--release` leaves a tombstone, and the first cut of this shipped the bug it fixes.** Deleting the row
made the number invisible to the scan, so a released id was immediately handed back out — against the
same rule the inbox ids follow (TRIAGE_LOOP.md § 1: ids are never reused, because whatever argued the
withdrawal refers to the number). The tombstone must *contain* the number literally, since a scan is all
that reads it.

Configuration (environment, all optional):
    INDEX        the document that indexes your tracked items  (default ./ROADMAP.md)
    SCAN_DIR     directory tree scanned for used ids           (default INDEX's parent)
    ITEM_PREFIX  tracked-item id prefix                        (default RM, so `RM47`)
    ANCHOR       heading under which a reservation is inserted  (default: the first `## ` heading)

Python, and named `.py` for it: run it or pass it to `python3`, never to `bash`. See the extension
gotcha in TRIAGE_LOOP.md § 5 — under bash this docstring is executed as commands.

Usage:
    item-next.py                  # reserve the next id and print it
    item-next.py --dry-run        # print what it would reserve, write nothing
    item-next.py --note "..."     # a few words on what it is for, into the placeholder
    item-next.py --list           # standing reservations
    item-next.py --release RM47   # drop a reservation you are not going to use
    item-next.py --help           # this text, and no reservation

**An unknown flag is refused, not ignored, and that is load-bearing rather than tidy.** The reserve
path is also the no-flag path, so anything this tool did not recognise used to fall through to
`allocate()` — `--help` spent a real id that way, and a typo spent a second one while the first was
being repaired. Both keep a tombstone row, since ids are never reused. A tool whose default path
mutates owes its argument surface a refusal.
"""

import os
import pathlib
import re
import sys

try:
    import fcntl as _FCNTL
except ImportError:  # pragma: no cover - non-POSIX; the degradation is announced, never silent
    _FCNTL = None

ITEM_PREFIX = os.environ.get("ITEM_PREFIX", "RM")
# Local to just-vcf: default derived from this script's location, not $PWD.
INDEX = pathlib.Path(os.environ.get("INDEX", pathlib.Path(__file__).resolve().parent.parent / "docs" / "ROADMAP.md")).resolve()
SCAN_DIR = pathlib.Path(os.environ.get("SCAN_DIR", INDEX.parent)).resolve()
#: Where a reservation is inserted. Default is the first `## ` heading, which is the least surprising
#: place in an unknown document; set it if your index has a natural open-items section. It matters
#: because a row appended at end-of-file lands under whatever prose closes the document and reads as
#: part of it — the same furniture hazard § 5 records for the inbox documents.
ANCHOR = os.environ.get("ANCHOR", "")

#: The reservation row. The marker is what a later scan counts, so it must contain the number.
RESERVED_RE = re.compile(rf"^- 🔷 \*\*reserved\*\* \*\*{re.escape(ITEM_PREFIX)}(\d+)\*\*")

#: A released reservation's row; `_tombstone` writes it.
TOMBSTONE_RE = re.compile(rf"^- ✖ \*\*{re.escape(ITEM_PREFIX)}(\d+)\*\*")

#: The lines that own an id: an item's heading, a reservation, a tombstone. Anchored at line start, so
#: an id quoted in a body (ours or a sibling's) is never read as taken.
OWNER_RES = (
    re.compile(rf"^#{{1,6}}\s+{re.escape(ITEM_PREFIX)}(\d+)\b"),
    RESERVED_RE,
    TOMBSTONE_RE,
)

LOCK_UNAVAILABLE = (
    "advisory locking is unavailable for {where} ({why}), so this reservation is NOT excluded from a "
    "concurrent one. Check that {prefix}{n} is unused before writing the entry."
)


def used_numbers() -> dict[int, str]:
    """`{number: first file it was seen in}` for every id a line owns anywhere under `SCAN_DIR`.

    A reservation placeholder counts as used, which is the whole mechanism.
    """
    seen: dict[int, str] = {}
    for path in sorted(SCAN_DIR.rglob("*.md")):
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            # Unreadable is not absent. A file this cannot open may hold the very number about to be
            # handed out, so it is reported rather than skipped silently.
            print(f"warning: could not read {path}, so its ids are unknown", file=sys.stderr)
            continue
        for line in text.splitlines():
            for owner in OWNER_RES:
                found = owner.match(line)
                if found:
                    seen.setdefault(int(found.group(1)), str(path))
    return seen


def next_free(used: dict[int, str]) -> int:
    """The next number above the current maximum.

    **Above the maximum, never filling a gap.** A gap can only mean an item was withdrawn, and its
    number stays referenced in whatever argued the withdrawal — so reusing it would make two different
    items answer to one id in the record. Same rule as the inbox ids.
    """
    return max(used, default=0) + 1


def reservations(text: str) -> dict[int, str]:
    """`{number: the row}` for every reservation currently standing in the index."""
    out: dict[int, str] = {}
    for line in text.splitlines():
        found = RESERVED_RE.match(line)
        if found:
            out[int(found.group(1))] = line
    return out


def _row(number: int, note: str | None) -> str:
    what = note.strip() if note and note.strip() else "claimed, entry not yet written"
    return (
        f"- 🔷 **reserved** **{ITEM_PREFIX}{number}** — {what}. Claimed by an allocator run so a "
        f"concurrent session cannot take the same id; replace this line with the item's real entry, "
        f"or drop it with `item-next.py --release {ITEM_PREFIX}{number}`."
    )


def _tombstone(number: int) -> str:
    """What a released reservation leaves behind: a row that still *contains* the number.

    It has to contain it literally, because `used_numbers` is a scan — a tombstone the scanner cannot
    see is the same as a deleted row, which is the bug this replaced.
    """
    return (
        f"- ✖ **{ITEM_PREFIX}{number}** — reserved and released without an item being written. The id "
        f"is spent rather than free: ids are never reused, so a later reference to "
        f"{ITEM_PREFIX}{number} can only mean this."
    )


def _insert(text: str, row: str) -> str:
    """Put `row` under `ANCHOR` (or the first `## ` heading): before that section's first list item, or
    at the section's end when it has none.

    **The walk stops at the next heading.** It used to run on to the first list item anywhere below, so
    under an empty `## Items` the row landed inside the bullet list of whichever item came first with
    one, and read as part of that item.

    The end-of-file fallback is announced rather than silent: a row that lands under a document's
    closing prose reads as part of it, and finding that out later is worse than being told now.
    """
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(ANCHOR) if ANCHOR else line.startswith("## "):
            at = i + 1
            while at < len(lines) and not lines[at].startswith(("- ", "#")):
                at += 1
            if at == len(lines) or lines[at].startswith("#"):
                # No list in the section: end it with the row, set off by one blank line on each side,
                # in place of the blank run that ended it.
                end = at
                while at > i + 1 and not lines[at - 1].strip():
                    at -= 1
                lines[at:end] = ["", row] + ([""] if end < len(lines) else [])
            else:
                lines.insert(at, row)
            return "\n".join(lines) + "\n"
    wanted = f"{ANCHOR!r} heading" if ANCHOR else "'## ' heading"
    print(
        f"warning: no {wanted} in {INDEX.name}; appending at the end of the file instead",
        file=sys.stderr,
    )
    return text.rstrip("\n") + "\n" + row + "\n"


class _Unlocked:
    """The degradation, carried as a value so the caller reports it rather than assuming a lock."""

    def __init__(self, why: str) -> None:
        self.why = why


def _lock(fd: int) -> "_Unlocked | None":
    """Take the exclusive lock, or say why not. Never returns a false "we have it"."""
    try:
        _FCNTL.flock(fd, _FCNTL.LOCK_EX)
    except OSError as exc:
        # A network filesystem may refuse or emulate `flock`. Blocking rather than `LOCK_NB` on
        # purpose: the critical section is one read plus one append, so a wait is milliseconds, and a
        # refusal would send the caller off to allocate by hand — which is the failure mode itself.
        return _Unlocked(str(exc))
    return None


def allocate(note: str | None, dry_run: bool) -> int:
    """Reserve and return the next id, holding the index's directory for the read-and-append."""
    if dry_run:
        number = next_free(used_numbers())
        print(f"--dry-run: would reserve {ITEM_PREFIX}{number} (nothing written)", file=sys.stderr)
        return number
    if _FCNTL is None:
        number = next_free(used_numbers())
        print("warning: " + LOCK_UNAVAILABLE.format(
            where=INDEX.parent, why="this platform has no fcntl", prefix=ITEM_PREFIX, n=number,
        ), file=sys.stderr)
        INDEX.write_text(_insert(INDEX.read_text(), _row(number, note)))
        return number

    fd = os.open(INDEX.parent, os.O_RDONLY)
    try:
        degraded = _lock(fd)
        # **Everything below is inside the critical section, and the scan has to be too**: reading the
        # maximum outside the lock and appending inside it is the same race with a smaller window.
        number = next_free(used_numbers())
        if degraded is not None:
            print("warning: " + LOCK_UNAVAILABLE.format(
                where=INDEX.parent, why=degraded.why, prefix=ITEM_PREFIX, n=number,
            ), file=sys.stderr)
        INDEX.write_text(_insert(INDEX.read_text(), _row(number, note)))
        return number
    finally:
        # Closing releases the lock; an explicit LOCK_UN first would be a second way to do one thing.
        os.close(fd)


def release(number: int) -> int:
    """Drop a standing reservation, leaving a tombstone so the id is never handed out again."""
    fd = os.open(INDEX.parent, os.O_RDONLY) if _FCNTL is not None else None
    try:
        if fd is not None:
            _lock(fd)
        text = INDEX.read_text()
        standing = reservations(text)
        if number not in standing:
            print(
                f"{ITEM_PREFIX}{number} is not a standing reservation in {INDEX.name}",
                file=sys.stderr,
            )
            return 1
        kept = [
            _tombstone(number) if line == standing[number] else line for line in text.splitlines()
        ]
        INDEX.write_text("\n".join(kept) + "\n")
        print(
            f"released {ITEM_PREFIX}{number} — the id stays spent, so the next allocation moves past it"
        )
        return 0
    finally:
        if fd is not None:
            os.close(fd)


KNOWN_FLAGS = frozenset({"--dry-run", "--note", "--list", "--release", "--help", "-h"})


def main() -> int:
    argv = sys.argv[1:]
    # The token after --note is its value, not a flag: a note may legitimately open with a dash.
    at_note = argv.index("--note") if "--note" in argv else None
    note_value = at_note + 1 if at_note is not None else None
    flags = {a for i, a in enumerate(argv) if a.startswith("-") and i != note_value}
    if flags & {"--help", "-h"}:
        print(__doc__.strip())
        return 0
    unknown = sorted(flags - KNOWN_FLAGS)
    if unknown:
        print(
            f"unknown flag(s): {' '.join(unknown)} — refusing rather than reserving a number; "
            f"see item-next.py --help",
            file=sys.stderr,
        )
        return 2
    note = None
    if "--note" in argv:
        at = argv.index("--note")
        if at + 1 >= len(argv):
            print("--note needs a value", file=sys.stderr)
            return 2
        note = argv[at + 1]

    if not INDEX.is_file():
        print(
            f"{INDEX} does not exist. Point INDEX at the document that indexes your tracked items.",
            file=sys.stderr,
        )
        return 2

    if "--release" in flags:
        rest = [a for a in argv if not a.startswith("--") and a != note]
        if len(rest) != 1:
            print(f"usage: item-next.py --release {ITEM_PREFIX}47", file=sys.stderr)
            return 2
        return release(int(rest[0].removeprefix(ITEM_PREFIX)))

    if "--list" in flags:
        standing = reservations(INDEX.read_text())
        if not standing:
            print("no standing reservations")
            return 0
        for number in sorted(standing):
            print(standing[number])
        return 0

    number = allocate(note, dry_run="--dry-run" in flags)
    used = used_numbers()
    print(f"highest in use: {ITEM_PREFIX}{max(used, default=0)}", file=sys.stderr)
    print(f"{ITEM_PREFIX}{number}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
