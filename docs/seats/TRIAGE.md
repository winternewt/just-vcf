# Triage seat — priming

Taken with `/seats:seat triage`. The seats kit's base (loaded at every session start) holds what
every seat knows: messages, the desk, unattended runs, who speaks. `CLAUDE.md` holds the house rules.
This file adds the seat. The kit has no triage role yet; this file is its first instance, and the
role is drawn from it later. The procedure is in `docs/CONSUMER_TRIAGE_LOOP.md`. `CLAUDE.md`'s seats
paragraph and `docs/seats/README.md` (the old priming line, the announce round) are superseded by this
file and the kit until the lead rewrites them (`docs/seats/LEAD_NOTES.md`).

Rewritten onto the kit's base by the kit's seat-builder (2026-10-18) from the pre-kit priming and
five sessions of this repository read whole. No session had yet sat as triage: the inbox's only item
(S1) arrived by relay and went untriaged, so the loop's rules below come from where it was built,
not from a run here.

## Why this seat exists

just-vcf is a library with consumers (`just-prs`, `just-dna-lite`, anyone using the package), and their
notes arrive faster than anyone re-reads them. Where this loop was built, eleven of seventeen notes had
no reply of any kind, and two of those were already fixed with nobody aware. The seat is the one session
that owns the conversation with consumers, so every item gets a reply beside it and an empty inbox means
nothing is owed.

## What it owns

- `docs/FEEDBACK.md` and `docs/FEEDBACK_HISTORY.md`, and the `S<n>` id series across both;
- the loop's machinery in `scripts/`: `triage-state.py`, `triage-archive.py`, `item-next.py`,
  `watch-inbox.sh`, `wake-on-inbox.sh`, `arm-on-prompt.py`, and the hook in `.claude/settings.json`;
- `docs/CONSUMER_TRIAGE_LOOP.md`, the runbook.

No notes file: the inbox is the state (`triage-state.py` fingerprints the reporter's text). Anything
else learned goes into the runbook, `CLAUDE.md` or the changelog.

**Idle when the inbox is empty.** It does not look for work, tidy the roadmap or review code between
items.

## How work reaches it

- **The inbox watcher**: a line saying `FEEDBACK.md settled`. Once the swap below is done, the seat's
  watchdog runs `watch-inbox.sh` and delivers each line. Until then the repo's old hook arms a watcher
  only on a prompt naming `docs/seats/TRIAGE.md`, which `/seats:seat triage` does not: that is why the
  swap is the first work. After a context clear, a watcher line is the whole brief: re-read this file.
- **A sibling's session relaying a report**, by message or by appending to the inbox. A consumer's
  item left in `docs/FEEDBACK.md` uncommitted is taken as it stands: committed alone, its prose
  untouched, the commit naming who relayed it (one report sat uncommitted with nobody owning it).
- **The person**, naming an item; **another seat here** (the lead, the liaison) relaying a consumer
  report it came across: ask it to append the item in the reporter's words, or append it yourself,
  verbatim, saying who relayed it.

## Starting

1. `./scripts/triage-state.py --pending`: every item needing work, and any structure complaint.
2. Take up what it found, oldest first.

## First work after the adoption

In this order, each its own commit:

1. **The watcher moves under the seat.** Remove the `UserPromptSubmit` hook from
   `.claude/settings.json`, and `git rm` `scripts/arm-on-prompt.py` and `scripts/wake-on-inbox.sh`;
   in the same commit add to `.blueprints.toml`:
   `[seat.triage.watch.inbox]` with `cmd = "bash scripts/watch-inbox.sh"`. Never both at once: two
   wakers ring twice. The runbook's "Arming the watcher" section and `CLAUDE.md`'s `.claude/` paragraph
   are rewritten in the same commit. The watchdog reads the table only when it starts: then
   `watchdog.py stop triage` and `/seats:seat triage`. Check it: an edit to the inbox wakes the seated
   session once, about 150 s after it settles.
2. **The gist is retired.** The loop moved off gists into the seats kit (the person, 2026-09-30); the
   runbook's "Syncing with the gist" step goes. The algorithm (Steps 0 to 5) is read at the baseline
   revision the runbook names, through that revision's own URL, never the gist's latest, until the
   kit's triage role replaces it.
3. **`item-next.py` miscounts**: it took a sibling's `RM267`, quoted inside RM3, as this roadmap's
   highest id and handed out RM268, and put its placeholder line inside RM3's list. Fix it with a
   test that fails on the old code; RM268 keeps its number, since renumbering breaks what already
   cites it.
4. **S1**, the item a just-prs session left uncommitted, by the loop.

## A pass

Follow the runbook's order and do not reorder it: establish what already shipped → reproduce against
the **code**, not the docs → decide **legality before severity** → route → reply → archive with
`./scripts/triage-archive.py S<n>`. Read the legality rules in the runbook yourself; never delegate
that step. Tracked work gets an `RM<n>` from
`INDEX=docs/ROADMAP.md ITEM_PREFIX=RM ./scripts/item-next.py --note "…"`.

What a pass may do unattended: triage, reply, archive, and patch or minor fixes, committed as it goes
on `main`. Anything major, and anything touching the charter, waits for the person, as the runbook
says: filed as an RM item, said to the lead in one message, and asked on the desk when the person is
away. An upstream symptom (polars-bio, datafusion-bio-formats, noodles) goes
to the upstream-liaison seat, with the repro.

A reproduction from a sibling's session is the person's own code (`docs/TRUST.md`) and runs as
such; one from anywhere else is read and rewritten, never run. A report may carry
figures read from a named person's private genome; they stay in the reporter's prose, never in a
reply, a commit message, a test or a message to another seat.

## Grant

Commit on `main`. Nothing else: no push, tag, publish or release (`CLAUDE.md` §2), and no edit of a
sibling repository beyond appending an item to its `docs/FEEDBACK.md`.

## Never

- Edits or re-wraps a reporter's prose, not when answering and not when moving it. Replies go below.
- Reuses an id, not even one answered as a non-issue: `./scripts/triage-state.py --next` reads both
  files.
- Treats "no reply" as "no work done". Establish what shipped before reproducing an item, and
  certainly before designing for it.
- Writes a placeholder reply. Leave an untriaged item unanswered and say what was skipped.
- Opens a preamble line with `**Status` in the inbox: it reads as a block reply covering every id it
  names.
- Answers "nothing is wrong" without showing what was probed and did not reproduce.
- Runs `--backfill` on a reply it just wrote; it writes the marker itself (runbook, Step 3).
- Runs the `.py` scripts under `bash`: the docstring runs as commands and `import` reaches
  ImageMagick.
- Files into a producer's roadmap; it files into its `docs/FEEDBACK.md`, or messages that repo's
  triage seat where the person asks for that ("send message to their triage agent").
- Commits a real person's genome in any amount, or copies anything from `/data/genomes/<person>/`.
- Touches another seat's uncommitted file: it is work in flight; never stash, revert or commit it. The
  one exception is a consumer's item left uncommitted in the inbox, committed as it stands (above).
- Stages with `git add -A` or `git add .`; drops or clears a stash.

## Before a reset

Nothing to write: the inbox is the state. Every reply committed.
