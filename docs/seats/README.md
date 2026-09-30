# Seats: standing roles a session can take in this repo

A **seat** is a fixed role a Claude Code session takes for its whole run. Every seat loads the same
`CLAUDE.md`, so the house rules are shared. What differs is the seat's own priming file: who it is,
what it owns, what it never does and how it starts. A session becomes a seat by being pointed at that
file, and it can be reset by starting a fresh session and pointing it there again.

| Seat | Priming | Notes | Inbox | What it does |
|---|---|---|---|---|
| **triage** | [TRIAGE.md](TRIAGE.md) | none: its state is `docs/FEEDBACK.md` | `docs/FEEDBACK.md` | Answers consumer feedback: triages, replies in the document, routes, archives. Keeps the loop's scripts in sync with the gist. Never designs a major change unasked. |
| **upstream liaison** | [UPSTREAM_LIAISON.md](UPSTREAM_LIAISON.md) | none: its state is `docs/POLARS_BIO_PENDING_FIXES.md` | none | Owns what we know about polars-bio and the layers under it: reproduces, locates, decides verdicts, drafts filings. Never files without the user's go. |

## Priming a seat

Start a new session in this checkout and give it one line: `prime from docs/seats/<SEAT>.md`. The seat
reads its priming file, then its notes file if it has one, then follows the startup steps in the
priming file. For the triage seat that line also arms the inbox watcher (see [TRIAGE.md](TRIAGE.md)).

## What a priming file holds

- **Why the seat exists**, in a sentence or two, so the next session can tell a rule from a habit.
- **What it owns**: files, id series, scripts. One writer per file; a seat that does not own a file
  reports to the seat that does.
- **How work reaches it**: a message, a watcher line, a file in an inbox.
- **Numbered startup steps.**
- **What it never does**, written out in full. Prohibitions from `CLAUDE.md` that bear on the seat are
  repeated there, because a prohibition behind a link does not get read.
- **What to write before a reset.**

It points at the procedure, it does not restate it. The triage algorithm lives in the loop runbook and
the upstream procedure in the ledger's header; a second copy drifts from the first.

## Notes files, and when a seat has none

A seat whose working state is not already in a tracked document gets `<SEAT>_NOTES.md` beside its
priming file: what is in flight, what it waits on, what it learned and has not yet written anywhere
durable. It is rewritten to the present before a reset, not appended to, and committed by path. It
holds state, not rules; a lesson still true next week moves into the priming file or `CLAUDE.md`.

**Both seats here have none, on purpose.** The triage loop derives its state from `docs/FEEDBACK.md` by
fingerprinting the reporter's text, and the liaison's state is its ledger. A notes file for either would
be a second copy of state that already has one home, which is the side-car the loop exists to remove.
Add one only for state that genuinely has nowhere else to live.

## Inboxes

A seat's inbox is where work for it lands when nobody can message it directly. The triage seat's inbox
is the repo's consumer inbox, `docs/FEEDBACK.md`, with its own header saying how to add an item. A seat
that needs a private one gets a lowercase folder `docs/seats/<seat>_inbox/` with a `README.md` giving the
shape of one request and a `LOG.md` with one line per request taken up (UTC time, file, what it became).

## When more than one seat runs at once

Seats run in the same checkout, so the working directory does not say which seat a session is.

- **Find peers with `ListAgents` each session.** Session names change on restart; never reuse one from
  memory.
- **Announce once after priming**, to each peer working in this checkout:
  `SEAT ANNOUNCE: I am the <seat> seat, session <name>, primed <UTC time>. Send me: <what now comes to
  this seat>. Answer with SEAT REPLY.` A seat that receives one answers once with `SEAT REPLY` in the same
  shape without the last sentence, and never answers a reply. That ends the round.
- **Two sessions primed as one seat is a conflict.** Stop and tell the user; neither acts on the seat's
  work until the user says which keeps it. Two triage seats would answer the same item twice.
- **Another seat's uncommitted file is work in flight.** Commit your own paths by name; never stash,
  revert or commit theirs.

## What stays true for every seat

Repeated from `CLAUDE.md` because these are the ones a seat is most likely to meet:

- Never `git add -A` or `git add .`; stage explicit paths and read `git status` first. Commit as you go;
  never `git push`, publish or release unless told. No `--amend`, rebase, `reset --hard` or force.
  Never `git stash drop` or `git stash clear`.
- Never `uv pip install`; never bare `python`. The stdlib triage scripts run directly, never under `bash`.
- Never edit or commit a sibling repo beyond appending an item to its `docs/FEEDBACK.md`.
- Never resolve a contradiction between two rules by inference: run the questionnaire in `CLAUDE.md` §1.
