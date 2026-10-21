# Seats: standing roles a session can take in this repo

A **seat** is a fixed role a Claude Code session takes for its whole run. The seats here run on the
seats kit (`agent-seat-blueprints`): its plugin is pinned in this clone, and its base, loaded at every
session start, holds what every seat knows: messages, the desk, unattended runs, who speaks, one
writer per file. Every seat also loads `CLAUDE.md`, so the house rules are shared. What differs is the
seat's own priming file: why it exists, what it owns, how work reaches it, its grant and what it
never does.

Take a seat with `/seats:seat <seat>` in a session started in this checkout. A session that has not
taken one does the task it was given and takes no seat's work.

| Seat | Priming | Notes | Reached by | What it does |
|---|---|---|---|---|
| **triage** | [TRIAGE.md](TRIAGE.md) | none: its state is `docs/FEEDBACK.md` | the inbox watcher; a relayed consumer report | Answers consumer feedback: triages, replies in the document, routes, archives. Owns the runbook. Never designs a major change unasked. |
| **upstream-liaison** | [UPSTREAM_LIAISON.md](UPSTREAM_LIAISON.md) | none: its state is `docs/POLARS_BIO_PENDING_FIXES.md` | an upstream symptom from the lead or triage; its own checks at each start | Owns what we know about polars-bio and the layers under it: reproduces, locates, decides verdicts, drafts filings. Never files without the user's go. |
| **lead** | [LEAD.md](LEAD.md), on the kit's `lead` role | [LEAD_NOTES.md](LEAD_NOTES.md) | the user; an RM item from triage; a mitigation from the liaison | The user's own coding sessions: the roadmap, the design, the build. Hands upstream symptoms to the liaison and consumer reports to triage. |
| **tooling** | [TOOLING.md](TOOLING.md), on the kit's `plumber` role | none: its state is its issues | a tool's fault from any seat; a kit release or reply | Owns the machinery: the watcher commands, the leak gate's install, the kit's pin (it alone moves it). The triage loop's scripts are the kit's; it files their faults and the kit's other gaps on the kit's issues. Never judges what its tools carry. |

`.blueprints.toml` declares the four seats and who may wake each one. `team.toml` lists the people on
the team by chat account. `docs/TRUST.md` is the trust register every seat checks before outside code
runs here.

## Notes files, and when a seat has none

A seat whose working state is not already in a tracked document gets `<SEAT>_NOTES.md` beside its
priming file: what is in flight, what it waits on, what it learned and has not yet written anywhere
durable. It is rewritten to the present before a reset, not appended to, and committed by path. It
holds state, not rules; a lesson still true next week moves into the priming file or `CLAUDE.md`.

The lead has one. Triage, the liaison and the tooling seat have none, on purpose. The triage loop
derives its state from `docs/FEEDBACK.md` by fingerprinting the reporter's text, the liaison's state is
its ledger, and the tooling seat's is its issues. A notes file for any of them would be a second copy of
state that already has one home.

## What stays true for every seat

Repeated from `CLAUDE.md` because these are the ones a seat is most likely to meet:

- Never `git add -A` or `git add .`; stage explicit paths and read `git status` first. Commit as you go;
  never `git push`, publish or release unless told (the lead's and the tooling seat's grants tell them to push `main`
  as a fast-forward). No `--amend`, rebase, `reset --hard` or force.
  Never `git stash drop` or `git stash clear`.
- Another seat's uncommitted file is work in flight: commit your own paths by name; never stash,
  revert or commit theirs. The one exception is a consumer's item left uncommitted in the inbox, which
  triage commits as it stands.
- Never `uv pip install`; never bare `python`. The kit's triage scripts run through the kit's `bin/py.sh` (`CLAUDE.md` §8), never under `bash`.
- Never edit or commit a sibling repo beyond appending an item to its `docs/FEEDBACK.md`.
- Never run outside code without checking its standing in `docs/TRUST.md`.
- Never resolve a contradiction between two rules by inference: run the questionnaire in `CLAUDE.md` §1.
