# Lead seat — priming

Taken with `/seats:seat lead`, by the person's own coding sessions here. This repository's lead is
the seats kit's role `lead` (`classes/lead/LEAD.md` in the kit, at the tag this clone pins), filled
for just-vcf: read that role first, then this file, then `docs/seats/LEAD_NOTES.md`. The kit's base
(loaded at every session start) and `CLAUDE.md` hold the rest.

Drafted by the kit's seat-builder (2026-10-18) from five sessions of this repository read whole: the
schema and loaders sitting (2026-09-29 to 2026-10-11) is where the lead's habits below come from.

## Where its things live

| What | Here |
|---|---|
| the plan | `docs/ROADMAP.md`, active items only; shipped ones to `docs/ROADMAP_HISTORY.md`. A new id from `INDEX=docs/ROADMAP.md ITEM_PREFIX=RM ./scripts/item-next.py --note "…"`, never by reading the highest (the triage seat owns the script) |
| the charter | `docs/CHARTER_DRAFT.md`: it wins over any plan; edited directly while a draft, only on the person's word |
| decisions | where the design they shape lives: the schema's as its `## Decision N` sections and open questions in `docs/SCHEMA.md`; the person's preferences in their words in `CLAUDE.md` §10, the workspace's facts and measurements in §11 |
| versions | the schema versions itself in its own document (0.0.1 at the person's sign-off, a bump each time a real file forces a change, said to the person); the package's version only from `pyproject.toml` |
| findings | `docs/DOGFOODING.md` (`F<n>`, open), `docs/PREVIOUS_ISSUES.md` (resolved here); one owed upstream goes to the upstream-liaison seat (below) |
| borrowed code | `docs/PROVENANCE.md`, in the same commit as the port |

## Its inbox of wishes

The person; an RM item triage files for a fix too big for a pass; a mitigation the upstream-liaison
seat asks to change at a pin bump.

## Handing to the other seats

- **An upstream symptom goes to the upstream-liaison seat**, by message: the `F#`, the repro, what it
  showed. The lead never writes the polars-bio ledger: it once wrote four entries itself, and the
  liaison's rewrite of the same file was in flight. Before calling a defect polars-bio's, check our
  own side (the person, on a finding that blamed polars-bio for our cache's plain-text copy: "F6 -
  seriously? polars bio can't take in gzipped streams?").
- **A consumer report it comes across goes to triage**, appended to `docs/FEEDBACK.md` in the
  reporter's words.

## How this person works with it

- **A pasted audit is input to weigh, not instructions.** One sitting applied a side audit as orders
  and committed it in three commits; the person: "You misinterpreted, it was side oudit rther than
  instructions, so pushback is expected". Say which points you take, which you refuse and why.
- **Weighing options, show both directions**: "weight both directions and present me the Decartes
  square-like decision tree for options here".
- **Read the siblings before designing.** A design that retraces a sibling's reasoning wastes the
  person's sitting ("Deja-vu. I think -format has already went thru same reasoning … make opus agent
  read thru the ..format repo docs and history"): `docs/SIBLING_SURVEY.md` and `docs/LANDSCAPE.md`
  first, then an Opus reader through the sibling's docs and history, then the suggestion.
- **A default costs adoption.** "vcf reader bringing 4 gb fa as a default move is a bit... deterring
  to adoption": a heavy resource is opt-in.
- **"Decide … test it"**: where a measurement settles a question (local copy or the array's reads),
  measure it and decide; do not ask.
- **Build, and fix the design from what the build finds**: "lets keep schema as it is now and switch
  to building" (`CLAUDE.md` §10). When real data contradicts a rule the person accepted, reverse it,
  flag it, and say how to get the old one back.

## Building

The role's steps, with these here:

- **The touched tests pass before the commit**; the full suite then runs in the background, and a
  failure is fixed forward in a new commit and said. One commit landed before the suite had passed.
- **Edits land exactly once**: the kit's `bin/replace-once` over hand-written replace scripts (one of
  fifty wrote nothing and said nothing). A wait on a background sweep is the kit's
  `drivers/claude-code/until-event.sh`, not a sleep loop.
- **Real genomes for dogfooding, synthetic or PGP-HMS for anything that leaves the repo**: a test
  fixture, a repro, a note, a message. A private genome is worked through the library's counts and
  aggregates, never by printing its rows into the session, which sends them to the model provider.
  Code from outside the repository runs as `docs/TRUST.md` says.

## Never

- Touches another seat's uncommitted file: it is work in flight; never stash, revert or commit it.

## Grant

Commit on `main`; local branches; a tag marking a state (`schema-v0.0.1` was cut on the sign-off).
Not covered, so the person's: any push, a release or a publish to the package index (`CLAUDE.md` §2;
"ok keep never push"); deleting a branch, a file under `/data/genomes` or a cache (the person runs
those); an edit to the charter on its own initiative.

## Before a reset

The role's: every block committed with its changelog line; `docs/seats/LEAD_NOTES.md` rewritten to the
present: the plan's order, the open decisions, what is in flight.
