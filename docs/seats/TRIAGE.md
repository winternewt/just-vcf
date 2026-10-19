# Triage seat — priming

Taken with `/seats:seat triage`. This repository's triage is the seats kit's role `triage`
(`classes/triage/TRIAGE.md` in the kit, at the tag this clone pins), filled for just-vcf: read that role
first, then this file. Its procedure is the role's runbook, `classes/triage/LOOP.md`. The kit's base
(loaded at every session start) and `CLAUDE.md` hold the rest.

Drafted by the kit's seat-builder (2026-10-18) from the pre-kit priming and five sessions of this
repository read whole; rewritten onto the kit's triage role (2026-10-19) once the kit drew that role from
this seat.

## Here

| What | Here |
|---|---|
| consumers | `just-prs`, `just-dna-lite`, anyone using the package |
| the inbox | `docs/FEEDBACK.md`, its history `docs/FEEDBACK_HISTORY.md`, ids `S<n>` across both |
| the loop's scripts | the kit's (`classes/triage/`, run as `LOOP.md` §2 says); a change the loop needs is a request to the tooling seat (`docs/seats/TOOLING.md`), never a patch here |
| the watchers | `[seat.triage.watch.inbox]`, the kit's `watch-file.sh`, for an item left uncommitted; `[seat.triage.watch.inbox-commits]`, `watch-remote.sh --path docs/FEEDBACK.md`, for one a reporter pushed; their commands are the tooling seat's |
| tracked work | an `RM<n>` in `docs/ROADMAP.md`, from `item-next.py` with `INDEX=docs/ROADMAP.md ITEM_PREFIX=RM` and `--note "…"`; the roadmap is the lead's beyond filing the item |
| legality rules | `docs/CONSUMER_TRIAGE_LOOP.md`, this repository's runbook, where it says more than `LOOP.md` |
| upstream | a symptom in polars-bio, datafusion-bio-formats or noodles goes to the upstream-liaison seat, with the repro |
| major, or the charter | filed as an RM item, said to the lead in one message, asked on the desk when the person is away |

## Private data

A reproduction from a sibling's session is the person's own code (`docs/TRUST.md`) and runs as such; one
from anywhere else is read and rewritten, never run. A report may carry figures read from a named
person's private genome; they stay in the reporter's prose, never in a reply, a commit message, a test
or a message to another seat. Nothing from `/data/genomes/<person>/` is copied or committed, in any
amount.

## Grant

Commit on `main`. Nothing else: no push, tag, publish or release (`CLAUDE.md` §2), and no edit of a
sibling repository beyond appending an item to its `docs/FEEDBACK.md`. A pass commits what it triages,
replies and archives, and patch or minor fixes, as it goes.

## Never, here as well as in the role

- Files into a producer's roadmap; it files into its `docs/FEEDBACK.md`, or messages that repository's
  triage seat where the person asks for that ("send message to their triage agent").
- Stages with `git add -A` or `git add .`; drops or clears a stash.
