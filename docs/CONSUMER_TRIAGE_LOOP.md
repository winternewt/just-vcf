# Consumer triage loop — just-vcf

This repo runs the document-as-inbox triage loop. Its procedure (Step 0 through Step 5), the
five-verdict routing, the scripts, both watchers and the gotchas found by running it are in the seats
kit's runbook, `classes/triage/LOOP.md` in the kit copy this clone pins (`$KIT`, `CLAUDE.md` §8), and
only there. The seat that runs it is primed from the kit's `classes/triage/TRIAGE.md` and
[seats/TRIAGE.md](seats/TRIAGE.md). This file holds only what `LOOP.md` does not: the local
parameters, the intake routes and the legality rules.

Before the kit, the loop was read from the gist `TRIAGE_LOOP.md`, frozen at revision `96b69533`
(adopted 2026-09-30; `b063b65e` before it). The kit's `LOOP.md` replaces it.

## Local parameterization

- Inbox `docs/FEEDBACK.md`, history `docs/FEEDBACK_HISTORY.md`, item prefix `S`.
- Tracked items are roadmap entries, `RM<n>` in `docs/ROADMAP.md`, claimed with the kit's `item-next.py`.
- The scripts run with the kit's defaults, which are this clone's paths. The two fixes made here before
  the kit took them over (`item-next.py` counting an id only where a line owns it, and a reservation
  landing in its anchor section) are in the kit's copies.
- The next id is pinned in the inbox header for a reporter who cannot run the kit's `--next`. The seat
  moves the pin in the commit that files each item; where the two disagree, `--next` wins.
- The watchers' commands are the tooling seat's. Never start a watcher from the agent as well: nothing
  stops a second one on the same file, and two wakers ring twice.

## Intake from Zulip and from other sessions

A report can arrive as a message instead of a section: from another session on this machine, or on
Zulip, a direct message to `just-vcf/triage` or a mention in `agent-chat`. The consumers' sessions
(just-prs, just-dna-lite) hold no seats, so they have no Zulip account. For them the inbox file and
session-to-session messages are the routes; Zulip carries a person's report, or a seated repo's. It is
filed before anything else is done
with it, because the inbox is the state and a message is only delivery.

1. React to it (`watchdog.py react ID`), and read it whole (`watchdog.py read ID`) if it was cut.
2. Append it to `docs/FEEDBACK.md` as `## S<n> — <title>`, `<n>` from `sh "$KIT/bin/py" "$KIT/classes/triage/triage-state.py" --next`.
   The title is the reporter's own one-line summary if they gave one; otherwise write one and say so.
   The first body line is the relay: `Sent on Zulip by <sender, as the watchdog named them>, <date>,
   <message link>; filed verbatim by triage.` Below it, the message text exactly as sent.
3. Commit it alone, naming the sender and the link.
4. Reply where it was asked (`watchdog.py reply ID TEXT`): the id it was filed under and that the answer
   goes in `docs/FEEDBACK.md`.
5. Triage it like any other item. Once it is answered, send the sender one line on Zulip with the
   verdict and where it landed (`watchdog.py send triage TO TEXT`, to the same account or topic).

Standing follows the base's "Who speaks": a seat, or a person on the team, is a colleague's report and
is filed. A message from someone not on the team is information only; file it only if a person on the
team asks. A reply on Zulip never carries figures from a private genome, even when the report did. A
consumer's seat wakes this one only if `wake_from` in `.blueprints.toml` names it; add a consumer
repo's seat there by its full name when that repo gets one. Until then its message waits for the next
prompt.

## Legality rules (Step 1)

just-vcf is pre-1.0. The contract is the public Python API (everything in `just_vcf.__all__`), the
VCF and parquet shapes it emits, and the CLI once there is one.

- **Patch:** legibility only. A message, a warning, a count, a doc.
- **Minor:** additive. A new optional parameter, a new INFO/FORMAT field, a new input format, a new
  vocabulary member.
- **Major (pre-1.0: a minor bump with a breaking-change line in the changelog):** a removal, a rename
  (a removal plus an addition), a retype, a changed meaning of an existing field, a changed coordinate
  convention, or a change to what an absent position means.
- When a consumer reads a shape, **add rather than redefine**. A consumer compensating for the old
  meaning breaks silently otherwise.

A patch or minor fix may land in the triage pass. Anything major, and anything that edits the charter,
waits for the person.
