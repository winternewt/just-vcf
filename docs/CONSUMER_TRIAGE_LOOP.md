# Consumer triage loop — just-vcf

This repo runs the document-as-inbox triage loop. Its algorithm (Step 0 through Step 5), the
five-verdict routing and the gotchas found by running it are in `TRIAGE_LOOP.md` at gist revision
`96b6953362bb4854e05db04334086c43c9eeb170`, and only there:
<https://gist.githubusercontent.com/winternewt/54b94bda01812be937b892146d1bb254/raw/96b6953362bb4854e05db04334086c43c9eeb170/TRIAGE_LOOP.md>.
Read it through that revision's own URL, never the gist's latest. The loop moved off gists into the
seats kit (the person, 2026-09-30), so nothing is synced from the gist any more; this revision stands
until the kit's triage role replaces it. This file holds only what is local to just-vcf.

The session that runs this loop is the **triage seat**, primed from
[seats/TRIAGE.md](seats/TRIAGE.md). That file says what the seat owns, how it starts and what it never
does; this one holds the local rules the pass reads.

**Frozen at gist revision `96b6953362bb4854e05db04334086c43c9eeb170`, adopted 2026-09-30.** Before it,
`b063b65ef777cf06002e6fa14953d6503692d1f1`, adopted 2026-09-28. The scripts' docstrings cite
`TRIAGE_LOOP.md` sections; they mean that revision.

## Local parameterization

- Inbox `docs/FEEDBACK.md`, history `docs/FEEDBACK_HISTORY.md`, item prefix `S`.
- Tracked items are roadmap entries, `RM<n>` in `docs/ROADMAP.md`, claimed with `scripts/item-next.py`.
- The three scripts live in `scripts/`. The defaults of `triage-state.py`, `triage-archive.py` and
  `item-next.py` are derived from the script's own location, not from the working directory. The
  watcher is the kit's `drivers/git/watch-file.sh`; the gist's `watch-inbox.sh` was retired for it. `item-next.py` also counts
  an id as taken only where a heading, a reservation or a tombstone owns it, never where prose quotes it
  (a sibling's `RM267` in RM3 once set the maximum), and puts a reservation inside its anchor section,
  never in a later item's list. The gist's `wake-on-inbox.sh` and `arm-on-prompt.py` were dropped when the
  watcher moved under the seat.
- The watcher guards branch `main`.

## Arming the watcher

Nothing to arm by hand. `.blueprints.toml` declares `[seat.triage.watch.inbox]` with
`cmd = "sh kit:drivers/git/watch-file.sh docs/FEEDBACK.md"`, and the triage seat's
watchdog runs it from the repo root for as long as a session holds the seat: started at punch-in,
restarted with backoff if it exits, stopped with the seat. When the inbox differs from HEAD and has sat
unchanged for 150 s, a watcher line arrives, `docs/FEEDBACK.md: uncommitted, +N -M lines against HEAD:
git diff HEAD -- docs/FEEDBACK.md`, and wakes an idle seat. Each state is said once, across a restart
too, so an item a sibling left while no one held the seat is said at the next punch-in. Off `main` it
pauses silently. The line names no item; the pass's `--pending` does. The watchdog reads the table only when it starts, so a
change to it takes `watchdog.py stop triage` and `/seats:seat triage` again.

Never start a watcher from the agent as well: nothing stops a second one on the same file, and two
wakers ring twice.

## Intake from Zulip and from other sessions

A report can arrive as a message instead of a section: from another session on this machine, or on
Zulip, a direct message to `just-vcf/triage` or a mention in `agent-chat`. The consumers' sessions
(just-prs, just-dna-lite) hold no seats, so they have no Zulip account. For them the inbox file and
session-to-session messages are the routes; Zulip carries a person's report, or a seated repo's. It is
filed before anything else is done
with it, because the inbox is the state and a message is only delivery.

1. React to it (`watchdog.py react ID`), and read it whole (`watchdog.py read ID`) if it was cut.
2. Append it to `docs/FEEDBACK.md` as `## S<n> — <title>`, `<n>` from `./scripts/triage-state.py --next`.
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

## What a pass may do unattended

Triage, reply and archive in this repo, committing as it goes on `main`. Fixes that fall under
"patch" or "minor" may land in the same pass. Anything "major", and anything that edits the charter,
waits for the user.
