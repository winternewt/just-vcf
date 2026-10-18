# Consumer triage loop — just-vcf

This repo runs the document-as-inbox triage loop published at
<https://gist.github.com/winternewt/54b94bda01812be937b892146d1bb254>. The gist's `TRIAGE_LOOP.md` is
the full pattern: the algorithm (Step 0 through Step 5), the five-verdict routing, and every gotcha
found by running it. This file holds only what is local to just-vcf.

The session that runs this loop is the **triage seat**, primed from
[seats/TRIAGE.md](seats/TRIAGE.md). That file says what the seat owns, how it starts and what it never
does; this one holds the local rules the pass reads.

**Adopted through gist revision `96b6953362bb4854e05db04334086c43c9eeb170`, 2026-09-30.** That line is
the baseline. Update it whenever an adoption lands. This revision brought the Claude Code arming,
`wake-on-inbox.sh` and `arm-on-prompt.py`; the four existing scripts were identical to it apart from our
path defaults. Before it, `b063b65ef777cf06002e6fa14953d6503692d1f1`, adopted 2026-09-28.

## Local parameterization

- Inbox `docs/FEEDBACK.md`, history `docs/FEEDBACK_HISTORY.md`, item prefix `S`.
- Tracked items are roadmap entries, `RM<n>` in `docs/ROADMAP.md`, claimed with `scripts/item-next.py`.
- The four scripts live in `scripts/`. The defaults of `triage-state.py`, `triage-archive.py`,
  `item-next.py` and `watch-inbox.sh` are derived from the script's own location, not from the working
  directory, and the watcher's `RUNBOOK` default is `docs/seats/TRIAGE.md`. Those are the only changes
  from the gist copies. The gist's `wake-on-inbox.sh` and `arm-on-prompt.py` were dropped when the
  watcher moved under the seat.
- The watcher guards branch `main`.

## Arming the watcher

Nothing to arm by hand. `.blueprints.toml` declares `[seat.triage.watch.inbox]` with
`cmd = "bash scripts/watch-inbox.sh"`, and the triage seat's watchdog runs it from the repo root for as
long as a session holds the seat: started at punch-in, restarted with backoff if it exits, stopped with
the seat. Each settle arrives as a watcher line, `FEEDBACK.md settled: S<n>(new) … — triage per
docs/seats/TRIAGE.md`, and wakes an idle seat. The watchdog reads the table only when it starts, so a
change to it takes `watchdog.py stop triage` and `/seats:seat triage` again.

The watcher is a singleton per inbox: a second one started for the same file replaces the first. Never
start it from the agent as well; two wakers ring twice.

## Syncing with the gist

```sh
GIST=54b94bda01812be937b892146d1bb254
curl -sf "https://api.github.com/gists/$GIST" \
  | python3 -c 'import json,sys; print(json.load(sys.stdin)["history"][0]["version"])'
```

The same sha as the baseline means nothing arrived. A different one means fetch the files pinned to
that revision (`<gist>/raw/<version>/<name>`; the unversioned raw URL is CDN-stale for a while after a
write), diff the code with comments stripped, and adopt behaviour changes only. Run the ledger over
both documents afterwards. If a fingerprint moved, the gist changed `fingerprint()`, and the markers
need restamping, said so in the commit message.

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
