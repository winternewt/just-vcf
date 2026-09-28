# Consumer triage loop — just-vcf

This repo runs the document-as-inbox triage loop published at
<https://gist.github.com/winternewt/54b94bda01812be937b892146d1bb254>. The gist's `TRIAGE_LOOP.md` is
the full pattern: the algorithm (Step 0 through Step 5), the five-verdict routing, and every gotcha
found by running it. This file holds only what is local to just-vcf.

**Adopted through gist revision `b063b65ef777cf06002e6fa14953d6503692d1f1`, 2026-09-28.** That line is
the baseline. Update it whenever an adoption lands.

## Local parameterization

- Inbox `docs/FEEDBACK.md`, history `docs/FEEDBACK_HISTORY.md`, item prefix `S`.
- Tracked items are roadmap entries, `RM<n>` in `docs/ROADMAP.md`, claimed with `scripts/item-next.py`.
- All four scripts live in `scripts/`. Their defaults are derived from the script's own location, not
  from the working directory. That is the only change from the gist copies and must not be overwritten
  when syncing.
- The watcher guards branch `main`.

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
