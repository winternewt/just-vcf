# Lead seat — notes

State as of 2026-10-18; no session has held the seat yet. Rewritten to the present before each reset.

## Open

1. **The data cache holds whole private genomes** (RM268): eviction is owed.
2. **The leak gate**: the person said yes on 2026-10-18, with the real names of the people under
   `/data/genomes/` only, not their folder names (those are already in tracked prose). The lists are
   git-ignored (4ff8393). Waiting on the person to write `.leak_words`; then the lead sets
   `core.hooksPath` to the pinned kit copy's `.githooks/` and tells triage that its commits are
   gated. That path is versioned, so after a kit move it keeps the old gate until it is re-pointed;
   the kit owes the real install (agent-seat-blueprints#83, filed on the person's word).

## Tools, settled at the first seating (2026-10-18)

- **Kept**: `scripts/pgp_leech.py`, `pgp_participants.py`, `sniff_corpus.py`, the repo's data tools;
  nothing in the kit overlaps them. The lead runs no watcher of its own.
- **Taken up**: the kit's `bin/mutate`, run from the pinned kit copy, for every new check (the house
  rule as a tool). Nothing to install.
- **Declined for now**: `bin/merge-changelog`. Every seat here commits on `main` in one clone, so the
  changelog never meets a three-way merge; the driver line would also pin a path into the versioned
  plugin cache, which moves with every kit release. Take it up when build branches start merging
  back with changelog lines of their own: one `.gitattributes` line plus the local `git config` line
  in its docstring.

## In flight

Read from the last sittings: RM2 (loaders) and RM3 (REF design) open; the schema at 0.0.3; RM269 (a
person's karyotype finding, from S1) filed by triage for the lead's schema call. Check the roadmap
against what exists before taking these as current (the role, step 1).
