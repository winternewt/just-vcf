# Lead seat — notes

State as of 2026-10-18; no session has held the seat yet. Rewritten to the present before each reset.

## Open

1. **The data cache holds whole private genomes** (RM268): eviction is owed.
2. **The leak gate** (the person's call): the kit's pre-commit gate and `bin/leak-check` read a
   git-ignored `.leak_words`; here there is none, so a commit or a kit report is never checked for a
   named person. A list of the names under `/data/genomes/` would gate both; the gate then refuses
   every commit while the list is missing.

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
