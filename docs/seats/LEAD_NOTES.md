# Lead seat — notes

State as of 2026-10-18; no session has held the seat yet. Rewritten to the present before each reset.

## Open

1. **The data cache holds whole private genomes** (RM268): eviction is owed.
2. **The leak gate is on** (2026-10-18, the person's yes): `core.hooksPath` points at the v0.1.67
   kit copy's `.githooks/`. `.leak_words` (git-ignored, 4ff8393) holds the VCF sample names of the
   private genomes, read with `vcf_samples` and written without printing them; DeepVariant's
   `default` was dropped. Folder names are left out, as they are in tracked prose already. The path
   is versioned: after a kit move, re-point it, or it keeps the old gate. The kit owes the real
   install (agent-seat-blueprints#83, filed on the person's word). Not covered: a name in a file
   name rather than a VCF header.

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
