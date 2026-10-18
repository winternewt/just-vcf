# Lead seat — notes

State as of 2026-10-19, after the first seating (2026-10-18 to 19). Rewritten to the present before
each reset.

## Open

1. **The data cache holds whole private genomes** (RM268): eviction is owed.
2. **The leak gate is on** (2026-10-18, the person's yes), installed by the kit's `bin/leak-gate`
   (v0.1.68): a hook in `.git/hooks/` that runs the pinned kit copy's gate, which `move-to-tag`
   re-points; `leak-gate status` says which copy. `.leak_words` (git-ignored, 4ff8393) holds the VCF sample names of the
   private genomes, read with `vcf_samples` and written without printing them; DeepVariant's
   `default` was dropped. Folder names are left out, as they are in tracked prose already. The kit
   built the install from agent-seat-blueprints#83, filed on the person's word. Not covered: a name in a file
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

The roadmap was read against the code on 2026-10-18 (da410cc): nothing open had shipped whole; the
status lines say what exists. Schema at 0.0.4.

- **RM270** (split multi-allelic rows): B, the flags, shipped (1559215). C, withholding the uncertain
  `0`s slot by slot, waits on desk #45 (default: B only). D, the re-join, goes with RM6. The detector's
  false-positive count on the two private gVCFs was running at the commit; its result goes into RM270.
- **Desk #43**: the order of work (default: the roadmap pass first), which this sitting followed.
- **Open, not started**: RM269 (karyotype finding and caller ploidy model, a schema call), RM2's next
  loaders, RM3, RM4, RM6, RM268 (cache eviction).
