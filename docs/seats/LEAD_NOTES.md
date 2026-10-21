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

- **RM270** (split multi-allelic rows) shipped 2026-10-19 and is in the roadmap history: B, the
  flags (1559215, e82e4af, the join placement fixed in 5f08965), and D, `rejoin_split` (2e84f1f), the
  person's choice on desk #45. Checked on a real genome split the bcftools way: all 88,713 sites with an
  ALT come back as the unsplit file has them. Not run: a file split by bcftools itself, Hail or GATK.
- **F8** (a sites-only VCF does not load, our spill, older than RM270): open in docs/DOGFOODING.md.
- **Desk #43**: the order of work; closed unanswered 2026-10-21, so its default (the roadmap pass first) stands, as followed.
- **Desk #57**: start RM269 now (survey, critic, counts on the local genomes, then the design as its own
  desk question)? Filed 2026-10-21 after the chat question went unanswered; its default is to start, so
  the lead works on under it. No schema change ships before that second question is answered.
- **Open, not started**: RM2's next loaders, RM3, RM4, RM6, RM268 (cache eviction).
