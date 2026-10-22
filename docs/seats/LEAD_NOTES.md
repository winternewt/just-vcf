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
- **F8** resolved 2026-10-22 (desk #61's default): a sites-only VCF is refused with `NoSampleError`.
- **Desk #43**: the order of work; closed unanswered 2026-10-21, so its default (the roadmap pass first) stands, as followed.
- **RM269**: the person answered desk #61 in chat on 2026-10-22, "269c" = option 3, the full label
  (XX/XY/XXY/X0), which includes the facts and the X-copy estimate. Order: an RM4 slice (build from the
  header's contig lengths, as INFERRED) → facts on the facet → X evidence (diploid-written calls only;
  a haploid GT is the caller's assertion, its own line) → the Y region set (time-boxed to an hour; if
  nothing citable exists for GRCh38, back to the desk) → calibration on PGP-HMS with self-reported sex.
- **Builder run** (the person, 2026-10-22: "use opus subagent builders + worktrees pattern on decided items
  until I cancel"): Opus builders in worktrees, the lead reviews, merges, bumps the schema and pushes.
  Wave 1: RM3 (REF lookup, FASTA tier), RM4 + assembly label and seqcol, RM2 VCF side (source detector,
  dialects), RM2 array loaders. Next waves: RM269 build, RM6 writer, then the rest of Part 1 of
  docs/MVP_READINESS.md. Part 2 waits for the interview.
  Merged: RM4 (ee4d5d1, schema 0.0.5), RM2 arrays (53600af; 18 of 39 formats load). Wave 2 running: RM6
  writer, RM5 ancestry, RM269 evidence (handed the depth result: X/auto 0.42–0.50 vs 0.89–1.00, X-degenerate
  Y ALT DP/auto 0.56–0.82 vs 0.35–0.40). Open from the array builder: `dropped_fields` is `[]` for the new
  table formats and null for the rsid tables; aligning them changes values on files that already load.
- **Open, not started**: RM2's next loaders, RM3, RM4, RM6, RM268 (cache eviction).
