# Lead seat — notes

State as of 2026-10-19, after the first seating (2026-10-18 to 19). Rewritten to the present before
each reset.

## Open

1. **The data cache holds whole private genomes** (RM268, shipped 731d38c): stale versions are swept on each
   step, and `clear_cache(source=...)` removes every copy of a deleted genome. Entries written before the
   index (REWRITE_VERSION 7) are swept as stale.
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
  Merged next: RM2 VCF side (bd5c02e, 0.0.6), RM3 (db67f83, 0.0.7). Running: RM6 writer, RM5, RM269
  evidence, RM2 rest (CGI native, fbvs key=value). Held for the person or later, from RM3: a
  caller-asserted `assembly=` in `resolve_ref` (a label would stand in for evidence: declined for now);
  the `[reference]` extra's pysam backs nothing since RM3 (removing it is a surface change); UCSC
  `chr*_random`/`chrUn_*` names do not meet Ensembl accessions (1% of GRCh38 rows unread); one PGP
  VCF under 150 MB raised ValueError in `load` (not yet identified).
- **Desk #72** (who is who in felis, ksuhaster, newton; default: unknown, never paired) and **#73** (may
  Privé's ancestry files be redistributed; default: local only), filed 2026-10-22 when the person went
  unattended. `docs/ANCESTRY_PANELS.md` stays out of CLAUDE.md's read-first list (the lead's call: it
  matters only to ancestry work). Merged since: RM5 (4c6c9c1), RM6 writer (a1f3317), F10 + SV END
  (ec95596), RM6 normalize (9469c6b).
- **Wave 3** (2026-10-22): merged F9 (2ff9b00; strand finding renumbered F13), RM269 evidence layer
  (7384428: XY/XXY or unknown; XX and XO wait for a coverage facet), RM271 VRS ids (b41bda5; 19,000
  alleles equal vrs-python; `vrs_allele_ids` follows SCHEMA's `[ref] + alts` layout, not one per ALT).
  The lead decided RM268's outline (cache index, stale-version sweep, `clear_cache`, no LRU). Running:
  RM2 CGI loaders, RM6 follow-ups, RM268, RM6 liftover, RM4 chip detection + F13. Open from RM269:
  freeBayes males read an X het share of 0.38–0.40 (inside the unknown band; depth says one X);
  two Helix exomes read `wgs` (RM4); a starling upload fails in polars-bio ("invalid BGZF ISIZE").
- **Merged since wave 3:** RM6 follow-ups (16277e8: split tag read back, dedup on what a row states, ALT
  `-` anchored, `<*>` dropped on split), RM268 (731d38c). Running: RM2 CGI loaders, RM6 liftover, RM4 chip
  detection + F13. Open: should `ids` leave the dedup key (one upload keeps 9,213 same-allele pairs
  because of it); `resolve_ref`, `normalize` and `mint_vrs_ids` rewrite their spill, so running one twice
  on a facet breaks the first result's frame.
- **Wave 4** (2026-10-22): merged RM6 liftover (047f0f2, 0.0.8), filtered genomes (3dfeaad, 0.0.9), spill
  safety (59f2a58), RM2 rest: CGI native + key=value (aa095fd; 25 of 39 formats), RM4 chip + F13 (9d18bcb,
  0.0.10). Running: cache safety for load/liftover/CGI (also: the CGI and key=value work dirs were outside
  RM268's index, so `clear_cache` could not remove them), open loader findings F1–F4/F7 + Omni2.5 + the
  starling BGZF failure, RM6 remainder (tabix, gVCF output, minus-strand lift indels). Open for later:
  vcf_cgi's `absence_means` (vcfBeta's silence is a reference call: `blocks` understates it, a new enum
  member is a schema change); the empty allele is `""` in CGI and `-` in Genos; just-prs lists AncestryDNA v2
  as GSA v3, but the 16 stated V2.0 exports overlap OmniExpress (0.64–0.67) far more than GSA (0.26–0.28):
  a FEEDBACK item for just-prs once checked.
- **Wave 5** (2026-10-22/23): merged cache safety for load/liftover/CGI (6a68c5b), loader findings (602f27c:
  F2, F7, F15 fixed; F1 needs a dependency, desk #83; F3 and F4 surfaced), RM6 remainder (8ca6f40: tabix,
  gVCF output, minus-strand lift indels), follow-ups (6ce13b0: lifted rows realigned, DamagedBgzfError,
  F4 split against the reference; F16 resolved, F17 open). Running: F3 second facet + chip by marker
  overlap. Desk waiting: #72, #73, #80, #83. Part 2 of docs/MVP_READINESS.md waits for the interview.
- **Review of 76b0a2d..main** (2026-10-23, read-only Opus): queued for one fix builder once the consistency
  bundle merges (it touches the empty allele): (1) CGI native and key=value facets lack `variant_classes` /
  `content_filter` (merge seam; CGI's unanchored `""` alleles must count as indels first); (2) an fbvs merge
  record without INFO/SF lands in no facet, uncounted; (3) CHANGELOG still says F3 is open; (4) a non-damage
  BGZF read failure walks the file twice; (6) the appended-export facet takes the first export's build as
  stated; (8) CGI native `wgs` only by `#TITLE`, no `classify` fallback. Lead's decision on (5): the stale
  sweep removes only entries of an older version, never a newer one (two installs share the user cache).
  (7) REWRITE_VERSION: no stale copy reachable; bump on the next rewrite-rule change.
- **Builder run paused, 2026-10-23** (decided work done): merged consistency decisions (36dbecc, 0.0.11)
  and the review fixes (2aa5436). 1,046 passed. Desk #86 asks what next (interview Part 2, CLI, 0.1 prep,
  or stop; default stop). Also open: ruff is not a dev dependency (about 200 findings in the tree); the
  local `.env` samples dir has no PGP-HMS, so the opt-in PGP tests skip unless run with
  `JUST_VCF_SAMPLES_DIR=/mnt/krypton/Genomes/genomes`.
- **Standing rule (the person, 2026-10-23, approving a side agent's note):** a read-only merge audit after
  every six merged builder branches, because a hand merge once dropped a per-facet field from two loaders
  while every test passed. Last audits: 2aa5436 (after wave 5) and one running over 2aa5436..5b6ec3b. A
  guard test (every format's example loaded, no facet column null without a stated reason) is being built.
  Audit 2 (2aa5436..5b6ec3b) found: the CLI printed unredacted header text (a path in `##reference=` or a
  version line), `chip_absence` mixing facets, `links` dropped by `Person` and the CLI's facet pick,
  `rejoin_split` dropping SV columns, liftover swaps leaving per-ALT copy numbers, stale version notes.
  All seven fixes are with one builder, plus conftest defaulting the opt-in variables to empty (the
  auditor's own first test run read `.env` and could have opened local genomes for a few seconds; its
  output was never read and nothing left the machine).
- **Open, not started**: CLI (MVP row 18; not decided), the rest of Part 2 of docs/MVP_READINESS.md (the
  interview).
- **2026-10-24, the 7-day quota at 100% (resets 23:59 CEST):** no new builders until then. Two runs were
  left to finish: the consensus calibration builder (rates by confidence band, commit 78d6ede on its branch,
  its real-data check written) and merge audit 4 over 01fca0b..ad2f185 (report to the scratchpad's
  audit4.md). `tests/test_facet_gate.py` is gone from main; `tests/test_load_seams.py` is the one guard.
- **Audit 4 (01fca0b..ad2f185), done 2026-10-24:** no hunk lost in a hand merge. Fixed in 1a45d09: a
  consensus facet stood as an independent source in `Person`/`check_pedigree`; the consensus ploidy test
  built a pre-F22 atlas row and guarded nothing; SCHEMA.md history (orphan heading, stale bump wording).
  Calibration merged at 86daede; the rebuilt atlas is in main's `data/output/`. Next audit after six more
  merges.
- **2026-10-24 00:18 CEST, quota reset; building resumed** on layer 7 (82%), the one layer under 90%: one
  builder builds desk #94–#96 as non-default options (defaults unchanged, so an answer is a one-line
  change); one measures why the 60+ phred band is 2x overconfident and fixes it if no decision is needed.
- **Flaky test, 2026-10-24:** `tests/test_ancestry.py::test_prive_proportions_of_each_group[EUR]` failed once
  in a full run and passed alone and in the next full run. Its inputs are seeded (seed 11), so the
  nondeterminism is elsewhere (solver, thread order, or another test's state). Unexplained; a builder's.
