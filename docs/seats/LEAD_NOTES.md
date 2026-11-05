# Lead seat — notes

State as of 2026-10-19, after the first seating (2026-10-18 to 19). Rewritten to the present before
each reset.

## Open

0. **One or two builders at a time** (desk #103's default, again from 2026-11-01; the person's "Turn of
   parallelism until Oct 31" of 2026-10-26 has lapsed).
0. **Merge audits run on two outside model families** (the person, 2026-10-26: "For audit runs, switch
   pasing to using grok4.7 via cursor + sol6.1 via codex subagents"), from audit 13 on (audit 12 was already
   running on Opus): Grok 4.7 through Cursor's CLI (the `cursor-grok` agent, `--mode ask`) and GPT-6.1 Sol
   through Codex (`codex exec -m gpt-6.1-sol -s read-only -p secrets-denied`), each on a `git archive HEAD`
   copy of the tracked files after the kit's `repo-prep`, never the working tree. The lead checks every
   finding against the code before a fix builder gets it. Both runs together count as the one item at a time.
0. **0.1 release prep is decided** (desk #86, the person, 2026-10-26 23:51Z: "A: 1, 2, 3 / 0.1 release after 2,
   pt 2 interview goes to 0.2 and doesn't block 2 and 3"). Option 2, the CLI, shipped as f00bb1f; the Part 2
   interview was held 2026-10-23. Option 3 is the prep: CHANGELOG's Unreleased cut as 0.1.0, a README pass,
   and a check that the built sdist and wheel carry no private-derived document. The publish, and the
   `v0.1.0` tag that would mark it, stay the person's. Brought up to date 2026-11-05: the
   Unreleased entries folded into 0.1.0 with schema 0.0.22 (d9c3c43), the README names coverage, karyotype
   and the consensus (c8f891f), `tests/test_package.py` green on that code.
0. **Audit 14 done** (2026-11-03, b58f033..bf4b4f4, Sol via Codex and Grok via Cursor): six findings, all
   confirmed. Four merged as 040fcea (karyotype per GT slot and array letters, rejoin leaves an unstated split
   site as written with F34 filed, normalize keeps the matched duplicate, query's variant_class). The fifth,
   somalier, is on branch `fix/audit14` (eb575f8, pushed, unmerged): it feeds somalier just-vcf's site
   genotypes, which ends its independence; blocked: desk #178, default no change (the path errors on WGS
   facets today: `--unknown` is a `relate` flag, not `extract`'s). Also unchecked there: somalier's hg38
   sites say `chr1`, we write `1`. Next audit after 6 more merges.
0. **Shipped since, kept out of Open**: F23 as schema 0.0.22 (2026-10-27); audit 12's fixes and the loosened
   own-sample test, merge 6eed88a (2026-10-26). The roadmap's stale "awaiting the lead's merge" lines (RM274,
   RM276, RM5) corrected 2026-11-03.
0. **Done 2026-11-05: RM287 merged as 08a94c5 (merge 5 of 6 before audit 15)**: `calls.imputation_r2` (minimac R2,
   else Beagle DR2) and `imputation_empirical_r2` (ER2), schema 0.0.23 with its migration step, not read by the
   grade; `state_imputed_from` and the `imputed_from` link (`links.STATED_RELATIONS`). Calls for veto, in
   docs/SCHEMA.md: two columns by quantity, not by key; R2 and DR2 that differ on a row withheld; the tool in
   `provenance.step` and the panel a `derivation_inputs` entry; `derived_from` set on a facet read from a file.
   **The lead's follow-up, for the audit-15 fix builder**: `Person.of` now treats a chip and its stated
   imputation as one callset and never runs `same_person` between them, so a wrong `imputed_from` statement
   (a chip from another person) goes unchecked; the more complete behaviour checks the stated pair once and
   keeps them one callset against everyone else. Found, not fixed: `source._source` does not recognise Beagle's
   `##source` (caller null on Beagle files); a store already holding the imputed facet ignores a later
   statement (documented). **The changelog has no Unreleased section**: every merge until 0.1.0 is published
   adds its lines at the top of the 0.1.0 section, and the README is checked against each merge. Every roadmap
   item carries a `**State:**` line since ea6a199, read by the plan watcher (8c0314f); the `## Items` heading it
   flags stays, since `item-next` anchors its rows on it (kit #281).
0. **Done 2026-11-05: RM276 closed by the lead's call, moved to history** (no builder, not a merge): a migrated
   VCF facet's `method_steps` and `gl_source` stay null, told apart from "not a VCF" by the manifest's
   `migrations` record and `loaded()`'s warning; an own value (a `not_recorded` enum member, a marker step) is
   rejected in docs/SCHEMA.md "Store migration", for the person's veto.
0. **Building 2026-11-05: RM282's decompose** (Opus builder in a worktree; merge 6 of 6, then audit 15): the
   lead's call is that a symbolic `<DEL>`/`<DUP>`/`<INV>` meets another caller's breakends as truvari's default
   `decompose` does, the row saying so, and a lone breakend near a symbolic SV reads `other_alleles`. Filed the
   same day: desk #208 (a copy-number allele against `<DEL>`/`<DUP>`, default strict) and #209 (a callset
   recording its BAM, default no change). At merge, point RM282's State at #208.
0. **Done 2026-11-03, nothing running**: RM5's consensus merged as 6b5503c (merge 4 of 6): one panel's methods
   are one witness, averaged; only panels the caller declares disjoint multiply; its calls for veto in
   docs/SCHEMA.md "Combining ancestry estimates" (disjointness never inferred, every reading weighed alike,
   arithmetic mean). No decided buildable item is left that needs no person's data. RM6's lift across a
   gap of one chain merged as 855ab49 (merge 3 of 6; `GAP_LIMIT` 50 and "the gap does not follow the person",
   the gap's target bases written as context, are its calls for veto in docs/SCHEMA.md "Lifting" 2b). RM278's thresholds file and quantized
   runs merged (merge 2 of 6), its veto list in docs/SCHEMA.md "Coverage facets" (a contig's last base its own
   bin, unlike mosdepth; a headerless five-column BED ending in a whole number refused). RM269's person-level roll-up
   merged as 6bebc0c (merge 1 of 6 before audit 15); its builder's nine calls beyond the pinned ones (an X or
   Y line from a withheld callset still contradicts, a torn callset withholds the person, `x_copies` weighed
   too, the member saying most stands for a callset, and five more) are in its commit df943ad and in
   docs/LAYERS.md layer 4 for the person's veto; measured over estimated stays open.
0. **Imputation** (the person asked 2026-11-03 whether it belongs here or downstream): desk #180, default
   downstream, nothing built here. The read-back half is RM287 (per-site R2/DR2, an `imputed_from` link, the
   way out checked against a server's input rules), blocked on #180, yes after a present.
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
- **2026-10-24 evening, paused on the person's word** ("pause after this item I restart"). All eight layers at or
  past 90% (docs/LAYERS.md); desk #100–#103 filed. Audit 5 (1a45d09..cc08843) done and its four findings fixed
  and merged. Open follow-ups from that fix, none started: HQ declared Integer with an unparseable value
  (`30.5,40`, `x,40`) still fails the whole 4.2 load (polars-bio sometimes names the field, so a retry dropping
  it is possible: a design choice); GQ declared Float under 4.2 truncates silently; PSO/PSQ declared
  `Number=1` with two values fail; FT declared Integer under 4.2 fails; CGI's `_number` nulls a non-integer
  score silently. Next audit after six more merges. The CGI var builder's worktree is kept for its
  `data/output` runs (copied to main's `data/output/`).
- **Flaky Privé test:** fixed in 71fe6fd (the inner join's row order); likely cause, failure message not kept.
