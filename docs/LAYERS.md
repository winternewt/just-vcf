# Layers of the extended shape: progress and blockers

The person's goal of 2026-10-23: build unattended, todos first, then take repeated swings at layers 1–8
until each is at 90% or every layer is blocked by a decision only the person can make. Each blocker is
written here, under its layer, the moment it is met, and the work moves to another layer. The lead
seat keeps this file; percentages are the lead's estimate against what the layer would need to be
usable for one person's facets, and say what they rest on.

| # | Layer | Estimate | Rests on |
|---|---|---|---|
| 1 | Reading facets in | 90% | 32 of 39 formats load (plus SV/CNV, PSL); the seven left are decisions (#83 Deflate64, #93 Oxford .gen), chip manifests (Geno 2.0, headerless Final Report), no positions (Promethease) or not facets (CGI summary, deCODEme info) |
| 2 | Provenance per facet | 90% | caller and version, vendor, callset type, chip, content filter, ploidy lines, subject, derivation inputs, reference, aligner, the method chain (`method_steps`), coverage facets |
| 3 | Identity (comparing facets) | 92% | REF, VRS, assembly inference incl. GRCh36/CHM13, scaffold aliases (F19: a PGP GRCh38 VCF reads 100%); `_fix` patches need a patch store |
| 4 | A person | 90% | subject; Person with links and relatedness as its default check; a supplied GA4GH pedigree checked against measured bands (agrees/disagrees/cannot_tell, never corrected); exome variant-only sites read as unknown; the durable store; somalier waits on desk #88 |
| 5 | Agreement measured | 90% | `concordance` (span-guarded VRS pairing, haploid vs a/a a match with ploidy flagged); the method-error atlas incl. the CGI pairs by VRS id and by each side's confidence band; the consensus posterior (RM280) |
| 6 | Facet graph and liftover | 90% | links frame, liftover, normalize, rejoin_split |
| 7 | Graded evidence across sources | 90% | `genotype_evidence` with the PL/GL source and CGI's allele scores (`calls.hq`, so a CGI `var` facet is banded too); the consensus posterior weighs each row by its confidence band's measured rate (held-out SNVs from 10 phred up: 0.85 observed/predicted with the banded var facet, 1.01 with GIAB's easy/difficult partition as regions; docs/CONSENSUS.md); desk #94–#96 built as options; what is left is decisions |
| 8 | Export | 92% | VCF/gVCF with tabix, VCF 4.5 PSL/PSO/PSQ, CisPhasedBlock export naming what withheld it, VCF Zarr |

## Blockers waiting for the person

Each names its layer, the decision, the options, the lead's recommendation, and what the work does
meanwhile. The desk carries the ones filed there.

- **Layer 1, desk #83:** how to open Deflate64 zips (an optional extra, an external program, or refuse).
  Meanwhile they are refused with an error that says why.
- **Layer 1, F17:** Living DNA indel REF needs the Sirius chip's site alleles, and no manifest is public
  or on disk. Needs a source the person may have, or stays withheld.
- **Layer 4, desk #72:** who is who in felis, ksuhaster and newton. Meanwhile no file in those folders is
  assumed to be another's person; the atlas pairs only files that share one alignment chain (felis's
  three callers on one BAM) or one folder's caller versions (newton's DeepVariant releases), each checked
  by `same_person` (the lead's call: the four open facts concern other pairings).
- **Layer 5/7, desk #80:** merge duplicate rows that differ only in record id. Meanwhile both stay.
- **Layer 3, desk #89:** mint SVs as VRS 2 CopyNumber/Adjacency objects, or keep withholding (RM282).
  Meanwhile withheld (`sv_no_vrs_model`).
- **Layer 8, desk #90:** offer CisPhasedBlocks per window, so one symbolic allele does not withhold a whole
  chromosome's block (RM284). Meanwhile whole blocks.
- **Layer 4, desk #73:** may Privé's ancestry files be redistributed as a fixture. Meanwhile local only.
- **Layer 8, a reading to veto:** when a 4.4+ record's first slot carries a PSL name, the loader reads it
  as phased (`phased[0] = true`); VCF 4.5 writes an unphased allele as `.`. The alternative is null, which
  withholds every such block as `phase_unknown`. Built as the first; reversible (SCHEMA Decision 2).
- **Layer 3, the lead's call (not now):** lifting onto GRCh38 `_fix`/`_alt` patch scaffolds needs a store
  holding patches (an NCBI or UCSC full analysis set); 40,099 HG002 rows stay `target_contig_unmatched`.
  Deferred: patch scaffolds carry almost no personal-genome calls.
- **Layer 4, desk #88:** the somalier binary has no row in docs/TRUST.md, so the wrapper (`just_vcf.somalier`)
  has never been run; only the person sets a source's standing. Meanwhile relatedness uses just-vcf's own
  KING-robust estimate, which is the default check.
- **Layer 3, desk #92:** a build the file does not state (CASAVA, SOAPsnp, Big Y, Nirvana): let the caller
  assert it (`assembly_label=`, basis `asserted`, checked against the reference). Meanwhile those facets get
  no REF check and no VRS ids.
- **Layer 1, desk #93:** where an Oxford .gen file's chromosome comes from. Meanwhile refused.
- **Layer 4, for the person in chat (not filed: the question itself is derived from a genome):** whether
  one alias's two files are one person. Meanwhile that pair is left out of the atlas.
- **Layer 7, desk #94:** grade a gVCF reference block by the lower of its GQ floor and its likelihoods
  (recommended), by its likelihoods (now), or by GQ. All three are built (`block_grade=`, `BlockGrade`);
  the answer is a one-line default change. Meanwhile likelihoods. Open beside it: whether `gq_on_reference`
  keeps covering calls only (today `genotype_evidence` grades a GQ-only block by GQ under it, the consensus does not).
- **Layer 7, desk #95:** use GQ as the confidence of reference calls by default for single-sample callers
  (96.6% agreement measured). Meanwhile off, with the option built and tested through `consensus`.
- **Layer 7, desk #96:** a population prior from a public frequency catalogue. The machinery is built
  (`AlleleFrequencyPrior`: a caller's table, Hardy-Weinberg priors, the prior named per row in
  `support.prior`); which catalogue, and whether frequencies of exactly 0 or 1 are allowed, are the
  person's. Meanwhile the flat prior, named on every row.
- **Layer 5, the lead's calls on RM280 (built as recommended, listed for veto):** a pair's rate split by
  triangulation where three methods were compared, else half; a GP row read as its call; an unmeasured
  method gets the worst measured rate in its class; support counts in a sidecar frame; facets from the same
  reads (freeBayes and VarScan2) treated as independent, which overstates agreement (documented).
- **Layer 5/7, the lead's calls on band calibration (built, listed for veto):** a row's whole error is its
  band's measured rate, not its own phred floored at the method's (the data ruled the floor out both ways:
  CGI rows claiming 20-30 phred err about 4%, freeBayes rows past 100 phred err 0.05%); with one pair,
  band rates are fixed only up to a constant moved between the two methods, split evenly where both sides
  sit in their most confident band; a band needs 1,000 alleles to be used (`MIN_BAND_ALLELES`). The 60+
  band's 2x overconfidence was measured (docs/CONSENSUS.md): freeBayes 0.9.15 borrowing 0.9.20's bands,
  whose confidence scale differs (0.9.15's GLs floor at -10); CGI and the array were ruled out, the even
  split is second order. Whether to borrow bands across versions at all is the person's call, three
  options there; meanwhile borrowed, with a warning naming the facet.
- **Layer 7, which regions weigh an allele (desk #100):** today an allele takes the highest rate among
  every GIAB class holding it, and the classes nest, so the 10-30 phred bands read about 0.6 observed/predicted.
  GIAB's easy/difficult partition (`giab_stratifications(..., classes=["easy", "difficult"])`) brings it to
  0.86-1.01. The rule is written in SCHEMA, so changing the default is the person's. Meanwhile all classes.
- **Layer 5/7, one caller version's bands differ between people (desk #101):** freeBayes 0.9.15's
  bands run about 2x higher for huAF3C63 than for hu1B2689, and the atlas keys a pair by caller and version only,
  so a rebuilt var pair and the vcfBeta pair of one person collide and would count that person twice if pooled.
  Extends desk #99. Meanwhile the var atlas stays a separate file (`data/output/method_atlas_cgi_var.parquet`).
- **Layer 5/7, CGI's sample-level VQLOW (desk #102):** `carried_pass` reads FILTER only; CGI puts its
  low-quality mark in `sample_filters`. Measured 2026-10-25 on both public `var` pairs: every CGI SNV row at
  10-20 phred is a VQLOW hom-ALT call (12,888 and 19,711 rows), and the VQLOW hets sit at 20-40 (18,051 and
  26,333), so reading it removes CGI's whole 10-20 band and about 31k and 46k SNVs, and moves every CGI rate.
  Meanwhile FILTER only, as documented.
- **Layer 7, CGI's lowest band against freeBayes's top band (measured, no decision):** the 11-13% (52% on
  huAF3C63) disagreement there is genotype selection: CGI's 10-20 band is its VQLOW hom-ALT calls, freeBayes's
  60-100 band mostly its het calls (78% and 97%), so the cell is mostly "CGI hom-ALT, freeBayes het" (638 of 643
  mismatches). The fitted band rate is the band's average error (4.70% against a 5.51% mismatch share), so the
  rule stays; the bands-by-genotype and partner-conditional alternatives are rejected in docs/CONSENSUS.md.
- **Layer 7/8, the lead's calls on a VCF's HQ (built, listed for veto):** CGI's vcfBeta is read by its HQ, the key VCF
  reserves, not its EHQ (equal on all 4.2M scored rows of the public hu1B2689 file, so nothing moves today); a score on
  an uncalled allele is dropped, as the native CGI loader drops `?` slots; HQ grades only where GP, GL and GQ do not;
  a row with more than two GT alleles is written without HQ (counted in a warning) rather than declaring `Number=.`.
- **Layer 6/8, the lead's calls on hq through normalize and Zarr (built, listed for veto):** `rejoin_split` joins
  `hq` per slot, as it joins GT (the old keep-when-equal rule lost `hq` at every split heterozygous site); rows that
  disagree on one slot null only that slot; `call_HQ` sits on the `ploidy` dimension, unconfirmed against the VCF
  Zarr 0.4 spec (not on this machine), a one-string rename if wrong.
- **Layer 4, sex chromosomes (RM269), desk #135 and #136:** from a variant file alone the estimate still names XY and
  XXY though it cannot rule out an extra Y (#135, default: keep, caveat documented); a coverage facet now counts X and
  Y copies past two. Chip files read X from their own hets; their Y probes per vendor (#136, answered 2026-10-26: on by
  default for AncestryDNA, whose Y calls are measured, off for 23andMe and every other vendor; `array_y` forces either). The lead's calls, for veto: three X with two Y
  is `OTHER_KARYOTYPE`; the array two-X expectation is 3/4 (Ne_X = 3/4 Ne_A), not PLINK's 0.2/0.8; copy counts are new
  nullable `x_copy_count`/`y_copy_count`, not new `XCopies` members; a table that states no build answers `no_build` as loaded and is estimated on `resolve_ref`'s output, the result naming the build's source; at the person level (`person_karyotype`, 2026-11-03) no callset outranks another, so a measured genome and an estimated array that differ withhold the answer as `facets_disagree`.
- **Layer 1, the lead's calls on VCF loading (built, listed for veto):** GT declared Integer, Float or Character is
  read as the genotype with a warning, not refused; a VCF with no GT anywhere loads with the genotype unknown on every
  row and a warning, not refused as no calls; Float-declared integer keys (GQ, DP, PS, END, ...) are read from their
  text in the copy polars-bio reads, since no 32-bit float value can be certified whole (`REWRITE_VERSION` 15; unlike
  the GT retype, no warning names it); a PS declared but not carried, or negative, withholds the phase of every phased
  record it would have set, haploid calls included, while a record whose PSL names a set keeps it; one PGP Final
  Report (no Chr/MapInfo columns) waits for its chip manifest, HumanOmni2.5-8v1_A, which is not on disk.
- **Layer 1, more of the lead's calls on loading (built, listed for veto):** a key=value or Nirvana genotype that is
  not allele indices (`0/x`) is withheld whole (`gt`, `gt_alleles`, `phased` null, the text in `source_gt`), not read
  as `[0, null]`, which would look like a written half-call; under VCF 4.3+ a carried key left out of the copy is
  declared `Number=.` with its reserved Type there so its values are skipped unread (`REWRITE_VERSION` 16). GRCh36:
  UCSC's public hg19ToHg18 chain was downloaded to derive the Y region set (md5 checked); GRCh36 arrays are
  estimated, GRCh36 genomes answer `uncalibrated_build`; one test pins two counts read off the carried chain file.
- **Layer 1, the three loader gaps closed (4167e0e), one call for veto:** in a VCF 4.4+ file, where `source_gt`
  is withheld (P4), `0/300` gives the same row as a half-call `0/.` (`gt [0, null]`) and only the past-ALT warning
  tells them apart; the alternative is withholding such a genotype whole in 4.4+ VCFs. `READER_DEFINITIONS` mirrors the
  pinned noodles fork and needs re-extracting when the pin moves. A Flag written with a value (`DB=1`) still fails a
  4.3+ read. A single `0/x` GT fails a whole VCF load upstream (P22).
- **Layers 4/8, the lead's calls on coverage and SV export (built, listed for veto):** a coverage facet is paired
  with a callset by build, since no frame records which alignment a callset came from (whether a callset should
  carry its alignment's digest is open); several usable coverage facets on one build pair with none and are all
  named; the coverage spill is keyed on polars-bio's version as well as `COVERAGE_VERSION`. SVLEN is written from
  END on length-type ALTs (desk #148 asks whether to store it as written); a CIPOS/CIEND offset of -1 or -2 collides
  with VCF Zarr's integer missing/fill sentinels, so such an interval is written missing and counted.
- **Layer 4, the lead's calls on store migration (built, listed for veto):** a store is never migrated on read;
  `migrate_store` deletes the old version's files once the verified swap commits (`keep_previous=True` keeps them);
  on a migrated VCF facet a null `method_steps` or `gl_source` means "not read", told apart from "not a VCF" or "no
  likelihoods" only by the manifest and the `loaded()` warnings, which now carry the migration caveats too.
- **Layer 1, the lead's calls on whole-genome cost (built, listed for veto):** `reference_support` and
  `genotype_evidence_at` collect the covering records in memory instead of joining the lazy frames (a whole-genome
  gVCF answers 1M sites within 6 GB; grading stays lazy). `load` writes its spill in canonical order once, in parts
  of whole contigs up to 1M rows, so the frames read as one scan with no sort (F27; 1M lookups on hu416394's gVCF
  249 s to 4 s, load peak unchanged at 7.8 GB). The layout change is keyed by a new `loading.FRAMES_LAYOUT`, not by
  bumping `cache.SPILL_LAYOUT`, so the resolve, lift, VRS and coverage spills are not retired with it (alternative
  rejected: one `SPILL_LAYOUT` bump, which recomputes all of them). Left: `rejoin_split`'s spill still sorts per
  contig on read; hu80855C's 138M-block load not re-measured.
- **Layer 5, the lead's calls on breakends in concordance (built, listed for veto):** truvari 5.5's `bnd_match` is
  mirrored as its code runs, including its reading of CIPOS's upper offset on both sides (`bounds()` sets the key
  before testing it, an apparent truvari bug; not reported upstream, the person's call). Mates are found by the
  ends each record states, Manta's within each other's CIPOS (F28), not by MATEID, which the schema drops; a pair
  of mates is one row. Single breakends are counted as `single_breakend`, as truvari's default leaves them out.
  The `breakend` exclusion reason is gone (pre-1.0 break, in the changelog). Left for the person: whether a
  symbolic `<DEL>`/`<DUP>`/`<INV>` should meet the breakends another caller writes for it (truvari's `decompose`).
- **Layer 3, the lead's calls on breakends in a lift (built, listed for veto):** the mate lifts through the chain and
  the brackets turn by each end's own strand per VCF 4.3 §5.4 (liftoverSV lifts the mate but drops the strand, so
  never turns a bracket; Picard, CrossMap and bcftools +liftover lift no mate). The mate is named as the target
  collection names it. `cipos`/`ciend` mirror without the SV lift's swap, since a breakend's POS stays its own.
  `breakend_not_lifted` narrowed to ALTs the lift cannot read (pre-1.0 break); new `breakend_*` reasons for the rest.
  Left: F30 (`write_vcf`'s `contig_names` does not reach a mate inside ALT); an imprecise breakend landing on a minus
  strand is not moved to the left of its CIPOS interval.
- **Layer 4, the lead's calls on audit 10 (built, listed for veto):** every step's spill key now covers its step
  versions and its set of file groups, and a generation found with other groups is rewritten (a lift generation
  missing `mate_contigs` had crashed `liftover`). The cost is one-time: each step's next run writes its result once
  more (a WGS load sorts again, `alignment_coverage` reads its BAM again), the old copy kept until `clear_cache`.
  Alternative rejected: keying the lift alone, which leaves the other six steps one forgotten bump from the same
  crash. Superseded generations are not deleted on a run (an older session may still read one within the grace);
  the docstrings now say so. F30 resolved (merge 0bc11c1); F31 logged, the lift's naming rule, not decided.
- **Layer 1, the lead's call on a genotype table's variant classes (not built, listed for veto):** they stay null.
  An array writes no REF or ALT at load, so its `AA` is hom-ref or hom-alt until `resolve_ref`. Rejected: counting
  after `resolve_ref` (about an afternoon), which counts the alleles the person carries in a column that, for a VCF,
  counts the ALT alleles the file writes; and counting probe shapes at load (letters as SNVs, I/D codes as indels),
  a new meaning or a new column. Nothing that gates absence reads the field for an array.
- **Layer 1, the lead's call on a cached load's time (not built yet, listed for veto):** a cached load still reads the
  VCF twice through polars-bio (warning counts, and the staged spill whose bytes key the generation): 88 s of
  hu416394's, 804 s of the Isaac gVCF's. The generation stays keyed by the staged bytes; the warning counts are to be
  cached by the copy's digest, which saves the first pass. Rejected: keying the generation on its inputs (source
  digest, options, a loader version), where one forgotten bump serves rows an older loader wrote. *Built
  2026-10-26:* a cached load of hu416394's gVCF 86 s to 48 s. The stats key carries a digest of the loader
  modules' source as well as `STATS_VERSION`, so a forgotten bump costs one pass, never stale counts (accepted
  cost: any edit to `loaders/vcf.py` or `loaders/text.py` re-counts each file once); a second session waits on the
  entry's lock rather than counting twice.
- **Layer 2, the lead's call on F23's bounds (built 2026-10-27, listed for veto):** desk #151 decided that a source
  REF of only `N` states no base (`ref_check = unstated`, lifted like `match`). Built only where the reference
  reads a base: an `N` REF over an `N` in the reference (an assembly gap) stays `match`, as before, because the
  letters do agree and VRS ids minted over gap bases rest on it. A public CNV upload holds 36 such rows. Rejected:
  `unstated` over gaps too, which takes `match` away from rows whose letters agree. An SNV whose REF is `N` keeps
  its ALT unminted (`ref_unstated`); no file in hand has one. Since merge audit 13 (2026-10-28) the same ALT is
  also left unlifted, unwritten and unanswered by `query` (`schema.unstated_sequence_alt`). Since 2026-10-31 it is
  not compared by `concordance` (excluded as `ref_unstated`), states nothing in `consensus`, is graded
  `grade_basis = ref_unstated` in genotype evidence (a new member of a returned frame's enum, not of the stored
  schema; rejected: `no_call`, which says the source called nothing, where it wrote `0/0`), and is left out of
  the sex-chromosome ALT counts. Variant classes are counted at load, before any REF is read, so such a row still
  counts as an SNV there; no file in hand has one.
- **Tooling, desk #91:** ruff as a dev dependency (about 200 findings in the tree).
- **Housekeeping for the person:** a builder left `stash@{0}: a94f-lift-check` (an old copy of `liftover.py`, it
  says safe to drop); dropping stashes is the person's.
