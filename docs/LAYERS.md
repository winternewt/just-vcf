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
| 7 | Graded evidence across sources | 88% | `genotype_evidence` with the PL/GL source; the consensus posterior weighs each row by its confidence band's measured rate (held-out SNVs: 1.15x predicted against observed, from 1.9x; the 60+ band's 2x overconfidence was freeBayes 0.9.15 borrowing another facet's bands, and closes to 1.02 with its own, docs/CONSENSUS.md); desk #94–#96 built as options |
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
- **Tooling, desk #91:** ruff as a dev dependency (about 200 findings in the tree).
- **Housekeeping for the person:** a builder left `stash@{0}: a94f-lift-check` (an old copy of `liftover.py`, it
  says safe to drop); dropping stashes is the person's.
