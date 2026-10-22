# Roadmap

Active items only, newest decisions first. Shipped items move to
[ROADMAP_HISTORY.md](ROADMAP_HISTORY.md). Claim a new id with
`ITEM_PREFIX=RM sh "${BLUEPRINTS_KIT:?}/bin/py.sh" "$BLUEPRINTS_KIT/classes/triage/item-next.py" --note "…"`; `CLAUDE.md` §8 says where `$BLUEPRINTS_KIT` comes from.

## Items

## RM1 — Final shape and plan
Severity: blocker · Status: charter drafted, not signed · Owner: user

The shape as the user described it on 2026-09-28 is written into the draft
`docs/CHARTER_DRAFT.md`, committed so work can build against it and edited directly until signed off. The frame is the product (Arrow/Polars, CompoSeq-style). There is a minimal
many-in/one-out shape and an extended shape of per-person facets (assemblies, callers, phasings,
coverage) with graded evidence, provenance, concordance and a graph. It is GA4GH first: VRS through
gtars, seqcol through refget. The user keeps it as a draft; planning of RM2–RM6 follows from it.
Surveys: docs/SIBLING_SURVEY.md, docs/LANDSCAPE.md.

## RM2 — Loading consumer formats
Severity: high · Status: in progress (detection shipped 2026-09-29; 23andMe-layout and plain VCF loaders shipped 2026-09-30; on 2026-10-22 the vendor-table loader, the VCF `source` detector and the VCF dialects `bcf`, `vcf_no_meta`, `vcf_cgi` and `vcf_freebayes_varscan2` (from the merged VCF; not yet from the `-orig` table): 22 of the 39 declared formats load, the 9 rsid-table layouts, `vcf` and its four dialects, AncestryDNA, FTDNA (both exports), the RSID/RESULT table (MyHeritage), Color, deCODEme, Genos and the Illumina Final Report with plus-strand alleles; `freebayes_varscan2_keyvalue` and CGI `var`/`masterVarBeta` are next; Geno 2.0 and Final Reports without positions or plus-strand alleles wait for chip manifests) · Owner: lead

**Shipped:** content-based detection. `just_vcf.sniff` and the declarative registry in
`just_vcf.formats` (docs/FORMATS.md). Over the PGP-HMS participants, 8 of about 1,740 files stay unknown,
all prose or empty.

**Next: per-format loaders on polars-bio**, emitting the frames in docs/SCHEMA.md (schema v0.0.1), driven by each `FormatSpec`. Start from the
contract gaps listed at the end of docs/FORMATS.md (versioned columns, section headers, 0-based CGI
conversion, half-calls, chromosome outside the row) and from polars-bio pending fixes P1/P2. Port from
just-prs `arrays.py`, `io_utils.py` and `vcf.py`; leave out the observed-allele-as-REF encoding (charter
invariant 3). Living DNA writes indels as run-together alleles (`CCA`); Color positions are 0-based.

The multi-facet participants in docs/PGP_HMS.md are the convergence tests: one person's facets, loaded
independently, must agree where they overlap.

**Shipped 2026-10-22: the `source` detector** (`just_vcf/source.py`), naming the caller or vendor from
VCF headers, because `Family` cannot say "array or imputed content inside a VCF" and that decides
whether absent means unknown. It fills `facets.vendor`, `caller`, `caller_version`, `phasing_method`
(Eagle), the new `ploidy_statements` (RM269 half 1) and, where the header states array or imputed
content, `callset_type`. Every row of the table below is a rule, plus `##bcftools_callVersion` and
`##DRAGENVersion`; `##workflow_type=Illumina_WGS_CLIA`, `##source_version` and `##content` (Isaac and
starling) were added from the PGP-HMS headers. Decisions:

- `callset_type` takes only what the header states (`array`, `imputed`); rejected: inferring it from
  the rows, which is RM4's. RM4's inference should not override a stated value.
- `vendor` is the party the header names as delivering the file; a caller's maker is not a vendor
  (DRAGEN and CASAVA lines do not make Illumina the vendor; `##workflow_type=Illumina_*` does). Rejected:
  a vendor per caller, which would call a Sequencing.com delivery an Illumina one.
- Post-processing tools (GATK SelectVariants, ApplyVQSR, VariantFiltration, bcftools view and concat)
  are not caller evidence, and GenotypeGVCFs is HaplotypeCaller's own pipeline. Rejected: GATK as the
  caller whenever a GATK line exists, which names GATK for a DeepVariant file run through SelectVariants.
- Two callers, or two versions of one, leave the field null with one warning naming each line
  (charter 8: no winner picked). A `vcf-merge` header is not a conflict; the merged file's facets
  carry one caller each (Decision 8).
- `ploidy_statements` keeps lines verbatim with no vocabulary, as RM269's critique recommends; an
  option's path value is written `<path>`.

**Open from `vcf_cgi` (2026-10-22):** a vcfBeta lists its no-called spans and half-calls but, in the
rows counted, never a reference-called record, which suggests an absent position was called
reference. The loader keeps `absence_means = blocks` (absent is not reference) until CGI's own
documentation or a comparison with the same genome's `var` file settles it; `weak_reference` or a
stronger value would be the change.

Header evidence found while declaring the VCF formats:

| Header line | Caller / vendor |
|---|---|
| `##source=freeBayes…` + `##source=VarScan2` + `vcf-merge` | PGP freeBayes+VarScan2 merge |
| `##center=Complete Genomics`, `##source=CGAPipeline_*` | Complete Genomics |
| `##source=IsaacVariantCaller`, `##source=starling` | Illumina Isaac/starling (Veritas) |
| `##source=…CASAVA…makeConsensusVcf.pl` | Illumina CASAVA |
| `##workflow_type=Illumina_GenotypingToVCF` | Illumina array as VCF |
| `##DRAGENCommandLine`, `##source=DRAGEN_SNV` | DRAGEN |
| `##GATKCommandLine*`, `##source=HaplotypeCaller…`, `##UnifiedGenotyper=` | GATK |
| `##DeepVariant_version=` | DeepVariant |
| `##source=Platypus_Version_*` | Platypus |
| `##source=GenerateSVCandidates` | Manta |
| `##source=aengine` | FamilyTreeDNA Big Y |
| `##source=loimpute`, `##generatedby=Gencove` | Gencove low-pass imputation |
| `##source=Minimac3`, `##eagleVersion=` | Minimac3 + Eagle imputation |
| `##Genes_for_Good_version=` | Genes for Good |
| `##source=open_humans_data_processing.twenty_three_and_me` | Open Humans 23andMe conversion |
| `##source=intersection_dbSNP_23andme` | DNA.Land 23andMe conversion |
| `##pipeline=helix-*` | Helix |
| `##dataAnalysisProvider=Sequencing.com` | Sequencing.com |
| `##AnnotationEngineVersion=` + FORMAT `CGIANN_*` | Genos exome |

Not identifiable from the header: 23andMe's exome VCF, the Dante- and Nebula-named uploads (generic GATK).

## RM3 — REF lookup and padding
Severity: high · Status: open (checked 2026-10-18: nothing built; every row's `ref_check` is `not_run`) · Owner: lead

Port just-prs `reference_allele.py` (reference-panel tier, FASTA tier, FASTA build check). New work:
ALT derivation from genotype against REF, indel anchor bases, chr-prefixed FASTA support, emitting
reference-supported positions.

**Measured 2026-09-30: gtars RefgetStore as the sequence source** (gtars 0.10, Ensembl GRCh38
primary-assembly FASTA, 3.15 GB, 194 sequences):

- Build: `RefgetStore.on_disk` + `add_sequence_collection_from_fasta`, 25 s. Encoded store 1.1 GB; the
  25 primary contigs are 1.16 GB of it (2-bit, or 3-bit where a contig holds `N`). chr1 is 93 MB and
  gzips only to 67.5 MB, so compressing a hosted store buys little.
- Local lookups: about 8 µs each, for 1 base or a 100-bp window. The OS cache was warm (the store had
  just been written). 200 lookups checked against the FASTA: 0 mismatches.
- Remote (`open_remote` against a local HTTP range server, so no real network latency): opening
  fetches 57 KB of metadata. A single `get_substring` is one HTTP range request of about 20 bytes,
  1.2 ms on localhost; over the internet each one costs a round trip. `get_substrings` with many ranges
  downloads the whole contig once (chr1: one 93 MB GET) and serves every later lookup locally (3 ms
  for 500). Both modes wrote the contig into the cache directory, even after `disable_persistence()`.
- No public hosted store was found (gtars docs use placeholder URLs; searched 2026-09-30). A remote
  tier means hosting the store's static files ourselves.
- just-dna-format has no sequence access except seqrepo REST, and mints SNV ids offline from a committed
  GRCh38 accession table: a substitution's VRS id does not encode REF. Their RM267 shows a variant
  table used as the REF/position source producing wrong events that a REF check passes, which rules out
  a reference panel as the REF value here.

## RM4 — Input type inference
Severity: medium · Status: in progress (2026-10-22: a VCF facet gets `callset_type` (`wgs`, `wes`, `panel`, `imputed`, or null) and the `absence_means` that follows, and `assembly_label` from its contig lengths; left: chip generation, array build, gVCF `END` expansion) · Owner: lead

Port just-dna-lite `detect_callset_scope` (WGS vs targeted), gVCF detection, and just-prs chip
generation detection. New work: WES vs WGS vs panel, and expanding gVCF `END` blocks.

**Shipped 2026-10-22** (`just_vcf.callset`, `just_vcf.assembly`; rules in docs/SCHEMA.md Decision 5):
the WGS rule is just-dna-lite's, unchanged; exome and panel bands were added with unknown between
every pair of classes. Decisions taken, each with what was rejected:

- gVCF stays out of `callset_type`: a gVCF can be a genome or an exome, so its blocks set
  `absence_means = blocks` and its calls set the type. Rejected: a `gvcf` member, which would change
  the `pl.Enum` dtype and answer a question the blocks frame already answers.
- A region extract reads as `panel`. Rejected: a separate "extract" outcome, which the calls cannot
  tell from a panel; either way absence outside the regions means nothing was looked at.
- Sequencing classes need a read-evidence FORMAT key, and rows whose ALTs are all symbolic do not
  count as sites. Rejected: reading GP/DS as imputation (DRAGEN writes GP), and counting SV records
  (a 5-sample DELLY merge in PGP-HMS had 1.04M SV sites at 93.5% breadth and read as `wgs`).
- `sorted_name_length_pairs_digest` is its own column; `seqcol_digest` stays null without bases.
  Rejected: filling `seqcol_digest` from a names-and-lengths lookup (the public seqcol server holds
  two collections for each of UCSC hg19's and hg38's coordinate systems).
- The classifier and the inference stay internal (`just_vcf.callset`, `just_vcf.assembly`, not in
  `__all__`); the facet columns are the contract. Rejected: exporting them before a consumer asks.

**Left:** just-prs chip generation detection; array build inference (just-prs has none, it defaults
to GRCh37, which is not ported); contigs named by RefSeq accession (`NC_000001.11`), which the label
inference does not read yet; gVCF `END` expansion; an SNV-only or indel-only split of a genome
(`bcftools view -v snps`) reads as `wgs`, so its absences read as weak reference for the other
variant class too. Unmarked imputation output with read-evidence keys would read as `wgs`.

## RM5 — Population inference
Severity: medium · Status: open · Owner: lead

Port the just-prs `ancestry/` package (PCA projection, KNN, proportions, consensus) behind the
`[ancestry]` extra, without the `PRSCatalog` wrappers or the Hugging Face pull.

## RM6 — Normalization and VCF writing
Severity: high · Status: open (checked 2026-10-18: no writer yet) · Owner: lead

Reimplement chromosome naming and quality filters (just-dna-lite's versions are AGPL and not the user's,
so they are a behaviour reference, not a source) and port coordinate liftover from just-prs.
New work: left-align, split multi-allelic sites, strand resolution (A/T and C/G), dedup, sort, and a
writer that keeps FORMAT/GT (just-dna-lite's writer drops them).

When it splits multi-allelic sites, the writer must not use the convention RM270 found in pre-split
input, where a split row's `0` can mean another ALT: it writes `.` for the other ALT, as vt, Glow,
vcflib and `bcftools --multi-overlaps .` do (S2, from just-dna-lite). `rejoin_split` (RM270) re-joins
such input today.


## RM268 — Cache eviction
Severity: medium · Status: open · Owner: lead

Nothing ever removes anything from `JUST_VCF_CACHE_DIR`. It holds two kinds of entry:

- `vcf/<key>/`: rewritten copies of inputs polars-bio cannot read as they are
  (`loaders/vcf.py` `readable_path`), plus a `structural/` copy for a file polars-bio refused. The key
  hashes the source digest, the archive member and `REWRITE_VERSION`.
- `frames/<facet_id>-<SCHEMA_VERSION>/`: the calls and blocks parquet spill from `load`, and since
  RM270 the calls spill of a `rejoin_split` facet, under the derived `facet_id`. Nothing in the cache
  ties that id to its parent, so an eviction that keys on a source file does not see it.

Because the version is part of each key, every bump of `REWRITE_VERSION` or `SCHEMA_VERSION` orphans
every earlier entry, and nothing collects them. Measured 2026-10-03: after the 2026-09-30 dogfooding
runs, the cache under `data/interim/cache` had reached 32 GB in 36 `vcf/` entries. Most of it was
orphans: the two largest were plain `input.vcf` copies of 14.5 GB and 7.4 GB, from before the rewrite
switched to BGZF `input.vcf.gz`. The user deleted it by hand.

Constraints on a fix:

- A copy has to outlive the lazy frames that read it, so evicting the entry a live `Loaded` points at
  breaks a later `collect`. Orphans by version are always safe to remove; current entries need an age
  or size policy.
- The copies are whole private genomes. Deleting a source under `/data/genomes` leaves its copy in the
  cache, so eviction is also how a removed genome actually goes away.

Options, not decided: sweep stale-version entries on each load; an explicit `clear_cache()` in the
API (and in the CLI once it exists); a size cap with least-recently-used eviction.

## RM269 — Karyotype finding and the caller's ploidy model
Severity: high · Status: in progress (the person chose the full label, desk #61, 2026-10-22; the evidence layer, X copies and an XY/XXY estimate per WGS facet built 2026-10-22, checked on 27 PGP-HMS VCFs; XX and XO wait for a coverage facet, as do XXX, XYY and arrays) · Owner: lead

Motivating report: S1 (just-prs, 2026-10-11). Callers, reference panels and scoring tools each pick a
convention for male X and Y and rarely record it: a caller with a diploid-everywhere model (DeepVariant
1.1.0) writes a man's non-PAR X as `1/1` and `0/1`, a sex-aware one (DRAGEN) writes or filters it
differently, and arrays split by vendor (23andMe one letter; AncestryDNA the same letter twice, per the
`ancestrydna` format notes; that format has no loader yet, RM2). A consumer reading `gt_alleles`
cannot tell one copy from two.

Reproduced 2026-10-18 against the code, on a synthetic VCF: `##DeepVariant_version=1.1.0` in the
header still leaves `facets.caller` and `caller_version` null, and a male-X `1` loads as one allele
beside `1/1` as two, exactly as written (SCHEMA Decision 2). Nothing in the frames can say the karyotype;
SCHEMA.md's XYY stress case already lists that as open.

Two halves:

1. **Caller, version and ploidy model on the facet.** `caller` and `caller_version` exist in `facets`
   and are never filled; RM2's planned `source` detector fills them. This item adds what RM2's table
   lacks: `##bcftools_callVersion`, and the caller's ploidy model as a new nullable facet field (sex-aware,
   diploid everywhere, or an explicit ploidy), read from the header where it says so and null where it
   does not. Never guessed from the caller's name alone. Checked 2026-10-21 on the 32 local private
   VCFs (header keys and option names only): DRAGEN writes `##referenceSexKaryotype=` (a two-letter
   value, the karyotype it called with) beside `##DRAGENCommandLine`; GATK writes the ploidy in its
   recorded command lines, `sample_ploidy=2` in all 7 GATK files (the 6 GenotypeGVCFs ones also
   `--sample-ploidy 2`), which says diploid everywhere unless a per-contig ploidy is also given; the 15
   DeepVariant files carry only `##DeepVariant_version`, so their model is unknown from the header.
   The 2 CrossMap/Picard liftovers and 7 files with a bare `##source` say nothing. PGP-HMS files
   (freeBayes, VarScan2) are not on disk now and are unchecked.
2. **A person-level karyotype finding with its own provenance.** Evidence without depth: the non-PAR X
   heterozygous fraction and the presence of PASS Y calls, which separate XX, XY, XXY and X0 on a WGS
   facet. XYY needs depth, so it waits for a coverage facet. Weak or contradictory evidence gives
   unknown, never the nearest class.
   Counted 2026-10-21 on the local single-sample VCFs with PASS calls (16 files, several of them
   liftovers of one callset; aggregates only, and kept here without naming a file, since a file's
   class would name a person's sex). The non-PAR X heterozygous share of PASS ALT calls falls in two
   clusters, 0.009–0.115 (11 files) and 0.580–0.626 (5 files), with nothing between: that half of the
   evidence works. **PASS Y calls do not**: every DeepVariant and GATK file in both clusters has
   2,666–10,759 PASS Y ALT calls outside the PARs, so their presence cannot tell XX from XY or XXY.
   Y evidence needs a region rule (the male-specific X-degenerate region, without the X-transposed
   and ampliconic parts; which tools restrict it so is for the survey to check) or depth, and until it has one the
   finding cannot separate XXY from XX or X0 from XY. The build is needed for the PAR bounds, and
   `load` fills `assembly_label` only where the header states it (none of these 16 did); the count
   used the header's chrX length. A sex-aware caller (DRAGEN, per its documentation) writes one
   allele on non-PAR X and Y for an XY sample, so the X share must count a haploid call as
   hemizygous evidence, not as missing. The vocabulary is GA4GH Phenopackets' `KaryotypicSex`, used
   verbatim (charter invariant 2). PAR bounds are per assembly and are domain constants.

Critique (2026-10-22, an adversarial read of this entry; tool behaviour marked unverified is the
critic's recollection, not checked here):

- Without a working Y signal no class can be named honestly: a high X share fits XX, XXY or XXX, a low
  one XY, X0 or XYY, so "never the nearest class" leaves every answer unknown. A Y region rule (or
  depth) is a precondition for any class, not a refinement.
- A haploid GT is the caller's (or the array vendor's) sex decision, not evidence: counting it fills a
  value from the source that checks it. Keep what the caller asserted as its own line, count het/hom
  only on calls written diploid, and read a caller's conflict filter (DRAGEN `PloidyConflict`) as Y
  evidence.
- The two clusters come from a few unlabelled people; a cut-off fitted to them is a number read off a
  data dump. Calibrate on PGP-HMS (self-reported sex on the profile, unverified; a second caller
  family) or cite a published rule (plink `--check-sex` F 0.2/0.8, unverified), with a minimum site count.
- "Ploidy model" interprets tools: GATK `sample_ploidy=2` is a default that per-interval options may
  override (unverified), DRAGEN's karyotype line may be given or estimated (unverified), DeepVariant
  1.6's `--haploid_contigs` does not reach the header (unverified, consistent with the 15 local files).
  Better: the header's ploidy lines verbatim, plus observed counts of one-, two- and n-allele GTs on
  autosomes, non-PAR X, non-PAR Y and MT.
- Depth is already in the frames (`dp`, `ad`, `blocks.dp_min`): X-het allele balance is the signal for
  contamination and mosaic 45,X/46,XX, which het share alone misses.
- There is no person in the schema. Evidence belongs per facet (a struct on `facets` or a small findings
  frame) and a verdict to a function returning a result with a reason (`no_build`, `y_unusable`,
  `caller_asserted_only`, `too_few_sites`, `intermediate`, `facets_disagree`, `not_wgs`); liftovers of
  one callset are one piece of evidence, not several (`derived_from`).
- `KaryotypicSex` belongs to a Phenopackets Individual and has no place for an inference: use its values
  inside a just-vcf-named inferred field, and say the result is not cytogenetic. An XX male gives XX
  correctly (karyotype, not sex). just-prs needs the X copy count for dosage; that may be the thing to
  expose first.
- LANDSCAPE names somalier for sex checks: say why this is built here, or take its sites. Arrays may do
  better than planned (female Y probes are no-calls; AncestryDNA doubles X letters).

Constraints: list length stays the number of haplotypes the source wrote (Decision 2). The finding sits
beside the facets and is never a rewrite of `calls`: no GT is collapsed to haploid by the loader. An
array facet without Y probes, or a targeted panel, says unknown.

Legality: minor. Filling `caller`/`caller_version`, a new nullable facet field and a new person-level
finding are all additive. Where the finding lives (a new frame or columns on an existing one) is a
schema decision for the lead. just-prs builds a stopgap now (karyotype from its normalized parquet,
caller and version as parquet metadata) and moves to this when it exists.

**Depth measured 2026-10-22 (lead), on the same 16 local private single-sample WGS VCFs** (aggregates
only, no file named). Medians of the calls' own FORMAT DP, each over the median DP of PASS calls on
chr20:

- non-PAR X, PASS calls: 0.42–0.50 in the low-X-het cluster (11 files), 0.89–1.00 in the high one (5).
  Nothing between.
- X-degenerate Y (GIAB v3.6 chrY nonPAR minus XTR minus ampliconic), PASS calls carrying an ALT:
  0.56–0.82 in the low-X-het cluster, 0.35–0.40 in the high one. Nothing between. Every row with DP in
  that region, PASS or not, overlaps less cleanly (0.68–1.06 against 0.35–0.57), so the ALT-call
  median is the one used.
- One sex-aware caller's file writes X haploid; its X depth ratio (0.48) agrees with its haploid calls,
  the independent check the critique asked for.
- Liftovers that dropped FORMAT carry no DP and give no depth line.

So X copies have two independent lines (X depth and the diploid X het share) that agree on all 16
files, and Y presence has one (the X-degenerate Y depth). Y-region ALT calls with DP were in the
thousands per file. Rerun through `sex_chromosome_evidence` on the 17 local single-sample VCFs it can
load: 12 facets read one X (X depth 0.419–0.500, Y depth 0.560–0.824) and 5 two X (het share
0.580–0.626, X depth 0.892–1.000, Y depth 0.351–0.400); no facet fell in an unknown band.

**PGP-HMS calibration, 2026-10-22** (consented-public, CC0, so per-file values are kept). 27 uploaded
single-sample VCFs with DP, read from the Krypton copy; labels are the participants' self-reported
`Sex/Gender` in `PGP_Participant_Survey.csv` of the HU-PGP-DB-2021-06-18 Arvados snapshot (14 files
labelled male, 4 female, 9 unlabelled; unverified). The profile pages could not be read: the site's
certificate names another host. Shipped rules, file by file (`x` and `y` are the depth ratios):

| Participant / upload | Label | Build, caller | X evidence | Y ratio | Answer |
|---|---|---|---|---|---|
| huE0F0CE/4084 | female | GRCh38 DRAGEN | share 0.593, x 0.926 | 0.222 | two X, Y not shown |
| huF7A4DE/3936 | female | GRCh38 GATK HC | share 0.618, x 0.938 | 0.125 | two X, Y not shown |
| hu21731A/2058, huAA16BD/825 | female | GRCh37 exomes | share 0.588, 0.467 | - | `not_wgs` |
| hu6ABACE/3885 | male | GRCh37 DRAGEN | share 0.041, x 0.446 | 0.607 | XY |
| huE18D82/3910 | male | GRCh37 DRAGEN | share 0.043, x 0.474 | 0.842 | XY |
| hu0486D6/1048 | male | GRCh37 Isaac | haploid X, x 0.560 | 0.560 | XY |
| huA2692E/2816 | male | GRCh37 Isaac | share 0.058, x 0.552 | 0.862 | XY |
| hu3BDC4B/3479, /3481 | male | GRCh37 GATK HC (bwa, bowtie2) | share 0.016, 0.012; x 0.484, 0.500 | 0.484, 0.433 | one X, Y not shown |
| hu2B76CC/2104, huF06AD0/1506 | male | GRCh37 freeBayes (no AD) | share 0.396, 0.378; x 0.500 | 1.450, 1.188 | share intermediate |
| hu80855C/1126, /1128 | male | GRCh37 Isaac | haploid X, x 0.552 | 162 sites | one X, Y too few sites |
| hu832966/3674, hu6ABACE/3884 | male | GRCh38 Helix | haploid X, x 0.500, 0.517 | 286, 251 sites | one X, Y too few sites |
| huFD37C8/3255 | male | GRCh37 exome | share 0.088 | - | `not_wgs` |
| huE2F2FC/2771 | male | starling | - | - | load failed: polars-bio "invalid BGZF ISIZE" |
| huBEF5E2, hu24385B, hu1C1368, huC39675, huBE0518, hu416394, hu2BA77F | unlabelled | GATK, DRAGEN, Isaac, bcftools | one X each | 0.548–1.970 | XY |
| hu876E3D/4118 | unlabelled | GRCh37 DRAGEN | haploid X, x 0.471 | 0.471 | one X, Y not shown |
| huCDFB7C/3206 | unlabelled | GRCh37 GATK UG | 33 diploid X calls | 0 sites | too few sites |

What it changed:

- **No Y threshold says absent.** Self-reported women read 0.125 and 0.222, the local two-X genomes
  0.351–0.400, but one man's GATK genome read 0.433 and 0.484 (two alignments of the same reads). The
  first build of this layer had an absent band at 0.44 and below and called him `XO`. Y presence is
  now only ever `True` (ratio 0.52 or more, which every labelled or clustered two-X genome stays 0.12
  under), and `XX` and `XO`, which need a Y shown absent, are not emitted.
- **A sex-aware caller's few diploid X calls are no het share.** huBE0518's DRAGEN wrote 103,154
  haploid X ALT calls and 1,315 diploid ones, all hets; read as a share they said two X against a
  one-X depth. The share now counts only where most X ALT calls were written diploid. Its header's
  `##referenceSexKaryotype=XY` agrees with the depth-based answer, an independent check.
- **freeBayes males read an intermediate het share** (0.378, 0.396; those files have no AD). The share
  fell in the unknown band, so they answer unknown although their X depth (0.500) says one X: a share
  in the band is not overruled by the other line. Whether freeBayes needs its own band is open.
- Isaac males sit at X depth 0.548–0.560, above the local 0.42–0.50 and under the 0.60 edge.
- The X-degenerate Y het share, also measured, does not transfer: 0.40–1.00 in local one-X facets,
  0.06–0.12 in hu3BDC4B's GATK files. It is reported, not used.
- Two Helix uploads labelled exome read `callset_type = wgs` (RM4), with only 3,211 and 4,214 chr20
  ALT calls; their Y line had too few sites anyway.

**Built 2026-10-22: the evidence layer** (`just_vcf.karyotype`, exported `sex_chromosome_evidence`,
`SexChromosomeEvidence`, `KaryotypicSex`, `XCopies`; docs/SCHEMA.md "Sex-chromosome evidence"). Per
facet: GT ploidy counts per group, the X het share on diploid-written calls with the AD allele balance,
the X and X-degenerate Y depth ratios against chr20, and the facet's `ploidy_statements` carried
verbatim. The estimate gives `x_copies` and `y_present` separately, then `XY`, `XXY` or
`UNKNOWN_KARYOTYPE` with a reason (`no_build`, `unsupported_build`, `not_wgs`, `too_few_sites`,
`intermediate`, `lines_disagree`, `caller_asserted_only`, `no_depth`, `y_unusable`; `y_reason` may also
be `y_not_shown`). `XX` and `XO` are not emitted (above).

Thresholds: each unknown band is the central half of the measured gap, so every class edge keeps a
quarter of the gap as margin. X het share: one X at most 0.23, two at least 0.46 (gap 0.115–0.580).
X depth ratio: one X at most 0.60, two from 0.79 to 1.25 (gap 0.50–0.89; above 1.25 fits three copies
better than two, so it is unknown, not two). Y depth ratio: present at 0.52 or more, otherwise not
shown, never absent (PGP-HMS, above). Minimum sites: 1,000 diploid X ALT calls for the share, 1,000 X and 1,000 chr20 ALT
calls with DP for the X ratio, 300 X-degenerate Y ALT calls with DP for the Y ratio. Where the het share
and the X depth both decide, they must agree. The share counts only where most X ALT calls were written
diploid; otherwise (a sex-aware caller) depth alone decides.

Decisions, each with what was rejected:

- Computed on demand, one result per facet. Rejected: a struct column on `facets` (a calls scan at
  every load, a schema bump that orphans the cache, thresholds frozen into data) and a person-level
  verdict (no person in the schema; the critique's advice).
- A `Loaded` with several facets (a freeBayes+VarScan2 merge) gets one answer per facet, in `facets`
  order. Rejected: a `facet_id` argument, which makes the common case a lookup.
- "Not failed" is `filters` null or `["PASS"]`, and the same for FORMAT `FT`. Rejected: literal PASS
  only, which leaves a GATK file whose every FILTER is `.` with no rows.
- The depth denominator is chr20, because the thresholds were measured against it. Rejected: all
  autosomes, which would move every ratio off its calibration.
- Rows carrying a `split_signal` stay out of the het share (a `0` there can be another ALT, RM270).
- The allele balance is reported, not used: no calibration says where a mosaic begins.
- X copies stay known when the karyotype is not (`y_unusable`, `no_depth`): what just-prs needs for
  dosage (S1).
- T2T-CHM13 answers `unsupported_build`, not `no_build`: it has a label but no region set here.
- The header's ploidy lines are carried and the caller's haploid genotypes counted, never as evidence.
- Y presence is never `False`. Rejected: an absent band under a threshold, which named a man `XO`
  (above); and dropping the Y line altogether, which would lose `XY` and `XXY` where Y depth is high
  and the women measured sit far below.

**Left:** `XX` and `XO` need a Y shown absent, and XXX, XYY and the other classes a copy count:
both wait for a coverage facet (depth at called sites cannot show a region empty, and does not
separate one Y from two). An array's female Y probes are no-calls, which may be the absence signal
arrays have and genomes lack. Arrays (no `gt`
before RM3, no depth) and exomes answer unknown; an array's own lines (female Y probes no-calls,
AncestryDNA's doubled X letters) are not built. Liftovers of one callset are not yet weighed as one
piece of evidence (`facets_disagree` belongs to the person-level layer). GRCh37's last X-degenerate
interval may reach about 99 kb into NCBI's PAR2 (docs/PROVENANCE.md). The somalier question in the
critique is not answered.
