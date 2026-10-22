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
Severity: high · Status: in progress (2026-10-22: the FASTA tier shipped as `resolve_ref` on a gtars refget store, for VCF and array rows; indel anchor bases and reference-supported positions are left) · Owner: lead

Port just-prs `reference_allele.py` (reference-panel tier, FASTA tier, FASTA build check). New work:
ALT derivation from genotype against REF, indel anchor bases, chr-prefixed FASTA support, emitting
reference-supported positions.

**Shipped 2026-10-22: `resolve_ref` (`just_vcf.reference`).** The rules are in docs/SCHEMA.md
(Decision 2, "How RM3 fills these slots"). Choices made, each with the alternative rejected:

- **Same facet, filled, not a derived facet.** Every column it writes is one the schema defines as
  "not read yet", and the result is fixed by the facet and its assembly's sequence. A derived facet
  (`relation = ref_checked`) would leave every loaded facet `not_run` for ever beside a checked twin
  with a second id. The consequence: a resolved and an unresolved frame of one `facet_id` can both
  exist in a session, so a consumer must not join a resolved `facets` with an unresolved `calls`;
  `ref_seqcol_digest` (null before) and `ref` tell them apart. The lead may still prefer the derived
  facet; nothing else depends on the choice.
- **The reference is a gtars `RefgetStore` on disk** (`JUST_VCF_REFERENCE_DIR`, default
  `reference/` under the cache dir), over pysam faidx: it is a core dependency, exact, and gives the
  seqcol digest the facet records. A FASTA is added once; a fingerprint alias (path, size, mtime)
  keeps it from being digested again.
- **Assembly matching per contig, by stated length and md5, else by label.** The reference's label
  comes from its contig 1 length (the just-prs fingerprint idea), or a label the caller asserts when
  adding a FASTA, which its lengths must not contradict. Rejected: matching by label alone (hg19 and
  GRCh37 share it and differ on MT), and refusing a facet over any one unmatched contig (an hg19 VCF
  would lose every autosome over chrM).
- **Array ALTs sorted by allele string**, not in genotype order, so two exports of one chip agree.
- **A store holding several collections** picks the one matching the most rows; a tie is refused
  unless `collection=` names one.

**Dogfooded 2026-10-22** against Ensembl GRCh37 and GRCh38 primary assemblies in one store (built
from the gzipped FASTAs in 103 s and 25 s, 2.2 GB; each assembly named by its own contig 1 length),
on PGP-HMS files read from the Krypton share, aggregates only:

- Six plain VCFs (73–142 MB). Five resolved: 34.1M rows read, all but 2 matching their source REF, and
  both of those against an IUPAC code in the reference. One stated lengths for 194 GRCh38 contigs in
  UCSC names; its `chr*_random` and `chrUn_*` contigs (44,594 rows, 1%) stayed unread, because
  Ensembl names them by accession (`KI270721.1`) and `canonical_chrom` does not map the two. One VCF
  with neither `##contig` lines nor a build hint was refused. A sixth failed to load (not RM3's).
- Four build-37 23andMe exports (1.0–1.4M rows): 97.5–99.4% of rows read; the rest are I/D codes,
  no-calls and MT. 14.7–19.4% of the rows read carry no REF allele, the range expected for
  homozygous non-reference calls on a chip, not for a wrong build or strand. Four build-36 exports
  were refused: no GRCh36 reference in the store.
- `resolve_ref` on the 1.4M-row export: 11.4 s (load 1.3 s), of which 1.3 s reading bases.

Order with `rejoin_split`: re-join first, then `resolve_ref`. The other way round, the re-joined rows
carry no `ref` while the facet keeps `ref_seqcol_digest`, and `resolve_ref` refuses it as already read.

Left: indel anchor bases (no loaded format lacks them: 23andMe I/D codes need the chip manifest, and
Living DNA's run-together `CCA` cannot be split without the site's alleles); emitting
reference-supported positions; an opt-in test against a real assembly; mapping UCSC `_random`/`chrUn_` names to the accessions Ensembl uses (by the accession inside the name, checked by length). A VCF with neither
`##contig` lengths nor a build hint, and an array without a stated build, are refused: nothing can
show which assembly they are on.

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


## RM271 — VRS allele ids
Severity: high · Status: in progress (2026-10-22: `mint_vrs_ids` and `vrs_withheld` shipped for VCF and array rows, equal to vrs-python on 19,000 real alleles; left: the gtars RLE gap filed upstream, an opt-in test on a real assembly) · Owner: lead

Charter invariant 4: a call's identity is the VRS Allele id of its fully justified allele on a
refget-addressed sequence, and that id is the join key between facets on the same assembly. Today
`calls.vrs_allele_ids` is filled nowhere, though gtars (a core dependency) mints VRS ids and normalizes
alleles. It became possible with RM3: `resolve_ref` reads REF from a refget store and records the
collection on `facets.ref_seqcol_digest`, so each allele can be placed on a sequence with a refget
accession. An id that cannot be minted conformantly (an allele that cannot be justified, a contig with
no refget accession, an unresolved REF, an array I/D code) is withheld with a counted reason in
`vrs_errors`, never approximated. vrs-python is a test oracle only (CLAUDE.md §5).

**Shipped 2026-10-22: `mint_vrs_ids(loaded, reference)` and `vrs_withheld(calls)` (`just_vcf.vrs`).**
The rules are in docs/SCHEMA.md ("When an id is minted", "Why an id is missing"). Choices made, each
with the alternative rejected:

- **A step after `resolve_ref`, not inside it.** It refuses a facet whose REF was never read, and
  reads the collection `ref_seqcol_digest` names from the same store. Rejected: minting inside
  `resolve_ref`, which would make every REF read pay for ids, and leave no way to re-mint when the
  rule changes (`MINT_VERSION`) without reading REF again.
- **Every allele waits for `ref`, substitutions too.** Rejected: SCHEMA's earlier no-bases path for
  substitutions from a table of contig accessions per assembly. No table exists, and a contig
  `resolve_ref` could not match has no accession shown to be the facet's.
- **`[ref] + alts`, the REF slot minted** as VRS's reference allele (a `ReferenceLengthExpression`
  over REF's span), as SCHEMA and gnomAD have it. Rejected: one id per ALT only, which the brief
  named; it would retype a signed-off column.
- **The state choice and canonical JSON are ours; rolling and sha512t24u are gtars'.** gtars 0.10
  rolls an indel correctly but always returns a literal sequence, and its Python API digests only
  `LiteralSequenceExpression`; vrs-python gives a `ReferenceLengthExpression` for every deletion, every
  repeat insertion and the reference allele (VRS 2.x), so gtars' own ids differ from vrs-python's
  there. Rejected: withholding those alleles until gtars has it (most indels and every REF slot), and
  vrs-python in core (CLAUDE.md §5). Owed upstream: a VRS 2.x state in gtars' normalization and an RLE
  digest in its Python API; then `_justify` and the JSON shrink to calls. Not filed.
- **Rows walked as runs of one contig** (`rle` over the spilled input), reading only the columns
  minting needs; the ids go to their own files and sit beside the spilled rows by position. Rejected:
  a filtered scan per contig with a re-sort (resolve_ref's shape), one scan of the file per contig (170
  on the GRCh38 genome).
- **The oracle is `scripts/vrs_oracle.py`**: a vrs-python data proxy over the same gtars store, so
  both sides read the same sequence; vrs-python is a dev dependency.

**Dogfooded 2026-10-22** on three PGP-HMS VCFs, aggregates only (`scripts/vrs_conformance.py`), each
sample drawn with a fixed seed from the minted alleles and checked against vrs-python on the same
store:

- A GRCh38 genome, 4,651,432 rows: 9,256,478 of 9,345,929 alleles minted (363,409 deletions, 412,401
  insertions); 89,451 withheld as `ref_unchecked`, all on contigs `resolve_ref` left unread (UCSC
  `_random`/`chrUn_` names, RM3). 8,000 checked, 0 mismatch.
- A GRCh37 indel-only genome split, 956,064 rows, 1,039,271 ALTs: every allele minted; 9,000 checked
  (deletions, insertions, REF), 0 mismatch.
- A GRCh37 SNV-only genome split, 3,461,640 rows: all but 4 alleles minted (the 4 on rows whose source
  REF the reference contradicts); 2,000 checked, 0 mismatch.
- Time on the 4.65M-row genome: 59 s wall and 103–112 s CPU with the cache on a fast local disk; 112–208
  s wall with the cache under `data/interim/` on the NTFS disk, which `dd` measured at 2.5 MB/s during
  the runs (other sessions were writing to it; the workspace note says about 350 MB/s). Peak RSS of the
  minting process 1.5–2.0 GB; one contig's alleles and, for its indels, its sequence are in memory.

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
Severity: high · Status: in progress (the person chose the full label, desk #61, 2026-10-22) · Owner: lead

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
