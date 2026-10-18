# Roadmap

Active items only, newest decisions first. Shipped items move to
[ROADMAP_HISTORY.md](ROADMAP_HISTORY.md). Claim a new id with
`INDEX=docs/ROADMAP.md ITEM_PREFIX=RM ./scripts/item-next.py --note "…"`.

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
Severity: high · Status: in progress (detection shipped 2026-09-29; 23andMe-layout and plain VCF loaders shipped 2026-09-30; checked 2026-10-18: 10 of the 39 declared formats load, the 9 rsid-table layouts and `vcf`; the VCF dialects `vcf_cgi`, `vcf_freebayes_varscan2`, `vcf_no_meta` and `bcf`, CGI, AncestryDNA, FTDNA and the vendor CSVs are next) · Owner: lead

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

**Also wanted: a `source` detector** naming the caller or vendor from VCF headers, because `Family`
cannot say "array or imputed content inside a VCF" and that decides whether absent means unknown.
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
Severity: medium · Status: open (checked 2026-10-18: array and imputed tables set `callset_type`, and Illumina gVCF blocks go to `blocks`; a VCF facet's `callset_type` is still null, and WGS/WES/panel is not inferred) · Owner: lead

Port just-dna-lite `detect_callset_scope` (WGS vs targeted), gVCF detection, and just-prs chip
generation detection. New work: WES vs WGS vs panel, and expanding gVCF `END` blocks.

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


## RM268 — Cache eviction
Severity: medium · Status: open · Owner: lead

Nothing ever removes anything from `JUST_VCF_CACHE_DIR`. It holds two kinds of entry:

- `vcf/<key>/`: rewritten copies of inputs polars-bio cannot read as they are
  (`loaders/vcf.py` `readable_path`), plus a `structural/` copy for a file polars-bio refused. The key
  hashes the source digest, the archive member and `REWRITE_VERSION`.
- `frames/<facet_id>-<SCHEMA_VERSION>/`: the calls and blocks parquet spill from `load`.

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
Severity: high · Status: open · Owner: lead

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
   does not. Which headers say so is not checked yet: DRAGEN's and GATK's recorded command lines are
   the candidates, and whether any DeepVariant header records haploid contigs is open. Never guessed
   from the caller's name alone.
2. **A person-level karyotype finding with its own provenance.** Evidence without depth: the non-PAR X
   heterozygous fraction and the presence of PASS Y calls, which separate XX, XY, XXY and X0 on a WGS
   facet. XYY needs depth, so it waits for a coverage facet. Weak or contradictory evidence gives
   unknown, never the nearest class. The vocabulary is GA4GH Phenopackets' `KaryotypicSex`, used
   verbatim (charter invariant 2). PAR bounds are per assembly and are domain constants.

Constraints: list length stays the number of haplotypes the source wrote (Decision 2). The finding sits
beside the facets and is never a rewrite of `calls`: no GT is collapsed to haploid by the loader. An
array facet without Y probes, or a targeted panel, says unknown.

Legality: minor. Filling `caller`/`caller_version`, a new nullable facet field and a new person-level
finding are all additive. Where the finding lives (a new frame or columns on an existing one) is a
schema decision for the lead. just-prs builds a stopgap now (karyotype from its normalized parquet,
caller and version as parquet metadata) and moves to this when it exists.

## RM270 — Pre-split multi-allelic input
Severity: high · Status: in progress (B, the flags, shipped 2026-10-19 as schema 0.0.4; C, withholding, waits on desk #45; D goes with RM6) · Owner: lead

Motivating report: S2 (just-dna-lite, 2026-10-18). Many VCFs arrive already split one row per ALT, in
the `bcftools norm -m-` convention, and their headers do not say so. In a split row the `0` of a
genotype can stand for "another ALT", not REF: a `1/2` at `A>G,T` becomes `1/0` at `A>G` and `0/1` at
`A>T`, and a `1/1` becomes `1/1` plus a `0/0` row for the uncarried ALT.

Reproduced 2026-10-18 against the code, on a synthetic VCF: the loader keeps one `calls` row per record
(SCHEMA Decision 1, A), so the split rows land as written, and `gt_alleles` maps each `0` onto
`source_ref`. The split `1/2` loads as `["G","A"]` and `["A","T"]`, two calls with a reference copy the
person does not have; the uncarried ALT's `0/0` row loads as `["C","C"]`, a hom-ref call. `warnings` is
empty. The unsplit record of the same site loads as `["G","T"]`.

**Prior art** (read from source and docs, 2026-10-18). bcftools `norm -m-` writes the other ALT as `0`
by default, keeping GT slots, and `.` with `--multi-overlaps .` (1.17+); its `--old-rec-tag` INFO
(1.12+) is opt-in. Hail's `split_multi_hts` downcodes to `0` and sorts unphased calls, so its split
`1/2` is `0/1` + `0/1`, and its `a_index`/`was_split` are row fields an exported VCF lacks. GATK
`LeftAlignAndTrimVariants --split-multi-allelics` writes `0`, except that every sample becomes `./.` when
any sample is `1/2`. vt `decompose`, Glow and vcflib `vcfbreakmulti` write `.`, and vt and Glow add INFO
`OLD_MULTIALLELIC`. The VCF spec says nothing about split records. `bcftools norm -m+` re-joins slot by
slot (`1/0` + `0/1` gives `1/2`). Rows written with `.` already load as half-calls and are safe.

**Detection**, three signals per `pos` within a contig, reported separately and never folded into one
flag: a **marker** (INFO `OLD_MULTIALLELIC`, or a value of the `--old-rec-tag` shape
`CHROM|POS|REF|ALT|IDX`); the **slot pattern** (two or more non-reference rows, each GT slot non-zero in
at most one); and an **all-zero row of any ploidy beside a non-reference row**. On the facet, a marker
is Kleene: true when seen, null otherwise, never false, because a split file can show no signal at all
(Hail's sorted output, split indels that left-alignment moved to another POS). Measured 2026-10-18 on
four unsplit private WGS VCFs (DeepVariant ×2, DRAGEN, GATK HaplotypeCaller), 22.2 million calls: no
signal fired. Only the DRAGEN file had positions with two records (10), none matching. The rate on a
real split file is not measured (no bcftools here); a synthetic split VCF is caught at both its sites.
The two unsplit DeepVariant gVCFs gave no signal either: no two call rows share a position once their
reference blocks go to `blocks`, which never take part.

**What the frame says about a `0` it cannot place.** No option changes `gt`, the row grain, or
`source_gt` where it is kept (null on VCF 4.4+, P4). Not decided:

| | What the frame says | Legality | Cost |
|---|---|---|---|
| A. Warn only | a warning per signal, with a count | patch | every flagged row stays wrong for a reader that skips warnings |
| B. A + flags | `facets.split_marker` (Kleene) and a nullable row column `split_signal` naming the signal | minor | values stay wrong, but every flagged row can be found |
| C. B + slot-aware withhold | in a flagged group, a `0` becomes null in `gt_alleles` only in a slot where another row of the group is non-zero; a slot zero in every row keeps its REF | values read today change (pre-1.0, changelog says so) | a false positive loses a true REF copy, never states a false one; amends Decision 2: `gt` and `gt_alleles` agree except where withheld, and a null means uncalled or withheld |
| D. Re-join as a derived facet | each flagged group merged into one site slot by slot, as `bcftools -m+` does (`relation = normalized`, charter 5) | minor | the most work; GT only, since bcftools subsets PL/AD to REF + one ALT; phasing across slots may be lost |

`gl`, `gp` and `ad` on a flagged row stay as written: they compare REF with one ALT and inherit the
ambiguity, and the row flag is what tells a reader so. The lead recommends B and C now, D inside RM6,
put to the person on the desk 2026-10-19 (desk #45; on silence, B is built and C waits for a yes).

**RM6's own split output** must not write this convention: it writes `.` for the other ALT, as vt,
Glow, vcflib and `bcftools --multi-overlaps .` do.

Legality: detection, a warning and a new nullable field are additive (minor). Nulling `gt_alleles` for
detected rows changes the meaning of values already emitted.
