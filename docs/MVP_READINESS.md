# MVP readiness and the extended shape

Written 2026-10-22 by the lead seat, to plan a long build run after an interview with the person. Part 1
checks the charter's minimal shape and its invariants against the code at f77a324: status lines come
from `docs/ROADMAP.md`'s dated checks and from greps of `src/`, not from memory. Part 2 lists candidates
for the extended shape, each with the question the interview has to settle before the run starts, so the
run never has to stop and ask. Nothing in Part 2 is on the roadmap until the interview accepts it; ids
are claimed then, not before.

**In one line (updated 2026-10-22, after the first builder waves):** the readers, REF lookup and the
VCF writer exist; VRS ids, normalization and the sex-chromosome finding are being built. just-vcf 0.1.0
is not on PyPI.

## Part 1 — the minimal shape against the code

Status: **shipped**, **partial**, or **not started**.

| # | Item | Charter | Status | Evidence | Gap |
|---|---|---|---|---|---|
| 1 | Recognise a file from its content | minimal shape | shipped | `sniff`, 39 declared formats, compression peeled by magic bytes | none for MVP |
| 2 | Load VCF and gVCF, one sample | minimal shape | shipped | `load`, `sample=`; gVCF blocks to `blocks`; spill to parquet | dialects `vcf_cgi`, `vcf_freebayes_varscan2`, `vcf_no_meta`, `bcf` do not load (RM2) |
| 3 | Load consumer arrays | minimal shape | partial | 22 of 39 formats load (RM2, 2026-10-22): the rsid tables, AncestryDNA, FTDNA, MyHeritage-layout, Color, deCODEme, Genos, Illumina Final Report, and the VCF dialects | CGI native files and PGP's key=value table (a builder is on them); Geno 2.0 and Forward/Design-only Final Reports need chip manifests |
| 4 | Refuse what it cannot hold | inv. 9 | shipped | `UnsupportedFormatError`, `MultiSampleError`, `NoSampleError` (F8) | none |
| 5 | Assembly identified by content (seqcol digest) | inv. 3 | partial | RM4: `sorted_name_length_pairs_digest` (the header's coordinate system); RM3: `ref_seqcol_digest` names the collection REF was read from | `seqcol_digest` stays null by design: no header carries per-sequence digests (SCHEMA.md) |
| 6 | Assembly label | inv. 3, 9 | partial | RM4 infers it from `##contig` lengths (basis `inferred`) | arrays that state no build: F9 (a builder infers it from the rows against the reference) |
| 7 | Allele identity (VRS ids) | inv. 4 | in progress | RM271 (a builder, on gtars after RM3) | — |
| 8 | REF from the reference | inv. 10 | shipped | `resolve_ref` (RM3) on a gtars refget store; arrays get REF and ALT from the reference; 34.1M PGP VCF rows checked, 2 mismatches (both IUPAC) | arrays without a stated build (F9); UCSC `chrUn`/`_random` names (1% of GRCh38 rows) |
| 9 | Callset type and what absence means | inv. 11 | shipped | RM4: `wgs`/`wes`/`panel`/`imputed` from the rows, `array`/`imputed` from the header (RM2), `absence_means` follows | SNV-only genome splits read as `wgs` (open question in RM4) |
| 10 | Graded genotype evidence | inv. 7 | partial | `gl`, `gp`, `gq`, `dp`, `ad`, `evidence_basis` carried as written | no calibrated vector across sources; arrays carry none |
| 11 | Provenance per facet | inv. 6 | partial | source path, sha, format, vendor, `caller`, `caller_version`, `phasing_method`, `ploidy_statements` from the header (RM2); `derived_from` for derived facets | aligner and post-processing steps are not recorded |
| 12 | Pre-split multi-allelic input | inv. 8 | shipped | `split_signal`, `rejoin_split` (RM270) | not run on files split by bcftools itself, Hail or GATK |
| 13 | Sex chromosomes and caller ploidy | inv. 9 | in progress | RM269, full label chosen (desk #61); paused at step 1 of 5 (f77a324) | the whole build |
| 14 | Normalize: chrom names, left-align, split, dedup, sort, strand | minimal shape | in progress | RM6 second slice (a builder): split, left-align, dedup as a derived facet | strand for arrays; liftover |
| 15 | Write VCF back out | minimal shape (export) | shipped | `write_vcf` (RM6, a1f3317): VCF 4.3, BGZF, deterministic bytes; 6.4M-call round trip with 0 key mismatches; a tabix index beside BGZF output | gVCF blocks, SVLEN, PSL phase sets |
| 16 | Liftover with a `liftover_to` link | inv. 5 | not started | RM6 (port from just-prs) | `bcftools +liftover` per LANDSCAPE |
| 17 | Deterministic output | inv. 12 | shipped | canonical order per contig; `write_vcf` gives the same bytes for the same frames (tested) | — |
| 18 | CLI | (CLAUDE.md §4) | not started | no `[project.scripts]` | Typer |
| 19 | Cache eviction | — | shipped | RM268: a manifest per entry, the stale-version sweep, `clear_cache(source, older_than, dry_run)` | no size cap (rejected: it can evict under a live frame); no CLI yet |
| 20 | Population inference | — | partial | `just_vcf.ancestry` (RM5): reference-PCA projection, neighbour vote, proportions, consensus, Privé's 21-group proportions with standard errors; gated on the facet's assembly and `absence_means`; run on PGP-HMS genomes and arrays against Privé | GRCh37 only for Privé; no reference-PCA panel on this machine; no liftover (RM6) |
| 21 | Charter signed | RM1 | not started | `CHARTER_DRAFT.md`, the person's call | — |
| 22 | First PyPI release | — | not started | 0.1.0 in `pyproject.toml`, absent from PyPI | a release is the person's |

What an MVP needs, at the least: rows 3, 5–9 and 14–15, so a consumer can hand in an array or a VCF
and get back a REF-checked, build-known, VRS-keyed frame and a VCF. just-prs and just-dna-lite each
carry their own copies of rows 6, 8, 14 and 16 today (docs/SIBLING_SURVEY.md); those are what they
would delete first.

Rough sizes, for one lead working through them: row 8 (REF lookup with a FASTA, arrays included) 2
days; rows 5 and 6 together 1 day; row 7 1 day once 5 and 8 exist; row 9 1 day; rows 14–15 3 days;
row 16 2 days; row 3's remaining vendor loaders half a day each; RM269 (row 13) 3–4 days.

## Part 2 — the extended shape: candidates for the interview

Sources: the charter's extended-shape paragraph (facets across assemblies, callers and phasings;
coverage as a facet; agreement measured; a graph), LANDSCAPE's "genuinely absent" list, and the
dreams report's convergence (docs/DREAMS.md, untracked, which says of itself "inspiration, not a
plan"). Labels are local (E1, E2, …), not roadmap ids. Most items need Part 1's rows 5, 7 and 8 first,
and several need a person (a subject) in the schema, which does not exist: frames are per facet today.

Calibration data for several items is on the Krypton array: 775 PGP-HMS participants and 326 VCFs,
consented-public, some with several facets per person.

| # | Candidate | Needs first | Serves | Size | Question for the interview | Options | Recommended |
|---|---|---|---|---|---|---|---|
| E1 | A person (subject) and a store where facets accumulate | rows 5, 7 | everything below; just-prs, just-dna-lite | 3 days | What is a person in the schema? | a) a subject id the caller supplies, nothing inferred; b) inferred from relatedness between facets; c) none yet, facets only | a: inference is E6's job, and a supplied id cannot be wrong in a way the library hides |
| E2 | Concordance between facets on graded evidence, by region and method | E1, row 7 | inv. 8; method atlas (E8) | 3 days | Concordance on what? | a) hard genotypes only, as every tool does; b) graded (GL/GP) where both facets carry it, hard otherwise, said per row; c) graded only | b: the charter wants graded, the data often has none |
| E3 | Coverage as a facet (mosdepth summaries, polars-bio depth over a BAM) | E1 | inv. 11 (absence with a cause); RM269's Y evidence | 2 days | Which coverage inputs? | a) mosdepth/samtools summaries only; b) plus polars-bio `depth()` over a BAM the person has; c) neither for now | b: the person's 2026-09-28 note names polars-bio depth as in scope |
| E4 | The facet graph: `liftover_to`, `derived_from` and other relations as edges | row 16, E1 | inv. 5 | 1 day | Edges as a frame or as columns? | a) a small `links` frame; b) columns on `facets` only | a: one facet can have several parents (a merge) |
| E5 | Phasing as VRS CisPhasedBlock | row 7 | charter (GA4GH phasing) | 2 days | Mint CisPhasedBlock objects, or keep phase sets as columns and export blocks on demand? | a) columns, export on demand; b) a phased-block frame | a: lazy and additive |
| E6 | Families: GA4GH Pedigree, relatedness | E1 | charter (pedigrees) | 3 days | Relatedness from our own sites or from somalier as a subprocess? | a) somalier; b) our own; c) pedigree as supplied data only | c first, then a |
| E7 | A consensus posterior per site, as a derived facet | E1, E2, row 7 | dreams theme B; just-prs scoring | 4 days | What does the posterior use as error priors? | a) flat priors per caller; b) priors from the method atlas (E8); c) not yet | c until E8 exists, then b |
| E8 | A method-error atlas from multi-facet PGP-HMS people | E2 | inv. 8 (reliability by method and region) | 3 days | Which regions define "region"? | a) GIAB stratifications; b) our own segdup/repeat labels; c) per chromosome only | a, if their licence lets us carry the BEDs |
| E9 | Chip manifests: arrays resolve I/D codes and strand, and know which probes they have | row 8 | row 3; dreams theme F | 3 days | Where do manifests come from? | a) vendor manifests where public; b) learned from paired PGP-HMS array+WGS people; c) both | c |
| E10 | T2T-CHM13 facets | rows 5, 8 | charter (assemblies) | 1 day | Support CHM13 now or after GRCh37/38 are done? | a) now; b) after | b |
| E11 | Structural-variant facets | row 7 | charter (callers) | 3 days | VRS has no SV model in 2.x: withhold ids for SVs, or use a just-vcf shape? | a) withhold; b) own shape, documented | a: inv. 2 and 4 |
| E12 | Privacy as structure: agents get counts and answers, never rows | none | dreams theme E; CLAUDE.md's one threat | 2 days | A query API that returns only aggregates? | a) yes, as the default public face; b) a helper beside `load`; c) not in the library | b |
| E13 | A durable per-person store on disk (parquet or VCF Zarr) | E1 | dreams theme A | 2 days | Which format? | a) parquet per facet in the cache; b) VCF Zarr (no gVCF blocks in its spec); c) both | a |
| E14 | The format registry and broken-file corpus as a public artifact | none | dreams theme C | 1 day | Publish the corpus of shapes (no data) as a spec? | a) yes, in this repo; b) a separate repo; c) no | a |

Order the run would take if every recommendation were accepted: Part 1 rows 5–9, then 13 (RM269),
then 14–16, then E1, E4, E13, E2, E3, E8, E7, with E5, E6, E9–E12 and E14 as the interview ranks
them. Each step lands as its own commit with its test and changelog line; each Part 2 item gets its
`RMn` when it is accepted.
