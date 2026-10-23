# Live samples for the extended shape

Who's who among the genomes on this machine, for live tests of the charter's extended shape (one person,
many facets). Surveyed 2026-10-22 by a read-only Opus agent for the lead seat. Private genomes appear only
by their folder alias, which CLAUDE.md §11 already names; no sample name, genotype or genotype-derived fact is
recorded here, and a library number is a sequencing-run label, not a person. Tests that use the private sets
run locally and skip without them; anything shown or shipped uses PGP-HMS.

The opt-in tests read their directories from the shell only: `tests/conftest.py` defaults `JUST_VCF_SAMPLES_DIR`, `JUST_VCF_ANCESTRY_DIR`, `JUST_VCF_MANIFESTS_DIR`, `JUST_VCF_TEST_REFERENCE_DIR` and `JUST_VCF_TEST_CHAIN_DIR` to empty, whatever the repo `.env` says. Name the ones a run needs on its command line:

```bash
uv run pytest -q                                                  # plain: every opt-in test skips
JUST_VCF_SAMPLES_DIR=/data/genomes uv run pytest -q tests/test_person.py
```


Sources: the Krypton copy `/mnt/krypton/Genomes/genomes/<alias>/` (the fuller one), the local
`/data/genomes/<alias>/`, the a3-bench inputs copy of newton, `../just-prs/data/input/vcf/`. Every fact comes from
headers (VCF `##` lines, BAM `@PG`/`@SQ`), magic bytes, sizes and record counts. File names are masked. Sample
names, genotype rows and anything genotype-derived were not read into this report. Assembly comes from the
declared chr1 length (248956422 = GRCh38, 249250621 = GRCh37). Anything marked *inferred* is my reading, not a
header fact.

Layout facts first:
- `newton_n1` and `aadr` are **empty on Krypton** (0 files). `/data/genomes/_backups` is empty too.
- **`/data/genomes` no longer holds `pgp-hms` or `angellos`.** Both exist only on Krypton. CLAUDE.md §11 still
  says `/data/genomes/pgp-hms/`, so that line is now stale.
- Index numbers (`#k`) follow the sorted Krypton listing. Index files (`.tbi`, `.bai`), `.md5` files and HTML/JSON
  reports are left out of the tables. A `.bam.<x>` sidecar of about 0 MB with 10 columns is a flagstat file
  (*inferred*: `1000genomes_hg00157` names that extension outright).

---

## 1. Private aliases

### felis (40 files, 295 GB on Krypton; local copy holds the library-12 DeepVariant pair only)

The folder holds **three sequencing libraries, numbered 12, 13 and 14** in the file names, from one sequencing run
(the same 2025-09-19 run prefix). Library 13 is **a different person** from library 12 (BAM metadata), so the two must never be paired as one
person's facets. Who library 14 is, is unknown.

| # | kind | caller / tool + version | mode | assembly | derived from | size |
|---|---|---|---|---|---|---|
| 17, 18 | FASTQ R1/R2, lib 12 | - | reads | - | - | 31 + 32 GB |
| 38, 39 | FASTQ R1/R2, lib 13 | - | reads | - | - | 38 + 40 GB |
| 7, 8 | FASTQ R1/R2, lib 14 | - | reads | - | - | 41 + 43 GB |
| 1 | BAM, lib 12 (`aligned/`) | bwa-mem2 2.2.1, sambamba 1.0 | alignment, sorted + markdup | GRCh38 Ensembl primary (194 contigs, no alts) | lib-12 FASTQ | 39 GB |
| 4 | BAM, lib 13 | bwa-mem2 2.2.1, sambamba 1.0 | alignment | GRCh38 Ensembl primary | lib-13 FASTQ | 48 GB |
| 23 | gVCF, lib 12 (`variants/deepvariant/`) | DeepVariant 1.9.0 | gVCF, bgzf, **extension `.g.vcf.bgz`** | GRCh38 Ensembl (194) | BAM #1 | 599 MB, 50.9M records |
| 26 | VCF, lib 12 | DeepVariant 1.9.0 | VCF, bgzf, `.vcf.bgz` | GRCh38 Ensembl (194) | BAM #1 | 119 MB, 6.58M records, 4.73M PASS |
| 29 | zip of #23, #26, indexes, report | - | archive | - | a copy of #23 and #26 | 720 MB |
| 30 | VCF (`vcfs/`) | DeepVariant 1.9.0 | VCF | GRCh38 **UCSC-style with alts** (455 contigs: 261 alt, 169 random/Un) | **a second alignment** to a different hg38; that BAM is not on disk | 129 MB, 7.13M records, 4.60M PASS |
| 33 | VCF `dv_norm` | DV 1.9.0 + bcftools 1.23.1 `norm -m -both --old-rec-tag ORIG` | normalized | GRCh38 (455) | #30 (*inferred* from the file name `dv.vcf` in the norm command) | 126 MB, 7.21M records |
| 35 | VCF `variants_norm` | the same, an earlier norm run | normalized | GRCh38 (455) | #30, the same record count as #33 | 124 MB, 7.21M |
| 36, 37 | `variants_norm_clinvar`, `variants_norm_snpeff` | bcftools annotate (ClinVar), SnpEff 5.4a (GRCh38.mane.1.2) | annotation only | GRCh38 | #35 | 126 / 203 MB, 7.21M each |
| 34 | VCF `gatk_norm` | GATK 4.6.2.0 HaplotypeCaller `--dragen-mode true`, per-chromosome, bcftools concat + norm | VCF | GRCh38 (455) | the second alignment | 162 MB, 4.97M records, **PASS = 0 (every FILTER is `.`)** |
| 31 | VCF `clair3_norm` | Clair3 1.0.10 `--platform=ilmn`, bcftools norm | VCF | GRCh38, **header has 3366 contigs** (alts + decoys + HLA, the DRAGEN-style hg38) | the second alignment | 75 MB, 4.88M records |
| 32 | VCF, 3-caller consensus | bcftools view `-f PASS,.` then `isec -n=3 -w1` over the Clair3, GATK and DV normalized PASS files | consensus | GRCh38 (3366-contig header) | #31 + #34 + #33 | 70 MB, 4.54M records |
| 12-15 | ClinVar CSV/HTML reports | - | annotation | - | a call set | under 1 MB |

The `vcfs/` series does not name which library it came from. The `--bam_fn` / `--input` basenames carry no library
number.

**Independent vs derived.** Library 12 has reads, a BAM and a DeepVariant 1.9.0 gVCF+VCF pair, all one chain.
The `vcfs/` series is a separate chain: another alignment to a different hg38 build (UCSC naming with alts, not
Ensembl primary), then three callers on the same reads, then normalization and a consensus. Library 13 has reads
and a BAM but no calls. Library 14 has reads only.

**What it can test:**
- Caller concordance on identical reads: DV vs GATK-DRAGEN-mode vs Clair3, with the person's own `isec` consensus
  as a comparison point.
- The same caller version on two references: DV 1.9.0 on Ensembl primary (#26) vs on UCSC hg38 with alts (#30).
  This is exactly the "hg38 with alts vs no-alt" digest case of charter invariant 3.
- gVCF vs VCF (#23 vs #26).
- Reads to a coverage facet (BAMs #1 and #4).
- Contig-header mismatches inside one pipeline: 455 vs 3366 declared contigs.
- FILTER handling: GATK with an all-`.` FILTER column vs DV's PASS/RefCall.
- Pre-normalized multi-allelic splitting (`--old-rec-tag ORIG`).
- Annotation-only derivatives that must not count as new facets.
- Library 13 (and 14, once called) is another person sequenced in the same run: a two-person case once a
  subject exists in the schema.

### newton (6 files on Krypton; local copy holds the DeepVariant pair)

| # | kind | caller / tool + version | mode | assembly | derived from | size |
|---|---|---|---|---|---|---|
| 1 | BAM (`aligned/`) | bwa-mem2 2.2.1, samtools 1.11, markdup | alignment | GRCh38 Ensembl primary (194) | FASTQ (not on disk) | 85 GB |
| 4 | gVCF (`variants/deepvariant/`) | DeepVariant 1.1.0 | gVCF, bgzf | GRCh38 Ensembl (194) | BAM #1 (the same library tag in the name) | 458 MB |
| 5 | VCF | DeepVariant 1.1.0 | VCF, bgzf | GRCh38 Ensembl (194) | BAM #1 | 116 MB, 6,380,925 records |
| a3-bench copy | VCF, named `GRCh38.ens106` | DeepVariant 1.3.0 | VCF, bgzf | GRCh38 Ensembl (194) | unknown BAM | 125 MB, 6,138,868 records |
| just-prs `newton_winter.vcf` (and a byte-identical `.vcf.gz`) | VCF | DeepVariant 1.8.0 | VCF, **BGZF named `.vcf`** | GRCh38 Ensembl (194) | unknown BAM | 128 MB, 6,139,024 records |

**newton_n1 is empty on Krypton.** I found no liftover of newton anywhere I looked: Krypton, `/data/genomes`,
`/data/a3-bench`, just-prs, or a `find /data -maxdepth 6` for the alias. The person said to expect one. It may live
on another machine or under another alias. `anton`'s GRCh38 input is not newton's: it is also DeepVariant 1.1.0,
but its record count differs (6,079,744 vs 6,380,925).

**What it can test:** **DeepVariant version drift on one person: 1.1.0, 1.3.0 and 1.8.0.** The 1.3.0 and 1.8.0
counts differ by only 156 records, so those two are probably re-calls of the same BAM (*inferred*). 1.1.0 differs
by about 240k. Also gVCF vs VCF (1.1.0), reads to a coverage facet (the 85 GB BAM), and BGZF named `.vcf`.

### anton (14 files on Krypton, 2.6 GB; local copy holds the 7 VCFs)

| # | kind | caller / tool + version | mode | assembly | derived from | size |
|---|---|---|---|---|---|---|
| 12 | VCF `hg38` (`variants/`) | DeepVariant 1.1.0 | VCF, plain text | GRCh38 Ensembl (194) | - | 483 MB, 6,079,744 records |
| 1 | VCF, crossmap variant A | CrossMap 0.6.3, chain `GRCh38_to_GRCh37`, target GRCh37.75 | lifted | GRCh37 (297 contigs) | #12 | 445 MB, 5,608,690 |
| 3 | VCF, crossmap variant B | CrossMap 0.6.3, the same chain | lifted | GRCh37 | #12 | 444 MB, 5,594,866 |
| 2, 4 | CrossMap `.unmap` | - | unmapped records | GRCh38 header | #12 | 44-46 MB |
| 9 | VCF `hg37` | Picard-style LiftoverVcf (*inferred*: `.idx`, `rejected`, `picard.txt` note), reference header GRCh37.75 | lifted | GRCh37 (297) | #12 | 449 MB, 5,593,760 |
| 5 | VCF `hg37 recover` | the same with swapped REF/ALT recovery (*inferred* from `recover` in the name) | lifted | GRCh37 | #12 | 450 MB, 5,607,006 |
| 7, 8 | rejected VCFs | - | rejects | GRCh38 header | #12 | 42-44 MB, 472,738 / 485,984 |
| 6, 10 | `.idx` | - | index | - | - | 17 MB |
| 11 | `.7z` archive of an hg37 VCF | - | archive | - | - | 282 MB |
| 13, 14 | `crossmap.txt`, `picard.txt` notes | - | text | - | - | small |

Lifted plus rejected adds up exactly to the source: 5,593,760 + 485,984 and 5,607,006 + 472,738 both equal 6,079,744.

**What it can test:** **liftover vs liftover on one source**: CrossMap vs Picard, each in two settings, with
rejects and unmapped kept. That is the `liftover_to` relation and its rejects. There is no direct GRCh37 call, so
"lifted vs direct" cannot be tested here.

### olga (9 files on Krypton; local copy holds 5 VCFs)

| # | kind | caller / tool + version | mode | assembly | derived from | size |
|---|---|---|---|---|---|---|
| 9 | VCF | GATK 4.1.8.1 HaplotypeCaller (gVCF, FPGA build "ACC1"), then GATK 3.8-1 MGI GenotypeGVCFs, VQSR and SelectVariants | VCF, bgzf | GRCh38 UCSC no-alt (195: 169 random/Un, 1 EBV) | - | 221 MB, 4,737,368 |
| 6 | VCF `liftover` | Picard-style liftover (*inferred*: `.log`, `.idx`, `reject`), target GRCh37.75 | lifted | GRCh37 (297) | #9 | 211 MB, 4,520,991 |
| 2 | VCF `liftover swap` | the same with REF/ALT swap recovery (*inferred*) | lifted | GRCh37 | #9 | 211 MB, 4,532,668 |
| 5, 1 | reject VCFs | - | rejects | GRCh38 | #9 | 10 MB each, 216,377 / 204,700 |
| 3, 4, 7 | idx and logs | - | - | - | - | small |
| 8 | zip of #1-#7 | - | archive | - | a copy | 445 MB |

Lifted plus rejected equals 4,737,368 in both runs. **What it can test:** liftover with and without swap
recovery, rejects as data, the GRCh38 `chr` naming going to GRCh37 Ensembl naming, and a GATK3+4 VQSR pipeline's
FILTER vocabulary.

### alina (2 files; local holds the VCF)

| # | kind | caller | mode | assembly | derived | size |
|---|---|---|---|---|---|---|
| 2 | VCF | the same MGI pipeline as olga (GATK 4.1.8.1 HC, GATK 3.8 MGI GenotypeGVCFs/VQSR) | VCF, **plain gzip, not bgzf** | GRCh38 UCSC no-alt (195) | - | 219 MB, 4,768,231 |
| 1 | log | - | - | - | - | small |

One facet. Useful as a non-BGZF `.vcf.gz` case, and as a second person through the same lab pipeline as olga.

### sergey (2 files; both local)

| # | kind | caller | mode | assembly | derived | size |
|---|---|---|---|---|---|---|
| 1 | VCF `hard-filtered` | DRAGEN 4.2.4 (SW 05.121.676.4.2.4), hash table 4.2.1, ML model 7; sources `DRAGEN_SNV`, `DRAGEN_TARGETED` | VCF, bgzf | GRCh38 DRAGEN-style hg38 (3366 contigs: alts, decoys, HLA) | - | 462 MB, 4,833,398 |
| 2 | VCF `bqsr.hc` | GATK 3.8-1 MGI-6.1 HaplotypeCaller | VCF, bgzf | GRCh38 GCA_000001405.15 analysis set with decoys (2580 contigs, no alts) | - | 225 MB, 4,835,656 |

**What it can test:** **two vendors' callers and two hg38 contig sets for one person.** DRAGEN vs GATK3
concordance. Two different seqcol digests that should both resolve to GRCh38 primary sequences. DRAGEN-specific
FILTER and INFO vocabulary (29 FILTERs, EVENTTYPE). There are no reads. Whether both came from one sequencing run
or two is not stated (*unknown*).

### ksuhaster (50 files on Krypton, 678 GB; local holds the two DeepVariant pairs)

| # | kind | caller | mode | assembly | derived | size |
|---|---|---|---|---|---|---|
| 41, 42 | FASTQ pair (plain gzip) | - | reads | - | - | 56 + 58 GB |
| 43-50 | 4 FASTQ pairs (bgzf) | - | reads | - | - | 19-28 GB each |
| 1, 6, 11, 16, 21 | 5 BAMs (`aligned/`), read group with ID only | bwa-mem2 2.2.1 | alignment | GRCh38 Ensembl (194) | the FASTQ pairs | 71, 27, 34, 34, 35 GB |
| 5, 10, 15, 20, 25 | 5 BAMs with full read groups (ID/LB/PL/PU/SM) | bwa-mem2 2.2.1 + samtools 1.15.1 reheader | alignment | GRCh38 Ensembl | each from its sibling BAM | 73, 29, 36, 36, 37 GB |
| 30, 32 | gVCF + VCF (`deepvariant/`) | DeepVariant 1.4.0 | gVCF / VCF | GRCh38 Ensembl | one of the BAMs | 712 MB / 109 MB, 6,218,137 or 6,606,810 records (which count goes with which pair was not recorded) |
| 35, 37 | gVCF + VCF | DeepVariant 1.4.0 | gVCF / VCF | GRCh38 Ensembl | another BAM | 977 MB / 114 MB |
| 40 | zip of #30, #32 | - | archive | - | a copy | 824 MB |
| 27 | text manifest mentioning Illumina | - | - | - | - | small |

Five libraries, one with deeper sequencing and four shallower. Whether they are one person or several, and which
BAM fed which DeepVariant run, is **unknown** from metadata. **Ask the person.** If they are one person, this is
the best "same caller, different sequencing runs" set. If they are several, it is a small cohort on one pipeline.
It also tests a gVCF vs VCF pair twice over, and reads to coverage.

### advocate (10 files, 101 GB; local holds the DeepVariant pair)

FASTQ R1/R2 (33 + 32 GB) → BAM (bwa-mem2 2.2.1, 42 GB, GRCh38 Ensembl) → **DeepVariant 1.3.0** gVCF (903 MB) +
VCF (113 MB, 6,436,476 records). One linear chain. It tests gVCF vs VCF and reads to coverage, and sits as a
DV-1.3.0 data point next to newton's 1.3.0.

### john (14 files, 634 GB; local holds the DeepVariant pair)

Six FASTQs: 4 per-lane files (64-69 GB each) plus 2 merged files whose sizes equal the sum of two lanes each
(129 / 138 GB; merge *inferred* from the sizes). BAM (bwa-mem2 2.2.1, 147 GB, GRCh38 Ensembl). **DeepVariant
1.3.0** gVCF (283 MB) + VCF (117 MB, 6,079,801 records). One chain. It tests deep-coverage reads to a coverage
facet, and gVCF vs VCF.

### angellos (2 files, Krypton only)

| # | kind | vendor | build | size |
|---|---|---|---|---|
| 1 | array export in 23andMe format (`genotec`) | header mentions 23andMe; Genotek service | build 37 | 24 MB, 959,840 rows |
| 2 | a `.md` note | - | - | small |

Only an array, and the only consumer array among the private aliases. It does not pair with any WGS here.

### gennet (11 files; local holds 7)

PLINK 1.90 VCF exports, **60 samples** (HapMap CEU per CLAUDE.md §11). The source declares PLINK-style contig
lengths (chr1 = 247,195,691, the maximum position rather than the true length). From that number I *infer*
NCBI36. Picard-style liftover to GRCh37 (297 contigs), with and without swap recovery, each with rejects, `.idx`
and log. A multi-sample panel, not one person: it tests a liftover of a multi-sample file and PLINK's
non-assembly contig lengths.

### nanopore (4 files, Krypton only)

ONT FASTQ (75 GB, named by a 13-digit vendor id, MinKNOW read-header keys, no basecaller model in the header),
two vendor summary files (Dante extension), and a 1.2 MB tar.gz. **Reads only, no calls, no alignment.** Which
person it belongs to is unknown.

### oxford (3 files), 1000genomes_hg00157 (9 files), hg00157 (2 files)

- **oxford:** header-only stub BAMs (about 3 KB; minimap2 2.22, GRCh38 Ensembl) plus flagstat. Probably an
  aborted ONT alignment (*inferred*).
- **1000genomes_hg00157:** stub BAMs (bwa-mem2 2.2.1, about 0-3 KB), flagstat files, an HTML/JSON QC report.
- **hg00157:** a public 1000 Genomes FASTQ pair, about 1 GB each (low coverage or a subset, *inferred* from size).
  Being public, it is the safe fixture candidate here, but it has no calls yet.

### alphagenome (21 files, 367 GB)

AlphaGenome variant-effect prediction TSVs (284 GB, 88 GB and 21 GB tabix-indexed `.tsv.gz` inside zips), PDFs
and notes. These are **annotations of some person's variants, not genome facets.** Out of scope for facet tests.

### Reference panels (one line each)

- **prive:** two 0.9 GB CSV.gz reference tables (ancestry PCA/reference, *inferred*).
- **aadr:** empty on Krypton.
- **dataverse_66p1:** AADR v66 EIGENSTRAT/PLINK-style `.snp`/`.ind`/`.geno` sets (binary 2-7 GB), an mtDNA FASTA,
  READMEs.

---

## 2. Five best private multi-facet sets for live tests

1. **felis (library 12 plus the `vcfs/` chain).** Reads, two alignments to two hg38 builds, DV 1.9.0 on both,
   GATK DRAGEN-mode, Clair3, a consensus, gVCF, annotation-only derivatives, and a second person
   (library 13). This covers more extended-shape axes than any other alias: caller concordance on the same reads,
   the same caller on two references, two seqcol digests, FILTER-`.` handling, and coverage from BAM. Missing: a
   second assembly and an array.
2. **newton.** Three DeepVariant versions (1.1.0, 1.3.0, 1.8.0) for one person, plus an 85 GB BAM and a gVCF.
   This is the version-drift test. One caveat: newton_n1 is empty, and no newton liftover was found.
3. **anton.** One GRCh38 call set lifted by two tools (CrossMap 0.6.3 and Picard-style) in two settings each, with
   rejects and unmapped records. This is the `liftover_to` test, and the record counts reconcile exactly. There is
   no direct GRCh37 call.
4. **sergey.** DRAGEN 4.2.4 vs GATK 3.8 for one person, on two different hg38 contig sets (3366 vs 2580). This tests
   cross-vendor concordance and digest identity. There are no reads.
5. **ksuhaster.** Five sequencing libraries, five BAM pairs (raw and reheadered) and two DeepVariant 1.4.0 runs.
   If it is one person, this is the "same caller, independent sequencing runs" test, plus coverage. Needs the
   person to confirm who is who first.

Runner-up: **olga** (GATK3/4 VQSR source plus Picard-style liftover with and without swap, with rejects). It
complements anton with a different caller. **advocate** and **john** are clean single chains (FASTQ→BAM→DV 1.3.0
gVCF+VCF): good for reads → coverage, not for concordance.

---

## 3. PGP-HMS: the 10 richest participants by facet diversity

Source: `/mnt/krypton/Genomes/genomes/pgp-hms/` (`catalog/sniff_report.tsv`, `participants/manifest.jsonl`,
`catalog/arv_collections_*.json`), plus VCF headers read for classification. "Arvados BAM" and "CGI raw/coverage"
are listed in Arvados but **not fetched**. "PGP fb+VS" is PGP's hg19 freeBayes+VarScan2 merge.

| # | hu | facets (count) | assemblies | callers / producers | extended-shape items |
|---|---|---|---|---|---|
| 1 | **hu832966** | PGP fb+VS VCF; Arvados BAM 32 GB (not fetched); Genos exome (gVCF-banded VCF); Helix exome+ VCF; Gencove low-pass + loimpute imputed VCF; Genes for Good array VCF; 23andMe ×2 (v36 and v37 builds); AncestryDNA; Color genotypes and variant set; Genographic | GRCh36, GRCh37, GRCh38 (Helix) | freeBayes 0.9.15 + VarScan2; Genos pipeline 1.0.0; helix-v3.0.0; loimpute (Gencove); GfG 1.2 | array vs WGS vs exome vs imputed; 3 assemblies; graded evidence (imputed dosages vs called); reads (once fetched) |
| 2 | **hu6ABACE** | Complete Genomics var/masterVar/vcfBeta; DRAGEN 3.2.8 WGS (snp, indel and cnv split files, hs37d5); 2012 GATK UnifiedGenotyper exome; Helix exome+; 23andMe v36; a tar.gz VCF (unclassified) | GRCh36, GRCh37, GRCh38 | CGI; DRAGEN 05.011.281.3.2.8; GATK UG + SnpEff 2.0.5; helix-v3.0.0 | **four independent sequencing assays** (CGI, Illumina DRAGEN, exome, Helix); CGI vs DRAGEN concordance; a split snp/indel/cnv set |
| 3 | **hu2FEC01** | CGI var; 2012 GATK UG exome ×2; FTDNA Big Y (aengine chrY VCF ×2 + a BAM upload); 23andMe v36; AncestryDNA; FTDNA ×4; plink-23 | **GRCh36, GRCh37, GRCh38**, the only one spanning all three | CGI; GATK UG; aengine v2.1.3 | targeted chrY (Big Y) as a facet whose absent sites are not hom-ref; 3 assemblies; reads |
| 4 | **huA2692E** | PGP fb+VS; Veritas gVCF (bz2, 1.3 GB, classified Isaac/Starling-style); Arvados BAM 40 GB; Genos exome; Genes for Good filtered + **imputed and phased**; AncestryDNA; a 178 GB FGC+Veritas+Genos upload (skipped, too large) | GRCh37 | freeBayes 0.9.21 + VarScan2; Veritas (Isaac/Starling, *inferred*); Genos | **two callers on the same Veritas reads** (Veritas native gVCF vs PGP fb+VS, *inferred* same reads); gVCF; statistical phasing |
| 5 | **hu416394** | Sequencing.com 30x WGS (bcftools 1.11 mpileup/call, Canvas 1.39 CNV, Manta-style SV); BAM upload (bwa 0.7.17, GRCh37, 1.7 GB); 23andMe; AncestryDNA; Color; Genes for Good; Full Genomes 57 GB and Sequencing.com 64 GB (skipped, too large) | GRCh37, GRCh38 | bcftools; Canvas; GenerateSVCandidates 1.5 | WGS hg38 vs arrays hg37; SNV + CNV + SV facets; an explicit "missing = not sequenced" header |
| 6 | **hu2B76CC** | PGP fb+VS; Arvados BAM 24 GB; Full-Genomes-style freeBayes WGS (0.9.20, hg19, two annotated variants); chrY UnifiedGenotyper (GATK 3.4); AncestryDNA; FTDNA ×4 | GRCh37 | freeBayes 0.9.15 (PGP) vs freeBayes 0.9.20 (FGC) | **the same caller family at two versions on two sequencings**; chrY-only facet |
| 7 | **huF06AD0** | CGI ×2 deliveries (**software 2.0.4 vs 2.5.0**, arvados + evidence); freeBayes 0.9.20 WGS on hs37d5 (SnpEff + VEP 78); 2012 GATK UG exome ×2; 23andMe v36 ×2; FGC/Novogene 5.2 GB (skipped) | GRCh36, GRCh37 | CGI 2.0.4 / 2.5.0; freeBayes; GATK UG | **CGI version drift**; CGI vs Illumina freeBayes; exome vs WGS |
| 8 | **huE18D82** | PGP fb+VS; DRAGEN upload; Arvados BAM 44 GB; 23andMe | GRCh37 (DRAGEN build not stated in the parsed header) | freeBayes + VarScan2 vs DRAGEN | Illumina caller concordance (fb+VS vs DRAGEN); reads; an array |
| 9 | **hu1B2689** (already in PGP_HMS.md) | CGI var/masterVar/vcfBeta/summary; PGP fb+VS; Arvados BAM 25 GB; 23andMe | GRCh37 | CGI vs freeBayes + VarScan2 | **CGI vs Illumina WGS** (independent technologies) plus an array, one assembly |
| 10 | **hu250634** | CGI ×2 deliveries (**2.5.0.33 vs 2.5.0.37**); Illumina HumanOmni2.5-8 array exported as VCF (`Illumina_GenotypingToVCF` v1.4); 23andMe | GRCh37 | CGI | CGI minor-version drift; array-as-VCF vs CGI |

Close behind:
- **hu5BB600:** FASTQ and BAM uploads, Gencove imputed, Genes for Good imputed, 23andMe ×2.
- **huAF3C63:** CGI, PGP fb+VS, BAM, 23andMe v36.
- **hu1C1368:** DRAGEN, 23andMe, AncestryDNA, Living DNA.
- **huB1488D / hu826751:** CGI, CGI raw coverage (not fetched), 23andMe, AncestryDNA.
- **hu2BA77F:** Isaac/Starling-style **GRCh38** gVCF.
- **huF7A4DE:** GATK HC GRCh38 plus 23andMe and FTDNA.

**CGI version drift is common:** 21 participants have two or more `var` files with different CGI software
versions. Examples: hu3A8D13 1.11 → 2.0.4, hu6E4515 1.12 → 2.0.3, huF06AD0 2.0.4 → 2.5.0, and **hu9385BA** with
three (2.0.1, 2.5.0.33, 2.5.0.37). All CGI files here declare NCBI build 37. Whether a pair is a re-sequencing or a
re-analysis was not checked; the GS sample ids would settle it.

PGP-wide counts of VCFs by producer:

| producer | VCFs |
|---|---|
| PGP fb+VS | 115 |
| 2012 GATK UG exome | 26 |
| unclassified | 18 |
| Genes for Good imputed | 15 |
| Genes for Good array | 12 |
| Isaac/Starling-style | 10 |
| DRAGEN | 10 |
| Gencove imputed | 8 |
| GATK HC | 6 |
| Genos | 5 |
| GATK UG | 3 |
| Big Y | 3 |
| Illumina array VCF | 2 |
| FGC freeBayes | 2 |
| Helix | 2 |
| Sequencing.com | 1 |
| 23andMe-as-VCF | 1 |
| other imputed | 1 |

---

## 4. Gaps: facets nobody has

- **T2T-CHM13:** none, private or PGP. Every file is NCBI36, GRCh37 or GRCh38.
- **Read-backed or HiFi phasing:** none. The only phased data is Genes for Good's statistically imputed and
  phased VCFs (PGP). No PacBio anywhere.
- **Long-read calls:** none. `nanopore` has raw ONT FASTQ only; `oxford` has stub BAMs. No Clair3-ONT, no Dorado,
  no Sniffles output.
- **Trio / pedigree:** none. felis has a second person (library 13) and an unknown third (library 14), with no
  pedigree file and calls only for library 12. The pipeline repo's DeepTrio input is a unit-test chr20
  config, not real data.
- **Coverage summaries (mosdepth, samtools depth):** none on disk. There are only BAMs (felis, newton, advocate,
  john, ksuhaster) for polars-bio depth, and PGP's CGI `coverageRefScore` (not fetched).
- **Direct GRCh37 and direct GRCh38 calls for the same private person:** none. Private GRCh37 exists only as
  liftovers (anton, olga, gennet). PGP has cross-assay pairs instead: hu832966 (Helix hg38 + PGP hg19), hu416394
  (Sequencing.com hg38 + hg19 arrays), hu6ABACE (Helix hg38 + DRAGEN hs37d5 + CGI). Lifting one and comparing it to
  the other works across assays, not on one call set.
- **A newton liftover and newton_n1:** not found; newton_n1 is empty on Krypton.
- **A private array paired with private WGS:** none. angellos is array-only, and no WGS alias has an array.
  PGP covers this well (hu832966, hu1B2689, huE18D82).
- **CNV/SV facets:** only PGP hu416394 (Canvas + SV) and hu6ABACE (DRAGEN cnv). Nothing private.

Open questions for the person:
1. Who are ksuhaster's five libraries?
2. What is felis library 14?
3. Which library fed felis `vcfs/`?
4. Where are newton_n1 and the newton liftovers?
