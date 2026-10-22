# PGP-HMS as a test-data source

Assessment notes from 2026-09-28/29. The Harvard Personal Genome Project is the source of consented,
public test data for just-vcf, and especially of **people with several independent facets** (a
consumer array export, a PGP-sequenced VCF, a Complete Genomics `var` file). Importing each facet
correctly should make them converge wherever they overlap, which makes them the strongest test setups.

Tools: `scripts/pgp_leech.py` (listing, one data type or a per-type sample) and
`scripts/pgp_participants.py` (everything per participant). Data lives under
`$JUST_VCF_SAMPLES_DIR/pgp-hms/`. Catalog snapshots taken during the assessment are in
`/data/genomes/pgp-hms/catalog/`, not in git, because the labels carry participants' real names.

## Licence and manners

- CC0. PGP Researcher FAQ: "All genomic and health data available on the Harvard Personal Genome
  Project page are under a CC0 license allowing unencumbered and free use by anyone around the world."
  (https://github.com/PGPHarvard/FAQ/wiki/Researcher-FAQ)
- `my.pgp-hms.org/robots.txt` asks for `Crawl-delay: 5`. Both scripts wait 5 s between requests to one
  host.
- Labels and file names can carry real names. Check a file before it becomes an `assets/` fixture.

## Where the data is

| Source | What | Access |
|---|---|---|
| Listing `my.pgp-hms.org/public_genetic_data` | 3,960 rows: 3,171 participant uploads + 789 PGP-produced (531 Complete Genomics, 258 microbiome) | server-rendered HTML, complete (no pagination, no default filter) |
| `public_genetic_data.json` | same rows | lacks participant id, size and locator; PGP rows have `data_type: null` |
| Participant uploads | `/user_file/download/<id>` | see "Download behaviour" |
| PGP Complete Genomics rows | 381 link to `evidence.pgp-hms.org/genome_download.php` (GET-Evidence); 150 have no link at all | GET-Evidence works over https; the 150 are not publicly reachable |
| Arvados cluster `ac2it` | 2,571 collections visible anonymously, in 700 projects | API with the anonymous token the cluster publishes at `/arvados/v1/config`; files via keep-web `https://<pdh with + → ->.collections.ac2it.arvadosapi.com/_/<path>`, public without a token |
| DB snapshots `HU-PGP-DB-*` (20, latest 2021-06-18) | `uploaded_data.tsv` is a copy of the listing; other tables are surveys, conditions, specimens | keep-web |

Rough sizes: participant uploads 7.0 TB; Arvados per-participant datasets 5.85 TB (BAM 4.78 TB, VCF 53 GB,
CGI ASM 0.13 TB, other 0.88 TB). Arvados also holds 1,980 copies of the uploads (2.28 TB, which does not
reconcile with the listing's 7.0 TB; not investigated).

## Arvados per-participant datasets

Projects are named `huXXXXXX datasets` and `huXXXXXX reports` (339 participants).

- **VCF collections (109, 108 participants).** All four naming styles (Veritas-like `…-EXT`, 14-digit ids,
  7-character ids, `WGC…D`) come from one PGP pipeline, run 2016–2018 on hg19: freeBayes
  (v0.9.15/0.9.21) and VarScan2 per chromosome, merged with `vcf-merge`. **Which caller made a row is in
  INFO/SF (0 = freeBayes, 1 = VarScan2), not in the sample columns.** Checked 2026-09-29 over the
  arvados merges: 13 of 110 carry two populated sample columns; 48 declare one sample; 50 declare two
  (`unknown`, `Sample1`) but fill only one per row, so where both callers found a site VarScan2's
  genotype is gone. An earlier note here said every file kept one column per caller; that was wrong.
  Each collection also has a `…-orig` key=value table of both callers' fields (spec
  `freebayes_varscan2_keyvalue`), which still holds both.
- **Complete Genomics deliveries (68 `GS…-DID`, 55 "Raw sequencing data" folders, some "copy of").**
  Native CGI format: `var`, `masterVarBeta`, CGI's own `vcfBeta`, `summary`, `dbSNPAnnotated`, gene and
  CNV files, and in the raw folders per-chromosome `coverageRefScore`/`evidence*` files (about 35 GB per
  person). That coverage could serve as a callability facet later.
- **BAMs (109)**, including 26 GB `WGC….bam.tgz`. Not fetched.
- **Reports (222):** GET-Evidence and ClinVar HTML. Not fetched.

## Data types are free text

The upload form offers 16 standard types plus "other"; choosing "other" stores whatever the participant
typed. 360 rows carry 243 free-text strings, confirmed by `?data_type=other` returning exactly 360.

Classification of those 360 rows (Haiku pass plus two rule passes; every rule override is marked in the
note column) is in `catalog/freetext_classified.tsv`:

| content | rows |
|---|---|
| genotype table | 122 (AncestryDNA 76 under ~40 spellings) |
| variant calls | 77 (Dante, Genes for Good incl. imputed, Full Genomes, Genos, DNA.Land imputed VCF/BCF, Nebula, Sequencing.com, Helix gVCF) |
| other document | 41 |
| reads or alignments | 34 |
| derived report | 33 |
| unknown, needs content sniffing | 21 |
| image or media | 17 |
| health or phenotype | 15 |

**Sniffing the 21 unsettled rows (2026-09-29).** Nine are on disk. Seven of those are variant data:
four VCFs (one of them a Genos gVCF with sorted meta lines), Gencove's Oxford `.gen` imputation zip,
and a SelfDecode GRCh38 table. The other two are not variant data: an iWork document and a gzipped
FASTQ in a zip. The twelve not on disk: five were over the 2.5 GB cap, one is not public, and six have no
manifest entry yet (not planned, or not reached by the fetch still running; not checked which). Reading their heads over HTTP would
settle them without downloading them.

**The label is unreliable both ways.** A `23andMe`-typed upload can be a VCF, and real genotype files sit
under free text. Standard types hide reads as well: a "Veritas Genetics" person is one `<id> - VCF` row
plus one `<id>.chrN.bam - BAM` row per chromosome. The fetchers filter by label and by name for now.
Content-based format detection is the real answer and belongs in the library (RM2).

## Download behaviour (all seen live)

- `/user_file/download/<id>` answers three ways: served directly (older uploads, 200 with
  Content-Disposition); a 302 to an Arvados collection's `/_/` page listing its files; or a 302 straight
  to a file in the collection.
- **401 from the file host**: listed but not shared (29 recent uploads). Recorded as `not_public`.
- **The listing's size column is sometimes wrong.** Five Veritas `- VCF` uploads were listed at 1 byte or
  20 bytes short, and were byte-identical to the person's Arvados VCF.
- Two old uploads return an HTML page instead of a file.
- GET-Evidence: links are `http://` and redirect to https; there is no Content-Length; the body arrives
  with `Content-Encoding: gzip` over a bz2 payload. Integrity is checked by decompressing the whole file.
  Some rows return 404.
- The participant cell can read `huF06AD0, PGP113` (the original enrolment number); parsed as a
  separate field.
- Arvados list queries with a long owner filter need POST with `_method=GET` (the GET gets 414).

## What is on disk

`/data/genomes/pgp-hms/participants/<hu>/<source>-<id>/<file>` with `participants/manifest.jsonl` (one
record per file: outcome, path, reason, and what a duplicate duplicates).

- Sources: `arvados-<collection uuid>` (VCF, CGI `var`/`masterVarBeta`/`vcfBeta`/`summary`),
  `upload-<download id>`, `evidence-<genome id>` (GET-Evidence CGI).
- Duplicates within a person are fetched once: same name and size; an upload of the same size as a
  dataset file; an upload labelled `<id> - VCF` matching dataset file `<id>.vcf.gz`; a GET-Evidence file
  whose header names a file the person already has.
- 2026-09-29 03:20, second run in progress: 632 participants planned, ~1,500 files present, 212 GB.
  Rerunning `uv run python scripts/pgp_participants.py` resumes; present files are skipped offline.
- Not fetched: files over 2.5 GB (19), BAMs, reports, per-chromosome CGI evidence/coverage.
- **One damaged upload: huE2F2FC download 2771** (an Illumina starling genome VCF, 886,854,046 bytes,
  the size the listing states). Walked block by block, one 163,477-byte span at offsets
  473,789,118–473,952,595 does not inflate ("invalid distance code"), and every block before and after
  it is intact, so `load` fails in polars-bio with "invalid BGZF ISIZE", which is right. The companion
  index (download 2772) points at nine block starts inside that span, none of them a valid block in this
  copy, while 2,000 of 2,000 checked outside it are: the index was built on bytes this copy does not have.
  So the damage came after indexing, in the upload, at PGP or in the download; a re-fetch to tell which
  failed on 2026-10-22 (my.pgp-hms.org served a certificate for another host, to httpx and curl alike).
  The sibling download 2774 (the same genome's variant VCF) is intact and loads, 4,763,188 calls.
- The by-data-type layout of the first runs was migrated into `participants/`;
  `manifest.legacy-by-type.jsonl` is its old manifest.

## Good first test setups

- **hu1B2689**: 23andMe export + 4 Complete Genomics files + PGP VCF.
- **huF5F075**: Complete Genomics + PGP VCF.
- **hu832966**: 23andMe, AncestryDNA, Color, Gencove (low-pass + imputation), Genos exome, Helix, plus
  PGP VCF.
- **huA2692E**: Veritas gVCF, Genes for Good, Genos exome, AncestryDNA, Full Genomes, plus PGP VCF.
- Any PGP VCF on its own: freeBayes versus VarScan2 in one file.
