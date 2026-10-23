# The first method-error atlas (RM279)

Built 2026-10-23 with `scripts/method_atlas.py` over eight same-person facet pairs, on GIAB v3.6
stratifications. The full table (768 rows: every region class, variant class and scope) is
`data/output/method_atlas.parquet` and `.tsv`, git-ignored, with `method_atlas.pairs.tsv` listing each
pair's methods and checks. Everything here is an aggregate over a whole genome: no position, sample
name or file is named, and the private pairs are described by their methods only. Columns and rules
are in docs/SCHEMA.md ("A method-error atlas").

## The pairs

| Set | Pairs | Assembly | What differs |
|---|---|---|---|
| A same-person set on the same reads (one GRCh38 alignment with alts, each output normalized by `bcftools norm -m -both`) | 3: DeepVariant 1.9.0, GATK HaplotypeCaller 4.6.2.0 (DRAGEN mode), Clair3 1.0.10, each against each | GRCh38, REF read | the caller only |
| One person called by three DeepVariant releases | 3: 1.1.0, 1.3.0, 1.8.0, each against each | GRCh38, REF read | the release; the alignment behind 1.3.0 and 1.8.0 is not stated |
| Two PGP-HMS participants (public) | 2: Complete Genomics vcfBeta (CGA 2.5.0.24, 2.5.0.44) against freeBayes 0.9.15 / 0.9.20 on Illumina reads | GRCh37, compared by label | the sequencing technology and the caller |

`same_person` (KING-robust, RM274) answered True for all eight, so none was skipped. The
DeepVariant 1.3.0 / 1.8.0 pair has almost no one-sided allele (below 0.01% of SNVs), which fits the
reading in docs/LIVE_SAMPLES.md that those two are re-calls of one alignment; its discordance is then
the release's own.

## Headline rates

Scope `carried_pass`: alleles at least one side's genotype carries, every stating row PASS or
unfiltered. Discordant is `mismatch + ploidy_differs` over the alleles both sides called; one-sided
is over all carried alleles. The 95% Wilson intervals are within ±0.1 percentage point for every
`all` and `easy` cell of the private pairs and are in the full table.

**SNV genotype discordance, %**

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 0.52 | 0.05 | 2.08 | 3.42 | 2.19 | 0.86 | 4.05 | 3.49 | 4.93 | 1.12 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 2.67 | 0.40 | 9.71 | 19.11 | 12.76 | 2.94 | 8.81 | 9.19 | 14.16 | 4.64 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 2.25 | 0.32 | 8.68 | 18.95 | 11.19 | 2.47 | 7.50 | 7.67 | 12.14 | 3.90 |
| DeepVariant 1.1.0 / 1.3.0 | 0.48 | 0.02 | 1.95 | 3.04 | 2.40 | 1.18 | 7.58 | 3.86 | 2.92 | 0.92 |
| DeepVariant 1.1.0 / 1.8.0 | 0.52 | 0.02 | 2.13 | 2.66 | 2.02 | 1.90 | 8.19 | 5.56 | 3.39 | 1.09 |
| DeepVariant 1.3.0 / 1.8.0 | 0.37 | 0.02 | 1.48 | 2.15 | 1.67 | 1.21 | 4.10 | 3.27 | 2.32 | 0.71 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | 2.71 | 2.38 | 4.45 | 7.19 | 6.15 | 3.71 | 13.23 | 5.19 | 4.97 | 2.55 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | 2.98 | 2.48 | 5.85 | 8.27 | 6.05 | 5.53 | 26.81 | 6.27 | 4.87 | 4.03 |

**Indel genotype discordance, %** (private pairs; the PGP pairs' indels are not comparable, below)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 3.58 | 0.67 | 4.29 | 8.09 | 6.74 | 1.42 | 5.64 | 5.73 | 8.69 | 4.91 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 7.84 | 2.03 | 9.25 | 20.98 | 16.62 | 4.45 | 10.67 | 12.12 | 19.64 | 12.28 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 7.13 | 1.44 | 8.51 | 21.07 | 16.00 | 4.13 | 10.35 | 9.89 | 16.55 | 10.92 |
| DeepVariant 1.1.0 / 1.3.0 | 7.13 | 0.28 | 8.86 | 10.24 | 9.41 | 2.42 | 14.65 | 8.99 | 7.70 | 7.63 |
| DeepVariant 1.1.0 / 1.8.0 | 7.95 | 0.53 | 9.82 | 10.51 | 9.82 | 3.29 | 14.39 | 13.25 | 10.21 | 9.03 |
| DeepVariant 1.3.0 / 1.8.0 | 3.68 | 0.31 | 4.48 | 5.87 | 5.44 | 1.75 | 6.06 | 6.62 | 5.61 | 4.29 |

**One-sided SNVs, % of carried** (private pairs)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | TR 201–10k bp |
|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 2.38 | 0.27 | 8.71 | 14.60 | 9.04 | 20.90 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 6.77 | 0.68 | 21.68 | 30.37 | 25.89 | 51.32 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 2.05 | 0.19 | 7.77 | 11.78 | 5.71 | 28.48 |
| DeepVariant 1.1.0 / 1.3.0 | 1.33 | 0.03 | 5.26 | 12.05 | 4.25 | 6.10 |
| DeepVariant 1.1.0 / 1.8.0 | 1.23 | 0.03 | 4.93 | 11.03 | 4.10 | 7.15 |
| DeepVariant 1.3.0 / 1.8.0 | 0.00 | 0.00 | 0.01 | 0.01 | 0.01 | 0.03 |

Mean `p_same` over the graded SNVs (`gl_flat_prior`, `all`): 0.985 for DeepVariant 1.9.0 / GATK,
0.980 to 0.984 for the DeepVariant releases, 0.967 and 0.970 for the PGP pairs. Clair3's rows carry
no likelihood vector the comparison could use, so its pairs have no graded allele.

## Reading it

- **Easy regions agree.** On the same reads, any two of the three callers disagree on fewer than 0.4%
  of SNV genotypes in GIAB's easy regions, and DeepVariant releases on 0.02%. Low mappability,
  segmental duplications, long tandem repeats and homopolymers of 12 or more bases carry most of
  every pair's disagreement, at several to over a hundred times the easy rate.
- **GATK is the odd one out on the same reads.** DeepVariant and Clair3 agree far better with each
  other (0.52% SNV, 3.58% indel) than either does with GATK in DRAGEN mode (2.25 to 2.67%, 7.1 to
  7.8%), and GATK's disagreement concentrates in low mappability (19%). Every GATK row here is
  unfiltered as written (its FILTER is `.` everywhere), so `carried_pass` keeps all of them while it
  drops DeepVariant's `RefCall` rows.
- **DeepVariant 1.1.0 indels differ from the later releases** (7.1 to 8.0% against 3.7% between 1.3.0
  and 1.8.0), mostly in homopolymers of 12 or more bases (14.5%).
- **One pair per row.** Every row pools one pair (`pairs` = 1). The interval is sampling error over
  the alleles compared, not the spread between people or runs; it says nothing about how a method
  behaves on another genome.

## What the PGP rows cannot say yet

The two CGI vcfBeta facets state a chrY whose md5 differs from the store's GRCh37 collection, so
`resolve_ref` refuses them and the pairs were compared by their equal `assembly_label` (the script
records `assembly_check = label`). Without REF read there are no VRS ids, so alleles meet only by
position and spelling. CGI and freeBayes spell indels and complex alleles differently, and 99% of the
PGP indels come out one-sided: that is representation, not method error, and the indel and one-sided
columns of those pairs are left out above. Their SNV discordance (alleles both called at one
position) stands. Reading CGI's REF needs a GRCh37 collection matching CGI's own reference, or a
per-contig resolve that leaves chrY unread; either is a library question, not the atlas's.

X and Y classes are in the full table only: for a private pair they would state an aggregate about
the person's sex chromosomes.
