# The first method-error atlas (RM279)

Built 2026-10-23 with `scripts/method_atlas.py` over eight same-person facet pairs, on GIAB v3.6
stratifications. The two PGP-HMS pairs were rebuilt the same day once `resolve_ref` read CGI's REF
(F20), with VRS ids minted on both sides. The full table (768 rows: every region class, variant class
and scope) is `data/output/method_atlas.parquet` and `.tsv`, git-ignored, with
`method_atlas.pairs.tsv` listing each pair's methods and checks. Everything here is an aggregate over a whole genome: no position, sample
name or file is named, and the private pairs are described by their methods only. Columns and rules
are in docs/SCHEMA.md ("A method-error atlas").

## The pairs

| Set | Pairs | Assembly | What differs |
|---|---|---|---|
| A same-person set on the same reads (one GRCh38 alignment with alts, each output normalized by `bcftools norm -m -both`) | 3: DeepVariant 1.9.0, GATK HaplotypeCaller 4.6.2.0 (DRAGEN mode), Clair3 1.0.10, each against each | GRCh38, REF read | the caller only |
| One person called by three DeepVariant releases | 3: 1.1.0, 1.3.0, 1.8.0, each against each | GRCh38, REF read, not normalized | the release; the alignment behind 1.3.0 and 1.8.0 is not stated |
| Two PGP-HMS participants (public) | 2: Complete Genomics vcfBeta (CGA 2.5.0.24, 2.5.0.44) against freeBayes 0.9.15 / 0.9.20 on Illumina reads | GRCh37, REF read (CGI's chrY left unread, F20), both sides minted, alleles paired by VRS id | the sequencing technology and the caller |

`same_person` (KING-robust, RM274) answered True for all eight, so none was skipped. The
DeepVariant 1.3.0 / 1.8.0 pair has almost no one-sided allele (below 0.01% of SNVs), which fits the
reading in docs/LIVE_SAMPLES.md that those two are re-calls of one alignment; its discordance is then
the release's own.

How alleles met differs by set. The six private pairs were built before minting was in the script:
their alleles met by position and spelling, which holds where both files spell an allele alike (the
same-reads set went through one `bcftools norm`; the DeepVariant releases, though not normalized,
share a caller's spelling). Their counts stand from that build, and `method_atlas.pairs.tsv` leaves
their `pairing` empty. The PGP pairs met by VRS id (`pairing = vrs`), which CGI against freeBayes
needs: the two spell indels and complex alleles differently.

## Headline rates

Scope `carried_pass`: alleles at least one side's genotype carries, every stating row PASS or
unfiltered. Discordant is `mismatch + ploidy_differs` over the alleles both sides called; one-sided
is over all carried alleles. The 95% Wilson intervals are within ±0.1 percentage point for every
`all` and `easy` cell of the private pairs and are in the full table. The PGP pairs have their own
tables below, because most of their `discordant` count is a ploidy convention on X.

**SNV genotype discordance, %** (private pairs)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 0.52 | 0.05 | 2.08 | 3.42 | 2.19 | 0.86 | 4.05 | 3.49 | 4.93 | 1.12 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 2.67 | 0.40 | 9.71 | 19.11 | 12.76 | 2.94 | 8.81 | 9.19 | 14.16 | 4.64 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 2.25 | 0.32 | 8.68 | 18.95 | 11.19 | 2.47 | 7.50 | 7.67 | 12.14 | 3.90 |
| DeepVariant 1.1.0 / 1.3.0 | 0.48 | 0.02 | 1.95 | 3.04 | 2.40 | 1.18 | 7.58 | 3.86 | 2.92 | 0.92 |
| DeepVariant 1.1.0 / 1.8.0 | 0.52 | 0.02 | 2.13 | 2.66 | 2.02 | 1.90 | 8.19 | 5.56 | 3.39 | 1.09 |
| DeepVariant 1.3.0 / 1.8.0 | 0.37 | 0.02 | 1.48 | 2.15 | 1.67 | 1.21 | 4.10 | 3.27 | 2.32 | 0.71 |

**Indel genotype discordance, %** (private pairs)

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
0.980 to 0.984 for the DeepVariant releases, 0.966 and 0.969 for the PGP pairs. Clair3's rows carry
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

## The PGP pairs: CGI against freeBayes

The first build compared these two pairs by `assembly_label` only, with alleles meeting by position
and spelling: `resolve_ref` refused both CGI vcfBeta facets over a chrY md5 that differs from the
store's GRCh37. Since F20 a differing X, Y or MT leaves that contig unread instead, so all four facets
now have REF read from the store's Ensembl GRCh37 collection (every read row matches its source REF
but 3, against IUPAC codes; CGI's chrY, 1,744 and 1,710 rows, stays unread), and `mint_vrs_ids`
mints them (99.96% and over 99.99% of alleles on one pair, 99.96% and 95.1% on the other, whose
freeBayes file has 256,986 rows on `hs37d5` and `NC_007605`, contigs the store does not have). Both pairs then met allele by allele through VRS
ids (`assembly_check = content`, `pairing = vrs`).

Two conventions sit inside these pairs' counts and are not method error:

- **X outside the PARs is haploid in CGI and diploid in freeBayes** (every genotype there: CGI writes
  one allele, freeBayes two). Every allele both call there is `ploidy_differs`: 72,462 and 64,818
  SNVs, 65% and 67% of the pairs' discordant SNVs. The tables here therefore give **mismatch over
  called both**: `ploidy_differs` is left out of the count and stays in `called_both`. The atlas's
  `discordance_rate` (3.40% and 3.07% for SNVs) is in the full table.
- **A complex freeBayes row can hold the allele under another spelling.** freeBayes writes a
  haplotype as one multi-allelic row: REF `CGCCC`, ALT `TGCCT,TGCCC`, GT `TGCCT/TGCCT`. Its ALT
  `TGCCC` trims to `C>T`, so its VRS id is CGI's `C>T` (GT `T/T`), and the two meet. `concordance`
  then counts the copies of `TGCCC` in the freeBayes genotype, which are zero, although each
  `TGCCT` haplotype carries the same `C>T`: a `mismatch` that is spelling. Counted on the VRS-paired
  mismatches by looking for a carried sibling ALT in such a row (REF longer than one base): 20,434
  SNVs in the CGA 2.5.0.44 pair whose sibling carries exactly the allele's change, 23,014 with any
  carried sibling; 2,601 and 2,950 in the 2.5.0.24 pair; and up to 5,982 and 4,627 indels (for an
  indel the sibling is not checked base by base, so that is an upper bound). The script that counted
  this is not part of the repo. Fixing it is a rule for `concordance`, not the atlas: the `site` path
  already refuses to read a row over another REF span, and the `vrs` path has no such guard.

**Before and after, region `all`, scope `carried_pass`, %.** "By position" is the first build (label
check, alleles by position and spelling); the REF it did not read changes nothing there, since
alleles fall back to the source REF. "Less the artefact" takes the complex-row mismatches (any
carried sibling) out of both counts.

| Methods (a / b) | Class | By position: carried | mismatch | one-sided | By VRS id: carried | called both | mismatch | less the artefact | one-sided | of which silent side reads reference |
|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | SNV | 4,441,897 | 0.49 | 25.97 | 4,416,242 | 3,281,233 | 1.19 | 0.49 | 24.93 | 55 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | SNV | 4,685,178 | 0.93 | 31.98 | 4,651,492 | 3,160,612 | 1.02 | 0.93 | 31.35 | 56 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | indel | 1,440,360 | 44.70 | 99.13 | 987,570 | 410,440 | 5.31 | 3.91 | 52.87 | 60 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | indel | 1,193,130 | 53.12 | 99.03 | 819,248 | 342,076 | 8.05 | 6.79 | 52.93 | 60 |

Pairing by VRS id is what the indels needed. The two spellings of one indel used to count as two
one-sided alleles: in the CGA 2.5.0.44 pair carried indels fall by a third (1.44 million to 0.99
million) and the indels both facets call rise from 11,891 to 410,440; in both pairs one-sided indels
drop from 99% to 53%, 16 to 18% in GIAB's easy regions. SNVs barely move: with the complex-row
mismatches taken out, the mismatch rate is the first build's to two decimals (an SNV's VRS id and
its position and spelling name one allele), which says the new pairing added no other error.
`uncalled` grows with it (indels 420 to 43,543 in the CGA 2.5.0.24 pair, as CGI's half-calls now meet
freeBayes), so `called_both` is less than the carried alleles minus the one-sided ones. MNVs are in
the full table; complex rows are a larger share of them, and their mean `p_same` falls (0.93 to 0.66
in the CGA 2.5.0.44 pair) for the same reason.

The artefact is counted for region `all` only. The VRS-paired SNV rates per region carry it in every
class, so the SNV table below is the first build's, by position, which has no artefact and agrees
with the VRS pairing once it is taken out (above). The indel table is the VRS pairing's, complex rows
still in.

**SNV mismatch, % of called both** (by position, the first build; `ploidy_differs` left out)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | 0.49 | 0.29 | 1.60 | 2.70 | 1.63 | 1.25 | 10.75 | 3.12 | 2.53 | 0.67 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | 0.93 | 0.55 | 3.08 | 3.70 | 1.89 | 3.25 | 24.69 | 4.35 | 2.86 | 2.35 |

**Indel mismatch, % of called both** (by VRS id, complex rows still in; `ploidy_differs` left out)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | 5.31 | 1.54 | 7.90 | 6.72 | 6.10 | 4.37 | 14.61 | 12.65 | 9.78 | 4.75 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | 8.05 | 3.03 | 12.52 | 9.42 | 7.60 | 8.09 | 26.57 | 15.31 | 9.72 | 8.88 |

**One-sided, % of carried** (by VRS id)

| Methods (a / b) | Class | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | SNV | 24.93 | 9.81 | 59.49 | 76.53 | 72.47 | 29.74 | 76.20 | 61.86 | 85.13 | 42.75 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | SNV | 31.35 | 14.50 | 61.59 | 74.51 | 66.77 | 52.58 | 89.72 | 65.28 | 81.42 | 48.44 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | indel | 52.87 | 16.13 | 62.27 | 82.66 | 78.69 | 29.93 | 67.09 | 84.76 | 87.74 | 59.76 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | indel | 52.93 | 17.96 | 63.26 | 81.30 | 74.85 | 38.82 | 68.95 | 85.98 | 84.68 | 55.53 |

The 95% Wilson intervals of these `all` and `easy` cells are within ±0.2 percentage point.

Reading them:

- **Easy regions mostly agree across technologies.** Leaving X's ploidy out, CGI and freeBayes
  disagree on 0.29 to 0.55% of SNV genotypes in easy regions, around the same-reads callers' rate
  (0.05 to 0.40%), and on 1.5 to 3.0% of indel genotypes with complex rows still in. Homopolymers of
  12 or more bases are again the worst class, at 11 to 27%.
- **One-sided alleles are the larger difference.** A quarter to a third of carried SNVs, and half the
  indels, are stated by one side only, rising to 75 to 90% in low mappability, segmental duplications
  and long tandem repeats. Most are freeBayes's (`only_b` is 6 times `only_a` for SNVs), and in
  55 to 60% of them the silent side reads as reference (`one_sided_reference`: a reference block over
  the position, or a whole-genome facet's absence rule), a stated disagreement, where the rest is
  silence. In the CGA 2.5.0.24 pair, 262,836 freeBayes alleles (before the carried and PASS
  filters) lie on `hs37d5` and `NC_007605`, which the CGI facet never declared; they are one-sided
  with `absence = unknown` and fall in `all` only.

X and Y classes are in the full table only: for a private pair they would state an aggregate about
the person's sex chromosomes.
