# The first method-error atlas (RM279)

Built 2026-10-23 with `scripts/method_atlas.py` over eight same-person facet pairs, on GIAB v3.6
stratifications. The two PGP-HMS pairs were rebuilt the same day once `resolve_ref` read CGI's REF
(F20), with VRS ids minted on both sides, and again once `concordance` read complex rows by the
allele's change (F21) and a haploid call against the same allele's homozygote as a match (F22). The
six private pairs were rebuilt the same way on 2026-10-26: REF read, VRS ids minted on both sides,
alleles met by VRS id, with both rules and the band cells. The
full table (768 rows: every region class, variant class and scope) is `data/output/method_atlas.parquet`
and `.tsv`, git-ignored, with `method_atlas.pairs.tsv` listing each pair's methods and checks; the
last CGI build is `data/output/method_atlas_cgi.*` beside it (192 pooled rows; rebuilt 2026-10-24 with the cells of both sides' confidence bands beside them, 10,552 rows, the pooled ones unchanged), all git-ignored. The private pairs' rebuild (17,896 rows with the band cells) was written to `data/output/method_atlas_private.*` on the builder's branch and not kept; the counts here are taken from it. The full table's private rows are the first build's, by position. Everything here is an aggregate over a whole genome: no position, sample
name or file is named, and the private pairs are described by their methods only. Columns and rules
are in docs/SCHEMA.md ("A method-error atlas").

## The pairs

| Set | Pairs | Assembly | What differs |
|---|---|---|---|
| A same-person set on the same reads (one GRCh38 alignment with alts, each output normalized by `bcftools norm -m -both`) | 3: DeepVariant 1.9.0, GATK HaplotypeCaller 4.6.2.0 (DRAGEN mode), Clair3 1.0.10, each against each | GRCh38, REF read, both sides minted, alleles paired by VRS id | the caller only |
| One person called by three DeepVariant releases | 3: 1.1.0, 1.3.0, 1.8.0, each against each | GRCh38, REF read, not normalized, both sides minted, alleles paired by VRS id | the release; the alignment behind 1.3.0 and 1.8.0 is not stated |
| Two PGP-HMS participants (public) | 2: Complete Genomics vcfBeta (CGA 2.5.0.24, 2.5.0.44) against freeBayes 0.9.15 / 0.9.20 on Illumina reads | GRCh37, REF read (CGI's chrY left unread, F20), both sides minted, alleles paired by VRS id | the sequencing technology and the caller |

`same_person` (KING-robust, RM274) answered True for all eight, so none was skipped. The
DeepVariant 1.3.0 / 1.8.0 pair has almost no one-sided allele (below 0.01% of SNVs), which fits the
reading in docs/LIVE_SAMPLES.md that those two are re-calls of one alignment; its discordance is then
the release's own.

All eight pairs now meet allele by allele through VRS ids (`assembly_check = content`,
`pairing = vrs`). The six private pairs were first built before minting was in the script, their
alleles meeting by position and spelling; they were rebuilt with minting on 2026-10-26, and what that
changed is in "The private pairs met by VRS id" below. The PGP pairs need it most: CGI and freeBayes
spell indels and complex alleles differently.

## Headline rates

Scope `carried_pass`: alleles at least one side's genotype carries, every stating row PASS or
unfiltered. Discordant is `mismatch` over the alleles both sides called (`match + mismatch`); a
haploid call against the same allele's diploid homozygote is a match, its ploidy counted apart in
`ploidy_differs` (F22). One-sided is over all carried alleles. The 95% Wilson intervals are within
±0.1 percentage point for every `all` and `easy` cell of the private pairs and are in the full table.
The private pairs' rates below are the 2026-10-26 rebuild (REF read, alleles met by VRS id, F21 and
F22). Their pairs have no allele whose two genotypes differ in ploidy (`ploidy_differs` is 0 in every
row, before and after), so the F22 rule leaves them as they are. The PGP pairs have their own tables
below.

**SNV genotype discordance, %** (private pairs)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 0.52 | 0.05 | 2.08 | 3.42 | 2.19 | 0.86 | 4.04 | 3.49 | 4.93 | 1.12 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 2.67 | 0.40 | 9.70 | 19.11 | 12.76 | 2.94 | 8.79 | 9.19 | 14.16 | 4.64 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 2.25 | 0.32 | 8.68 | 18.95 | 11.19 | 2.46 | 7.49 | 7.67 | 12.14 | 3.90 |
| DeepVariant 1.1.0 / 1.3.0 | 0.48 | 0.02 | 1.93 | 3.04 | 2.40 | 1.16 | 7.41 | 3.77 | 2.92 | 0.92 |
| DeepVariant 1.1.0 / 1.8.0 | 0.52 | 0.02 | 2.11 | 2.66 | 2.02 | 1.88 | 8.01 | 5.48 | 3.39 | 1.09 |
| DeepVariant 1.3.0 / 1.8.0 | 0.38 | 0.02 | 1.48 | 2.15 | 1.67 | 1.21 | 4.10 | 3.28 | 2.33 | 0.71 |

**Indel genotype discordance, %** (private pairs)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 3.56 | 0.67 | 4.26 | 8.03 | 6.72 | 1.39 | 5.60 | 5.69 | 8.63 | 4.86 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 7.82 | 2.03 | 9.22 | 20.93 | 16.58 | 4.43 | 10.63 | 12.07 | 19.58 | 12.23 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 7.12 | 1.44 | 8.50 | 21.07 | 16.00 | 4.13 | 10.35 | 9.88 | 16.55 | 10.92 |
| DeepVariant 1.1.0 / 1.3.0 | 7.96 | 0.31 | 9.85 | 11.30 | 10.55 | 2.53 | 16.45 | 9.43 | 7.84 | 8.24 |
| DeepVariant 1.1.0 / 1.8.0 | 8.69 | 0.57 | 10.69 | 11.50 | 10.81 | 3.38 | 15.95 | 13.66 | 10.40 | 9.59 |
| DeepVariant 1.3.0 / 1.8.0 | 3.78 | 0.31 | 4.60 | 6.17 | 5.76 | 1.77 | 6.27 | 6.69 | 5.66 | 4.36 |

**One-sided SNVs, % of carried** (private pairs)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | TR 201–10k bp |
|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | 2.38 | 0.27 | 8.71 | 14.60 | 9.04 | 20.90 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | 6.77 | 0.68 | 21.68 | 30.37 | 25.89 | 51.32 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | 2.05 | 0.19 | 7.77 | 11.78 | 5.71 | 28.48 |
| DeepVariant 1.1.0 / 1.3.0 | 1.32 | 0.03 | 5.22 | 12.04 | 4.24 | 6.07 |
| DeepVariant 1.1.0 / 1.8.0 | 1.22 | 0.03 | 4.88 | 11.02 | 4.09 | 7.10 |
| DeepVariant 1.3.0 / 1.8.0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 |

Mean `p_same` over the graded SNVs (`gl_flat_prior`, `all`): 0.985 for DeepVariant 1.9.0 / GATK,
0.980 to 0.984 for the DeepVariant releases, 0.985 and 0.990 for the PGP pairs. Clair3's rows carry
no likelihood vector the comparison could use, so its pairs have no graded allele (its rows still
state a GQ, so they are banded below).

## Reading it

- **Easy regions agree.** On the same reads, any two of the three callers disagree on fewer than 0.4%
  of SNV genotypes in GIAB's easy regions, and DeepVariant releases on 0.02%. Low mappability,
  segmental duplications, long tandem repeats and homopolymers of 12 or more bases carry most of
  every pair's disagreement, at several to over a hundred times the easy rate.
- **GATK is the odd one out on the same reads.** DeepVariant and Clair3 agree far better with each
  other (0.52% SNV, 3.56% indel) than either does with GATK in DRAGEN mode (2.25 to 2.67%, 7.1 to
  7.8%), and GATK's disagreement concentrates in low mappability (19%). Every GATK row here is
  unfiltered as written (its FILTER is `.` everywhere), so `carried_pass` keeps all of them while it
  drops DeepVariant's `RefCall` rows.
- **DeepVariant 1.1.0 indels differ from the later releases** (8.0 to 8.7% against 3.8% between 1.3.0
  and 1.8.0), mostly in homopolymers of 12 or more bases (16.0 to 16.5%). Part of it was hidden by
  spelling until the alleles met by VRS id (below).
- **One pair per row.** Every row pools one pair (`pairs` = 1). The interval is sampling error over
  the alleles compared, not the spread between people or runs; it says nothing about how a method
  behaves on another genome.

## The private pairs met by VRS id (2026-10-26)

The six private pairs were rebuilt with `scripts/method_atlas.py` as the PGP pairs were: each facet's
REF read by `resolve_ref` from the store's GRCh38 collection, VRS ids minted on both sides, alleles met
by VRS id, with the F21 and F22 rules and the band cells. All six facets had their REF read on every
row, and every source REF matched the store's (no mismatch on any facet, no contig left unread).
`mint_vrs_ids` minted every allele, except GATK's 18,248 spanning-deletion `*` alleles, which name an
allele called at another record and have no id of their own. `same_person` answered True for all six
again, and every pair is `assembly_check = content`, `pairing = vrs`.

**Before and after, region `all`, scope `carried_pass`.** "By position" is the first build (REF read,
alleles met by position and spelling, before F21 and F22); "by VRS id" is the rebuild.

| Methods (a / b) | Class | Carried, by position → by VRS id | Called both | Discordance % | One-sided % | `unreadable` (new) |
|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | SNV | 3,789,408 → 3,789,408 | 3,699,288 → 3,699,277 | 0.52 → 0.52 | 2.38 → 2.38 | 11 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | SNV | 4,023,711 → 4,023,711 | 3,751,266 → 3,751,251 | 2.67 → 2.67 | 6.77 → 6.77 | 15 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | SNV | 3,757,154 → 3,757,154 | 3,680,178 → 3,680,169 | 2.25 → 2.25 | 2.05 → 2.05 | 9 |
| DeepVariant 1.1.0 / 1.3.0 | SNV | 3,828,921 → 3,828,588 | 3,777,971 → 3,777,862 | 0.48 → 0.48 | 1.33 → 1.32 | 198 |
| DeepVariant 1.1.0 / 1.8.0 | SNV | 3,808,266 → 3,807,925 | 3,761,355 → 3,761,243 | 0.52 → 0.52 | 1.23 → 1.22 | 214 |
| DeepVariant 1.3.0 / 1.8.0 | SNV | 3,788,298 → 3,788,257 | 3,788,161 → 3,788,191 | 0.37 → 0.38 | 0.00 → 0.00 | 32 |
| Clair3 1.0.10 / DeepVariant 1.9.0 | indel | 968,605 → 968,507 | 920,446 → 920,318 | 3.58 → 3.56 | 4.97 → 4.97 | 49 |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | indel | 1,040,463 → 1,040,344 | 924,455 → 924,335 | 7.84 → 7.82 | 11.15 → 11.15 | 36 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | indel | 989,611 → 989,609 | 924,240 → 924,214 | 7.13 → 7.12 | 6.61 → 6.61 | 26 |
| DeepVariant 1.1.0 / 1.3.0 | indel | 993,577 → 979,597 | 905,877 → 920,421 | 7.13 → 7.96 | 8.83 → 6.01 | 292 |
| DeepVariant 1.1.0 / 1.8.0 | indel | 1,005,119 → 991,042 | 912,934 → 928,041 | 7.95 → 8.69 | 9.17 → 6.33 | 299 |
| DeepVariant 1.3.0 / 1.8.0 | indel | 949,609 → 947,519 | 943,696 → 945,550 | 3.68 → 3.78 | 0.62 → 0.20 | 37 |

The 95% Wilson intervals of the rebuilt `all` and `easy` cells are within ±0.07 percentage point for
discordance and ±0.1 for one-sided.

- **The same-reads set barely moves.** One `bcftools norm` had already spelled the three callers'
  alleles alike, so VRS ids pair no allele that position missed: carried SNVs are the same to the
  allele, carried indels fall by 2 to 119. What changes is F21: 9 to 49 alleles per class meet a
  genotype holding another allele of its row over a different span and are now `unreadable`, outside
  `called_both`. Every such SNV was a mismatch before (SNV mismatches fall by exactly that count);
  indel mismatches fall by 26 to 261. No rate moves by more than 0.02 point.
- **DeepVariant 1.1.0 spells about 14,000 indels unlike the later releases.** The releases are not
  normalized, and by position each such indel counted twice, once one-sided on each side. Met by VRS
  id, 1.1.0's one-sided indels against 1.3.0 fall from 42,803 to 28,817 and 1.3.0's from 44,897 to
  30,067 (against 1.8.0 alike), carried indels fall by 13,980 and 14,077, and `called_both` rises by
  14,544 and 15,107. Those newly met indels disagree more often than the rest: mismatches rise by
  8,636 and 8,045, so indel discordance rises from 7.1 to 8.0% (1.3.0) and from 8.0 to 8.7% (1.8.0),
  while one-sided falls from 8.8 to 6.0% and 9.2 to 6.3%. In homopolymers of 12 or more bases
  discordance rises from 14.65 to 16.45% and 14.39 to 15.95%. Part of what the first build counted as
  1.1.0 missing an indel is 1.1.0 writing it elsewhere, and often with another genotype.
- **DeepVariant 1.3.0 and 1.8.0 differ in spelling too, at a smaller scale.** One-sided indels fall
  from 2,843 and 3,070 to 866 and 1,066 (0.62 to 0.20% of carried); the alleles both call rise by
  1,854, and discordance from 3.68 to 3.78%. The pair's SNVs hardly change (one-sided 137 to 34).
- **SNVs of the release pairs**: against 1.1.0, carried falls by 333 and 341, one-sided SNVs by 422
  and 443, and 198 and 214 SNVs are now `unreadable`; the rates do not move.

**By confidence band.** SNV discordance, region `all`, scope `carried_pass`, by one side's band with
the other side's bands summed (% of the alleles both called in that band; a cell is left empty where
the band holds under 0.1% of the pair's alleles both called, and no SNV of these pairs states no
confidence):

| Pair | Side | 0–10 | 10–20 | 20–30 | 30–40 | 40–60 | 60–100 | 100+ |
|---|---|---|---|---|---|---|---|---|
| Clair3 1.0.10 / DeepVariant 1.9.0 | Clair3 | 19.76 | 3.48 | 0.12 | 0.01 | 0.00 | 0.01 | |
| Clair3 1.0.10 / DeepVariant 1.9.0 | DeepVariant | 14.74 | 3.36 | 0.24 | 0.01 | 0.00 | | |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | Clair3 | 40.11 | 24.39 | 0.64 | 0.07 | 0.03 | 0.01 | |
| Clair3 1.0.10 / GATK HC 4.6.2.0 | GATK | 47.85 | 36.41 | 13.40 | 1.09 | 0.15 | 0.19 | 0.85 |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | DeepVariant | 46.36 | 26.57 | 2.48 | 0.05 | 0.01 | | |
| DeepVariant 1.9.0 / GATK HC 4.6.2.0 | GATK | 49.10 | 35.93 | 12.90 | 1.11 | 0.13 | 0.11 | 0.64 |
| DeepVariant 1.1.0 / 1.3.0 | 1.1.0 | 15.65 | 4.12 | 1.38 | 0.32 | 0.01 | 0.00 | |
| DeepVariant 1.1.0 / 1.3.0 | 1.3.0 | 11.28 | 1.64 | 0.17 | 0.05 | 0.01 | 0.01 | |
| DeepVariant 1.1.0 / 1.8.0 | 1.1.0 | 18.61 | 5.02 | 1.52 | 0.30 | 0.01 | 0.00 | |
| DeepVariant 1.1.0 / 1.8.0 | 1.8.0 | 10.60 | 2.36 | 0.25 | 0.10 | 0.20 | | |
| DeepVariant 1.3.0 / 1.8.0 | 1.3.0 | 9.61 | 1.09 | 0.07 | 0.01 | 0.00 | 0.00 | |
| DeepVariant 1.3.0 / 1.8.0 | 1.8.0 | 7.73 | 1.49 | 0.17 | 0.05 | 0.07 | | |

- **Below 20 phred the callers disagree on 1 to 49% of their SNVs; from 30 phred up on 1.1% or less,
  mostly under 0.1%.** Every caller here writes a GQ, so every pair is banded, Clair3's included.
- **GATK's top band is not its surest.** GATK in DRAGEN mode puts 9% of its SNVs both called at 100
  phred or more, where it disagrees with DeepVariant on 0.64% and with Clair3 on 0.85% (2,147 of
  333,936 and 2,847 of 334,753), against 0.11 to 0.19% at 60–100 and 0.13 to 0.15% at 40–60. The rise
  is measured here, not diagnosed.
- **The DeepVariant releases use their bands differently**: of the SNVs both called, 1.1.0 puts 69%
  at 40–60 and 19% at 60–100, 1.3.0 puts 42% at 40–60 and 0.2% above, and 1.8.0 puts 73% at 20–40,
  21% at 40–60 and almost none above 60. A band is one release's scale, not one across releases.

## The PGP pairs: CGI against freeBayes

The first build compared these two pairs by `assembly_label` only, with alleles meeting by position
and spelling: `resolve_ref` refused both CGI vcfBeta facets over a chrY md5 that differs from the
store's GRCh37. Since F20 a differing X, Y or MT leaves that contig unread instead, so all four facets
now have REF read from the store's Ensembl GRCh37 collection (every read row matches its source REF
but 3, against IUPAC codes; CGI's chrY, 1,744 and 1,710 rows, stays unread), and `mint_vrs_ids`
mints them (99.96% and over 99.99% of alleles on one pair, 99.96% and 95.1% on the other, whose
freeBayes file has 256,986 rows on `hs37d5` and `NC_007605`, contigs the store does not have). Both pairs then met allele by allele through VRS
ids (`assembly_check = content`, `pairing = vrs`).

Two conventions sat inside these pairs' counts at that build, and neither is method error. Both are
now rules of `concordance`:

- **X outside the PARs is haploid in CGI and diploid in freeBayes** (every genotype there: CGI writes
  one allele, freeBayes two). Every allele both called there was `ploidy_differs` and counted as
  discordant: 72,462 and 64,818 SNVs, 65% and 67% of the pairs' discordant SNVs. Since F22 a haploid
  call against the same allele's homozygote is a `match`, with `ploidy_differs` set apart: of the
  71,715 and 64,382 flagged SNVs in `x_nonpar`, 71,193 and 61,920 now match, and 522 and 2,462 are
  mismatches, a haploid call against a freeBayes het on a male X, which is a real disagreement and
  stays in the rate.
- **A complex freeBayes row can hold the allele under another spelling.** freeBayes writes a
  haplotype as one multi-allelic row: REF `CGCCC`, ALT `TGCCT,TGCCC`, GT `TGCCT/TGCCT`. Its ALT
  `TGCCC` trims to `C>T`, so its VRS id is CGI's `C>T` (GT `T/T`), and the two meet. `concordance`
  counted the copies of the spelling `TGCCC`, zero, although each `TGCCT` haplotype carries the same
  `C>T`. Counted at that build on the VRS-paired mismatches with a script outside the repo: 20,434
  SNVs whose sibling carries exactly the allele's change in the CGA 2.5.0.44 pair, 23,014 with any
  carried sibling; 2,601 and 2,950 in the 2.5.0.24 pair. Since F21 each genotype is read for the
  allele's change in its own row; a sibling of another length over the allele's span cannot be read
  and is `unreadable` (2,315 and 201 SNVs), outside `called_both`.

**Before and after, region `all`, scope `carried_pass`.** "By position" is the first build (label
check, alleles by position and spelling). "By VRS id" is the build with REF read and alleles met by
VRS id, as it stood before F21 and F22; its discordance counted `ploidy_differs`. "Now" is the same
pairing with both rules.

| Methods (a / b) | Class | By position: mismatch % of called both, ploidy left out | By VRS id: discordance % | mismatch % of called both, ploidy left out | Now: called both | discordance % | `ploidy_differs` | `unreadable` | mean `p_same` (gl), before → now |
|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | SNV | 0.49 | 3.40 | 1.19 | 3,278,624 | 0.54 | 72,421 | 2,315 | 0.969 → 0.990 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | SNV | 0.93 | 3.07 | 1.02 | 3,160,432 | 1.02 | 64,818 | 201 | 0.966 → 0.985 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | indel | 44.70 | 7.95 | 5.31 | 405,102 | 4.98 | 10,784 | 5,658 | 0.863 → 0.879 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | indel | 53.12 | 10.41 | 8.05 | 338,839 | 7.92 | 8,063 | 3,391 | 0.865 → 0.882 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | MNV | | 12.37 | 10.77 | 87,290 | 3.86 | 1,328 | 4,701 | 0.661 → 0.706 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | MNV | | 12.23 | 10.92 | 51,206 | 5.76 | 713 | 3,609 | 0.776 → 0.808 |

The 95% Wilson interval of every SNV cell here is within ±0.02 percentage point, of every indel cell
within ±0.1. In the CGA 2.5.0.44 pair the SNV rate falls from 3.40% to 0.54%, near the first build's
0.49%: outside `x_nonpar` its mismatches fall by 21,873, the complex rows (20,434 counted exactly). In the CGA 2.5.0.24
pair it falls from 3.07% to 1.02%, which equals the old mismatch share only by coincidence: outside
`x_nonpar` its mismatches fall by 2,460 (the 2,601 counted), and the 2,462 haploid-against-het X
mismatches, hidden inside `ploidy_differs` before, are now counted. Carried MNVs rise by a quarter
(264,478 to 332,824 in the CGA 2.5.0.44 pair): the atlas now counts an allele as carried where
`concordance` read a copy of it under another spelling, which a VRS-met MNV whose CGI row holds it
only inside another ALT was not before. SNV and indel one-sided rates do not move; the MNV one falls
with the larger carried count (71.0% to 66.9% in the CGA 2.5.0.44 pair).

Pairing by VRS id is what the indels needed. The two spellings of one indel used to count as two
one-sided alleles: in the CGA 2.5.0.44 pair carried indels fall by a third (1.44 million to 0.99
million) and the indels both facets call rise from 11,891 to 405,102; in both pairs one-sided indels
drop from 99% to 53%, 16 to 18% in GIAB's easy regions. `uncalled` grows with it (indels 420 to
43,389 in the CGA 2.5.0.24 pair, as CGI's half-calls now meet freeBayes), so `called_both` is less
than the carried alleles minus the one-sided ones.

**SNV discordance, % of called both** (by VRS id, both rules)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | 0.54 | 0.30 | 1.85 | 3.01 | 2.01 | 1.43 | 11.34 | 3.72 | 2.92 | 0.76 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | 1.02 | 0.62 | 3.35 | 4.21 | 2.27 | 3.47 | 25.00 | 4.84 | 3.27 | 2.42 |

**Indel discordance, % of called both** (by VRS id, both rules)

| Methods (a / b) | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | 4.98 | 1.47 | 7.43 | 6.62 | 5.90 | 3.94 | 14.08 | 12.12 | 9.43 | 4.47 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | 7.92 | 3.10 | 12.27 | 9.80 | 7.76 | 7.81 | 26.52 | 15.01 | 9.55 | 8.67 |

**One-sided, % of carried** (by VRS id)

| Methods (a / b) | Class | all | easy | difficult | low mappability | segdups | homopolymer 7–11 | homopolymer ≥12 | TR 51–200 bp | TR 201–10k bp | GC <25 or >65 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | SNV | 24.93 | 9.81 | 59.49 | 76.53 | 72.47 | 29.74 | 76.20 | 61.86 | 85.13 | 42.75 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | SNV | 31.35 | 14.50 | 61.59 | 74.51 | 66.77 | 52.58 | 89.72 | 65.28 | 81.42 | 48.44 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | indel | 52.87 | 16.13 | 62.27 | 82.66 | 78.69 | 29.93 | 67.09 | 84.76 | 87.74 | 59.76 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | indel | 52.93 | 17.96 | 63.26 | 81.30 | 74.85 | 38.82 | 68.95 | 85.98 | 84.68 | 55.53 |

The 95% Wilson intervals of these `all` and `easy` cells are within ±0.2 percentage point.

Reading them:

- **Easy regions mostly agree across technologies.** CGI and freeBayes disagree on 0.30 to 0.62% of
  SNV genotypes in easy regions, around the same-reads callers' rate (0.05 to 0.40%), and on 1.5 to
  3.1% of indel genotypes. Homopolymers of 12 or more bases are again the worst class, at 11 to 27%.
- **One-sided alleles are the larger difference.** A quarter to a third of carried SNVs, and half the
  indels, are stated by one side only, rising to 75 to 90% in low mappability, segmental duplications
  and long tandem repeats. Most are freeBayes's (`only_b` is 6 times `only_a` for SNVs), and in
  55 to 60% of them the silent side reads as reference (`one_sided_reference`: a reference block over
  the position, or a whole-genome facet's absence rule), a stated disagreement, where the rest is
  silence. In the CGA 2.5.0.24 pair, 262,836 freeBayes alleles (before the carried and PASS
  filters) lie on `hs37d5` and `NC_007605`, which the CGI facet never declared; they are one-sided
  with `absence = unknown` and fall in `all` only.

**By confidence band.** The 2026-10-24 rebuild splits each pair by both sides' own confidence in
the allele (`genotype_evidence`'s confidence of each side's row, docs/SCHEMA.md). SNV discordance,
region `all`, scope `carried_pass`, by one side's band with the other side's pooled (% of the alleles
both called in that band; CGI's few rows below 10 phred and those with no confidence left out):

| Pair | Side | 10–20 | 20–30 | 30–40 | 40–60 | 60–100 | 100+ |
|---|---|---|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 | CGI | 6.17 | 5.35 | 3.29 | 2.11 | 1.18 | 0.37 |
| CGA 2.5.0.44 / freeBayes 0.9.15 | freeBayes | 6.28 | 3.66 | 1.13 | 0.73 | 0.30 | |
| CGA 2.5.0.24 / freeBayes 0.9.20 | CGI | 7.25 | 6.79 | 4.78 | 3.17 | 2.19 | 0.73 |
| CGA 2.5.0.24 / freeBayes 0.9.20 | freeBayes | 9.18 | 2.47 | 0.78 | 0.47 | 0.29 | 0.16 |

freeBayes 0.9.15 has no SNV row at 100 phred or more, so its top band is empty: PGP's merge floors
every GL entry at -10, which caps its confidence at 97 phred, and most of its rows sit at the cap.
The two freeBayes facets' bands are therefore not one scale: in `easy`, against CGI at 100 phred or
more, 0.9.15's rows at 40–60 phred disagree on 0.52% and 0.9.20's on 0.105% (docs/CONSENSUS.md,
"Where the overconfidence at 60 phred came from"). 90 and 91% of CGI's
SNVs both called sit at 100 phred or more, where the pair disagrees on 0.37 and 0.73%, against 0.54
and 1.02% overall; at CGI's rows below 60 phred it disagrees three to eleven times more often. Much
of a method's average rate is its unsure rows. `method_error_rates` fits each pair's cells as one
rate per band of each side (docs/CONSENSUS.md has the rates).

X and Y classes are in the full table only: for a private pair they would state an aggregate about
the person's sex chromosomes.

## The CGI `var` pairs, banded (2026-10-24)

Until schema 0.0.19 a CGI `var` facet stated no confidence: its rows graded `gt_only`, every allele
fell in band `none`, and `method_error_rates` could fit no band for either side of a `var` pair (the
fit needs two graded bands on each). Since `calls.hq` carries CGI's allele scores and
`genotype_evidence` grades a call by the lowest of them (docs/SCHEMA.md, "CGI's allele scores"),
both public `var` pairs were built again with `scripts/method_atlas.py`, alleles met by VRS id:

| Pair | Called both (SNV, `carried_pass`, `all`) | Discordance | One-sided |
|---|---|---|---|
| hu1B2689: CGA 2.5.0.44 `var` / freeBayes 0.9.15 | 3,278,880 | 0.50% | 24.9% |
| huAF3C63: CGA 2.5.0.33 `var` / freeBayes 0.9.15 | 3,141,288 | 1.77% | 29.3% |

hu1B2689's `var` file is the same assembly as its vcfBeta (above: 0.54%, 3,278,624 called both).
Neither file writes a varScoreVAF that differs from its varScoreEAF, and none past u16 (no load warning
in either). The output is `data/output/method_atlas_cgi_var.parquet` and `.tsv` beside the others,
git-ignored.

**SNV discordance by one side's band**, the other side pooled (% of the alleles both called in that
band, region `all`; CGI's few rows below 10 phred are counted in `none`, which holds the 1,059 and
1,301 alleles whose CGI row states no score):

| Pair | Side | 0–10 | 10–20 | 20–30 | 30–40 | 40–60 | 60–100 | 100+ |
|---|---|---|---|---|---|---|---|---|
| hu1B2689 `var` | CGI | | 5.51 | 4.76 | 3.02 | 1.89 | 1.08 | 0.35 |
| hu1B2689 `var` | freeBayes | 8.78 | 6.08 | 3.51 | 1.11 | 0.72 | 0.26 | |
| huAF3C63 `var` | CGI | | 5.14 | 4.11 | 3.06 | 2.28 | 1.80 | 1.66 |
| huAF3C63 `var` | freeBayes | 17.11 | 10.79 | 7.29 | 1.95 | 1.42 | 0.63 | |

**A `var` facet's bands are its vcfBeta's.** On hu1B2689 the `var` pair's CGI bands (5.51, 4.76,
3.02, 1.89, 1.08, 0.35) sit within a fifth of the vcfBeta pair's (6.17, 5.35, 3.29, 2.11, 1.18,
0.37), and freeBayes's side of the two pairs is the same within a tenth. The two facets write one
confidence scale: CGI's vcfBeta GQ is the lowest `var` haplotype score at 42,894 of 42,895 loci of a
chromosome, and the `hq` confidence of the `var` file's 4.02M scored calls (p10/p50/p90 68/374/1,037)
matches the vcfBeta's GL confidence (357 and 1,022 at p50 and p90). That is what makes it sound for
a `var` and a vcfBeta facet of one CGA version to share an atlas key (`caller`, `caller_version`):
the 60+ diagnosis (docs/CONSENSUS.md) found two facets of one freeBayes version that did not share a
scale; these two do. A band also says which genotype a row called: CGI's 10–20 SNV band is its VQLOW
hom-ALT calls, freeBayes 0.9.15's 60–100 band mostly its hets, so the cell of the two is mostly a
genotype disagreement (docs/CONSENSUS.md, "What CGI's lowest band holds"). Two such facets of one person are still one person's pair twice, and never both
go into one atlas (the check below holds hu1B2689's vcfBeta pair out when it adds the `var` one).

**The second person's fit is broken by the even split.** The one-pair fit fixes the two sides' band
rates up to a constant and splits the top cell (CGI 100+ against freeBayes 60–100) evenly. On
huAF3C63 that cell disagrees on 0.64%, so CGI's 100+ band takes 0.32%, and CGA 2.5.0.33's 40–60 and
60–100 bands come out below it, on the Wilson floor (0.004% and 0.002%; 1.80, 1.60, 0.62% at 10 to
40). The disagreement there is mostly freeBayes's: its own side runs about twice hu1B2689's at every
band (16.7, 10.5, 7.0, 1.66, 1.13, 0.32% fitted at 0–10 to 60–100, against 8.1, 5.6, 3.1, 0.76,
0.41, 0.081%). So one freeBayes version does not keep one error per band across people: PGP's merge
for huAF3C63 is a noisier run, by its discordance with CGI (1.77% against 0.50%) as by its bands. A
consensus weighing CGA 2.5.0.33 rows with these bands would take its 40–100 phred rows as nearly
never wrong. The even split is first order on this pair, where docs/CONSENSUS.md found it second order
on hu1B2689's; the rule stands, and the person's options are in docs/CONSENSUS.md.

## More PGP pairs: a third CGI–freeBayes person and CGI's own run-to-run drift (2026-11-08)

Built with `scripts/method_atlas.py`, alleles met by VRS id, every pair `same_person`. The outputs are
`data/output/method_atlas_pgp_timing.*` and `data/output/method_atlas_cgi_drift.*`, git-ignored; each
pair's counts are in the `.pairs/` folder beside them. Columns: SNVs both called (`carried_pass`,
region `all`), their discordance, the one-sided share of all alleles, the easy-region SNV discordance,
and the indel discordance (`all`).

**A third CGI against freeBayes 0.9.15 person.** Only three PGP-HMS people have both a CGI `var` and
PGP's freeBayes+VarScan2 VCF. huF5F075, CGA 2.5.0.48 `var` against freeBayes 0.9.15: 3,126,771 SNVs
both called, discordance 0.36%, one-sided 23.6%; easy SNVs 0.106%, indels 4.10%. That is below both
hu1B2689's 0.50% and huAF3C63's 1.77% above.

**Two CGI sequencing runs of one person.** 16 PGP-HMS people have two CGI `var` files from different
assemblies with different sample ids, so two runs, not one run analysed twice. Fifteen are built here,
along with hu250634 in the timing build; huF06AD0's 2.0.4.14/2.5.0.24 pair is among them. Nine share one
key, 2.5.0.33 against 2.5.0.37; the drift build's pooled row for it holds eight (`pairs` = 8), and
hu250634's is in the timing build:

| Person | Versions | SNVs both called | Discordance | One-sided | Easy SNV | Indel |
|---|---|---|---|---|---|---|
| hu34D5B9 | 2.0.3.5 vs 2.0.4.14 | 3,277,190 | 0.08% | 6.7% | 0.036% | 3.41% |
| huF06AD0 | 2.0.4.14 vs 2.5.0.24 | 3,214,858 | 0.17% | 9.1% | 0.065% | 6.90% |
| hu066C78 | 2.0.4.18 vs 2.5.0.33 | 3,252,442 | 0.15% | 7.9% | 0.054% | 7.12% |
| hu4C20BC | 2.5.0.24 vs 2.5.0.33 | 3,136,293 | 0.20% | 8.7% | 0.084% | 9.10% |
| hu03E3D2 | 2.5.0.33 vs 2.5.0.36 | 3,243,773 | 0.17% | 7.1% | 0.057% | 9.56% |
| hu24A473 | 2.5.0.33 vs 2.5.0.36 | 3,116,775 | 0.63% | 10.9% | 0.485% | 10.19% |
| huD2B804 | 2.5.0.33 vs 2.5.0.36 | 2,833,189 | 3.66% | 18.0% | 3.505% | 13.64% |
| hu1725A5 | 2.5.0.33 vs 2.5.0.37 | 2,973,980 | 1.82% | 14.1% | 1.649% | 11.06% |
| hu250634 | 2.5.0.33 vs 2.5.0.37 | 3,100,963 | 0.91% | 10.8% | 0.753% | 10.26% |
| hu48AC54 | 2.5.0.33 vs 2.5.0.37 | 3,198,281 | 0.39% | 8.9% | 0.260% | 9.97% |
| hu524B5B | 2.5.0.33 vs 2.5.0.37 | 3,124,160 | 0.93% | 11.1% | 0.780% | 10.22% |
| hu87BAAD | 2.5.0.33 vs 2.5.0.37 | 2,997,604 | 1.69% | 13.8% | 1.516% | 10.91% |
| hu9385BA | 2.5.0.33 vs 2.5.0.37 | 3,117,324 | 0.91% | 10.9% | 0.757% | 10.23% |
| huBB058E | 2.5.0.33 vs 2.5.0.37 | 3,222,653 | 0.26% | 8.1% | 0.131% | 9.96% |
| huD58ABC | 2.5.0.33 vs 2.5.0.37 | 2,922,809 | 2.03% | 15.9% | 1.841% | 11.86% |
| huD649F1 | 2.5.0.33 vs 2.5.0.37 | 2,666,437 | 1.27% | 21.6% | 1.094% | 9.50% |

- **One key, a fourteen-fold spread between people.** Across the nine 2.5.0.33/2.5.0.37 pairs, easy-region
  SNV discordance runs from 0.131% to 1.841% (median 0.78%), and the three 2.5.0.33/2.5.0.36 pairs from
  0.057% to 3.505%. A pooled row's interval is sampling error over its alleles, so it says nothing of
  this spread. The measurement went to desk #101 (whether one person's rates may stand for another's);
  pooling stays as today until it is answered.
- **Pairs with a 2.0 or 2.5.0.24 run agree best** (0.036% to 0.084% easy SNV); the 2.5.0.33 against
  2.5.0.36 or 2.5.0.37 pairs run from 0.057% to 3.505%. Whether some runs were weaker sequencing or the
  later versions changed the calling is not settled here: the version and the run change together in
  every pair.
- **Two CGI runs disagree more in easy regions than CGI does with freeBayes on huF5F075** (0.106%) for
  all nine 2.5.0.33/2.5.0.37 people. One-sided shares, 7 to 22% for SNVs, are lower than CGI
  against freeBayes (24 to 29%).
