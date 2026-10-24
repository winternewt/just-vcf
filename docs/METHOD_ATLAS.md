# The first method-error atlas (RM279)

Built 2026-10-23 with `scripts/method_atlas.py` over eight same-person facet pairs, on GIAB v3.6
stratifications. The two PGP-HMS pairs were rebuilt the same day once `resolve_ref` read CGI's REF
(F20), with VRS ids minted on both sides, and again once `concordance` read complex rows by the
allele's change (F21) and a haploid call against the same allele's homozygote as a match (F22). The
full table (768 rows: every region class, variant class and scope) is `data/output/method_atlas.parquet`
and `.tsv`, git-ignored, with `method_atlas.pairs.tsv` listing each pair's methods and checks; the
last CGI build is `data/output/method_atlas_cgi.*` beside it (192 pooled rows; rebuilt 2026-10-24 with the cells of both sides' confidence bands beside them, 10,552 rows, the pooled ones unchanged). Everything here is an aggregate over a whole genome: no position, sample
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
unfiltered. Discordant is `mismatch` over the alleles both sides called (`match + mismatch`); a
haploid call against the same allele's diploid homozygote is a match, its ploidy counted apart in
`ploidy_differs` (F22). One-sided is over all carried alleles. The 95% Wilson intervals are within
±0.1 percentage point for every `all` and `easy` cell of the private pairs and are in the full table.
The private pairs have no allele whose two genotypes differ in ploidy (`ploidy_differs` is 0 in every
row of their build), so the F22 rule leaves their rates as they are; they met by position, and the
F21 reading changes only an allele whose genotype holds another ALT of one row that carries its
change, which a caller writing one allele per site barely has. They were not rebuilt. The PGP pairs
have their own tables below.

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
0.980 to 0.984 for the DeepVariant releases, 0.985 and 0.990 for the PGP pairs. Clair3's rows carry
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
scale; these two do. Two such facets of one person are still one person's pair twice, and never both
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
