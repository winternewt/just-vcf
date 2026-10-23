# The first consensus check (RM280)

Run 2026-10-23 with `scripts/consensus_check.py` on one public PGP-HMS participant, hu1B2689, over
the 22 autosomes on GRCh37: four facets, `consensus` over all of them and once without freeBayes,
with GIAB v3.6 regions. Everything here is an aggregate over the whole autosomal genome; no position
is named. The rules are in docs/SCHEMA.md ("A consensus posterior").

## The facets and their rates

| Facet | Method | Absence | Rate basis |
|---|---|---|---|
| Complete Genomics vcfBeta | CGAPipeline 2.5.0.44, REF not read (its chrY md5 differs from the store's GRCh37) | `reference` | `caller`: CGAPipeline 2.5.0.24 |
| freeBayes, PGP's merge | freeBayes 0.9.15 | `weak_reference` | `caller`: freeBayes 0.9.20 |
| VarScan2, the same merge (kin of freeBayes: the same reads) | VarScan2 | `unknown` | `default` |
| 23andMe export | array, REF read | `unknown` | `default` |

The atlas (docs/METHOD_ATLAS.md) holds this person's CGI/freeBayes pair. It was **held out**: the
rates come from the other PGP pair (huF06AD0, CGA 2.5.0.24 against freeBayes 0.9.20) and the private
same-reads set, so the check measures them on a genome they were not drawn from. With one pair behind
each PGP method there is no triangle, and CGI and freeBayes each get half of that pair's genotype
mismatch: 0.47% for an SNV overall, 0.28% in GIAB's easy regions, 1.9% in low mappability.
VarScan2 and the array were never measured and get the worst rate any method has in each class
(2.2% overall, GATK on the same reads; 0.32% easy; 17% low mappability). The silence rates are
directional: freeBayes left 6.2% of CGI's carried SNVs one-sided, CGI 28.7% of freeBayes's. The PGP
indel and MNV rates (27% and 2.6%, 98% silence) are representation, not method error
(docs/METHOD_ATLAS.md), so they mostly mute those classes; the headline here is SNVs.

**The atlas's PGP SNV discordance is mostly X ploidy.** Of the 88,096 discordant SNVs of this
person's pair (2.71%), 72,068 are `ploidy_differs`, 71,362 of them on X outside the PARs: CGI writes
a male's X haploid, freeBayes diploid. Genotype mismatch alone is 0.50% (16,028), and 0.94% for the
other pair, whose 2.98% is 69% ploidy. `method_error_rates` leaves `ploidy_differs` out for that
reason. The first run of this check took the atlas's `discordant` as it is and was about seven times
underconfident; that rate was the cause, not the combination.

The consensus has 4.44M SNV, 0.35M MNV and 1.43M indel rows. Of the SNV rows, 86.2% have no
contradicting facet, 13.6% one and 0.2% two. At most 13% of those contradictions are two stated
genotypes of the allele disagreeing (78k mismatches with the sources, below); the rest come from a
facet's silence read as reference or its row of another allele. 14.5% of SNV rows had a facet
speaking with a default rate.

## Agreement with each source

The share of SNV alleles either side carries, both called, where the consensus's most probable
genotype is the source's (`concordance(consensus, source)`, `match / (match + mismatch)`):

| Source | Both called | Agreement |
|---|---|---|
| CGI vcfBeta | 3.35M | 99.88% |
| 23andMe | 448k | 99.85% |
| VarScan2 | 196k | 98.58% |
| freeBayes | 4.10M | 98.27% |

Where CGI and freeBayes disagree, their rates are equal, so the consensus settles it by their own
likelihoods and by the facets silent or present there; CGI's GL is the more decisive of the two, and
the consensus sits nearer CGI. That is a reading of the evidence as given, not a judgement that CGI
is right; the per-facet calls stay in the parents.

## Calibration against the held-out facet

The consensus of CGI, VarScan2 and the array, without freeBayes, against freeBayes on the 3.21M SNV
alleles both called and either carries. Predicted is the probability freeBayes disagrees with the
consensus's top genotype, under the consensus's posterior and freeBayes's own overall rate (0.47%);
observed is how often it did, with the 95% Wilson interval:

| Confidence (phred) | Alleles | Predicted | Observed (95% CI) |
|---|---|---|---|
| below 10 | 39,071 | 37.4% | 37.4% (36.9–37.9) |
| 10 to 20 | 458,783 | 2.26% | 1.16% (1.13–1.19) |
| 20 to 30 | 2,279,595 | 0.81% | 0.28% (0.28–0.29) |
| 30 to 40 | 27,863 | 0.55% | 0.16% (0.12–0.21) |
| 40 to 60 | 404,120 | 0.47% | 0.21% (0.19–0.22) |

Where the consensus is unsure it is calibrated. Above 10 phred it is **underconfident by two to three
times**: freeBayes disagrees with it on 0.2 to 1.2% of the alleles, where the posterior predicts 0.5
to 2.3%. The 39k alleles below 10 phred hold 54% of the disagreements (14,612 of 27,247), and a
single rate per method and region class spreads its error over every allele of the class, so it
caps a facet's confident likelihoods at the class average. Two smaller causes sit inside the
prediction too: freeBayes's own rate is taken overall, not per class, and the prediction counts
freeBayes's error as the half-pair attribution gives it. Rates conditioned on the facet's own
confidence (an atlas row per GQ or GL band) are the obvious next measure; none is built.
