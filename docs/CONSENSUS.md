# The first consensus check (RM280)

Run 2026-10-23 with `scripts/consensus_check.py` on one public PGP-HMS participant, hu1B2689, over
the 22 autosomes on GRCh37: four facets, `consensus` over all of them and once without freeBayes.
Everything here is an aggregate over the whole autosomal genome; no position is named. The rules are in
docs/SCHEMA.md ("A consensus posterior").

## The facets and their rates

| Facet | Method | Rate basis |
|---|---|---|
| Complete Genomics vcfBeta | CGAPipeline 2.5.0.44, REF not read (its chrY md5 differs from the store's GRCh37) | `caller`: CGAPipeline 2.5.0.24 |
| freeBayes, PGP's merge | freeBayes 0.9.15 | `caller`: freeBayes 0.9.20 |
| VarScan2, the same merge (kin of freeBayes: the same reads) | VarScan2 | `default` |
| 23andMe export | array, REF read | `default` |

The atlas (docs/METHOD_ATLAS.md) has this person's CGI/freeBayes pair. It was **held out**: the rates
come from the other PGP pair (huF06AD0, CGA 2.5.0.24 against freeBayes 0.9.20) and the private
same-reads set, so the check measures the rates on a genome they were not drawn from. With one pair
behind each PGP method there is no triangle, so CGI and freeBayes each get half of that pair's
discordance: 1.49% for an SNV genotype overall, 1.24% in GIAB's easy regions, 4.1% in low
mappability. VarScan2 and the array were never measured and get the worst rate any method has per
class (2.2% overall, from GATK on the same reads). The silence rates are directional: freeBayes left
6.2% of CGI's carried SNVs one-sided, CGI 28.7% of freeBayes's. The PGP indel and MNV rates (27% and
2.6% genotype, 98% silence) are representation, not method error (docs/METHOD_ATLAS.md), so they mostly
mute those classes rather than weigh them; the headline here is SNVs.

The consensus has 4.44M SNV, 0.35M MNV and 1.43M indel rows. Of the SNV rows, 86.2% have no
contradicting facet, 13.6% one and 0.2% two; most of the one-contradiction rows are a call one WGS
facet made where the other's silence reads as reference. 14.5% of SNV rows had a facet speaking with a
default rate.

## Agreement with each source

The share of SNV alleles either side carries, both called, where the consensus's most probable
genotype is the source's (`concordance(consensus, source)`, `match / (match + mismatch)`):

| Source | Both called | Agreement |
|---|---|---|
| CGI vcfBeta | 3.35M | 99.88% |
| 23andMe | 448k | 99.86% |
| VarScan2 | 196k | 98.57% |
| freeBayes | 4.10M | 98.26% |

CGI and freeBayes differ on 2.7% of the SNV genotypes they both call (the atlas, for this
pair). With equal rates the consensus settles those by their likelihoods and the facets silent there,
and CGI's GL is the more decisive of the two, so the consensus sits nearer CGI. That is a reading of
the evidence as given, not a judgement that CGI is right; the per-facet calls stay in the parents.

## Calibration against the held-out facet

The consensus of CGI, VarScan2 and the array, without freeBayes, against freeBayes on the SNV alleles
both called and either carries (3.21M). Predicted is the probability freeBayes disagrees with the
consensus's top genotype under the consensus's posterior and freeBayes's own rate (1.49%); observed is
how often it did, with the 95% Wilson interval:

| Confidence (phred) | Alleles | Predicted | Observed (95% CI) |
|---|---|---|---|
| below 10 | 39,758 | 38.0% | 37.1% (36.6–37.5) |
| 10 to 20 | 2,735,070 | 3.09% | 0.43% (0.42–0.44) |
| 20 to 30 | 30,155 | 1.64% | 0.18% (0.14–0.23) |
| 30 to 40 | 94,858 | 1.50% | 0.20% (0.17–0.23) |
| 40 to 60 | 309,591 | 1.50% | 0.21% (0.19–0.22) |

Where the consensus is unsure it is calibrated. Everywhere else it is **underconfident by about seven
times**: freeBayes disagrees with it on 0.2 to 0.4% of the alleles, where one method's rate alone
predicts 1.5%. The atlas's 2.7 to 3.0% CGI/freeBayes discordance is concentrated in a small set of
sites (the 40k below 10 phred hold 54% of the disagreements), and a single per-method rate spreads it
evenly: it caps a facet's confident likelihoods at its average error, so 85% of the alleles, those
one WGS facet carries alone or both agree on, sit at 10 to 20 phred whatever their own GL says. The
region classes narrow this only a little (easy 1.24% against 1.49% overall). Rates conditioned on the
facet's own confidence (an atlas row per GQ or GL band), or a rate that tempers the vector less, are
the obvious next measures; neither is built.
