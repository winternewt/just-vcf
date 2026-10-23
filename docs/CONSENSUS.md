# The consensus checks (RM280, and its calibration)

Run with `scripts/consensus_check.py` on one public PGP-HMS participant, hu1B2689, over the 22
autosomes on GRCh37: four facets, `consensus` over all of them and once without freeBayes, with GIAB
v3.6 regions. Everything here is an aggregate over the whole autosomal genome; no position is named.
The rules are in docs/SCHEMA.md ("A consensus posterior"). The first check (2026-10-23, RM280) ran on
the atlas as it stood before F21 and F22; the second (2026-10-24) on the atlas rebuilt with both
rules, once with each method's average rate ("pooled") and once with its rates by confidence band
("banded", `CONSENSUS_VERSION` 2).

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
mismatch: 0.51% for an SNV overall, 0.31% in GIAB's easy regions, 2.1% in low mappability (0.47, 0.28
and 1.9% in the first check, before F22 counted a haploid call against a het on a male X). VarScan2
and the array were never measured and get the worst rate any method has in each class (2.2% overall,
GATK on the same reads; 0.32% easy; 17% low mappability). The silence rates are directional:
freeBayes left 6.1% of CGI's carried SNVs one-sided, CGI 28.2% of freeBayes's. The PGP indel and MNV
rates are representation as much as method error (docs/METHOD_ATLAS.md), so the headline here is
SNVs.

**By confidence band** (SNV, region `all`; fitted from the other pair's cells of both sides' bands,
docs/SCHEMA.md "Rates by confidence band"):

| Band (phred) | CGI | freeBayes |
|---|---|---|
| 0–10 | (366 alleles, too few) | 27.9% |
| 10–20 | 4.2% | 8.7% |
| 20–30 | 4.3% | 2.0% |
| 30–40 | 2.6% | 0.41% |
| 40–60 | 1.3% | 0.18% |
| 60–100 | 0.60% | 0.084% |
| 100+ | 0.049% | 0.049% |
| average | 0.51% | 0.51% |

The averages hide a factor of ten either way: 89% of CGI's SNVs both called sit at 100 phred or more
and err about ten times less often than CGI's average, and its rows below 60 phred err three to eight
times more. freeBayes's rows are spread wider, and its own confidence tracks its error more closely.
Neither method's GL is calibrated as a probability: a CGI row at 20 to 30 phred states at most 1%
error and is wrong on about 4%.

## Agreement with each source

The share of SNV alleles either side carries, both called, where the consensus's most probable
genotype is the source's (`concordance(consensus, source)`, `match / (match + mismatch)`):

| Source | Both called | First check | Pooled | Banded |
|---|---|---|---|---|
| CGI vcfBeta | 3.35M | 99.88% | 99.88% | 99.79% |
| 23andMe | 448k | 99.85% | 99.85% | 99.85% |
| VarScan2 | 196k | 98.58% | 98.58% | 98.58% |
| freeBayes | 4.10M | 98.27% | 98.27% | 98.33% |

With pooled rates CGI and freeBayes weigh alike, so where they disagreed CGI's more decisive GL
settled it and the consensus sat nearer CGI. With banded rates a freeBayes row at 40 phred outweighs a
CGI row at 30: CGI's mismatches with the consensus rise by 3,081 and freeBayes's fall by 2,510. That is a reading of
the evidence as measured, not a judgement that either is right; the per-facet calls stay in the
parents. The consensus has 4.44M SNV, 0.35M MNV and 1.43M indel rows. Of the SNV rows, 86.2% have no
contradicting facet, 13.6% one and 0.2% two; 14.5% had a facet speaking with a default rate.

## Calibration against the held-out facet

The consensus of CGI, VarScan2 and the array, without freeBayes, against freeBayes on the SNV
alleles both called and either carries (3.21M in the first check, 3.21M now). Predicted is the
probability freeBayes disagrees with the consensus's top genotype, under the consensus's posterior,
with freeBayes erring at its own rate at that allele by the rule the consensus weighs a facet by (the
highest among the region classes holding it, of its row's confidence band where the atlas has one);
observed is how often it did, with the 95% Wilson interval. The first check took freeBayes's one
overall rate instead.

**First check** (RM280, the atlas before F21 and F22, one overall rate for freeBayes):

| Confidence (phred) | Alleles | Predicted | Observed (95% CI) |
|---|---|---|---|
| below 10 | 39,071 | 37.4% | 37.4% (36.9–37.9) |
| 10 to 20 | 458,783 | 2.26% | 1.16% (1.13–1.19) |
| 20 to 30 | 2,279,595 | 0.81% | 0.28% (0.28–0.29) |
| 30 to 40 | 27,863 | 0.55% | 0.16% (0.12–0.21) |
| 40 to 60 | 404,120 | 0.47% | 0.21% (0.19–0.22) |

**Pooled rates** (the atlas with F21 and F22, each method's average):

| Confidence (phred) | Alleles | Predicted | Observed (95% CI) |
|---|---|---|---|
| below 10 | 39,085 | 40.8% | 37.4% (36.9–37.9) |
| 10 to 20 | 459,012 | 3.64% | 1.16% (1.13–1.19) |
| 20 to 30 | 2,284,581 | 0.73% | 0.28% (0.28–0.29) |
| 30 to 40 | 22,654 | 1.74% | 0.14% (0.10–0.20) |
| 40 to 60 | 404,100 | 0.36% | 0.21% (0.19–0.22) |
| all | 3,209,432 | 1.60% | 0.85% |

**Banded rates** (the same atlas, rates by confidence band):

| Confidence (phred) | Alleles | Predicted | Observed (95% CI) | Ratio |
|---|---|---|---|---|
| below 10 | 29,465 | 47.2% | 48.8% (48.2–49.4) | 1.03 |
| 10 to 20 | 156,810 | 4.57% | 2.43% (2.36–2.51) | 0.53 |
| 20 to 30 | 391,701 | 1.41% | 0.64% (0.62–0.67) | 0.46 |
| 30 to 40 | 2,228,593 | 0.193% | 0.256% (0.250–0.263) | 1.33 |
| 40 to 60 | 107,045 | 0.228% | 0.192% (0.168–0.221) | 0.85 |
| 60 and up | 295,818 | 0.115% | 0.210% (0.194–0.227) | 1.82 |
| all | 3,209,432 | 0.98% | 0.85% | 0.87 |

Over all alleles the pooled consensus predicts 51,225 disagreements and sees 27,243 (1.9 times too
many predicted); the banded one predicts 31,462 against 27,245 (1.15 times). The 2.23M alleles CGI
alone carries at 100 phred or more moved from the 20–30 band, where the pooled rate capped them, to
30–40, where predicted and observed are within a third of each other. Where the consensus is unsure
it stays calibrated. Two bands are still off, in opposite directions:

- **10 to 30 phred, about twice underconfident** (548k alleles). These are mostly CGI rows of middle
  confidence, alone or against an array. Their band rates (4.2 and 4.3%) come from one other
  person's CGA 2.5.0.24 against freeBayes 0.9.20, whose SNV discordance (1.02%) is nearly twice this
  person's pair's (0.54%): the same rate on another genome and another release.
- **60 phred and up, about twice overconfident** (296k alleles, most where CGI and the array agree).
  The prediction there is almost all freeBayes's own rate (the consensus's own error is below 1e-6),
  0.05 to 0.18% by its band, and freeBayes disagrees on 0.21%. The array and CGI are not
  independent of freeBayes's errors where all three meet, which the independence the posterior
  assumes cannot see, and the observed rate barely moves from 40 to over 60 phred (0.19 to 0.21%):
  a floor of disagreement no confidence band of either side explains.

Both rest on one pair per PGP method. A second CGI/freeBayes pair, from another participant, would
let the band rates be pooled over people and show how far they move from one genome to the next.

## Options left to the person

The checks above ran with every default. Three choices that move the consensus are the person's
(desk #94–#96), each built as an option whose default is what ran here, so an answer is a one-line
default change; the rules are in docs/SCHEMA.md.

- **#94, how a gVCF reference block speaks** (`block_grade=`): by its likelihoods (the default), by its
  GQ, or by the lower of the two. The checks have not been run under the other two.
- **#95, GQ as the confidence of a reference call** (`gq_on_reference=`): off by default, as VCF
  defines GQ (conditional on the site being variant); on, a row that lists the allele and calls
  hom-ref with a GQ and no likelihoods is graded by that GQ. It reaches calls, not blocks, which
  `block_grade` grades.
- **#96, a population prior** (`prior=AlleleFrequencyPrior(...)`): flat by default, as here. With a
  caller's allele-frequency table, an allele it holds gets Hardy-Weinberg genotype priors from its
  frequency and the support row names the table. Which catalogue, if any, is the person's choice;
  none is downloaded or packaged.
