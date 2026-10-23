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
  person's pair's (0.54%): the same rate on another genome and another release. It stays when
  freeBayes and CGI take their own versions' bands (0.60 and 0.56 below), so it is not the borrowing
  of the next point.
- **60 phred and up, about twice overconfident** (296k alleles, all where CGI and the array agree).
  The prediction there is almost all freeBayes's own rate (the consensus's own error is below 1e-6),
  and freeBayes disagrees on 0.21% against 0.115% predicted. The cause is freeBayes 0.9.15's band
  rates, all borrowed from freeBayes 0.9.20 on another person, whose confidence scale is not 0.9.15's:
  measured in the next section, where freeBayes's own bands, fitted on other chromosomes, close the
  gap (1.02).

## Where the overconfidence at 60 phred came from (2026-10-24)

Two causes were named for the 60+ band and neither was measured: freeBayes 0.9.15 borrowing
freeBayes 0.9.20's band rates, or CGI and the array sharing freeBayes's errors where all three meet,
which the posterior's independence cannot see. A third was the even split of the one pair's top
cell between its two methods. All three were put to the data with `scripts/consensus_check.py`, which
now writes `<out>.calibration_detail.tsv`: the same sums per confidence bin, band of the held-out
facet's row, and which other facets wrote a row at the allele. Every run is over the 22 autosomes of
hu1B2689, against main's atlas of 2026-10-24.

**It is the borrowed bands.** With this person's CGI/freeBayes pair held out, freeBayes 0.9.15 has
no atlas row of its own, and every one of its rates is freeBayes 0.9.20's (`rate_basis = caller`),
measured on huF06AD0's upload. Split by the band of freeBayes's own row, the gap sits in its middle
bands and nowhere else; its top band is calibrated:

| Consensus bin | freeBayes band | Alleles | Predicted | Observed |
|---|---|---|---|---|
| 60 and up | 60–100 | 237,126 | 0.027% | 0.032% |
| 60 and up | 40–60 | 45,195 | 0.093% | 0.496% |
| 60 and up | 30–40 | 10,536 | 0.325% | 0.674% |
| 30 to 40 | 60–100 | 1,780,226 | 0.086% | 0.056% |
| 30 to 40 | 40–60 | 335,821 | 0.156% | 0.514% |
| 30 to 40 | 30–40 | 84,964 | 0.384% | 0.872% |
| 10 to 20 | 0–10 | 2,248 | 29.5% | 7.1% |

freeBayes 0.9.15's rows at 40 to 60 phred err about five times as often as 0.9.20's did, and its
rows below 10 phred about four times less often. The atlas's own cells say the same with no
consensus in between. Region `easy`, CGI at 100 phred or more, discordance with freeBayes by
freeBayes's band:

| Pair | 30–40 | 40–60 | 60–100 | 100+ |
|---|---|---|---|---|
| CGA 2.5.0.44 / freeBayes 0.9.15 (hu1B2689) | 0.72% | 0.52% | 0.064% | (none) |
| CGA 2.5.0.24 / freeBayes 0.9.20 (huF06AD0) | 0.36% | 0.105% | 0.032% | 0.021% |

CGI's own error at 100 phred in hu1B2689 can be no more than 0.064% (the 60–100 cell), so at least
0.45 of the 0.52% is freeBayes 0.9.15's. The two facets write their confidence differently: the PGP
merge's freeBayes 0.9.15 floors every GL entry at -10 (97% of its rows on chromosome 22 hold an entry
at exactly -10.0), so its confidence stops at 97 phred, the median row sits at that cap, and it has
no 100+ band; the 0.9.20 upload's GLs reach -1,450 and its confidence 4,127 phred. A band is a range of
a facet's stated confidence, and the same range means a different error on the two. The version and
the run (pipeline, depth) differ together here, so this measures that the bands of one freeBayes facet
do not predict another's, not that 0.9.15 differs from 0.9.20 as a program.

**freeBayes's own bands close it.** The test, without circularity: this person's pair built again
on the odd autosomes only (`scripts/method_atlas.py --contigs 1 3 … 21`), put in the atlas in place of
the whole-genome pair, and the consensus checked on the even autosomes. freeBayes 0.9.15 then speaks
with its own version's bands (`method`), measured on other chromosomes of the same facet:

| Confidence (phred) | Held out, borrowed bands: alleles | Predicted | Observed | Ratio | Own bands, odd/even split: alleles | Predicted | Observed (95% CI) | Ratio |
|---|---|---|---|---|---|---|---|---|
| below 10 | 29,465 | 47.2% | 48.8% | 1.04 | 12,022 | 49.3% | 54.9% (54.0–55.8) | 1.11 |
| 10 to 20 | 156,810 | 4.57% | 2.43% | 0.53 | 66,936 | 4.64% | 2.78% (2.66–2.91) | 0.60 |
| 20 to 30 | 391,701 | 1.41% | 0.64% | 0.46 | 198,342 | 1.13% | 0.63% (0.60–0.67) | 0.56 |
| 30 to 40 | 2,228,593 | 0.193% | 0.256% | 1.33 | 1,100,069 | 0.283% | 0.254% (0.244–0.263) | 0.90 |
| 40 to 60 | 107,045 | 0.228% | 0.192% | 0.84 | 82,105 | 0.343% | 0.285% (0.251–0.324) | 0.83 |
| 60 and up | 295,818 | 0.115% | 0.210% | 1.82 | 153,444 | 0.209% | 0.212% (0.191–0.237) | 1.02 |
| all | 3,209,432 | 0.98% | 0.85% | 0.87 | 1,612,918 | 0.93% | 0.81% | 0.87 |

Within the 60+ bin the odd/even run predicts 0.494% for freeBayes's 40–60 rows and sees 0.489%, and
0.042% against 0.036% for its 60–100 rows. The fitted bands it used (SNV, region `all`, on the odd
autosomes) against the borrowed ones:

| Band (phred) | freeBayes 0.9.20, borrowed | freeBayes 0.9.15, own (odd autosomes) |
|---|---|---|
| 0–10 | 27.9% | 9.4% |
| 10–20 | 8.7% | 10.0% |
| 20–30 | 2.0% | 5.3% |
| 30–40 | 0.41% | 0.88% |
| 40–60 | 0.18% | 0.45% |
| 60–100 | 0.084% | 0.093% |
| 100+ | 0.049% | (no rows) |

The same run on the whole genome with the pair left in (in-sample, so circular for freeBayes) gives
1.23 at 60 and up and 1.05 at 30 to 40, every freeBayes band holding 2,000 alleles within the 60+
bin inside a fifth of its prediction.

**CGI and the array are not the cause.** Held out in turn, with freeBayes and the array in the
consensus, CGI disagrees with them on 0.002% of 156,486 alleles at 40 to 60 phred and on 0.010% of
245,666 at 60 and up, against 0.073% and 0.029% predicted: CGI's top band is underconfident, not
over. Where freeBayes and the array agree, the error they share is at most that 0.010%, too little
to make a 0.1-point gap; and the gap is there without the array too (30 to 40 phred, where CGI alone
speaks beside freeBayes: 2.19M alleles, 0.184% predicted, 0.257% observed). VarScan2 never wrote a
row where freeBayes was compared (the PGP merge keeps one caller's genotype per site), so it is in
none of these bins.

**The even split is a second-order error.** The one-pair fit puts half the top cell's disagreement
on each method; CGI's near-zero rate where freeBayes and the array agree says CGI's share is smaller,
so the split hands freeBayes's error to CGI. Moving all of it would raise freeBayes 0.9.20's bands by
at most the constant moved, 0.049 points (0.019 in `easy`), taking the 60+ prediction from 0.115% to
about 0.13%, not to 0.21%; it would also lower CGI's top band, which the leave-CGI-out run says is
the right direction. It is not what made the 60+ band twice overconfident, and the rule stays.

**A second freeBayes 0.9.15 facet** (huAF3C63, PGP's merge again, against a CGI `var` file) gives no
band fit: the `var` facet states no confidence, so only freeBayes's side splits by band (10.7, 7.1,
1.8, 1.3 and 0.48% at 10–20 to 60–100 in `easy`, each holding CGI's unbanded share of a pair whose
SNVs disagree on 1.6% in `easy`). It cannot say whether one freeBayes version keeps one scale across
people.

**What changed in the library.** Nothing in the rule: a facet whose version the atlas lacks still
borrows its caller's bands. `consensus` now names such facets in a warning, since a borrowed band was
measured on another facet's confidence scale. Whether to borrow at all is the person's call (below);
`CONSENSUS_VERSION` stays 2, no row changes.

**For the person: bands of a caller measured at another version.**

1. **Borrow, and warn (as built).** Held out: 1.15 times as many disagreements predicted as seen in
   all; freeBayes-like facets run up to five times overconfident in their middle bands and four
   times underconfident in their lowest. Costs nothing where the atlas holds the facet's own version.
2. **Never borrow bands across versions.** A facet without its own version's bands speaks with its
   caller's average, as before bands. On this held-out check both CGI and freeBayes lose their bands
   and the consensus falls back to the pooled calibration (51,225 predicted against 27,243 seen, 1.9
   times). It gives up CGI, whose bands carry over across its versions within a third
   (4.19/4.26/2.62/1.25/0.60% against 5.39/4.74/2.71/1.61/0.73% at 10 to 100 phred) except its top
   band (0.049 against 0.098%).
3. **Borrow only between versions shown to agree.** Bands of another version are used where the atlas
   holds two versions of the caller and their bands agree within their intervals; otherwise the
   average. It needs two measured versions per caller and a tolerance, and on this check it equals
   option 2.

Recommended: 1, and measure each facet's own version wherever the person has a partner facet for it,
which is the normal use (an atlas with the person's own pairs: in-sample 1.14 overall). The odd/even
split above is the honest way to check that use, and it is calibrated at 30 phred and up.
