# The consensus checks (RM280, and its calibration)

Run with `scripts/consensus_check.py` on one public PGP-HMS participant, hu1B2689, over the 22
autosomes on GRCh37: four facets, `consensus` over all of them and once without freeBayes, with GIAB
v3.6 regions. Everything here is an aggregate over the whole autosomal genome; no position is named.
The rules are in docs/SCHEMA.md ("A consensus posterior"). The first check (2026-10-23, RM280) ran on
the atlas as it stood before F21 and F22; the second (2026-10-24) on the atlas rebuilt with both
rules, once with each method's average rate ("pooled") and once with its rates by confidence band
("banded", run at `CONSENSUS_VERSION` 2).

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
  of the next point. Measured since (2026-10-24, "Where the underconfidence at 10 to 30 phred came
  from", below): it is mostly the region rule, an allele weighted by the highest rate among the
  nested GIAB classes holding it.
- **60 phred and up, about twice overconfident** (296k alleles, all where CGI and the array agree).
  The prediction there is almost all freeBayes's own rate (the consensus's own error is below 1e-6),
  and freeBayes disagrees on 0.21% against 0.115% predicted. The cause is freeBayes 0.9.15's band
  rates, all borrowed from freeBayes 0.9.20 on another person, whose confidence scale is not 0.9.15's:
  measured in the next section, where freeBayes's own bands, fitted on other chromosomes, close the
  gap (1.02).

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
people. Since schema 0.0.19 the `var` facet is banded and the pair fits (docs/METHOD_ATLAS.md, "The CGI
`var` pairs"): it does not keep one; huAF3C63's freeBayes bands run about twice hu1B2689's.

**What changed in the library.** Nothing in the rule: a facet whose version the atlas lacks still
borrows its caller's bands. `consensus` now names such facets in a warning, since a borrowed band was
measured on another facet's confidence scale. Whether to borrow at all is the person's call (below).
The warning changed no row, so it left `CONSENSUS_VERSION` at the 2 it was then (later changes have
moved it since; `consensus.py` holds the current value).

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

## The CGI `var` facet banded (2026-10-24)

Schema 0.0.19 reads CGI's allele scores into `calls.hq`, and `genotype_evidence` grades a `var` row
by the lowest of them (docs/SCHEMA.md, "CGI's allele scores"), so a `var` facet speaks with its
confidence band where it spoke with its method's average. The check, as the odd/even split above:
hu1B2689's facets with the **`var` file** in place of the vcfBeta (CGI `var`, freeBayes, VarScan2,
23andMe; freeBayes held out), the atlas being main's without hu1B2689's CGI/freeBayes pair plus the
`var` pair built on the odd autosomes, the consensus checked on the even ones. "Before" is the code
before 0.0.19 (a copy of it, run on the same files): the `var` pair then fitted no band on either
side, so CGI spoke with its average and freeBayes 0.9.15 borrowed 0.9.20's bands. Both runs also
held huAF3C63's `var` pair whole-genome (unbanded before, banded after), which gives freeBayes 0.9.15 a
second partner; a third run leaves it out.

| Confidence (phred) | Before: alleles | Predicted | Observed | Ratio | After: alleles | Predicted | Observed | Ratio | After, without huAF3C63: ratio |
|---|---|---|---|---|---|---|---|---|---|
| below 10 | 11,698 | 64.8% | 96.9% | 1.50 | 11,593 | 65.0% | 97.1% | 1.49 | 1.49 |
| 10 to 20 | 504,960 | 1.89% | 0.79% | 0.42 | 66,802 | 4.30% | 3.26% | 0.76 | 0.82 |
| 20 to 30 | 882,447 | 1.07% | 0.30% | 0.28 | 198,706 | 1.25% | 0.77% | 0.62 | 0.76 |
| 30 to 40 | 58,361 | 0.257% | 0.192% | 0.75 | 1,101,077 | 0.448% | 0.255% | 0.57 | 0.92 |
| 40 to 60 | 156,093 | 0.131% | 0.211% | 1.62 | 231,135 | 0.434% | 0.242% | 0.56 | 0.85 |
| 60 and up | | | | | 4,246 | 0.575% | 0.141% | 0.25 | 1.02 |
| 10 and up | 1,601,861 | 1.21% | 0.44% | 0.37 | 1,601,966 | 0.71% | 0.44% | 0.63 | 0.85 |

- **Banded, the `var` facet calibrates as the vcfBeta does.** Without huAF3C63 the `var` run is
  0.82, 0.76, 0.92, 0.85, 1.02 from 10 phred up, against the vcfBeta's 0.60, 0.56, 0.90, 0.83, 1.02
  (above): the same shape, the two facets sharing one confidence scale (docs/METHOD_ATLAS.md). Before,
  1.39M alleles sat at 10 to 30 phred, mostly CGI rows weighted by CGI's average rate, predicted at
  two and a half to three and a half times the disagreement seen.
- **huAF3C63 makes it worse, through freeBayes.** Its pair gives freeBayes 0.9.15 a second partner
  whose bands are about twice hu1B2689's (docs/METHOD_ATLAS.md), and a method's band rate is the mean
  over its partners, so freeBayes's predicted error rises with it, and at 30 phred and up, where the
  prediction is mostly freeBayes's own, the run sees 0.57 and 0.56 of it. A
  version's bands measured on one person's facet do not carry to another's facet of that version.
- **Below 10 phred is no calibration signal.** Those 11.6k alleles are where no facet in the
  consensus spoke (`speaking` 0: the `var` file has a row at the position but over another span, or a
  half-call), so the posterior is flat, its top genotype the tie's zero copies, and freeBayes,
  carrying the allele, disagrees by construction. The rows "10 and up" leave them out.

## Where the underconfidence at 10 to 30 phred came from (2026-10-24)

With each facet's own bands the consensus still saw about 0.6 of the disagreements it predicted at
10 to 30 phred (0.60, 0.56 above). `scripts/consensus_check.py` now splits its detail by GIAB region
(`region`: `easy`, `difficult`, or neither) and runs two variants of the region rule (`--no-regions`,
every facet weighted by its `all` rates; `--classes easy --classes difficult`, by GIAB's partition of
the genome alone), all on the vcfBeta odd/even split above:

| Rule | 10–20 | 20–30 | 30–40 | 40–60 | 60+ | 10+ | `easy` | `difficult` |
|---|---|---|---|---|---|---|---|---|
| highest rate among all GIAB classes holding the allele (the rule) | 0.60 | 0.56 | 0.90 | 0.83 | 1.02 | 0.71 | 0.99 | 0.80 |
| `all` rates only (`--no-regions`) | 0.85 | 0.80 | 0.87 | 0.84 | | 0.86 | 0.80 | 1.17 |
| `easy` / `difficult` only (`--classes`) | 0.80 | 0.71 | 0.98 | 0.89 | 1.00 | 0.86 | 1.06 | 0.91 |
| the same, the `var` facet banded | 1.00 | 1.04 | 1.00 | 0.91 | 1.00 | 1.01 | | |

(observed over predicted; the region columns are over every bin, the one below 10 phred included.)

**It is the region rule.** 224k of the 265k alleles at 10 to 30 phred lie in GIAB's `difficult`
regions, almost all CGI rows speaking alone, and there the rule weighs each by the highest rate among
the classes holding it. The classes nest and overlap: `difficult` is the union of the rest, and an
allele in a mild class (a 7–11 bp homopolymer, a short tandem repeat) also sits in `difficult`, whose rate is the
union's average, raised by low mappability and segmental duplications. The highest of overlapping
averages overstates the error of most alleles it is applied to: `difficult` alleles were predicted
1.25 times the disagreement seen, the 10 to 30 phred ones almost twice. `easy`, which nests nothing,
was right (0.99). Taking `all` alone swings the other way, `easy` over-predicted and `difficult`
under (0.80, 1.17), and agrees overall only by cancelling. GIAB's partition, each allele in exactly
one of `easy` and `difficult`, makes `easy` calibrated in every bin (1.00, 1.00, 0.99, 0.89, 1.00 from
10 phred up) and with the banded `var` facet every bin from 10 phred up lies within a tenth of 1.

**What is left** is the vcfBeta run's `difficult` 10 to 30 phred alleles under the partition (0.77,
0.67): mostly CGI rows of 100+ phred alone, at the band's `difficult` rate, which the even split of the
pair's top cell sets at 0.39% (the `var` pair's, where the partition calibrates, at 0.30%). The cells
also say the additive fit does not hold for CGI's lowest band: CGI rows at 10 to 20 phred disagree
with freeBayes on 0.6 to 1.4% where freeBayes is at 30 to 60 phred and on 11 to 12% where freeBayes
is at its 60–100 cap, so no one rate fits the band; it is a second, smaller cause, and not chased here.

**For the person: which regions weigh an allele.** The library rule (`consensus(regions=)`: the
highest rate among the classes of the `RegionSet` holding the allele) is unchanged; the choice is
which `RegionSet` a caller passes, and what the docs and checks recommend.

1. **Every GIAB class (as run until now).** Conservative: no allele is weighted below any class it
   sits in. Costs calibration in difficult regions, 0.71 overall from 10 phred up, about 0.6 at 10 to
   30 phred.
2. **GIAB's partition, `easy` and `difficult` (recommended).** `giab_stratifications(..., classes=
   ["easy", "difficult"])`. Calibrated in `easy` in every bin, 0.86 (vcfBeta) and 1.01 (`var`) from
   10 phred up. Gives up the finer classes' rates, so a low-mappability allele is weighted as an
   average difficult one; within `difficult` the remaining gap is the top cell's even split.
3. **A finer partition built from the classes** (each allele in one class: low mappability first,
   then segmental duplications, and so on, by a stated precedence). Keeps the hard classes' rates
   without the overlap; not built, and each class then holds fewer alleles per band
   (`MIN_BAND_ALLELES`).
4. **`all` only.** Calibrated overall by cancellation; wrong in both directions by region.
