# Relatedness between facets

`just_vcf.relatedness` (RM274) measures how alike two facets' genotypes are and names the relationship
band that puts them in. It is what `Person.of` (RM273) runs to check that facets a caller files under
one subject are one person's. This page says how it works, what it was measured on, and what it cannot
do. Nothing here identifies anyone: it compares genotypes, and a band is a statement about two facets,
not about who their people are.

## How

- **Sites.** A panel of autosomal, strand-unambiguous SNVs on one assembly. The package carries one per
  build, GRCh37 and GRCh38, cut from somalier's published site VCFs (17,766 common gnomAD SNVs each) to
  15,289: the 382 sites off the autosomes and the 2,095 A/T and C/G sites dropped, every REF checked
  against the refget store (docs/PROVENANCE.md). Any other `SitePanel` can be passed; it is refused if
  it breaks those rules.
- **Genotypes.** Each facet is read at the panel's sites by `just_vcf.site_genotypes`, the reading
  population inference uses: a row's literal alleles against the panel's two, so an array export (no
  REF) reads like a VCF; a no-call, a half call, a haploid call or a row of a pre-split site is unknown,
  never 0; so is a row whose source REF is only `N` (`ref_check = unstated`, a copy-number caller's
  REF) with a sequence ALT, whose letters were spelled against the `N` (`N` → `G`, `1/1`, names no
  base beside `G`): the site is unknown, `ref_unstated`, and still a row there, so the absence rule
  does not reach it; an absent site is homozygous reference only where the facet's `absence_means` says so
  (`weak_reference`, inside a `ref` block for `blocks`, outside every block for `reference`). An
  array's silence stays unknown.
- **Kinship.** KING-robust (Manichaikul et al. 2010), the between-family form PLINK 2's `--make-king`
  uses, over the sites both facets have a genotype at. No allele frequencies, so it holds across
  ancestries better than an estimator that needs them. One person is 0.5, first degree 0.25, second
  0.125, third 0.0625.
- **Bands and uncertainty.** KING's floors: 0.354 one person, 0.177 first degree, 0.0884 second, 0.0442
  third, below that unrelated. A delete-one-block jackknife over 50 contiguous blocks of panel sites
  gives a standard error; a band is named only when kinship ± 2 SE lies inside it, otherwise the
  answer is `unknown`, reason `between_bands`. `same_person` is True above the one-person floor, False
  below it, None across it.
- **First degree, split by IBS0.** Kinship puts a parent and child and full siblings at the same 0.25.
  What tells them apart is IBS0, the share of shared sites where the two are homozygous for opposite
  alleles: a parent and child share an allele at every site, so only genotyping error gives them any,
  while full siblings share no allele identical by descent over about a quarter of the genome. The
  rule is KING's (Manichaikul et al. 2010, Bioinformatics 26:2867, section 2.3: "IBS between a
  parent–offspring pair is always 1 or 2 at any SNP in the absence of genotyping errors"). Neither the
  paper nor KING's manual prints a cutoff; the numbers are the fixed defaults in KING 2.3.2's source
  (`KingCore.cpp`, the first-degree split `--build` makes before it refits the cutoff to the sample;
  `relationship.cpp` falls back to the same 0.005): IBS0 at most 0.005 is parent and offspring, above
  0.008 full siblings, and between the two KING breaks the tie by kinship above 0.24. Here a
  first-degree band reads `parent_child` when the IBS0 share's interval lies wholly at or below 0.005,
  `full_siblings` when wholly above 0.008, and stays `first_degree` otherwise. The interval is
  Wilson's score interval at two standard errors, so a share of 0 needs about 800 sites before it
  rules out 0.005 (at 700 it reaches 0.0057). For comparison, GWASTools' `ibdAssignRelatednessKing`
  takes 0.003 as its parent-offspring cutoff. Rejected: KING's tie-break by kinship (both kinds expect
  0.25, so it is a guess); a cutoff fitted to the sample's own pairs, as KING does when it can (a
  threshold fitted to local data, and a person has few pairs); a jackknife SE for the share, as the
  kinship has (a parent and child with no IBS0 in any block have an SE of zero, and would read
  `parent_child` from a few dozen sites). somalier's pairs are split the same way, from its own IBS0
  count over its `n` sites.
- **Unknown, with a reason.** `assembly_unknown` (a facet has no `assembly_label`; run `resolve_ref`
  first), `assembly_mismatch` (positions do not compare across builds; nothing is lifted here),
  `no_panel`, `too_few_sites` (under 500 shared, by default), `no_heterozygous_sites`, `between_bands`,
  `variant_sites_only` (a facet is homozygous reference at none of the shared sites, so its own
  genotype chose them; kinship and SE are kept, the band is not named; measured below).
- **VRS cross-check.** Where both facets carry minted VRS ids for a site's ALT, the ids are compared
  and a difference counted (`n_vrs_differ`). On one assembly there should be none.
- **somalier.** `just_vcf.somalier.somalier_relate` runs somalier (MIT) when its binary is on PATH and
  reads its `pairs.tsv` into the same frame. Not run on this machine: docs/TRUST.md has no row for it.
- **Pedigrees.** `read_pedigree` reads a GA4GH Pedigree as supplied. `expected_kinship` reads what it
  says each pair's kinship should be: the pedigree kinship coefficient through its parent links
  (KIN:003 and the terms under it, KIN:032 the other way), stated degrees where a term fixes one (full
  and half siblings, twins, grandparents, parental siblings, cousins), founders unrelated where nothing
  links them (basis `assumed`). Which term allows which kinship is just-vcf's own reading of KIN, listed
  in the `just_vcf.pedigree` docstring; a term that allows two (a sibling of either kind, a cousin)
  allows both bands. A first-degree value is expected as `parent_child` for a listed parent link,
  `full_siblings` for two shared listed parents or a stated sibling term, and `first_degree`, either
  kind, when the pedigree reaches it some other way (a twin's child and the other twin).
  `check_pedigree` sets that beside the measured relatedness of the facets whose
  `subject` is each individual's id: `agrees` when every band the measured interval (kinship ± 2 SE)
  reaches is allowed, `disagrees` when none is, otherwise `cannot_tell` with a reason: `between_bands`
  when the interval reaches allowed and other bands, `kind_not_told` when it lies wholly in the first
  degree and the IBS0 count cannot tell parent from sibling where the pedigree expects one. A disagreement
  is a row and, through `people(..., pedigree=...)`, a warning naming facet ids. It never corrects the
  pedigree or relabels a facet.

## Measured on PGP-HMS

2026-10-23, 22 facets of six participants (hu832966, huF31554, huA4F281, hu1B2689, huE18D82 and
their files on Krypton), all pairs through `pairwise_relatedness` with the packaged panel: 210 pairs in
400 s once loaded. Counts and aggregates only.

| Pairs | Band | Pairs | Kinship | Largest SE | IBS0 | Shared sites |
|---|---|---|---|---|---|---|
| one participant's facets | `same_person` | 28 | 0.437 to 0.500, median 0.492 | 0.029 | at most 0.016 | 531 to 15,244 |
| two participants | `unrelated` | 102 | −0.083 to 0.010, median −0.022 | 0.043 | 0.073 to 0.100 | 1,402 to 15,234 |
| two participants | `unknown`, `between_bands` | 12 | −0.95 to 0.062 | 1.54 | | 512 to 7,003 |
| either | `unknown`, `too_few_sites` | 47 | | | | 5 to 498 |
| either | `unknown`, `assembly_mismatch` | 20 | | | | |
| one participant's | `unknown`, `no_heterozygous_sites` | 1 | | | | 531 |

The table was measured before `variant_sites_only` existed and has not been re-run. Counted after
(2026-10-23, per facet, the five participants' facets that are not array exports): only two facets
have no homozygous-reference genotype at any panel site, the Genos exome (7,003 sites) and the
VarScan2 part of hu832966's merge (270 sites, `too_few_sites` already). Arrays write every letter
pair, so they always carry hom-ref. The rule therefore moves the exome's pairs alone: those with its
own participant's facets from `same_person`, those with other participants' from `between_bands`,
both to `unknown`, `variant_sites_only`, with `same_person` unchanged.

What the table says:

- **Every pair of one participant's facets that could be compared read `same_person`**, and none of
  two participants did: no false call either way. The facets cover 23andMe and AncestryDNA exports,
  an array exported as VCF, a Genos exome gVCF, Gencove's imputed VCF, Complete Genomics' `var` and
  VCF, a DRAGEN-style upload and PGP's freeBayes+VarScan2 merge. Both halves of F3's file (huF31554)
  read as the participant's 23andMe export does: kinship 0.4996 to 0.5.
- **Kinship below 0.5 for one person is genotype disagreement**, largest where one side is PGP's
  freeBayes+VarScan2 merge (0.437 to 0.489) or an imputed facet (0.459 to 0.475); arrays against arrays,
  CGI and DRAGEN stay at 0.495 or above.
- **`between_bands` across participants** are pairs whose interval reaches over 0.0442 (the
  third-degree floor): either few sites (one facet a 500-site VarScan2 part), or the Genos exome, whose
  kinship against every other participant sits slightly high (0.017 to 0.033, the others median
  −0.022). `same_person` still says False for all of them. The exome's offset is settled below; its
  pairs now read `variant_sites_only`.
- **The Genos exome's offset is which sites it has, chosen by its own genotype** (measured
  2026-10-23, hu832966). Its gVCF's 259,070 blocks are all no-call, so a panel site has a genotype
  only where a variant row stands: 7,003 of 15,289 sites, 4,495 heterozygous, 2,508 homozygous ALT,
  none homozygous reference (heterozygous share 0.64, arrays 0.42). Against three other
  participants' arrays it read −0.010, 0.017 and 0.033. The same participant's 23andMe export read
  −0.030, −0.012 and −0.014 against them on the full panel, and −0.021, 0.012 and 0.033 when cut to
  the exome's sites: the array reproduces the offset on those sites, so it is the selection, not
  the exome's genotypes. Two other participants' arrays cut to the same sites do not move (−0.010
  to −0.018, −0.047 to −0.076, −0.002 to −0.013, within their SEs of 0.02 to 0.05), because the
  sites were not chosen by their genotypes. A pair where either facet is homozygous reference at
  none of the shared sites now reads `unknown`, `variant_sites_only`, with its kinship and SE kept.
  `same_person` still reads the kinship (the exome against its own array is 0.5), and the pedigree
  check reads such a pair `cannot_tell`. A facet with a few hom-ref rows among many variant rows
  (a VCF that writes some 0/0 records) is chosen the same way and is not caught by this; none was
  measured.
- **`too_few_sites`** are the merges' caller-split parts (VarScan2's indel-heavy part shares under 300
  panel sites with anything). **`assembly_mismatch`** is huA4F281's headerless export, which
  `resolve_ref` placed on GRCh36 (NCBI36), against GRCh37 facets; huF31554's export was placed on GRCh37.
- **Not measured:** first- and second-degree pairs, and so the first-degree split. No PGP-HMS
  pedigree is on this machine and none of the participants measured here is related to another; the
  synthetic family in `tests/test_relatedness.py` is the only check of those bands. No VRS ids were
  minted in this run, so the VRS cross-check is tested only on synthetic data.

## Limits

- A monozygotic twin reads as the same person; genotypes cannot tell them apart.
- KING's IBS0 cutoffs assume array-grade genotypes, and facets of one person measured here disagree
  more than that: the same-person pairs above reach an IBS0 of 0.016 (one person's kinship is lowest
  where one side is PGP's freeBayes+VarScan2 merge or an imputed facet). A parent and child read from
  facets that disagree as much would show IBS0 past 0.008 and read `full_siblings`, and a pedigree's
  parent link against them `disagrees`. The test fixture shows it: a child's genotypes with one
  haplotype's allele flipped at 8% of sites (IBS0 0.013) read as the parent's full sibling. A facet's
  own error rate is not estimated, so nothing here can tell such a pair from real siblings; `ibs0`
  beside the band is what to read.
- Facets on different builds are not compared. Lifting one first (`liftover`) makes them comparable.
- The jackknife captures the sampling of sites, not the spread of real relationships (two second-degree
  relatives share anywhere from about 0.09 to 0.17), so a band near a floor can be `unknown` for a
  genuine relative. The pedigree check compares a pedigree's expected value with the measured band and
  so shares this: a real cousin can read `unrelated` (and `disagrees`), a real half sibling `third_degree`.
  The fixture's first cousins at 9,000 sites read `agrees` or `cannot_tell`, never `disagrees`, which is
  what the test asks; no PGP-HMS pedigree is on this machine to measure it on.
- A pedigree is checked pair by pair; its internal consistency (a grandparent stated where the parent
  links say otherwise) is not validated. That is an OWL reasoner's work over KIN. A stated term takes
  precedence over the parent links for its own pair.
