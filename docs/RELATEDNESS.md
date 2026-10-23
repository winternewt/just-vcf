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
  never 0; an absent site is homozygous reference only where the facet's `absence_means` says so
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
- **Unknown, with a reason.** `assembly_unknown` (a facet has no `assembly_label`; run `resolve_ref`
  first), `assembly_mismatch` (positions do not compare across builds; nothing is lifted here),
  `no_panel`, `too_few_sites` (under 500 shared, by default), `no_heterozygous_sites`, `between_bands`.
- **VRS cross-check.** Where both facets carry minted VRS ids for a site's ALT, the ids are compared
  and a difference counted (`n_vrs_differ`). On one assembly there should be none.
- **somalier.** `just_vcf.somalier.somalier_relate` runs somalier (MIT) when its binary is on PATH and
  reads its `pairs.tsv` into the same frame. Not run on this machine: docs/TRUST.md has no row for it.
- **Pedigrees.** `read_pedigree` reads a GA4GH Pedigree as supplied; nothing compares it with the
  measured bands yet.

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
  −0.022). `same_person` still says False for all of them. Whether the exome's offset comes from its
  callable regions or its genotype quality is not settled.
- **`too_few_sites`** are the merges' caller-split parts (VarScan2's indel-heavy part shares under 300
  panel sites with anything). **`assembly_mismatch`** is huA4F281's headerless export, which
  `resolve_ref` placed on GRCh36 (NCBI36), against GRCh37 facets; huF31554's export was placed on GRCh37.
- **Not measured:** first- and second-degree pairs. No PGP-HMS pedigree is on this machine; the
  synthetic family in `tests/test_relatedness.py` is the only check of those bands. No VRS ids were
  minted in this run, so the VRS cross-check is tested only on synthetic data.

## Limits

- A monozygotic twin reads as the same person; genotypes cannot tell them apart.
- First degree is one band: parent and child (IBS0 near 0) and full siblings (IBS0 above 0) are told
  apart only by reading `ibs0` beside the band.
- Facets on different builds are not compared. Lifting one first (`liftover`) makes them comparable.
- The jackknife captures the sampling of sites, not the spread of real relationships (two second-degree
  relatives share anywhere from about 0.09 to 0.17), so a band near a floor can be `unknown` for a
  genuine relative.
