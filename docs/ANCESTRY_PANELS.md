# Ancestry reference panels

`just_vcf.ancestry` (RM5) places a facet against a reference panel the caller passes as a local path.
Nothing is downloaded and no panel ships with the package. This page says which panels work, where
they come from and what each one assumes. Every estimate is a placement relative to the panel's
groups, with its uncertainty; none is a statement of anyone's identity.

A facet is matched to a panel only on the facet's own `assembly_label`. A facet without one, or on
an assembly the panel does not have, gets an unknown estimate (`assembly_unknown`,
`assembly_mismatch`); nothing is lifted over here.

## Privé's worldwide reference (works here, GRCh37)

Privé 2022, "Using the UK Biobank as a global reference of worldwide populations" (Bioinformatics),
as used by the bigsnpr `ancestry` vignette. Two public files on figshare:

| File | figshare | Columns |
|---|---|---|
| `ref_freqs.csv.gz` | file 31620968 | `chr, pos, a0, a1, rsid`, then the frequency of `a1` in 21 groups: `Africa (West)`, `Africa (South)`, `Africa (East)`, `Africa (North)`, `Middle East`, `Ashkenazi`, `Italy`, `Europe (South East)`, `Europe (North East)`, `Finland`, `Scandinavia`, `United Kingdom`, `Ireland`, `Europe (South West)`, `South America`, `Sri Lanka`, `Pakistan`, `Bangladesh`, `Asia (East)`, `Japan`, `Philippines` |
| `projection.csv.gz` | file 31620953 | `chr, pos, a0, a1, rsid, PC1..PC16` (the loadings) |

A copy of both is on the Krypton array under `genomes/prive/` (checked 2026-10-22 by their header
lines; just-prs's `prive.py` reads the same two files). Build once, then point
`prive_proportions` at the output directory:

```python
from pathlib import Path
from just_vcf.ancestry import build_prive_reference, prive_proportions

build_prive_reference(Path("ref_freqs.csv.gz"), Path("projection.csv.gz"), Path("panels"))
estimate = prive_proportions(loaded, Path("panels"))
```

The build joins the two files on `chr, pos, a0, a1` into `prive_reference.parquet` (5,816,590 sites,
1.4 GB, 140 s from Krypton on 2026-10-22). The method projects the genome onto the 16 PCs with the
vignette's per-PC shrinkage correction and solves for group proportions by simplex-constrained least
squares; the standard errors come from a bootstrap over 200 contiguous blocks of sites.

`a0` is the GRCh37 reference allele, which matters because an absent site in a whole-genome facet is
read as homozygous `a0`. Checked 2026-10-22 on a PGP-HMS DRAGEN 3.2.8 SNV genome (hu6ABACE, GRCh37,
3,863,397 records): where its records met the reference's positions, its REF was `a0` at 2,389,449
of 2,389,956 and `a1` at 252.

Measured the same day on PGP-HMS data (counts only):

| Facet | `absence_means` | Sites used of 5,816,590 | observed / absent as reference / unusable | Largest standard error | Time |
|---|---|---|---|---|---|
| DRAGEN SNV genome (hu6ABACE) | `weak_reference` | 5,814,017 | 2,387,383 / 3,426,634 / 2,573 | 0.016 | 13 s after an 87 s load |
| 23andMe build 37 export (hu832966) | `unknown` | 888,688 | 888,688 / 0 / 6,761 | 0.016 | 6 s |

Caveats: GRCh37 only; South America is the only American group; Middle Eastern weight has no
continental counterpart and is dropped from a consensus. The two files are data from a GPL-licensed
project's publication; just-vcf reads them and redistributes nothing.

## Reference-PCA panels in just-prs's layout (none built on this machine)

Three files per panel and assembly (`ancestry_model_<panel>_<assembly>.parquet`,
`ancestry_refpcs_<panel>_<assembly>.parquet`, `ancestry_model_<panel>_<assembly>_meta.json`), read by
`infer_ancestry(loaded, panel_dir, panel)`; `build_pca_panel` writes them from a genotype matrix.

just-prs builds them offline with plink2 (QC, LD pruning) and numpy:

- `1000g` and `hgdp_1kg`, GRCh38, from the PGS Catalog's processed reference panels
  (`pgsc_1000G_v1`, `pgsc_HGDP+1kGP_v1`, EBI FTP). just-prs publishes them to the Hugging Face dataset
  `just-dna-seq/pgs-catalog` under `data/ancestry/`; a copy fetched from there by hand reads here
  unchanged, but none was on this machine on 2026-10-22, so that path is untested against a real
  panel.
- `aadr_ho`, GRCh37, from AADR Human Origins present-day West Eurasians (an EIGENSTRAT reader and
  selection in just-prs `ancestry/aadr.py`, `eigenstrat.py`, not ported). The AADR files are on Krypton
  under `genomes/dataverse_66p1/`.

A panel's `ref` column has to be the assembly's reference allele: an absent site in a whole-genome
facet is read as homozygous for it.

## Testing against a panel

`JUST_VCF_ANCESTRY_DIR` names one directory holding any of the above. `tests/test_ancestry_panels.py`
draws people from the Privé groups' own frequencies, and rebuilds each PCA panel's reference samples
from its stored factors, and checks each lands where it was drawn from. Both skip without the variable.
