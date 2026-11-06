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

## Reference-PCA panels in just-prs's layout

Three files per panel and assembly (`ancestry_model_<panel>_<assembly>.parquet`,
`ancestry_refpcs_<panel>_<assembly>.parquet`, `ancestry_model_<panel>_<assembly>_meta.json`), read by
`infer_ancestry(loaded, panel_dir, panel)`; `build_pca_panel` writes them from a genotype matrix. It
checks a builder's `extra_meta` (no key of its own replaced, JSON only) before writing anything, then
writes the three under temporary names and moves them into place, the meta file last, so a rebuild
that is refused or fails leaves the panel already there whole.

just-prs builds them offline with plink2 (QC, LD pruning) and numpy:

- `1000g` and `hgdp_1kg`, GRCh38, from the PGS Catalog's processed reference panels
  (`pgsc_1000G_v1`, `pgsc_HGDP+1kGP_v1`, EBI FTP). just-prs publishes them to the Hugging Face dataset
  `just-dna-seq/pgs-catalog` under `data/ancestry/`; a copy fetched from there by hand reads here
  unchanged, but none was on this machine on 2026-10-22, so that path is untested against a real
  panel.
- `aadr_ho`, GRCh37, from AADR Human Origins present-day West Eurasians. Built here by
  `build_aadr_panel` (next section).

A panel's `ref` column has to be the assembly's reference allele: an absent site in a whole-genome
facet is read as homozygous for it.

`infer_ancestry` reports the neighbour vote and its margin, the centroid proportions and the projected
PC coordinates, and by default their uncertainty from a leave-one-chromosome-out jackknife: the facet
projected again with each shared autosome left out, giving a standard error per PC coordinate
(`pc_errors`) and per proportion (`proportion_errors`), and the share of those replicates, one per
shared autosome, whose vote agrees (`vote_stability`). A facet that reads sites on fewer than five autosomes gets none, with
`jackknife_reason = too_few_blocks`; `jackknife=False` skips it. The method, what each replicate
recomputes and the minimum are in docs/SCHEMA.md, "The PCA path's uncertainty". Its tests use
synthetic panels (`tests/test_ancestry.py`); no run on a real panel is recorded yet.

## AADR Human Origins (`aadr_ho`, GRCh37): built here

AADR, the Allen Ancient DNA Resource (Reich lab, https://dataverse.harvard.edu/dataverse/reich_lab),
publishes ancient and present-day individuals in EIGENSTRAT files. Release v66.p1 is on Krypton under
`genomes/dataverse_66p1/`, with AADR's `v66.p1__files.md5sum`. The builder reads four of them, from
whatever path the caller gives:

| File | What is read |
|---|---|
| `v66.p1_compatibility_HO.aadr.patch.PUB.geno` | genotypes, packed `TGENO` (1.9 GB; only the chosen individuals' records are read) |
| `….snp` | 276,725 sites: position (hg19) and two alleles, the first one counted |
| `….ind` | 27,594 individuals: genetic id and group label |
| `v66.p1_compatibility_HO.aadr.PUB.anno` | per individual: date, data type, assessment, political entity, coordinates |

`compatibility_HO` is AADR's subset of the Human Origins sites chosen so that array and shotgun data
can be analysed together (its README), which the panel does. The plain `HO` files (584,131 sites) have
the same layout but were not built here. The `.geno` may also be the older site-major packed `GENO` or
text EIGENSTRAT; those two were tested on synthetic files only.

```python
from pathlib import Path
from just_vcf.ancestry import build_aadr_panel, infer_ancestry

aadr = Path("/mnt/krypton/Genomes/genomes/dataverse_66p1")
build_aadr_panel(aadr / "v66.p1_compatibility_HO.aadr.patch.PUB", aadr / "v66.p1_compatibility_HO.aadr.PUB.anno", Path("panels"))
estimate = infer_ancestry(loaded, Path("panels"), "aadr_ho")
```

It needs a refget store holding one GRCh37 collection (`reference=`, default the library's own store,
`JUST_VCF_REFERENCE_DIR`), and refuses without one: every site's REF is read from it.

**Who goes in.** Each rule in order; the first one an individual fails is the reason it is counted
under:

1. listed in the `.anno`. An `.anno` that lists one Genetic ID twice is refused, naming how many
   repeat, since the join would read that person's genotypes twice (v66.p1's has no repeat);
2. present-day (`Date mean in BP` 0): the ancient individuals are out, as in just-prs;
3. diploid data (`Data type` `HumanOrigins` or `Shotgun.diploid`): pseudo-haploid shotgun and capture
   data, one read per site, would look homozygous everywhere, and the reference genomes (human, chimp,
   gorilla, ancestral) are not people;
4. `ASSESSMENT` `Pass`, `PROVISIONAL_PASS` or `MERGE_PASS`;
5. a group label without `Ignore` at its start and without a QC qualifier (`group_label_flag`). A
   label is a name followed by any number of `-` qualifiers, and every one is read. AADR's README names
   two marks of individuals to filter out of primary analyses, `Ignore_` and outliers; an outlier is
   `-o` or any `-o…` (`-oPCA`, `-oRelative`, `-oAfrica`, `-olowEEF`). The other QC qualifiers,
   `QC_QUALIFIERS`, were read off v66.p1's `.ind`: `-QCremove`, `-Dup`, `-Discovery` (the Human Origins
   ascertainment samples), `-NotPermittedToPublish` and `-DontUseHighlyProblematic`. A qualifier of
   digits only is a numbered sub-group (`Kumyk-1`, `Uzbek-2`, `Hazara_1-1`) and is kept, under its label
   as written. A qualifier that is neither is counted apart (`group label qualifier not recognised`),
   with a warning naming the labels; in v66.p1 those (`-A`, `-B`, `-C`, `-lowEastEU`, `-highEastEU`,
   `-Africa`, `-3CentralEuropean`) are all on ancient individuals, so rule 2 has already left them out;
6. a latitude and longitude;
7. inside the box, by default just-prs's West Eurasia: latitude 30 to 75, longitude -25 to 90
   (`box=` takes another);
8. a region for its political entity (below);
9. once per person: where AADR lists one person under two genetic ids (`Individual ID` shared, usually
   a Human Origins genotype and a shotgun genome), the `HumanOrigins` one is kept, else the first.

`population` is the `.ind` group label as written (172 of them in the default build). A numbered
sub-group stays a population of its own: nothing in the panel code needs the suffix gone, and AADR
lists `Kumyk-1` and `Kumyk-2` apart. Two things follow, measured on v66.p1. Several numbered labels
are the Human Origins genotype of a person whose shotgun genome is under the bare label (`Hazara_1-1`
and `Hazara_1`, `Italian_North-1` and `Italian_North`), and rule 9 keeps the Human Origins one, so
11 people moved from the bare label to the numbered one (`Hazara_1` emptied into `Hazara_1-1`; two of
`Italian_North`'s 40 are now `Italian_North-1`). And a numbered sub-group holds one to five people, so it cannot lead a 20-neighbour vote at
`resolution="population"`. `superpop` is a region:
the political entity's United Nations M49 region, Europe as one region and Sub-Saharan Africa as one,
with Russia split at 60°E into Europe and Siberia, which M49 does not do (`AADR_REGIONS`,
`region_of`). Where AADR's samples were collected away from their ancestry's region (the USA, Barbados:
1000 Genomes' CEU, ASW, ACB) there is no region and they are left out. A region is geography, not
genetics: Kalmyks, west of the Urals, are in Europe, Iran is in Southern Asia, and the box's slice of China is
Xinjiang. For a consensus, a region folds onto a continental group through `continental_of` on its
name: Europe to EUR, Northern and Sub-Saharan Africa to AFR, Latin America to AMR. The Asian regions,
Oceania and Siberia fold onto nothing and are dropped, as Middle Eastern groups are: their samples are
not what 1000 Genomes' EAS or SAS sampled.

**Which sites.** Each rule in order: an autosome; one base each, A, C, G or T; not A/T or C/G (Human
Origins has none by design, so this drops nothing here); one record at the position; a contig the
GRCh37 store holds; one of the two alleles the reference base there, the alleles swapped where it is
the second (REF is read from the store, never taken from the file or a genotype); minor allele
frequency at least 0.02 over the chosen individuals. No LD pruning, as in just-prs. The panel then
keeps 20 components (just-prs's choice, for fine structure within Europe); `knn_k` stays 20.

**Measured 2026-10-26** from Krypton, GRCh37 store `/data/just-dna-cache/reference/store`:

| | |
|---|---|
| Individuals | 2,567 of 27,594. Out: 19,120 ancient, 244 not diploid, 1 assessment, 370 flagged labels, 0 unrecognised qualifiers, 541 without coordinates, 4,176 outside the box, 575 second representations. 22 of the 2,567 are in a numbered sub-group (16 labels) |
| Regions | Europe 1,394, Western Asia 347, Southern Asia 307, Siberia 216, Central Asia 166, Northern Africa 105, Eastern Asia 32 |
| Sites | 184,882 of 276,725; every one of the 276,725 had the GRCh37 base as its first allele (none swapped, none dropped for it); 91,843 below the frequency floor |
| Files | sites 116 MB, reference samples 1.6 MB |
| Build | 17 to 24 s and peak RSS 6.1 GiB with the shipped code, the `.geno` read warm; 50 to 54 s on a first, cold read with an earlier cut that held a float copy of the matrix (peak 13.4 GiB) |
| Leave-one-out | 0.942 by region. Per region: Central Asia 163/166, Europe 1,351/1,394, Northern Africa 103/105, Siberia 195/216, Southern Asia 290/307, Western Asia 317/347, Eastern Asia 0/32 (its Xibo, Uyghur and Kyrgyz samples are voted Central Asian) |
| hu6ABACE (PGP-HMS, DRAGEN SNV genome, GRCh37) | 184,881 of 184,882 sites used: 85,981 observed, 98,900 absent read as reference, 1 unusable. Region Europe, all 20 neighbours; population GBR, margin 0.55. Load and placement 41 s |

The build before numbered sub-groups were kept (rule 5 then excluded every `-` qualifier) chose 2,556
individuals and 184,833 sites, leave-one-out 0.944; hu6ABACE was placed the same way.

The regions overlap where the populations form a cline. `tests/test_ancestry_panels.py` rebuilds each
region's first reference sample from the panel's factors and places it. On this panel two land across
a border (2026-10-26, the build above): the first Central Asian sample, a Karakalpak, is placed in
Europe with 0.35 of its neighbours Central Asian and 0.60 European, and Eastern Asia's first, a Xibo,
in Central Asia, 0.40 Eastern Asian and 0.60 Central Asian. The check accepts that (the person's choice,
2026-10-26): a sample is placed well when it is in its own region, or when the vote is two-way (its own
region and the leader hold at least 0.9 of it together) and its own region holds at least a third. The
same panel with its region labels shuffled fails the check for six of its seven regions (their votes
spread over four to six regions), and the test asserts that it fails for at least one, so the looser
rule still catches a panel that places its samples among unrelated regions.

AADR is a published research dataset. The builder reads it from the caller's path and nothing derived
from it is in the repository: the tests write synthetic EIGENSTRAT files with made-up individuals.

## Combining estimates from several panels

`consensus(estimates)` folds each estimate's vote and proportions onto the five 1000 Genomes groups.
Panels that share reference samples are one witness: their methods see the same people at the same
sites, so agreeing with each other is not more evidence. Inside a witness the answer is the mean of
its readings, with their spread per group beside it; witnesses are multiplied as independent only
where the caller declares their panels' reference samples disjoint. Nothing here can tell that from
the panels, so by default every panel given is one witness:

```python
fused = consensus([infer_ancestry(loaded, panels, "aadr_ho"), prive_proportions(loaded, panels)])
# one witness, combination "mean"
fused = consensus([...], disjoint_panels=[("aadr_ho", "prive")])
# two witnesses, combination "product_of_means": only if the caller knows no person is in both
```

Which of the panels above share people: `hgdp_1kg` holds `1000g`'s (HGDP+1kGP is the two projects
together); `aadr_ho` holds 1000 Genomes people (hu6ABACE's nearest population in it is `GBR`; read from its
labels, not checked by sample id), so it shares with both. Privé's groups are drawn from UK Biobank participants (Privé 2022); whether any of
them is in another panel here has not been checked, so no pair of these panels is known disjoint. The rule and
what was rejected for it are in docs/SCHEMA.md, "Combining ancestry estimates".

## Testing against a panel

`JUST_VCF_ANCESTRY_DIR` names one directory holding any of the above, and may hold AADR's
`compatibility_HO` release at its top level or one directory down (`dataverse_66p1/`):
`tests/test_ancestry_aadr.py` then builds `aadr_ho` from it into a temporary directory and places
PGP-HMS hu6ABACE, read from `JUST_VCF_SAMPLES_DIR`, with the GRCh37 store `JUST_VCF_TEST_REFERENCE_DIR`
names; it skips without any of the three. `tests/test_ancestry_panels.py`
draws people from the Privé groups' own frequencies, and rebuilds each PCA panel's reference samples
from its stored factors, and checks each lands in its own group or on that group's border with one
other (above), and that a shuffled copy of the panel fails the same check. Both skip without the
variable; the check's rule itself is tested always, on a synthetic panel.
