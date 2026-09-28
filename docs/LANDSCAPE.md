# Landscape: what exists, what we take, what is missing

Surveyed 2026-09-28 by two research agents against primary sources (repos, docs, PyPI and GitHub API
metadata). gtars was re-checked first-hand afterwards. Versions and dates are as of that day. Re-check
before depending on anything here.

**Verdict: nothing makes just-vcf redundant.** The substrate exists (polars-bio), and so do most of the
algorithms, as subprocesses or Rust crates. What does not exist is a lazy Arrow model of one person's
genome across facets, keyed by VRS and seqcol.

## Take rather than make

| Piece | Tool | How | Why |
|---|---|---|---|
| Readers, interval ops, VCF export, depth | **polars-bio** 0.36 (Apache-2.0) + datafusion-bio-formats/-functions | core dependency | VCF/BCF/VCF-Zarr/BGEN/PGEN scans, `sink_vcf`, overlap/merge/subtract, and mosdepth-compatible `depth()` from BAM for the coverage facet. Lacks normalization, gVCF block semantics, liftover, VRS and array readers |
| VRS ids, allele normalization, seqcol | **gtars** 0.10 (BSD-2, Rust + Python) and **refget** 0.12 (BSD-2) | dependency; bulk path as a Rust UDF beside polars-bio | `gtars-vrs` does digests, normalization and parallel VCF→VRS, tested against the VRS spec and vrs-python vectors. `gtars-refget` builds a refget store from FASTA and computes seqcol digests and comparisons |
| VRS conformance, Cat-VRS and CisPhasedBlock models | **vrs-python** 2.3.3, **ga4gh.cat_vrs** 0.7.2 (Apache-2.0) | dev/test only | Reference implementation as test oracle. About 6k variants/s, too slow for the hot path (vrs-python issue #534) |
| Allele-aware liftover with provenance | **bcftools +liftover** (MIT, freeseek/score) | subprocess, then read back | Only tool that handles indel REF/ALT swaps. `--write-src` and SWAP/FLIP tags give what a `liftover_to` edge needs. T2T works through UCSC hs1 chains |
| Left-align, split multi-allelic sites | **bcftools norm** | subprocess | No Arrow-native equivalent exists outside Spark (Glow) |
| N-caller same-sample merge and regional concordance | **aardvark** 1.0 (MIT, PacBio) | subprocess | Haplotype-aware merge with per-row `SOURCES`, plus agree/disagree region BEDs. Scores hard GT only |
| Consumer array sniffing, build markers, chip clusters | **snps** 2.12 (BSD-3) | vendor selected modules | Covers 23andMe, Ancestry, FTDNA, LivingDNA, MyHeritage and more. Its pandas frame, REF handling, symbolic I/D alleles and REST remap are not taken |
| Ancestry, sex, relatedness | **somalier** 0.3.5 (MIT) | subprocess | Maintained. Light pure-Python alternatives: ezancestry, fraposa-pgsc. just-prs `ancestry/` is the in-house option |
| Coordinate-only lifts (array rows) | **liftover** (MIT) or agct | dependency | Positions only |
| SVs, phasing comparison | **truvari** 5.5 (MIT), **whatshap compare** | subprocess | For later facets |
| Durable per-facet store | **VCF Zarr** via bio2zarr | option | polars-bio reads it natively. The spec has no gVCF reference blocks |

Design references, not dependencies:
- **Hail VDS** for reference blocks plus local alleles (LA/LGT/LPL).
- **Glow** for liftover status columns.
- **gnomAD's `VRS_*` parallel arrays** as the column convention for VRS in tables.

## Genuinely absent

- A per-person store of facets across assemblies, callers and phasings, with a graph between them.
- Arrow-native left-align normalization and Arrow-native VRS computation (gtars is the closest, Rust,
  with `list[dict]` Python output).
- Concordance on graded evidence (GL/PL/GP). Every benchmarking tool scores hard genotypes.
- gVCF reference blocks in Arrow.
- `liftover_to` edges, CisPhasedBlock or GA4GH Pedigree anywhere in the tabular stack. Pedigree has no
  Python tooling at all (pedigree-tools is Node.js).
- Resolving 23andMe I/D codes to real alleles (needs chip manifests), and a Genotek reader.
- A callset-type classifier (array / WES / WGS / gVCF).

## Considered and set aside

- **Stale or dead:** biobear/exon, GLnexus, hap.py, GA4GH benchmarking-tools, SURVIVOR/Jasmine, peddy,
  rye.
- **Own storage or runtime:** TileDB-VCF (own storage, Python ≤3.12), GenomicsDB, Hail (JVM/Spark),
  Glow (Spark).
- **GPL, subprocess at most:** CrossMap, vcfdist, transanno.
- **Wrong scope:** leviosam2 lifts alignments, not VCF. For us it is provenance metadata on a facet.
  anyvar is a multi-person registry service.
- **Converters that drop I/D rows or leave REF provisional:** bcftools `--tsv2vcf`, PLINK 1.9
  `--23file`. PLINK2 `--23file` is not implemented.
