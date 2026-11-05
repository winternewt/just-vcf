# just-vcf

Turn consumer and personal genomics files into VCF, and say how sure each step is.

just-vcf reads what another tool already called (a 23andMe or AncestryDNA raw data download, a
genotype table, a VCF or gVCF from a sequencing pipeline) and loads it into lazy
[polars](https://pola.rs) frames on top of [polars-bio](https://github.com/biodatageeks/polars-bio).
From there it can read REF from a reference, normalize, lift over to another assembly, mint GA4GH VRS
ids, and write VCF, gVCF or VCF Zarr. It does not call variants, does not impute genotypes, and its
output is not clinical-grade.

**Pre-1.0.** The Python API, the frames and the command line may change between minor releases. Every
break is recorded in the changelog.

## What it does

- **Loads** the file's format, recognised from its content (compression and archive layers peeled by
  magic bytes, never by extension). Array and table exports: 23andMe (every chip version, and a copy
  re-saved as CSV), AncestryDNA, FamilyTreeDNA (Family Finder and Big Y), Living DNA, Genes for Good,
  deCODEme, SelfDecode, Genos, Color, Illumina GenomeStudio final reports, PLINK `--recode 23`, and
  bare rsid tables. Sequencing outputs: VCF, BCF and gVCF, Complete Genomics `var` and `masterVar`,
  Illumina CASAVA SNP and indel tables, Nirvana JSON, SOAPsnp and GET-Evidence GFF. `LOADABLE_FORMATS`
  lists them all.
- **Describes each facet** it loads: the genome build (stated by the file or inferred from it, and
  unknown when it cannot tell, never a default), the callset type (whole genome, exome, panel, array
  or imputed), the vendor and caller, and for arrays the chip. A chip's absent site is unknown, not
  reference. An imputation server's output keeps its per-site quality (minimac's R2 and ER2, Beagle's
  DR2), and `state_imputed_from` records which chip it was made from, so the two count once.
- **Reads REF from a reference** (`resolve_ref`) kept in a local [refget](https://ga4gh.github.io/refget/)
  store; a FASTA is added to it once. The source's own REF is checked against it, never trusted, and
  never written from the genotype.
- **Normalizes** (`normalize`): splits multi-allelic sites, left-aligns and trims indels.
- **Lifts over** (`liftover`) through a UCSC chain file, checked against both assemblies first, with
  a report row for every row that did not lift.
- **Mints VRS allele ids** (`mint_vrs_ids`) with [gtars](https://github.com/databio/gtars).
- **Writes** VCF 4.3, or 4.5 where a row holds a phase set (`write_vcf`, BGZF with a tabix index),
  gVCF with the reference blocks (`blocks=True`), and VCF Zarr (`write_vcf_zarr`).
- **Compares** facets: `concordance` allele by allele, `relatedness` (KING-robust) and `same_person`,
  and `people` to group one person's files and check they agree. `query` answers counts and presence
  questions about a genome without returning a row.
- **Reads coverage** as a facet of its own: `read_coverage` takes what mosdepth or `samtools depth` /
  `samtools coverage` already wrote, and `alignment_coverage` measures depth over a BAM or CRAM with
  polars-bio's `depth()`, in fixed windows.
- **Estimates sex-chromosome karyotype**: `sex_chromosome_evidence` per facet, and `person_karyotype`
  for one answer per person across their files, where a liftover of a file counts once and files that
  disagree withhold the answer.

## Install

```bash
pip install just-vcf
# or
uv add just-vcf
```

Python 3.13 or later. Ancestry estimation (`just_vcf.ancestry`: PCA projection onto a reference
panel, and `consensus` across several estimates, where one panel's methods count as one witness)
needs numpy, which comes with the `ancestry` extra:

```bash
pip install "just-vcf[ancestry]"
```

## Use it from Python

```python
import just_vcf

loaded = just_vcf.load("23andme_export.txt")
print(loaded.facets.select("format_id", "assembly_label", "callset_type", "chip").collect())

store = just_vcf.ReferenceStore.open()   # the local store, under the user cache directory
store.add_fasta("GRCh37.fa")             # read once; later calls find it there

loaded = just_vcf.resolve_ref(loaded, store)
loaded = just_vcf.normalize(loaded, store)
loaded = just_vcf.mint_vrs_ids(loaded, store)
written = just_vcf.write_vcf(loaded, "genome.vcf.gz", resolved_only=True)
print(written.rows_written, written.rows_skipped)
```

`load` returns the facets, calls, reference blocks and links as lazy frames, and `loaded.warnings`
counts every row it could not read, by reason. Rows that cannot be valid VCF (an array call whose REF
could not be read, say) raise `UnwritableRowsError` unless `resolved_only=True` writes the rest.

## Command line

The command line prints counts and aggregates only, never rows, sample names or paths:

```bash
just-vcf sniff FILE                      # format, family, container layers, stated build
just-vcf load FILE [--sample NAME]       # each facet: assembly, callset type, caller, chip, row counts
just-vcf check FILE --reference REF      # REF read from the reference; ref_check and VRS counts
just-vcf convert FILE OUT.vcf.gz [--reference REF] [--normalize] [--lift CHAIN --target-reference FASTA] [--vrs] [--gvcf]
```

Every command takes `--json`. Exit codes: 0 done, 1 crash, 2 usage, 3 unsupported format, 4 several
samples (or `--sample` not in the file), 5 no sample, 6 reference mismatch, 7 damaged BGZF, 8 table
layout, 9 rows that cannot be VCF, 10 unsorted rows, 11 several facets without `--facet`.

## Configuration

| Variable | Default | Meaning |
|---|---|---|
| `JUST_VCF_CACHE_DIR` | the platform's user cache | rewritten inputs, spills and counts the lazy frames read; `just_vcf.clear_cache` removes them |
| `JUST_VCF_REFERENCE_DIR` | `reference/` under the cache | the refget store REF is read from |
| `JUST_VCF_MANIFESTS_DIR` | unset | Illumina manifest zips, used to resolve an array's `I`/`D` indel codes; unset, those rows stay unread, with a warning |

The command line also reads them from a `.env` file in the working directory or the nearest one above it.

## License

MIT.
