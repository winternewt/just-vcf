# just-vcf

A [polars-bio](https://github.com/biodatageeks/polars-bio) based library for converting consumer and
personal genomics files into VCF. It covers 23andMe-style array exports and other freeform genotype
tables that polars-bio's own VCF readers don't handle.

Early stage. Nothing is implemented yet beyond the package skeleton. The plan is in
[docs/ROADMAP.md](docs/ROADMAP.md).

```bash
uv sync                       # core
uv sync --extra reference     # + FASTA REF lookup (pysam; not on Windows)
uv sync --all-extras          # + liftover, ancestry
uv run pytest
```

The command line prints counts and aggregates only, never rows, sample names or paths:

```bash
uv run just-vcf sniff FILE                      # format, family, container layers, stated build
uv run just-vcf load FILE [--sample NAME]       # each facet: assembly, callset type, caller, chip, row counts
uv run just-vcf check FILE --reference REF      # REF read from the reference; ref_check and VRS counts
uv run just-vcf convert FILE OUT.vcf.gz [--reference REF] [--normalize] [--lift CHAIN --target-reference FASTA] [--vrs] [--gvcf]
```

Every command takes `--json`. Exit codes: 0 done, 1 crash, 2 usage, 3 unsupported format, 4 several
samples (or `--sample` not in the file), 5 no sample, 6 reference mismatch, 7 damaged BGZF, 8 table
layout, 9 rows that cannot be VCF, 10 unsorted rows, 11 several facets without `--facet`.

MIT licensed.
