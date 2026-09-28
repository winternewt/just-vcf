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

MIT licensed.
