# Survey of just-prs and just-dna-lite (2026-09-28)

What the two sibling repos already implement for the five starting features, with pointers for porting.
just-dna-lite reads VCF only and imports the consumer-format code from just-prs, so **just-prs is the
main source**. Paths are as of just-prs `54cd4b0` and just-dna-lite at the same date. Check
`git blame` in just-dna-lite before copying anything: only `annotation/restoration.py` is entirely the
user's (see PROVENANCE.md).

## 1. Loading

- just-prs `just_prs/arrays.py`:
  - `_read_array_text`: zip (prefers a .txt/.csv member) and gzip by magic bytes.
  - `detect_array_format`: vendor strings in the first 4 KB, then delimiter and column count. Covers
    23andMe, AncestryDNA, MyHeritage/FTDNA. LivingDNA falls through to 23andMe.
  - `normalize_array`: maps chr 23/24/25/26 to X/Y/X/MT and drops no-calls. **It writes the observed
    allele as REF** (het ref=a1 alt=a2; hom ref=alt=a1). Do not port that part. Build is assumed
    GRCh37; I/D codes pass through as letters; hemizygous calls become homozygous.
  - `detect_chip_generation`: marker-count ranges per vendor and version.
- just-prs `vcf.py`: `scan_vcf` registers under a content-hash name behind a lock, because concurrent
  `pb.scan_vcf` panics ("table already exists"), and adds +1 to the 0-based `start`.
  `detect_genome_build` uses header hints, then votes on contig lengths (at least 3 votes).
- just-prs `normalize.py`: `normalize_vcf`, `genotype_expr` (handles polars-bio's `|`-joined
  multi-allelic ALT), `_expand_pass_filters` (FILTER `.` decodes as `""`). `io_utils.open_maybe_compressed`
  detects compression by magic bytes.
- just-dna-lite `io.py` (not the user's code): `_VCF_INFO_BLOCKLIST={"ALLELE_ID"}` works around a
  polars-bio Rust panic on DRAGEN VCFs; there is also INFO/FORMAT name-collision dedup.

## 2. REF lookup

- just-prs `reference_allele.py` `resolve_reference_alleles`: the reference panel `.pvar` via DuckDB,
  then a pysam faidx lookup for SNVs only; indels stay unresolved. `_verify_fasta_build` checks chr1
  length and rejects chr-prefixed FASTAs. Its docstring says it was meant to be extractable.
- Missing everywhere: ALT derived by comparing the genotype to REF, indel anchor bases, emitting hom-ref.
- Open bug in just-prs: arrays drop no-calls, then chip-scope restoration refills those positions as
  hom-ref (`docs/refcall-resolution-plan.md` §D).

## 3. Input type

- just-dna-lite `annotation/restoration.py` (the user's code): `detect_callset_scope` separates WGS from
  targeted with thresholds `MIN_WGS_SITES=1_000_000` and `MIN_WGS_BREADTH=0.75`, measured and documented
  in the file. `infer_genotype_input_mode` detects gVCF from `<NON_REF>`/`RefCall`; run it on the
  post-filter data.
- just-prs `prs.py` `_infer_genotype_input_mode` (the same idea) and `chip_coverage.py` (GSA manifests).
- Missing: WES vs WGS vs panel, and gVCF `END` block expansion.

## 4. Population

- just-prs `ancestry/`: PCA projection (vendored FRAPOSA, MIT) plus KNN (k=20), UNKNOWN below 20% site
  coverage, simplex least-squares proportions, Bayesian consensus across panels. Wrappers in
  `prs_catalog.py` pull a ~250 MB model from Hugging Face `just-dna-seq/pgs-catalog`. For arrays,
  `assume_homref_absent=False` must be set.
- Tests: `test_ancestry.py`; methodology doc `docs/sample-ancestry-methodology.md` (it names the wrong
  Hugging Face repo).

## 5. Normalization

- just-dna-lite `annotation/chromosomes.py` (not the user's code): chr prefix handling. The gotcha
  learned there: a plain `^chr` strip turns chrM into `M` and drops mtDNA from joins.
- just-prs `liftover.py`: coordinates only (pyliftover plus UCSC chains); strand-flipped records are
  dropped, not complemented.
- just-dna-lite `vcf_export_logic.py` writes via `pb.write_vcf` but drops FORMAT/GT, so it is not a
  genotype writer.
- Missing everywhere: left-align, splitting multi-allelic sites, strand resolution, dedup, sort.

External tools that fill these gaps are in LANDSCAPE.md.
