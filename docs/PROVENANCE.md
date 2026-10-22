# Provenance of borrowed code

Every piece ported from a sibling repo gets a row: where it came from (repo, path, commit), what was
cut, and what was changed. Update this in the same commit as the port.

Licensing: just-prs is MIT. just-dna-lite is AGPL. From just-dna-lite, only
`just-dna-pipelines/src/just_dna_pipelines/annotation/restoration.py` may be copied. Every line of it
is the user's (git blame, 2026-09-28), and the user relicenses it to MIT here. `annotation/chromosomes.py`
and `io.py` are entirely antonkulaga's, and `QualityFilters` in `module_config.py` is mostly his, so
those are reimplemented from their behaviour, never copied, unless he gives permission. Check `git blame`
before porting anything else from that repo.

| just-vcf module | Source | Commit | Cut | Changed |
|---|---|---|---|---|
| `just_vcf.loaders.vcf.scan` | just-prs `src/just_prs/vcf.py` `scan_vcf` / `_vcf_registration_key` | c424bec | eliot, `chunk_size`/`concurrent_fetches` overrides, `format_fields` | table name prefixed `just_vcf_`; the 0-based `start` it converts is checked against the file's first record on every load |
| `just_vcf.callset.spread`, `classify` | just-dna-lite `just-dna-pipelines/src/just_dna_pipelines/annotation/restoration.py` `detect_callset_scope`, `_coverage_breadth`, `MIN_WGS_SITES`, `MIN_WGS_BREADTH`, `restoration_max_flank_bp` default (`configs.py`, a constant only) | 18d28d3 | the two-way `CallsetScope` (WGS / targeted), the restoration context, `infer_genotype_input_mode` (its `<NON_REF>` test was already in `loaders.vcf._is_block_record`, stricter; its `RefCall` test is not ported, because a DeepVariant `RefCall` row is a call here and does not change what an absent position means) | targeted split into `wes` and `panel` with unknown bands between every class; symbolic-only rows do not count as sites; a read-evidence FORMAT key is required and an imputation marker wins; the result carries a reason and its measurements; breadth computed per contig in one lazy aggregation instead of a sorted eager frame |
| `just_vcf.assembly.infer_assembly` | just-prs `src/just_prs/vcf.py` `detect_genome_build` (the contig-length vote) and its GRCh37/GRCh38 length tables | fce5e26 | the header-hint pass (sniff's build hints already do it), the bare majority vote | T2T-CHM13 added; lengths checked against UCSC `chrom.sizes`; a conflicting length, or lengths matching two assemblies, gives no label instead of the majority; the mitochondrion only contradicts. just-prs's array build "detection" is a `GRCh37` default, so nothing of it was ported: an array table's label stays what its header states, or null |

## Public reference data

Not code, but taken from outside and carried in the package; each row says what was derived.

| just-vcf constant | Source | Files | Derived |
|---|---|---|---|
| `just_vcf.karyotype.PAR_BOUNDS`, `X_DEGENERATE_Y` (RM269) | GIAB genome stratifications v3.6, `XY` group (NIST, public), https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/ , downloaded 2026-10-22 | `GRCh37_*` and `GRCh38_*`: `chrX_PAR`, `chrY_PAR`, `chrY_nonPAR`, `chrY_XTR`, `chrY_ampliconic` (`.bed.gz`, carried decompressed in `assets/giab_stratifications_v3.6_xy/` with GIAB's README for each build, renamed `.txt`). GIAB built the XTR and ampliconic files from SexChrLab's SexChrCoordinates (commit 4ad1faf, per the README) | just-vcf's derivation, not a GIAB product. PAR: each arm from the first to the last base of GIAB's PAR intervals (GIAB leaves gaps out). X-degenerate Y: `chrY_nonPAR` minus `chrY_XTR` minus `chrY_ampliconic`, kept between PAR1 and PAR2, which drops the single base GRCh38's `chrY_nonPAR` lists before PAR1. BED (0-based, half-open) becomes 1-based closed `pos`. `tests/test_karyotype.py` re-derives both from the carried files. Two places where GIAB's PAR is not NCBI's published figure, kept as GIAB has them: GRCh37 Y PAR1 ends at 2,749,806 (NCBI: 2,649,520) and PAR2 starts at 59,133,471 (NCBI: 59,034,050); several GIAB PAR starts are one base after NCBI's 1-based start (GRCh38 X PAR2 at 155,701,384, NCBI 155,701,383). NCBI's figures are as commonly cited, not re-checked here. If they hold, GRCh37's last X-degenerate interval (58,967,657–59,133,470) reaches about 99 kb into NCBI's PAR2, about 1% of the set; it is kept as derived and listed as open in RM269 |
