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
| `just_vcf.reference` | just-prs `just-prs/src/just_prs/reference_allele.py` | 46d103d | not ported: reimplemented on gtars `RefgetStore`. Dropped: the reference-panel tier (a variant table as the REF source; RM3 rules it out), pysam faidx, the `^chr` strip, the SNV-only gate, eliot | the chr1-length build fingerprint became `CHR1_LENGTHS`, extended by per-contig length and md5 checks against the facet's header |
