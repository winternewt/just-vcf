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
