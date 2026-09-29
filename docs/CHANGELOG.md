# Changelog

What shipped, newest first. Cross-repo integration changes made on our side are recorded here too.

## Unreleased

- `just_vcf.sniff(path)` recognises an input file from its content: it peels gzip/BGZF, bz2, xz, zip
  and tar layers by magic bytes, reads a bounded head and matches a declarative format registry
  (`just_vcf.formats`). The answer says which negative it is when nothing matched: unknown, unreadable
  or not variant data. The build is reported only when the header states it. First formats: 23andMe,
  VCF and BCF. See docs/FORMATS.md.
- `scripts/pgp_participants.py` gathers PGP-HMS data per participant: PGP's own VCFs and Complete
  Genomics deliveries from the public Arvados cluster, plus the person's variant uploads, under
  `participants/<hu>/<source>-<id>/`. Same-person duplicates are fetched once, and the older by-type
  layout is migrated in place.
- `scripts/pgp_leech.py` downloads consented-public PGP-HMS genotype files into
  `$JUST_VCF_SAMPLES_DIR/pgp-hms/<data type>/<participant>/<download id>/`, with a sorted, resumable
  `manifest.jsonl`. It honours the site's 5 s crawl delay.
- gtars and refget join the core, for VRS allele ids and refget seqcol assembly identity.
- Package skeleton: uv library on hatchling, lean core (polars-bio, pydantic, platformdirs, typer,
  python-dotenv) with `[reference]`, `[liftover]` and `[ancestry]` extras.
- House rules (`CLAUDE.md`, charter, doc set) adopted from the repo ruleset gist, and the consumer
  triage loop adopted from its gist.
