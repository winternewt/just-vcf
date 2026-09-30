# Changelog

What shipped, newest first. Cross-repo integration changes made on our side are recorded here too.

## Unreleased

- `just_vcf.load(path)` loads the 23andMe-layout genotype tables (23andMe tab and comma, Living DNA,
  SelfDecode, PLINK `--recode 23`, Genes for Good measured and imputed, and the preamble-less rsid
  tables) into the `facets`, `calls` and `blocks` frames. The genotype text is kept verbatim, no REF is
  written, and rows that do not fit are dropped and counted in `Loaded.warnings`, one line per reason.
  Other formats raise `UnsupportedFormatError`. Over the 503 PGP-HMS files of these formats: no
  errors, 432 million rows. Schema 0.0.2 loosens two "never null" rules the real files broke.
- `just_vcf.open_payload(sniffed)` reopens the payload `sniff` matched, replaying its compression and
  archive layers (nested zips and tar members included), so a loader reads exactly the bytes that were
  recognised.
- `just_vcf.schema` holds the frame contract of docs/SCHEMA.md v0.0.1: the `FACETS`, `CALLS` and
  `BLOCKS` schemas, the closed vocabularies as enums, canonical contig names (`chrM`/`M` → `MT`, array
  codes 23–26, `25` kept as `XY`) and the canonical row order. A test fails when the schema tables in
  the doc and the code disagree.
- `just_vcf.sniff(path)` recognises an input file from its content: it peels gzip/BGZF, bz2, xz, zip
  and tar layers by magic bytes, reads a bounded head and matches a declarative format registry
  (`just_vcf.formats`). The answer says which negative it is when nothing matched: unknown, unreadable
  or not variant data. The build is reported only when the header states it. First formats: 23andMe,
  VCF and BCF. See docs/FORMATS.md.
- VCF dialects that parse differently get their own ids: Complete Genomics `vcf_cgi`, PGP's merged
  `vcf_freebayes_varscan2`, `vcf_no_meta` (a `#CHROM` line without meta lines) and the
  `freebayes_varscan2_keyvalue` table. `##fileformat` may now come after other meta lines, and build
  hints read CGI and common reference file names. gVCF stays `vcf`: whether a file is a gVCF is decided
  from its rows later, since headers mislead both ways.
- Other genotype-bearing tables: Illumina GenomeStudio Final Reports (with and without a column
  line), Oxford `.gen` imputation output, Color variant sets, Illumina Nirvana JSON, GET-Evidence GFF
  (build 36 and b37), FTDNA Big Y, CASAVA SNP and indel tables and SOAPsnp consensus. `sniff` also
  names genetic-genealogy files that hold no positioned calls (Y-STR, Y-SNP panels, the rCRS/RSRS
  table, Genographic reports, ancestry segments) and GenBank ASN.1 submissions.
- Vendor genotype tables: AncestryDNA, Genos, Color genotypes, deCODEme (and its build-stating info
  file), National Geographic Genographic Geno 2.0 (with or without its column line), FamilyTreeDNA in comma and 2010-era space-separated
  layouts, and a neutral `RSID,CHROMOSOME,POSITION,RESULT` fallback for unbranded exports.
- 23andMe-shaped genotype tables: Genes for Good (unphased and imputed), 23andMe (including 2026
  wording and a spreadsheet re-save), Living DNA, SelfDecode, PLINK's 23andMe-style output and the
  Promethease tab-genotype report, plus two shape-only fallbacks for header-stripped and headerless
  rsid tables. `sniff` now normalises CRLF and bare-CR line endings before matching.
- Complete Genomics native files are recognised by `#TYPE` and their `>` header line: `cgi_var`,
  `cgi_master_var` and `cgi_summary`, format versions 1.6 to 2.5. `FormatSpec.header_prefix` names a
  marker on the column-name line.
- `scripts/sniff_corpus.py` sniffs every file under a directory and writes a TSV report with an
  aggregate by outcome and format; it is how the PGP-HMS corpus is checked against the registry.
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
