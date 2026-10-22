# Changelog

What shipped, newest first. Cross-repo integration changes made on our side are recorded here too.

## Unreleased

- VRS allele ids (RM271): `mint_vrs_ids(loaded, reference)` fills `calls.vrs_allele_ids` after
  `resolve_ref`, one VRS 2.x id per slot of `[ref] + alts`, equal to vrs-python's on every kind of
  allele tested and on 19,000 real PGP-HMS alleles; `vrs_withheld(calls)` says why an id is missing.
  Substitutions now wait for REF to be read, like indels (SCHEMA.md; no version bump, no column
  changed). vrs-python joins the dev dependencies as the test oracle.

- Schema 0.0.7 (RM3): `facets.ref_seqcol_digest` is new; `resolve_ref` fills `ref`, `ref_check`, and an
  array row's `alts` and `gt` from a reference.

- Schema 0.0.6 (RM2's VCF side): `facets.ploidy_statements` is new, and the caller, its version, the
  vendor and the phasing method are filled from VCF headers. `REWRITE_VERSION` 6 and the schema bump
  orphan earlier cache entries (RM268).

- Schema 0.0.5 (RM4): `facets.sorted_name_length_pairs_digest` is new, and VCF facets now fill
  `assembly_label` (basis `inferred`), `callset_type` and `absence_means` where the evidence settles them.
  Every earlier frame spill in the cache is orphaned by the version bump (RM268).

- A VCF facet now gets `callset_type` (RM4): `wgs`, `wes`, `panel` or `imputed` from the spread of
  its calls and its header, ported from just-dna-lite `detect_callset_scope` (docs/PROVENANCE.md),
  and null wherever the evidence falls between classes. `absence_means` follows: `weak_reference`
  for `wgs` without blocks, `blocks` wherever blocks exist, `unknown` otherwise. Before, every VCF
  facet without blocks said `unknown`, so a whole-genome VCF's facet changes value here (minor: no
  column changes). On PGP-HMS uploads: four Genos exome VCFs read `wes` (44,020–44,796 sites, 13.6%
  breadth, against the 75% a genome needs); DRAGEN genomes split SNV-only read `wgs` (3.46M and
  3.84M sites, 92–93% breadth) and `weak_reference`; two Gencove files (`##source=loimpute`) and a
  minimac upload (INFO R2, ER2) read `imputed`; a five-sample DELLY SV merge, an array-derived VCF of
  542,266 sites, an indel-only DRAGEN split (958,046 sites), a smaller SNV-only upload (831,137
  sites) and an all-sites file of the sex chromosomes stay unknown.

- A VCF facet's `assembly_label` is inferred from its `##contig` lengths when the header states no
  build (RM4): three or more primary chromosomes must match one of GRCh37, GRCh38 or T2T-CHM13 and
  none may match another or none, else it stays null. A stated label is never replaced; a header whose
  lengths disagree with it, or match two assemblies, gets a warning with the counts. New facet column
  `sorted_name_length_pairs_digest`: the seqcol attribute digest of the header's contig names and
  lengths, checked against refget's compliance vectors. `seqcol_digest` stays null, because its
  top-level digest needs per-sequence digests no header carries (docs/SCHEMA.md). Minor: one new
  nullable column; facet ids do not change.

- `load` reads eight more genotype tables (RM2), through a new loader driven by each format's spec:
  `ancestrydna`, `ftdna`, `ftdna_space`, `rsid_result_csv` (where MyHeritage exports land),
  `color_genotypes`, `decodeme`, `genos` and `illumina_final_report`. 18 of the 39 declared formats
  now load. Columns are found by name; a file whose columns lack what its format needs is refused
  with the new `TableLayoutError` (exported), for instance a Final Report with only Forward or Design
  alleles, which need the chip manifest to reach the plus strand. `load(sample=…)` also picks one
  sample of a table with a sample column (a Final Report, Color), and refuses several samples with
  no choice, as for a VCF. No REF is written; an absent site means unknown. AncestryDNA's doubled
  letters on haploid sites stay two alleles, the vendor's convention. Color's 0-based positions are
  shifted by one and its build read from the rows; deCODEme's `-` strand rows are complemented in
  `gt_alleles` (`source_gt` keeps the letters) and its build comes from the info member of its zip;
  Genos insertion sites sit at the base before them with `end_pos` null until RM3. Rows on
  chromosome 0, a column line repeated inside the file and unreadable genotypes are each one counted
  warning. Columns not carried are listed in `facets.dropped_fields` by their file names.
  `FormatSpec` gains `HeaderMode.SECTION` and a `section` field (the Final Report's `[Data]`). Minor:
  new formats, a new exception and a new spec field; nothing that loaded before changes. Measured on
  the PGP-HMS copy: 41 AncestryDNA, 103 FTDNA (7 of them two exports joined), 2 FTDNA space-separated,
  4 deCODEme (build 36 from the zip), 4 Genos, 1 Color and 1 RSID/RESULT table load with no error; of
  the 2 Final Reports, the GSA one loads (665,604 calls) and the Forward/Design one is refused.

- A CGI no-called span written with a haploid `.` (male X and Y) now goes to `blocks` like `./.`.
  polars-bio 0.36 returns a haploid `.` GT as null, so these 139,310 records of a PGP-HMS vcfBeta had
  stayed in `calls` with no genotype at all. The same polars-bio behaviour makes a haploid `.` in any
  VCF load with `gt`, `gt_alleles` and `source_gt` null rather than one null haplotype (open).

- `load` reads PGP-HMS's freeBayes+VarScan2 merges (`vcf_freebayes_varscan2`, RM2) as two facets in
  one `Loaded`, one per caller (SCHEMA Decision 8), each with `relation = caller_split` and its caller
  and version. A row goes to a caller when INFO/SF lists its index. In the one-column layout most PGP
  merges have, the VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the
  shared sites whose VarScan2 genotype the merge dropped. A `Loaded` can now hold several facets:
  `facets` has a row each, and `calls` and `blocks` run facet by facet in `facet_id` order. Not done:
  reading VarScan2's genotypes back from the `-orig` key=value table, which has no loader.
- `load` reads Complete Genomics' VCF (`vcf_cgi`, RM2). Its `<CGA_NOCALL>` spans and its ALT `.`
  records with no allele called go to `blocks` as `no_call` (to END, or over REF); its `<CGA_CNVWIN>`
  copy-number windows are left out with a count; half-calls keep their called allele. GL is read in
  phred units (gl = GL / 10) wherever the header names CGAPipeline: in every one of 126,604 rows with
  GQ in a PGP-HMS vcfBeta, GQ was the gap between the two best GL values. CGI's AD (`Number=2`, per GT
  slot) is not carried. The rules follow the declared alleles, so a CGI VCF that lost its `##center`
  line and sniffs as `vcf` loads the same way.
- A VCF's `ad`, `gl` and `gp` are now withheld where they cannot be indexed over the row's alleles:
  a row whose AD is not one value per allele, or whose PL, GL, GP or PP is not one per genotype of its
  ploidy, gets null there, with one warning per key and a count. A key declared with a fixed count
  (VarScan2's AD, `Number=1`, the ALT depth; CGI's AD, `Number=2`, reads per GT slot) is not carried
  and is listed in `dropped_fields`. Found loading the PGP-HMS `vcf_no_meta` file, where 141,814
  VarScan2 rows had put their ALT depth into `ad` as if it were per allele. Files whose vectors fit
  are unchanged.
- `load` reads `vcf_no_meta`, a VCF body under a bare `#CHROM` line (RM2). The loader synthesises
  `##fileformat` and declares the carried keys the rows use; `vcf_version` stays null. For every VCF,
  an undeclared INFO `END` or `SVLEN` is now declared too, so a block whose header forgot END keeps
  its extent instead of ending at POS, and `dropped_fields` also lists the keys the first 1,000 rows use
  without a declaration (polars-bio drops those unread, P1). Cached rewrites are made once more
  (`REWRITE_VERSION` 6). Dogfooded on the one PGP-HMS `vcf_no_meta` file (a freeBayes+VarScan2 merge
  stripped of its meta lines): 7,123,137 calls in 236 s over the SMB share; its INFO/SF split by
  caller is not applied (no `##sourceFiles` line says which index is which caller).
- `load` reads BCF (`bcf`, RM2). polars-bio reads the records; just-vcf reads the header out of the
  binary layout itself, and the first record's contig and position to check polars-bio's coordinates,
  as it does for text. polars-bio takes BCF only from a path named `.bcf`, so a BCF under another name
  or inside an archive is read from a copy in `JUST_VCF_CACHE_DIR`. Rows equal those of the same file
  as text VCF, except `source_gt`, which stays null: a BCF wrote no genotype text. Tested on BCFs
  written by pysam; no real BCF is on disk.
- A loaded VCF's facet now says who made it, from its header alone (RM2's `source` detector,
  `just_vcf/source.py`): `vendor` (the party the header names as delivering the file), `caller` and
  `caller_version` (the program that wrote the genotypes: GATK HaplotypeCaller or UnifiedGenotyper,
  DRAGEN, DeepVariant, freeBayes, VarScan2, CGAPipeline, Isaac/starling, CASAVA, Platypus, Manta,
  aengine, bcftools call, loimpute, Minimac), every row of RM2's header-evidence table plus
  `##bcftools_callVersion`. Tools that only touched the file (SelectVariants, ApplyVQSR, bcftools view)
  are not caller evidence. Two callers in one header leave `caller` null with one warning naming both;
  a `vcf-merge` of callers is not a conflict. Where the header says the content is array or imputed
  genotypes (Illumina GenotypingToVCF, the Open Humans and DNA.Land 23andMe conversions, Gencove,
  Minimac, Genes for Good), `callset_type` is `array` or `imputed`; nothing else sets it yet (RM4).
  `phasing_method` is filled from `##eagleVersion`. New facet column `ploidy_statements` (RM269 half
  1): the header's ploidy lines verbatim (`##referenceSexKaryotype=XY`, GATK `sample_ploidy=2`,
  `bcftools call --ploidy GRCh38`), with no vocabulary. Schema change: one new nullable column, minor.
  Checked against the 233 PGP-HMS pipeline VCF headers (CGI and freeBayes+VarScan2): every one names
  its caller, none conflicts.

- `resolve_ref(loaded, reference)` (RM3, charter 10) reads REF from the facet's assembly sequence,
  held in a gtars refget store (`ReferenceStore`; a FASTA, plain or gzipped, is added to the store at
  `JUST_VCF_REFERENCE_DIR` once). A VCF row gets `ref` over its `source_ref` span and `ref_check`
  `match` or `mismatch`; an array row gets `ref` from the reference, `alts` from its genotype's
  non-REF alleles (sorted) and `gt` into `[ref] + alts`. REF is never taken from the genotype. Before
  any base is read, each contig is matched with the reference by its stated length and md5, or, for
  arrays, by assembly label; a primary contig of another length refuses the facet
  (`ReferenceMismatchError`), and any other unmatched contig (hg19's chrM against GRCh37's MT) is left
  unread, counted in the warnings. I/D codes, no-calls and spans past a contig's end stay unread too.
  `chr`-prefixed and plain names meet through `canonical_chrom`. The facet keeps its id and gains
  `facets.ref_seqcol_digest` (schema: one new nullable column, version bump at merge). Minor.

- A sites-only VCF (no FORMAT column, no sample) is refused by `load` with the new `NoSampleError`
  (exported), instead of failing inside the parquet spill with a polars conversion error (F8). Minor:
  a new exception; nothing that loaded before changes, facet ids included.

- The kit's scripts are named by their new names since kit v0.1.86 (`bin/py.sh`, `bin/move-to-tag.py`,
  `bin/claim.sh`) in `CLAUDE.md`, the seat texts, the inbox headers and the tooling seat's issues watcher.
  The old names still work as symlinks in the kit. The pin is at v0.1.87.

- `rejoin_split(loaded)` (RM270 D, the person's choice on desk #45) re-joins the rows a pre-split VCF
  wrote one per ALT into one site each, as a derived facet (`relation = split_rejoined`; the parent is
  untouched). It works from the rows alone, with no FASTA and no reads, so it gives back the genotype
  (`gt`, `gt_alleles`) and withholds the likelihood vectors (`gl`, `gp`, `ad`, `gq`), which split rows
  cannot rebuild. ALTs come out sorted by allele string. A position it cannot read as one site stays
  as written, its rows still flagged, with a warning per reason. Tested against the unsplit form of the
  same synthetic sites, and on a real genome: its 89,855 multi-allelic records split the bcftools way
  (185,402 rows), loaded and re-joined, gave back every one of the 88,713 sites carrying an ALT,
  matching the unsplit file; the 1,142 split hom-ref or no-call sites stay as two rows, correctly
  unflagged. 23 s on that 6.4M-call genome; with nothing flagged it passes every row through. Not run:
  a file split by bcftools itself, or by Hail or GATK. RM270 moves to the roadmap history; RM6 keeps
  its rule that its own split output writes `.` for the other ALT.

- A filter on a loaded VCF's `calls` that had split-flagged rows ran about 20 times slower than a full
  collect (169 s against 8.5 s on a real genome split one row per ALT), because the flagged positions
  were joined below the per-contig sort, which kept the filter from reaching the parquet scan. They
  join above it now, and not at all when nothing is flagged. A file with nothing flagged was not
  affected (1.7 s against 1.5 s for the same filter).
- The triage commands name `$BLUEPRINTS_KIT`, which the seats kit's session-start hook exports from kit
  v0.1.78 (kit #96), instead of a `$KIT` each session had to find by hand.
- The triage runbook, `docs/CONSUMER_TRIAGE_LOOP.md`, now points at the seats kit's `LOOP.md` (v0.1.77)
  instead of the frozen gist, and keeps only the local parameters, the intake routes and the legality
  rules. The inbox header pins the next id again (S3), for a reporter who cannot run the kit's `--next`.
- The triage seat has a second watcher, the kit's `watch-remote.sh --path docs/FEEDBACK.md`, for a report
  committed and pushed from elsewhere, which the inbox watcher cannot see once the file matches HEAD.
- The triage loop runs the seats kit's scripts (`classes/triage/` from kit v0.1.77) instead of copies in
  `scripts/`. `triage-state.py`, `triage-archive.py`, `item-next.py` and `tests/test_item_next.py` are
  gone; the kit's versions gave the same ledger, the same next `S` id and the same next `RM` id on this
  clone's documents. The commands now read `sh "$KIT/bin/py" "$KIT/classes/triage/<script>"`, with
  `$KIT` the pinned kit copy (`CLAUDE.md` §8). The kit's `item-next.py` also counts an id mentioned in a
  changelog as taken.
- Schema 0.0.4 (RM270, from S2): a VCF row that looks like one record of a multi-allelic site written one
  row per ALT gets `calls.split_signal` (`marker`, `slot_pattern` or `all_zero_beside_alt`), the
  file gets `facets.split_marker` when vt's `OLD_MULTIALLELIC` or a bcftools `--old-rec-tag` key tagged
  a record, and each signal is one warning with its count. In such a row a `0` can mean another ALT,
  and `gt_alleles` still reads it as REF: no value changed. Whether to withhold those `0`s is on the
  person's desk (#45). Measured: no signal on four unsplit WGS VCFs (22.2 million calls) or on two
  unsplit gVCFs. Not run: a real split file (only synthetic ones). The patterns are read from the
  parquet spill, about 1.2 s on a 30 s load of 6.4 million calls; a second parse of the VCF had cost
  12 to 17 s. A cached copy is rewritten once (`REWRITE_VERSION` 5),
  since the copy made for a refused header now keeps the marker keys.

- The clone is on the seats kit v0.1.69, whose `watch-file.sh` asks the remote for its default branch
  (kit #88), so the inbox watcher's command no longer names `--branch main`.
- The inbox watcher is the kit's `drivers/git/watch-file.sh` (kit #76); `scripts/watch-inbox.sh` is
  gone. Its line points at the uncommitted diff rather than naming pending ids, which the pass's
  `triage-state.py --pending` lists anyway, and an item a sibling leaves while no session holds the
  triage seat is now said at the next punch-in, where the old watcher stayed quiet about any edit that
  predated it.
- The pinned polars-bio stack (0.36.0, with datafusion-bio-formats, the noodles and bigtools forks,
  datafusion-bio-functions and OpenDAL at their locked revisions) was audited once for code that sends
  data off the machine, and stamped "no network exfiltration code" in `docs/TRUST.md`. Every network
  request starts from a `scheme://` path the caller passed; just-vcf passes local paths only.
- Triage answered S2 (just-dna-lite: a VCF split one row per ALT reads wrong, since a split row's `0`
  can mean another ALT, not REF). Reproduced: the loader reads such a file the same way, with no
  warning. Filed as RM270, pre-split multi-allelic input; it also binds RM6's split output.
- A report can reach the triage seat as a message: from another session on this machine, or on Zulip
  (`just-vcf/triage`) from a person or a seated repo. The seat files it into `docs/FEEDBACK.md` verbatim
  and answers there. Appending to the inbox stays the main route. `docs/FEEDBACK.md` says how,
  and the runbook has the intake steps.
- Triage answered S1, the first consumer item (just-prs: male X/Y ploidy conventions differ by caller,
  vendor and panel, and nothing records which). Nothing had shipped for it; it is filed as RM269, a
  karyotype finding and the caller's ploidy model, beside RM2's planned source detector.
- `scripts/item-next.py` counts an id as taken only where a heading, a reservation or a tombstone row
  owns it, so a sibling's id quoted in prose no longer sets the maximum (just-dna-format's `RM267`,
  quoted in RM3, made it hand out RM268; RM268 keeps its number). A reservation now lands in its anchor
  section instead of the first bullet list below it. `tests/test_item_next.py` covers both.
- The triage loop no longer syncs from the gist: its algorithm is read at the frozen revision
  `96b69533` through that revision's URL until the seats kit's triage role replaces it.
- The inbox watcher runs under the triage seat (`[seat.triage.watch.inbox]` in `.blueprints.toml`),
  started and stopped with the seat by its watchdog; a settle with nothing pending is dropped, as before. The `UserPromptSubmit` hook, `.claude/settings.json`
  that held it, `scripts/arm-on-prompt.py` and `scripts/wake-on-inbox.sh` are gone.
- The repo runs on the seats kit (`agent-seat-blueprints`, plugin pinned at v0.1.60, local scope): three
  seats, triage, upstream-liaison and a new lead, are taken with `/seats:seat <seat>` and declared in
  `.blueprints.toml`, and `team.toml` lists the people. `docs/TRUST.md`, the trust register, is a
  threat model per source with exfiltration of a private genome as the one threat; `CLAUDE.md` §2
  states its two prohibitions. `docs/seats/README.md`
  drops the old priming line and announce round, which the kit now covers.
- The charter is committed as `docs/CHARTER_DRAFT.md`. While it is a draft it is edited directly, with
  no amendment procedure; it becomes `docs/CHARTER.md` when the user signs it off.
- A triage seat (`docs/seats/TRIAGE.md`) owns the consumer inbox, and `docs/seats/README.md` sets out how
  seats work here. The triage loop is synced to gist revision `96b69533`: `wake-on-inbox.sh` and
  `arm-on-prompt.py` join `scripts/`, and a `UserPromptSubmit` hook in `.claude/settings.json` arms the
  watcher when a prompt cites the seat file.
- The seats pattern and the UPPER_SNAKE doc-name rule went out to the onboarding gist
  (`9b04fb4c`, revision `12382b8e`, 2026-09-30): §3.5, a new §3.8, a triage-seat note in §4, and the
  template.
- Every markdown file under `docs/` is named in UPPER_SNAKE: `dogfooding.md`, `previous_issues.md`,
  `polars-bio-pending-fixes.md` and `seats/upstream-liaison.md` were renamed, and links followed.
- VCF `calls` and `blocks` are spilled to parquet in the cache (`JUST_VCF_CACHE_DIR/frames/`) and read
  back one contig at a time in canonical order, so a sorted load needs the memory of one contig, not of
  the file: 3.6 GB instead of 11.2 GB for a 6.4M-call VCF. `load` now makes that pass itself.
- The loader's cached copy of a VCF polars-bio cannot read in place (bz2, zip, tar, a repaired header) is
  written as BGZF instead of plain text: about a quarter of the size and as fast to read.
- Schema 0.0.3: `gq` is `u16`; Illumina gVCF blocks (`ALT .` with `INFO/END`) go to `blocks`; GP's
  scale is read from values that settle it, whatever the header version says. Found by loading the 124
  PGP-HMS VCFs: 1.5 billion rows, no `end_pos` before `pos`, no genotype index out of range. A cached
  rewrite from older loader code is no longer reused.
- `just_vcf.load` reads single-sample and multi-sample (one sample chosen) VCF and gVCF into `calls` and
  `blocks`: `<NON_REF>` written as `<*>`, records whose ALTs are all `<*>` and whose GT asserts reference
  or nothing go to `blocks`, `end_pos` follows VCF 4.5's END rule, `dp`/`ad` come from FORMAT by
  namespace, GP's scale is checked against the file version, and REF goes to `source_ref` only.
  Headers polars-bio would misread are repaired in a cached copy (`JUST_VCF_CACHE_DIR`). The module is
  `just_vcf.loading` (was `just_vcf.load`, which the function of that name shadowed). `vcf_samples(path)`
  lists a VCF's samples; several samples and no choice raise `MultiSampleError`.
- The VCF loader parses only the FORMAT keys and the sample it carries, and when polars-bio refuses a
  file (a raw backslash in a quoted header value, a malformed value in a field the frames do not carry)
  it reads it again with a header of only the carried keys, and says so in the warnings.
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
