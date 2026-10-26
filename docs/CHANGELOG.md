# Changelog

What shipped, newest first. Cross-repo integration changes made on our side are recorded here too.

## Unreleased

- `absence_support(facet, variant_class, force=False)` answers what `absence_is_reference` answers,
  with the gate that decided it (`AbsenceReason`): why a whole genome's silence reads as reference
  (`whole_genome_silence`, `omitted_called_reference`, `content_filter_forced`), or why an absence says
  nothing (`outside_blocks`, `array`, `targeted`, `imputed`, `callset_type_unknown`,
  `class_not_counted`, `class_not_called`, `content_filter_held`, `class_cut_by_filter`).
  `absence_is_reference` now returns its `reads_reference`; its answers do not change
  (docs/SCHEMA.md Decision 4).

- `concordance` compares structural-variant calls across two facets (RM282). A symbolic allele
  (`<DEL>`, `<DUP:TANDEM>`, `<INV>`, `<INS>`, `<CNV>`, `<CNn>`) pairs with at most one of the other
  facet's by type and extent, under truvari bench's matching and defaults (`SvMatching`: `refdist`
  500, `pctsize` 0.7, `pctovl` 0.0, `typeignore` off; truvari 5.5's `var_match` and `--pick single`).
  The type is the ALT's, not SVTYPE (DRAGEN writes `SVTYPE=CNV` on a `<DEL>`); CIPOS and CIEND widen
  the distance test, which truvari does only for breakends; an `<INS>` has no size in the frame, so
  its size test is withheld and counted, never passed. Pairs compare genotypes as sequence alleles do
  (`match`, `mismatch`, graded `p_same` from PL). A one-sided SV is `only_a`/`only_b` with `absence`
  `other_alleles` (an SV of the other facet within `refdist`), `reference` only where that facet is a
  whole-genome callset that called SVs and was not filtered, else `unknown`: an array, an exome or a
  small-variant callset never turns an SV into a disagreement. `sites` gains `sv_kind`, `b_alt`,
  `a_end_pos`, `b_end_pos`, `size_similarity`, `reciprocal_overlap`, `ci_needed`; `key` gains
  `extent`; `summary` groups by `sv_kind`. **Values change**: symbolic alleles were counted in
  `excluded` as `symbolic` and are now rows of `sites` (and of RM279's atlas, class `sv`); breakends
  are still left out, now as `breakend` (docs/SCHEMA.md, "Concordance between two facets").

- `resolve_ref`'s row sample for the build inference and the strand check (F9, F13;
  `assembly.sample_snv_rows`) keeps only rows over one base (`end_pos == pos`). It sampled a Genos
  insertion site, `G/G` inserted after a `G` at the anchor with `end_pos` null, as an SNV at the anchor
  base, and would have sampled a Genos or GET-Evidence row called over two bases whose alleles are
  single letters. Every other array format already wrote `end_pos = pos` exactly where the sample
  took a row, and CGI's and CASAVA's insertion points have a source REF, so only those two change.
  `RESOLVE_VERSION` 12: the strand is not in the spill key, so a cached result could hold rows a
  mixed strand withheld under the old sample. **Values change** where such rows tipped a facet's
  inferred `assembly_label` or `strand_check` (docs/SCHEMA.md, "Changes since sign-off").

- A VCF or a coverage facet whose header lengths are NCBI36's (hg18) gets `assembly_label = GRCh36`,
  basis `inferred` (RM4). `assembly.ASSEMBLY_LENGTHS` held GRCh37, GRCh38 and T2T-CHM13 only, so a
  build-36 header's lengths named no assembly and the label stayed null. The lengths come from NCBI's
  assembly report for GCF_000001405.12, MT from hg18's chrom.sizes (the report lists none), both
  carried in `assets/` (docs/PROVENANCE.md). No primary length is shared with another assembly, and
  the one shared mitochondrion length never supports a label, so a header with only chrM or two
  chromosomes stays unknown. **Values change**: null to `GRCh36` for those facets
  (docs/SCHEMA.md, "Changes since sign-off").

- `liftover` lifts symbolic SVs (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<CNn>`, `<INS>`) and copy-number
  segments by both ends (RM6). POS and `end_pos` are lifted as two points; the row is kept where both
  land on one contig and strand, in order, with a length within 5% of the source's
  (`liftover.SV_LENGTH_TOLERANCE`, liftoverSV's rule), and otherwise left behind under the new reasons
  `sv_end_unlifted`, `sv_ends_apart`, `sv_ends_reversed` and `sv_length_changed`. REF is read from the
  target at the new POS; on a minus strand the SV is anchored on the target's base before the event's
  image and `cipos`/`ciend` swap and mirror. Breakends stay behind as `breakend_not_lifted`, where they
  were `allele_not_sequence`. **Fixed**: a copy-number segment (`ALT .` with END and CN) lifted with
  its `end_pos` cut to its one-base REF span. **Values change**: lifted facets get a new `facet_id`
  (`LIFT_VERSION` 9). Additive otherwise: minor (docs/SCHEMA.md, "Lifting symbolic SVs").

- `liftover` lifts the `blocks` frame (RM6). A gVCF reference or no-call block is cut wherever the
  chain's mapping changes: the pieces inside one chain block with one target lift, with `kind`, `gq`,
  `dp_min` and `gl` unchanged, and the rest stays behind, counted in bases per reason in the new
  `Lifted.block_reasons` and in the warnings. A block crossing a chain gap becomes several blocks,
  never one stretched over the gap; minus-strand and other-contig pieces lift as rows do. The frame
  keeps no block REF, so `write_vcf(..., blocks=True)` on the lifted facet reads it from the target.
  The lifted facet had an empty `blocks` frame; its `absence_means` stays `unknown`. **Values change**:
  lifted facets get a new `facet_id` (`LIFT_VERSION` 8). Additive otherwise: minor (docs/SCHEMA.md,
  "Lifting blocks").

- A facet records the rule that read its REF, and the store refuses a held facet given again in
  another REF state (schema 0.0.20). `resolve_ref` keeps a facet's `facet_id`, so `Store.add` of the
  facet resolved under a newer `RESOLVE_VERSION`, or of the loaded facet once its resolved twin was
  stored (or the reverse), reported it `present` and kept the held rows. `facets.resolve_version` is
  new (the `RESOLVE_VERSION`, null until REF is read, inherited by derived facets), and `Store.add`
  compares it and `ref_seqcol_digest` with the held facet's, raising `StoreConflictError` (new, a
  `StoreError`) before anything is written. A held facet is never rewritten: the other state goes into
  a new store. `migrate_store` adds the column to a 0.0.19 store, null and named as not recorded on
  facets whose REF was read (docs/SCHEMA.md, "The durable store" and "Store migration").

- `alignment_coverage` names a BAM or CRAM with no index beside it by the sha256 of its bytes, now in
  `CoverageFacet.source_sha256`, its `facet_id` and its spill key. Its header and byte size alone named
  it, so a file rewritten in place with the same header and size but other reads was served the old
  file's windows under the old `facet_id`. The cost is one more read of an unindexed file per call. An
  indexed file is named as before. **Values change** for an unindexed file: its `facet_id`, and
  `source_sha256` from null (docs/SCHEMA.md, "Coverage facets").

- `migrate_store` called after a migration that was cut off between committing the new manifest and
  removing the old version's files returns that migration, read from the manifest's `migrations`,
  and finishes the cleanup. It returned `from_version = to_version`, `migrated` null, which reads as
  nothing applied though the store had just been migrated (docs/SCHEMA.md, "Store migration").

- `check_pedigree` gives a pair whose kinship lies wholly in the first degree, but whose IBS0 count
  cannot tell parent from sibling where the pedigree expects one, the reason `kind_not_told`
  (`person.KIND_NOT_TOLD`). It read `between_bands`, which says the kinship's interval crosses a band
  boundary: kinship 0.25 ± 0.01 with no IBS0 over 300 sites against a parent link was reported as
  between bands. The reasons a check row may carry are now one vocabulary,
  `person.PEDIGREE_CHECK_REASONS`, checked where the row is built. **Values change** for a reader of
  `reason`: those rows read `kind_not_told`; `agreement` is unchanged (docs/SCHEMA.md, "Relatedness
  and pedigrees").

- `loaders.text.GT_INDEX_MAX` is an int (255), read from `gt`'s dtype. It was `pl.UInt8.max()`, an
  expression annotated `int`; the warning for an index `gt` cannot hold printed "(255)" only because
  polars renders that literal bare. No value changes.

- A store written under an older schema version can be migrated (RM276): `migrate_store(path)`
  rewrites it to this version by one step per schema bump since the store shipped under 0.0.14
  (`just_vcf.store_migration.MIGRATIONS`, five steps to 0.0.19), and returns a `StoreMigration`
  naming the steps, the facets, and every value it filled unknown (`FilledUnknown`, with the facet
  ids and the reason: a migrated VCF facet's `method_steps` and `gl_source` are null for not read,
  a CGI facet's `hq` and hapLink phase likewise; load the source again to read them). Each facet is
  written beside the old files and checked against them before a marked swap under the store's lock;
  a migration cut off before the mark leaves the store as it was, one cut off after finishes on the
  next call. The old version's files are removed once the new manifest is committed, or kept with
  `keep_previous=True`. Nothing migrates on read: `loaded()` and `add` on an older store are still
  refused, now naming `migrate_store`; a store from a later just-vcf is refused with "upgrade
  just-vcf". The manifest gains `migrations` (default empty, `store_layout` unchanged), and
  `loaded()` repeats each filled-unknown set as a warning (docs/SCHEMA.md, "Store migration").

- `resolve_ref` reads the anchor base of an insertion point and REF over a Genos row's span (RM3,
  `RESOLVE_VERSION` 10). An insertion point (CGI `var` and `masterVarBeta` loci of zero length, a
  CASAVA indel table's `-/T`, whose facet states no build and is still refused, a Genos
  `BEGIN == END` row) sits at the base before it; that base, read
  from the reference, is now `ref` and is written before every allele (`["", "GG"]` becomes
  `["G", "GGG"]`), `end_pos = pos`, `ref_check` null (the file wrote no REF base). Before a contig's
  first base, VCF's rule: the base after, appended, at POS 1. A point where nothing was called stays
  unread, counted, and so does one whose every called haplotype is the empty allele (CGI `""/""` and
  `""/?`, Genos `-/-`; `RESOLVE_VERSION` 11): it says no insertion is there, and anchoring it wrote
  `REF=X ALT=. 0/0`, a hom-ref call at the base before that the source never made, beside the
  source's own call there. It keeps its loaded `["", ""]`, counted under its own warning; Genos's is
  no longer counted among the array rows too. A Genos row gets REF over its stated span whatever its length, and its alleles
  derived against it as an array row's (`AT/GC`, `T/-`); before, only its single-base rows were read,
  and a one-base insertion site was read as a call at the anchor base: `G/G` inserted after a `G`
  read hom-ref there (checked on the old code), another base an SNV. `normalize` aligns the anchored
  rows (`NORMALIZE_VERSION` 9), and anchors a deletion at a contig's first base on the base after
  its span, where it left it unanchored. `genotype_evidence` names an anchored point's alleles with
  `ref`. **Breaking for a reader of a resolved CGI or CASAVA facet**: an insertion point's `alts`,
  `gt_alleles` and `source_alleles` change spelling, and `ref_check` goes from `not_run` to null
  (docs/SCHEMA.md, "How RM3 fills these slots" and "Changes since sign-off"). On PGP-HMS hu1B2689's
  `var`, 263,689 of 263,690 anchored insertions with an ALT equal CGI's own vcfBeta record at the same
  POS, REF and ALT set alike (RM3 in docs/ROADMAP.md).

- `write_vcf_zarr` writes `cipos`, `ciend` and `alt_copy_numbers` (RM276), which it counted as not
  written: `variant_CIPOS` and `variant_CIEND` (`<i4`, two offsets per ALT, on `CIPOS_dim` and
  `CIEND_dim`) and `variant_CN` (`<f4` on `alt_alleles`), each where any row has one, read back by
  polars-bio's `scan_vcf_zarr`. The warning "fields VCF Zarr export does not write" is gone. An
  interval holding an offset of -1 or -2, VCF Zarr's integer missing and fill values, is written
  missing on that row with a warning, since a reader would take the value for no value
  (docs/SCHEMA.md, "VCF Zarr export").

- `write_vcf` writes INFO `SVLEN` beside END on a length-type SV allele (`<DEL>`, `<DUP>`, `<INV>`,
  `<CNV>`, `<CNn>`; RM6), which VCF 4.5 requires and tools that read SVLEN rather than END need. The
  frame keeps no SVLEN, so the value is `end_pos - pos`, END restated: `.` on an `<INS>`, a breakend or
  a sequence ALT; VCF 4.3's `Number=.`, negative for a deletion; 4.5's `Number=A`, unsigned (a facet
  with a PSL set). A record with several length-type ALTs gets it only where VCF 4.5 says they share
  one length (copy-number alleles, or `<DEL>`/`<DUP>`/`<CNV>` beside FORMAT CN); elsewhere it is left
  out with a warning. A written file's reload reads the same `end_pos` (docs/SCHEMA.md, "Structural
  variants").

- A `Person` carries its coverage facets (RM278): `read_coverage` and `alignment_coverage` take
  `subject=`, written on the new `CoverageFacet.subject` (not in `facet_id`), and `Person.of` and
  `people` take `coverage=[...]`, grouped by that label and refused (`SubjectError`) without one or
  under another subject. `person_sex_chromosome_evidence(person, pairs=None)` is new: each facet's
  estimate with the person's coverage facet on its build, and a new `coverage_pairing` field saying
  which was used and why (`only_one_on_build`, `only_usable_on_build`, `stated`, …). Of two usable
  coverage facets on one build none is used and both are named, since nothing records which
  alignment a callset came from; `pairs` states one. `sex_chromosome_evidence` is unchanged and
  leaves `coverage_pairing` null. Additive.

- `alignment_coverage` reads a BAM or CRAM once (RM278): its windows are spilled to the cache under
  `coverage/` (a new entry kind, `coverage_spill`, with a manifest, the stale-version sweep and
  `clear_cache`, as RM268's), and a second call with the same file, window and read filters reads the
  spill instead of running polars-bio's `depth()` again. The key is the file's identity as `facet_id`
  takes it (header, index, size) with the parameters, `COVERAGE_VERSION` and the polars-bio version, so
  no spill made with other parameters or by another version is served; a stated `assembly_label` reuses
  the spill. The returned facet is the same, its Int16 warning included; `intervals` is now a scan of
  the spill's parquet. `just_vcf.coverage.POLARS_BIO_VERSION` is new.

- **Breaking, pre-1.0:** relatedness splits the first degree into `parent_child` and `full_siblings`
  by IBS0 (RM274), two new `RelationBand` members. A parent and child are never homozygous for
  opposite alleles absent genotyping error; full siblings can be, over the quarter of the genome they
  share nothing of by descent. The
  rule is KING's (Manichaikul et al. 2010, section 2.3) and the cutoffs KING 2.3.2's fixed defaults
  (IBS0 at most 0.005 parent and offspring, above 0.008 full siblings); a band is named only when the
  IBS0 share's Wilson interval (two SEs) lies wholly on one side, so KING's grey zone and too few
  sites stay `first_degree`, which now means the kind is not told. A consumer selecting
  `band == "first_degree"` gets fewer rows: the value narrowed, it was not added beside (docs/SCHEMA.md,
  "Changes since sign-off"). somalier's pairs are split the same way. `expected_kinship` expects
  `parent_child` for a listed parent link and `full_siblings` for two shared listed parents or a stated
  sibling term, and `check_pedigree` reads a parent link measured as siblings (or the reverse) as
  `disagrees`, where it read `agrees`. `same_person` is unchanged. KING's cutoffs assume array-grade
  genotypes: a parent and child read from facets as noisy as PGP's merges can read `full_siblings`
  (docs/RELATEDNESS.md, Limits).

- `normalize` no longer fails on a contig whose one multi-allelic site to split is an ALT beside `<*>`
  (`2 60 . A C,<*> 0/1`), which DeepVariant writes on every variant record, so a contig with one variant
  met it: the split's pieces were a frame of one row, where polars 1.44's `list.slice` at a column
  offset inside `when/then` raised `ComputeError` ("could not extract number from any-value of dtype:
  'Null'") on the `cipos`/`ciend` slice. Those two values are read with `list.get` now. The output of
  every input that normalized before is unchanged (`NORMALIZE_VERSION` stays 8).

- The VCF, key=value and Nirvana loaders count a genotype whose allele index the row has but `gt`
  cannot hold (`0/300` over 300 ALTs: `gt` is a `u8` list), with its own warning (`… genotype(s) hold
  an allele index past what gt can hold (255), though the row has that allele; gt and gt_alleles have
  a null there, as an uncalled slot does`). It read `gt [0, null]`, the frame of the half-call `0/.`,
  counted nowhere, since it is not past the ALTs. A genotype past the ALTs is counted as that alone.
  The frame does not change; widening `gt` to `u16`, a retype, is left to the person.

- A VCF 4.3+ refused because a key the version reserves as a Flag (INFO DB, SOMATIC) is written with a
  value (`DB=1`) now names the key: `HeaderRefused` says `INFO/DB written with a value: VCF 4.3 reserves
  DB as a Flag…`, with `key` set, where polars-bio's error ("invalid flag") named none. Where no record
  read ahead writes one, it says so and lists the version's Flags. The file is still refused, also when
  its own line declares the key otherwise (`Number=1,Type=Integer`), which the docs had said it was not:
  the reader refuses that line and the copy leaves it out, so the key is read as its Flag.

- A VCF 4.3+ key its version reserves and the loader never carries (FORMAT MQ, NQ, INFO DP, AN, and
  the rest of the reader's table) is skipped unread in the retry whatever its values, declared in the
  file or not: the structural copy left its line out, polars-bio's reader parsed it under the
  version's own definition, and FORMAT MQ `60.5` or INFO DP `30.5` under 4.3 failed the load with
  `HeaderRefused`, naming no key. The copy now declares every reserved key it does not declare
  otherwise `Number=.` with the reserved Type, GT and Flags excepted (`REWRITE_VERSION` 17). A Flag
  written with a value (`DB=1`) still fails the read.

- A VCF whose INFO holds a value polars-bio cannot parse reaches the loader's retry, as a FORMAT one
  does, in every version: INFO DP `Type=Integer` holding `30.5` beside a carried SVTYPE or a split
  marker failed the load with polars' raw `ComputeError` from the spill, since the pass the retry wraps
  read no INFO column. It now reads every INFO column the spill reads, and the file loads with DP not
  carried, under VCF 4.2 as under 4.3+ (the entry above).

- The VCF loader counts a genotype indexing past the row's ALT list, with the warning the key=value
  table and Nirvana JSON give (`… genotype(s) index an allele past the row's ALT list; gt_alleles has a
  null there`), counted from the GT text's digits: `0/300` loaded as `gt [0, null]` and `0/2` over one
  ALT as `gt [0, 2]`, both unwarned. The frame does not change.

- The key=value table and Nirvana JSON withhold a genotype that is not allele indices whole: `gt`,
  `gt_alleles` and `phased` are null and `source_gt` keeps the text, where the readable slots were read
  (`0/x` was `gt [0, null]`, the frame of the half-call `0/.`, and `N/A` was `[null, null]`, the frame
  of a written no-call). The warning says so (`… kept with gt, gt_alleles and phased null, and the
  text in source_gt`). A consumer reading these facets sees `gt` null on such rows; no column changes.

- A VCF 4.3+ key its version refuses under the declared Type (DP, GQ or PS `Type=Float`, CN `Float`
  under 4.3, INFO END `Float`) is skipped unread whatever its values, as Decision 9 says: `30` loaded
  with the key not carried, but `30.5` failed the whole load with `HeaderRefused` ("invalid digit found
  in string"), because polars-bio's reader parsed the key, left undeclared in the copy, under the
  version's own definition. The structural copy now declares such a key `Number=.` with the reserved
  Type, which the reader never parses for a key not requested (`REWRITE_VERSION` 16).

- Nirvana JSON counts a genotype indexing past the row's ALT list from the genotype text, as the
  key=value table does: `1/300` reads `gt [1, null]` (a `u8` list) and was counted from `gt`, so it
  went unwarned. Each genotype is counted under one reason there too (`1:2` is unreadable only).

- The key=value table counts each genotype under one warning reason: one not readable as allele
  indices (`1:2`) was counted there and again as indexing past the row's ALT list, from its digits; a
  readable genotype alone is checked against the ALTs now. A genotype with an index too long for 64
  bits (`0/99999999999999999999`) is counted as past the ALTs, where it failed the load with polars'
  `InvalidOperationError`.

- A BCF whose header declares GT with a Type other than String is refused with `HeaderRefused` naming
  GT and its Type, where polars-bio's own refusal named nothing under 4.3 ("invalid record"). The text
  VCF's warning that such a GT is read through a copy declaring it String no longer applies to a BCF,
  which has no copy.

- The tooling seat signs what it posts on this repository's issues with `Seat: just-vcf/tooling`, so its
  issues watcher no longer wakes it for its own comments and closes (kit v0.1.94, #165). The pin is at
  v0.1.94.

- The key=value table's and Nirvana JSON's warning for a genotype that is not allele indices says what
  the row holds: a null in `gt` and `gt_alleles` for each slot that is not an index, and the text in
  `source_gt`. It said "gt and gt_alleles null", but `gt` is a list of one slot per haplotype the text
  separates (`N/A` reads `[null, null]`), and a slot that is an index is read (`0/x` reads `[0, null]`).
  Wording only. **Superseded within this release** by the entry above that withholds such a genotype
  whole: `N/A` and `0/x` both read `gt` null now, with the text in `source_gt`, and the warning says so.
  This entry records the wording in between; neither frame it describes ships.

- A VCF read through a copy that declares its `Float` integer keys `Type=String` (`retyped_whole`)
  says so: one load warning names each key so read and its declared Type (`FORMAT/DP (Type=Float)`),
  as the GT retype's warning does. A key not carried is named by its own warning instead, and a BCF,
  read as stored, gets none. Additive.

- An array facet on GRCh36 is estimated (RM269), where it answered `uncalibrated_build` (and before
  GRCh36 had a region set, `unsupported_build`): `karyotype.ARRAY_ESTIMATE_BUILDS` adds GRCh36 to the
  builds the array rule reads, since its bands are not measured thresholds and the build only places
  the regions. A whole genome on GRCh36 still answers `uncalibrated_build`. Checked on PGP-HMS: the
  250 GRCh36 23andMe exports read as the GRCh37 ones do (73 two-letter X, all self-reported women,
  `XX` with `array_y`; 173 one-letter X, the vendor's call), and of the 69 tables `resolve_ref` places
  on GRCh36, 3 give one X and none reads against its label. A value change; no column changes.

- `sex_chromosome_evidence` places the sex chromosomes on GRCh36 (NCBI36, hg18), where it answered
  `unsupported_build` (RM269): `PAR_BOUNDS["GRCh36"]` are NCBI36's PARs as Ensembl 54 records them
  (checked on hg18's own sequence, X and Y identical inside each PAR), and `X_DEGENERATE_Y["GRCh36"]`
  is the GRCh37 set lifted through UCSC's hg19ToHg18 chain, kept between those PARs (docs/PROVENANCE.md;
  the sources are carried in `assets/` and the test re-derives both). A GRCh36 facet is grouped and
  measured, and answers `uncalibrated_build`, as T2T-CHM13 does: no NCBI36 genome was among the
  calibration files. `unsupported_build` is left for labels with no region set (NCBI35, a label a
  caller wrote). A coverage facet on GRCh36 or T2T-CHM13 is now set aside as `uncalibrated_build`,
  where it said `unsupported_build` although both have region sets. A value change; no column changes.

- `sex_chromosome_evidence` on an array facet tells a written no-call from an unreadable row (RM269):
  the check that the export writes its no-calls counted every autosomal row not called, so one row
  whose genotype no loader could read (`gt_alleles` null) let an export that drops its no-calls read
  its Y line with `array_y=True` (700 called Y rows of 700 read as a Y). It now counts the rows written
  as no-calls (`--`, a letter null), on `ArrayLines`' new `autosome_no_calls`, and both the Y ratio
  and `no_calls_not_written` read that count. Additive field.

- An integer key declared `Float` (GQ, DP, MIN_DP, LEN, PS, AD, PSO, PSQ, HQ, INFO END, SVLEN, CIPOS,
  CIEND) is read from its text, where polars-bio's 32-bit float rounded the fraction away and the value
  read as whole: DP `30.0000001` read 30, END `5000000.1` read 5000000, PS `10000000.5` read set
  `10000000`. The copy polars-bio reads declares such a key `Type=String` (`REWRITE_VERSION` 15), and
  whether a value is whole is decided by its digits, never by a float, for `String` keys too
  (`30.00000000000000001` was 30 there): a fraction is left out and counted. No magnitude threshold
  could do it, since `30.0000001` rounds well below 2^22 (docs/SCHEMA.md, Decision 9). A whole number
  past 2^24 now reads exactly (PS or END 123456789, which was left out as inexact). Exponent forms
  read by their digits as before (`1e2` is 100, `1.25e1` is left out). A BCF's Float is its stored
  32-bit value and reads as before.

- A VCF that declares PS but whose PS the loader does not carry (a Type its version refuses, such as
  `PS Type=String` under 4.3; a PS the retry left out) withholds the phase of every phased record,
  with one warning counting them, where `0|1:100` and `1|0:200` read as phased in one set across the
  chromosome, the implicit set of records that wrote no PS: the values are skipped unread, so which
  records wrote one is unknown. A record whose PSL names a set keeps it, here and where a PS value is
  left out, since PS is ignored there (it was withheld). A file that declares no PS is unchanged.

- A negative PS (`0|1:-5`, `0|1:-1`) withholds the record's phase (`phased` null where it was true),
  with a warning counting the values (`… PS value(s) negative …`), where the record landed in its
  chromosome's implicit phase set: VCF's PS is a non-negative integer, so the value names no set.
  PS `0` reads as before.

- A Nirvana JSON file whose sample entries carry no genotype loads with `gt` and `gt_alleles` null,
  where it failed the spill with polars' `InvalidOperationError`, the VCF loader's failure below; and
  `normalize` splits a contig whose every record has no GT, where its `<*>` check raised
  `InvalidOperationError` on the same Null-typed gather. Every `list.gather` indexed by a column that
  can be null is now guarded on that column, also in the CASAVA SNP (`ad`) and Color (`gt_alleles`)
  readers. Those two build eager frames, which carry the Null type in their own schema, and the
  frame's cast to the calls schema absorbed it: no input fails through `load` without their guard,
  and their sink tests, which pass on the unguarded readers too, are regression shape tests only. What
  the guard changes there is the reader's own frame, typed instead of `List(Null)`, which two tests
  check on `read_site_table`. The key=value table and CASAVA indels never reach one: a row without a
  genotype is dropped before `gt` is read.

- The forge gate keeps checking every session's GitHub writes (issues, comments, pull request bodies)
  against `.leak_words`, not only a seated session's: `[repo] forge_gate_unseated = true` in
  `.blueprints.toml`, since kit v0.1.93 (#162) gates unseated sessions only when asked. The pin is at
  v0.1.93.

- A VCF 4.3, 4.4 or 4.5 file declaring GT `Type=Integer`, `Float` or `Character` loads with its
  genotypes, where the retry after polars-bio's refusal left GT out and the load failed with
  `InvalidOperationError` (or, past the fix below, would have had no genotype on any row). GT is never
  left out: the copy polars-bio reads declares it `Type=String` and the rows are read as written, in
  every version, with a warning naming the declared Type. A 4.2 file, which loaded already, now gets
  that warning too. `REWRITE_VERSION` 14.

- A VCF whose records carry no GT (a FORMAT of `DP` alone, whether GT is declared or not) loads with
  `gt` and `gt_alleles` null, the genotype unknown, where a run of such records failed the spill with
  polars' `InvalidOperationError` (`list.gather` over a batch of all-null indices makes a list of a
  Null inner type, which the parquet sink cannot write). A file that declares no GT and whose records
  read ahead use none says so in a warning (`no GT: …`).

- `SexChromosomeEvidence` says where its build came from (RM269): `assembly_label_basis` (the
  facet's column), `assembly_label_source` (`stated`, `contig_lengths`, `rows` for a build
  `resolve_ref` inferred from an array's rows, `liftover`; null for a label no rule accounts for) and
  `ref_seqcol_digest`. All three are new nullable fields, carried and deciding nothing. The estimate
  already read a build `resolve_ref` inferred: an array table that states no build (FTDNA, an rsid
  table, a Final Report) answers `no_build` as loaded and is estimated on `resolve_ref`'s output, and
  `no_build`'s detail now says so. Nothing infers or defaults a build in the estimate.

- A VCF 4.2 file declaring FT, PSL or INFO SVTYPE with a Type other than String (`Integer`, `Float`,
  `Character`) loads, where the first value (`PASS`, `DEL`) failed the whole read: the key is not
  carried, with a warning that names it and its Type, and it is in `dropped_fields`, as a 4.3+ header
  declaring them so already was. `REWRITE_VERSION` 13.

- A PSO or PSQ declared `Number=1` loads, where a diploid row's two values failed the load with
  `HeaderRefused`: VCF defines PSL, PSO and PSQ as one value per GT allele, and the copy polars-bio
  reads now declares a `Number=1` one `Number=.`, as it does `Number=P`, with each row's count
  checked against GT and PSL. A PSL declared `Number=1` reads its values too, where its records' phase
  was withheld. `Number=2` reads as before. `REWRITE_VERSION` 12.

- A carried integer key declared `Float` or `String` (a VCF 4.2 header, or MIN_DP, LEN, PSO or PSQ in
  a version that holds no definition of them) keeps only whole numbers, where it was truncated or
  lost without a word: GQ `30.7` read as 30, AD `3.5,4` as `[3, 4]`, END `200.5` as 200, and GQ `x`
  as null. A whole number reads as it is (`30.0` is 30); a fraction, text that is no number, and a
  whole number a 32-bit float cannot hold exactly are left out, one warning per key and reason with
  its count. The last mattered: polars-bio reads a Float as a 32-bit float, so a PS or END of
  123456789 came back 123456792. A PS left out withholds a phased record's phase rather than putting
  it in the chromosome's implicit set. A PS, LEN or INFO END declared `String` with text that is no
  integer no longer fails the load with `InvalidOperationError`. Covers GQ, DP, MIN_DP, LEN, PS, AD,
  PSO, PSQ and INFO END, SVLEN, CIPOS, CIEND.

- A VCF with a carried key declared as VCF reserves it but holding a value polars-bio cannot parse
  under that declaration (HQ `Integer` with `30.5,40` or `x,40`, AD with `x,4`, INFO SVLEN with `x`)
  loads without that key, where it failed with `HeaderRefused`: when the structural retry is refused
  too and polars-bio's error names the key, the loader reads again with it left out, warns once
  naming it and polars-bio's reason, and lists it in `dropped_fields`. Each key is left out once at
  most. A single-value key (GQ `30.5` under `Integer`) still fails the load: polars-bio's error names
  no key there, and the loader does not guess.

- `sex_chromosome_evidence` reads array facets (RM269), which answered `not_wgs` before. From the
  letters as written (no REF needed): X copies from the het share of two-letter base calls on non-PAR
  X over the autosomes' (one X at most 0.1875, two or more from 0.5625, the central half of the gap
  between 0 and the 3/4 that two X carry under the neutral expectation), and a Y line from the
  X-degenerate Y probes' called share over the autosomes' (absent at most 0.25, present from 0.75).
  An X written with one letter is the vendor's call and answers `caller_asserted_only`. The Y line is
  read only with the new `array_y=True`, since a vendor may no-call Y for every sample it calls female;
  without it `y_reason` is the new `vendor_asserted` and the karyotype unknown. An export that writes
  no no-call on the autosomes has no Y line (`no_calls_not_written`): its rows are not its probes.
  The counts are on the result's new nullable `array` field (`just_vcf.karyotype.ArrayLines`).
  PGP-HMS check in docs/ROADMAP.md RM269. Minor.

- `sex_chromosome_evidence(loaded, coverage=…)` counts copies past two from the coverage facet
  (RM269): each copy reads half chr20's depth, so X reads one to three copies and Y none to two, each
  0.125 either side of its value, and names `XXX`, `XYY`, `XXYY`, `XXXY`, and `OTHER_KARYOTYPE` for
  three X with two Y. New nullable result fields `x_copy_count` and `y_copy_count`, and a reason
  `copies_uncounted`. Values that change: a coverage X ratio of 1.125–1.25 (was two X), one under
  0.375 (was one X) and a coverage Y ratio of 0.625–0.875 or over 1.125 (was a Y present) now answer
  unknown, `intermediate`; the calls' X depth above 1.25 allows three or more X where it read unknown,
  so `x_copies` is `two_or_more` there and the karyotype `copies_uncounted` unless coverage counts it.
  No XXX or XYY genome was measured: the bands are copy number alone. From calls alone `XY` and `XXY`
  are still named, though they cannot tell XY from XYY. Minor.

- docs/CONSENSUS.md measures why CGI's 10–20 phred band disagrees with freeBayes's capped top band on
  13% (hu1B2689) and 52% (huAF3C63) and with its 30–60 rows on about 1%: each side's band says which
  genotype it called (CGI's 10–20 SNVs are all VQLOW hom-ALT calls, freeBayes's 60–100 rows mostly
  hets), so the cell is mostly a genotype disagreement, and the fitted band rate is still the band's
  average error. No rule or default changed. `scripts/band_cells.py` counts the cells by genotype,
  VQLOW, region, depth and an optional third facet, aggregates only.

- A CGI `var` or masterVarBeta score that is not a whole number is left out of `calls.hq` and
  counted, where it was lost without a word: the score was cast to an integer leniently, so `30.0`
  and `30.5` alike became no score, silently. Now a whole number written as a float reads as it is
  (`30.0` is 30), a number with a fraction is left out (never truncated), text that is no number is
  left out, each reason with its own warning and count: the rule a VCF's HQ follows (below). A `?`
  slot's score is still not read, so not counted.

- The CGI `var` and masterVarBeta warning that counts haplotypes whose varScoreVAF differs from the
  varScoreEAF read into `hq` counts only haplotypes whose EAF is read. It counted a half-call's `?`
  slot, a partly called haplotype and a score past u16 too, whose EAF never reaches `hq`.

- An HQ value that is not a whole number is left out of `calls.hq` and counted, where it was
  truncated or lost without a word: a VCF 4.2 header declaring HQ `Type=Float` read `30.5,40.5` as
  `[30, 40]`, and one declaring it `String` read `30,x` as `[30, null]` silently. Now a number with
  a fraction is left out (never truncated or rounded), text that is no number is left out, and each
  reason has its own warning with its count. A whole number written as a float (`30.0`) reads as 30.

- A VCF whose HQ the loader does not carry loads again: HQ declared with another Number and a row
  writing more values than it allows (`Number=1` with `30,40`), or declared `Type=Float` under VCF
  4.3+, failed with `HeaderRefused` since HQ became carried, where it loaded before. The structural
  copy the loader retries on now leaves out every carried key declared with a Number it is not read
  under, which covers AD, PL, GL, GP and PP declared `Number=1` too. And a carried FORMAT or INFO key
  declared with a Type its VCF 4.3, 4.4 or 4.5 version reserves otherwise (GQ, DP, PS, END, SVLEN
  and the rest) is not carried, with a warning and the key in `dropped_fields`, as 4.3's CN already
  was; those files failed to load before. `REWRITE_VERSION` 11.

- `write_vcf_zarr` writes `hq` as `call_HQ`, shaped (variants, samples, ploidy) like `call_genotype`:
  `-1` on a slot with no score, `-2` past a row's ploidy (a haploid `[51]` is stored `51, -2`), and a
  row of three or more alleles is written whole, which `write_vcf`'s `Number=2` HQ cannot hold.
  polars-bio's `scan_vcf_zarr` reads it back as the `HQ` genotype field. An `hq` that is not one value
  per GT slot is left out, counted in the warning `write_vcf` gives. Before, the export dropped `hq`.

- `normalize` no longer leaves a split row's `hq` whole: a GT slot it writes `.` (the slot held
  another ALT) loses its haplotype quality too, and a row left with no score has `hq` null, so `1/2`
  with HQ `51,40` splits into `[51, null]` and `[null, 40]`. Before, both rows kept `[51, 40]`, and
  `write_vcf` wrote a score the loader then refused on the uncalled allele. `rejoin_split` now joins
  `hq` per GT slot, as it joins GT (one score per slot from the rows, null where they give two or
  none, or where the slot is uncalled), so a normalized site re-joins with its parent's scores, and a
  vt-split file's HQ joins where it was dropped. Grades are unchanged (a split row with a `.` slot is
  not graded by `hq`). `NORMALIZE_VERSION` 8, `REJOIN_VERSION` 5: derived `facet_id`s change.

- `write_vcf` writes `hq` as FORMAT `HQ` (`Number=2`, as VCF reserves it), so a CGI `var` or
  masterVarBeta facet keeps its allele scores in its VCF export, graded alike on both sides, and a
  VCF with HQ round-trips. A haploid row is written `51,.`, which `load` reads back as `[51]`; a row
  with more than two GT alleles cannot be written under `Number=2` and has its `hq` left out, counted
  in a warning.

- The VCF loaders read FORMAT `HQ` (haplotype qualities) into `calls.hq`, one value per GT allele, so
  a VCF row with HQ and no GP, GL or GQ that grades it is graded by its lowest HQ (basis `hq`), a
  reference genotype included. A haploid row's `Number=2` pair with its second value missing (`51,.`,
  as CGI's vcfBeta writes every haploid row) reads `[51]`; any other count withholds the row's `hq`,
  and a value on an uncalled allele or outside 0..65535 is left out, each counted in one warning. HQ
  declared with a Number other than `2`, `P` or `.` is not carried and stays in `dropped_fields`.
  CGI's vcfBeta EHQ is not read (equal to HQ on all 4.21M rows of one PGP-HMS genome). No schema
  version change; `CONSENSUS_VERSION` 4, `REWRITE_VERSION` 10.

- Each of the four seats now hears the desk's answers to its own questions as a watcher line
  (`[seat.<seat>.watch.desk]`, the kit's `classes/desk/watch-answers.py`, kit v0.1.88, D133), so an
  answer no longer waits on a secretary's relay. The pin is at v0.1.88.

- docs/METHOD_ATLAS.md and docs/CONSENSUS.md: the two public CGI `var` pairs rebuilt, now banded
  (hu1B2689, CGA 2.5.0.44; huAF3C63, CGA 2.5.0.33). A `var` facet's bands match its vcfBeta's (one
  confidence scale), and on the odd/even check a banded `var` facet calibrates as the vcfBeta does
  (from 0.37 to 0.85 of predicted disagreements seen, 10 phred and up). huAF3C63's freeBayes 0.9.15
  bands run twice hu1B2689's, and its pair's even split leaves CGA 2.5.0.33's 40–100 phred bands on the
  floor. The 10–30 phred underconfidence is measured: the region rule, the highest rate among nested
  GIAB classes, over-predicts in `difficult` regions; with GIAB's `easy`/`difficult` partition alone
  `easy` is calibrated in every bin and, with the `var` facet, every bin from 10 phred up is within a
  tenth of 1 (1.01 overall). Which regions to pass is the person's call (four options there); the
  library rule is unchanged. `scripts/consensus_check.py` gains `--no-regions`, `--classes` and a
  `region` column in its calibration detail.

- `build_prive_reference` writes the reference's rows in the frequency file's order. The join chose an
  order of its own, which moved from one build to the next; the bootstrap's blocks are runs of rows, so
  `prive_proportions`' `proportion_errors` moved with it (an EUR person's error for the group it holds none
  of took values such as 0.0015 and 0.0027 from build to build, close to the test's `0 < e`: the likely cause of the flaky
  `test_prive_proportions_of_each_group[EUR]`, whose failure message was not kept). Proportions were not
  affected. A reference built before this keeps its old order; rebuilding it gives stable errors.

- Schema 0.0.19: `calls.hq`, each haplotype's phred-scaled score parallel to `gt`, filled from
  Complete Genomics `var` and masterVarBeta files (`varScoreEAF`, the lowest of a `var` haplotype's
  sub-rows; none on a `?` or partly called slot). `genotype_evidence` grades a fully called genotype
  whose every slot is scored by the lowest score, basis `hq` (`GradeBasis.HQ`), reference genotypes
  included and no vector. CGI `var` rows were `gt_only` (band `none`) and can now be banded in the
  method-error atlas; CGI's own VCF writes GQ as exactly that lowest score (42,894 of 42,895 loci on
  one chromosome of hu1B2689). A load warns where a file's VAF and EAF scores differ, and counts a
  score past u16. `varScoreEAF` leaves CGI's `dropped_fields`. `consensus` weighs such rows by their
  band's rate (`CONSENSUS_VERSION` 3); `liftover` carries `hq` (`LIFT_VERSION` 7) and `rejoin_split`
  keeps it where split rows agree (`REJOIN_VERSION` 4).

- `consensus(prior=)` takes an `AlleleFrequencyPrior`, a caller-supplied allele-frequency table keyed
  as consensus rows are (`chrom`, `pos`, `ref`, `alt`, `alt_frequency`; `ALLELE_FREQUENCIES`). An
  allele the table holds gets Hardy-Weinberg genotype priors over its row's ploidy; every other allele
  keeps the flat prior. `Consensus.support` gains `prior`, naming `flat` or the table on each row. The
  table is validated (frequency in [0, 1], no allele twice, no nulls, its assembly the facets'), and
  its name and digest enter the consensus `facet_id` and `derivation_inputs`. Which catalogue, if any,
  is the person's call (desk #96); nothing is downloaded or packaged, and the default stays flat, with
  default outputs and `facet_id`s unchanged.

- `genotype_evidence`, `genotype_evidence_at` and `consensus` take `block_grade=` (`BlockGrade`), the
  rule a gVCF `ref` block is graded by: `likelihoods` (the default, unchanged), `gq`, or
  `lower_of_gq_and_likelihoods`, which takes whichever of GQ and the likelihoods' confidence in hom-ref
  is lower. Which one should be the default is the person's call (desk #94). A block lacking what its
  rule reads is `block`, reference asserted with no number, never graded by another rule's quantity.
  Default outputs, consensus `facet_id`s included, are unchanged. `gq_on_reference` (desk #95) was
  already built end to end for calls; tests now cover it in `genotype_evidence_at` and `consensus`.

- docs/CONSENSUS.md: the consensus's 2x overconfidence at 60 phred and up, measured. It was freeBayes
  0.9.15 speaking with freeBayes 0.9.20's band rates, measured on another person's facet whose
  confidence scale differs (PGP's 0.9.15 floors GL at -10 and caps at 97 phred; its 40-60 phred rows
  err about five times as often as 0.9.20's). With its own bands fitted on the odd autosomes, the
  even autosomes' 60+ band is at 1.02 (0.209% predicted, 0.212% seen). CGI and the array sharing
  freeBayes's errors is ruled out (where freeBayes and the array agree, CGI disagrees on 0.010%), and
  the even split of a pair's top cell is second order. `scripts/consensus_check.py` writes
  `<out>.calibration_detail.tsv`, the calibration by band of the held-out facet's row and by which
  facets wrote a row there. Whether to borrow bands across versions is left to the person (three
  options in docs/CONSENSUS.md).

- `consensus` warns, naming each facet weighted by confidence-band rates its caller was measured at
  in other versions (`rate_basis = caller` on a band row). The rates and rows are unchanged.

- `Person` and `check_pedigree` no longer check a merged or consensus facet against the facets it was
  made from, nor stand two of them for one source: computed from its parents, it agreed with them by
  construction and counted as evidence that they are one person (merge audit 4).

- The stale-version sweep in `load` and the derived steps no longer removes an older version's cache
  entry that was used in the last `JUST_VCF_CACHE_GRACE_DAYS` (default 7; `0` restores sweeping at
  once). Two sessions on one cache dir at two `SCHEMA_VERSION`s showed why: the newer one's sweep
  removed a generation the older one's live frame read, and its next `collect` failed with
  `FileNotFoundError`. The sweep judges each entry again under its lock, logs how many older entries
  it kept, and reports them in the new `CacheClearReport.recently_used` (an added field, empty from
  `clear_cache`). An entry with no manifest is judged by the newest modification time in it. A frame
  kept open past the grace period without its step running again is still not protected.

- docs/CONSENSUS.md: the consensus calibrated again on hu1B2689 against held-out freeBayes, on the
  atlas with F21 and F22. With rates by confidence band it predicts 31,462 disagreements and sees
  27,245 (pooled rates: 51,225 against 27,243); 10 to 30 phred is still about twice underconfident and
  60 phred and up about twice overconfident. `scripts/consensus_check.py` predicts with the held-out
  facet's own rate per allele (its region classes and band). docs/METHOD_ATLAS.md gains the PGP
  pairs' discordance by confidence band; its pooled PGP rates already were F22's.

- The method-error atlas counts each pair per cell of the two sides' confidence bands (`a_band`,
  `b_band`: `genotype_evidence`'s confidence of each side's row, `0-10` to `100+` phred, `none`),
  beside the pooled rows (`all`, `all`). `method_error_rates` gives a rate per band of a method's own
  confidence (`confidence_band`) from those cells, and `consensus` weighs a row in a measured band by
  that band's rate instead of its method's average, so a confident row of a method is trusted as far
  as such rows were measured to be right, and an unsure one less (`CONSENSUS_VERSION` 2;
  `Consensus.support` gains `banded`). A consensus given an atlas with region rows and no `regions`
  warns. An atlas without band columns reads as before.

- docs/METHOD_ATLAS.md: the two PGP-HMS CGI pairs rebuilt with F21 and F22. SNV discordance falls
  from 3.40% to 0.54% and from 3.07% to 1.02%, indel from 7.95% to 4.98% and from 10.41% to 7.92%;
  the per-region SNV table is now the VRS pairing's, with no complex-row artefact left to take out.

- **Breaking:** `concordance` reads a haploid call against the same allele's diploid homozygote (a
  male's X outside the PARs, `1` against `1/1`) as `match`, with a new `sites.ploidy_differs` flag
  (null on a one-sided allele) that `summary` groups by (F22). `Agreement.PLOIDY_DIFFERS` is removed:
  a haploid call against a het is a `mismatch`. The atlas's `ploidy_differs` now counts the flagged
  alleles of `called_both`, which is `match + mismatch`, and `discordant` is `mismatch`, so the
  convention leaves the discordance rate. The graded comparison meets a haploid `k` with a diploid
  `2k`.

- `concordance` reads each side's genotype for the allele's change in its own row (F21). A complex
  multi-allelic row met by VRS id (freeBayes REF `CGCCC`, ALT `TGCCT,TGCCC`, GT `TGCCT/TGCCT` against
  `C>T` `T/T`) now reads two copies instead of zero, where the other ALT holds the same change; one
  that cannot be read (another length over the allele's span) is never zero copies. New `agreement`
  value `unreadable` for a pair that leaves undecided, and a matching `unreadable` count in the
  atlas, outside `called_both`; the graded distribution and the `site` pairing read alleles the same
  way, and the atlas counts an allele carried where `concordance` read a copy under another spelling.
  An added enum value and column: minor.

- `scripts/method_atlas.py` mints VRS ids on both facets of a pair once REF is read, so alleles pair
  by VRS id (`--no-mint` pairs by position and spelling, as before); `method_atlas.pairs.tsv` gains a
  `pairing` column. The atlas's two PGP-HMS pairs (CGI vcfBeta against freeBayes) were rebuilt that
  way, with CGI's REF read since F20: one-sided indels fall from 99% to 53%, and docs/METHOD_ATLAS.md
  now gives their indel and one-sided rates, with two conventions counted apart (X haploid against
  diploid, and freeBayes's complex rows holding an allele under another spelling).

- `resolve_ref` refuses a facet only when an autosome's stated length or md5 differs from the
  reference (F20). A differing X, Y or MT leaves that contig's rows unread (`ref` null, `ref_check`
  `not_run`) and the matching contigs are read; the warning names each mismatching contig with its
  row count and the stated against the reference value. A PGP-HMS CGI vcfBeta, refused before over its
  chrY md5, now gets 99.93% of its rows read, every one matching. `RESOLVE_VERSION` 9; `SCHEMA_VERSION`
  unchanged. `concordance` across two collections and `liftover` keep their own rule, which refuses on
  a differing X or Y too.

- The reference conflict warnings fire again for files whose names are withheld: `##reference=` and the call
  steps' reference options are compared by their base names as written, and only the stored `reference` is
  written `<file>`. Before, two different unknown names both became `<file>` first and agreed, so neither
  "more than one `##reference`" nor "differs from the caller's" fired. The three warnings now name no file
  ("the header names 2 different ##reference files"), known reference names included.

- `facets.caller_version` written as a path or a file name (`##DeepVariant_version=/home/<name>/…`) is
  `<path>`, as a method step's version already was, so the frame never holds the path; a caller split's
  facets get the same. Two versions are still compared as written, so two different paths stay a conflict,
  and the warning quotes each as `<path>`.

- `ftdna_bigy` is in the `variant_calls` family, not `genotype_table`: Big Y is targeted chrY sequencing
  with the reference base stated, as CASAVA's and SOAPsnp's tables are, so `chip_absence` now refuses it as
  not an array export instead of reading its absent sites against a chip. A Big Y facet's `family` changes.

- The first consensus check (RM280, docs/CONSENSUS.md): PGP-HMS hu1B2689's CGI, freeBayes, VarScan2 and
  23andMe facets on the autosomes, rates held out from this person's atlas pair. The consensus's top SNV
  genotype agrees with CGI on 99.88%, the array 99.85%, VarScan2 98.58%, freeBayes 98.27%; against
  freeBayes held out it is calibrated below 10 phred and two to three times underconfident above.
  `scripts/consensus_check.py` runs it from a plan, aggregates only.

- `method_error_rates` takes a pair's genotype disagreement as `mismatch / (match + mismatch)`, leaving
  `ploidy_differs` out: 69 to 82% of the PGP pairs' SNV `discordant` is a male X written haploid by CGI
  and diploid by freeBayes, a convention and not a genotype error. With it in, the consensus check
  was about seven times underconfident. `METHOD_RATES.called_both` is now `compared`.

- `consensus(person, atlas=...)` (RM280, `just_vcf.consensus`): a derived facet (`relation = consensus`, a
  link from every parent) with one row per ALT allele the parents state on one assembly, its `gp` the
  posterior over the copies of it. Each parent's reading (its row, its vector tempered by its rate, a
  site row, a `ref` block, or silence where `absence_is_reference` allows) is weighted by its method's
  error rate; `method_error_rates` turns the atlas's pair rates into per-method genotype and silence
  rates (triangulated where three methods were compared, else half a pair's), by caller and version,
  else caller, else a stated default (the worst measured), and `regions` picks the highest rate among
  the classes holding an allele. `Consensus.support` counts per row the parents that supported,
  contradicted or were silent and those weighted by a default; the parents keep their calls. Refuses
  two assemblies, a facet with its own derivative, and two subjects. No column changes;
  `SCHEMA_VERSION` unchanged.

- `links.CONSENSUS` (`consensus`): a documented multi-parent lineage relation, for a facet whose rows are
  computed from all its parents (RM280). `merge` now means the parents' rows put together.
  `link_disagreements` checks a `consensus` facet as it checks a `merge` one.

- Schema 0.0.18: `facets.gl_source` says whether `gl` was read from PL or GL; `genotype_evidence` reports it
  on rows graded `gl`.

- `genotype_evidence` and `genotype_evidence_at` gain `gl_source`: on a row graded `gl`, the key its likelihoods were read from (`gl`, `pl`, `gl_over_pl`), from `facets.gl_source`; null on every other basis. `grade_basis` is unchanged: GL and PL stay one basis, `gl`.

- `facets.gl_source` names the FORMAT key a VCF facet's `gl` is read from: `gl`, `pl`, or `gl_over_pl` where the header declares both (GL per row, PL where a row's GL is missing). Null where the facet carries no likelihoods. A new `facets` column; `SCHEMA_VERSION` 0.0.18.

- `method_atlas(pairs, regions)` (RM279, `just_vcf.atlas`): discordance and one-sided rates per
  (each side's caller and version, region class, variant class, scope) over pairs of one person's
  facets, with 95% Wilson intervals and the mean `p_same` per basis; counts and rates only, never a
  per-site row. `giab_stratifications` reads the GIAB v3.6 BEDs (easy/difficult, low mappability,
  segdups, homopolymers and tandem repeats by length, GC extremes, X and Y) into a `RegionSet` merged
  once and reused per pair; each pair is compared once and its alleles placed in their classes, and
  only alleles a genotype carries count. `pair_counts`/`pool_counts` split the run from the pooling;
  `scripts/method_atlas.py` builds one from a plan of pairs, each gated by `same_person`. No column
  changes; `SCHEMA_VERSION` unchanged.

- The first method-error atlas (RM279), over eight same-person pairs (three callers on the same reads,
  three DeepVariant releases, two PGP-HMS CGI-against-freeBayes pairs): headline rates in
  docs/METHOD_ATLAS.md. `scripts/method_atlas.py` takes one facet of a file split per caller
  (`a_caller`/`b_caller` in the plan) and compares a pair by label where `resolve_ref` refuses a facet
  (CGI vcfBeta's chrY md5 differs from the store's GRCh37), saying so per pair.

- `UnsupportedCoverageError` names the kind of input (`the coverage input is an empty file`, `the
  BAM header lists no reference sequence …`, `the alignment input is not a BAM or CRAM …`), never
  the file's name, as `DamagedBgzfError` carries no path.

- The store's `manifest.json` keeps each facet's load warnings redacted as `query` returns them
  (the facets' sample, source path and member, and anything shaped like a path), where it kept them
  as the load gave them; `Store.loaded().warnings` reads them back so.

- `reference`, and a reference or chain named in a method step's options, keep a file's base name
  only where every word of it is a known reference's (`hg38.fa`, `hs37d5.fa.gz`,
  `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, `GRCh38_to_GRCh37.chain.gz`, DRAGEN's
  `hg38_alt_aware`); any other is written `<file>`, since a base name like `ALICESMITH-hg38.fa` names
  someone. A value change on facets whose header names such a file.

- A method step's `tool` or `version` written as a path or file name is `<path>`, as an option's
  value is: `##liftOverProgram=</home/<name>/bin/CrossMap,…>` put a home directory in
  `method_steps`. A value change on such facets only.

- `chip_absence` reads only a genotype table as an array export: a headerless VCF, a CASAVA or
  SOAPsnp site table, a Nirvana, Color or GET-Evidence export and a CGI `var` or `masterVarBeta` file
  state no VCF version and load with no callset type, so they were read as arrays, and a site such a
  file has no row at could come back as a chip's `no_call`. They are refused as not an array export.

- A lifted facet's `reference` is null: it kept the parent's, which named the source assembly's
  file, while its rows are against the target (`ref_seqcol_digest`). `aligner`, `aligner_version`
  and `method_steps` stay the parent's, since the reads were aligned to the source. `LIFT_VERSION` 6:
  lifted facets get new `facet_id`s.

- `write_vcf_zarr` writes a row with a slot in a VCF 4.5 PSL phase set unphased, counted in a
  warning, as it does a hapLink's: VCF Zarr has one phase flag per call and no PSL, so `true` put the
  row in the chromosome's implicit set. The store's `fileformat` follows `write_vcf` (VCF 4.5 where a
  row holds a PSL set) instead of always saying 4.3.

- Relatedness names no band when a facet is homozygous reference at none of the shared sites
  (`variant_sites_only`): its own genotype chose them. This is the Genos exome's offset from RM274's
  note (a gVCF of no-call blocks only), measured: the participant's array cut to the exome's sites
  reproduces it. Kinship and SE are kept, `same_person` still reads them, and the pedigree check reads
  such a pair `cannot_tell`. Values change: those pairs read `unknown` where they read a band.

- A supplied GA4GH Pedigree is checked against measured relatedness (RM274, layer 4).
  `expected_kinship(pedigree)` gives each pair of individuals its expected kinship and bands, through
  parent links, stated KIN degrees, or the founders-unrelated assumption, each row saying which.
  `check_pedigree(pedigree, *people_or_loaded)` matches individuals to facets by subject and reports
  `agrees`, `disagrees` or `cannot_tell` with a reason per pair of facets; `people(..., pedigree=...)`
  warns on disagreements naming facet ids. Nothing in the pedigree is changed. `pairwise_relatedness`
  takes `pairs=` to measure only the pairs asked for. No schema version change: two computed frames.

- Schema 0.0.17: `phase_set_kinds` gains `cgi_haplink`; CGI hapLinks load as phase sets, and masterVarBeta read
  counts fill `ad`.

- CGI `var` files of format 1.6 and 1.7 (CGA 1.11, 1.12: one `totalScore`, no filter column) are
  tested: they load as 2.5 does, with `sample_filters` null and `totalScore` in `dropped_fields`
  (RM2). No code change; the loader already read columns by name.

- CGI's `hapLink` is read as phase sets (RM2): a new `PhaseSetKind.CGI_HAPLINK` (`cgi_haplink`); a
  linked haplotype of a `var` or `masterVarBeta` row has `phased` true and the id in `phase_sets`,
  and `cis_phased_blocks` gives one VRS CisPhasedBlock per id and chromosome, the same from both
  files. `write_vcf` has no VCF field for a hapLink: it writes the haplotypes one phases unphased,
  with a warning, since `|` without the set would read back as the chromosome's implicit set;
  `write_vcf_zarr` writes such a row's `call_genotype_phased` false, for the same reason. Adding the enum
  member changes the `phase_set_kinds` dtype: schema 0.0.17.

- CGI `masterVarBeta` rows carry `ad` (RM2): REF's count is `referenceAlleleReadCount`, each ALT's the
  read count of the slot that carries it, indexed over `source_alleles = [source_ref] + alts`. A row
  with an empty count, a partly called allele or counts that contradict each other keeps `ad` null,
  counted in a load warning by reason (`loaders.cgi.AD_WITHHELD_REASONS`). The read count columns
  leave `dropped_fields`. A value change on masterVarBeta facets; no column changes.

- Schema 0.0.16: the method chain on the facet (`method_steps`, `reference`, `aligner`, `aligner_version`);
  the CLI and `query` pass `reference` and the aligner through their path and name redaction.

- The two-process tests in `tests/test_rerun.py` fail within ten minutes instead of hanging the suite.
  A child put its result (about 57 KiB) on a `multiprocessing.Queue`; once the parent's loop failed
  (its 120 s deadline passes under a heavy load), nothing read the pipe, the second child waited in
  its queue's feeder thread forever, and pytest's exit joined it forever. Children now write their
  results to files, are daemonic and are killed when a test fails, and a missed deadline prints each
  live child's stack. Test-only: the library starts no processes.

- Graded genotype evidence across sources (layer 7, charter 7): `genotype_evidence(loaded)` gives every
  `calls` and `blocks` row a probability vector over its genotypes in VCF order where the source
  supports one (GP as given; GL/PL normalized under a flat prior, named in `prior`; GQ alone on its
  called genotype, the rest spread evenly), a phred-like `confidence` (null where nothing quantifies it,
  never 0 for unknown) and its `grade_basis` (`gp`, `gl`, `gq`, `dosage`, `gt_only`, `block`, `no_call`);
  `genotype_evidence_at(loaded, sites, variant_class=…)` adds a facet's absence reading at asked
  positions (`absent_reference`, `absent_unknown`, through `absence_is_reference`) and, with
  `coverage=`, callability beside it. Computed, not stored: no column changes and `SCHEMA_VERSION`
  stays. GQ grades no reference genotype unless `gq_on_reference=True` (VCF defines it conditional on
  the site being variant). docs/SCHEMA.md, "Graded genotype evidence".

- Schema 0.0.15: `calls.phase_set_ordinals`, `phase_set_qualities` (VCF 4.4+ PSO/PSQ); PSL phase sets load and
  write back under VCF 4.5.

- A BCF facet now takes the build its header states (`##reference=…hs37d5.fa` gives GRCh37, basis
  `stated`), as the same header in a text VCF always did. `sniff` stops at a BCF's magic bytes, so the
  loader reads the header's build lines itself, by the VCF spec's hints. Before, a BCF's
  `assembly_label` came only from its contig lengths. Found by the facet guard (`tests/test_load_seams.py`).

- A guard against merge seams in the loaders (`tests/test_load_seams.py`): every loadable format's
  spec example, a BCF written from the `vcf` one and the appended-export table are loaded from a zip
  through `load`, and any `facets` column left null fails unless the test states why null is right for
  that format, caller, family or every facet. A reason has to hold for every facet it covers, so a
  family-wide excuse cannot hide one loader that drops a column another fills. A new `FACETS` column
  fails for every format until its reasons are written. The same loads check that each load path
  returns `links` of the `LINKS` schema that agree with its facets, and that `subject` reaches every
  facet. One gap it found, a BCF's header build hints never read (no `assembly_label` where the
  same header as text gets GRCh37), is fixed in the entry above.

- The CLI printed header text unredacted. `sniff` printed `build_evidence`, the whole matched header line
  (a `##reference=file:///home/<name>/…`); it prints `build_hint` instead, the statement's key and the
  build it names (`{"kind": "##reference", "label": "GRCh38"}`). `load`, `convert` and `check` printed the
  facet's free-text fields as read (a `##DeepVariant_version` path landed in `caller_version`); every
  one `query` redacts (`query.PROVENANCE_TEXT`), and `dropped_fields`, now goes through `query`'s
  redaction (`query.FacetRedactor`) and then the CLI's labels. The CLI's warnings and refusals also
  lose anything shaped like a path, by the same rule (`query.redact_paths`), which now also catches
  `file:///` URLs, `//host` shares and one-component `~/` and `C:\` paths, in `query`'s answers too.
  Breaking for a script reading `build_evidence` from `sniff --json`.

- `chip_absence` and `learn_probe_sets` refuse a `Loaded` of more than one facet (an appended export)
  with `ValueError`, as the derived steps do: they read the first facet's chip and assembly against
  every facet's rows. `scripts/learn_probe_sets.py` skips such a file and counts it.

- `Person` carries `links` (RM275's frame), which it dropped: every link its facets' loads carried,
  once each, in `canonical_links` order, so a lifted facet grouped with its parent keeps its
  `liftover_to` link. The CLI's `--facet` narrowing kept no links either, so a step after it lost the
  load's kin links; it keeps them now. No column changes.

- `rejoin_split` keeps the SV columns RM282 added, which it wrote null: `sv_type` and `copy_number`
  where the split rows agree, and `cipos`, `ciend` and `alt_copy_numbers` concatenated in the joined
  `alts` order (null when a row does not state them per ALT, or two rows disagree on one ALT).
  `REJOIN_VERSION` 2: re-joined facets get new `facet_id`s.

- `liftover` withholds a swapped row's `alt_copy_numbers` (INFO CN, one per ALT), as it withholds its
  `dosage`: the swapped ALT's slot holds the old REF, whose copy number the source did not state; a
  warning counts them. `LIFT_VERSION` 3. `NORMALIZE_VERSION` 6, for RM282's change to what `normalize`
  emits (per-ALT SV columns split with their ALT, a moved row's intervals cleared), which shipped
  without a bump. Lifted and normalized facets get new `facet_id`s.

- The suite's opt-in variables (`JUST_VCF_SAMPLES_DIR`, `JUST_VCF_ANCESTRY_DIR`, `JUST_VCF_MANIFESTS_DIR`,
  `JUST_VCF_TEST_REFERENCE_DIR`, `JUST_VCF_TEST_CHAIN_DIR`) are empty unless the shell sets them:
  `tests/conftest.py` no longer reads them from the repo `.env`, and defaults each to `""` so the CLI's
  `load_dotenv` cannot fill them either. A plain `uv run pytest` never opens a private genome; an opt-in
  run names the variables on its command line (CLAUDE.md §4, docs/LIVE_SAMPLES.md). No library change.

- A withheld CisPhasedBlock names the first slot that withheld it: three new columns of
  `cis_phased_blocks`, `withheld_chrom`, `withheld_pos` and `withheld_allele` (RM284), so the one
  symbolic ALT that withholds a chromosome's implicit set can be found. Null when the block has an id
  or is `too_few_members`. Offering blocks per window instead stays the person's call.

- VCF 4.4+ phase-set lists are read and written (layer 8). The VCF loader reads FORMAT `PSL` into
  `phase_sets`/`phase_set_kinds` (kind `psl`) per allele, and `PSO` and `PSQ` into two new nullable
  `calls` columns, `phase_set_ordinals` and `phase_set_qualities` (schema 0.0.15; `liftover` and
  `rejoin_split` carry them, `LIFT_VERSION` 5 and `REJOIN_VERSION` 3). A
  PSL name on a VCF 4.4+ record's first slot makes that slot phased, which P4 otherwise leaves
  unknown. polars-bio refuses `Number=P` (P10), so such a header is read through the rewritten copy
  declaring `Number=.` (`REWRITE_VERSION` 9); before this, a file with PSL loaded with PSL dropped.
  `write_vcf` writes VCF 4.5 when a row holds a PSL set (PSL, PSO, PSQ, a leading GT phase indicator
  where it differs from the implicit one, CN declared Float), and 4.3 otherwise, as before; a PSL facet
  round-trips through `load`, and htslib reads it. `cis_phased_blocks` keeps a PSL record's
  implicit-set slots out of its PSL sets.

- A facet records its method chain (layer 2, charter 6): four new nullable `facets` columns from the
  `source` detector. `method_steps` lists every step the header records, in the order they ran, each
  with its tool, version, role (`align`, `call`, `filter`, `normalize`, `merge`, `lift`, `annotate`,
  `phase`, `other`), key options (paths written `<path>`, a reference or chain by its base name, no
  sample names) and the meta key it came from: GATK in its three header styles, DRAGEN, bcftools,
  vcftools, freeBayes, starling, Clair3, SnpEff, VEP, gvcftools, CrossMap, Eagle, WhatsHap, SHAPEIT4
  and HapCUT2. Dated steps (GATK, DRAGEN, bcftools) are put in date order, because htsjdk sorts meta
  lines by text. `reference` is `##reference=`'s base name or the caller's; `aligner` and
  `aligner_version` are DRAGEN's when it mapped the reads. `phasing_method` now also names a caller
  that declares FORMAT `PS` in its own words (DRAGEN, GATK HaplotypeCaller, DeepVariant, CGI) and the
  phasers above; two phasing steps leave it null with a warning. `##clair3_version` names Clair3 as
  the caller. A `##reference=` that differs from the caller's with no liftover recorded is a warning.
  The CLI summary prints `reference`, `aligner` and each step's role, tool and version. Schema
  0.0.16.

- `chip_absence` separates a written no-call from a call (RM281): a site the export writes as `--`
  (AncestryDNA: `0 0`) is `no_call`, where it was `in_export`, with or without a probe set; `in_export`
  now means the export has a genotype there. `AbsenceCause` keeps its three members; the meaning of
  `in_export` narrows, a value change for any caller reading it on no-call sites.

- `resolve_ref` resolves `I`/`D` indel codes only on a chip whose calibration passed (RM281,
  `just_vcf.probes.ID_CALIBRATION`): `23andme_v5` and `ancestrydna_v2`, where 98.7% and 98.4% of the
  resolved homozygous rows read hom-ref. `23andme_v3` and `23andme_v4` (76%, the same sites hom-alt in
  every export; on v3 a same-person CGI genome reads all of them hom-ref, so the codes run the other way
  there), `23andme_v1_v2` (6 rows), Final Reports and unknown chips keep their code rows as
  loaded, with one warning counting them as `id_orientation_unverified` and naming the evidence. A value
  change on those facets when the manifests are configured; `RESOLVE_VERSION` 8.

- Illumina Annotation Engine (Nirvana) JSON loads (`nirvana_json`, RM2, layer 1): each position object
  is the VCF record it annotated, so rows get `source_ref`, `alts`, `gt`, `phased`, `gq`, `dp`, `ad`,
  `qual` and `filters`; the annotations are listed in `dropped_fields`. Several samples load one at a
  time with `sample=`. 32 of 39 declared formats load. No column changes.

- Six more formats load (RM2, layer 1), through the new `loaders.site_tables`: CASAVA SNP and indel
  tables, SOAPsnp consensus rows, FamilyTreeDNA Big Y CSVs, Color's `variant_set.csv` and GET-Evidence genome
  GFFs. Each writes a reference base beside the genotype, so rows get `source_ref`, `alts` and a `gt`
  over `[source_ref] + alts` at load; CASAVA rows also `ad` and `dp` (indels `ad`, `gq`, `qual`, the
  padding dashes dropped, insertion points at the base before them), SOAPsnp `gq` and `dp`;
  GET-Evidence's `REF` rows are `ref` blocks (`absence_means = blocks`). Big Y's Known SNP rows, which
  have no position, are dropped with a count. 31 of 39 declared formats load. No column changes.
  docs/FORMATS.md says why each remaining format does not load, and which are not facets.

- In an array facet whose strand is `mixed` (F13), `resolve_ref` reads the rows it withheld by their
  probe's two alleles on the plus strand, from the Illumina manifests under `JUST_VCF_MANIFESTS_DIR`
  (RM281): matched by rsid, read only where the pair holds the reference base and is not `A/T` or
  `C/G`; a row whose alleles are the pair's complement is a minus-strand row, and its `gt_alleles` are
  complemented (`source_gt` keeps the text). Homozygous rows of such probes, withheld until now, get
  `alts` and `gt`. Forward facets are unchanged. One warning counts every outcome.

- An array's absence gets a cause (RM281): `chip_absence(loaded, sites)` says per site whether the
  export has a row there (`in_export`), the chip has a probe there and the export reports nothing
  (`no_call`), or the chip has no probe there (`not_on_chip`); null with a reason when the facet has no
  probe set. `probe_set(chip, assembly_label)` is the chip's own manifest for an `illumina:<manifest>`
  chip (configured), or the sites at least half of the chip's exports carry, learned by `learn_probe_sets` into the
  cache's `probes/learned/` (`scripts/learn_probe_sets.py` walks a directory of exports). Platform
  families named by marker overlap get none. New public names: `chip_absence`, `ChipAbsence`,
  `AbsenceCause`, `probe_set`, `ProbeSet`, `ProbeSetKind`, `learn_probe_sets`. No column changes.

- `resolve_ref` resolves an array's `I`/`D` indel codes (23andMe, AncestryDNA, Final Reports) against
  the indel probes of the Illumina manifests under `JUST_VCF_MANIFESTS_DIR` (RM281, `just_vcf.probes`):
  a row is matched with a probe by rsid, the probe's site is located in the reference by its flanks,
  and the row gets `ref` (from the reference), `alts`, `gt`, `end_pos` and `gt_alleles` spelled as the
  site's two sequences, REF starting at the row's own position; `source_gt` keeps the codes. Each
  manifest is read once into a probe table under the cache's `probes/`. Rows no probe resolves stay as
  loaded, each reason counted; without the manifests every code row stays unread, with a warning that
  says so. A value change on array facets with I/D rows when the manifests are configured. Tests see
  `JUST_VCF_MANIFESTS_DIR` empty; the opt-in manifest test takes the new `manifests_dir` fixture.

- A GRCh37 or GRCh38 scaffold now meets the reference store under any name NCBI's assembly reports
  give it (F19): UCSC's `chr1_KI270706v1_random`, GenBank's and Ensembl's `KI270706.1`, RefSeq's
  `NT_187361.1`, GRCh37's `chr4_gl000194_random`. `ReferenceStore.contigs_by_chrom` lists a contig
  under its report names only where its length is the report's and the collection holds no contig of
  that name; stated lengths and chain sizes are checked against it as before. `resolve_ref`,
  `normalize`, `mint_vrs_ids` and `liftover`'s target contigs all go through it, and `resolve_ref`
  says how many rows it read that way. `chrom` keeps the source's spelling. Values change on files
  that already load: such rows now get `ref`, alignment, VRS ids and lifts. `RESOLVE_VERSION` 7,
  `NORMALIZE_VERSION` 7 and `LIFT_VERSION` 4 (derived `facet_id`s change), `MINT_VERSION` 3. Measured:
  a PGP-HMS GRCh38 VCF goes from 44,594 rows unread to 0; GIAB HG002 lifted CHM13 → GRCh38 from
  75,001 `target_contig_unmatched` to 40,099 (all on `_fix` patches, which the Ensembl store lacks).

- GRCh37's and GRCh38's assembled-molecule accessions map to canonical `chrom`, as T2T-CHM13's do
  (`CM000663.2`, `NC_000001.11` → `1`; `J01415.2`, `NC_012920.1` → `MT`), from NCBI's GRCh37.p13 and
  GRCh38.p14 assembly reports, now carried in `assets/` and derived into the package table
  `contig_names.tsv` (`just_vcf.contig_names`, F19). A header named by accession is labelled with its
  assembly. Values change on files that already load: such a file's `chrom` and `contig_map` now say
  `1`…`MT` instead of the accession.

- `write_vcf_zarr(loaded, path)` exports one facet's calls as a VCF Zarr 0.4 store, which polars-bio
  reads back with `scan_vcf_zarr` (RM276). Rows follow `write_vcf`'s rules (REF from the reference
  or the file, never the genotype; unwritable rows refused or counted). The VCF Zarr spec has no
  blocks, so a facet's blocks stay in the parquet store and are counted. Written with the standard
  library (gzip chunks, no timestamp), so no extra and the same bytes for the same frame.

- A durable store per person (RM276): `Store.open(path, subject=...)` creates or opens a directory
  the caller names, `store.add(loaded)` writes each facet's four frames as parquet with a manifest
  (schema version, subject, facet ids, added times in UTC), and `store.loaded()` / `store.person()`
  read everything back lazily. Facets accumulate: a facet already held is a no-op and adding one never
  rewrites another; each add is atomic per facet and holds a lock, so concurrent adds are safe. Another
  schema version or subject is refused, never mixed. No column changes; the schema version is not
  bumped.

- `Person.of` and `people` check a subject's facets with `same_person` (RM274) when no check is named;
  `check=None` runs none, as no check did before. A pair of facets too small or on two assemblies reads
  `could_not_tell`. A behaviour change on RM273's API, no column changes.

- GA4GH Pedigree read as supplied data (RM274): `read_pedigree` takes the standard's Phenopacket
  `GA4GHPedigree` message as JSON (both protobuf JSON spellings of a field) into a `Pedigree` with an
  individuals frame and a relationships frame, KIN terms and sex values verbatim. Nothing is inferred
  from it or into it; a document without its required parts raises `PedigreeError`.

- somalier as an optional second opinion on relatedness (RM274): `just_vcf.somalier.somalier_relate`
  writes each facet at the panel's positions with `write_vcf`, runs `somalier extract` and `relate`,
  and reads `pairs.tsv` into the relatedness frame (`method = somalier`, its `relatedness` kept
  verbatim, half of it as `kinship`). Without a `somalier` binary on PATH it raises
  `SomalierNotFoundError`. Not run on this machine: docs/TRUST.md has no row for its binary yet. The
  module's names (`somalier_relate`, `read_somalier_pairs`, `somalier_binary`,
  `SomalierNotFoundError`) are module-level, under `just_vcf.somalier`, not in `just_vcf.__all__`: an
  optional runner of an outside binary, imported only by a caller who wants it.

- Relatedness between facets (RM274): `relatedness(a, b)` gives KING-robust kinship (between-family
  form), IBS0, IBS2, the shared site counts, a jackknife standard error and a band (`same_person`,
  `first_degree`, `second_degree`, `third_degree`, `unrelated`, or `unknown` with a reason);
  `pairwise_relatedness` does every pair of facets across several `Loaded` as a frame, and
  `same_person(a, b)` answers True, False or None. Sites come from a panel of 15,289 autosomal,
  strand-unambiguous SNVs per build (GRCh37, GRCh38) derived from somalier's site VCFs and shipped in
  the package (`relatedness_sites.json`); genotypes are read as ancestry reads them, now from the core
  module `just_vcf.site_genotypes`, so arrays and whole genomes compare alike. New API, no column
  changes.

- Schema 0.0.12 (RM273): `facets.subject`; `Person.of` groups one subject's facets and checks the claim.

- `Person.of(*loaded, check=...)` and `people(*loaded, check=...)` group facets by the subject the
  caller supplied, with their calls and blocks, and check the claim with a pluggable relatedness test
  (RM273). `Person.checks` keeps each pair's answer (`same`, `different`, `could_not_tell`, `not_run`)
  and `Person.one_person` their Kleene AND; a `different` pair is a warning naming both facets, and the
  facets stay where the caller put them. A file loaded under two subjects is allowed and warned. No
  check is wired by default yet: RM274's relatedness takes that place.

- `load(..., subject=...)` labels a file's facets with the person they belong to, in a new nullable
  `facets.subject` column (RM273). The label is the caller's: nothing fills it from the sample name or
  the file, and `facet_id` does not digest it. Every derived step (`resolve_ref`, `normalize`,
  `liftover`, `rejoin_split`, `mint_vrs_ids`) carries its parent's. Schema change, landed in 0.0.12.

- `query(source)` (RM283), the aggregate-only helper beside `load`, for agents and apps: counts per
  facet by contig, variant class (ALT alleles), filter set, genotype class, REF check, evidence basis
  and block kind (with the bases blocks cover); GQ/DP/QUAL/MIN_DP summaries; provenance without the
  sample, paths, archive member or header lines; warnings grouped by reason with paths and the sample
  replaced; and `carries(...)`, whether a facet carries one allele asked by VRS id or position, as a
  Kleene `Presence` with its basis (stated, site genotyped, ref or no-call block, absence read through
  `absence_is_reference`, unknown). Not-carried is answered for SNVs only. No row, sample name or
  unasked position appears in any answer (docs/TRUST.md). No column changes.

- The facet graph as a fourth frame, `Loaded.links` (RM275, SCHEMA.md Decision 12): one row per link,
  `from_facet_id`, `to_facet_id`, `relation` and a `provenance` struct (the step, its rule's version,
  the files it read). `rejoin_split`, `normalize` and `liftover` add a lineage link and carry their
  parent's links, so a facet several steps from its file keeps every step; `resolve_ref` and
  `mint_vrs_ids` pass them through. A load of several facets links them: a caller split's facets in
  every pair, an appended export from the first export. A plain load has none. `facets.derived_from`
  and `relation` stay for the one-parent case, and `link_disagreements` lists every place the two
  disagree. `LINKS` and `link_disagreements` are exported. No column of the other frames changes; the
  frame landed in schema 0.0.13.

- A command line, `just-vcf` (RM272): `sniff`, `load` (a summary of each facet), `check` (REF read
  from a reference, `ref_check` and withheld VRS ids counted) and `convert` (to VCF or gVCF through
  `write_vcf`, with `--reference`, `--normalize`, `--lift` and `--target-reference`, `--vrs`, `--gvcf`,
  `--resolved-only`). It prints counts and aggregates only, as text or `--json`; `--sample` is never
  echoed, library warnings and refusals are redacted of the input's path, file name, sample and member,
  and logging is off while a command runs (polars-bio's import configures the root logger, which would
  print each library warning unredacted). Refusals exit 3 to 11, one code per kind (`cli.ExitCode`); a
  crash exits 1 with its message withheld. No column changes; typer and python-dotenv were already core.

- docs/FORMAT_REGISTRY.md (RM285): every declared format as a language-neutral spec (detection rules,
  columns, coordinates, genotype encoding, build statements, whether it loads, a synthetic example)
  and the broken and odd shapes met in real files, F1-F17 and the dialects behind FORMATS.md's
  contract gaps, each with a synthetic example and what just-vcf does with it. Generated by
  `scripts/format_registry.py`; a test fails when the doc and the code disagree, and every shape's
  example is sniffed in the suite. No code or frame changes.

- Structural-variant and copy-number VCFs load as SV facets (RM282). Five nullable `calls` columns:
  `sv_type` (INFO SVTYPE as written), `cipos` and `ciend` (two offsets per ALT; VCF 4.2's one pair
  for the record repeated per ALT), `copy_number` (FORMAT CN) and `alt_copy_numbers` (INFO CN, per
  ALT); they landed in schema 0.0.14. Canvas's `<CN0>`,
  `<CN2>`… end at INFO END, and its copy-neutral segments (ALT `.` with END and CN) are calls, not
  reference blocks (F18), so a Canvas facet's absence is `unknown`, not `blocks`. `vrs_withheld`
  gives `sv_no_vrs_model` for an SV allele or breakend (it gave `symbolic`, which now names only
  symbolic alleles that are no variant, CGI's `<CGA_NOCALL>`). `write_vcf` writes SVTYPE, CIPOS,
  CIEND and CN back (CN as the Integer VCF 4.3 reserves; a fractional one is left out and counted),
  and the rows load back equal. A VCF 4.3 header declaring CN Float, which polars-bio refuses, loads
  without CN, with a warning. `normalize` splits the per-ALT columns with their ALT and clears the
  intervals of a row it moves, as `liftover` does. `REWRITE_VERSION` 8. Values change on SV and CNV
  files that already load. Breaking for a consumer matching `vrs_withheld == "symbolic"` on SVs.

- `cis_phased_blocks(loaded)` exports the phase sets as VRS 2.x CisPhasedBlocks on demand (RM284):
  one row per haplotype of each PS set (per GT slot) and per PSL set, the member Allele ids and the
  block's `ga4gh:CPB.` id, equal to vrs-python's on synthetic facets and on a phased Genes for Good
  upload from PGP-HMS. Members are the non-reference alleles unless `reference_alleles=True`. A
  block that cannot be stated whole is withheld with a reason (`phase_unknown`, `uncalled`, the
  member allele's `vrs_withheld` reason, `too_few_members`), never shortened. No column changes.

- Sex-chromosome evidence on T2T-CHM13 (RM286): `PAR_BOUNDS` and `X_DEGENERATE_Y` gain T2T-CHM13 (GIAB
  v3.6 CHM13 PAR and XTR, the T2T consortium's X-degenerate class), so a T2T facet's X and Y are grouped
  and its het share and depth ratios measured; the estimate is withheld with the new reason
  `uncalibrated_build`, since the thresholds were measured on GRCh37/38 genomes only. It answered
  `unsupported_build` before. `karyotype.REGION_SET` is now a dict by build (`region_set` names the T2T
  source too), and `REGION_BUILDS` joins `ESTIMATE_BUILDS`. A value change; no column changes.

- T2T-CHM13 names (RM286): a VCF whose contigs are T2T-CHM13v2.0's GenBank (`CP068277.2`…) or RefSeq
  (`NC_060925.1`…) accessions maps them to canonical `chrom` and gets its `assembly_label` from their
  lengths; those names kept their source spelling before. A `##reference` naming `chm13`, `t2t`, UCSC's
  `hs1` or the assembly's accession, or a `##contig` `assembly=` of `chm13…`, `t2t…` or `hs1`, states
  `T2T-CHM13`. The `T2T-CHM13` label cannot vouch for MT (the CHM13 mitochondrion and rCRS are both
  16,569 bases), so an array labelled T2T-CHM13 keeps MT unread. Value changes on files that already
  load; no column changes.

- `sex_chromosome_evidence(loaded, coverage=…)` takes a coverage facet (RM269 with RM278): the
  median read depth over non-PAR X and the X-degenerate Y against chr20 is a third X line and a Y line
  that can show a Y absent, so `XX` and `XO` are now named where the lines agree. Without coverage
  nothing changes; the result gains a nullable `coverage` field (`just_vcf.karyotype.CoverageRatios`,
  module-level like the result's other nested models) and two reasons
  (`coverage_no_build`, `coverage_build_mismatch`, `coverage_uneven` for an exome-shaped depth).

- Coverage as a facet (RM278): `read_coverage` reads mosdepth's per-base, regions, quantized and
  summary outputs and `samtools depth` / `samtools coverage` into a `Coverage` (a `CoverageFacet` and
  a lazy `intervals` frame of depth mean, bounds and described bases per interval);
  `alignment_coverage` runs polars-bio's `depth()` over a BAM or CRAM and folds its per-base depth
  into 1,000-base windows while it streams (13.3 GB peak, 12.7 minutes on a 27 GB BAM);
  `callable_sites` and `callable_at` say per position whether the reads reached a `min_depth`, `True`, `False` or unknown with a reason. A headerless
  three-column file that could be a BED3 is refused unless `kind=` states it. Nothing in `facets`,
  `calls` or `blocks` changes, so `SCHEMA_VERSION` stays (docs/SCHEMA.md, "Coverage facets").

- `concordance(a, b, *, region=None)` compares two facets of one person on one assembly, allele by
  allele (RM277). Alleles meet by VRS id where both facets minted one, else by position, REF and ALT,
  else against the other facet's one row at the position. Each row states the agreement (match,
  mismatch, ploidy differs, uncalled, one-sided), what a silent facet's absence means, and, where both
  rows carry likelihoods or both posteriors, `p_same`, the probability the two hold the same copies of
  the allele (likelihoods under a named flat prior). A summary counts alleles by agreement and variant
  class. Disagreement is kept, never resolved. Facets on different assemblies are refused. New public
  names: `concordance`, `Concordance`, `Agreement`, `PairKey`, `AbsenceReading`, `ComparisonBasis`.
  No column changes; `SCHEMA_VERSION` is unchanged.

- A PGP-HMS freeBayes+VarScan2 merge record whose INFO/SF names no caller (SF absent, or an index no
  `##sourceFiles` line declares) is counted in a load warning; it lands in no facet, as before, but
  was dropped without a word. Counting it reads the file a second time, beside the first part's stats
  pass. No column changes.

- A Complete Genomics `var` or `masterVarBeta` whose header does not say it is a whole genome (28 of
  250 PGP-HMS files have no `#TITLE`) gets its `callset_type` from its rows, by the rule a VCF's
  follows; it was null. `absence_means` stays `blocks`. A value change on files that already load.

- Complete Genomics' native `var` and `masterVarBeta` facets and the PGP key=value table's facets fill
  `variant_classes`, `content_filter` and `content_filter_basis` as a VCF facet does, from the rows and
  the file name; they were null. CGI's unanchored deletions and insertions (`""`) count as indels, so
  a whole CGI genome shows no content filter rather than reading `snv_only`. A value change on files
  that already load. No column changes.

- An appended export (F3) no longer takes the build the first export's preamble states: its
  `assembly_label` is null until `resolve_ref` infers it from its own rows (F9), where it read the
  first export's label as `stated`. A value change on files that already load. No column changes.

- A read of an intact BGZF input that fails for another reason walks the file's blocks once, where it
  walked them twice (once at the facet, once at the source) before re-raising the reader's error.

- The cache sweep removes only entries an older version wrote, comparing version keys numerically
  (`0.0.9` is behind `0.0.11`). An entry a newer just-vcf wrote is kept, so two installs sharing the
  platform user cache no longer delete each other's live spills on every load; before, any version
  other than the running one was swept.

- Schema 0.0.11: `absence_means` gains `reference` (a Complete Genomics VCF typed `wgs`). Values change on
  files that already load: the empty allele is `""` everywhere (Genos and upload `-` included), and
  `dropped_fields` is `[]` rather than null where nothing was dropped.

- A zero-length allele has one spelling, `""` (VRS's empty literal), in `ref`, `alts`,
  `source_alleles` and `gt_alleles`. An upload's ALT `-` (F12) and Genos's `-` load as `""`, where they
  were `"-"`; `source_gt` keeps what the file wrote. Every reader agrees: `""` is an `indel` in
  `variant_classes` and a site in the callset count; `normalize` anchors CGI's native deletions on the
  base before, as it does F12's (`NORMALIZE_VERSION` 5); `mint_vrs_ids` mints a resolved unanchored
  deletion with the id of its anchored spelling, where it withheld it as `not_sequence`
  (`MINT_VERSION` 2); `write_vcf` refuses an empty REF (a CGI insertion point) as `ref_unresolved`,
  where it wrote an empty REF field. Values change on files that already load.

- `normalize` aligns and trims the rows `resolve_ref` gave a REF with no source REF to check, in a
  facet read from a file: a Living DNA genotype split against the reference that spells an MNV or a
  deletion with extra shared bases (`CTA`/`CGA`) is trimmed (`T>G` one base on). Before, such rows passed
  through as unchecked. Array SNVs are already minimal and do not change. Normalized facets get new ids
  (`NORMALIZE_VERSION` 4).

- A Complete Genomics VCF reads its silence as reference calls: `AbsenceMeans` gains `reference`
  ("every position the file omits was called reference by its producer"), set on a VCF whose caller
  line names CGAPipeline when its rows show a whole genome, and `absence_is_reference` answers True for
  it, by the same class and content-filter gates as `weak_reference`. Settled by comparing a PGP-HMS
  vcfBeta with the same assembly's `var` (0 of 317.7M non-reference bases missed). A smaller CGI VCF
  keeps `blocks`, with a warning. Ancestry reads an omitted site of such a facet as reference outside
  its no-call blocks. A value change on files that already load (`blocks` becomes `reference`), and a
  new member of the `absence_means` enum (schema 0.0.11). CGI's native files keep
  `blocks`.

- `facets.dropped_fields` means the same in every loader: `[]` when the loader read the file's fields
  and dropped none, null only for a format with no notion of fields. The 23andMe-layout tables (all
  nine) and a Genos file without its `##Columns=` line wrote null; they now write `[]`. A value change
  on files that already load.

- `resolve_ref` splits Living DNA's run-together indel genotypes (`CCA`) against the reference where
  exactly one split fits (both pieces start with one base, one equals the reference over its span),
  and sets REF and ALT where only one allele is the reference. On the real file: 58 of 83 rows split, 4
  with REF and ALT; 25 stay unsplit, counted (F4 resolved; the rest is F17). Values change for
  `living_dna` facets (`gt_alleles`, and on those 4 rows `ref`, `alts`, `gt`, `end_pos`);
  `RESOLVE_VERSION` 4.

- A damaged BGZF input makes `load` raise `DamagedBgzfError` (exported) with the facet's id, the
  compressed byte offset of the first block that does not read back and what is wrong with it, in
  place of polars-bio's raw `ComputeError`; no path in it (F16). The damaged PGP-HMS starling upload
  now names offset 473,789,118. `just_vcf.bgzf.find_damage` is the block walker behind it.

- A VCF re-saved by a spreadsheet loads when a field is longer than 128 KiB: the re-save is undone
  without Python's `csv`, whose reader refused such a field with `csv.Error` (F7 follow-up).

- `normalize` left-aligns and trims the rows of a lifted facet (`relation = liftover_to`), whose REF
  `liftover` read from the target: an indel that a minus-strand block left at a repeat's right end
  moves to its left end, as `bcftools norm -f` on the target puts it. Before, lifted rows passed through
  as unchecked. Normalized facets get new ids (`NORMALIZE_VERSION` 3).

- A 23andMe-layout table that ends in 1,000 or more rows of one other width gets a warning that it
  looks like a second export appended, which is not read (F3 surfaced). Since resolved: an appended
  export of one allele per column loads as its own facet (the F3 entry below); a run of any other
  width keeps this warning.

- An FTDNA or MyHeritage-style export, or a Final Report that names no manifest, gets `facets.chip` by
  marker overlap with Illumina's manifests (RM4): `illumina_omniexpress` or `illumina_gsa`, basis
  `inferred`, only when one manifest holds at least 90% of the export's sampled rsids, 30 points above
  the other, and the export at least 40% of the manifest's; otherwise unknown with the shares in a
  warning. A sample of each manifest's rsids (1 in 61, 187 KB) ships as `just_vcf/chip_markers.json`,
  derived by `scripts/chip_markers.py`. On PGP-HMS, 57 of 103 FTDNA exports read `illumina_omniexpress`;
  the 43 X-chromosome exports, the 2010 Affymetrix exports and the one rsid CSV read unknown. No column
  changes.

- A 23andMe-layout table that ends in 1,000 or more rows of one allele per column (an AncestryDNA-shaped
  export appended after the first) loads as two facets in one `Loaded`: the first export as before,
  and the appended one read by its own layout (`0` is no call, chromosomes 23-26 are X, Y, XY, MT) with
  `relation = appended_export`, a new just-vcf relation (F3, resolved). The two are never merged. A run
  of any other width still gets the warning that it is not read. No column changes.

- A Final Report refused for want of positions or plus-strand alleles now names the manifest that
  would supply them, as its `Content` line states it (`HumanOmni2.5-8v1_A.bpm`); `chip` gains
  `final_report_content` for that line.

- A VCF re-saved by a spreadsheet (quoted fields with doubled quotes, tab-padded meta lines, bare-CR
  line ends) loads under `vcf`: the re-save is undone with `csv` and polars-bio reads a rewritten copy,
  with a warning (F7, now resolved). `vcf_samples` reads such a file too.

- A VCF whose column line is written `##CHROM` (two Veritas SV/CNV uploads, F7) loads: that line is
  read as the column line where no `#CHROM` exists, and `###` meta lines are kept without being read as
  declarations. An ALT of letters outside `ACGTN` (`CNV`) is no VCF allele: it gets no variant class
  and no callset site, and is counted in a warning (F15). Values change for files with such ALTs:
  `variant_classes` no longer counts them as indels or SNVs.

- A 23andMe-layout row with a doubled tab right after the id (`rs15842⇥⇥1⇥948921⇥CT`) is read without
  the empty field and counted in a warning, instead of being dropped (F2). Two 2012–2014 PGP-HMS exports
  gain 35 and 14 rows; an empty field anywhere else still drops the row.

- `liftover` lifts indels on minus-strand chain blocks (RM6): an indel whose alleles share their first
  base is reverse-complemented, that base dropped from every allele and the target's base before the
  lifted span put in front, so the row moves one base left; REF is read from the target over the moved
  span. Positions agree with pyliftover's image of the base after the source span. An indel that would
  land at a target contig's first base stays behind as `indel_minus_strand` (its meaning narrows); one
  whose alleles share no first base lifts complemented in place. Re-anchored indels are not
  left-aligned on the target. `LIFT_VERSION` 2: lifted facets get new ids.

- `write_vcf(…, blocks=True, reference=…)` writes a gVCF (RM6): the `blocks` frame as reference blocks
  (ALT `<*>`, or `<NON_REF>` with `block_allele=BlockAllele.NON_REF`; INFO `END`; FORMAT GT, GQ,
  MIN_DP and GL as the frame holds them), merged with the calls by position, contig by contig. A
  block's REF base is read from the collection `resolve_ref` read; a facet whose REF was not read is
  refused, and a block the collection cannot vouch for counts as `block_ref_unread`. `load` reads the
  output back to equal frames. `VcfWritten.blocks_written` is new, `blocks_not_written` counts only
  blocks left out, and `rows_written` counts calls only. Exports `BlockAllele`.

- `write_vcf` writes a tabix index, `<path>.tbi`, beside BGZF output (RM6), built from the offsets it
  wrote rather than by a second pass: `index=None` (the default) indexes every BGZF write, `index=False`
  skips it, and a plain-text write has none. A write without an index removes a `<path>.tbi` left
  from an earlier one. `VcfWritten.index_path` names it. htslib (pysam, opt-in test) answers region
  queries through it as the frame does, and builds the same linear index for the same file.
  `BgzfWriter` records its blocks' compressed offsets (`block_starts`, `virtual_offset`).

- Schema 0.0.10: `facets.chip`, `chip_basis` (array chip, stated or inferred) and `strand_check` (set by
  `resolve_ref`; a `mixed` facet withholds the alleles that read differently on the two strands).

- `resolve_ref` checks the strand of an array facet's genotypes (F13): `facets.strand_check` is
  `forward`, `mixed` (more than 2% of sampled hets carry the reference base's complement and not the
  base: part of the probes is written on the minus strand), or null (not an array, or fewer than 200
  hets). In a `mixed` facet the rows whose genotype reads differently on the two strands (homozygous
  calls, half-calls, hets without the reference base) keep `ref` and get `alts` and `gt` null, counted
  per reason; nothing is complemented. Before, such a facet's flipped hom-ref rows read as hom-alt.
  Values change on the one PGP-HMS export known to be mixed (an FTDNA file); exports `StrandCheck`.

- An array facet now says which genotyping chip it came from (RM4): `facets.chip` and
  `facets.chip_basis`. Stated where the header names it (AncestryDNA's array version, an Illumina
  Final Report's `Content` manifest as `illumina:<manifest>`, Living DNA's genotype chip); inferred for
  23andMe exports, and preamble-less tables carrying 23andMe's internal ids, when the row count and the
  internal-id profile both point at one generation (`23andme_v1_v2`, `_v3`, `_v4`, `_v5`); unknown
  otherwise, with a warning naming the features. FamilyTreeDNA and MyHeritage exports state nothing
  and stay unknown. Bands measured on 341 PGP-HMS 23andMe exports and checked against Illumina's GSA
  and OmniExpress manifests (docs/ROADMAP.md RM4). New columns are additive, nullable; exports
  `ChipBasis`. A 23andMe export too small for any band now carries one `chip unknown` warning.

- Loading a file again no longer touches the files an earlier `Loaded` of it reads, in this process
  or another. `load`'s spill (VCF, Complete Genomics, key=value) is a generation per content key under
  `frames/<facet_id>-<SCHEMA_VERSION>/`, like the derived steps': a reload sinks as before, finds the
  same bytes and keeps the first files. Before, it sank over `calls.parquet` and `blocks.parquet` in
  place. The manifest carries `SPILL_LAYOUT` 2, so the next `load` sweeps every load spill written
  before. The Complete Genomics and key=value loaders' working copy of the decompressed payload
  (`cgi/`, `keyvalue/`, a whole genome) is now an indexed entry (kind `work_copy`) under a lock: two
  loads of one file no longer write the copy at once, the copy goes even when a load fails, and
  `clear_cache(source)` removes a copy a killed load left, with the file's spills. Their spills had
  no manifest at all before, so `clear_cache(source)` missed them.

- `liftover` run again no longer breaks the `Lifted` an earlier run returned. It spills like the other
  derived steps: a generation per content key under `frames/<facet_id>-<SCHEMA_VERSION>/` (kind
  `lift_spill`, with a manifest naming the root source, so `clear_cache(source)` removes it), reused
  with its warnings and reasons when the same input comes again, a second one for other rows. Before,
  each run removed and rewrote `frames/<facet_id>-<SCHEMA_VERSION>-lift/`, which the index did not
  see at all; those old directories are swept by the next `liftover`.

- Schema 0.0.9: `facets.variant_classes`, `content_filter`, `content_filter_basis` (a genome file that
  kept one variant class is still `wgs`, but absence reads as reference only for the classes it carries;
  `absence_is_reference` answers it for consumers).

- A VCF facet now says which variant classes its rows carry and whether it was cut down to one:
  `facets.variant_classes` counts the ALT alleles per class (`snv`, `mnv`, `indel`, `sv`), and
  `content_filter` (`snv_only`, `indel_only`, `sv_only`) with `content_filter_basis` (`rows`,
  `stated`, `file_name`) names the cut. The rows decide first, then a bcftools or GATK SelectVariants
  command line in the header, then the file name's words (`.snp.`, `.indel.`, `.sv.`); a weaker piece
  that disagrees is a warning. The `.snp.` half of a split genome stays `callset_type = wgs` with
  `absence_means = weak_reference`, and the new public `absence_is_reference(facet_row, variant_class,
  force=False)` answers `True` or `None` (never `False`): `None` for a class the facet does not carry
  or that a filter cut, and for the kept class of a filtered facet unless `force=True`. Consumers that
  restore hom-ref from absence (just-prs, just-dna-lite) can gate on it; before, an SNV-only genome's
  absences read as weak reference for indels too. Three additive facet columns; the schema version
  is left for the merge.

- Running `rejoin_split`, `resolve_ref`, `normalize` or `mint_vrs_ids` again on a facet no longer
  breaks the frame an earlier run returned. Each run stages its spill and commits it as a generation
  named by the content it came from: the same input again reuses it (warnings included) and writes
  nothing, other rows under the same `facet_id` get a generation of their own, and two processes on
  one facet wait on the entry's lock and share one generation. Before, a second run removed and
  rewrote the files, and the first frame's next `collect` failed (`Invalid thrift`) when the rows
  differed or another process was mid-write. The derived spills' manifests carry `SPILL_LAYOUT` 2, so
  the next run of each step sweeps the entries written the old way. `load`'s own spill is unchanged.

- `normalize` drops `<*>` (and `<NON_REF>`) from the ALTs it splits: DeepVariant's variant records,
  `A>G,<*>`, become plain biallelic rows (no split flag) with `ad` and `gl` sliced exactly, and are
  then left-aligned like any other. Sites with two or more concrete ALTs beside `<*>` split as before,
  each row `split_signal = marker`. A site whose genotype calls `<*>`, and sites holding another
  symbolic allele, stay as written, counted. Before, every site holding `<*>` was left multi-allelic.

- `normalize` judges a duplicate on what the row states: the columns a move clears (`source_ref`,
  `ref_check`, `source_gt`, the VRS ids) no longer make two rows differ, so an indel the source wrote
  twice, once already left-aligned, with the same genotype and evidence, is one row (the unmoved one)
  rather than a reported conflict. Rows of one allele with other genotypes, evidence, filters or ids
  are all kept, with the warning, as before.

- A deletion written as ALT `-` (F12, about 8,080 rows in each of five PGP-HMS uploads) is now read as
  an unanchored deletion of the REF bases: the loader keeps it verbatim and counts it, `normalize`
  anchors it on the base before once REF has matched and then left-aligns it, and `write_vcf` refuses
  any ALT that is no VCF allele as `alt_not_an_allele` (a new `UnwritableRowsError` reason). Before,
  `write_vcf` wrote `-` and htslib refused the file. `NORMALIZE_VERSION` is now 2, which changes every
  normalized `facet_id`; the one bump also covers the dedup and `<*>` changes of the same day. No
  schema version change.

- `load` reads `write_vcf`'s own `INFO/JUST_VCF_SPLIT` back as a split marker: those rows get
  `split_signal = marker` and the facet's `split_marker` names the key, so a normalized facet written
  and reloaded keeps its split rows flagged and `rejoin_split` joins them again. Before, the tag was
  dropped as an unknown INFO key and the loader's patterns flagged some of the rows as `slot_pattern`.
  `writing.SPLIT_INFO_KEY` now lives in `loaders.vcf` (still importable from `writing`). No schema
  version change.

- PL now loads as the Float32 nearest `-PL/10`, divided in Float64 and rounded once, so a PL file and a
  GL file stating one likelihood compare equal on `gl` (F11). Complete Genomics' phred-unit GL is
  divided the same way. This changes `gl` by one Float32 step on some entries of files that already
  load. No schema version change (a rule change, no column).

- Cache eviction (RM268). Every cache entry now has a manifest beside it (`<name>.json`: kind, version
  keys, the root source's sha256 and path, the facet it holds and its parent, size, created and last
  used, ISO-8601 UTC) and a lock file its writer holds while it writes. `load` and each derived step
  (`rejoin_split`, `resolve_ref`, `normalize`, `mint_vrs_ids`) first remove the entries of their own
  kind an older version wrote, logging the bytes freed. New `clear_cache(source=None, older_than=None,
  dry_run=False)` removes every entry derived from a source (named by path, or by sha256 once the file
  is gone; derived spills included) or unused for longer than a duration, skips an entry being
  written, never removes an entry older than the index, and never touches the reference store; it
  returns a `CacheClearReport`. A manifest a newer version wrote is kept as unknown origin, never a
  reason for `load` to fail. Removing an entry a live `Loaded` reads breaks its next `collect`.
  `REWRITE_VERSION` 6 → 7, so every copy without a manifest is known stale and swept: the first load
  of each file after upgrading rewrites its copy once. Rules in docs/SCHEMA.md ("The cache under the
  lazy frames"). Measured on three PGP-HMS VCFs in a scratch cache: 4 entries, 48.1 MB; an old-version
  orphan of 34.8 MB was swept by the next load; `clear_cache(source)` removed one file's 2 entries,
  10.2 MB. Minor: new exports (`clear_cache`, `CacheClearReport`, `CacheEntryInfo`).

- Liftover as a derived facet (RM6, third slice): `liftover(loaded, chain, reference)` lifts a facet
  `resolve_ref` has read to another assembly through a UCSC chain file (`read_chain`, plain or gzipped)
  and returns a `Lifted`: the derived facet (`relation = liftover_to`, the Cat-VRS term) and a report
  with every parent row's target or its reason for staying behind. The chain's contig sizes are checked
  against both collections first (a primary contig of another size refuses it). A row's REF span must
  lie in one ungapped block of one chain; REF is read from the target, never carried over; minus-strand
  blocks reverse-complement the alleles; a target REF equal to an ALT swaps the two and re-indexes the
  genotype, one equal to no allele leaves the row unlifted. Indels on minus-strand blocks, symbolic
  alleles, pre-split rows and blocks are not lifted yet, each counted. Schema 0.0.8:
  `facets.derivation_inputs` names the chain by file name and sha256 (additive; spills from 0.0.7 are
  not reused). Minor: new exports (`liftover`, `read_chain`, `Lifted`, `Chain`).

- New `sex_chromosome_evidence(loaded)` (RM269): for each facet, what its calls say about X and Y, and
  an estimate of the karyotype. It measures GT ploidy as written per group (autosomes, non-PAR X,
  non-PAR Y, PAR, MT), the het share of non-PAR X calls written diploid with their AD allele balance,
  and the median DP of ALT calls on non-PAR X and on the X-degenerate Y (from GIAB v3.6
  stratifications, carried in `assets/`) over chr20's. It answers `x_copies` (one, two or more) and
  `y_present` (true, or unknown: call-site depth cannot show a Y absent), and from them a GA4GH
  Phenopackets `KaryotypicSex` value (`XY`, `XXY`, or `UNKNOWN_KARYOTYPE` with a reason), on GRCh37
  and GRCh38 whole-genome facets only; `XX` and `XO` wait for a coverage facet. X thresholds sit in the
  middle of the gaps the lead measured on 16 local genomes; checked on 27 PGP-HMS VCFs with
  self-reported sex, where every labelled man with a decided answer read `XY` and no woman read a Y
  (table in docs/ROADMAP.md RM269). A karyotype, not a sex, estimated from variant calls; not a
  cytogenetic result. Haploid genotypes and `ploidy_statements` are carried as the caller's assertion
  and never counted. Computed, not stored: no schema change. Minor: new exports
  (`sex_chromosome_evidence`, `SexChromosomeEvidence`, `KaryotypicSex`, `XCopies`).

- Normalization as a derived facet (RM6, second slice): `normalize(loaded, reference)` splits
  multi-allelic sites one row per ALT (the other ALT written `.`, each row `split_signal = marker`),
  left-aligns and trims indels whose REF matched the reference, drops exact duplicates (counted), keeps
  conflicting rows side by side with a warning, and sorts. `ad` is sliced, `gl` sliced only where the
  genotype holds no other ALT, `gp` withheld on split rows. Relation `normalized`; rules in
  docs/SCHEMA.md ("Normalizing"); no schema version change. Agrees with `bcftools norm -m-any
  --multi-overlaps . -f` on a synthetic file (sites, genotypes, AD, the GL it keeps), except that it
  drops the exact duplicates bcftools keeps. `rejoin_split` now leaves a group holding a row with no
  `source_ref` as written instead of failing on it. On six PGP-HMS uploads, chromosome 22: 235,911 of
  236,425 rows as bcftools writes them, the rest pre-split input left for `rejoin_split`. New opt-in
  test variable `JUST_VCF_TEST_REFERENCE_DIR` (`.env.template`). Found on the way: F11 (PL to `gl` one
  float step off) and F12 (uploads writing a deletion as ALT `-`).

- A symbolic SV's `end_pos` now reads INFO END: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>` and their subtypes
  end at END where the record gives one, else at POS + |SVLEN|, as SCHEMA.md's `end_pos` rule now says.
  Before, END was ignored there, so an SV with END and no SVLEN (a caller that writes none, and
  `write_vcf`'s own output) got `end_pos = pos`. An SVLEN declared `Number=1` is read too. This changes
  values on files that already load. END and SVLEN that disagree, and an END that is not read (`<INS>`,
  breakends, Delly's `<TRA>`), are counted in warnings. The per-row Python call that read SVLEN is now a
  Polars expression. No schema version change (a rule change, no column).

- GP within [0, 1] that does not sum to 1 on every row is now read as probabilities under a VCF 4.3+
  header, which defines GP so, with the rows counted in a warning (F10). Before, it was withheld as
  ambiguous, so just-vcf lost `gp` reloading its own `write_vcf` output of a phred GP. This changes
  values on files that already load: such 4.3+ files now get `gp`, and `source_alleles` on those rows.
  A 4.2 or older header, or none, still withholds it; a header without a version no longer counts as
  4.3 for the warnings. No schema version change (a rule change, no column).

- F9 closed after dogfooding on 153 PGP-HMS array exports: the row inference's floor is 60%, not 70%
  (FTDNA's X-only exports fit their build at 66.5–84.3%, every wrong build stays at or below 33.4%),
  and a VCF facet is not scored by its rows, so its refusal reads as before. F13 opened: an FTDNA
  export with part of its probes on the minus strand.

- A reference whose contig 1 is 247,249,719 bases is labelled `GRCh36`, so build-36 exports (23andMe
  v36, deCODEme) resolve against UCSC hg18 (where to get it: RM3). Like `GRCh37`, the label leaves an
  array's MT unread.

- `resolve_ref` infers an array's build from its own rows when the facet states none (F9): the share
  of sampled SNV rows whose genotype carries each collection's reference base names the build, with
  basis `inferred` and the scores in a warning, or the facet is refused with the scores. No column
  changed.

- VRS allele ids (RM271): `mint_vrs_ids(loaded, reference)` fills `calls.vrs_allele_ids` after
  `resolve_ref`, one VRS 2.x id per slot of `[ref] + alts`, equal to vrs-python's on every kind of
  allele tested and on 19,000 real PGP-HMS alleles; `vrs_withheld(calls)` says why an id is missing.
  Substitutions now wait for REF to be read, like indels (SCHEMA.md; no version bump, no column
  changed). vrs-python joins the dev dependencies as the test oracle.

- Schema 0.0.7 (RM3): `facets.ref_seqcol_digest` is new; `resolve_ref` fills `ref`, `ref_check`, and an
  array row's `alts` and `gt` from a reference.

- The CGI convergence tests (opt-in, PGP-HMS): one assembly's `var` and `masterVarBeta` give the same
  blocks and sites, and its `var` and `vcfBeta` agree on 4,138,184 sites, with 27,505 partial calls
  vcfBeta writes `.`, 356,513 loci vcfBeta writes as no-calls, and 1,665 (0.04%) left where vcfBeta
  joins neighbouring loci. vcfBeta covers every base `var` does not call reference: its silence is a
  reference call. `vcf_cgi` still says `absence_means = blocks`; the value that should replace it is
  open (RM2). The comparison found masterVarBeta writing two no-called stretches as `??`; it is now
  read as one `?`, as `var` joins them, so both files give the same blocks.

- Complete Genomics' native callsets load (RM2): `cgi_var` and `cgi_master_var` become one facet each,
  caller `CGAPipeline` at the file's `#SOFTWARE_VERSION`, `sample` its `#ASSEMBLY_ID`, `absence_means
  = blocks`, `callset_type = wgs` where `#TITLE` says whole-genome sequencing. Uniform loci become
  `ref` and `no_call` blocks; a variant locus becomes one call, its haplotypes joined from `var`'s
  sub-rows or read from `masterVarBeta`'s two allele columns. Indels stay unpadded (CGI writes no
  anchor base and REF only comes from the sequence): a deleted allele is `""`, and an insertion point
  sits at the base before it with `end_pos` null and `source_ref = ""`. Partial calls (`?A`) are kept
  verbatim; `phased`, `gl` and `gq` stay null (the allele scores are no likelihood, and `hapLink` is
  not read). No column changes, so `SCHEMA_VERSION` is unchanged.

- `just_vcf.table_samples(path)` lists the sample ids of a genotype table with a sample column (an
  Illumina Final Report of several samples, Color) in the order they first appear, the table
  counterpart of `vcf_samples`; `[]` for a table format without one. The refusal of a several-sample
  table now names it.

- `freebayes_varscan2_keyvalue` loads (RM2, SCHEMA Decision 8): PGP's `-orig` key=value table
  becomes one `caller_split` facet per caller, rows where that caller's genotype is not `NULL`, with
  `gt` over `[REF] + ALT`, the genotype text in `source_gt` and the caller's total depth in `dp`. Its
  VarScan2 facet keeps the genotype at sites both callers found, which a one-column merge gives to
  freeBayes. On the one pair on disk the table is a subset of its merge (2.04M of 5.30M rows, every
  chromosome, the gene-dense ones most): its freeBayes rows all match the merge's, and its VarScan2
  rows agree with the merge's wherever both have one. `callset_type` is inferred from the rows as for a
  VCF (RM4); that table's spread fits no class, so its facets say `absence_means = unknown`.

- Schema 0.0.6 (RM2's VCF side): `facets.ploidy_statements` is new, and the caller, its version, the
  vendor and the phasing method are filled from VCF headers. `REWRITE_VERSION` 6 and the schema bump
  orphan earlier cache entries (RM268).

- Privé's worldwide reference (RM5): `just_vcf.ancestry.build_prive_reference` joins the two public
  figshare files (group frequencies and 16 PC loadings, GRCh37) into one parquet, and
  `prive_proportions(loaded, reference_dir)` gives a facet's proportions over the 21 groups with
  block-bootstrap standard errors (200 blocks, 100 resamples, seeded). Ported from just-prs
  `ancestry/prive.py` (docs/PROVENANCE.md); a facet not on GRCh37 is unknown with a reason rather than
  lifted. On PGP-HMS data: a DRAGEN GRCh37 genome used 5,814,017 of 5,816,590 sites (3,426,634 of them
  absences read as reference) in 13 s, a 23andMe build 37 export 888,688; largest standard error 0.016
  in both. Which panels work, and how they were checked: docs/ANCESTRY_PANELS.md. New opt-in test
  variable `JUST_VCF_ANCESTRY_DIR` (`.env.template`). No schema change.

- Population inference behind the `[ancestry]` extra (RM5): `just_vcf.ancestry.infer_ancestry(loaded,
  panel_dir, panel)` places a facet on a reference-PCA panel (FRAPOSA's shrinkage-corrected projection,
  then a 20-neighbour vote and PC-centroid proportions), and `consensus` combines estimates over the five
  continental groups. Ported from just-prs `ancestry/` without its `PRSCatalog` wrappers or Hugging Face
  pull; a panel is a local directory, and one just-prs built reads unchanged (docs/PROVENANCE.md). The
  answer is an `AncestryEstimate`: the leading group, the vote, its margin, proportions and how many
  panel sites it rests on, split into observed and absent-read-as-reference; or None with a reason
  (`assembly_unknown`, `assembly_mismatch`, `panel_missing`, `too_few_sites`). The facet's own
  `assembly_label` picks the panel's assembly, and its `absence_means` decides whether an absent site
  is homozygous reference (`weak_reference`, or inside a `ref` block) or unknown; a no-call is never 0.
  `just_vcf` itself does not import it, so the core still imports without numpy. No schema change.

- `write_vcf(loaded, path)` writes a facet back out as a single-sample VCF 4.3, BGZF by default (RM6,
  first slice; new exports `write_vcf`, `VcfWritten`, `ContigNames`, `UnwritableRowsError`,
  `UnsortedRowsError`). It keeps GT with its per-haplotype phasing, GQ, DP, AD, GL, GP, DS, FT and PS,
  listing per row only the keys that row has; declares every FILTER, INFO and FORMAT key it writes;
  names the facet, its source digest and this package's version in the header, never the source
  path. REF is `ref`, else what the file wrote, never the genotype: a row with neither (an array before
  RM3), a REF mismatch, or no ALTs is refused with a count per reason, or skipped and counted with
  `resolved_only=True`. Rows RM270 flagged as split are written as held, tagged `INFO/JUST_VCF_SPLIT`
  and counted; blocks are counted and not written (a gVCF block needs a REF base the blocks frame does
  not hold before RM3). Contigs as stored, `chr`-prefixed (`MT` → `chrM`), or as the source spelled
  them. Two writes of the same frames give the same bytes. Round trips through `load`: a PGP-HMS VCF
  upload (956,064 calls) and a local DeepVariant genome (6,436,476 calls) came back with 0 key
  mismatches; htslib reads both outputs record for record, and the genome's output, reloaded and
  written again, gives the same body byte for byte. The PGP-HMS file's GP did not come back (F10).
  Minor: new functions and types only.

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
  `facets.ref_seqcol_digest` (schema 0.0.7: one new nullable column). Minor.

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
