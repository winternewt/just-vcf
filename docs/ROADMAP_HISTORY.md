# Roadmap history
Shipped roadmap items, moved here from [ROADMAP.md](ROADMAP.md) with their rationale. Nothing is
deleted; it is relocated.

## RM300 — Consensus reads no row of a facet for a class it carries none of
Severity: medium · Status: shipped, found 2026-11-23 by RM299's builder, not yet shown failing. Merge 70748f1: shown failing first (the array's row read as a contradiction moved a het insertion's posterior from 25 to 5 phred); such a row is now silent, `class_not_carried`, and CONSENSUS_VERSION is 6. Closed 2026-11-23.)

`consensus` has its own reading of a facet's one row at an allele's position over the same REF
(docstring point 2: "zero copies, `site`"), parallel to the one RM299 fixed in `concordance`. It very
likely reads an array's or an SNV-only file's row at the base an insertion is anchored on as evidence
of no copy, which a consensus would then weigh against the caller that saw the insertion. Where
`class_gate` says the facet carries no allele of the class, its row says nothing of that class.

## RM299 — A facet that carries no indels reads no copies of one
Severity: medium · Status: shipped, found 2026-11-23 by RM279's last pairs. Merge ec9925b; the RM279 pairs rebuilt the same day: the array pair's 15,869 indels read as called by both are one-sided, the DRAGEN SNV file's 4,009 are 1,318 both files wrote (docs/METHOD_ATLAS.md). Closed 2026-11-23. Left: consensus's parallel reading (RM300), a genotype table's classes (RM301), site_genotypes' (RM302).)

`concordance` reads the other facet's row at an allele's position when the row's REF is the allele's
(`single` in `_other_side`), and `_holds` then reads a haplotype equal to REF as holding no copy. That
is right for a whole-genome caller, which would have written the insertion in that row had it seen one.
It is wrong for a facet that carries no indels at all: an array export with no indel probes, or a
whole genome cut to its SNVs (`content_filter` `snv_only`). Its SNV row at the base an insertion is
anchored on says nothing of the insertion, yet the allele is counted as called by both and discordant.
On RM279's CGI `var` against an Omni2.5 array exported as VCF, all 15,869 indel alleles read so were
discordant, and none could have matched: the array file holds no indel row. CGI against a DRAGEN SNV
file read 4,009 such alleles, 78.7% discordant. Where the facet's absence says nothing of a class
(`absence_is_reference`: the class not in `variant_classes`, or cut by `content_filter`), its row
should say nothing of that class either: the allele stays one-sided, as it already is where that facet
has no row, and the reading is counted with its reason. Then the RM279 pairs' indel cells are rebuilt.

## RM279 — A method-error atlas
Severity: medium · Status: shipped (`method_atlas` and a first atlas over eight pairs shipped 2026-10-23, merge 4047968, docs/METHOD_ATLAS.md; the PGP CGI pairs rebuilt with REF read, minted and met by VRS id, merge 93dd5f0, with complex rows read by the allele's change (F21) and X's haploid calls matched to homozygotes (F22), merge c6da26c; counts per confidence band, merge 86daede, 2026-10-24; CGI `var` pairs banded by `calls.hq`, merge 8334138, 2026-10-24; the six private pairs rebuilt with minting through `scripts/method_atlas.py` (REF read, met by VRS id, banded), 2026-10-26, merged 2026-10-26. PGP-HMS pairs that pool more than one genome built 2026-11-08 (docs/METHOD_ATLAS.md, "More PGP pairs"): a third CGI `var` against freeBayes 0.9.15 person, the last on disk, and 16 people's two CGI sequencing runs, nine of them under one key (2.5.0.33 against 2.5.0.37) with easy-region SNV discordance from 0.13% to 1.84% between people, sent to desk #101. Single-person pairs of other methods (DRAGEN against freeBayes, CGI against DRAGEN, CGI against an Omni2.5 array as VCF) are on disk and planned in `data/interim/rm279/plan_pgp_methods.json`, not built: each adds a method with one genome, not a pooled row. Settled by default: pooling stays by caller and version, the `var` atlas a separate file (desk #101); the private folders are never paired (desk #72); CGI's sample-level VQLOW is not read by `carried_pass` (desk #102). The last single-person PGP pairs built 2026-11-23, three with rows, docs/METHOD_ATLAS.md, their indel cells left to RM299. Closed 2026-11-23: every planned pair built)

The person chose GIAB stratifications as regions (E8): error and discordance rates by caller and region from multi-facet PGP-HMS people (docs/LIVE_SAMPLES.md). Needs RM277.

Measured 2026-10-26, the six private pairs met by VRS id (docs/METHOD_ATLAS.md): every row of the six facets had its REF read and matched, every allele minted but GATK's 18,248 spanning deletions; the same-reads pairs' rates move by at most 0.02 point (F21's `unreadable`, 9 to 49 alleles per class); DeepVariant 1.1.0 spells about 14,000 indels unlike 1.3.0 and 1.8.0, which now meet, so its indel discordance rises from 7.1 to 8.0% and 8.0 to 8.7% and its one-sided indels fall from 8.8 to 6.0% and 9.2 to 6.3%; 1.3.0 against 1.8.0, 3.68 to 3.78% and 0.62 to 0.20%. SNV rates move by at most 0.01 point.

## RM280 — A consensus posterior
Severity: low · Status: shipped 2026-10-23, closed 2026-11-23 (`consensus` and `method_error_rates` shipped 2026-10-23, merge c1b0e1f, docs/SCHEMA.md "A consensus posterior", checked on PGP-HMS in docs/CONSENSUS.md; rates by the facet's own confidence band, merge 86daede, 2026-10-24, with the 60+ band traced to borrowed bands (9bf578e) and CGI `var` rows banded by `calls.hq` (8334138); desk #94–#96 built as options with every default unchanged, merge 960b6c6 (`BlockGrade`, GQ on reference, `AlleleFrequencyPrior`). Closed 2026-11-23: every question put to the person ended on its default: #94 blocks graded by their likelihoods, #95 GQ off for reference calls, #96 withdrawn 2026-11-16 (the prior an option), #99 bands borrowed across caller versions with a warning, #100 every GIAB class with the highest rate winning (the easy/difficult partition an argument), #101 pooling by caller and version, #102 FILTER only. The lead's calls listed for veto in docs/LAYERS.md stay listed there)

The person chose after the atlas (E7): a per-site posterior across one person's facets as a derived facet, with RM279's rates as priors; never before them.

## RM284 — Phasing as CisPhasedBlock
Severity: low · Status: shipped 2026-10-23, closed 2026-11-22 (`cis_phased_blocks` shipped 2026-10-23, merge 4129dc7; PSL, PSO and PSQ read and written back under VCF 4.5, with withheld blocks naming the allele that withheld them, merge 8b72a46, schema 0.0.15. Closed 2026-11-22: desk #90 ended on its default, so blocks stay whole phased stretches and one symbolic allele withholds its block; still outside the plan: a 4.4+ slot-0 block's phase on polars-bio keeping the leading phase indicator (P4 in docs/POLARS_BIO_PENDING_FIXES.md, not filed, the upstream-liaison seat's); the reading of a 4.4+ slot 0 that a PSL names as phased is listed for the person's veto in docs/LAYERS.md, layer 8)

The person chose columns with an export on demand (E5): phase sets stay columns; VRS CisPhasedBlock objects are minted only when asked.

**Built:** `cis_phased_blocks` (docs/SCHEMA.md, "CisPhasedBlock export"). The VCF loader reads PSL,
PSO and PSQ since 2026-10-23, and `write_vcf` writes them back under VCF 4.5, so PSL sets are tested on
a loaded facet against vrs-python. Left open: in a VCF 4.4+ file a slot-0 block with an ALT is
`phase_unknown` until polars-bio keeps the leading phase indicator (P4), unless PSL names that slot's set.
A withheld block names the first allele that withheld it (`withheld_chrom`, `withheld_pos`,
`withheld_allele`). On a whole imputed genome
an implicit set is a chromosome, so one symbolic ALT withholds the haplotype's block (a phased Genes
for Good upload: 45 of 46 withheld `symbolic`); whether a block per window or per contiguous run of
sequence alleles is wanted is the person's call.

## RM271 — VRS allele ids

Severity: high · Status: shipped 2026-10-22, closed 2026-11-21 (`mint_vrs_ids` and `vrs_withheld` shipped 2026-10-22 for VCF and array rows, merge b41bda5, equal to vrs-python on 19,000 real alleles. The opt-in test in `tests/test_vrs.py` that mints on the real store `JUST_VCF_TEST_REFERENCE_DIR` names and compares every slot with vrs-python (`scripts/vrs_oracle.py`) shipped 2026-10-26, merge 19a9753; measured below, 0 mismatch. Closed 2026-11-21: desk #77's term ended on its default, so the gtars gap G1 is not filed; it stays recorded in docs/GTARS_PENDING_FIXES.md and just-vcf keeps writing the VRS state itself)

Charter invariant 4: a call's identity is the VRS Allele id of its fully justified allele on a
refget-addressed sequence, and that id is the join key between facets on the same assembly. Today
`calls.vrs_allele_ids` is filled nowhere, though gtars (a core dependency) mints VRS ids and normalizes
alleles. It became possible with RM3: `resolve_ref` reads REF from a refget store and records the
collection on `facets.ref_seqcol_digest`, so each allele can be placed on a sequence with a refget
accession. An id that cannot be minted conformantly (an allele that cannot be justified, a contig with
no refget accession, an unresolved REF, an array I/D code) is withheld with a counted reason in
`vrs_errors`, never approximated. vrs-python is a test oracle only (CLAUDE.md §5).
**Shipped 2026-10-22: `mint_vrs_ids(loaded, reference)` and `vrs_withheld(calls)` (`just_vcf.vrs`).**
The rules are in docs/SCHEMA.md ("When an id is minted", "Why an id is missing"). Choices made, each
with the alternative rejected:
- **A step after `resolve_ref`, not inside it.** It refuses a facet whose REF was never read, and
  reads the collection `ref_seqcol_digest` names from the same store. Rejected: minting inside
  `resolve_ref`, which would make every REF read pay for ids, and leave no way to re-mint when the
  rule changes (`MINT_VERSION`) without reading REF again.
- **Every allele waits for `ref`, substitutions too.** Rejected: SCHEMA's earlier no-bases path for
  substitutions from a table of contig accessions per assembly. No table exists, and a contig
  `resolve_ref` could not match has no accession shown to be the facet's.
- **`[ref] + alts`, the REF slot minted** as VRS's reference allele (a `ReferenceLengthExpression`
  over REF's span), as SCHEMA and gnomAD have it. Rejected: one id per ALT only, which the brief
  named; it would retype a signed-off column.
- **The state choice and canonical JSON are ours; rolling and sha512t24u are gtars'.** gtars 0.10
  rolls an indel correctly but always returns a literal sequence, and its Python API digests only
  `LiteralSequenceExpression`; vrs-python gives a `ReferenceLengthExpression` for every deletion, every
  repeat insertion and the reference allele (VRS 2.x), so gtars' own ids differ from vrs-python's
  there. Rejected: withholding those alleles until gtars has it (most indels and every REF slot), and
  vrs-python in core (CLAUDE.md §5). Owed upstream: a VRS 2.x state in gtars' normalization and an RLE
  digest in its Python API; then `_justify` and the JSON shrink to calls. Not filed.
- **Rows walked as runs of one contig** (`rle` over the spilled input), reading only the columns
  minting needs; the ids go to their own files and sit beside the spilled rows by position. Rejected:
  a filtered scan per contig with a re-sort (resolve_ref's shape), one scan of the file per contig (170
  on the GRCh38 genome). *Amended 2026-10-26:* the runs are minted a group at a time, consecutive runs
  up to `SPILL_PART_ROWS` rows (`loading.bounded_groups`), one slice, one pass and one ids part per
  group, since a pass per run was a fixed cost per contig (F24, docs/PREVIOUS_ISSUES.md).
- **The oracle is `scripts/vrs_oracle.py`**: a vrs-python data proxy over the same gtars store, so
  both sides read the same sequence; vrs-python is a dev dependency.
**Dogfooded 2026-10-22** on three PGP-HMS VCFs, aggregates only (`scripts/vrs_conformance.py`), each
sample drawn with a fixed seed from the minted alleles and checked against vrs-python on the same
store:
- A GRCh38 genome, 4,651,432 rows: 9,256,478 of 9,345,929 alleles minted (363,409 deletions, 412,401
  insertions); 89,451 withheld as `ref_unchecked`, all on contigs `resolve_ref` left unread (UCSC
  `_random`/`chrUn_` names, RM3). 8,000 checked, 0 mismatch.
- A GRCh37 indel-only genome split, 956,064 rows, 1,039,271 ALTs: every allele minted; 9,000 checked
  (deletions, insertions, REF), 0 mismatch.
- A GRCh37 SNV-only genome split, 3,461,640 rows: all but 4 alleles minted (the 4 on rows whose source
  REF the reference contradicts); 2,000 checked, 0 mismatch.
- Time on the 4.65M-row genome: 59 s wall and 103–112 s CPU with the cache on a fast local disk; 112–208
  s wall with the cache under `data/interim/` on the NTFS disk, which `dd` measured at 2.5 MB/s during
  the runs (other sessions were writing to it; the workspace note says about 350 MB/s). Peak RSS of the
  minting process 1.5–2.0 GB; one contig's alleles and, for its indels, its sequence are in memory.
- *Measured 2026-10-26* (F24's amendment, library only, VRS spill fresh after a cached `load` and
  `resolve_ref`, every output hashing the same before and after): hu416394's bcftools gVCF, 9,651,143
  calls on 2,057 contigs, minted in 97.5–108.8 s a pass per contig, 70.5–75.0 s in 12 groups, at a
  sampled peak of 4.7 GB either way; its Manta file, 9,373 rows on 942 contigs, 18.6–18.8 s → 11.4 s. What
  remains is per allele (24.4 million gtars digest calls, the indels' justification) and 25 whole-contig
  fetches; docs/PREVIOUS_ISSUES.md, F24.
**Measured 2026-10-26**, the opt-in tests in `tests/test_vrs.py` on the store at
`/data/just-dna-cache/reference/store` (GRCh37 and GRCh38) and PGP-HMS uploads from the Krypton copy,
aggregates only; every expected id computed by vrs-python at test time:
- Synthetic alleles on real chromosome 22 (a 400 kb window, its repeats found at test time) and MT:
  SNVs, MNVs, insertions, deletions, delins, one and two units into homopolymer and di-, tri- and
  tetranucleotide repeats, part-unit deletions, a deletion anchored on the last leading N,
  multi-allelic rows, a lowercase ALT, a monomorphic row, MT written `MT` (GRCh37) and `chrM`
  (GRCh38), and a seeded random batch over the window and inside its repeats. GRCh37: 700 rows, 1,394
  alleles equal (697 REF, 120 SNV, 122 MNV, 273 insertions, 182 deletions); GRCh38: 694 rows, 1,382
  equal (691 REF, 118 SNV, 118 MNV, 278 insertions, 177 deletions). 11 alleles withheld on each, every
  one with SCHEMA's reason: `spanning_deletion`, `sv_no_vrs_model` (`<DEL>`, a breakend),
  `unspecified`, `not_sequence`, `ref_mismatch` (2), `ref_unchecked` (2 for REF past MT's end, 2 on
  a contig the collection lacks, `hs37d5` and `chrEBV`).
- A 23andMe-shaped build 37 export on real chromosome 22: hom-ref and het ids equal; `DI` is
  `array_code`, `--` is `ref_unresolved`.
- PGP-HMS, the smallest upload per build the store places, the first 3,000 chromosome 22 rows and up
  to 500 MT rows: GRCh37, 3,295 rows, 6,000 alleles equal (3,000 REF, 2,940 SNV, 27 insertions, 33
  deletions), the 590 MT alleles `ref_unchecked` (`resolve_ref` left that file's MT unread); GRCh38,
  3,032 rows, 6,113 equal (3,032 REF, 2,596 SNV, 233 insertions, 252 deletions), none withheld.
- Sensitivity: a deletion's `repeatSubunitLength` taken from the expanded REF instead of the
  deleted length fails both synthetic runs.
- Not reachable on the real store: `refused` and `no_contig`. After `resolve_ref`, REF was read from
  the very sequence minting rolls on, and neither collection holds a name twice; `_justify` is driven
  directly for `refused` in the synthetic suite.
## RM298 — The remaining readers take a leading BOM as their neighbours do
Severity: low · Status: shipped 2026-11-18 (merge 26ee0c9)

**Built:** a VCF opening with a BOM loads as the file without it, plain, gzip or BGZF, re-saved or not: the loader
reads past it and such a file always goes through the rewritten copy (`Header.bom`, `REWRITE_VERSION` 20; the copy
19 wrote put the header among the records). `rsid_table.column_line_ok` and `liftover.read_chain` decode with
`utf-8-sig`. Nirvana already read past it (`json.loads` on bytes); a test pins it. **Not upstream:** VCF 4.3 and 4.5
forbid a BOM, so polars-bio refusing one ("empty input") is conformant.

## RM297 — Coverage detection reads past a BOM
Severity: low · Status: shipped 2026-11-18 (merge ee8cb8a)
**Built:** `coverage._head_lines` drops a BOM from the first line, gzip or not, so six kinds detect as without it
(mosdepth summary, `samtools coverage`, `samtools depth` with `#CHROM`, thresholds, per-base and quantized BEDs
opening with a comment); a thresholds file with a BOM reads with `kind=` stated too.
## RM296 — An EIGENSTRAT `.ind` is read as its `.snp` is
Severity: low · Status: shipped 2026-11-18 (merge b5e1f94)

**Built:** `ancestry.eigenstrat._fields` reads both files' lines with `scan_lines`: a leading UTF-8 BOM dropped, a
line empty or whitespace alone skipped, fields split on runs of whitespace, a width error naming the line in the file.
`read_ind` and `read_snp` only pick their columns. A gzip `.ind` now reads as a `.snp` did (pinned by a test); a
bare-CR `.ind` is refused, not misread; invalid UTF-8 raises polars' `ComputeError`.
## RM295 — A blank line in an EIGENSTRAT `.snp` is no site
Severity: low · Status: shipped 2026-11-18 (merge b5e1f94)

**Built:** `read_snp` skips a blank or whitespace-only line, as `read_ind` did, so a `.snp` ending in one no longer
fails the panel build on the `.geno` header's site count; `AadrBuild.n_snps`/`n_individuals` and `read_genotypes` say
they count rows, not lines.

## RM294 — A text line holding 0x1F is read
Severity: low · Status: shipped 2026-11-18 (merge 47edcd6; full suite 2819 passed on polars 2.1.0 and 1.44.2)
**Built:** `loaders.text.lines` and `eigenstrat.read_snp` read lines with `scan_lines`; `lines()` decodes invalid UTF-8
with replacement first, the same U+FFFD per bad sequence as `utf8-lossy` gave. **Left:** `loaders/table.py` and
`consensus.py` still use 0x1F as a join key (a collision needs a field holding it).
## RM293 — Coverage's read and its row count agree on quotes
Severity: low · Status: shipped 2026-11-18 (merge 47edcd6)

**Built:** `coverage._read_fields` passes `quote_char=None`, so a double quote is text to the read as to its count of
cut and padded rows. A field written `"chr1"` keeps its quotes; a quoted number in a typed column fails the read.
## RM291 — Coverage rows cut or padded are counted
Severity: low · Status: shipped 2026-11-17 (merge 77cbd0b; full suite 2799 passed on polars 2.1.0 and on 1.44.2)

**Built:** `coverage._read_fields` counts each row's fields from the same bytes it reads (`scan_lines`,
`str.count_matches`), summed across blocks, and the facet's warnings say per reason how many rows were wider than the
format (fields past it not read), short of it (missing fields null), or blank (all-null rows); a samtools `depth` file of
several alignment files says the first depth column was read and how many there were. About 12–14% slower on a 90 MB
per-base BED. **Left:** the read honours `"` quoting while the count splits on raw tabs, so a quoted field holding a tab
would count as a wide row (no coverage tool writes quotes; `quote_char=None` would align them).
## RM292 — `read_bed` skips track lines
Severity: low · Status: shipped 2026-11-17 (merge 77cbd0b)

**Built:** `atlas.read_bed` reads three named fields with `schema=`, `truncate_ragged_lines` and `coverage._RAGGED`, no
quote character, and drops lines matching `^(track|browser)(\s|$)` wherever they stand; a BED of comments alone reads as
no intervals on both polars majors. GIAB reads are unchanged.

## RM289 — Coverage BED rows of mixed width on polars 2
Severity: low · Status: shipped 2026-11-17 (merge aac313d; full suite 2740 passed on polars 2.1.0, the touched tests on 1.44.2 too)
**Built:** coverage's two reads go through `_read_fields`: `schema=`, `truncate_ragged_lines=True`, and on polars 2.x
`extra_columns="ignore", missing_columns="insert"` (`coverage._RAGGED`), so fields past a format's columns are dropped
and missing ones null on both versions, at the parser's speed. A file of no rows (a header or comments alone) is a facet
of no intervals on both. **Left:** the cut and padded rows are not counted (RM291); `atlas.read_bed` fails on a `track`
line (RM292); `loaders.text.lines` and `eigenstrat.read_snp` refuse a line holding `\x1f`, alike on both versions.
## RM290 — An unclosed file in the VCF loader
Severity: low · Status: shipped 2026-11-17 (merge aac313d)

**Built:** the re-save undoing stream (`_Undone`) owns its text wrapper and detaches it on close, so a payload read to
its end is no longer closed under its opener when the wrapper is collected (the `ResourceWarning`, a closed handle as
well as a warning); `_vcf_text` closes the stream before the payload on every exit. The full suite shows no
`ResourceWarning`. `-W error::ResourceWarning` does not turn a collection-time warning into a test failure, so the test
asserts the handle and the warnings directly.

## RM288 — polars-bio 0.37 and polars 2.x
Severity: medium · Status: shipped 2026-11-16 (merge 0adf89c; full suite 2706 passed on polars 2.1.0, 2705 on 1.44.2)
polars-bio 0.37.0 (2026-11-07) allows polars `>=1.40,<2.2`, which removed the reason the lock stayed on polars 1.44.2.
Built with no question put to the person (CLAUDE.md §10, 2026-11-16: a choice with only pros is not a decision).

**Built:** the floor is `polars-bio>=0.37.0` with no direct polars pin, so a consumer on polars 1.4x co-installs; the
lock is polars-bio 0.37.0 and polars 2.1.0. The key=value and CGI `var` loaders split lines with
`loaders.text.tab_fields` (polars 2 refuses a headerless ragged `scan_csv`); VRS minting's horizontal `concat` passes
`strict` only to polars 1.x. CLAUDE.md §5 says, per polars 1.44 quirk, what 2.1.0 does (three fixed, two still so, one
default changed; the workarounds stay while 1.x is admitted) and names three polars 2 changes met. P16 is fixed
upstream, but the CLI keeps logging off for Python's last-resort handler, pinned by a subprocess test. **Left:** drop
the 1.x workarounds and `vrs._EQUAL_HEIGHTS` once the floor is polars 2; `coverage.py`'s `read_csv` with `new_columns`
would refuse a BED of mixed row widths on polars 2 (left as is); a `ResourceWarning` from `loaders/vcf.py` (an
unclosed file, on both versions).

## RM276 — A durable per-person store
Severity: medium · Status: shipped 2026-11-05 (the lead kept a migrated VCF facet's null `method_steps` and `gl_source`, told apart by the manifest's record and `loaded()`'s warning, docs/SCHEMA.md "Store migration"; was: in progress, decided by the person in the Part 2 interview, 2026-10-23; `Store` and `write_vcf_zarr` shipped 2026-10-23, merge 6830fb7; `call_HQ` in the Zarr export 2026-10-24, merge fbcb3a9; `variant_CIPOS`, `variant_CIEND` and `variant_CN` in the Zarr export 2026-10-26, an interval holding -1 or -2, the spec's integer sentinels, written missing and counted; 2026-10-26, `migrate_store`: a store of an older schema version rewritten to this one by a step per bump since 0.0.14 (`just_vcf.store_migration`), written beside, checked and swapped, what it cannot carry named as filled unknown, docs/SCHEMA.md "Store migration", merge 401cbf0. Nothing buildable left. For the lead: whether `method_steps` and `gl_source` null on a migrated VCF facet (not read) should get its own value instead of sharing null with "not a VCF" / "no likelihoods", which only the manifest and `loaded()`'s warnings tell apart now) · Owner: lead
The person chose both (E13): a directory per person holding parquet per facet, the links frame and a manifest, read lazily; and a VCF Zarr export (its spec has no gVCF blocks, so blocks stay parquet-only). Facets accumulate; nothing is overwritten (invariant 5).
**Built:** `Store` (docs/SCHEMA.md, "The durable store"): parquet per facet, links per facet, a
manifest, a lock, lazy reads into `Loaded` and `Person`; and `write_vcf_zarr`, one facet's calls as
VCF Zarr 0.4, blocks counted and left in the store. **Left:** a migration between schema versions (a
store of another version is refused, not migrated); `cipos`, `ciend` and `alt_copy_numbers` in the
Zarr export; whether `Person` should carry the links (`store.loaded().links` has them).
*Corrected 2026-10-26:* `cipos`, `ciend` and `alt_copy_numbers` are written (docs/SCHEMA.md, "VCF
Zarr export").
*Corrected 2026-10-26:* `Person` carries `links` (merge b8bef07, 2026-10-23, "Person and the CLI keep
links"), and `Store.person()` builds it from `loaded()`, so its links are the store's; that question is
closed.
*Corrected 2026-10-26:* the migration is built: `migrate_store(path)` rewrites an older store by the
steps in `just_vcf.store_migration` and reports what it filled unknown (docs/SCHEMA.md, "Store
migration").
## RM277 — Concordance between facets
Severity: medium · Status: shipped 2026-10-23 (merge 2446c53; F21 and F22 followed the same day, merge c6da26c; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose graded where possible (E2): GL/GP compared where both facets carry them, hard genotypes otherwise, and every row says which. Kept as data, never resolved to a winner (invariant 8). Needs RM273 and VRS ids (RM271).

## RM273 — A person (subject) in the schema
Severity: high · Status: shipped 2026-10-23 (subject and Person, merge 93745c2; RM274's `same_person` as the default check, merge 6dd3591; was: in progress, subject and Person shipped 2026-10-23, the default check RM274's `same_person`) · Owner: lead
The person chose "Both: supplied, checked" (E1): the caller supplies a subject id when loading; nothing is inferred and facets are never regrouped by the library. The library checks the claim with relatedness (RM274) and warns on a mismatch. Most of the extended shape sits on this: RM275, RM276, RM277, RM280.
## RM285 — The format registry as a spec
Severity: low · Status: shipped 2026-10-23 (merge d0edd91; was: built 2026-10-23, docs/FORMAT_REGISTRY.md generated by scripts/format_registry.py, to move here on merge) · Owner: lead

The person chose this repo (E14): the declared formats and the broken shapes (F1–F17) documented as a readable spec in docs/ with synthetic examples; no genome data.

## RM272 — CLI
Severity: high · Status: shipped 2026-10-23 (merge 5b8700f; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead
Typer, registered in `[project.scripts]`, loading `.env` with python-dotenv (CLAUDE.md §4): `just-vcf sniff`, `load` (a summary of the facets), `convert` (to VCF/gVCF via `write_vcf`, optionally through `resolve_ref`, `normalize`, `liftover`, `mint_vrs_ids`) and `check`. Output is counts and aggregates only, never rows: the CLI is a privacy boundary like E12.
## RM275 — The facet graph as a links frame
Severity: medium · Status: shipped 2026-10-23 (merge c503f57, schema 0.0.13; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose a small links frame (E4): (from, to, relation, provenance), so a facet can have several parents (a merge, a consensus). `derived_from` and `relation` stay on `facets` for the one-parent case. Relations verbatim from GA4GH where it has them (`liftover_to`), just-vcf's own otherwise, documented.

## RM283 — An aggregate-only query helper
Severity: medium · Status: shipped 2026-10-23 (merge 95d8b61; was: built on a builder branch 2026-10-23, `just_vcf.query`, waiting for the lead's merge; decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead
The person chose a helper beside `load` (E12): an API that answers counts, shares and presence questions and never returns rows, for agents and apps; `load` stays for code that may read rows.
## RM268 — Cache eviction
Severity: medium · Status: shipped 2026-10-22: the cache index, the stale-version sweep, `clear_cache`; no size cap · Owner: lead

Nothing ever removes anything from `JUST_VCF_CACHE_DIR`. It holds two kinds of entry:

- `vcf/<key>/`: rewritten copies of inputs polars-bio cannot read as they are
  (`loaders/vcf.py` `readable_path`), plus a `structural/` copy for a file polars-bio refused. The key
  hashes the source digest, the archive member and `REWRITE_VERSION`.
- `frames/<facet_id>-<SCHEMA_VERSION>/`: the calls and blocks parquet spill from `load`, and since
  RM270 the calls spill of a `rejoin_split` facet, under the derived `facet_id`. Nothing in the cache
  ties that id to its parent, so an eviction that keys on a source file does not see it.
Because the version is part of each key, every bump of `REWRITE_VERSION` or `SCHEMA_VERSION` orphans
every earlier entry, and nothing collects them. Measured 2026-10-03: after the 2026-09-30 dogfooding
runs, the cache under `data/interim/cache` had reached 32 GB in 36 `vcf/` entries. Most of it was
orphans: the two largest were plain `input.vcf` copies of 14.5 GB and 7.4 GB, from before the rewrite
switched to BGZF `input.vcf.gz`. The user deleted it by hand.
Constraints on a fix:
- A copy has to outlive the lazy frames that read it, so evicting the entry a live `Loaded` points at
  breaks a later `collect`. Orphans by version are always safe to remove; current entries need an age
  or size policy.
- The copies are whole private genomes. Deleting a source under `/data/genomes` leaves its copy in the
  cache, so eviction is also how a removed genome actually goes away.
Options, not decided: sweep stale-version entries on each load; an explicit `clear_cache()` in the
API (and in the CLI once it exists); a size cap with least-recently-used eviction.

**Shipped 2026-10-22** (the lead's outline, 2026-10-22). Rules in docs/SCHEMA.md ("The cache under the
lazy frames"); code in `just_vcf/cache.py`.

1. **The index.** Every entry has `<name>.json` beside it: kind (`rewrite_copy`, `frames_spill`,
   `rejoin_spill`, `resolve_spill`, `vrs_spill`, `normalize_spill`), version keys, the root source's
   sha256 and resolved path, `facet_id` and `parent_facet_id`, size, `created` and `last_used`
   (ISO-8601 UTC). A derived spill takes its root source from the facet row, which every derived step
   carries down (`source_sha256`, `source_path`), so `rejoin_split`'s spill is no longer cut off from
   its file.
2. **Stale versions are swept** by the step that owns the kind, before it writes, with the bytes freed
   logged. An entry from before the index is judged by its name: a frames name carries
   `SCHEMA_VERSION` (and the resolve or mint version); a `vcf/` key hashes its version, so
   `REWRITE_VERSION` went to 7 with the index and a copy without a manifest is stale by construction.
3. **`clear_cache(source=None, older_than=None, dry_run=False)`** returns a `CacheClearReport`
   (`removed`, `being_written`, `unknown_origin`, `bytes_removed`; no source paths, so it can be shown).
   Writers hold an OS file lock (`flock`; `msvcrt` on Windows, unexercised) beside the entry, and
   `clear_cache` skips an entry whose lock is held. A waiting writer re-checks that the lock file is
   still the one at its path, since a remover unlinks it.
4. **No size cap.** Rejected: LRU under a cap. `last_used` records writes and reuses, not reads, so a
   cap could evict the files a live frame scans with nobody asking.

Decisions taken while building, each with the rejected alternative:

- The manifest sits beside the entry, not inside it: `resolve_ref`, `normalize` and `mint_vrs_ids`
  remove their spill directory before writing it, which would take a manifest or lock inside with it.
- Bump `REWRITE_VERSION` rather than keep unindexed copies until someone deletes them by hand: the
  first load of each file after the upgrade rewrites its copy once, and in exchange no whole-genome
  copy stays behind unattributed.
- `source` resolves as: an existing file is hashed and matches by sha256 or by recorded path; a
  64-character hex string that is no file is the sha256; any other path matches the recorded path (a
  deleted file). Rejected: sha256 only, since the person who deleted a genome knows its path, not its
  digest.
- `older_than` reads `last_used`, not `created`: a copy reused every day is not old.
- `clear_cache()` with neither filter clears every indexed entry, with `dry_run` as the valve.
  Rejected: refusing without a filter, which would make clearing the cache take
  `older_than=timedelta(0)`.
- The manifest is written when the lock is taken as well as when the write ends, so an entry whose
  writer died still names its source and can be cleared by it.

Measured 2026-10-22 on three PGP-HMS VCFs (one that needs a rewrite copy) in a scratch cache under the
builder's worktree: 4 entries, 48,099,333 bytes (3 load spills, 1 copy). A copy of one spill under an
older schema name, without a manifest, 34,780,557 bytes, was swept by the next `load`.
`clear_cache(source)` on the file with the copy removed its 2 entries, 10,202,328 bytes, dry run and
real run agreeing; `older_than=timedelta(days=1)` removed nothing; `clear_cache()` removed the last 2.

Open: the derived steps still remove and rewrite their spill in place, so running one twice on the
same facet breaks the frame the first run returned. That predates RM268 and is not changed here.
Closed later the same day: a derived spill is now one generation per content key (docs/SCHEMA.md,
"The cache under the lazy frames"). `liftover`'s spill and `load`'s own followed the same day, and the
Complete Genomics and key=value loaders' working copies came under the index (`work_copy`).
Amended 2026-10-23 (the lead's decision, after two sessions at two `SCHEMA_VERSION`s on one cache dir
broke the older one's live frame): "orphans by version are always safe to remove" was wrong while an
older install still runs. The sweep now keeps an older entry used within `JUST_VCF_CACHE_GRACE_DAYS`
(7 days) or whose lock is held (docs/SCHEMA.md, "The cache under the lazy frames").

## RM270 — Pre-split multi-allelic input
Severity: high · Status: shipped 2026-10-19: B, the flags (schema 0.0.4), and D, `rejoin_split`; C not built, the person chose D (desk #45) · Owner: lead
Motivating report: S2 (just-dna-lite, 2026-10-18). Many VCFs arrive already split one row per ALT, in
the `bcftools norm -m-` convention, and their headers do not say so. In a split row the `0` of a
genotype can stand for "another ALT", not REF: a `1/2` at `A>G,T` becomes `1/0` at `A>G` and `0/1` at
`A>T`, and a `1/1` becomes `1/1` plus a `0/0` row for the uncarried ALT.

Reproduced 2026-10-18 against the code, on a synthetic VCF: the loader keeps one `calls` row per record
(SCHEMA Decision 1, A), so the split rows land as written, and `gt_alleles` maps each `0` onto
`source_ref`. The split `1/2` loads as `["G","A"]` and `["A","T"]`, two calls with a reference copy the
person does not have; the uncarried ALT's `0/0` row loads as `["C","C"]`, a hom-ref call. `warnings` is
empty. The unsplit record of the same site loads as `["G","T"]`.

**Prior art** (read from source and docs, 2026-10-18). bcftools `norm -m-` writes the other ALT as `0`
by default, keeping GT slots, and `.` with `--multi-overlaps .` (1.17+); its `--old-rec-tag` INFO
(1.12+) is opt-in. Hail's `split_multi_hts` downcodes to `0` and sorts unphased calls, so its split
`1/2` is `0/1` + `0/1`, and its `a_index`/`was_split` are row fields an exported VCF lacks. GATK
`LeftAlignAndTrimVariants --split-multi-allelics` writes `0`, except that every sample becomes `./.` when
any sample is `1/2`. vt `decompose`, Glow and vcflib `vcfbreakmulti` write `.`, and vt and Glow add INFO
`OLD_MULTIALLELIC`. The VCF spec says nothing about split records. `bcftools norm -m+` re-joins slot by
slot (`1/0` + `0/1` gives `1/2`). Rows written with `.` already load as half-calls and are safe.
**Detection**, three signals per `pos` within a contig, reported separately and never folded into one
flag: a **marker** (INFO `OLD_MULTIALLELIC`, or a value of the `--old-rec-tag` shape
`CHROM|POS|REF|ALT|IDX`); the **slot pattern** (two or more non-reference rows, each GT slot non-zero in
at most one); and an **all-zero row of any ploidy beside a non-reference row**. On the facet, a marker
is Kleene: true when seen, null otherwise, never false, because a split file can show no signal at all
(Hail's sorted output, split indels that left-alignment moved to another POS). Measured 2026-10-18 on
four unsplit private WGS VCFs (DeepVariant ×2, DRAGEN, GATK HaplotypeCaller), 22.2 million calls: no
signal fired. Only the DRAGEN file had positions with two records (10), none matching. The rate on a
real split file is not measured (no bcftools here); a synthetic split VCF is caught at both its sites.
The two unsplit DeepVariant gVCFs gave no signal either: no two call rows share a position once their
reference blocks go to `blocks`, which never take part.
**What the frame says about a `0` it cannot place.** No option changes `gt`, the row grain, or
`source_gt` where it is kept (null on VCF 4.4+, P4). Not decided:

| | What the frame says | Legality | Cost |
|---|---|---|---|
| A. Warn only | a warning per signal, with a count | patch | every flagged row stays wrong for a reader that skips warnings |
| B. A + flags | `facets.split_marker` (Kleene) and a nullable row column `split_signal` naming the signal | minor | values stay wrong, but every flagged row can be found |
| C. B + slot-aware withhold | in a flagged group, a `0` becomes null in `gt_alleles` only in a slot where another row of the group is non-zero; a slot zero in every row keeps its REF | values read today change (pre-1.0, changelog says so) | a false positive loses a true REF copy, never states a false one; amends Decision 2: `gt` and `gt_alleles` agree except where withheld, and a null means uncalled or withheld |
| D. Re-join as a derived facet | each flagged group merged into one site slot by slot, as `bcftools -m+` does (`relation = normalized`, charter 5) | minor | the most work; GT only, since bcftools subsets PL/AD to REF + one ALT; phasing across slots may be lost |
`gl`, `gp` and `ad` on a flagged row stay as written: they compare REF with one ALT and inherit the
ambiguity, and the row flag is what tells a reader so. The lead recommended B and C now and D inside
RM6 (desk #45). **Decided 2026-10-19**, the person: "4 if it can be implemented w/o bam, else 3". It
can: the re-join needs only the split rows, as `bcftools norm -m+` does; reads would only rebuild the
full likelihood vectors, which is variant calling (charter). So D is built now, GT only, beside B;
C is not built.

**RM6's own split output** must not write this convention: it writes `.` for the other ALT, as vt,
Glow, vcflib and `bcftools --multi-overlaps .` do.
Legality: detection, a warning and a new nullable field are additive (minor). Nulling `gt_alleles` for
detected rows changes the meaning of values already emitted.

RM6 carries point 3: its own split output writes `.` for the other ALT.

## RM8 — Output frame schema
Severity: high · Status: signed off 2026-09-30 as schema v0.0.1 · Owner: user decides, agent drafts
Before any per-format loader is written, draft a proposal for the one frame every loader emits, and get
the user's sign-off. The charter (docs/CHARTER.md) is still an uncommitted draft, so the schema cannot
be derived from it alone.
The proposal must settle:
- **Columns and dtypes** of the minimal shape (many formats in, one frame out), and which of them the
  extended shape (facets per assembly, caller, phasing) adds.
- **Coordinate naming.** 1-based `POS` at every public boundary; name each column for the convention it
  holds so 0-based sources (CGI, Color) cannot leak through unconverted.
- **Graded genotype evidence.** The user's words: a fuzzy scale from 0 (unknown) to PHRED-like
  (confident), probably an n-vector of probabilities or log-likelihoods per position; never a binary
  called/not-called. Say how array calls, imputed dosages (Oxford `.gen`, Genes for Good imputed),
  half-calls (`1|.`) and CGI partial no-calls (`GC?`) map onto it.
- **Absence semantics.** Where an absent site means unknown (arrays, imputed content, VCFs of array
  data) versus weak evidence of reference (whole-genome callsets); likely a per-facet property.
- **Provenance per row:** source file, format id, caller or vendor (the `source` detector in RM2),
  build as stated versus inferred.
- **GA4GH fit:** VRS allele ids as the join key (withheld when non-conformant), seqcol for assembly
  identity.
- **Multi-caller files:** PGP merges carry the caller in INFO/SF, not per sample column
  (docs/PGP_HMS.md).
Inputs: the "Contract gaps" list at the end of docs/FORMATS.md, RM2, docs/POLARS_BIO_PENDING_FIXES.md
(P1: undeclared INFO/FORMAT keys are dropped), and the charter draft's invariants.

Signed off 2026-09-30 as docs/SCHEMA.md, schema v0.0.1. The user's words: "lets keep schema as it is
now and switch to building - we'll find out if some pieces are missing in progress". Gaps found while
building the loaders are fixed in SCHEMA.md as they come up, each with a version bump.

## RM7 — PGP-HMS sample leecher
Severity: medium · Status: shipped 2026-09-28 · Owner: Claude

A script in `scripts/` that downloads consented-public 23andMe exports from
<https://my.pgp-hms.org/public_genetic_data?data_type=23andMe> into `data/input/`, so tests and
dogfooding have real array files that are allowed to be used. httpx is the dev dependency for it.

Shipped as `scripts/pgp_leech.py`.
Download links answer two ways: older uploads are served directly, newer ones redirect to an Arvados
collection page. The listing's participant cell can carry the original enrolment number
(`huF06AD0, PGP113`), which the live test caught. Rows over 200 MB are skipped by default: the
23andMe listing totals about 101 GB, and 19 of its files are multi-GB uploads.
