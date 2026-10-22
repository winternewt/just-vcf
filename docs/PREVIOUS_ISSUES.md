# Previous issues — resolved here

Findings from [DOGFOODING.md](DOGFOODING.md) resolved in this repo, each with its resolution and a
code pointer. Check here before re-investigating something that looks fixed.

## F7 — Two VCF shapes `sniff` accepts but the loader cannot read

Severity: low · Found 2026-09-30, same sweep.

A VCF re-saved by a spreadsheet (tab-padded, quoted `"##` lines) and Veritas SV/CNV VCFs whose every
header line carries an extra `#` (`###contig`, `##CHROM`) are recognised as `vcf` but fail with "no
#CHROM line". The Veritas shape could be read by accepting `##CHROM` as the column line; the spreadsheet
one needs its own format id, as `23andme_csv` has. Not fixed: three files.

**Veritas shape fixed 2026-10-22.** Not every header line has the extra `#`: `##fileformat` and `##ALT`
are written right, `###FORMAT`, `###contig` and `###reference` carry a third, and the column line is
`##CHROM`. Where no `#CHROM` line exists, a meta line that is VCF's column line with one `#` too many
is now the column line, and the file is read through a rewritten copy; the `###` lines are kept as
the meta lines VCF says they are (a key starting with `#`), so no contig and no FORMAT is declared from
them, and GT and END are declared by the loader as for any undeclared key. Both files load (5,614 SV
calls, all `0/1`; 2,632 CNV rows, all `0/1`, whose ALT `CNV` is F15). Code: `loaders/vcf.py`
`read_header`, `HASHED_COLUMN_LINE`; test `test_a_column_line_written_with_two_hashes_is_read`.


**Spreadsheet shape fixed 2026-10-22; F7 resolved.** The re-save is CSV quoting applied to a tab file:
a field holding a comma or a quote is wrapped in quotes with its quotes doubled (115,123 of 115,220
rows, their INFO holding SnpEff's `EFF=…,…` or their sample field an AD of `12,0`; the 97 unquoted
rows are `./.`), meta lines are padded with tabs to ten fields, and lines end in a bare CR, which is
why `read_header` found no `#CHROM` line in 64 MiB. Undoing that is exact, so the rejected option was
the separate format id: the delimiter did not change, as it did for `23andme_csv`, and the quoting
and CRs are a container's mangling, like CRLF, which `sniff` already reads through. The loader looks
at the payload's first 64 KiB (a line starting `"##`, or CRs with no LF); a re-saved file is read
through `csv` with universal newlines, trailing empty fields dropped from `##` lines, and polars-bio
reads the rewritten copy; a warning says so, and that a number whose text the spreadsheet rewrote
(GATK's `6980.90` as `6980.9`) cannot be told. The upload loads all 115,220 rows (a GATK
UnifiedGenotyper exome, `wes`, GRCh37 from its 25 contigs), its genotype tally equal to the raw
text's. Code: `loaders/vcf.py` `is_resaved`, `undo_resave`; `loading.py` `_vcf_text`; test
`test_a_vcf_resaved_by_a_spreadsheet_loads_as_the_text_it_was_saved_from`.

The third failure in the sweep, a corrupt upload (`gzip -t` fails too), keeps its correct error.

**Follow-up 2026-10-22: a field over 128 KiB.** The re-save was undone with Python's `csv`, whose
reader refuses a field over `csv.field_size_limit()` (131,072 bytes by default) with `csv.Error`, and
a VCF field has no such limit (SnpEff's EFF on a busy site, quoted for its commas). The limit is
process-wide, and the reader is consumed lazily while the copy is written, so it cannot be raised for
this reader alone. The re-save is now undone without `csv`: lines by universal newlines, cut at tabs,
a field wrapped in quotes unwrapped with its doubled quotes made single. A VCF field holds no tab and
no line break, so nothing carries across fields or lines. Rejected: raising the process-wide limit,
which widens every other `csv` reader in the caller's process and removes the guard against an
unclosed quote swallowing the rest of a file. Code: `loaders/vcf.py` `undo_resave`, `_unquote`; tests
`test_a_resaved_vcf_with_a_field_past_the_csv_limit_loads` (about 200 KiB of EFF) and
`test_undo_resave_reads_quotes_as_csv_does` (the same output as `csv` on the shapes it reads).

## F15 — An ALT that is no VCF allele counted as an indel

Severity: low · Found 2026-10-22 loading the Veritas CNV upload F7's fix made readable.

The file writes every ALT as the letters `CNV`, without angle brackets. The class rule took any run of
letters for a sequence allele, so the facet's `variant_classes` said 2,632 indels and the callset
inference counted the rows as small-variant sites. VCF allows only `ACGTN` in a sequence allele.

Resolved 2026-10-22: a sequence ALT is `ACGTN` bases only, in `variant_class_expr` and in the
callset's site count, so `CNV` (or an ambiguity code such as `R`) has no class and is no site; the
loader keeps it verbatim and counts the records in one warning. Code: `just_vcf/callset.py`
`variant_class_expr` and `spread`, `loaders/vcf.py` `scan_stats`; test
`test_an_alt_that_is_no_vcf_allele_has_no_variant_class`.

## F2 — 23andMe rows with a doubled tab are dropped

Severity: low · Found 2026-09-30 by loading the 503 PGP-HMS 23andMe-layout files.

Two 2012–2014 23andMe exports carry rows with an empty field between the rsid and the chromosome
(`rs15842⇥⇥1⇥948921⇥CT`): 14 and 35 rows out of about 960,000. The loader drops them and says so in
`Loaded.warnings`. Not fixed: accepting "exactly one empty field after the id" is a guess about what the
vendor meant, and the rows are few. Revisit if another vendor shows the same shape.

Resolved 2026-10-22. The rows were checked rather than the vendor's intent guessed: in both exports
every such row (49 distinct ids, all rsids of four or five digits) has its one empty field right after
the id, and each of the 49 sits at the same chromosome and position in other PGP-HMS 23andMe exports
(the first 120 checked), with none placed elsewhere. It is not tab-stop alignment: most ids of the
same length are written with one tab. So a row of exactly five fields whose second field is empty is
read without that field, and counted in its own warning line; an empty field anywhere else still drops
the row (the stray first row of the F3 file, empty between position and genotype, stays dropped). Both
exports now load 960,613 rows, the same count. Code: `just_vcf/loaders/rsid_table.py` `read_calls`;
test `test_a_doubled_tab_after_the_id_is_read_past_and_counted`.

## F12 — Some VCF uploads write a deletion as ALT `-`

Severity: low · Found 2026-10-22 running `normalize` against `bcftools norm` on PGP-HMS uploads (RM6).

Five of 25 PGP-HMS VCF uploads of 5–60 MB hold about 8,080 rows each whose ALT is `-`, the dash
notation for a deletion that VCF does not have. The loader keeps it verbatim as an allele, `normalize`
leaves those rows as they are and counts them, and `write_vcf` writes them back out, so the export is a
file htslib refuses ("Non-ACGTN alternate allele"). Nothing in the frame says the allele is not one.

Not fixed here. Reading `-` as a deletion needs the anchor base and its REF span, which the row does
not state, so it is not a spelling to rewrite; the options are for the loader (RM2: flag or refuse
such rows, counted) and the writer (refuse them as unwritable, as it does a row without REF).

Resolved 2026-10-22 (RM6). A survey of the five uploads settled the reading: every such row has one
ALT, `-`, and a REF of one to eight or more bases (about 4,680 of one base, 1,800 of two), and REF
matches the GRCh37 reference at POS on all 40,394 of them, so REF is the deleted span and only the
anchor base is missing. Almost all are `0/0`. Now:

- the loader keeps `-` verbatim and counts the rows in one warning that names them unanchored
  deletions;
- `normalize` reads `-` as the empty allele once `resolve_ref` has matched REF, and the left-align loop
  anchors it on the base before (`CT>-` at 30 becomes `GCT>G` at 29), then aligns it like any
  deletion. At a contig's first base there is no base before; such a row stays as written, counted
  (and the retry loop now stops there instead of widening its window for ever);
- `write_vcf` refuses any ALT that is no VCF allele as `alt_not_an_allele`, counted, or skips it under
  `resolved_only=True`.

Anchoring in `resolve_ref` was set aside, since it fills slots of the same facet and moves nothing;
refusing the rows in the loader was set aside, since their REF matches everywhere it was checked. On
the five uploads: 40,392 rows anchored; the other 2 are rows the loader flagged as one record of a
pre-split site, which `normalize` leaves for `rejoin_split` and `write_vcf` skips as
`alt_not_an_allele`. htslib reads every row of the five exports. Code: `loaders/vcf.py`
`DASH_DELETION`, `normalize.py` `_align`, `writing.py` `_unwritable`; tests:
`tests/test_load_vcf.py` `test_a_dash_alt_is_kept_verbatim_and_counted`, `tests/test_normalize.py`
`test_a_dash_deletion_is_anchored_on_the_base_before`, `tests/test_write_vcf.py`
`test_an_alt_that_is_no_vcf_allele_is_never_written`.

## F11 — PL to `gl` lands one float step off the nearest value

Severity: low · Found 2026-10-22 by comparing `normalize` with `bcftools norm` (RM6).

The loader turns PL into `gl` as `PL / -10` in Float32, and Polars computes that division so that PL 21
gives -2.1000001, while the nearest Float32 to -2.1 is -2.0999999. A file that writes the same
likelihood as GL (`-2.1`) therefore loads one step away from one that writes PL 21. Deterministic, so no
digest or test of ours moves, and the size is far below any PHRED rounding; but two facets of one
genome, one with PL and one with GL, do not compare equal on `gl`, and a test comparing against an
outside tool needs a tolerance (`test_agrees_with_bcftools_norm`). Not fixed here: it is the loader's
PL rule (RM2). The repair is to divide in Float64 and cast the result to Float32 once.

Resolved 2026-10-22 as proposed: PL is divided by -10 in Float64 and rounded to Float32 once, and so is
Complete Genomics' phred-unit GL (divided by 10). PL `p` now loads as the Float32 nearest `-p/10`, the
value a GL file's `-p/10` is parsed to; checked for every PL from 0 to 299. Documenting a tolerance
instead was set aside: the repair costs one cast. This changes `gl` values on files that already load,
by one Float32 step on some entries. The comparison with bcftools now holds exactly. Code:
`loaders/vcf.py` `to_rows` (`gl`); tests: `tests/test_load_vcf.py`
`test_pl_and_the_same_likelihood_as_gl_load_equal`, `tests/test_normalize.py`
`test_agrees_with_bcftools_norm`.


## F14 — Ensembl's GRCh38_to_GRCh37 chain lists every block twice

Severity: high · Found 2026-10-22 lifting anton's GRCh38 callset with the chain CrossMap had used (RM6).

The Ensembl chain (`GRCh38_to_GRCh37.chain.gz`, from ftp.ensembl.org `assembly_mapping`) holds 22,364
chains, and every block of the main ones appears twice, under two chain ids, with the same source span,
target and strand. The first `liftover` counted blocks per base, so 5,717,715 of 6,079,744 rows read as
lying in more than one block (`multiple_targets`) and none was lifted. UCSC's hg38ToHg19 chain has no
overlapping block at all, so the synthetic tests, built like it, could not show it.

Resolved the same day in `liftover._segments`: the source axis is cut where the set of *targets*
changes, a target being a contig, strand and offset, so two blocks sending every shared base to the same
place count once. A position under blocks with different targets is still `multiple_targets`. With the
fix the same chain lifts 5,578,476 rows, and 108,776 stay ambiguous. Test:
`test_a_chain_listing_its_blocks_twice_lifts_as_once` (fails on the first version, checked).

## F13 — An FTDNA export reports some probes on the minus strand

Severity: medium · Found 2026-10-22 by F9's row inference on the 103 PGP-HMS FTDNA exports (counts
only).

One autosomal FTDNA export (about 709,000 rows) fits GRCh36 at 68.2% of sampled SNV rows, against
80–83% for the other 59 autosomal FTDNA exports. Of its heterozygous rows, only 81% carry the reference
base (99.96–100% in every other export with autosomes), and 13.0% of all sampled rows carry the complement of
the reference base (0.0–3.0% elsewhere). So a share of its probes, about one
in eight, is written on the minus strand. `resolve_ref` reads REF for 99.8% of its rows, and on the
flipped rows derives `alts` and `gt` against the wrong strand: a flipped hom-ref row reads as hom-alt.
The existing warning counts them among the rows that "carry no REF allele", which does not single
them out.

Not fixed. Candidate repairs, each wrong for now: complementing a row whose alleles miss REF (an A/T
or C/G probe is its own complement, and a true hom-alt call would be flipped into a hom-ref one);
scoring both strands (a complemented het covers all four bases, so a wrong build would fit as well as
the right one); a per-chip strand table (the right fix, but FTDNA ships none in the export, and the
chip is not identified yet: RM4's chip generation detection).

Resolved 2026-10-22 by withholding, not by complementing. `resolve_ref` now checks the strand of every
array facet it reads, on F9's row sample against the collection it reads, through the heterozygous
rows: on a two-allele site a het carrying the reference base is safe on either strand, and a het
carrying the base's complement and not the base is impossible on the forward strand. More than 2% of
at least 200 such hets makes the facet `strand_check = mixed`; then every row that is not strand-safe
(homozygous calls, half-calls, hets without the base) keeps `ref` and gets `alts` and `gt` null, each
reason counted in the warnings, and the safe hets keep their derived alleles. Rule: docs/SCHEMA.md
Decision 2 ("The strand of an array's genotypes"). Code: `strand.py` `check_strand`,
`withheld_reason`; `reference.py` `_check_strand`, `_derive_array_alleles`. Tests:
`tests/test_reference.py` from `test_a_mixed_strand_array_keeps_ref_and_withholds_the_alleles_it_cannot_read`
(each fails on the code before this change).

On F13's export: `mixed`, 2,351 of about 12,600 sampled hets (18.65%) carry the complement and not the base,
none carry neither. On the other FTDNA exports of 600,000 rows or more that ran before the scan was
stopped (10 of them, on OmniExpress like this one): 0.00–0.02%, all `forward`. On 32 more, two of
each format and detected chip (23andMe v1/v2, v3, v4, v5 and unknown, its CSV re-save, AncestryDNA
V1.0 and V2.0, deCODEme, FTDNA space-separated, Genes for Good measured and imputed, the GSA Final
Report, Living DNA, SelfDecode, the rsid CSV and the rsid tables), 1,171–13,794 hets sampled each:
every one `forward`, 0 or 1 hets carrying the complement and none carrying neither (at most 0.01%).
The 2% threshold sits between 0.02% and 18.65%.

Why each other repair is wrong:

- **A per-chip strand table.** The chip of this export is not detected (FTDNA states none and RM4 has
  no rule for it), but its probe set is OmniExpress's: 98–99% of its sampled rsids are in Illumina's
  Infinium OmniExpress-24 v1.2 manifest, as for every FTDNA export of that size. That manifest's three
  strand columns do not separate its minus-strand probes. Over all hets of the export (aggregates,
  matched by rsid): proven-forward or strand-invariant hets split RefStrand `+` 83,306 / `-` 83,552,
  proven-minus hets 19,337 / 19,467; IlmnStrand TOP 85,247 / BOT 81,611 against 19,821 / 18,983;
  SourceStrand TOP 83,213 / BOT 83,645 against 19,205 / 19,599. The flipped probes are on every
  autosome (16–24% of hets each) and not in runs (mean run 1.54 hets in file order). So no public
  table says which probes FTDNA wrote flipped in this export, and complementing by any of these columns
  would flip about half of the right rows and half of the wrong ones.
- **Complementing the hets that prove a minus probe.** It is a deduction for those rows, but the
  homozygous rows of the same probes cannot be found, so the facet would mix corrected and uncorrected
  rows with nothing in the row to tell them apart. Withholding is uniform: every non-safe row is null.
- **Complementing every row whose alleles miss REF, or scoring both strands**: as above (an A/T or
  C/G probe is its own complement; a complemented het fits any build).
- **Withholding nothing and warning only**: a flipped hom-ref row would still read as a confident
  hom-alt call downstream, which is the error F13 is about.


## F10 — A GP that does not sum to 1 is lost when just-vcf reloads its own export

Severity: low · Found 2026-10-22 by round-tripping a PGP-HMS VCF upload through `write_vcf` and `load`
(RM6).

The upload writes GP phred-scaled (values above 1), so the loader converts it and every row gets `gp`.
Phred values rounded to integers do not convert to posteriors that sum to 1: 57,259 of its 956,057 GP
rows are more than 0.02 short. `write_vcf` writes `gp` as probabilities under a VCF 4.3 header, where
the spec defines GP as a probability. Reloading that output, the loader reads the scale from the values
alone: all within [0, 1], and not all summing to 1, so it calls GP ambiguous and drops it from every
row. Everything else in the file came back byte for byte.

Not the writer's to fix: normalizing the posteriors would change the evidence, and writing them back as
phred under a 4.2 header would not read back bit for bit.

Resolved 2026-10-22 in the loader's GP rule: values the numbers cannot settle (all within [0, 1], not
all summing to 1) are read as probabilities under a VCF 4.3+ header, which defines GP so, and the rows
that do not sum to 1 are counted in a warning. A phred GP whose every value is at most 1 would mean
genotypes all but equally likely, which is not what any caller writes. Under 4.2 or earlier, or a header
with no version, they stay withheld as before; values above 1 are phred and values summing to 1 are
probabilities, whatever the header. The FORMAT line's Description was considered as a second header
signal and set aside: it is free text, not a declaration. On the upload that found it: the reload keeps
`gp` on all 956,057 rows, equal to the first load's (largest difference 0.0), with one warning counting
the 57,259 rows that do not sum to 1. Code: `loaders/vcf.py` `_gp`; tests:
`tests/test_load_vcf.py` `test_gp_from_phred_survives_its_own_export` and the 4.3, 4.5 and no-version
cases of `test_gp_scale_is_read_from_decisive_values`.
## F9 — An array export that states no build cannot get its REF

Severity: high · Found 2026-10-22 by running `resolve_ref` (RM3) on the table formats RM2 just made
loadable, on public PGP-HMS uploads (counts only).

`resolve_ref` worked where the facet carries an assembly: an AncestryDNA export (stated GRCh37) had REF
read for 681,662 of 701,478 rows (97.2%), a Color file for all of its rows. Two FTDNA exports were
refused with `ReferenceMismatchError` ("the facet states neither contig lengths nor an assembly"): an
FTDNA file names no build, and an array has no `##contig` lines, so RM4's length inference has nothing
to read. A deCODEme export (GRCh36, stated in its zip) was refused because the store holds only GRCh37
and GRCh38. The refusals are correct (never assume a build), so the gap is evidence, not the check.

Not fixed. Candidate repairs: (a) infer an array's build from its own rows against each collection in
the store, as the share of SNV rows whose called alleles include the reference base at that position
(the right build scores near 100%, a wrong one near chance), recorded as `inferred` with the scores; (b)
add NCBI36 to the reference store for the GRCh36 exports. A caller-asserted `assembly=` argument was
considered and set aside: it lets a label stand in for evidence.


Resolved 2026-10-22, both repairs. (a) `resolve_ref` infers the build of a facet that states neither
a label nor contig lengths from its own rows, against every collection it considers, and records it as
`assembly_label` with basis `inferred` and the scores in one warning line; rows that do not settle it
refuse the facet with the scores (rule: docs/SCHEMA.md Decision 2; choices and rejected alternatives:
ROADMAP RM3). (b) A collection whose contig 1 is 247,249,719 bases is labelled `GRCh36`; UCSC hg18 is
the FASTA to add (where and why not Ensembl 54: RM3). Code: `assembly.py` `sample_snv_rows`,
`infer_assembly_from_rows`, `MIN_SHARE`/`MIN_MARGIN`; `reference.py` `_infer_from_rows`, `_bases_at`,
`CHR1_LENGTHS`. Tests: `tests/test_reference.py` `test_an_array_without_a_build_gets_it_from_its_rows`
and the refusals after it, `test_a_grch36_reference_is_named_by_its_contig_1_and_cannot_vouch_for_mt`;
`tests/test_assembly.py` from `test_the_sample_keeps_scorable_array_snvs_only`.

Dogfooded on 153 loadable PGP-HMS array exports from the Krypton share, against a store holding Ensembl
GRCh37 and GRCh38 and UCSC hg18 (aggregates only). Share of sampled SNV rows carrying each collection's
reference base, the right build against the wrong ones, and the share of rows whose REF `resolve_ref`
then read:

| Exports | Files | Inferred | Right build | Wrong builds | REF read |
|---|---|---|---|---|---|
| FTDNA, autosomal + X | 60 | 36 GRCh36, 24 GRCh37 | 0.682–0.825 | 0.311–0.327 | 95.4–99.9% |
| FTDNA, X only (about 2,000 rows) | 43 | 30 GRCh36, 13 GRCh37 | 0.665–0.843 | 0.230–0.334 | 94.4–99.9% |
| FTDNA, space-separated | 2 | GRCh36 | 0.825–0.827 | 0.313–0.318 | 99.4–99.6% |
| deCODEme (build stated in the zip) | 4 | GRCh36, as stated | 0.811–0.814 | 0.307–0.316 | 96.0–99.9% |
| 23andMe, stated build 36 | 13 | GRCh36, as stated | 0.793–0.810 | 0.310–0.329 | 98.3–99.4% |
| 23andMe, stated build 37 | 13 | GRCh37, as stated | 0.801–0.906 | 0.283–0.324 | 95.5–99.1% |
| AncestryDNA, Genes for Good, Living DNA, 23andMe CSV (stated) | 10 | as stated | 0.806–0.892 | 0.283–0.327 | 96.2–99.99% |
| rsid tables, PLINK 23, rsid CSV, Final Report (no build) | 8 | 2 GRCh36, 6 GRCh37 | 0.796–0.913 | 0.263–0.326 | 96.6–99.4% |

For the stated-build exports the inference ran as if no build were stated: it named the stated build
in all 40. hu832966's two 23andMe exports infer GRCh36 and GRCh37, each its stated build; REF is read
for 99.2% and 97.9% of their rows. Before this, every FTDNA export, every unlabelled table and every
build-36 export was refused. One Final Report did not load (`TableLayoutError`, not this finding's).

Strand: forward only. The share of rows carrying the complement of the right build's base was
0.0001–0.030 in 152 files, so complementing would gain nothing there and would let a wrong build fit
(a complemented het covers all four bases). The 153rd, an FTDNA export at 0.13, has part of its probes
on the minus strand: F13. The floor was first 0.70; it refused 27 FTDNA X-only exports and F13's file,
whose right build stood 0.36–0.46 above every other, so it is 0.60 and the 20-point margin carries the
separation. Cost: with the store warm, sampling and reading the bases of three collections takes about
0.3 s per facet.


## F8 — A sites-only VCF does not load

Severity: medium · Found 2026-10-19 while testing `rejoin_split` (RM270).

A VCF with no FORMAT column and no sample (a sites-only file, as annotation sources and some
callers' site lists are) fails in `load` with polars' `InvalidOperationError: to_arrow() conversion
failed: cannot convert (Null) to (... LargeUtf8 ...)`, raised by the parquet spill (`_load_vcf`'s
`sink_parquet`). It failed the same way before RM270 (checked at c746325), so it is older than the
split flags. polars-bio itself reads the file: `pb.scan_vcf` returns the two records with their eight
site columns. So the fault is on our side, most likely a column the loader adds as an untyped null
for the absent sample. Repro: a two-record VCF 4.2 with one `##contig` line and an eight-column
header. A sites-only file has no genotype, so it is a `calls` frame of nulls at best.

Resolved 2026-10-22: `load` refuses it with `NoSampleError`, before anything is read or spilled. just-vcf
holds a person's calls (the charter's first line), and a file with no sample column has none to hold;
a `calls` frame of null genotypes would look like a person with nothing called. `vcf_samples` still
reads the file and returns `[]`, and `sample=` on it is refused as a missing sample, as before. Code:
`loading.py` `_load_vcf`; test: `tests/test_load_vcf.py` `test_sites_only_vcf_is_refused`.

## F6 — A bz2 or zip VCF was copied out as plain text first

Found 2026-09-30 by loading the 124 PGP-HMS VCFs. polars-bio reads plain, gzip and BGZF VCF directly
(79 of the 124); a VCF inside bz2, xz, zip or tar is decompressed into the loader cache first. That copy
was plain text: a 1.27 GB bz2 Veritas gVCF became a 14 GB file.

Resolved 2026-09-30: the copy is written as BGZF by `just_vcf.bgzf.BgzfWriter` (stdlib `zlib` on a thread
pool, no new dependency). On a 2 GB slice of a real VCF: BGZF at level 1 is 0.50 GB and takes 17 s to
write; polars-bio counts its 18.5 million records in 9.5 s, against 9.1 s for plain text and 24.8 s for a
plain gzip copy. Code: `loaders/vcf.py` `readable_path`.

## F5 — The canonical sort cost most of a VCF load's memory

Found 2026-09-30 loading a DeepVariant VCF of 6.4 million calls: writing `calls` sorted peaked at
11.2 GB, against 3.5 GB unsorted. The user picked "keep file order within a contig, sort only contig
blocks and ties". Measuring it showed that reordering contig blocks is not cheap either: polars-bio does
not use the `.tbi` index, so a per-contig filter is a full scan (4.3 s each, against 6.2 s for the file).

Resolved 2026-09-30 without changing SCHEMA.md Decision 11: `load` spills the VCF frames unsorted to
parquet in the cache, in one streaming pass, and returns scans that read them back one contig at a time
in canonical order (`schema.canonical_order_per_contig`), using parquet row-group statistics to skip to
each contig. Same file: 77 s and 3.6 GB peak, row order identical to the global sort. The spill was
0.14 GB. Code: `loading.py` `_load_vcf`; test: `tests/test_schema.py`.
