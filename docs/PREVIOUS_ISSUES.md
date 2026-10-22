# Previous issues — resolved here

Findings from [DOGFOODING.md](DOGFOODING.md) resolved in this repo, each with its resolution and a
code pointer. Check here before re-investigating something that looks fixed.

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
