# Previous issues — resolved here

Findings from [DOGFOODING.md](DOGFOODING.md) resolved in this repo, each with its resolution and a
code pointer. Check here before re-investigating something that looks fixed.


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
