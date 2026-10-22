# Dogfooding findings — open

Quirks, bugs and gaps found by using just-vcf on real genomes. Each carries a stable `F#` id and
moves to [PREVIOUS_ISSUES.md](PREVIOUS_ISSUES.md) when resolved here, or to
`<UPSTREAM>_PENDING_FIXES.md` when it needs an upstream change first.

## F1 — Deflate64 zips cannot be opened

Severity: low · Found 2026-09-29 by `scripts/sniff_corpus.py` over the PGP-HMS participants.

One PGP upload (a Veritas VCF of 2.9 GB, zipped by Windows Explorer) uses zip method 9, Deflate64.
The standard library cannot read it, so `sniff` reports `unreadable` and names the method. Reading it
needs a third-party inflater (the `zipfile-deflate64` package, or `7z` as a subprocess) behind an extra.
Not fixed: one file, and the same person's VCF may exist unzipped elsewhere.

## F2 — 23andMe rows with a doubled tab are dropped

Severity: low · Found 2026-09-30 by loading the 503 PGP-HMS 23andMe-layout files.

Two 2012–2014 23andMe exports carry rows with an empty field between the rsid and the chromosome
(`rs15842⇥⇥1⇥948921⇥CT`): 14 and 35 rows out of about 960,000. The loader drops them and says so in
`Loaded.warnings`. Not fixed: accepting "exactly one empty field after the id" is a guess about what the
vendor meant, and the rows are few. Revisit if another vendor shows the same shape.

## F3 — A headerless file holding two exports loses the second one

Severity: low · Found 2026-09-30, same sweep.

One upload (`rsid_table_headerless`) is a 23andMe-shaped export followed by an AncestryDNA-shaped one
(5 columns). The loader keeps the 643,487 four-column rows and drops the 677,437 five-column rows with
a warning, as FORMATS.md predicted. Not fixed: splitting one file into two facets needs a second format
per file, which the loader does not have.

## F4 — Living DNA indel genotypes cannot be split

Severity: medium · Found 2026-09-30, same sweep.

Living DNA writes an indel's two allele sequences back to back (`CCA`, and up to 52 letters in one
file). 83 rows in the one real file keep the text in `source_gt` with `gt_alleles` null. Splitting needs
the site's alleles, from a chip manifest or from REF once RM3 reads the sequence.

## F7 — Two VCF shapes `sniff` accepts but the loader cannot read

Severity: low · Found 2026-09-30, same sweep.

A VCF re-saved by a spreadsheet (tab-padded, quoted `"##` lines) and Veritas SV/CNV VCFs whose every
header line carries an extra `#` (`###contig`, `##CHROM`) are recognised as `vcf` but fail with "no
#CHROM line". The Veritas shape could be read by accepting `##CHROM` as the column line; the spreadsheet
one needs its own format id, as `23andme_csv` has. Not fixed: three files.

A third failure in the sweep is a corrupt upload (`gzip -t` fails too), and its error is correct.

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

## F10 — A GP that does not sum to 1 is lost when just-vcf reloads its own export

Severity: low · Found 2026-10-22 by round-tripping a PGP-HMS VCF upload through `write_vcf` and `load`
(RM6).

The upload writes GP phred-scaled (values above 1), so the loader converts it and every row gets `gp`.
Phred values rounded to integers do not convert to posteriors that sum to 1: 57,259 of its 956,057 GP
rows are more than 0.02 short. `write_vcf` writes `gp` as probabilities under a VCF 4.3 header, where
the spec defines GP as a probability. Reloading that output, the loader reads the scale from the values
alone: all within [0, 1], and not all summing to 1, so it calls GP ambiguous and drops it from every
row. Everything else in the file came back byte for byte.

Not fixed, and not the writer's to fix. Normalizing the posteriors would change the evidence, and
writing them back as phred under a 4.2 header would not read back bit for bit. The repair belongs in
the loader's GP rule (RM2): values within [0, 1] under a 4.3+ header are probabilities, because the
header settles what the values cannot, and a phred GP whose every value is at most 1 would mean
genotypes all but equally likely, which is not what any caller writes.

## F11 — PL to `gl` lands one float step off the nearest value

Severity: low · Found 2026-10-22 by comparing `normalize` with `bcftools norm` (RM6).

The loader turns PL into `gl` as `PL / -10` in Float32, and Polars computes that division so that PL 21
gives -2.1000001, while the nearest Float32 to -2.1 is -2.0999999. A file that writes the same
likelihood as GL (`-2.1`) therefore loads one step away from one that writes PL 21. Deterministic, so no
digest or test of ours moves, and the size is far below any PHRED rounding; but two facets of one
genome, one with PL and one with GL, do not compare equal on `gl`, and a test comparing against an
outside tool needs a tolerance (`test_agrees_with_bcftools_norm`). Not fixed here: it is the loader's
PL rule (RM2). The repair is to divide in Float64 and cast the result to Float32 once.
