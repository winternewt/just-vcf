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

**Veritas shape fixed 2026-10-22.** Not every header line has the extra `#`: `##fileformat` and `##ALT`
are written right, `###FORMAT`, `###contig` and `###reference` carry a third, and the column line is
`##CHROM`. Where no `#CHROM` line exists, a meta line that is VCF's column line with one `#` too many
is now the column line, and the file is read through a rewritten copy; the `###` lines are kept as
the meta lines VCF says they are (a key starting with `#`), so no contig and no FORMAT is declared from
them, and GT and END are declared by the loader as for any undeclared key. Both files load (5,614 SV
calls, all `0/1`; 2,632 CNV rows, all `0/1`, whose ALT `CNV` is F15). Code: `loaders/vcf.py`
`read_header`, `HASHED_COLUMN_LINE`; test `test_a_column_line_written_with_two_hashes_is_read`.

A third failure in the sweep is a corrupt upload (`gzip -t` fails too), and its error is correct.

