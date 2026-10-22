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

**Measured and surfaced 2026-10-22.** Both halves are one person: of the 164,378 ids both halves call,
164,275 have the same genotype (99.94%), and 166,047 of 166,528 shared ids sit at the same position.
The first line is a stray 5-field row (an empty field between position and genotype, not F2's place),
then 643,486 four-field rows, then 677,436 five-field rows to the end (chromosomes 1–26, `0` for no
call). The load now says so: a file that ends in at least 1,000 rows of one other width gets a warning
that it looks like a second export appended, not read. Still not loaded, and why each repair is wrong:

- *Read the 5-field rows as a second facet in this loader:* the library cannot load that half even as
  its own file (`ancestrydna` needs its header line), so it needs a new headerless two-allele format id
  and a part option in the facet id, for one file whose halves state no vendor. Calling the second half
  AncestryDNA would be the shape standing in for a vendor statement.
- *Split the file into two inputs and load each:* that is routing around the loader (CLAUDE.md §7).

What would close it: a headerless two-allele-column layout in `formats.rsid_tables`, and `load`
returning one facet per layout block of a file (the shape `Loaded` already has for caller splits).
Code of the warning: `loaders/rsid_table.py` `_appended_export`; test
`test_an_appended_export_in_another_layout_is_named_not_read`.

## F4 — Living DNA indel genotypes cannot be split

Severity: medium · Found 2026-09-30, same sweep.

Living DNA writes an indel's two allele sequences back to back (`CCA`, and up to 52 letters in one
file). 83 rows in the one real file keep the text in `source_gt` with `gt_alleles` null. Splitting needs
the site's alleles, from a chip manifest or from REF once RM3 reads the sequence.
