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

## F11 — PL to `gl` lands one float step off the nearest value

Severity: low · Found 2026-10-22 by comparing `normalize` with `bcftools norm` (RM6).

The loader turns PL into `gl` as `PL / -10` in Float32, and Polars computes that division so that PL 21
gives -2.1000001, while the nearest Float32 to -2.1 is -2.0999999. A file that writes the same
likelihood as GL (`-2.1`) therefore loads one step away from one that writes PL 21. Deterministic, so no
digest or test of ours moves, and the size is far below any PHRED rounding; but two facets of one
genome, one with PL and one with GL, do not compare equal on `gl`, and a test comparing against an
outside tool needs a tolerance (`test_agrees_with_bcftools_norm`). Not fixed here: it is the loader's
PL rule (RM2). The repair is to divide in Float64 and cast the result to Float32 once.

## F12 — Some VCF uploads write a deletion as ALT `-`

Severity: low · Found 2026-10-22 running `normalize` against `bcftools norm` on PGP-HMS uploads (RM6).

Five of 25 PGP-HMS VCF uploads of 5–60 MB hold about 8,080 rows each whose ALT is `-`, the dash
notation for a deletion that VCF does not have. The loader keeps it verbatim as an allele, `normalize`
leaves those rows as they are and counts them, and `write_vcf` writes them back out, so the export is a
file htslib refuses ("Non-ACGTN alternate allele"). Nothing in the frame says the allele is not one.

Not fixed here. Reading `-` as a deletion needs the anchor base and its REF span, which the row does
not state, so it is not a spelling to rewrite; the options are for the loader (RM2: flag or refuse
such rows, counted) and the writer (refuse them as unwritable, as it does a row without REF).

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
