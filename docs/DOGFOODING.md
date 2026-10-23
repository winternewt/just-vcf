# Dogfooding findings — open

Quirks, bugs and gaps found by using just-vcf on real genomes. Each carries a stable `F#` id and
moves to [PREVIOUS_ISSUES.md](PREVIOUS_ISSUES.md) when resolved here, or to
`<UPSTREAM>_PENDING_FIXES.md` when it needs an upstream change first.

## F20 — A CGI vcfBeta's chrY md5 refuses the whole facet

Severity: medium · Found 2026-10-23 by RM279 (the method-error atlas) comparing PGP-HMS CGI vcfBeta
facets against freeBayes.

`resolve_ref` refuses a facet when any primary contig (1–22, X, Y) states a length or md5 the
reference does not have. Checked with the library on the smallest PGP-HMS vcfBeta against the shared
store: the header states 25 contigs (1–22, X, Y, M), every autosome and X agree with Ensembl GRCh37
in length and md5, and Y agrees in length (59,373,566) but not in md5. So the whole facet is refused:
no `ref`, no `ref_check`, no VRS ids, and the atlas compares the CGI pairs by label only
(docs/METHOD_ATLAS.md). Assemblies of one build legitimately differ on X, Y and MT (PAR masking in an
analysis set, rCRS against the Yoruba mitochondrion); an autosome of another sequence is another
assembly. The lead's call (2026-10-23): when no autosome mismatches, read every matching contig and
leave X, Y or MT rows that mismatch unread, counted per contig; an autosome mismatch still refuses.

## F1 — Deflate64 zips cannot be opened

Severity: low · Found 2026-09-29 by `scripts/sniff_corpus.py` over the PGP-HMS participants.

One PGP upload (a Veritas VCF of 2.9 GB, zipped by Windows Explorer) uses zip method 9, Deflate64.
The standard library cannot read it, so `sniff` reports `unreadable` and names the method. Reading it
needs a third-party inflater (the `zipfile-deflate64` package, or `7z` as a subprocess) behind an extra.
Not fixed: one file, and the same person's VCF may exist unzipped elsewhere.

**Corrected 2026-10-22.** The member is not a Veritas VCF: read through Info-ZIP `unzip -p` (built
with Deflate64), its header is PGP's own freeBayes+VarScan2 `vcf-merge` on hg19, 2.93 GB unzipped
from 678 MB, a format the loader already reads (`vcf_freebayes_varscan2`). The person has nothing
else on disk but a 23andMe export, so this is their only sequencing callset here. `sniff` already
answers right (`unreadable`, naming the method). Reading it needs a Deflate64 inflater, and each
option is a dependency decision, which is the lead's (§5):

- `zipfile-deflate64` or `inflate64` (a C extension, the one py7zr uses) behind a new extra, opened
  where `sniff` and `open_payload` hit method 9; a new dependency, licence to check;
- `unzip` or `7z` as a subprocess when one is on `PATH`, streaming the member: no Python dependency,
  but a system one, and a subprocess in a library;
- leave it: one file in 1,851.

## F17 — Living DNA indel calls need the chip's site alleles for REF

Severity: low · Found 2026-10-22 resolving F4 on the one real Living DNA file.

`resolve_ref` now splits Living DNA's run-together indel genotypes against the reference (F4, in
[PREVIOUS_ISSUES.md](PREVIOUS_ISSUES.md)). Of the 83 rows, 4 get REF and ALT. 54 are split with REF
unread: 38 homozygous calls, where the text shows one of the site's two alleles and so cannot say
whether the site is `CA>C` (the call is hom-ref) or `C>CA` (hom-alt), and 16 hets whose two alleles
both equal the reference over their own spans (an insertion and a deletion fit alike; every het
deletion is one). 25 stay unsplit: 24 fit the reference at several splits (repeats), 1 at none. These
rows cannot be written to a VCF or minted. What would settle them is the site's own two alleles, from
the Sirius chip's manifest (LANDSCAPE: chip manifests), which the library does not read. Why each
repair from the text alone is wrong:

- *Read a homozygous call whose allele equals the reference as hom-ref:* it picks one of two sites the
  text fits, which is the genotype the row is about.
- *Take the shortest allele as REF on a het:* it calls every ambiguous het an insertion.
- *Split a repeat at the anchored point nearest the middle:* it chooses among splits the reference
  fits equally.

Checked 2026-10-23 for RM281 (chip manifests), on the same file, counts only: no public manifest
describes the Sirius chip. It is not a GSA: 14,360 of its 691,603 rows carry Affymetrix Axiom ids
(`AX-…`, and one `AFFX-…`), and 405,273 rows sit at a GSA-24v3-0 A1 probe position. Of the 83 indel
rows, 2 carry an rsid that names a GSA indel probe; the other 81 match none by name, and 2 by position
one base off. So the GSA manifest settles at most 2 rows, and those only by assuming Living DNA's
probe tests GSA's allele. A learned probe set cannot settle it either: PGP-HMS holds one Living DNA
export, and learning its site alleles from the same person's WGS would make the array agree with that
WGS by construction (a value filled from the source that checks it). Stays open until a second
source of the Sirius chip's site alleles appears.

