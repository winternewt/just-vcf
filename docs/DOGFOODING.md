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

## F25 — UCSC scaffold names written without `chr` are not matched to the store

Severity: low · Found 2026-10-26 dogfooding RM6's SV lift on hu416394's Manta file.

Sequencing.com's GRCh38 header names the scaffolds `1_KI270706v1_random` and `Un_KI270302v1`: UCSC's
names with `chr` cut off. F19's alias table carries the names NCBI's assembly reports give
(`chr1_KI270706v1_random`, `KI270706.1`, `NT_187361.1`), and this spelling is none of them, so 275
rows (96 on `_random` scaffolds, 179 on unplaced `Un_` ones) stay `ref_unread` and cannot lift, though
the store holds the sequences. Repairs: add the `chr`-less spelling of each UCSC name to the alias
table, which is the name-pattern rule F19 rejected (its reasons are in PREVIOUS_ISSUES.md);
or match a header contig by its stated length when exactly one scaffold of the
assembly has it, which is a new matching rule. Low: 3% of this file, and only its breakends and SVs
touching scaffolds.

*Updated 2026-10-26 (breakends in a lift):* the same names hold breakends back on both ends: 163
breakend records on these scaffolds stay `ref_unread`, and 20 records elsewhere whose mate lies on one
stay `breakend_mate_not_lifted` (the mate's contig is read with `canonical_chrom`, which does not know
this spelling either, so the chain's `chr`-named scaffold is not found; whether the chain would place
those mates once the name matched was not checked).

## F31 — A lift can name one target contig two ways in one facet

Severity: low · Found 2026-10-26 writing hu1C1368's PGP-HMS Manta file lifted to GRCh38 (F30's repair).

The lifted facet has rows on `chrUn_KI270742v1`, the name UCSC's chain gives the scaffold, while a
breakend lifted onto that scaffold names its mate `KI270742.1`, the target collection's `Contig.name`
(RM6's "Lifting breakends" names mates that way), which `canonical_chrom` keeps as it is. So one
contig carries two canonical names in one facet. `write_vcf` copes since F30: both names are declared
under `stored` and `chr`, merged into one `##contig` line under `source`. What does not cope:
`concordance` matching a mate against another facet's `chrom` misses a mate on that contig, and if rows
ever stood on both names, a `source` body would show one contig in two separate runs, which tabix
refuses. It predates F30; the lift's rows take the chain's name and its mates the collection's.
Repairs: name lifted rows by the collection too, which renames every lifted scaffold row; or map the
chain's name to the collection's through F19's alias table before either is written, which is the
table's job. Not fixed: which name a lifted facet should carry is a naming rule for the whole lift.
