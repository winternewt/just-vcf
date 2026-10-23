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


## F21 — A VRS-met allele is read as zero copies in a complex multi-allelic row

Severity: medium · Found 2026-10-23 by RM279 rebuilding the PGP-HMS CGI against freeBayes pairs with
alleles met by VRS id.

freeBayes writes a haplotype as one multi-allelic row: REF `CGCCC`, ALT `TGCCT,TGCCC`, GT
`TGCCT/TGCCT`. Its ALT `TGCCC` trims to `C>T`, so its VRS id is CGI's `C>T` (GT `T/T`) and the two
meet (`key = vrs`). `concordance` then counts the copies of the spelling `TGCCC` in the freeBayes
genotype, zero, although each `TGCCT` haplotype carries the same `C>T`: a `mismatch` that is only
representation. 20,434 to 23,014 SNV mismatches on the CGA 2.5.0.44 pair, 2,601 to 2,950 on the
2.5.0.24 one, up to 5,982 and 4,627 indels (docs/METHOD_ATLAS.md). The `site` path guards the same
case by reading only a row over the same REF span; the `vrs` path has no guard. The graded
comparison has the same fault: its distribution counts the copies of the spelling too.

## F22 — A haploid call against the same allele's diploid homozygous call counts as discordant

Severity: medium · Found 2026-10-23 by RM279 on the PGP-HMS CGI against freeBayes pairs.

Outside the PARs CGI writes a male's X call haploid (`1`) and freeBayes diploid (`1/1`). Every such
pair is `ploidy_differs`, which the atlas counts in `discordant`: 72,462 and 64,818 SNVs, 65% and
67% of the two pairs' discordant SNVs. Both say the person carries the allele on every copy of X
they have; the difference is the writer's ploidy convention, not a method error, and it swamps the
rate the atlas exists to give.
