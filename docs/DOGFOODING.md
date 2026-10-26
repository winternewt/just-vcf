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

## F23 — A CNV caller's `N` REF reads as a mismatch, so no CNV facet lifts

Severity: medium · Found 2026-10-26 dogfooding RM6's SV lift on the PGP-HMS CNV files.

DRAGEN's and Canvas's copy-number VCFs write REF `N` on every record: the caller never reads the
reference. `resolve_ref` compares that `N` with the base it reads and sets `ref_check = mismatch`
on every row, and `liftover` leaves every `mismatch` row behind as `source_ref_mismatch`. Measured
with the library: hu6ABACE's DRAGEN 3.2.8 `.cnv.` (586 rows, GRCh37), hu1C1368's DRAGEN 3.4.5
`.cnv.` (556 rows, GRCh37) and hu416394's Canvas `.cnv.` (182 rows, 92 copy-number calls and 90
copy-neutral segments, GRCh38): 0 match, 100% mismatch, 0 rows lifted through `hg19ToHg38` or
`hg38ToHg19`. So the copy-number half of RM6's SV lift has never run on a real file; the 31 `N` rows
of hu416394's Manta file sit on decoys and are unread anyway.

The IUPAC note in the resolve warning does not cover it: `_iupac_explains` counts an ambiguity code on
the reference side (a gap base) against a concrete source base, not a source `N` against a concrete
reference base. The SV lift itself needs nothing new: on a plus strand `_place_sv` compares the
target's base with `ref`, the base read from the source assembly, not with `source_ref`. Only the
`ref_check` gate stops these rows. Each candidate repair needs a decision:

- *Read a source `N` as matching any base (`match`):* the row moves from honestly unverified to
  apparently verified, a check that never ran reported as passed (§5, Kleene).
- *A new `RefCheck` value for a source REF that states no base (`unstated` or similar), lifted like
  `match`:* the honest reading, but `RefCheck` is a `pl.Enum` in the wire schema, so a new value is a
  schema change, and `mint_vrs_ids`, `write_vcf(resolved_only=True)` and `normalize` each have to say
  what they do with it.
- *`not_run` for an all-`N` source REF, with `ref` read:* no schema change, but `not_run` today means
  the comparison did not happen for want of a reference, and the two causes would share one value.
- *Special-case symbolic rows in `liftover`:* makes SVs the one row kind whose unstated REF passes,
  and leaves `resolve_ref`'s counts calling 1,324 public rows mismatched.

Blocks: lifting any CNV facet (every one in hand); also any QC that reads `ref_check` on them.

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

## F26 — Loading a WGS gVCF holds gigabytes

Severity: low · Found 2026-10-26 dogfooding RM6's block lift on the PGP-HMS gVCFs.

`load` of hu416394's bcftools gVCF (a 397 MB member of a zip, 9,651,143 calls and 9,473,303 blocks)
peaks at 7.8 GB on a fresh cache (151 s) and at 7.7 GB again with its frames already cached (89 s);
`resolve_ref`, `liftover` and `write_vcf` after it stay within that peak (sampled per phase). The
whole pipeline over hu80855C's Isaac gVCF (1.9 GB, 45,242,155 calls and 138,019,217 blocks) peaked at
46 GB by the end of `load` (1,008 s) and 49 GB on the reload, on a 62 GB machine; that run's own
counting over the frames is inside the figure, so `load`'s share there is not measured. Not looked
into: where in `load` the memory goes (the per-contig spill is meant to bound it, PREVIOUS_ISSUES has
the 6.4-million-call DeepVariant file at 3.6 GB), and why a load whose frames are cached still takes
89 s at the same peak. Blocks nothing yet; a 30x Isaac gVCF is near the edge of this machine.

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
