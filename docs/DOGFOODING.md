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

## F33 — `write_vcf(resolved_only=True)` keeps a phase set whose member it skipped

Severity: low · Found 2026-10-31 by the builder making phase blocks follow F23's rule (merge c94251f).

`resolved_only=True` skips each unwritable row and writes the rest, phase sets included. When a skipped
row sat inside a phase set (`PS`, or a VCF 4.5 `PSL`), the set's other rows are written with the same
set and phase, so on reload `cis_phased_blocks` builds that haplotype's block from fewer alleles than the
source had, and a block the source facet withholds (for the skipped row's allele) comes back stated. It
predates F23 and holds for every unwritable reason, `ref_mismatch` and `ref_unstated` alike. Repairs,
none made: write the set's remaining rows unphased, which drops phase the source did state for them; cut
the set at the skipped row into two sets, which states a break the source did not; or refuse a set with
a skipped member and count it, which loses rows `resolved_only` promises to write. Which of the three is
a judgement on what a written set may claim.

## F34 — `rejoin_split` on a facet whose REF was read loses the joined rows' REF for good

Severity: low · Found 2026-11-03 by the builder fixing merge audit 14's finding 4 (`rejoin_split` and
the `unstated` rule).

`rejoin_split` builds each joined row from its split rows' `source_ref` alone: `ref` null, `ref_check =
not_run`, no VRS ids. On a facet straight from `load` that is what every row has, and `resolve_ref` reads
it next. On a facet already through `resolve_ref`, the facet keeps its `ref_seqcol_digest`, so a second
`resolve_ref` refuses it ("already read"), and the joined rows can never be read: measured on a synthetic
split site whose two rows were both `match`, the joined row came back `ref` null, `ref_check = not_run`,
and `mint_vrs_ids` withheld all three of its alleles as `ref_unchecked`, while the unsplit row beside it
kept `match` and its ids. `liftover` leaves such a row behind (`ref_unread`) and `normalize` does not
align it. Nothing is wrong in what is stated: `not_run` is the honest value for a REF no step compared.
The order "rejoin before `resolve_ref`" is written in `normalize`'s docstring, not enforced. (A group
with a row spelled against an `N`, `ref_check = unstated`, is no longer joined at all, audit 14: there
`not_run` would have made its sequence ALT readable again.)

Repairs, none made:

- Carry `ref` and `match` over when every row of the group is `match` and the joined REF is the
  longest row's own `ref`. A joined row's `ref` then comes from the split rows' reads, not a read of the
  joined span; that holds where the REFs are prefixes of one another (the join requires it), but the
  row's `ref_check = match` would describe a comparison made on other rows, which is a check reported
  for a representation it never saw. VRS ids would still need minting again.
- Let `resolve_ref` read the `not_run` rows of a facet already read, against its own collection. It
  breaks the rule that a facet is read once, the guard against reading one facet against two
  collections, and needs a rule for which rows a second pass may touch.
- Make `rejoin_split` refuse a facet whose REF was read, as `resolve_ref` refuses one twice, naming the
  order. The cheapest and the most honest, but a caller holding only a resolved facet (a store's) has no
  way back to the split rows' parent short of reloading the file.
