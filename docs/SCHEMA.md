# Frame schema: proposal (RM8)

Status: **signed off 2026-09-30 as schema v0.0.1; now v0.0.22.** Every per-format loader (RM2) emits the frames
described below. Pieces found missing while building are added here as they come up, and each change
bumps the schema version (`just_vcf.schema_version` in the frame metadata).

The charter (docs/CHARTER_DRAFT.md) is itself still a draft, so each decision below cites the charter
invariant it leans on. If the charter changes, those decisions are the ones to reopen.

Each decision lists its options with what each one costs, and a recommendation. The recommended
options together make up the "At a glance" section at the end.

## Facts the draft rests on

Probed against polars-bio 0.36 on 2026-09-29 with a synthetic two-row VCF (one multi-allelic SNV, one
gVCF block):

- `pb.scan_vcf` with default options returns **1-based closed** coordinates: `POS=100` gives
  `start=100, end=100`, and a block with `INFO/END=300` gives `start=200, end=300`. The same file through
  `pb.register_vcf` + `pb.sql` returns **0-based half-open** (`start=99`). `scan_vcf` takes
  `use_zero_based`. The loader passes it explicitly and checks it against a known row. It never relies
  on the default.
- Multi-allelic ALT is joined with `|` (`G|T`). FILTER `.` decodes as `""`.
- With one sample, FORMAT fields come out as top-level columns (`GT: str`, `GQ: i32`, `PL: list[i32]`).
  With two or more they come out as one `genotypes` struct of per-sample lists (`GT: list[str]`,
  `PL: list[list[i32]]`). `scan_vcf(samples=[…])` selects samples.
- For an indel, `end = pos + len(REF) - 1`: `ACG>A` at 100 gives `end=102`, `A>ATT` at 200 gives
  `end=200`.
- INFO/FORMAT keys that the header does not declare are dropped silently (polars-bio P1).

More probes on 2026-09-30, after reading VCF 4.5 (hts-specs `VCFv4.5.tex` as of 2026-02-24):

- A 4.5 reference block with FORMAT `LEN=14` at 100 gives `end=113`, which is correct.
- `<DEL>` with `SVLEN=30` at 200 gives `end=229`. The spec's computed END is `POS + SVLEN = 230`,
  because a symbolic SV's POS is the base *before* the event (polars-bio P3).
- A partially phased GT `|0/1/2` comes back as `0/1/2`. The leading phase indicator is dropped, so the
  first allele's phasing is lost (polars-bio P4).
- ALT `.` (a monomorphic reference record) decodes as `""`, the same as FILTER `.`. The `*` allele is
  kept (`C|*`).
- **A column's name does not say its namespace.** When INFO and FORMAT both declare `DP`, INFO/DP is
  `DP` and FORMAT/DP is `fmt_DP`. When only FORMAT declares it, FORMAT/DP is `DP`. So a bare `DP`
  column means the cohort's combined depth in one file and this sample's depth in another.

## Decision 1: row grain

What one row of the calls frame is.

| Option | Cost |
|---|---|
| **A. One row per site, as the source stated it.** `alts` is a list; allele-level values are lists parallel to it | Joins on allele identity need an `explode` first. It matches polars-bio, VCF and `sink_vcf` one to one, and the genotype and likelihood vector stay whole |
| B. One row per allele (split multi-allelics at load) | Joins on VRS ids are direct. But a `1/2` genotype and its six-entry likelihood vector belong to the site, not to either allele, so splitting either duplicates them or rewrites them into per-allele approximations. That is a transformation of the evidence, done before anyone asked for it |
| C. Both, stored | Two copies that can drift |

**Recommendation: A.** The per-allele view is a lazy function over the frame (`explode` on the parallel
lists), not a stored shape. Splitting and left-aligning stay RM6 operations that produce a *derived*
facet (charter 5), not something the loader does silently. Re-joining a file that arrived split one row
per ALT is one too: `rejoin_split` (RM270), under `split_signal` below.

## Decision 2: how the genotype is written

**REF.** Charter 10 allows only a REF read from the facet's assembly sequence. So `ref` stays null
until RM3 has read that sequence, for every input, VCFs included. What the file wrote goes in
`source_ref`, which is null for sources that write none (arrays). RM3 fills `ref` from the sequence and
records the comparison in `ref_check`:

- `not_run`: `source_ref` is there and has not been compared yet;
- `match` and `mismatch`: the outcome. A mismatch keeps both values side by side, and neither is
  overwritten;
- `unstated` (0.0.22, F23): `source_ref` is made only of `N` (any case) where the reference reads a
  base. The source stated no base: DRAGEN's and Canvas's copy-number VCFs write REF `N` on every record,
  since the caller never reads the reference. `ref` is read as for any row; nothing was compared, so it
  is neither `match` (a check that never ran reported as passed) nor `mismatch` (a contradiction the
  source never made). A real base beside an `N` (`AN`, `NA`) states that base and is compared as
  usual. Where the reference reads `N` too (a gap), the letters agree and the row is `match`, as
  before: the id minted there is the one vrs-python mints over `N` bases (RM271);
- null: there is no `source_ref` to compare (arrays).

What each reader does with `unstated` follows one rule: **a sequence ALT on an `unstated` row is not
read** (`schema.unstated_sequence_alt`). It was spelled against the `N`, not against `ref`, so set
beside the reference's base it would be an allele the source never wrote (`N` → `NAC` read as `C` →
`NAC`). A symbolic ALT (`<DEL>`, `<CN0>`, ...), no ALT, and the REF slot are read as on a `match` row,
with `ref` standing in for the base the source left out. Every copy-number file in hand has only
those (merge audit 13, 2026-10-28, found liftover, the writers and `query` reading the sequence ALT;
`concordance` and `consensus` followed on 2026-10-31, then the graded-evidence readers,
`sex_chromosome_evidence`, the panel-site reading behind relatedness and ancestry, and phase blocks):

- `liftover` lifts it as a `match` row (the target's base is compared with `ref`, read from the
  source assembly; on a plus strand an SV's padding base must still agree, `ref_mismatch_after_lift`),
  unless an ALT is a sequence allele: then it stays behind as `source_ref_unstated`. Lifted, that ALT
  would be compared with the target's base and could swap with it, re-indexing the genotype (`N` → `G`,
  `0/0`, over `A` where the target reads `G`, would become `G` → `A`, `1/1`), and the lifted row,
  whose `ref_check` is null, would mint an id `mint_vrs_ids` withholds on the parent.
  `mismatch` is `source_ref_mismatch`.
- `write_vcf` and `write_vcf_zarr` write it, REF the base `ref` holds, ALTs as written, unless an
  ALT is a sequence allele: that row is unwritable (`ref_unstated`; `resolved_only=True` skips it,
  counted), since `N` → `NAC` written as `C` → `NAC` would load back as a stated `match`. `mismatch`
  is unwritable as `ref_mismatch`. A written row's `gl`, `gp` and `ad` are left out as on any row
  whose `source_alleles` are not the written `[REF] + ALTs` (the source's REF slot is the `N`).
- `mint_vrs_ids` mints the REF slot (the reference's own allele) and withholds every sequence ALT as
  `ref_unstated`: `N` → `NAC` read as `C` → `NAC` would be a delins the source never wrote. A symbolic
  ALT keeps its own reason (`sv_no_vrs_model`).
- `normalize` leaves it as written, counted ("the source REF states no base"): trimming an ALT spelled
  against `N` against the reference's base would rewrite an allele.
- `query` never answers from it for a sequence allele, neither `stated` nor `site_genotyped`: its
  sequence ALTs were spelled against the `N`, and a REF the source did not state says nothing about
  which base the site holds. An allele asked at its POS, spelled by it or not, is `overlapping_call`
  (carried unknown); read as stated, `N` → `G` with GT `0/0` would answer `hom_ref` for a base the
  source never named. `counts` by `ref_check` lists it, as does the
  CLI's `check` (`ref_check` counts, `--json` included), which groups by the column.
- `resolve_ref` counts it apart in its warning ("… with a source REF of only N, which states no base
  to compare"), never among the mismatches, and never as an IUPAC ambiguity.
- `concordance` compares none of its sequence ALTs: each is counted in `excluded` as `ref_unstated`
  (2026-10-31; before, `N` → `G` was compared as `ref` → `G`, so it met another facet's `T` → `G` as
  a concordance the source never stated). Nor is its genotype read for another facet's allele at its
  position (the `site` key): that allele stays one-sided, `absence = other_alleles`. A symbolic ALT
  of the row, a copy-number call, meets by type and extent as over a stated base.
- `consensus` takes none of its sequence ALTs as a stated allele, and reads its genotype for no other
  facet's allele at its position: that facet is silent there as `ref_unstated` (`CONSENSUS_VERSION`
  5). The row still stands at its position, so the facet's silence is not read by a block or its
  absence rule there: a whole-genome facet would otherwise read zero copies where it wrote a call.
- `genotype_evidence`, `genotype_evidence_at` and `reference_support` grade no genotype on it:
  `grade_basis = ref_unstated`, `confidence_kind = unknown`, no `genotype_index`, alleles, vector or
  confidence, its GP or likelihoods included, which are over `[N] + alts` as much as its GT is
  (2026-10-31; before, `N` → `G` with GT `0/0` was hom-ref, and `reference_support` answered
  reference at a base the source never named, the case `query` already refused). The row still covers
  its position, so `reference_support` answers `call` there (`call_at_position`, or
  `call_spans_position`), never the facet's absence rule; among calls of one span it ranks as one
  stating nothing. `consensus` and the atlas's confidence bands read the same grade (`atlas` takes
  the lowest confidence a facet states at a position, and this row states none); `consensus` already
  read no such row, so its output does not move. Before `resolve_ref` the row's `ref_check` is
  `not_run` and it is graded as written: the rule needs the REF read.
- `sex_chromosome_evidence` (`karyotype.measure`) reads no GT slot on a sequence ALT of it: a row with
  such a slot feeds neither the het share, nor the allele balance, nor the depth line (2026-10-31;
  before, `N` → `G`, `0/1`, on X was a het of two bases the source never named). Its ploidy as written
  still counts, being the caller's statement of how many haplotypes it saw, not an allele. A slot on a
  symbolic ALT counts as any symbolic call does: a called `<DEL>` with GT `0/1` is one more diploid ALT
  call and a het, on `N` → `<DEL>,G` as on `N` → `<DEL>` (2026-11-03; before, the second ALT made the
  whole row count as none). `N` → `<DEL>,G` with `0/2` or `1/2` has a slot on `G` and counts no ALT.
  An array's letter measures read no letter of it either: `N` → `C,G`, `1/2`, from a chip VCF, is
  neither a two-letter base call nor a letter het (2026-11-03), while its call as written still counts
  (`called`, `no_calls`, one letter or two), as its ploidy does. On a facet not yet through
  `resolve_ref` the rule cannot apply (`not_run`).
- `site_genotypes.read_site_genotypes`, the panel-site reading `relatedness` and `ancestry` share,
  reads no genotype from it at a panel site: the row is not usable, and a site whose every row is
  such a row is unknown, `basis = ref_unstated` (2026-10-31; before, `N` → `G`, `1/1`, at a site whose
  effect allele is `G` read two copies, `observed`, as a stated `C` → `G` does, while its `0/1` and
  `0/0` were `unusable` only because the letter `N` is neither panel allele). It is still a row at the
  position, so the facet's absence rule is not read there; a stated row beside it at the same
  position is the one usable row. `ancestry` counts such sites apart (`SiteDosages.n_ref_unstated`,
  `AncestryEstimate.n_sites_ref_unstated`); `relatedness` reads only sites with a count, as it does for
  every unknown. A symbolic ALT is read as before (it is no panel allele, so `unusable`); before
  `resolve_ref` the row is read as written (`not_run`).
- `somalier_relate` hands somalier that same reading, one `0/0`, `0/1` or `1/1` record per panel site
  it gives a count and none where it reads unknown, and runs `somalier relate` without `--unknown`
  (2026-11-03; before, each facet's rows went through `write_vcf(resolved_only=True)` and a
  whole-genome facet was meant to be read with `--unknown`, which is a flag of `relate`, not of
  `extract`, and reads every sample's unknowns as homozygous reference: a site whose only row the writer
  withheld, `ref_unstated` among its reasons, had no record and would have read as zero copies).
- `cis_phased_blocks` takes no slot of it into a block, its REF slot included: the GT indexes
  `[N] + alts`, so a `0` names no base either, as `genotype_evidence` reads it (2026-10-31; before,
  `reference_alleles=True` took the REF slot's minted id, so `N` → `G`, `0|0`, joined both blocks
  stating the reference's base). A slot that would be a member withholds its block `ref_unstated`,
  `withheld_allele` naming it as `gt_alleles` holds it (`N` on the REF slot). By default a REF slot
  joins nothing and the sequence ALT is already withheld (`mint_vrs_ids`), so the default blocks do
  not move. A row with only a symbolic ALT is read as before: its REF slot joins with
  `reference_alleles=True`, and its ALT withholds the block under its own reason (`sv_no_vrs_model`).
  `write_vcf` reads the phase columns, not the blocks, so what it writes does not move.
- `facets.variant_classes` (`callset.variant_class_counts`) cannot keep the rule: it is counted at
  load, before any REF is read, when every row's `ref_check` is `not_run`, and `resolve_ref` keeps the
  count. `N` → `G` counts as the `snv` its lengths make it, a reading that holds whatever base the `N`
  stands for; what it cannot show is that the site varies (over a reference `G` it is no variant). A
  whole-genome facet whose only SNV-shaped rows were spelled against `N` would so read its silence as
  reference for SNVs. No file in hand writes such rows: every copy-number file writes symbolic ALTs
  or none, and its `sv` count is as stated. A recount at `resolve_ref` would change a load-time count
  in a resolve step, and is not made (`tests/test_unstated_ref.py` pins the count).
- No other reader: loaders, `rejoin_split` and `consensus` write `not_run` or null; `liftover` and
  `normalize` null it on a row they move or re-derive (`liftover` moves none with a sequence ALT); the
  facet guard (`tests/test_load_seams.py`) reads `facets` columns only.

A store written before 0.0.22 keeps its rows' `mismatch` on such a REF (the step to 0.0.22 never
rewrites a value; it names the facets holding those rows, "Store migration"), and the facet records
`RESOLVE_VERSION` 12, so the same facet resolved today (13) is another REF state of it.

RM3 never fills `ref` from `source_ref` or the reverse. That would be filling a value from the source
that checks it.

Array exports write no REF at all (charter 10: REF never comes from the genotype), so a VCF-style index
genotype (`0/1`) cannot be assigned to them until RM3 has read REF and derived the alts.

**How RM3 fills these slots (`resolve_ref`).** It returns the same facet, same `facet_id`, with the
slots above filled, and leaves the `Loaded` it was given unchanged. That is not an overwrite in the
charter's sense (invariant 5): every column it writes is one this schema defines as "not read yet"
until the sequence is read (`ref` null, `ref_check = not_run`, an array row's `alts` and `gt` null),
nothing the source wrote is touched, and the result is fixed by the facet and its assembly's
sequence, not by a choice. The same holds for an array's `assembly_label` when the facet stated none
(null is "unknown", and the rows against the sequence settle it, below); a stated label is never
replaced. The alternative, a derived facet (`relation = ref_checked`), was rejected:
every facet would then stay `not_run` for ever and its checked twin would be a second identity for the
same calls. A facet already read is refused rather than read again. Since the id stays, the facet
records the rule that read it: `facets.resolve_version` is the `RESOLVE_VERSION` (0.0.20), so a rule
change that rewrites rows is visible on the facet, and the durable store refuses to take the same
facet in another REF state as `present`. The rows also depend on the chip manifests an array's `I`/`D`
codes and a mixed strand's withheld rows are read with (`JUST_VCF_MANIFESTS_DIR`), which the rule's
version does not fix: `facets.resolve_inputs` (0.0.21) digests the ones that read rows, so the same
array resolved without manifests and then with them is two states too. Rejected: digesting `RESOLVE_VERSION` into the `facet_id`,
as `normalize`, `liftover` and `rejoin_split` digest theirs, which would make every resolved facet the
second identity rejected above.

- A row with `source_ref`: `ref` is the reference over the same span, `pos` to
  `pos + len(source_ref) - 1`. `ref_check` compares case-insensitively (VCF REF bases are), so a
  lowercase `source_ref` matches; a reference ambiguity code against a concrete source base (VCF 4.5
  writes `R` as `A`) is a `mismatch` and counted apart. A mismatch keeps both values.
- An array row whose alleles are all single bases with at least one called: `ref` is the base at
  `pos`; `alts` are the called alleles other than REF, **sorted by allele string** (an array states
  no ALT order, so two files of one chip agree); `gt` indexes `[ref] + alts`, with null for an
  uncalled haplotype. A hom-ref row gets `alts = []`. `ref_check` stays null: there is no source REF.
- A Living DNA row whose genotype is two allele sequences written back to back (`CCA`, bases only, F4;
  the vendor states that encoding, so only `living_dna` facets are split): of every split into two
  pieces that start with the same base (VCF's padding base), keep those where one piece equals the
  reference over its own span from `pos`. Exactly one such split gives `gt_alleles`, in the order
  written; several (a repeat: `AAAA` over `AAA…` fits `A/AAA` and `AA/AA`) or none leave the row as
  loaded. A split row gets `ref`, `alts = [the other allele]`, `gt` and `end_pos` only when exactly one
  of its two alleles equals the reference over its own span: a biallelic site's REF is one of its
  alleles, so the other is ALT (`C/CA` over `CG…`). Where both do (`C/CA` over `CA…`, and every
  deletion, whose one-base ALT always equals the reference: an insertion and a deletion fit alike), or
  the call is homozygous (`CA/CA` shows one of the site's two alleles, so it cannot say whether the
  site is `CA>C`, which makes it hom-ref, or `C>CA`, hom-alt), `ref`, `alts` and `gt` stay null with
  `gt_alleles` set. Not split at all in a facet whose strand is `mixed` (F13). `ref_check` stays null.
  `normalize` aligns and trims a split row that got `ref` like any row whose REF was read (a padded
  MNV `CTA`/`CGA` becomes `T>G` one base on). The split itself never leaves an insertion that could
  move left: its inserted bases would end in the padding base, and the text would then fit a second
  split.
- An array row whose genotype is `I`/`D` codes (23andMe, AncestryDNA, a Final Report: `I` the allele
  with the inserted bases, `D` the one without), on an assembly a configured Illumina manifest maps
  its probes to (RM281, `just_vcf.probes`): matched with a manifest indel probe by an rsid among its
  `ids`, on the same contig. The probe's source sequence names the inserted bases and about 60 bases
  of flank each side; the site is located in the reference by those flanks, within 200 bases of the
  probe's own span around its `MapInfo` (the flanks with the bases between them: the reference
  carries them, a deletion site, `I` is REF; the flanks joined: an insertion site, `D` is REF). Taken
  only when the flanks are found at exactly one place, on the strand the manifest's `SourceStrand`
  gives. The site is then spelled with REF starting at the row's own `pos`, shifted along a repeat
  where the export placed it elsewhere in one: the reference bases from `pos` up to the indel come
  before both alleles, and where there are none both alleles end with the base after the indel (VCF's
  padding base, on the right; `normalize` moves it left as it aligns). A site before `pos` that no
  shift brings there, or more than 64 bases after it, stays unread. The row gets `ref` (read from the
  reference, never from the manifest), `alts` (the called alleles other than REF, as for any array
  row), `gt`, `end_pos`, and `gt_alleles` spelled as the two sequences, so `gt` indexes them like any
  row's; `source_gt` keeps the codes. A row whose rsids name more than one located site, or none, is
  left as loaded, each reason counted. Not affected by `strand_check`: `I` and `D` read the same on
  either strand. `ref_check` stays null. **Gated per chip** (`probes.ID_CALIBRATION`): a vendor's `I`
  is matched with Illumina's by rsid only on a chip where that was measured to hold, the share of
  resolved homozygous rows that read hom-ref (`23andme_v5` 98.7%, `ancestrydna_v2` 98.4%, PGP-HMS,
  2026-10-23). Every other chip (`23andme_v3` and `23andme_v4` at 76%, with the same sites hom-alt in
  every export, and on v3 hom-ref in the same person's CGI genome; `23andme_v1_v2` on 6 rows; an `illumina:` Final Report, a platform family or an unknown
  chip, never measured) has its code rows left as loaded, counted `id_orientation_unverified`.
- **An insertion point** (a row over no reference base: `source_ref = ""` in CGI's native files and a
  CASAVA indel table; a Genos row with `end_pos` null) gets VCF's anchor base. Every loader places the
  point at the base before it (1-based `pos = begin` of a zero-length 0-based locus), so the anchor is
  the reference base at `pos`: `ref` is that base, it is written before every allele in `gt_alleles`,
  `alts` and `source_alleles`, in place (`gt` and the vectors stay indexed; a partly called `?A` becomes
  `G?A` and stays partial), and `end_pos = pos`. A point before a contig's first base (`pos` 0, a
  zero-length locus at `begin` 0) has no base before it: VCF's rule takes the base after, appended to
  every allele, at `pos` 1, the one case where `resolve_ref` writes `pos` (0 was the loader's
  placement of a zero-length locus there, not a position read from the file). The source wrote no REF
  base, so a CGI or CASAVA point's `ref_check` becomes null, never `match`; `source_ref` stays `""`,
  what the file wrote, and the alleles as written are the anchored ones less their first (last, at
  POS 1) base. A Genos point is then a spanned row (below). A point where no haplotype was called
  (`gt_alleles` all null: CGI's zero-length no-calls, 346,474 in one PGP-HMS masterVarBeta) is left
  as loaded, counted: it has no allele to anchor, and REF with no ALT would state a no-call at the
  anchor base. So is a point whose every called haplotype is the empty allele (CGI `""/""` and `""/?`,
  Genos `-/-`), counted apart (`RESOLVE_VERSION` 11): it says nothing is inserted between two bases,
  which is evidence about the gap, not a genotype at the base before. Anchored it read `REF=X`, no
  ALT, `0/0`, a hom-ref call the source never made, written to VCF beside the source's own call at
  that base (`G . 0/0` beside `G A 1/1`). No frame carries a statement about a gap without one: a
  `calls` row at `pos` is a genotype at that base, and a block spans bases. The row keeps the loaded
  `["", ""]`, which is the statement as the file made it, with `ref` and `end_pos` null. A point
  with an inserted allele on any haplotype (`""/GG`) is anchored as before: there the `0` is the
  haplotype without the insertion, as VCF writes it. Anchoring places the insertion where the file
  did; `normalize` left-aligns it.
  Measured on a PGP-HMS assembly's `var` against CGI's own vcfBeta of it (2026-10-26, GRCh37):
  263,689 of the 263,690 anchored points carrying an ALT have a vcfBeta record at the same POS with
  the same REF and the same ALT set; the one other has none there. Shifted one base either way, at
  most 18 agree, so `pos = begin` is the anchor.
- **A Genos row** (`reference.SPANNED_GENOTYPES`): Genos states each row's span (0-based half-open
  BEGIN/END) and spells its alleles on the plus strand, `-` the empty allele, and writes no REF. `ref`
  is read over `pos`..`end_pos`, whatever its length; `alts` are the called alleles other than REF,
  sorted by string, and `gt` indexes `[ref] + alts`, as for an array row. So `AT/GC` over a two-base
  span reads REF and gets an MNV, `T/-` a deletion with the empty ALT that `normalize` anchors, and an
  insertion site `-/A`, anchored above, `G/GA`. `ref_check` stays null.
- Left unread, each reason counted once: I/D codes no manifest probe resolves (above), whole
  no-calls (the probe's REF span is unknown), genotypes that could not be split, spans past the
  contig's end or at POS 0, insertion points where nothing was called, where nothing inserted was
  called (every called haplotype the empty allele), or past their contig's end, and
  every row of a contig that could not be matched with the reference (below), which includes the
  insertion points and Genos rows on it, left unanchored.

**Matching the reference to the facet.** Per contig the calls use, read through `canonical_chrom` on
both sides (so `chr1`/`1` and `chrM`/`MT` meet; array code `25` reads X). A GRCh37 or GRCh38 scaffold
also meets the reference under the other names NCBI's assembly reports give it (UCSC's
`chr1_KI270706v1_random`, GenBank's and Ensembl's `KI270706.1`, RefSeq's `NT_187361.1`; F19,
`ReferenceStore.contigs_by_chrom`), only where the reference contig's length is the report's and the
reference holds no contig of that name itself; the row keeps its source spelling in `chrom`. NCBI36's
report names its scaffolds by RefSeq accession, and gives a UCSC name to two of them only, the alt
haplotypes `chr6_cox_hap1` (`NT_113891.1`, 4,731,698 bases) and `chr22_h2_hap1` (`NT_113959.1`, 63,661),
which meet each other's name. It adds no other name for the rest (those GRCh37 kept, such as
`NT_113878.1`, already meet GRCh37's names):

- the facet stated the contig's length (`header_contigs`): the reference's length must equal it, and
  its md5 too where the header gives one;
- it did not (arrays): the facet's `assembly_label` must equal the reference's. The reference's label
  comes from its contig 1 length (GRCh36, GRCh37, GRCh38, T2T-CHM13), or from one the caller asserted
  when adding the FASTA. A reference holding two sequences named contig 1 (`chr1` and `1`) has no
  label from its lengths, since which is contig 1 is unknown. An asserted label is refused when any
  sequence named contig 1 has a length that names another assembly, the doubled contig 1 included;
  it is kept where every contig 1 length names it or names nothing. `GRCh37` cannot vouch for MT, since it covers hg19's chrM (16,571 bases) as
  well as rCRS (16,569), and `GRCh36` neither (NCBI36 and hg18 carry the 16,571-base chrM), so an
  array's MT stays unread. `T2T-CHM13` neither: `chm13v2.0.fa` carries the CHM13 cell line's own
  mitochondrion and the T2T analysis set's `chm13v2.0_maskedY_rCRS.fa.gz` carries rCRS, both 16,569
  bases.
- it stated neither lengths nor a label (an FTDNA export, a headerless table): the label is inferred
  from the rows (F9), and then matched as above. Per collection in the store, a deterministic sample
  of the facet's array SNV rows (at most 2,000 per chromosome on 1–22, X and Y, every k-th row by
  position; MT never) is scored: the share whose called alleles include the collection's base at that
  position, on the forward strand only. An SNV row here has no source REF, a span of one base
  (`end_pos = pos`) and single-base called alleles. A Genos insertion site, whose alleles are the
  inserted bases at the anchor before them (`end_pos` null), and a row called over more than one base
  are not SNVs at `pos` and are left out, as is an insertion point with an empty source REF. A row on a contig the collection lacks is not scored on it; a
  row past the contig's end counts against it. The best collection names its label when it was scored
  on at least 1,000 rows, reaches 60%, and stands 20 points above both the best collection of any
  other label (an unlabelled collection counts as another) and chance (the mean number of distinct
  called alleles over four). The facet then gets `assembly_label` with basis `inferred`, and the
  scores go in one warning line. Rows that do not settle it refuse the facet, with the scores in the
  error. Only array rows are scored: for a VCF without lengths or a build hint, `ref_check` would then
  confirm the very comparison that chose its build, so it is still refused.

**The strand of an array's genotypes (F13).** An array derives `alts` and `gt` on the reading that its
genotypes are on the forward strand, which every export claims. `resolve_ref` checks that claim on the
same row sample, against the collection it reads, through the heterozygous rows (two distinct called
bases): on a two-allele site a het carrying the reference base is safe on either strand (its other
allele is the complement, which reads the same flipped, or another base, which flipped would make three
alleles), and a het carrying the complement of the reference base and not the base is impossible on the
forward strand. With at least 200 such hets sampled, the facet's `strand_check` is `mixed` when more
than 2% of them carry the complement and not the base, else `forward`; with fewer it stays null (not
checked), the alleles are derived as before, and a warning says so. A VCF is not checked: its REF is
the source's own. In a `mixed` facet every array row keeps `ref`, and the rows that are not
strand-safe (homozygous calls, half-calls, and hets without the reference base) get `alts` and `gt`
null, each reason counted in the warnings: `AA` at a reference `A` is a forward hom-ref call or a
minus-strand hom-alt call of an A/T probe, and the rows cannot say which. No genotype is complemented
by the facet's rows alone, not even on a het that proves its probe minus: complementing exactly the
minus probes needs to know each probe's alleles.

A probe's own alleles decide it where a configured Illumina manifest gives them (RM281,
`just_vcf.probes.snv_catalog`): a withheld row of a `mixed` facet is matched with a manifest SNV probe
by an rsid among its `ids`, on the same contig, and the probe's two alleles are put on the plus strand
(the `SNP` column as written where `RefStrand` is `+`, complemented where it is `-`). The row is read
only when that pair holds the row's reference base (so the pair is on the facet's own plus strand,
whatever build the manifest maps to) and is not its own complement (`A/T`, `C/G`). Called alleles all
in the pair: forward, derived as written. All in the pair's complement: minus, `gt_alleles`
complemented, then derived; `source_gt` keeps what the export wrote. Otherwise, or with no probe (or
an rsid whose probes give two pairs), withheld as before; each outcome is counted in one warning.
`TT` at a probe whose plus pair is `A/G` is a minus-strand `AA` whichever column the vendor flipped
by, so this reads what F13's strand columns could not. Forward facets are not touched: reading them
by the pairs would rewrite rows of every export whose rsids a manifest shares, on evidence the strand
check already found sound.

An autosome (1–22) whose stated length or md5 differs refuses the whole facet
(`ReferenceMismatchError`): it is another assembly. Any other contig that does not match leaves its
rows unread, and the matching contigs are read: X, Y and MT, where assemblies of one build
legitimately differ (Complete Genomics' GRCh37 chrY has another md5 than Ensembl's at the same length,
F20; hg19's chrM against GRCh37's MT), an alt contig the reference lacks, a name the reference holds
twice. A length or md5 mismatch is counted per contig in the warnings, with what was stated and what
the reference holds. A facet of which nothing matches is refused. (`concordance` across two collections and
`liftover` keep their own rule, which refuses on X and Y too.) When a store holds several collections, the
one that matches the most rows is used; a tie is refused unless the caller names one.

| Option | Cost |
|---|---|
| A. VCF indices only (`gt: list[u8]`) | Impossible for arrays until RM3. The loader would have to either block on a FASTA or invent a REF |
| B. Literal alleles only (`gt_alleles: list[str]`, one entry per haplotype) | Works for every input at load time. But likelihood vectors are indexed by VCF genotype order, so they cannot be read without `[ref] + alts` next to them |
| **C. Both: `gt_alleles` always, `gt` when there is an allele list to index** | Two columns that must agree, enforced by a test. `gt` indexes the file's own allele list, `[source_ref] + alts`, whenever the source wrote a REF. So a VCF has `gt` at load time. Otherwise it indexes `[ref] + alts` once RM3 has read `ref` and derived `alts`, and until then it is null. RM3 never writes REF from `gt_alleles` |

**Recommendation: C.**

**`source_alleles`: what the vectors are indexed over.** A likelihood or posterior vector needs an
allele order, and REF is not always known when the row is loaded. An Oxford `.gen` row names `alleleA`
and `alleleB` and gives three posteriors, but does not say which allele is the reference. So each row
carries `source_alleles: list[str]`, the allele list in the order the source indexed its vectors:

- a VCF with full vectors (PL, GL, GP, AD): `[source_ref] + alts`;
- a VCF with only local-allele vectors (`LAA` with `LPL`, `LAD`, `LGP`, and no full vector):
  `[source_ref] + alts[LAA]`. Where both a full and a local vector exist, the full one is used, because
  the spec requires them to carry the same information;
- Oxford `.gen`: `[alleleA, alleleB]`;
- an array hard call with no vector: null.

`gl`, `gp` and `ad` (Decision 3) are indexed over `source_alleles` and are never reordered in place.
When RM3 reads REF and it turns out to be `alleleB`, `ref` and `alts` are filled and the vectors stay as
the source wrote them. A consumer that wants `[ref] + alts` order derives it, and that derivation gets
a test.

How the awkward cases land under C:

- **Half-call** (`1|.`, `0/.`): a null inside the list, e.g. `gt_alleles = ["G", null]`. The list
  length is the number of haplotypes the source wrote. It is not an asserted ploidy.
- **Whole no-call** (`./.`, `--`, `00`): every entry null. The row is kept when the source wrote it,
  because the source said it looked and could not tell, which is not the same as not looking.
- **Hemizygous X/Y in arrays**: kept exactly as written (`["A"]` or `["A", "A"]`). An array writing
  `AA` on male X is not saying "homozygous". The facet's ploidy handling is a known open issue, not
  something the loader corrects.
- **23andMe `I`/`D`**: kept verbatim as symbolic alleles at load. `resolve_ref` spells them as
  sequences where a configured chip manifest's indel probe resolves them (above, RM281); until then,
  and where none does, no VRS id is minted for them.
- **CGI partial no-call inside an allele** (`GC?`): kept verbatim in `gt_alleles` (`"GC?"`), with no VRS
  id. The `?` marks the bases CGI could not call. Because `gt` indexes `[source_ref] + alts`, a partial
  allele is in `alts` too (the native `var` and `masterVarBeta` files; CGI's VCF writes such a
  haplotype `.`).
- **No anchor base** (CGI native files): CGI writes an indel over its own span, without VCF's padding
  base, and REF only comes from the sequence. So an allele can be the empty string: a deletion of the
  whole span is `""` (`source_ref = "T"`, `alts = [""]` at `pos = begin + 1`). A locus of zero length
  (an insertion point) has `source_ref = ""`, sits at the base before it (`pos = begin`, where VCF puts
  the anchor) and has `end_pos` null until `resolve_ref` reads that base and writes it into every allele
  ("How RM3 fills these slots", above), as for Genos insertion sites. An empty
  string in `alts`, `source_alleles` or `gt_alleles` occurs only in such unanchored rows, and in
  `source_ref` on every insertion point, which keeps what the file wrote once anchored (its `ref` then
  holds the anchor); that is how RM3 finds the rows it has to pad. An insertion point CGI could not
  call stays a row with every haplotype null: it covers no base, so it cannot be a block, and it has
  no allele to anchor. One that calls only the empty allele stays unanchored too (above).
- **The empty allele has one spelling: `""`** (`schema.EMPTY_ALLELE`), VRS's zero-length literal, in
  every loader: CGI's unpadded indels above, an upload's ALT `-` (below), and Genos's `-` (an insertion
  site's empty REF allele, or a deletion's ALT; `source_gt` keeps the `-/A` the file wrote). Every reader
  takes `""` as a sequence allele of length zero: it is an `indel` in `variant_classes`, a site in RM4's
  count, the insertion point `resolve_ref` anchors, the deletion `normalize` anchors, and the allele
  `mint_vrs_ids` justifies (a resolved `CT>""` at 31 gets the id of `GCT>G` at 30). `write_vcf` refuses
  it, REF or ALT, until it is anchored (VCF has
  no empty allele), and `liftover` leaves it unlifted (`allele_not_sequence`): run `normalize` first.
  A `-` that a format uses for a no-call (FTDNA's `-A`, a Final Report's `-`) is no allele and stays
  a null haplotype.
- **The `*` allele** (VCF: "allele missing due to overlapping deletion") is kept as `"*"` in `alts` and
  `gt_alleles`. It says the haplotype is covered by a deletion called at another record. It is not
  a reference call, and it gets no VRS id. `0/2` over `ALT=C,*` must never be read as a het with a
  reference haplotype and nothing else.
- **A deletion written as ALT `-`** (some uploads, F12: REF holds the deleted bases at POS, with no
  anchor base) is read as the empty allele `""` in `alts`, `source_alleles` and `gt_alleles`, counted
  in a load warning. It is an unanchored deletion, which VCF has no allele for, and the anchor needs
  the reference base before POS, so the loader cannot anchor it. `normalize` anchors it once
  `resolve_ref` has matched its REF; until then `write_vcf` refuses it.
- **A monomorphic record** (ALT `.`) has `alts = []`, which is a stated "no alternate allele". A null
  `alts` means the source did not say, which only happens for non-VCF sources (an array row before
  RM3). ALT is a required VCF column, so polars-bio's `""` always means `.` and becomes `[]`.
- **`<NON_REF>`** is written as `<*>` in `alts`. VCF 4.5 makes them aliases.
- **`source_gt`**: the file's genotype text, byte for byte, on every row that had one. It costs one
  string column. In exchange, every derived column can be checked against what the file said, and no
  representation choice here ever loses information. It cannot be polars-bio's `GT` column, because
  polars-bio rewrites `|0/1/2` as `0/1/2` (P4). The loader has to get the raw text itself. `phased[0]`
  is read from that text, never from polars-bio's GT.

Phasing is per haplotype, not per row. Since VCF 4.4 each allele in GT carries its own phase
indicator (`|0/1/2`: the first allele is phased and the other two are not), and `PSL` gives each allele
its own phase set. So:

- `phased: list[bool]`, parallel to `gt_alleles`. An element is null where the source did not say
  (an array export).
- `phase_sets: list[str]`, parallel to `gt_alleles`: the raw value the file wrote, with no prefix. A
  `PS` integer is written as its decimal string on each phased haplotype; a `PSL` entry is taken as is.
  Unphased haplotypes are null.
- `phase_set_kinds: list[Enum{ps, psl, cgi_haplink}]`, parallel to `phase_sets`: which field the value
  came from. VCF 4.5 says PS and PSL "are not interoperable", so a PS `5` and a PSL `5` are different
  sets, and a CGI hapLink `5` is a third.
- **CGI's `hapLink`** (native `var` and `masterVarBeta`, kind `cgi_haplink`, just-vcf's own term) is
  an id on one haplotype of a locus, shared by the haplotypes of other loci on the same molecule: a
  linked slot has `phased = true` and the id verbatim in `phase_sets` (an integer before format 2.5,
  mostly `Phased_a_b_c` in 2.5); an unlinked slot has both null, since the file says nothing about it,
  and a row with no linked slot has all three columns null. Like a PSL name it is a molecule, not a GT
  slot; unlike PSL it is keyed per chromosome, `(facet_id, chrom, cgi_haplink, value)`: on a PGP-HMS
  assembly no id leaves its chromosome, and CGI's own VCF writes the same phasing as PS sets of two
  ids each (25,696 sets, 51,392 ids). A slot written `?` or partly called keeps its link.
- **The key of a phase set depends on its kind.** A PS is scoped to its chromosome ("on the same
  chromosome and have the same PS value"), so its key is `(facet_id, chrom, ps, value)`. `PS=100` on
  chr1 and on chr2 are two sets. A PSL name is global to the file ("records with different CHROM but
  the same phase-set name are considered part of the same phase set"), so its key is
  `(facet_id, psl, value)`, and one set can span chromosomes.
- A phased haplotype with no PS or PSL has `phased = true` and `phase_set` null (kind `ps`). It belongs
  to its chromosome's one implicit set (VCF 4.5: "All phased genotypes that do not contain a PS
  subfield are assumed to belong to the same phased set"), key `(facet_id, chrom, ps, null)`.
- **A PS the loader cannot read withholds the phase** (`phased` null where it was true, the phase-set
  lists null), never the implicit set: that set is for a record that wrote no PS, and one whose PS was
  not read did write one. The cases, each counted in a warning that says the phase was withheld: a PS
  value left out as not a whole number (Decision 9), and a negative PS (`-5`, `-1`), since VCF's PS is
  a non-negative integer (2026-10-25; the list arithmetic read `-1` as the implicit set). PS `0` is a
  set like any other. And a PS the header declares but the loader does not carry (a Type its version
  refuses, `PS Type=String` under 4.3; a key the retry left out, Decision 9): its values are skipped
  unread, so which records wrote one is unknown, and every phased record has its phase withheld
  (2026-10-25; it read `0|1:100` and `1|0:200` as one set across the chromosome). That includes a
  haploid call, phased by VCF's reading, which may have written a PS too. A record whose PSL names a
  set keeps it, since PS is ignored there. A file that declares no PS keeps VCF's rule: its phased
  records are in the implicit set.
- **How the VCF loader reads PSL** (VCF 4.4+, `Number=P`: one value per allele in GT). On a record
  whose PSL names a set, each slot is read from its own entry: a phased slot with a name is in that
  set (kind `psl`); a phased slot with `.` is in the implicit set above (kind `ps`, set null), since it
  has no PS subfield either; an unphased slot is in none, and a name written there is not read
  (counted: VCF 4.5 says an unphased allele has `.`). A PSL that is not one value per GT allele
  withholds the record's phase (`phased` null on every slot, counted), since the set each phased
  allele is in cannot be read. PS is ignored on a record whose PSL names a set (counted; VCF 4.5
  forbids both).
- **The first slot of a VCF 4.4+ record** has `phased` null, because polars-bio drops the leading
  indicator (P4). A PSL name there says it is phased: VCF 4.5 gives an unphased allele `.`, so the
  loader reads `phased[0] = true` from the name. Rejected: leaving it null, which withholds every PSL
  block whose first-slot allele is an ALT (`phase_unknown`), although the file stated the set. A `.`
  there says nothing (a phased allele may have no name), so it stays null.
- **PSO and PSQ** are carried beside PSL, on the slots in a PSL set: `phase_set_ordinals` (the
  traversal order along a derivative chromosome, which VCF 4.5 asks for when a set holds a
  symbolic SV or a breakend) and `phase_set_qualities` (PSQ, the phasing evidence of each slot,
  Decision 3). A list that is not one value per PSL entry is dropped and counted. PQ, the PS
  counterpart of PSQ, is not carried.
- polars-bio refuses a header that declares `Number=P` (P10), so the loader reads such a file
  through its rewritten copy with those keys declared `Number=.` (`REWRITE_VERSION` 9), with a warning;
  the count is checked against GT here. A PSL, PSO or PSQ declared `Number=1` is read the same way
  (`REWRITE_VERSION` 12, 2026-10-25): VCF 4.4 and 4.5 define all three as one value per GT allele, and
  a diploid row's two values do not parse under a single-value declaration, which failed the load.
  The per-row counts decide, as for `P`, so a haploid row and a diploid one both read. Another fixed
  count (`Number=2`) reads as a list already and is read as it is. Rejected: not carrying a key so
  declared, which would drop the values of every file whose rows are haploid (they load today), and
  judging the declaration by the first rows, since a later diploid row would still fail the load.

The VRS CisPhasedBlock is built from `phase_sets` at export, by `cis_phased_blocks` ("CisPhasedBlock
export" below). It is not stored per row.

**Genotype text that is not allele indices** (2026-10-25; the key=value table and Nirvana JSON, whose
genotypes the loader reads from text itself). A genotype is readable when each slot is an index or `.`,
separated by `/` or `|`. One that is not (`N/A`, `x`, `-`, `1:2`, `0/x`) is withheld whole: `gt`,
`gt_alleles` and `phased` are null, `source_gt` keeps the text, and a warning counts the rows. Its
readable slots are not read either: `0/x` as `[0, null]` would be the frame of `0/.`, a half-call the
file did not write, and the frame could no longer tell a slot the file left uncalled from one the
loader could not read. Withheld, the two stay apart: `gt` null with `source_gt` set is unreadable, a
list of nulls is a written no-call, and `gt` null with `source_gt` null is no genotype. The array
loaders withhold a genotype they cannot read the same way (`gt_alleles` null, the text kept).
Rejected: keeping the readable slots, which is the honest reading only if the text means what its
readable part seems to, and a text the loader cannot read is the evidence that it may not. The VCF
loader never reaches the case: polars-bio parses GT itself and refuses the whole file on `0/x`
(`HeaderRefused`, "Error reading FORMAT field 'GT'"), and since GT is never left out, the load fails.
Each genotype is counted under one reason at most: one not readable is counted as that
alone, and only a readable one is checked against the row's ALT list. So `1:2` or `x/5` on a one-ALT
row is unreadable, never also past the ALTs, though its digits name an allele the row lacks. The
indices are read from the text's digits, never from `gt`: `gt` is a `u8` list, so `1/300` on a row
with fewer than 300 ALTs reads `[1, null]` there and is still counted past the ALTs, and an index too
long for 64 bits is past them too (it failed the key=value table's load). The VCF loader counts the
same, with the same warning, from the GT text polars-bio hands back (2026-10-26; `0/300` loaded as
`gt [0, null]` unwarned), and its frame is the other loaders' (`gt [0, null]`, `gt_alleles` null
there). In a VCF 4.4+ file, which has no `source_gt` (P4), that row looks like the half-call `0/.`, and
only the warning's count says otherwise.

An index the row has but `gt` cannot hold (2026-10-26) is a third reason, counted apart in all three
loaders: on a row with 300 ALTs `0/300` names an allele that exists, yet `gt` reads `[0, null]` and
`gt_alleles` `[REF, null]`, the frame of `0/.` again, and it was counted nowhere, since it is not past
the ALTs. Its warning says `gt` cannot hold the index (255) though the row has the allele. A genotype
past the ALTs (`0/301`, `300/301` over 300 ALTs) is counted as that alone, so each genotype still
falls under one reason at most. Not done: widening `gt` to `u16`, which would hold
the index; that retypes a column every consumer reads, which §8 of CLAUDE.md makes a major change (a
retype, like a rename), allowed before 1.0 when the changelog says so. It is the person's call, and
the warning stands until it is made.

## Decision 3: graded evidence (charter 7)

The user's words: "a bayesian based fuzzy 0 (unknown) .. PHRED-like (confident genotype). … each
position is probably an n-vector of probabilities or log or something, tbd; No binary call-no call
definitely."

Sources give different quantities, and they are not interchangeable. PL and GL are **likelihoods**
(the probability of the reads given each genotype). GP and imputed dosages are **posteriors**, and
those already contain a prior (a population panel). Converting a posterior into a likelihood is wrong,
and so is the reverse without the prior. So:

| Option | Cost |
|---|---|
| A. One vector column, everything converted into it | Mixes likelihoods with posteriors under one name. That is exactly the kind of approximation the charter forbids |
| **B. Separate columns per quantity** | `gl: list[f32]`: log10 likelihoods (from PL as `-PL/10`, computed in Float64 and rounded to Float32 once, so PL and GL stating one likelihood load equal; or GL directly). `gp: list[f32]`: posterior probabilities (GP, Oxford `.gen` triples). Both are in VCF genotype order over `source_alleles` (Decision 2). Which of them is non-null says which quantities the source gave, and a row can have more than one. A consumer has to handle two columns, and the schema is honest about what each one is |
| C. Keep native columns (`PL`, `GL`, `GP`, `DS`) as-is | No normalisation at all. Every consumer then re-learns every caller's dialect |

**Recommendation: B.** Scalars that callers commonly give sit beside it: `gq: u16`, `dp: u32`,
`ad: list[u32]` (parallel to `source_alleles`), all nullable.

Sources for `gp`: `GP`, and `PP` converted as `10^(-PP/10)`. Both are posteriors, so this is a change
of scale, not a change of quantity. One trap: before VCF 4.3, `GP` was phred-scaled, and 4.3 made it a
probability. Real files do not follow the version: imputation tools write probabilities under 4.1 and
4.2 headers (9 of the 124 PGP-HMS VCFs). So the values decide when they can: any value above 1 means
phred (`gp = 10^(-GP/10)`), and values all within [0, 1] that sum to 1 on every row mean probabilities.
The warnings name a header that disagreed. Values that settle neither (within [0, 1] but not summing to
1) are left to the header: under VCF 4.3 or later, which defines GP as a probability, they are read as
probabilities and the rows that do not sum to 1 are counted in a warning (posteriors converted from
rounded phred, as just-vcf's own 4.3 export writes them, F10). Under 4.2 or earlier, or a header that
states no version, they leave `gp` null, with a warning. The FORMAT line's Description is not read: it
is free text, not a declaration. The loader never guesses the scale.

`dp` and `ad` come from FORMAT (this sample), never from INFO (all samples combined). INFO/DP and
FORMAT/DP are different quantities with the same name. The loader resolves each key by the header's
namespace and never by the column name polars-bio gives it.

A vector is carried only where it can be indexed: `ad` needs one value per allele, `gl` and `gp` one
per genotype of the row's ploidy (checked for haploid and diploid rows). A row whose vector has
another length gets null there, counted in the warnings. A key declared with a fixed count is a
different quantity under the same name and is not carried at all (VarScan2's AD, `Number=1`, is the
ALT depth; CGI's AD, `Number=2`, counts reads per GT slot); it is listed in `dropped_fields`.

CGI's native `masterVarBeta` states its counts per allele, not per slot, so it does fill `ad`, over
`source_alleles = [source_ref] + alts`: REF gets `referenceAlleleReadCount`, each ALT the
`allele<n>ReadCount` of the slot that carries it, and `dp` is `totalReadCount`. On a PGP-HMS file a
slot carrying the reference always counts what `referenceAlleleReadCount` does and both slots of a
hom ALT count alike; their sum passes `totalReadCount` (the same reads twice), so slots are never
summed, and since CGI counts a read for every allele it supports, `ad` can add up past `dp`. A row
that calls no haplotype (an uncalled insertion point, whose counts the file leaves empty), whose
counts are empty (a half-called reference), contradict that rule, or belong to a partly called
allele (`?A`) has `ad` and `source_alleles` null, counted in a warning by reason.

A dosage (VCF `DS`, the expected count of the second allele) is one number. It cannot be turned into a
posterior triple, because two unknowns do not follow from one number, so it gets its own column:
`dosage: f32`, with `gp` left null unless the source also gave GP.

`evidence_basis` says only where the row's evidence came from, not which quantities it has. `gl`, `gp`
and `dosage` already say that. It is an enum:

- `assay`: the source wrote a record at this position (a variant call, a monomorphic record, an array
  genotype, an imputed row);
- `ref_block`: materialised from a reference block. A gVCF block that carries PL is `ref_block` and still
  has `gl`;
- `absence`: materialised from a whole-genome callset's silence (charter 11).

A hard call is a row whose `gl`, `gp` and `dosage` are all null. It has no value of its own in the
enum.

**What an array hard call gets** is the user's open "tbd". Options:

| Option | Cost |
|---|---|
| **A. No vector: `gl`, `gp` and `dosage` null** | Honest now. Consumers that need a number have none until C lands |
| B. A fixed pseudo-PHRED per chip (e.g. "Q20") | A made-up number that then looks measured. Rejected by the no-fabricated-values rule |
| C. A per-facet, per-region error rate measured by concordance with other facets of the same person (charter 8), turned into a vector later | The right end state, but it needs a second facet and the concordance machinery, so it cannot be ready at load time |

**Recommendation: A now, C as the path**, recorded as a facet-level measurement and not written into
rows by the loader.

## Decision 4: reference blocks and absence (charter 11)

gVCF `<NON_REF>` blocks and CGI `var` ref/no-call blocks are intervals. A WGS gVCF has hundreds of
millions of them.

| Option | Cost |
|---|---|
| A. Blocks as rows in the calls frame, with a `kind` column | Satisfies "one frame out" literally. Every variant query then scans the blocks too, and a block has no alts, no VRS id and no genotype in the usual sense |
| **B. Two frames per facet: `calls` and `blocks`** (the Hail VDS layout) | A second frame to know about. `blocks` is small in columns: `facet_id, chrom, pos, end_pos, kind (ref / no_call), gq, dp_min, gl`. `gl` is the block's own likelihoods over `[REF, <*>]` (a gVCF block carries PL, e.g. `0,60,900`), because that *is* the graded reference evidence (charter 7). It is null only where the source gave none (CGI blocks) |
| C. Expand blocks to one row per base | Billions of rows. Rejected |

**Which VCF records are blocks:** a record goes to `blocks` only when every ALT is `<*>` (`<NON_REF>`
included), or when it has no ALT (`.`) and an `INFO/END` extent, which is how Illumina's gVCF writes a
reference block. Either way its GT must assert reference or nothing; a GT calling `<*>` stays a call. A record with `ALT=TC,<*>` is a variant call that also carries the `<*>` allele, and it
stays in `calls` with `<*>` in `alts`. A monomorphic record (ALT `.`) stays in `calls` too.
An `ALT .` record to `END` that carries a copy number (FORMAT or INFO `CN`) is a copy-number segment,
not a block (RM282): Canvas writes its copy-neutral segments so (`GT 0/0`, `CN 2`, 90 of 182 records
of a PGP-HMS Canvas file, which before RM282 loaded as 90 `ref` blocks and `absence_means = blocks`).
It states the region's copy number and nothing about the small variants in it, so it is a `calls`
row with `alts = []`, `end_pos = END` and `copy_number`, and the facet's absence stays `unknown`.

Complete Genomics' native `var` and `masterVarBeta` files state their blocks directly: a `ref` locus is
a `ref` block, and `no-call`, `no-ref` (an assembly gap) and `PAR-called-in-X` (chrY stretches whose
calls are on chrX) are `no_call`, as is a variant locus whose every haplotype is no-called. Each block
is `pos = begin + 1`, `end_pos = end`. A locus of zero length is never a block (above). The files tile
the genome, so `absence_means = blocks`.

Complete Genomics' VCF (`vcf_cgi`, or any VCF declaring its alleles) adds two block kinds, both
`no_call`: a `<CGA_NOCALL>` record spans POS to INFO `END`, and, where the header declares
`CGA_NOCALL`, a record with ALT `.`, no `END` and no allele called (`./.`, `.|.`, or a haploid `.`,
which polars-bio 0.36 returns as a null GT) spans its REF. Its
`<CGA_CNVWIN>` records are 2 kb copy-number windows, not genotypes, and are left out with a count. A
half-called reference span (`0|.`, ALT `.`) stays a call.

What a CGI VCF leaves out, CGI called reference: on a PGP-HMS genome its vcfBeta covers every base the
same assembly's `var` file does not call reference (0 bases missed of 317.7M, 2026-10-22), and it
writes no reference-called record. So a VCF whose caller line names CGAPipeline (`vcf_cgi`, and an
upload that lost its `##center` line) gets `absence_means = reference` when its rows show a whole
genome (`callset_type = wgs`, RM4). Otherwise (a region cut from one, or too few rows to tell) it keeps
`blocks`, with a warning: claiming reference for the regions a cut took away would read a deletion of
data as a reference call. The native `var` and `masterVarBeta` files keep `blocks`: they state their
reference calls as `ref` blocks, so nothing they omit was called.

A 4.5 block states its length per sample in FORMAT `LEN`. INFO `END` is deprecated and only a fallback
when `LEN` is missing. Blocks keep the source's `MIN_DP` as `dp_min`, and not the block's average `DP`.
A block with an average depth of 25 can still hold an uncovered base.

**Writing blocks back out** (`write_vcf(…, blocks=True)`, RM6). The frame holds no REF base, so a
block record takes the base at `pos` from the collection `resolve_ref` read the facet's REF from
(charter 10); a facet whose REF was not read is refused, and a block on a contig that collection
cannot vouch for (not in it, named twice, a stated length that differs) or past its end is counted
`block_ref_unread`, never written with a made-up base. The record is `REF <*>` (or GATK's
`<NON_REF>`), INFO `END = end_pos`, and FORMAT GT, GQ, MIN_DP (`dp_min`) and GL (`gl`) where the row
has them. The frame keeps no block GT, only whether it called anything (`kind`), so GT is written
`0/0` for `ref` and `./.` for `no_call`, at the ploidy the block's GL implies (k + 1 values for
ploidy k), else the ploidy of the contig's calls, else of the facet's; that ploidy is the writer's
choice, and `load` drops it again. Blocks and calls are merged by position one contig at a time,
calls first where both start at one position. `load` reads the file back to equal `blocks` and
`calls` frames.

**Recommendation: B.** The charter's "one frame out" is amended to name the three frames (facets,
calls, blocks): a gVCF cannot fit in one calls frame. Every format lands in the same `calls` schema, and
a format that has blocks also fills `blocks`.

**Absence is a facet property, not rows.** The facet carries `absence_means`: `unknown` (arrays,
imputed content, VCFs of array data), `weak_reference` (whole-genome callsets without blocks),
`blocks` (gVCF, CGI's native files: only positions inside a `ref` block are reference), or `reference`
(CGI's VCF, above: every position the file omits was called reference by its producer). A position a
block covers is not omitted; its block says what it is. `reference` is just-vcf's own term (charter
2), and it is stronger than `weak_reference`, whose silence is only a caller not reporting a variant
there. Nothing is materialised
until a consumer asks for a position. When it does, the answer comes out as a `calls`-shaped row with
`evidence_basis = absence` or `ref_block`.

**Absence is per variant class.** `weak_reference` and `reference` hold only for the classes in the facet's
`variant_classes`; for any other class, absence is unknown. A small-variant caller that wrote no
symbolic ALT says nothing about SVs. And a whole genome cut down to one class keeps
`callset_type = wgs` (the assay did not change) but gets a `content_filter` (`snv_only`, `indel_only`,
`sv_only`). Its absences say nothing about the classes cut away, even where a residue of them is
present (`bcftools view -v snps` keeps a multi-allelic site whose other ALT is an indel), and the cut
holds them back for the kept class too.

`absence_is_reference(facet_row, variant_class, force=False)` is the rule in code, for just-vcf and
for consumers that restore hom-ref themselves (just-prs, just-dna-lite). A position is absent when no
`calls` row and no `blocks` row covers it. It returns `True` only when
`absence_means` is `weak_reference` or `reference`, the class is present in `variant_classes`, and the facet has no
`content_filter`. `force=True` lifts the hold for the filter's kept class only. Every other case is
`None`, never `False`, including `blocks` (a gVCF decides per position, in its `blocks` frame). Any
hom-ref restoration from absence that just-vcf grows must ask it, and refuse a facet with a
`content_filter` unless called with `force=True`.

`absence_support(facet_row, variant_class, force=False)` gives the same answer with the gate that
decided it (`AbsenceReason`, just-vcf's own vocabulary), so a refusal says why. It reads as reference
under `whole_genome_silence` (`weak_reference`), `omitted_called_reference` (`reference`) and
`content_filter_forced` (either, on a filtered genome, for the kept class, with `force=True`). It is
unknown under `outside_blocks` (`blocks`); under `array`, `targeted` (an exome or a panel), `imputed`,
`callset_type_unknown` (a null type) or `absence_not_carried` (a whole genome: a derived facet whose
producer withheld what its absence means, as `liftover` and `consensus` do), from the facet's
`callset_type`, when `absence_means` is `unknown`; and
under `class_not_counted`, `class_not_called`, `content_filter_held` (the kept class, without
`force`) and `class_cut_by_filter`, which a class the filter cut gets even where the facet wrote none
of it (the cut is why). `absence_is_reference` is its `reads_reference`, so the two cannot drift.
`reference_support` ("Graded genotype evidence", below) is the hom-ref restoration just-vcf grows: per
asked position it answers from a call, a block or this rule, and never reads an absence as reference
where `absence_support` does not.

**On an array, an absence has a cause (RM281).** An array's absence stays `unknown`, but
`chip_absence(loaded, sites)` says why, per asked site (`chrom` canonical, `pos` 1-based), from the
probe set of the facet's `chip` on its `assembly_label` (`just_vcf.probes.probe_set`):
`absence_cause` is `in_export` (the export has a genotype there: an allele called in any of its rows
at the site, or a text it could not split), `no_call` (the chip probed the site and the export gives
no genotype: a no-call it wrote, `--` or AncestryDNA's `0 0`, every allele null; or, from the probe
set, a site the export has no row at, a no-call it dropped or a probe left out of this release) or
`not_on_chip` (no probe there). A written no-call is `no_call` with or without a probe set, since the
row shows the chip looked. The sites with no row get null when there is no probe set (chip or
assembly unknown, or nothing configured or learned), and every site when the facet is not an array
(an array is a `genotype_table` facet whose `callset_type` is `array` or null: a VCF, a site table, an
annotation export and a CGI callset are not, though they state no VCF version either); the answer
carries the reason. Neither cause is ever evidence of reference. A probe set is one of two kinds
(`ProbeSetKind`): `manifest`, the chip's own Illumina manifest when a Final Report names it
(`illumina:<manifest>`) and it is configured (`JUST_VCF_MANIFESTS_DIR`), its probes on that build
(`Chr 0` probes are unplaced and left out); or `learned`, the positions at least half of the exports
of that chip carry, written into the cache by `learn_probe_sets` from three exports or more, for a
vendor that publishes no manifest (23andMe and AncestryDNA write a row per probe, no-calls included).
Not the union: one mislabelled export would put its whole coordinate system on the chip. A
platform family named by marker overlap (`illumina_gsa`, `illumina_omniexpress`) has none: another
version's manifest would put probes on the chip that it may not carry. A learned set is as complete as
the exports it was learned from; `of_exports` and each site's `exports` count say how many. The cause
is computed on request, not stored: it changes when a probe set is learned or configured, and the facet
does not.

`content_filter` comes from the strongest evidence that settles it, and each weaker piece that
disagrees is a warning naming both sides; the strongest is kept:

1. **The rows** (`content_filter_basis = rows`): at least `MIN_FILTER_ALLELES` (10,000) ALT alleles of
   one class, and the other classes together at most `MAX_FILTER_RESIDUE` (1%) of all ALT alleles. At
   or above `MIN_UNFILTERED_RESIDUE` (5%) the rows show no filter, and they then overrule a header line
   or a file name that says otherwise (a `bcftools concat` of the two halves keeps the first half's
   `view -v snps` line). Between the two, or with fewer alleles, the rows do not settle it.
2. **A header line** (`stated`): a bcftools `view`/`filter` command line with `-v`/`--types`,
   `-V`/`--exclude-types` or an `-i`/`-e` expression that is a bare `TYPE` comparison and no soft
   filter (`-s`), or a GATK SelectVariants line with one select type. `##source=DRAGEN_SNV` is not
   one: it names DRAGEN's small-variant caller, which writes indels under the same line.
3. **The file name** (`file_name`), only when nothing above speaks: one class word between
   separators (`.snp.`, `.snps.`, `.snv.`, `.indel.`, `.indels.`, `.sv.`, `.cnv.`). `hard-filtered`
   is quality filtering and names no class. Warnings quote the word, never the name.

A genotype table writes no REF/ALT to count, so its `variant_classes` is null; its absence is
`unknown` anyway. The other loaders that read a REF and ALT count them as the VCF loader does: the
PGP key=value table, and CGI's native `var` and `masterVarBeta`, whose unanchored alleles (`""`, a
deletion's ALT or an insertion point's REF) are indels. Their content filter rests on the rows and
the file name; neither format has a header line that records a filtering command. Counting `""`
matters there: CGI writes every indel unanchored, so a genome counted without them would read
`snv_only`.

Measured 2026-10-22 on PGP-HMS uploads: the other classes are 0% of eleven `.snp.`/`.snps.` and
`.indel.`/`.indels.` halves (five with no header line to say so, Isaac and CASAVA; three cut by GATK
SelectVariants, three by `bcftools filter -i 'TYPE=…'`), 0.16% of a `bcftools view -v snps` DRAGEN
half and 0.48% of its `-v indels` half;
indels are 20.5% of an unsplit DRAGEN genome's ALT alleles (`##source=DRAGEN_SNV`) and 25.8% of a
local DeepVariant genome's. DRAGEN's `.cnv.` files and Manta's `.sv.` file hold under 10,000 alleles,
so their filter rests on the header or the name.

## Decision 5: the facet table and provenance (charter 5, 6)

Provenance that is constant across a file lives once in a `facets` frame. Each row of `calls` and
`blocks` carries only `facet_id`, plus the provenance that actually varies per row (`source_ref` and
`ref_check`, Decision 2).

`facets` columns, all nullable unless marked:

| Column | Meaning |
|---|---|
| `facet_id` (required) | deterministic: a digest of the source bytes, the archive member, the sample and the load options. Loading the same file twice gives the same id (charter 12) |
| `source_sha256`, `source_path`, `member` | what was read. The path is a label, not identity |
| `format_id`, `family` | from `sniff` |
| `vendor`, `caller`, `caller_version` | the RM2 `source` detector, from header lines only. Null = not identifiable, never guessed. `vendor` is the party the header names as delivering the file (`##center`, `##generatedby`), falling back to the format's own vendor (`vcf_cgi`); a caller's maker is not a vendor. `caller` is the program that wrote the genotypes; tools that only touched the file afterwards (GATK SelectVariants, ApplyVQSR, bcftools view) are not evidence. Two callers in one header leave both fields null, with a warning; so do two versions of one caller for `caller_version`, compared as written. A `caller_version` written as a path or a file name (`##DeepVariant_version=/home/<name>/…`) is `<path>`, as a step's version in `method_steps` is, and so is each such version a conflict warning quotes |
| `ploidy_statements` | the header's ploidy lines verbatim, in header order (RM269 half 1): a meta line whose key names a ploidy or karyotype (DRAGEN `##referenceSexKaryotype=XY`, CGI `##source_MAX_PLOIDY`), and each ploidy option of a recorded command line, prefixed with the line it came from (`GATKCommandLine HaplotypeCaller: --sample-ploidy 2`, `bcftools_callCommand: --ploidy GRCh38`). An option whose value is a path writes `<path>` for it. No vocabulary: reading a model out of these is the karyotype work's job. `[]` = a VCF header with none; null = not a VCF. just-vcf's own shape |
| `sample` | the sample column name, verbatim; for a CGI native file its `#ASSEMBLY_ID`, which CGI's VCF uses as its sample column. **Identifying**: PGP uploads carry real names here. Logs, warnings and error messages name the `facet_id`, never the sample |
| `assembly_label`, `assembly_label_basis` | e.g. `GRCh37`; basis is `stated` (header), `inferred` or null. A stated label is kept when the lengths disagree, with a warning (RM4, below). `inferred` has three sources: a VCF's contig lengths at load (RM4), an array's own rows scored against the reference by `resolve_ref` (F9, Decision 2), and a lifted facet's target collection (`liftover`, "Lifting"). A consumer tells them apart by `relation = liftover_to` for the third, then by `header_contigs`, which only the first has |
| `seqcol_digest` | refget Sequence Collection digest (charter 3). Null until computed; the label never stands in for it. No VCF header can give it (below) |
| `sorted_name_length_pairs_digest` | the seqcol `sorted_name_length_pairs` attribute digest of the header's `##contig` names (verbatim) and lengths: the coordinate system, not the collection. Null when a contig has no length or there are none |
| `ref_seqcol_digest` | the Sequence Collection `ref` was read from (RM3, `resolve_ref`). Null until REF is read. It is kept apart from `seqcol_digest`: the reference that checks the facet is not evidence of the facet's own collection, and filling one from the other would make that check agree with itself |
| `resolve_version` | the `RESOLVE_VERSION` of the `resolve_ref` that read REF (u32). Null until REF is read, and on a facet a store migrated from before 0.0.20 whose REF was read (not recorded). A derived facet keeps its parent's; a consensus its parents' where they share one, else null. `resolve_ref` keeps the `facet_id`, so this column, `ref_seqcol_digest` and `resolve_inputs` are what tell two REF states of one facet apart (the durable store refuses an add that meets a held facet in another). just-vcf's own shape |
| `resolve_inputs` | what else `resolve_ref` read rows with: a 16-hex digest over the chip manifests whose indel probes read an array's `I`/`D` codes and whose SNV allele pairs read a mixed strand's withheld rows (`probes.catalog_digest` of each, so the manifests' bytes and `PROBES_VERSION`). Null until REF is read, where no manifest read a row (a VCF, an array with no code row and no strand withheld, or none configured), and on an array or consensus facet a store migrated from before 0.0.21 whose REF was read (not recorded). A derived facet keeps its parent's. A consensus has a digest of its parents' REF states instead (each parent's `facet_id`, `ref_seqcol_digest`, `resolve_version` and `resolve_inputs`), always set, since its `resolve_version` is null wherever they differ. With `ref_seqcol_digest` and `resolve_version` it is the facet's REF state (`store.REF_STATE`). just-vcf's own shape |
| `callset_type` | `wgs`, `wes`, `panel`, `array`, `imputed`, or null (RM4, below). A VCF gets `array` or `imputed` where its header says what the content is (Illumina GenotypingToVCF, a 23andMe conversion, Gencove, Minimac, Genes for Good: the `source` detector, RM2); otherwise RM4 infers it from the rows |
| `absence_means` | Decision 4 |
| `variant_classes` | `list[struct{variant_class, alleles}]`: the ALT alleles of each class the rows carry, counted at load, in the order `snv`, `mnv`, `indel`, `sv`, a class with none left out. `snv` and `mnv` are sequence ALTs as long as their REF (one base, more); `indel` a sequence ALT of another length (the empty allele `""` included: an unanchored deletion, F12's `-`), whatever its length; `sv` a symbolic ALT or a breakend. `<*>`, `*`, `<CGA_NOCALL>` and ALT `.` are no class, and neither is an ALT of letters outside `ACGTN`, which is no VCF allele (Veritas's `CNV`, F15; counted in a load warning). `[]` = a facet whose rows carry REF and ALT and no variant allele; null = not counted (a genotype table, which writes no REF or ALT). just-vcf's own shape |
| `content_filter` | `snv_only`, `indel_only`, `sv_only`: the one class the rows were cut down to, or null when no evidence settles one. Decision 4 says what it does to absence |
| `content_filter_basis` | `rows`, `stated` (a header line records the filtering command) or `file_name` (only the name's words say so); null with a null `content_filter` |
| `phasing_method` | the tool and version of the one `phase` step in `method_steps` (`Eagle 2.3.5+htslib-1.4`, `DRAGEN 05.021.408.3.4.11`, `WhatsHap 1.7`). Null when the header records no phasing step, or more than one (with a warning). A caller counts as phasing its own calls only where the header declares FORMAT `PS` in that caller's own words (`method_steps`, below) |
| `reference` | the reference the header says the rows are against: the base name of `##reference=` (`hg38.fa`, `human_g1k_v37.fasta`; DRAGEN's `…/<dir>/reference.bin` gives `<dir>`), written `<file>` unless every word of it is a known reference's (`source._known_reference`: an assembly word such as `hg38`, `GRCh37`, `hs37d5`, `GCA_000001405`, `chm13`, beside the bundles' descriptors and extensions), since a base name can carry a person's or a lab's name (`ALICESMITH-hg38.fa`); or, with no such line, the one reference the call steps' command lines name (`--reference`, `reference_sequence=`, `--ref-dir`, freeBayes `-f`). Two references leave it null with a warning. A `##reference=` that differs from the caller's command line is kept, with a warning unless a `lift` step is recorded: a lift the header does not name (Picard LiftoverVcf writes no line) shows up this way. Both checks compare the base names as written, before either is written `<file>`, so two different unknown files are two references; their warnings name no file. A label, never identity: `seqcol_digest` is that. just-vcf's own |
| `aligner`, `aligner_version` | the tool of the one `align` step: today only DRAGEN, where its command line says `--enable-map-align true`. Null when no header line records an alignment, which is most files: the callers that read a BAM do not copy its `@PG` lines |
| `method_steps` | `list[struct{tool, version, role, options, line}]`: every step the header records, in the order the steps ran (charter 6). `role` is one of `align`, `call`, `filter`, `normalize`, `merge`, `lift`, `annotate`, `phase`, `other` (`schema.STEP_ROLES`, an open vocabulary, so a string). `options` are the step's key options as written (`-m -both`, `selectTypeToInclude=[SNP]`, `mode=SNP`), from a list kept per tool: inputs, outputs, sample names and defaults are not kept, a path or file name is written `<path>`, a reference or chain keeps its base name where it is a known reference's, as `reference` does, else `<file>`. A `tool` or `version` written as a path (`/home/<name>/bin/CrossMap`) is `<path>` too. `line` is the meta key it came from. One call step per caller name; command lines of the same caller add their options to it. Order is header order, except that steps whose line states a run date (GATK, DRAGEN, bcftools) are sorted by it among the places they hold, since htsjdk writes meta lines sorted by text. These are the steps the header records: a facet just-vcf derived (`normalize`, `liftover`, `rejoin_split`) keeps its parent's list, and what just-vcf did is in `links.provenance` (Decision 12). Each facet of a caller split (Decision 8) carries the whole file's chain, both call steps and the merge; `caller` says which half it is. `[]` = a VCF header that records no step; null = not a VCF. just-vcf's own shape (charter 2) |
| `derived_from`, `relation` | parent facet and how: `liftover_to` (Cat-VRS, written by `liftover`), or a documented just-vcf relation (`caller_split`, `merge`, `normalized`, `split_rejoined`, `appended_export`). `appended_export`, like `caller_split`, has no parent: a second export appended to a file after the first, in another layout (F3), loaded as its own facet beside the first's; `derived_from` stays null, since neither export was made from the other, and the two share `source_sha256` and `member` (below, "Two exports in one file"). A facet with several parents (`merge`) has `derived_from` null and its parents in `links` (Decision 12), which also links a caller split's facets and an appended export to the first; Decision 12 states how the two agree |
| `contig_map` | `list[struct{source, chrom}]`: how this file's contig names were mapped (Decision 7) |
| `header_contigs` | `list[struct{name, length, md5}]` from `##contig` lines. It is evidence for the assembly (RM4) and for the seqcol digest, not the digest itself |
| `vcf_version` | from `##fileformat`. It decides how `GP` is scaled (Decision 3) |
| `sample_meta` | the file's `##SAMPLE` and `##PEDIGREE` lines for this sample, verbatim. Pedigree becomes GA4GH Pedigree later (charter) |
| `dropped_fields` | INFO/FORMAT keys present in the file but not carried (Decision 9). For a delimited table, the file's own column names that no frame column carries (a Final Report's `GenCall Score`), or, where the file names none, the format's names for them; for CGI's native files and the key=value table, their columns and keys. `[]` = the loader read the file's fields and dropped none (every 23andMe-layout table: its four fields are all carried). Null = the format has no notion of fields; no loaded format is one today. A derived facet keeps its parent's |
| `split_marker` | the declared split-marker INFO key(s), comma-joined, when any record was tagged as split off a multi-allelic site (vt's `OLD_MULTIALLELIC`, `write_vcf`'s `JUST_VCF_SPLIT`, a bcftools `--old-rec-tag` key). Null = no record was tagged, which is not "not split": see `split_signal` |
| `gl_source` | which FORMAT key a VCF facet's `gl` is read from, by the header's declarations: `gl`, `pl`, or `gl_over_pl` where both are declared (a row's `gl` is its GL, and its PL only where the GL is missing or does not fit the genotypes; which key a given row took is not kept). A key declared with a fixed count or a refused Type is not read, so it does not count; one the rows use undeclared, which the loader declares, does. CGAPipeline's phred-unit GL is `gl`: the value names the key, the scale is in a load warning. Null = the facet carries no likelihoods: not a VCF, or a header with neither key, and then `gl` is null on every row. A derived facet keeps its parent's. just-vcf's own shape |
| `derivation_inputs` | `list[struct{role, name, sha256}]`: the files a derived facet was made with besides its parent and its reference collection (`ref_seqcol_digest`). A liftover writes its chain (`role = chain`, the file's name, not its path, and the sha256 of its bytes). Null on a facet read from a file. just-vcf's own shape |
| `chip`, `chip_basis` | the genotyping chip an array export came from (RM4, `just_vcf.chip`), and how it is known: `stated` (the header names it: AncestryDNA's array version, a Final Report's `Content` manifest, Living DNA's genotype chip) or `inferred` (a 23andMe export's row count and internal-id profile both point at one generation; a preamble-less table in its layout is read so when it carries 23andMe's `i…` ids; an FTDNA or MyHeritage export, or a Final Report that names no manifest, whose rsids one Illumina manifest clearly holds, by marker overlap). An open vocabulary: `23andme_v1_v2` (the two are not told apart), `23andme_v3`, `23andme_v4`, `23andme_v5`, `ancestrydna_v1`, `ancestrydna_v2`, `living_dna_sirius`, `illumina:<manifest>` as a Final Report names it, and the platform families `illumina_gsa` and `illumina_omniexpress`, which only marker overlap names: a stated value carries the manifest, an inferred one the family, since one manifest per family is compared with, and a GSA v2 Final Report matches the v3 manifest at 0.97: the overlap does not tell versions apart. Null = unknown (a statement with no entry is quoted in a warning; an overlap that names no chip says why in one) or not an array. just-vcf's own shape |
| `strand_check` | `forward` or `mixed`: what `resolve_ref` found about the strand an array's genotypes are written on (F13, Decision 2). Null = not checked: not an array, REF not read yet, or too few heterozygous rows. Never read as forward |
| `subject` | the person the facet belongs to, as the caller labels them (RM273): supplied at `load(..., subject=...)`, never read from the file, never a sample name, never inferred. Null = not supplied. Every derived step's facet inherits its parent's. It is a claim, not identity: `facet_id` does not digest it, so one file loaded under two subjects gives one `facet_id` twice. It identifies a person if the caller makes it so: logs, warnings and errors name the `facet_id`, never the subject. Grouping by it and the check of the claim are under "Facets of one person" below. just-vcf's own shape |

No load timestamp: it would break determinism, and the file digest already says which bytes were read.

Where this table says when a column is null, `tests/test_load_seams.py` holds it per format: it loads
every format's example and fails on a null column without a stated reason. A column added here gets
its reasons there, or the loaders fill it.

**Assembly identity from a VCF header (RM4, checked 2026-10-22).**

- `seqcol_digest` cannot be computed from a header. The top-level digest is taken over the `names`
  and `sequences` attributes (refget 0.12's schema: `inherent: [names, sequences]`), and `sequences`
  holds each contig's refget sha512t24u digest, which needs the bases. A `##contig` `md5` is another
  digest of the same bases and does not convert. Filling it from the label, or from a table keyed by
  names and lengths, would be a lookup standing in for content (charter 3): the public seqcol server
  holds two different collections for UCSC hg38's names and lengths, and two for hg19's. So it stays
  null until RM3 reads bases.
- `sorted_name_length_pairs_digest` is conformant on its own: it is an attribute the seqcol spec
  defines, computed with refget's own functions (checked against the six compliance vectors refget
  ships, and against the public server seqcolapi.databio.org, where the digests computed from UCSC
  hg19's and hg38's `chrom.sizes` equal this attribute of two registered collections each). It gets
  its own column rather than `seqcol_digest`'s, because it identifies a coordinate system, which
  several collections share.
  Names are verbatim, so `chr1` and `1` headers give different digests, as seqcol intends.
- `assembly_label` is inferred from primary-chromosome lengths against GRCh36 (NCBI36, hg18), GRCh37,
  GRCh38 and T2T-CHM13 (domain constants, `just_vcf.assembly`; GRCh36's from NCBI's assembly report for
  GCF_000001405.12, carried in `assets/ncbi_assembly_reports/`): at least three must match one
  assembly and none may match another or none. No primary length is shared between any two of them.
  The mitochondrion only contradicts (GRCh37, GRCh38 and T2T-CHM13 share 16,569; hg19's 16,571 counts
  as GRCh37, and is GRCh36's only one, as hg18 and Ensembl 54 carry it: NCBI's report lists no
  mitochondrion). So a header with only chrM, or with two GRCh36 chromosomes, leaves it null, never
  the assembly with the nearer lengths. Fewer, or a mixture, leaves it null. A stated label (sniff's
  header hints) is never replaced; a disagreement is a warning.

**`callset_type` and `absence_means` for a VCF (RM4).** From the spread of the facet's calls on the
primary assembly (`just_vcf.callset`, ported from just-dna-lite `detect_callset_scope`): `wgs` needs
at least 1,000,000 sites and 75% of the span within 10 kb of a call; `wes` 10,000–400,000 sites on
all 22 autosomes with at most 60% breadth; `panel` 50–5,000 sites reaching at most 50 Mb. A region
extract of a larger callset reads as `panel`, which is what its absences then mean. `imputed` when the
header names an imputation tool or declares an imputation INFO key. Every sequencing class needs a
read-evidence FORMAT key (AD, DP, GQ, PL, GL, MIN_DP) in the header. A site is a row with a sequence
ALT (bases `ACGTN`) or none; rows whose ALTs are all symbolic do not count, so a genome-wide SV callset is not `wgs`.
Anything between the bands is null. A CGAPipeline VCF typed `wgs` gets `absence_means = reference`
(Decision 4). A gVCF is not a callset type: its blocks set
`absence_means = blocks`, whatever its scope. Otherwise `wgs` gives `weak_reference` and every other
type, null included, `unknown`. A lifted facet is `unknown` whatever its type ("Lifting"). Measured on PGP-HMS uploads: four Genos exomes at 44,020–44,796 sites
and 13.6% breadth, about 400 Mb reached.

CGI's native `var` and `masterVarBeta` are typed `wgs` where `#TITLE` says the file is a whole
genome. Where the header is silent (28 of 250 PGP-HMS files have no `#TITLE`), the same rule reads
the rows; the read-evidence gate counts as met, since the format is CGAPipeline's own output from
reads. Either way `absence_means` stays `blocks`: the files state their reference calls. Rejected:
typing `wgs` from the format alone, which a region cut from a native file would then claim too.

There is no `ref_basis` column. A non-null `ref` always came from the facet's assembly sequence
(charter 10), so the basis would have one value. Which collection it came from is `ref_seqcol_digest`.

## Decision 6: the `calls` columns

With Decisions 1–5 as recommended:

| Column | Dtype | Null means |
|---|---|---|
| `facet_id` | `str` | never null |
| `chrom` | `str` | never null (Decision 7) |
| `pos` | `u32` | never null. 1-based, VCF `POS` |
| `end_pos` | `u32` | the length of the REF span is unknown (a 23andMe `I`/`D` row, a Living DNA run-together indel before RM3, an insertion point before `resolve_ref` reads its anchor). 1-based, inclusive. Rule below |
| `ids` | `list[str]` | the source had no id. rsids, 23andMe `i` ids and VCF `ID` split on `;` |
| `ref` | `str` | the assembly sequence has not been read yet (charter 10). Only ever filled from the sequence |
| `source_ref` | `str` | the source wrote no REF (arrays). What the file wrote, as written |
| `ref_check` | `Enum` | there is no `source_ref` to compare, or it is the empty REF of an insertion point whose anchor `resolve_ref` read (no base to compare). `not_run`, `match`, `mismatch` or `unstated` (a source REF of only `N` over bases the reference reads, Decision 2) |
| `alts` | `list[str]` | the source did not say (a non-VCF row before RM3), or RM3 withheld it: an array row of a `strand_check = mixed` facet whose genotype reads differently on the two strands keeps `ref` and has `alts` and `gt` null (Decision 2). `[]` = a stated monomorphic record (ALT `.`). Symbolic alleles and `*` kept verbatim, except `<NON_REF>`, which is written `<*>` |
| `source_alleles` | `list[str]` | the source gave no vector to index (Decision 2) |
| `vrs_allele_ids` | `list[str]` | parallel to `[ref] + alts`, gnomAD's convention: the REF slot first, then one per ALT in `alts` order. An element is null when it was withheld (charter 4); the whole list is null when nothing was minted. Minting rule below ("When an id is minted") |
| `vrs_errors` | `list[str]` | parallel to `vrs_allele_ids`: what minting itself reported for an allele it tried and could not mint (`no_contig`, `refused`). Null where minting succeeded or never ran. Only facts that minting alone learns go here; see "Why an id is missing" below |
| `gt_alleles` | `list[str]` | the genotype text cannot be split into haplotypes without the site's alleles (Living DNA writes an indel's two allele sequences back to back, `CCA`; `resolve_ref` splits those where one split fits the reference, Decision 2), or is not readable at all (an array's, the key=value table's or Nirvana's, Decision 2); `source_gt` keeps the text. Entries null per uncalled haplotype |
| `gt` | `list[u8]` | there is no allele list to index yet, or RM3 withheld it with `alts` (a mixed-strand array, Decision 2), or the genotype text is not allele indices (the key=value table, Nirvana JSON; `source_gt` keeps it, Decision 2). Indices into `[source_ref] + alts` when the source wrote a REF, otherwise into `[ref] + alts` after RM3 (Decision 2) |
| `source_gt` | `str` | the source wrote no genotype text (a BCF encodes GT in binary). The file's raw text, not polars-bio's rewrite of it |
| `phased` | `list[bool]` | the source says nothing about phase for the whole row (every array export), or `gt` is withheld as genotype text that is not allele indices (the key=value table, Nirvana JSON, Decision 2). Otherwise parallel to `gt_alleles`, with an element null where the source did not say for that haplotype |
| `phase_sets` | `list[str]` | parallel to `gt_alleles`. The raw PS, PSL or CGI hapLink value. Null for an unphased haplotype, and for a phased one in its chromosome's implicit set. Keyed by kind: PS and hapLink per chromosome, PSL per file (Decision 2) |
| `phase_set_kinds` | `list[Enum]` | parallel to `phase_sets`: `ps`, `psl` or `cgi_haplink`. Null where the haplotype is unphased |
| `phase_set_ordinals` | `list[u32]` | the row has no PSL set or its PSO gave no value. Parallel to `gt_alleles`: FORMAT `PSO`, the order a haplotype's alleles are met along a derivative chromosome, on the slots in a PSL set; null elsewhere and per `.` (Decision 2) |
| `phase_set_qualities` | `list[u16]` | the row has no PSL set or its PSQ gave no value. Parallel to `gt_alleles`: FORMAT `PSQ`, the phred quality of each slot's PSL phasing (encoded like PQ), on the slots in a PSL set; null elsewhere, per `.` (unknown quality) and out of 0..65535 |
| `gl` | `list[f32]` | no likelihoods given |
| `gp` | `list[f32]` | no posteriors given |
| `dosage` | `f32` | no dosage given |
| `evidence_basis` | `Enum` | never null. Origin only: `assay`, `ref_block`, `absence` (Decision 3) |
| `gq` | `u16` | not given, or outside 0..65535 (counted in the warnings) |
| `hq` | `list[u16]` | no slot has a score. Parallel to `gt`: a phred-scaled confidence in each haplotype's call (VCF's reserved FORMAT `HQ`, haplotype qualities). From CGI's native files, `varScoreEAF` per haplotype, the lowest of a `var` haplotype's sub-rows (rule in "CGI's allele scores", below); from a VCF, FORMAT `HQ` per GT allele (rule in "A VCF's haplotype qualities", below). An element null on an uncalled (`.`, `?`) or partly called slot, where the source gave none, or outside 0..65535 (counted). `write_vcf` writes it as FORMAT `HQ` (`Number=2`; a haploid row as `51,.`; a row of more than two GT alleles left out, counted), `write_vcf_zarr` as `call_HQ` on the ploidy dimension |
| `dp` | `u32` | not given |
| `ad` | `list[u32]` | not given |
| `qual` | `f32` | not given. It is the confidence in what ALT asserts: that the variant is there on a variant record, and that no variant is there on a monomorphic record (`alts = []`, where VCF defines it as `-10·log10 P(variant)`). The same threshold therefore works on both |
| `filters` | `list[str]` | the source said `.` or has no filter column. `["PASS"]` is literal. Site-level |
| `sample_filters` | `list[str]` | FORMAT `FT`: whether *this sample's* genotype passed, separate from the site. Null = not given |
| `split_signal` | `Enum` | no sign that the row is one record of a multi-allelic site split one row per ALT. Otherwise why it looks so: `marker`, `slot_pattern` or `all_zero_beside_alt` (rule below). In such a row a `0` can mean another ALT, not REF |
| `sv_type` | `str` | the source wrote no INFO `SVTYPE`. Verbatim, an open vocabulary (`DEL`, `BND`, DRAGEN's `CNV` on a `<DEL>`, Canvas's `LOH`, BreakDancer's `ITX`). Never derived from ALT, which holds the allele type already (structural variants, below) |
| `cipos` | `list[i32]` | the source wrote no INFO `CIPOS`, or not one pair per ALT. Two offsets from `pos` per ALT, VCF 4.5's shape; an element null where the source wrote `.` (structural variants, below) |
| `ciend` | `list[i32]` | as `cipos`, from INFO `CIEND`: offsets from the end. Never filled from `cipos`, though VCF 4.5 says a missing CIEND matches it |
| `copy_number` | `f32` | the source wrote no FORMAT `CN`: this sample's copy number of the segment |
| `alt_copy_numbers` | `list[f32]` | the source wrote no INFO `CN`, or not one value per ALT. Parallel to `alts`: VCF 4.5's copy number of each ALT allele; a VCF 4.2 INFO CN (`Number=1`, the copy number of the segment) is the one-ALT list |

**CGI's allele scores (`hq`, 0.0.19).** A Complete Genomics `var` or masterVarBeta file scores each
haplotype's call: `varScoreVAF` and `varScoreEAF` per `var` sub-row, `allele<n>VarScoreVAF/EAF` per
masterVarBeta slot, phred-scaled, empty on blocks and on `?` stretches. The loader reads EAF, the
equal-allele-fraction (diploid) model's score, which cgatools takes for germline comparison
(`calldiff --diploid` uses it instead of VAF); VAF is only compared, and a file where any
haplotype whose EAF is read into `hq` has two scores that differ says how many in a load warning (a
`?`, partly called or out-of-range haplotype's EAF is not read, so its VAF is not compared). In both PGP-HMS `var` files measured
(CGA 2.5.0.33 and 2.5.0.44, chromosome 22, 2026-10-24) the two are equal on every scored row, so the
choice moves nothing there. Per haplotype:

- a `var` haplotype of several sub-rows takes the lowest of their scores, cgatools' rule for the
  masterVarBeta it builds ("for multi-call loci, this is the minimum score across all calls",
  cgatools 1.8.0 command-line reference). 52% and 56% of the multi-row haplotypes in those files carry
  differing scores, and every called sub-row has one;
- a slot written `?`, or only partly called (`?A`), gets none: there is no one sequence for the score
  to speak for. masterVarBeta writes a score on a half-call's `?` slot too, which is not read;
- a score is a whole phred number, the rule a VCF's HQ follows (below): `30.0` reads as 30, a
  number with a fraction (`30.5`) is left out, never truncated or rounded, and text that is no finite
  number (`x`, `nan`) is left out too, each reason counted per value on a called haplotype in its own
  warning. A `var` haplotype one of whose sub-rows is left out so has a sub-row without a score, and
  gets none. Both PGP-HMS `var` files write only digits there (8,249,694 and 8,154,926 EAF fields,
  the VAF fields alike, counted 2026-10-25), so nothing moves on them;
- format 1.x's one `totalScore` is neither score, so those files leave `hq` null.

CGI's own VCF of an assembly writes `GQ` as the lowest of these per locus: on chromosome 22 of a
PGP-HMS assembly (hu1B2689), vcfBeta's GQ equals the lowest `var` haplotype score at 42,894 of the
42,895 fully called loci where both write one. Rejected: writing that lowest score into `gq`, which
VCF defines for a variant site's genotype and which `genotype_evidence` does not read for a reference
genotype, where CGI scores reference haplotypes as it does any other; and keeping the scores out of
the frame, which left every CGI `var` row unquantified (`gt_only`), so a `var` facet could never be
banded in the method-error atlas.

**A VCF's haplotype qualities (`hq`, 2026-10-24).** VCF reserves FORMAT `HQ` (`Number=2`,
`Integer`): "haplotype qualities, two comma separated phred qualities", one per haplotype of a diploid
genotype. The loader reads it into `hq`, value `i` for GT allele `i`, and never fills it from GQ or any
other key:

- **The Number.** `2`, VCF 4.4's `P` (one per GT allele, read through the `Number=.` copy, P10) and `.`
  are read. Any other declaration (`1`, `A`, `R`, `G`, another count) names another quantity under the
  key: not carried, listed in `dropped_fields`, with the load warning that names keys declared so,
  and a row writing more values than that Number allows (`51,40` under `Number=1`) does not fail the
  load. An HQ declared with a Type other than Integer in a VCF 4.3+ header is not carried either
  (Decision 9). An undeclared HQ the rows use is declared `Number=2`, as for any carried key.
- **Per row.** Values are read where there is one per GT allele. VCF's fixed count makes a writer pad
  a haploid row to two: `51,.`, as CGI's vcfBeta does on every haploid row (90,251 of 90,251 with HQ
  in a PGP-HMS vcfBeta, hu1B2689, whole genome, 2026-10-24), so values past the GT's alleles are allowed
  where they are missing, and the row reads `[51]`. Any other count (a haploid `51,40`, a diploid `30`,
  a polyploid row under `Number=2` with fewer values than alleles) leaves the row's `hq` null, since
  which allele each value speaks for is unknown; so does a record with no GT. Counted in one warning.
- **Per value.** `.` is no value. A value on an uncalled allele (`./1:20,30`) speaks for no sequence
  and is not read, the rule CGI's native `?` slots follow; counted. A value outside 0..65535 is left out,
  counted. A VCF 4.2 header can declare HQ `Float` or `String` (polars-bio holds only 4.3+ to the
  reserved Integer, Decision 9): a whole number is read as it is (`30.0` is 30), a number with a
  fraction (`30.5`) is left out, never truncated or rounded, since that would be a score the file did
  not write, and text that is no finite number (`x`, `nan`) is left out too, each reason counted in
  its own warning. A row left with no value is null (`.` and `.,.` alike).
- **Blocks** hold no `hq`; an HQ on a reference-block record is not read.

CGI's vcfBeta (`vcf_cgi`) writes HQ beside `EHQ`, "haplotype quality, equal allele fraction
assumption"; the native loader reads the equal-allele-fraction score (`varScoreEAF`, above). The VCF
loader reads HQ, the key VCF reserves, in every dialect, and EHQ stays in `dropped_fields`. In the
hu1B2689 vcfBeta the two are equal on all 4,210,685 rows with HQ (2026-10-24), so the choice moves
nothing there. A row with likelihoods, or with a GQ on a genotype carrying an ALT, is graded by those
first (`genotype_evidence`, below): HQ grades a VCF row only where neither does.

**`split_signal`.** A file split one row per ALT (bcftools `norm -m-`, Hail, GATK) can write the
other ALT of a genotype as `0`, so `1/2` at `A>G,T` becomes `1/0` at `A>G` and `0/1` at `A>T`, and
`gt_alleles` reads each `0` as REF (RM270). The loader flags such rows and changes no value. Within a
contig, over the rows that become `calls` (blocks never take part), first match wins:

1. `marker`: the row carries INFO `OLD_MULTIALLELIC` (vt, Glow) or `JUST_VCF_SPLIT` (`write_vcf`'s
   own tag, any value), or a value with five `|` fields under a key declared with bcftools'
   `--old-rec-tag` Description. A record bcftools only realigned has four and is not flagged. So a
   facet `normalize` split keeps `split_signal = marker` through `write_vcf` and `load`; a row RM270
   flagged by a pattern comes back as `marker` too, since the file now carries the tag.
2. `slot_pattern`: two or more non-reference rows at one `pos`, each GT slot non-zero in at most one.
3. `all_zero_beside_alt`: a fully called all-`0` row, of any ploidy, beside a non-reference row at
   the same `pos`. Every row of the position gets the signal.

A signal is evidence, not proof: an SNV and an indel on two haplotypes at one POS (`0|1`, `1|0`) fit
the slot pattern. Out of reach: Hail's split output, which sorts unphased calls (`0/1` + `0/1`); split
indels that left-alignment moved to another POS; and a split file where no position kept two rows.
Rows written with `.` for the other ALT (vt, Glow, vcflib, `bcftools --multi-overlaps .`) are already
half-calls. Each signal is one warning line with its count. The patterns are read from the loader's
parquet spill of `calls`, three columns, so they cost no second parse of the VCF.

**Re-joining: `rejoin_split(loaded)`** returns a derived facet (`relation = split_rejoined`,
`derived_from` the parent, which is never changed). Per position, the rows that carry a
`split_signal` (only those) become one site, from the rows alone, as `bcftools norm -m+` does:

- `source_ref` is the longest REF; every other must be a prefix of it, and each ALT is extended by the
  REF bases it lacks (`A>G` beside `AT>A` becomes `AT>GT`). No FASTA is read.
- `alts` are the distinct extended alleles, **sorted by allele string**: split rows do not say the
  source's ALT order, so a consumer compares ALTs as a set, or through `gt_alleles`.
- Each GT slot takes the one allele the rows place there; `0` when no row places an ALT and a row says
  `0`; null when every row is null there. bcftools' `0` and vt's `.` read alike.
- `hq` is joined per slot too: each slot takes the one score the rows give there, null where they give
  none or two different ones, or where the joined GT slot is uncalled; null when no slot keeps one. So
  a score copied whole to every row (bcftools copies a `Number=2` field) and one kept on a row's own
  slot only (vt's `.`, `normalize`) read alike (`REJOIN_VERSION` 5; before, `hq` was kept only when
  every row's whole list agreed, which `normalize`'s rows never do at a heterozygous site).
- `ids` and `filters` are unions (sorted; `PASS` only when no row failed a filter); `phased`, the phase
  sets, `dp`, `qual`, `sv_type` and `copy_number` are kept when every row agrees and null otherwise.
- The per-ALT columns, `cipos`, `ciend` (two values per ALT) and `alt_copy_numbers` (one), are each
  row's values concatenated in the joined `alts` order; null when a row does not state them for each
  of its ALTs, or two rows giving one ALT disagree.
- Withheld: `gl`, `gp`, `ad`, `gq`, `dosage`, `source_alleles` (split rows hold REF against one ALT,
  and the full vectors would need the reads, which is calling genotypes), and `source_gt` (no file
  wrote the joined text).

A position stays as written, its rows still flagged, when it cannot be read as one site: two ALTs in
one slot, rows of different ploidy, a symbolic allele, REFs that are not prefixes of one another, an
index past a row's ALTs, a row with no genotype, or a row with no `source_ref` (one `normalize` moved).
Each reason is one warning with a count. Not done:
grouping by the bcftools marker's original record, which would also catch split indels that
left-alignment moved to another POS.

**Normalizing: `normalize(loaded, reference)`** (RM6) returns a derived facet (`relation = normalized`,
`derived_from` the parent, which is never changed; its `facet_id` digests the parent's, the
collection read and the rule's version). It takes a facet `resolve_ref` has read, and reads more bases
from the same collection (`ref_seqcol_digest`) of the same store. In order:

1. **Split.** A row with two or more ALTs becomes one row per ALT, as `bcftools norm -m-any
   --multi-overlaps .`, vt and Glow write them: a GT slot holding another ALT is null (`.`), never
   `0`, so `1/2` at `A>G,T` becomes `1/.` and `./1`, and `1/1` gives `./.` on the uncarried ALT. Each
   split row has `split_signal = marker`: in the frame it is tagged as split, as vt's output is, and
   `rejoin_split` can join it again. `source_gt` is null (no file wrote the text). `hq` keeps a slot's
   score only where its GT slot is kept (REF or this row's ALT): a score speaks for the allele in its
   slot, and the loaders give an uncalled slot none, so `1/2` with HQ `51,40` gives `[51, null]` and
   `[null, 40]`, and a row left with no score has `hq` null (`NORMALIZE_VERSION` 8). Vectors:
   - `ad` is `[REF, ALT]`, sliced from the parent's: read counts per allele do not depend on the other
     alleles;
   - `gl` is sliced to the genotypes over `[REF, ALT]` (haploid `0, k`; diploid `0/0, 0/k, k/k`) only
     where no GT slot holds another ALT. Where one does, no entry of the slice is the genotype the
     caller chose, and the slice alone would point at one it rejected, so it is withheld;
   - `gp` and `dosage` are withheld: a posterior over part of the genotypes does not sum to 1, and
     rescaling it would change the values;
   - `source_alleles` is the row's `[REF, ALT]` where any vector is kept, null otherwise.
   `<*>` beside concrete ALTs (DeepVariant writes it on every variant record of a gVCF, and
   `<NON_REF>`, loaded as `<*>`, in its VCF) names no allele, so it gets no row: the concrete ALTs are
   split as above, indexed over the parent's own allele list, so `ad` and `gl` are sliced exactly as
   for any other ALT and withheld where the slice would not hold the called genotype. A site left with
   one concrete ALT is an ordinary biallelic row, with `split_signal` null; `split_signal = marker` is
   set only where two or more concrete ALTs were split. Left as written, counted: a site whose
   genotype calls `<*>` (a haplotype carries some other allele, which no split row could hold), and a
   site holding any other symbolic allele (`<DEL>`, a breakend). bcftools' own `0/.` row for `<*>` was
   set aside: the loader would read it back as a reference block. `*` is split like any allele.
2. **Left-align and trim**, only on biallelic sequence rows whose `ref` was read from the facet's own
   collection and not contradicted: `ref_check = match`, or any row with `ref` set in a facet whose
   relation read every row's REF from that collection (`liftover_to`: REF is read from the target over
   the lifted span, and a row is kept only where it equals one of its alleles; its `ref_check` is null
   because the check it would record compared the source's REF on the source assembly), or, in a
   facet read from a file (no `derived_from`), a row with `ref` set and no `source_ref`, or an empty
   one, and `ref_check` null: `resolve_ref` is the only writer of `ref` there, and it read that REF
   with nothing from the source to compare it with (an array row, a Genos row, a Living DNA genotype
   split against the reference, F4; an insertion point's anchor, RM3). Drop a common last
   base, and when an allele runs empty, prepend the reference base before both and step `pos` back;
   repeat; then drop common first bases while both alleles keep one. At position 1 an allele keeps
   its last base. On a row this changed, `ref` is the new REF (read from the reference), `gt_alleles`
   and `source_alleles` follow the new spelling (same order, so the vectors stay indexed), and
   `source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids` and `vrs_errors` are null: they described
   the representation the source wrote, at the position it wrote it. A `ref_check = match` here would be
   bases read from the sequence compared with that same sequence. Rows whose REF was not shown to
   match pass through, counted, and so do rows whose ALT is no VCF allele. A deletion written without
   its anchor base, ALT the empty allele (an upload's `-`, F12; CGI's native deletions; Genos's `-`), is
   anchored on the base before when its REF was read (`CT>""` at 30 becomes `GCT>G` at 29) and then
   aligned like any deletion. At a contig's first base there is no base before, and VCF's rule anchors
   it on the base after its REF span (`CT>""` at 1 becomes `CTA>A` at 1), counted; only a deletion of
   the whole contig, with no base after either, stays as written, counted. An insertion point
   `resolve_ref` anchored (`source_ref = ""`, `ref_check` null, `ref` the anchor) is aligned as a row
   whose REF was read with nothing from the source to check (NORMALIZE_VERSION 9).
3. **Dedup and sort.** Rows that state the same thing are kept once, counted: equal in every column
   except the ones a move clears (`source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids`,
   `vrs_errors`, `cipos`, `ciend`), which say how the source wrote a row, not what it states. So one indel written twice
   by the source, once already left-aligned and once further right, with the same genotype and
   evidence, is one row, and the one kept is the row that did not move, with its source REF, check and
   GT text. Rows that state the same allele at one position and differ in anything else (genotype,
   phasing, the likelihoods, depths, quality, filters, ids) are all kept, with a warning (charter 8).
   Then canonical order (Decision 11).

Rows the loader flagged as one record of a pre-split site (RM270) pass through as written: run
`rejoin_split`, then `resolve_ref`, then `normalize`. Blocks pass through unchanged.

**Lifting: `liftover(loaded, chain, reference)`** (RM6) returns a `Lifted`: a derived facet
(`relation = liftover_to`, the Cat-VRS term verbatim; `derived_from` the parent, which is never
changed; its `facet_id` digests the parent's, the chain's sha256, the target collection and the rule's
version), and a `report` with one row per parent row, in the parent's order: the parent's `calls`
columns, then `lift_chrom`, `lift_pos` (1-based), `lift_strand`, `lift_swapped` and `lift_reason` (null
when lifted). The report is CrossMap's `.unmap` and bcftools' `--write-src` in one frame; it is a
result beside the facet, not a stored column. It takes a facet `resolve_ref` has read and a store
holding that collection and the target.

1. **The chain is checked against both assemblies.** Each chain header states its contig's size on
   both sides. Through `canonical_chrom` (and a scaffold's report names, as for REF above), the source sizes must equal the collection the parent's REF
   was read from and the target sizes the target collection. A primary contig (1-22, X, Y) that differs
   refuses the lift (`ReferenceMismatchError`: the chain is for another assembly); any other contig that
   differs leaves its rows unlifted. So UCSC's hg38ToHg19 chain lifts no MT row into GRCh37: its chrM is
   hg19's 16,571-base sequence, not rCRS. Without `target=`, the target is the one collection whose
   sizes match the chain's target side on every primary contig it shares and on at least one.
2. **Positions.** A row's REF span (`pos` to `pos + len(ref) - 1`) must lie inside one ungapped block,
   and inside no block with another target (blocks that send the same bases to the same contig,
   strand and offset count once: Ensembl's chains list them twice, F14). The span's start then moves by the block's offset; on a
   minus-strand block the forward start is the image of the span's last base.
3. **Alleles.** On a minus-strand block every allele (`alts`, `gt_alleles`, `source_alleles`) is
   reverse-complemented, in place, so the vectors stay indexed. `ref` is read from the target over the
   lifted span. Where it equals the parent's REF (carried across, complemented on a minus block), the
   alleles stay; where it equals an ALT, that ALT and the old REF swap places in the allele list and the
   genotype indices follow (`0/1` becomes `1/0`); where it equals none, the row is not lifted. Array rows
   (no `source_ref`) take their ALTs from the genotype against the new REF, sorted, as `resolve_ref`
   derives them.

   **A minus-strand indel is re-anchored** when its alleles all start with one base (VCF's padding
   base): reverse-complemented, that base ends every allele instead, so it is dropped from each and the
   target's base before the lifted span is put in front, in all four allele lists alike. The row moves
   one base left (`pos` is the lifted span's start minus one; the span keeps its length) and REF is
   read over the moved span, its first base being the new padding base; the rest of it is compared
   with the alleles as above. So `A>AC` at source 100, on a block reversing the contig, becomes
   `B>BG` at the image of source 101, `B` being the target's base there.
   The new padding base is the target's, wherever it lies, even past the chain block: REF is read from
   the target, and the event is placed by the bases the block does cover. A minus-strand indel whose
   alleles share no first base (`A>T,AG`, an SNV beside an insertion; `AT>G`) needs no padding base,
   and lifts complemented in place like an MNV. A re-anchored indel is not left-aligned on the target
   by `liftover`: the shift that was leftmost on the source is rightmost on the target. Run
   `normalize` on the lifted facet to left-align it there ("Normalizing", step 2).
4. **What a lifted row carries.** `source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids` and
   `vrs_errors` are null, and on a sequence row `cipos` and `ciend` too: they described the source's
   representation on the source assembly (the intervals are offsets a minus-strand block would
   reverse). A symbolic SV keeps them, mirrored on a minus strand ("Lifting symbolic SVs" below). `gt`
   indexes `[ref] + alts`. `gl`, `gp`, `ad` and the order of `source_alleles` are the parent's (Decision
   2: never reordered). A swapped row's `dosage` is null, since it counted the other allele, and so
   are its `alt_copy_numbers`: the swapped ALT's slot holds the old REF, whose copy number the source
   did not state. A row
   lifted to another contig than its source's has `phased`, `phase_sets` and `phase_set_kinds` null: a
   PS set is keyed per contig, and two sets would otherwise merge. `ids`, `filters`, `qual`, `gq`, `dp`
   and `evidence_basis` are kept. Re-mint ids on the target with `mint_vrs_ids`.
5. **The derived facet.** `ref_seqcol_digest` is the target collection; `assembly_label` is the
   target's label, basis `inferred` (from the collection's lengths, which the chain's target side
   matched); `seqcol_digest`, `sorted_name_length_pairs_digest` and `header_contigs` are null (they were
   the parent's header); `reference` is null too, since it named the source assembly's file (the
   target is `ref_seqcol_digest`), while `aligner`, `aligner_version` and `method_steps` stay the
   parent's: the reads were aligned to the source, and what just-vcf did is in `links.provenance`;
   `contig_map` maps the target collection's names of the contigs the rows use;
   `derivation_inputs` names the chain; `absence_means` is `unknown`, because the rows the lift dropped
   and the regions no chain covers are absent from it without being reference (the XYY stress case
   below).

Not lifted, each counted in the warnings and named in the report's `lift_reason` (an open vocabulary,
`UNLIFTED_REASONS`): `ref_unread`, `source_ref_mismatch` (`ref_check = mismatch`; an `unstated` row
lifts as a `match` one, Decision 2), `source_ref_unstated` (`ref_check = unstated` with a sequence ALT,
spelled against the `N`, `LIFT_VERSION` 13), `allele_not_sequence` (the empty allele of
an unanchored deletion, which `normalize` anchors first; `*`, `<*>` and symbolic alleles the SV rule
below does not read, such as Delly's `<TRA>` or a symbolic SV beside a sequence ALT; an SV or segment
whose REF is longer than its padding base), `split_flagged` (run
`rejoin_split` first), `contig_not_in_chain`, `chain_contig_mismatch`, `no_chain_block`,
`split_across_blocks`, `multiple_targets`, `target_contig_unmatched`, `indel_minus_strand` (a
minus-strand indel or SV lifted to the target contig's first base, with no base before it to
re-anchor on), `ref_mismatch_after_lift`, for SVs `sv_end_unlifted`, `sv_ends_apart`,
`sv_ends_reversed` and `sv_length_changed`, and for breakends `breakend_not_lifted` (an ALT with
brackets that is none of the four forms, a breakend beside another kind of allele, or a REF longer
than one base), `breakend_mate_not_lifted`, `breakend_ends_apart`, `breakend_ends_reversed` and
`breakend_distance_changed`.

**Lifting symbolic SVs** (2026-10-26, `LIFT_VERSION` 9). A row whose every ALT is a symbolic SV the
loader gives an extent (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<CNn>`, subtypes included) or `<INS>`, and
a copy-number segment (no ALT, `end_pos` past REF: Canvas's copy-neutral records), is lifted by its two
ends rather than its REF span, with one base of REF:

1. POS and `end_pos` are lifted as two single positions, each by step 2's rule (one block, one
   target). POS's failure is the row's reason as for any row; `end_pos`'s is `sv_end_unlifted`.
2. Both must land on one target contig and one strand (`sv_ends_apart`), in their order
   (`sv_ends_reversed`), and the lifted length `end_pos - pos` may differ from the source's by at most
   `SV_LENGTH_TOLERANCE` (5%) of it (`sv_length_changed`). This is liftoverSV's rule and default
   (`-p/--percent 0.05`: dropped on positions on different chromosomes, reversed order between the
   lifted positions, or a length change above 5% of SVLEN). Rejected: Picard LiftoverVcf's, which
   lifts the whole POS..END interval with `LIFTOVER_MIN_MATCH` 1.0, so one base of chain gap inside a
   200 kb deletion drops it; CrossMap's, which rewrites END from the lifted REF span and never lifts it
   (an SV comes out one base long); and bcftools +liftover's `--lift-end`, which lifts END under the
   same contig, strand and order checks but on failure keeps the record with END missing, an SV with
   no extent. A point (`<INS>`, length 0) stays a point.
3. On a plus-strand block the row lands at `image(POS)`..`image(END)`. On a minus-strand block the
   event's image runs the other way. A padded SV (a symbolic ALT, whose event is `[POS+1, END]`) lands
   with POS at `image(END) - 1`, the target's base before the event's image, and `end_pos` at
   `image(POS) - 1`: the length is kept, and the anchor is the target's base, as for a re-anchored
   indel (step 3); one that would be anchored before the contig's first base is `indel_minus_strand`.
   A segment (`[POS, END]` inclusive, as a block) lands at `image(END)`..`image(POS)`.
4. REF, one base, is read from the target at the new POS (charter 10). On a plus-strand block it is
   the source's anchor carried over, so it must equal the source's REF (`ref_mismatch_after_lift`
   otherwise, as liftoverSV and Picard drop a changed REF); on a minus-strand block the anchor is
   another base and nothing is compared, as a re-anchored indel's old padding base is not. The ALTs
   stay as written (Picard and bcftools leave a symbolic allele as it is on either strand);
   `gt_alleles` and `source_alleles` take the new REF where they held the old; `gt`, `sv_type`,
   `copy_number`, `alt_copy_numbers` and `dosage` carry. Nothing swaps.
5. `cipos` and `ciend` carry on a plus-strand block. On a minus-strand block the new POS is the old
   END's image and the new END the old POS's, so each pair `[lo, hi]` mirrors to `[-hi, -lo]` and the
   two swap; where the source wrote no CIEND, VCF 4.5's rule (assume CIPOS) supplies the old END's
   interval, so both are set. Where it wrote CIEND and no CIPOS, its POS is precise and its END is
   not: the new END, the old POS's image, gets `ciend = [0, 0]` per ALT, since a null there would read
   as the new CIPOS by the same rule (`LIFT_VERSION` 10). Neither written stays neither.
6. Breakends (`N[2:123[`, `]2:123]N`, `.N`) are not this rule's: they lift by their own position and
   their mate's ("Lifting breakends" below).

Tested on the synthetic chain against the construction (each reason, both strands, a segment, the
5% boundary, `cipos`/`ciend` on both strands, a write and reload), and opt-in against pyliftover's
placing of both ends of 550 `<DEL>`s of 50 bases to 200 kb over the GRCh37 autosomes through
`hg19ToHg38` (`tests/test_liftover_chain.py`).

**Lifting breakends** (2026-10-26, `LIFT_VERSION` 11). A row whose every ALT is a breakend, one of
VCF's four forms (`t[p[`, `t]p]`, `]p]t`, `[p[t`, brackets paired, bases `ACGTN`) or a single
breakend (`t.`, `.t`), with a REF of one base, is lifted by its own position and each ALT's mate
(`just_vcf.breakends` reads the ALT; `concordance` reads it there too):

1. POS lifts as one base, by step 2's rule. The target's base there is the new REF, and must equal
   the source's REF, complemented on a minus-strand block (`ref_mismatch_after_lift`): the base at a
   breakend's POS is part of the joined sequence, not a padding base before an event as a symbolic
   SV's is, so it is compared on both strands, as an SNV's is.
2. Each ALT's mate position `p` lifts as one base through the chain blocks of the contig the ALT names
   (read with `canonical_chrom`), whichever contig that is, except a contig of the `##assembly` file
   (step 7). A mate with no single target (in no block,
   in several, on a contig the chain does not cover or whose size it states differently, or landing
   on a contig the target collection lacks) leaves the row behind as `breakend_mate_not_lifted`; a
   warning counts the mates by their own position's reason.
3. A mate on the row's own source contig is held to liftoverSV's three checks on its two ends: they
   land on one target contig (`breakend_ends_apart`); where both land on one strand, their order is
   kept on a plus strand and reversed on a minus one (`breakend_ends_reversed`); and their distance
   differs from the source's by at most `SV_LENGTH_TOLERANCE` (5%) of it (`breakend_distance_changed`).
   Ends on two strands have no order to keep. A mate on another source contig is checked by nothing
   more (liftoverSV's translocation).
4. The ALT is rewritten: the mate's contig as the target collection names it (`Contig.name`, so
   `chr2` or `2` as the target FASTA did), its position the image of `p`, and the form turned by both
   ends' strands (the table below). `t` is reverse-complemented when POS lands on a minus-strand
   block, so REF's base, now its complement, moves to the other end of `t` and any inserted bases are
   reverse-complemented on the other side of it; a flip of the mate alone leaves `t` as written.
5. `cipos` is about POS and mirrors (`[lo, hi]` to `[-hi, -lo]`) when POS lands on a minus-strand
   block; `ciend`, which `concordance` reads as the mate's interval (truvari's reading), mirrors per
   ALT when that ALT's mate does. Unlike a symbolic SV's, they never swap: POS stays POS. A breakend
   VCF placed at the left-most position of its interval is not moved back there on a minus strand.
6. `end_pos` is the new POS; `gt_alleles` and `source_alleles` take the new REF and ALTs where they
   held the old; `gt`, `sv_type`, `dosage` and the rest carry, the phase cleared on another contig as
   for any row. The lifted facet's `contig_map` lists every contig an ALT names, rows or none.
7. **A mate in angle brackets** (`C[<ctg1>:7[`, `]<ctg1>:214]A`; 2026-10-26, `LIFT_VERSION` 12, the
   lead's call) is not lifted, and does not hold the row back. VCF 4.3 §5.4.2 gives a large insertion
   as a join into "a contig from the assembly file", and §1.6.1 says an angle-bracketed ID points "to
   a contig in the assembly file (cf. the ##assembly line in the header)", a file of its own, not the
   reference. No chain covers it, and it is the same sequence on any reference, so step 2 places
   nothing: the contig, the position and the mate's side stay as written, and only POS lifts. The
   form then turns by POS's strand alone (the strand table below with the mate on a plus strand:
   `t[p[` on a minus POS becomes `[p[t`, `]p]t` becomes `t]p]`), `t` reverse-complemented as
   in step 4, `ciend` unmirrored. This is consistent with §5.4's four sentences: read through them, the
   lifted ALT on the target and the source ALT on the source describe one molecule, in one direction
   or the other, given the contig's sequence on both sides (checked per form and strand, and with
   inserted bases, in `tests/test_liftover_bnd.py`); a lift and its inverse give the row back. So the
   fallback the lead allowed, leaving such a row behind as `breakend_mate_in_assembly`, is not needed
   and no such reason exists. The contig is not in the lifted facet's `contig_map`, and the writers
   never declare it. A record whose own CHROM is such a contig (§5.4.2's circular insertion) is not
   on the reference either and stays behind as `contig_not_in_chain`.

Written back out, a breakend's mate is a contig of the file like CHROM (F30): `write_vcf` reads the
mate's contig canonically and writes it under the `contig_names` style's one canonical → written
mapping, the one CHROM is written under, and declares every contig a written ALT names in a
`##contig` line, rows or none (VCF 4.3 §1.4.7). `stored` writes the mate canonical (a lifted `chr2`
becomes `2`), `chr` writes `chr2`, and `source` writes it as `contig_map` spells it, so a lifted
facet's ALTs come back as the lift wrote them. Under `stored` and `chr`, a mate the style has no other
name for keeps its canonical name, as such a CHROM does. Under `source` (2026-10-26, the lead's call)
a mate is written as the source spelled that contig: as the facet spells it (`contig_map`, then
`header_contigs`) where it knows the contig, even where the ALT spelled it otherwise (CHROM `2` beside
an ALT's `chr2` writes `2`, so the file agrees with its own CHROM column), and as the ALT spelled it
where the facet does not (a chr-named file with no `##contig` lines writes `A[chr2:500[` and declares
`chr2`, where it wrote `A[2:500[`). A contig only ALTs name, spelled two ways by them (`chr3`, `3`), is
refused under `source`, as a contig its CHROM or header spells two ways is: which spelling the source
meant cannot be told. Single breakends and symbolic alleles are written as held.
So is a mate in angle brackets (`]<ctg1>:7]T`, VCF 4.3 §5.4.2's large insertion): `<ctg1>` is a
contig of the `##assembly` file, not of the reference (§1.6.1), so it gets no `##contig` line under
any style (2026-10-26; htslib warned "Invalid contig name" on the line it got before). Only the
written rows' mates are declared: a row `resolved_only` leaves out declares none.

Mate records stay consistent without pairing them: each lifts the same two points by the same
function, so the record at `p` lands where its mate's ALT now points, whatever the strands (checked on
676 lifted records of synthetic translocations through `hg19ToHg38`, every one). Manta's mates stated
a few bases apart (F28) are two points, and could straddle a chain edge; each is lifted as it says.
The frame keeps no MATEID, and none is needed.

*The strand table.* VCF 4.3 §5.4: `t[p[` "piece extending to the right of p is joined after t";
`t]p]` "reverse comp piece extending left of p is joined after t"; `]p]t` "piece extending to the
left of p is joined before t"; `[p[t` "reverse comp piece extending right of p is joined before t".
So each form states two sides: the side of POS whose sequence takes part in the join (left of POS
when `t` comes first, right when it comes last) and the side of `p` (right of `p` with `[`, left with
`]`): `t[p[` (L, R), `t]p]` (L, L), `]p]t` (R, L), `[p[t` (R, R). A minus-strand block reverses the
sequence around an end, so what lay left of it lies right of its image: the record's side flips when
POS lands on one, the mate's when `p` does, and the lifted form is the one stating the new pair.

| Source form | POS +, mate + | POS +, mate − | POS −, mate + | POS −, mate − |
|---|---|---|---|---|
| `t[p[` | `t[p[` | `t]p]` | `[p[t` | `]p]t` |
| `t]p]` | `t]p]` | `t[p[` | `]p]t` | `[p[t` |
| `]p]t` | `]p]t` | `[p[t` | `t]p]` | `t[p[` |
| `[p[t` | `[p[t` | `]p]t` | `t[p[` | `t]p]` |

A flip at POS alone moves `t` across the brackets and keeps their direction; a flip at the mate alone
turns the brackets and keeps `t` where it is; both give the full reversal (`t[p[` and `]p]t`, `t]p]`
and `[p[t`), which is how a mate pair's two records already relate. A single breakend's dot moves
with POS's side: `t.` becomes `.t'` on a minus block, `t'` being `t` reverse-complemented. The table
is checked against the sentences themselves: the tests read the novel sequence each ALT asserts from
the source and the target, and require the two to be one molecule, read one way or the other.

*Against the tools that lift VCF* (read 2026-10-26, none run):

- **liftoverSV** 0.3.1_beta (GPL-3, commit e62dcbf; `share/python3/liftoverSV/core/liftover_engine.py`,
  `LiftoverEngine._lift_square_bracketed_alt_notation`) is the one that lifts a breakend's mate: it
  lifts `p` through the chain, rewrites `contig:pos`, takes the anchor base from the target, and for
  a mate on the record's own contig (its `normalize_svtype` gives no `TRA`) drops the record when the
  ends land on two contigs, in reversed order, or at a distance changed by more than 5%. Followed
  here: the mate's lift, the rewrite, and the three checks. Not followed, each where it contradicts
  the spec or this rule: its `ChainLifter.lift` (`io_tools/chain_lifter.py`) keeps pyliftover's
  contig and position and drops the strand, so a breakend on a minus-strand block keeps its brackets
  and `t`, which §5.4 contradicts (the table above); its order check ignores strands, so a pair lifted
  correctly through one inverted block, whose order must reverse, would be dropped (here order is
  checked only where both ends share a strand, reversed expected on minus); and it takes the target's
  anchor base without comparing, which this lift compares, as it does for a symbolic SV on a plus
  strand.
- **bcftools +liftover** (MIT, `freeseek/score` `liftover.c`, `LIFTOVER_VERSION` 2026-02-06,
  `process`): `bcf_is_symbolic` tests only `<` and `*`, so a breakend is read as a sequence allele and,
  on a minus-strand block, `reverse_complement` turns the whole ALT string: its IUPAC table (`rev_nt`)
  maps `[` to `]`, so the brackets turn as the both-minus column does, but the mate's `contig:pos` is
  reversed with them (`A[2:100[` becomes `]001:2]T`) and never lifted.
- **Picard LiftoverVcf** (`picard/util/LiftoverUtils.java`, `reverseComplement(Allele)`) returns a
  symbolic allele unchanged, and htsjdk counts breakends as symbolic (`Allele.wouldBeSymbolicAllele`:
  `wouldBeBreakpoint`, `wouldBeSingleBreakend`): the ALT is kept verbatim on either strand, the mate
  unlifted.
- **CrossMap** (`cmmodule/mapvcf.py`, `crossmap_vcf_file`, commit 2fefee4): `is_dna` fails on a
  bracket, so the ALT is appended verbatim, the mate unlifted.

Tested on a synthetic chain with plus- and minus-strand blocks (`tests/test_liftover_bnd.py`): each
form under each pair of strands, a mate on another contig, a mate in a gap and on a contig the chain
does not cover, each same-contig check, single breakends on both strands, inserted bases, two ALTs,
the intervals, a write and reload, and a lift through the inverse chain giving back every record.
Opt-in (`tests/test_liftover_chain.py`): 800 synthetic translocation records over the GRCh37
autosomes through `hg19ToHg38`, each end where pyliftover puts it, the form by the table, and
`hg38ToHg19` giving all 676 lifted records back.

Measured 2026-10-26 on the PGP-HMS Manta files (aggregates only). hu1C1368 (GRCh37 to 38): 1,216 of
1,228 breakend records lift (4 with POS and 4 mates on minus-strand blocks); behind, 10
`breakend_distance_changed`, 1 `breakend_mate_not_lifted`, 1 `no_chain_block`. Lifted back through
`hg38ToHg19`, 1,215 return and 1,213 are their source record exactly; the 2 that are not have one end
that pyliftover also sends elsewhere on the way back (the two chains are not each other's inverse
there). `concordance` of the returned facet against the source: 607 adjacencies `match`, 7 `only_b`
(both records left behind), 1 `only_a`. hu416394 (GRCh38 to 37): 1,043 of 3,350 lift (7 POS and 7
mates on minus-strand blocks); behind, 1,310 `ref_unread` (1,147 on decoys, 163 on F25's scaffolds),
748 `breakend_mate_not_lifted` (691 mates on decoys, which UCSC's chain does not cover, 20 on F25's
scaffolds, 37 in chain gaps), 233 `no_chain_block`, 12 `breakend_distance_changed`, 4
`breakend_ends_reversed`. Back through `hg19ToHg38`, 1,038 return, 1,030 exactly, the 8 others again
ends pyliftover sends elsewhere; `concordance`: 518 `match`, 4 `only_a`, 1,157 `only_b`. About 30 s per
file for load, `resolve_ref` and both lifts.

**Lifting blocks** (2026-10-26, `LIFT_VERSION` 8). A `blocks` row's span (`pos` to `end_pos`) is cut
on the axis step 2 judges a row's span on: the chain's blocks, cut wherever the set of targets
changes, so at every gap on either side, every block edge and every change of offset. Each piece
inside one chain block with one target lifts by that block's offset, as a row's start does (on a
minus-strand block the piece's forward start is the image of its last base); a block has no alleles,
so nothing is complemented. A block that maps whole to one interval of the same length lifts as one
block; one that crosses a gap or an edge becomes one block per piece, never one block stretched over
the gap, and never merged with its neighbour even where the two pieces touch on the target (a source
deletion leaves two adjacent target intervals whose bases were not adjacent). A piece on another
contig or a minus-strand block lifts like a row (rows there are lifted too). `kind`, `gq`, `dp_min`
and `gl` carry unchanged onto every piece: they are the block's own evidence, and `gl` is over
`[REF, <*>]`, which a lift does not reorder. The frame keeps no block REF, so there is nothing to
re-read or carry: `write_vcf(..., blocks=True)` reads each block's first base from the facet's
`ref_seqcol_digest`, which on a lifted facet is the target (charter 10). A piece in no chain block, in
more than one target, on a target contig the collection lacks, or on a source contig the chain does
not cover or covers at another size stays behind with the row reason of the same name, and is counted
in bases: `Lifted.block_reasons` (bases per reason) and a warning line per reason with the pieces and
bases. Lifted pieces are sorted into canonical order per target contig. Two source regions the chain
sends to one place give overlapping lifted blocks; both are kept and counted (charter 8), and every
reader of blocks finds both where they cover a position (`reference_support`, `genotype_evidence_at`,
`concordance`, `consensus`; one answers where one must, by `reference_support`'s order). The lifted
facet's `absence_means` stays `unknown`: a position no lifted block or row covers says nothing, which
is what `blocks` would say too (`absence_is_reference` answers `None` for both), and `unknown` keeps
one rule for every lifted facet. Checked against pyliftover base by base on synthetic blocks over the
GRCh37 autosomes through UCSC's `hg19ToHg38` (opt-in, `tests/test_liftover_chain.py`).

**`end_pos`.** It is named apart from polars-bio's `end` because that value is wrong for symbolic SV
alleles (P3), so the two must not share a name. The rule follows VCF 4.5's computed END:

1. Take `abs(SVLEN)` for each allele (the spec says to read a negative SVLEN as positive). An SVLEN
   declared `Number=1` is read as the one-ALT list it is; one that is not one value per ALT cannot be
   matched to its alleles and is not read.
2. For each allele compute an end:
   - `<DEL>`, `<DUP>`, `<INV>`, `<CNV>` and their subtypes, and the copy-number alleles `<CN0>`,
     `<CN2>`… (Canvas, the 1000 Genomes SV release; RM282): the record's INFO `END` where it gives one,
     else `pos + abs(SVLEN)`. POS is the base before the event, so `<DEL>` at 200 with `SVLEN=30` ends
     at 230. END is the record's stated end (VCF 4.5 deprecates it but defines it as this very
     maximum, and a missing SVLEN is to be inferred from it); where the two disagree, END is followed
     and the records are counted in a warning. SVLEN first, as FORMAT `LEN` comes before END for a
     block, was set aside: END is a coordinate the writer stated, while `pos + abs(SVLEN)` depends on
     POS being the base before the event, which a caller counting POS as the first base breaks;
   - `<INS>` and breakend alleles: the anchor, `pos`, whatever END says;
   - `<*>`: `pos + LEN - 1`, because a block includes POS, else END;
   - a sequence allele: `pos + len(source_ref) - 1`.
3. `end_pos` is the largest of these over the site.

END is read only for those SV types, `<*>`, `<CGA_NOCALL>` and ALT `.` blocks. Another symbolic
allele's END is not known to lie on this contig (Delly's `<TRA>` writes the partner's position on CHR2
there), so it keeps its REF span. A record whose END is not read and differs from the REF span is
counted in a warning.

Where REF is not known yet (an array `I`/`D` row) `end_pos` is null. A single-position array call has
`end_pos = pos`.

**Structural variants (RM282).** An SV or CNV VCF (Manta, DELLY, Canvas, DRAGEN's `.sv`/`.cnv`,
Sniffles) loads as any VCF does: its symbolic ALTs (`<DEL>`, `<DUP:TANDEM>`, `<INV>`, `<INS>`,
`<CNV>`, `<CN0>`…) and breakends (`A]2:7000]`, `.A`) are kept verbatim in `alts`, each counted as an
`sv` allele in `variant_classes`, its extent is `end_pos` (above), and what the caller wrote about the
event goes to five columns, each just-vcf's own shape (VRS has no Allele for these; charter 2):

- `sv_type` is INFO `SVTYPE` as written. VCF 4.5 deprecates it as redundant with ALT, but callers
  write a type ALT does not say: DRAGEN's copy-number caller writes `SVTYPE=CNV` on `<DEL>` and `<DUP>`
  (an abundance claim, VCF 4.5's `SVCLAIM=D`, not an adjacency), Canvas `LOH` on `<CN0>,<CN2>`.
  Rejected: deriving it from ALT where SVTYPE is absent, which would mix stated and derived values in
  one column.
- `cipos` and `ciend` hold two offsets per ALT. VCF 4.2 declared one pair for the record
  (`Number=2`), and every ALT of a record shares its POS, so that pair is repeated for each ALT; a
  record with no ALT (a copy-neutral segment) keeps its one pair; any other count is withheld and
  counted. Rejected: storing the values as written, which leaves every consumer to read the facet's
  `vcf_version` to know what a list of two means on a two-ALT record. A missing `ciend` stays null:
  VCF 4.5 says to assume CIPOS, and a consumer applies that rule, the frame does not fill it.
- `copy_number` is FORMAT `CN`, the sample's total; `alt_copy_numbers` is INFO `CN`, per ALT under
  VCF 4.5 (its own example: `<CNV>,<CNV>` with `CN=1,2` and FORMAT `CN` 3). The two are kept apart:
  neither is filled from the other, nor from a `<CNn>` allele's name.
- A VCF 4.3 header declaring `CN` with a Type other than the Integer 4.3 reserves (it is Float from
  4.4 on) is refused by polars-bio's reader, so that key is not carried, with a warning and the key in
  `dropped_fields`; an undeclared `CN` the rows use is declared as the file's version reserves it.
  `write_vcf` writes VCF 4.3, so it declares `CN` Integer and leaves out, counted, a copy number that
  is not a whole number.
- The frame keeps no SVLEN: it goes into `end_pos` and nowhere else, so an `<INS>`'s inserted length
  and a source SVLEN that disagreed with END are not kept. `write_vcf` still writes INFO `SVLEN`
  beside END (2026-10-26), since VCF 4.5 requires it on a symbolic SV allele: `end_pos - pos` on each
  length-type allele (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<CNn>`), `.` on the record's other ALTs. It
  is END restated, never a value read from the source. Under VCF 4.3 it is `Number=.` and negative for
  a deletion, as 4.3's examples write it; under 4.5 `Number=A` and unsigned. A record with several
  length-type ALTs gets it only where VCF 4.5 says they share one length (copy-number alleles alone,
  or `<DEL>`/`<DUP>`/`<CNV>` beside FORMAT CN); otherwise its ALTs' own lengths are not in the frame,
  and SVLEN is left out and counted rather than written as one ALT's length for another. Rejected:
  a `svlen` column here, which would add an `<INS>`'s length and each ALT's own length on a record of
  several, but is a schema change for the person to decide; and SVLEN from `end_pos` on every ALT, which writes a length a `<DEL>`
  beside a longer `<DUP>` does not have.

Each of these is null on a row that is no SV. `normalize` leaves symbolic sites as written; on a
sequence site it splits, each row keeps its own ALT's pair and copy number, and a row it moves loses
`cipos` and `ciend` (offsets from the POS the source wrote), as a lifted sequence row does. `liftover`
lifts a symbolic site by both ends and keeps its intervals ("Lifting symbolic SVs"). `concordance`
compares SV calls across facets by type and extent, and breakends by both ends and the join, as
truvari does ("Concordance between two facets", RM282).

Measured 2026-10-23 on the PGP-HMS SV and CNV files (aggregates): a Manta `.sv.` file, 9,373 calls
(4,366 sequence-written indels, 3,350 breakends, 848 `<DEL>`, 558 `<DUP:TANDEM>`, 251 `<INS>`); a
Canvas `.cnv.` file, 182 records, 92 copy-number calls and 90 copy-neutral segments; a DRAGEN
`.cnv.` file, 586 calls (450 `<DEL>`, 136 `<DUP>`, every one with `SVTYPE=CNV` and a FORMAT `CN`);
two Veritas files (F7's `##CHROM` shape, which loads since 2026-10-22), 5,614 SV calls and 2,632 CNV
rows whose ALT is the bare word `CNV` (F15: no class, no extent, `end_pos = pos`).

**When an id is minted (RM271, `mint_vrs_ids`).** An id is the VRS 2.x Allele id of the fully
justified allele, on the refget sequence (`SQ.…`) of the contig in the collection REF was read from
(`facets.ref_seqcol_digest`). Every slot of `[ref] + alts` gets one: the REF slot is the reference
allele, which VRS writes as a `ReferenceLengthExpression` over REF's span (length and
`repeatSubunitLength` both `len(ref)`), as vrs-python does for REF = ALT. An ALT is normalized as the
spec says: trimmed of shared bases; a substitution keeps a `LiteralSequenceExpression` on the trimmed
span; an insertion or deletion is rolled over its whole region of ambiguity and becomes a
`ReferenceLengthExpression` when it is a deletion or an insertion that continues a repeat, and a
`LiteralSequenceExpression` otherwise. Inside VRS objects coordinates are 0-based interbase; `pos`
stays 1-based, and `reference._span` is the one conversion (charter 13).

An allele is minted only when its row's `ref` was read from the sequence (`resolve_ref`) and the
source did not contradict it (`ref_check` is `match`, or null on an array row). On a row whose source
REF states no base (`unstated`) only the REF slot is minted; its sequence ALTs are `ref_unstated`.
Substitutions too: `mint_vrs_ids` refuses a facet whose REF was never read. The rule the draft had, minting a
substitution from a table of contig accessions per assembly without reading bases, was not built:
no such table exists, and a contig `resolve_ref` could not match has no accession that was shown to
be the facet's own. Array rows mint once `resolve_ref` has derived their `alts`. Minting has been
compared with vrs-python on the same store for every kind of allele (tests) and on real PGP-HMS VCFs
(docs/ROADMAP.md, RM271).

gtars 0.10 rolls the allele (`normalize_allele`) and computes sha512t24u, but its Python surface
digests a `LiteralSequenceExpression` only, and its normalization never returns a
`ReferenceLengthExpression`. So the choice of state and the canonical JSON (keys sorted, no
whitespace, inherent properties only: an RLE's `sequence` is not digested) are just-vcf's, in one
place (`just_vcf.vrs`), and tested against vrs-python.

**Why an id is missing** is not a stored column. Most of the reasons are already facts in the frame,
and a stored copy of them drifts: RM3 fills `ref`, and a stored reason still says "unchecked". So the
reason is computed by one function, `vrs_withheld(calls)`, which returns the lazy `calls` frame with a
`vrs_withheld` column: one reason per slot of `[ref] + alts`, parallel to `vrs_allele_ids`, null where
an id exists. A row with neither `ref` nor `alts` (an array row whose REF was not read) has one slot.
It checks in this fixed order and the first hit wins:

1. **What the allele is**, read from the allele string: `unspecified` (`<*>`), `spanning_deletion`
   (`*`), `sv_no_vrs_model` (an allele `variant_classes` counts as `sv`: `<DEL>`, `<CNV:TR>`, `<CN0>`,
   breakends, `.`-anchored single breakends; RM282), `symbolic` (a symbolic allele that names no
   variant: Complete Genomics' `<CGA_NOCALL>` where it stays a call), `partial_call`
   (CGI `?`), `array_code` (23andMe `I`/`D` where the source wrote no REF), `not_sequence` (anything
   else outside `ACGTN`, such as an ambiguity code in ALT).
2. **Pipeline state**, read from the row: `array_code` (an array row with an `I`/`D` genotype, whose
   REF stays unread), `telomere` (`pos` 0), `ref_unresolved` (`ref` and `source_ref` both null),
   `ref_unchecked` (a source REF, but `ref` not read: a contig the reference could not match, or a
   span past the contig's end), `ref_mismatch` (`ref_check = mismatch`); then, per ALT slot,
   `ref_unstated` (a sequence ALT on a row whose `ref_check` is `unstated`: spelled against an `N`).
3. **What minting reported**, the `vrs_errors` element: `no_contig` (the collection has no single
   sequence for the contig) or `refused` (gtars refused the allele; its message is counted in the
   warnings). After `resolve_ref` neither is expected: `ref` was read from that very sequence.
4. Otherwise `not_run`: nothing forbids minting and it has not happened yet.

`mint_vrs_ids` counts the withheld slots per reason in the facet's warnings, one line per reason.

The reasons form an open vocabulary (a `frozenset[str]` with a validator, stored as `str`), not a
`pl.Enum`: every new kind of facet brings reasons, and an Enum would change the frame's dtype each
time. A new reason is an edit to one function, with no schema change. When a frame is exported for
readers without the library, the function's result is written out as a column at that point.

`u32` holds every human contig, T2T chr1 included (248,387,328). Closed vocabularies are `pl.Enum`.
`format_id` stays `str`, because `FORMATS` grows additively.

The schema carries a version in the Arrow/parquet metadata (`just_vcf.schema_version`), bumped by the
same rule as the package: a new nullable column is minor, and a rename or retype is major.

## Decision 7: contig names

| Option | Cost |
|---|---|
| A. Verbatim per source (`chr1`, `1`, `23`) | Every cross-facet join has to know every spelling |
| **B. One canonical spelling in `chrom`, the mapping recorded per facet in `contig_map`** | A mapping table to maintain. Canonical: `1`–`22`, `X`, `Y`, `MT`. `chrM`/`M`/`MT` are mapped explicitly, never with a `^chr` strip. Array codes `23`→`X`, `24`→`Y`, `26`→`MT`. T2T-CHM13v2.0's GenBank and RefSeq accessions (`CP068277.2`→`1`, `NC_060948.1`→`Y`, `CP068254.1`→`MT`; RM286), and GRCh37.p13's and GRCh38.p14's assembled molecules' (`CM000663.1`, `NC_000001.11`→`1`; F19) and NCBI36's RefSeq ones (`NC_000001.9`→`1`; NCBI36 has no GenBank accessions; RM4) map by the NCBI assembly reports, versions verbatim. Other contigs (alts, decoys, `chrUn_*`, scaffolds by any name) keep their source name |
| C. The seqcol's own names | Needs the seqcol at load time, and it is not always known then |

**Recommendation: B.** Array code `25` (the X pseudo-autosomal region) is kept as its own contig, `XY`
(PLINK's name for it), and not mapped to `X`, because that would drop the PAR fact. `contig_map`
records `25 → XY`. It ranks after `MT` (Decision 11).

## Decision 8: PGP freeBayes + VarScan2 merges

PGP's merged VCFs say which caller made a row in INFO/SF, and in about half the files, where both
callers found a site, VarScan2's genotype is gone (docs/PGP_HMS.md). The `…-orig` key=value table in
the same collection still has both.

| Option | Cost |
|---|---|
| A. One facet, caller `freeBayes+VarScan2`, SF passed through | Charter 8 wants disagreement between methods measured. Two callers in one facet hide it |
| **B. Two facets split on SF (`relation = caller_split`), sourced from the `-orig` table where present** | The loader needs to read the key=value table (it is already a declared format). Where only the merged VCF exists, its split loses VarScan2 on shared sites, and the facet has to say so |
| C. Merged VCF only, split on SF | Always loses VarScan2 on shared sites, even when the `-orig` table would have kept it |

**Recommendation: B.**

**Built 2026-10-22, from the merged VCF only.** `load` returns both facets in one `Loaded`, in
`facet_id` order, each with `relation = caller_split` and its own `caller` (read from the merge's
`##sourceFiles` names, `.FB.` and `.VS.`; null where they do not say). A row belongs to a caller when
INFO/SF lists its index; a record whose SF lists no caller's index (SF absent, or naming an index no
`##sourceFiles` line declares) belongs to no facet, and a load warning counts it. With a column per caller, each facet reads its own column. With one data
column, the column holds the first caller's genotype wherever that caller found the site, so the
VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the shared sites it
lost. `derived_from` stays null: the merged file is not loaded as a facet of its own (that is option
A), so there is no parent row to point at.

**The `-orig` key=value table, built 2026-10-22.** It loads the same way, as two `caller_split`
facets (`caller_index` 0 and 1 in the facet id, as for the merge), each holding the rows where its
caller's genotype is not `NULL`. It is not a second copy of the merge: on the one pair on disk it holds
2.04M of the merge's 5.30M rows, on every chromosome, the gene-dense ones most, and every row carries
gene annotation, so it reads as a gene-annotated subset. Within it, its freeBayes rows are the merge's,
and its VarScan2 facet keeps the genotypes the one-column merge lost (1.75M shared sites there). So
the two loads complement each other: the merge's VarScan2 facet has VarScan2's lone sites genome-wide,
the table's has VarScan2 at shared sites inside the subset. Nothing links the two loads: their facet ids
differ, and `derived_from` stays null on both, since neither file is made from the other. A table
facet's `callset_type` is inferred from its rows as for a VCF (RM4); the one on disk fits no class, so
it says `absence_means = unknown`.

**Two exports in one file (F3, built 2026-10-22).** One PGP-HMS upload is a headerless 23andMe-shaped
export followed by an export with one allele per column (AncestryDNA's data layout, without its
preamble). Both are one person: of the 164,378 ids both call, 164,275 have the same genotype
(99.94%). `load` returns them as two facets in one `Loaded`, as for the caller splits above: the
first under the file's options, the second with `relation = appended_export` and `part = 1` in its
facet id. Each is read by its own layout's rules (the second: `0` in an allele column is no call,
`source_gt` the two fields joined, chromosome codes 23-26 read as X, Y, XY, MT), and they are never
merged into one facet (charter 8). The second export is recognised by the rows: the file ends in at
least 1,000 rows of one other width, after other rows, and their first five are all in the two-allele
layout. A run of any other layout is not read, with a warning. Both facets keep the file's
`format_id`: the format names the file `sniff` read. Rejected, each with its cost:

- a registered headerless two-allele format id for the second facet: a new public id for a layout no
  file on disk has on its own (a headerless file of that shape sniffs `unknown`), added only to name
  half of one file; the facet's `format_id` would then name a file nobody gave;
- `format_id = ancestrydna` for it: the shape would stand in for a vendor statement the file never
  makes;
- `derived_from` pointing the second facet at the first: that column is lineage, and charter 5
  prefers a direct facet over a derived one, which would demote one export below the other for no
  reason in the data;
- reading the second export's chip with the 23andMe rule: that rule reads 23andMe's internal ids,
  which the second layout does not number. Its `chip` stays null;
- giving the second facet the build the file's preamble states (`assembly_label_basis = stated`): the
  preamble is the first export's, and nothing in the file states the second's build. Its label stays
  null until `resolve_ref` scores its own rows against the reference (F9, basis `inferred`).

## Decision 9: INFO and FORMAT fields the schema does not name

| Option | Cost |
|---|---|
| A. Drop them | Lean frame. Loses caller-specific fields (DRAGEN, CGI `CGIANN_*`) |
| **B. Drop by default and list them in `facets.dropped_fields`; opt in with `keep_fields=[…]` into an `extra: struct` column** | `extra`'s struct shape differs per facet, so it is not part of the one shared schema. That is fine while it stays opt-in and the name says so |
| C. Always keep everything as `list[struct{key, value: str}]` | One schema, but millions of rows times dozens of stringly keys on WGS |

Kept keys are named with their namespace (`info.DP`, `format.DP`), because the bare token is ambiguous
(see the probe above). Some fields VCF 4.5 reserves are not genotype calls, and they land here too:
base modifications (`M5mC` and the rest of `M[0-9]+[ACGTUN]`), copy number (`CN`, which is a float
since 4.4), and tandem repeat counts (`RUC`, a float with an interval that can be open above). A
methylation or copy-number facet is plausible later. Nothing here rules it out, and nothing here
designs it.

**Recommendation: B.** Under P1, a key missing from the header never reaches the loader at all. The
loader builds complete headers for the formats it synthesises (`vcf_no_meta`, the key=value table).

**A carried key the reader cannot follow is not carried** (2026-10-24), in three cases, each with a load
warning that names the key and with the key in `dropped_fields`:

- declared with a Number it is not read under: a vector over the alleles or genotypes with a fixed
  count (`AD`, `PL`, `GL`, `GP`, `PP`), or HQ with any Number but `2`, `P` or `.`;
- declared with a Type other than the one its VCF version reserves, in a 4.3, 4.4 or 4.5 header
  (`HQ` `Type=Float`; `CN` `Type=Float` under 4.3, which reserves it Integer). polars-bio's reader
  holds a definition per version and refuses the whole header for such a line; it checks only the
  Type, and nothing under 4.2. A key it has no definition for in that version (`DS`, `MIN_DP`, `LEN`
  before 4.5, `PSL`/`PSO`/`PSQ` before 4.4) is read under any Type, except a text key below;
- a text key (`FT`, `PSL`, INFO `SVTYPE`) declared with a Type other than String, in any version
  (2026-10-25): its values (`PASS`, a set name, `DEL`) do not parse as an Integer, a Float or a
  Character, and under VCF 4.2, where polars-bio refuses no Type, the first one failed the whole read. The warning names the key
  and its Type.

polars-bio parses every key its header declares, requested or not, so a value that does not parse
under such a declaration (`51,40` under HQ `Number=1`) would fail the load. The loader's retry reads a
copy whose header leaves all three out, and the rows' values under them are skipped unread, whatever
they are (2026-10-25). Under VCF 4.3, 4.4 and 4.5 leaving the line out is not enough: polars-bio's
reader applies its own definition of the version's reserved keys to one the header does not declare,
and parses a single value (DP, GQ, PS, CN, INFO END) as the record is read, requested or not, so DP
`Type=Float` under 4.3 loaded with `30` and failed the whole load with `30.5` ("invalid digit found in
string", naming no key). So where the version holds a definition of a key left out, the copy declares
it `Number=.` with the Type the version reserves (`skipped_declarations`, `REWRITE_VERSION` 16): the
reader checks only a line's Type, and it parses a list value only when the value is requested, which
a key not carried never is. Rejected: declaring it `Type=String`, which the reader refuses as it
refused the file's own line; removing the key's slot from every record, which rewrites the data
polars-bio reads and costs a pass over each sample column; and writing the copy as VCF 4.2, whose
reader holds no definitions, which would read the whole file under another version (GT's leading
phase indicator, P4).

The same holds for a reserved key the loader never carries (2026-10-26): the structural copy keeps
only the carried keys' lines, so FORMAT MQ or INFO DP, declared by the file or used undeclared, was
read under the version's definition, and MQ `60.5` under 4.3's Integer failed the retry. So the copy
declares every key the version reserves and the copy does not declare otherwise `Number=.` with the
reserved Type (`READER_DEFINITIONS`, the pinned noodles fork's own table; `REWRITE_VERSION` 17), every
one, not only those the rows read ahead use, since a later record may use another and a declared key
no record uses costs nothing. Two are not declared: GT, which the reader parses as a genotype whatever
its line says, and a Flag (INFO DB, SOMATIC), which the reader holds under no Number but 0. So under
4.3+ the reader holds a reserved Flag as one whatever the file declares, and a Flag written with a
value (`DB=1`) fails every read with `HeaderRefused`, even where the file's own line makes the value
valid (`##INFO=<ID=DB,Number=1,Type=Integer>`): the reader refuses that line, the copy leaves it out
like any key not carried, and the reader reads `DB` as its Flag ("invalid flag"). polars-bio's error
names no key, so the refusal names it from the records read ahead (2026-10-26): `INFO/DB written with a
value`, with `HeaderRefused.key` set where it is one key; where no record read ahead writes a Flag with
a value, it says so and lists the version's Flags. Under 4.2, whose reader holds no definitions, the
copy reads such a key as text and the file loads. Rejected: dropping the value from the copy's records,
which rewrites the data polars-bio reads, as for the keys above.

GT is never left out (2026-10-25): a facet without genotypes is not a calls facet, and GT's values are
genotype text whatever its line says. Declared `Integer`, `Float` or `Character`, it is read as written
through a copy whose header declares it `Type=String`, in every version, with a warning naming the Type.
Under 4.2 polars-bio read it so already; under 4.3, 4.4 and 4.5 it refused the header, and the retry
that left GT out had no row with a genotype. A BCF is read as stored, with no copy, and polars-bio
refuses a BCF header whose GT is not String in every version, so such a BCF is refused with
`HeaderRefused` naming GT and its Type.

One more case is found only by reading (2026-10-25): a carried key declared as its version reserves it,
with a value that does not parse under that declaration (HQ `Integer` with `30.5,40` or `x,40` under
VCF 4.2). Where the structural copy is refused too and polars-bio's error names the key, the copy is
written again without it and read again; the key is not carried, with a warning that names it and
polars-bio's reason, and it is in `dropped_fields`. Each key is left out once at most, and GT never, so
the retries end. polars-bio names the key only for a list-valued field (HQ, AD, PL, INFO SVLEN…): a
single-value one (GQ `30.5` under `Integer`, an INFO END) fails inside the reader with no key named,
and the load is refused as before, since which key failed would be a guess.

An INFO value takes the same path as a FORMAT one, in every version (2026-10-26). polars-bio parses a
record's INFO only when an INFO column is read, and the pass the retry wraps read none where the
carried INFO keys were ones it counts nothing of (SVTYPE, a split marker), so the first read of INFO
was the spill's, after the retry: INFO DP `Type=Integer` holding `30.5` beside SVTYPE failed the load
with polars' raw `ComputeError`. That pass now reads every INFO column the spill reads, in the same
pass over the file. The value is in a key not carried, so the structural copy, which leaves DP's line
out, loads the file with DP in `dropped_fields`; one in a carried key is left out by name or refused
as above.

**A carried integer key declared `Float` or `String` keeps only whole numbers** (2026-10-25). VCF 4.2
reserves no Type, and a later version none for a key it holds no definition of (MIN_DP; LEN before
4.5; PSO and PSQ before 4.4), so such a header can declare GQ, DP, MIN_DP, LEN, PS, AD, PSO, PSQ, HQ or
INFO END, SVLEN, CIPOS or CIEND `Float` or `String`. A whole number is read as it is (`30.0` is 30); a
number with a fraction (`30.7`) is left out, never truncated or rounded; text that is no finite number
(`x`, `nan`) is left out; and a whole number past what the value was read through holds exactly
(2^53, a 64-bit float) is left out too. Each key and reason is one warning with its count. A PS left out
withholds the phase of a phased record (`phased` null where it was true): a phased record without PS
is in its chromosome's implicit set, and saying so would be a claim the file did not make. A PS not
carried at all (any of the cases above) withholds the phase of every phased record (Decision 2). PL and
PP declared `Float` are read as the numbers they are; a fraction there is a likelihood, not a count cut
short.

**Whole is decided from the digits the file wrote, never from a float** (2026-10-25). polars-bio reads
a `Float` as a 32-bit float, and no value it hands back can be certified whole, at any magnitude: a
32-bit float has 24 significant bits, so `30.0000001` (30 + 1e-7, under half the step of 1.9e-6 at
30) and `5000000.1` (under half the step of 0.5 above 2^22) both come back whole, and above 2^24 the
whole numbers themselves are rounded (PS 123456789 came back 123456792). A threshold such as 2^22,
from which the step is 0.5 or more, catches the large values only: below it a fraction finer than half
the step still rounds away (`30.0000001`), and a file can write as many digits as it likes. 64 bits
only moves the problem (`30.00000000000000001` is 30.0 there).
So every copy polars-bio reads declares these keys `Type=String` where the file says `Float` and its
version lets polars-bio read the line (`retyped_whole`; `REWRITE_VERSION` 15), and the text is judged
by its digits (`_whole_value`): whole when the fraction digits, past the exponent's shift, are all
zeros (`30.000`, `1e2`, `120e-1`), its number then the 64-bit parse, exact below 2^53. Rejected: a
magnitude threshold on the 32-bit value (2^22, 2^24), which keeps `30.0000001` as 30, a value the file
did not write; and leaving every Float-declared integer key out, which drops the whole numbers most such
files write. A BCF stores a `Float` as a 32-bit float, so its value is the file's own; whole where it
has no fraction, and left out from 2^24 as before. A key its version refuses as `Float` (PS under 4.3)
is not carried (above).
One load warning names the keys a copy so declares, each with its declared Type (`FORMAT/DP
(Type=Float)`), as GT's retype has one; a key not carried is named by its own warning instead, and a
BCF, whose header is not rewritten, gets none.

## Decision 10: multi-sample files

The charter is about one person. A trio VCF or a PGP two-column merge has more than one sample column.

| Option | Cost |
|---|---|
| **A. One facet per sample. With more than one sample, the caller names which (`samples=[…]`); without that the loader refuses and lists the sample names** | One extra argument for multi-sample files. Nothing is loaded by accident |
| B. Load every sample as its own facet by default | Easy to end up holding a relative's genome without asking for it |
| C. Keep samples as a struct of lists (polars-bio's native shape) | Breaks the one-schema rule |

**Recommendation: A**, except for PGP merges (Decision 8), where the columns are one person split by
caller.

## Decision 11: canonical order (charter 12)

Frames are emitted sorted by `facet_id`, then contig rank, then `pos`, `end_pos`, `source_ref`, `alts`. Contig rank
follows the seqcol order when it is known. Otherwise it is `1`–`22`, `X`, `Y`, `MT`, `XY`, then the remaining
names sorted as strings. Nulls sort last. Every loader gets a test that loads a shuffled input and
compares it with the sorted one.

## Decision 12: the `links` frame (charter 5, RM275)

The facet graph's edges are a fourth lazy frame, `Loaded.links`, one row per link between two facets.
`facets.derived_from` holds one parent; a merge or a consensus has several, and two facets read from one
source are related without either being made from the other. The person chose a small links frame
over a list column on `facets` (2026-10-23, E4). just-vcf's own shape (charter 2), except for the
relation terms GA4GH defines.

| Column | Meaning |
|---|---|
| `from_facet_id` (required) | the facet the link starts at: for lineage, the facet the step read |
| `to_facet_id` (required) | the facet it ends at: for lineage, the facet the step made |
| `relation` (required) | one of the relations below. An open vocabulary (`pl.String`); `just_vcf.links.RELATIONS` lists the documented values |
| `provenance` | `struct{step, step_version, derivation_inputs}`: the public function that made the link (`load`, `rejoin_split`, `normalize`, `liftover`), the version of its rule as a string (the `REJOIN_VERSION`, `NORMALIZE_VERSION` or `LIFT_VERSION` its `facet_id` digests; null for `load`, which has no versioned rule and digests none), and the files this step read besides its parent and its reference collection, in `facets.derivation_inputs`' shape (a liftover's chain; `[]` when it read none) |

**Direction.** A lineage link points from the parent to the facet made from it, as Cat-VRS reads
`liftover_to`: from the defining allele to the one related to it (the term verbatim, charter 2). So
`facets.derived_from = P` on facet `C` is the link `P → C` with `C`'s relation.

**Relations.**

| Relation | Kind | Written by | Meaning |
|---|---|---|---|
| `liftover_to` | lineage | `liftover` | Cat-VRS's term, verbatim: `to` is `from` lifted to another assembly |
| `normalized` | lineage | `normalize` | just-vcf's: `to` is `from` split, left-aligned, trimmed and deduplicated |
| `split_rejoined` | lineage | `rejoin_split` | just-vcf's: `to` is `from` with its pre-split sites joined again |
| `merge` | lineage, several parents | none yet | just-vcf's: `to` was made from every `from` linked to it by putting their rows together (a merge) |
| `consensus` | lineage, several parents | `consensus` | just-vcf's: `to` holds a posterior per allele computed from every `from` linked to it (RM280); the parents keep their own calls |
| `caller_split` | kin | `load` | just-vcf's: two facets split by caller from one merge (Decision 8), the merge itself no facet |
| `appended_export` | kin | `load` | just-vcf's: `to` is the export appended to a file after `from` (Decision 8, "Two exports in one file") |

GA4GH defines no term for the others; Cat-VRS's relations relate alleles across assemblies and
molecules (`liftover_to`, `transcribed_to`, `translated_to`), and only the first fits a facet.

**Kin links.** A caller split's facets are linked in every pair, from the lower `facet_id` to the
higher: the merge they came from is not a facet (Decision 8, option A), so there is no hub to link
them to. An appended export is linked from the first export to the appended one, file order. Neither
is lineage: `derived_from` stays null on both ends.

**How `links` and `facets` agree** (`just_vcf.links.link_disagreements`, tested on every step):

- `derived_from` set: exactly one lineage link ends at the facet, from `derived_from`, with the
  facet's `relation`;
- `derived_from` null and `relation` `merge` or `consensus`: two or more links of that relation end at
  the facet, and no other lineage link;
- `derived_from` null and a kin `relation`: no lineage link ends at the facet, and a link of that
  relation names it;
- both null: no lineage link ends at the facet.

`facets.derivation_inputs` and a link's `provenance.derivation_inputs` can differ, so the rule above
does not compare them: as the steps write it today, the facet's is carried from parent to child
(`rejoin_split` and `normalize` copy the parent's row, so a facet normalized after a liftover still
names the chain), while the link's names only the files its one step read. Whether the facet's
column should be per step too is open.

**What a `Loaded` holds.** A plain load: no links, except the kin links of a load that returns two
facets or more. A derived step returns its parent's links and the one it adds, so a facet three steps
from its file still says how each step got there when only the last `Loaded` is kept; its links can
name facets its own `facets` frame does not hold. `resolve_ref` and `mint_vrs_ids` keep the facet's
`facet_id` and pass its links through. Links are sorted by `from_facet_id`, `to_facet_id`, `relation`,
each key once (charter 12).

Rejected:

- a `parents: list[str]` column on `facets`: one row per facet cannot carry a provenance per parent,
  and a kin pair would be written twice;
- links only into the facets of the `Loaded` at hand: a chain `load → normalize → liftover` would lose
  its first step's link the moment the middle `Loaded` is dropped;
- a kin link from a star's hub: a caller split's merge is no facet, and inventing one would be the
  option Decision 8 rejected;
- the reference collection in `provenance.derivation_inputs`: it is a seqcol digest, not a file, and
  the `to` facet's `ref_seqcol_digest` already names it;
- just-vcf's own package version in `provenance`: a new release would change the links of facets
  whose `facet_id` it does not change.

## Other VCF 4.5 points the loaders have to respect

None of these changes a column. Each is a rule a loader follows, and each gets a test.

- **REF may be a reduced IUPAC code.** The spec requires an ambiguous reference base to be written as
  the alphabetically first concrete base (`R` becomes `A`). When RM3 checks a source REF against a FASTA
  that holds `R`, a mismatch there is expected and is not a coordinate error.
- **POS 0 and N+1 mark telomeres.** `pos` allows 0. No VRS id is minted for either.
- **For a symbolic SV, POS is the base before the event.** For `<*>` the block includes POS. The
  `end_pos` rule above follows both.
- **GQ is conditional** on the site being variant. It is not the confidence of a reference call.
- **Percent encoding.** `%3B`, `%3D`, `%2C` and the rest can appear in INFO values and must be decoded.
  Whether polars-bio decodes them has not been checked yet.
- **Ploidy is whatever GT says.** GT can have one allele (haploid Y, male X, MT) or three and more, and
  the genotype order for `gl`/`gp` follows the spec's general formula for any ploidy. Without GT, the
  spec says to assume diploid.
- **`ID` is a `;` list**, split into `ids`.

## Stress case: XYY with a Y deletion the reference does not have

A person with 47,XYY (Jacobs syndrome) carries a deletion on Y. Their Y is shorter than the
reference's, and the deleted region lies where GRCh38 Y has no sequence at all (an unassembled gap,
all `N`), while T2T-CHM13 has it. Facets: a GRCh38 short-read VCF, a T2T HiFi VCF, and the T2T facet
lifted to GRCh38. How the draft holds up:

- **Ploidy.** A caller run with a "male means haploid Y" model writes one allele per Y site, although
  this person has two Ys. `gt_alleles` keeps what the file wrote, because list length is never an
  asserted ploidy (Decision 2), so nothing is silently "corrected". **Open:** nothing in the draft can
  *say* the karyotype. It is a fact about the person, measured from evidence (Y depth near 2x, X depth
  near 1x). It is not a fact about any one facet. A person-level finding with its own provenance does
  not exist yet. Since RM269 each facet's own evidence can be measured (below); it reads this person's
  X as one copy and Y as present, so `XY` from the VCF alone: the Y depth of a variant-call VCF does
  not separate one Y from two (`y_copy_count` stays null). With a coverage facet of the same person
  the X-degenerate Y reads near 1.0 of chr20's depth, two copies, and the answer is `XYY`.
- **The deleted region, absent from the GRCh38 VCF.** The GRCh38 facet is WGS without blocks, so
  `absence_means = weak_reference`, and a question about a site in the deletion gets "weak evidence of
  reference". For a deletion on both Y copies, that is wrong. **Open:** absence has to be conditioned
  on a coverage facet (charter: coverage is a facet). Zero depth there turns "weak reference" into "no
  evidence". The draft names the coverage facet but does not yet say how absence reads it.
- **The same deletion on T2T.** Called against T2T, it has a REF from real sequence and can get a VRS
  id (or, if it is written as `<DEL>`, `symbolic`). Lifting it to GRCh38 fails, because the target is
  a gap. The derived facet has no row for it. **Open:** a lift that drops a row has to be counted on
  the derived facet (the Glow-style liftover status in LANDSCAPE). Otherwise "absent in the lifted
  facet" reads as "absent in the person".
- **A GRCh38 REF that is `N`.** A site inside the gap has no base to compare, and `ref_check`
  compares against `N`. Where the source REF is `N` too, the id is minted: gtars and vrs-python both
  justify an allele over `N` bases and agree on the id (tested, RM271). Were gtars to refuse, the
  refusal would land in `vrs_errors` as `refused`.

None of this changes a column today. It marks where the person-level layer and the coverage facet have
to connect.

## Facets of one person (RM273)

The person chose "supplied, checked" (Part 2 interview, E1): the caller says whose a file is, the
library checks the claim and warns, and it never infers a subject or moves a facet to another one.

- **Supplied.** `load(..., subject=...)` writes the label into `facets.subject` on every facet the file
  loads as (Decision 5). Derived steps copy their parent's row, so the label carries down. Nothing
  else writes it.
- **Grouped.** `Person.of(*loaded, check=...)` gathers facets that all carry one subject: their
  `facets`, `calls` and `blocks`, concatenated in `facet_id` order, which keeps the canonical order
  (Decision 11), and `links`, every link the given loads carried, once each in `canonical_links`
  order (Decision 12). `people(*loaded, check=...)` returns one `Person` per subject, sorted by subject. A
  facet without a subject, or two subjects handed to `Person.of`, is refused with `SubjectError`; one
  `facet_id` given twice is refused, since its rows would count twice.
- **Checked.** `check` is a relatedness test over two facets, each handed in as a `Loaded` of that one
  facet, answering `True` (one person), `False` (two) or `None` (could not tell). Facets of one source
  (the same `source_sha256`, `member` and `sample`: a file, the facets derived from it, or the facets
  one file loads as) are one person's data by construction and are not checked against each other.
  Every other pair of sources is checked once, each source standing as its facet read from the file
  (else its lowest `facet_id`). `Person.checks` holds one row per pair (`facet_id`, `other_facet_id`,
  `outcome`: `same`, `different`, `could_not_tell`, `not_run`), and `Person.one_person` is their Kleene
  AND: `False` when any pair is `different`, `True` only when every pair is `same`, `None` otherwise,
  one source alone included. `different` and `could_not_tell` pairs are warnings naming both facet ids;
  the facets stay where the caller put them. The default check is `just_vcf.relatedness.same_person`
  (RM274): KING-robust kinship over the package's site panel, `True` when its interval lies above the
  same-person floor, `False` below it, `None` across it or when the pair cannot be compared (no shared
  assembly, too few shared sites). `check=None` runs none and leaves every pair `not_run`, which is not
  a pass.
- **Coverage facets** (RM278). `read_coverage(..., subject=...)` and `alignment_coverage(...,
  subject=...)` write the label on the `CoverageFacet` (`subject`, nullable, not digested into its
  `facet_id` and not in the cache spill's key). `Person.of(*loaded, coverage=[...])` and
  `people(*loaded, coverage=[...])` carry them on `Person.coverage`, in `facet_id` order. One without
  a subject, under another subject than the facets, or (in `people`) under a subject no facet carries
  is refused with `SubjectError`; one given twice with `ValueError`. The label is their only check: a
  depth file holds no genotype, so no relatedness test reads one and `Person.checks` has no row for it.
- **One source under two subjects** (a file loaded twice under two labels) is allowed: `people` warns
  in each person concerned, naming the facet ids. Two samples of one VCF under two subjects (a trio)
  are different sources and draw no warning.
- **The label is not identity.** `facet_id` does not digest it, so it cannot split one file's data into
  two ids, and relabelling changes no `facet_id`, its own or a derived one's. A derived step's cache
  key digests the facet row, subject included, so one step run on one file under two labels writes two
  cache generations: disk, never correctness (rejected: stripping the subject from that key in
  `just_vcf.cache`, a shared module, for a cost only relabelling pays). The subject identifies a person if
  the caller makes it so; warnings and errors name facet ids, never the subject. Rejected: a subject
  inferred from relatedness (the person's option b), which E6's relatedness may suggest but never
  assign; a subject read from the sample name, which is a file's label, often a real name, and the
  same for different people across files.

## Relatedness and pedigrees: computed or supplied, not stored (RM274)

Neither is a column of `facets`, `calls` or `blocks`. Relatedness is computed on demand from two facets'
`calls` (and `blocks`, for a gVCF's absences) and returned; a pedigree is a document the caller
supplies, read into frames of its own. docs/RELATEDNESS.md has the method.

- **`RELATEDNESS_FRAME`** (`pairwise_relatedness`, `just_vcf.somalier`; one `Relatedness` model per
  row): `facet_a`, `facet_b`, `method` (`king_robust` or `somalier`, an open vocabulary), `band` (a
  `pl.Enum`: `same_person`, `first_degree`, `parent_child`, `full_siblings`, `second_degree`,
  `third_degree`, `unrelated`, `unknown`; a first-degree kinship reads `parent_child` or
  `full_siblings` where the IBS0 share's interval lies wholly on one side of KING's cutoffs, and
  `first_degree` where it does not, the kind not told),
  `reason` (why `unknown`: `assembly_unknown`, `assembly_mismatch`, `no_panel`, `too_few_sites`,
  `no_heterozygous_sites`, `between_bands`, `variant_sites_only`; an open vocabulary), `assembly_label`, `panel`, `kinship`
  and `kinship_se` (KING's scale, one person 0.5), `ibs0` and `ibs2` (shares of the shared sites), the
  counts `n_sites`, `n_het_a`, `n_het_b`, `n_het_both`, `n_ibs0`, `n_ibs2`, `n_vrs_compared`,
  `n_vrs_differ`, and `somalier_relatedness` (somalier's own value, verbatim). A pair not compared at
  all (no shared assembly, no panel) has null counts: not measured is not zero. The bands are
  just-vcf's own terms (KING's inference criteria); GA4GH's KIN ontology names stated relationships,
  not degrees measured from genotypes.
- **`PEDIGREE_INDIVIDUALS`** (`pedigree_id`, `individual_id`, `sex`) and **`PEDIGREE_RELATIONSHIPS`**
  (`pedigree_id`, `individual_id`, `relation_id`, `relation_label`, `relative_id`): a GA4GH Pedigree as
  the caller supplied it, the KIN term and sex value as written. `individual_id` is `relation` of
  `relative_id`, the standard's direction.
- **`EXPECTED_KINSHIP`** (`expected_kinship(pedigree)`): one row per pair of the pedigree's
  individuals, in its order: `pedigree_id`, `individual_id`, `relative_id`, `relation_id` and
  `relation_label` (the first statement relating the two, verbatim, null when none does), `basis`
  (`stated`, `derived` through parent links, or `assumed`: no biological path, founders taken as
  unrelated), `expected_kinship` (null unless the pedigree allows one value), `expected_bands` (every
  band the allowed values fall in, least related first; a first-degree value is `parent_child` for a
  listed parent link, `full_siblings` for two shared listed parents or a stated sibling term, and
  `first_degree` reached any other way) and `expected_reason` (why there is none:
  `term_without_degree`, `more_than_two_parents`, `statements_conflict`; an open vocabulary). Which KIN
  term allows which kinship is just-vcf's own reading (KIN has no kinship values); the
  `just_vcf.pedigree` docstring lists it.
- **`PEDIGREE_CHECKS`** (`check_pedigree(pedigree, *people_or_loaded)`): `EXPECTED_KINSHIP`'s columns,
  then `facet_a` (the individual's) and `facet_b` (the relative's), the measured `band`, `kinship`,
  `kinship_se`, `ibs0`, `n_sites` and `relatedness_reason`, and `agreement` (`agrees`: every band the
  measured interval reaches is allowed; `disagrees`: none is; `cannot_tell`), with `reason` for
  `cannot_tell` (`person.PEDIGREE_CHECK_REASONS`: `no_facet`, a pedigree reason, a relatedness
  reason, `between_bands` when the interval reaches allowed and other bands, or `kind_not_told`
  when it lies wholly in the first degree and the IBS0 count does not tell the kind the pedigree
  needs; an open vocabulary, checked where the row is built). The two first-degree kinds are bands of
  their own here, and `first_degree` on either side stands for both: a parent link measured
  `full_siblings` disagrees, one measured `first_degree` (the kind not told, a thin IBS0 interval or
  no count) is `cannot_tell`, `kind_not_told`; it read `between_bands` until 2026-10-26, which said the
  kinship crossed a band boundary when it did not. A pedigree individual is matched to facets by its
  `individual_id` as `subject` (RM273), each source once as `Person.checks` stands it; one row per pair
  of those facets, one row with null facets when an individual has none. The pedigree is never changed
  and nothing is relabelled (charter invariant 8). `people(..., pedigree=...)` runs it and warns each
  person in a disagreeing pair, naming facet ids and bands, never a subject.

## Sex-chromosome evidence: computed per facet, not stored (RM269)

`sex_chromosome_evidence(loaded)` returns one `SexChromosomeEvidence` per row of `loaded.facets`. Like
`vrs_withheld()`, it is a function over the frames, not a column, so it adds nothing to the schema and
needs no version bump. Rejected: a struct column on `facets`, which would cost a full scan of `calls`
at every load, orphan every cache spill with its version bump (RM268) and freeze today's thresholds
into stored data. There is no person in the schema yet, so nothing person-level is stored either;
weighing facets against each other (liftovers of one callset are one piece of evidence) is for that
later layer.

What it reads, over `assay` rows no filter failed (`filters` and `sample_filters` null or `["PASS"]`;
null is VCF's "filters not applied"):

- `ploidy`: rows by the number of haplotypes their genotype was written with, per group (autosomes,
  non-PAR X, non-PAR Y, PAR, MT). What the source wrote, not a measured ploidy.
- `x_heterozygosity`: on non-PAR X, the het share of fully called two-allele genotypes with an ALT
  (rows with a `split_signal` left out), the count of one-allele ALT genotypes (the caller's
  assertion, never counted as evidence), and the median het allele balance from `ad`.
  `y_heterozygosity`: the same on the X-degenerate Y, reported only.
- `depth`: the median `dp` of ALT-carrying rows on non-PAR X and on the X-degenerate Y, each over the
  same median on chr20.
- A row whose source REF states no base and which spells a sequence ALT (`ref_check = unstated`,
  Decision 2) carries no ALT for the two lines above (2026-10-31): its genotype was spelled against
  the `N`. Its ploidy as written still counts.
- `caller` and `ploidy_statements`, the facet's own columns, carried through unread.

It answers `x_copies` (`one`, `two_or_more`) and `y_present` separately, so a facet whose Y line is
unusable still says how many X copies it holds, and `inferred_karyotypic_sex`, a GA4GH Phenopackets
`KaryotypicSex` value (charter invariant 2: the vocabulary is the standard's, the field is
just-vcf's, since Phenopackets has no place for an inference). From calls alone only `XY`, `XXY` and
`UNKNOWN_KARYOTYPE` are estimated, the last always with a reason: `y_present` is `True` or unknown,
never `False`, because depth at called sites cannot show a Y absent (a man's GATK genome read below the
present band, close to the two-X range), so `XX` and `XO` wait for a coverage facet. The PAR bounds and the X-degenerate
Y come from GIAB's stratifications for GRCh37 and GRCh38, for T2T-CHM13 from GIAB's (PAR, XTR) with
the T2T consortium's X-degenerate class, and for GRCh36 from NCBI36's PARs as Ensembl 54 records them
with GRCh37's X-degenerate Y lifted through UCSC's hg19ToHg18 chain (docs/PROVENANCE.md); a label
with none of these (NCBI35, or one a caller wrote) answers unknown (`unsupported_build`). The
thresholds were measured on GRCh37 and GRCh38 genomes only, so a T2T-CHM13 facet, or a GRCh36 one
that is not an array, is grouped and measured, and its estimate withheld (`uncalibrated_build`, RM286, RM269) until callsets of known
karyotype on it are measured. A coverage facet on such a build is set aside with the same reason.
Rules and thresholds: `just_vcf/karyotype.py` and docs/ROADMAP.md RM269.

`sex_chromosome_evidence(loaded, coverage=…)` takes a coverage facet (RM278) of the same person and
build, paired by the caller. Its `coverage` field (`just_vcf.karyotype.CoverageRatios`, module-level
like the result's other nested models) holds the length-weighted median read depth over non-PAR X,
over the X-degenerate Y and over chr20 (zero-depth intervals left out on chr20 and X, where they are
assembly gaps; kept on the Y, where zero is the signal), from intervals wholly inside each region and
described at every base, and the two ratios. It is a third X line and a second Y line, and the only
line that counts past two: copy number puts each copy at half chr20's depth, so it reads one, two or
three X (0.5, 1.0, 1.5) and none, one or two Y (0, 0.5, 1.0), each within 0.125, and anything
between or beyond those bands is unknown (`intermediate`). Each line allows a range of counts (the het
share one, or two or more; the calls' X depth one, two, or above 1.25 three or more; the calls' Y
depth one or more) and the answer is what all of them allow: none in common is `lines_disagree`, as
is coverage showing no Y where the calls' depth shows one. Where it is one number it is
`x_copy_count` or `y_copy_count` (nullable ints; `x_copies` and `y_present` keep their meaning).
The counts name `XO`, `XY`, `XYY`, `XX`, `XXY`, `XXYY`, `XXX`, `XXXY`, and `OTHER_KARYOTYPE` for
three X with two Y, which Phenopackets does not list. Without a count, the earlier reading stands:
two or more X read two and a Y the calls show reads one, so the calls alone still name `XY` and `XXY`
(they cannot tell XY from XYY); three or more X that nothing counts are `copies_uncounted`. No XXX or
XYY genome was measured: the bands are copy number alone. A coverage facet without a build, or on
another one, is set aside with `coverage.reason` (`coverage_no_build`, `coverage_build_mismatch`),
and the answer is the calls' alone. So is a coverage facet whose chr20 depth is too uneven for a whole
genome (`coverage_uneven`: under 75% of its bases within half to twice the median, an exome's or a
panel's shape). The field defaults to null, so nothing changes without coverage.

**At the person level**, `person_sex_chromosome_evidence(person, pairs=None, array_y=None)` answers
for each facet of a `Person` (RM273) with one of the person's coverage facets paired to it, and the
result's `coverage_pairing` (`just_vcf.karyotype.CoveragePairing`: `coverage_facet_id`, `basis`,
`candidates`, `detail`; null from the facet-level function) says which and why. Each answer is
exactly what `sex_chromosome_evidence` gives that facet with that coverage facet. Pairing is by build,
since nothing in the frames records which alignment a callset was called from: the person's one
coverage facet on the facet's `assembly_label` is used (`only_one_on_build`); of several there, the
one `measure_coverage` finds usable (`only_usable_on_build`, the others named with their reasons);
two or more usable is `several_on_build`, none usable `none_usable_on_build`, none there
`none_on_build`, none at all `no_coverage`, and in those cases no coverage facet is used and the
answer is the calls' alone. `pairs` (a facet's `facet_id` to a coverage facet's, or None) states the
pairing instead (`stated`). The vocabulary is `just_vcf.karyotype.COVERAGE_PAIRINGS`. Rejected: picking
one of several usable coverage facets (the deepest, the first), a guess that would hide another
person's BAM labelled with this subject, which `lines_disagree` otherwise shows; pairing on the BAM's
read-group sample against the VCF's sample column, a file's label that is often a real name and the
same for different people (the reason the subject is never read from it); a blend of several coverage
lines.

An array facet (`callset_type = array`) gets its own lines, on the result's `array` field
(`just_vcf.karyotype.ArrayLines`, null for every other facet), read from `gt_alleles` as the export
wrote them, so no REF and no `gt` are needed: rows, called rows, written no-calls, one-letter rows, two-letter base
calls and their hets, per autosomes, non-PAR X and the X-degenerate Y. X copies come from the het
share of two-letter base calls on non-PAR X over the autosomes' het share (`x_het_ratio`): one X at
most 0.1875, two or more from 0.5625, the central half of the gap between 0 (one X, genotyping errors
only) and 3/4 (two X, the neutral expectation of X diversity against an autosome's). An X written
with one letter on most called rows is the vendor's decision (23andMe), counted, never read
(`caller_asserted_only`). Y comes from the X-degenerate Y probes' called share over the autosomes'
(`y_call_ratio`): absent at most 0.25, present from 0.75. Whether it is read is decided per facet by
its `vendor`: a vendor that no-calls Y for the samples it calls female makes the line its own sex
call, and the file cannot tell that mask from a chip that saw no Y. By default
(`sex_chromosome_evidence(…, array_y=None)`) the line is read for a vendor in
`just_vcf.karyotype.ARRAY_Y_VENDORS`, AncestryDNA alone, whose PGP-HMS exports call what the chip
measured (its women call 0.8–9.9% of the Y probes, 23andMe's none), and for no other vendor nor a facet
whose `vendor` is null (the person's decision, desk #136, 2026-10-26). `array_y=True` reads it for
every array facet, `array_y=False` for none; a person with an AncestryDNA and a 23andMe facet gets the
first's Y read and not the second's under the default. Not read, `y_reason` is `vendor_asserted`, the
`detail` says whether the vendor's default or the caller turned it off, and the karyotype is unknown. An export with no no-call row
on the autosomes drops its no-calls, so its rows are not its probes and its Y share says nothing:
`y_call_ratio` is null and `y_reason` `no_calls_not_written`. A no-call row is one written so (`--`,
a letter null, `autosome_no_calls`); a row whose genotype could not be read (`gt_alleles` null, a
code no loader reads) is not called but says nothing of the export's no-calls (2026-10-25: one such
row read a no-call-dropping export's Y). Arrays never count copies, so
they name at most `XO`, `XY`, `XX` and `XXY`. They are estimated on GRCh36 as on GRCh37 and GRCh38
(`just_vcf.karyotype.ARRAY_ESTIMATE_BUILDS`), where a genome is not: their bands are no measured
threshold, and GRCh36's 23andMe exports read as the GRCh37 ones do (RM269). Computed, not stored: no schema change.

The build is the facet's `assembly_label` as the frames carry it, whatever in the library filled it:
the source's own statement at load, the header's `##contig` lengths at load (RM4), an array's rows
scored against the reference by `resolve_ref` (F9), or a liftover's target. The estimate never infers
or defaults one, and a facet without one answers `no_build`. An array table that states no build
(FTDNA, a bare rsid table, a Final Report) therefore answers `no_build` as `load` returns it and is
estimated on `resolve_ref`'s output, which carries the build its rows gave. The result says where its
build came from: `assembly_label_basis` (the facet's column verbatim), `assembly_label_source`
(`stated`, `contig_lengths`, `rows`, `liftover`; `just_vcf.karyotype.ASSEMBLY_LABEL_SOURCES`), read
off the facet's columns by Decision 5's rule (`relation = liftover_to`, then a length in
`header_contigs`, then `ref_seqcol_digest` on a facet that is no VCF), and null for a label no rule
accounts for (one a caller wrote into the frame, which is still used), and `ref_seqcol_digest`, the
collection REF was read from. These three are carried and decide nothing. Rejected: a `reference`
argument that runs F9's inference inside the estimate (a REF read hidden in a karyotype call, and a
second path to a build `resolve_ref` already records), and widening `assembly_label_basis` to name the
three inferred sources (a dtype change and a schema bump for what the columns already tell apart).

## CisPhasedBlock export: computed on demand, not stored (RM284)

The person chose columns with an export on demand: phase sets stay the three `calls` columns of
Decision 2, and `cis_phased_blocks(loaded, reference_alleles=False)` returns a lazy frame of VRS 2.x
`CisPhasedBlock`s when asked. It adds no column to the three frames and needs no version bump. Rejected:
a stored block table, which every load would pay for and every change to an allele id would orphan.

| Column | Dtype | Null means |
|---|---|---|
| `facet_id` | `str` | never null |
| `chrom` | `str` | a PSL set, which is keyed per file and may span chromosomes |
| `phase_set_kind` | `Enum{ps, psl, cgi_haplink}` | never null |
| `phase_set` | `str` | the chromosome's implicit set (phased, no PS) |
| `haplotype` | `u32` | a PSL set or a CGI hapLink, each one haplotype itself (a hapLink block keeps its `chrom`). For a PS set, the GT slot: slot `h` of every record in the set is one molecule |
| `sites` | `u32` | never null. Records whose slot is in the set on this haplotype |
| `members` | `list[str]` | the block is withheld. Otherwise the member Alleles' VRS ids, sorted, each once |
| `cpb_id` | `str` | the block is withheld. Otherwise `ga4gh:CPB.` + the sha512t24u of `{"members":[<Allele digests, sorted>],"type":"CisPhasedBlock"}`, as vrs-python identifies a block of inline Alleles |
| `withheld` | `str` | the block has an id. Otherwise why not (below) |
| `withheld_chrom` | `str` | the block has an id, or is `too_few_members`. Otherwise the contig of the first slot that withheld it |
| `withheld_pos` | `u32` | as `withheld_chrom`: that slot's record position, 1-based |
| `withheld_allele` | `str` | as `withheld_chrom`, and null for an uncalled slot: that slot's allele as `gt_alleles` holds it |

- **Members** are the non-reference alleles the haplotype carries. A `0` on a record split one row per
  ALT (RM270) or under an overlapping deletion says "not this ALT", not "the reference base", so REF
  slots join only with `reference_alleles=True`. `*` never joins: its deletion is a member at its own
  record.
- **Withheld, never shortened** (charter 4). A block missing one of the haplotype's alleles would carry
  an id that looks like the haplotype's. So the first slot in row order that cannot be stated gives the
  reason for the whole block: `phase_unknown` (the slot's `phased` is null, as on slot 0 of every VCF
  4.4+ record, P4, and its allele would be a member: a REF slot by default is not), `uncalled`,
  `not_indexed` (no `gt` index), or the allele's own `vrs_withheld` reason (`not_run` before
  `mint_vrs_ids`, `symbolic`, …; `ref_unstated` on every slot of a row whose source REF states no
  base and which spells a sequence ALT, its REF slot included, Decision 2). Then `too_few_members`:
  VRS requires two members. The vocabulary is
  `phasing.CPB_WITHHELD_REASONS`, open like `vrs.WITHHELD_REASONS`. On a whole imputed genome this
  bites: an implicit set is a chromosome, and one symbolic ALT on it withholds the haplotype's block
  (a phased Genes for Good upload: 45 of 46 blocks withheld `symbolic`).
- **A CGI hapLink** is one molecule: one block per id and chromosome, `haplotype` null, whichever
  slot each allele sits in, so `var` and `masterVarBeta` (whose slot orders differ) give the same
  blocks. A linked haplotype CGI wrote `?` withholds its block `uncalled`. A locus with every haplotype
  no-called is a block, not a `calls` row, so whatever links it carried are not in `calls`.
- **Haploid calls** (no separator, so phased in VCF's reading) sit in their chromosome's implicit set
  even in a file nobody phased. On one copy (male X outside the PARs, Y, MT) that is true by ploidy, so
  they make a block. Rejected: leaving out haploid rows, which would drop a true statement to guard
  against a caller that writes a two-copy region haploid.
- **Which allele withheld it**: `withheld_chrom`, `withheld_pos` and `withheld_allele` name the slot
  that gave the reason, the first in row order: its record's contig and position, and the allele
  as `gt_alleles` holds it (null for an uncalled slot). On a whole-chromosome implicit set this is
  how a reader finds the one symbolic ALT that withheld the block. Null when the block has an id or
  is `too_few_members`, which no single allele causes.
- **A PSL record's slots are split by kind**: its named slots are in their PSL sets, its `ps` slots
  (phased, PSL `.`) in the chromosome's implicit set, and only those make the record part of it.
- **Order**: by the set's first row in `calls`, then kind, set name and haplotype.
- An IRI is not an Allele to vrs-python 2.3: a block whose members are written as `ga4gh:VA.…` strings
  digests to another id there. Resolve the members to Alleles before identifying it again.

## Coverage facets: read depth beside the frames (RM278)

A coverage facet says how deep the reads were, interval by interval: what was callable (charter, "a
facet too"). `read_coverage(path)` reads a depth tool's output; `alignment_coverage(path)` runs
polars-bio's `depth()` over a BAM or CRAM, the substrate's own alignment processing, used as is. Both
return a `Coverage`: a `CoverageFacet` model and a lazy `intervals` frame. It is not a fourth frame of
`Loaded` and changes no column of the three, so `SCHEMA_VERSION` stays.

`intervals` (`COVERAGE_INTERVALS`), sorted by contig rank, `pos`, `end_pos`:

| Column | Meaning |
|---|---|
| `chrom` | canonical (Decision 7) |
| `pos`, `end_pos` | 1-based, inclusive, like `calls` and `blocks`. A BED's `start` is `pos - 1`; `samtools depth` and `samtools coverage` positions are taken as written |
| `depth_mean` | mean depth over the interval's described bases; null for a quantized bin |
| `depth_min`, `depth_max` | the lowest and highest depth at any described base. Null where the source bounds nothing (a regions mean); `depth_max` null for an open bin (`4:inf`) |
| `bases` | how many of the interval's bases the source described; fewer than its length means a window only partly described |
| `label` | a quantized bin's or a named region's label, verbatim |

`CoverageFacet`: `facet_id` (a digest of the source, the format, the options and `COVERAGE_VERSION`),
`format_id` (`per_base_bed`, `regions_bed`, `quantized_bed`, `mosdepth_summary`, `samtools_depth`,
`samtools_coverage`, `alignment_depth`: named by shape, since mosdepth's per-base BED and bedtools'
`genomecov -bg` are one shape), `source_path`, `source_sha256`, `assembly_label` and its basis,
`contig_lengths`, `resolution` (`base`, `window`, `bin`, `region`), `window`, `options` (the read
filters `depth()` ran with), `zero_filled`, `dropped_fields`, `warnings`, `subject` (the caller's
label, "Facets of one person" above). just-vcf's own shape
(charter 2): GA4GH has no coverage model.

Rules:

- **A position no interval covers is unknown, never zero.** A file may be a region cut, and `samtools
  depth` leaves zero-depth bases out unless run with `-a`. The exception is `alignment_depth`, whose
  BAM header names every contig: its windows tile them and an uncovered base is zero
  (`zero_filled = True`).
- **Per-base sources become windows while they stream** (`window`, default 1,000): the folded sums,
  bounds and described bases are exact, and a whole genome is about 3.1 million rows whatever its depth.
  `window=None` keeps every source row, read whole, for small files. Regions, bins and summaries keep
  their own intervals.
- **The kind is read from the content**, compression from the magic bytes. A `samtools depth` with
  several depth columns is refused (one alignment file per facet), as are mosdepth's thresholds file and
  anything whose rows are not positions. A headerless three-column file whose third column always
  exceeds the second is refused too: a BED3 of intervals and a `samtools depth` near a contig's start
  (chrM, thousands deep) look the same. `kind=` states the format where the content cannot.
- **The build** is stated by the caller, or inferred from contig lengths where the source gives them
  (an alignment header, mosdepth's summary). Otherwise it is null.
- **An indexed alignment file is not hashed**: `source_sha256` is null and `facet_id` digests its
  header, its index and its size. Rejected: the sha256 of the bytes, which would read a 30x BAM a
  second time over the network for a name. **One with no index is hashed** (2026-10-26):
  `source_sha256` is its bytes' sha256, which `facet_id` and the spill key digest beside its header
  and size. Its header and size alone named it, so a BAM rewritten in place with the same header and
  byte size but other reads was served the first one's windows under the first one's `facet_id`. The
  cost is one more read of the file on every call, spill or not; an index is what a pipeline leaves
  beside its BAMs. Rejected: not spilling an unindexed file, which reads it in full by `depth()` on
  every call and still gives two such files one `facet_id`.
- **An alignment file is read once.** `alignment_coverage` spills its windows to the cache
  (`coverage/<key>/`, kind `coverage_spill`, with a manifest, swept and cleared as the loader's spills
  are, docs/ROADMAP_HISTORY.md RM268), and a later call with the same inputs reads the spill. The key
  digests the file's header, index and size (its sha256 too, with no index), `window`, `min_mapping_quality`, `filter_flag`,
  `COVERAGE_VERSION` and the polars-bio version, so a spill made with other parameters, or by another
  `COVERAGE_VERSION` or polars-bio, is never served. `assembly_label` is left out: no window depends on
  it, and the facet built on a reused spill still carries the build the call states. The spill trusts
  the file's identity exactly as far as `facet_id` does: two indexed files with one header, index and
  byte size would share it. `clear_cache(source=bam)` matches the spill by path, but hashes an existing file to
  look for digests too, which is one full read of the BAM. Rejected: keying by the just-vcf release,
  which would read every BAM again at each upgrade where `COVERAGE_VERSION` already says when the
  windows change.
- **polars-bio 0.36 writes depth as Int16 by a wrapping cast** from a 32-bit count. A negative depth is
  left out of its window (`bases` drops) and counted in a warning; a depth past 65,535 wraps to a small
  positive value that cannot be told apart.

**Callable.** `callable_sites(coverage, sites, min_depth=…)` adds `callable` and `callable_reason` to a
frame of `chrom, pos`; `callable_at` asks for one position and returns a `Callability` (truthy only
when callable). `True` when the covering interval's `depth_min` reaches `min_depth` and every base was
described; `False` when its `depth_max` is under it, likewise; otherwise null with a reason:
`no_interval`, `no_bounds` (a mean or a renamed bin bounds no base), `straddles`, `partly_described`.
`min_depth` has no default: what depth makes a base callable is the caller's choice. It is evidence
beside `absence_means` (Decision 4) and changes no answer `absence_is_reference` gives.

Rejected, each with its cost:

- A fourth frame on `Loaded`, filled by `load`: changes the load contract, the sniff registry and the
  cache spill of every format, for evidence that belongs to an alignment, not to a callset. It can
  move there when a person-level layer pairs facets.
- Coverage rows in `blocks`: a block is a caller's reference or no-call statement; reads reaching a
  base are not a call, and `absence_means = blocks` would then read depth as reference.
- Keeping `depth()`'s run-length blocks as they come: 498 million runs on a 27 GB BAM (RM278), where
  windows with exact bounds answer the same questions. `alignment_coverage` reads `depth()`'s
  per-base mode, which peaks lower than its block mode (docs/ROADMAP.md RM278).

## Graded genotype evidence: computed on demand, not stored (layer 7, charter 7)

The frames keep each quantity as the source wrote it (Decision 3). `just_vcf.evidence` reads them into
one graded value a consumer can rank and compare across sources, and stores nothing: it is a function
over `calls` and `blocks`, like `sex_chromosome_evidence`, so no column changes and `SCHEMA_VERSION`
stays. The shape is just-vcf's own (charter 2): VRS 2.x has no genotype model and GA4GH defines no
graded-evidence value. Rejected: a stored column, which would freeze today's rules (the prior, the GQ
reading) into every cache spill and need a version bump each time they move.

`genotype_evidence(loaded)` gives one row per `calls` row and per `blocks` row (`GENOTYPE_EVIDENCE`):
`facet_id`, `chrom`, `pos`, `end_pos`, `evidence_basis` (the origin, Decision 3: `assay` or
`ref_block`), `block_kind`, and

| Column | Meaning |
|---|---|
| `grade_basis` | which quantity the value rests on (below). Never null |
| `gl_source` | on a row graded `gl`, the facet's `gl_source`: `gl`, `pl`, or `gl_over_pl` where the header declares both (the row took its GL, else its PL, and the frame does not keep which). Null on every row of another basis, which rests on neither key |
| `confidence_kind` | `quantified` (a number), `asserted` (a genotype stated without a quantity) or `unknown` (none stated). It follows from the basis, and is a column so that a null `confidence` is never read as a low one |
| `genotype_alleles` | the alleles `genotype_probabilities` and `genotype_index` index: `source_alleles` for a vector (Decision 2: never reordered; for a VCF row that is `[source_ref] + alts`), `[source_ref or ref] + alts` for a genotype, which is what `gt` indexes (`ref` where `source_ref` is the empty REF of an insertion point `resolve_ref` anchored, whose alleles then start with the anchor), and `[null, "<*>"]` for a block, whose REF base the frame does not keep |
| `genotype_probabilities` | `list[f64]`, one probability per genotype in VCF genotype order at the row's ploidy; null where the source supports no vector |
| `prior` | how the vector was made: `source` (GP as given, under whatever prior the source used), `flat` (GL/PL normalized with every genotype equally likely beforehand, just-vcf's choice), `even_remainder` (GQ alone, below) |
| `genotype_index` | the genotype `confidence` speaks for: the top one of a vector (the first on a tie), or the called one (`gt`, haploid or diploid). 0 for reference asserted by a block or by absence |
| `confidence` | phred-like: `-10·log10(1 - p)` of that genotype, from 0 (no confidence; a flat three-genotype vector gives 1.76) up. Null where nothing quantifies it, never 0 for unknown. Infinite where the source gave every other genotype probability zero (a rounded GP `0,1,0`); it is computed in log space, so a PL of 3,000 gives 3,000 and not infinity by rounding |

**The bases**, the first that applies to a row:

- `ref_unstated` (0.0.22, F23): a row whose source REF states no base and which spells a sequence
  ALT (`ref_check = unstated`, Decision 2). Its genotype and any vector were spelled against the `N`,
  so nothing is graded: unknown, no `genotype_index`, alleles, vector or number. A row with only
  symbolic ALTs, or none, is graded by the bases below.
- `gp`: a usable GP vector, as written.
- `gl`: a usable likelihood vector, normalized under a flat prior. GL and PL are one quantity on two
  scales and load into one `gl` column (Decision 3), so there is one basis for both; the row's
  `gl_source` says which key the facet read (`facets.gl_source`). Rejected: a `pl` basis beside `gl`,
  which would split one quantity in two, change what `gl` means to a consumer already filtering on it,
  and have no answer for a facet declaring both keys, whose rows do not keep which key each took.
- `gq`: a fully called genotype carrying an ALT, with a GQ. VCF defines GQ as conditional on the site
  being variant ("Other VCF 4.5 points" above), so by default GQ grades no reference genotype and no
  reference block; `gq_on_reference=True` lets it, as single-sample callers compute it. The called
  genotype gets `1 - 10^(-GQ/10)` and the rest is spread evenly over the other genotypes, which the
  source did not rank (`prior = even_remainder`); `confidence` is GQ. No vector beyond diploid, or for a
  block, whose ploidy the frame does not keep without its likelihoods.
- `hq`: a fully called genotype whose every slot has a haplotype score (`hq`, CGI's allele scores,
  "CGI's allele scores" above), of any ploidy, reference genotypes included: CGI scores a reference
  haplotype's call as it does an ALT's, so the `gq` basis's restriction to genotypes carrying an ALT
  does not apply. `confidence` is the lowest score, since the genotype is wrong where any haplotype is
  (CGI's own VCF writes GQ so, above). No vector and no `prior`: one number per haplotype ranks no
  other genotype, and Decision 3 rejects a made-up distribution. `genotype_alleles` and
  `genotype_index` name the called genotype. A half-call or a partly called slot has no score, so the
  row stays `gt_only`. A VCF's FORMAT HQ ("A VCF's haplotype qualities" above) grades the same way,
  and keeps its place in this order: a VCF row with usable GP or GL is graded by them and one with a
  GQ on a genotype carrying an ALT by GQ, whatever its HQ; HQ grades the rest, which by default
  includes a reference genotype with a GQ (`gq` grades none, unless `gq_on_reference=True`, which
  puts GQ first there too). Rejected: HQ before GQ, which would rank a per-haplotype score above the
  caller's own score of the genotype it called.
- `dosage`: only a dosage. One number gives no vector (Decision 3), so no vector and no number.
- `gt_only`: a genotype with at least one allele called and no quantity: asserted, unquantified. A
  half-call (`0/.`) is here even with a GQ: there is no one genotype for GQ to grade.
- `block`: a `ref` block with no quantity the block rule reads (CGI's native blocks, or a block lacking
  what `block_grade` asks for): reference asserted, unquantified.
- `no_call`: the source wrote that it looked and called nothing (`./.`, a `no_call` block) and gave no
  vector: unknown. A no-call that carries likelihoods is graded by them (a gVCF `./.` block with PL
  `0,0,0` is `gl` at 1.76).
- `absent_reference` and `absent_unknown`: only from `genotype_evidence_at` (below).

A vector is usable when it has one value per genotype at the row's ploidy (GT's, else diploid as the
spec says, ploidy 1 to 4), at least two genotypes, no missing value, and for GP no negative value and a
positive sum. Otherwise the row falls to the next basis; nothing is patched.

**Asked positions.** `genotype_evidence_at(loaded, sites, variant_class=…)` answers a one-facet
`loaded` at each `chrom, pos` (`GENOTYPE_EVIDENCE_AT`): one row per record that covers the site (a
`calls` row from `pos` to `end_pos`, or its `pos` alone where the span is unknown; a `blocks` row over
its span), graded as above, with the record's span as `record_pos`, `record_end_pos`. A site inside a
called deletion is not absent. Every block over a site is found, overlapping ones included: a
gVCF's and CGI's blocks tile the genome (none of 10.8M CGI `var` blocks or 14.6M DeepVariant blocks
overlaps its predecessor), so the block over a site is the last one starting at or before it, an
as-of join; not every gVCF's quite do (927 of the 9.24M blocks on the primary contigs of PGP-HMS
hu416394's bcftools gVCF reach the next block's start, measured 2026-10-26), and a lifted facet's
do not (two source regions a chain sends to one place give overlapping blocks, "Lifting blocks"),
and a block that reaches its successor goes through a range join instead (`intervals.covering`). A site's rows come by record start, origin and end, then by what each states
(block kind, basis, confidence, alleles, genotype), so records over one span come out in one order
whatever the file's. Only the matched records are graded, so a panel's answer does not grade a
genome's every row. A site nothing covers gets one row, `evidence_basis = absence`:
`absent_reference` (asserted, no number, `genotype_index` 0) where `absence_is_reference(facet,
variant_class, force=force)` is True, else `absent_unknown`; `force` lifts a content filter's hold for
the class it kept, as there. Decision 3 rejected a made-up number for an
unquantified call, and that holds for absence too. `absence_means` comes along on every row, so a
CGI `reference` silence and a `weak_reference` one stay apart. `coverage=` with `min_depth=` adds
`callable`, `callable_reason` and the depths of `callable_sites` beside the answer; they change no
basis, because reads reaching a base are not a call (coverage facets, above).

Measured 2026-10-23 on PGP-HMS facets, each asked at the packaged relatedness panel's 15,289 SNVs:

| Facet | Bases of its rows | Panel sites |
|---|---|---|
| 23andMe, AncestryDNA, Living DNA, Genes for Good (array and its 23andMe-layout imputed export, which writes no probabilities) | `gt_only`, `no_call` 0.0–1.9% | 26–51% `gt_only`, rest `absent_unknown` |
| Genos exome VCF (GATK, gVCF-banded) | `gl` 44,784 (confidence p10/p50/p90 102/581/2,473), `no_call` blocks 259,070 | 46% `gl`, 25% `no_call` block, 29% `absent_unknown` |
| CGI `var` (whole genome, hu1B2689) | since 0.0.19 `hq` 4.02M (68/374/1,037 at p10/p50/p90), `gt_only` 205,177 (half-calls and partly called slots), `no_call` 356k calls; `block` 7.7M, `no_call` 3.1M blocks. Before 0.0.19 every call row was `gt_only` (4.2M) | 60% `gt_only` before 0.0.19 (not asked again since), 38% `block`, 2% `no_call`; none absent |
| CGI VCF of the same genome | `gl` 4.16M (3/357/1,022 at p1/p50/p90), `gt_only` 48,655, `no_call` blocks | 60% `gl`, **5,763 `absent_reference`: the same 5,763 sites its `var` file calls with a `ref` block** |
| DRAGEN `.snp.` half (GP) | `gp` 3.86M (41/48/175) | 62% `gp`, 38% `absent_unknown` (its `content_filter` holds absence back unless `force`) |
| DRAGEN `.indel.` half | `gp` 956,057 | all `absent_unknown`: an SNV question, an indel-only cut |
| Isaac/Starling GRCh38 VCF | `gq` 4.63M (45/137/235) | 62% `gq`, 38% `absent_reference` |
| PGP freeBayes + VarScan2, the two caller facets | freeBayes `gl` 4.85M (42/97/97), `gt_only` 114,598; VarScan2 `gq` 340,150 (16/42/141) | freeBayes 59% `gl`, 40% `absent_reference`; VarScan2 2% `gq`, rest `absent_unknown` |

No PGP facet carried both GQ and likelihoods on a reference genotype, so the GQ reading was checked on
a local DeepVariant gVCF (aggregates only): on 1.32M hom-ref calls GQ is within 1.5 of the confidence
their PL gives (clipped at 99, the highest GQ DeepVariant wrote on them) for 96.6%, which is what `gq_on_reference=True`
assumes. On its 13.6M reference blocks only 26.9% agree, median difference 49: a block's GQ is the
band's floor (at most 50 there), its PL is not, so a block graded by its likelihoods (median 291) states more than its
weakest base supports. Which of the two should grade a block is the person's call (desk #94); until
then a block is graded by its likelihoods, as written.

**The block rule** (`block_grade=`, `BlockGrade`, on `genotype_evidence`, `genotype_evidence_at` and
`consensus`) is built for that answer to be a one-line default change. It governs `ref` blocks only;
a `no_call` block is graded by its likelihoods, or is `no_call`, under every rule.

| Rule | A `ref` block is graded | Lacking what it reads |
|---|---|---|
| `likelihoods` (the default) | by its likelihoods (`gl`); without usable ones, by its GQ where `gq_on_reference=True` | `block` |
| `gq` | by its GQ (`gq`, `prior = even_remainder`), whatever `gq_on_reference` says; the vector is spread over as many genotypes as its likelihoods give where it has usable ones, else none | no GQ: `block` |
| `lower_of_gq_and_likelihoods` | by whichever states less confidence in hom-ref: GQ against `-10·log10(1 − p_0)` of its likelihoods under a flat prior, the likelihoods on a tie; `grade_basis` says which (`gq` or `gl`) | either missing: `block` |

Genotype 0 is compared, not the likelihoods' top genotype: a block's GQ speaks for hom-ref, and where
its likelihoods put another genotype first they state less confidence in hom-ref than any GQ but a few
phred. A block lacking what its rule reads is `block`, reference asserted with the number withheld,
never graded by another rule's quantity (CLAUDE.md, no silent fallback).

**Reference-supported positions (RM3).** `reference_support(loaded, sites, variant_class=…)`
materialises what a one-facet `loaded` says about the reference allele at each asked position
(`chrom` canonical, `pos` 1-based), from calls, blocks and absence at once: one row per distinct
position (`REFERENCE_SUPPORT`), in canonical order (contig rank, then position) whatever order the
sites came in, as a lazy frame. It takes `gq_on_reference`, `block_grade` and `force` as
`genotype_evidence_at` does, and grades a covering record exactly as `genotype_evidence` does: the
grade columns are that frame's (`grade_basis`, `gl_source`, `confidence_kind`, `genotype_alleles`,
`genotype_probabilities`, `prior`, `genotype_index`, `confidence`). What it adds is the one answer per
position, what that answer rests on, and why:

| `support` | When | Grade | `reason` |
|---|---|---|---|
| `call` | a `calls` row covers the position (from `pos` to `end_pos`, both inclusive) | as `genotype_evidence`: `genotype_index` 0 states hom-ref, any other an ALT on some haplotype, null none; on a row with a `split_signal` a `0` may be another ALT (Decision 2), so `rejoin_split` first | `call_at_position`, `call_spans_position` |
| `block` | a `blocks` row covers it (both ends inclusive) | a `ref` block by `block_grade`, `genotype_index` 0; a `no_call` block by its likelihoods, else `no_call` | `ref_block`, `no_call_block` |
| `absence` | nothing covers it, and `absence_support` reads the silence as reference | `absent_reference`, asserted, `genotype_index` 0, no number | `whole_genome_silence`, `omitted_called_reference`, `content_filter_forced` |
| `none` | nothing covers it, and the silence says nothing | `absent_unknown`, unknown, no genotype | the `AbsenceReason` (Decision 4): `outside_blocks`, `array`, `targeted`, `imputed`, `callset_type_unknown`, `absence_not_carried`, `class_not_counted`, `class_not_called`, `content_filter_held`, `class_cut_by_filter` |

`absence_means` comes along on every row, so a producer's reference call by omission (`reference`) and
a caller's silence (`weak_reference`) stay apart on an `absence` row. That is the absence's own grade,
below every record's: Decision 3 makes up no number for an unquantified statement, and a silence is
one. A number would have to come from the method's own false-negative rate per region (RM279's atlas),
not from a constant; until the person asks for one, none is given. `records` counts the records
covering the position (0 for `absence` and `none`). Where several do, one answers: a call starting at
the position, else the call starting nearest before it, else the block starting nearest before it
(a lifted facet's blocks can overlap), the shortest first, a `no_call` block before a `ref` one; among calls of one span, one stating an ALT before one
stating nothing before one stating hom-ref, then a quantified one before an unquantified one, the
least `confidence` first, then by allele list (and, past that, the genotype, the basis and the
vector), so no tie depends on the order a join emits rows. A GATK gVCF does start a
`ref` block inside a deletion's span, and a reference statement must not hide a variant stated at the
same base. Rejected: one row per covering record, as `genotype_evidence_at` gives, which leaves every
consumer restoring hom-ref to pick among them, each by its own rule. The asked positions and the facet
row are read up front; every join over `calls` and `blocks` is pre-filtered to the asked contigs, and
only the covering records are graded (`genotype_evidence_at`'s matching, shared).

## The cache under the lazy frames (RM268)

The frames `load` and the derived steps return are lazy scans over files in `cache_dir()`
(`JUST_VCF_CACHE_DIR`). Not part of the schema, and no version bump; the rules are here because a
frame stops working when its files go.

- Entries are directories under `vcf/` (rewrite copies, kind `rewrite_copy`), `keyvalue/` and `cgi/`
  (the decompressed payload the key=value and Complete Genomics loaders read while they spill, kind
  `work_copy`, emptied when the load ends, failed or not), `frames/` (spills: `frames_spill`,
  `rejoin_spill`, `resolve_spill`, `vrs_spill`, `normalize_spill`, `lift_spill`), `coverage/`
  (`coverage_spill`) and `stats/` (a VCF load's warning counts, `scan_stats`, below). Beside each is
  `<name>.json`, its manifest, and `<name>.lock`, an OS file lock the writer holds while it writes.
  The manifest names the root source (sha256 and resolved path), which every derived facet carries
  down from its parent, so a genome's derived spills go with its copies.
- An entry of an older version is removed when the step that owns its kind next runs (`load` for
  copies and load spills). Older means no version key ahead of the running version's and at least
  one behind, compared numerically (`0.0.9` is behind `0.0.11`); a key the running version has and
  the entry lacks counts as behind (the entry predates the key). An entry with any key ahead, a key
  the running version does not know, or a value that is not a number is a newer install's, or may
  be: kept. Two installs of different versions share the platform user cache, and the older one must
  not remove the newer one's live spills; the newer one retires the older one's. Rejected: sweeping
  every version other than the running one, which made two installs delete each other's spills on
  every load.
- The older install may still be running, with frames over its entries: two sessions on one cache
  dir at two `SCHEMA_VERSION`s hit this on 2026-10-23, when the newer one's sweep removed a generation
  the older one's live frame read and its next `collect` failed (`FileNotFoundError` on
  `calls-00000.parquet`). So an older entry goes only when (a) no process holds its lock, and (b) its
  `last_used` is older than the grace period, `JUST_VCF_CACHE_GRACE_DAYS` (7 days unset; `0` sweeps
  at once, as before; a value that is not a non-negative number is refused, not defaulted). Both are
  judged again under the entry's lock, so a step that ran between the sweep's scan and its removal
  keeps its entry. A step that writes or reuses a generation moves `last_used` (installs from before
  this rule too), so a frame stays safe for the grace period after its step last ran. Nothing marks a
  `collect`: a frame kept open longer than that, whose step never ran again, can still lose its
  files. An entry with no manifest takes the newest modification time in it as its last use. The
  sweep logs how many older entries it kept, and `CacheClearReport.recently_used` lists them. A
  writer that starts waiting on the lock while the sweep holds it finds the lock file gone and writes
  the entry afresh.
  An entry with no manifest predates the index: removed only by that sweep, where its name shows an old
  version or where it must be old (a `vcf/` entry, below `REWRITE_VERSION` 7; a `keyvalue/` or `cgi/`
  one, from before those copies were indexed).
- A manifest this version cannot read (a newer version's, with a kind or field it does not know) makes
  its entry unknown origin: kept, never swept, and it does not stop `load`. Kinds are checked where
  an entry is written.
- A spill entry, `load`'s and each derived step's (`rejoin_split`, `resolve_ref`, `normalize`, `mint_vrs_ids`, `liftover`), holds one
  directory per generation. A run writes into a staging directory inside the entry, under the entry's
  lock, and renames it to `<key>/` when it is complete, `<key>` being a sha256 of the content the
  result was computed from (the staged input's bytes, or for `rejoin_split` and `load` the result's,
  and the facet row) and of the generation's shape (the step's versions, as its manifest records them,
  and the set of file groups it holds), and `result.json` in it listing its files by group and the
  step's row warnings. A run whose key is
  already there reuses that generation and writes nothing; a run on other rows under the same
  `facet_id` (a filtered `Loaded`) commits a second generation beside the first. No run rewrites or
  removes a complete generation, so a frame an earlier run returned, in this process or another,
  keeps working. A staging directory left by a writer that died, and a generation missing a file or
  holding other groups than its step's, are removed by the next run. A step that adds or drops a group,
  or changes a version, so writes a generation of its own beside the old one and never reads the old
  shape: a lift generation 45126a2 wrote, without the `mate_contigs` group 0d3599a added at the same
  `LIFT_VERSION`, made 0d3599a's `liftover` fail with a `KeyError` while the key covered the input
  alone. Rejected: keying reuse on the facet id and versions alone, which would
  hand a filtered `Loaded` the full facet's rows; reading a missing group as empty, which would hand
  back a `contig_map` silently missing the mates' contigs; and replacing the files in place, which broke the
  earlier frame (`Invalid thrift` on its next `collect`). The manifest's `SPILL_LAYOUT` key (2) retires,
  through the sweep, entries written before generations. A frame's files are parts of whole contigs in
  canonical order, consecutive contigs grouped up to `loading.SPILL_PART_ROWS` rows
  (`loading.canonical_batches`), so the frame is one scan with no sort however many contigs it has:
  `load`'s and `rejoin_split`'s (F27), `resolve_ref`'s, `normalize`'s and `liftover`'s (F24; a lift's
  report in parts of whole source contigs, in the parent's order within each). Each of these steps
  names its layout in its versions (`FRAMES_LAYOUT`, `REJOIN_LAYOUT`, `RESOLVE_LAYOUT`,
  `NORMALIZE_LAYOUT`, `LIFT_LAYOUT`), so a generation an older version wrote with a file per contig is
  never read, and the sweep retires its entry. `mint_vrs_ids`'s generation holds its staged input
  (`calls`) as it came and the `ids` beside it by position, one part per group of consecutive runs of
  one contig in row order, by the same bound (`loading.bounded_groups`); `VRS_LAYOUT` names that layout,
  so a generation written a pass per run, its id parts cut every million rows, is not read (it would
  read the same, but is written once more under the key). `liftover`'s generation holds five groups,
  the lifted `calls`, the `report`, the lifted `blocks`, the `block_report` and the `mate_contigs`
  lifted breakends name; its key is the bytes of
  the staged parent `calls` and `blocks` (the blocks since `LIFT_VERSION` 8, when they began to be
  lifted; before, only their count), and takes the chain's file name and the two collections' labels,
  which its warnings name. `Lifted.reasons` and `Lifted.block_reasons` are read back from the two
  reports, so a reused generation answers as the run that wrote it. Its spill, once
  `frames/<facet_id>-<SCHEMA_VERSION>-lift/` with no manifest, removed and rewritten on every run, is
  now `frames/<facet_id>-<SCHEMA_VERSION>/` like the others; the old name is swept by `liftover`.
- `load` sinks `calls.parquet` and `blocks.parquet` into its stage on every run, as before, and keys
  them by their own bytes: the same file sinks to the same bytes, so a reload, in this process or
  another, finds the generation and touches nothing. Measured 2026-10-22 on a PGP-HMS VCF (4,065,482
  calls, 145 MB of spill): two loads in two processes wrote byte-identical files, and hashing them
  takes 0.09 s (`content_key` over the two files, warm page cache) against a load of over a minute.
  A reload never skips the sink. Rejected: reusing a complete spill by
  `facet_id` without sinking, which would hand a load the rows an older loader wrote whenever a fix
  to a loader does not bump `SCHEMA_VERSION`; and writing under a temporary name and renaming over
  the old files, which swaps the files under a frame that is mid-`collect`.
- A rewrite copy under `vcf/` is read only while `load` sinks its spill; no frame `load` returns reads
  it, and a complete copy (its `.complete` marker present) is never written again, so it needs no
  generations.
- A VCF load's warning counts and the facet fields they feed (`contig_map`, `absence_means`, the GP
  scale, the ploidy) come from one polars-bio pass over the records (`loaders.vcf.scan_stats`), separate
  from the sink. Its result is stored as `stats/<key>/stats.json` (kind `scan_stats`), so a reload reads
  the records once, to sink them (F26). `<key>` is a sha256 of everything the pass depends on
  (`loaders.vcf.stats_key`): the bytes polars-bio reads (the copy's, hashed on every load, or the
  source's when it is read as it is) and the file's name; the source's digest and member, which fix
  the header the loader parsed; the keys and sample column asked for, the registration options, the
  header's VCF version, the keys left out of the copy and those the retry left out; the part (a
  merge's caller index, `SF` selection); `READER_DEFINITIONS`; `STATS_VERSION`; the polars and
  polars-bio versions; and a digest of the source of `loaders/vcf.py` and `loaders/text.py`, so an edit
  that forgets the `STATS_VERSION` bump costs a pass, never a stale count. The manifest carries
  `STATS_VERSION`, `POLARS_BIO_VERSION`, `POLARS_VERSION` and that code digest, `STATS_CODE`, and
  `load` sweeps older ones like its other kinds. A digest has no order, so an entry whose code digest
  differs is swept unless one of the three numbered versions is ahead (`loaders.vcf.STATS_DIGESTS`);
  that may sweep a newer install's entry of the same numbered versions, and costs it one pass. The entry is read and written under its lock, through the pass: a second session
  loading the same file waits, then reads the first one's result; the file is written under a
  temporary name and moved into place. A file that does not read whole, or answers another key, is
  counted again and rewritten. Only a pass that succeeded is stored: a read polars-bio refuses is met
  again on every load, and retried as the first load did. A stored result means the spill's read of
  the same bytes, under the same request and polars-bio, is the read that already succeeded, so the
  sink, still under `_damage_named`, is now where a cached load would meet a damaged block. Rejected
  (F26): keying the load's generation by its inputs instead of its staged bytes, which would make a
  reload a lookup but hand back the rows an older loader wrote whenever a fix forgets its bump.
- `clear_cache` removes what the caller names, even under a live `Loaded`; that frame's next
  `collect` then fails. It skips an entry whose lock is held, and touches nothing outside `vcf/`,
  `keyvalue/`, `cgi/`, `frames/`, `coverage/` and `stats/`: the reference store under `reference/` is
  not a cache of a genome. A source's cached stats carry its digest and go with its copies.
- Rejected: a size cap with least-recently-used eviction. `last_used` is when a step wrote or reused an
  entry, not when a frame last read it, so a cap could evict under a live frame without anyone asking.

## The durable store (RM276)

The person chose both (Part 2 interview, E13): a directory per person holding parquet per facet, the
links and a manifest, read lazily; and a VCF Zarr export of one facet's calls (below). The store
adds no column and needs no version bump; it holds the four frames as they are. just-vcf's own
layout (charter 2): GA4GH defines no on-disk store for a person's callsets.

```
<dir>/manifest.json                    store_layout, schema_version, subject, created; per facet its
                                       facet_id, added and the warnings its load gave
<dir>/facets/<facet_id>/facets.parquet, calls.parquet, blocks.parquet, links.parquet
<dir>/.lock                            the OS file lock every add and migration holds (just_vcf.cache's)
```

The manifest also lists the store's migrations (`migrations`, below), empty until one runs.

The manifest's warnings are redacted as `query` returns them (`FacetRedactor`: the facets' sample,
source path and member, and anything shaped like a path): it is the store's one plain-text file, and
a load warning may quote a sample or a file name.

- **Where.** `Store.open(path, subject=...)` creates the store in the directory the caller names,
  absent or empty, and never anywhere else: a store holds whole genomes. Opening an existing one
  with `subject` given checks it is the store's.
- **One person, one schema version.** The subject is named at creation. `add` refuses a facet with
  no subject or another one (`SubjectError`, naming facet ids, never the subject), and a `Loaded` of
  another `SCHEMA_VERSION` (`StoreVersionError`). A store of another version opens, so its manifest
  can be read, but its frames are not read and nothing is added: versions are never mixed, and nothing
  is migrated on read. The refusal names both versions and what to do: `migrate_store(path)` for an
  older store ("Store migration", below), upgrading just-vcf for a later one. A later `store_layout`
  is refused.
- **Facets accumulate (charter 5).** A facet's directory is written once and never rewritten or
  removed, except by a migration, which rewrites every facet's files with the same facets and rows;
  adding a facet the store holds, by `facet_id`, is a no-op that writes nothing, when its REF state is
  the held one's (`store.REF_STATE`: `ref_seqcol_digest`, `resolve_version` and `resolve_inputs`). One
  in another state (held as loaded and given resolved, or the reverse, or read by another
  `RESOLVE_VERSION`, from another collection or with other chip manifests) is refused with `StoreConflictError`, naming the facet and both states, before
  anything is written: `resolve_ref` keeps the `facet_id`, and a `present` there kept the held rows
  while the caller thought the new ones stored. A held facet's files are never rewritten, so the other
  state goes into a new store. Rejected: replacing the held facet (charter 5; a replace is the
  person's call, not an add's), and comparing the rows (a full read of both on every re-add). The links
  are kept per facet, so adding one never touches a file another facet owns: each facet keeps the
  links of the `Loaded` it came with that end at it or at one of its lineage ancestors, so a facet
  added without its parent still names it. `loaded()` reads them all back once each, in Decision 12's
  order. Rejected: one `links.parquet` for the store, which every add would rewrite.
- **Atomic per facet, safe under concurrent adds.** An add holds the store's lock from its check of
  what is present to its last manifest write: each facet is sunk into a staging directory, synced,
  renamed to `facets/<facet_id>/`, and committed by rewriting the manifest (temporary name, then
  `os.replace`). Since the lock spans the whole add, a staging directory or a facet directory the
  manifest does not name is a dead writer's when found under the lock: the next add removes it.
  `open` reports one without taking the lock, so it may be an add still running. Rejected:
  sinking outside the lock and locking only the commit, which lets concurrent adds overlap but leaves
  no safe moment to remove a dead writer's staging directory.
- **Lazy reads.** `loaded()` returns one `Loaded` over every facet, scans of its files in `facet_id`
  order (the canonical order, Decision 11), the links each once, and the facets' warnings each once,
  then one line per set of values a migration filled unknown (below); an empty store gives empty frames. `person(check=...)` is `Person.of` over it; an empty store has
  no person (`StoreError`). A frame read from a store keeps working after later adds.
- **Deterministic (charter 12).** The same `Loaded` sinks to the same parquet bytes under one polars
  release. The manifest holds times and differs between two stores of one input.

### Store migration

`migrate_store(path, keep_previous=False)` rewrites a store written under an older `SCHEMA_VERSION`
to this one and returns a `StoreMigration`: the steps applied, the facets rewritten, and every value
it filled unknown (`FilledUnknown`: the step, the frame, the columns, the facet ids, the reason). The
same record goes into the manifest's `migrations`, and `loaded()` repeats each `FilledUnknown` as a
warning, so a null a migration left is never read as one a load wrote. Nothing migrates on read.

The steps are `just_vcf.store_migration.MIGRATIONS`, one per version bump since the store shipped
under 0.0.14 (merge 6830fb7, between 0.0.14's merge and 0.0.15's), each a pure function per frame,
applied in order from the store's version. A bump that changed no stored frame would get an identity
step; none since 0.0.14 is one. `blocks` and `links` have not changed since the store shipped, so no
step touches them. Each step states its dtypes as they were at its version, and the suite checks the
chain runs from 0.0.14 to `SCHEMA_VERSION` and ends in exactly `schema`'s frames, so a bump without
a step fails it.

| Step | Frames | What it does | Filled unknown, on which facets |
|---|---|---|---|
| 0.0.14 → 0.0.15 | calls | adds `phase_set_ordinals`, `phase_set_qualities` (null) after `phase_set_kinds` | facets whose `dropped_fields` list FORMAT `PSL`, `PSO` or `PSQ`: their PSL sets were left out, so `phased`, `phase_sets`, `phase_set_kinds` lack them and the two new columns are null for not read |
| 0.0.15 → 0.0.16 | facets | adds `reference`, `aligner`, `aligner_version`, `method_steps` (null) after `phasing_method` | VCF facets (`format_id` a VCF loader's): `method_steps` null means not read, not "not a VCF"; the other three null mean not read; `phasing_method` keeps the earlier rule's value |
| 0.0.16 → 0.0.17 | calls | casts `phase_set_kinds` to the enum with `cgi_haplink` | CGI `var` and masterVarBeta facets: `phased`, `phase_sets`, `phase_set_kinds` null for hapLink not read; masterVarBeta facets: `ad`, `source_alleles` for read counts not read |
| 0.0.17 → 0.0.18 | facets | adds `gl_source` (null) after `split_marker` | VCF facets: null means not recorded, not "no likelihoods" |
| 0.0.18 → 0.0.19 | calls | adds `hq` (null) after `gq` | CGI native facets, and facets whose `dropped_fields` list FORMAT `HQ`: `hq` null for not read, so their rows grade `gt_only` |
| 0.0.19 → 0.0.20 | facets | adds `resolve_version` (null) after `ref_seqcol_digest` | facets whose REF was read (`ref_seqcol_digest` set): null means not recorded, not "REF not read", so the store refuses the same facet resolved today as another state of it |
| 0.0.20 → 0.0.21 | facets | adds `resolve_inputs` (null) after `resolve_version` | array facets (`family = genotype_table`) whose REF was read: null means not recorded, not "no manifest", so the store refuses the same facet resolved today with manifests as another state of it; consensus facets: null means the parents' REF states were not recorded, so the store refuses the same consensus built today |
| 0.0.21 → 0.0.22 | calls | casts `ref_check` to the enum with `unstated` | facets whose REF was read and that hold a row with a source REF of only `N` and `ref_check = mismatch`: kept as `mismatch` where today's rule reads `unstated`, so `liftover` leaves those rows behind; the facet records `RESOLVE_VERSION` 12, so the store refuses it resolved today. A caveat may name rows as well as facets (`Caveat.rows`), and this is the first that does |

A step never invents a value: an added column is null, and where that null (or a value an older loader
left) reads differently from what a load of the same source gives today, the facets it applies to are
named, with the reason. Loading the source again is the only way to get those values. Rejected:
deriving `gl_source` from whether a facet's rows have `gl` (it would still guess GL against PL), and
writing `method_steps = []` for a VCF facet (`[]` says the header was read and recorded no step).
A value changed by a rule under one version ("No version change" below) is not a bump, so it gets no
step: a store keeps the values its load gave until the source is loaded again.

- **Refused, with nothing written** (`StoreVersionError`): a store of a later version (no downgrade:
  upgrade just-vcf), and one of a version no store was written under (before 0.0.14, or not a version).
  A store at this version is left as it is: nothing applied, `migrated` null.
- **Written beside, checked, then swapped.** Under the store's lock, every facet is written into
  `.migrating-<pid>-<random>/facets/<facet_id>/`, with the same sink as `add`, and checked against its
  old files: this version's schema exactly, the same row count, and the old columns' values unchanged
  (a digest of every row, the new files cast back to the old dtypes). Only then are the new manifest
  and a `READY` marker written. The swap: the old manifest is copied and the old `facets/` moved to
  `.previous-<old version>/`, the new `facets/` moved in, and the new manifest moved over the old, which
  commits it; then `.previous-*` is removed (kept with `keep_previous=True`) and the migration's
  directory last. A reader meanwhile reads the manifest's version, which is the old one until the
  commit, and is refused as before.
- **Cut off, finished or undone on the next call.** A migration directory without `READY` is a
  migration cut off before every facet was checked; the store is as it was, and the next
  `migrate_store` removes it. One with `READY` is finished by the next `migrate_store`, each move
  skipped where it is already done; one cut off after its commit, with only the cleanup left, is
  cleaned up and its record, read from the committed manifest's `migrations`, is what that call
  returns, so the call does not read as nothing applied. `open` reports a `.migrating-*` or `.previous-*` directory with a
  warning. Rejected: renaming the old `facets/` aside and relying on a crash never landing between the
  two renames; parking old facet directories under `facets/`, where an add would take them for a dead
  writer's and remove them.

### VCF Zarr export

`write_vcf_zarr(loaded, path)` writes one facet's `calls` as a VCF Zarr 0.4 store (Zarr format 2),
the version polars-bio's `scan_vcf_zarr` reads. It is an export at the edge, beside `write_vcf`, and
follows its rules for the rows: REF is `ref`, else `source_ref`, never the genotype; a row VCF
cannot hold is refused, counted per reason, unless `resolved_only=True` (an array facet whose REF was
not read writes no row until `resolve_ref` has run). The path must not exist; the store is written
beside it and renamed into place.

- **Blocks are not representable.** The VCF Zarr spec has no reference or no-call blocks, so a
  facet's `blocks` stay in the parquet store: counted in `blocks_not_written`, in a warning and in
  the root's attributes, never written as calls.
- **Arrays.** `contig_id` (and `contig_length` where the header gave lengths), `variant_contig`,
  `variant_position` (1-based), `variant_length` (`end_pos - pos + 1`, so a symbolic ALT keeps its
  span), `variant_allele` (REF then ALTs), `variant_id` with `variant_id_mask`, `variant_quality`,
  `filter_id` (`PASS` first, then the others sorted) and `variant_filter`; `sample_id`,
  `call_genotype` and `call_genotype_phased` (true where every separator after the first slot is
  `|`, as `write_vcf` writes GT; false, counted in a warning, on a row phased by a CGI hapLink or with
  a slot in a VCF 4.5 PSL set: the spec has one phase flag per call, read as PS's set or the
  chromosome's implicit one, and no field for either, so `true` would put the row in the wrong set); `call_HQ`
  where any row has `hq`, on `call_genotype`'s `ploidy` dimension (VCF 4.4's `Number=P` reading of HQ,
  one value per GT allele): `-1` on a slot with no score, `-2` past the row's ploidy, so a haploid
  `[51]` is stored `51, -2` and a row of more than two alleles is written whole, where `write_vcf`'s
  `Number=2` HQ cannot hold it; an `hq` not one per GT slot (no loader writes one) is left out,
  counted; and, where any row has them, `call_GQ`, `call_DP`, `call_AD`,
  `call_GL`, `call_GP`, `call_DS`, `call_PS` (PS-kind sets as integers), `call_FT`, `call_CN`,
  `variant_SVTYPE`, `variant_CIPOS`, `variant_CIEND`, `variant_CN` and `variant_JUST_VCF_SPLIT`. AD, GL
  and GP are written only where they are indexed over the written REF and ALTs. `variant_CIPOS` and
  `variant_CIEND` (2026-10-26) hold `cipos` and `ciend` as the frame does, two offsets per ALT in ALT
  order, `<i4` on a dimension of their own (`CIPOS_dim`, `CIEND_dim`; the spec lets a `Number=.` field
  name any). An offset of -1 or -2 is a real value there but the spec's integer missing and fill, which
  a reader can only take as "no value" (the spec does not say so; it fixes the sentinels for every
  integer field): such a row's interval is written missing, counted in a warning. Rejected: writing it
  anyway, which stores a value the reader reads as another; and a float array, which keeps every value
  but declares an Integer field Float. `variant_CN` holds `alt_copy_numbers`, `<f4` on `alt_alleles`
  (the spec's `Number=A` dimension), a fractional value as it is, where `write_vcf`'s VCF 4.3 leaves it
  out. polars-bio reads all three back as lists, sentinels as stored.
- **Breakend mates** (2026-10-26, F32) follow `write_vcf`'s rule ("Lifting breakends" above): the
  mate's contig inside a breakend ALT in `variant_allele` is read canonically and written under the
  `contig_names` style CHROM is written under, and every contig a written row's ALT names is listed
  in `contig_id`, rows or none (no `contig_length` for one the header does not state, written
  missing). Under `source` a mate keeps the source's spelling, and a contig only ALTs name, spelled
  two ways by them (`chr3`, `3`), is refused with a `ValueError` before anything is written, as
  `write_vcf` refuses it, where the store wrote both as held. A mate in angle brackets (`<ctg1>`, a
  contig of the `##assembly` file) is written as held and listed nowhere, and a row `resolved_only`
  leaves out lists no mate. The writer's own helpers do it (`writing._mate_chroms`,
  `writing._contig_names`, `breakends.mate_renamed`), so the two exports cannot drift apart.
- **Missing and padding** are the spec's: `-1` and `-2` for integers, the NaNs `0x7F800001` and
  `0x7F800002` for floats, `"."` and `""` for strings. The last chunk is padded whole, as Zarr v2
  stores it. polars-bio hands the integer sentinels back as they are.
- **The root's attributes** name the facet (`facet_id`, `source_sha256`, `format_id`,
  `derived_from`, the schema version), never its `source_path`. `vcf_meta_information` gives the
  `fileformat` `write_vcf` would write for the facet: VCF 4.5 where a row holds a PSL set, else 4.3.
- **Written with the standard library alone**: JSON metadata, `vlen-utf8` strings, `gzip` chunks
  with no timestamp. So it needs no extra, and two writes of one frame give the same bytes (charter
  12). polars-bio's reader takes `gzip` or no compressor and refuses `zlib` (checked 2026-10-23).
  Rejected: bio2zarr, which converts VCF files, so the frame would go out through `write_vcf` and
  back in through text, with cyvcf2 and numpy as dependencies; and the `zarr` package behind an
  extra, for a format whose metadata is a few JSON files.
- **Cost**: the columns written are held in memory and encoded in Python, chunk by chunk. Measured
  2026-10-23 on a PGP-HMS participant's largest VCF facet: 37,559,140 rows in 378 s, a 454 MB
  store, which polars-bio counts back in 2.4 s. The same person's 11 facets (60 M calls, 27 M
  blocks) make a 1.26 GB parquet store, each add taking seconds after its load.

## Concordance between two facets (RM277)

`concordance(a, b, *, region=None, reference=None)` (`just_vcf.concordance`) compares two facets of
one person and returns a `Concordance`: a `sites` frame, a `summary`, the alleles it left out
(`excluded`, per facet and reason) and its warnings. It is a measurement over two facets, not a facet:
nothing is written into either facet's frames, no winner is picked where they disagree (charter 8),
and no derived facet comes out. RM280's consensus is the place that would be. Its frames are
just-vcf's own shape (charter 2: GA4GH has no concordance model). It adds no column to `facets`,
`calls` or `blocks`, so `SCHEMA_VERSION` is unchanged.

**One assembly.** Both facets need their REF read (`resolve_ref`). Equal `ref_seqcol_digest`s compare
every contig. Two collections compare the contigs whose refget sequence is the same in both (Ensembl
`1` and UCSC `chr1` of GRCh38 are one sequence; a VCF's alt contigs that `resolve_ref` could not
match stay unread either way), and the rest are counted as `contig_not_shared`. A primary contig
(1-22, X, Y) whose sequence differs, or no primary contig shared, refuses the pair
(`ReferenceMismatchError`): lift one facet first. Neither facet read and one equal `assembly_label`
on both is compared with a warning, since a label is not the assembly (charter 3); anything else is
refused. Rejected: comparing by label alone. GRCh37 and hg19 share it and differ in MT.

**The unit is one ALT allele**, Decision 1's per-allele view. Sequence alleles meet as below;
symbolic structural-variant alleles meet by type and extent ("Structural variants", below, RM282),
breakends by their adjacency ("Breakends", below). `<*>`, `*`, single breakends (`single_breakend`),
breakend ALTs that do not parse (`breakend_unparsed`), ALTs that are no VCF allele, rows whose REF
was never read, and array rows whose `alts` are still null are counted in `excluded`, and so is a
sequence ALT of a row whose source REF states no base (`ref_unstated`, `ref_check = unstated`,
Decision 2): spelled against the `N`, it reads as neither a match nor a mismatch. A sequence allele meets the other
facet's in this order, each row saying which in `key`:

1. `vrs`: both rows minted an id for it (the ALT's slot of `vrs_allele_ids`, never the REF slot: a
   reference allele's id depends on how long the row's REF is). A VRS id is the fully justified
   allele, so an indel written left-aligned in one file and at the right end of its repeat in the
   other meets here without `normalize` (charter 4).
2. `position`: the rest, by `(chrom, pos, REF, ALT)` as written, `REF` being `ref` or else
   `source_ref`. Unnormalized indels that were not minted can miss each other here; mint both first.
3. `site`: an allele the other facet does not state, where it wrote exactly one row at the
   position over the same REF span: the other genotype is read for copies of this allele (an array's
   hom-ref row, a DeepVariant `RefCall`, a monomorphic record, the other ALT of a multi-allelic
   site). A row over another span is not read: a haplotype caller's complex allele (`AC>GT` beside
   another facet's `A>G`) holds no allele of this span, so its genotype would read as zero copies, a
   mismatch that is only representation. Several rows there and none with the allele (a split file)
   are not read against, counted in a warning. Nor is a row whose source REF states no base and
   which spells a sequence ALT (Decision 2): its haplotypes were spelled against the `N`. Either way
   the allele stays one-sided, `absence = other_alleles`.

An allele stated by several rows of one facet (`normalize` keeps rows that differ, charter 8) pairs
with each, counted in a warning.

**Genotypes are compared as copies of the allele** in `gt_alleles`, unordered: phase sets are per
facet and cannot be compared across facets, so `0|1` against `1|0` is a `match`. Copies count only
this allele, so a split file's `0` for another ALT (bcftools' default) and a `.` (vt, `normalize`)
read correctly: neither is this allele, and `.` is unknown.

Each side's genotype is read against its own row's REF and ALT, so a pair met by VRS id is read in the
row the id came from. A haplotype caller writes a haplotype as one multi-allelic row (REF `CGCCC`, ALT
`TGCCT,TGCCC`, GT `TGCCT/TGCCT`): its `TGCCC` trims to `C>T` and meets another facet's `C>T`, and each
`TGCCT` holds that change too, so counting only the spelling `TGCCC` would read zero copies, a
mismatch that is only representation (F21). A haplotype that is neither REF nor the ALT is therefore
read for the allele's change, both trimmed against REF (common prefix, then suffix):

- a change whose span does not meet the allele's span holds REF there: not a copy;
- of one length with REF and the ALT: a copy when its bases over the allele's span are the ALT's
  (`TGCCT` holds `C>T` at the first base; `CGCCT` does not);
- a length change against another length change or a substitution is another allele of the site
  (a deletion and an insertion of one repeat unit in one row, GT `1/2`): not a copy;
- anything else, a change of another length that also replaces bases over the allele's span, or a
  base other than ACGT there, cannot be read. Such a haplotype is never read as zero copies: it is
  open like a null one, and a pair left undecided by it is `unreadable`, counted in a warning.

The `site` pairing reads the other row the same way, against this allele's REF (the row has the same
REF span). The graded comparison counts copies by the same reading of each of the row's alleles, and
a row with an allele that cannot be read gives no distribution. `agreement`:

| Value | When |
|---|---|
| `match` | both fully called with the same copies; where the ploidies differ, each holds the allele on every haplotype or each on none (a haploid `1` against `1/1`) |
| `mismatch` | both fully called and the copies differ (of two ploidies: not both homozygous alike, so `1` against `0/1`); or with a null haplotype, the called ones already contradict (`1/.` against `0/0`, `1` against `0/.`) |
| `uncalled` | a haplotype is null on either side and the rest do not contradict; a whole no-call (`./.`) is the source saying it could not tell, never a mismatch |
| `unreadable` | a haplotype holds another allele of its row that cannot be read for this one (above) and the rest do not contradict: not decided, and not the source's no-call |
| `only_a`, `only_b` | only one facet states the allele and the other gave no single row to read against |

**Ploidy is a flag beside the agreement.** `ploidy_differs` is True where both genotypes were written
with different numbers of haplotypes, False where with the same, null where either side has no
genotype (a one-sided allele). Outside the PARs one caller writes a male's X call haploid (`1`) and
another diploid (`1/1`): both say every copy of X the person has holds the allele, so they agree, and
the difference is the writer's convention, which the flag keeps (F22). A haploid call against a
heterozygote disagrees: the haploid side says no copy lacks the allele. Rejected: an agreement value
of its own, which a rate over called alleles has to count as either agreeing or not, and which hid
the 1-against-`0/1` disagreements among the convention's matches.

**Graded where possible** (charter 7). Each row's vector, in VCF genotype order over its own
`source_alleles` (Decision 3), becomes a distribution over the copies of the allele: haploid when it
has one entry per allele, diploid when it has n(n+1)/2, summed over the genotypes with each number of
copies, so a multi-allelic row compares with a biallelic one. `p_same = Σ_k p_a(k) · p_b(k)` is the
probability the two facets hold the same number of copies, the facets' errors taken as independent.
Copies are compared as a share of the ploidy, so a haploid `k` meets a diploid `2k`, as the hard
comparison reads them: `p_same = p_a(0) p_b(0) + p_a(1) p_b(2)` for a haploid against a diploid row.
`basis` says what it rests on:

- `gl_flat_prior`: both rows carry likelihoods. Decision 3 forbids turning a likelihood into a
  posterior without a prior, so the prior is named: flat over the genotypes the vector lists
  (`p = 10^gl / Σ 10^gl`, in Float64). Preferred over `gp`: the likelihoods are the assay's evidence,
  and two facets' posteriors can share a population prior that makes them agree for reasons other
  than the reads.
- `gp`: both rows carry posteriors, renormalized to sum to 1 (they are rounded, F10).
- `hard`: anything else, `p_same` null. A likelihood is never compared with a posterior: they are
  different quantities. A vector with a null entry, a length that fits neither ploidy, or an allele
  list that does not hold the ALT exactly once is not used. A `site` pairing is `hard`: the other
  row's vector indexes other alleles.

Rejected: a distance in log space (no probability to read), and comparing the most likely genotypes
only, which is the hard comparison again.

**One-sided alleles** carry `absence`, what the silent facet says, checked in this order: a row there
without the allele (`other_alleles`, also on `site` pairings), a `ref` or `no_call` block over the
position (`ref_block`, `no_call_block`), else the facet's rule: `reference` when
`absence_is_reference(facet, variant_class)` is True, `unknown` otherwise and on a contig the facet
neither declared in its header nor wrote a row on. Never a negative (Decision 4). `absence` is null
where both facets state the allele. A one-sided allele stays one-sided whatever the absence says:
reading `reference` or a block as `0/0` would turn graded absence into a call.

**`sites`**, one row per allele either facet states, in canonical order (contig rank, `pos`, `ref`,
`alt`, then each facet's row order):

| Column | Dtype | Meaning |
|---|---|---|
| `chrom`, `pos`, `ref`, `alt` | `str`, `u32`, `str`, `str` | the allele as the facet stating it wrote it, `a` first |
| `variant_class` | `Enum` | `snv`, `mnv`, `indel` or `sv` (Decision 5's rule) |
| `vrs_allele_id` | `str` | either facet's id for the allele; null where neither minted one |
| `key` | `Enum` | `vrs`, `position`, `site`, `extent` (an SV pair), `breakend` (a pair of breakend adjacencies); null on a one-sided allele |
| `agreement` | `Enum` | above |
| `absence` | `Enum` | above; null where both facets state the allele |
| `basis` | `Enum` | `gl_flat_prior`, `gp`, `hard`; null on a one-sided allele |
| `p_same` | `f64` | null unless `basis` is graded |
| `a_pos`, `b_pos` | `u32` | where each facet wrote the allele (they differ for an indel met by VRS id); on a `site` pairing the other row's position; null on a silent side |
| `a_gt_alleles`, `b_gt_alleles` | `list[str]` | each facet's genotype (on a `site` pairing, the other row's) |
| `a_copies`, `b_copies` | `u8` | copies of the allele; null where a haplotype is uncalled or unreadable |
| `ploidy_differs` | `bool` | the two genotypes were written with different numbers of haplotypes; null on a one-sided allele |
| `a_filters`, `b_filters`, `a_gq`, `b_gq` | `list[str]`, `u16` | each row's FILTER and GQ, so a consumer can restrict to PASS |
| `sv_kind` | `str` | an SV allele's type as matched (below), `BND` for a breakend; null on a sequence allele |
| `b_alt` | `str` | on an `extent` or `breakend` pair, `b`'s ALT as written (`<DUP>` beside `a`'s `<DUP:TANDEM>`); null otherwise |
| `a_end_pos`, `b_end_pos` | `u32` | an SV allele's `end_pos` on each side stating it (a breakend's own position); null on a sequence allele |
| `size_similarity`, `reciprocal_overlap` | `f64` | truvari's measures of an `extent` pair (below); `size_similarity` null where either size is unknown; null on a `breakend` pair |
| `ci_needed` | `bool` | on an `extent` or `breakend` pair, True where the two met only once CIPOS/CIEND widened them; null otherwise |

The last seven columns came with RM282 (2026-10-26), after the rest: added, nothing redefined.

**`summary`** (`summarize(sites)`) counts alleles per `(agreement, ploidy_differs, key, absence,
variant_class, sv_kind, basis)` with the mean `p_same`; its counts sum to the height of `sites`. Run it on a filtered `sites`
(PASS only, alleles either facet carries) to summarize that subset.

**`region`** keeps the alleles whose `pos` lies in one of its intervals: a frame with `chrom`, `pos`
and `end_pos`, 1-based and closed as `blocks` writes them (a BED's 0-based start is `pos - 1`), any
contig spelling. Both facets are cut before they meet, so an indel met by VRS id whose two spellings
straddle an edge counts on each side alone. RM279's atlas does not cut per stratum: it compares once
and places each allele in its strata ("A method-error atlas"). Rejected: carrying the region's own
columns (a stratum name) into `sites`, which overlapping strata would duplicate.

**Structural variants (RM282)** meet by type and extent, as truvari bench matches them (truvari 5.5,
MIT; its wiki page "bench", "Matching Parameters", and the code those defaults run:
`VariantParams.DEFAULTS` and `VariantRecord.var_match`, read 2026-10-26). A symbolic allele
(`<DEL>`, `<DUP:TANDEM>`, `<INV>`, `<INS>`, `<CNV>`, `<CN0>`) is not paired by VRS id or spelling,
since the frame withholds its id (`sv_no_vrs_model`) and two callers rarely write one END. The
thresholds are `SvMatching`'s, named and defaulted as truvari's:

| Parameter | Default | Rule |
|---|---|---|
| `refdist` | 500 | the two extents lie within this many bases: truvari's `overlaps(start - refdist, end + refdist, other start, other end)`, strict, so a call 500 bases past the other's end still meets it, 501 does not. Looser than the wiki's "within from the base call's start/end": the code is followed |
| `pctsize` | 0.7 | `min(size) / max(size)`; a zero size counts as one base, two zeros are equal (truvari's `sizesim`) |
| `pctovl` | 0.0 | shared bases over the longer extent (truvari's `reciprocal_overlap`), so 0.0 takes any pair within `refdist` |
| `typeignore` | False | True pairs across types (truvari's `--typeignore`) |

Coordinates: a call's extent is truvari's half-open `[pos - 1, end_pos)`, the anchor base and the
event, converted from the frame's 1-based closed positions in one place (`concordance._sv_side`).
Size is `end_pos - pos` for the types whose `end_pos` is the event's end (`<DEL>`, `<DUP>`, `<INV>`,
`<CNV>`, `<CNn>`), as `write_vcf` writes SVLEN and as truvari reads a stated SVLEN; truvari's
fallback where SVLEN is missing, `end - start` on its half-open extent, is one base more. Where the record has several ALTs, its `end_pos` is the site's, the largest of
its alleles' ends.

What differs from truvari, each with its reason:

- **Type from the ALT, not SVTYPE.** truvari reads INFO `SVTYPE` first. The frame's `sv_type` is the
  record's claim, and DRAGEN's copy-number caller writes `SVTYPE=CNV` on both its `<DEL>` and its
  `<DUP>`, Canvas `LOH` on `<CN0>,<CN2>`: read first, it would pair a loss with a gain and never meet
  another caller's `<DEL>`. The type is the symbolic ALT's name before any `:`, so `<DUP:TANDEM>` is
  `DUP` (truvari also drops subtypes); a copy-number allele keeps its number (`CN0`), so a loss of
  both copies does not meet a loss of one.
- **No size where the frame keeps none.** truvari takes SVLEN, else `end - start`. The frame keeps
  no SVLEN (above, the person's open call), so an `<INS>`'s length, and any type's whose `end_pos` is
  its anchor (Delly's `<TRA>`), is unknown: `size_similarity` is null, the size threshold is not
  checked, and the pair is counted in a warning. truvari without SVLEN sizes an insertion by its
  `end - start`, one base, so any two insertions are equal in size (similarity 1) on no evidence.
  Such pairs meet by type and distance only, which is weaker than truvari's match with SVLEN. truvari
  also widens each insertion by half its size before the overlap; with no size that cannot be done,
  so under a `pctovl` above 0 two insertions meet only at one anchor.
- **CIPOS and CIEND widen the distance test.** truvari reads them only for breakends (`bnd_match`,
  below). Here each call's extent is widened to the earliest its start and the latest its end can be: the
  start by `|CIPOS[0]|`, the end by `|CIEND[1]|`, or by `|CIPOS[1]|` where CIEND is missing (VCF 4.5:
  assume CIPOS; the frame leaves that to its consumer, which this is). Only the `refdist` test uses
  the widened extents; size and overlap use the stated ones. A `.` widens nothing,
  except that a `.` for the end's upper bound falls back to CIPOS's, as a missing CIEND does.
  `ci_needed` says when a pair met only by the widening.
- **No sequence similarity.** truvari's `pctseq` (0.7) runs only when both calls carry sequence; no
  symbolic allele carries one in the frame, so it is never run, as with `--pctseq 0`. A deletion
  written out base by base is an `indel`, compared allele by allele above, never against a symbolic
  call: matching across the two representations is not built.
- **A symbolic allele never meets a breakend.** truvari's default `decompose` turns a `<DEL>`,
  `<DUP>` or `<INV>` into the breakends it implies to meet another caller's breakends
  (`VariantRecord.cpx_match`). Not built: whether a `<DEL>` and the breakend pair one caller writes
  for it are one call is a judgement the roadmap has not made (RM282, open). Breakends meet
  breakends ("Breakends", below).
- **No size filter.** truvari leaves out calls under `sizemin` (50). A symbolic allele is an SV by
  its writer's statement whatever its length, and every one is compared.

**One pair per allele.** Each allele is taken at most once (truvari's `--pick single`): the passing
pairs are taken best first, by score (size similarity plus reciprocal overlap, truvari's TruScore
without its sequence term; a withheld size counts 0), then the nearer start, then the nearer end, then
each side's row order, which truvari leaves to its sort and is the tie-break that makes the pick
deterministic. A sequence allele stated twice keeps every pairing, since its pairing is identity; an
extent match is a judgement of nearness, and taking a call twice would count one event as two.

**Genotypes of a pair** compare as copies of each side's own ALT, with the agreement, ploidy flag and
graded `p_same` above. A haplotype holding another allele of its record (`<CN2>` beside `<CN0>`, GT
`1/2`) is not a copy: a symbolic allele has no bases for the sequence reading above to trim.

**A one-sided SV** reads the silent facet's SVs first: an SV of its own whose extent lies within
`refdist` of this one's (widened as above), of any type or size, makes `other_alleles`. Otherwise the
facet's rule for the `sv` class (`absence_is_reference`): `reference` only for a whole-genome
callset without blocks that counted `sv` alleles and was not filtered to one class, `unknown`
otherwise. So an array, an exome, a small-variant WGS callset, a `.sv.` or `.cnv.` file (filtered to
SVs, so it says nothing even about its own class without `force`) leave an SV the other facet called
`unknown`, never `mismatch` and never `reference`. Blocks are not read: a block states the bases a
small-variant caller saw, and an SV's anchor base lies outside the event. The `site` reading of a
sequence allele is not done either: a row at the anchor says nothing about the event after it.

Measured 2026-10-26 on hu416394's Manta, Canvas and bcftools gVCF facets (one sequencing run;
aggregates only): Manta against Canvas meets no pair at the defaults, every Canvas allele being a
`<CNn>`; with `typeignore` 35 meet, 10 of them `<INS>` beside copy-number alleles, unsized. Against the
gVCF every SV is one-sided and `unknown`. Whether a copy-number allele should meet a `<DEL>` or `<DUP>`
as a loss or a gain is open (RM282).

**Breakends (RM282)** meet as truvari bench's `bnd_match` meets them (truvari 5.5, MIT;
`truvari/variant_record.py`, `VariantRecord.bnd_match`, `bnd_position` and `bnd_direction_strand`;
`truvari/comparisons.py`, `overlaps`; read from the 5.5.0 sdist 2026-10-26). An ALT is one of VCF
4.2-4.5's four forms, `t[p[`, `t]p]`, `]p]t` or `[p[t` (`t` the bases at POS, `p` the mate's
`contig:position`). Two breakend records meet when all of these hold:

1. **One contig**, and their positions' windows overlap strictly: each window is `pos - bnddist - c`
   to `pos + bnddist + c` (`SvMatching.bnddist`, truvari's default 100; `c` the record's
   `|CIPOS[1]|`), and truvari's `overlaps` is `max(starts) < min(ends)`. So two positions meet when
   they lie less than `2 * bnddist + ca + cb` apart: 199 bases apart at the default, not 200. As
   with `refdist` above, the code is followed, not the parameter's description ("Maximum allowed
   distance for breakend (BND) comparisons").
2. **The mates parsed from ALT lie on one contig**, their windows overlapping the same way, widened
   by `|CIEND[1]|`. With no CIEND the mate is widened by nothing: truvari does not fall back to CIPOS
   here, unlike the symbolic rule above.
3. **The same one of the four forms**: truvari's `(direction, strand)`, `[` or `]` and whether the
   bases come first, which is exactly the form.

Index 1 of CIPOS and CIEND, the upper offset, is truvari's as it runs: `bounds()` sets `key = 'CI' +
key` before `idx = 0 if key == 'POS' else 1`, so index 0 is never read. The code is followed, so a
lower offset alone (`CIPOS=-30,5`) widens by 5. `ci_needed` is True where a pair met only by the
widening. `bnddist` 0 meets only by an interval, as in truvari.

**The unit is the adjacency, not the record.** A junction is usually written twice, once from each
side (Manta's `MantaBND:…:0` and `:1`). Each record states two ends, a `(contig, position, side)`
each, the side being where the sequence that takes part in the join lies: left of POS for `t[p[`
and `t]p]`, right for `]p]t` and `[p[t`; right of `p` for `t[p[` and `[p[t`, left for the other two.
VCF's mate pairs (`t[p[` beside `]p]t`, `t]p]` beside `t]p]`, `[p[t` beside `[p[t`) state one pair
of ends. Records stating one pair of ends are one adjacency and one row. Records whose ends differ by
no more than their intervals are one adjacency too: each one's mate contig and side are the other's
own, and its mate position lies in the other's `[POS + CIPOS[0], POS + CIPOS[1]]`. Manta writes a
mate's POS at one end of the junction's microhomology, which CIPOS spans, and the ALT pointing at it
from the other record at the other end; on the PGP-HMS Manta files 448 of 3,350 and 342 of 1,228
records pair only that way, every one inside the interval and none ambiguous (F28). Only a
one-to-one pairing folds; a record that could be the mate of several stays alone, counted in a
warning. MATEID and EVENT are not read: the frame does not carry them (`dropped_fields`), and the
ends are what VCF requires mates to agree on. A record whose mate is missing (not written, cut by
`region`, on a contig the comparison leaves out) is its adjacency alone. An EVENT of several
adjacencies (an inversion's two) is several rows. An adjacency that joins an end to itself (a
foldback, `t]p]` pointing at its own POS) has both mates at one end; more records than its ends
hold are counted in a warning, as are mates whose genotypes differ.

Two adjacencies meet when any record of one meets any record of the other (truvari compares records,
and a pair of mates is two). The row shows the pair of records that met best, by truvari's score
`max(0, 1 - mean(|start distance|, |mate distance|) / (2 * bnddist))` times 100 (0 at `bnddist` 0),
then the nearer start and mate (its `MatchResult` order), then the record at the adjacency's first
end, so `a_pos` and `b_pos` lie on `chrom`. Each adjacency is taken once, best first (`--pick
single`), ties to each side's record order. Key `breakend`, `sv_kind` `BND`; `size_similarity` and
`reciprocal_overlap` are null, and `a_end_pos`/`b_end_pos` are the positions (a breakend's `end_pos`).
Genotypes compare as for symbolic alleles. A one-sided adjacency reads the silent facet's breakends
whose position window meets one of its records' (truvari's first test, any mate or form;
`other_alleles`), else the facet's absence rule for `sv` on a contig it declared, at both ends. A
symbolic SV of the other facet near a breakend does not make `other_alleles`.

What differs from truvari, each with its reason:

- **Mate contigs by canonical name.** truvari compares the ALT's contig text; here it is read with
  `canonical_chrom`, as `chrom` is, so Ensembl `2` and UCSC `chr2` are one contig, and a mate on a
  contig the two collections do not share is `contig_not_shared`. A mate in angle brackets
  (`C[<ctg1>:7[`, VCF 4.3 §5.4.2) is the exception: it names a contig of the `##assembly` file, not of
  the reference (§1.6.1), so it has no canonical name and is on neither collection. It is compared as
  truvari compares it, by its text as written (5.5's `bnd_position` reads `<ctg1>`, brackets and all,
  with `[^\[\]:]+`) and its position, on facets of one collection or of two, never `contig_not_shared`.
  Whether two facets' `<ctg1>` hold one sequence is not checked: their assembly files are not read.
  Left alone, such an adjacency reads `unknown`, never `reference`: a facet could state a join into
  another file's assembly contig only with that file, so its silence there says nothing (2026-10-26).
- **A contig name may hold `:`.** truvari's `[^\[\]:]+` refuses `HLA-A*01:01:01:01`; here the name
  runs to the last `:` before the digits.
- **Brackets must pair.** truvari reads `A[2:100]` as a right-direct join (it tests the first and last
  characters only); VCF has no such form, and it is `breakend_unparsed`.
- **Single breakends are counted, not dropped silently.** truvari's default `no_single_bnd` leaves
  `.A` and `A.` out of the comparison; here too, as `single_breakend` in `excluded`. A single
  breakend has no mate to compare, and truvari's `bnd_position` cannot read one either.
- **Adjacencies, not records.** truvari bench counts each record, so a mate pair met by a mate pair
  is two true positives; here it is one row.

Measured 2026-10-26 (aggregates only): hu416394's Manta file (GRCh38, 3,350 breakend records) holds
1,675 adjacencies and hu1C1368's (GRCh37, 1,228) 614; each against itself meets every adjacency
(`match`, none `ci_needed`). Against itself with the records at odd positions dropped, hu1C1368's
meets the 467 adjacencies that kept a record and leaves the 147 that kept none `only_a`, none
`only_b`. Manta against Canvas leaves all 1,675 `only_a`, `unknown`: Canvas writes no breakends. No
person here has two breakend facets, so no two callers' breakends were compared.

## A method-error atlas (RM279)

`method_atlas(pairs, regions)` (`just_vcf.atlas`) turns `concordance` over pairs of one person's
facets into rates per method and region (charter 8: a facet's reliability is stated by method and by
region). It returns one frame of counts and rates and never a per-site row: a consumer gets how often
DeepVariant and GATK disagree in homopolymers of 7 to 11 bases, not where. It is just-vcf's own shape
(charter 2; GA4GH has no model for it) and adds no column to `facets`, `calls` or `blocks`, so
`SCHEMA_VERSION` is unchanged.

**Regions** are a `RegionSet`: named classes of 1-based closed intervals on one assembly (`chrom`,
`pos`, `end_pos`, any contig spelling), each merged once when the set is built and reused for every
pair. `giab_stratifications(directory, assembly_label)` reads the GIAB v3.6 BEDs for GRCh37 or
GRCh38 found by file name under `directory`, and refuses a missing one by name rather than leaving
its class out:

| Class | GIAB v3.6 file (`GRCh38_…`) |
|---|---|
| `easy`, `difficult` | `notinalldifficultregions`, `alldifficultregions` |
| `low_mappability`, `segdups` | `lowmappabilityall`, `segdups` |
| `homopolymer_4to6`, `homopolymer_7to11`, `homopolymer_ge12` | `SimpleRepeat_homopolymer_<n>_slop5` |
| `tandem_repeat_le50bp`, `…_51to200bp`, `…_201to10000bp`, `…_ge10001bp` | `AllTandemRepeats_<n>_slop5` |
| `gc_lt25_or_gt65` | `gclt25orgt65_slop50` |
| `x_par`, `x_nonpar`, `y_nonpar` | `chrX_PAR`, `chrX_nonPAR`, `chrY_nonPAR` |

`all` is the class with no region: every allele the pair compared. Classes overlap (a homopolymer is
also `difficult`), so rows of one pair never sum across classes. A pair whose facets' `assembly_label`
is not the set's is refused: lift one first.

**One comparison per pair.** `concordance(a, b)` runs once and each allele is placed in every class
whose intervals hold its `pos` (the spelling of the facet that states it, `a` first). Cutting each
facet per class with `region=` would give the same counts except for an indel the two facets spelled
at different positions across an edge, which the cut splits into two one-sided alleles; placing keeps
it one paired allele. `by_contig=True` runs `concordance` per contig and sums the counts, the same
counts with one contig's sites in memory at a time. Rejected: `region=` per class, which repeats the
whole comparison for each of sixteen classes.

**The unit is a carried allele.** A row counts an allele when at least one side's genotype holds a
copy of it (its `gt_alleles` contain the ALT, an uncalled haplotype beside it included, or its
`copies` are above zero, where another ALT of a complex row holds the allele). An allele a
row states and no genotype carries (a DeepVariant `RefCall` `0/0`, a split row's `0` for its ALT,
`./.` on both sides) says nothing about either method's error and is left out of every count. Scope
`carried` counts all of them; `carried_pass` keeps those whose every stating row passed its filters
(FILTER `PASS`, or none written; a silent side has no row and passes).

**The two sides are put in method order** (caller, then caller version, unknown first) and `only_a`
and `only_b` follow them, so a pair counts the same whichever facet comes first and pairs of the same
two methods pool.

**`atlas`**, one row per (`region_set`, `assembly_label`, `a_caller`, `a_caller_version`,
`b_caller`, `b_caller_version`, `region_class`, `variant_class`, `scope`, `a_band`, `b_band`), sorted
by them (region classes in the table's order, bands in `CONFIDENCE_BANDS`' order after `all`). The
pooled rows have `a_band = b_band = all`, and every variant class but `sv` has one in every class, zero
counts included. Beside them, one row per cell of the two sides' **confidence bands** where alleles
fall: each side's band is that of `genotype_evidence`'s `confidence` of its row at the allele (the
lowest where it wrote several rows at the position), `0-10`, `10-20`, `20-30`, `30-40`, `40-60`,
`60-100` or `100+` phred, `none` where the row states no confidence or the side wrote no row. A
side's bands, summed over the other's, give the pooled row. Counts written before the bands (no band
columns) pool as pooled rows. The cells say how often two methods disagree where each was sure and
where it was not, which a method's average spreads over all its rows (docs/CONSENSUS.md: fitted on
one PGP pair, a CGI row at 100 phred or more errs on 0.05% of SNVs, against CGI's 0.51% average):

| Column | Dtype | Meaning |
|---|---|---|
| `pairs` | `u32` | pairs pooled into the row |
| `alleles` | `u64` | carried alleles compared |
| `match`, `mismatch`, `uncalled`, `unreadable`, `only_a`, `only_b` | `u64` | `concordance`'s `agreement` per allele |
| `ploidy_differs` | `u64` | of `called_both`, the alleles whose two genotypes were written with different ploidies (`concordance`'s flag): counted apart, inside `match` or `mismatch` as they agree |
| `called_both` | `u64` | `match + mismatch`: both sides called a genotype |
| `discordant` | `u64` | `mismatch` |
| `discordance_rate`, `discordance_ci_low`, `discordance_ci_high` | `f64` | `discordant / called_both` and its 95% Wilson score interval |
| `one_sided_rate`, `one_sided_ci_low`, `one_sided_ci_high` | `f64` | `(only_a + only_b) / alleles` and its interval |
| `one_sided_reference` | `u64` | one-sided alleles the silent side reads as reference (`absence` `reference` or `ref_block`): a stated disagreement, where the rest are silence |
| `graded_gl`, `p_same_gl_mean` | `u64`, `f64` | alleles with `basis = gl_flat_prior` and their mean `p_same` |
| `graded_gp`, `p_same_gp_mean` | `u64`, `f64` | the same for `basis = gp` |

A rate whose denominator is zero is null, never zero (charter 9). `p_same` is pooled per basis: a mean
over likelihoods and posteriors together would compare the two (Decision 3). The interval is the
Wilson score interval without continuity correction (Newcombe 1998, method 3), which stays inside
[0, 1] and is not degenerate at zero discordant alleles; it measures sampling error over the alleles
compared, not the error of treating a few pairs as a method's whole behaviour, so `pairs` says how
far a row generalizes.

`pair_counts(a, b, regions)` returns one pair's counts (the `COUNTS` columns: the keys, `pairs`, and
the summed counts with `p_same_gl_sum`, `p_same_gp_sum`), and `pool_counts(frames)` sums them per key
and computes the rates from the sums; `method_atlas` is the two together. Counts written per pair can
be pooled later and in any grouping. `scripts/method_atlas.py` builds an atlas from a plan of pairs,
skipping a pair `same_person` does not answer True for; docs/METHOD_ATLAS.md holds the first one's
headline rates.

## A consensus posterior (RM280)

`consensus(person, *, atlas, regions=None, prior="flat", …)` (`just_vcf.consensus`) reads every facet
of one person at the alleles they state and gives each allele a posterior over the copies of it the
person carries, weighing each facet by its method's error rate from a method-error atlas. The result is
one more facet (charter 5), derived from all of them: `relation = consensus` (Decision 12), one
`consensus` link from each parent, `derived_from` null. The parents' calls stay as they are, and each
row says how many facets supported its most probable genotype, contradicted it or were silent
(charter 8: the consensus resolves nothing away; it is a reading beside the others). The person chose
it after the atlas (Part 2 interview, E7). Its shape is just-vcf's own (charter 2): GA4GH has no
model for a cross-callset posterior. It adds no column to any frame, so `SCHEMA_VERSION` is unchanged.

**Parents.** A `Person` (refused when its check read two people, warned when not every pair read as
one) or `Loaded`s of one subject: two facets or more, none given twice, no consensus facet among
them, and none made from another given one (its lineage, through `links` and `derived_from`): a facet
and its `normalize`d child are the same reads twice. Kin facets (a PGP merge's freeBayes and VarScan2)
are allowed. All on one assembly: equal `ref_seqcol_digest`s, or, where any facet's REF was not read,
one equal `assembly_label` with a warning (charter 3). Anything else is `ReferenceMismatchError`:
lift first.

**The unit is one ALT allele**, as in concordance and the atlas. A row is `(chrom, pos, REF, ALT)`
with REF as written, upper-cased, and its vector is over the copies of that ALT: 0, 1, 2 diploid, 0, 1
haploid. The rows are every sequence ALT a parent states in a row the consensus reads (it passed its
filters under `pass_only`, and it is graded or carries a genotype), on every contig. A multi-allelic
site gives one row per ALT, each a marginal: copies of this ALT against anything else. Indels meet by
spelling only (VRS ids are not used): two parents that spell one indel differently give two rows,
each with the other parent silent or reading `other_span`. `normalize` both, or compare SNVs.

**What a facet says** about an allele, the first that applies:

| Reading | When | Its distribution over copies `q` |
|---|---|---|
| `stated` | a row states the allele | its vector (`genotype_evidence`'s `gl`, under a flat prior, or `gq`) summed by copies; else its genotype, a point mass. A `gp` row is read as its genotype unless `gp_as = likelihood`; an `hq` row (CGI's allele scores, a VCF's FORMAT HQ) has no vector and is read as its genotype, its confidence placing it in a band (`CONSENSUS_VERSION` 3; 4 for a VCF's HQ) |
| `site` | one row at the position over the same REF span, without the allele (a hom-ref record, another ALT) | zero copies, a point mass: its vector indexes other alleles |
| `ref_block` | a `ref` block covers the position | by `block_grade` (above, "The block rule"): under the default, the block's likelihoods over `[REF, <*>]`, `<*>` read as this ALT; graded by its GQ, zero copies at `1 − 10^(−GQ/10)` and the rest spread evenly over the other copies, at the ploidy its likelihoods give, else the site's; lacking what the rule reads, zero copies |
| `absence` | nothing covers it and `absence_is_reference(facet, class, force=force)` is True on a contig the facet declared or wrote | zero copies |

`gq_on_reference` reaches a row stating the allele, as `genotype_evidence` grades it; it does not
reach a block, which `block_grade` alone grades (so the default reads a block with a GQ and no
likelihoods as zero copies, unquantified, whatever `gq_on_reference` says).

Silent otherwise, each kind counted in the warnings per facet: `filtered` (`pass_only` and the row
failed a filter), `no_call` (no haplotype called and no vector), `covered` (inside another record's
span, a deletion), `several_rows` and `other_span` (concordance's rule: such a row holds no allele of
this span), `ref_unstated` (the one row there has a source REF that states no base and spells a
sequence ALT, Decision 2: its genotype was spelled against the `N`), `no_call_block`, `unknown`
(silence that says nothing: an array, a content filter, a class the facet did not call),
`ploidy_differs` (below).

A row whose source REF states no base (`ref_check = unstated`) and which spells a sequence ALT
states no allele of the consensus (`CONSENSUS_VERSION` 5): `N` → `G` read beside the reference's `T`
would be `T` → `G`, an allele the source never wrote. It still counts as a row at its position, so
the facet there is `ref_unstated`, never `ref_block` or `absence`: a whole-genome facet's rule would
read zero copies at a position where it wrote a call it could not spell.

**Error rates** (`method_error_rates(atlas)`). The atlas states how often two methods disagree, not
which erred. With independent rare errors a pair's genotype disagreement `d = mismatch / (match +
mismatch)` is about the sum of its two methods' rates. A haploid call against the same allele's
diploid homozygote (a male X written two ways) is a `match` in the atlas since F22, so two conventions
are no genotype error; a haploid call against a het stays a mismatch, a real disagreement. A method's
`genotype_error` is `(d_xy + d_xz − d_yz) / 2` where it forms a triangle with two
partners also compared with each other (the same-reads set), averaged over its triangles and kept at
or below its lowest pair discordance (`triangle`); otherwise half of each pair's discordance, averaged
over its partners (`half_pair`). Its `silence_error` is the share of the alleles a partner carries
that it left one-sided, the lowest over its partners: an upper bound, since a partner's false calls
count against it too. No rate goes below the 95% Wilson upper bound of zero events in the alleles it
rests on. Rows of `scope` (`carried_pass`, or `carried` when `pass_only=False`) only, assemblies pooled.
These are the `confidence_band = all` rows of `METHOD_RATES`, from the pooled atlas rows only.

**Rates by confidence band** (`confidence_band` a band, `genotype_attribution = band`), from the
atlas's cells: per pair and class, each cell's discordance `d_ij` is fitted, weighted by its alleles,
as `e_x(i) + e_y(j)`, one rate per band of each side (least squares, by alternating means). One pair
fixes these rates only up to a constant moved from one side to the other, as it fixes the averages
only up to their split; the split taken is even in the cell where both sides are at their most
confident band holding `MIN_BAND_ALLELES` (1,000) alleles both called, the `half_pair` rule applied
where neither row's own confidence explains a disagreement. A band is given where it holds
`MIN_BAND_ALLELES` and both sides have two such bands (a side that states no confidence cannot be told
apart from its partner band by band); floored at the Wilson bound of its alleles; averaged over
partners. A band row carries the method's average silence rate. The `all` rate stays the pooled rule,
so a method's bands weighted by their alleles need not average to it. Rejected: one side's band less
its partner's average, which credits the partner with its unsure rows' errors at alleles where it
was sure (it put freeBayes's confident bands at the Wilson floor, about 4e-6); and the row's own
`1 - p` plus a method residual, which the PGP cells refute: CGI's rows at 20 to 30 phred state at
most 1% error and disagree with freeBayes on 5.4 and 6.8% of SNVs in the two pairs.

A facet takes, per region class and variant class, its caller's rates at its version (`method`),
else its caller's pooled over versions without self pairs (`caller`), else the stated default
(`default`): `default_rate=` when given, otherwise the worst rate any measured method has in that
class. An unmeasured method is taken as no better than the worst one measured, not as average. Every
row counts the facets that spoke with a default (`defaulted`), and `Consensus.rates` lists each
parent's rates and their basis. With `regions` (a `RegionSet`) an allele is weighted by the highest
rate among the classes holding its position, `all` where none does; without it, by `all`, with a
warning where the atlas carries rates in other region classes (the atlas holds class names, not their
intervals). A facet's band rates are its method's, else its caller's (`rate_basis`), never a
default: a band no atlas row measured is absent from `Consensus.rates`, and a row in it takes the
average. Among the classes holding an allele, a class without the row's band gives its average. A
facet weighted by band rates its caller was measured at in other versions is named in a warning: a
band's measured error belongs to the confidence scale of the facets it was measured on, and two
facets of one caller can write that scale differently (docs/CONSENSUS.md, freeBayes 0.9.15 against
0.9.20).

**The posterior.** A stated reading whose row's own confidence falls in a band its method has a rate
for gets that rate `e_b`, the measured share of such rows whose genotype was wrong: `L(top) = 1 − e_b`
for the copies its distribution `q` puts first (the fewer on a tie), and `e_b` shared among the other
copies in proportion to `(1 − e_b)·q(k) + e_b·(1 − q(k))/P`, which ranks them as `q` does and is never
zero. The band rate is the whole error of such a row, so the row's own `1 − p` is not added to it.
Any other reading (a site reading, a stated one without a measured band) gets the genotype error, a
block or an absence the silence error, and becomes the likelihood `L(k) = (1 − e)·q(k) + e·(1 −
q(k))/P` over `k = 0..P` copies: a genotype gets `1 − e` and spreads `e` evenly over the rest, and a
graded vector is tempered so that no facet says more than its method's measured error allows (a PL
of 3,000 against two facets that disagree is still one facet). The posterior is `prior(k) · ∏ L_f(k)`, normalized, in log
space in Float64; the prior is flat over the copies (`prior = "flat"`, the default) unless the caller
supplies a population prior (below). `P` is the ploidy most speaking readings state (diploid on a tie or where none states one); a reading of
another ploidy is silent (`ploidy_differs`). The facets are taken as independent: two callers on
the same reads are not, so where they agree the posterior is overconfident, and the atlas cannot
correct that.

**A population prior** (`prior=AlleleFrequencyPrior(name=, assembly_label=, frequencies=)`, desk
#96). Which catalogue to use, and whether to use one, is the person's decision: nothing is downloaded
or packaged, and the default stays flat. The table (`ALLELE_FREQUENCIES`) is keyed as consensus rows
are, `chrom`, `pos`, `ref`, `alt` on the facets' assembly, with `alt_frequency` the population
frequency of that ALT; VRS ids are not used, since consensus rows carry none. On the way in contigs
are made canonical as the frames' are (`chr1` is `1`, `chrM` is `MT`) and REF and ALT upper-cased; a
missing column, a null, a `pos` below 1, a frequency outside [0, 1] (NaN included) or an allele named
twice is refused, never patched, and so is the name `flat` (the default's). A table on another
`assembly_label` than the facets' is refused (lift one first), as `regions` is. An allele the table
holds gets Hardy-Weinberg genotype priors over its row's ploidy from its frequency `f`: diploid
`(1 − f)², 2f(1 − f), f²`, haploid `1 − f, f`. That is the marginal of copies of this ALT against
anything else, so a multi-allelic site's rows each take their own ALT's frequency. Every other allele
keeps the flat prior, and `Consensus.support.prior` names which each row took. A frequency of 0 or 1
makes the genotypes it rules out impossible whatever the facets say (the confidence is then infinite):
a catalogue that lists an allele at 0 is a strong statement to multiply in. The warnings count the
alleles the table held, and say when it held none (most often a table on other contig names or
another assembly). With `gp_as = likelihood` a warning says the population prior then counts more than
once. The table's name and the `sha256` of its normalized CSV enter the `facet_id` and the facet's
`derivation_inputs` (`role = prior`); the flat prior adds nothing, so a default consensus keeps the id
it had.

**The facet.** `facets`: `caller = "just-vcf consensus"`, `caller_version = CONSENSUS_VERSION`, one
`merge` step in `method_steps`, `derivation_inputs` one `error_rates` entry whose `sha256` digests
`Consensus.rates` (and a `prior` entry under a population prior), `absence_means = unknown`, `variant_classes` counted from its rows, the parents'
`assembly_label`, `ref_seqcol_digest`, `resolve_version`, `callset_type` and `subject` where they share one,
and `resolve_inputs` a digest of every parent's REF state (its `facet_id`, `ref_seqcol_digest`,
`resolve_version` and `resolve_inputs`, in `facet_id` order). `facet_id`
digests the sorted parent ids, `CONSENSUS_VERSION`, the options and the rates. A parent keeps its id
when `resolve_ref` reads it, so the id does not change with the parents' REF states, and
`resolve_version` is null where they differ: `resolve_inputs` is what lets the durable store tell two
builds over parents in other states apart. Rejected: digesting the parents' states into the
`facet_id`, which would change the id of every consensus already made. `calls`, one row per
allele in canonical order: `alts = [alt]`, `source_alleles = [REF, ALT]`, `source_ref` the REF,
`ref` the same where a parent read it from the reference (charter 10), `ref_check = not_run`, `gp` the
posterior (VCF genotype order over `[REF, ALT]`), `gt` and `gt_alleles` its most probable genotype
(the fewer copies on a tie), `evidence_basis = assay`, the rest null. `genotype_evidence` grades
every row `gp`. `blocks` is empty.

**`Consensus.support`**, one row per `calls` row in its order: `chrom`, `pos`, `ref`, `alt`,
`variant_class`, `ploidy`, `copies` (the most probable), `confidence` (`-10·log10(1 − p)` of it from
the log-space posterior: `gp` is Float32 and rounds a remainder below about 1e-45 to zero, which
`genotype_evidence` then reads as infinite), `facets`, `supported`, `contradicted`, `silent` (they
sum to `facets`), `defaulted`, `gp_as_likelihood` and `banded` (the speaking facets read with a
default rate, with a posterior multiplied in, and with their row's band rate), and `prior` (`flat`, or
the population prior's name where its table held the allele). A sidecar, not `calls` columns: Decision 6's columns hold what a
source states about a site, and these say how the consensus came about. Rejected: columns in `calls`,
which every other facet would carry null and which would bump `SCHEMA_VERSION`; and the
per-facet readings as rows, which copy the parents' calls into a second place.

The first check on a public genome is docs/CONSENSUS.md.

Rejected: a winner per site (charter 8); a joint posterior over a site's union of alleles, which the
parents' allele lists do not share and which no other just-vcf measurement uses; multiplying GP as a
likelihood by default (Decision 3: a posterior already holds its source's prior, which would count
once per facet); one rate for every reading, which would let a whole genome's silence (one-sided
rates of 1 to 30%, docs/METHOD_ATLAS.md) weigh like a stated genotype (0.05 to 3%).

## Open questions

Q1 (CGI `GC?`), Q2 (array code 25), Q3 (withheld VRS ids), Q4 (sample names), Q5 (`end_pos`) and Q6
(PSL keys) were settled on 2026-09-30 and are folded in above.


## Changes since sign-off

- **An enum member, landed in 0.0.22** (2026-10-27, F23: a source REF that states no base):
  `RefCheck.UNSTATED` (`unstated`), a source REF made only of `N` where the reference reads a base,
  the REF DRAGEN's and Canvas's copy-number VCFs write on every record. Adding a member changes the
  `pl.Enum` dtype of `calls.ref_check`, a schema change; pre-1.0, it breaks a consumer that matches on
  the three values. Rule changes (Decision 2): `resolve_ref` writes `unstated` where it wrote
  `mismatch` (`RESOLVE_VERSION` 13); `liftover`, `write_vcf` and `write_vcf_zarr` take such a row as
  they take a `match` when no ALT is a sequence allele, and leave it behind (`source_ref_unstated`,
  `LIFT_VERSION` 13) or refuse it (`ref_unstated`) when one is (merge audit 13, 2026-10-28);
  `mint_vrs_ids` mints its REF slot and withholds a sequence ALT as `ref_unstated` (`MINT_VERSION` 4);
  `normalize` leaves it as written, counted; `query` answers from it neither `stated` for a sequence
  ALT nor `site_genotyped`; `concordance` leaves a sequence ALT out (`excluded`, `ref_unstated`) and
  reads no such row for another facet's allele at its position, and `consensus` states no allele from
  it and is silent there as `ref_unstated` (`CONSENSUS_VERSION` 5; 2026-10-31); `genotype_evidence`,
  `genotype_evidence_at` and `reference_support` grade no genotype on it, `grade_basis = ref_unstated`
  (a new `GradeBasis` member in those computed frames, none stored; 2026-10-31), and
  `sex_chromosome_evidence` counts no ALT on it (2026-10-31), and `read_site_genotypes`, behind
  relatedness and ancestry, reads no genotype from it at a panel site, `basis = ref_unstated` (a new
  `SiteBasis` member in that computed frame, counted by ancestry as `n_sites_ref_unstated`, none
  stored; 2026-10-31); `cis_phased_blocks` takes no slot of it into a block, REF included, and
  withholds a block it would join as `ref_unstated`, which moves only `reference_alleles=True`
  (2026-10-31); `facets.variant_classes`, counted at
  load before any REF is read, cannot apply the rule and counts such an ALT by its lengths. **Values change**: the 1,324 rows of the three PGP-HMS copy-number files
  F23 measured were all `mismatch` and none lifted; all are `unstated`, and 1,224 lift. A store of
  0.0.21 migrates with the enum widened and its values kept, its facets holding such rows named.

- **No version change** (2026-10-26, F32: a breakend's mate in the Zarr export): no column changes.
  `write_vcf_zarr` names a breakend mate inside `variant_allele` as `write_vcf` names it and lists
  every mate contig in `contig_id` ("VCF Zarr export" above). **Values change** in the store only:
  under `stored` a mate the source spelled `chr2` is stored `2`, under `chr` a `2` is stored `chr2`,
  and `contig_id` gains an entry per mate contig no row stands on.
- **No version change** (2026-10-26, assembly mates lifted and compared): no column changes.
  `liftover` lifts a breakend whose mate is in angle brackets (`C[<ctg1>:7[`, a contig of the
  `##assembly` file) by its own position, keeping the mate as written and turning the form by POS's
  strand alone ("Lifting breakends", step 7), where it left the row behind as
  `breakend_mate_not_lifted`. Lifted facets get a new `facet_id` (`LIFT_VERSION` 12). `concordance`
  compares such a mate by its text as written and no longer leaves it out as `contig_not_shared` when
  the two collections differ ("Breakends"). **Values change**: such rows are lifted, and compared.
- **No version change** (2026-10-26, F30's `source` mates): no column changes. Under
  `contig_names="source"`, `write_vcf` writes a breakend mate on a contig the facet does not know as
  the ALT spelled it, where it wrote the canonical name, and refuses one ALTs spelled two ways; a
  contig the facet knows keeps the facet's spelling ("Lifting breakends" above). **Values change** in
  the written file only, under `source`: a chr-named source's mate on a contig with no row and no
  header line comes back `chr2`, where it was written `2`.
- **No version change** (2026-10-26, F30's assembly mates): no column changes. `write_vcf` declares
  no `##contig` line for a breakend mate in angle brackets (`]<ctg1>:7]T`, a contig of the
  `##assembly` file) and writes that ALT as held under every style, and declares only the mates of
  the rows it writes ("Lifting breakends" above). **Values change** in the written header only: the
  `##contig=<ID=<ctg1>>` line, and a line for a mate only a left-out row named, are gone.
- **No version change** (2026-10-26, F30: a breakend's mate written as CHROM is): no column
  changes. `write_vcf` writes the mate contig inside a breakend ALT under the `contig_names` style
  CHROM is written under, and declares every contig a written ALT names ("Lifting breakends" above).
  **Values change** in the written file only: under `stored` a mate the source spelled `chr2` is
  written `2`, and the header gains a `##contig` line per mate contig no row stands on.
- **No version change** (2026-10-26, RM6: breakends lifted): no column changes. `liftover` lifts
  breakends by their own position and each mate's ("Lifting breakends" above), where it left every
  one behind as `breakend_not_lifted`. `UNLIFTED_REASONS` gains `breakend_mate_not_lifted`,
  `breakend_ends_apart`, `breakend_ends_reversed` and `breakend_distance_changed` (an open
  vocabulary); `breakend_not_lifted` stays, narrowed to the breakends the lift cannot read (an ALT
  with brackets that is none of the four forms, a breakend beside another kind of allele, a REF of
  more than one base). Lifted facets get a new `facet_id` (`LIFT_VERSION` 11), and a lifted facet's
  `contig_map` also lists the contigs only a breakend ALT names. **Values change**: lifted breakend
  rows appear, and the reason of those left behind. A narrowed meaning is a break of an open
  vocabulary, pre-1.0.

- **No version change** (2026-10-26, RM282 breakends in concordance): no column of `facets`, `calls`
  or `blocks` changes, and `sites` keeps its columns. `concordance` compares breakends ("Breakends",
  above): `key` adds `breakend`, `sv_kind` takes `BND`, `SvMatching` adds `bnddist` (100).
  **Values change**: a breakend, counted in `excluded` as `breakend` before, is now a row of
  `sites`, one per adjacency, so `sites`, `summary` and RM279's atlas (class `sv`) count it.
  `EXCLUDED_REASONS` drops `breakend`, which nothing emits any more, and adds `single_breakend` and
  `breakend_unparsed`; `contig_not_shared` also covers a breakend whose mate lies on such a contig.
  The removal is a break of an open vocabulary, pre-1.0.

- **No version change** (2026-10-26, RM3's reference-supported positions): no column of `facets`,
  `calls` or `blocks` changes. `reference_support` and `absence_support` are computed on demand, as
  `genotype_evidence_at` is ("Graded genotype evidence"; Decision 4), and `REFERENCE_SUPPORT` is
  their answer's shape, outside the three frames. `absence_is_reference` answers as before. Additive:
  minor.

- **No version change** (2026-10-26, F9's row sample): no column changes. The rows `resolve_ref`
  scores to infer an array's build and check its strand leave out rows not over one base, so a Genos
  insertion site (`G/G` inserted after a `G`) is no longer scored as an SNV at its anchor
  (`RESOLVE_VERSION` 12). A facet's inferred `assembly_label`, its `strand_check`, and the rows a
  mixed strand withholds can change where such rows tipped the scores. Values only.

- **No version change** (2026-10-26, RM4's GRCh36 lengths): no column changes. `infer_assembly`
  knows NCBI36's primary-chromosome lengths, so a VCF, or a coverage facet, whose `##contig` or BAM
  header lengths are build 36's gets `assembly_label = GRCh36`, basis `inferred`, where it was null.
  A value that was unknown becomes known; no known value changes. Additive: minor.

- **No version change** (2026-10-26, RM4's NCBI36 accessions): no column changes. NCBI36's
  assembled molecules' RefSeq accessions (`NC_000001.9` … `NC_000024.8`) map to canonical `chrom`
  (Decision 7), from NCBI's NCBI36 assembly report, as GRCh37's and GRCh38's do. **Values change**
  on files that already load: a contig named by one reads as `1`…`Y` in `calls.chrom` and
  `contig_map`, where it kept its accession, and a header naming its chromosomes so gets
  `assembly_label = GRCh36` from their lengths, where it was null. The `chrom` change is a known
  value replaced, not an unknown filled, so it is a rule change on existing files, not an additive
  one; the load spill is keyed by its rows' bytes, so a cached load of such a file is written again,
  never read back with the accession.

- **No version change** (2026-10-26, RM6: symbolic SVs lifted): no column changes. `liftover` lifts
  symbolic SVs and copy-number segments by both ends ("Lifting symbolic SVs" above), where it left
  the SVs behind as `allele_not_sequence` and lifted a segment with its `end_pos` cut to its REF span.
  `UNLIFTED_REASONS` gains `breakend_not_lifted`, `sv_end_unlifted`, `sv_ends_apart`,
  `sv_ends_reversed` and `sv_length_changed` (an open vocabulary). Lifted facets get a new `facet_id`
  (`LIFT_VERSION` 9). **Values change**: a lifted segment's `end_pos`, and breakends' reason.

- **No version change** (2026-10-26, RM6: blocks lifted): no column changes. `liftover` lifts the
  `blocks` frame ("Lifting blocks" above), where the lifted facet had none; `Lifted` gains
  `block_reasons`. Lifted facets get a new `facet_id` (`LIFT_VERSION` 8). Additive: minor.

- **No version change** (2026-10-26, RM276's store migration): no column of `facets`, `calls`,
  `blocks` or `links` changes. `migrate_store` rewrites a store of an older version to this one by a
  step per bump since 0.0.14 ("Store migration" above); the manifest gains `migrations` (empty by
  default, so `store_layout` stays 1 and every existing manifest reads as before), and `loaded()`'s
  warnings gain a line per value a migration filled unknown. A store of another version is still
  refused on read, now with a message naming the call that migrates it, or saying it is from a later
  just-vcf. Additive: minor. From now on a version bump that changes a stored frame comes with its
  step in `just_vcf.store_migration`, or the suite fails.

- **No version change** (2026-10-26, RM3's insertion points and Genos spans): no column's shape
  changes. `resolve_ref` reads an insertion point's anchor base (CGI's native files, a CASAVA indel
  table, Genos) and REF over a Genos row's stated span ("How RM3 fills these slots" above;
  `RESOLVE_VERSION` 10); `normalize` aligns the anchored rows and anchors a deletion at a contig's
  first base on the base after (`NORMALIZE_VERSION` 9). Rows that stayed unread now get `ref`,
  `alts`, `gt` and `end_pos`, except a point whose every called haplotype is the empty allele, which
  stays unread (`RESOLVE_VERSION` 11, "How RM3 fills these slots"). **Not additive for a reader of a resolved CGI or CASAVA facet's alleles**:
  an insertion point's `alts`, `gt_alleles` and `source_alleles` are now spelled with the anchor base
  (`["", "GG"]` read `["G", "GGG"]`), its `ref_check` goes from `not_run` to null, and one before a
  contig's first base moves from `pos` 0 to 1. The loaded frame is unchanged. Allowed before 1.0 and
  said so in the changelog; rejected: leaving the alleles unanchored beside `ref` = the anchor, which
  every reader of `[ref] + alts` (VRS, `write_vcf`, concordance) would read as an SNV.

- **Two enum members, outside the three frames** (2026-10-26, RM274's first-degree split):
  `RelationBand.PARENT_CHILD` (`parent_child`) and `RelationBand.FULL_SIBLINGS` (`full_siblings`),
  just-vcf's own terms, in the `band` `pl.Enum` of `RELATEDNESS_FRAME` and `PEDIGREE_CHECKS` and in
  `EXPECTED_KINSHIP`'s `expected_bands` ("Relatedness and pedigrees" above). No column of `facets`,
  `calls` or `blocks` changes, so by RM274's precedent no `SCHEMA_VERSION` change; the lead's call at
  merge. Legality: the vocabulary grows without a removal or a rename, but the change is **not
  additive for a reader of `band`**. A pair that read `first_degree` now reads `parent_child` or
  `full_siblings` wherever its IBS0 share is told, so `first_degree` narrows from "first degree" to
  "first degree, the kind not told", and a consumer selecting `band == "first_degree"` gets fewer rows.
  That is a redefinition of a value a consumer reads (§8 says add, not redefine), allowed before 1.0
  and said so in the changelog; no sibling repo reads `band` (just-prs and just-dna-lite checked
  2026-10-26). The pedigree check's rule changes with it: a parent link measured `full_siblings`
  (or the reverse) is `disagrees`, where it was `agrees`. Rejected: a separate nullable column (say
  `first_degree_kind`) beside an unchanged `band`, strictly additive but a second column every reader
  must join to the first to know what a pair is.

- **No version change** (2026-10-26, RM278's person-level coverage): no column of `facets`, `calls` or
  `blocks` changes. Three nullable or defaulted fields outside the frames: `CoverageFacet.subject`,
  `Person.coverage` (empty by default) and `SexChromosomeEvidence.coverage_pairing` (null from
  `sex_chromosome_evidence`, which is unchanged). Additive: minor.

- **No version change** (2026-10-26, an allele index `gt` cannot hold): the VCF, key=value and Nirvana
  loaders count a genotype whose index names an allele of the row past 255 (`0/300` over 300 ALTs)
  under its own warning; the frame is as before (`gt [0, null]`). Widening `gt` to `u16` is left to
  the person (Decision 2, "Genotype text that is not allele indices").

- **No version change** (2026-10-25, genotype text that is not allele indices): the key=value table
  and Nirvana JSON withhold such a genotype whole (`gt`, `gt_alleles`, `phased` null, `source_gt`
  kept), where they read its readable slots (`0/x` was `[0, null]`). No column's shape changes; `gt`
  gains a stated reason to be null (Decision 2, "Genotype text that is not allele indices").

- **No version change** (2026-10-24, layer 7, `hq` through a split): `normalize` writes a split row's
  `hq` null on each GT slot it writes `.` (the slot held another ALT), and the row's `hq` null when no
  slot keeps a score, so `hq` holds one score per called GT slot, as every loader writes it. Before,
  the parent's whole list stayed on every split row (`NORMALIZE_VERSION` 8). `rejoin_split` joins
  `hq` per slot, as it joins GT, so `normalize`'s rows and other tools' (bcftools copies HQ whole, vt
  keeps it on a row's own slot) give the parent's scores back (`REJOIN_VERSION` 5). Grades do not
  move: `genotype_evidence` reads `hq` only on a fully called genotype, and a split row with a `.` slot
  is not one.

- **No version change** (2026-10-24, layer 7, a VCF's haplotype qualities): the VCF loaders read
  FORMAT `HQ` into `calls.hq` ("A VCF's haplotype qualities"), a column 0.0.19 defined as VCF's HQ
  and filled from CGI's native files only. No column's shape or meaning changes, so `SCHEMA_VERSION`
  stays; a VCF facet with HQ now has `hq` where it was null, and `FORMAT/HQ` leaves its
  `dropped_fields` (an HQ declared with a Number other than `2`, `P` or `.` stays there). Rules
  changed: `genotype_evidence` grades such a row by its lowest HQ where GP, GL and GQ do not
  (basis `hq`, unchanged), so a VCF reference genotype with HQ is `quantified` where it was
  `gt_only`; `consensus` weighs it by its band (`CONSENSUS_VERSION` 4). `REWRITE_VERSION` 10: the
  rewritten copy declares and keeps HQ.
  `write_vcf` writes `hq` as FORMAT `HQ`, declared `Number=2`, so a CGI `var` facet's scores survive
  its VCF export and a VCF with HQ round-trips; a haploid row is written `51,.`, which the loader
  reads back as `[51]`, and a row of more than two GT alleles has its `hq` left out, counted.
  `write_vcf_zarr` writes `hq` as `call_HQ` on the ploidy dimension, a row of any ploidy whole.

- **A calls column, landed in 0.0.19** (2026-10-24, layer 7, CGI `var` confidence): `calls.hq`
  (`list[u16]`, parallel to `gt`), each haplotype's score, filled from CGI's native files'
  `varScoreEAF` ("CGI's allele scores", Decision 6); null from every other loader. Additive. Rules
  changed: `genotype_evidence` grades a fully called genotype whose every slot is scored by the lowest
  score (basis `hq`), so a CGI `var` or masterVarBeta row is `quantified` where it was `gt_only`, and
  the method-error atlas bands it; `consensus` weighs such a row by its band's rate where the atlas
  measured one (`CONSENSUS_VERSION` 3). `liftover` carries `hq` (`LIFT_VERSION` 7) and `rejoin_split`
  keeps it where the split rows agree (`REJOIN_VERSION` 4). `varScoreEAF` and
  `allele<n>VarScoreEAF` leave CGI's `dropped_fields`.

- **No version change** (2026-10-24, layer 7, the options for desk #94–#96): no column of `facets`,
  `calls` or `blocks` changes, and no default output does. `BlockGrade` (`block_grade=` on
  `genotype_evidence`, `genotype_evidence_at` and `consensus`) names the rule a `ref` block is graded
  by: `likelihoods` (the default, as before), `gq`, or `lower_of_gq_and_likelihoods` ("Graded
  genotype evidence", the block rule). A non-default rule enters the consensus `facet_id`; the default
  leaves it as it was. `consensus(prior=)` also takes an `AlleleFrequencyPrior`, a caller's
  allele-frequency table (`ALLELE_FREQUENCIES`, "A consensus posterior"), and `CONSENSUS_SUPPORT`
  gains `prior` (`flat`, or the table's name), its last column: a sidecar outside the three frames,
  as `banded` was, so by that precedent no version change; whether it takes one is the lead's call at
  merge. Every other support column, and every default row, is as before.

- **No version change** (2026-10-23, calibrating the consensus, layers 5 and 7): no column of
  `facets`, `calls` or `blocks` changes. The atlas gains the key columns `a_band` and `b_band` and its
  rows per cell of the two sides' confidence bands; `METHOD_RATES` and `FACET_RATES` gain
  `confidence_band`; `CONSENSUS_SUPPORT` gains `banded`. Rules changed: a consensus row stated in a
  measured band of its method's confidence is weighted by that band's rate (`CONSENSUS_VERSION` 2),
  and a consensus given an atlas with region rows and no regions warns.

- **No version change** (2026-10-23, the lead's three follow-ups to that audit): no column changes.
  Rules changed: `reference`'s two conflict checks compare base names as written, before either is
  `<file>`, and their warnings name no file; a `caller_version` written as a path is `<path>`, as a
  step's version is, compared as written first (Decision 5); `ftdna_bigy`'s `family` is
  `variant_calls`, so `chip_absence` refuses a Big Y facet as not an array. Values change on facets
  these reach.

- **A facet column, landed in 0.0.18** (2026-10-23, the lead's decision on GL and PL): `facets.gl_source` (`schema.GlSource`: `gl`, `pl`, `gl_over_pl`) names the FORMAT key `gl` is read from (Decision 5). Adding a column changes the `facets` schema; `SCHEMA_VERSION` went to 0.0.18 at merge. No rule changes: `gl` loads as before.

- **No version change** (2026-10-23, the audit of 5b6ec3b..01fca0b): no column changes. Rules
  changed: a lifted facet's `reference` is null ("Lifting", step 5; `LIFT_VERSION` 6); `reference`
  and a method step's reference and chain options keep a base name only when it is a known
  reference's, else `<file>`, and a step's `tool` or `version` written as a path is `<path>`
  (Decision 5); `chip_absence` reads only a `genotype_table` facet as an array; `write_vcf_zarr`
  writes a PSL-phased row unphased and names the `fileformat` `write_vcf` would. Values change on
  facets these reach.

- **No version change** (2026-10-23, RM274's pedigree check, layer 4): two computed frames outside the
  three, `EXPECTED_KINSHIP` and `PEDIGREE_CHECKS` ("Relatedness and pedigrees" above). No column of
  `facets`, `calls` or `blocks` changes; a pedigree individual is now matched to a `subject` by id.

- **An enum member, landed in 0.0.17** (2026-10-23, RM2's CGI hapLink): `PhaseSetKind.CGI_HAPLINK`
  (`cgi_haplink`), CGI's `hapLink` as a phase set, just-vcf's own term. Adding a member changes the
  `pl.Enum` dtype of `calls.phase_set_kinds` and of `cis_phased_blocks`' `phase_set_kind`, which is a
  schema change. Rule changes (Decision 2): CGI's native rows fill `phased`, `phase_sets` and
  `phase_set_kinds` from `hapLink`, where they were null, and `hapLink` leaves `dropped_fields`;
  `cis_phased_blocks` gives one block per hapLink and chromosome; `write_vcf` writes a haplotype
  phased by a hapLink `/`, with no PS (VCF has no field for a hapLink, and `|` without its set would
  read back as the chromosome's implicit set), and `write_vcf_zarr` its `call_genotype_phased` false.
  Values change on CGI files that already load.

- **No version change** (2026-10-23, RM2's CGI read counts): no column changes. Rule changed
  (Decision 3): a `masterVarBeta` row fills `ad` and `source_alleles` from its per-allele read counts,
  where they were null. Values change on files that already load; `var` rows are unchanged.

- **No version change** (2026-10-23, the audit of RM282's derived steps): no column changes. Rule
  changed in `liftover`: a swapped row's `alt_copy_numbers` are null, as its `dosage` is
  (`LIFT_VERSION` 3). `normalize` gets a new `facet_id` (`NORMALIZE_VERSION` 6) for RM282's rule,
  which changed its output without a bump: per-ALT SV columns split with their ALT, a moved row's
  intervals cleared. Rule changed in `rejoin_split`: a joined site keeps `sv_type` and `copy_number`
  where its rows agree and concatenates `cipos`, `ciend` and `alt_copy_numbers` in its `alts` order,
  where it wrote them null (`REJOIN_VERSION` 2).

- **Additive, landed in 0.0.15** (2026-10-23, PSL/PSO/PSQ phase sets, layer 8): two
  nullable `calls` columns, `phase_set_ordinals` and `phase_set_qualities` (Decision 2, "How the VCF
  loader reads PSL"), and three columns of the computed `cis_phased_blocks` frame, `withheld_chrom`,
  `withheld_pos` and `withheld_allele`. A new nullable column is a minor bump by the rule above,
  0.0.15; `liftover` and `rejoin_split` carry the two columns (`LIFT_VERSION` 5, `REJOIN_VERSION` 3). Values change
  on VCF 4.4+ files that declare PSL: they loaded with PSL dropped (or not at all, Number=P, P10), and
  now carry its sets, and `phased[0]` where PSL names the first slot's set.

- **No version change** (2026-10-23, RM274): no column of `facets`, `calls` or `blocks` changes.
  Added outside the three frames: `RELATEDNESS_FRAME` and the two pedigree frames ("Relatedness and
  pedigrees" above). `Person.of`'s default check is now `same_person`.

- **Four facet columns, landed in 0.0.16** (2026-10-23, layer 2, the method chain):
  `reference`, `aligner`, `aligner_version` and `method_steps` (Decision 5), all nullable and additive,
  and `phasing_method` now filled from the chain's one `phase` step. A new nullable column is a minor
  bump by the rule above, 0.0.16. `facet_id` digests none of them.

- **No version change** (2026-10-23, RM281's written no-calls): no column changes, no enum member
  added. Rule changed (Decision 4): `chip_absence` gives a site the export writes as a no-call
  (`--`, `0 0`) `no_call`, where it gave `in_export`; `in_export` now means a genotype. Computed on
  request, so nothing stored changes.

- **No version change** (2026-10-23, RM281's I/D calibration): no column changes. Rule changed
  (Decision 2): `I`/`D` codes are resolved only on a chip `probes.ID_CALIBRATION` passes
  (`23andme_v5`, `ancestrydna_v2`); every other chip's code rows stay as loaded, counted
  `id_orientation_unverified`. Values change on `23andme_v3`, `23andme_v4`, `23andme_v1_v2`, Final
  Report and unknown-chip facets when the manifests are configured. `RESOLVE_VERSION` 8.

- **One facet column, landed in 0.0.12** (2026-10-23, RM273): `facets.subject`, the caller's
  label for the person a facet belongs to, supplied at `load(..., subject=...)` and inherited by every
  derived facet. Nullable and additive. `facet_id` does not digest the subject (Decision 5).

- **No version change** (2026-10-23, RM276, the durable store): `Store` holds the four frames per
  facet as parquet with a manifest, and `write_vcf_zarr` exports one facet's calls ("The durable
  store"). No column of any frame changes; `SCHEMA_VERSION` is not bumped.

- **A fourth frame, landed in 0.0.13** (2026-10-23, RM275, the `links` frame): `Loaded.links`
  (Decision 12), and no column of `facets`, `calls` or `blocks` changes. Rules added: every derived step adds a lineage link and
  carries its parent's; a load of several facets links them (`caller_split`, `appended_export`);
  `link_disagreements` states how `links` agrees with `derived_from` and `relation`.

- **Additive, landed in 0.0.14** (2026-10-23, RM282 structural-variant facets): five
  nullable `calls` columns, `sv_type`, `cipos`, `ciend`, `copy_number` and `alt_copy_numbers`
  (Decision 6, "Structural variants"). Rules changed: `<CN0>`, `<CN2>`… end at INFO `END` (the `end_pos`
  rule); an `ALT .` record with a copy number is a call, not a block (Decision 4); `vrs_withheld` gives
  `sv_no_vrs_model` for an SV allele where it gave `symbolic` (`symbolic` now names only symbolic
  alleles that are no variant). Values change on SV and CNV files that already load.

- **No version change** (2026-10-23, RM284): no column changes. A computed frame,
  `cis_phased_blocks`, exports the phase sets as VRS CisPhasedBlocks ("CisPhasedBlock export").

- **No version change** (2026-10-23, RM281's allele pairs): no column changes. Rule changed (Decision
  2, "The strand of an array's genotypes"): in a `strand_check = mixed` facet, a withheld row whose
  rsid names a configured manifest's SNV probe, with a pair that holds the reference base and is not
  `A/T` or `C/G`, is read by that pair: forward as written, or minus with `gt_alleles` complemented.
  Values change on mixed facets when the manifests are configured. `RESOLVE_VERSION` 6.

- **No version change** (2026-10-23, RM281's chip manifests): no column changes. Rule added (Decision
  2): `resolve_ref` resolves an array's `I`/`D` rows against the indel probes of the Illumina manifests
  under `JUST_VCF_MANIFESTS_DIR`, filling `ref`, `alts`, `gt`, `end_pos` and replacing `gt_alleles`'
  codes with the site's sequences (`source_gt` keeps them). Values change on array facets with such
  rows when the manifests are configured; without them nothing changes but the warning's wording.
  `RESOLVE_VERSION` 5.

- **No version change** (2026-10-23, RM269 with coverage): `sex_chromosome_evidence` takes an
  optional coverage facet and its result a nullable `coverage` field; it may now answer `XX` and `XO`.
  Not stored, so no frame changes.

- **No version change** (2026-10-23, RM278 coverage facets): no column of `facets`, `calls` or
  `blocks` changes. `read_coverage`, `alignment_coverage` and `callable_sites` return their own
  `Coverage` shape beside the frames ("Coverage facets", above).

- **No version change** (2026-10-26, RM282 SV concordance): no column of `facets`, `calls` or
  `blocks` changes. `concordance`'s own frames gain what an SV pair needs: `sites` adds `sv_kind`,
  `b_alt`, `a_end_pos`, `b_end_pos`, `size_similarity`, `reciprocal_overlap` and `ci_needed`, `key`
  adds `extent`, `summary` groups by `sv_kind` too, and `excluded` counts breakends as `breakend`.
  **Values change**: a symbolic allele, counted in `excluded` as `symbolic` before, is now a row of
  `sites`, so `sites`, `summary` and RM279's atlas (class `sv`) count it. Additive otherwise.

- **No version change** (2026-10-23, RM277 concordance): no column changes to `facets`, `calls` or
  `blocks`. `concordance` returns its own frames ("Concordance between two facets"), a measurement
  beside the facets, not a facet.

- **No version change** (2026-10-23, CGI native `callset_type`): no column changes. Rule changed
  (Decision 5, RM4): a `var` or `masterVarBeta` whose header has no whole-genome `#TITLE` is typed from
  its rows, where it was null. `absence_means` stays `blocks`. Values change on files that already
  load.

- **No version change** (2026-10-23, classes of CGI and key=value facets): no column changes. Rule
  stated (Decision 4): CGI's native files and the PGP key=value table fill `variant_classes`,
  `content_filter` and `content_filter_basis` as a VCF does, from the rows and the file name. They
  wrote null. Values change on files that already load; `absence_is_reference` is unchanged for them
  (their absence is `blocks` and `unknown`), but a consumer gating by class now has the classes.

- **No version change** (2026-10-23, the appended export's build): no column changes. Rule changed
  (Decision 8, "Two exports in one file"): an appended export's facet gets no `assembly_label` from
  the first export's preamble; null until `resolve_ref` infers it from the rows (F9). Values change on
  files that already load, where the first export states a build.

- **No version change** (2026-10-23, the empty allele): no column changes. Rule stated: a zero-length
  allele is `""` in every allele column (Decision 2). Values change on files that already load: an
  upload's ALT `-` (F12) and Genos's `-` were `"-"` and are now `""` (`source_gt` keeps the text);
  `""` is an `indel` in `variant_classes` and a site in RM4's count; `normalize` anchors CGI's native
  deletions as it anchored F12's (`NORMALIZE_VERSION` 5); `mint_vrs_ids` mints a resolved `""` ALT,
  which it withheld as `not_sequence` (`MINT_VERSION` 2); `write_vcf` refuses an empty REF as
  `ref_unresolved`, where it wrote it.

- **No version change** (2026-10-23, `normalize` on resolved array rows): no column changes. Rule
  changed in `normalize` ("Normalizing", step 2): in a facet read from a file, a row with `ref` set and
  no `source_ref` (an array row, a split Living DNA genotype) is aligned and trimmed, where it passed
  through as unchecked. Normalized facets get a new `facet_id` (`NORMALIZE_VERSION` 4).

- **An enum member, landed in 0.0.11** (2026-10-23, RM2's CGI question): `AbsenceMeans.REFERENCE`
  (`reference`), "every position the file omits was called reference by its producer", just-vcf's own
  term. Adding a member changes the `pl.Enum` dtype of `facets.absence_means`, which is a schema
  change, which landed in 0.0.11. Rule change: a CGAPipeline VCF typed `wgs` gets `reference` where it got `blocks`, and
  `absence_is_reference` answers True for it (Decision 4). Values change on files that already load.

- **No version change** (2026-10-23, `dropped_fields`): no column changes. Rule stated: `[]` means the
  loader read the file's fields and dropped none, null that the format has no notion of fields
  (Decision 5). 23andMe-layout tables and a Genos file without `##Columns=` wrote null; they now write
  `[]`. Values change on files that already load.

- **No version change** (2026-10-22, F4): no column changes. Rule added in `resolve_ref` (Decision 2):
  a Living DNA run-together indel genotype is split against the reference, so `gt_alleles` is set on
  rows where it was null, and `ref`, `alts`, `gt`, `end_pos` where one allele alone is the reference.
  A split row with `gt_alleles` set and `ref` null is a new combination for that format (an `I`/`D` row
  already had it). `RESOLVE_VERSION` 4.

- **No version change** (2026-10-22, RM6: lifted rows normalized): no column changes. Rule changed in
  `normalize` ("Normalizing", step 2): rows of a `liftover_to` facet are left-aligned and trimmed on
  the target, where they passed through as unchecked. A lifted row that moves follows the moved-row
  rules; a swapped one keeps no `gl`, `gp` or `ad` when it moves, since its `source_alleles` are not
  its new `[ref, alt]` (counted, as for any row). Normalized facets get a new `facet_id`
  (`NORMALIZE_VERSION` 3).

- **No version change** (2026-10-22, RM4's marker overlap): no column changes. `chip` gains two values,
  `illumina_gsa` and `illumina_omniexpress` (an open vocabulary stored as `str`), and is now inferred for
  FTDNA (`ftdna`, `ftdna_space`) and MyHeritage-style (`rsid_result_csv`) exports and for Final Reports
  that name no manifest, where it was always null. Values change on such facets.

- **No version change** (2026-10-22, F3): no column changes. Rule added: a second export appended to a
  23andMe-layout file in the two-allele layout is its own facet, `relation = appended_export`, a new
  value of the open `relation` vocabulary (Decision 8, "Two exports in one file").

- **Rule change, still 0.0.10** (2026-10-22, F15): a sequence ALT is `ACGTN` bases only, in the
  variant-class rule (`variant_classes`) and in the callset's site count (RM4). An ALT of other letters
  (`CNV`, an ambiguity code) gets no class and is no site, counted in a load warning. Values change on
  files with such ALTs; no column changed, so no version bump.

- **No version change** (2026-10-22, RM6's minus-strand indels): no column changes. Rule changed in
  `liftover` ("Lifting", step 3): minus-strand indels are re-anchored and lifted, where they were left
  behind; `indel_minus_strand` now names only those that land at a target contig's first base. Lifted
  facets get a new `facet_id` (`LIFT_VERSION` 2).

- **No version change** (2026-10-22, RM6's gVCF output): no column changes. Rule added: writing
  blocks back out (Decision 4).

- **Three facet columns, landed in 0.0.10** (2026-10-22, RM4's chip detection and F13):
  `facets.chip`, `facets.chip_basis` and `facets.strand_check`, all nullable and additive, which by
  the rule above is a minor bump, 0.0.10. Rule change in `resolve_ref`:
  in a `mixed` facet, array rows that are not strand-safe keep `ref` with `alts` and `gt` null
  (Decision 2), a new null combination on rows that were filled before.

- **0.0.22** (2026-10-27, F23): `RefCheck.UNSTATED`, a member of `calls.ref_check`'s enum (the entry
  at the top of this list).

- **0.0.21** (2026-10-26, the manifests that read REF): `facets.resolve_inputs` (`str` or null), a
  digest of the chip manifests `resolve_ref` read an array's rows with, inherited by derived facets;
  on a consensus, a digest of its parents' REF states. Additive. A store of 0.0.20 migrates with it
  null, named as not recorded on its resolved array facets and its consensus facets. Rule change in
  the store: `REF_STATE` gains it, so the same array resolved without and then with manifests, or a
  consensus rebuilt over parents in other REF states, is refused (`StoreConflictError`) where it was
  reported `present`.

- **0.0.20** (2026-10-26, the rule that read REF): `facets.resolve_version` (`u32` or null), the
  `RESOLVE_VERSION` that read the facet's REF, inherited by derived facets. Additive. A store of
  0.0.19 migrates with it null, named as not recorded on its resolved facets. Rule change in the
  store: re-adding a held facet in another REF state is refused (`StoreConflictError`) where it was
  reported `present`.

- **0.0.19** (2026-10-24, CGI's allele scores): `calls.hq` (`list[u16]` or null). Additive.

- **0.0.18** (2026-10-23, the PL/GL source): `facets.gl_source` (`gl`, `pl`, `gl_over_pl` or null). Additive.

- **0.0.17** (2026-10-23, the CGI leftovers): `PhaseSetKind.CGI_HAPLINK` (the `phase_set_kinds` dtype gains a
  member); CGI masterVarBeta fills `ad`; CGI `var` 1.x layouts tested.

- **0.0.16** (2026-10-23, layer 2's method chain): `facets.method_steps`, `reference`, `aligner`,
  `aligner_version`. Additive; `phasing_method` now comes from the phase step.

- **0.0.15** (2026-10-23, PSL phase sets): `calls.phase_set_ordinals` (PSO) and `phase_set_qualities` (PSQ);
  PSL read into `phase_sets` with kind `psl`. Additive.

- **0.0.14** (2026-10-23, RM282): `calls.sv_type`, `cipos`, `ciend`, `copy_number`, `alt_copy_numbers`.
  Additive. Values change: Canvas copy-neutral `ALT .` records are calls, not blocks (F18); `symbolic` in
  VRS withholding now covers only CGI no-calls, SVs read `sv_no_vrs_model`.

- **0.0.13** (2026-10-23, RM275): a fourth frame, `links` (`Loaded.links`, `schema.LINKS`). Additive.

- **0.0.12** (2026-10-23, RM273): `facets.subject`, the caller's subject label. Additive.

- **0.0.11** (2026-10-23): `AbsenceMeans.REFERENCE` (a CGI VCF's silence is its producer's reference call;
  the `absence_means` enum dtype gains a member). Values change: the empty allele is `""` in every loader,
  and `dropped_fields` is `[]` wherever a loader read the fields and dropped none.

- **0.0.10** (2026-10-22, RM4's chip detection and F13): `facets.chip`, `chip_basis` and `strand_check`.
  Additive. In a `mixed` facet an array row may have `ref` set while `alts` and `gt` are null.

- **0.0.9** (2026-10-22, RM4's filtered genomes): `facets.variant_classes`,
  `content_filter` and `content_filter_basis`, and the per-class reading of absence (Decision 4) with
  `absence_is_reference`. Additive; no value of an existing column changes.

- **0.0.8** (2026-10-22, RM6's third slice): `facets.derivation_inputs`, the files a derived facet
  was made with (a liftover's chain, by name and sha256). Additive. `liftover` emits the `liftover_to`
  relation the facets table already listed; rules under "Lifting".

- **0.0.7, continued** (2026-10-22, RM271; no version bump: no column or dtype changed):
  `mint_vrs_ids` fills `vrs_allele_ids` and `vrs_errors`, and `vrs_withheld` computes why an id is
  missing. Substitutions now wait for `resolve_ref` like every allele; `no_seqcol` left the reasons
  (a read `ref` implies its collection); `not_sequence`, `no_contig` and `refused` joined them;
  `telomere` covers POS 0 only (a span past the contig's end reads as `ref_unchecked`).

- **No version change** (2026-10-22, RM6's second slice): `normalize` emits the `normalized` relation
  the facets table already listed, and no column changes. Rules added: the normalizing rules above;
  on a row `normalize` moved or trimmed, `ref_check` and `source_ref` are null although `ref` is set
  (Decision 6's "null means no `source_ref`" covers it: the derived row has none).

- **Rule changes, still 0.0.7** (2026-10-22, F10 and the SV read side): GP values within [0, 1] that do
  not sum to 1 are probabilities under a 4.3+ header (Decision 3); a symbolic SV's `end_pos` reads INFO
  END (the `end_pos` rule). Values change on files that already load; no column changed, so no version
  bump.
- **Rule changes, still 0.0.7** (2026-10-22, RM6 follow-ups): PL becomes `gl` as `-PL/10` in Float64,
  rounded to Float32 once (F11); `JUST_VCF_SPLIT` is a split marker (`split_signal`); ALT `-` is an
  unanchored deletion that `normalize` anchors (F12); `normalize` drops `<*>` from split rows and
  judges a duplicate without the columns a move clears. Values change; no column changed, so no
  version bump (`NORMALIZE_VERSION` is 2).

- **0.0.7** (2026-10-22, RM3): `facets.ref_seqcol_digest` names the collection
  `resolve_ref` read REF from. Additive. `resolve_ref` fills `ref`, `ref_check`, and an array row's
  `alts` and `gt` (Decision 2).

- **0.0.6** (2026-10-22, RM2's VCF side): `facets.ploidy_statements` (the header's ploidy lines
  verbatim); `vendor`, `caller`, `caller_version`, `phasing_method` now filled from the header; BCF rows
  carry no `source_gt`; CGI haploid no-calls are blocks; an AD/PL/GL/GP vector that does not fit its row
  is withheld; a freeBayes+VarScan2 merge loads as two `caller_split` facets. Additive.

- **0.0.5** (2026-10-22, RM4): `facets.sorted_name_length_pairs_digest`, the header's coordinate
  system as a seqcol attribute digest. Additive. `assembly_label` (basis `inferred`), `callset_type`
  and `absence_means = weak_reference` are now filled for VCF facets where the evidence settles them.

- **0.0.4** (2026-10-19, RM270, from S2): `calls.split_signal` and `facets.split_marker` flag rows
  that look like one multi-allelic site written one row per ALT. Both additive; no value changed.

- **0.0.3** (2026-09-30, found by loading the 124 PGP-HMS VCFs): `gq` widens to `u16` (22 files had GQ
  above 255); an `ALT .` record with `INFO/END` is a block (Illumina gVCF; a Veritas file had 183 million
  such records land in `calls`), and its `end_pos` is `END`; GP's scale is read from decisive values rather
  than from the header version.
- **0.0.2** (2026-09-30, found by the first loader over 503 PGP-HMS array files): `gt_alleles` may be
  null when the genotype text cannot be split without the site's alleles, and `phased` is null as a
  whole when the source says nothing about phase. Both loosen a "never null"; no column changed.

## At a glance

- Three frames: `facets` (one row per facet), `calls` (one row per site per facet), `blocks` (reference
  and no-call intervals, for formats that have them). A fourth, `links`, holds the facet graph's edges
  (from, to, relation, provenance), so a facet can have several parents (Decision 12).
- Coordinates: `pos` and `end_pos` are 1-based and inclusive. VRS objects keep 0-based interbase coordinates
  internally, and nothing else does.
- REF: `ref` only from the assembly sequence; `source_ref` is what the file wrote; `ref_check` compares
  them. VRS ids wait for `ref` to be read, and are withheld on a mismatch.
- Genotype: literal `gt_alleles` always; `gt` indices into the file's own allele list; `source_gt` is
  the raw text.
  Vectors are indexed over `source_alleles` and never reordered.
- Evidence: `gl` (likelihoods), `gp` (posteriors) and `dosage` never mixed, and a row can have more
  than one. `evidence_basis` is the origin only (`assay`, `ref_block`, `absence`). Array calls have no
  vector until concordance can measure one.
- Identity: `vrs_allele_ids` parallel to `[ref] + alts`, minted by `mint_vrs_ids` once REF is read.
  `vrs_errors` stores only what minting reported; why an id is missing is computed by
  `vrs_withheld()`. `seqcol_digest` per facet. SV alleles have no id (`sv_no_vrs_model`); what the
  caller wrote about the event is `sv_type`, `cipos`, `ciend`, `copy_number`, `alt_copy_numbers`.
- Absence: a facet property (`absence_means`), turned into a row only on request.
