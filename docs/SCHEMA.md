# Frame schema: proposal (RM8)

Status: **signed off 2026-09-30 as schema v0.0.1; now v0.0.11.** Every per-format loader (RM2) emits the frames
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
- null: there is no `source_ref` to compare (arrays).

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
same calls. A facet already read is refused rather than read again.

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
- Left unread, each reason counted once: I/D codes, whole no-calls (the probe's REF span is unknown),
  genotypes that could not be split, spans past the contig's end or at POS 0, and every row of a
  contig that could not be matched with the reference (below).

**Matching the reference to the facet.** Per contig the calls use, read through `canonical_chrom` on
both sides (so `chr1`/`1` and `chrM`/`MT` meet; array code `25` reads X):

- the facet stated the contig's length (`header_contigs`): the reference's length must equal it, and
  its md5 too where the header gives one;
- it did not (arrays): the facet's `assembly_label` must equal the reference's. The reference's label
  comes from its contig 1 length (GRCh36, GRCh37, GRCh38, T2T-CHM13), or from one the caller asserted
  when adding the FASTA. `GRCh37` cannot vouch for MT, since it covers hg19's chrM (16,571 bases) as
  well as rCRS (16,569), and `GRCh36` neither (NCBI36 and hg18 carry the 16,571-base chrM), so an
  array's MT stays unread.
- it stated neither lengths nor a label (an FTDNA export, a headerless table): the label is inferred
  from the rows (F9), and then matched as above. Per collection in the store, a deterministic sample
  of the facet's array SNV rows (at most 2,000 per chromosome on 1–22, X and Y, every k-th row by
  position; MT never) is scored: the share whose called alleles include the collection's base at that
  position, on the forward strand only. A row on a contig the collection lacks is not scored on it; a
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
minus-strand hom-alt call of an A/T probe, and the rows cannot say which. No genotype is complemented,
not even on a het that proves its probe minus: complementing exactly the minus probes needs the
vendor's own strand table.

A primary contig (1–22, X, Y) whose stated length or md5 differs refuses the whole facet
(`ReferenceMismatchError`): it is another assembly. Any other contig that does not match (hg19's chrM
against GRCh37's MT, an alt contig the reference lacks, a name the reference holds twice) leaves its
rows unread. A facet of which nothing matches is refused. When a store holds several collections, the
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
- **23andMe `I`/`D`**: kept verbatim as symbolic alleles. Resolving them needs chip manifests
  (LANDSCAPE). No VRS id is minted for them.
- **CGI partial no-call inside an allele** (`GC?`): kept verbatim in `gt_alleles` (`"GC?"`), with no VRS
  id. The `?` marks the bases CGI could not call. Because `gt` indexes `[source_ref] + alts`, a partial
  allele is in `alts` too (the native `var` and `masterVarBeta` files; CGI's VCF writes such a
  haplotype `.`).
- **No anchor base** (CGI native files): CGI writes an indel over its own span, without VCF's padding
  base, and REF only comes from the sequence. So an allele can be the empty string: a deletion of the
  whole span is `""` (`source_ref = "T"`, `alts = [""]` at `pos = begin + 1`). A locus of zero length
  (an insertion point) has `source_ref = ""`, sits at the base before it (`pos = begin`, where VCF puts
  the anchor) and has `end_pos` null until RM3 reads that base, as Genos insertion sites do. An empty
  string in `source_ref`, `alts`, `source_alleles` or `gt_alleles` occurs only in such unanchored rows,
  which is how RM3 finds the rows it has to pad. An insertion point CGI could not call stays a row with
  every haplotype null: it covers no base, so it cannot be a block.
- **The empty allele has one spelling: `""`** (`schema.EMPTY_ALLELE`), VRS's zero-length literal, in
  every loader: CGI's unpadded indels above, an upload's ALT `-` (below), and Genos's `-` (an insertion
  site's empty REF allele, or a deletion's ALT; `source_gt` keeps the `-/A` the file wrote). Every reader
  takes `""` as a sequence allele of length zero: it is an `indel` in `variant_classes`, a site in RM4's
  count, the deletion `normalize` anchors, and the allele `mint_vrs_ids` justifies (a resolved `CT>""` at
  31 gets the id of `GCT>G` at 30). `write_vcf` refuses it, REF or ALT, until it is anchored (VCF has
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
- `phase_set_kinds: list[Enum{ps, psl}]`, parallel to `phase_sets`: which field the value came from.
  VCF 4.5 says PS and PSL "are not interoperable", so a PS `5` and a PSL `5` are different sets.
- **The key of a phase set depends on its kind.** A PS is scoped to its chromosome ("on the same
  chromosome and have the same PS value"), so its key is `(facet_id, chrom, ps, value)`. `PS=100` on
  chr1 and on chr2 are two sets. A PSL name is global to the file ("records with different CHROM but
  the same phase-set name are considered part of the same phase set"), so its key is
  `(facet_id, psl, value)`, and one set can span chromosomes.
- A phased haplotype with no PS or PSL has `phased = true` and `phase_set` null (kind `ps`). It belongs
  to its chromosome's one implicit set (VCF 4.5: "All phased genotypes that do not contain a PS
  subfield are assumed to belong to the same phased set"), key `(facet_id, chrom, ps, null)`.

The VRS CisPhasedBlock is built from `phase_sets` at export. It is not stored per row.

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
| `vendor`, `caller`, `caller_version` | the RM2 `source` detector, from header lines only. Null = not identifiable, never guessed. `vendor` is the party the header names as delivering the file (`##center`, `##generatedby`), falling back to the format's own vendor (`vcf_cgi`); a caller's maker is not a vendor. `caller` is the program that wrote the genotypes; tools that only touched the file afterwards (GATK SelectVariants, ApplyVQSR, bcftools view) are not evidence. Two callers in one header leave both fields null, with a warning; so do two versions of one caller for `caller_version` |
| `ploidy_statements` | the header's ploidy lines verbatim, in header order (RM269 half 1): a meta line whose key names a ploidy or karyotype (DRAGEN `##referenceSexKaryotype=XY`, CGI `##source_MAX_PLOIDY`), and each ploidy option of a recorded command line, prefixed with the line it came from (`GATKCommandLine HaplotypeCaller: --sample-ploidy 2`, `bcftools_callCommand: --ploidy GRCh38`). An option whose value is a path writes `<path>` for it. No vocabulary: reading a model out of these is the karyotype work's job. `[]` = a VCF header with none; null = not a VCF. just-vcf's own shape |
| `sample` | the sample column name, verbatim; for a CGI native file its `#ASSEMBLY_ID`, which CGI's VCF uses as its sample column. **Identifying**: PGP uploads carry real names here. Logs, warnings and error messages name the `facet_id`, never the sample |
| `assembly_label`, `assembly_label_basis` | e.g. `GRCh37`; basis is `stated` (header), `inferred` or null. A stated label is kept when the lengths disagree, with a warning (RM4, below). `inferred` has three sources: a VCF's contig lengths at load (RM4), an array's own rows scored against the reference by `resolve_ref` (F9, Decision 2), and a lifted facet's target collection (`liftover`, "Lifting"). A consumer tells them apart by `relation = liftover_to` for the third, then by `header_contigs`, which only the first has |
| `seqcol_digest` | refget Sequence Collection digest (charter 3). Null until computed; the label never stands in for it. No VCF header can give it (below) |
| `sorted_name_length_pairs_digest` | the seqcol `sorted_name_length_pairs` attribute digest of the header's `##contig` names (verbatim) and lengths: the coordinate system, not the collection. Null when a contig has no length or there are none |
| `ref_seqcol_digest` | the Sequence Collection `ref` was read from (RM3, `resolve_ref`). Null until REF is read. It is kept apart from `seqcol_digest`: the reference that checks the facet is not evidence of the facet's own collection, and filling one from the other would make that check agree with itself |
| `callset_type` | `wgs`, `wes`, `panel`, `array`, `imputed`, or null (RM4, below). A VCF gets `array` or `imputed` where its header says what the content is (Illumina GenotypingToVCF, a 23andMe conversion, Gencove, Minimac, Genes for Good: the `source` detector, RM2); otherwise RM4 infers it from the rows |
| `absence_means` | Decision 4 |
| `variant_classes` | `list[struct{variant_class, alleles}]`: the ALT alleles of each class the rows carry, counted at load, in the order `snv`, `mnv`, `indel`, `sv`, a class with none left out. `snv` and `mnv` are sequence ALTs as long as their REF (one base, more); `indel` a sequence ALT of another length (the empty allele `""` included: an unanchored deletion, F12's `-`), whatever its length; `sv` a symbolic ALT or a breakend. `<*>`, `*`, `<CGA_NOCALL>` and ALT `.` are no class, and neither is an ALT of letters outside `ACGTN`, which is no VCF allele (Veritas's `CNV`, F15; counted in a load warning). `[]` = a facet whose rows carry REF and ALT and no variant allele; null = not counted (a genotype table, which writes no REF or ALT). just-vcf's own shape |
| `content_filter` | `snv_only`, `indel_only`, `sv_only`: the one class the rows were cut down to, or null when no evidence settles one. Decision 4 says what it does to absence |
| `content_filter_basis` | `rows`, `stated` (a header line records the filtering command) or `file_name` (only the name's words say so); null with a null `content_filter` |
| `phasing_method` | null unless stated |
| `derived_from`, `relation` | parent facet and how: `liftover_to` (Cat-VRS, written by `liftover`), or a documented just-vcf relation (`caller_split`, `merge`, `normalized`, `split_rejoined`, `appended_export`). `appended_export`, like `caller_split`, has no parent: a second export appended to a file after the first, in another layout (F3), loaded as its own facet beside the first's; `derived_from` stays null, since neither export was made from the other, and the two share `source_sha256` and `member` (below, "Two exports in one file"). A facet with several parents (`merge`) has `derived_from` null and its parents in `links` (Decision 12), which also links a caller split's facets and an appended export to the first; Decision 12 states how the two agree |
| `contig_map` | `list[struct{source, chrom}]`: how this file's contig names were mapped (Decision 7) |
| `header_contigs` | `list[struct{name, length, md5}]` from `##contig` lines. It is evidence for the assembly (RM4) and for the seqcol digest, not the digest itself |
| `vcf_version` | from `##fileformat`. It decides how `GP` is scaled (Decision 3) |
| `sample_meta` | the file's `##SAMPLE` and `##PEDIGREE` lines for this sample, verbatim. Pedigree becomes GA4GH Pedigree later (charter) |
| `dropped_fields` | INFO/FORMAT keys present in the file but not carried (Decision 9). For a delimited table, the file's own column names that no frame column carries (a Final Report's `GenCall Score`), or, where the file names none, the format's names for them; for CGI's native files and the key=value table, their columns and keys. `[]` = the loader read the file's fields and dropped none (every 23andMe-layout table: its four fields are all carried). Null = the format has no notion of fields; no loaded format is one today. A derived facet keeps its parent's |
| `split_marker` | the declared split-marker INFO key(s), comma-joined, when any record was tagged as split off a multi-allelic site (vt's `OLD_MULTIALLELIC`, `write_vcf`'s `JUST_VCF_SPLIT`, a bcftools `--old-rec-tag` key). Null = no record was tagged, which is not "not split": see `split_signal` |
| `derivation_inputs` | `list[struct{role, name, sha256}]`: the files a derived facet was made with besides its parent and its reference collection (`ref_seqcol_digest`). A liftover writes its chain (`role = chain`, the file's name, not its path, and the sha256 of its bytes). Null on a facet read from a file. just-vcf's own shape |
| `chip`, `chip_basis` | the genotyping chip an array export came from (RM4, `just_vcf.chip`), and how it is known: `stated` (the header names it: AncestryDNA's array version, a Final Report's `Content` manifest, Living DNA's genotype chip) or `inferred` (a 23andMe export's row count and internal-id profile both point at one generation; a preamble-less table in its layout is read so when it carries 23andMe's `i…` ids; an FTDNA or MyHeritage export, or a Final Report that names no manifest, whose rsids one Illumina manifest clearly holds, by marker overlap). An open vocabulary: `23andme_v1_v2` (the two are not told apart), `23andme_v3`, `23andme_v4`, `23andme_v5`, `ancestrydna_v1`, `ancestrydna_v2`, `living_dna_sirius`, `illumina:<manifest>` as a Final Report names it, and the platform families `illumina_gsa` and `illumina_omniexpress`, which only marker overlap names: a stated value carries the manifest, an inferred one the family, since one manifest per family is compared with, and a GSA v2 Final Report matches the v3 manifest at 0.97: the overlap does not tell versions apart. Null = unknown (a statement with no entry is quoted in a warning; an overlap that names no chip says why in one) or not an array. just-vcf's own shape |
| `strand_check` | `forward` or `mixed`: what `resolve_ref` found about the strand an array's genotypes are written on (F13, Decision 2). Null = not checked: not an array, REF not read yet, or too few heterozygous rows. Never read as forward |

No load timestamp: it would break determinism, and the file digest already says which bytes were read.

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
- `assembly_label` is inferred from primary-chromosome lengths against GRCh37, GRCh38 and T2T-CHM13
  (domain constants, `just_vcf.assembly`): at least three must match one assembly and none may
  match another or none. The mitochondrion only contradicts (all three share 16,569; hg19's 16,571
  counts as GRCh37). Fewer, or a mixture, leaves it null. A stated label (sniff's header hints) is
  never replaced; a disagreement is a warning.

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
| `end_pos` | `u32` | the length of the REF span is unknown (a 23andMe `I`/`D` row, a Living DNA run-together indel before RM3). 1-based, inclusive. Rule below |
| `ids` | `list[str]` | the source had no id. rsids, 23andMe `i` ids and VCF `ID` split on `;` |
| `ref` | `str` | the assembly sequence has not been read yet (charter 10). Only ever filled from the sequence |
| `source_ref` | `str` | the source wrote no REF (arrays). What the file wrote, as written |
| `ref_check` | `Enum` | there is no `source_ref` to compare. `not_run`, `match` or `mismatch` (Decision 2) |
| `alts` | `list[str]` | the source did not say (a non-VCF row before RM3), or RM3 withheld it: an array row of a `strand_check = mixed` facet whose genotype reads differently on the two strands keeps `ref` and has `alts` and `gt` null (Decision 2). `[]` = a stated monomorphic record (ALT `.`). Symbolic alleles and `*` kept verbatim, except `<NON_REF>`, which is written `<*>` |
| `source_alleles` | `list[str]` | the source gave no vector to index (Decision 2) |
| `vrs_allele_ids` | `list[str]` | parallel to `[ref] + alts`, gnomAD's convention: the REF slot first, then one per ALT in `alts` order. An element is null when it was withheld (charter 4); the whole list is null when nothing was minted. Minting rule below ("When an id is minted") |
| `vrs_errors` | `list[str]` | parallel to `vrs_allele_ids`: what minting itself reported for an allele it tried and could not mint (`no_contig`, `refused`). Null where minting succeeded or never ran. Only facts that minting alone learns go here; see "Why an id is missing" below |
| `gt_alleles` | `list[str]` | the genotype text cannot be split into haplotypes without the site's alleles (Living DNA writes an indel's two allele sequences back to back, `CCA`; `resolve_ref` splits those where one split fits the reference, Decision 2); `source_gt` keeps the text. Entries null per uncalled haplotype |
| `gt` | `list[u8]` | there is no allele list to index yet, or RM3 withheld it with `alts` (a mixed-strand array, Decision 2). Indices into `[source_ref] + alts` when the source wrote a REF, otherwise into `[ref] + alts` after RM3 (Decision 2) |
| `source_gt` | `str` | the source wrote no genotype text (a BCF encodes GT in binary). The file's raw text, not polars-bio's rewrite of it |
| `phased` | `list[bool]` | the source says nothing about phase for the whole row (every array export). Otherwise parallel to `gt_alleles`, with an element null where the source did not say for that haplotype |
| `phase_sets` | `list[str]` | parallel to `gt_alleles`. The raw PS or PSL value. Null for an unphased haplotype, and for a phased one in its chromosome's implicit set. Keyed by kind: PS per chromosome, PSL per file (Decision 2) |
| `phase_set_kinds` | `list[Enum]` | parallel to `phase_sets`: `ps` or `psl`. Null where the haplotype is unphased |
| `gl` | `list[f32]` | no likelihoods given |
| `gp` | `list[f32]` | no posteriors given |
| `dosage` | `f32` | no dosage given |
| `evidence_basis` | `Enum` | never null. Origin only: `assay`, `ref_block`, `absence` (Decision 3) |
| `gq` | `u16` | not given, or outside 0..65535 (counted in the warnings) |
| `dp` | `u32` | not given |
| `ad` | `list[u32]` | not given |
| `qual` | `f32` | not given. It is the confidence in what ALT asserts: that the variant is there on a variant record, and that no variant is there on a monomorphic record (`alts = []`, where VCF defines it as `-10·log10 P(variant)`). The same threshold therefore works on both |
| `filters` | `list[str]` | the source said `.` or has no filter column. `["PASS"]` is literal. Site-level |
| `sample_filters` | `list[str]` | FORMAT `FT`: whether *this sample's* genotype passed, separate from the site. Null = not given |
| `split_signal` | `Enum` | no sign that the row is one record of a multi-allelic site split one row per ALT. Otherwise why it looks so: `marker`, `slot_pattern` or `all_zero_beside_alt` (rule below). In such a row a `0` can mean another ALT, not REF |

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
- `ids` and `filters` are unions (sorted; `PASS` only when no row failed a filter); `phased`, the phase
  sets, `dp` and `qual` are kept when every row agrees and null otherwise.
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
   `rejoin_split` can join it again. `source_gt` is null (no file wrote the text). Vectors:
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
   facet read from a file (no `derived_from`), a row with `ref` set and no `source_ref`: `resolve_ref`
   is the only writer of `ref` there, and it read that REF with nothing from the source to compare it
   with (an array row, a Living DNA genotype split against the reference; F4). Drop a common last
   base, and when an allele runs empty, prepend the reference base before both and step `pos` back;
   repeat; then drop common first bases while both alleles keep one. At position 1 an allele keeps
   its last base. On a row this changed, `ref` is the new REF (read from the reference), `gt_alleles`
   and `source_alleles` follow the new spelling (same order, so the vectors stay indexed), and
   `source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids` and `vrs_errors` are null: they described
   the representation the source wrote, at the position it wrote it. A `ref_check = match` here would be
   bases read from the sequence compared with that same sequence. Rows whose REF was not shown to
   match pass through, counted, and so do rows whose ALT is no VCF allele. A deletion written without
   its anchor base, ALT the empty allele (an upload's `-`, F12; CGI's native deletions), is anchored on
   the base before when matched (`CT>""` at 30 becomes `GCT>G` at 29) and then aligned like any
   deletion. At a contig's first base there is no base before, and the row stays as written, counted.
3. **Dedup and sort.** Rows that state the same thing are kept once, counted: equal in every column
   except the ones a move clears (`source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids`,
   `vrs_errors`), which say how the source wrote a row, not what it states. So one indel written twice
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
   both sides. Through `canonical_chrom`, the source sizes must equal the collection the parent's REF
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
   `vrs_errors` are null: they described the source's representation on the source assembly. `gt`
   indexes `[ref] + alts`. `gl`, `gp`, `ad` and the order of `source_alleles` are the parent's (Decision
   2: never reordered). A swapped row's `dosage` is null, since it counted the other allele. A row
   lifted to another contig than its source's has `phased`, `phase_sets` and `phase_set_kinds` null: a
   PS set is keyed per contig, and two sets would otherwise merge. `ids`, `filters`, `qual`, `gq`, `dp`
   and `evidence_basis` are kept. Re-mint ids on the target with `mint_vrs_ids`.
5. **The derived facet.** `ref_seqcol_digest` is the target collection; `assembly_label` is the
   target's label, basis `inferred` (from the collection's lengths, which the chain's target side
   matched); `seqcol_digest`, `sorted_name_length_pairs_digest` and `header_contigs` are null (they were
   the parent's header); `contig_map` maps the target collection's names of the contigs the rows use;
   `derivation_inputs` names the chain; `absence_means` is `unknown`, because the rows the lift dropped
   and the regions no chain covers are absent from it without being reference (the XYY stress case
   below).

Not lifted, each counted in the warnings and named in the report's `lift_reason` (an open vocabulary,
`UNLIFTED_REASONS`): `ref_unread`, `source_ref_mismatch`, `allele_not_sequence` (the empty allele of an unanchored deletion, which `normalize` anchors first; symbolic alleles, `*`,
`<*>`, breakends: an SV's extent is not lifted), `split_flagged` (run `rejoin_split` first),
`contig_not_in_chain`, `chain_contig_mismatch`, `no_chain_block`, `split_across_blocks`,
`multiple_targets`, `target_contig_unmatched`, `indel_minus_strand` (a minus-strand indel lifted to
the target contig's first base, with no base before it to re-anchor on) and `ref_mismatch_after_lift`. The `blocks` frame is not lifted; the
derived facet has none, and the warnings count the parent's.

**`end_pos`.** It is named apart from polars-bio's `end` because that value is wrong for symbolic SV
alleles (P3), so the two must not share a name. The rule follows VCF 4.5's computed END:

1. Take `abs(SVLEN)` for each allele (the spec says to read a negative SVLEN as positive). An SVLEN
   declared `Number=1` is read as the one-ALT list it is; one that is not one value per ALT cannot be
   matched to its alleles and is not read.
2. For each allele compute an end:
   - `<DEL>`, `<DUP>`, `<INV>`, `<CNV>` and their subtypes: the record's INFO `END` where it gives one,
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
source did not contradict it (`ref_check` is `match`, or null on an array row). Substitutions too:
`mint_vrs_ids` refuses a facet whose REF was never read. The rule the draft had, minting a
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
   (`*`), `symbolic` (`<DEL>`, `<CNV:TR>`, breakends, `.`-anchored single breakends), `partial_call`
   (CGI `?`), `array_code` (23andMe `I`/`D` where the source wrote no REF), `not_sequence` (anything
   else outside `ACGTN`, such as an ambiguity code in ALT).
2. **Pipeline state**, read from the row: `array_code` (an array row with an `I`/`D` genotype, whose
   REF stays unread), `telomere` (`pos` 0), `ref_unresolved` (`ref` and `source_ref` both null),
   `ref_unchecked` (a source REF, but `ref` not read: a contig the reference could not match, or a
   span past the contig's end), `ref_mismatch` (`ref_check = mismatch`).
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
| **B. One canonical spelling in `chrom`, the mapping recorded per facet in `contig_map`** | A mapping table to maintain. Canonical: `1`–`22`, `X`, `Y`, `MT`. `chrM`/`M`/`MT` are mapped explicitly, never with a `^chr` strip. Array codes `23`→`X`, `24`→`Y`, `26`→`MT`. Other contigs (alts, decoys, `chrUn_*`) keep their source name |
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
| `merge` | lineage, several parents | none yet | just-vcf's: `to` was made from every `from` linked to it (a merge, a consensus) |
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
- `derived_from` null and `relation = merge`: two or more `merge` links end at the facet;
- `derived_from` null and a kin `relation`: no lineage link ends at the facet, and a link of that
  relation names it;
- both null: no lineage link ends at the facet.

`facets.derivation_inputs` and a link's `provenance.derivation_inputs` differ by design: the facet's
is carried from parent to child (a facet normalized after a liftover still names the chain), the
link's names the files that one step read.

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
  X as one copy and Y as present, so `XY`: the Y depth of a variant-call VCF does not
  separate one Y from two, and XYY stays open until a coverage facet.
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
- `caller` and `ploidy_statements`, the facet's own columns, carried through unread.

It answers `x_copies` (`one`, `two_or_more`) and `y_present` separately, so a facet whose Y line is
unusable still says how many X copies it holds, and `inferred_karyotypic_sex`, a GA4GH Phenopackets
`KaryotypicSex` value (charter invariant 2: the vocabulary is the standard's, the field is
just-vcf's, since Phenopackets has no place for an inference). Only `XY`, `XXY` and
`UNKNOWN_KARYOTYPE` are estimated, the last always with a reason: `y_present` is `True` or unknown,
never `False`, because depth at called sites cannot show a Y absent (a man's GATK genome read below the
present band, close to the two-X range), so `XX` and `XO` wait for a coverage facet. The PAR bounds and the X-degenerate
Y come from GIAB's stratifications for GRCh37 and GRCh38 (docs/PROVENANCE.md); other builds answer
unknown. Rules and thresholds: `just_vcf/karyotype.py` and docs/ROADMAP.md RM269.

## The cache under the lazy frames (RM268)

The frames `load` and the derived steps return are lazy scans over files in `cache_dir()`
(`JUST_VCF_CACHE_DIR`). Not part of the schema, and no version bump; the rules are here because a
frame stops working when its files go.

- Entries are directories under `vcf/` (rewrite copies, kind `rewrite_copy`), `keyvalue/` and `cgi/`
  (the decompressed payload the key=value and Complete Genomics loaders read while they spill, kind
  `work_copy`, emptied when the load ends, failed or not) and `frames/` (spills: `frames_spill`,
  `rejoin_spill`, `resolve_spill`, `vrs_spill`, `normalize_spill`, `lift_spill`). Beside each is
  `<name>.json`, its manifest, and `<name>.lock`, an OS file lock the writer holds while it writes.
  The manifest names the root source (sha256 and resolved path), which every derived facet carries
  down from its parent, so a genome's derived spills go with its copies.
- An entry of an older version is removed when the step that owns its kind next runs (`load` for
  copies and load spills). Older means no version key ahead of the running version's and at least
  one behind, compared numerically (`0.0.9` is behind `0.0.11`); a key the running version has and
  the entry lacks counts as behind (the entry predates the key). An entry with any key ahead, a key
  the running version does not know, or a value that is not a number is a newer install's, or may
  be: kept. Two installs of different versions share the platform user cache, and the older one must
  not remove the newer one's live spills; the newer one retires the older one's. No current frame
  reads an old version, so this is safe. Rejected: sweeping every version other than the running
  one, which made two installs delete each other's spills on every load.
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
  and the facet row) and `result.json` in it listing its files and the step's row warnings. A run whose key is
  already there reuses that generation and writes nothing; a run on other rows under the same
  `facet_id` (a filtered `Loaded`) commits a second generation beside the first. No run rewrites or
  removes a complete generation, so a frame an earlier run returned, in this process or another,
  keeps working. A staging directory left by a writer that died, and a generation missing a file,
  are removed by the next run. Rejected: keying reuse on the facet id and versions alone, which would
  hand a filtered `Loaded` the full facet's rows; and replacing the files in place, which broke the
  earlier frame (`Invalid thrift` on its next `collect`). The manifest's `SPILL_LAYOUT` key (2) retires,
  through the sweep, entries written before generations. `liftover`'s generation holds two groups,
  the lifted `calls` and the `report`; its key also takes the chain's file name, the parent's block
  count and the two collections' labels, which its warnings name, and `Lifted.reasons` is read back
  from the report, so a reused generation answers as the run that wrote it. Its spill, once
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
- `clear_cache` removes what the caller names, even under a live `Loaded`; that frame's next
  `collect` then fails. It skips an entry whose lock is held, and touches nothing outside `vcf/`,
  `keyvalue/`, `cgi/` and `frames/`: the reference store under `reference/` is not a cache of a genome.
- Rejected: a size cap with least-recently-used eviction. `last_used` is when a step wrote or reused an
  entry, not when a frame last read it, so a cap could evict under a live frame without anyone asking.

## Open questions

Q1 (CGI `GC?`), Q2 (array code 25), Q3 (withheld VRS ids), Q4 (sample names), Q5 (`end_pos`) and Q6
(PSL keys) were settled on 2026-09-30 and are folded in above.


## Changes since sign-off

- **No version change** (2026-10-23, RM275, the `links` frame): a fourth frame, `Loaded.links`
  (Decision 12), and no column of `facets`, `calls` or `blocks` changes. `SCHEMA_VERSION` is not
  bumped here: the lead bumps it when merging. Rules added: every derived step adds a lineage link and
  carries its parent's; a load of several facets links them (`caller_split`, `appended_export`);
  `link_disagreements` states how `links` agrees with `derived_from` and `relation`.

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

- **An enum member, version bump pending** (2026-10-23, RM2's CGI question): `AbsenceMeans.REFERENCE`
  (`reference`), "every position the file omits was called reference by its producer", just-vcf's own
  term. Adding a member changes the `pl.Enum` dtype of `facets.absence_means`, which is a schema
  change; `SCHEMA_VERSION` was deliberately left at 0.0.10 on this builder branch, and the lead bumps it
  at merge. Rule change: a CGAPipeline VCF typed `wgs` gets `reference` where it got `blocks`, and
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

- **Three facet columns, version bump pending** (2026-10-22, RM4's chip detection and F13):
  `facets.chip`, `facets.chip_basis` and `facets.strand_check`, all nullable and additive, which by
  the rule above is a minor bump. `SCHEMA_VERSION` was deliberately left at 0.0.9 on this branch,
  which runs beside other builders' branches; the lead bumps it at merge. Rule change in `resolve_ref`:
  in a `mixed` facet, array rows that are not strand-safe keep `ref` with `alts` and `gt` null
  (Decision 2), a new null combination on rows that were filled before.

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
  `vrs_withheld()`. `seqcol_digest` per facet.
- Absence: a facet property (`absence_means`), turned into a row only on request.
