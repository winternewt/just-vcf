# Frame schema: proposal (RM8)

Status: **draft, awaiting the user's sign-off.** Nothing here is implemented. Every per-format loader
(RM2) emits the frames described below, so the loaders wait for this document to be settled.

The charter (docs/CHARTER.md) is itself still a draft, so each decision below cites the charter
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
facet (charter 5), not something the loader does silently.

## Decision 2: how the genotype is written

Array exports have no REF when they are loaded (charter 10: REF never comes from the genotype), so a
VCF-style index genotype (`0/1`) cannot be assigned until RM3 has resolved REF.

| Option | Cost |
|---|---|
| A. VCF indices only (`gt: list[u8]` into `[ref] + alts`) | Impossible for arrays until RM3. The loader would have to either block on a FASTA or invent a REF |
| B. Literal alleles only (`gt_alleles: list[str]`, one entry per haplotype) | Works for every input at load time. But likelihood vectors are indexed by VCF genotype order, so they cannot be read without `[ref] + alts` next to them |
| **C. Both: `gt_alleles` always, `gt` when REF is known** | Two columns that must agree: `gt_alleles == ([ref] + alts)[gt]`, enforced by a test. `gt` is null while `ref` is unresolved. Filling it is RM3's job, and it never writes REF from `gt_alleles` |

**Recommendation: C.**

**`source_alleles`: what the vectors are indexed over.** A likelihood or posterior vector needs an
allele order, and REF is not always known when the row is loaded. An Oxford `.gen` row names `alleleA`
and `alleleB` and gives three posteriors, but does not say which allele is the reference. So each row
carries `source_alleles: list[str]`, the allele list in the order the source indexed its vectors. For a
VCF that is `[ref] + alts`. For `.gen` it is `[alleleA, alleleB]`. For an array hard call with no
vector it is null. `gl`, `gp` and `ad` (Decision 3) are indexed over `source_alleles` and are never
reordered in place. When RM3 resolves REF and it turns out to be `alleleB`, `ref` and `alts` are filled
and the vectors stay as the source wrote them. A consumer that wants `[ref] + alts` order derives it,
and that derivation gets a test.

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
- **CGI partial no-call inside an allele** (`GC?`): see open question Q1.
- **The `*` allele** (VCF: "allele missing due to overlapping deletion") is kept as `"*"` in `alts` and
  `gt_alleles`. It says the haplotype is covered by a deletion called at another record. It is not
  a reference call, and it gets no VRS id. `0/2` over `ALT=C,*` must never be read as a het with a
  reference haplotype and nothing else.
- **A monomorphic record** (ALT `.`) has `alts = []`, which is a stated "no alternate allele". A null
  `alts` means the source did not say, which only happens for non-VCF sources (an array row before
  RM3). ALT is a required VCF column, so polars-bio's `""` always means `.` and becomes `[]`.
- **`<NON_REF>`** is written as `<*>` in `alts`. VCF 4.5 makes them aliases.
- **`source_gt`**: the source's genotype text, verbatim, on every row that had one. It costs one string
  column. In exchange, every derived column can be checked against what the file said, and no
  representation choice here ever loses information.

Phasing is per haplotype, not per row. Since VCF 4.4 each allele in GT carries its own phase
indicator (`|0/1/2`: the first allele is phased and the other two are not), and `PSL` gives each allele
its own phase set. So:

- `phased: list[bool]`, parallel to `gt_alleles`. An element is null where the source did not say
  (an array export).
- `phase_sets: list[str]`, parallel to `gt_alleles`. It comes from `PSL` directly, since PSL names are
  global to the file. A `PS` is scoped to its chromosome ("on the same chromosome and have the same PS
  value"), so it is stored as `PS:<chrom>:<n>` on every phased haplotype. A phased haplotype with no PS
  or PSL belongs to the chromosome's one implicit set (VCF 4.5: "All phased genotypes that do not
  contain a PS subfield are assumed to belong to the same phased set") and gets `PS:<chrom>:implicit`.
  Unphased haplotypes are null. PS and PSL names never share a set (VCF 4.5: "PS and PSL are not
  interoperable"), and the prefix keeps them apart.

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
| **B. Two vector columns plus a basis** | `gl: list[f32]`: log10 likelihoods (from PL as `-PL/10`, or GL directly). `gp: list[f32]`: posterior probabilities (GP, Oxford `.gen` triples). Both are in VCF genotype order over `source_alleles` (Decision 2). `evidence_basis` says which source field each came from. A consumer has to handle two columns, and the schema is honest about what each one is |
| C. Keep native columns (`PL`, `GL`, `GP`, `DS`) as-is | No normalisation at all. Every consumer then re-learns every caller's dialect |

**Recommendation: B.** Scalars that callers commonly give sit beside it: `gq: u8`, `dp: u32`,
`ad: list[u32]` (parallel to `source_alleles`), all nullable.

Sources for `gp`: `GP`, and `PP` converted as `10^(-PP/10)`. Both are posteriors, so this is a change
of scale, not a change of quantity. One trap: before VCF 4.3, `GP` was phred-scaled. The loader reads
the `##fileformat` version and checks the values against it (a probability vector sums to about 1). If
the two disagree, `gp` stays null and the facet says why. The loader never guesses the scale.

Local alleles (VCF 4.5 `LAA` with `LPL`, `LAD`, `LGP`) fit without a new column. `source_alleles` is
then `[ref] + alts[LAA]`, which is exactly the set the local vectors are indexed over.

`dp` and `ad` come from FORMAT (this sample), never from INFO (all samples combined). INFO/DP and
FORMAT/DP are different quantities with the same name. The loader resolves each key by the header's
namespace and never by the column name polars-bio gives it.

A dosage (VCF `DS`, the expected count of the second allele) is one number. It cannot be turned into a
posterior triple, because two unknowns do not follow from one number, so it gets its own column:
`dosage: f32`, with `gp` left null unless the source also gave GP.

`evidence_basis` is an enum: `likelihood` (PL/GL), `posterior` (GP, `.gen`), `dosage`, `hard_call` (a
genotype with no vector, which is every consumer array), `ref_block` (restored from a gVCF or CGI
reference block), `absence` (restored from a whole-genome callset's silence, charter 11).

**What an array hard call gets** is the user's open "tbd". Options:

| Option | Cost |
|---|---|
| **A. No vector. `evidence_basis = hard_call`, `gl`/`gp` null** | Honest now. Consumers that need a number have none until C lands |
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
| **B. Two frames per facet: `calls` and `blocks`** (the Hail VDS layout) | A second frame to know about. `blocks` is small in columns: `facet_id, chrom, pos, end, kind (ref / no_call), gq, dp_min, gl`. `gl` is the block's own likelihoods over `[REF, <*>]` (a gVCF block carries PL, e.g. `0,60,900`), because that *is* the graded reference evidence (charter 7). It is null only where the source gave none (CGI blocks) |
| C. Expand blocks to one row per base | Billions of rows. Rejected |

A 4.5 block states its length per sample in FORMAT `LEN`. INFO `END` is deprecated and only a fallback
when `LEN` is missing. Blocks keep the source's `MIN_DP` as `dp_min`, and not the block's average `DP`.
A block with an average depth of 25 can still hold an uncovered base.

**Recommendation: B.** "One frame out" is then read as "one schema out": every format lands in the same
`calls` schema, and a format that has blocks also fills `blocks`.

**Absence is a facet property, not rows.** The facet carries `absence_means`: `unknown` (arrays,
imputed content, VCFs of array data), `weak_reference` (whole-genome callsets without blocks), or
`blocks` (gVCF, CGI: only positions inside a `ref` block are reference). Nothing is materialised
until a consumer asks for a position. When it does, the answer comes out as a `calls`-shaped row with
`evidence_basis = absence` or `ref_block`.

## Decision 5: the facet table and provenance (charter 5, 6)

Provenance that is constant across a file lives once in a `facets` frame. Each row of `calls` and
`blocks` carries only `facet_id`, plus the provenance that actually varies per row (`ref_basis`,
below).

`facets` columns, all nullable unless marked:

| Column | Meaning |
|---|---|
| `facet_id` (required) | deterministic: a digest of the source bytes, the archive member, the sample and the load options. Loading the same file twice gives the same id (charter 12) |
| `source_sha256`, `source_path`, `member` | what was read. The path is a label, not identity |
| `format_id`, `family` | from `sniff` |
| `vendor`, `caller`, `caller_version` | the RM2 `source` detector. Null = not identifiable, never guessed |
| `sample` | the sample column name as written in the file. It can be a real name; see Q4 |
| `assembly_label`, `assembly_label_basis` | e.g. `GRCh37`; basis is `stated` (header), `inferred` (contig lengths) or null |
| `seqcol_digest` | refget Sequence Collection digest (charter 3). Null until computed; the label never stands in for it |
| `callset_type` | `wgs`, `wes`, `panel`, `array`, `imputed`, or null (RM4) |
| `absence_means` | Decision 4 |
| `phasing_method` | null unless stated |
| `derived_from`, `relation` | parent facet and how: `liftover_to` (Cat-VRS), or a documented just-vcf relation (`caller_split`, `merge`, `normalized`) |
| `contig_map` | `list[struct{source, chrom}]`: how this file's contig names were mapped (Decision 7) |
| `header_contigs` | `list[struct{name, length, md5}]` from `##contig` lines. It is evidence for the assembly (RM4) and for the seqcol digest, not the digest itself |
| `vcf_version` | from `##fileformat`. It decides how `GP` is scaled (Decision 3) |
| `sample_meta` | the file's `##SAMPLE` and `##PEDIGREE` lines for this sample, verbatim. Pedigree becomes GA4GH Pedigree later (charter) |
| `dropped_fields` | INFO/FORMAT keys present in the file but not carried (Decision 9) |

No load timestamp: it would break determinism, and the file digest already says which bytes were read.

Row-level: `ref_basis` is an enum: `source` (the file stated REF), `fasta`, `panel`, `unresolved`. A
source-stated REF is later *checked* against the FASTA by RM3, and never filled from it.

## Decision 6: the `calls` columns

With Decisions 1–5 as recommended:

| Column | Dtype | Null means |
|---|---|---|
| `facet_id` | `str` | never null |
| `chrom` | `str` | never null (Decision 7) |
| `pos` | `u32` | never null. 1-based, VCF `POS` |
| `end` | `u32` | the length of the REF span is unknown (a 23andMe `I`/`D` row, a Living DNA run-together indel before RM3). 1-based, inclusive. It is VCF 4.5's computed END: the largest of `pos + len(ref) - 1`, `pos + SVLEN` for a symbolic SV allele, and `pos + LEN - 1` for `<*>`. It is `pos` for a single-position array call. The loader computes it and does not take polars-bio's `end` for SVs (P3) |
| `ids` | `list[str]` | the source had no id. rsids, 23andMe `i` ids and VCF `ID` split on `;` |
| `ref` | `str` | unresolved (charter 10) |
| `ref_basis` | `Enum` | never null |
| `alts` | `list[str]` | the source did not say (a non-VCF row before RM3). `[]` = a stated monomorphic record (ALT `.`). Symbolic alleles and `*` kept verbatim, except `<NON_REF>`, which is written `<*>` |
| `source_alleles` | `list[str]` | the source gave no vector to index (Decision 2) |
| `vrs_allele_ids` | `list[str]` | parallel to `[ref] + alts`, gnomAD's convention. An element is null when it was withheld (charter 4); the whole list is null when nothing could be minted |
| `gt_alleles` | `list[str]` | never null on a genotyped row; entries null per uncalled haplotype |
| `gt` | `list[u8]` | REF unresolved. Indices into `[ref] + alts` |
| `source_gt` | `str` | the source wrote no genotype text |
| `phased` | `list[bool]` | parallel to `gt_alleles`. An element is null where the source did not say |
| `phase_sets` | `list[str]` | parallel to `gt_alleles`. An element is null for an unphased haplotype. `PS:<chrom>:<n>`, `PS:<chrom>:implicit` or a PSL name (Decision 2) |
| `gl` | `list[f32]` | no likelihoods given |
| `gp` | `list[f32]` | no posteriors given |
| `dosage` | `f32` | no dosage given |
| `evidence_basis` | `Enum` | never null |
| `gq` | `u8` | not given |
| `dp` | `u32` | not given |
| `ad` | `list[u32]` | not given |
| `qual` | `f32` | not given. It is the confidence in what ALT asserts: that the variant is there on a variant record, and that no variant is there on a monomorphic record (`alts = []`, where VCF defines it as `-10·log10 P(variant)`). The same threshold therefore works on both |
| `filters` | `list[str]` | the source said `.` or has no filter column. `["PASS"]` is literal. Site-level |
| `sample_filters` | `list[str]` | FORMAT `FT`: whether *this sample's* genotype passed, separate from the site. Null = not given |

`end` is not named for its convention, although charter 13 asks that. The reason: a 1-based
inclusive end and a 0-based half-open end are the same number, so only the start differs between the
two systems, and `pos` already names its convention. It is still a naming decision, so Q5 puts it to
the user.

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

**Recommendation: B.** One catch: array code `25` is the X pseudo-autosomal region. Mapping it to `X`
drops the PAR fact. Q2 asks what to do with it.

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

Frames are emitted sorted by `facet_id`, then contig rank, then `pos`, `end`, `ref`, `alts`. Contig rank
follows the seqcol order when it is known. Otherwise it is `1`–`22`, `X`, `Y`, `MT`, then the remaining
names sorted as strings. Nulls sort last. Every loader gets a test that loads a shuffled input and
compares it with the sorted one.

## Other VCF 4.5 points the loaders have to respect

None of these changes a column. Each is a rule a loader follows, and each gets a test.

- **REF may be a reduced IUPAC code.** The spec requires an ambiguous reference base to be written as
  the alphabetically first concrete base (`R` becomes `A`). When RM3 checks a source REF against a FASTA
  that holds `R`, a mismatch there is expected and is not a coordinate error.
- **POS 0 and N+1 mark telomeres.** `pos` allows 0. No VRS id is minted for either.
- **For a symbolic SV, POS is the base before the event.** For `<*>` the block includes POS. The `end`
  rule above follows both.
- **GQ is conditional** on the site being variant. It is not the confidence of a reference call.
- **Percent encoding.** `%3B`, `%3D`, `%2C` and the rest can appear in INFO values and must be decoded.
  Whether polars-bio decodes them has not been checked yet.
- **Ploidy is whatever GT says.** GT can have one allele (haploid Y, male X, MT) or three and more, and
  the genotype order for `gl`/`gp` follows the spec's general formula for any ploidy. Without GT, the
  spec says to assume diploid.
- **`ID` is a `;` list**, split into `ids`.

## Open questions

- **Q1. CGI partial no-calls** (`GC?`: some bases of an allele called, the rest not). Options: the
  haplotype becomes null and `source_gt` keeps the text (simple, loses the partial); or the allele is
  kept verbatim with `?` and no VRS id (keeps it, but `?` is not a base any consumer knows). Leaning
  towards the first.
- **Q2. Array code `25` (X PAR).** Map to `X` and record it in `contig_map` (the position then says
  whether it is PAR), or keep `XY` as its own contig label.
- **Q3. Withheld VRS ids.** Should a reason travel per row (`vrs_withheld: Enum`: symbolic, ref
  unresolved, no seqcol, not justifiable), or only as aggregate counts on the facet? Per row costs a
  column. Aggregates follow the house rule on repeated warnings.
- **Q4. Sample names.** PGP uploads can carry real names in the sample column. Keep them verbatim in
  `facets.sample` (it is what the file says), or store a digest and let the caller opt in to the name.
- **Q5. `end` or `end_pos`.** Keep polars-bio's `end` (the value is the same under both coordinate
  systems), or rename it `end_pos` so that every coordinate column says its convention (charter 13).
  Leaning towards `end`: one fewer rename at the polars-bio edge.

## At a glance (all recommendations taken)

- Three frames: `facets` (one row per facet), `calls` (one row per site per facet), `blocks` (reference
  and no-call intervals, for formats that have them).
- Coordinates: `pos` and `end` are 1-based and inclusive. VRS objects keep 0-based interbase coordinates
  internally, and nothing else does.
- Genotype: literal `gt_alleles` always; `gt` indices once REF is known; `source_gt` verbatim.
  Vectors are indexed over `source_alleles` and never reordered.
- Evidence: `gl` (likelihoods), `gp` (posteriors) and `dosage` never mixed; `evidence_basis` says which. Array
  calls are `hard_call` with no vector until concordance can measure one.
- Identity: `vrs_allele_ids` parallel to `[ref] + alts`, null where withheld; `seqcol_digest` per facet.
- Absence: a facet property (`absence_means`), turned into a row only on request.
