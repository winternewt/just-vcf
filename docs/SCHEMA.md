# Frame schema: proposal (RM8)

Status: **signed off 2026-09-30 as schema v0.0.1; now v0.0.7.** Every per-format loader (RM2) emits the frames
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
sequence, not by a choice. The alternative, a derived facet (`relation = ref_checked`), was rejected:
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
- Left unread, each reason counted once: I/D codes, whole no-calls (the probe's REF span is unknown),
  genotypes that could not be split, spans past the contig's end or at POS 0, and every row of a
  contig that could not be matched with the reference (below).

**Matching the reference to the facet.** Per contig the calls use, read through `canonical_chrom` on
both sides (so `chr1`/`1` and `chrM`/`MT` meet; array code `25` reads X):

- the facet stated the contig's length (`header_contigs`): the reference's length must equal it, and
  its md5 too where the header gives one;
- it did not (arrays): the facet's `assembly_label` must equal the reference's. The reference's label
  comes from its contig 1 length (GRCh37, GRCh38, T2T-CHM13), or from one the caller asserted when
  adding the FASTA. `GRCh37` cannot vouch for MT, since it covers hg19's chrM (16,571 bases) as well
  as rCRS (16,569), so an array's MT stays unread.

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
  id. The `?` marks the bases CGI could not call.
- **The `*` allele** (VCF: "allele missing due to overlapping deletion") is kept as `"*"` in `alts` and
  `gt_alleles`. It says the haplotype is covered by a deletion called at another record. It is not
  a reference call, and it gets no VRS id. `0/2` over `ALT=C,*` must never be read as a het with a
  reference haplotype and nothing else.
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
| **B. Separate columns per quantity** | `gl: list[f32]`: log10 likelihoods (from PL as `-PL/10`, or GL directly). `gp: list[f32]`: posterior probabilities (GP, Oxford `.gen` triples). Both are in VCF genotype order over `source_alleles` (Decision 2). Which of them is non-null says which quantities the source gave, and a row can have more than one. A consumer has to handle two columns, and the schema is honest about what each one is |
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

Complete Genomics' VCF (`vcf_cgi`, or any VCF declaring its alleles) adds two block kinds, both
`no_call`: a `<CGA_NOCALL>` record spans POS to INFO `END`, and, where the header declares
`CGA_NOCALL`, a record with ALT `.`, no `END` and no allele called (`./.`, `.|.`, or a haploid `.`,
which polars-bio 0.36 returns as a null GT) spans its REF. Its
`<CGA_CNVWIN>` records are 2 kb copy-number windows, not genotypes, and are left out with a count. A
half-called reference span (`0|.`, ALT `.`) stays a call.

A 4.5 block states its length per sample in FORMAT `LEN`. INFO `END` is deprecated and only a fallback
when `LEN` is missing. Blocks keep the source's `MIN_DP` as `dp_min`, and not the block's average `DP`.
A block with an average depth of 25 can still hold an uncovered base.

**Recommendation: B.** The charter's "one frame out" is amended to name the three frames (facets,
calls, blocks): a gVCF cannot fit in one calls frame. Every format lands in the same `calls` schema, and
a format that has blocks also fills `blocks`.

**Absence is a facet property, not rows.** The facet carries `absence_means`: `unknown` (arrays,
imputed content, VCFs of array data), `weak_reference` (whole-genome callsets without blocks), or
`blocks` (gVCF, CGI: only positions inside a `ref` block are reference). Nothing is materialised
until a consumer asks for a position. When it does, the answer comes out as a `calls`-shaped row with
`evidence_basis = absence` or `ref_block`.

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
| `sample` | the sample column name, verbatim. **Identifying**: PGP uploads carry real names here. Logs, warnings and error messages name the `facet_id`, never the sample |
| `assembly_label`, `assembly_label_basis` | e.g. `GRCh37`; basis is `stated` (header), `inferred` (contig lengths) or null. A stated label is kept when the lengths disagree, with a warning (RM4, below) |
| `seqcol_digest` | refget Sequence Collection digest (charter 3). Null until computed; the label never stands in for it. No VCF header can give it (below) |
| `sorted_name_length_pairs_digest` | the seqcol `sorted_name_length_pairs` attribute digest of the header's `##contig` names (verbatim) and lengths: the coordinate system, not the collection. Null when a contig has no length or there are none |
| `ref_seqcol_digest` | the Sequence Collection `ref` was read from (RM3, `resolve_ref`). Null until REF is read. It is kept apart from `seqcol_digest`: the reference that checks the facet is not evidence of the facet's own collection, and filling one from the other would make that check agree with itself |
| `callset_type` | `wgs`, `wes`, `panel`, `array`, `imputed`, or null (RM4, below). A VCF gets `array` or `imputed` where its header says what the content is (Illumina GenotypingToVCF, a 23andMe conversion, Gencove, Minimac, Genes for Good: the `source` detector, RM2); otherwise RM4 infers it from the rows |
| `absence_means` | Decision 4 |
| `phasing_method` | null unless stated |
| `derived_from`, `relation` | parent facet and how: `liftover_to` (Cat-VRS), or a documented just-vcf relation (`caller_split`, `merge`, `normalized`, `split_rejoined`) |
| `contig_map` | `list[struct{source, chrom}]`: how this file's contig names were mapped (Decision 7) |
| `header_contigs` | `list[struct{name, length, md5}]` from `##contig` lines. It is evidence for the assembly (RM4) and for the seqcol digest, not the digest itself |
| `vcf_version` | from `##fileformat`. It decides how `GP` is scaled (Decision 3) |
| `sample_meta` | the file's `##SAMPLE` and `##PEDIGREE` lines for this sample, verbatim. Pedigree becomes GA4GH Pedigree later (charter) |
| `dropped_fields` | INFO/FORMAT keys present in the file but not carried (Decision 9). For a delimited table, the file's own column names that no frame column carries (a Final Report's `GenCall Score`) |
| `split_marker` | the declared split-marker INFO key(s), comma-joined, when any record was tagged as split off a multi-allelic site (vt's `OLD_MULTIALLELIC`, a bcftools `--old-rec-tag` key). Null = no record was tagged, which is not "not split": see `split_signal` |

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
ALT or none; rows whose ALTs are all symbolic do not count, so a genome-wide SV callset is not `wgs`.
Anything between the bands is null. A gVCF is not a callset type: its blocks set
`absence_means = blocks`, whatever its scope. Otherwise `wgs` gives `weak_reference` and every other
type, null included, `unknown`. Measured on PGP-HMS uploads: four Genos exomes at 44,020–44,796 sites
and 13.6% breadth, about 400 Mb reached.

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
| `alts` | `list[str]` | the source did not say (a non-VCF row before RM3). `[]` = a stated monomorphic record (ALT `.`). Symbolic alleles and `*` kept verbatim, except `<NON_REF>`, which is written `<*>` |
| `source_alleles` | `list[str]` | the source gave no vector to index (Decision 2) |
| `vrs_allele_ids` | `list[str]` | parallel to `[ref] + alts`, gnomAD's convention. An element is null when it was withheld (charter 4); the whole list is null when nothing could be minted. Minting rule below ("When an id is minted") |
| `vrs_errors` | `list[str]` | parallel to `vrs_allele_ids`: what minting itself reported for an allele it tried and could not mint (`not_justifiable`, `no_contig`, …). Null where minting succeeded or never ran. Only facts that minting alone learns go here; see "Why an id is missing" below |
| `gt_alleles` | `list[str]` | the genotype text cannot be split into haplotypes without the site's alleles (Living DNA writes an indel's two allele sequences back to back, `CCA`); `source_gt` keeps the text. Entries null per uncalled haplotype |
| `gt` | `list[u8]` | there is no allele list to index yet. Indices into `[source_ref] + alts` when the source wrote a REF, otherwise into `[ref] + alts` after RM3 (Decision 2) |
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

1. `marker`: the row carries INFO `OLD_MULTIALLELIC` (vt, Glow), or a value with five `|` fields under
   a key declared with bcftools' `--old-rec-tag` Description. A record bcftools only realigned has
   four and is not flagged.
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
   A site holding a symbolic allele (`<*>`, `<DEL>`, a breakend) is left as written: a split `<*>` row
   would read back as a reference block. `*` is split like any allele.
2. **Left-align and trim**, only on biallelic sequence rows with `ref_check = match`: drop a common last
   base, and when an allele runs empty, prepend the reference base before both and step `pos` back;
   repeat; then drop common first bases while both alleles keep one. At position 1 an allele keeps
   its last base. On a row this changed, `ref` is the new REF (read from the reference), `gt_alleles`
   and `source_alleles` follow the new spelling (same order, so the vectors stay indexed), and
   `source_ref`, `ref_check`, `source_gt`, `vrs_allele_ids` and `vrs_errors` are null: they described
   the representation the source wrote, at the position it wrote it. A `ref_check = match` here would be
   bases read from the sequence compared with that same sequence. Rows whose REF was not shown to
   match pass through, counted, and so do rows whose ALT is no VCF allele (a `-` for a deletion).
3. **Dedup and sort.** Rows equal in every column are kept once, counted. Rows that state the same
   allele at one position and differ in anything else are all kept, with a warning (charter 8). Then
   canonical order (Decision 11). Because a moved row loses `source_ref`, `ref_check` and `source_gt`,
   one indel written twice by the source, once already left-aligned and once further right, stays as
   two rows and is reported as a conflict, not dropped as a duplicate.

Rows the loader flagged as one record of a pre-split site (RM270) pass through as written: run
`rejoin_split`, then `resolve_ref`, then `normalize`. Blocks pass through unchanged.

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

**When an id is minted.** A VRS Allele id for a substitution does not encode REF. It needs only the
contig's refget accession, the position and the allele, so it is minted without reading any bases.
An indel id needs the sequence around it (VRS normalises it against that sequence), so it waits for
bases. Per allele:

- **Substitution** (REF and allele the same length): minted as soon as the facet's assembly is
  identified and the contig has an accession, from a small table of contig accessions per assembly.
  No sequence is read. `ref_check = not_run` does not block it. `ref_check = mismatch` withholds it,
  because a mismatch is evidence that the position or the build is wrong.
- **Indel or MNV with a length change**: minted only after bases have been read and
  `ref_check = match`.
- **Array rows**: minted once RM3 has derived `alts`, under the same two rules.

So loading a VCF needs no reference download. The table of contig accessions is derived by gtars from
each assembly's sequence collection and checked by a test, never typed in by hand. Where the bases come
from when they are needed (a local store, or a hosted one read lazily) is RM3's decision; the
measurements are in docs/ROADMAP.md, RM3.

**Why an id is missing** is not a stored column. Most of the reasons are already facts in the frame,
and a stored copy of them drifts: RM3 fills `ref`, and a stored reason still says "unchecked". So the
reason is computed by one function, `vrs_withheld()`, a lazy Polars expression over the frame. It
returns one reason per allele, parallel to `vrs_allele_ids`, and null where an id exists. It checks in
this fixed order and the first hit wins:

1. **What the allele is**, read from the allele string: `symbolic` (`<DEL>`, `<CNV:TR>`, breakends),
   `unspecified` (`<*>`), `spanning_deletion` (`*`), `array_code` (23andMe `I`/`D`), `partial_call`
   (CGI `?`).
2. **Pipeline state**, read from other columns: `ref_unresolved` (`ref` and `source_ref` both null),
   `ref_unchecked` (`ref_check = not_run`, indels only; substitutions mint without the check),
   `ref_mismatch` (`ref_check = mismatch`), `no_seqcol`
   (`facets.seqcol_digest` null), `telomere` (`pos` 0, or past the contig's length).
3. **What minting reported**: the `vrs_errors` element.
4. Otherwise `not_run`: nothing forbids minting and it has not happened yet.

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
INFO/SF lists its index. With a column per caller, each facet reads its own column. With one data
column, the column holds the first caller's genotype wherever that caller found the site, so the
VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the shared sites it
lost. `derived_from` stays null: the merged file is not loaded as a facet of its own (that is option
A), so there is no parent row to point at. Not built: sourcing from the `-orig` key=value table,
which has no loader yet.

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
  not exist yet.
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
- **A GRCh38 REF that is `N`.** A site inside the gap has no base to compare, and `ref_check` would
  compare against `N`. Whether gtars mints on an `N` reference has not been checked. If it refuses,
  the refusal lands in `vrs_errors`.

None of this changes a column today. It marks where the person-level layer and the coverage facet have
to connect.

## Open questions

Q1 (CGI `GC?`), Q2 (array code 25), Q3 (withheld VRS ids), Q4 (sample names), Q5 (`end_pos`) and Q6
(PSL keys) were settled on 2026-09-30 and are folded in above.


## Changes since sign-off

- **No version change** (2026-10-22, RM6's second slice): `normalize` emits the `normalized` relation
  the facets table already listed, and no column changes. Rules added: the normalizing rules above;
  on a row `normalize` moved or trimmed, `ref_check` and `source_ref` are null although `ref` is set
  (Decision 6's "null means no `source_ref`" covers it: the derived row has none).

- **Rule changes, still 0.0.7** (2026-10-22, F10 and the SV read side): GP values within [0, 1] that do
  not sum to 1 are probabilities under a 4.3+ header (Decision 3); a symbolic SV's `end_pos` reads INFO
  END (the `end_pos` rule). Values change on files that already load; no column changed, so no version
  bump.

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
  and no-call intervals, for formats that have them).
- Coordinates: `pos` and `end_pos` are 1-based and inclusive. VRS objects keep 0-based interbase coordinates
  internally, and nothing else does.
- REF: `ref` only from the assembly sequence; `source_ref` is what the file wrote; `ref_check` compares
  them. Substitution ids mint without bases and are withheld on a mismatch; indel ids wait for a match.
- Genotype: literal `gt_alleles` always; `gt` indices into the file's own allele list; `source_gt` is
  the raw text.
  Vectors are indexed over `source_alleles` and never reordered.
- Evidence: `gl` (likelihoods), `gp` (posteriors) and `dosage` never mixed, and a row can have more
  than one. `evidence_basis` is the origin only (`assay`, `ref_block`, `absence`). Array calls have no
  vector until concordance can measure one.
- Identity: `vrs_allele_ids` parallel to `[ref] + alts`. `vrs_errors` stores only what minting
  reported; why an id is missing is computed by `vrs_withheld()`. `seqcol_digest` per facet.
- Absence: a facet property (`absence_means`), turned into a row only on request.
