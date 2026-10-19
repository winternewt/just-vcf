# Roadmap history

Shipped roadmap items, moved here from [ROADMAP.md](ROADMAP.md) with their rationale. Nothing is
deleted; it is relocated.

## RM270 — Pre-split multi-allelic input
Severity: high · Status: shipped 2026-10-19: B, the flags (schema 0.0.4), and D, `rejoin_split`; C not built, the person chose D (desk #45) · Owner: lead

Motivating report: S2 (just-dna-lite, 2026-10-18). Many VCFs arrive already split one row per ALT, in
the `bcftools norm -m-` convention, and their headers do not say so. In a split row the `0` of a
genotype can stand for "another ALT", not REF: a `1/2` at `A>G,T` becomes `1/0` at `A>G` and `0/1` at
`A>T`, and a `1/1` becomes `1/1` plus a `0/0` row for the uncarried ALT.

Reproduced 2026-10-18 against the code, on a synthetic VCF: the loader keeps one `calls` row per record
(SCHEMA Decision 1, A), so the split rows land as written, and `gt_alleles` maps each `0` onto
`source_ref`. The split `1/2` loads as `["G","A"]` and `["A","T"]`, two calls with a reference copy the
person does not have; the uncarried ALT's `0/0` row loads as `["C","C"]`, a hom-ref call. `warnings` is
empty. The unsplit record of the same site loads as `["G","T"]`.

**Prior art** (read from source and docs, 2026-10-18). bcftools `norm -m-` writes the other ALT as `0`
by default, keeping GT slots, and `.` with `--multi-overlaps .` (1.17+); its `--old-rec-tag` INFO
(1.12+) is opt-in. Hail's `split_multi_hts` downcodes to `0` and sorts unphased calls, so its split
`1/2` is `0/1` + `0/1`, and its `a_index`/`was_split` are row fields an exported VCF lacks. GATK
`LeftAlignAndTrimVariants --split-multi-allelics` writes `0`, except that every sample becomes `./.` when
any sample is `1/2`. vt `decompose`, Glow and vcflib `vcfbreakmulti` write `.`, and vt and Glow add INFO
`OLD_MULTIALLELIC`. The VCF spec says nothing about split records. `bcftools norm -m+` re-joins slot by
slot (`1/0` + `0/1` gives `1/2`). Rows written with `.` already load as half-calls and are safe.

**Detection**, three signals per `pos` within a contig, reported separately and never folded into one
flag: a **marker** (INFO `OLD_MULTIALLELIC`, or a value of the `--old-rec-tag` shape
`CHROM|POS|REF|ALT|IDX`); the **slot pattern** (two or more non-reference rows, each GT slot non-zero in
at most one); and an **all-zero row of any ploidy beside a non-reference row**. On the facet, a marker
is Kleene: true when seen, null otherwise, never false, because a split file can show no signal at all
(Hail's sorted output, split indels that left-alignment moved to another POS). Measured 2026-10-18 on
four unsplit private WGS VCFs (DeepVariant ×2, DRAGEN, GATK HaplotypeCaller), 22.2 million calls: no
signal fired. Only the DRAGEN file had positions with two records (10), none matching. The rate on a
real split file is not measured (no bcftools here); a synthetic split VCF is caught at both its sites.
The two unsplit DeepVariant gVCFs gave no signal either: no two call rows share a position once their
reference blocks go to `blocks`, which never take part.

**What the frame says about a `0` it cannot place.** No option changes `gt`, the row grain, or
`source_gt` where it is kept (null on VCF 4.4+, P4). Not decided:

| | What the frame says | Legality | Cost |
|---|---|---|---|
| A. Warn only | a warning per signal, with a count | patch | every flagged row stays wrong for a reader that skips warnings |
| B. A + flags | `facets.split_marker` (Kleene) and a nullable row column `split_signal` naming the signal | minor | values stay wrong, but every flagged row can be found |
| C. B + slot-aware withhold | in a flagged group, a `0` becomes null in `gt_alleles` only in a slot where another row of the group is non-zero; a slot zero in every row keeps its REF | values read today change (pre-1.0, changelog says so) | a false positive loses a true REF copy, never states a false one; amends Decision 2: `gt` and `gt_alleles` agree except where withheld, and a null means uncalled or withheld |
| D. Re-join as a derived facet | each flagged group merged into one site slot by slot, as `bcftools -m+` does (`relation = normalized`, charter 5) | minor | the most work; GT only, since bcftools subsets PL/AD to REF + one ALT; phasing across slots may be lost |

`gl`, `gp` and `ad` on a flagged row stay as written: they compare REF with one ALT and inherit the
ambiguity, and the row flag is what tells a reader so. The lead recommended B and C now and D inside
RM6 (desk #45). **Decided 2026-10-19**, the person: "4 if it can be implemented w/o bam, else 3". It
can: the re-join needs only the split rows, as `bcftools norm -m+` does; reads would only rebuild the
full likelihood vectors, which is variant calling (charter). So D is built now, GT only, beside B;
C is not built.

**RM6's own split output** must not write this convention: it writes `.` for the other ALT, as vt,
Glow, vcflib and `bcftools --multi-overlaps .` do.

Legality: detection, a warning and a new nullable field are additive (minor). Nulling `gt_alleles` for
detected rows changes the meaning of values already emitted.

RM6 carries point 3: its own split output writes `.` for the other ALT.

## RM8 — Output frame schema
Severity: high · Status: signed off 2026-09-30 as schema v0.0.1 · Owner: user decides, agent drafts

Before any per-format loader is written, draft a proposal for the one frame every loader emits, and get
the user's sign-off. The charter (docs/CHARTER.md) is still an uncommitted draft, so the schema cannot
be derived from it alone.

The proposal must settle:
- **Columns and dtypes** of the minimal shape (many formats in, one frame out), and which of them the
  extended shape (facets per assembly, caller, phasing) adds.
- **Coordinate naming.** 1-based `POS` at every public boundary; name each column for the convention it
  holds so 0-based sources (CGI, Color) cannot leak through unconverted.
- **Graded genotype evidence.** The user's words: a fuzzy scale from 0 (unknown) to PHRED-like
  (confident), probably an n-vector of probabilities or log-likelihoods per position; never a binary
  called/not-called. Say how array calls, imputed dosages (Oxford `.gen`, Genes for Good imputed),
  half-calls (`1|.`) and CGI partial no-calls (`GC?`) map onto it.
- **Absence semantics.** Where an absent site means unknown (arrays, imputed content, VCFs of array
  data) versus weak evidence of reference (whole-genome callsets); likely a per-facet property.
- **Provenance per row:** source file, format id, caller or vendor (the `source` detector in RM2),
  build as stated versus inferred.
- **GA4GH fit:** VRS allele ids as the join key (withheld when non-conformant), seqcol for assembly
  identity.
- **Multi-caller files:** PGP merges carry the caller in INFO/SF, not per sample column
  (docs/PGP_HMS.md).

Inputs: the "Contract gaps" list at the end of docs/FORMATS.md, RM2, docs/POLARS_BIO_PENDING_FIXES.md
(P1: undeclared INFO/FORMAT keys are dropped), and the charter draft's invariants.

Signed off 2026-09-30 as docs/SCHEMA.md, schema v0.0.1. The user's words: "lets keep schema as it is
now and switch to building - we'll find out if some pieces are missing in progress". Gaps found while
building the loaders are fixed in SCHEMA.md as they come up, each with a version bump.


## RM7 — PGP-HMS sample leecher
Severity: medium · Status: shipped 2026-09-28 · Owner: Claude

A script in `scripts/` that downloads consented-public 23andMe exports from
<https://my.pgp-hms.org/public_genetic_data?data_type=23andMe> into `data/input/`, so tests and
dogfooding have real array files that are allowed to be used. httpx is the dev dependency for it.

Shipped as `scripts/pgp_leech.py`.
Download links answer two ways: older uploads are served directly, newer ones redirect to an Arvados
collection page. The listing's participant cell can carry the original enrolment number
(`huF06AD0, PGP113`), which the live test caught. Rows over 200 MB are skipped by default: the
23andMe listing totals about 101 GB, and 19 of its files are multi-GB uploads.
