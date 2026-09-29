# Roadmap

Active items only, newest decisions first. Shipped items move to
[ROADMAP_HISTORY.md](ROADMAP_HISTORY.md). Claim a new id with
`INDEX=docs/ROADMAP.md ITEM_PREFIX=RM ./scripts/item-next.py --note "…"`.

## Items

## RM1 — Final shape and plan
Severity: blocker · Status: charter drafted, not signed · Owner: user

The shape as the user described it on 2026-09-28 is written into the uncommitted draft
`docs/CHARTER.md`. The frame is the product (Arrow/Polars, CompoSeq-style). There is a minimal
many-in/one-out shape and an extended shape of per-person facets (assemblies, callers, phasings,
coverage) with graded evidence, provenance, concordance and a graph. It is GA4GH first: VRS through
gtars, seqcol through refget. The user keeps it as a draft; planning of RM2–RM6 follows from it.
Surveys: docs/SIBLING_SURVEY.md, docs/LANDSCAPE.md.

## RM2 — Loading consumer formats
Severity: high · Status: in progress (detection shipped 2026-09-29, loaders next) · Owner: unassigned

**Shipped:** content-based detection. `just_vcf.sniff` and the declarative registry in
`just_vcf.formats` (docs/FORMATS.md). Over the PGP-HMS participants, 8 of about 1,740 files stay unknown,
all prose or empty.

**Next: per-format loaders on polars-bio**, blocked on RM8 (the frame schema), one frame out, driven by each `FormatSpec`. Start from the
contract gaps listed at the end of docs/FORMATS.md (versioned columns, section headers, 0-based CGI
conversion, half-calls, chromosome outside the row) and from polars-bio pending fixes P1/P2. Port from
just-prs `arrays.py`, `io_utils.py` and `vcf.py`; leave out the observed-allele-as-REF encoding (charter
invariant 3). Living DNA writes indels as run-together alleles (`CCA`); Color positions are 0-based.

The multi-facet participants in docs/PGP_HMS.md are the convergence tests: one person's facets, loaded
independently, must agree where they overlap.

**Also wanted: a `source` detector** naming the caller or vendor from VCF headers, because `Family`
cannot say "array or imputed content inside a VCF" and that decides whether absent means unknown.
Header evidence found while declaring the VCF formats:

| Header line | Caller / vendor |
|---|---|
| `##source=freeBayes…` + `##source=VarScan2` + `vcf-merge` | PGP freeBayes+VarScan2 merge |
| `##center=Complete Genomics`, `##source=CGAPipeline_*` | Complete Genomics |
| `##source=IsaacVariantCaller`, `##source=starling` | Illumina Isaac/starling (Veritas) |
| `##source=…CASAVA…makeConsensusVcf.pl` | Illumina CASAVA |
| `##workflow_type=Illumina_GenotypingToVCF` | Illumina array as VCF |
| `##DRAGENCommandLine`, `##source=DRAGEN_SNV` | DRAGEN |
| `##GATKCommandLine*`, `##source=HaplotypeCaller…`, `##UnifiedGenotyper=` | GATK |
| `##DeepVariant_version=` | DeepVariant |
| `##source=Platypus_Version_*` | Platypus |
| `##source=GenerateSVCandidates` | Manta |
| `##source=aengine` | FamilyTreeDNA Big Y |
| `##source=loimpute`, `##generatedby=Gencove` | Gencove low-pass imputation |
| `##source=Minimac3`, `##eagleVersion=` | Minimac3 + Eagle imputation |
| `##Genes_for_Good_version=` | Genes for Good |
| `##source=open_humans_data_processing.twenty_three_and_me` | Open Humans 23andMe conversion |
| `##source=intersection_dbSNP_23andme` | DNA.Land 23andMe conversion |
| `##pipeline=helix-*` | Helix |
| `##dataAnalysisProvider=Sequencing.com` | Sequencing.com |
| `##AnnotationEngineVersion=` + FORMAT `CGIANN_*` | Genos exome |

Not identifiable from the header: 23andMe's exome VCF, the Dante- and Nebula-named uploads (generic GATK).

## RM3 — REF lookup and padding
Severity: high · Status: open · Owner: unassigned

Port just-prs `reference_allele.py` (reference-panel tier, FASTA tier, FASTA build check). New work:
ALT derivation from genotype against REF, indel anchor bases, chr-prefixed FASTA support, emitting
reference-supported positions.

## RM4 — Input type inference
Severity: medium · Status: open · Owner: unassigned

Port just-dna-lite `detect_callset_scope` (WGS vs targeted), gVCF detection, and just-prs chip
generation detection. New work: WES vs WGS vs panel, and expanding gVCF `END` blocks.

## RM5 — Population inference
Severity: medium · Status: open · Owner: unassigned

Port the just-prs `ancestry/` package (PCA projection, KNN, proportions, consensus) behind the
`[ancestry]` extra, without the `PRSCatalog` wrappers or the Hugging Face pull.

## RM6 — Normalization and VCF writing
Severity: high · Status: open · Owner: unassigned

Reimplement chromosome naming and quality filters (just-dna-lite's versions are AGPL and not the user's,
so they are a behaviour reference, not a source) and port coordinate liftover from just-prs.
New work: left-align, split multi-allelic sites, strand resolution (A/T and C/G), dedup, sort, and a
writer that keeps FORMAT/GT (just-dna-lite's writer drops them).

## RM8 — Output frame schema
Severity: high · Status: drafted 2026-09-29 in docs/SCHEMA.md, awaiting sign-off; blocks the RM2 loaders · Owner: user decides, agent drafts

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

Inputs: the "Contract gaps" list at the end of docs/FORMATS.md, RM2, docs/polars-bio-pending-fixes.md
(P1: undeclared INFO/FORMAT keys are dropped), and the charter draft's invariants.
