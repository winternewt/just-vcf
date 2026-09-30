# Roadmap history

Shipped roadmap items, moved here from [ROADMAP.md](ROADMAP.md) with their rationale. Nothing is
deleted; it is relocated.

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
