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
Severity: high · Status: open · Owner: unassigned

Port from just-prs `arrays.py`, `io_utils.py` and `vcf.py`: 23andMe, AncestryDNA, MyHeritage and FTDNA
parsing; zip/gzip by magic bytes; VCF build detection. New work: a Living DNA parser (check whether Genotek exports differ from 23andMe; the local one carries a 23andMe header),
build detection for arrays, I/D indel codes, hemizygous calls. Leave out the observed-allele-as-REF
encoding (charter invariant 3).

Detect the format from content, never from the label. PGP-HMS shows labels are wrong in both
directions (docs/PGP_HMS.md). The 21 free-text uploads its classification could not settle, and the
multi-facet participants, are the first test inputs.

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
