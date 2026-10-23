# Layers of the extended shape: progress and blockers

The person's goal of 2026-10-23: build unattended, todos first, then take repeated swings at layers 1–8
until each is at 90% or every layer is blocked by a decision only the person can make. Each blocker is
written here, under its layer, the moment it is met, and the work moves to another layer. The lead
seat keeps this file; percentages are the lead's estimate against what the layer would need to be
usable for one person's facets, and say what they rest on.

| # | Layer | Estimate | Rests on |
|---|---|---|---|
| 1 | Reading facets in | 80% | 25 of 39 formats, SV/CNV, PSL phase sets, BCF stated builds; remaining formats being swung at; F1 and F17 blocked |
| 2 | Provenance per facet | 75% | as before plus scaffold aliases; the method chain (aligner, post-processing) is being built; coverage building |
| 3 | Identity (comparing facets) | 92% | REF, VRS, assembly inference incl. GRCh36/CHM13, scaffold aliases (F19: a PGP GRCh38 VCF reads 100%); `_fix` patches need a patch store |
| 4 | A person | 65% | subject, Person with links, the durable store; relatedness (default check) building |
| 5 | Agreement measured | 15% | concordance building; atlas and posterior wait for it |
| 6 | Facet graph and liftover | 90% | links frame, liftover, normalize, rejoin_split |
| 7 | Graded evidence across sources | 30% | values carried as written; calibration waits for the atlas |
| 8 | Export | 92% | VCF/gVCF with tabix, VCF 4.5 PSL/PSO/PSQ, CisPhasedBlock export naming what withheld it, VCF Zarr |

## Blockers waiting for the person

Each names its layer, the decision, the options, the lead's recommendation, and what the work does
meanwhile. The desk carries the ones filed there.

- **Layer 1, desk #83:** how to open Deflate64 zips (an optional extra, an external program, or refuse).
  Meanwhile they are refused with an error that says why.
- **Layer 1, F17:** Living DNA indel REF needs the Sirius chip's site alleles, and no manifest is public
  or on disk. Needs a source the person may have, or stays withheld.
- **Layer 4, desk #72:** who is who in felis, ksuhaster and newton. Meanwhile those folders are never
  treated as one person's facets.
- **Layer 5/7, desk #80:** merge duplicate rows that differ only in record id. Meanwhile both stay.
- **Layer 3, open:** mint SVs as VRS 2 CopyNumber/Adjacency objects, or keep withholding (RM282).
  Meanwhile withheld (`sv_no_vrs_model`).
- **Layer 8, open:** offer CisPhasedBlocks per window, so one symbolic allele does not withhold a whole
  chromosome's block (RM284). Meanwhile whole blocks.
- **Layer 4, desk #73:** may Privé's ancestry files be redistributed as a fixture. Meanwhile local only.
- **Layer 8, a reading to veto:** when a 4.4+ record's first slot carries a PSL name, the loader reads it
  as phased (`phased[0] = true`); VCF 4.5 writes an unphased allele as `.`. The alternative is null, which
  withholds every such block as `phase_unknown`. Built as the first; reversible (SCHEMA Decision 2).
- **Layer 3, open:** lifting onto GRCh38 `_fix`/`_alt` patch scaffolds needs a store holding patches (an
  NCBI or UCSC full analysis set). Meanwhile 40,099 HG002 rows stay `target_contig_unmatched`.
- **Tooling, open:** ruff as a dev dependency (about 200 findings in the tree).
