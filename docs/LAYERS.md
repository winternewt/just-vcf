# Layers of the extended shape: progress and blockers

The person's goal of 2026-10-23: build unattended, todos first, then take repeated swings at layers 1–8
until each is at 90% or every layer is blocked by a decision only the person can make. Each blocker is
written here, under its layer, the moment it is met, and the work moves to another layer. The lead
seat keeps this file; percentages are the lead's estimate against what the layer would need to be
usable for one person's facets, and say what they rest on.

| # | Layer | Estimate | Rests on |
|---|---|---|---|
| 1 | Reading facets in | 75% | 25 of 39 formats load, plus SV/CNV; chip detection; F1 and F17 open |
| 2 | Provenance per facet | 70% | caller, version, vendor, callset type, chip, content filter, ploidy lines, subject, derivation inputs; no aligner or post-processing steps; coverage facet building |
| 3 | Identity (comparing facets) | 85% | REF from a reference, VRS ids equal to vrs-python, assembly inference incl. GRCh36 and CHM13; UCSC scaffold aliases building (F19); seqcol digest withheld by design |
| 4 | A person | 45% | subject and Person; relatedness and the store building |
| 5 | Agreement measured | 15% | concordance building; atlas and posterior wait for it |
| 6 | Facet graph and liftover | 90% | links frame, liftover, normalize, rejoin_split |
| 7 | Graded evidence across sources | 30% | values carried as written; calibration waits for the atlas |
| 8 | Export | 85% | VCF and gVCF with tabix, CisPhasedBlock export; VCF Zarr building with the store; PSL phase sets unread |

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
- **Tooling, open:** ruff as a dev dependency (about 200 findings in the tree).
