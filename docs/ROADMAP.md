# Roadmap

Active items only, newest decisions first. Shipped items move to
[ROADMAP_HISTORY.md](ROADMAP_HISTORY.md). Claim a new id with
`ITEM_PREFIX=RM sh "${BLUEPRINTS_KIT:?}/bin/py.sh" "$BLUEPRINTS_KIT/classes/triage/item-next.py" --note "…"`; `CLAUDE.md` §8 says where `$BLUEPRINTS_KIT` comes from.

## Items

## RM272 — CLI
Severity: high · Status: open (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

Typer, registered in `[project.scripts]`, loading `.env` with python-dotenv (CLAUDE.md §4): `just-vcf sniff`, `load` (a summary of the facets), `convert` (to VCF/gVCF via `write_vcf`, optionally through `resolve_ref`, `normalize`, `liftover`, `mint_vrs_ids`) and `check`. Output is counts and aggregates only, never rows: the CLI is a privacy boundary like E12.

## RM273 — A person (subject) in the schema
Severity: high · Status: in progress (subject and Person shipped 2026-10-23; the default check is RM274's `same_person`) · Owner: lead

The person chose "Both: supplied, checked" (E1): the caller supplies a subject id when loading; nothing is inferred and facets are never regrouped by the library. The library checks the claim with relatedness (RM274) and warns on a mismatch. Most of the extended shape sits on this: RM275, RM276, RM277, RM280.

## RM274 — Relatedness and families
Severity: medium · Status: in progress (KING-robust relatedness, `same_person` as RM273's default check, the somalier wrapper and GA4GH Pedigree reading shipped 2026-10-23; left: somalier run for real once docs/TRUST.md has its row, pedigree compared with measured bands, a first-degree split by IBS0) · Owner: lead

The person chose both (E6): somalier (MIT) as an optional subprocess, skipped cleanly when absent, and just-vcf's own relatedness on the frames (IBS/KING-style over a validated site panel, so arrays work too). Pedigrees as supplied GA4GH Pedigree data. Feeds RM273's check.

## RM275 — The facet graph as a links frame
Severity: medium · Status: open (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose a small links frame (E4): (from, to, relation, provenance), so a facet can have several parents (a merge, a consensus). `derived_from` and `relation` stay on `facets` for the one-parent case. Relations verbatim from GA4GH where it has them (`liftover_to`), just-vcf's own otherwise, documented.

## RM276 — A durable per-person store
Severity: medium · Status: built 2026-10-23, awaiting merge (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose both (E13): a directory per person holding parquet per facet, the links frame and a manifest, read lazily; and a VCF Zarr export (its spec has no gVCF blocks, so blocks stay parquet-only). Facets accumulate; nothing is overwritten (invariant 5).

**Built:** `Store` (docs/SCHEMA.md, "The durable store"): parquet per facet, links per facet, a
manifest, a lock, lazy reads into `Loaded` and `Person`; and `write_vcf_zarr`, one facet's calls as
VCF Zarr 0.4, blocks counted and left in the store. **Left:** a migration between schema versions (a
store of another version is refused, not migrated); `cipos`, `ciend` and `alt_copy_numbers` in the
Zarr export; whether `Person` should carry the links (`store.loaded().links` has them).

## RM277 — Concordance between facets
Severity: medium · Status: open (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose graded where possible (E2): GL/GP compared where both facets carry them, hard genotypes otherwise, and every row says which. Kept as data, never resolved to a winner (invariant 8). Needs RM273 and VRS ids (RM271).

## RM278 — Coverage as a facet
Severity: medium · Status: in progress (built 2026-10-23: readers, BAM depth, callable) · Owner: lead

The person chose summaries and BAM depth (E3): mosdepth/samtools summaries, and polars-bio's `depth()` over a BAM the person has (scope note of 2026-09-28). Absence with a cause. RM269's XX and XO labels wait for this (the person chose to wait, 2026-10-23).

**Built 2026-10-23** (`just_vcf.coverage`; docs/SCHEMA.md "Coverage facets"): `read_coverage` for
mosdepth's four outputs and `samtools depth`/`coverage`, told by content; `alignment_coverage` over a
BAM or CRAM through polars-bio `depth()`, folded into 1,000-base windows with exact mean, min, max and
described bases, tiling the header's contigs with zeros; `callable_sites`/`callable_at` with
`min_depth` required. Decisions, each with what was rejected: a `Coverage` beside `Loaded`, not a
fourth frame (the load contract and the cache of every format would change for evidence about an
alignment); windows, not `depth()`'s raw runs (498 million on the 27 GB BAM below); per-base mode, not block mode (block mode's peak grows with the reads, below); a mean
bounds no base, so regions and summaries answer callable unknown (rejected: reading a mean of 30 as
"every base above 4"); `min_depth` with no default (rejected: GATK CallableLoci's 4, a convention the
caller should choose); a BAM's `facet_id` from its header, index and size (rejected: hashing the bytes,
a second full read); a headerless three-column file whose rows read as intervals is refused unless `kind=` states it (rejected: taking it as `samtools depth`, which reads a BED3 as depths). mosdepth's thresholds file is refused for now.

Substrate findings, polars-bio 0.36 (for the upstream-liaison seat to file; read in the mirror under
`data/interim/`, `bio-function-pileup` and `polars-bio/src/pileup.rs`):

- `depth()`'s output column `coverage` is Int16, written by `prev_cov as i16` from a 32-bit count, so a
  base deeper than 32,767 comes out negative and one past 65,535 as a wrong positive. Shown by a test
  (`test_a_depth_past_int16_is_left_out_and_counted`, 32,770 reads over five bases). Mitigated here:
  negative depths are left out of their window and counted in a warning.
- `DepthTableProvider::scan` drops pushed-down filters and `depth()` takes no region, and `PileupExec`
  holds a dense `Vec<i32>` per touched contig (its full length) for the whole file before emitting
  anything. Peak memory therefore grows with the genome, not with a window; it cannot be bounded from
  just-vcf.
- In block mode `merge_dense_results` builds every contig's run-length batches into `pending_batches`
  before emitting the first, so the peak adds all the runs (498 million on the 27 GB BAM below) to
  the dense arrays. Per-base mode emits contig by contig. `alignment_coverage` uses per-base mode.
- Block mode loses depth at a contig's end: on a synthetic BAM whose reads reach a contig's last
  base it emits no run for that contig at all, and on the 27 GB BAM the last window of 60 contigs
  read a zero there that per-base mode did not (the only differences across 2.94 million covered
  windows). Shown by `test_a_read_on_a_contigs_last_base_is_counted`, which fails on block mode.
- An indexed full scan reads the BAM twice: the dedicated no-coordinate partition scans "from the
  first alignment record and filter[s] to no-coor rows" (`bio-format-bam/src/physical_exec.rs`), the
  comment says, rather than seek past the placed reads. Measured as bytes read equal to twice the
  file (below).

**Measured 2026-10-23** on two private GRCh38 bwa-mem2 BAMs on the Krypton share (aggregates only;
the machine shared with about eleven other sessions, load average 17–60):

| BAM | depth mode | time | peak RSS | bytes read | windows |
|---|---|---|---|---|---|
| 27.4 GB | blocks | 18.1 min | 22.6 GB | not recorded | 3,099,848 |
| 27.4 GB | per-base | 12.7 min | 13.3 GB | 49.6 GB | 3,099,848 |
| 48.3 GB | blocks | 38.3 min | 33.6 GB | 96.7 GB | 3,099,848 |
| 42.1 GB | per-base | 20.4 min | 16.3 GB | 83.1 GB | 3,099,848 |
| 39.2 GB | per-base | 16.7 min | 13.6 GB | 78.5 GB | 3,099,848 |

The four come from three person folders (the 48.3 and 39.2 GB BAMs are two libraries in one folder,
which docs/LIVE_SAMPLES.md says are different people).

Each BAM had a few windows (18 on the first) holding bases deeper than 32,767, reported as wrapped
(the Int16 finding above).

Left: a cache spill of the windows keyed by `facet_id` (RM268's pattern), so a BAM is read once; a
quantized output from one depth pass; mosdepth's thresholds file; a person-level layer that pairs a
coverage facet with callsets.

## RM279 — A method-error atlas
Severity: medium · Status: open (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose GIAB stratifications as regions (E8): error and discordance rates by caller and region from multi-facet PGP-HMS people (docs/LIVE_SAMPLES.md). Needs RM277.

## RM280 — A consensus posterior
Severity: low · Status: open (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose after the atlas (E7): a per-site posterior across one person's facets as a derived facet, with RM279's rates as priors; never before them.

## RM281 — Chip manifests
Severity: medium · Status: in progress (2026-10-23: I/D codes resolved against Illumina's indel probes, a mixed-strand array's withheld rows read by their probe's allele pair, absence with a cause from manifest or learned probe sets; left: WGS-learned I/D orientation, F17, the Forward/Design-only Final Report) · Owner: lead

The person chose both sources (E9): public vendor manifests (Illumina GSA/OmniExpress are on disk) and probe sets learned from PGP people with an array and a WGS. Resolves I/D codes, strand per probe, and which sites a chip has (absence with a cause). Closes F17.

**Shipped 2026-10-23** (`just_vcf.probes`; rules in docs/SCHEMA.md Decision 2 and Decision 4). The
manifests are read at run time from `JUST_VCF_MANIFESTS_DIR`, each once into a probe table in the cache.

- **I/D codes** (`resolve_ref`): matched with a manifest indel probe by rsid; the site located in the
  reference by the probe's flanks (GSA-24v3-0 A1: 9,618 of its 10,118 indel probes found at exactly one
  place against GRCh37, on the strand `SourceStrand` gives for all but 4); spelled with REF at the
  export's own position, read from the reference. Measured on PGP-HMS GRCh37 exports (counts only):

  | Chip | Exports | I/D rows each | Resolved | Resolved homozygous rows read hom-ref |
  |---|---|---|---|---|
  | `23andme_v5` | 2 | 4,815–5,166 | 46.6–48.6% | 98.7% |
  | `ancestrydna_v2` | 1 | 7,631 | 43.0% | 98.4% |
  | `23andme_v3` | 2 | 1,032–1,038 | 17.1–17.6% | 76% |
  | `23andme_v4` | 2 | 4,840–4,866 | 3.6–3.7% | 76% |
  | `23andme_v1_v2` (build 37 re-export) | 1 | 121 | 6.6% | 6 of 6 |
  | `ancestrydna_v1` | 2 | 0 | | |
  | `living_dna_sirius` | 1 | 0 codes (83 concatenated indels, F17) | | |

  Most unresolved rows carry no rsid a GSA indel probe names (in two other v4 exports, 94% of about
  15,200 I/D rows each carry internal `i…` ids); about 2% have flanks not found, under 1% a site no shift brings to the export's position.
  Build-36 exports get none: the manifests map to GRCh37. **Open:** on 23andMe v3 and v4 the same sites
  read hom-alt in every export (43 of 179 v3 sites and 40 of 156 v4 sites in all of four exports each,
  against 21 of 2,313 on v5), with hardly a het. REF agrees with Ensembl's GRCh38 REF on long or short
  for all 154 biallelic v3 sites checked, so the reference read is right; whether those sites are ones
  where the reference carries the rare allele or where 23andMe's I/D mean the other way round needs the
  same person's WGS (left, below). **Checked on one v3 person** (2026-10-23, hu0E7AAF, the export's
  codes resolved with the gate lifted in memory, against the same person's CGI `var` on GRCh37): of
  171 resolved homozygous sites, the 44 the export reads hom-alt and the 127 it reads hom-ref all lie
  inside a CGI reference block with no variant within 2 bases. The join does find CGI indels: the
  export's only 2 resolved het sites both meet a het indel in the `var` (58,864 indel calls on the
  contigs read); a small control, but not an empty one. So on those probes 23andMe's I/D is the
  other way round from Illumina's, not a rare reference allele; v3 stays withheld (now with that
  evidence in `ID_CALIBRATION`), and a per-probe flip waits for more people. v4 is not checked yet:
  loading one CGI `var` from Krypton took about 25 minutes.
- **Strand per probe** (`resolve_ref`, mixed facets only): a withheld row is read by its probe's two
  alleles on the plus strand (`SNP` as written for `RefStrand` `+`, complemented for `-`; the GRCh37 base
  at `MapInfo` is in that pair for all but 2 of 20,000 sampled SNV probes of each manifest). Where the
  pair holds the row's reference base and is not `A/T` or `C/G`, the row reads forward or minus by its
  alleles alone. This is not the repair F13 rejected: F13 showed that no strand column predicts which
  probes FTDNA flipped; the pair decides each row whatever the vendor flipped by.
- **Absence with a cause** (`chip_absence`, `probe_set`, `learn_probe_sets`): `in_export`, `no_call`
  or `not_on_chip` per asked site (a written `--` is `no_call`, not `in_export`: the lead's call,
  2026-10-23; rejected: a fourth member for a written no-call, since both say the chip looked and saw
  nothing, and the row is there for a caller who needs to know it was written), from the chip's own manifest (`illumina:<manifest>`, configured) or a
  set learned from its exports. Learned from every loadable 23andMe and AncestryDNA export under PGP-HMS
  (sites carried by at least half of a chip's exports; every file counts, so a participant's copies of
  one upload count more than once):

  | Chip, assembly | Exports | Sites on the chip | Positions any export carries |
  |---|---|---|---|
  | `23andme_v1_v2`, GRCh36 | 81 | 576,342 | 1,175,160 |
  | `23andme_v1_v2`, GRCh37 | 10 | 574,476 | 594,753 |
  | `23andme_v3`, GRCh36 | 164 | 966,581 | 1,003,300 |
  | `23andme_v3`, GRCh37 | 107 | 960,586 | 1,001,149 |
  | `23andme_v4`, GRCh37 | 77 | 610,326 | 614,508 |
  | `23andme_v5`, GRCh37 | 28 | 635,966 | 654,063 |
  | `ancestrydna_v1`, GRCh37 | 25 | 701,478 | 944,173 |
  | `ancestrydna_v2`, GRCh37 | 16 | 797,008 | 828,658 |

  AncestryDNA's `V2.0` covers two contents (its exports have 664,429 rows; 797,008 sites are in half
  of them or more), and `V1.0` likewise (425,316 sites in every export): the stated version does not
  pin the probe set. A learned set is the cache's, never shipped: at 0.6–1.0 million positions per
  chip it is several MB, past the few hundred KB a packaged constant may take.

Decisions, each with what was rejected:

- The manifests' indel probes serve every array on their build, matched by rsid. Rejected: only the
  chip's own manifest, which resolves nothing on PGP-HMS (no export's chip has its manifest on disk).
  Risk kept: a vendor's probe for a multi-allelic rsid may test another allele than Illumina's; the
  hom-ref share above is the check, and v3/v4's is the open question.
- I/D codes resolve only on a chip whose calibration passed (`probes.ID_CALIBRATION`, the lead's
  call 2026-10-23): `23andme_v5` and `ancestrydna_v2`; the rest withheld, counted
  `id_orientation_unverified`. Rejected: resolving every chip (v3 and v4 would ship a quarter of their
  homozygous rows possibly backwards), and flipping v3/v4's codes (nothing yet says which way is
  wrong). An `illumina:` Final Report is withheld too, though its own manifest defines its `I`: no
  Final Report with I/D rows has been measured, and the catalog matches rsids across every configured
  manifest, not only the chip's own.
- A site is located by its flanks, not by `MapInfo`, which is not the anchor (the gap is one base
  before it for 5,594 GSA probes, at or after it for the rest). REF comes from the reference; the
  manifest gives the inserted bases and the locator only.
- REF starts at the export's own position, the site shifted along a repeat to get there, padded on the
  right where nothing lies between. Rejected: moving `pos` (it is what the source wrote) and the empty
  allele (VCF has none, and an empty REF reads as an insertion after `pos`).
- `gt_alleles` of a resolved row becomes the two sequences, `source_gt` keeping the codes. Rejected:
  leaving `I`/`D` there, which would break `gt` indexing `[ref] + alts`.
- Strand by pair only in mixed facets. Rejected: reading forward facets by it too, which would rewrite
  rows of every export whose rsids a manifest shares, on evidence the strand check found sound.
- A learned site needs half of the chip's exports. Rejected: the union, which took 578,794 positions
  from one of 81 `23andme_v1_v2` build-36 exports into the set.
- A platform family (`illumina_gsa`, `illumina_omniexpress`) has no probe set. Rejected: the family's
  manifest, which would put another version's probes on the chip.

Left: WGS-learned I/D orientation for the codes no manifest names (23andMe v4's internal ids) and the
v3/v4 hom-alt question, from PGP people with a 23andMe export and a CGI `var` (its reference blocks say
hom-ref outright); F17 (DOGFOODING: no public Sirius manifest, one export); the Final Report with only
Forward/Design alleles (its manifest, HumanOmni2.5-8v1_A, is not on disk); exact-manifest absence on
real data (no PGP-HMS export's chip has its own manifest on disk).

## RM282 — Structural-variant facets
Severity: low · Status: loading built 2026-10-23, awaiting merge (decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose load with ids withheld (E11): SV calls load as rows with type and `end_pos`, VRS ids withheld with a reason (VRS 2.x has no SV model; charter: never a non-conforming GA4GH value).

**Built:** the SV columns `sv_type`, `cipos`, `ciend`, `copy_number`, `alt_copy_numbers`, `<CNn>`
extents, copy-number segments kept out of `blocks` (F18), `sv_no_vrs_model`, and the `write_vcf`
round trip (docs/SCHEMA.md, "Structural variants"). All five PGP-HMS SV and CNV files load. **Left:**
comparing SV calls across facets (truvari, LANDSCAPE); no DELLY or Sniffles file exists to check
against (synthetic tests only); and the question whether VRS 2's `CopyNumberCount`/`CopyNumberChange`
and `Adjacency` should be minted for CNVs and breakends, which the reason's name presumes away.

## RM283 — An aggregate-only query helper
Severity: medium · Status: built on a builder branch 2026-10-23 (`just_vcf.query`), waiting for the lead's merge; decided by the person in the Part 2 interview, 2026-10-23 · Owner: lead

The person chose a helper beside `load` (E12): an API that answers counts, shares and presence questions and never returns rows, for agents and apps; `load` stays for code that may read rows.

## RM284 — Phasing as CisPhasedBlock
Severity: low · Status: built 2026-10-23 on a builder branch, awaiting the lead's merge · Owner: lead

The person chose columns with an export on demand (E5): phase sets stay columns; VRS CisPhasedBlock objects are minted only when asked.

**Built:** `cis_phased_blocks` (docs/SCHEMA.md, "CisPhasedBlock export"). The VCF loader reads PSL,
PSO and PSQ since 2026-10-23, and `write_vcf` writes them back under VCF 4.5, so PSL sets are tested on
a loaded facet against vrs-python. Left open: in a VCF 4.4+ file a slot-0 block with an ALT is
`phase_unknown` until polars-bio keeps the leading phase indicator (P4), unless PSL names that slot's set.
A withheld block names the first allele that withheld it (`withheld_chrom`, `withheld_pos`,
`withheld_allele`). On a whole imputed genome
an implicit set is a chromosome, so one symbolic ALT withholds the haplotype's block (a phased Genes
for Good upload: 45 of 46 withheld `symbolic`); whether a block per window or per contiguous run of
sequence alleles is wanted is the person's call.

## RM285 — The format registry as a spec
Severity: low · Status: built 2026-10-23 (docs/FORMAT_REGISTRY.md, generated by scripts/format_registry.py; moves to ROADMAP_HISTORY.md on merge) · Owner: lead

The person chose this repo (E14): the declared formats and the broken shapes (F1–F17) documented as a readable spec in docs/ with synthetic examples; no genome data.

## RM286 — T2T-CHM13
Severity: low · Status: built 2026-10-23 except the karyotype calibration (below) · Owner: lead

The person chose after GRCh37/38 (E10): CHM13 in the reference store, the assembly table, liftover chains and the sex-chromosome regions.

**Built 2026-10-23.** Most of the assembly side was already there: `ASSEMBLY_LENGTHS` and
`CHR1_LENGTHS` carried T2T-CHM13 v2.0's lengths, so a `chr1`… header was labelled and a store could
name a CHM13 collection. Added:

- T2T-CHM13v2.0's GenBank (`CP068277.2`…) and RefSeq (`NC_060925.1`…) accessions map to canonical
  `chrom` (`schema.CHM13_ACCESSIONS`, from the NCBI assembly report carried in `assets/`), so an
  accession-named header is labelled too. VCF build hints state T2T-CHM13 from `chm13`, `t2t`, `hs1`
  and the assembly accession. The `T2T-CHM13` label cannot vouch for MT: `chm13v2.0.fa` carries the
  CHM13 mitochondrion, the analysis set's `_maskedY_rCRS` variant rCRS, both 16,569 bases.
- `PAR_BOUNDS` and `X_DEGENERATE_Y` for T2T-CHM13 (docs/PROVENANCE.md). A T2T facet's sex chromosomes
  are grouped and measured; the estimate answers `uncalibrated_build`.
- A store with all three assemblies, built as a new directory from the public FASTAs:
  `/data/just-dna-cache/reference/store-chm13` (T2T-CHM13 `YsDEAoUZA7_uLIT_YwqS3SUWV03NtbOY`, GRCh38
  `oLfPx0NOBKKXMIngGeQ4YewtU4Ge_wKz`, GRCh37 `fOP7cqLAXIU31eAMrGsK0WN_Htp__AcH`; the last two equal
  `store`'s). UCSC's `hg38ToHs1` and `hs1ToHg38` chains are in `/data/just-dna-cache/liftover/`.
  `tests/test_chm13.py` runs against them with `JUST_VCF_TEST_REFERENCE_DIR=…/store-chm13` and
  `JUST_VCF_TEST_CHAIN_DIR=/data/just-dna-cache/liftover` (about 2.5 minutes), and skips otherwise.

Decisions, each with what was rejected:

- The CHM13 X-degenerate Y is the T2T consortium's `X-DEG` class less GIAB's XTR and PAR. Rejected:
  GIAB's nonPAR minus XTR as for GRCh37/38. GIAB's CHM13 set has no ampliconic file, and on a gapless Y
  that difference keeps the ampliconic, heterochromatic, centromeric and satellite classes.
- T2T-CHM13 is measured but not estimated (`uncalibrated_build`). Rejected: the GRCh37/38 bands
  applied unmeasured; its Y region set differs and no T2T genome was among the calibration files.
- Accessions are matched with their versions, verbatim. Rejected: version-insensitive matching, which
  would let a later accession version (another sequence) pass as this one.
- Opt-in tests skip on a store without a T2T collection, through the existing
  `JUST_VCF_TEST_REFERENCE_DIR`. Rejected: a second store variable.

**Dogfood, 2026-10-23: GIAB's HG002 Q100 v5.0q small-variant benchmark on CHM13v2.0** (public;
docs/PROVENANCE.md), with the library only. Load: 5,829,374 rows, labelled `T2T-CHM13` (inferred from
its `chr1`… lengths; the file states no build), `wgs`. `resolve_ref` against store-chm13: 5,829,374 of
5,829,374 match. Sex-chromosome evidence: `uncalibrated_build`; the file could not have been
estimated anyway, since dipcall writes non-PAR X and Y haploid (119,112 and 430 rows, all one allele)
and no DP. Lifted to GRCh38 through `hs1ToHg38`: 4,502,840 lifted (1,644,705 with REF and an ALT
swapped, 85,844 on minus-strand blocks), and not lifted: 917,871 `no_chain_block`, 167,650
`allele_not_sequence` (`*`), 142,505 `split_across_blocks`, 75,001 `target_contig_unmatched`, 23,507
`ref_mismatch_after_lift`.

The 75,001 are a finding: 34,178 of the chain's 856,832 blocks target UCSC's scaffold names
(`chr1_KI270706v1_random`, `chrUn_…`, and `_fix` patches), which the Ensembl GRCh38 primary assembly
names `KI270706.1` or does not hold (patches). The rows are counted, not lost silently, but they cannot
lift onto an Ensembl-named GRCh38 store. A UCSC↔GenBank scaffold alias table would fix the first kind.
It does since 2026-10-23 (F19, docs/PREVIOUS_ISSUES.md): 40,099 left, all on `_fix` patches.

**Left:**
- Calibrating the karyotype thresholds on T2T callsets of known karyotype with FORMAT DP (the HG002
  benchmark has none); then `ESTIMATE_BUILDS` takes T2T-CHM13.
- Making store-chm13 the default test store, or adding CHM13 to `store`, is the lead's call.

## RM1 — Final shape and plan
Severity: blocker · Status: charter drafted, not signed · Owner: user

The shape as the user described it on 2026-09-28 is written into the draft
`docs/CHARTER_DRAFT.md`, committed so work can build against it and edited directly until signed off. The frame is the product (Arrow/Polars, CompoSeq-style). There is a minimal
many-in/one-out shape and an extended shape of per-person facets (assemblies, callers, phasings,
coverage) with graded evidence, provenance, concordance and a graph. It is GA4GH first: VRS through
gtars, seqcol through refget. The user keeps it as a draft; planning of RM2–RM6 follows from it.
Surveys: docs/SIBLING_SURVEY.md, docs/LANDSCAPE.md.

## RM2 — Loading consumer formats
Severity: high · Status: in progress (on 2026-10-23, 32 of the 39 declared formats load: the earlier 25, plus the CASAVA SNP and indel tables, SOAPsnp, FTDNA Big Y, Color's variant set, GET-Evidence GFF and Nirvana JSON; `oxford_gen` waits for a decision (docs/LAYERS.md), Geno 2.0 and Final Reports without positions wait for chip manifests, `promethease_tab_genotype` has no positions, `cgi_summary` and `decodeme_info` are not facets) · Owner: lead

**Shipped 2026-10-23: the site tables** (`loaders/site_tables.py`): CASAVA SNP and indel tables,
SOAPsnp, Big Y, Color's variant set and GET-Evidence GFFs; then Nirvana JSON (`loaders/nirvana.py`),
32 of 39 formats now. What still does not load, and why, is
in docs/FORMATS.md ("Not loadable yet"); the `oxford_gen` chromosome question waits for the
person. Builder's choices, each with what was rejected:

- A one-letter call is one `gt_alleles` entry (CASAVA `Diff`, a plain SOAPsnp base, GET-Evidence
  `alleles G`). Rejected: doubling it to a homozygous pair, which asserts a diploid call on male X and
  Y, where the same files write one letter too; the ploidy is the karyotype work's.
- CASAVA `dp` is `Total Bases`, and `ad` the base counts of `[source_ref] + alts` (the four counts sum
  to `Bases Used`). Rejected: `dp = Bases Used`, which is the filtered depth GATK calls the AD sum.
- `Sample Name` is not `facets.sample`: it is free text, in real files a person's name, not a sample
  column. Rejected: carrying it as the identifying sample, which the spec says to drop.
- No build is inferred for CASAVA, SOAPsnp or Big Y, and `caller` stays null for all five: neither is
  stated. Rejected: CASAVA or SOAPsnp as the caller from the layout, which the file does not name.

**Shipped:** content-based detection. `just_vcf.sniff` and the declarative registry in
`just_vcf.formats` (docs/FORMATS.md). Over the PGP-HMS participants, 8 of about 1,740 files stay unknown,
all prose or empty.

**Next: per-format loaders on polars-bio**, emitting the frames in docs/SCHEMA.md (schema v0.0.1), driven by each `FormatSpec`. Start from the
contract gaps listed at the end of docs/FORMATS.md (versioned columns, section headers, 0-based CGI
conversion, half-calls, chromosome outside the row) and from polars-bio pending fixes P1/P2. Port from
just-prs `arrays.py`, `io_utils.py` and `vcf.py`; leave out the observed-allele-as-REF encoding (charter
invariant 3). Living DNA writes indels as run-together alleles (`CCA`); Color positions are 0-based.

The multi-facet participants in docs/PGP_HMS.md are the convergence tests: one person's facets, loaded
independently, must agree where they overlap.

**Shipped 2026-10-22: the `source` detector** (`just_vcf/source.py`), naming the caller or vendor from
VCF headers, because `Family` cannot say "array or imputed content inside a VCF" and that decides
whether absent means unknown. It fills `facets.vendor`, `caller`, `caller_version`, `phasing_method`
(Eagle), the new `ploidy_statements` (RM269 half 1) and, where the header states array or imputed
content, `callset_type`. Every row of the table below is a rule, plus `##bcftools_callVersion` and
`##DRAGENVersion`; `##workflow_type=Illumina_WGS_CLIA`, `##source_version` and `##content` (Isaac and
starling) were added from the PGP-HMS headers. Decisions:

- `callset_type` takes only what the header states (`array`, `imputed`); rejected: inferring it from
  the rows, which is RM4's. RM4's inference should not override a stated value.
- `vendor` is the party the header names as delivering the file; a caller's maker is not a vendor
  (DRAGEN and CASAVA lines do not make Illumina the vendor; `##workflow_type=Illumina_*` does). Rejected:
  a vendor per caller, which would call a Sequencing.com delivery an Illumina one.
- Post-processing tools (GATK SelectVariants, ApplyVQSR, VariantFiltration, bcftools view and concat)
  are not caller evidence, and GenotypeGVCFs is HaplotypeCaller's own pipeline. Rejected: GATK as the
  caller whenever a GATK line exists, which names GATK for a DeepVariant file run through SelectVariants.
- Two callers, or two versions of one, leave the field null with one warning naming each line
  (charter 8: no winner picked). A `vcf-merge` header is not a conflict; the merged file's facets
  carry one caller each (Decision 8).
- `ploidy_statements` keeps lines verbatim with no vocabulary, as RM269's critique recommends; an
  option's path value is written `<path>`.

**Shipped 2026-10-23: the method chain** (layer 2, charter 6), `source.method_chain`: `facets.method_steps`,
`reference`, `aligner`, `aligner_version`, and `phasing_method` from the chain. Measured over 392
PGP-HMS VCF headers and 52 local ones (the Krypton copy of `/data/genomes`, ancestry panels left out):
`reference` 95% and 38% (DeepVariant writes no `##reference`), at least one step 80% and 83%, a step
after the call 76% and 48%, `phasing_method` 34% and 38%, `aligner` 2% and 2% (DRAGEN only; no header
on disk copies a BAM's `@PG` lines). Decisions, each with what was rejected:

- Step order is header order, except dated steps (GATK, DRAGEN, bcftools), sorted by their stated
  date among the places they hold: every GATK pipeline header on disk lists ApplyRecalibration before
  HaplotypeCaller, because htsjdk sorts meta lines by text. Rejected: header order alone, which is
  wrong for every GATK header; and a sort of all steps by date, which has nowhere to put the undated.
- `reference` is `##reference=`, the coordinate system the rows are on, with the caller's reference
  only as a fallback. Rejected: the caller's reference first, which for a file lifted by a tool that
  records nothing (Picard) names the reference the rows are no longer on. When the two differ and no
  `lift` step is recorded, a warning says so (four local headers, the Picard lifts).
- A caller phased its own calls only where FORMAT `PS` is declared in the caller's own words
  (`_CALLER_PHASE_SETS`). Rejected: any `PS` declaration, which credits the caller with a phaser's work.
  HaplotypeCaller's PGT/PID alone is not a phase step: no `calls` column carries it.
- Phasers with no header on disk get a rule only where their source code shows the line they write
  (WhatsHap, SHAPEIT4, HapCUT2). Beagle, SHAPEIT5 and DRAGEN's `--vc-enable-phasing` without a `PS`
  declaration have no rule: no header on disk, no source read.
- `options` come from a list kept per tool, never the whole command: GATK writes forty defaults, and a
  command carries inputs, outputs and sample names. Rejected: every option with paths redacted, which
  still keeps `--vc-sample-name`.
- One call step per caller name, at its first line, collecting its command lines' options; every
  other line is its own step. Rejected: one step per caller line, which repeats a caller run per
  chromosome.
- The role vocabulary is open (`schema.STEP_ROLES`, a string column): `other` holds a conversion or a
  gvcftools pass, and a new role is additive.

**Open from `vcf_cgi` (2026-10-22):** a vcfBeta lists its no-called spans and half-calls but, in the
rows counted, never a reference-called record, which suggests an absent position was called
reference. The loader keeps `absence_means = blocks` (absent is not reference) until CGI's own
documentation or a comparison with the same genome's `var` file settles it; `weak_reference` or a
stronger value would be the change. **Settled by the comparison, the value not yet chosen
(2026-10-22):** on the smallest PGP-HMS assembly, vcfBeta covers every base its `var` does not call
reference (0 bases missed of 317.7M) and touches 904,739 of the 2.78G bases `var` calls reference,
the anchor bases of its indels and no-calls. So a position absent from a vcfBeta is one CGI called
reference, as explicitly as a `var` ref block. The test is `test_real_var_converges_with_its_vcf_beta`.
**Value chosen 2026-10-23 (the lead):** a new member, `AbsenceMeans.REFERENCE` (`reference`), just-vcf's
own (SCHEMA Decision 4); the native files keep `blocks`. Builder's choices, each with what was rejected:

- Keyed on the caller line (`##source=CGAPipeline_*`, `callset.REFERENCE_BY_OMISSION_CALLERS`), not on
  `format_id = vcf_cgi`: an upload that lost `##center` sniffs as `vcf` and is still CGI's VCF, as the
  GL-scale rule already reads it. Rejected: the format id alone, which would give the same file two
  readings of its silence.
- Only where RM4 types the rows `wgs`. A region cut from a CGI VCF (RM4's `panel`, or nothing it can
  name) keeps `blocks`, with a warning. Rejected: the caller alone, which would claim reference for the
  regions a cut took away; RM4's `content_filter` covers only cuts by class.
- `absence_is_reference` applies the same class and content-filter gates to `reference` as to
  `weak_reference`. Rejected: True for every class, which would read the absence of SVs, which
  CGAPipeline's small-variant VCF never writes, as their absence.
- Ancestry's absent-site rule (`ancestry/sites.py`) reads `reference` as reference outside any block,
  where it read the member as unknown.

**CGI native files, built 2026-10-22** (`loaders/cgi.py`). Decisions, each with what was rejected:

- Indels stay unpadded: a deleted allele is `""`, an insertion point sits at `pos = begin` with
  `end_pos` null (the Genos rule). Rejected: an anchor base, which needs the sequence (RM3); and
  `pos = begin + 1` with `end_pos = pos - 1` for an insertion, which keeps one interbase rule but puts
  an `end_pos` below `pos` in the frame.
- An uncalled insertion point (a zero-length no-call, 346,474 in one genome) is a `calls` row with
  every haplotype null. Rejected: dropping it with a count, which loses "looked and could not tell"
  (Decision 2); and a block on the anchor base, as vcfBeta writes it, which overlaps a `ref` block.
- `alts` sorted by string, slots as each file wrote them. Neither file's slot order is a phase (the
  two files of one assembly order the slots differently).
- `hapLink` is a phase set of its own kind, `cgi_haplink` (since 2026-10-23): an id per haplotype,
  one molecule, keyed per chromosome; a linked slot is `phased` with the id, an unlinked one null.
  CGI's documentation was not at hand; the reading rests on a PGP-HMS assembly, where no id leaves
  its chromosome, both haplotypes of a locus never share one, and CGI's own VCF writes the same
  phasing as 25,696 PS sets of two ids each (51,392 ids; 2,708,191 linked loci, 2,707,622 `|`
  records). `cis_phased_blocks` gives one block per id. Rejected: kind `psl`, which would name a VCF
  field the file does not have and key the set per file; kind `ps` with the slot as the haplotype,
  which is wrong because the slot is no phase.
- `gq`, `gl` null; the allele scores in `dropped_fields`. Rejected: `gq` as the lower
  allele score, the convention CGI's VCF writer uses, which is a derivation, not a value the file
  states.
- `varQuality` values kept verbatim (VQHIGH). Rejected: mapping VQHIGH to PASS without a CGI document
  that says so.
- `callset_type = wgs` only where `#TITLE` says whole-genome sequencing. Rejected: `wgs` from the
  format id alone.
- `masterVarBeta`'s `totalReadCount` fills `dp`. Since 2026-10-23 its read counts fill `ad` over
  `[source_ref] + alts`: REF from `referenceAlleleReadCount`, each ALT from the slot that carries it
  (measured: a reference slot always counts what `referenceAlleleReadCount` does, both slots of a hom
  ALT alike). Withheld, counted by reason, where a count is empty, an allele is partly called, or the
  counts contradict that. Rejected: summing a hom ALT's two slots, which counts the same reads twice
  (the sum passes `totalReadCount`); `ad` per slot, which is CGI's VCF AD and not VCF's.

The 1.x layout (formats 1.6 and 1.7: `totalScore`, no filter column) loads as it is, since the loader
takes columns by name; two PGP-HMS files have it (GET-Evidence collections), and a whole 1.6 file
loads as a `wgs` genome with no row dropped (tested, opt-in).

**Left for CGI:** memory (a whole `var` peaks at 8.9 GB, the haplotype sort);
GET-Evidence copies named `<assembly>.tsv.bz2`, not checked.

**Answered 2026-10-22: the HumanOmni2.5 Final Report that RM4's chip sweep could not load.** Its
columns are SNP Name, Sample ID, Allele1/Allele2 - Forward, Allele1/Allele2 - Design and GC Score:
no chromosome, no position, no plus-strand allele (2,379,855 rows of one sample; the header counts
38). No reading of it is not a guess. A position cannot come from the row, and an rsid lookup would
be another source standing in for the file's. Forward is dbSNP's orientation of the probe when the
manifest was built and Design the probe's own strand, so neither is the genomic plus strand. Loading
them as `gt_alleles` and leaving F13's strand check to catch the minus probes was rejected: a
consumer that never runs `resolve_ref` reads every minus-strand probe as the complement call, and the
strand check only withholds rows, it does not place them. The manifest the header names
(`HumanOmni2.5-8v1_A`) gives both the positions and the strand, which is RM4's manifest work. The
`TableLayoutError` now names that manifest, from the `Content` line.

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
| `##clair3_version=` | Clair3 |
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
Severity: high · Status: in progress (2026-10-22: the FASTA tier shipped as `resolve_ref` on a gtars refget store, for VCF and array rows, an array without a stated build inferred from its rows (F9); indel anchor bases and reference-supported positions are left) · Owner: lead

Port just-prs `reference_allele.py` (reference-panel tier, FASTA tier, FASTA build check). New work:
ALT derivation from genotype against REF, indel anchor bases, chr-prefixed FASTA support, emitting
reference-supported positions.

**Shipped 2026-10-22: `resolve_ref` (`just_vcf.reference`).** The rules are in docs/SCHEMA.md
(Decision 2, "How RM3 fills these slots"). Choices made, each with the alternative rejected:

- **Same facet, filled, not a derived facet.** Every column it writes is one the schema defines as
  "not read yet", and the result is fixed by the facet and its assembly's sequence. A derived facet
  (`relation = ref_checked`) would leave every loaded facet `not_run` for ever beside a checked twin
  with a second id. The consequence: a resolved and an unresolved frame of one `facet_id` can both
  exist in a session, so a consumer must not join a resolved `facets` with an unresolved `calls`;
  `ref_seqcol_digest` (null before) and `ref` tell them apart. The lead may still prefer the derived
  facet; nothing else depends on the choice.
- **The reference is a gtars `RefgetStore` on disk** (`JUST_VCF_REFERENCE_DIR`, default
  `reference/` under the cache dir), over pysam faidx: it is a core dependency, exact, and gives the
  seqcol digest the facet records. A FASTA is added once; a fingerprint alias (path, size, mtime)
  keeps it from being digested again.
- **Assembly matching per contig, by stated length and md5, else by label.** The reference's label
  comes from its contig 1 length (the just-prs fingerprint idea), or a label the caller asserts when
  adding a FASTA, which its lengths must not contradict. Rejected: matching by label alone (hg19 and
  GRCh37 share it and differ on MT), and refusing a facet over any one unmatched contig (an hg19 VCF
  would lose every autosome over chrM).
- **Array ALTs sorted by allele string**, not in genotype order, so two exports of one chip agree.
- **A store holding several collections** picks the one matching the most rows; a tie is refused
  unless `collection=` names one.

**Dogfooded 2026-10-22** against Ensembl GRCh37 and GRCh38 primary assemblies in one store (built
from the gzipped FASTAs in 103 s and 25 s, 2.2 GB; each assembly named by its own contig 1 length),
on PGP-HMS files read from the Krypton share, aggregates only:

- Six plain VCFs (73–142 MB). Five resolved: 34.1M rows read, all but 2 matching their source REF, and
  both of those against an IUPAC code in the reference. One stated lengths for 194 GRCh38 contigs in
  UCSC names; its `chr*_random` and `chrUn_*` contigs (44,594 rows, 1%) stayed unread, because
  Ensembl names them by accession (`KI270721.1`) and `canonical_chrom` does not map the two. One VCF
  with neither `##contig` lines nor a build hint was refused. A sixth failed to load (not RM3's).
- Four build-37 23andMe exports (1.0–1.4M rows): 97.5–99.4% of rows read; the rest are I/D codes,
  no-calls and MT. 14.7–19.4% of the rows read carry no REF allele, the range expected for
  homozygous non-reference calls on a chip, not for a wrong build or strand. Four build-36 exports
  were refused: no GRCh36 reference in the store.
- `resolve_ref` on the 1.4M-row export: 11.4 s (load 1.3 s), of which 1.3 s reading bases.

Order with `rejoin_split`: re-join first, then `resolve_ref`. The other way round, the re-joined rows
carry no `ref` while the facet keeps `ref_seqcol_digest`, and `resolve_ref` refuses it as already read.

Left: indel anchor bases (no loaded format lacks them: 23andMe I/D codes need the chip manifest, and
Living DNA's run-together `CCA` is split against the reference since 2026-10-22, F4, and gets REF on 4
of 83 rows; the rest need the site's alleles, F17. Checked 2026-10-23: CGI's native files and Genos do
write indels without the anchor base, as the empty allele `""`. A CGI deletion's REF is read and
`normalize` anchors it; a CGI or Genos insertion point's anchor base is not read, and neither is REF
for Genos's multi-base rows, whose span the file states, so `normalize` cannot reach them yet); emitting
reference-supported positions; an opt-in test against a real assembly. UCSC's scaffold names meet
Ensembl's accessions since 2026-10-23, through NCBI's assembly reports rather than the accession inside
the name (F19, docs/PREVIOUS_ISSUES.md). A VCF with neither
`##contig` lengths nor a build hint is refused: nothing can show which assembly it is on. An array
without a stated build has it inferred from its rows (F9, below).

**Shipped 2026-10-22: an array's build from its own rows (F9).** `resolve_ref` scores a facet that
states neither a label nor contig lengths against every collection it considers
(`just_vcf.assembly.infer_assembly_from_rows`, called from `reference._infer_from_rows`); the rule
is in docs/SCHEMA.md Decision 2. Choices, each with the alternative rejected:

- **The scores go in a warning line**, and the label into the existing `assembly_label` with basis
  `inferred`. Rejected: a `facets` column for the scores, which is a schema change for a value only a
  reader of the inference needs; and a new `AssemblyLabelBasis` member (`rows`), which retypes the
  `pl.Enum` dtype, when `header_contigs` already tells the two `inferred` sources apart.
- **Inside `resolve_ref`, not at load.** Scoring needs the reference's bases, and loading reads no
  reference. Rejected: a separate public `infer_build(loaded, store)` step before a consumer asks for
  one; the function stays internal, as RM4's inference does.
- **Margin against another label, not another collection.** hg19 and Ensembl GRCh37 score alike, and
  are one assembly. An unlabelled collection counts as a rival (it could be another build). With one
  collection, chance is the rival: the mean number of distinct called alleles over four.
- **Arrays only.** A VCF without lengths or a build hint could be scored by its source REF, and would
  then pass a `ref_check` chosen by that same comparison. It stays refused.
- **MT is never scored**: the GRCh37 and GRCh36 labels each cover two mitochondrial sequences (the
  Yoruba chrM of hg19, hg18 and NCBI36, and rCRS).
- **A sample, not every row**: at most 2,000 rows per chromosome, every k-th by position, so the
  same file always scores the same rows.

**GRCh36 in the store (2026-10-22).** A reference with contig 1 of 247,249,719 bases is labelled
`GRCh36`, and that label, like `GRCh37`'s, cannot vouch for MT (NCBI36 and hg18 carry the Yoruba
chrM, 16,571 bases; an array labelled build 36 does not say which MT it used). Where to get it: UCSC
`https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.fa.gz` (943 MB, md5
`5cf361b2f92d01691ee7d555ffca721f` per the `md5sum.txt` beside it, checked), then
`ReferenceStore.open(…).add_fasta("hg18.fa.gz")`: 176 s, 49 sequences (`chr*_random` and the
haplotypes included). Not Ensembl release 54's `Homo_sapiens.NCBI36.54.dna.chromosome.*.fa.gz`: its
Y is a sub-range (`chromosome:NCBI36:Y:2709521:57443437`, the PARs cut off), so every Y position in it
is offset by 2,709,520, and the release publishes no checksums.

**Measured 2026-09-30: gtars RefgetStore as the sequence source** (gtars 0.10, Ensembl GRCh38
primary-assembly FASTA, 3.15 GB, 194 sequences):

- Build: `RefgetStore.on_disk` + `add_sequence_collection_from_fasta`, 25 s. Encoded store 1.1 GB; the
  25 primary contigs are 1.16 GB of it (2-bit, or 3-bit where a contig holds `N`). chr1 is 93 MB and
  gzips only to 67.5 MB, so compressing a hosted store buys little.
- Local lookups: about 8 µs each, for 1 base or a 100-bp window. The OS cache was warm (the store had
  just been written). 200 lookups checked against the FASTA: 0 mismatches.
- Remote (`open_remote` against a local HTTP range server, so no real network latency): opening
  fetches 57 KB of metadata. A single `get_substring` is one HTTP range request of about 20 bytes,
  1.2 ms on localhost; over the internet each one costs a round trip. `get_substrings` with many ranges
  downloads the whole contig once (chr1: one 93 MB GET) and serves every later lookup locally (3 ms
  for 500). Both modes wrote the contig into the cache directory, even after `disable_persistence()`.
- No public hosted store was found (gtars docs use placeholder URLs; searched 2026-09-30). A remote
  tier means hosting the store's static files ourselves.
- just-dna-format has no sequence access except seqrepo REST, and mints SNV ids offline from a committed
  GRCh38 accession table: a substitution's VRS id does not encode REF. Their RM267 shows a variant
  table used as the REF/position source producing wrong events that a REF check passes, which rules out
  a reference panel as the REF value here.

## RM4 — Input type inference
Severity: medium · Status: in progress (2026-10-22: a VCF facet gets `callset_type` (`wgs`, `wes`, `panel`, `imputed`, or null) and the `absence_means` that follows, and `assembly_label` from its contig lengths; an array without a stated build gets one from its rows in `resolve_ref` (F9, RM3); a facet says which variant classes it carries and whether it was cut to one (`content_filter`); chip from the header or, for 23andMe, the rows (`facets.chip`), and for FTDNA, MyHeritage and unstated Final Reports by marker overlap with Illumina's manifests; left: gVCF `END` expansion) · Owner: lead

Port just-dna-lite `detect_callset_scope` (WGS vs targeted), gVCF detection, and just-prs chip
generation detection. New work: WES vs WGS vs panel, and expanding gVCF `END` blocks.

**Shipped 2026-10-22** (`just_vcf.callset`, `just_vcf.assembly`; rules in docs/SCHEMA.md Decision 5):
the WGS rule is just-dna-lite's, unchanged; exome and panel bands were added with unknown between
every pair of classes. Decisions taken, each with what was rejected:

- gVCF stays out of `callset_type`: a gVCF can be a genome or an exome, so its blocks set
  `absence_means = blocks` and its calls set the type. Rejected: a `gvcf` member, which would change
  the `pl.Enum` dtype and answer a question the blocks frame already answers.
- A region extract reads as `panel`. Rejected: a separate "extract" outcome, which the calls cannot
  tell from a panel; either way absence outside the regions means nothing was looked at.
- Sequencing classes need a read-evidence FORMAT key, and rows whose ALTs are all symbolic do not
  count as sites. Rejected: reading GP/DS as imputation (DRAGEN writes GP), and counting SV records
  (a 5-sample DELLY merge in PGP-HMS had 1.04M SV sites at 93.5% breadth and read as `wgs`).
- `sorted_name_length_pairs_digest` is its own column; `seqcol_digest` stays null without bases.
  Rejected: filling `seqcol_digest` from a names-and-lengths lookup (the public seqcol server holds
  two collections for each of UCSC hg19's and hg38's coordinate systems).
- The classifier and the inference stay internal (`just_vcf.callset`, `just_vcf.assembly`, not in
  `__all__`); the facet columns are the contract. Rejected: exporting them before a consumer asks.

**Shipped 2026-10-22: chip detection** (`just_vcf.chip`; `facets.chip`, `chip_basis`; rules in
docs/SCHEMA.md Decision 5). Stated where the header names a chip: AncestryDNA's `array version`
(V1.0, V2.0), a Final Report's `Content` manifest (kept as `illumina:<manifest>`), Living DNA's
`Genotype chip`. Inferred for 23andMe exports from two features of the rows, each read alone: the row
count and the internal-id profile (share of rows with an `i…` id, and the share of those numbered
below `i1000000`). Measured on 341 PGP-HMS 23andMe exports, by generation:

| Generation | Exports | Rows | Internal ids, share of rows | Ids below i1000000 | rsids in GSA-24v3 / OmniExpress-24 v1.2 |
|---|---|---|---|---|---|
| v1/v2 | 61 | 571,676–595,526 | 0.30–0.56% | none | 0.20 / 0.55–0.56 |
| v4 | 77 | 598,897–612,647 | 7.71–9.72% | few | 0.20–0.22 / 0.55–0.56 |
| v5 | 28 | 631,455–643,535 | 2.69–2.80% | nearly all | 0.92–0.94 / 0.26–0.27 |
| v3 | 166 (2 read unknown, below) | 960,520–1,003,773 | 1.02–1.37% | few | 0.20 / 0.70–0.72 |
| none (1.41M rows, v5-style ids; 1,764–2,459-row extracts) | 9 | | | | |

The overlap columns are an independent check, not part of the rule: only the cluster named v5 is on
the GSA, consistent with 23andMe's v5 being its GSA chip, and the v3 cluster carries the most OmniExpress content.
The gaps the bands sit in are thin in two places: rows 595,526 (v1/v2) against 598,897 (v4), split at
597,000; internal-id share 0.56% (v1/v2) against 1.02% (v3), split at 0.8%. Measured on public files
only; a v1/v2 or v4 export outside them reads unknown, never the neighbour.

Decisions, each with what was rejected:

- Two features must agree. Rejected: the row count alone (just-prs), which cannot tell v1/v2 from v4
  (just-prs put both in one band and named it v4), and its fallback that names any unknown count
  "GSA v3 (estimated)".
- v1 and v2 are one value, `23andme_v1_v2`. Rejected: splitting them on the highest internal id
  (`i4001526` or `i4990307` in build-36 exports), which reflects the download date as much as the chip.
- FTDNA and MyHeritage get no row-count rule. Every FTDNA export of 548,000 rows or more in PGP-HMS is
  either the 2010 Affymetrix layout (548,011 rows, `ftdna_space` and three comma exports) or on
  OmniExpress (98–99% of sampled rsids in its manifest, 680,544–720,449 rows); none is on a later chip,
  so nothing here says what one looks like. Rejected: naming them by row count. *Superseded the same
  day by marker overlap (below), which reads the probe set rather than its size.*
- AncestryDNA is taken as stated. Its V2.0 exports are not on the GSA, where just-prs's chip table puts
  them: 0.26–0.28 of their sampled rsids are in the GSA-24v3 manifest and 0.64–0.67 in OmniExpress-24's
  (V1.0: 0.24–0.26 and 0.64–0.99).
- The chip is a facet column, not a computed function like `vrs_withheld`: it is read once at load,
  from the header and the rows, and does not change as the facet is processed.
- An open vocabulary (`CHIPS` + `is_chip`), stored as `str`. Rejected: a `pl.Enum`, whose dtype would
  change with every vendor release.

Dogfooded through `load` on the 522 PGP-HMS array uploads that sniff as a loadable table (one Final
Report did not load: `TableLayoutError`, its manifest `HumanOmni2.5-8v1_A` therefore unread):

| Chip | Basis | Exports |
|---|---|---|
| `23andme_v3` | inferred | 164 23andMe, 1 rsid table |
| `23andme_v4` | inferred | 77 23andMe |
| `23andme_v1_v2` | inferred | 61 23andMe, 1 headerless rsid table |
| `23andme_v5` | inferred | 28 23andMe, 1 headerless rsid table |
| `ancestrydna_v1` | stated | 24 AncestryDNA |
| `ancestrydna_v2` | stated | 16 AncestryDNA |
| `illumina:GSA-24v2-0_A1` | stated | 1 Final Report |
| `living_dna_sirius` | stated | 1 Living DNA |
| unknown | | 12 23andMe-shaped: 4 of 1.41M rows and 5 extracts (no band), 2 v3-sized whose internal ids are a v3/v5 mix (the features disagree), 1 spreadsheet re-save of 70,000 rows; 103 FTDNA, 2 FTDNA space-separated, 4 deCODEme, 15 Genes for Good, 4 Genos, 2 rsid tables without 23andMe ids, and one each of Color, PLINK, SelfDecode and the rsid CSV (no rule) |

**Shipped 2026-10-22: chip by marker overlap** (`just_vcf.chip.overlap_chip`; SCHEMA.md Decision 5,
`chip`). For `ftdna`, `ftdna_space`, `rsid_result_csv` and a Final Report whose header names no
manifest, the export's rsids are compared with a sample of the rsids of Illumina's GSA-24v3-0 A1 and
OmniExpress-24 v1.2 A1 manifests (every rsid whose number is divisible by 61, applied alike to both
sides; `just_vcf/chip_markers.json`, 187 KB, derived by `scripts/chip_markers.py`, sources in
docs/PROVENANCE.md). Per manifest: coverage, the share of the export's sampled rsids in the manifest,
and recall, the share of the manifest's sampled rsids in the export. A chip family is named, basis
`inferred`, when the export has at least 1,000 sampled rsids, the best manifest's coverage is at least
0.90 and 0.30 above the other's, and its recall at least 0.40; otherwise the chip is unknown and a
warning gives both shares per manifest.

Measured through `load` on every loadable array export under the Krypton copy's
`pgp-hms/participants/`, duplicates as uploaded (660 facets; 606 with 61,000 rsids or more, where the
sampled shares are within 0.011 of the full sets'). That corpus is larger than the 341 23andMe exports
the generation bands above were measured on, so the 23andMe counts below are not the same measurement
as that table's. Coverage / recall, sampled:

| Exports | GSA-24v3 | OmniExpress-24 v1.2 |
|---|---|---|
| 57 FTDNA, 680,544–720,449 rows | 0.24–0.25 / 0.26–0.28 | 0.98–0.99 / 0.94–0.99 |
| 3 FTDNA + 2 `ftdna_space`, 548,011 rows (2010 Affymetrix) | 0.15 / 0.13 | 0.33 / 0.26 |
| 43 FTDNA, 15,789–18,091 rows (X chromosome only; 245–270 sampled rsids) | 0.22–0.52 / ≤0.01 | 0.22–0.99 / ≤0.02 |
| 1 rsid CSV (`rsid_result_csv`, vendor unstated), 631,369 rows | 0.84 / 0.84 | 0.36 / 0.32 |
| 1 Final Report, stated `GSA-24v2-0_A1` (half its probes named by position) | 0.97 / 0.52 | 0.41 / 0.19 |
| 28 23andMe v5 (a custom GSA) | 0.93 / 0.92 | 0.28 / 0.24 |
| 271 23andMe v3 | 0.20 / 0.31 | 0.71–0.72 / 0.96–0.98 |
| 25 AncestryDNA V1; 16 V2 | 0.25–0.26 / 0.28; 0.26–0.28 / 0.27–0.30 | 0.64–0.99 / 0.60–0.97; 0.64–0.67 / 0.59–0.63 |
| F3's appended export (`appended_export` facet) | 0.28 / 0.30 | 0.66 / 0.63 |

Outcome on the formats in scope, through `load`: 57 FTDNA `illumina_omniexpress` (inferred); 43 FTDNA
unknown, too few sampled rsids; 3 FTDNA, 2 `ftdna_space` and the rsid CSV unknown, no manifest clearly
holds them; the Final Report keeps its stated chip, and the Omni2.5 one still does not load. No
Final Report on disk names no manifest, so that branch is tested on synthetic rows only. No file on
disk states MyHeritage (the rsid CSV states no vendor).

Decisions, each with what was rejected:

- The value names the family (`illumina_omniexpress`, `illumina_gsa`), not the manifest compared with.
  Rejected: `illumina:InfiniumOmniExpress-24v1-2_A1`, the Final Report spelling: the GSA Final Report
  above states GSA v2 and reaches 0.97 against the v3 manifest, so the overlap does not tell a
  manifest's version from its neighbour's; naming the one compared with would be a statement the
  export never made.
- Coverage, margin and recall all have to hold. Rejected: coverage alone, which an extract of a chip
  passes (39 of the 43 FTDNA X exports have 0.96–0.99 of their rsids in OmniExpress's manifest, and
  hold 2% of it); and a
  recall floor of 0.8, which would refuse a report like the GSA one that names half its probes by
  position. The coverage floor sits between 0.72 (23andMe v3 against OmniExpress, a custom chip around
  it) and 0.98 (the FTDNA exports named); the rsid CSV at 0.84 and SelfDecode at 0.80 against GSA sit
  in that gap and read unknown.
- A sampled marker list ships in the package. Rejected: reading the manifests at run time from a
  configured path (every consumer would have to fetch 150 MB of Illumina files for one column); and
  the full lists, 1.3M rsids, above what a source tree should carry. 1 in 61 keeps about 10,000 rsids
  of each manifest; the sampled shares stay within 0.011 of the full ones.
- 23andMe stays on its own rule. Rejected: the overlap for every array format: 23andMe's v5 would read
  `illumina_gsa`, hiding the generation the vendor rule names, and AncestryDNA states its version.
  The F3 appended export is not read either (no vendor stated; its shares match AncestryDNA V2's and no
  rule names that).

**Left:** more manifests (other OmniExpress and GSA versions, an Affymetrix one that might place the
2010 FTDNA exports, whatever the rsid CSV's chip is), which could name more of the unknowns above; an FTDNA
X-chromosome export, which reads like an extract of its chip, could take the chip of the same
person's autosomal export, but nothing links the two files yet; GRCh36 in
`ASSEMBLY_LENGTHS`, so a build-36 VCF's `##contig` lengths name it (the reference side knows GRCh36's
contig 1 since F9; the header side reads those lengths as matching no assembly); contigs named by RefSeq accession (`NC_000001.11`), which the label
inference does not read yet; gVCF `END` expansion. Unmarked imputation output with read-evidence keys
would read as `wgs`.

**Answered 2026-10-22: an SNV-only or indel-only split of a genome.** The person: "it can [be
treated as WGS], but with a caveat that ref restoration is additionally gated (needs force)". It
stays `wgs`; `facets.variant_classes` and `content_filter` say what it carries, and
`absence_is_reference` holds its absences back (SCHEMA.md Decision 4). Decisions, each with what was
rejected:

- The assay class stays `wgs` and the cut is its own column. Rejected: a new `callset_type` member
  (`wgs_snv`), which mixes what was sequenced with what was kept and changes the `pl.Enum` dtype.
- The rows decide first, with a residue allowed. Rejected: zero tolerance, which `bcftools view -v
  snps` defeats (it keeps a multi-allelic site with an indel ALT); and the file name first, which a
  renamed or concatenated file defeats.
- `##source=DRAGEN_SNV` is not a filter statement. Rejected: reading it as `snv_only`, as first
  proposed; a PGP-HMS DRAGEN `hard-filtered` upload carries indels under that very line.
- `variant_classes` counts alleles per class rather than listing the classes, so a residue stays
  visible. Rejected: a bare list of names, where one stray indel and a million look the same.
- `sv` is a symbolic ALT or a breakend; a deletion written out base by base is an `indel` at any
  length. Rejected: a 50 bp cut-off, which would split one caller's indels across two classes.
- The helper answers at facet level and says `None` for `blocks`. Rejected: answering per position
  here, which belongs with the materialising function that does not exist yet.

Dogfooded 2026-10-22 on PGP-HMS (hu1C1368, hu6ABACE, huBE0518) and one local DeepVariant genome:
each `.snp.` half is `wgs` + `snv_only` (rows), `absence_is_reference` `None` for `snv` until forced
and `None` for `indel` either way; the unsplit DRAGEN and DeepVariant genomes have no filter and
answer `True` for `snv` and `indel`, `None` for `sv`; a `.snp.` half reached under a name saying
`.indel.` keeps `snv_only` from the rows, with one warning quoting the word.

Left from it: an indel-only genome half reads as `callset_type` null, so its absences are unknown
anyway. Four PGP-HMS indel halves of 0.32M–1.05M ALT alleles did; the two DRAGEN ones reach 91% of
the span but hold 953,506 and 958,046 sites, just under `MIN_WGS_SITES`. Whether an indel half
should read as `wgs` (and so be gated by its filter instead) is open. A genotype table's classes are
not counted.

## RM5 — Population inference
Severity: medium · Status: in progress (2026-10-22: PCA projection, neighbour vote, centroid proportions, consensus and Privé's 21-group proportions shipped as `just_vcf.ancestry`; panels in docs/ANCESTRY_PANELS.md) · Owner: lead

Port the just-prs `ancestry/` package (PCA projection, KNN, proportions, consensus) behind the
`[ancestry]` extra, without the `PRSCatalog` wrappers or the Hugging Face pull.

Left: no reference-PCA panel exists on this machine, so `infer_ancestry` has run only on synthetic
panels (a just-prs `1000g`/`hgdp_1kg` build fetched by hand would be the first real one); the AADR
Human Origins builder (just-prs `aadr.py`, `eigenstrat.py`) is not ported; there is no uncertainty for
the PCA path beyond the vote and its margin (a site bootstrap does not apply to OADP as it does to
Privé's sums); the consensus multiplies methods from one panel as if independent; a facet on
another assembly waits for RM6's liftover.

## RM6 — Normalization and VCF writing
Severity: high · Status: in progress (2026-10-22: the VCF writer shipped, `write_vcf`, with a tabix index and gVCF blocks, normalization, `normalize`: split, left-align, trim, dedup, sort, and liftover, `liftover`: SNVs, MNVs and indels through a chain, run GRCh38 to GRCh37 (a GRCh36 chain reads, untested on data); left: strand, blocks in a lift, chromosome-name normalization of input, quality filters) · Owner: lead

Reimplement chromosome naming and quality filters (just-dna-lite's versions are AGPL and not the user's,
so they are a behaviour reference, not a source) and port coordinate liftover from just-prs.
New work: left-align, split multi-allelic sites, strand resolution (A/T and C/G), dedup, sort, and a
writer that keeps FORMAT/GT (just-dna-lite's writer drops them).

When it splits multi-allelic sites, the writer must not use the convention RM270 found in pre-split
input, where a split row's `0` can mean another ALT: it writes `.` for the other ALT, as vt, Glow,
vcflib and `bcftools --multi-overlaps .` do (S2, from just-dna-lite). `rejoin_split` (RM270) re-joins
such input today.

**Shipped 2026-10-22: the writer** (`just_vcf.writing.write_vcf`; rules in its module docstring).
One facet's `calls` out as a single-sample VCF 4.3, BGZF by default, through the repo's `BgzfWriter`.
Round-tripped through `load` on a PGP-HMS upload (956,064 calls) and a local DeepVariant genome
(6,436,476 calls) with 0 key mismatches, and checked by htslib (pysam, opt-in test). Decisions, each
with what was rejected:

- Our own text writer, not polars-bio's `sink_vcf`. Rejected: `sink_vcf`, which drops an undeclared
  column silently (P13), wants polars-bio's genotype shape rather than `gt` + `phased`, and leaves float
  and GT formatting, so byte determinism, to another library.
- VCF 4.3. Rejected: 4.2, where GP is phred and `gp` holds probabilities; 4.4+, which our own loader
  reads with the first haplotype's phase withheld (P4).
- `gl` is written as GL. Rejected: PL by rounding `-10·gl`, since the frame does not remember whether
  the source wrote PL, and a GL source would lose precision. `gl`, `gp` and `ad` are written only where
  `source_alleles` is the written `[REF] + ALTs`; elsewhere they are left out and counted.
- FORMAT is listed per row, keys with a value only. Rejected: one FORMAT string per file with `.` for
  missing, which turns a GT the source left out into an uncalled `.` on reload.
- REF is `ref`, else `source_ref`. Unwritable rows (no REF, a REF mismatch, no ALTs, an ALT that is
  no VCF allele such as F12's `-`, a genotype not indexed or indexed past its ALTs) raise `UnwritableRowsError` with a count per reason, unless
  `resolved_only=True`. Rejected: writing them with a placeholder REF, or quietly leaving them out.
- Split-flagged rows (RM270) are written as held, tagged `INFO/JUST_VCF_SPLIT=<signal>`, with a
  warning and the count. Rejected: refusing until `rejoin_split` runs, which would block exporting a
  pre-split file as it came; and writing them untagged, which loses the marker signal (the marker INFO
  key is not in the frame).
- Blocks are counted in the result and the warnings, not written, unless `blocks=True` (the gVCF slice
  below): a gVCF block record needs the REF base at its start, and the `blocks` frame holds none.
- Contig names: `stored` (canonical), `chr` (a listed map, `MT` → `chrM`, nothing by pattern; other
  contigs unchanged), or `source` (the facet's `contig_map`, and its `header_contigs` for a declared
  contig with no rows; refused where the source spelled one contig two ways). `##contig` lines are sorted in the rows' canonical order, so a tool that checks
  records against the header's order (GATK) accepts the file. Rejected: the source header's order,
  which can put `chrM` first while the rows put it after `Y`.
- The header names `facet_id`, `source_sha256`, `format_id`, a derived facet's parent and relation, the
  assembly label and the package version. It leaves out `source_path` and the source's `##SAMPLE` and
  `##PEDIGREE` lines, which can carry a real name. The sample column is the facet's `sample`, else
  the `facet_id` for a source that names none, or the caller's `sample=`.
- Rows must arrive grouped by contig with positions ascending, as `load` emits them; otherwise
  `UnsortedRowsError`, and nothing is left at the target path. Rejected: sorting inside the writer,
  which holds the whole frame in memory on a genome.

**Shipped 2026-10-22: normalization** (`just_vcf.normalize.normalize`; rules in docs/SCHEMA.md,
"Normalizing"). A derived facet, `relation = normalized`, from a facet `resolve_ref` has read, with
bases read from the same collection. Checked against `bcftools norm -m-any --multi-overlaps . -f`
(htslib's, bundled in pysam) on a synthetic file of 64 records: the same rows, genotypes and AD, and
the same GL wherever we keep one. Dogfooded on six PGP-HMS VCF uploads (GRCh37, 176 thousand to
7.7 million rows, read from the Krypton share), chromosome 22 against bcftools: 235,911 of 236,425
rows identical in position, REF, ALT and genotype; all 514 that differ, in two files, are rows the
loader flagged as one record of a pre-split site, which `normalize` leaves for `rejoin_split`.
`normalize` took 1–34 s per genome. Splits per genome ranged from 22 to 92,403 sites, trimmed rows
up to 55,350; one upload holds 9,213 positions with two records of one allele (they
differ in `ids`), kept side by side. Found on the way: F11, F12. Decisions, each with what was
rejected:

- Our own Polars pass, not a bcftools subprocess: the frame is the product, bcftools would need the
  facet written out and read back, and it is not a core dependency. Splitting and rewriting are Polars
  expressions; only the left-align loop runs in Python, on `(pos, REF, ALT)` scalars of the candidate
  rows (building rows as Python dicts took 78 s on an 831-thousand-row genome, against 8 s now).
  Reference bases are fetched in one batch per contig, 100 bases left of each indel, and again 16 times
  wider for the rows a repeat carried further.
- The other ALT is written `.`, not `0` (RM270's finding). Rejected: bcftools' default `0`, which reads
  as a REF haplotype the person does not have.
- Split rows carry `split_signal = marker`, as vt's output does once loaded, so `rejoin_split` can
  join them. Rejected: a new `SplitSignal` value (a new Enum category changes the column's dtype, a
  schema bump), and no flag at all (the row would not say it is one part of a site). In a written file
  the tag is `write_vcf`'s `JUST_VCF_SPLIT`, and since 2026-10-22 the loader reads it back as a marker
  (any value), so a written normalized facet keeps its flag. Rejected for reading it back: restoring
  the signal the value names (`slot_pattern`, `all_zero_beside_alt`), since the loader finds those
  patterns again from the rows and the file now carries a marker key. vt's `OLD_MULTIALLELIC` and
  bcftools' `--old-rec-tag` were rejected for writing because their values are the original record,
  which the frame does not keep per row.
- `gl` on a split row only where no GT slot holds another ALT. Rejected: bcftools' slice on every row,
  which on `1/2` leaves a `1/.` row whose likelihoods point at `0/0`. `gp` withheld on split rows.
  Rejected: bcftools' unrescaled slice (it does not sum to 1, and the loader drops such GP, F10) and a
  rescaled one (other values than the caller's).
- On a moved or trimmed row, `source_ref`, `ref_check`, `source_gt` and the VRS ids are null and `ref`
  carries the new REF. Rejected: `ref_check = match` on the new bases, which would compare the sequence
  with itself; and keeping the old `source_ref`, which no longer describes the row's span.
- Sites with a symbolic allele (`<DEL>`, breakends) are left as written. Since 2026-10-22 `<*>` beside
  concrete ALTs (every variant record of a DeepVariant gVCF) is dropped from the split rows: the
  concrete ALTs are split, the vectors sliced over the parent's alleles where exact, and a site left
  with one ALT is an ordinary biallelic row with no split flag. A site whose genotype calls `<*>` stays
  as written, counted. Dogfooded on a local DeepVariant gVCF: every one of its 6,079,801
  variant records carried `<*>`, none remains after `normalize` (173 s), 82,113 sites with two
  concrete ALTs split into 164,466 flagged rows, no genotype called `<*>`, and htslib reads all
  6,162,154 rows of the export. Rejected: bcftools' own `0/.` row for `<*>`, which our loader would read back as
  a reference block; dropping `<*>` where the genotype calls it, which loses a haplotype's stated
  non-reference allele; and flagging every row of a `<*>` site, which would make each DeepVariant
  variant record look pre-split, so a second `normalize` would hold it and `rejoin_split` would try to
  join a lone row.
- Only `ref_check = match` rows are aligned or trimmed. Rejected: trimming unchecked rows, which needs
  no bases but rewrites a REF nobody has shown to be right.
- A deletion written as ALT `-` (F12) is anchored by `normalize`, on the base before, once REF has
  matched; the loader keeps it verbatim and counts it, and `write_vcf` refuses it as
  `alt_not_an_allele`. Rejected: anchoring in `resolve_ref`, which fills slots of the same facet and
  touches nothing the source wrote, while anchoring moves `pos` and rewrites REF and ALT; rewriting
  it in the loader, which has no reference base to anchor on; and the loader refusing such rows,
  which would drop calls whose REF matches the reference everywhere it was checked.
- Pre-split input (RM270 flags) passes through. Rejected: aligning those rows, which would move one
  row of a site away from the other and keep `rejoin_split` from finding it.
- Exact duplicates are dropped; rows stating one allele with other values are all kept (charter 8).
  Rejected: bcftools `-d`, which keeps the first. Since 2026-10-22 a duplicate is judged on what the
  row states: the columns a move clears (`source_ref`, `ref_check`, `source_gt`, the VRS ids) are left
  out of the comparison, so a moved indel and its unmoved twin are one row, the unmoved one kept.
  Every other column still counts, `ids`, `filters` and `qual` included, so the upload whose 9,213
  positions hold two records of one allele that differ in `ids` keeps both. Rejected: keying on
  `(chrom, pos, ref, alts, gt, phased)` and the vectors alone, which would drop a row whose filters,
  quality or ids differ and so lose what the source said.

**Shipped 2026-10-22: liftover** (`just_vcf.liftover.liftover`, `read_chain`; rules in docs/SCHEMA.md,
"Lifting"). A derived facet, `relation = liftover_to`, from a facet `resolve_ref` has read, through a
UCSC chain file, with REF re-read from the target collection of the same store. It returns `Lifted`:
the facet, and a report of every parent row with its target position, strand and swap, or the reason
it stayed behind. Decisions, each with what was rejected:

- Our own chain reader and block lookup in Polars (the axis cut at every block edge, two as-of joins
  per row), in the core, no extra. Rejected: pyliftover per position, as just-prs does: one Python call
  per row, and it does not say which block it used, so a span crossing a gap is only seen by lifting
  both ends. pyliftover stays the test oracle.
- The chain's own contig sizes are the assembly check (charter 3), both sides, against the store's
  collections; a primary contig of another size refuses the lift, any other contig leaves its rows
  unlifted. Rejected: trusting the chain's file name or the facet's label.
- A span in more than one chain's block is unlifted (`multiple_targets`), not given the best-scoring
  chain's target. Rejected: pyliftover's and just-prs's first hit, which picks a target the data does
  not choose (charter 8).
- A target REF equal to an ALT swaps it with the old REF in place and re-indexes `gt`; one equal to no
  allele leaves the row unlifted (`ref_mismatch_after_lift`). Rejected: bcftools' SWAP=-1, which adds
  the target REF as a new allele and shifts the genotype: a REF matching none of the alleles is also
  how a wrong mapping looks, and nothing in the row tells the two apart.
- A lifted row's `source_ref`, `ref_check` and `source_gt` are null, as on a row `normalize` moved.
  Rejected: writing the parent's REF, carried across, as `source_ref` with `ref_check = match`; it
  would let `normalize` re-align lifted indels, but no file wrote that REF, which the column means.
  Since 2026-10-22 `normalize` re-aligns lifted rows anyway: a `liftover_to` facet read every row's REF
  from its target collection, which is what `normalize` needs (SCHEMA "Normalizing", step 2). The gate
  is the facet's relation (`REREADS_REF`). Rejected: a row-level gate (`ref` set and `ref_check` not
  `mismatch`), which would also align array rows and a normalized facet's moved rows; equivalent on
  lifted facets, but wider than the decision, and it would quietly start aligning any future row
  kind whose `ref_check` is null for another reason. Since 2026-10-23 (the lead's decision) array rows
  are aligned too, by a gate of the same kind: a facet read from a file (`derived_from` null), where
  `resolve_ref` is the only writer of `ref`, aligns its rows with `ref` set and no `source_ref` (a split
  Living DNA genotype, F4). Rejected again: the row-level gate, which would also take in a derived
  facet's rows that another step wrote `ref` on; and `ref_check = match` alone, which left every array
  indel as written.
- Indels on minus-strand blocks were unlifted (`indel_minus_strand`) in this slice; since 2026-10-22
  they are re-anchored, below.
- The derived facet's `absence_means` is `unknown`. Rejected: the parent's `weak_reference`, which
  would read every dropped row and every uncovered region as reference.
- Its label is the target collection's, basis `inferred` (from that collection's lengths, which the
  chain's target side matched). Rejected: `stated`, since no header stated it; and a new basis value,
  a dtype change for one case.
- The chain is recorded on the facet in a new column, `derivation_inputs` (schema 0.0.8), by file name
  and sha256. Rejected: only inside the `facet_id` digest, which says a chain was used but not which;
  and the chain's path, which is a label and can name a person.
- A row lifted to another contig than its source's loses its phase sets (a PS set is keyed per
  contig); VRS ids are not carried (re-mint with `mint_vrs_ids`, tested on lifted rows).

Compared on anton (docs/LIVE_SAMPLES.md: one GRCh38 DeepVariant callset, 6,079,744 rows, lifted to
GRCh37 by CrossMap 0.6.3 in two settings and by Picard with and without swap recovery; counts only).
With the Ensembl chain CrossMap used, `liftover` lifted 5,578,476 rows in 281 s (plus 175 s for
`resolve_ref`), swapped 13,098 and lifted 20,424 on minus-strand blocks. Every position it lifted is
in CrossMap's first output; 5,561,478 rows have the same REF, ALT and genotype there. Against Picard
with swap recovery, 5,578,457 of its rows are the same in position, REF, ALT and genotype, and 19 are
ours only. The tools lifted 28,000 to 30,000 positions this did not; the slice withholds 108,776 rows
under chains with different targets, 2,865 minus-strand indels and 342 REFs matching no allele (which
rows the tools' extra positions came from cannot be told: their output keeps no source position).
Unlifted: 271,958 in no block, 90,071 on contigs the chain does not name (unplaced GRCh38 contigs),
26,244 whose target contig the GRCh37 collection lacks, 1,012 across a block edge. With UCSC's
hg38ToHg19 chain it lifts 5,724,577 (UCSC's chain is netted, one target per base, so nothing is
ambiguous; its chrM is hg19's, so MT is `target_contig_unmatched`). Found on the way: F14 (resolved).
No bcftools on this machine, so the opt-in bcftools comparison was not run.

**Shipped 2026-10-22: a tabix index beside the output** (`just_vcf.tabix`, `write_vcf(index=…)`). The
writer hands each line's contig, interval and byte offsets to an indexer as it writes them, and
`BgzfWriter` records where each compressed block starts, so the `.tbi` costs no second pass. Checked
by htslib (pysam, opt-in): region queries through our index return the rows the frame says overlap
each window (305 windows, every bin level, symbolic deletions found by their END), and htslib's own
index of the same file has the same linear index and the same per-contig counts. Decisions, each with
what was rejected:

- Our own indexer, in the core. Rejected: `pysam.tabix_index` after the write, which needs the
  `[reference]` extra (pysam, not on Windows) and reads the whole file again.
- An interval runs to INFO `END` where the line writes one, else over REF, which is how htslib's VCF
  preset reads the line back. Rejected: always the REF span, which would leave a gVCF block or a
  symbolic deletion out of a query over its middle.
- A window no record overlaps takes the next window's offset, and the end of a BGZF block is written as
  the next block's start, both as htslib writes them (found by comparing with its index). Rejected:
  the previous window's offset, valid but looser and not byte-equal to htslib's.
- Not done: htslib's folding of small bins into their parents (`compress_binning`), which only makes
  the index smaller; and CSI, for contigs past 2^29 bases, which the writer refuses to index (pass
  `index=False`).
- On by default for BGZF output. Rejected: off by default, which leaves the default export unqueryable
  by region for a cost of one small file.

**Shipped 2026-10-22: gVCF output** (`write_vcf(…, blocks=True, reference=…)`; rules in docs/SCHEMA.md,
Decision 4, "Writing blocks back out"). The `blocks` frame goes out as reference blocks interleaved
with the calls, REF read from the collection `resolve_ref` read, and `load` reads the file back to
equal `blocks` and `calls` frames (synthetic gVCF with diploid, haploid and no-call blocks; htslib
parses it and its index finds a block from inside its span). Decisions, each with what was rejected:

- REF from the reference store, never from the frame or the source: the loader keeps no block REF.
  A facet whose REF was not read is refused, and a block the collection cannot vouch for is counted
  `block_ref_unread` (refused, or skipped under `resolved_only`). Rejected: `N` as REF, which htslib
  accepts but which states a base nobody read.
- ALT `<*>` by default, `<NON_REF>` on request (`BlockAllele`). Rejected: `<NON_REF>` only, GATK's
  spelling, where VCF 4.3 defines `<*>`; `load` reads both alike.
- GT `0/0` for a `ref` block and `./.` for `no_call`, at the ploidy the block's GL implies, else the
  contig's calls'. The frame keeps no block GT and `load` reads only whether one called anything.
  Rejected: always diploid, which writes a haploid X block of a male sample as `0/0` beside its
  haploid calls.
- FORMAT GT, GQ, MIN_DP and GL, the values the frame holds. Rejected: PL by rounding, for the reason
  calls write GL; and DP, which the frame does not hold for a block (it keeps `MIN_DP`, Decision 4).
- Blocks and calls merged one contig at a time, by position, calls first at a shared position, after
  both frames' order is checked as a stream. Rejected: a global sort, which the writer already refuses
  to do, and `merge_sorted` over the whole genome, which needs one key across contigs and a REF
  lookup inside the lazy plan.

Dogfooded on a local DeepVariant gVCF (counts only): 6,436,476 calls and 79,329,981 blocks written in
461 s (947 MB, a 3.2 MB index), none refused. Reloaded, the `blocks` frame is equal to the input's,
and every `calls` column is equal row by row except `ref` and `ref_check`, which a reload leaves
unread until `resolve_ref` runs again. htslib's whole-contig fetch through the index returns the
frames' count of records on all 182 contigs.

**Shipped 2026-10-22: minus-strand indels in a lift** (rules in docs/SCHEMA.md, "Lifting", step 3).
An indel on a minus-strand block whose alleles share their first base (VCF's padding base) is
reverse-complemented, that base dropped from every allele (`alts`, `gt_alleles`, `source_alleles`)
and the target's base before the lifted span put in front; REF is read over the span moved one base
left and compared as for any row, so a swap still works. Tested on a synthetic chain against
hand-computed alleles (an insertion, a two-base deletion, an indel at a block's edge) and against
pyliftover's image of the source base after the span. `LIFT_VERSION` 2.

On anton (the same callset and Ensembl chain as above, counts only): 2,894 indels land on minus-strand
blocks; 2,855 lift (2,844 re-anchored, 11 with no shared first base complemented in place), 10 now
fail as `ref_mismatch_after_lift`, 29 stay `target_contig_unmatched`, and none is `indel_minus_strand`
(before: 2,865 of them). Against the tools' output at the 2,959 rows at those lifted positions:
Picard with swap recovery has 794 with the same position, REF and ALT, and 2,151 more with an indel
of the same lengths up to 50 bases to the left, which is what its left-alignment does to a shift
that is rightmost on the target; CrossMap (variant A) has 105 identical, so it does not write these
indels this way, and which representation it writes was not worked out. The lift takes 93 s. Not
checked row by row: that each left-shifted Picard row is the same allele as ours.

Decisions, each with what was rejected:

- The new padding base is the target's base before the lifted span, wherever it lies, past the
  block's edge included. Rejected: requiring it inside the block, which would leave behind indels the
  block places exactly; the base is read from the target, so the record states the target's own
  sequence there.
- An indel whose alleles share no first base (`A>T,AG`, `AT>G`) lifts complemented in place, as an MNV
  does: no allele is empty, so it needs no padding base. Rejected: leaving it unlifted, which drops a
  representable row; and padding it on the left anyway, which changes a span no allele needs changed.
- An indel that would land at a target contig's first base stays behind as `indel_minus_strand` (the
  reason's meaning narrows; the name is kept). Rejected: anchoring on the base after the span, as VCF
  allows at position 1, which no chain in hand needs and no test covers.
- Not left-aligned on the target: the shift that was leftmost on the source is rightmost there.
  Rejected: aligning inside `liftover`, which is `normalize`'s job; `normalize` does not realign
  lifted rows today (`ref_check` is null), which stays open below.

**Left in RM6:** SVLEN for symbolic
SVs (the writer writes END, which the loader does not read for symbolic SVs yet, so their `end_pos`
does not survive a reload); PSL phase sets (dropped and counted; 4.3 has no PSL, and no loader reads PSL yet, so no frame holds a `psl` set to write: reading it comes first); strand resolution;
in a liftover, left-aligning lifted indels (`normalize` skips rows with no `ref_check`), blocks,
symbolic SVs, anchoring across a small gap of one chain, and T2T (UCSC hs1 chains).


## RM271 — VRS allele ids
Severity: high · Status: in progress (2026-10-22: `mint_vrs_ids` and `vrs_withheld` shipped for VCF and array rows, equal to vrs-python on 19,000 real alleles; left: the gtars RLE gap filed upstream, an opt-in test on a real assembly) · Owner: lead

Charter invariant 4: a call's identity is the VRS Allele id of its fully justified allele on a
refget-addressed sequence, and that id is the join key between facets on the same assembly. Today
`calls.vrs_allele_ids` is filled nowhere, though gtars (a core dependency) mints VRS ids and normalizes
alleles. It became possible with RM3: `resolve_ref` reads REF from a refget store and records the
collection on `facets.ref_seqcol_digest`, so each allele can be placed on a sequence with a refget
accession. An id that cannot be minted conformantly (an allele that cannot be justified, a contig with
no refget accession, an unresolved REF, an array I/D code) is withheld with a counted reason in
`vrs_errors`, never approximated. vrs-python is a test oracle only (CLAUDE.md §5).

**Shipped 2026-10-22: `mint_vrs_ids(loaded, reference)` and `vrs_withheld(calls)` (`just_vcf.vrs`).**
The rules are in docs/SCHEMA.md ("When an id is minted", "Why an id is missing"). Choices made, each
with the alternative rejected:

- **A step after `resolve_ref`, not inside it.** It refuses a facet whose REF was never read, and
  reads the collection `ref_seqcol_digest` names from the same store. Rejected: minting inside
  `resolve_ref`, which would make every REF read pay for ids, and leave no way to re-mint when the
  rule changes (`MINT_VERSION`) without reading REF again.
- **Every allele waits for `ref`, substitutions too.** Rejected: SCHEMA's earlier no-bases path for
  substitutions from a table of contig accessions per assembly. No table exists, and a contig
  `resolve_ref` could not match has no accession shown to be the facet's.
- **`[ref] + alts`, the REF slot minted** as VRS's reference allele (a `ReferenceLengthExpression`
  over REF's span), as SCHEMA and gnomAD have it. Rejected: one id per ALT only, which the brief
  named; it would retype a signed-off column.
- **The state choice and canonical JSON are ours; rolling and sha512t24u are gtars'.** gtars 0.10
  rolls an indel correctly but always returns a literal sequence, and its Python API digests only
  `LiteralSequenceExpression`; vrs-python gives a `ReferenceLengthExpression` for every deletion, every
  repeat insertion and the reference allele (VRS 2.x), so gtars' own ids differ from vrs-python's
  there. Rejected: withholding those alleles until gtars has it (most indels and every REF slot), and
  vrs-python in core (CLAUDE.md §5). Owed upstream: a VRS 2.x state in gtars' normalization and an RLE
  digest in its Python API; then `_justify` and the JSON shrink to calls. Not filed.
- **Rows walked as runs of one contig** (`rle` over the spilled input), reading only the columns
  minting needs; the ids go to their own files and sit beside the spilled rows by position. Rejected:
  a filtered scan per contig with a re-sort (resolve_ref's shape), one scan of the file per contig (170
  on the GRCh38 genome).
- **The oracle is `scripts/vrs_oracle.py`**: a vrs-python data proxy over the same gtars store, so
  both sides read the same sequence; vrs-python is a dev dependency.

**Dogfooded 2026-10-22** on three PGP-HMS VCFs, aggregates only (`scripts/vrs_conformance.py`), each
sample drawn with a fixed seed from the minted alleles and checked against vrs-python on the same
store:

- A GRCh38 genome, 4,651,432 rows: 9,256,478 of 9,345,929 alleles minted (363,409 deletions, 412,401
  insertions); 89,451 withheld as `ref_unchecked`, all on contigs `resolve_ref` left unread (UCSC
  `_random`/`chrUn_` names, RM3). 8,000 checked, 0 mismatch.
- A GRCh37 indel-only genome split, 956,064 rows, 1,039,271 ALTs: every allele minted; 9,000 checked
  (deletions, insertions, REF), 0 mismatch.
- A GRCh37 SNV-only genome split, 3,461,640 rows: all but 4 alleles minted (the 4 on rows whose source
  REF the reference contradicts); 2,000 checked, 0 mismatch.
- Time on the 4.65M-row genome: 59 s wall and 103–112 s CPU with the cache on a fast local disk; 112–208
  s wall with the cache under `data/interim/` on the NTFS disk, which `dd` measured at 2.5 MB/s during
  the runs (other sessions were writing to it; the workspace note says about 350 MB/s). Peak RSS of the
  minting process 1.5–2.0 GB; one contig's alleles and, for its indels, its sequence are in memory.

## RM269 — Karyotype finding and the caller's ploidy model
Severity: high · Status: in progress (the person chose the full label, desk #61, 2026-10-22; the evidence layer, X copies and an XY/XXY estimate per WGS facet built 2026-10-22, checked on 27 PGP-HMS VCFs; XX and XO from a coverage facet built 2026-10-23 (RM278); XXX, XYY and arrays wait) · Owner: lead

Motivating report: S1 (just-prs, 2026-10-11). Callers, reference panels and scoring tools each pick a
convention for male X and Y and rarely record it: a caller with a diploid-everywhere model (DeepVariant
1.1.0) writes a man's non-PAR X as `1/1` and `0/1`, a sex-aware one (DRAGEN) writes or filters it
differently, and arrays split by vendor (23andMe one letter; AncestryDNA the same letter twice, per the
`ancestrydna` format notes; that format has no loader yet, RM2). A consumer reading `gt_alleles`
cannot tell one copy from two.

Reproduced 2026-10-18 against the code, on a synthetic VCF: `##DeepVariant_version=1.1.0` in the
header still leaves `facets.caller` and `caller_version` null, and a male-X `1` loads as one allele
beside `1/1` as two, exactly as written (SCHEMA Decision 2). Nothing in the frames can say the karyotype;
SCHEMA.md's XYY stress case already lists that as open.

Two halves:

1. **Caller, version and ploidy model on the facet.** `caller` and `caller_version` exist in `facets`
   and are never filled; RM2's planned `source` detector fills them. This item adds what RM2's table
   lacks: `##bcftools_callVersion`, and the caller's ploidy model as a new nullable facet field (sex-aware,
   diploid everywhere, or an explicit ploidy), read from the header where it says so and null where it
   does not. Never guessed from the caller's name alone. Checked 2026-10-21 on the 32 local private
   VCFs (header keys and option names only): DRAGEN writes `##referenceSexKaryotype=` (a two-letter
   value, the karyotype it called with) beside `##DRAGENCommandLine`; GATK writes the ploidy in its
   recorded command lines, `sample_ploidy=2` in all 7 GATK files (the 6 GenotypeGVCFs ones also
   `--sample-ploidy 2`), which says diploid everywhere unless a per-contig ploidy is also given; the 15
   DeepVariant files carry only `##DeepVariant_version`, so their model is unknown from the header.
   The 2 CrossMap/Picard liftovers and 7 files with a bare `##source` say nothing. PGP-HMS files
   (freeBayes, VarScan2) are not on disk now and are unchecked.
2. **A person-level karyotype finding with its own provenance.** Evidence without depth: the non-PAR X
   heterozygous fraction and the presence of PASS Y calls, which separate XX, XY, XXY and X0 on a WGS
   facet. XYY needs depth, so it waits for a coverage facet. Weak or contradictory evidence gives
   unknown, never the nearest class.
   Counted 2026-10-21 on the local single-sample VCFs with PASS calls (16 files, several of them
   liftovers of one callset; aggregates only, and kept here without naming a file, since a file's
   class would name a person's sex). The non-PAR X heterozygous share of PASS ALT calls falls in two
   clusters, 0.009–0.115 (11 files) and 0.580–0.626 (5 files), with nothing between: that half of the
   evidence works. **PASS Y calls do not**: every DeepVariant and GATK file in both clusters has
   2,666–10,759 PASS Y ALT calls outside the PARs, so their presence cannot tell XX from XY or XXY.
   Y evidence needs a region rule (the male-specific X-degenerate region, without the X-transposed
   and ampliconic parts; which tools restrict it so is for the survey to check) or depth, and until it has one the
   finding cannot separate XXY from XX or X0 from XY. The build is needed for the PAR bounds, and
   `load` fills `assembly_label` only where the header states it (none of these 16 did); the count
   used the header's chrX length. A sex-aware caller (DRAGEN, per its documentation) writes one
   allele on non-PAR X and Y for an XY sample, so the X share must count a haploid call as
   hemizygous evidence, not as missing. The vocabulary is GA4GH Phenopackets' `KaryotypicSex`, used
   verbatim (charter invariant 2). PAR bounds are per assembly and are domain constants.

Critique (2026-10-22, an adversarial read of this entry; tool behaviour marked unverified is the
critic's recollection, not checked here):

- Without a working Y signal no class can be named honestly: a high X share fits XX, XXY or XXX, a low
  one XY, X0 or XYY, so "never the nearest class" leaves every answer unknown. A Y region rule (or
  depth) is a precondition for any class, not a refinement.
- A haploid GT is the caller's (or the array vendor's) sex decision, not evidence: counting it fills a
  value from the source that checks it. Keep what the caller asserted as its own line, count het/hom
  only on calls written diploid, and read a caller's conflict filter (DRAGEN `PloidyConflict`) as Y
  evidence.
- The two clusters come from a few unlabelled people; a cut-off fitted to them is a number read off a
  data dump. Calibrate on PGP-HMS (self-reported sex on the profile, unverified; a second caller
  family) or cite a published rule (plink `--check-sex` F 0.2/0.8, unverified), with a minimum site count.
- "Ploidy model" interprets tools: GATK `sample_ploidy=2` is a default that per-interval options may
  override (unverified), DRAGEN's karyotype line may be given or estimated (unverified), DeepVariant
  1.6's `--haploid_contigs` does not reach the header (unverified, consistent with the 15 local files).
  Better: the header's ploidy lines verbatim, plus observed counts of one-, two- and n-allele GTs on
  autosomes, non-PAR X, non-PAR Y and MT.
- Depth is already in the frames (`dp`, `ad`, `blocks.dp_min`): X-het allele balance is the signal for
  contamination and mosaic 45,X/46,XX, which het share alone misses.
- There is no person in the schema. Evidence belongs per facet (a struct on `facets` or a small findings
  frame) and a verdict to a function returning a result with a reason (`no_build`, `y_unusable`,
  `caller_asserted_only`, `too_few_sites`, `intermediate`, `facets_disagree`, `not_wgs`); liftovers of
  one callset are one piece of evidence, not several (`derived_from`).
- `KaryotypicSex` belongs to a Phenopackets Individual and has no place for an inference: use its values
  inside a just-vcf-named inferred field, and say the result is not cytogenetic. An XX male gives XX
  correctly (karyotype, not sex). just-prs needs the X copy count for dosage; that may be the thing to
  expose first.
- LANDSCAPE names somalier for sex checks: say why this is built here, or take its sites. Arrays may do
  better than planned (female Y probes are no-calls; AncestryDNA doubles X letters).

Constraints: list length stays the number of haplotypes the source wrote (Decision 2). The finding sits
beside the facets and is never a rewrite of `calls`: no GT is collapsed to haploid by the loader. An
array facet without Y probes, or a targeted panel, says unknown.

Legality: minor. Filling `caller`/`caller_version`, a new nullable facet field and a new person-level
finding are all additive. Where the finding lives (a new frame or columns on an existing one) is a
schema decision for the lead. just-prs builds a stopgap now (karyotype from its normalized parquet,
caller and version as parquet metadata) and moves to this when it exists.

**Depth measured 2026-10-22 (lead), on the same 16 local private single-sample WGS VCFs** (aggregates
only, no file named). Medians of the calls' own FORMAT DP, each over the median DP of PASS calls on
chr20:

- non-PAR X, PASS calls: 0.42–0.50 in the low-X-het cluster (11 files), 0.89–1.00 in the high one (5).
  Nothing between.
- X-degenerate Y (GIAB v3.6 chrY nonPAR minus XTR minus ampliconic), PASS calls carrying an ALT:
  0.56–0.82 in the low-X-het cluster, 0.35–0.40 in the high one. Nothing between. Every row with DP in
  that region, PASS or not, overlaps less cleanly (0.68–1.06 against 0.35–0.57), so the ALT-call
  median is the one used.
- One sex-aware caller's file writes X haploid; its X depth ratio (0.48) agrees with its haploid calls,
  the independent check the critique asked for.
- Liftovers that dropped FORMAT carry no DP and give no depth line.

So X copies have two independent lines (X depth and the diploid X het share) that agree on all 16
files, and Y presence has one (the X-degenerate Y depth). Y-region ALT calls with DP were in the
thousands per file. Rerun through `sex_chromosome_evidence` on the 17 local single-sample VCFs it can
load: 12 facets read one X (X depth 0.419–0.500, Y depth 0.560–0.824) and 5 two X (het share
0.580–0.626, X depth 0.892–1.000, Y depth 0.351–0.400); no facet fell in an unknown band.

**PGP-HMS calibration, 2026-10-22** (consented-public, CC0, so per-file values are kept). 27 uploaded
single-sample VCFs with DP, read from the Krypton copy; labels are the participants' self-reported
`Sex/Gender` in `PGP_Participant_Survey.csv` of the HU-PGP-DB-2021-06-18 Arvados snapshot (14 files
labelled male, 4 female, 9 unlabelled; unverified). The profile pages could not be read: the site's
certificate names another host. Shipped rules, file by file (`x` and `y` are the depth ratios):

| Participant / upload | Label | Build, caller | X evidence | Y ratio | Answer |
|---|---|---|---|---|---|
| huE0F0CE/4084 | female | GRCh38 DRAGEN | share 0.593, x 0.926 | 0.222 | two X, Y not shown |
| huF7A4DE/3936 | female | GRCh38 GATK HC | share 0.618, x 0.938 | 0.125 | two X, Y not shown |
| hu21731A/2058, huAA16BD/825 | female | GRCh37 exomes | share 0.588, 0.467 | - | `not_wgs` |
| hu6ABACE/3885 | male | GRCh37 DRAGEN | share 0.041, x 0.446 | 0.607 | XY |
| huE18D82/3910 | male | GRCh37 DRAGEN | share 0.043, x 0.474 | 0.842 | XY |
| hu0486D6/1048 | male | GRCh37 Isaac | haploid X, x 0.560 | 0.560 | XY |
| huA2692E/2816 | male | GRCh37 Isaac | share 0.058, x 0.552 | 0.862 | XY |
| hu3BDC4B/3479, /3481 | male | GRCh37 GATK HC (bwa, bowtie2) | share 0.016, 0.012; x 0.484, 0.500 | 0.484, 0.433 | one X, Y not shown |
| hu2B76CC/2104, huF06AD0/1506 | male | GRCh37 freeBayes (no AD) | share 0.396, 0.378; x 0.500 | 1.450, 1.188 | share intermediate |
| hu80855C/1126, /1128 | male | GRCh37 Isaac | haploid X, x 0.552 | 162 sites | one X, Y too few sites |
| hu832966/3674, hu6ABACE/3884 | male | GRCh38 Helix | haploid X, x 0.500, 0.517 | 286, 251 sites | one X, Y too few sites |
| huFD37C8/3255 | male | GRCh37 exome | share 0.088 | - | `not_wgs` |
| huE2F2FC/2771 | male | starling | - | - | load failed: polars-bio "invalid BGZF ISIZE" |
| huBEF5E2, hu24385B, hu1C1368, huC39675, huBE0518, hu416394, hu2BA77F | unlabelled | GATK, DRAGEN, Isaac, bcftools | one X each | 0.548–1.970 | XY |
| hu876E3D/4118 | unlabelled | GRCh37 DRAGEN | haploid X, x 0.471 | 0.471 | one X, Y not shown |
| huCDFB7C/3206 | unlabelled | GRCh37 GATK UG | 33 diploid X calls | 0 sites | too few sites |

What it changed:

- **No Y threshold says absent.** Self-reported women read 0.125 and 0.222, the local two-X genomes
  0.351–0.400, but one man's GATK genome read 0.433 and 0.484 (two alignments of the same reads). The
  first build of this layer had an absent band at 0.44 and below and called him `XO`. Y presence is
  now only ever `True` (ratio 0.52 or more, which every labelled or clustered two-X genome stays 0.12
  under), and `XX` and `XO`, which need a Y shown absent, are not emitted.
- **A sex-aware caller's few diploid X calls are no het share.** huBE0518's DRAGEN wrote 103,154
  haploid X ALT calls and 1,315 diploid ones, all hets; read as a share they said two X against a
  one-X depth. The share now counts only where most X ALT calls were written diploid. Its header's
  `##referenceSexKaryotype=XY` agrees with the depth-based answer, an independent check.
- **freeBayes males read an intermediate het share** (0.378, 0.396; those files have no AD). The share
  fell in the unknown band, so they answer unknown although their X depth (0.500) says one X: a share
  in the band is not overruled by the other line. Whether freeBayes needs its own band is open.
- Isaac males sit at X depth 0.548–0.560, above the local 0.42–0.50 and under the 0.60 edge.
- The X-degenerate Y het share, also measured, does not transfer: 0.40–1.00 in local one-X facets,
  0.06–0.12 in hu3BDC4B's GATK files. It is reported, not used.
- Two Helix uploads labelled exome read `callset_type = wgs` (RM4), with only 3,211 and 4,214 chr20
  ALT calls; their Y line had too few sites anyway.

**Built 2026-10-22: the evidence layer** (`just_vcf.karyotype`, exported `sex_chromosome_evidence`,
`SexChromosomeEvidence`, `KaryotypicSex`, `XCopies`; docs/SCHEMA.md "Sex-chromosome evidence"). Per
facet: GT ploidy counts per group, the X het share on diploid-written calls with the AD allele balance,
the X and X-degenerate Y depth ratios against chr20, and the facet's `ploidy_statements` carried
verbatim. The estimate gives `x_copies` and `y_present` separately, then `XY`, `XXY` or
`UNKNOWN_KARYOTYPE` with a reason (`no_build`, `unsupported_build`, `not_wgs`, `too_few_sites`,
`intermediate`, `lines_disagree`, `caller_asserted_only`, `no_depth`, `y_unusable`; `y_reason` may also
be `y_not_shown`). `XX` and `XO` are not emitted (above).

Thresholds: each unknown band is the central half of the measured gap, so every class edge keeps a
quarter of the gap as margin. X het share: one X at most 0.23, two at least 0.46 (gap 0.115–0.580).
X depth ratio: one X at most 0.60, two from 0.79 to 1.25 (gap 0.50–0.89; above 1.25 fits three copies
better than two, so it is unknown, not two). Y depth ratio: present at 0.52 or more, otherwise not
shown, never absent (PGP-HMS, above). Minimum sites: 1,000 diploid X ALT calls for the share, 1,000 X and 1,000 chr20 ALT
calls with DP for the X ratio, 300 X-degenerate Y ALT calls with DP for the Y ratio. Where the het share
and the X depth both decide, they must agree. The share counts only where most X ALT calls were written
diploid; otherwise (a sex-aware caller) depth alone decides.

Decisions, each with what was rejected:

- Computed on demand, one result per facet. Rejected: a struct column on `facets` (a calls scan at
  every load, a schema bump that orphans the cache, thresholds frozen into data) and a person-level
  verdict (no person in the schema; the critique's advice).
- A `Loaded` with several facets (a freeBayes+VarScan2 merge) gets one answer per facet, in `facets`
  order. Rejected: a `facet_id` argument, which makes the common case a lookup.
- "Not failed" is `filters` null or `["PASS"]`, and the same for FORMAT `FT`. Rejected: literal PASS
  only, which leaves a GATK file whose every FILTER is `.` with no rows.
- The depth denominator is chr20, because the thresholds were measured against it. Rejected: all
  autosomes, which would move every ratio off its calibration.
- Rows carrying a `split_signal` stay out of the het share (a `0` there can be another ALT, RM270).
- The allele balance is reported, not used: no calibration says where a mosaic begins.
- X copies stay known when the karyotype is not (`y_unusable`, `no_depth`): what just-prs needs for
  dosage (S1).
- T2T-CHM13 answered `unsupported_build`, not `no_build`: it had a label but no region set. Since RM286
  it has region sets and answers `uncalibrated_build`, measured but not estimated.
- The header's ploidy lines are carried and the caller's haploid genotypes counted, never as evidence.
- Y presence is never `False`. Rejected: an absent band under a threshold, which named a man `XO`
  (above); and dropping the Y line altogether, which would lose `XY` and `XXY` where Y depth is high
  and the women measured sit far below.

**Built 2026-10-23: coverage as an X and a Y line** (RM278). `sex_chromosome_evidence(loaded,
coverage=…)` reads a coverage facet's length-weighted median depth over non-PAR X, the X-degenerate Y
and chr20 (`CoverageRatios`, on the result's `coverage` field). X: a third line, one X at a ratio of
at most 0.625, two from 0.875 to 1.25, and all measured X lines must agree. Y: present at 0.375 or
more, absent at 0.125 or less, unknown between; absent with one X names `XO`, with two `XX`. Each band
is the central half of the copy-number gap (one copy 0.5, two 1.0; a Y 0.5, none 0), the rule the
call-based bands followed. Decisions, each with what was rejected:

Calibration, 2026-10-23, on the four private BAMs of RM278's table (aggregates only, pooled, no file
named or classed): every X and Y coverage ratio lay within 0.034 of its nearest copy-number
expectation (0, 0.5 or 1.0), so every file fell inside a class band with at least 0.09 to spare; each
file's X and Y pair was one of the usual two (one X with a Y, two X without); chr20's even share was
0.972 or more. The X-degenerate Y coverage of the files without a Y read near 0, against 0.35–0.40
from the calls' DP on the local two-X VCFs: the call-site line sees the mismapped stretches only.
Four files from three person folders confirm the bands' placement, not their width.

- Coverage shows a Y absent where call-site depth could not: the coverage median runs over every base
  of the region, the calls' DP only over the sites a caller found something at, which in a genome
  without a Y are the mismapped stretches. Rejected: an absent band on the calls' Y ratio (it named a
  man `XO`, above).
- Coverage showing no Y against calls showing one is `lines_disagree`. Rejected: letting coverage
  win, which would hide a mispaired coverage facet (another person's BAM).
- Zero-depth intervals are left out on chr20 and X (assembly gaps; each span holds a centromere) and
  kept on the X-degenerate Y, which GIAB cuts free of gaps and where zero is the signal. Rejected:
  counting them everywhere, which lowers both medians by the gaps' share.
- A coverage facet whose chr20 depth is uneven (under 75% of its bases within half to twice the
  median: an exome's or a panel's shape) is set aside as `coverage_uneven`. Rejected: trusting the
  paired facet's `callset_type`, which says what the calls are, not what the reads are.
- The coverage facet is paired by the caller and used only on its own build. Rejected: placing it on
  the facet's build, which would read GRCh37 coordinates against GRCh38 regions.

**Left:** XXX, XYY and the other classes need a copy count past two, which the bands above do not
read (a Y ratio near 1.0 still reads as one Y present). An array's female Y probes are no-calls, which may be the absence signal
arrays have and genomes lack. Arrays (no `gt`
before RM3, no depth) and exomes answer unknown; an array's own lines (female Y probes no-calls,
AncestryDNA's doubled X letters) are not built. Liftovers of one callset are not yet weighed as one
piece of evidence (`facets_disagree` belongs to the person-level layer). GRCh37's last X-degenerate
interval may reach about 99 kb into NCBI's PAR2 (docs/PROVENANCE.md). The somalier question in the
critique is not answered.
