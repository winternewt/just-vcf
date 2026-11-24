# Roadmap history
Shipped roadmap items, moved here from [ROADMAP.md](ROADMAP.md) with their rationale. Nothing is
deleted; it is relocated.

## RM269 — Karyotype finding and the caller's ploidy model
Severity: high · Status: shipped 2026-11-24 (the person chose the full label, desk #61, 2026-10-22. Shipped: half 1 as the `source` detector's caller, version and `ploidy_statements` (RM2, merge bd5c02e, 2026-10-22); the evidence layer, X copies and an XY/XXY estimate per WGS facet, merge 7384428, 2026-10-22, checked on 27 PGP-HMS VCFs; XX and XO from a coverage facet, merge dd11011 (RM278), 2026-10-23; XXX, XYY, XXYY and XXXY from a coverage facet's copy counts, and arrays' X from their hets with Y behind `array_y`, merge 5ccbffc, 2026-10-25, checked on 270 PGP-HMS exports; the build's source on the result, with tables that state no build estimated on `resolve_ref`'s output, merge e73b320, 2026-10-25, checked on 109 PGP-HMS tables; GRCh36's region set, its arrays estimated and its genomes measured only, merge 3a80fc0, 2026-10-25, checked on 250 PGP-HMS 23andMe exports and 69 resolved tables; an array's Y read by default per vendor, on for AncestryDNA alone (`ARRAY_Y_VENDORS`, desk #136 answered by the person, 2026-10-26); a person's facets weighed against each other, liftovers of one callset answering once, a merge or a consensus set aside and independent callsets that differ withheld as `facets_disagree` (`person_karyotype`, 2026-11-03, synthetic facets only, no real genome run). Decided 2026-11-24, desk #135's term ended on its default: a VCF alone still names XY and XXY, though it cannot rule out an extra Y, as docs/SCHEMA.md and `karyotype` already say. Waiting on the person: the lead's calls listed for veto in docs/LAYERS.md, layer 4 (`OTHER_KARYOTYPE`, the 3/4 array two-X expectation against PLINK's shape, `x_copy_count`/`y_copy_count`, GRCh36 genomes left `uncalibrated_build`). Buildable now: an answer to the critique's somalier question. Open, no desk number: whether GRCh37's last X-degenerate interval reaches about 99 kb into NCBI's PAR2 (docs/PROVENANCE.md); whether a measured answer (a genome's calls and depth, a coverage facet's counts) should outrank an estimated one (an array's letters) at the person level, where the lead's call, open for veto in docs/SCHEMA.md, ranks none. Exomes answer unknown) Closed 2026-11-24: desk #135's term ended on its default (a VCF alone still names XY and XXY, though it cannot rule out an extra Y), and the open GRCh37 PAR question is settled on hg19's own sequence: X and Y are identical over every base of NCBI's PARs and nothing past them, so GRCh37's PARs are now Ensembl's (NCBI's), GIAB's 100,286-base PAR1 tail (SRY among it) joins the X-degenerate Y and the 99,421 bases GIAB put before PAR2 leave it (merge f97e51d; suite 2942 passed, 78 skipped, 8 xfailed on both polars majors). Left as they are: the lead's calls in docs/LAYERS.md layer 4 and the person-level ranking in docs/SCHEMA.md, standing until the person vetoes one; GRCh38's one-base PAR start offsets, noted in docs/PROVENANCE.md and unchecked, since Ensembl's GRCh38 masks Y's PARs; `person_karyotype` still checked on synthetic facets only.

Motivating report: S1 (just-prs, 2026-10-11). Callers, reference panels and scoring tools each pick a
convention for male X and Y and rarely record it: a caller with a diploid-everywhere model (DeepVariant
1.1.0) writes a man's non-PAR X as `1/1` and `0/1`, a sex-aware one (DRAGEN) writes or filters it
differently, and arrays split by vendor (23andMe one letter; AncestryDNA the same letter twice, per the
`ancestrydna` format notes; that format has no loader yet, RM2). A consumer reading `gt_alleles`
cannot tell one copy from two.
*Corrected 2026-10-26:* AncestryDNA loads since 2026-10-22 (`loaders/table.py`, RM2's array side,
merge 53600af).

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
   *Corrected 2026-10-26:* `caller` and `caller_version` are filled since 2026-10-22 by RM2's `source`
   detector (b3caebc, merge bd5c02e), which also reads `##bcftools_callVersion`; the ploidy model became
   `facets.ploidy_statements`, the header's lines kept verbatim with no vocabulary, as the critique below
   recommends.
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

*Answered 2026-11-06 (the lead), the somalier point:* built here, and somalier's sites not taken.
somalier's own sex rule (`relate.nim`, `add_parents_and_check_sex`, read as text from its source)
is binary: male where the X het count over the X hom-alt count is below 0.05, female above 0.4, with
more than 10 X sites with depth, and nothing in between; it reads no Y and has no class for XXY, X0,
XYY or a mosaic, which are the cases this item exists for (desk #61's full label). It reads only its
own site list (published for hg19 and hg38; its size, recalled as about 17,000 sites, is not checked
here), where this
estimate reads every non-PAR X call of a genome, an array's letters, and a coverage facet's depth, on
GRCh36, GRCh37, GRCh38 and T2T-CHM13. Its het-over-hom ratio is the same quantity as this item's het
share, so taking its sites would only shrink the sample. somalier stays what it is here, the
relatedness tool (`just_vcf.somalier`), not the karyotype's source.

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

**Built 2026-10-25: copy counts past two from coverage.** The coverage line counts: each copy reads
half chr20's depth, so X reads one, two or three copies (0.5, 1.0, 1.5) and Y none, one or two (0,
0.5, 1.0), each band 0.125 either side of its value, the central half of every gap unknown, the rule
the earlier bands followed. Each line now allows a range of counts and the answer is their overlap
(`x_copy_count`, `y_copy_count`, nullable, beside `x_copies` and `y_present`, which keep their
meaning). New classes: `XXX`, `XYY`, `XXYY`, `XXXY`, and `OTHER_KARYOTYPE` for three X with two Y.
No genome of those classes was measured; the bands are copy number alone. Decisions, each with what
was rejected:

- The coverage bands narrowed to fit: two X was 0.875–1.25, now 0.875–1.125; one X gained a floor
  at 0.375; one Y was 0.375 and up, now 0.375–0.625. A ratio in the new gaps reads unknown where it
  named a class. None of the four calibration BAMs is near a new edge (each within 0.034 of its
  copy-number value). Rejected: keeping 1.25 as two's top, which leaves no unknown band before three.
- The calls' X depth above 1.25 allows three or more, where it read unknown, so a three-X genome's
  own calls do not veto a coverage three. It names no count itself: with no coverage line the answer
  is `copies_uncounted` (and `x_copies = two_or_more`, which is true). Rejected: leaving it
  `intermediate`, which makes XXX unreachable for any facet carrying DP.
- The calls alone still name `XY` and `XXY`, reading their Y "one or more" as one and a het share's
  "two or more" as two, as they did before. They cannot tell XY from XYY or XXY from XXYY (men's
  call-site Y depth spans 0.43–1.97 across callers). Rejected for now: answering unknown there, which
  turns every VCF-only answer unknown; it is the person's choice (below).
- Three X with two Y is `OTHER_KARYOTYPE`, the standard's own term for a karyotype it does not list.
  Rejected: unknown, which would say the evidence was weak when the counts are not.

**Built 2026-10-25: array facets.** A facet with `callset_type = array` gets lines from its letters
as written (`ArrayLines` on the result's `array` field; no REF or `gt` needed). X: the het share of
two-letter base calls on non-PAR X over the autosomes' het share. One X expects 0 (genotyping errors
only), two X 3/4 (Ne_X = 3/4 Ne_A, the neutral expectation), and the bands are the gap's central
half: one X at most 0.1875, two or more from 0.5625. An X written with one letter on most called rows
is the vendor's call (`caller_asserted_only`), never read, as a sex-aware caller's haploid GT is not.
Y: the called share of the X-degenerate Y probes over the autosomes' called share, absent at most
0.25, present from 0.75 (none 0, a Y 1; the gap's central half). Arrays cannot count copies.

Checked 2026-10-25 on PGP-HMS through `load` and the estimate (the 229 GRCh37 23andMe exports, the 41
AncestryDNA exports and 109 other array tables; labels as in the VCF table above, the survey's
`Sex/Gender` and `Anatomical sex at birth` and `demographics.tsv`'s `Gender`, unverified; aggregates
only):

| Vendor | X written | X het ratio | Y call ratio | Files |
|---|---|---|---|---|
| 23andMe | two letters | 0.72–1.52 (91; v5's autosomal het share is about half the others') | exactly 0 (91) | 91 |
| 23andMe | one letter | not read | 0.76–0.99 | 137 |
| AncestryDNA | two letters, women | 0.83–0.95 (8) | 0.008–0.099 | 8 |
| AncestryDNA | two letters, men | 0.0002–0.0059 (17) | 0.746–0.999 | 17 |
| AncestryDNA | two letters, unlabelled | 9 at 0.0002–0.0069, 7 at 0.76–0.94 | 0.0008–1.0 | 16 |

- **23andMe's sex chromosomes follow one switch, the vendor's.** Every one of its 228 exports with
  a usable Y line has either a two-letter X and every X-degenerate Y probe no-called, or a one-letter
  X and the Y called. No two-letter X read a low het ratio. Exactly 0 over 709–2,841 probes, where
  AncestryDNA's women call a few, reads as a mask. So for 23andMe the Y line and
  the letter count are the vendor's sex call, and only the two-letter X's hets are independent.
- **AncestryDNA calls what it measures.** Its women call a few Y probes (0.8–9.9%), its men's X carry
  a genotyping error's hets (0.01–0.2% of two-letter calls): neither is a mask.
- **The Y line is read only with `array_y=True`** (default off, `y_reason = vendor_asserted`): nothing
  in a file tells a mask from a chip that saw no Y. With it on: AncestryDNA, every labelled woman `XX`
  (8), 16 of 17 men `XY`, one man unknown (Y 0.746, in the band); 23andMe, 73 of 76 women `XX` (3
  have a one-letter X, the vendor's male call, and answer `caller_asserted_only`), 106 of 108 men
  `caller_asserted_only`, one too few Y probes, and one self-reported man `XX` (two-letter X at 0.91,
  Y masked: vendor and hets agree against the self-report). No AncestryDNA file read against its label.
  Off, every array answers unknown, `x_copies` filled wherever the X is written with two letters.
  *Corrected 2026-10-26:* the default is per vendor since desk #136: read for AncestryDNA
  (`ARRAY_Y_VENDORS`), off for every other vendor (the superseded decision below).
- 23andMe on GRCh36 (250 exports) answered `unsupported_build` (it has a region set since, below);
  FTDNA (103), FTDNA space (2), rsid tables (3) and one Final Report state no build at load and answer `no_build` as loaded. The
  estimate reads the facets frame, not what `load` alone wrote: on `resolve_ref`'s output, which
  carries the build F9 inferred from the rows, they are estimated (below, "Built 2026-10-25: the
  build's source"). One Final Report did not load (`TableLayoutError`).

Decisions, each with what was rejected:

- The two-X expectation is 3/4, not 1. Rejected: 1 (X as diverse as an autosome), whose bands 0.25
  and 0.75 match PLINK `--check-sex`'s F cut-offs 0.2 and 0.8 (F from allele-frequency expected het,
  which one export cannot supply); it puts 3 of the 106 two-X exports (23andMe v4, 0.72) in the
  unknown band. The 3/4 is population genetics, not read off these files; the person may prefer the
  PLINK-shaped edges.
- The Y line behind an argument, default off, for every vendor. Rejected: reading it by default
  (23andMe's is a mask), and a per-vendor list of trusted vendors (one measured property of
  PGP-HMS exports, which a vendor's next release can change silently).
  *Superseded 2026-10-26:* the person chose the per-vendor default on desk #136 ("yes, per-vendor
  decision"). `array_y` became `bool | None`, default None: read by default for a vendor in
  `just_vcf.karyotype.ARRAY_Y_VENDORS` (AncestryDNA alone, the one vendor checked and found to call
  what it measures), off for 23andMe, for every other vendor and for a facet whose `vendor` is null;
  `True` and `False` apply to every array facet as before, so an explicit caller keeps its meaning.
  Off, `y_reason` stays `vendor_asserted` and the detail says whether the default or the caller turned
  it off; no reason value was added, since explicit `array_y=False` callers already read
  `vendor_asserted`. The risk the rejection named stands: a vendor release can change the property
  silently, and `array_y=False` is the override. Run 2026-10-26 under the default through `load` and
  `sex_chromosome_evidence` on every PGP-HMS upload that sniffs as AncestryDNA or 23andMe (Krypton
  copy; aggregates only): the 41 AncestryDNA exports answer 15 `XX`, 23 `XY` and 3 unknown with the Y
  in the band (`intermediate`), which fits the table above (8 women and 7 unlabelled two-letter X
  exports give the 15 `XX`; 17 men and 9 unlabelled one-X exports give 23 `XY` and the 3 in the band,
  one of them the man at 0.746); the 480 23andMe exports (230 GRCh37, 250 GRCh36) read no Y: 474
  `vendor_asserted` (164 of them two-letter X, `x_copies = two_or_more`; 310 one-letter X,
  `caller_asserted_only`) and 6 too few sites. None failed to load.
- An array's absent site is unknown (charter 11); only probes the export writes count, and a probe
  it writes as a no-call counts as not called. So the Y share needs an export that writes its
  no-calls: one with no no-call row on the autosomes (a VCF made from an array, a table that drops
  them) has only called Y rows, where a woman's cross-hybridising calls would read as a Y, and answers
  `no_calls_not_written`. Rejected: the chip's probe set as the denominator (`just_vcf.probes`), which
  exists only for learned or configured chips.

**Built 2026-10-25: the build's source.** The estimate already read a build `resolve_ref` inferred:
it reads the facets frame, and `resolve_ref` writes the label F9's row scoring gave
(`assembly_label_basis = inferred`) into it. What was missing was saying so. The result now carries
`assembly_label_basis` (verbatim), `assembly_label_source` (`stated`, `contig_lengths`, `rows`,
`liftover`, read off the facet's columns by SCHEMA.md Decision 5's rule; null for a label no rule
accounts for) and `ref_seqcol_digest`, carried and deciding nothing; `no_build`'s detail points at
`resolve_ref`. Nothing in the estimate infers or defaults a build. Decisions, each with what was
rejected:

- The estimate takes `resolve_ref`'s output; it does not resolve. Rejected: a `reference` argument
  running F9 inside the estimate (a hidden REF read of about 13 s per FTDNA table, and a second path to
  a build the frames already record).
- The finer source is derived on the result. Rejected: widening `AssemblyLabelBasis` to name the three
  inferred sources, a dtype change and a schema bump for what the columns already tell apart. The
  person may prefer it stored.

Checked 2026-10-25 on the same 109 PGP-HMS tables, through `load`, `resolve_ref` against
`/data/just-dna-cache/reference/store-ncbi36` (GRCh36, GRCh37, GRCh38) and the estimate (aggregates
only): all 109 resolved, every one with `assembly_label_source = rows`. 69 read GRCh36 (FTDNA 66,
FTDNA space 2, rsid table 1) and answered `unsupported_build` (below: GRCh36), F9's split of these
files. Of the 40 on GRCh37, 32 FTDNA exports answer `too_few_sites`: 19 carry under 1,000 X rows (autosomal files) and 13
under 1,000 autosomal rows (X-only files, no autosomal het share to divide by). The other 8 give
`x_copies`: 5 FTDNA with too few X-degenerate Y probes, and 3 (2 rsid tables, the Final Report)
`vendor_asserted`, named with `array_y`. Against `store-chm13`, which has no GRCh36, the 69 were
refused (`no_fit`), as they must be. FTDNA's autosomal and X files of one person are two facets; a
person-level layer could pair them.

**Built 2026-10-25: GRCh36's region set.** GRCh36 (NCBI36, hg18) had no PAR bounds and no
X-degenerate Y, so its facets answered `unsupported_build`. Now:

- PAR bounds: NCBI36's published PARs as Ensembl 54, the last Ensembl release on NCBI36, records
  them (`assembly_exception`, type PAR): X and Y 1-2,709,520, X 154,584,238-154,913,754 against Y
  57,443,438-57,772,954. GIAB publishes no NCBI36 set and UCSC no hg18 `par` table. Checked on hg18's
  own sequence: X and Y identical at every base inside each PAR, chance agreement in the kilobase
  outside.
- X-degenerate Y: GRCh37's set lifted through UCSC's hg19ToHg18 chain, interval by interval from its
  first to its last base on chrY (the rule UCSC liftOver applies to a region; 117 bases that land on
  hg18 chrX dropped), kept between the GRCh36 PARs: 13 intervals, 9.79 Mb. Lifted back through
  hg18ToHg19, every aligned base lands inside the GRCh37 set.
- The sources are carried in `assets/` and `tests/test_karyotype.py` re-derives both constants;
  the whole derivation and its checks are in docs/PROVENANCE.md.

A GRCh36 whole genome is now grouped and measured and answers `uncalibrated_build`, as T2T-CHM13
does: no NCBI36 genome was among the calibration files. A GRCh36 array is estimated
(`ARRAY_ESTIMATE_BUILDS`): its bands are not measured thresholds, the build only places the regions,
and its exports read as the GRCh37 ones do (below). `unsupported_build` is left for labels with no region
set (NCBI35, a label a caller wrote into the frame). A coverage facet on GRCh36 or T2T-CHM13 is set
aside as `uncalibrated_build`, where it said `unsupported_build`.

Decisions, each with what was rejected:

- The clip at PAR2. GRCh37's last X-degenerate interval lifts to 57,377,045-57,542,858, across NCBI36's
  PAR2 start; its first 66 kb are kept, the 99 kb inside PAR2 dropped. Rejected: dropping the whole
  interval (66 kb of male-specific sequence lost), and keeping the lift as it lands (a coverage line
  would read PAR2, where a woman's X reads land too, as Y).
- The GRCh37 set is not changed here, although the lift corroborates its open question (below).
- Ensembl 54's PARs, not GIAB's GRCh37 PARs lifted: those put Y's PAR1 end at 2,809,806 and PAR2
  start at 57,542,859, about 100 kb past where hg18's X and Y stop being identical (GIAB's GRCh37 Y
  PARs are about 100 kb longer than NCBI's, docs/PROVENANCE.md).
- Arrays on GRCh36 estimated, genomes not. Rejected: holding the arrays back with the genomes (their
  rule has nothing measured on a build to transfer), and estimating GRCh36 genomes (the depth bands
  were measured on GRCh37/38 genomes, PGP-HMS has no GRCh36 genome to check them on, and the lifted
  region set is the same sequence but the alignments are not). The person may prefer either.

Checked 2026-10-25 on PGP-HMS through the library (`load`, for tables `resolve_ref` against
`store-ncbi36`, then the estimate; labels as above, unverified; aggregates only):

- **The 250 23andMe exports on GRCh36** (chip v1/v2 81, v3 164, unknown 5) follow the vendor's one
  switch as the GRCh37 ones do. 73 write X with two letters, every one labelled a woman: X het ratio
  0.83-1.02, every X-degenerate Y probe no-called (Y call ratio exactly 0 over 708 or more probes).
  173 write X with one letter (165 men, 1 woman, 1 conflicting label, 6 unlabelled): Y call ratio
  0.84-1.00. Answers: the 73 give `x_copies = two_or_more` and, with `array_y`, `XX`; the 173 answer
  `caller_asserted_only`; 4 have too few X calls to read (`too_few_sites`). No two-letter X read a low
  het ratio, and no export read against its label but the one woman whose X the vendor wrote with one
  letter (its call, not read).
- **The 69 tables `resolve_ref` places on GRCh36** (FTDNA 66, FTDNA space 2, rsid table 1): FTDNA
  writes no Y, so no Y line anywhere but the rsid table. 64 answer `too_few_sites`, FTDNA's split
  files (35 with under 1,000 X rows, 29 with under 1,000 autosomal rows); 3 FTDNA give `x_copies = one`
  (het ratio 0.006; one man, two unlabelled); 2 write X with one letter (`caller_asserted_only`), one
  of them the rsid table, whose Y calls 0.94 of its probes.

**Left:** whether a VCF alone may name `XY`/`XXY` when it cannot count the Y, or must answer unknown
until a coverage facet counts it, is the person's decision. The array two-X expectation (3/4 or
PLINK's shape) is the person's too (the AncestryDNA `array_y` default was answered on desk #136,
2026-10-26: per vendor, above), and so is whether GRCh36 genomes stay `uncalibrated_build`. Exomes answer unknown. Liftovers of one callset are weighed as one
piece of evidence since 2026-11-03 (below, "Built 2026-11-03"). GRCh37's last X-degenerate
interval may reach about 99 kb into NCBI's PAR2 (docs/PROVENANCE.md); the GRCh36 lift corroborates
it, since that stretch lifts inside NCBI36's PAR2, where hg18's X and Y are identical. The somalier
question in the critique is not answered.

**Built 2026-11-03: the person level** (`person_karyotype`, `PersonKaryotype`; docs/SCHEMA.md
"Weighed against each other"). The lead pinned the rules; the facet-level answers are not touched.

- A callset is the person's facets joined by one-parent lineage (`liftover_to`, `normalized`,
  `split_rejoined`, from `links` and `derived_from`) to one file as loaded; it answers once, by the
  member that says the most, the root first, else the nearest member, and its members listed. Members
  that differ withhold it (`members_disagree`).
- A merge or a consensus, and anything made from one, is set aside with a reason, never counted.
- No callset outranks another: different answers withhold the person's (`facets_disagree`), each
  callset's own answer listed; withheld answers count for neither side; only withheld callsets give
  `all_withheld`.

Decisions beyond the pin, each with what was rejected:

- A callset withheld on its karyotype still states its X and Y lines, and those can contradict:
  23andMe's two-letter X (`x_copies = two_or_more`, its Y not read) against a genome's `XY` is
  `facets_disagree`. Rejected: comparing named karyotypes alone, which would read that pair as the
  genome's `XY`.
- A callset whose members disagree withholds the person's answer too (`members_disagree`), even where
  other callsets agree. Rejected: letting it fall silent, which resolves its contradiction away.
- `x_copies` is weighed on its own by the same rule, so the person keeps an X count where the
  karyotype is withheld (just-prs's dosage, S1). Rejected: a person-level Y roll-up as a field (Y lines
  still veto the karyotype).
- Kin facets (`caller_split`, `appended_export`) are separate callsets, as the pin's lineage set
  implies.

Checked on synthetic facets only (tests/test_karyotype.py); no real genome was run through it.

## RM304 — `query` reads no row of a facet for a class it carries none of
Severity: low · Status: shipped, found 2026-11-23 by RM302's builder. Merge 357e5a1: a row of a non-array facet whose SNVs `class_gate` says it does not carry, or an array's row of another probe class, answers `class_not_carried` with `carried` None (new `PresenceBasis` member); rows that state the allele and `ref` blocks are read as before. Left, noted, not filed: an indel-only facet's leftover SNV row that states the asked SNV still answers `stated`, as concordance reads it, where RM302's panel read withholds it. Closed 2026-11-23.

`just_vcf.query` answers `carried=False` for an SNV from a row at its position (`SITE_GENOTYPED`)
without asking the class rules: an indel-only facet's hom-ref insertion row there, or a chip's `DD`
insertion probe (a one-base REF with no ALT), answers as if it had genotyped the SNV. The rule: a
non-array facet's row answers unknown where `class_gate(facet, "snv")` is in `CLASS_UNSAID`; an
array's row answers only for its own probe's class (`callset.reads_per_probe`, `row_class_expr`).
The block read needs no change, as RM303's builder found: `_presence` answers `carried=False` from a
`ref` block for an SNV only, which a one-base GFF block speaks for, and a caller's gVCF block stays
read as RM303 settled.

## RM303 — An array's row speaks only for its own probe's class
Severity: low · Status: shipped, found 2026-11-23 by RM301's builder. Merge 7ce6b7f: an array's row (`callset_type` array) is read only for alleles of its own probe's class (`row_class_expr`: an I/D code or a length change is indel, a hom-ref row over a one-base REF an SNV probe, a row whose REF was not read is read for every class), in concordance, consensus and `site_genotypes`, so a hom-ref deletion probe now counts against a deletion; no `content_filter` is inferred for an array; a GET-Evidence GFF's one-base `ref` blocks speak only for an SNV at their base (`PROBE_BLOCK_FORMATS`), a caller's gVCF blocks read as before. CONSENSUS_VERSION 9. Closed 2026-11-23.

RM299 to RM301 gate a facet's row by the classes the whole facet carries. For a whole-genome caller
that is the right unit: a caller that wrote an SNV row at a base would have written an insertion
anchored there too, had it seen one. An array is a set of probes, each for its own allele, and the
facet-level gate misses three cases. A chip with even one resolved non-REF indel counts `indel`, so
its SNV probe rows at insertion anchors read as zero copies again. A chip whose resolved `I`/`D` rows
are all hom-ref counts no indel, so a deletion probe that looked and saw none says nothing of the
deletion. And the VCF loader infers a `content_filter` from an array export's rows (a chip is mostly
SNVs, so it reads `snv_only`), a whole genome's notion that an array does not fit: `class_gate` then
silences a resolved indel probe row. The rule: an array facet's row is read only for alleles of the
row's own class (an SNV probe for SNVs at its position, an `I`/`D` probe for indels over its span),
in `concordance`, `consensus` and `site_genotypes` alike, and no `content_filter` is inferred for a
facet whose callset type is an array.

Blocks too. A GET-Evidence GFF turns each `REF` row into a one-base `ref` block beside `SNP` rows, so
its block at an insertion's anchor is read as zero copies of the insertion (`ref_block`), in
concordance, consensus and `site_genotypes`. Such a block is probe-shaped and speaks only for SNVs at
its base. A facet-wide gate was built and left unmerged (`rm302` branch, 4187b91): it also withheld a
gVCF's caller-written blocks when the file was cut to one class, which are true reference calls for
the cut class as well, and every block of a gVCF with no variant row. Read the GFF's blocks by their
probe's class, and leave a caller's blocks as they are read now.

## RM302 — Ancestry and relatedness read no row of a facet for a class it carries none of
Severity: low · Status: shipped, found 2026-11-23 by RM300's builder. Merge f14ce8e of 28e119c: `site_genotypes` reads no row of a facet whose SNVs `class_gate` says it does not carry (`SiteBasis.CLASS_NOT_CARRIED`, counted by ancestry). The block part, 4187b91 on branch rm302, was left unmerged: a facet-wide gate on blocks also withholds a cut gVCF's caller-written blocks, which are true reference calls for the cut class, and every block of a gVCF with no variant row; the GET-Evidence GFF's one-base blocks it was right for go to RM303. Closed 2026-11-23.

`site_genotypes` (step 1, which ancestry and relatedness share) reads any row at a panel SNV's
position whose two haplotypes are among the panel's two alleles. An indel-only facet's hom-ref indel
row there has haplotypes REF and REF, and would be read as zero copies of the panel's ALT; only the
absence reading (step 2) asks `class_gate` today. Narrow: it needs an indel-only file with `0/0` rows
at panel sites. The `ref` block readings in concordance and consensus are not class-gated either; a
filtered gVCF should leave a gap where the cut row was, so that is likely moot, and the test says so
either way.

## RM301 — A genotype table's variant classes are counted once its REF is read
Severity: low · Status: shipped, found 2026-11-23 by RM299's builder. Merge 8aa3a42: `resolve_ref` counts a table's classes with the VCF loader's counter, null where no row's ALT is known; RESOLVE_VERSION 15, so a store holding a table resolved by 14 refuses the re-resolved facet as another REF state. Closed 2026-11-23. Left: an array's row read by its own probe's class (RM303).

A 23andMe-style text export keeps `variant_classes` null even after `resolve_ref`, so RM299's rule
cannot tell whether it carries indels, and its SNV rows still read as zero copies of an insertion
anchored on them (the same array exported as VCF had its classes counted, and is fixed). Count the
classes for a genotype table once its alleles are known, as the VCF loader counts them; where they
cannot be known, they stay null.

## RM300 — Consensus reads no row of a facet for a class it carries none of
Severity: medium · Status: shipped, found 2026-11-23 by RM299's builder. Merge 70748f1: shown failing first (the array's row read as a contradiction moved a het insertion's posterior from 25 to 5 phred); such a row is now silent, `class_not_carried`, and CONSENSUS_VERSION is 6. Closed 2026-11-23.

`consensus` has its own reading of a facet's one row at an allele's position over the same REF
(docstring point 2: "zero copies, `site`"), parallel to the one RM299 fixed in `concordance`. It very
likely reads an array's or an SNV-only file's row at the base an insertion is anchored on as evidence
of no copy, which a consensus would then weigh against the caller that saw the insertion. Where
`class_gate` says the facet carries no allele of the class, its row says nothing of that class.

## RM299 — A facet that carries no indels reads no copies of one
Severity: medium · Status: shipped, found 2026-11-23 by RM279's last pairs. Merge ec9925b; the RM279 pairs rebuilt the same day: the array pair's 15,869 indels read as called by both are one-sided, the DRAGEN SNV file's 4,009 are 1,318 both files wrote (docs/METHOD_ATLAS.md). Closed 2026-11-23. Left: consensus's parallel reading (RM300), a genotype table's classes (RM301), site_genotypes' (RM302).

`concordance` reads the other facet's row at an allele's position when the row's REF is the allele's
(`single` in `_other_side`), and `_holds` then reads a haplotype equal to REF as holding no copy. That
is right for a whole-genome caller, which would have written the insertion in that row had it seen one.
It is wrong for a facet that carries no indels at all: an array export with no indel probes, or a
whole genome cut to its SNVs (`content_filter` `snv_only`). Its SNV row at the base an insertion is
anchored on says nothing of the insertion, yet the allele is counted as called by both and discordant.
On RM279's CGI `var` against an Omni2.5 array exported as VCF, all 15,869 indel alleles read so were
discordant, and none could have matched: the array file holds no indel row. CGI against a DRAGEN SNV
file read 4,009 such alleles, 78.7% discordant. Where the facet's absence says nothing of a class
(`absence_is_reference`: the class not in `variant_classes`, or cut by `content_filter`), its row
should say nothing of that class either: the allele stays one-sided, as it already is where that facet
has no row, and the reading is counted with its reason. Then the RM279 pairs' indel cells are rebuilt.

## RM279 — A method-error atlas
Severity: medium · Status: shipped (`method_atlas` and a first atlas over eight pairs shipped 2026-10-23, merge 4047968, docs/METHOD_ATLAS.md; the PGP CGI pairs rebuilt with REF read, minted and met by VRS id, merge 93dd5f0, with complex rows read by the allele's change (F21) and X's haploid calls matched to homozygotes (F22), merge c6da26c; counts per confidence band, merge 86daede, 2026-10-24; CGI `var` pairs banded by `calls.hq`, merge 8334138, 2026-10-24; the six private pairs rebuilt with minting through `scripts/method_atlas.py` (REF read, met by VRS id, banded), 2026-10-26, merged 2026-10-26. PGP-HMS pairs that pool more than one genome built 2026-11-08 (docs/METHOD_ATLAS.md, "More PGP pairs"): a third CGI `var` against freeBayes 0.9.15 person, the last on disk, and 16 people's two CGI sequencing runs, nine of them under one key (2.5.0.33 against 2.5.0.37) with easy-region SNV discordance from 0.13% to 1.84% between people, sent to desk #101. Single-person pairs of other methods (DRAGEN against freeBayes, CGI against DRAGEN, CGI against an Omni2.5 array as VCF) are on disk and planned in `data/interim/rm279/plan_pgp_methods.json`, not built: each adds a method with one genome, not a pooled row. Settled by default: pooling stays by caller and version, the `var` atlas a separate file (desk #101); the private folders are never paired (desk #72); CGI's sample-level VQLOW is not read by `carried_pass` (desk #102). The last single-person PGP pairs built 2026-11-23, three with rows, docs/METHOD_ATLAS.md, their indel cells left to RM299. Closed 2026-11-23: every planned pair built)

The person chose GIAB stratifications as regions (E8): error and discordance rates by caller and region from multi-facet PGP-HMS people (docs/LIVE_SAMPLES.md). Needs RM277.

Measured 2026-10-26, the six private pairs met by VRS id (docs/METHOD_ATLAS.md): every row of the six facets had its REF read and matched, every allele minted but GATK's 18,248 spanning deletions; the same-reads pairs' rates move by at most 0.02 point (F21's `unreadable`, 9 to 49 alleles per class); DeepVariant 1.1.0 spells about 14,000 indels unlike 1.3.0 and 1.8.0, which now meet, so its indel discordance rises from 7.1 to 8.0% and 8.0 to 8.7% and its one-sided indels fall from 8.8 to 6.0% and 9.2 to 6.3%; 1.3.0 against 1.8.0, 3.68 to 3.78% and 0.62 to 0.20%. SNV rates move by at most 0.01 point.

## RM280 — A consensus posterior
Severity: low · Status: shipped 2026-10-23, closed 2026-11-23 (`consensus` and `method_error_rates` shipped 2026-10-23, merge c1b0e1f, docs/SCHEMA.md "A consensus posterior", checked on PGP-HMS in docs/CONSENSUS.md; rates by the facet's own confidence band, merge 86daede, 2026-10-24, with the 60+ band traced to borrowed bands (9bf578e) and CGI `var` rows banded by `calls.hq` (8334138); desk #94–#96 built as options with every default unchanged, merge 960b6c6 (`BlockGrade`, GQ on reference, `AlleleFrequencyPrior`). Closed 2026-11-23: every question put to the person ended on its default: #94 blocks graded by their likelihoods, #95 GQ off for reference calls, #96 withdrawn 2026-11-16 (the prior an option), #99 bands borrowed across caller versions with a warning, #100 every GIAB class with the highest rate winning (the easy/difficult partition an argument), #101 pooling by caller and version, #102 FILTER only. The lead's calls listed for veto in docs/LAYERS.md stay listed there)

The person chose after the atlas (E7): a per-site posterior across one person's facets as a derived facet, with RM279's rates as priors; never before them.

## RM284 — Phasing as CisPhasedBlock
Severity: low · Status: shipped 2026-10-23, closed 2026-11-22 (`cis_phased_blocks` shipped 2026-10-23, merge 4129dc7; PSL, PSO and PSQ read and written back under VCF 4.5, with withheld blocks naming the allele that withheld them, merge 8b72a46, schema 0.0.15. Closed 2026-11-22: desk #90 ended on its default, so blocks stay whole phased stretches and one symbolic allele withholds its block; still outside the plan: a 4.4+ slot-0 block's phase on polars-bio keeping the leading phase indicator (P4 in docs/POLARS_BIO_PENDING_FIXES.md, not filed, the upstream-liaison seat's); the reading of a 4.4+ slot 0 that a PSL names as phased is listed for the person's veto in docs/LAYERS.md, layer 8)

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

## RM271 — VRS allele ids

Severity: high · Status: shipped 2026-10-22, closed 2026-11-21 (`mint_vrs_ids` and `vrs_withheld` shipped 2026-10-22 for VCF and array rows, merge b41bda5, equal to vrs-python on 19,000 real alleles. The opt-in test in `tests/test_vrs.py` that mints on the real store `JUST_VCF_TEST_REFERENCE_DIR` names and compares every slot with vrs-python (`scripts/vrs_oracle.py`) shipped 2026-10-26, merge 19a9753; measured below, 0 mismatch. Closed 2026-11-21: desk #77's term ended on its default, so the gtars gap G1 is not filed; it stays recorded in docs/GTARS_PENDING_FIXES.md and just-vcf keeps writing the VRS state itself)

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
  on the GRCh38 genome). *Amended 2026-10-26:* the runs are minted a group at a time, consecutive runs
  up to `SPILL_PART_ROWS` rows (`loading.bounded_groups`), one slice, one pass and one ids part per
  group, since a pass per run was a fixed cost per contig (F24, docs/PREVIOUS_ISSUES.md).
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
- *Measured 2026-10-26* (F24's amendment, library only, VRS spill fresh after a cached `load` and
  `resolve_ref`, every output hashing the same before and after): hu416394's bcftools gVCF, 9,651,143
  calls on 2,057 contigs, minted in 97.5–108.8 s a pass per contig, 70.5–75.0 s in 12 groups, at a
  sampled peak of 4.7 GB either way; its Manta file, 9,373 rows on 942 contigs, 18.6–18.8 s → 11.4 s. What
  remains is per allele (24.4 million gtars digest calls, the indels' justification) and 25 whole-contig
  fetches; docs/PREVIOUS_ISSUES.md, F24.
**Measured 2026-10-26**, the opt-in tests in `tests/test_vrs.py` on the store at
`/data/just-dna-cache/reference/store` (GRCh37 and GRCh38) and PGP-HMS uploads from the Krypton copy,
aggregates only; every expected id computed by vrs-python at test time:
- Synthetic alleles on real chromosome 22 (a 400 kb window, its repeats found at test time) and MT:
  SNVs, MNVs, insertions, deletions, delins, one and two units into homopolymer and di-, tri- and
  tetranucleotide repeats, part-unit deletions, a deletion anchored on the last leading N,
  multi-allelic rows, a lowercase ALT, a monomorphic row, MT written `MT` (GRCh37) and `chrM`
  (GRCh38), and a seeded random batch over the window and inside its repeats. GRCh37: 700 rows, 1,394
  alleles equal (697 REF, 120 SNV, 122 MNV, 273 insertions, 182 deletions); GRCh38: 694 rows, 1,382
  equal (691 REF, 118 SNV, 118 MNV, 278 insertions, 177 deletions). 11 alleles withheld on each, every
  one with SCHEMA's reason: `spanning_deletion`, `sv_no_vrs_model` (`<DEL>`, a breakend),
  `unspecified`, `not_sequence`, `ref_mismatch` (2), `ref_unchecked` (2 for REF past MT's end, 2 on
  a contig the collection lacks, `hs37d5` and `chrEBV`).
- A 23andMe-shaped build 37 export on real chromosome 22: hom-ref and het ids equal; `DI` is
  `array_code`, `--` is `ref_unresolved`.
- PGP-HMS, the smallest upload per build the store places, the first 3,000 chromosome 22 rows and up
  to 500 MT rows: GRCh37, 3,295 rows, 6,000 alleles equal (3,000 REF, 2,940 SNV, 27 insertions, 33
  deletions), the 590 MT alleles `ref_unchecked` (`resolve_ref` left that file's MT unread); GRCh38,
  3,032 rows, 6,113 equal (3,032 REF, 2,596 SNV, 233 insertions, 252 deletions), none withheld.
- Sensitivity: a deletion's `repeatSubunitLength` taken from the expanded REF instead of the
  deleted length fails both synthetic runs.
- Not reachable on the real store: `refused` and `no_contig`. After `resolve_ref`, REF was read from
  the very sequence minting rolls on, and neither collection holds a name twice; `_justify` is driven
  directly for `refused` in the synthetic suite.
## RM298 — The remaining readers take a leading BOM as their neighbours do
Severity: low · Status: shipped 2026-11-18 (merge 26ee0c9)

**Built:** a VCF opening with a BOM loads as the file without it, plain, gzip or BGZF, re-saved or not: the loader
reads past it and such a file always goes through the rewritten copy (`Header.bom`, `REWRITE_VERSION` 20; the copy
19 wrote put the header among the records). `rsid_table.column_line_ok` and `liftover.read_chain` decode with
`utf-8-sig`. Nirvana already read past it (`json.loads` on bytes); a test pins it. **Not upstream:** VCF 4.3 and 4.5
forbid a BOM, so polars-bio refusing one ("empty input") is conformant.

## RM297 — Coverage detection reads past a BOM
Severity: low · Status: shipped 2026-11-18 (merge ee8cb8a)
**Built:** `coverage._head_lines` drops a BOM from the first line, gzip or not, so six kinds detect as without it
(mosdepth summary, `samtools coverage`, `samtools depth` with `#CHROM`, thresholds, per-base and quantized BEDs
opening with a comment); a thresholds file with a BOM reads with `kind=` stated too.
## RM296 — An EIGENSTRAT `.ind` is read as its `.snp` is
Severity: low · Status: shipped 2026-11-18 (merge b5e1f94)

**Built:** `ancestry.eigenstrat._fields` reads both files' lines with `scan_lines`: a leading UTF-8 BOM dropped, a
line empty or whitespace alone skipped, fields split on runs of whitespace, a width error naming the line in the file.
`read_ind` and `read_snp` only pick their columns. A gzip `.ind` now reads as a `.snp` did (pinned by a test); a
bare-CR `.ind` is refused, not misread; invalid UTF-8 raises polars' `ComputeError`.
## RM295 — A blank line in an EIGENSTRAT `.snp` is no site
Severity: low · Status: shipped 2026-11-18 (merge b5e1f94)

**Built:** `read_snp` skips a blank or whitespace-only line, as `read_ind` did, so a `.snp` ending in one no longer
fails the panel build on the `.geno` header's site count; `AadrBuild.n_snps`/`n_individuals` and `read_genotypes` say
they count rows, not lines.

## RM294 — A text line holding 0x1F is read
Severity: low · Status: shipped 2026-11-18 (merge 47edcd6; full suite 2819 passed on polars 2.1.0 and 1.44.2)
**Built:** `loaders.text.lines` and `eigenstrat.read_snp` read lines with `scan_lines`; `lines()` decodes invalid UTF-8
with replacement first, the same U+FFFD per bad sequence as `utf8-lossy` gave. **Left:** `loaders/table.py` and
`consensus.py` still use 0x1F as a join key (a collision needs a field holding it).
## RM293 — Coverage's read and its row count agree on quotes
Severity: low · Status: shipped 2026-11-18 (merge 47edcd6)

**Built:** `coverage._read_fields` passes `quote_char=None`, so a double quote is text to the read as to its count of
cut and padded rows. A field written `"chr1"` keeps its quotes; a quoted number in a typed column fails the read.
## RM291 — Coverage rows cut or padded are counted
Severity: low · Status: shipped 2026-11-17 (merge 77cbd0b; full suite 2799 passed on polars 2.1.0 and on 1.44.2)

**Built:** `coverage._read_fields` counts each row's fields from the same bytes it reads (`scan_lines`,
`str.count_matches`), summed across blocks, and the facet's warnings say per reason how many rows were wider than the
format (fields past it not read), short of it (missing fields null), or blank (all-null rows); a samtools `depth` file of
several alignment files says the first depth column was read and how many there were. About 12–14% slower on a 90 MB
per-base BED. **Left:** the read honours `"` quoting while the count splits on raw tabs, so a quoted field holding a tab
would count as a wide row (no coverage tool writes quotes; `quote_char=None` would align them).
## RM292 — `read_bed` skips track lines
Severity: low · Status: shipped 2026-11-17 (merge 77cbd0b)

**Built:** `atlas.read_bed` reads three named fields with `schema=`, `truncate_ragged_lines` and `coverage._RAGGED`, no
quote character, and drops lines matching `^(track|browser)(\s|$)` wherever they stand; a BED of comments alone reads as
no intervals on both polars majors. GIAB reads are unchanged.

## RM289 — Coverage BED rows of mixed width on polars 2
Severity: low · Status: shipped 2026-11-17 (merge aac313d; full suite 2740 passed on polars 2.1.0, the touched tests on 1.44.2 too)
**Built:** coverage's two reads go through `_read_fields`: `schema=`, `truncate_ragged_lines=True`, and on polars 2.x
`extra_columns="ignore", missing_columns="insert"` (`coverage._RAGGED`), so fields past a format's columns are dropped
and missing ones null on both versions, at the parser's speed. A file of no rows (a header or comments alone) is a facet
of no intervals on both. **Left:** the cut and padded rows are not counted (RM291); `atlas.read_bed` fails on a `track`
line (RM292); `loaders.text.lines` and `eigenstrat.read_snp` refuse a line holding `\x1f`, alike on both versions.
## RM290 — An unclosed file in the VCF loader
Severity: low · Status: shipped 2026-11-17 (merge aac313d)

**Built:** the re-save undoing stream (`_Undone`) owns its text wrapper and detaches it on close, so a payload read to
its end is no longer closed under its opener when the wrapper is collected (the `ResourceWarning`, a closed handle as
well as a warning); `_vcf_text` closes the stream before the payload on every exit. The full suite shows no
`ResourceWarning`. `-W error::ResourceWarning` does not turn a collection-time warning into a test failure, so the test
asserts the handle and the warnings directly.

## RM288 — polars-bio 0.37 and polars 2.x
Severity: medium · Status: shipped 2026-11-16 (merge 0adf89c; full suite 2706 passed on polars 2.1.0, 2705 on 1.44.2)
polars-bio 0.37.0 (2026-11-07) allows polars `>=1.40,<2.2`, which removed the reason the lock stayed on polars 1.44.2.
Built with no question put to the person (CLAUDE.md §10, 2026-11-16: a choice with only pros is not a decision).

**Built:** the floor is `polars-bio>=0.37.0` with no direct polars pin, so a consumer on polars 1.4x co-installs; the
lock is polars-bio 0.37.0 and polars 2.1.0. The key=value and CGI `var` loaders split lines with
`loaders.text.tab_fields` (polars 2 refuses a headerless ragged `scan_csv`); VRS minting's horizontal `concat` passes
`strict` only to polars 1.x. CLAUDE.md §5 says, per polars 1.44 quirk, what 2.1.0 does (three fixed, two still so, one
default changed; the workarounds stay while 1.x is admitted) and names three polars 2 changes met. P16 is fixed
upstream, but the CLI keeps logging off for Python's last-resort handler, pinned by a subprocess test. **Left:** drop
the 1.x workarounds and `vrs._EQUAL_HEIGHTS` once the floor is polars 2; `coverage.py`'s `read_csv` with `new_columns`
would refuse a BED of mixed row widths on polars 2 (left as is); a `ResourceWarning` from `loaders/vcf.py` (an
unclosed file, on both versions).

## RM276 — A durable per-person store
Severity: medium · Status: shipped 2026-11-05 (the lead kept a migrated VCF facet's null `method_steps` and `gl_source`, told apart by the manifest's record and `loaded()`'s warning, docs/SCHEMA.md "Store migration"; was: in progress, decided by the person in the Part 2 interview, 2026-10-23; `Store` and `write_vcf_zarr` shipped 2026-10-23, merge 6830fb7; `call_HQ` in the Zarr export 2026-10-24, merge fbcb3a9; `variant_CIPOS`, `variant_CIEND` and `variant_CN` in the Zarr export 2026-10-26, an interval holding -1 or -2, the spec's integer sentinels, written missing and counted; 2026-10-26, `migrate_store`: a store of an older schema version rewritten to this one by a step per bump since 0.0.14 (`just_vcf.store_migration`), written beside, checked and swapped, what it cannot carry named as filled unknown, docs/SCHEMA.md "Store migration", merge 401cbf0. Nothing buildable left. For the lead: whether `method_steps` and `gl_source` null on a migrated VCF facet (not read) should get its own value instead of sharing null with "not a VCF" / "no likelihoods", which only the manifest and `loaded()`'s warnings tell apart now) · Owner: lead
The person chose both (E13): a directory per person holding parquet per facet, the links frame and a manifest, read lazily; and a VCF Zarr export (its spec has no gVCF blocks, so blocks stay parquet-only). Facets accumulate; nothing is overwritten (invariant 5).
**Built:** `Store` (docs/SCHEMA.md, "The durable store"): parquet per facet, links per facet, a
manifest, a lock, lazy reads into `Loaded` and `Person`; and `write_vcf_zarr`, one facet's calls as
VCF Zarr 0.4, blocks counted and left in the store. **Left:** a migration between schema versions (a
store of another version is refused, not migrated); `cipos`, `ciend` and `alt_copy_numbers` in the
Zarr export; whether `Person` should carry the links (`store.loaded().links` has them).
*Corrected 2026-10-26:* `cipos`, `ciend` and `alt_copy_numbers` are written (docs/SCHEMA.md, "VCF
Zarr export").
*Corrected 2026-10-26:* `Person` carries `links` (merge b8bef07, 2026-10-23, "Person and the CLI keep
links"), and `Store.person()` builds it from `loaded()`, so its links are the store's; that question is
closed.
*Corrected 2026-10-26:* the migration is built: `migrate_store(path)` rewrites an older store by the
steps in `just_vcf.store_migration` and reports what it filled unknown (docs/SCHEMA.md, "Store
migration").
## RM277 — Concordance between facets
Severity: medium · Status: shipped 2026-10-23 (merge 2446c53; F21 and F22 followed the same day, merge c6da26c; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose graded where possible (E2): GL/GP compared where both facets carry them, hard genotypes otherwise, and every row says which. Kept as data, never resolved to a winner (invariant 8). Needs RM273 and VRS ids (RM271).

## RM273 — A person (subject) in the schema
Severity: high · Status: shipped 2026-10-23 (subject and Person, merge 93745c2; RM274's `same_person` as the default check, merge 6dd3591; was: in progress, subject and Person shipped 2026-10-23, the default check RM274's `same_person`) · Owner: lead
The person chose "Both: supplied, checked" (E1): the caller supplies a subject id when loading; nothing is inferred and facets are never regrouped by the library. The library checks the claim with relatedness (RM274) and warns on a mismatch. Most of the extended shape sits on this: RM275, RM276, RM277, RM280.
## RM285 — The format registry as a spec
Severity: low · Status: shipped 2026-10-23 (merge d0edd91; was: built 2026-10-23, docs/FORMAT_REGISTRY.md generated by scripts/format_registry.py, to move here on merge) · Owner: lead

The person chose this repo (E14): the declared formats and the broken shapes (F1–F17) documented as a readable spec in docs/ with synthetic examples; no genome data.

## RM272 — CLI
Severity: high · Status: shipped 2026-10-23 (merge 5b8700f; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead
Typer, registered in `[project.scripts]`, loading `.env` with python-dotenv (CLAUDE.md §4): `just-vcf sniff`, `load` (a summary of the facets), `convert` (to VCF/gVCF via `write_vcf`, optionally through `resolve_ref`, `normalize`, `liftover`, `mint_vrs_ids`) and `check`. Output is counts and aggregates only, never rows: the CLI is a privacy boundary like E12.
## RM275 — The facet graph as a links frame
Severity: medium · Status: shipped 2026-10-23 (merge c503f57, schema 0.0.13; was: open, decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead

The person chose a small links frame (E4): (from, to, relation, provenance), so a facet can have several parents (a merge, a consensus). `derived_from` and `relation` stay on `facets` for the one-parent case. Relations verbatim from GA4GH where it has them (`liftover_to`), just-vcf's own otherwise, documented.

## RM283 — An aggregate-only query helper
Severity: medium · Status: shipped 2026-10-23 (merge 95d8b61; was: built on a builder branch 2026-10-23, `just_vcf.query`, waiting for the lead's merge; decided by the person in the Part 2 interview, 2026-10-23) · Owner: lead
The person chose a helper beside `load` (E12): an API that answers counts, shares and presence questions and never returns rows, for agents and apps; `load` stays for code that may read rows.
## RM268 — Cache eviction
Severity: medium · Status: shipped 2026-10-22: the cache index, the stale-version sweep, `clear_cache`; no size cap · Owner: lead

Nothing ever removes anything from `JUST_VCF_CACHE_DIR`. It holds two kinds of entry:

- `vcf/<key>/`: rewritten copies of inputs polars-bio cannot read as they are
  (`loaders/vcf.py` `readable_path`), plus a `structural/` copy for a file polars-bio refused. The key
  hashes the source digest, the archive member and `REWRITE_VERSION`.
- `frames/<facet_id>-<SCHEMA_VERSION>/`: the calls and blocks parquet spill from `load`, and since
  RM270 the calls spill of a `rejoin_split` facet, under the derived `facet_id`. Nothing in the cache
  ties that id to its parent, so an eviction that keys on a source file does not see it.
Because the version is part of each key, every bump of `REWRITE_VERSION` or `SCHEMA_VERSION` orphans
every earlier entry, and nothing collects them. Measured 2026-10-03: after the 2026-09-30 dogfooding
runs, the cache under `data/interim/cache` had reached 32 GB in 36 `vcf/` entries. Most of it was
orphans: the two largest were plain `input.vcf` copies of 14.5 GB and 7.4 GB, from before the rewrite
switched to BGZF `input.vcf.gz`. The user deleted it by hand.
Constraints on a fix:
- A copy has to outlive the lazy frames that read it, so evicting the entry a live `Loaded` points at
  breaks a later `collect`. Orphans by version are always safe to remove; current entries need an age
  or size policy.
- The copies are whole private genomes. Deleting a source under `/data/genomes` leaves its copy in the
  cache, so eviction is also how a removed genome actually goes away.
Options, not decided: sweep stale-version entries on each load; an explicit `clear_cache()` in the
API (and in the CLI once it exists); a size cap with least-recently-used eviction.

**Shipped 2026-10-22** (the lead's outline, 2026-10-22). Rules in docs/SCHEMA.md ("The cache under the
lazy frames"); code in `just_vcf/cache.py`.

1. **The index.** Every entry has `<name>.json` beside it: kind (`rewrite_copy`, `frames_spill`,
   `rejoin_spill`, `resolve_spill`, `vrs_spill`, `normalize_spill`), version keys, the root source's
   sha256 and resolved path, `facet_id` and `parent_facet_id`, size, `created` and `last_used`
   (ISO-8601 UTC). A derived spill takes its root source from the facet row, which every derived step
   carries down (`source_sha256`, `source_path`), so `rejoin_split`'s spill is no longer cut off from
   its file.
2. **Stale versions are swept** by the step that owns the kind, before it writes, with the bytes freed
   logged. An entry from before the index is judged by its name: a frames name carries
   `SCHEMA_VERSION` (and the resolve or mint version); a `vcf/` key hashes its version, so
   `REWRITE_VERSION` went to 7 with the index and a copy without a manifest is stale by construction.
3. **`clear_cache(source=None, older_than=None, dry_run=False)`** returns a `CacheClearReport`
   (`removed`, `being_written`, `unknown_origin`, `bytes_removed`; no source paths, so it can be shown).
   Writers hold an OS file lock (`flock`; `msvcrt` on Windows, unexercised) beside the entry, and
   `clear_cache` skips an entry whose lock is held. A waiting writer re-checks that the lock file is
   still the one at its path, since a remover unlinks it.
4. **No size cap.** Rejected: LRU under a cap. `last_used` records writes and reuses, not reads, so a
   cap could evict the files a live frame scans with nobody asking.

Decisions taken while building, each with the rejected alternative:

- The manifest sits beside the entry, not inside it: `resolve_ref`, `normalize` and `mint_vrs_ids`
  remove their spill directory before writing it, which would take a manifest or lock inside with it.
- Bump `REWRITE_VERSION` rather than keep unindexed copies until someone deletes them by hand: the
  first load of each file after the upgrade rewrites its copy once, and in exchange no whole-genome
  copy stays behind unattributed.
- `source` resolves as: an existing file is hashed and matches by sha256 or by recorded path; a
  64-character hex string that is no file is the sha256; any other path matches the recorded path (a
  deleted file). Rejected: sha256 only, since the person who deleted a genome knows its path, not its
  digest.
- `older_than` reads `last_used`, not `created`: a copy reused every day is not old.
- `clear_cache()` with neither filter clears every indexed entry, with `dry_run` as the valve.
  Rejected: refusing without a filter, which would make clearing the cache take
  `older_than=timedelta(0)`.
- The manifest is written when the lock is taken as well as when the write ends, so an entry whose
  writer died still names its source and can be cleared by it.

Measured 2026-10-22 on three PGP-HMS VCFs (one that needs a rewrite copy) in a scratch cache under the
builder's worktree: 4 entries, 48,099,333 bytes (3 load spills, 1 copy). A copy of one spill under an
older schema name, without a manifest, 34,780,557 bytes, was swept by the next `load`.
`clear_cache(source)` on the file with the copy removed its 2 entries, 10,202,328 bytes, dry run and
real run agreeing; `older_than=timedelta(days=1)` removed nothing; `clear_cache()` removed the last 2.

Open: the derived steps still remove and rewrite their spill in place, so running one twice on the
same facet breaks the frame the first run returned. That predates RM268 and is not changed here.
Closed later the same day: a derived spill is now one generation per content key (docs/SCHEMA.md,
"The cache under the lazy frames"). `liftover`'s spill and `load`'s own followed the same day, and the
Complete Genomics and key=value loaders' working copies came under the index (`work_copy`).
Amended 2026-10-23 (the lead's decision, after two sessions at two `SCHEMA_VERSION`s on one cache dir
broke the older one's live frame): "orphans by version are always safe to remove" was wrong while an
older install still runs. The sweep now keeps an older entry used within `JUST_VCF_CACHE_GRACE_DAYS`
(7 days) or whose lock is held (docs/SCHEMA.md, "The cache under the lazy frames").

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
