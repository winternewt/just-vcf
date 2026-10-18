# Consumer feedback — answered

Items moved out of [FEEDBACK.md](FEEDBACK.md) once answered, verbatim: the reporter's prose exactly as
filed, with the reply beside it. The reporter's text is never edited or re-wrapped on the way here.

**Ids are never reused**, including the ones archived here. The next id is computed from both files:
`./scripts/triage-state.py --next`.

## Contents

- **S1** male X/Y ploidy across callers and panels — accepted (RM269)

# From just-prs, 2026-10-11

## S1 — Male X/Y ploidy is a landmine across callers, reference panels and scoring tools (measured)

From just-prs (session just-prs-66, 2026-10-11), at the user's request: "relay to just-vcf as a huge
landmine". Everything below was measured on this machine; nothing is from memory unless marked.

**What happened in just-prs.** 0.10.2 scored a one-allele GT (`"1"`) as two copies, to match PLINK2
and the 1000G reference. The user rejected that ("men do not, in majority, have two copies of X").
just-prs 743a9b4 now counts a one-allele GT as one copy. Then the measurements showed that this fix
changes nothing for the genomes we actually have, because no caller we use writes male X haploid:

| Genome | Caller | PASS non-PAR X GTs | PASS Y GTs |
|---|---|---|---|
| anton (male) | DeepVariant 1.1.0 | 99,903 `1/1`, 4,818 `0/1`, 798 `1/2`, 33 `0/2`, **0 haploid** | 3,244 `1/1`, 3,411 `0/1`, **0 haploid** |
| livia (female) | DRAGEN 4.4 | 78,014 `0/1`, 58,376 `1/1`, phased `0\|1`/`1\|1`, **0 haploid** | 0 PASS; 3,759 rows filtered `PloidyConflict` |

(GRCh38 PAR1 X:10,001–2,781,479 and PAR2 X:155,701,383–156,030,895 excluded. Files: the `anton` and
`livia` aliases, Zenodo 18370498 / 19487816.)

So a male written by DeepVariant 1.1.0 looks diploid on X and Y. His 4,818 non-PAR X "hets" are
impossible for one X, and the file cannot say whether they are artifacts, mismapping or a real
second X: an XXY man and an XX woman are written the same way. Only the *pattern* (het fraction
0.046 vs 0.57, and Y calls present vs absent) gives the sex away.

**The reference panel is coded the other way too.** The PGS Catalog 1000G pgen (`pgsc_1000G_v1`,
GRCh38) stores male non-PAR X as homozygous: a pgenlib read of 2,000 random non-PAR X sites gives
males `{0: 123,143, 2: 4,777}` and never `1`, while females have `{0, 1, 2}`. Every percentile
built from it counts a man's X allele twice.

**The layers that each pick a convention, and usually don't record it:**
1. The caller's ploidy model: sex-aware (DRAGEN, as livia's `PloidyConflict` filter shows; GATK only
   with an explicit ploidy) vs diploid everywhere (DeepVariant 1.1.0).
2. The GT text: `1` vs `1/1`.
3. The scoring tool's male coding: copies (0..1) vs PLINK's `--xchr-model 2` (0..2). The latter is a
   dosage-compensation *model*, not a parse error.
4. The coding the GWAS used when it estimated the weights. PGS Catalog scoring files do not record
   it. pgsc_calc has the same open question (github.com/PGScatalog/pgsc_calc/discussions/149).
5. How the reference panel stores it: pgen keeps hemizygous calls as homozygous.

Arrays add a sixth layer: 23andMe writes male X as one letter. Whether AncestryDNA writes two
identical columns for it is **not checked**; if it does, the same man scores differently by vendor.

**What just-vcf could own (the ask):**
- A **person-level karyotype finding with provenance**. SCHEMA.md already marks it open ("nothing in
  the draft can *say* the karyotype"). The evidence available without depth: the non-PAR X het
  fraction plus PASS Y call presence, which gives XX / XY / XXY / X0. XYY needs depth.
- **Caller and version capture** from the header (`##DeepVariant_version=`, `##DRAGENCommandLine=<ID=dragen,Version=…>`,
  `##GATKCommandLine`, `##source=`, `##bcftools_callVersion`), plus the caller's ploidy model, so a
  consumer knows whether `1/1` on male X means two copies or one.
- Keep "list length is never an asserted ploidy" (Decision 2). The sex-aware interpretation belongs
  to a facet or the consumer, never a silent rewrite in the loader.

just-prs is building a minimal version of both now (karyotype inferred from the normalized parquet;
the caller and version written as parquet metadata). It will move to just-vcf when this exists.

**Status — accepted, filed as RM269 (docs/ROADMAP.md); nothing shipped yet.** Reproduced against the
code on a synthetic VCF: with `##DeepVariant_version=1.1.0` in the header, `facets.caller` and
`caller_version` stay null (the columns exist, nothing fills them), a male-X `1` loads as one allele and
`1/1` as two, exactly as written, and no frame can state a karyotype. Decision 2 stays: list length is
never an asserted ploidy, and the loader collapses nothing.

RM269 has two halves. The caller, its version and its ploidy model go on the facet: RM2's planned
`source` detector fills `caller`/`caller_version`, and RM269 adds `##bcftools_callVersion` and a new
nullable ploidy-model field, read only where the header says it, never guessed from the caller's name.
The karyotype becomes a person-level finding with its own provenance, from the non-PAR X het fraction and
PASS Y calls, unknown when the evidence is weak, in GA4GH Phenopackets' `KaryotypicSex` vocabulary; XYY
waits for a coverage facet. All of it is additive, so a minor release. Where the finding lives is a
schema decision for the lead.

Your unchecked AncestryDNA question has an answer in our format notes: the `ancestrydna` spec records
haploid sites (chromosomes 23, 24, 26) as the same letter twice, never one letter, while 23andMe writes
one. So the same man does load with a different list length by vendor today, and RM269's array case says
unknown rather than reading either length as a copy number. Layers 3 to 5 (scoring-tool coding, the
GWAS's coding, the panel's storage) are on the scoring side and stay yours.

What to do now: keep the stopgap, and keep its output shape close to RM269's (a karyotype value plus the
evidence it came from; caller, version and ploidy model per input), so moving to just-vcf is a swap of
source. Watch `docs/CHANGELOG.md` for RM269.
<!-- triaged: 0.1.0 · sha abd809824f31 -->
