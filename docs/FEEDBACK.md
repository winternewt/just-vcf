# Consumer feedback — open inbox

The intake for repos and people that **consume** just-vcf: today `just-prs` and `just-dna-lite`, plus
anyone using the published package.

**An empty inbox means nothing is owed.** Once an item is answered it moves to
[FEEDBACK_HISTORY.md](FEEDBACK_HISTORY.md) verbatim, reply and all.

The runbook for triaging this file is [CONSUMER_TRIAGE_LOOP.md](CONSUMER_TRIAGE_LOOP.md).

## How to add an item

Append a top-level section at the end of this file:

```
## S<n> — <what happened, in one line>

What you were doing, what you expected, what you saw. Paths, inputs and versions.
```

Take `<n>` from `./scripts/triage-state.py --next`, which reads this file and the history file
together. Ids are never reused. Your prose is never edited, and the reply goes below it.

Next id: S1.

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
