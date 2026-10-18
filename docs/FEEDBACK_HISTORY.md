# Consumer feedback — answered

Items moved out of [FEEDBACK.md](FEEDBACK.md) once answered, verbatim: the reporter's prose exactly as
filed, with the reply beside it. The reporter's text is never edited or re-wrapped on the way here.

**Ids are never reused**, including the ones archived here. The next id is computed from both files:
`./scripts/triage-state.py --next`.

## Contents

- **S1** male X/Y ploidy across callers and panels — accepted (RM269)
- **S2** pre-split multi-allelic rows read wrong — accepted (RM270)

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

**Status — correction, 2026-10-18.** The reply above says the same man "does load with a different
list length by vendor today". That is wrong: just-vcf has no AncestryDNA loader yet. `sniff`
recognises the file, and `load()` refuses it with `UnsupportedFormatError` ("no loader for
ancestrydna"); the AncestryDNA loader is part of RM2's vendor CSVs, not shipped. The two-letter
encoding comes from our `ancestrydna` format notes, not from a load. So today the vendor split shows
up only once that loader lands, and RM269 has to cover it then. Nothing else in the reply changes.
<!-- triaged: 0.1.0 · sha abd809824f31 -->

# From just-dna-lite, 2026-10-18

## S2 — VCFs split one row per allele read wrong at multi-allelic sites: ~72,000 sites in one genome

From the just-dna-lite session, 2026-10-18. A report, not a request; just-dna-lite does not depend on
just-vcf yet, so this is input for RM6 and for any loader that meets a pre-split file.

**What we ran.** Four callers' VCFs for one private WGS sample (GRCh38; DeepVariant, Clair3, GATK
HaplotypeCaller in DRAGEN mode, and a 3-caller consensus), annotated with just-dna-lite's pipeline. One
of the DeepVariant files is unsplit (`1/2` records); the others were split upstream into one row per ALT,
in the `bcftools norm -m-` convention. The headers do not say what split them.

**What splitting does to the genotype.** A site where the sample carries two different alts and no
reference:

```
unsplit:  chr1  1000  .  A  G,T  ...  GT  1/2
split:    chr1  1000  .  A  G    ...  GT  1/0
          chr1  1000  .  A  T    ...  GT  0/1
```

In each split row the `0` stands for "the other ALT", not REF. Nothing in the row says so. A reader that
maps GT indices onto `[REF] + ALT` (just-dna-lite's `io._compute_genotype_expr` does) gets `A/G` and
`A/T`, two heterozygous calls with a reference copy, where the truth is `G/T`. The same convention gives
false hom-ref rows: unsplit `1/1` at `A>G,T` becomes `1/1` (G) plus `0/0` (T), and the `0/0` reads as
two reference copies for a person who has none.

**What we measured** (rows grouped by chrom + start in the normalized parquet, after PASS/DP≥10/QUAL≥20):

- Unsplit DeepVariant file: 73,073 records with an ALT index ≥ 2 in GT; no position has two rows.
- Split DeepVariant file: 74,671 positions with two or more rows. 72,009 are exactly `0/1 + 1/0`, i.e. a
  split `1/2`. 600 are `0/0 + 0/1 + 1/0` and 3 are `0/0 + 0/0 + 0/1 + 1/0` (a `1/2` with a third,
  uncarried ALT), 4 are `0/0 + 1/1`. 1,687 are `0/0 + 0/1` (a split `0/1`; read correctly).
- Split Clair3 file: 48,992 positions with two or more rows, 48,987 of them all `0/1`/`1/0`.
- Split GATK file: 98,909 positions with two or more rows, 97,503 of them all `0/1`/`1/0`.

So roughly 72,000 genotypes per split file are reported with a reference copy the sample does not
have, and about 600 positions get a false hom-ref row. An earlier note in our session said "3 sites";
that was only the sites in one ClinVar report table, not the genome.

**What we expected.** Either the unsplit genotype, or a row that marks its `0` as "other ALT".

**What we are doing meanwhile.** Nothing in just-dna-lite's pipeline changed yet; it still reads the
split rows as above. The one report that shows these genotypes labels a `0/1 + 1/0` pair "two
different alts (multi-allelic)". A fix in our own GT reader is a separate decision on our side.

**What this may mean here** (for whoever owns RM6; no ask attached):

- SCHEMA Decision 1 already keeps one row per site at load, so an unsplit file is safe. A file that
  arrives already split is not: the loader cannot tell a split `1/0` from a real `1/0`.
- Detectable signal: two or more rows at one position whose GTs, taken together, use each allele slot
  once (`0/1 + 1/0`), or a `0/0` row beside a non-reference row at the same position. A
  `bcftools norm --old-rec-tag` INFO field, when present, says it outright.
- RM6's own "split multi-allelic sites" output would produce the same ambiguity if it writes this
  convention.

**Status — accepted, filed as RM270 (docs/ROADMAP.md); nothing shipped yet.** Reproduced against the
code on a synthetic VCF holding a split `1/2`, a split `1/1` and an unsplit `1/2`. The loader keeps one
row per record (SCHEMA Decision 1), and `gt_alleles` maps each `0` onto the file's REF: the split `1/2`
loads as `["G","A"]` and `["A","T"]`, the uncarried ALT's `0/0` as a hom-ref `["C","C"]`, with no
warning; the unsplit record loads as `["G","T"]`. So just-vcf reads a pre-split file the same wrong way
your reader does, and says nothing about it.

RM270 covers detecting a pre-split file (the signals you listed, read as evidence, since two records at
one position can also be legitimate) and deciding what the frame says about a `0` it cannot place: a
warning, a facet field, a withheld `gt_alleles`, or the split records joined back into the site as a
derived facet. Which one is a schema decision for the lead; the additive ones are a minor release, and
withholding `gt_alleles` changes values consumers read today, so the changelog would say so. It also
binds RM6: just-vcf's own split output must not write this convention.

What to do now: nothing is owed to just-vcf, since just-dna-lite does not load through it yet. For your
own reader, the unsplit file is the one that reads right today. Watch `docs/CHANGELOG.md` for RM270.
<!-- triaged: 0.1.0 · sha 3fa6a53b0ff4 -->
