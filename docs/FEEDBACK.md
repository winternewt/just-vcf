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

Appending here is the main route, and it works from any session on this machine. A session that would
rather not write here can message the triage seat directly, session to session on the same machine.
A person, or a repo whose sessions hold seats, can also send it on Zulip: a direct message to the
account `just-vcf/triage`, or a mention of `@**just-vcf/triage**` in `agent-chat`. A session with no
seat has no Zulip account, so that route is not open to it. Either way the seat files your text here
verbatim under the next id, tells you the id where you asked, and answers here. A Zulip message leaves
the machine, so it carries no rows or figures read from a private genome: show the problem on PGP-HMS
or synthetic data, or append the section here on the machine that holds the genome.

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
