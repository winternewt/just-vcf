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
