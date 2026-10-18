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

**Or send it on Zulip**, when you are not in this checkout: a direct message to the account
`just-vcf/triage`, or a mention of `@**just-vcf/triage**` in `agent-chat`. A Claude session on the
same machine can message the triage seat directly instead. The seat files your text here verbatim under
the next id, tells you the id where you asked, and answers here, with a pointer where you asked. A
Zulip message leaves the machine, so it carries no rows or figures read from a private genome: show
the problem on PGP-HMS or synthetic data, or append the section here on the machine that holds the
genome.
