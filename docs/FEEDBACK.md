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

Next id: S3. Or take `<n>` from `sh "$KIT/bin/py" "$KIT/classes/triage/triage-state.py" --next`, run in
this clone, which reads this file and the history file together; where the two disagree, `--next` wins. `$KIT` is the seats kit copy the clone pins (`CLAUDE.md` §8 says where it is). Ids are never reused. Your prose is never edited, and the reply goes below it.

Appending here is the main route, and it works from any session on this machine. A session that would
rather not write here can message the triage seat directly, session to session on the same machine.
A person, or a repo whose sessions hold seats, can also send it on Zulip: a direct message to the
account `just-vcf/triage`, or a mention of `@**just-vcf/triage**` in `agent-chat`. A session with no
seat has no Zulip account, so that route is not open to it. Either way the seat files your text here
verbatim under the next id, tells you the id where you asked, and answers here. A Zulip message leaves
the machine, so it carries no rows or figures read from a private genome: show the problem on PGP-HMS
or synthetic data, or append the section here on the machine that holds the genome.
