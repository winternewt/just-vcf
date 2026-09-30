# Seat: upstream liaison

The seat that owns just-vcf's relationship with the code it stands on but does not own: polars-bio and
the two layers under it, datafusion-bio-formats and noodles. Any session that finds an upstream defect,
or needs to know whether one is known, sits here for the duration.

This file describes the shape of the work. The work itself lives in the per-upstream ledger,
`docs/polars-bio-pending-fixes.md` today; a second upstream gets its own `docs/<upstream>-pending-fixes.md`.

## What the seat answers

1. **Is this ours or theirs?** A symptom seen through polars-bio can live in our loader, in polars-bio,
   in datafusion-bio-formats, in noodles, or in the input file. Only the last four belong in the ledger,
   and the input-file case is not a defect at all.
2. **Is it known?** Reported, fixed, intended, or declared a non-goal upstream.
3. **What exactly is wrong, and where?** A repro, a `repo/path:line`, and the spec line it breaks.
4. **Should it be filed, where, and has anyone answered?**

## Inputs

- Our own findings, from tests and dogfooding (`docs/dogfooding.md`).
- The consumer siblings that use polars-bio (`../just-prs`, `../just-dna-lite`, `../just-dna-format`).
  Their workarounds are unfiled bug reports. Harvest them read-only; never edit a sibling.
- An upstream release. Every entry is re-checked when the pinned polars-bio version moves.

## The knowledge base

A local, greppable mirror under `data/interim/` (git-ignored). The ledger's header lists what is there,
the pinned revisions and how to refresh it. The rules for using it:

- **Mirror the pinned revision, not `main`.** Check the installed wheel against the clone (`diff` the
  Python package) and check out the tag or rev the lockfiles pin, all the way down the stack. A bug
  fixed on `main` is still present in what we run.
- **Tracker dumps hold titles and bodies only.** A maintainer's "by design" or "won't fix" is in the
  comments. Read them (`gh api repos/<owner>/<repo>/issues/<n>/comments`) before quoting a stance.
- **Search every layer's tracker, not only the one the symptom surfaced in.**
- **Docs include `openspec/` designs and their Non-Goals, the CHANGELOG and the docstrings.** A
  docstring is a claim to test, not a fact: `register_vcf` documents 1-based and returns 0-based (P6b).
- **Our own ledger is a claim to test too.** P2 said the spec allowed a late `##fileformat`; the spec
  says the opposite. Quote the spec line, never paraphrase it from memory.

## Per-item procedure

1. **Dedupe** against the ledger. A duplicate adds evidence to the existing entry; it does not get a
   new id.
2. **Reproduce** on the pinned version with the smallest synthetic input. Record the command and the
   output. No real genome in a repro, a note or an issue.
3. **Search** the trackers, docs and non-goals of every layer.
4. **Locate** the root cause as `repo/path:line` in the layer that owns it. The layer decides where to
   file: the one that makes the wrong decision, not the one that prints the error. (P9: noodles returns
   `None` by design; datafusion-bio-formats turns that into a failure, so it is filed there.)
5. **Check the spec** for anything a format rule settles, and cite the line.
6. **Decide the verdict**, one of:
   - *real bug* — file at the owning layer;
   - *intended, undocumented* — file a docs request, or nothing;
   - *spec-correct* — the input is invalid; repairing it is our job; do not file;
   - *fixed upstream* — name the release; the entry closes when our pin reaches it;
   - *cannot reproduce* — say what was tried.
7. **Write the entry** (template below) and update the verdicts table at the top of the ledger.

## Filing

Filing publishes, so it is the user's call, per batch. Draft first, then ask.

- One issue per defect per layer. Two defects in one function may share an issue (P6b and P7).
- Reformulate; do not forward. State it as a fact about the thing they own: symptom, minimal repro,
  `path:line`, spec line, the smallest fix we can see, and the related ids found in the search.
- Say we are relaying where a sibling found it first, and name the motivating case.
- Record the issue URL in the entry the moment it exists, and add a line to `docs/CHANGELOG.md` when
  filing changes how we integrate.

## Reply tracking

Each filed entry carries `Filed: <url> (<date>)` and `Last checked: <date>, <state>`. On each pass:

- read new comments; a maintainer stance goes into the entry verbatim with its date;
- a merged fix gets the release it lands in; when our pin reaches that release, re-run the repro,
  remove our mitigation in the same commit as the pin bump, and move the entry to
  `docs/previous_issues.md` with the resolution;
- a *won't fix* keeps the mitigation and turns the entry into a permanent note.

## The entry template

```
## Pn — <symptom, in one line>

Seen with <package> <version> (<date>) on <kind of input>. <What happens.>

Mitigation: <what our code does, path>, and what it tells the user.

*Upstream (<date>).* <Trackers searched and hits, with stances.> <Root cause: repo/path:line.>
<Spec: file:line.> **<Verdict>. File at <repo>** / **Do not file.**
Filed: <url> (<date>) · Last checked: <date>, <state>
```

## Hard rules

- Never edit or commit an upstream or sibling repo. Reading and cloning into `data/interim/` is fine.
- Never file, comment or react upstream without the user's explicit go for that batch.
- Never mark an item "not a bug" without the spec line or the maintainer comment that says so.
- Never drop a mitigation because an issue closed; drop it when the pinned release passes the repro.
- Delegated harvesting (subagents across siblings) returns findings to this seat. The seat writes the
  ledger; a harvester never edits it, so there is one writer and no merge.
