# Trust register

just-vcf's tests and dogfooding read `JUST_VCF_SAMPLES_DIR`, which on the person's machine is the
whole genomes folder: the consented-public PGP-HMS downloads, and beside them the genomes of named
people, all private (`CLAUDE.md` §11). Code that runs here can read every one of them, so this
repository keeps a register (the seats kit's base, "Code from elsewhere, where a repo holds private
things").

**The rule:** before code from outside this repository's own history runs here in any way (tests, an
install, a build, a script, a reproduction), check its source's standing below.

- **trusted**: it runs.
- **audit**: the whole incoming code is read first. Every file outside what it claims to do is
  named, and so is every place that runs by itself: build scripts and hooks, test fixtures, CI, agent
  settings and instructions, and every network call, subprocess or file write outside the repository
  and the cache. The person approves before anything executes, unless the source's row says otherwise.
- **untrusted**: read, never run. A fix from it is rewritten from the reading.

A source that is not listed is untrusted. Only the person sets a standing; a session never raises
one itself.

**Outside code run on its own runs against public genomes only**: a consumer's reproduction, an
upstream clone's tests, a build of the stack from source. It runs with `JUST_VCF_SAMPLES_DIR` pointed
at the PGP-HMS folder (consented-public, CC0), never the whole genomes folder, unless the person says
otherwise for that run. A private genome is never needed to show an upstream or consumer defect: a
repro is synthetic (the liaison's rule) or PGP-HMS. This repository's own tests and dogfooding, with
the packages `uv.lock` pins, keep reading the whole folder.

## Sources

| Source | Standing | Why | Decided |
|---|---|---|---|
| polars-bio and what its wheel builds in (datafusion-bio-formats, the noodles fork), at the release `uv.lock` pins from the package index | **the pin, with one audit owed** (proposed) | Pinned in the lock, so outside the rule like any pinned package: the wheel keeps running in the tests while the audit is owed. The person, 2026-10-18: "maybe one-time audit of upstream is due". One audit of the source the pinned release was built from (the mirror under `data/interim/`, at the pinned tags): what runs at import, every network call, subprocess and file write outside the paths it is given. The upstream-liaison seat runs it as its first item, its read-in-full pass done by an Opus subagent; the seat's verdict, recorded under "Audits", stands once the person confirms it. A pin move is read by the liaison's own re-check of the ledger, not audited again. | the person's merge of the first adoption pull request |
| The same stack built from source: a clone under `data/interim/`, a branch, a pull request, a git dependency | **audit** (proposed) | Code no release vouches for, run to test a fix or a mitigation: read, and the person approves, each time. | the person's merge of the first adoption pull request |
| A reproduction in a consumer's inbox item (`docs/FEEDBACK.md`) | **audit** (proposed) | Written by another repository's session: read whole, and the person approves, before it runs, against public genomes only. A repro rewritten by the seat from the reading is this repository's own code. | the person's merge of the first adoption pull request |
| The person's sibling repositories (`../just-prs`, `../just-dna-lite`, `../just-dna-format`), at their own committed history | **trusted** (proposed) | The person's own code, already run on this machine. Their uncommitted working trees are not covered. | the person's merge of the first adoption pull request |
| The seats kit, at the tag this clone pins | **trusted** (proposed) | The person's own repository, pinned by tag. Its hooks and watchdog run in every session here. Each new release adopted on its announcement. | the person's merge of the first adoption pull request |

## Outside the rule

- **This repository's own history.**
- **Packages pinned in `uv.lock` from the package index**: the pin is their standing (the row above
  adds one audit for the polars-bio stack). A new or changed dependency is still read before it lands.
- **The person's user-level Claude Code setup**: it loads in sessions here, but this repository does
  not pull it in.

## Not listed, so untrusted

Code from an upstream tracker's issue or comment (the tracker dumps under `data/interim/` included),
from another person's fork or pull request before upstream merges it, or pasted into a chat. Any file
another session leaves in the working tree uncommitted.

## Audits

None yet. The polars-bio stack's one audit is recorded here when it is done: the release, what was
read in full, what was not, the verdict and the person's confirmation. Reading the source does not
vouch for the compiled extension in the wheel from the package index: the record says so.
