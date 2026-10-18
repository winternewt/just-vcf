# Trust register

What code from outside this repository may do here, as a threat model per source rather than one rule
for all (the seats kit's base, "Code from elsewhere, where a repo holds private things"). The person,
2026-10-18: "there should be a threat model for each lib, this is mature security aproach, not this
random alarmism". Only the person sets a standing; a session never raises one itself. Each row below
is proposed and set by the person's merge of the first adoption pull request.

## What is worth protecting

- **The private genomes**: `/data/genomes/<person>/`, named people's, and their decompressed copies in
  the loader's cache (`JUST_VCF_CACHE_DIR`, by default the platform's user cache; RM268). The PGP-HMS downloads beside them are
  consented-public (CC0).
- Nothing else here is secret: no tokens or keys in `.env`; the forge and messenger credentials belong
  to the machine, not this repository.

**The one threat is exfiltration**: a private genome, or anything read from one, leaving the
machine. The person: "vcf lib fears data exfiltration, that's basically it. Local perimiter is
considered safe and secure". Code reading those files locally is the library's job, and the tests and
dogfooding read them on purpose; what any code does inside the machine is out of this model.

## The setting

A private workstation, one person, a trusted local perimeter: nothing here is walled off from
anything else, and nothing needs to be. The controls are on what leaves it: the client's own
exfiltration guard on agent actions, and the audits below for code that could send data out without
an agent asking.

## Per source

| Source | What it runs here | Realistic threat | Control | Residual, accepted |
|---|---|---|---|---|
| **polars-bio**, with what it builds in (datafusion-bio-formats, the noodles fork); from the package index, pinned in `uv.lock` | A compiled extension and its Python layer, in every test and dogfooding run; reads every file it is given | Code that sends data out: telemetry, an upload, a remote read reaching beyond the URL it was given. It reads remote inputs by URL by design, so network code exists; the question is only where it can be pointed. Second: a compromised release at an upgrade | **One audit**, the upstream-liaison seat's first item, of the source the pinned release was built from: every network path and what drives it, import-time code, build scripts. A clean verdict is stamped "no network exfiltration code" and the source is trusted at any release, from a wheel or built from source. Upgrades are deliberate pin moves, each read by the liaison's re-check | The compiled wheel is not the source that was read; a later release is not re-audited |
| **Other packages** pinned in `uv.lock` | Library code in the same runs | A compromised release at an upgrade | The lock: an upgrade is a deliberate change, and a new dependency is a decision (`CLAUDE.md` §5) | Accepted as for any pinned package |
| **The person's sibling repositories** (`../just-prs`, `../just-dna-lite`, `../just-dna-format`) and the reproductions their sessions write into `docs/FEEDBACK.md` | Ported code; a repro run once by triage | None that sends data out: the person's own code, inside the perimeter | Trusted | Accepted |
| **The seats kit**, at the tag this clone pins | Hooks and a watchdog in every session; it talks to the messenger and the forge | The one real exit for genome data here is **an agent putting it into a message, an issue, a commit or a pull request**, which the kit carries out | Each seat's priming forbids genome data in any of those; the client's guard; the kit is the person's own, pinned by tag, each release adopted on its announcement | Accepted |
| **The model provider**, through the sessions themselves | Whatever an agent reads into its context is sent to the model | A genome's rows printed into a session | Work through the library's aggregates and counts, never by printing a private genome's rows into a session; PGP-HMS for anything shown | Counts and summaries do reach the provider |

## Not run

Code from an upstream tracker's issue or comment (the tracker dumps under `data/interim/` included),
from another person's fork or pull request before upstream merges it, or pasted into a chat: read,
and a fix rewritten from the reading. A source not listed here is in this group until the person adds
its row.

## Audits

None yet. Each audit is recorded here: the release read, what was read in full and what was not, every
network path found and what drives it, the verdict, and the stamp when it holds.
