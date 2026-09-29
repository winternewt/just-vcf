# Agent Guidelines — just-vcf

`just-vcf` is a published Python library (distribution `just-vcf`, import `just_vcf`) that converts
consumer and personal genomics files into VCF on top of polars-bio. That covers 23andMe-style array
exports, other freeform genotype tables, and VCF/gVCF normalization. It is not a variant caller, not an
annotation pipeline and not an app. Its consumers are `just-prs` and `just-dna-lite`, plus anyone using
the package. The core carries a lean dependency tier (§5).

`AGENTS.md` is a symlink to this file. If the two ever differ, that is a bug —
`ln -sf CLAUDE.md AGENTS.md`.

---

## Read these first, in this order

Obligatory. Read them yourself. **Do not delegate a document you are about to judge a design
against** — a subagent returns a summary, and a summary of a rule drops the qualifier the
decision turned on. Delegation is for finding things, never for deciding them.

1. **docs/CHARTER.md** — the durable charter. What this project is, what it will never do, the
   invariants every release upholds. When a plan conflicts with it, the charter wins. An audit that has
   not read it is incomplete. It is self-contained and names no other document; navigation lives here.
2. **docs/ROADMAP.md** — active-only, forward-only. One `## RMn — name` section per *open* item with a
   severity/status/owner line. Claim new ids with `scripts/item-next.py`, never by grepping (§8).
3. **docs/CHANGELOG.md** — what actually shipped, newest first.
4. **docs/dogfooding.md** — open findings from using the shipped surface for real work. Read it before
   touching the tool or API surface.
5. **docs/FEEDBACK.md** — the consumer inbox. Empty means nothing is owed.
   `docs/FEEDBACK_HISTORY.md` holds the answered items; `docs/CONSUMER_TRIAGE_LOOP.md` is the runbook.
6. **docs/PROVENANCE.md** — where each borrowed piece of code came from in `just-prs` / `just-dna-lite`,
   and what was cut or changed on the way in. Read it before porting anything else.
7. **docs/SIBLING_SURVEY.md** and **docs/LANDSCAPE.md** — what just-prs/just-dna-lite already implement,
   and which external tools to take rather than build. Read both before designing a feature.
8. **docs/FORMATS.md** — how `sniff` recognises a file, how a format is declared, and why the
   registry order is what it is. Read it before adding or reordering a format.
9. **docs/PGP_HMS.md** — the public test-data source: where the files are, their quirks, what is on
   disk, and which participants make good multi-facet test setups.

Everything below is self-contained: no rule here requires following a link to know what you
must not do. Links carry positive detail only.

There are no `.claude/` skills, agents or workflows in this repo yet.

---

## 1. Adopting these guidelines: ask, never infer

**When two rules conflict — this file against a sibling repo's, this file against the user's
global preferences, a rule against what the code actually does — stop and run a
questionnaire. Do not pick the one that looks better, do not synthesize a compromise, and do
not silently follow the more specific file.** A contradiction between two live rules is
almost always a real difference in the repos' natures (a published library's dependency tier
against an app's, a wire contract against internal state), and inferring which nature applies
here is exactly the guess that produces a rule nobody agreed to.

How to run one:

1. **Survey first, ask second.** Read the conflicting rules in full and find out *why* each
   side adopted its version. A question that does not carry the reason is unanswerable.
2. **One question per contradiction, batched** — do not drip-feed. Each gets two to four
   concrete options, never an open prompt.
3. **Each option states its cost**, not just its label: what breaks, what it forces on other
   repos, which existing rule it contradicts.
4. **Recommend one and say so.** A questionnaire with no recommendation offloads work the
   survey already did.
5. **Record the answer where it will be read again** — the rule into the relevant section
   below, the reasoning into *Learned user preferences* in the user's own words. An answered
   contradiction that is not written down gets re-asked, which is worse than having guessed.

The one exception: a conflict with the charter is not a questionnaire. The charter wins.

---

## 2. Non-negotiables

Read the whole list before the first edit. The reason follows each one, because a rule
without its reason gets rationalised away at 2 a.m.

- **Never `uv pip install`.** Use `uv sync` / `uv add` / `uv add --dev`. `uv pip install`
  writes into the venv without touching `pyproject.toml` or the lockfile, so the next clean
  checkout silently lacks the dependency.
- **Never call bare `python` / `python3`.** Always `uv run python …`, `uv run pytest …`,
  `uv run <script>`. A bare interpreter bypasses the project environment. The stdlib-only triage
  scripts in `scripts/` are the one exception; they run directly.
- **Never hardcode a version string.** It comes from `pyproject.toml`, via
  `importlib.metadata.version("just-vcf")`. Two sources of truth drift, and the one you
  read is the wrong one.
- **Never rename a user-facing command to dodge a stale `uv run` wrapper.** Bump this package's
  version and re-run `uv sync` so uv rebuilds the entry points.
- **Never use a placeholder path or a fabricated example value** in committed code —
  `/my/custom/path/`, a dummy digest, `rs999999999`, `1e-328`. A fabricated value proves
  nothing and outlives the session that invented it.
- **Never mock the transformation under test.** See §6.
- **Never claim a test "would have caught" a bug without first running it against the buggy
  code and watching it fail.**
- **Never commit large data, and never commit a real person's genome in any amount.** No
  VCF/parquet/gz/BAM/FASTA/`.db`. Anything over ~5 MB that genuinely must travel goes through Git LFS.
  Everything under `/data/genomes` belongs to a named person: not even a 1000-line head of it goes into
  `assets/`. Tests point at `JUST_VCF_SAMPLES_DIR` and skip when it is absent. Consented-public data
  (PGP-HMS) is the fixture source once the leecher exists.
- **Git: commit as you go, never push.** Land each coherent block as its own commit, semantically
  grouped: a behaviour plus its test, doc and changelog line in one commit. No conventional-commit
  prefixes. Tagging and local branches are allowed. **Never** `git push`, publish to PyPI, or cut a
  release unless explicitly told. Never delete a branch without asking. No `--amend`, rebase,
  `reset --hard` or force-anything; a mistake is fixed forward. **Never `git stash drop` /
  `git stash clear`**, even on explicit request. **Never `git add -A` or `git add .`** — it sweeps in
  `.env` files, editor swap files and junk; stage explicit paths and read `git status` first.
- **Never silently fall back when primary data is missing.** If the reference FASTA is absent, do not
  quietly skip REF lookup and emit something that looks resolved. Refuse with an explicit error or
  mark every affected row unresolved. The user cannot see that the data source differs.
- **Never nest a `try`/`except` inside another.** It hides the real error. Let typed
  exceptions propagate; wrap only where there is a genuine recovery path.
- **Never fill a value from the same source that checks it.** A cross-check compares an
  independently authored value against a source; filling it *from* that source makes the
  check compare a convention against itself, and it agrees perfectly. Worse, the row moves
  from honestly unverified to apparently verified.
- **Never let a tool write a checked value from a lookup result.** Lookups report, the human
  decides, the linter checks. Where a tool refuses to apply a value, preserve the refusal
  verbatim — the refusal is the feature.
- **Never treat a determinism gate as a correctness gate.** A digest match or a reproducible
  output means *reproducible*, not *right*.
- **Never collapse "unknown" into a boolean, and never collapse distinct failures into
  `False`.** `None` is unknown / did not run; a timeout is not a negative result. See §5.
- **Never reduce a genotype to a binary called / not-called.** The charter makes genotype evidence
  graded. Whatever the representation turns out to be, a boolean "is called" flag is not it.
- **Never write REF from the genotype.** REF comes from a reference (FASTA or a reference panel) or is
  unresolved. just-prs `arrays.normalize_array` writes the observed allele as REF; do not port that.
- **Never route around a capability the product lacks while dogfooding.** See §7.
- **Never resolve a contradiction between two rules by inference.** Run §1.
- **Never edit or re-wrap a reporter's prose** in the feedback inbox — not when answering it,
  not when moving it. Replies are appended below it.
- **Never reuse a feedback id**, not even one answered as a non-issue. Compute the next id
  from the inbox **and** the history file: `./scripts/triage-state.py --next`.
- **Never treat "no reply" as "no work done."** Establish what already shipped before
  reproducing a reported item, and certainly before designing for it.
- **Never write a placeholder reply.** Leave an untriaged item unanswered and say what was skipped.
- **Never open a preamble line with `**Status`** in an inbox. It is read as a block reply
  covering every id it names, which marks the whole backlog answered.
- **Never file into a producer's roadmap.** File into its `docs/FEEDBACK.md`. See §8.
- **Never edit or commit a sibling repo** (`../just-prs`, `../just-dna-lite`, …) beyond appending an
  item to its `docs/FEEDBACK.md`. Borrowing means copying code *in*, recorded in `docs/PROVENANCE.md`.

---

## 3. Repository layout, data and assets

```
src/just_vcf/         source (src layout; hatchling packages = ["src/just_vcf"])
tests/                pytest suite
docs/                 all markdown except this file and README.md
assets/               fixtures that MUST travel with the repo — committed, synthetic or consented-public only
data/input/           input samples          ─┐
data/interim/         generated intermediates ├─ git-ignored, never travels
data/output/          results                ─┘
scripts/              triage-loop tools and operational one-offs (the PGP-HMS leecher), not importable
```

- `data/` is git-ignored by ignore-all + allowlist (`data/.gitignore`: `*`, `!.gitignore`). To commit a
  subtree, add explicit `!<dir>/` and `!<dir>/**` lines. The root `.gitignore` also ignores every
  genomic file type everywhere.
- Committed test data lives in `assets/`, not `data/`. If a test writes, it writes under `tmp_path` or
  a resolved cache dir (`platformdirs`), never into the project tree.
- Never hardcode a platform-specific cache path; resolve it at runtime.
- **LFS history gotcha:** a blob committed *before* `git lfs track` stays in every past
  commit even after the pointer replaces it at HEAD. Detect it, then hand the remediation to the user:

  ```bash
  git lfs ls-files
  git rev-list --objects --all \
    | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' \
    | awk '/^blob/ && $3 > 5000000 {print $3, $4}' | sort -rn
  ```

Sibling repos, read-only unless a task explicitly targets them:

| Path | What we take from it |
|---|---|
| `../just-prs` (MIT) | array loading, VCF build detection, REF lookup (`reference_allele.py`), liftover, ancestry |
| `../just-dna-lite` (AGPL) | `annotation/restoration.py` only: callset-scope (WGS vs targeted) detection, gVCF detection, hom-ref restoration. The user wrote all of it (git blame, 2026-09-28) and relicenses it to MIT here. Everything else there is AGPL and mostly antonkulaga's (`chromosomes.py`, `io.py`, most of `QualityFilters`): **never copy it**. Reimplement from the behaviour, or ask the user to get permission first |
| `../just-dna-format` | producer runs the same triage loop; a model for inbox and charter shape |

---

## 4. Build, run, test

```bash
uv sync                      # core
uv sync --all-extras         # + reference (pysam), liftover (pyliftover), ancestry (numpy)
uv run pytest -vvv           # the suite; -vvv when diagnosing
```

There is no CLI entry point yet. When one lands it is Typer, registered in `[project.scripts]`, and
it loads `.env` via `python-dotenv` at startup before reading configuration. New configurable values
are read from env vars with sensible defaults, documented in `.env.template`, and listed here:

| Variable | Default | Meaning |
|---|---|---|
| `JUST_VCF_SAMPLES_DIR` | unset | directory of real genomes for opt-in tests and dogfooding; tests skip without it |

**Timestamps: store ISO-8601 UTC, display local.** Never write a naive `YYYY-MM-DD HH:MM:SS`.

---

## 5. Coding standards

- **Type hints mandatory** (Python ≥3.10 syntax). **`pathlib.Path`** for every path.
- **Absolute imports only**, and **no inline imports** — every import at module top level.
  The sole exception is a guarded module-level `try/except ImportError` for an *optional*
  dependency, raising a clear message with install instructions. The core must import
  cleanly with only its hard dependencies present.
- **Dependency tier.** The core may import polars / polars-bio, gtars, refget, pydantic, platformdirs,
  typer and python-dotenv, plus the stdlib. gtars (Rust, BSD-2) mints VRS ids and normalizes alleles;
  refget computes and compares seqcol assembly digests. vrs-python is a test oracle only, never core. Everything else goes behind an extra with a guarded import:
  `[reference]` = pysam (not on Windows), `[liftover]` = pyliftover, `[ancestry]` = numpy. httpx is a
  dev dependency for `scripts/` only and is never imported from `src/`. Adding a hard core dependency
  is a questionnaire (§1), not an edit. duckdb, eliot, huggingface-hub, Dagster and Reflex from the
  sibling repos do not come along.
- **Typer for every CLI.** **Pydantic 2** at every boundary — config, external input, public result
  models. Internal structures may use a dataclass where that is simpler; say so here when a framework
  forces it.
- **Standard-library `logging` for diagnostics — never `print`.** `print` is only for CLI output the
  user asked to see. Code ported from the sibling repos uses eliot `start_action`; convert it to
  `logging` on the way in.
- **Polars over Pandas.** Prefer expressions over Python loops; lazyframes (`scan_*`) and streaming
  (`sink_*`) on large paths; pre-filter before joining. Avoid a full `collect()` where a lazy
  aggregation or `pl.collect_all()` streams the result. A WGS VCF is millions of rows and a gVCF is
  hundreds of millions.
- **`__all__` is curated.** This is a published library with an API contract, so `just_vcf/__init__.py`
  curates `__all__` and keeps it in sync with its imports. Anything not in `__all__` is internal.
- **Constrained vocabularies are criterion-based.** A vocabulary that must grow additively inside a
  major version (input formats, vendors, chip versions) is `frozenset[str]` + a validator. A vocabulary
  bound to a fixed schema/dtype is a `str`-subclassing enum. Never a bare `Literal` in a wire model.
- **Kleene by default for any result carrying uncertainty. `True` = positive, `False` = negative,
  `None` = unknown / did not run.** When the answer is unknown, **withhold** — never report, never
  negate. `None` is never `False`. An absent measurement gets an explicit sentinel, never the lowest
  bucket. **A check that could not run is not a check that passed.** Use Kleene semantics:
  `unknown AND false` is `false`.
- **When "negative" has kinds, a bool is the wrong return type.** Return a small result object with
  `__bool__` that carries the discriminant (an outcome enum, a reason, the underlying error). Build
  detection, format detection and callset-type inference are all of this kind: "not GRCh38", "could not
  tell" and "file unreadable" are different answers.
- **Deterministic ordering is load-bearing** wherever bytes are hashed, diffed or shipped — and a VCF
  is shipped. Never derive emitted rows from `set`/`dict` iteration or from polars `mode()`/`unique()`
  without an explicit stable sort or tie-break. Every new ordering gets a test.
- **Aggregate repeated warnings.** Group by *reason*, with a count, never one line per row.
  Distinguish "the source did not say" from "the source said something we cannot represent".
- **Derive a guard from the model's rule, not beside it.** Test the guard against the model case by
  case rather than asserting a message.
- **Heed terminal warnings, deprecations especially.** A deprecation in code you touched is a
  **blocker**: find the current upstream API, fix it, and update this file.
- **Refactor internals aggressively** — no dead code, no old API kept for nostalgia. The contract —
  public API, the VCF and parquet shapes we emit, the CLI surface — breaks only deliberately and
  versioned.

### Genomics rules

These were learned the expensive way in the sibling repos. Each is a prohibition, so each is here.

- **Coordinates are 1-based VCF `POS` at every public boundary.** polars-bio returns 0-based `start`
  from `register_vcf`-named tables, and just-prs adds 1 to compensate. Never mix the two in one frame;
  convert at the edge and name the column for the convention it holds.
- **Never run concurrent `pb.scan_vcf` on the same table name.** It panics with "table already exists".
  just-prs registers under a content-hash name behind a lock (`vcf.scan_vcf`).
- **polars-bio joins multi-allelic ALT with `|`** and decodes FILTER `.` as `""`. Handle both
  explicitly; a naive equality against an allele silently misses every multi-allelic site.
- **Never strip `chr` with a plain `^chr` regex.** `chrM` becomes `M` and mtDNA silently disappears
  from any join against `MT`. Map `chrM`/`M`/`MT` explicitly.
- **Detect compression by magic bytes, never by extension.** Real inputs include BGZF files named
  `.vcf`.
- **A chip's absent site is not hom-ref.** An array export drops no-calls, so a position missing from it
  is unknown, not reference. Only a whole-genome callset supports "absent means reference", and even
  there it is graded evidence (see the charter), not a certainty.
- **Genome build is detected, never assumed.** When detection cannot tell, the answer is unknown and
  says so. It never silently defaults to GRCh37.

---

## 6. Testing — layer 1

- **Real data + ground truth.** Exercise the actual code path against real fixtures and **compute
  expected values at runtime** from the fixture. Do not mock a data transformation.
- **Fixtures:** synthetic inline strings for parser shape, consented-public samples (PGP-HMS) in
  `assets/` for real-file behaviour, and `JUST_VCF_SAMPLES_DIR` for opt-in tests against local genomes.
  The last kind skips when the variable is unset or the file is missing.
- **Meaningful assertions.** Prefer relationships, aggregates and set equality over existence
  or count checks. `assert set(source_ids) == set(output_ids)`, not `assert len(df) > 0`.
- **What to validate:** counts and aggregates; joins (pre/post counts, key coverage, nulls introduced by
  outer joins); transformations (round-trip survival, key preservation); data quality (ranges,
  malformed entries, duplicates).
- **Hardcoding domain constants** — chromosome lengths, a documented threshold — is fine. **Hardcoding a
  row or unique count read off a data dump is not**; it drifts.
- **Deterministic coverage.** Fixed seeds and explicit filters; representative *and* edge cases.
- **Be resilient to changing data.** `pytest.skip()` for optional or live data; skip gracefully when
  offline. `pytest.mark.xfail` with a clear reason for a known upstream defect.
- **Avoid the AI test anti-patterns**: happy-path-only; counts derived from inspecting data;
  mocking the transformation under test; redundant checks; ignoring nulls, empties, boundaries and
  malformed input.
- **A test that means "no credential" must say so** — neutralize in an autouse fixture with
  `setenv(VAR, "")`, not `delenv`.
- **Suspect ordering whenever a test passes alone and fails in the suite.**

---

## 7. Dogfooding — layer 2

Tests prove the code does what it was told. Dogfooding asks: *is this usable, and what is missing?*
Here that means converting real genomes end to end: the `/data/genomes` samples and PGP-HMS downloads.

**Do not "verify the tool's answers" with a second independent implementation while dogfooding — that
is a test, and it belongs in the suite.** Use the tool, notice the friction, write down what was not
there.

- **A capability the tool LACKS is the result, not an obstacle to route around.** The moment you reach
  for an ad-hoc script, `bcftools` or a hand-written query to get past something the library cannot
  do, the exercise has stopped producing signal. Record the gap; if it blocks the work, build it into
  the library and carry on **with the library**.
- **Run the adversarial round.** Be a beta-tester trying to show the thing fails at something it
  advertises, then switch back and fix. **Attack claims, not gaps**, and **use real data**.
- **Pick the probe where the design generalized from one case.** One vendor's format → try a second.
  One build → try the other. Two alleles → try a multi-allelic site.
- **Turn the tool on the work you just did.**
- **Dogfood a finding before you report it.** A demonstrated failure beats a plausible-looking
  mechanistic claim.
- **Finish each probe as a committed reference example whose README names what it broke.**
- **Separate "fix it" from "surface it" before writing any code**, and say *why each candidate repair is
  wrong* for the ones you surface.

### The finding log

Findings carry stable `F#` IDs and **move** between files — never duplicated, except where a finding
is mitigated here but still owed upstream, which legitimately appears in two.

| File | Holds |
|---|---|
| `docs/dogfooding.md` | open quirks, bugs and UX gaps found by using the shipped surface |
| `docs/previous_issues.md` | findings resolved **here**, each with its resolution and a code pointer |
| `docs/<upstream>-pending-fixes.md` | findings needing an **upstream** change first (polars-bio is the likely one), each noting the mitigation already in place |

---

## 8. Docs and their lifecycle

- **All new markdown goes in `docs/`** — the only exceptions are this file and `README.md`.
- **`docs/` is the single ground truth.** This file duplicates from it only where a fact is needed to
  *orient* — and every prohibition lives here in full.
- **`docs/ROADMAP.md` is active-only.** Shipped items move to `docs/ROADMAP_HISTORY.md` with their
  rationale. Nothing is deleted; it is relocated. Claim a new `RMn` with
  `INDEX=docs/ROADMAP.md ITEM_PREFIX=RM ./scripts/item-next.py --note "…"`, never by grepping the
  highest id — two sessions reading the same number both write it.
- **`docs/CHANGELOG.md` records what shipped**, newest first — including cross-repo integration
  changes, so agents in sibling repos are not surprised.
- **`docs/PROVENANCE.md` records every borrowed piece**: source repo, path and commit, what was cut,
  what was changed. Update it in the same commit as the port.
- **Consumer → producer intake.** When something belongs in a repo we consume but do not own
  (polars-bio, or a sibling library), **do not edit or commit that repo** beyond appending your item.
  File into its **`docs/FEEDBACK.md`** where it has one; polars-bio has no local clone and no inbox, so
  use its GitHub issue tracker and record the reference here. Record the integration change in our
  `docs/CHANGELOG.md`.
  - **Reformulate; do not forward.** Restate it as a fact about the thing they own, say you are
    relaying, and name the motivating case.
  - **Ids are per file.** Take the number from *their* ledger.
  - **Our own record still has to stand alone.** Name the id we filed and what we asked for.
  - `docs/<upstream>-pending-fixes.md` is where such an item waits before it is filed and after.

### The consumer inbox

`docs/FEEDBACK.md` holds **open** consumer items (`## Sn — …`); `docs/FEEDBACK_HISTORY.md` holds
answered ones verbatim. **Empty inbox ⟹ nothing owed.** The runbook is
`docs/CONSUMER_TRIAGE_LOOP.md`; the tools are in `scripts/` and run directly (stdlib only):

```bash
./scripts/triage-state.py --pending   # what needs work
./scripts/triage-state.py --next      # next unclaimed id, over BOTH files
./scripts/triage-archive.py S1        # move an answered item, verifying the move
./scripts/watch-inbox.sh              # one line when the inbox settles
```

The three `.py` scripts are Python. Run them directly or with `python3`, **never with `bash`**: under
bash the docstring executes as commands and `import` reaches ImageMagick.

Triage order: establish what already shipped → reproduce against the **code**, not the docs →
decide **legality before severity** → route → reply → archive.

- **Legality sizes the release; severity only orders the queue inside it.** A new optional field or an
  added flag is minor; a removal, a promotion to required, a retype — **including a rename** — is major.
  Pure legibility is a patch. When changing a shape a consumer reads, **add rather than redefine**.
  Pre-1.0 we are allowed to break, but the changelog says so every time.
- **A non-issue verdict is not the cheap outcome.** "Nothing is wrong here" has to be shown.

- **Update this file and the affected `docs/` in the same change as the refactor**, not after.
  Policy is written first; code complies.
- **Run the commands yourself** rather than telling the user to run them — except where a command
  genuinely needs an interactive terminal.

### Prose style

Natural, human prose. Avoid AI tells — em-dash pile-ups, filler transitions, marketing voice.
Never hallucinate documentation or overpromise an unimplemented feature. **just-vcf must never be
described as calling variants, as imputing genotypes, or as producing clinical-grade output.** It
converts and normalizes what another tool already called, and says how sure it is.

---

## 9. Self-correction

When outdated API knowledge causes a real crash or logic failure, fix the code **and** update this
file (and the affected `docs/`) with the correct pattern, so the next agent does not repeat it. The
same applies when the user corrects a preference: it goes below, in their words, with the reason.

## 10. Learned user preferences

- 2026-09-28, license: MIT. "the most valuable wgs restoration code in -lite built by me anyway,
  license stays, but I can reuse it" — just-dna-lite stays AGPL; the parts the user wrote may be
  ported here under MIT.
- 2026-09-28, dependencies: lean core plus extras, "but drop leech, just dev" — the PGP-HMS leecher's
  httpx is a dev dependency, not an extra.
- 2026-09-28, git: "Autocommit-autopush as per current global memory." The global rule on record says
  commit as you go and never push. Checked 2026-09-28: no sibling repo overrides that (just-dna-format's
  triage permit keeps push out too). The user said "ok keep never push", so §2 stands.
- 2026-09-28, charter: yes, "but no heavy stance on no-calls. … a bayesian based fuzzy 0(unknown) ..
  PHRED-like(confident genotype). … each position is probably an n-vector of probabilities or log or
  something, tbd; No binary call-no call definitely."

- 2026-09-28, scope: "never reads meant 'we don't implement it', if it is here for free, let it be".
  An out-of-scope capability means we don't build it. When the substrate (polars-bio) already has it,
  use it. Example: polars-bio `depth()` over BAM for coverage.

- 2026-09-28, VRS: "this one is not a lightweight lib, so gtars it is". Core carries gtars + refget for
  VRS ids and seqcol, over the earlier pick of just-dna-format (SNV-only, GRCh38-only). A heavy
  dependency that covers GRCh37, T2T and indels beats a light one that doesn't.

## 11. Learned workspace facts

- Real genomes live under `/data/genomes/<person>/`, all private. What is there by kind:
  - `angellos/2024_genotec.txt` — a 23andMe-format array export, the only consumer array file locally.
  - `advocate/deepvariant/`, `john/deepvariant/` — DeepVariant `.g.vcf.gz` + `.vcf.gz` pairs with `.tbi`.
  - `anton/variants/` — GRCh37 CrossMap/Picard liftover outputs, with reject files.
  - `olga/BorysovaGenome.liftover38to37.zip` — a zipped VCF; `ksuhaster/deepvariant/*.zip` likewise.
  - `sergey/` — a DRAGEN `hard-filtered.vcf.gz` and a GATK HaplotypeCaller VCF.
  - `gennet/` — HapMap CEU lifted VCFs. `prive/`, `aadr/`, `dataverse_66p1/` are ancestry reference data.
- `/data/genomes/pgp-hms/participants/<hu>/` holds PGP-HMS data per person (`scripts/pgp_participants.py`):
  `arvados-<collection>/` for PGP's own sequencing (freeBayes+VarScan2 hg19 VCFs with **two sample columns**,
  one per caller; Complete Genomics `var`/`masterVarBeta`/`vcfBeta`) and `upload-<download id>/` for the
  person's own uploads. `manifest.jsonl` there records every file and every folded duplicate.
- `/data/genomes/pgp-hms/` holds PGP-HMS downloads from `scripts/pgp_leech.py`. Those are consented-public,
  unlike the rest of `/data/genomes`, and CC0 per the PGP Researcher FAQ, so a small one may become an
  `assets/` fixture. Data types are self-reported, so a file labelled 23andMe may be something else. Check the file first:
  participant uploads can carry real names in file names.
- `../just-prs/data/input/vcf/newton_winter.vcf` is BGZF despite the `.vcf` name.
- Public WGS samples used by the sibling suites are on Zenodo, records 18370498 and 19487816.
- PGP-HMS publishes consented 23andMe exports at
  `https://my.pgp-hms.org/public_genetic_data?data_type=23andMe`; the leecher (roadmap) fetches them.
