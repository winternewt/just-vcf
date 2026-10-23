# Trust register

What code from outside this repository may do here, as a threat model per source rather than one rule
for all (the seats kit's base, "Code from elsewhere, where a repo holds private things"). The person,
2026-10-18: "there should be a threat model for each lib, this is mature security aproach, not this
random alarmism". Only the person sets a standing; a session never raises one itself. Each row below
is proposed and set by the person's merge of the first adoption pull request.

## What is worth protecting

- **The private genomes**: `/data/genomes/<person>/`, named people's, and their decompressed copies in
  the loader's cache (`JUST_VCF_CACHE_DIR`, by default the platform's user cache). `clear_cache(source)` removes a
  source's copies and every spill derived from it, by path or, once the file is gone, by sha256 (RM268). The PGP-HMS downloads beside them are
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

## The surface an agent reads through

An agent asking about a private genome uses `just_vcf.query` (RM283), not `load`. `query(path)` loads
the file itself and answers with counts by contig, variant class, filter set, genotype class and block
kind; summary statistics of GQ, DP and QUAL; the facet's provenance; the load's warnings grouped by
reason; and, for one variant asked by VRS id or by position and alleles, whether the facet carries it,
graded (`carried` true, false or unknown, and what the answer rests on). Every answer is a frozen
pydantic model of counts, numbers, the schema's enums and bounded strings. None has a field for a row,
a sample name, a path or an archive member's name, and a position, allele or genotype is returned only
when it was asked for; the sample, the paths and anything shaped like a path are replaced in the text
it returns. There is no position dimension and no region filter, so counts cannot be narrowed onto a
row. Small counts are not suppressed: the data is one person's, a small count names no one else, and a
suppressed cell is recovered by subtraction from its total. This is the model-provider row's control
made into an API, a convention for the session and not a sandbox: `load` stays for code that reads
rows on the machine, and its frames are never printed into a session.

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

Each audit is recorded here: the release read, what was read in full and what was not, every network
path found and what drives it, the verdict, and the stamp when it holds.

### polars-bio 0.36.0 and what it builds in (2026-10-18, upstream-liaison)

**Stamp: no network exfiltration code.** Nothing in the stack sends data anywhere on its own: no
telemetry, no update check, no crash reporter, no socket outside the HTTP client, no subprocess, and no
remote write, upload or delete reachable from polars-bio. Every network request starts from a path the
caller passed with a `scheme://` prefix. A local path never builds a network client and never loads a
credential. just-vcf hands polars-bio local paths only.

**What was read.** The source each piece of the pinned build came from, at the revision polars-bio's
`Cargo.lock` and `uv.lock` pin, all under `data/interim/` (the ledger's mirror table lists them). Five
read-only Opus subagents each took one slice. The seat re-checked every claim in the list below against
the code before recording it.

| Piece | Revision | How it was read |
|---|---|---|
| polars-bio Python package | `a1f724f` = 0.36.0; byte-identical to the installed wheel's `.py` files | all 24 modules in full, plus `pyproject.toml`, `Makefile`, `requirements.txt` |
| polars-bio Rust extension | same | all of `src/*.rs` in full, plus `Cargo.toml` and `rust-toolchain.toml` |
| datafusion-bio-formats | `v1.13.1` | `bio-format-core` and `bio-format-vcf` in full; the other 15 linked format crates by grep, reading in full any file a hit landed in; `bio-format-ensembl-cache` skipped (not linked) |
| noodles, biodatageeks fork | `42a3c01` | the fork's 21 commits over `noodles-0.113.0` (`fe2b112`) in full; the rest is zaeleus upstream, read by grep |
| bigtools, biodatageeks fork | `0d7a572` | the fork's 17 commits over upstream `b1fd283` in full; manifest and remote-file module in full |
| datafusion-bio-functions | `v0.22.2` | manifests in full; the three linked crates (ranges, pileup, fastqc) by grep, with the files a hit landed in read; `bio-function-vep` not linked |
| Apache OpenDAL | `1755db1` (main, 279 commits after `v0.55.0`, not a release) | the compiled-in services (s3, gcs, azblob, azure-common, http) and layers (logging, retry, timeout, concurrent-limit) in full; `opendal-core` by grep |

**Not read**: the crates.io dependencies, ridden on the lock pin as any package is (the row above).
Among them are the ones that matter most to this question: `reqsign-*` 3.0.1 (the credential chains),
`reqwest` 0.12.28 and 0.13.4, `object_store` 0.13.2, `ctor`, and the native build scripts
(`hdf5-metno-src`, `libsqlite3-sys`, `libdeflate-sys`, `aws-lc-sys` and others). Neither was the
compiled `.so` itself; the source is what was read.

**Every network path, and what drives it.** All of them sit behind
`datafusion-bio-formats/datafusion/bio-format-core/src/object_storage.rs:447-466`
(`get_storage_type`): a path is remote only when it carries a recognised `scheme://`, and an unknown
scheme falls back to local. `RemoteObject::open` returns an error for a local path (`:956-961`) before
any operator or credential chain is built.

- Reads from the caller's URL through an OpenDAL operator, read-only: S3 `:810-827`, Azure Blob
  `:859-872`, HTTP `:880-912`, GCS `:931-948`. The destination is the URL's bucket or host. Two
  environment variables can redirect it: `AWS_ENDPOINT_URL` (`:813`) for S3, which OpenDAL extends to
  `AWS_ENDPOINT` and `AWS_S3_ENDPOINT` (`opendal/core/services/s3/src/backend.rs:790-799`), and
  `AZURE_ENDPOINT_URL` (`object_storage.rs:566`).
- A companion index probe for a remote BCF or BGEN file, at the same host with `.csi` or `.bgi` appended
  (`bio-format-vcf/src/table_provider.rs:1025-1045`). A text VCF's index is looked up only for a local
  path (`:1089`).
- **One request to a host the caller did not name.** On an `s3://` path with neither `AWS_REGION` nor
  `AWS_DEFAULT_REGION` set, region discovery sends an unsigned HEAD to `https://s3.amazonaws.com/<bucket>`
  (`object_storage.rs:728-737`, through `opendal/core/services/s3/src/backend.rs:649-656`), even when
  `AWS_ENDPOINT_URL` points elsewhere. It carries the bucket name and nothing read from a file.
- Credentials, lazily at the first signed request, through `reqsign` (not read). Depending on the
  service and on `allow_anonymous`, that can contact the cloud's instance-metadata endpoint, STS or an
  OAuth token service. polars-bio's Rust `*ReadOptions::default()` sets `allow_anonymous: false`
  (`polars-bio/src/option.rs:353` and the others like it); the Python layer passes `True` by default.
- OpenDAL's HTTP client is `reqwest::Client::new()` (`opendal/core/core/src/raw/http_util/client.rs:50`).
  reqwest's defaults (not read) honour the `*_PROXY` variables and follow redirects, so a server's
  3xx can move a request to another host.

**What runs by itself.** On import: the Python package configures the root logger to stderr
(`polars_bio/logging.py:3-7`), sets `POLARS_FORCE_NEW_STREAMING` (`__init__.py:53-61`), and builds
one session context (`context.py:140`). The extension's module init forwards Rust logs to Python logging
(`src/lib.rs:1528`). OpenDAL's `#[ctor]` (`opendal/core/src/lib.rs:236-240`) fills an in-memory service
registry at load. None of it does I/O. No crate in the read set has a `build.rs` that is compiled in, and
none has `[build-dependencies]`.

**Compiled out, checked in `Cargo.lock`.** bigtools' remote reader (`attohttpc`, feature `remote`):
`bio-format-bbi` takes `default-features = false, features = ["read"]`, and `attohttpc` is absent from
the lock. noodles' `htsget` and `refget` clients (reqwest): not in the lock. Nothing in noodles-cram
fetches a reference over the network, as htslib does by default.

**Residual, accepted.** The not-read crates above. A later release is not re-audited, only re-checked
at its pin move (the row above). For a remote path only: the bucket-name HEAD to AWS, the credential
chains, and redirect following. None of these fires for the local files just-vcf reads.
