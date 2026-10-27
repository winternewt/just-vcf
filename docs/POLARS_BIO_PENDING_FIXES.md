# polars-bio: pending fixes

Behaviour of polars-bio that needs an upstream change first. Each entry notes the mitigation already in
place here. polars-bio has no inbox; filing goes to its GitHub issue tracker, and the issue reference
is recorded here once filed. The mirror under `data/interim/` (below) is a read-only clone at the
pinned revisions, for locating code paths; it is never edited and never filed into.

No watcher follows upstream (decided 2026-10-18, at the seat's first seating). A polars-bio release
changes nothing we run until `uv.lock` moves, and the seat's start check compares the lock against the
pin below, so a pin move still triggers the re-check pass. The kit's `watch-remote.sh` would also
report every polars-bio trunk commit, and `watch-issues.py` every updated issue in the repo; both are
noise while nothing here is filed. Revisit at the first filing: a `watch-issues.py` filtered to the
`Filed:` numbers would then replace the reply check at each start.

polars-bio reads VCF through two layers it does not own in-repo: `datafusion-bio-formats` (the table
provider, same org) and a biodatageeks fork of `noodles` (the parser). An entry names the layer the
defect lives in, because that is where it has to be filed.

## Upstream status sweep (2026-09-30)

Every entry below carries an *Upstream* block from one sweep: all issues and PRs, open and closed, of
polars-bio (153 issues, 320 PRs), datafusion-bio-formats (80, 178) and zaeleus/noodles (291 issues),
plus the docs, docstrings, `openspec/` designs and source of each. None of the nine is reported in any
of those trackers, and no upstream non-goal covers any of them. Nothing was filed.

P10–P14 and the harvest section at the end come from the same day's read-only sweep of the consumer
siblings (`../just-prs`, `../just-dna-lite`, `../just-dna-format`) for polars-bio notes and workarounds,
each re-checked on 0.36.0 the same way.

The mirror is under `data/interim/` (git-ignored), for grepping:

| Path | What | Pinned at |
|---|---|---|
| `polars-bio/` | full clone, docs in `docs/` | `a1f724f` (2026-09-21) = released 0.36.0; `diff` against the installed wheel is clean |
| `datafusion-bio-formats/` | full clone | tag `v1.13.1`, the one polars-bio 0.36.0 pins |
| `noodles/` | biodatageeks fork, full history; remote `upstream` = zaeleus/noodles | `42a3c01`, the rev datafusion-bio-formats pins; 21 commits on `noodles-0.113.0` (`fe2b112`); moving it means re-extracting `READER_DEFINITIONS` (P22) |
| `opendal/` | apache/opendal, full clone (remote reads: s3, gcs, azblob, http) | `1755db1`, the rev datafusion-bio-formats pins; main, 279 commits after `v0.55.0`, not a release |
| `bigtools/` | biodatageeks fork; remote `upstream` = jackh726/bigtools | `0d7a572`, the rev polars-bio's `Cargo.lock` pins; 17 commits on upstream `b1fd283` |
| `datafusion-bio-functions/` | full clone, LFS files left as pointers | tag `v0.22.2`, the one polars-bio 0.36.0 pins |
| `hts-specs/` | `VCFv4.2.tex`, `VCFv4.3.tex`, `VCFv4.5.tex` | samtools/hts-specs master |
| `polars-bio-issues/*.json` | tracker dumps; `q.py '<regex>' [file…]` greps title and body | 2026-09-30 |

Refresh a dump with `gh issue list -R biodatageeks/polars-bio --state all --limit 1000 --json
number,title,state,body,labels,closedAt,createdAt,url > issues.json` (likewise `gh pr list`, and the
other two repos). The dumps hold bodies only; a maintainer's "by design" usually lives in the comments,
so read `gh api repos/<owner>/<repo>/issues/<n>/comments` for any hit before relying on it.

## Verdicts at a glance (2026-09-30)

| Item | Verdict | File at |
|---|---|---|
| P1 undeclared keys dropped | real bug, `physical_exec.rs:1849` `continue` | datafusion-bio-formats |
| P2 `##fileformat` not first | **our doc was wrong** — spec says "must be the first line" (`VCFv4.5.tex:142`); refusal is correct | don't file |
| P3 SV `end` one short | real bug in noodles `variant_end` (`record.rs:104`); spec example confirms END = POS+\|SVLEN\|; negative SVLEN also errors | zaeleus/noodles |
| P4 leading `\|` dropped | noodles parses it; dbf discards it (`if !first`, `:1861`/`:2497`) | datafusion-bio-formats |
| P5 `fmt_DP` vs `DP` | intended (dbf #160, pb #354), undocumented for users; related open: pb #394, dbf #258 | polars-bio (docs) |
| P6a "table already exists" | reproduced; pb has a uuid lease fix but only for Cool/MSA (`io.py:5010`) | polars-bio |
| P6b (new) | `register_vcf` returns 0-based despite 1-based docstring + global; reproduced 99 vs 100 | polars-bio |
| P7 no `format_fields`/`samples` | never existed; `register_bcf` got them in #435 — ask parity | polars-bio |
| P8 backslash in header | spec-correct strictness | low priority, don't file |
| P9 POS 0 | noodles `None` by design (noodles #343); dbf turns it into the failure | datafusion-bio-formats |
| P10 FORMAT `Number=LA/LR/LG/P/M` | spec-legal 4.5 header refuses the file; pinned noodles fork parses none; zaeleus master parses only `P` (since 0.91.0) | zaeleus/noodles (LA/LR/LG/M); dbf fork bump (P) |
| P11 FILTER `.` read as `""` | intended-ish, undocumented; "not applied" becomes a value | polars-bio (docs), low |
| P12 multi-allelic ALT joined with `\|` | intended (dbf #252), documented only in the dbf README | polars-bio (docs), low |
| P13 `write_vcf` drops undeclared columns | write-side twin of P1; answered by recipe in pb #351; related open pb #474 | polars-bio |
| P14 fork after a polars-bio call hangs the child | inherent to native thread pools (plain Polars hangs too) | polars-bio (docs), low |
| P15 haploid no-call GT `.` read as null (2026-10-22) | noodles maps `.` to missing before its GT parser; dbf writes null, same as no GT | datafusion-bio-formats |
| P16 `import polars_bio` configures the root logger (2026-10-23) | `basicConfig()` plus a root level of WARNING at import; already open as pb #291 | polars-bio (comment on #291) |
| P17 `depth()` wraps past 32,767 (2026-10-24) | an `i32` count cast `as i16` into an Int16 column | datafusion-bio-functions |
| P18 `depth()` block mode drops a run that reaches a contig's last base, 1-based only (2026-10-24) | the end marker of a 1-based read lands one past the dense array | datafusion-bio-functions |
| P19 `depth()` block mode builds every contig's batches before emitting (2026-10-24) | `merge_dense_results` fills one queue for the whole file | datafusion-bio-functions |
| P20 an indexed BAM `depth()` reads the file twice (2026-10-24) | seen on real 30x BAMs, not reproduced on a 12 MB synthetic one | unconfirmed |
| P21 a scalar INFO/FORMAT value that fails to parse is refused without its key (2026-10-25) | noodles parses scalars inside its iterator and drops the key with the error; dbf wraps it bare, while its list paths name the key | datafusion-bio-formats |
| P22 an undeclared key's bad value fails the scan though the key is dropped; an unparsable GT fails it too (2026-10-26) | the first: noodles reads an undeclared reserved key by its reserved type, and dbf parses every key before it drops the undeclared ones; the second is dbf's stated choice | datafusion-bio-formats (the first, with P21); the second by design |
| P23 `replace` with a mapping holds every column of a streaming sink (2026-10-26) | polars, not polars-bio: `replace` and `replace_strict` run as a columnar function the streaming engine does not split; already open as pola-rs/polars #24393 | pola-rs/polars (optional comment on #24393) |

## Verdicts by kind (2026-09-30)

The same items, split by what they ask of upstream. An item sits in the first table that fits it. The
*Ours* column is the work we carry whatever upstream does, in `src/just_vcf/loaders/vcf.py` unless
another path is named.

### Bugs: upstream contradicts the spec or its own docs

| Item | Defect | Owner | Ours |
|---|---|---|---|
| P1 | undeclared FORMAT/INFO keys dropped without a word; the spec only says "should" declare | datafusion-bio-formats | declare them in a rewritten copy, warn (done) |
| P3 | symbolic SV `end` = POS+SVLEN−1; negative SVLEN errors | zaeleus/noodles | compute `end` ourselves (done) |
| P4 | leading GT phase indicator discarded after noodles parsed it | datafusion-bio-formats | withhold `phased[0]` and `source_gt` on 4.4+, warn (done) |
| P6a | concurrent collects of one VCF panic on the stem-named table | polars-bio | register once under a unique name (done) |
| P6b | `register_vcf` returns 0-based starts; docstring and global setting say 1-based | polars-bio | convert `start` and check the first record (done) |
| P9 | a spec-legal POS 0 record fails the whole scan | datafusion-bio-formats | re-read without POS 0 rows, count them (done) |
| P10 | spec-legal FORMAT `Number=LA/LR/LG/M` refuses the whole file | zaeleus/noodles | `Number=P` read through a copy declaring `Number=.`, warned (`8b72a46`, text VCF only); **none yet** for LA/LR/LG/M or for BCF |
| P15 | a haploid no-call GT `.` reads as null, the same as a sample with no GT | datafusion-bio-formats | a no-call block for CGI VCFs only (`518f796`); elsewhere `gt`, `gt_alleles`, `source_gt` stay null |
| P16 | `import polars_bio` adds a stderr handler to the root logger and sets its level to WARNING; Python's logging HOWTO says a library adds only a `NullHandler` | polars-bio | the CLI turns logging off while a command runs (`cli.py:149-152`); a library user still gets the handler |
| P17 | `depth()` returns a depth over 32,767 wrapped negative; mosdepth and samtools count in 32 bits | datafusion-bio-functions | a negative depth is left out of its window and counted in a warning (`coverage.py`) |
| P18 | 1-based `depth()` in block mode emits nothing for a run of coverage that reaches a contig's last base; 0-based and per-base mode count it | datafusion-bio-functions | per-base mode, folded into windows here (`coverage.py`) |
| P22 | a key the scan drops for being undeclared (P1) still fails the whole scan when its value does not parse as the type the spec reserves for its name | datafusion-bio-formats | under 4.3+ the copy declares every undeclared reserved key `Number=.` and skips it (`c5daa78`, `READER_DEFINITIONS`); a Flag written with a value still fails, refused by name where read ahead (`b191114`) |

### Feature requests: works as designed, we want more

| Item | Ask | Owner | Ours |
|---|---|---|---|
| P7 | `format_fields` and `samples` on `register_vcf`, as `register_bcf` has | polars-bio | call private `SQL._register_variant` (done; breaks on upgrade) |
| P6b | `use_zero_based` on `register_vcf`, as `register_bcf` and `register_vcf_zarr` have (same issue as the bug fix) | polars-bio | as above |
| P10 | move the noodles fork past 0.91.0, which parses `Number=P` | datafusion-bio-formats | as above |
| P19 | block-mode `depth()` emits contig by contig, as per-base mode does (22.6 GB peak on a 27 GB BAM against per-base mode's 13.3 GB) | datafusion-bio-functions | per-base mode (`coverage.py`) |
| P13 | warn or refuse when `write_vcf`/`sink_vcf` drop an undeclared column; see open pb #474 | polars-bio | nothing until we write VCF through polars-bio |
| P21 | name the key in every INFO/FORMAT parse error, as the list-valued paths already do | datafusion-bio-formats | the retry without a named key (`loaders/vcf.py:262`, `:609`); a keyless refusal still refuses the file |
| P23 | `replace` and `replace_strict` streamed batch by batch, as an elementwise map is | pola-rs/polars | `schema.canonical_chrom_expr` maps per batch (`map_batches(..., is_elementwise=True)`, `85dcebd`); `tests/test_spill_streams.py` flips when polars streams it |
| P5 | document the `fmt_` rename, and that one sample gives flat columns while several give a `genotypes` struct | polars-bio (docs) | map columns by the header's own namespace (done) |
| P11 | document that FILTER `.` reads as `""` | polars-bio (docs) | `""` → null (done, `:455`) |
| P12 | document that multi-allelic ALT is joined with `\|` | polars-bio (docs) | split on `\|` (done, `:390`) |

### Ours: the input or the platform, nothing to ask upstream

| Item | Why it is ours | Ours |
|---|---|---|
| P2 | a late `##fileformat` makes the file invalid (`VCFv4.5.tex:142`) | reorder the header, warn (done) |
| P8 | an unescaped backslash in a quoted header value is invalid (`VCFv4.5.tex:174`); a lenient noodles mode would be a nice-to-have, not a filing | retry on a stripped header, warn (done) |
| P14 | no native thread pool survives `fork()`; plain Polars hangs the same way | never fork after polars-bio runs; use `spawn` (nothing in `src/` forks today) |

## P1 — INFO and FORMAT keys missing from the header are dropped silently

Seen with polars-bio 0.36, while declaring the VCF formats (2026-09-29). A key that a row uses but the
header does not define is dropped without a warning, GT included. A header repaired or synthesised by
just-vcf (the `vcf_no_meta` and `freebayes_varscan2_keyvalue` formats) must define every key the rows
use, or the genotypes vanish.

Mitigation: the VCF loader declares every carried FORMAT key, and INFO `END` and `SVLEN`, that the rows
use without a declaration, in a rewritten copy of the file, and says so in `Loaded.warnings`
(`loaders/vcf.py`, `readable_path`); the other undeclared keys the first rows use are listed in
`facets.dropped_fields`. `vcf_no_meta` loads this way (2026-10-22). Not filed.

*Upstream (2026-09-30).* Not reported anywhere. The drop is in datafusion-bio-formats
`bio-format-vcf/src/physical_exec.rs:1849`: the column index is built from the header, and a row key
not in it hits `else { continue; }`. The spec says the format of each INFO key "should be specified in
the meta-information" (`VCFv4.5.tex:437`), so an undeclared key is tolerated, not invalid, and dropping
it without a word is a defect. **File at datafusion-bio-formats:** warn with a per-key count, or read
the key as String.

## P2 — Headers that do not start with `##fileformat` are refused

Seen with polars-bio 0.36. A VCF whose meta lines are sorted (`##fileformat` on line 199), or that has
only a `#CHROM` line, is refused. The sniffer recognises both (`vcf`, `vcf_no_meta`); the loader will
need to reorder or synthesise the meta block before handing the file over.

Mitigation: the VCF loader rewrites the header with `##fileformat` first and warns.

*Upstream (2026-09-30).* **Not upstream's bug; do not file.** The spec requires the line first: "A
single `fileformat' line is always required, must be the first line in the file" (`VCFv4.5.tex:142`,
same in 4.3), and again "with the exception of ##fileformat which must come first" (`:2804`). The
refusal is noodles' `MissingFileFormat` (`noodles-vcf/src/header/parser.rs:263`), which reads the first
line as the file format and nothing else. Such files are invalid, and repairing them is our job. An
earlier version of this entry said the spec allows `##fileformat` later; it does not.

## P3 — `end` for a symbolic SV is one short

Seen with polars-bio 0.36 (2026-09-30). `<DEL>` with `SVLEN=30` at POS 200 gives `end=229`. VCF 4.5
computes END for a symbolic SV as `POS + SVLEN` (230), because POS is the base before the event.

Mitigation: the loader computes `end` itself (docs/SCHEMA.md, Decision 6). Not filed.

*Upstream (2026-09-30).* Not reported anywhere (zaeleus/noodles #343 and #346 touch SVLEN and END but
are different bugs). The arithmetic is noodles' `variant_end` (`noodles-vcf/src/variant/record.rs:75`):
for VCF ≥ 4.5 it takes `max(len(REF), max SVLEN, max LEN)` and returns `start + len - 1` (`:104`),
treating SVLEN as if it counted from POS. The spec's own example says otherwise: `2 321682 T <DEL>
END=321887;SVLEN=-205` (`VCFv4.2.tex:521`), so END = POS + |SVLEN|. VCF < 4.5 files take INFO/END
instead, so the bug only shows on 4.5. The same function converts SVLEN with `usize::try_from` (`:167`),
so a negative SVLEN, which the spec says to read as its absolute value, is an error there.
datafusion-bio-formats calls it unchanged (`physical_exec.rs:834`). **File at zaeleus/noodles**, the
fork's upstream; both points belong in one issue.

## P4 — The leading phase indicator of a GT is dropped

Seen with polars-bio 0.36 (2026-09-30). `|0/1/2` comes back as `0/1/2`. Since VCF 4.4 the first phase
indicator is optional and meaningful: `|0/1/2` says the first allele is phased. Dropping it loses that.

Mitigation: for VCF 4.4+ files the loader withholds `phased[0]` and `source_gt`, with a warning; older
files cannot carry a leading indicator. Reading raw GT text would restore both. Not filed.

*Upstream (2026-09-30).* Not reported anywhere. noodles parses it correctly: the first allele's phasing
is the explicit leading indicator when present, else the implicit one (`/` if any separator is `/`,
`|` otherwise; `noodles-vcf/src/record/samples/series/value/genotype.rs:33`, spec
`VCFv4.5.tex:614`). datafusion-bio-formats throws it away when it rebuilds the string: `if !first`
skips the first allele's phasing (`physical_exec.rs:1861` and `:2497`). It could write the leading
indicator exactly when it differs from the implicit one, which reproduces the source for every valid
GT. **File at datafusion-bio-formats.**

## P5 — Whether a column is INFO or FORMAT depends on what else the header declares

Seen with polars-bio 0.36 (2026-09-30), with one sample. When INFO and FORMAT both declare `DP`, the
columns are `DP` (INFO) and `fmt_DP` (FORMAT). When only FORMAT declares it, FORMAT/DP is `DP`. The
same column name means different quantities in different files.

Mitigation: the loader maps keys by the header's namespace, read by its own header parser
(docs/SCHEMA.md, Decision 3). Not filed.

*Upstream (2026-09-30).* **Intended; a documentation gap, not a bug.** The rename was added on purpose
in datafusion-bio-formats PR #160 and polars-bio PR #354 ("handle INFO/FORMAT column name collision in
single-sample VCFs"). It is documented nowhere a user reads: the `scan_vcf` docstring says
single-sample FORMAT fields are top-level columns such as `GT` and `DP` and does not mention `fmt_`.
`docs/developers.md:444` documents Arrow field metadata (`bio.vcf.field.field_type`) that tells the two
apart, but polars-bio PR #464 notes that Polars strips field metadata, so a frame cannot carry it; only
`describe_vcf` or a header parse can. Related open items: polars-bio #394 (per-sample columns, raises
INFO/FORMAT collisions for the multi-sample case; no maintainer position yet) and datafusion-bio-formats
#258 (an INFO field named like a core column breaks the schema). If filed, ask polars-bio for the rename
rule in the `scan_vcf` docstring.

The shape also depends on the sample count, checked on 0.36.0 with INFO/DP and FORMAT/DP declared: one
sample gives flat `DP` and `fmt_DP`, two samples give INFO `DP` plus a nested `genotypes` struct holding
FORMAT/DP. just-dna-lite (`just-dna-pipelines/.../io.py`, ~L255-280) still drops a FORMAT field whenever
INFO declares the same key, a workaround for a `DuplicateQualifiedField` error around polars-bio 0.19;
on 0.36.0 that workaround only loses data.

## P6 — `scan_vcf` registers a stem-named table on every run, and `register_vcf` ignores the coordinate setting

Seen with polars-bio 0.36 (2026-09-30). Two scans of one file in one plan, or `pl.collect_all` over two
frames from it, panic with `The table <stem> already exists`. A named `register_vcf` avoids it, but
returns 0-based half-open `start` where `scan_vcf(use_zero_based=False)` returns 1-based.

Mitigation: `loaders/vcf.py` `scan` registers once per source under a unique name (ported from
just-prs) and converts `start`, checking the first record against the file. Not filed.

*Upstream (2026-09-30).* Two defects, neither reported.

(a) The race. Reproduced on 0.36.0: `pl.collect_all` over two `scan_vcf` frames of one file panics
with `Failed to register VCF table: Execution("The table t already exists")`. polars-bio already has
the fix, for other formats: `_lazy_scan` registers under a per-collect `_pb_<fmt>_collect_<uuid>` lease
(`polars_bio/io.py:5010`), but only for Cool, Stockholm, A2M and A3M. VCF still re-registers under the
stem name. The lease came in with Cooler support (PR #446) and its description gives no reason for
leaving the other formats out. **File at polars-bio:** extend the lease to VCF (and the other
file-backed formats).

(b) The coordinates. `register_vcf` has no `use_zero_based` parameter, unlike `register_bcf` and
`register_vcf_zarr`. `SQL._register_variant` builds `VcfReadOptions` without `zero_based`
(`polars_bio/sql.py:211`), so the Rust default `zero_based=true` (`src/option.rs:389`) wins over the
global `datafusion.bio.coordinate_system_zero_based`, which is `false`. The docstring promises "the VCF
reader uses 1-based coordinates", and PR #263 made 1-based the global default. Reproduced: a record at
POS 100 comes back `start=99` from `register_vcf` and `start=100` from `scan_vcf`, with the global
option `false`. **File at polars-bio:** pass the resolved global setting and add `use_zero_based`.

just-prs carries the same two workarounds (a lock around a named `register_vcf` and `+1` on `start`,
`src/just_prs/vcf.py:37-63`, test `tests/test_vcf.py:98`); its commit `c424bec` (2026-07-17) is the
earliest dated record of both symptoms.

## P7 — `register_vcf` cannot restrict FORMAT fields or samples

Seen with polars-bio 0.36 (2026-09-30). `register_vcf` passes `format_fields=None` and `samples=None` to
its internal `SQL._register_variant`, so every FORMAT field of every sample is parsed. One malformed
field the loader does not even carry then fails the load: a Genos exome VCF's `CGIANN_1000GAF` holds a
value that is not a float.

Mitigation: `loaders/vcf.py` `scan` calls `polars_bio.sql.SQL._register_variant` directly with the
carried FORMAT keys and the chosen sample. That is a private API; a polars-bio upgrade can break it, and
the tests will show it. Filing upstream: ask for `format_fields` and `samples` on `register_vcf`. Not filed.

*Upstream (2026-09-30).* Not reported, and not a removal: `register_vcf` never had these parameters
(checked at `4b73db1^`). PR #435 ("native BCF APIs with VCF parity") gave them to `register_bcf` and
left `register_vcf` hardcoding `None`, so the ask is parity with `register_bcf`. Some issue text calls
`register_vcf(format_fields=…)` as if it existed (datafusion-bio-formats #88–#91, open polars-bio PR
#334 with `samples=`); those calls would fail with a `TypeError` on 0.36.0, which is worth citing.
**File at polars-bio.** P6(b) is the same function and can go in the same issue.

## P8 — A backslash inside a quoted header value refuses the whole file

Seen with polars-bio 0.36 (2026-09-30) on a GATK 3.4 UnifiedGenotyper VCF: its
`##GATKCommandLine.UnifiedGenotyper=<…,CommandLineOptions="…">` line holds a Windows path with raw
backslashes, and registration fails with `IO error: invalid record`. A minimal repro is one such line,
`CommandLineOptions="a=C:\Users\x\ref.fa"`; the same line with `/` reads fine. VCF 4.5 does say a
backslash in a quoted value must be escaped, so the file is strictly invalid, but it is how GATK wrote
them.

Mitigation: when registration is refused, the loader retries once on a cached copy whose header keeps
only `##fileformat`, INFO, FORMAT, FILTER, ALT and contig lines, and warns. The full header is still read
by the loader's own parser for the facet. Not filed.

*Upstream (2026-09-30).* **Spec-correct; not a bug.** noodles' `parse_escaped_string` accepts only
`\\` and `\"` and returns `InvalidEscapeSequence` for anything else
(`noodles-vcf/src/header/parser/record/value/map/field/value/string.rs:74`), as `VCFv4.5.tex:174`
requires. Not reported anywhere. The only fair ask is a lenient mode that keeps an unknown escape
literally, at zaeleus/noodles; low priority while our retry works.

## P9 — A record at POS 0 fails the scan

Seen with polars-bio 0.36 (2026-09-30) on Illumina `GenotypingToVCF` exports, which put unmapped probes
at `CHROM NA`, `POS 0`, `REF N`: collecting fails with `Execution error: Missing variant start`. VCF 4.5
allows POS 0 (a telomere).

Mitigation: on that error the loader reads a copy without the POS 0 records and says how many were left
out. Not filed.

*Upstream (2026-09-30).* The error text is datafusion-bio-formats' (`physical_exec.rs:1173` and
`:3063`): noodles returns `None` from `variant_start()` for POS 0, and the reader turns `None` into a
failure. The noodles side is by design. zaeleus/noodles #343 (open since 2025-06) covers POS 0; the
maintainer answered that the record parses fine and the caller should handle it, and a commenter cited
the telomere rule (`VCFv4.5.tex:396`). **File at datafusion-bio-formats:** map `None` to position 0, or
at least skip the record with a count, instead of failing the scan. Cite noodles #343.

*Re-checked (2026-10-18), after the stack audit read the noodles fork's delta.* The fork does change POS 0,
but only in the `RecordBuf` parser (`noodles-vcf/src/io/reader/record_buf/position.rs`, 21 commits over
`noodles-0.113.0`): there `0` is now invalid and `.` is missing. datafusion-bio-formats reads VCF
through the lazy `noodles_vcf::Record` (`physical_exec.rs:1173`), whose `variant_start`
(`noodles-vcf/src/record/fields.rs:22-36`) the fork did not touch and which still returns `None` for
`0`. The owning layer stays datafusion-bio-formats. A synthetic two-record 4.5 file with one POS 0 row
still fails `scan_vcf(...).collect()` on 0.36.0 with `Execution error: Missing variant start`.

## P10 — FORMAT `Number=LA`, `LR`, `LG`, `P` or `M` refuses the whole file

Found 2026-09-30 through just-dna-format (`docs/probes/METHYLATION_PACBIO.md:286-291`, quoting rastair,
which writes `Number=.` instead to dodge exactly this). VCF 4.4/4.5 add these FORMAT-only Number values
(`VCFv4.5.tex:191-199`; `M` for base modifications, `:203-227`). On 0.36.0 a header declaring any of
`LA`, `LR`, `LG` or `M` fails registration with `IO error: invalid record`, which does not point at the
header. `P` was not run but takes the same path.

Mitigation: none yet. The cached-copy header rewrite could turn these into `Number=.` and warn.

*Upstream (2026-09-30).* No issue in any tracker. The pinned noodles fork's `parse_number`
(`noodles-vcf/src/header/parser/record/value/map/format/number.rs:32`) accepts only `A`, `R`, `G`, `.`
and integers. The `Number` enum already has variants for all five (`header/record/value/map/format/number.rs`).
zaeleus/noodles master parses `P` for VCF ≥ 4.4 (commit `a8f7445`, noodles-vcf 0.91.0) and still not
`LA`, `LR`, `LG` or `M`. **File at zaeleus/noodles** for the four; `P` needs datafusion-bio-formats to
move its noodles fork past 0.91.0.

*Re-checked (2026-10-18), after the stack audit read the noodles fork's delta.* The fork relaxes header
validation (`noodles-vcf/src/header/parser/record/value.rs:230-234`, `:254-256`), but only the check of a
reserved key's declared `Number` against the spec's table, which runs after parsing. `parse_number`
(`number.rs:32-41`) is unchanged and still rejects `LA`, `LR`, `LG`, `P` and `M`. A synthetic 4.5 file
declaring `FORMAT=<ID=LAD,Number=LR,...>` still fails registration on 0.36.0 with `IO error: invalid
record`. Verdict unchanged.

*Mitigated for `P` (2026-10-24), from the lead seat's PSL work (merged as `8b72a46`).* polars-bio 0.36 still
refuses a `Number=P` header (the PSL builder re-confirmed it). The text-VCF loader now reads any FORMAT
key declared `Number=P` (PSL, PSO, PSQ in practice) through its rewritten copy, declared `Number=.`, and
checks each value's count against GT (`READ_AS_ANY_NUMBER`, `loaders/vcf.py:83`; `REWRITE_VERSION` 9),
with a warning naming the keys (`loading.py:712`). `LA`, `LR`, `LG` and `M` are not rewritten and still
refuse the file, and a BCF declaring `Number=P` still refuses: the BCF path reads no rewritten copy. The
filings are unchanged.

## P11 — FILTER `.` comes back as an empty string

Recorded in just-prs (`src/just_prs/normalize.py:50`) and in our CLAUDE.md §5; reproduced on 0.36.0: a
record with FILTER `.` gives `filter == ""`. The spec says `.` means filters have not been applied
(`VCFv4.5.tex:426`), which is unknown, not a value.

Mitigation: the loader has to map `""` to null. No valid record has an empty FILTER, so the mapping is
lossless.

*Upstream (2026-09-30).* No issue, no docs. noodles yields no items for `.` (`noodles-vcf/src/record/filters.rs:48`)
and datafusion-bio-formats joins that empty iterator into `""` (`physical_exec.rs:992`, also `:1263`,
`:1531`, `:3165`); its `unwrap_or(".")` never fires. **Intended in effect, undocumented.** A docs note
at polars-bio, low priority; asking for null would be a breaking change for them.

## P12 — Multi-allelic ALT is joined with `|`

Recorded in just-prs (`normalize.py:67`), just-dna-lite (CLAUDE.md) and our CLAUDE.md §5; reproduced on
0.36.0: `ALT=T,G` gives `alt == "T|G"`.

Mitigation: the loader splits `alt` on `|`.

*Upstream (2026-09-30).* **Intended.** datafusion-bio-formats joins with `'|'` (`physical_exec.rs:990`
and the other readers) and PR #252 pins the encoding in a test. It is documented only in the
datafusion-bio-formats README (`bio-format-vcf/README.md:245`, the `vcf_ac`/`vcf_af` UDFs); the
`scan_vcf` docstring and polars-bio docs do not say it. A docs note at polars-bio, low priority.

## P13 — `write_vcf` and `sink_vcf` drop columns the header metadata does not declare

Recorded in just-dna-lite (CLAUDE.md:561-563, `annotation/vcf_export_logic.py:57-105`); reproduced on
0.36.0: a frame with a `score` column and no source metadata is written with INFO `.`, no `##INFO` line
and no warning. Declaring `score` through `set_source_metadata` writes it. The write-side twin of P1.
Missing core columns and a non-UInt32 `start` do raise, loudly.

Mitigation: none needed until just-vcf writes VCF through polars-bio. When it does, declare every
carried field and check the written header against the frame.

*Upstream (2026-09-30).* polars-bio #351 (closed) got the `set_source_metadata` recipe as its answer;
#464 (merged) added the `INFO_<id>` convention, which also needs the id declared; #474 (open) proposes
an explicit column-to-id mapping. A column declared in neither INFO nor FORMAT metadata passes through
`src/write.rs:395-405` unmarked and the serializer skips it. **File at polars-bio:** warn or refuse with
the dropped column names; cite #351 and #474.

## P14 — `fork()` after a polars-bio call hangs the child

Recorded in just-dna-lite (`docs/GRANIAN_POLARS_FORK_DEADLOCK.md`, `docs/DAGSTER_GUIDE.md:179`);
reproduced on 0.36.0 with synthetic input: a child forked after a warm `scan_vcf` or `describe_vcf`
hangs on its next polars-bio call and ignores SIGTERM. A child forked after a plain Polars sort hangs on
the next sort too, and `import polars_bio` alone is harmless.

Mitigation: just-vcf must not use fork-based multiprocessing once polars-bio has run; use `spawn`.

*Upstream (2026-09-30).* Nothing in any tracker or doc. polars-bio builds Tokio runtimes in its Rust
entry points (`src/scan.rs:1101`, `:1176`), and neither Tokio nor Rayon survives `fork()`. **Inherent,
not a bug.** At most a "use spawn" note in polars-bio docs.

## P15 — A haploid no-call GT `.` comes back as null

Found by the lead seat on a PGP-HMS Complete Genomics vcfBeta (merged as `bd5c02e`): on male X and Y the
no-called records write GT `.`, and 139,310 of them loaded with no genotype at all. VCF says a no-call
writes `.` for each missing allele, "`./.` for a diploid genotype and `.` for haploid genotype"
(`VCFv4.5.tex:619`), so `.` is a called-as-unknown haploid genotype, not an absent field. polars-bio
0.36 returns it as null, the same value as a sample whose FORMAT has no GT, so the ploidy and the fact
that the caller wrote a no-call are both lost. `./.` and `.|1` survive.

Mitigation: for Complete Genomics VCFs only (a header declaring `CGA_NOCALL`), a record with ALT `.`, no
`END` and a null GT goes to `blocks` as a no-call span (`518f796`, `loaders/vcf.py`, `to_rows`). In any
other VCF a haploid `.` loads with `gt`, `gt_alleles` and `source_gt` null. Not filed.

*Upstream (2026-10-22).* Not reported in any of the three trackers (searched live for "haploid"; noodles
#328, about BCF genotypes, is unrelated). Reproduced on 0.36.0 with a synthetic 4.5 file: GT `.` reads
as null, exactly as a record whose FORMAT is `DP` alone. The information goes in noodles: the lazy
record's `parse_value` returns `None` for `.` before it checks for the GT key
(`noodles-vcf/src/record/samples/series/value.rs:24-29` at the fork's `42a3c01`; zaeleus master has
the same order), so the genotype parser, which reads `.` as one missing allele, never sees it. datafusion-bio-formats
then writes null for any GT value that is not a `Genotype` (`physical_exec.rs:1853-1877`, and
`:2489-2514` for several samples). It can tell the two apart without a noodles change: a sample's
`iter` zips the FORMAT keys with the values present (`noodles-vcf/src/record/samples/sample.rs:46-51`),
so a GT key yielded with `None` means the text was `.`, and writing `.` there reproduces the source. A
sample written as `.` alone would then read as a haploid no-call too. **File at
datafusion-bio-formats**, citing the noodles ordering; the same one-line reorder in the biodatageeks fork
or at zaeleus/noodles would fix it at the root.

## P16 — `import polars_bio` configures the root logger

Found by the lead seat building the CLI (merged as `5b8700f`): after `import polars_bio`, every
`logging` warning in the process, just-vcf's included, prints to stderr, whatever the application set
up. polars-bio's `polars_bio/logging.py:3-5`, run on import through `__init__.py:26`, calls
`logging.basicConfig()` and sets the root logger to WARNING; `set_loglevel` changes the root level
again (`:30`, `:34`, `:38`). Reproduced on 0.36.0: a root logger with no handlers at DEBUG comes out
of the import with a stderr `StreamHandler` at WARNING. Python's logging HOWTO ("Configuring Logging
for a Library") says a library adds no handler but a `NullHandler` and leaves the root logger to the
application.

Mitigation: the CLI turns logging off while a command runs and back on when it closes
(`src/just_vcf/cli.py:149-152`), because a library error's message must not reach stderr unredacted.
A program that imports just-vcf as a library gets polars-bio's handler. Not worked around in the
library: it would mean undoing another package's import side effect. Not filed.

*Upstream (2026-10-23).* Already reported: pb #291, "How do I get polars-bio logging to not stomp on my
colorlog logging?" (open since 2026-02-10, labelled "analysis needed"; the maintainer said they would look,
the reporter asked again on 2026-09-17). It gives no cause and no fix. Upstream master after 0.36.0
leaves `logging.py` unchanged (0.36.2 is the newest tag). **Comment on pb #291**, not a new issue:
name the two lines, and the fix: drop `basicConfig()` and the root level, give the `polars_bio`
logger a `NullHandler`, and let `set_loglevel` set that logger's level alone.

## P17 — `depth()` wraps a depth past 32,767

Found by the lead seat building alignment coverage (`src/just_vcf/coverage.py`). polars-bio's
`depth()` comes from datafusion-bio-functions' pileup crate (tag `v0.22.2`, which polars-bio 0.36.0
pins). It counts in `i32` and writes the column as `Int16` with a wrapping cast: `prev_cov as i16` in
block mode (`bio-function-pileup/src/coverage.rs:76`, `:144`, `:189`) and `self.cumulative_cov as i16`
in per-base mode (`:374`), into the `Int16` field of both output schemas (`schema.rs:33`, `:52`). The
type is documented (the crate README's tables, `polars_bio/pileup_op.py:98-100`), the wrap is not.
Reproduced on 0.36.0 with a synthetic BAM: 32,770 reads over five bases come back as depth -32,766, in
both modes. Amplicon panels, mitochondria and centromeric repeats pass 32,767 on real data; mosdepth
and `samtools depth` count in 32 bits.

Mitigation: a negative depth is left out of its window and the count of such bases is warned
(`coverage.py`, `test_a_depth_past_int16_is_left_out_and_counted`). A depth past 65,535 would wrap back
to a positive value and is not caught. Not filed.

*Upstream (2026-10-24).* Not reported (polars-bio, datafusion-bio-functions and datafusion-bio-formats
searched live for "depth", "coverage", "pileup" and "Int16"). The pileup source is unchanged on
datafusion-bio-functions main after `v0.22.2` (newest tag `v0.23.1`). **File at
datafusion-bio-functions:** an `Int32` (or `UInt32`) coverage column, or at least a saturating cast.

## P18 — `depth()` block mode drops a run that reaches a contig's last base (1-based)

Found by the lead seat (`test_coverage.py::test_a_read_on_a_contigs_last_base_is_counted`; on a real
BAM the last window of 60 contigs read zero). Reproduced on 0.36.0 with a synthetic BAM (contigs of 100,
50 and 40 bases): with reads on chrX 30–40 and 40–50 and one on all of chrY, 1-based block mode returns
only chr1's block, nothing for chrX or chrY; `use_zero_based=True` returns all three, and
`dense_mode="disable"` (the sparse path) returns all three in 1-based too. A contig whose reads stop
one base short is fine.

The cause: polars-bio builds the BAM provider in the caller's coordinate system (`src/pileup.rs:166-168`),
so 1-based reads arrive with 1-based starts. The dense array is sized for 0-based ones, the contig's
length plus one sentinel (`bio-function-pileup/src/events.rs:184`). A 1-based read ending on the last
base puts its `-1` at length + 1, which `apply_binary_cigar_to_depth` skips as out of range
(`cigar.rs:206-209`). The run never closes, and the block walk emits a block only when the depth
changes and has no flush at the end (`coverage.rs:165-201`), so the whole run is lost, back to the
last base where the depth was zero. Per-base mode does not need the closing delta. A related fix,
datafusion-bio-functions #204 and pb #427 (0-based blocks made half-open), touched only the block ends.

Mitigation: just-vcf uses per-base mode and folds it into windows (`coverage.py`). Not filed.

*Upstream (2026-10-24).* Not reported (same searches as P17); the pileup source is unchanged after
`v0.22.2`. **File at datafusion-bio-functions:** size the dense array for the coordinate system in use
(or shift 1-based starts to 0-based before accumulating), with the repro. Also worth naming: the
dense merge drops a contig whose batches fail to build without an error (`physical_exec.rs:420-430`,
`if let Ok(...)`).

## P19 — `depth()` block mode builds every contig's batches before emitting any

Found by the lead seat on 30x BAMs: block mode peaks at 22.6 GB on a 27 GB BAM and 33.6 GB on a 48 GB
one, per-base mode at 13.3 GB on the first. In the pileup crate's dense merge, block mode turns every
contig into record batches and queues them all before it emits the first
(`bio-function-pileup/src/physical_exec.rs:418-435`); per-base mode queues the contigs and builds each
one's batches as it is emitted (`:406-416`). Both keep every contig's dense array, four bytes per base,
until the input is read. Not reproduced here, since a synthetic BAM small enough to build is too small
to show it; the code path is plain from the source.

Mitigation: per-base mode (`coverage.py`, the comment at `:596-598`). Not filed.

*Upstream (2026-10-24).* Not reported; datafusion-bio-functions #2 (closed) added the batch-size
chunking that block mode now uses, but not streaming. **File at datafusion-bio-functions**, together
with P18 if it is filed: emit block-mode batches contig by contig, as per-base mode does.

*Filing P17–P19 (2026-10-27).* The person's answer on desk #104: file all three at
datafusion-bio-functions, but draft them first and wait for the person to vet the drafts. The three
drafts are on desk #104 (comment of 2026-10-27), with the code linked at main `ea47537`, where the
pileup crate (`datafusion/bio-function-pileup/`) is unchanged since `v0.22.2` but for its Cargo.toml.
The repros for P17 and P18 were rerun as written on 0.36.0 with synthetic BAMs, and the trackers were
searched again: no duplicates. Nothing is posted until the person says so; record `Filed:` here then.

## P20 — An indexed BAM `depth()` scan reads the file twice

Found by the lead seat on real 30x BAMs: the bytes read come to about twice the file's size when it has
an index. **Not reproduced**: on a 12 MB synthetic BAM (two 5 Mb contigs, 300,000 reads), `depth()` reads
1.00–1.01 times the file, with or without its `.bai`, in both modes. The likely place is the BAM
provider's index-driven partitioning in datafusion-bio-formats (`bio-format-core/src/partition_balancer.rs`
splits a long contig into ranges), which a small file may not trigger. Owner unknown until it is
reproduced; nothing to file yet. Mitigation: none needed for correctness.

## P21 — A scalar INFO/FORMAT value that fails to parse is refused without its key

Found by the lead seat on the VCF loader (merged `20fb9b0`). A scalar value that does not parse as its
declared type (FORMAT `GQ` `Number=1,Type=Integer` holding `30.5`, INFO `END=10.5`) fails the scan with
`Error reading FORMAT field: invalid digit found in string` or `Error reading INFO field: ...`, no key
named. The same fault in a list-valued key names it: `Error reading FORMAT field 'AD': ...`, `Error
reading INFO array field 'SVLEN': ...`. Reproduced on 0.36.0 with four one-record synthetic VCFs, one per
case.

The cause is in two layers. noodles parses a scalar value inside its field iterator and returns
`io::Result<(key, value)>`, attaching the key only on success (`noodles-vcf/src/record/info.rs:62`,
`record/samples/sample.rs:50`, fork `42a3c01`), and its message names no key. A list value is parsed
lazily, after the key is known. datafusion-bio-formats then wraps the iterator's error bare
(`bio-format-vcf/src/physical_exec.rs:736` for INFO, `:1846` and `:2474` for FORMAT, v1.13.1), while its
own helpers for the list paths name the key (`info_array_error`, `format_value_error`, `:671-685`).
Unchanged on its main after v1.14.0 (`:736`, `:1852`, `:2480`).

Mitigation: the loader reads again without a key polars-bio names in its refusal (`_UNPARSED_FIELD`,
`loaders/vcf.py:262`; `can_leave_out`, `:609`). A keyless refusal gives it nothing to leave out, so such
a file still refuses as a whole. Not filed.

*Upstream (2026-10-25).* Not reported (polars-bio, datafusion-bio-formats and both noodles trackers
searched live for the two messages and "invalid digit"; pb #380, closed, was a bare INFO key, not this).
**File at datafusion-bio-formats**, which can fix it alone: on the error path, find the failing key by
reading each selected key with `info.get(header, key)` / `sample.get(header, key)`, which costs nothing
on a good record, and name it as the list paths do. A noodles change (the key in the iterator's error)
would be the deeper fix, and could be mentioned. Goes well in one batch with P15, also
datafusion-bio-formats.

## P22 — An undeclared key's bad value fails the scan; so does one unparsable GT

Found by the lead seat fixing genotype reading, the noodles half traced by the seat-builder. Two
behaviours, reproduced on 0.36.0 with synthetic two-record VCFs (VCF 4.3, one good record first):

- **An undeclared key that the scan drops can still fail it.** FORMAT `GQ` used but not declared, with
  the value `30.5`, fails the load (`Error reading FORMAT field: invalid digit found in string`, no key
  named, as in P21), although an undeclared key is never a column (P1). An undeclared INFO `END=10.5`
  fails it the same way once any INFO key is declared; with none declared the INFO pass is skipped and
  the file loads. An undeclared key the spec does not reserve (`FOO=30.5`) loads, read as a String and
  dropped. The cause: noodles takes an undeclared key's Number and Type from the definition the file's
  VCF version reserves for that name (`.or_else(|| definition(header.file_format(), key))`,
  `noodles-vcf/src/record/samples/series/value.rs:36`, `record/info/field.rs:23`, fork `42a3c01`), which
  is reasonable; datafusion-bio-formats iterates every key of a sample or an INFO column, parsing each,
  and only then skips the ones it has no column for (`bio-format-vcf/src/physical_exec.rs:736`, `:1846`
  and `:2474` against the `field_to_index` lookup after them). It iterates once on purpose, to avoid a
  `get` per key (the comment at `:666-667`), so the unselected keys are parsed too.
- **One unparsable GT fails the whole scan.** `GT=0/x` gives `Error reading FORMAT field 'GT': invalid
  digit found in string`. noodles parses GT lazily; datafusion-bio-formats refuses the scan on the first
  allele that does not decode, and says why in its own comments: mapping undecodable data to null "would
  present corrupt data as merely absent" (`:668-685`). By design, not a defect. just-vcf's text loaders
  withhold such a genotype instead (`loaders/text.py:123`); the VCF loader cannot, since the refusal comes
  before any row and GT is never left out (`can_leave_out`, `loaders/vcf.py:609`).

Mitigation: under VCF 4.3+ the structural copy declares every key the reader reserves and the file does
not declare, carried or not (FORMAT `MQ`, INFO `DP`), `Number=.` with the reader's Type for that version,
so it is read as text and skipped (`7c37f1a`, `c5daa78`; `READER_DEFINITIONS`, `loaders/vcf.py:176`,
`REWRITE_VERSION` 17). Flags are the exception: noodles refuses a Flag declared with a Number other than
0, so `DB=1` still fails a 4.3+ read, declared or not: a file's own `DB` line with another Type is
refused for the reserved name and left out of the copy, and noodles reads `DB` as its reserved Flag. The
loader refuses such a file naming the key when the records it reads ahead show it, and otherwise lists
the version's Flags (`named_flag_refusal`, `b191114`). The stats pass reads every INFO column the frames
read, so an INFO failure reaches the retry rather than escaping from the spill (`5f976d3`). The GT case
still refuses the file. Not filed.

`READER_DEFINITIONS` is copied by hand from the noodles fork at `42a3c01`
(`header/record/value/map/{info,format}/definition/v4_*.rs`): **re-extract it whenever the noodles pin
moves.**

*Upstream (2026-10-26).* Neither is reported (the searches of P21). **File the first at
datafusion-bio-formats, with P21**: one change covers both. On an error from the key iterator, find the
failing key with `get` (cost only on a bad record); if it has no column, skip it and go on; if it has,
fail naming it. The GT case is their stated choice; at most a wish for an opt-in lenient read that nulls
an undecodable value and counts it. Not worth filing on its own.

*Filing P15, P21 and P22(a) (2026-10-27).* The person's answer on desk #140, which covers P21
and P22(a): draft first and wait for the person to vet, as for P17–P19. Two drafts are on desk #140
(comment of 2026-10-27): P21 with P22(a) as one issue, and P15 as a second, since #140's recommended
option files it with them (P15's own question, desk #69, is unanswered). Code is linked at dbf main
`6af1fa1` and the noodles fork `42a3c01`, both unchanged for these paths; the repros were rerun as
written on 0.36.0 with synthetic VCFs, and the trackers searched again: no duplicates. P22(b), the GT
`0/x` case, is named in the draft as not asked. Nothing is posted until the person says so.

## P23 — `replace` with a mapping holds every column of a streaming sink

Not polars-bio's: polars itself, which polars-bio and just-vcf both build on. Found by the lead seat as
F26 (`docs/PREVIOUS_ISSUES.md`): a WGS gVCF `load` held 8 GB for 2.6 GB of rows. In a streaming sink,
polars 1.44 runs `pl.col(x).replace(mapping)` as a `columnar-function` node behind an in-memory engine
fallback, so the stream splits there and every other column of the plan is held until the function
has seen its whole column.

Reproduced on 1.44.2 with synthetic data only (30M rows; a constant string column mapped, an integer
and a formatted string beside it; `sink_parquet`): peak RSS 0.84 GB without the replace, 2.74 GB with
it, 0.91 GB with the same mapping inside `map_batches(..., is_elementwise=True)`. The physical graph
(`show_graph(engine="streaming", plan_stage="physical")`) shows the fallback node for
`replace(mapping)`, `replace(old_list, new_list)` and `replace_strict(mapping, default=...)` alike, and
none for the `map_batches` form.

Mitigation: `schema.canonical_chrom_expr` maps one batch at a time (`85dcebd`).
`tests/test_spill_streams.py` asserts a plain `replace` still falls back, so it fails once polars
streams it and the plain expression can return. Any new `replace`/`replace_strict` in a plan that
should stream needs the same form (CLAUDE.md §5 says so).

*Upstream (2026-10-26).* Already open as **pola-rs/polars #24393** (2025-09-07, labels `bug`,
`needs triage`, no maintainer reply): "Streaming `replace_strict` does in-memory-map, unless
`eval`uated in a list"; two later comments say it persists in 1.42.0 and show a `when/then` chain as
a workaround, and it links #23040 (the same function's memory inside `list.eval`). The issue names
only `replace_strict`; plain `replace` goes through the same fallback, which is the one thing a
comment would add, with 1.44.2 and the `map_batches` form. Nothing to file new. Not put to the
person on its own: the mitigation holds, and the issue's fix covers ours if `replace` shares the
code path, which is likely but not read in the source. Comment only alongside another batch.

## Harvested and closed without an entry (2026-09-30)

Items the sibling sweep found that need nothing from upstream now.

| Source | Item | Verdict |
|---|---|---|
| just-dna-lite `just-dna-pipelines/.../io.py:113-116` | INFO `ALLELE_ID` (Number=R with `.`) panicked | fixed upstream (pb #312/#314, dbf #74/#178); passes on 0.36.0 |
| just-dna-lite `io.py:308-312` | predicate pushdown bug in `_translate_in_expr` (0.23) | fixed by pb #407; semi-join and `is_in` pass on 0.36.0 |
| just-prs `CLAUDE.md:593`, just-dna-lite `CLAUDE.md:559` | PGEN 512 MiB companion cap (pb #453) | fixed by pb #454 in 0.36.0; #453 still shows open |
| just-dna-lite `webui/.../compute/tasks.py:21` | `scan_vcf` plan cannot be serialized or pickled | inherent to Polars Python IO sources (needs cloudpickle); not polars-bio's |
| just-prs `CLAUDE.md:592` | "all columns in a record batch must have the same length" on multi-column aggregation | not reproduced on 4-row or 300k-row inputs; no version or query recorded |
| just-prs `sample_scores/publish.py:226,248` | "invalid BGZF" on some `.vcf.gz` | not reproduced on plain, nameless, multi-member or misnamed gzip; needs the failing file's first 18 bytes |
| just-prs `vcf.py:49-55` | `register_vcf(chunk_size=8, concurrent_fetches=1)` | no defect recorded; tuning |
| just-dna-format | — | does not use polars-bio; its one indirect item is P10 |

Two sibling workarounds are now harmful rather than stale: just-dna-lite's INFO/FORMAT key drop (see P5)
and its `ALLELE_ID` blocklist. Both belong in just-dna-lite's `docs/FEEDBACK.md`, not here.
