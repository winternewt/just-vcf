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
| `noodles/` | biodatageeks fork, full history; remote `upstream` = zaeleus/noodles | `42a3c01`, the rev datafusion-bio-formats pins; 21 commits on `noodles-0.113.0` (`fe2b112`) |
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
| P10 | spec-legal FORMAT `Number=LA/LR/LG/M` refuses the whole file | zaeleus/noodles | **none yet**: rewrite to `Number=.` and warn |

### Feature requests: works as designed, we want more

| Item | Ask | Owner | Ours |
|---|---|---|---|
| P7 | `format_fields` and `samples` on `register_vcf`, as `register_bcf` has | polars-bio | call private `SQL._register_variant` (done; breaks on upgrade) |
| P6b | `use_zero_based` on `register_vcf`, as `register_bcf` and `register_vcf_zarr` have (same issue as the bug fix) | polars-bio | as above |
| P10 | move the noodles fork past 0.91.0, which parses `Number=P` | datafusion-bio-formats | as above |
| P13 | warn or refuse when `write_vcf`/`sink_vcf` drop an undeclared column; see open pb #474 | polars-bio | nothing until we write VCF through polars-bio |
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

Mitigation: the VCF loader declares every carried FORMAT key the rows use without a declaration, in a
rewritten copy of the file, and says so in `Loaded.warnings` (`loaders/vcf.py`, `readable_path`). Not filed.

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
