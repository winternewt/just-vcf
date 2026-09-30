# polars-bio: pending fixes

Behaviour of polars-bio that needs an upstream change first. Each entry notes the mitigation already in
place here. polars-bio has no local clone or inbox; filing goes to its GitHub issue tracker, and the
issue reference is recorded here once filed.

polars-bio reads VCF through two layers it does not own in-repo: `datafusion-bio-formats` (the table
provider, same org) and a biodatageeks fork of `noodles` (the parser). An entry names the layer the
defect lives in, because that is where it has to be filed.

## Upstream status sweep (2026-09-30)

Every entry below carries an *Upstream* block from one sweep: all issues and PRs, open and closed, of
polars-bio (153 issues, 320 PRs), datafusion-bio-formats (80, 178) and zaeleus/noodles (291 issues),
plus the docs, docstrings, `openspec/` designs and source of each. None of the nine is reported in any
of those trackers, and no upstream non-goal covers any of them. Nothing was filed.

The mirror is under `data/interim/` (git-ignored), for grepping:

| Path | What | Pinned at |
|---|---|---|
| `polars-bio/` | full clone, docs in `docs/` | `a1f724f` (2026-09-21) = released 0.36.0; `diff` against the installed wheel is clean |
| `datafusion-bio-formats/` | full clone | tag `v1.13.1`, the one polars-bio 0.36.0 pins |
| `noodles/` | biodatageeks fork, one commit | `42a3c01`, the rev datafusion-bio-formats pins |
| `hts-specs/` | `VCFv4.2.tex`, `VCFv4.3.tex`, `VCFv4.5.tex` | samtools/hts-specs master |
| `polars-bio-issues/*.json` | tracker dumps; `q.py '<regex>' [file…]` greps title and body | 2026-09-30 |

Refresh a dump with `gh issue list -R biodatageeks/polars-bio --state all --limit 1000 --json
number,title,state,body,labels,closedAt,createdAt,url > issues.json` (likewise `gh pr list`, and the
other two repos). The dumps hold bodies only; a maintainer's "by design" usually lives in the comments,
so read `gh api repos/<owner>/<repo>/issues/<n>/comments` for any hit before relying on it.

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
