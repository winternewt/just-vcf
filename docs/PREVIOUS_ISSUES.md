# Previous issues — resolved here

Findings from [DOGFOODING.md](DOGFOODING.md) resolved in this repo, each with its resolution and a
code pointer. Check here before re-investigating something that looks fixed.

## F23 — A CNV caller's `N` REF reads as a mismatch, so no CNV facet lifts

Severity: medium · Found 2026-10-26 dogfooding RM6's SV lift on the PGP-HMS CNV files.

DRAGEN's and Canvas's copy-number VCFs write REF `N` on every record: the caller never reads the
reference. `resolve_ref` compares that `N` with the base it reads and sets `ref_check = mismatch`
on every row, and `liftover` leaves every `mismatch` row behind as `source_ref_mismatch`. Measured
with the library: hu6ABACE's DRAGEN 3.2.8 `.cnv.` (586 rows, GRCh37), hu1C1368's DRAGEN 3.4.5
`.cnv.` (556 rows, GRCh37) and hu416394's Canvas `.cnv.` (182 rows, 92 copy-number calls and 90
copy-neutral segments, GRCh38): 0 match, 100% mismatch, 0 rows lifted through `hg19ToHg38` or
`hg38ToHg19`. So the copy-number half of RM6's SV lift has never run on a real file; the 31 `N` rows
of hu416394's Manta file sit on decoys and are unread anyway.

The IUPAC note in the resolve warning does not cover it: `_iupac_explains` counts an ambiguity code on
the reference side (a gap base) against a concrete source base, not a source `N` against a concrete
reference base. The SV lift itself needs nothing new: on a plus strand `_place_sv` compares the
target's base with `ref`, the base read from the source assembly, not with `source_ref`. Only the
`ref_check` gate stops these rows. Each candidate repair needs a decision:

- *Read a source `N` as matching any base (`match`):* the row moves from honestly unverified to
  apparently verified, a check that never ran reported as passed (§5, Kleene).
- *A new `RefCheck` value for a source REF that states no base (`unstated` or similar), lifted like
  `match`:* the honest reading, but `RefCheck` is a `pl.Enum` in the wire schema, so a new value is a
  schema change, and `mint_vrs_ids`, `write_vcf(resolved_only=True)` and `normalize` each have to say
  what they do with it.
- *`not_run` for an all-`N` source REF, with `ref` read:* no schema change, but `not_run` today means
  the comparison did not happen for want of a reference, and the two causes would share one value.
- *Special-case symbolic rows in `liftover`:* makes SVs the one row kind whose unstated REF passes,
  and leaves `resolve_ref`'s counts calling 1,324 public rows mismatched.

Blocks: lifting any CNV facet (every one in hand); also any QC that reads `ref_check` on them.

**Resolved 2026-10-27** (desk #151, the person's choice: the second repair above). `RefCheck` gains
`unstated` (schema 0.0.22): `resolve_ref` writes it for a source REF made only of `N`, any case, where
the reference reads a base, counts those rows apart in its warning (never as mismatches, never as
IUPAC), and records `RESOLVE_VERSION` 13. `liftover` gates on `mismatch` only, so an `unstated` row
lifts as a `match` one, the target's base compared with `ref`. A partial `N` (`AN`) states a base and
stays `mismatch`; an `N` over a reference gap stays `match` (the letters agree, and RM271's ids over
`N` bases stand). Readers that carry the alleles as written take it as `match` (`write_vcf`,
`write_vcf_zarr`); readers that compute a new representation withhold (`mint_vrs_ids` mints the REF
slot and withholds a sequence ALT as `ref_unstated`; `normalize` leaves the row as written); `query`
never answers `site_genotyped` from a REF the source did not state. Measured with the library after:
hu6ABACE 586 of 586 rows `unstated`, 568 lifted through `hg19ToHg38`; hu1C1368 556 `unstated`, 538
lifted; hu416394 182 `unstated` (90 segments, 92 symbolic), 118 lifted through `hg38ToHg19`; the rest
behind for chain reasons, none as `source_ref_mismatch` or `ref_mismatch_after_lift`. Code:
`reference._resolve_rows` and `_count_vcf`, `vrs._slot_reason`, `normalize._align`'s reasons,
`query._presence`, `store_migration.MIGRATIONS` (0.0.21 → 0.0.22); tests in
`tests/test_unstated_ref.py` (synthetic, and an opt-in count over every PGP-HMS `.cnv.` upload);
docs/SCHEMA.md Decision 2 says what each reader does with the value.

**Amended 2026-10-28** (merge audit 13): `liftover`, `write_vcf`, `write_vcf_zarr` and `query` still
read an `unstated` row's sequence ALT, which minting already withheld. One rule now holds in every
reader, `schema.unstated_sequence_alt`: such a row stays behind in a lift (`source_ref_unstated`,
`LIFT_VERSION` 13), is unwritable (`ref_unstated`), and answers `overlapping_call` in `query`, never
`stated`. Symbolic and ALT-less rows are unaffected: the three files above lift 568, 538 and 118 rows
as before, none left behind as `source_ref_unstated`.

**Amended 2026-10-31** (the two readers audit 13 left): `concordance` and `consensus` keyed an allele
on `ref`, else `source_ref`, and never read `ref_check`, so `N` → `G` over `T` was compared as `T` →
`G`. `concordance` now counts such an ALT in `excluded` as `ref_unstated` and reads such a row for no
other facet's allele at its position (`_alleles`, `_one_sided`); `consensus` states no allele from it
and reads it as silent, `ref_unstated`, for another facet's allele there, never as `site`, `ref_block`
or `absence` (`_stated_alleles`, `_readings`, `CONSENSUS_VERSION` 5). A symbolic ALT on such a row
still meets by type and extent. Tests in `tests/test_unstated_ref.py`.

**Amended 2026-10-31, later** (the last readers the previous builder found): `evidence._calls_evidence`
indexed GT into `[N] + alts` without `ref_check`, so `genotype_evidence` and `reference_support` read
`N` → `G`, `0/0`, as hom-ref. Such a row now grades nothing (`GradeBasis.REF_UNSTATED`, unknown), and
still covers its position, so `reference_support` answers `call`, never absence. A symbolic ALT on
such a row is graded as before. `karyotype.measure` counted such a row's GT as an ALT call (`N` →
`G`, `0/1`, a het on X); it now carries no ALT there, its ploidy as written still counted. Tests in
`tests/test_unstated_ref.py`.

## F32 — `write_vcf_zarr` leaves a breakend's mate as held and declares none

Severity: low · Found 2026-10-26 in an audit of F30's repair; until now only the CHANGELOG said so.

F30's repair taught `write_vcf` to name a breakend mate's contig inside ALT as it names CHROM and to
declare every mate contig; `write_vcf_zarr` shares its contig mapping (`writing._contig_names`) but
not the mate handling. A one-row synthetic VCF, `1 100 A[2:500[` with only contig 1 declared, written
with `contig_names="chr"`: `contig_id` is `['chr1']`, and `scan_vcf_zarr` reads the row back as
`chr1`, ALT `A[2:500[`. The mate is on a contig the store's `contig_id` does not list, under a name
the style would not give it, the shape F30 was about. Repairs: write `variant_allele` through
`breakends.mate_renamed` with the same mapping and add the mates (`writing._mate_chroms`, angle
bracketed ones left out) to `contig_id`, as `write_vcf` does; or refuse breakends under `chr` and
`source` until it does. Not fixed here: it is the zarr export's slice, logged so it is not lost.

**Resolved 2026-10-26.** `write_vcf_zarr` takes the first repair: after the rows are read it finds the
written rows' mates with `writing._mate_chroms` (assembly mates and rows `resolved_only` leaves out
excluded), passes them to `writing._contig_names` beside the rows' contigs, and writes
`variant_allele` through `breakends.mate_renamed` with the mapping that gives CHROM its name. The
repro now stores `chr1`, ALT `A[chr2:500[`, and `contig_id` `['chr1', 'chr2']`. Refusing breakends
under `chr` and `source` was rejected: the writer's rule already existed, and sharing its helpers keeps
the two exports from drifting. Tested by writing every source shape and style of
`tests/test_write_breakends.py` to a store and reading it back with polars-bio, `contig_id` from the
bytes (`tests/test_vcf_zarr.py`); the expectations come from the same `written_name` rule. Code:
`vcf_zarr.write_vcf_zarr`; docs/SCHEMA.md, "VCF Zarr export".

## F26 — Loading a WGS gVCF holds gigabytes

Severity: low · Found 2026-10-26 dogfooding RM6's block lift on the PGP-HMS gVCFs.

`load` of hu416394's bcftools gVCF (a 397 MB member of a zip, 9,651,143 calls and 9,473,303 blocks)
peaks at 7.8 GB on a fresh cache (151 s) and at 7.7 GB again with its frames already cached (89 s);
`resolve_ref`, `liftover` and `write_vcf` after it stay within that peak (sampled per phase). The
whole pipeline over hu80855C's Isaac gVCF (1.9 GB, 45,242,155 calls and 138,019,217 blocks) peaked at
46 GB by the end of `load` (1,008 s) and 49 GB on the reload, on a 62 GB machine; that run's own
counting over the frames is inside the figure, so `load`'s share there is not measured. Not looked
into: where in `load` the memory goes (the per-contig spill is meant to bound it, PREVIOUS_ISSUES has
the 6.4-million-call DeepVariant file at 3.6 GB), and why a load whose frames are cached still takes
89 s at the same peak. Blocks nothing yet; a 30x Isaac gVCF is near the edge of this machine.

**Resolved 2026-10-26** for the memory; the cached load's time is explained, not changed. Profiled
phase by phase through the library (RSS sampled every 50 ms, counts and digests only, the same file
from the Krypton share), hu416394 before the repair, fresh / cached: `file_sha256` 4.3 / 0.3 s; the
BGZF copy out of the zip 13 s / reused; `scan_stats`, one polars-bio pass for the warning counts,
40 / 39 s at 0.7 GB; the spill's sink of calls and blocks 53 / 65 s at **8.1 GB**; the spill key's
hash of the staged files 0.1 s; `canonical_parts` 25 s / skipped (generation found) at 6.0 GB, which
was the sink's memory not yet handed back; `split_groups`, `classify` and the variant classes 1 to
2 s each near 5 GB. Sunk alone, the calls plan held 8.0 GB and the blocks plan 2.0 GB, while the raw
records polars-bio hands over (19,124,446 rows) sank at no growth. Bisecting the calls plan by
columns showed memory growing with every column sunk, `chrom` and `pos` alone 0.7 GB, and polars'
physical graph of the plan marked one node as an in-memory fallback: `replace(CONTIG_ALIASES)`, the
canonical `chrom`. polars 1.44 runs `replace` with a mapping as a columnar function, fed through a
multiplexer and zipped back to the other columns, so the multiplexer held every other column until
the function had the whole column. The blocks plan, which also canonicalises `chrom`, held less only
because its rows are narrower. Not polars-bio's: the same plan over a synthetic parquet (30M rows,
three columns) held 2.9 GB against 0.9 GB without the `replace`.

Repair: `schema.canonical_chrom_expr` maps each batch on its own (`map_batches`, elementwise), the
same values, and the plan streams; the post-spill passes (`split_groups`, `callset.spread`,
`callset.variant_class_counts`) collect in the streaming engine. Rows, order, dtypes, facet and
warnings are the same (both frames hashed row by row), and so are the staged bytes: the generation
key is the one the code before wrote, so a cached generation is found again. hu416394 after: 3.9 GB
fresh (the sink 0.8 GB; the peak is now `canonical_parts` sorting one part at a time), 2.1 GB cached
(87.7 s). hu80855C's Isaac gVCF, under a 40 GB cap: 9.4 GB fresh (1,035 s: the copy 173 s,
`scan_stats` 279 s, the sink 434 s, `canonical_parts` 106 s, peak in the passes after the spill,
measured before they streamed) and 6.8 GB cached (804 s, with them; the peak is `callset.spread`'s
unique sites).

Why a cached load takes as long: it reads the VCF through polars-bio twice, as a fresh one does.
`scan_stats` counts the warnings, and the spill's sink writes the rows unsorted only so that
`cache.Spill.key` can hash them and find the generation an earlier load committed; the staged files
are then thrown away. Nothing re-reads the source for a content key or a facet field: `file_sha256`
of the 397 MB zip is 0.3 s warm, and the copy is reused. The two passes are 38 + 45 s of hu416394's 88
s and 371 + 406 s of the Isaac file's 804 s. Shortening it is a design call, not taken here: keying a
load's generation by its inputs (the source digest, the options and a loader version bumped with
every change to the rows) would make a cached load a hash and a lookup, which is what `_spill_files`'s
docstring rejected, since a loader fix that forgets the bump would hand back the rows an older loader
wrote. A narrower option keeps the key and caches `scan_stats`'s counts by the copy's digest, saving
the first pass only. Code: `schema.canonical_chrom_expr`, `loaders/vcf.split_groups`,
`callset.spread`, `callset.variant_class_counts`; test in `tests/test_spill_streams.py` (no spilled
plan of the gVCF, CGI, key=value or Nirvana loaders has an in-memory fallback node; failed on the
gVCF, key=value and Nirvana plans before).

**Amended 2026-10-26:** the first pass is now cached, the narrower option above, as the lead chose
(docs/LAYERS.md). `scan_stats`'s result is stored as `stats/<key>/stats.json`, the key a sha256 of
every input of the pass (`loaders.vcf.stats_key`: the bytes polars-bio reads and the file's name, the
source's digest and member, the request and `SCAN_OPTIONS`, the header's version and left-out keys,
the retry's `Header.unparsed`, the part, `READER_DEFINITIONS`, `STATS_VERSION`, the polars and
polars-bio versions, and a digest of `loaders/vcf.py`'s and `loaders/text.py`'s source); the
generation is still keyed by its staged bytes. hu416394, frames cached, same conditions as above
(Krypton, library only, 40 GB cap): 85.0 / 86.7 s before, 47.4 / 48.4 s after, peak 2.0 GiB both; the
copy's digest costs 0.3 s per load. The first load on this version runs the pass and stores it
(88.9 s). Warnings, facet row and row counts hash the same. The Isaac gVCF was not measured again;
its pass was 371 s of 804. A cached load still sinks the rows to find its generation (45 s here, the
rest of the time). Code: `loaders/vcf.read_part`, `stats_key`, `cached_stats`; tests in
`tests/test_stats_cache.py` (a reload runs no pass and hands `to_rows` the same values and types;
each key input changed reruns it; a truncated, garbage or other-key file is counted again; a refused
read is not stored; a second loader waits on the lock and reads the first one's result; failed on main
but for the copy-bytes test, which also fails with the copy's digest left out of the key).

## F24 — Every source contig costs a fixed pass, so a header full of decoys is slow and heavy

Severity: low · Found 2026-10-26 dogfooding RM6's SV lift on hu416394's Manta file (Sequencing.com,
GRCh38, a 2,581-contig hs38DH-style header).

`resolve_ref` and `liftover` work contig by contig, each contig a filtered scan, a few `collect`s and a
parquet file in the spill. The file has 9,373 rows on 942 contigs, 918 of them decoys and scaffolds
with a handful of rows each (1,228 rows on decoys, which no assembly holds and which stay unread, as
they should). Profiled on a fresh cache: `resolve_ref` 13.4 s (943 parquet sinks), `liftover` 27.2 s
(942 `_lift_contig` passes, 13,217 `collect`s, 997 sinks, `_segments` rebuilt 72 times). The same
rows cut to the 25 primary contigs resolve in 0.4 s at a 1,774 MB peak; all of them take 14.3 s and
3,757 MB. For scale, hu1C1368's GRCh37 Manta file (7,246 rows on 24 contigs) resolves in 0.5 s and
lifts in 2.7 s under 430 MB. Nothing is wrong in the output; the cost is per contig, not per row.
Repairs, each a change to the spill design: batch the small contigs into one pass; send rows on a
contig the store or the chain does not hold straight to their reason without a pass of their own.
Neither is free: the per-contig parts are what keep a WGS gVCF within memory, and their layout is a
cache key (`SPILL_LAYOUT`).

*Updated 2026-10-26:* the first repair is in for `load`'s own spill (F27, PREVIOUS_ISSUES.md): it is
written in parts of whole contigs, consecutive in canonical order, up to `SPILL_PART_ROWS` rows, so
contigs with a few rows share a file (`loading.canonical_parts`). hu416394's bcftools gVCF, rows on
2,179 contigs, is 12 call parts and 11 block parts. `resolve_ref` and `liftover` still make a pass and
write a file per source contig, which is what this finding measured; it stays open for them.

**Resolved 2026-10-26** by both repairs, for `resolve_ref`, `normalize` and `liftover`. Each reads its
staged input a group of whole contigs at a time, consecutive in canonical order up to `SPILL_PART_ROWS`
rows unless one contig alone has more (`loading.canonical_batches`, the rule `canonical_parts` already
wrote `load`'s spill by, now a function both use): one filtered scan and one part file per group, where
it was one per contig with rows. A lift's lifted `calls` and `blocks` are then written in parts of whole
target contigs through `canonical_parts`, and its report in parts of whole source contigs, in the
parent's order within each. Within a group, the rows on contigs nothing can be done for go to their
answer together: `resolve_ref` passes the rows of every contig it cannot read through as they are;
`liftover` gives the rows of every contig the chain does not cover, and of every contig it covers at
another size, their reason in one call each; `normalize` splits, deduplicates and sorts the rows of
every contig the collection does not hold in one call (its conflict count now groups by contig, which
the per-contig call did without saying). A contig whose sequence is read, or that the chain covers, is
still worked on its own, against its one sequence or through its own chain blocks; a lift places the
breakend mates of the whole group once per mate contig (`liftover._mate_places`), where it placed them
once per pair of a source contig and a mate contig. Memory is bounded by the group as it was by the
contig.

Measured on the two PGP-HMS files of this finding (library only; counts, timings and digests; the
code before and after over one scratch cache per file, the load spill shared, each step's spill written
fresh; 40 GB cap; the 62 GB machine shared with other sessions, load average 6 to 16, so times are
rough). Every output hashes the same before and after, row for row and dtype for dtype: the resolved,
lifted and normalized calls, the lifted blocks, the lift report, the reasons and block reasons, the
warnings and the facets' `contig_map`.

- hu416394's Manta file (9,373 rows on 942 contigs; `hg38ToHg19`): `resolve_ref` 4.6–4.9 s with 942
  part files → 0.3 s with 1; `liftover` 17.0–21.9 s with 942 report files and 30 lifted ones → 3.0–3.3
  s with one each (`_place` called 462 times → 120, `collect`s 15,040 → about 5,200; most of what is
  left is `_segments`, the chain's axis built once per chain contig a row or a mate lies on, 96 of
  them); `normalize` 12.3–32.8 s with 942 files → 0.4 s with 1. Peaks within 100 MB of the load's.
- hu416394's bcftools gVCF (9,651,143 calls on 2,057 contigs, 9,473,303 blocks; `hg38ToHg19`), in
  one process after a cached `load` (7.7–7.9 GB): `resolve_ref` 110.5 s → 27.8–29.7 s, 2,057 part files
  → 12, 4.6 GB above the start → 4.7–5.0 GB; `liftover` 226.6 s → 65.7 s, 2,057 report files → 12,
  68 lifted call files → 10, 70 block files → 10, 5.9 GB above the start → 4.3 GB; `normalize`
  206.2 s → 71.9–75.7 s, 2,057 files → 12, 3.1 GB → 3.1–3.4 GB. The run's peak 8.7 GB → 8.3 GB.
  `liftover` alone in a fresh process, twice each: 169–171 s at 4.8 and 6.1 GB above the start before,
  64–68 s at 5.7 GB both times after, so the lift's memory is the old one's within the noise; a group of
  up to a million rows holds about what the largest contig (chr1) did.

Each step's spill names its layout in its versions, as `load`'s and `rejoin_split`'s do (F27):
`RESOLVE_LAYOUT`, `NORMALIZE_LAYOUT` and `LIFT_LAYOUT`, each 1. A generation an older version wrote, a
file per contig, is keyed without the key, so no run of this version finds it; the sweep retires its
entry once unused past the grace period. Code: `loading.canonical_batches`, `reference.resolve_ref`,
`normalize.normalize`, `liftover._lift_into`, `_lift_group`, `_mate_places`, `_lift_blocks_into`,
`_target_parts`; tests in `tests/test_step_spill_order.py` (every output against the steps run a contig
per group, at three bounds; the parts and the plans; an older layout's generation not read, and swept).

Rejected: sorting the input into parts with `canonical_parts` and working part by part, which reuses the
writer as it is but re-sorts the input, so the lift report would leave the parent's order wherever the
parent's frame was not in canonical order already, and adds a sort pass over the whole input; one pass
over every unreadable or unchained row of the file at once, which is not bounded by the group;
`_normalize_contig` over a whole group, which would align many contigs against one sequence;
bumping `RESOLVE_VERSION`, `NORMALIZE_VERSION` or `LIFT_VERSION`, which `facets.resolve_version` and
the derived `facet_id`s carry, so the facets would be renamed for the same rows; and bumping
`cache.SPILL_LAYOUT`, which would retire the VRS and coverage spills too. Left open: `mint_vrs_ids`
already writes its ids in parts of `_IDS_PER_PART` rows, but still collects a slice per contig run
(not measured here); `_segments` is a Python loop per chain contig.

*Amended 2026-10-26:* `mint_vrs_ids` measured and repaired the same way. It wrote its ids in few parts
already, but made a slice-`collect` and a `_mint_rows` pass per run of one contig, about twenty
`collect`s each, so each contig with rows paid a fixed pass. It now reads a group of whole runs at a
time, consecutive in row order, up to `SPILL_PART_ROWS` rows unless one run alone has more
(`loading.bounded_groups`, the rule `canonical_batches` now groups by too; on `resolve_ref`'s output,
which is in canonical order, the groups are `canonical_batches`'): one slice, one pass and one ids part
per group. A pass places each allele on its contig's sequence through an integer column naming the
group's contigs (`_contig`, the accession joined in only on the distinct locations), withholds
the alleles of a contig with no single sequence as `no_contig` in the same pass, and justifies each
contig's indels together, its sequence read once per group and let go before the next contig's (it was
read once per run). Measured on the two files above (library only, after a cached `load` and
`resolve_ref`, the VRS spill written fresh, counts, timings and digests; 12 and 24 GB caps; load
average 3 to 8, so times are rough). Every output hashes the same before and after, on every run: the
minted calls (ids, `vrs_errors`, every other column, row order, dtypes), `vrs_withheld` and the
warnings.

- hu416394's Manta file (9,373 rows on 942 contigs): 18.6–18.8 s with 942 passes and 1,886 `collect`
  calls (about 20,000 inside polars) → 11.4 s, three runs, with 1 pass and 4 `collect`s; one ids part
  before and after; the process's peak 840 MB → 850 MB. What is left is 24 whole-contig fetches for
  4,225 indels: 7.2 s in `get_substrings` and 3.7 s in `fetch`'s `isupper` guard (profiled).
- hu416394's bcftools gVCF (9,651,143 calls on 2,057 contigs, 876,559 indel alleles justified):
  97.5–108.8 s with 2,057 passes, 4,116 `collect`s and 8 ids parts → 70.5–75.0 s with 12 passes, 26
  `collect`s and 12 parts (polars' own collect time 21.2 s → 2.0 s, profiled). Peak RSS sampled during
  the mint 4.72–4.75 GB before, 4.71–4.72 GB after (3.0–3.5 GB above the start before, 3.9–4.0 GB
  after, from a lower start); a group of up to a million rows holds about what the largest contig did.
  What is left is per allele: 24.4 million `sha512t24u_digest` calls (18.5 s profiled), the Python
  justification of each indel (about 17 s profiled, gtars' roll 0.8 s of it), and 25 whole-contig
  fetches (11 s).

`VRS_LAYOUT` (1) is in the spill's versions and generation key: an old generation would read the same
(its ids are row-aligned whatever their parts), but is written once more under the key, so the manifest
says what it holds and the sweep retires the old entry. `MINT_VERSION` stays 3 (it names the cache
entry only, no derived id). Code: `vrs.mint_vrs_ids`, `_mint_rows`, `_identify`, `_location_json`,
`loading.bounded_groups`; tests in `tests/test_step_spill_order.py` (minted ids, reasons and warnings
against the per-run oracle at three bounds, on the canonical frame and on one interleaved by position;
a contig with no sequence beside stored ones; one pass and one part per group; each contig's sequence
read once per group; an old-layout generation not read, and swept; the last four fail on main).
Rejected: filtering the staged input by `canonical_batches`, which would hand ids back out of row
order for an input not in canonical order and need a sort to put them beside the rows; keeping
`_IDS_PER_PART` apart from the group bound, two rules for one shape. Left open, both measured above
and neither F24's pattern: a whole contig is fetched to justify its first indel, about 0.3 s for
chromosome 1 in `get_substrings`, so on the Manta file 24 fetches for 4,225 indels are nearly all that
remains (a window around each indel would need a widen-and-retry rule to be shown output-preserving,
as gtars rolls over an unbounded region); and one `sha512t24u_digest` call per location and per
allele, about 24 million on the gVCF, a Python loop over gtars that only a vectorised digest in gtars
would shorten. The `isupper` guard, 0.15 s more per chromosome 1, went afterwards: the store settles
case where a sequence enters it (a collection `add_fasta` added is uppercase from gtars' FASTA import,
recorded by its `just_vcf_fasta` alias), so `fetch` scans only a sequence with no such record. The
Manta file then mints in 7.45–7.55 s (11.27–11.29 s with the guard), every output hashing the same.

## F30 — `write_vcf`'s contig names do not reach a breakend's mate inside ALT

Severity: low · Found 2026-10-26 dogfooding RM6's breakend lift on hu1C1368's Manta file.

`write_vcf(..., contig_names="chr")` writes CHROM as `chr1` but leaves each breakend ALT's mate as
the frame holds it, here `2` (both store collections name contigs without `chr`): all 1,216 lifted
breakend ALTs, and all 1,228 of the loaded facet's, name a mate contig the written header does not
declare. With the default `stored` names, one lifted ALT of 1,216 names a contig no row stands on,
which the header leaves out too. A tool reading the mate against the header (or joining by name)
misses it; `load` and `concordance` read it back, since they canonicalise the mate. Not new with the
lift: a loaded facet whose source spelled `chr2` inside ALT writes it so under `stored`. Repairs:
rewrite the mate contig inside a breakend ALT by the same mode as CHROM, and declare every contig an
ALT names in the header (a writer change, with its own test of the round trip); or have the lift
write the mate by its canonical name, which only moves the mismatch to the `source` and `chr` modes
and gives up "named as the target collection names it". Surfaced, not fixed: it is the writer's rule,
outside the lift's slice.

**Resolved 2026-10-26** by the first repair. The mate's contig inside each of the four bracket forms
is read canonically (`breakends.breakend_parts`, as `load` reads CHROM) and written through the one
canonical → written mapping CHROM is written through (`breakends.mate_renamed`), and every contig a
written ALT names joins the contigs the header declares, rows or none (VCF 4.3 §1.4.7). A mate the
style has no other name for keeps its canonical name, as such a CHROM does. Single breakends, symbolic
alleles and brackets that are none of the four forms are written as held. Under `stored` a source's
`chr2` inside ALT is now written `2`; under `source` a lifted facet's ALTs come back verbatim, its
`contig_map` holding the target's spellings. Two canonical names the style writes as one get one
`##contig` line: the lifted hu1C1368 facet stands rows on `chrUn_KI270742v1` (the chain's name) and
names their mates `KI270742.1` (the FASTA's), both `KI270742.1` under `source`; `vcf_zarr` indexes
`contig_id` by written name for the same reason. hu1C1368's Manta file written and re-loaded, mates the
header declares / does not (counts only): loaded facet, `stored` 1,228 / 0 → 1,228 / 0, `chr` 0 / 1,228
→ 1,228 / 0, `source` 1,228 / 0 → 1,228 / 0; lifted to GRCh38, `stored` 1,215 / 1 → 1,216 / 0, `chr`
0 / 1,216 → 1,216 / 0, `source` 1,216 / 0 → 1,216 / 0. No version constant: `write_vcf` has none, its
output is keyed into no cache, and no frame column changes. The bracket search costs 0.1 s on a 4.8
million-row VCF without breakends; the mate read runs only on a facet that has one. Code:
`writing._mate_chroms`, `writing._line`, `writing._contig_names`; tests in
`tests/test_write_breakends.py` and `tests/test_liftover_bnd.py`
(`test_a_lifted_facet_written_with_chr_names_declares_every_mate`).

**Amended 2026-10-26** (an audit of the repair). "Every contig a written ALT names" took in two it
should not have. A mate in angle brackets (VCF 4.3 §5.4.2's `]<ctg1>:7]T`) names a contig of the
`##assembly` file, not of the reference; it was declared `##contig=<ID=<ctg1>>` under every style,
which htslib warns is an invalid contig name. It is now written as held and declared nowhere, by
rule (`breakends.assembly_mate`) rather than because no mapping named it. And a row
`resolved_only=True` left out still declared its mate; now only the written rows' mates are
declared. Tests: `test_an_assembly_mate_is_written_verbatim_and_declared_nowhere`,
`test_htslib_reads_an_assembly_mates_header_without_a_warning` (pysam, skipped without the
`[reference]` extra) and `test_a_skipped_rows_mate_is_not_declared`.

**Amended 2026-10-26** (the same audit; the lead's call). "Under `source` a lifted facet's ALTs come
back verbatim" held because a lifted facet's `contig_map` lists every mate contig. A loaded facet's
lists only the contigs with rows, so "a mate the style has no other name for keeps its canonical
name" made `source` rewrite a source's own spelling: a chr-named VCF with no `##contig` lines wrote
`chr1 100 A[chr2:500[` as `A[2:500[` with `##contig=<ID=2>`, and one declaring chr1 and chr2 wrote a
mate on chr3 as `A[3:500[`. Under `source` a mate is now written as the source spelled the contig:
the facet's spelling (`contig_map`, then `header_contigs`) where it has one, so an ALT's `chr2`
beside CHROM `2` still writes `2` and the file agrees with its own CHROM column; otherwise the ALT's
own spelling, declared under it. A contig only ALTs name, spelled two ways by them, is refused, as a
contig CHROM or the header spells two ways already was. `stored` and `chr` are unchanged. Code:
`writing._contig_names` (`mates`), `writing._mate_chroms` (now each canonical with its spelling).
Tests: `test_source_writes_a_mate_the_file_alone_names_as_it_spelled_it`,
`test_source_prefers_the_chrom_spelling_of_a_mate_with_rows`,
`test_source_refuses_a_contig_only_alts_name_spelled_two_ways`, and the chr-named source with its
header declaring only the rows' contigs and with none, under every style.

## F29 — Reading blocks over a facet with hundreds of contigs exhausts memory

Severity: medium · Found 2026-10-26 dogfooding RM282's breakends on the PGP-HMS Manta files.

`concordance` of hu416394's Manta file (GRCh38, its calls on 942 contigs: alts, decoys, HLA) against
itself with the records at odd positions dropped grew past 36 GB resident before it was stopped; the
same probe on hu1C1368's file (24 contigs) finishes in 3 s. The memory goes to `_block_reading`, which
reads the other facet's blocks under each one-sided sequence allele through `intervals.covering`, and
the file has no blocks at all. `covering` builds one range-join branch per contig of its probes, each
over the facet's lazily scanned blocks, and that facet's spill is 942 per-contig parquet files. Over
hu416394's calls as probes, `covering` takes 0.5 s for one contig and fails to allocate under a 30 GB
address-space cap for all 942. On hu1C1368's (24 contigs, 24 files) it takes 2.7 s for all of them,
0.24 s for 8 and 0.01 s with the blocks as an empty in-memory frame: the cost grows with branches
times files. Manta against Canvas (a 22-file spill) ran; a file against itself never reads blocks,
since every sequence allele meets its twin.

Surfaced, not fixed: `covering` is shared with `genotype_evidence_at` and `reference_support`, and
its plan was chosen for gVCF blocks. A candidate repair is to collect the intervals that reach their
successor once and run the per-contig range join only on contigs that have any; whether collecting is
safe for a lifted 30x gVCF's overlapping blocks (hundreds of thousands) is the measurement it needs
first. Not measured: whether a GRCh38 gVCF with alt contigs (its blocks spilled per contig too) meets
the same wall in `genotype_evidence_at`.

**Resolved 2026-10-26** by two changes on main that this branch had not yet met: `intervals.covering`
collects the as-of join and the overlapping blocks in one `pl.collect_all` and range-joins in memory
(merge ccf4829), and `load` writes its spill in canonical order in parts of whole contigs, so the 942
per-contig files became a few parts read in one scan (F27, merge c4def55). The same probe, hu416394's
Manta file against itself with the odd positions dropped, now runs end to end in 6 s under a 30 GB
cap (1,248 `match`, 427 `only_a`, as expected). The candidate repair above is what ccf4829 did; the
alt-contig gVCF case it left unmeasured is covered by RM3's measurement on hu416394's gVCF.
Code: `src/just_vcf/intervals.py` (`covering`), `src/just_vcf/loading.py` (`_canonical_parts`).

## F27 — Every read of a loaded VCF facet's frames pays one scan and sort per source contig

Severity: medium · Found 2026-10-26 measuring `reference_support` and `genotype_evidence_at` at
whole-genome scale (RM3) on hu416394's bcftools gVCF (Sequencing.com, GRCh38, 9,651,143 calls and
9,473,303 blocks; rows on 2,179 of the header's contigs).

`load` hands back `calls` and `blocks` as one filtered, sorted scan of the spill per contig that has
rows (`canonical_order_per_contig`), so a facet from a decoy-heavy header is a union of 2,179 parts,
and every query re-runs it. Measured with the library, one process, after `load`: the calls filtered
to the 24 primary contigs and to the rows a lookup can use (866,279 spanning calls plus those at the
asked positions) take 83 s at about 3 GB above the post-load baseline, however few rows match; the
blocks' as-of join 22 s; `select(pl.len())` over both frames 16 s at 5 GB. The same lookups over those
frames collected into memory take 1.2 s for 1,000 sites and 2.9 s for a million, and over the lifted
facet's frames (24 parts in its spill) 1.6 s and 4.1 s. A filter on a million positions costs 128 s
more on the loaded facet (211 s against 83 s), consistent with the position set being evaluated in
each part's scan; on the lifted facet it costs nothing visible.

Joins over these frames also held far more than the rows they keep: before `_covering` read each
frame once, 1,000 sites ran past 30 GB (before 6750e05, also with the sites on chr21 and chr22 only)
and past 55 GB after it, where the same queries over the frames held in memory peak at 10 GB. `_covering` no longer joins the
lazy frames, so the lookups are within bounds; the per-query read cost stays.

The part count is what costs: hu80855C's Isaac gVCF (25 parts, 45,242,155 calls and 138,019,217
blocks, ten times the rows) answers 1,000 sites in 20 s and a million in 24 s (RM3). Not looked into:
what in each part's scan costs. Repairs, each a change
to `load`'s spill: write the spill in canonical order once, so the frames are one scan with no sort
at read time (the reason it is not is the 11 GB global sort `loading.py` notes; a per-contig sort at
load time, written in order, would keep that bound); or batch the small contigs into one part, F24's
first repair. Both change `SPILL_LAYOUT`. Blocks nothing; a panel query on such a facet takes
minutes instead of seconds.

Resolved 2026-10-26 by the lead's decision: the spill is written in canonical order once, at load time,
sorted a part at a time so load's memory stays bound, and the frames come back as one scan with no
sort. `_spill_files` still sinks each frame unsorted in one streaming pass, then `canonical_parts`
writes it again as parts of whole contigs, consecutive in canonical contig order, up to
`SPILL_PART_ROWS` rows unless one contig alone has more, each part filtered from the staged file,
sorted (`canonical_order`) and written; the unsorted file is removed. The frames are
`pl.scan_parquet(parts)`. Contigs with a few rows share a part, F24's first repair for `load`'s spill.
A load spill's manifest and generation key carry `FRAMES_LAYOUT` 1, so a generation an older version
wrote is never found, and the sweep retires an older entry; `cache.SPILL_LAYOUT` stays, since the
derived steps' spills, which carry it too, did not change. Every loader that spills goes through it
(VCF, Nirvana, key=value, CGI).

Measured on hu416394 (RM3, "after F27's repair"): one scan of 12 call parts and 11 block parts, no
sort; `select(pl.len())` 15.7 s → 0.3 s; `reference_support` at 1M sites 249 s → 4.0 s; both frames
and both answer sets hash the same row for row; `load` peaks at 7.7–7.9 GB as before, its sort pass
15 s. Code: `loading.py` `_spill_files`, `canonical_parts`; tests: `tests/test_spill_order.py` (the
frames against the old per-contig read of the same staged bytes, the plan with no sort node, the
parts holding whole contigs, an older layout's generation not read and swept).

Rejected: one sorted file per frame, which needs the whole frame sorted in memory (the 11 GB global
sort) or a parquet writer that appends a part at a time (pyarrow's, not a core dependency), when the
parts already read as one scan; and bumping `cache.SPILL_LAYOUT`, which would mark every resolve,
lift, VRS, normalize and coverage spill older too, so the sweep would retire those unused for the grace
period and their next run recompute them, though their layout did not change. Left at first:
`rejoin_split`'s own spill (`rejoin.py`) was still read back one contig at a time with
`canonical_order_per_contig`, so a rejoined whole-genome facet still paid a sort per contig on each
read.
## F28 — Manta's breakend mates do not state each other's positions exactly

Severity: medium · Found 2026-10-26 dogfooding RM282's breakends on the PGP-HMS Manta files.

Folding mates by the exact ends they state left hu416394's 3,350 breakend records as 1,899
adjacencies, not 1,675, and hu1C1368's 1,228 as 785, not 614: 448 and 342 records stood alone, each
mate counted as a second junction. Every lone record had a partner a few bases off (most under 10)
pointing back at it, joined the same way. Manta places a breakend's POS at one end of the junction's
microhomology, which CIPOS spans, while the mate's ALT points at the other end. The frame carries no
MATEID (`dropped_fields`), so the file's own pairing is not there to read.

Resolved 2026-10-26: two records whose mate contig and side are each other's own, and whose mate
position lies within the other's `[POS + CIPOS[0], POS + CIPOS[1]]`, are one adjacency, where the
pairing is one-to-one (`concordance._fold_mates`). Every one of the 224 and 171 pairs lay inside the
intervals, none ambiguous, and the counts are now 1,675 and 614, half the records each.
Rejected: pairing by `ids` (`MantaBND:…:0` beside `:1`), which is Manta's naming and no other
caller's; and carrying MATEID into `calls`, a schema change for the person, which these files do not
need. Test: `tests/test_concordance_bnd.py::test_mates_written_apart_within_their_intervals_count_once`.

Closed 2026-10-26 for `rejoin_split` too: its spill goes through the same `loading.canonical_parts`
(renamed from `_canonical_parts`, now shared), and its manifest and generation key carry a
rejoin-only `REJOIN_LAYOUT` 1, so an older rejoin generation is never read and no other spill is
retired. `REJOIN_VERSION` stays at 5: the rows are the same, and it is part of the derived
`facet_id`, so a bump would rename the facet and orphan the resolve, normalize and VRS spills made
from it. The old read took its contigs from the parent's `contig_map`, the new one from the rows;
every loader and `liftover` build `contig_map` from the contigs with rows, so the two agree. With no
caller left in `src/`, `canonical_order_per_contig` moved into `tests/test_spill_order.py` as the
oracle of the old read. On hu416394 (RM6, dated 2026-10-26): the rejoined calls' plan 2,179 scans and
sorts → one scan, no sort; `select(pl.len())` 17.2 s → under 0.01 s; collected whole 80.9 s → 0.75 s;
the same rows by hash; `rejoin_split` 9.9 s → 20.4 s for its one sort pass. Tests:
`tests/test_rejoin_spill_order.py`. No other spill sorts on read: `resolve_ref`, `normalize` and
`liftover` sort each contig when they write it, `coverage` sorts before it writes, and `mint_vrs_ids`
keeps its input's order.

## F22 — A haploid call against the same allele's diploid homozygous call counts as discordant

Severity: medium · Found 2026-10-23 by RM279 on the PGP-HMS CGI against freeBayes pairs.

Outside the PARs CGI writes a male's X call haploid (`1`) and freeBayes diploid (`1/1`). Every such
pair is `ploidy_differs`, which the atlas counts in `discordant`: 72,462 and 64,818 SNVs, 65% and
67% of the two pairs' discordant SNVs. Both say the person carries the allele on every copy of X
they have; the difference is the writer's ploidy convention, not a method error, and it swamps the
rate the atlas exists to give.

Resolved 2026-10-23 by the lead's rule: a haploid call against the same allele's diploid homozygote
is a `match`, and the ploidy is a flag of its own, `sites.ploidy_differs` (True, False, or null on a
one-sided allele), counted apart in `summary` (a group key) and in the atlas (`ploidy_differs`, the
flagged alleles of `called_both`). Of two ploidies the genotypes agree when each holds the allele on
every haplotype or each on none; a haploid `1` against `0/1` stays a `mismatch`, and with a null
haplotype a pair is decided only where every completion agrees (`1` against `1/.` is `uncalled`,
against `0/.` a `mismatch`). The atlas's `called_both` is `match + mismatch` and `discordant` is
`mismatch`, so the convention leaves the discordance rate. The graded comparison compares copies as a
share of the ploidy (haploid `k` against diploid `2k`), as the hard one does. `Agreement` loses
`ploidy_differs`.
Measured on the rebuilt PGP pairs (`x_nonpar`, `carried_pass`): of 71,715 and 64,382 flagged SNVs,
71,193 and 61,920 match and 522 and 2,462 are haploid calls against a freeBayes het, which stay
discordant. The pairs' SNV discordance falls from 3.40% to 0.54% and from 3.07% to 1.02%, with F21.

Rejected: keeping `ploidy_differs` as an agreement value and leaving it out of the rate (a haploid
call against a het is a real disagreement, and the value hid it among the convention's matches);
rewriting the haploid call as `1/1` (the source wrote one allele, Decision 2: list length is never
corrected). Code: `concordance._agreement`, `_ploidy_differs`, `_distribution` (`_dose`),
`atlas._tally`, `pool_counts`; tests `test_genotypes_compare_as_copies_of_the_allele`,
`test_a_haploid_call_is_graded_against_a_diploid_one_by_the_share_of_copies`
(tests/test_concordance.py), `test_counts_follow_the_cases_per_region_and_scope` (tests/test_atlas.py).

## F21 — A VRS-met allele is read as zero copies in a complex multi-allelic row

Severity: medium · Found 2026-10-23 by RM279 rebuilding the PGP-HMS CGI against freeBayes pairs with
alleles met by VRS id.

freeBayes writes a haplotype as one multi-allelic row: REF `CGCCC`, ALT `TGCCT,TGCCC`, GT
`TGCCT/TGCCT`. Its ALT `TGCCC` trims to `C>T`, so its VRS id is CGI's `C>T` (GT `T/T`) and the two
meet (`key = vrs`). `concordance` then counts the copies of the spelling `TGCCC` in the freeBayes
genotype, zero, although each `TGCCT` haplotype carries the same `C>T`: a `mismatch` that is only
representation. 20,434 to 23,014 SNV mismatches on the CGA 2.5.0.44 pair, 2,601 to 2,950 on the
2.5.0.24 one, up to 5,982 and 4,627 indels (docs/METHOD_ATLAS.md). The `site` path guards the same
case by reading only a row over the same REF span; the `vrs` path has no guard. The graded
comparison has the same fault: its distribution counts the copies of the spelling too.

Resolved 2026-10-23: each side's genotype is read for the allele's change in its own row. A haplotype
that is neither REF nor the met ALT is trimmed against REF like the ALT: a change away from the
allele's span is not a copy; one of the same length is a copy when its bases over the span are the
ALT's (`TGCCT` holds `C>T`, `CGCCT` does not); a length change against another length change or a
substitution is another allele of the site (a `1/2` deletion and insertion of one repeat unit), not a
copy. What is left, a change of another length that also replaces bases over the span, cannot be
read: it is open like a null haplotype, never zero copies, and a pair it leaves undecided gets the new
`agreement` value `unreadable`, counted in a warning and in the atlas's `unreadable` column (never in
`called_both`). The `site` pairing and the graded distribution read alleles the same way, and the
atlas counts an allele as carried where `concordance` read a copy of it under another spelling.
Measured on the rebuilt PGP pairs (region `all`, `carried_pass`, docs/METHOD_ATLAS.md): SNV
mismatches outside X's non-PAR fall by 21,873 in the CGA 2.5.0.44 pair (20,434 counted exactly before)
and by 2,460 in the 2.5.0.24 pair (2,601 counted); 2,315 and 201 SNVs are `unreadable`; mean `p_same`
over graded MNVs rises from 0.661 to 0.706 and from 0.776 to 0.808.

Rejected: marking every VRS pair whose row holds a sibling ALT unreadable (it drops the 20,000 SNVs
the same-length rule reads exactly, and the common multi-allelic indel row); reading only the
`vrs` pairing (the `site` and `position` pairings read the same genotypes, and the rule gives the old
answer wherever a genotype holds only REF and the ALT). Code: `concordance._holds`, `_codes`, `_read`,
`Agreement.UNREADABLE`, `atlas._counts`; tests
`test_a_vrs_met_allele_is_read_in_terms_of_its_own_change_in_a_complex_row`,
`test_a_multi_allelic_indel_row_reads_its_other_length_as_another_allele` (tests/test_concordance.py),
`test_a_complex_row_s_sibling_carries_the_allele_and_an_unreadable_one_is_counted_apart`
(tests/test_atlas.py).

## F20 — A CGI vcfBeta's chrY md5 refuses the whole facet

Severity: medium · Found 2026-10-23 by RM279 (the method-error atlas) comparing PGP-HMS CGI vcfBeta
facets against freeBayes.

`resolve_ref` refused a facet when any primary contig (1–22, X, Y) stated a length or md5 the
reference does not have. Checked with the library on the smallest PGP-HMS vcfBeta against the shared
store: the header states 25 contigs (1–22, X, Y, M), every autosome, X and M agree with Ensembl
GRCh37 in length and md5, and Y agrees in length (59,373,566) but not in md5. So the whole facet was
refused: no `ref`, no `ref_check`, no VRS ids, and the atlas compared the CGI pairs by label only
(docs/METHOD_ATLAS.md). Assemblies of one build legitimately differ on X, Y and MT (PAR masking in an
analysis set, rCRS against the Yoruba mitochondrion); an autosome of another sequence is another
assembly.

Resolved 2026-10-23 by the lead's rule: only an autosome whose stated length or md5 differs refuses
the facet. A differing X, Y or MT leaves that contig's rows unread (`ref` null, `ref_check`
`not_run`), and every matching contig is read; the warning names each mismatching contig with its
row count and the stated against the reference length or md5. A facet with nothing else to read is
still refused. Measured with the library on the same vcfBeta: 4,166,235 rows, 4,163,187 read (99.93%),
every one `match`; the 3,048 Y rows stay `not_run`, X (90,020) and MT (36) are read in full, and
`mint_vrs_ids` then mints 8,316,880 allele ids. `RESOLVE_VERSION` 9.

Rejected: dropping Y from the refusal only (the lead's rule names X, Y and MT together, as the
contigs where assemblies of one build vary; no X mismatch has been seen on real data yet); matching Y by length alone when the md5 differs (the md5 is the evidence that the
bases differ, and reading them would check CGI's REF against a sequence it did not use); changing
`PRIMARY`, which `concordance` and `liftover` read for their own rules. Code: `reference.AUTOSOMES`,
`reference._plan`, `reference._contig_warnings`; tests in `tests/test_reference.py`
(`test_a_sex_chromosome_that_differs_is_left_unread_when_the_autosomes_match` and the two beside it).

## F19 — UCSC scaffold names do not meet Ensembl's

Severity: medium · Found twice on 2026-10-22/23: by RM3 (`resolve_ref` on PGP-HMS GRCh38 VCFs) and by
RM286 (lifting GIAB HG002 Q100 from CHM13 to GRCh38 with UCSC's chain).

UCSC names GRCh38's unplaced, unlocalised and patch scaffolds `chr1_KI270706v1_random`, `chrUn_…`,
`…_fix`, `…_alt`; Ensembl's FASTA (the reference store) names the same sequences by accession
(`KI270706.1`). `canonical_chrom` maps the primary chromosomes and M/MT, not these. So 44,594 rows (1%)
of a PGP GRCh38 VCF stay unread by `resolve_ref`, and 75,001 lifted HG002 rows land as
`target_contig_unmatched`. Not fixed: the repair is an alias table from the GRCh38 assembly report
(UCSC name ↔ GenBank/RefSeq accession, the same source RM286 used for CHM13), checked by length, never a
name-pattern guess.

Resolved 2026-10-23 by the alias table F19 asked for. NCBI's GRCh37.p13 and GRCh38.p14 assembly
reports are carried in `assets/` and derived into the package table `contig_names.tsv`
(`just_vcf.contig_names`). A reference contig that a report lists as a scaffold, at the report's
length, is also served under that scaffold's other report names, unless the collection holds a
contig of that name itself; the stated length or the chain's size is then checked against it as
before. The row keeps its source name in `chrom`. Measured with the library on the same files:

- the PGP-HMS GRCh38 VCF against the Ensembl store: 44,594 rows unread before, 0 after; all
  4,651,432 rows read, every one matching its source REF, 145 scaffolds read under their accession;
- GIAB HG002 Q100 lifted CHM13 → GRCh38 through `hs1ToHg38` onto store-chm13:
  `target_contig_unmatched` 75,001 before, 40,099 after; 33,088 more rows lifted (28,690 onto
  `_random` scaffolds, 4,398 onto `chrUn_`), 1,814 more `ref_mismatch_after_lift`. Every chain block
  onto a `_random` or `chrUn_` scaffold (22,760 of 34,178) now meets the store; the 11,418 onto `_fix`
  patches do not, and the 40,099 are theirs: the Ensembl primary assembly holds no patch, so that is
  the store, not a naming gap.

Rejected: putting the GenBank accession in `chrom` for a scaffold (Decision 7 keeps the source's
spelling for every non-primary contig, and it would change the rows of every UCSC-named file on load);
reading the accession out of the UCSC name (`chr1_KI270706v1_random` → `KI270706.1`), a name pattern
that GRCh37's lower-case names (`chr4_gl000194_random`) and its haplotypes (`chr6_apd_hap1`) already
break. Code: `contig_names.scaffold_names`, `ReferenceStore.contigs_by_chrom`; tests in
`tests/test_contig_names.py`.

## F18 — A Canvas CNV file loads its copy-neutral segments as reference blocks

Severity: medium · Found 2026-10-23 loading the PGP-HMS SV and CNV files for RM282.

Canvas writes each copy-neutral segment as ALT `.` with INFO END, `GT 0/0` and FORMAT `CN 2`, which
is the shape of an Illumina gVCF reference block. 90 of the file's 182 records loaded as `ref`
blocks, and the facet got `absence_means = blocks`, so a consumer reading the blocks would take a
segment's whole span as called reference for small variants, which a copy-number caller never said.
Its `<CN0>`, `<CN2>`… calls also had no extent: `end_pos = pos` on all 92, since END was read only
for `<DEL>`, `<DUP>`, `<INV>` and `<CNV>`.

Resolved 2026-10-23: an `ALT .` record to END that carries a copy number (FORMAT or INFO CN) is a
`calls` row with `alts = []` and `copy_number`, never a block; `<CNn>` alleles end at END. The file
now loads 182 calls, no blocks, every `end_pos` past `pos`, and `absence_means = unknown`. Code:
`loaders/vcf.py` `_is_block_record`, `_SV_LENGTH_TYPES`; test
`test_sv_columns_are_read_as_written`.

## F3 — A headerless file holding two exports loses the second one

Severity: low · Found 2026-09-30, same sweep.

One upload (`rsid_table_headerless`) is a 23andMe-shaped export followed by an AncestryDNA-shaped one
(5 columns). The loader keeps the 643,487 four-column rows and drops the 677,437 five-column rows with
a warning, as FORMATS.md predicted. Not fixed: splitting one file into two facets needs a second format
per file, which the loader does not have.

**Measured and surfaced 2026-10-22.** Both halves are one person: of the 164,378 ids both halves call,
164,275 have the same genotype (99.94%), and 166,047 of 166,528 shared ids sit at the same position.
The first line is a stray 5-field row (an empty field between position and genotype, not F2's place),
then 643,487 four-field rows, then 677,436 five-field rows to the end (chromosomes 1–26, `0` for no
call). The load now says so: a file that ends in at least 1,000 rows of one other width gets a warning
that it looks like a second export appended, not read. Still not loaded, and why each repair is wrong:

- *Read the 5-field rows as a second facet in this loader:* the library cannot load that half even as
  its own file (`ancestrydna` needs its header line), so it needs a new headerless two-allele format id
  and a part option in the facet id, for one file whose halves state no vendor. Calling the second half
  AncestryDNA would be the shape standing in for a vendor statement.
- *Split the file into two inputs and load each:* that is routing around the loader (CLAUDE.md §7).

What would close it: a headerless two-allele-column layout in `formats.rsid_tables`, and `load`
returning one facet per layout block of a file (the shape `Loaded` already has for caller splits).
Code of the warning: `loaders/rsid_table.py` `_appended_export`; test
`test_an_appended_export_in_another_layout_is_named_not_read`.

**Fixed 2026-10-22; F3 resolved.** The loader reads the second export as a facet of its own, in the
same `Loaded` as the first: a run of at least 1,000 rows of one other width that ends the file, after
other rows, whose first five rows are in the two-allele layout (id, chromosome, position, one allele
per column), is read by that layout's rules (`0` is no call, 23-26 are X, Y, XY, MT) and gets
`relation = appended_export`, `part = 1` in its facet id; both facets keep the file's `format_id`.
The repairs listed above were wrong for the reasons given there, but the format id turned out not to
be needed: the registry names the file, and the loader names the layout of its second half (the
alternatives and their costs are in SCHEMA.md Decision 8, "Two exports in one file"). On the upload:
643,487 rows in the first facet (chip `23andme_v5`, inferred as before) and 677,436 in the second
(chip null, 26 contigs, 4,463 whole no-calls); the stray first row is still dropped and counted.
Joined on id through the two facets' `calls`, 164,275 of the 164,378 ids both call have the same
genotype and 166,047 of 166,528 shared ids the same position, the figures measured above by hand.
Code: `loaders/rsid_table.py` `Appended`, `_appended_start`, `_read_appended`; `loading.py`
`_load_rsid_table`; tests `test_an_appended_export_in_another_layout_is_its_own_facet`,
`test_appended_exports_keep_their_disagreements`, `test_a_tail_that_is_short_or_in_no_known_layout_is_not_read`.

## F7 — Two VCF shapes `sniff` accepts but the loader cannot read

Severity: low · Found 2026-09-30, same sweep.

A VCF re-saved by a spreadsheet (tab-padded, quoted `"##` lines) and Veritas SV/CNV VCFs whose every
header line carries an extra `#` (`###contig`, `##CHROM`) are recognised as `vcf` but fail with "no
#CHROM line". The Veritas shape could be read by accepting `##CHROM` as the column line; the spreadsheet
one needs its own format id, as `23andme_csv` has. Not fixed: three files.

**Veritas shape fixed 2026-10-22.** Not every header line has the extra `#`: `##fileformat` and `##ALT`
are written right, `###FORMAT`, `###contig` and `###reference` carry a third, and the column line is
`##CHROM`. Where no `#CHROM` line exists, a meta line that is VCF's column line with one `#` too many
is now the column line, and the file is read through a rewritten copy; the `###` lines are kept as
the meta lines VCF says they are (a key starting with `#`), so no contig and no FORMAT is declared from
them, and GT and END are declared by the loader as for any undeclared key. Both files load (5,614 SV
calls, all `0/1`; 2,632 CNV rows, all `0/1`, whose ALT `CNV` is F15). Code: `loaders/vcf.py`
`read_header`, `HASHED_COLUMN_LINE`; test `test_a_column_line_written_with_two_hashes_is_read`.


**Spreadsheet shape fixed 2026-10-22; F7 resolved.** The re-save is CSV quoting applied to a tab file:
a field holding a comma or a quote is wrapped in quotes with its quotes doubled (115,123 of 115,220
rows, their INFO holding SnpEff's `EFF=…,…` or their sample field an AD of `12,0`; the 97 unquoted
rows are `./.`), meta lines are padded with tabs to ten fields, and lines end in a bare CR, which is
why `read_header` found no `#CHROM` line in 64 MiB. Undoing that is exact, so the rejected option was
the separate format id: the delimiter did not change, as it did for `23andme_csv`, and the quoting
and CRs are a container's mangling, like CRLF, which `sniff` already reads through. The loader looks
at the payload's first 64 KiB (a line starting `"##`, or CRs with no LF); a re-saved file is read
through `csv` with universal newlines, trailing empty fields dropped from `##` lines, and polars-bio
reads the rewritten copy; a warning says so, and that a number whose text the spreadsheet rewrote
(GATK's `6980.90` as `6980.9`) cannot be told. The upload loads all 115,220 rows (a GATK
UnifiedGenotyper exome, `wes`, GRCh37 from its 25 contigs), its genotype tally equal to the raw
text's. Code: `loaders/vcf.py` `is_resaved`, `undo_resave`; `loading.py` `_vcf_text`; test
`test_a_vcf_resaved_by_a_spreadsheet_loads_as_the_text_it_was_saved_from`.

The third failure in the sweep, a corrupt upload (`gzip -t` fails too), keeps its correct error.

**Follow-up 2026-10-22: a field over 128 KiB.** The re-save was undone with Python's `csv`, whose
reader refuses a field over `csv.field_size_limit()` (131,072 bytes by default) with `csv.Error`, and
a VCF field has no such limit (SnpEff's EFF on a busy site, quoted for its commas). The limit is
process-wide, and the reader is consumed lazily while the copy is written, so it cannot be raised for
this reader alone. The re-save is now undone without `csv`: lines by universal newlines, cut at tabs,
a field wrapped in quotes unwrapped with its doubled quotes made single. A VCF field holds no tab and
no line break, so nothing carries across fields or lines. Rejected: raising the process-wide limit,
which widens every other `csv` reader in the caller's process and removes the guard against an
unclosed quote swallowing the rest of a file. Code: `loaders/vcf.py` `undo_resave`, `_unquote`; tests
`test_a_resaved_vcf_with_a_field_past_the_csv_limit_loads` (about 200 KiB of EFF) and
`test_undo_resave_reads_quotes_as_csv_does` (the same output as `csv` on the shapes it reads).

## F4 — Living DNA indel genotypes cannot be split

Severity: medium · Found 2026-09-30, same sweep.

Living DNA writes an indel's two allele sequences back to back (`CCA`, and up to 52 letters in one
file). 83 rows in the one real file keep the text in `source_gt` with `gt_alleles` null. Splitting needs
the site's alleles, from a chip manifest or from REF once RM3 reads the sequence.

**Anatomy, 2026-10-22.** Of the 83 texts (3 to 52 letters), 59 are a sequence written twice
(`GTTGTT`), so homozygous if split in half, and 24 are not. If both alleles share their first base,
as VCF's anchored alleles do, 58 have exactly one split (`CCA` → `C`/`CA`, `GCCTG` → `GCCT`/`G`) and 25
have several (`AAAAGAAAAG` splits four ways). Splitting at the one anchored point was rejected: the
shared anchor is VCF's convention, Living DNA states no encoding beyond "both allele sequences
concatenated", and a quarter of the rows stay ambiguous even under it. The honest path is
`resolve_ref`: read the sequence at `pos` and keep only the splits in which one allele is the
reference over its own span. Even that leaves a repeat ambiguous (`C`/`CA` where the base after is
`A` fits an insertion and a deletion alike), so the split may be known while REF and ALT are not.
Not built: it needs the reference, which the loader does not read.

**Resolved 2026-10-22, the split.** `resolve_ref` splits a `living_dna` facet's run-together texts
(bases only, three or more) where exactly one split fits: both pieces start with the same base and one
equals the reference over its own span from `pos` (rule in SCHEMA.md Decision 2). On the real file,
against GRCh37: 58 of the 83 split (`gt_alleles` set), 24 fit at several splits and 1 at none, and those
25 stay as loaded, counted. The 58 are not quite the 58 with one anchored split above: 57 of those fit
the reference and one does not (a doubled text whose one split matches no reference span), and one of
the 25 with several anchored splits fits at exactly one. REF and ALT follow on 4 of the 58, all
heterozygous insertions (`C/CA` where the base after `C` is not `A`); on the other 54 the split is known
and REF is not: 38 homozygous calls (`CA/CA` cannot say whether the site is `CA>C` or `C>CA`, so
whether the call is reference), and 16 hets whose two alleles both equal the reference over their own
spans (every het deletion is one: its one-base ALT always does). Those 54 and the 25 unsplit rows need
the chip's own site alleles, which is F17. Rejected: reading a homozygous call as hom-ref when its
allele equals the reference, which picks one of two sites the text fits; splitting other formats' long
texts, whose encoding no vendor states; splitting in a `mixed`-strand facet, where a text may be written on the minus strand, and a split
fitting the forward reference would then fit by chance. The 4 rows with REF are not left-aligned by `normalize`
(their `ref_check` is null and the facet is not derived). Code: `reference.py` `_split_concatenated`,
`CONCATENATED_ALLELES`; tests `test_living_dna_indel_genotypes_split_against_the_reference`,
`test_only_living_dna_texts_are_split`, `test_a_mixed_strand_living_dna_file_splits_nothing`.

## F16 — A damaged BGZF input fails with polars-bio's raw error

Severity: low · Found 2026-10-22 on the damaged starling upload ([PGP_HMS.md](PGP_HMS.md)).

One block of the upload does not inflate. `load` failed with a `polars.exceptions.ComputeError`
whose text depends on how the block is broken ("invalid BGZF ISIZE" on the real file; a checksum
or libdeflater message on synthetic ones) and says neither which facet nor where in the file.

Resolved 2026-10-22: when a read fails (polars-bio's `ComputeError`, or zlib's or `gzip`'s from a copy
made in Python) and the input's outer layer is BGZF, `load` walks the file block by block
(`bgzf.find_damage`: header, `BC` size, inflate, CRC32, ISIZE) and raises `DamagedBgzfError` with the
facet's id, the first bad block's compressed offset and what is wrong with it; a failure whose blocks
all read is re-raised unchanged. The error carries no path. Damage that keeps `sniff` from reading the
first blocks names the source's sha256 instead of a facet. The check costs nothing on a file that
reads: it runs only after a failure, and every read of the records happens inside `load` (the stats
pass reads the whole file), so the lazy frames never meet the damage. On the real upload: offset
473,789,118, `inflate`, found by a 4-second walk, as the hand walk in PGP_HMS.md had it. Rejected: a
block check before every BGZF load, a full extra decompression of every input; and matching
polars-bio's message texts, which differ by kind of damage and by release. Not covered: damage under
another layer (a BGZF inside a zip, a plain gzip), which keeps the reader's own error, since only an
outer BGZF has block offsets to name. Code: `bgzf.py` `find_damage`; `loading.py` `DamagedBgzfError`,
`_damage_named`; tests in `tests/test_bgzf.py` (`test_load_names_the_facet_and_the_block_never_the_path`
and the walker's cases).

## F15 — An ALT that is no VCF allele counted as an indel

Severity: low · Found 2026-10-22 loading the Veritas CNV upload F7's fix made readable.

The file writes every ALT as the letters `CNV`, without angle brackets. The class rule took any run of
letters for a sequence allele, so the facet's `variant_classes` said 2,632 indels and the callset
inference counted the rows as small-variant sites. VCF allows only `ACGTN` in a sequence allele.

Resolved 2026-10-22: a sequence ALT is `ACGTN` bases only, in `variant_class_expr` and in the
callset's site count, so `CNV` (or an ambiguity code such as `R`) has no class and is no site; the
loader keeps it verbatim and counts the records in one warning. Code: `just_vcf/callset.py`
`variant_class_expr` and `spread`, `loaders/vcf.py` `scan_stats`; test
`test_an_alt_that_is_no_vcf_allele_has_no_variant_class`.

## F2 — 23andMe rows with a doubled tab are dropped

Severity: low · Found 2026-09-30 by loading the 503 PGP-HMS 23andMe-layout files.

Two 2012–2014 23andMe exports carry rows with an empty field between the rsid and the chromosome
(`rs15842⇥⇥1⇥948921⇥CT`): 14 and 35 rows out of about 960,000. The loader drops them and says so in
`Loaded.warnings`. Not fixed: accepting "exactly one empty field after the id" is a guess about what the
vendor meant, and the rows are few. Revisit if another vendor shows the same shape.

Resolved 2026-10-22. The rows were checked rather than the vendor's intent guessed: in both exports
every such row (49 distinct ids, all rsids of four or five digits) has its one empty field right after
the id, and each of the 49 sits at the same chromosome and position in other PGP-HMS 23andMe exports
(the first 120 checked), with none placed elsewhere. It is not tab-stop alignment: most ids of the
same length are written with one tab. So a row of exactly five fields whose second field is empty is
read without that field, and counted in its own warning line; an empty field anywhere else still drops
the row (the stray first row of the F3 file, empty between position and genotype, stays dropped). Both
exports now load 960,613 rows, the same count. Code: `just_vcf/loaders/rsid_table.py` `read_calls`;
test `test_a_doubled_tab_after_the_id_is_read_past_and_counted`.

## F12 — Some VCF uploads write a deletion as ALT `-`

Severity: low · Found 2026-10-22 running `normalize` against `bcftools norm` on PGP-HMS uploads (RM6).

Five of 25 PGP-HMS VCF uploads of 5–60 MB hold about 8,080 rows each whose ALT is `-`, the dash
notation for a deletion that VCF does not have. The loader keeps it verbatim as an allele, `normalize`
leaves those rows as they are and counts them, and `write_vcf` writes them back out, so the export is a
file htslib refuses ("Non-ACGTN alternate allele"). Nothing in the frame says the allele is not one.

Not fixed here. Reading `-` as a deletion needs the anchor base and its REF span, which the row does
not state, so it is not a spelling to rewrite; the options are for the loader (RM2: flag or refuse
such rows, counted) and the writer (refuse them as unwritable, as it does a row without REF).

Resolved 2026-10-22 (RM6). A survey of the five uploads settled the reading: every such row has one
ALT, `-`, and a REF of one to eight or more bases (about 4,680 of one base, 1,800 of two), and REF
matches the GRCh37 reference at POS on all 40,394 of them, so REF is the deleted span and only the
anchor base is missing. Almost all are `0/0`. Now:

- the loader keeps `-` verbatim and counts the rows in one warning that names them unanchored
  deletions;
- `normalize` reads `-` as the empty allele once `resolve_ref` has matched REF, and the left-align loop
  anchors it on the base before (`CT>-` at 30 becomes `GCT>G` at 29), then aligns it like any
  deletion. At a contig's first base there is no base before; such a row stays as written, counted
  (and the retry loop now stops there instead of widening its window for ever);
- `write_vcf` refuses any ALT that is no VCF allele as `alt_not_an_allele`, counted, or skips it under
  `resolved_only=True`.

Anchoring in `resolve_ref` was set aside, since it fills slots of the same facet and moves nothing;
refusing the rows in the loader was set aside, since their REF matches everywhere it was checked. On
the five uploads: 40,392 rows anchored; the other 2 are rows the loader flagged as one record of a
pre-split site, which `normalize` leaves for `rejoin_split` and `write_vcf` skips as
`alt_not_an_allele`. htslib reads every row of the five exports. Code: `loaders/vcf.py`
`DASH_DELETION`, `normalize.py` `_align`, `writing.py` `_unwritable`; tests:
`tests/test_load_vcf.py` `test_a_dash_alt_is_read_as_the_empty_allele_and_counted`, `tests/test_normalize.py`
`test_a_dash_deletion_is_anchored_on_the_base_before`, `tests/test_write_vcf.py`
`test_an_alt_that_is_no_vcf_allele_is_never_written`.

Changed 2026-10-23 (the lead's decision: the empty allele has one spelling). The loader no longer keeps
`-` verbatim: it reads it as `""`, the spelling CGI's native deletions already had, so `alts`,
`source_alleles` and `gt_alleles` hold `""` and the warning still counts the records. `normalize`
anchors any empty ALT, CGI's included; `mint_vrs_ids` mints a matched one with its anchored spelling's
id (`tests/test_vrs.py` `test_a_deletion_without_its_anchor_base_mints_the_anchored_id`).

## F11 — PL to `gl` lands one float step off the nearest value

Severity: low · Found 2026-10-22 by comparing `normalize` with `bcftools norm` (RM6).

The loader turns PL into `gl` as `PL / -10` in Float32, and Polars computes that division so that PL 21
gives -2.1000001, while the nearest Float32 to -2.1 is -2.0999999. A file that writes the same
likelihood as GL (`-2.1`) therefore loads one step away from one that writes PL 21. Deterministic, so no
digest or test of ours moves, and the size is far below any PHRED rounding; but two facets of one
genome, one with PL and one with GL, do not compare equal on `gl`, and a test comparing against an
outside tool needs a tolerance (`test_agrees_with_bcftools_norm`). Not fixed here: it is the loader's
PL rule (RM2). The repair is to divide in Float64 and cast the result to Float32 once.

Resolved 2026-10-22 as proposed: PL is divided by -10 in Float64 and rounded to Float32 once, and so is
Complete Genomics' phred-unit GL (divided by 10). PL `p` now loads as the Float32 nearest `-p/10`, the
value a GL file's `-p/10` is parsed to; checked for every PL from 0 to 299. Documenting a tolerance
instead was set aside: the repair costs one cast. This changes `gl` values on files that already load,
by one Float32 step on some entries. The comparison with bcftools now holds exactly. Code:
`loaders/vcf.py` `to_rows` (`gl`); tests: `tests/test_load_vcf.py`
`test_pl_and_the_same_likelihood_as_gl_load_equal`, `tests/test_normalize.py`
`test_agrees_with_bcftools_norm`.


## F14 — Ensembl's GRCh38_to_GRCh37 chain lists every block twice

Severity: high · Found 2026-10-22 lifting anton's GRCh38 callset with the chain CrossMap had used (RM6).

The Ensembl chain (`GRCh38_to_GRCh37.chain.gz`, from ftp.ensembl.org `assembly_mapping`) holds 22,364
chains, and every block of the main ones appears twice, under two chain ids, with the same source span,
target and strand. The first `liftover` counted blocks per base, so 5,717,715 of 6,079,744 rows read as
lying in more than one block (`multiple_targets`) and none was lifted. UCSC's hg38ToHg19 chain has no
overlapping block at all, so the synthetic tests, built like it, could not show it.

Resolved the same day in `liftover._segments`: the source axis is cut where the set of *targets*
changes, a target being a contig, strand and offset, so two blocks sending every shared base to the same
place count once. A position under blocks with different targets is still `multiple_targets`. With the
fix the same chain lifts 5,578,476 rows, and 108,776 stay ambiguous. Test:
`test_a_chain_listing_its_blocks_twice_lifts_as_once` (fails on the first version, checked).

## F13 — An FTDNA export reports some probes on the minus strand

Severity: medium · Found 2026-10-22 by F9's row inference on the 103 PGP-HMS FTDNA exports (counts
only).

One autosomal FTDNA export (about 709,000 rows) fits GRCh36 at 68.2% of sampled SNV rows, against
80–83% for the other 59 autosomal FTDNA exports. Of its heterozygous rows, only 81% carry the reference
base (99.96–100% in every other export with autosomes), and 13.0% of all sampled rows carry the complement of
the reference base (0.0–3.0% elsewhere). So a share of its probes, about one
in eight, is written on the minus strand. `resolve_ref` reads REF for 99.8% of its rows, and on the
flipped rows derives `alts` and `gt` against the wrong strand: a flipped hom-ref row reads as hom-alt.
The existing warning counts them among the rows that "carry no REF allele", which does not single
them out.

Not fixed. Candidate repairs, each wrong for now: complementing a row whose alleles miss REF (an A/T
or C/G probe is its own complement, and a true hom-alt call would be flipped into a hom-ref one);
scoring both strands (a complemented het covers all four bases, so a wrong build would fit as well as
the right one); a per-chip strand table (the right fix, but FTDNA ships none in the export, and the
chip is not identified yet: RM4's chip generation detection).

Resolved 2026-10-22 by withholding, not by complementing. `resolve_ref` now checks the strand of every
array facet it reads, on F9's row sample against the collection it reads, through the heterozygous
rows: on a two-allele site a het carrying the reference base is safe on either strand, and a het
carrying the base's complement and not the base is impossible on the forward strand. More than 2% of
at least 200 such hets makes the facet `strand_check = mixed`; then every row that is not strand-safe
(homozygous calls, half-calls, hets without the base) keeps `ref` and gets `alts` and `gt` null, each
reason counted in the warnings, and the safe hets keep their derived alleles. Rule: docs/SCHEMA.md
Decision 2 ("The strand of an array's genotypes"). Code: `strand.py` `check_strand`,
`withheld_reason`; `reference.py` `_check_strand`, `_derive_array_alleles`. Tests:
`tests/test_reference.py` from `test_a_mixed_strand_array_keeps_ref_and_withholds_the_alleles_it_cannot_read`
(each fails on the code before this change).

On F13's export: `mixed`, 2,351 of about 12,600 sampled hets (18.65%) carry the complement and not the base,
none carry neither. On the other FTDNA exports of 600,000 rows or more that ran before the scan was
stopped (10 of them, on OmniExpress like this one): 0.00–0.02%, all `forward`. On 32 more, two of
each format and detected chip (23andMe v1/v2, v3, v4, v5 and unknown, its CSV re-save, AncestryDNA
V1.0 and V2.0, deCODEme, FTDNA space-separated, Genes for Good measured and imputed, the GSA Final
Report, Living DNA, SelfDecode, the rsid CSV and the rsid tables), 1,171–13,794 hets sampled each:
every one `forward`, 0 or 1 hets carrying the complement and none carrying neither (at most 0.01%).
The 2% threshold sits between 0.02% and 18.65%.

Why each other repair is wrong:

- **A per-chip strand table.** The chip of this export is not detected (FTDNA states none and RM4 has
  no rule for it), but its probe set is OmniExpress's: 98–99% of its sampled rsids are in Illumina's
  Infinium OmniExpress-24 v1.2 manifest, as for every FTDNA export of that size. That manifest's three
  strand columns do not separate its minus-strand probes. Over all hets of the export (aggregates,
  matched by rsid): proven-forward or strand-invariant hets split RefStrand `+` 83,306 / `-` 83,552,
  proven-minus hets 19,337 / 19,467; IlmnStrand TOP 85,247 / BOT 81,611 against 19,821 / 18,983;
  SourceStrand TOP 83,213 / BOT 83,645 against 19,205 / 19,599. The flipped probes are on every
  autosome (16–24% of hets each) and not in runs (mean run 1.54 hets in file order). So no public
  table says which probes FTDNA wrote flipped in this export, and complementing by any of these columns
  would flip about half of the right rows and half of the wrong ones.
- **Complementing the hets that prove a minus probe.** It is a deduction for those rows, but the
  homozygous rows of the same probes cannot be found, so the facet would mix corrected and uncorrected
  rows with nothing in the row to tell them apart. Withholding is uniform: every non-safe row is null.
- **Complementing every row whose alleles miss REF, or scoring both strands**: as above (an A/T or
  C/G probe is its own complement; a complemented het fits any build).
- **Withholding nothing and warning only**: a flipped hom-ref row would still read as a confident
  hom-alt call downstream, which is the error F13 is about.

Addendum 2026-10-23 (RM281): the withheld rows of a mixed facet are now read by their probe's two
alleles on the plus strand, where a configured manifest names the rsid, the pair holds the reference
base and is not `A/T` or `C/G` (docs/SCHEMA.md Decision 2). That is not the per-chip strand table
rejected above: the columns still do not say which probes FTDNA flipped, but `TT` at a probe whose plus
pair is `A/G` is a minus-strand `AA` whichever they are, so homozygous rows of such probes are read
too. Not yet measured on F13's export (build 36, so matched by rsid; the pair is checked against the
build-36 reference base): the search through the PGP-HMS FTDNA exports had not reached it when RM281
was handed over.


## F10 — A GP that does not sum to 1 is lost when just-vcf reloads its own export

Severity: low · Found 2026-10-22 by round-tripping a PGP-HMS VCF upload through `write_vcf` and `load`
(RM6).

The upload writes GP phred-scaled (values above 1), so the loader converts it and every row gets `gp`.
Phred values rounded to integers do not convert to posteriors that sum to 1: 57,259 of its 956,057 GP
rows are more than 0.02 short. `write_vcf` writes `gp` as probabilities under a VCF 4.3 header, where
the spec defines GP as a probability. Reloading that output, the loader reads the scale from the values
alone: all within [0, 1], and not all summing to 1, so it calls GP ambiguous and drops it from every
row. Everything else in the file came back byte for byte.

Not the writer's to fix: normalizing the posteriors would change the evidence, and writing them back as
phred under a 4.2 header would not read back bit for bit.

Resolved 2026-10-22 in the loader's GP rule: values the numbers cannot settle (all within [0, 1], not
all summing to 1) are read as probabilities under a VCF 4.3+ header, which defines GP so, and the rows
that do not sum to 1 are counted in a warning. A phred GP whose every value is at most 1 would mean
genotypes all but equally likely, which is not what any caller writes. Under 4.2 or earlier, or a header
with no version, they stay withheld as before; values above 1 are phred and values summing to 1 are
probabilities, whatever the header. The FORMAT line's Description was considered as a second header
signal and set aside: it is free text, not a declaration. On the upload that found it: the reload keeps
`gp` on all 956,057 rows, equal to the first load's (largest difference 0.0), with one warning counting
the 57,259 rows that do not sum to 1. Code: `loaders/vcf.py` `_gp`; tests:
`tests/test_load_vcf.py` `test_gp_from_phred_survives_its_own_export` and the 4.3, 4.5 and no-version
cases of `test_gp_scale_is_read_from_decisive_values`.
## F9 — An array export that states no build cannot get its REF

Severity: high · Found 2026-10-22 by running `resolve_ref` (RM3) on the table formats RM2 just made
loadable, on public PGP-HMS uploads (counts only).

`resolve_ref` worked where the facet carries an assembly: an AncestryDNA export (stated GRCh37) had REF
read for 681,662 of 701,478 rows (97.2%), a Color file for all of its rows. Two FTDNA exports were
refused with `ReferenceMismatchError` ("the facet states neither contig lengths nor an assembly"): an
FTDNA file names no build, and an array has no `##contig` lines, so RM4's length inference has nothing
to read. A deCODEme export (GRCh36, stated in its zip) was refused because the store holds only GRCh37
and GRCh38. The refusals are correct (never assume a build), so the gap is evidence, not the check.

Not fixed. Candidate repairs: (a) infer an array's build from its own rows against each collection in
the store, as the share of SNV rows whose called alleles include the reference base at that position
(the right build scores near 100%, a wrong one near chance), recorded as `inferred` with the scores; (b)
add NCBI36 to the reference store for the GRCh36 exports. A caller-asserted `assembly=` argument was
considered and set aside: it lets a label stand in for evidence.


Resolved 2026-10-22, both repairs. (a) `resolve_ref` infers the build of a facet that states neither
a label nor contig lengths from its own rows, against every collection it considers, and records it as
`assembly_label` with basis `inferred` and the scores in one warning line; rows that do not settle it
refuse the facet with the scores (rule: docs/SCHEMA.md Decision 2; choices and rejected alternatives:
ROADMAP RM3). (b) A collection whose contig 1 is 247,249,719 bases is labelled `GRCh36`; UCSC hg18 is
the FASTA to add (where and why not Ensembl 54: RM3). Code: `assembly.py` `sample_snv_rows`,
`infer_assembly_from_rows`, `MIN_SHARE`/`MIN_MARGIN`; `reference.py` `_infer_from_rows`, `_bases_at`,
`CHR1_LENGTHS`. Tests: `tests/test_reference.py` `test_an_array_without_a_build_gets_it_from_its_rows`
and the refusals after it, `test_a_grch36_reference_is_named_by_its_contig_1_and_cannot_vouch_for_mt`;
`tests/test_assembly.py` from `test_the_sample_keeps_scorable_array_snvs_only`.

Dogfooded on 153 loadable PGP-HMS array exports from the Krypton share, against a store holding Ensembl
GRCh37 and GRCh38 and UCSC hg18 (aggregates only). Share of sampled SNV rows carrying each collection's
reference base, the right build against the wrong ones, and the share of rows whose REF `resolve_ref`
then read:

| Exports | Files | Inferred | Right build | Wrong builds | REF read |
|---|---|---|---|---|---|
| FTDNA, autosomal + X | 60 | 36 GRCh36, 24 GRCh37 | 0.682–0.825 | 0.311–0.327 | 95.4–99.9% |
| FTDNA, X only (about 2,000 rows) | 43 | 30 GRCh36, 13 GRCh37 | 0.665–0.843 | 0.230–0.334 | 94.4–99.9% |
| FTDNA, space-separated | 2 | GRCh36 | 0.825–0.827 | 0.313–0.318 | 99.4–99.6% |
| deCODEme (build stated in the zip) | 4 | GRCh36, as stated | 0.811–0.814 | 0.307–0.316 | 96.0–99.9% |
| 23andMe, stated build 36 | 13 | GRCh36, as stated | 0.793–0.810 | 0.310–0.329 | 98.3–99.4% |
| 23andMe, stated build 37 | 13 | GRCh37, as stated | 0.801–0.906 | 0.283–0.324 | 95.5–99.1% |
| AncestryDNA, Genes for Good, Living DNA, 23andMe CSV (stated) | 10 | as stated | 0.806–0.892 | 0.283–0.327 | 96.2–99.99% |
| rsid tables, PLINK 23, rsid CSV, Final Report (no build) | 8 | 2 GRCh36, 6 GRCh37 | 0.796–0.913 | 0.263–0.326 | 96.6–99.4% |

For the stated-build exports the inference ran as if no build were stated: it named the stated build
in all 40. hu832966's two 23andMe exports infer GRCh36 and GRCh37, each its stated build; REF is read
for 99.2% and 97.9% of their rows. Before this, every FTDNA export, every unlabelled table and every
build-36 export was refused. One Final Report did not load (`TableLayoutError`, not this finding's).

Strand: forward only. The share of rows carrying the complement of the right build's base was
0.0001–0.030 in 152 files, so complementing would gain nothing there and would let a wrong build fit
(a complemented het covers all four bases). The 153rd, an FTDNA export at 0.13, has part of its probes
on the minus strand: F13. The floor was first 0.70; it refused 27 FTDNA X-only exports and F13's file,
whose right build stood 0.36–0.46 above every other, so it is 0.60 and the 20-point margin carries the
separation. Cost: with the store warm, sampling and reading the bases of three collections takes about
0.3 s per facet.


## F8 — A sites-only VCF does not load

Severity: medium · Found 2026-10-19 while testing `rejoin_split` (RM270).

A VCF with no FORMAT column and no sample (a sites-only file, as annotation sources and some
callers' site lists are) fails in `load` with polars' `InvalidOperationError: to_arrow() conversion
failed: cannot convert (Null) to (... LargeUtf8 ...)`, raised by the parquet spill (`_load_vcf`'s
`sink_parquet`). It failed the same way before RM270 (checked at c746325), so it is older than the
split flags. polars-bio itself reads the file: `pb.scan_vcf` returns the two records with their eight
site columns. So the fault is on our side, most likely a column the loader adds as an untyped null
for the absent sample. Repro: a two-record VCF 4.2 with one `##contig` line and an eight-column
header. A sites-only file has no genotype, so it is a `calls` frame of nulls at best.

Resolved 2026-10-22: `load` refuses it with `NoSampleError`, before anything is read or spilled. just-vcf
holds a person's calls (the charter's first line), and a file with no sample column has none to hold;
a `calls` frame of null genotypes would look like a person with nothing called. `vcf_samples` still
reads the file and returns `[]`, and `sample=` on it is refused as a missing sample, as before. Code:
`loading.py` `_load_vcf`; test: `tests/test_load_vcf.py` `test_sites_only_vcf_is_refused`.

## F6 — A bz2 or zip VCF was copied out as plain text first

Found 2026-09-30 by loading the 124 PGP-HMS VCFs. polars-bio reads plain, gzip and BGZF VCF directly
(79 of the 124); a VCF inside bz2, xz, zip or tar is decompressed into the loader cache first. That copy
was plain text: a 1.27 GB bz2 Veritas gVCF became a 14 GB file.

Resolved 2026-09-30: the copy is written as BGZF by `just_vcf.bgzf.BgzfWriter` (stdlib `zlib` on a thread
pool, no new dependency). On a 2 GB slice of a real VCF: BGZF at level 1 is 0.50 GB and takes 17 s to
write; polars-bio counts its 18.5 million records in 9.5 s, against 9.1 s for plain text and 24.8 s for a
plain gzip copy. Code: `loaders/vcf.py` `readable_path`.

## F5 — The canonical sort cost most of a VCF load's memory

Found 2026-09-30 loading a DeepVariant VCF of 6.4 million calls: writing `calls` sorted peaked at
11.2 GB, against 3.5 GB unsorted. The user picked "keep file order within a contig, sort only contig
blocks and ties". Measuring it showed that reordering contig blocks is not cheap either: polars-bio does
not use the `.tbi` index, so a per-contig filter is a full scan (4.3 s each, against 6.2 s for the file).

Resolved 2026-09-30 without changing SCHEMA.md Decision 11: `load` spills the VCF frames unsorted to
parquet in the cache, in one streaming pass, and returns scans that read them back one contig at a time
in canonical order (`schema.canonical_order_per_contig`), using parquet row-group statistics to skip to
each contig. Same file: 77 s and 3.6 GB peak, row order identical to the global sort. The spill was
0.14 GB. Code: `loading.py` `_load_vcf`; test: `tests/test_schema.py`.
