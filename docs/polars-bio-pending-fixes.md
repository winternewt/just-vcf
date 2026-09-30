# polars-bio: pending fixes

Behaviour of polars-bio that needs an upstream change first. Each entry notes the mitigation already in
place here. polars-bio has no local clone or inbox; filing goes to its GitHub issue tracker, and the
issue reference is recorded here once filed.

## P1 — INFO and FORMAT keys missing from the header are dropped silently

Seen with polars-bio 0.36, while declaring the VCF formats (2026-09-29). A key that a row uses but the
header does not define is dropped without a warning, GT included. A header repaired or synthesised by
just-vcf (the `vcf_no_meta` and `freebayes_varscan2_keyvalue` formats) must define every key the rows
use, or the genotypes vanish.

Mitigation: the VCF loader declares every carried FORMAT key the rows use without a declaration, in a
rewritten copy of the file, and says so in `Loaded.warnings` (`loaders/vcf.py`, `readable_path`). Not filed.

## P2 — Headers that do not start with `##fileformat` are refused

Seen with polars-bio 0.36. A VCF whose meta lines are sorted (`##fileformat` on line 199), or that has
only a `#CHROM` line, is refused. The sniffer recognises both (`vcf`, `vcf_no_meta`); the loader will
need to reorder or synthesise the meta block before handing the file over.

Mitigation: the VCF loader rewrites the header with `##fileformat` first and warns. Not filed; the VCF
spec allows `##fileformat` to come later, but it may not be
upstream's bug to fix.

## P3 — `end` for a symbolic SV is one short

Seen with polars-bio 0.36 (2026-09-30). `<DEL>` with `SVLEN=30` at POS 200 gives `end=229`. VCF 4.5
computes END for a symbolic SV as `POS + SVLEN` (230), because POS is the base before the event.

Mitigation: the loader computes `end` itself (docs/SCHEMA.md, Decision 6). Not filed.

## P4 — The leading phase indicator of a GT is dropped

Seen with polars-bio 0.36 (2026-09-30). `|0/1/2` comes back as `0/1/2`. Since VCF 4.4 the first phase
indicator is optional and meaningful: `|0/1/2` says the first allele is phased. Dropping it loses that.

Mitigation: for VCF 4.4+ files the loader withholds `phased[0]` and `source_gt`, with a warning; older
files cannot carry a leading indicator. Reading raw GT text would restore both. Not filed.

## P5 — Whether a column is INFO or FORMAT depends on what else the header declares

Seen with polars-bio 0.36 (2026-09-30), with one sample. When INFO and FORMAT both declare `DP`, the
columns are `DP` (INFO) and `fmt_DP` (FORMAT). When only FORMAT declares it, FORMAT/DP is `DP`. The
same column name means different quantities in different files.

Mitigation: the loader maps keys by the header's namespace (docs/SCHEMA.md, Decision 3). Not filed;
this may be intended upstream, but it is not documented.

## P6 — `scan_vcf` registers a stem-named table on every run

Seen with polars-bio 0.36 (2026-09-30). Two scans of one file in one plan, or `pl.collect_all` over two
frames from it, panic with `The table <stem> already exists`. A named `register_vcf` avoids it, but
returns 0-based half-open `start` where `scan_vcf(use_zero_based=False)` returns 1-based.

Mitigation: `loaders/vcf.py` `scan` registers once per source under a unique name (ported from
just-prs) and converts `start`, checking the first record against the file. Not filed.

## P7 — `register_vcf` cannot restrict FORMAT fields or samples

Seen with polars-bio 0.36 (2026-09-30). `register_vcf` passes `format_fields=None` and `samples=None` to
its internal `SQL._register_variant`, so every FORMAT field of every sample is parsed. One malformed
field the loader does not even carry then fails the load: a Genos exome VCF's `CGIANN_1000GAF` holds a
value that is not a float.

Mitigation: `loaders/vcf.py` `scan` calls `polars_bio.sql.SQL._register_variant` directly with the
carried FORMAT keys and the chosen sample. That is a private API; a polars-bio upgrade can break it, and
the tests will show it. Filing upstream: ask for `format_fields` and `samples` on `register_vcf`. Not filed.

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

## P9 — A record at POS 0 fails the scan

Seen with polars-bio 0.36 (2026-09-30) on Illumina `GenotypingToVCF` exports, which put unmapped probes
at `CHROM NA`, `POS 0`, `REF N`: collecting fails with `Execution error: Missing variant start`. VCF 4.5
allows POS 0 (a telomere).

Mitigation: on that error the loader reads a copy without the POS 0 records and says how many were left
out. Not filed.
