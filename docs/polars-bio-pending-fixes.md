# polars-bio: pending fixes

Behaviour of polars-bio that needs an upstream change first. Each entry notes the mitigation already in
place here. polars-bio has no local clone or inbox; filing goes to its GitHub issue tracker, and the
issue reference is recorded here once filed.

## P1 — INFO and FORMAT keys missing from the header are dropped silently

Seen with polars-bio 0.36, while declaring the VCF formats (2026-09-29). A key that a row uses but the
header does not define is dropped without a warning, GT included. A header repaired or synthesised by
just-vcf (the `vcf_no_meta` and `freebayes_varscan2_keyvalue` formats) must define every key the rows
use, or the genotypes vanish.

Mitigation: none yet; the loader step must build complete headers. Not filed.

## P2 — Headers that do not start with `##fileformat` are refused

Seen with polars-bio 0.36. A VCF whose meta lines are sorted (`##fileformat` on line 199), or that has
only a `#CHROM` line, is refused. The sniffer recognises both (`vcf`, `vcf_no_meta`); the loader will
need to reorder or synthesise the meta block before handing the file over.

Mitigation: none yet. Not filed; the VCF spec allows `##fileformat` to come later, but it may not be
upstream's bug to fix.

## P3 — `end` for a symbolic SV is one short

Seen with polars-bio 0.36 (2026-09-30). `<DEL>` with `SVLEN=30` at POS 200 gives `end=229`. VCF 4.5
computes END for a symbolic SV as `POS + SVLEN` (230), because POS is the base before the event.

Mitigation: the loader computes `end` itself (docs/SCHEMA.md, Decision 6). Not filed.

## P4 — The leading phase indicator of a GT is dropped

Seen with polars-bio 0.36 (2026-09-30). `|0/1/2` comes back as `0/1/2`. Since VCF 4.4 the first phase
indicator is optional and meaningful: `|0/1/2` says the first allele is phased. Dropping it loses that.

Mitigation: none yet. The loader has to read GT text where a leading indicator can occur. Not filed.

## P5 — Whether a column is INFO or FORMAT depends on what else the header declares

Seen with polars-bio 0.36 (2026-09-30), with one sample. When INFO and FORMAT both declare `DP`, the
columns are `DP` (INFO) and `fmt_DP` (FORMAT). When only FORMAT declares it, FORMAT/DP is `DP`. The
same column name means different quantities in different files.

Mitigation: the loader maps keys by the header's namespace (docs/SCHEMA.md, Decision 3). Not filed;
this may be intended upstream, but it is not documented.
