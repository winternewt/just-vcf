# Previous issues — resolved here

Findings from [dogfooding.md](dogfooding.md) resolved in this repo, each with its resolution and a
code pointer. Check here before re-investigating something that looks fixed.

None yet.

## F6 — A bz2 or zip VCF was copied out as plain text first

Found 2026-09-30 by loading the 124 PGP-HMS VCFs. polars-bio reads plain, gzip and BGZF VCF directly
(79 of the 124); a VCF inside bz2, xz, zip or tar is decompressed into the loader cache first. That copy
was plain text: a 1.27 GB bz2 Veritas gVCF became a 14 GB file.

Resolved 2026-09-30: the copy is written as BGZF by `just_vcf.bgzf.BgzfWriter` (stdlib `zlib` on a thread
pool, no new dependency). On a 2 GB slice of a real VCF: BGZF at level 1 is 0.50 GB and takes 17 s to
write; polars-bio counts its 18.5 million records in 9.5 s, against 9.1 s for plain text and 24.8 s for a
plain gzip copy. Code: `loaders/vcf.py` `readable_path`.
