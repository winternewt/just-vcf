# Roadmap history

Shipped roadmap items, moved here from [ROADMAP.md](ROADMAP.md) with their rationale. Nothing is
deleted; it is relocated.

## RM7 — PGP-HMS sample leecher
Severity: medium · Status: shipped 2026-09-28 · Owner: Claude

A script in `scripts/` that downloads consented-public 23andMe exports from
<https://my.pgp-hms.org/public_genetic_data?data_type=23andMe> into `data/input/`, so tests and
dogfooding have real array files that are allowed to be used. httpx is the dev dependency for it.

Shipped as `scripts/pgp_leech.py`.
Download links answer two ways: older uploads are served directly, newer ones redirect to an Arvados
collection page. The listing's participant cell can carry the original enrolment number
(`huF06AD0, PGP113`), which the live test caught. Rows over 200 MB are skipped by default: the
23andMe listing totals about 101 GB, and 19 of its files are multi-GB uploads.
