# Lead seat — notes

State as of 2026-10-18; no session has held the seat yet. Rewritten to the present before each reset.

## Open

1. **The data cache holds whole private genomes** (RM268): eviction is owed.
2. **The leak gate** (the person's call): the kit's pre-commit gate and `bin/leak-check` read a
   git-ignored `.leak_words`; here there is none, so a commit or a kit report is never checked for a
   named person. A list of the names under `/data/genomes/` would gate both; the gate then refuses
   every commit while the list is missing.

## In flight

Read from the last sittings: RM2 (loaders) and RM3 (REF design) open; the schema at 0.0.3; RM269 (a
person's karyotype finding, from S1) filed by triage for the lead's schema call. Check the roadmap
against what exists before taking these as current (the role, step 1).
