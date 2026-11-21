# gtars: pending fixes

Behaviour of gtars (databio/gtars, a core dependency: VRS allele ids, allele normalization, refget)
that needs an upstream change first. Each entry notes the mitigation already in place here. Filing
goes to gtars' GitHub issue tracker and needs the person's go per batch; the issue reference is
recorded here once filed. The upstream-liaison seat keeps this ledger, beside
`docs/POLARS_BIO_PENDING_FIXES.md`.

The mirror: `data/interim/gtars/`, a read-only clone of databio/gtars at the tag `uv.lock` pins
(`gtars-python-v0.10.0`, `90141b6`), for locating code paths; never edited, never filed into. Refresh
with `git -C data/interim/gtars fetch --tags` and a checkout of the new `gtars-python-v<version>` tag.
The seat's `upstream-index` watcher says a final gtars release on PyPI the day it is out (added
2026-11-16); the seat's start check still compares `uv.lock`'s gtars version against the pin above.
Read mark: gtars-python 0.12.0 (read 2026-11-20: G1 still true there, G2's one-base case reproduces
there); the lock is on 0.10.0. Bump note for 0.12: `tests/test_reference.py`, `test_reference_support.py`,
`test_unstated_ref.py` and `test_vrs.py` pass on it (105 passed, 5 skipped), but its collection aliases are
one-to-one and refuse a re-point (`OSError: alias conflict`), so `ReferenceStore.add_fasta(…,
assembly_label=X)` for a second collection of X raises there; on the lock it silently moves the label off
the first collection (synthetic two-FASTA probe, 2026-11-20). The lead closed it in `ed7edfa`: `add_fasta`
refuses a label another collection holds with a `ReferenceMismatchError`, on 0.10 and 0.12 alike, so the
bump has nothing left from this.

## Verdicts at a glance

| Item | Verdict | File at |
|---|---|---|
| G1 VRS ids differ from vrs-python on deletions, repeat insertions and reference alleles | normalization never returns a `ReferenceLengthExpression`; the Python API digests literal states only | not filed: desk #77 closed on its default (2026-11-21); databio/gtars if that changes |
| G2 `add_sequence` into an Encoded store reads back wrong bases (2026-10-26) | a record's ASCII bytes are kept as if already bit-packed; fixed for two bases and more in 0.11.0 (`03656182`); a one-base sequence still reads back `C` on 0.12.0 through Python `add_sequence`, whose Rust-only escape `ingest_sequence` Python lacks | not filed: rode with G1 on desk #77, closed on its default (2026-11-21) |

## G1 — VRS ids differ from vrs-python on deletions, repeat insertions and reference alleles

Found by the lead seat building RM271 (merged as `b41bda5`), on 19,000 real PGP-HMS alleles. VRS 2.x
normalization ends a deletion, a reference-derived (repeat) insertion and a reference allele in a
`ReferenceLengthExpression` state, not a literal sequence: vrs-python 2.3.3 does so at spec steps 2.a,
5.b and 5.c (`ga4gh/vrs/normalize.py:161-167`, `:199-203`, `:206-223`). gtars 0.10 always ends in a
literal sequence, so its ids for those alleles are not the VRS ids every other implementation mints.
SNVs and substitutions agree.

Mitigation: just-vcf writes the VRS state and the canonical JSON itself (`src/just_vcf/vrs.py`) and uses
gtars only for rolling and `sha512t24u`; its ids equal vrs-python's on the lead's 19,000 alleles.
RM271 in `docs/ROADMAP.md` records it, and `docs/LANDSCAPE.md`'s gtars row is corrected. Once gtars
has the state, `_justify` and the JSON writer shrink to calls. Not filed: desk #77 closed unanswered on
2026-11-21, so its default, not filed, is the decision.

*Upstream (2026-10-22).* Not reported (gtars' tracker searched live for "ReferenceLengthExpression" and
"vrs-python": nothing). Reproduced with a synthetic 32-base sequence and its own `SQ.` digest
(no real data): gtars 0.10.0 `normalize_allele` plus `vrs_id` against vrs-python 2.3.3 `normalize` plus
`ga4gh_identify` agree on an SNV and differ on a deletion in a T run, an insertion in that run, a
deletion with no repeat and a reference allele; vrs-python's state is a `ReferenceLengthExpression` in
all four. Where it is in gtars: the Rust model and digest know the state
(`gtars-vrs/src/models.rs:27-31`, `gtars-vrs/src/digest.rs:125-135`, checked by
`gtars-vrs/tests/vrs_compliance.rs:161-175`), but `normalize_ref` returns a start, an end and literal
bytes (`NormalizedAllele`, `gtars-vrs/src/normalize.rs:170-174`), the Python `vrs_id` and `vrs_digest` build only a literal
allele (`gtars-python/src/vrs/funcs.rs:10-41`), and the VCF batch path digests through the literal-only
fast writer (`gtars-vrs/src/vcf.rs:149-154`; `digest.rs:37-41` says it does not support the other
state). gtars main after 0.10.0 changes none of these files. **File at databio/gtars:** a VRS 2.x
normalization result (the RLE state where the spec ends in one), in `normalize_ref`, the Python API and
the VCF batch path; cite vrs-python's steps and the repro.

*Rechecked (2026-11-14).* Still true at gtars-python v0.11.2, the newest tag (2026-11-03): `NormalizedAllele` (`gtars-vrs/src/normalize.rs:170`) still carries literal bytes, the Python `vrs_id` and `vrs_digest` build only a literal allele (`gtars-python/src/vrs/funcs.rs:11`), and the VCF path digests through `allele_identifier_literal` (`gtars-vrs/src/vcf.rs:154`). Nothing on gtars' tracker. Desk #77 updated: still waiting, nothing filed by default. Again at gtars-python v0.12.0 (2026-11-20): the release touches only `ref_view_for` in `gtars-vrs/src/normalize.rs` (how a stored sequence is viewed); `NormalizedAllele` (`:170`), the literal-only Python `vrs_id` and `vrs_digest` (`gtars-python/src/vrs/funcs.rs:11`) and `allele_identifier_literal` in the VCF path (`gtars-vrs/src/vcf.rs:154`) are as they were.

## G2 — `add_sequence` into an Encoded store reads back wrong bases

Found by a just-vcf builder (the lead's relay, commit `00f442f` on its branch) as lowercase: a
`SequenceRecord` whose `sequence` was set by hand to `acgtacga` reads back as `CATCCATG` from a store in
`StorageMode.Encoded`, the mode `RefgetStore.in_memory()` starts in. Case is not the cause. Reproduced on
gtars 0.10.0 with synthetic sequences only: `ACGTACGA` reads back `CTTCCTTG` too, and so does the record
`digest_sequence(b"ACGTACGA")` makes, untouched, which is the usage `add_sequence`'s own docstring shows.
In `StorageMode.Raw` every case reads back as given (lowercase stays lowercase). The DNA lookup tables
map lowercase to the right codes (`gtars-refget/src/digest/alphabet.rs:129-141`, `:157-175`); the defect
is that 0.10.0's `add_sequence_record` stores a full record's bytes as they are
(`gtars-refget/src/store/readonly.rs:835-880`), while an Encoded store decodes what it holds as packed
bits. Its debug-only assert allows ASCII bytes in either mode, so nothing refuses them in a release
build. A store built from a FASTA is not affected: that import packs the bases itself.

Mitigation: none needed. just-vcf adds sequences only through `add_fasta` (the FASTA import); `src/`
never calls `add_sequence`, and the one test that does (`tests/test_reference.py`, from `00f442f`) sets
the store to Raw first. `ReferenceStore.fetch` still uppercases a sequence the store holds no FASTA
record for.

*Upstream (2026-10-26).* Not reported on gtars' tracker (searched live: "add_sequence", "encoded mode",
"SequenceRecord encoded"), and **already fixed on trunk**: `03656182` (2026-09-16, "Fix encoded sequence
ingestion and add byte reads") packs a record whose bytes look like ASCII before storing it in an
Encoded store (`pack_raw_ascii_for_encoded_mode`), and is in the `gtars-python-v0.11.0` tag; PyPI's
newest is still 0.10.0. Nothing to file. Read, not run: the new check takes bytes as ASCII only when
their length differs from the packed size, which for the 2- and 3-bit DNA alphabets is false at length 1,
so a one-base hand-added sequence would still be stored unpacked. Too small to file on a reading; check it
when `uv.lock` moves to 0.11 and file then if it reproduces.

*Rechecked (2026-11-20), gtars 0.12.0.* The one-base case reproduces, run in an overlay
(`uv run --with gtars==0.12.0`) on synthetic bases only: `digest_sequence(b"A")` (and `C`, `G`, `T`) added
with `add_sequence` to `RefgetStore.in_memory()` reads back `C` from `get_substring`; `AC` and `ACGTACGA`
read back as given. Upstream knows the ambiguity: 0.12.0's `pack_raw_ascii_for_encoded_mode`
(`gtars-refget/src/store/readonly.rs:1267-1272`) says a single-base sequence is ambiguous and points to
`ingest_sequence` (`:1308`, added with `03656182`), which takes raw bytes and packs them always. That
method is Rust only; the Python `add_sequence` (`gtars-python/src/refget/mod.rs:1766`), whose own
docstring shows `digest_sequence` then `add_sequence`, still goes through the ambiguous check. Not on
gtars' tracker (searched live: "single base", "add_sequence", "ingest_sequence", "one base"). Still no
mitigation needed here: `src/` never calls `add_sequence`. **File at databio/gtars:** a one-base record
from `digest_sequence` added with Python `add_sequence` to an Encoded store reads back wrong; expose
`ingest_sequence` to Python, or have `add_sequence` pack a record `digest_sequence` made (its bytes are
always raw). Rode with G1 on desk #77, which closed unanswered on 2026-11-21: not filed.
