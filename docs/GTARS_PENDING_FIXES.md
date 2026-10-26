# gtars: pending fixes

Behaviour of gtars (databio/gtars, a core dependency: VRS allele ids, allele normalization, refget)
that needs an upstream change first. Each entry notes the mitigation already in place here. Filing
goes to gtars' GitHub issue tracker and needs the person's go per batch; the issue reference is
recorded here once filed. The upstream-liaison seat keeps this ledger, beside
`docs/POLARS_BIO_PENDING_FIXES.md`.

The mirror: `data/interim/gtars/`, a read-only clone of databio/gtars at the tag `uv.lock` pins
(`gtars-python-v0.10.0`, `90141b6`), for locating code paths; never edited, never filed into. Refresh
with `git -C data/interim/gtars fetch --tags` and a checkout of the new `gtars-python-v<version>` tag.
Nothing watches upstream, for the reasons the polars-bio ledger gives; the seat's start check compares
`uv.lock`'s gtars version against the pin above.

## Verdicts at a glance

| Item | Verdict | File at |
|---|---|---|
| G1 VRS ids differ from vrs-python on deletions, repeat insertions and reference alleles | normalization never returns a `ReferenceLengthExpression`; the Python API digests literal states only | databio/gtars |
| G2 `add_sequence` into an Encoded store reads back wrong bases (2026-10-26) | a record's ASCII bytes are kept as if already bit-packed, any case; fixed upstream after 0.10.0 (`03656182`, in the `gtars-python-v0.11.0` tag, not on PyPI yet) | nothing to file; re-check at the 0.11 bump |

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
has the state, `_justify` and the JSON writer shrink to calls. Not filed.

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
