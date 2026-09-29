# Input formats and how they are recognised

`just_vcf.sniff(path)` recognises a file from its content. Labels, extensions and upload types are
never trusted: PGP-HMS shows them wrong in both directions (PGP_HMS.md).

## How sniffing works

1. **Layers.** Compression and archives are peeled by magic bytes: gzip and BGZF (told apart by the
   `BC` extra subfield), bz2, xz, zip and tar, nested up to four deep. zstd is reported unreadable.
   Zip members under `__MACOSX/`, AppleDouble `._*` forks and `.DS_Store` are skipped. A zip is first
   checked for being an office document (OOXML, ODF, iWork).
2. **Archives with several members** are sniffed member by member, largest first (ties by name). The
   answer leads with the largest detected member and carries every member's result in `members`; a
   warning says when members hold different formats.
3. **Head.** At most `HEAD_BYTE_CAP` (1 MiB) of decompressed payload is read, stopping early after 50
   lines that do not start with `#`. Text is decoded BOM-aware (UTF-16, UTF-8), falling back to
   latin-1. A NUL in the first 8 KiB makes it binary.
4. **Match.** Binary magic first (`BCF\2` is the only binary *format*; BAM, CRAM, indexes, PDF, images,
   OLE2 and friends are `not_variant_data` with a `kind`). Then the registry, in order; the first spec
   whose detection matches wins, and later matches are listed in `also_matched`. Last, text kinds that
   hold no variant data (FASTA, FASTQ, SAM, HTML, XML, GEDCOM).
5. **Build.** Only what the header states, through the winning spec's `build_hints`. No hint means
   `build=None`; hints naming different builds also mean `None`, with a warning. Contig-length voting
   belongs to the loader, not the sniffer.

`Sniffed` is truthy only for `detected`. The other outcomes are different negatives: `unknown`
(readable, nothing matched, or empty, or binary), `unreadable` (the container could not be opened,
with `reason`) and `not_variant_data` (recognised, holds no calls, with `kind`).

## Declaring a format

A format is a `FormatSpec` (`src/just_vcf/formats/spec.py`) in a family module under
`src/just_vcf/formats/`. The spec holds everything the sniffer and the loader need: detection
patterns and a data-line shape, comment prefix, delimiter, where the column names are, columns with
dtypes and alias spellings, coordinate convention, genotype encoding and no-call spellings, build
hints, loader notes, and a synthetic `example` head.

Format ids are a `frozenset[str]` (`just_vcf.formats.FORMATS`) that grows additively; `validate_format`
checks one. `Family` (what the content is) is a closed enum.

## Detection order

Each family module exports `SPECS` (a vendor signature) and `FALLBACKS` (data shape alone, such as a
headerless rsid table). The registry is every module's `SPECS` in module order, then every module's
`FALLBACKS`, so a signature always beats a shape. Inside a module the more specific spec comes first.

The ordering test (`tests/formats/test_registry.py`) sniffs every spec's `example` and requires it to
come back as that spec, so no earlier entry may capture a later one's example. A new spec that breaks
it is registered too late or detects too loosely.
