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

## Contract gaps the loader step starts from

Found while declaring formats; none changes what `sniff` reports.

- **Columns differ by format version.** A spec declares one layout (CGI declares 2.5). Earlier
  versions drop, add or rename columns; CGI 1.x has `totalScore` where 2.x has three score columns.
  The loader should read names from the file's header and look dtypes up by name.
- **An alias can hide a vocabulary change.** CGI `varQuality` (VQHIGH/VQLOW, ≤2.2) is aliased to
  `varFilter` (PASS/VQLOW/AMBIGUOUS, 2.5); the values need mapping, not just the name.
- **0-based conversion is more than `+1`.** CGI insertions have `begin == end` and need an anchor base
  from the reference; `reference` is `=` on ref and no-call blocks, so REF comes from the FASTA.
- **`no_call` is a flat tuple of spellings.** It cannot express a partial no-call inside an allele
  (CGI `GC?`), `half` zygosity, or 1.x `no-call-rc`. These live in `genotype_encoding` and `notes`.
- **Column names outside the table.** Genos names its columns in a `##Columns=` meta line; the spec
  uses `HeaderMode.NONE` with the columns in order.
- **One delimiter per spec.** 2010-era FTDNA Affymetrix exports have a comma header, space-separated
  rows, and quoted comma rows at the end. Quoting has no field either.
- **0-based single positions.** Color positions sit one below GRCh37; `Coordinates` has only 1-based
  and 0-based half-open, so the spec says half-open and a note explains.
- **Build stated elsewhere.** Color states it in a data column; deCODEme states it in a sibling file
  inside the zip, which the archive result does not carry over to the genotype member.
- **Optional columns.** `ColumnSpec` cannot say a column appears only in some versions (Promethease
  `Repute`).
- **A section line before the columns.** Illumina Final Reports put `[Data]` on its own line, and the
  Geno 2.0 export has no column line after it. `header_prefix` does not fit; a `HeaderMode.SECTION`
  with the section name would.
- **Tab or comma.** Final Reports use either, by export; the spec says `delimiter=None` and the loader
  reads it off the line after `[Data]`.
- **Chromosome outside the row.** Oxford `.gen` imputed rows carry it only in the zip member name; Big Y
  is always chrY.
- **Text kinds look at one line.** A headerless BED cannot be named `not_variant_data` without a
  multi-line shape check, so it stays `unknown`.
- **Array or imputed content inside a VCF.** Illumina GenotypingToVCF, 23andMe conversions, Genes for
  Good and Gencove VCFs parse as `vcf`, but an absent site means unknown. `Family` has no axis for it;
  a header-based `source` detector (the VCF agent's table is in the RM2 entry) would supply it.
- **Half-calls** (`1|.`, `0/.`) are not expressible in `no_call`.
- **Variable-width rows.** The key=value table has 12 fixed fields and an annotation tail of any length.
