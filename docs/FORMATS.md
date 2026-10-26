# Input formats and how they are recognised

`just_vcf.sniff(path)` recognises a file from its content. Labels, extensions and upload types are
never trusted: PGP-HMS shows them wrong in both directions (PGP_HMS.md).

This file is the design note. The per-format spec, every declared format with its detection rules,
columns, coordinates, genotype encoding, build statements and a synthetic example, plus the broken
shapes met in real files, is [FORMAT_REGISTRY.md](FORMAT_REGISTRY.md), generated from the specs by
`scripts/format_registry.py`. Regenerate it after changing a spec; a test fails until you do.

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

`open_payload(result)` reopens a detected result's payload for a loader, replaying the same layers
and members.

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

## Which formats load

`just_vcf.LOADABLE_FORMATS` lists them. `load` refuses any other format with `UnsupportedFormatError`.

| Loader | Formats |
|---|---|
| `loaders.vcf` | `vcf`, `vcf_cgi`, `vcf_freebayes_varscan2`, `vcf_no_meta`, `bcf` |
| `loaders.rsid_table` | `23andme`, `23andme_csv`, `living_dna`, `selfdecode`, `plink_23`, `genes_for_good`, `genes_for_good_imputed`, `rsid_table`, `rsid_table_headerless` |
| `loaders.table` | `ancestrydna`, `ftdna`, `ftdna_space`, `rsid_result_csv`, `color_genotypes`, `decodeme`, `genos`, `illumina_final_report` |
| `loaders.keyvalue` | `freebayes_varscan2_keyvalue`: one facet per caller (SCHEMA Decision 8) |
| `loaders.cgi` | `cgi_var`, `cgi_master_var`: Complete Genomics loci into calls and blocks |
| `loaders.site_tables` | `casava_snp_table`, `casava_indel_table`, `soapsnp_cns`, `ftdna_bigy`, `color_variant_set`, `get_evidence_gff`: tables that state REF at each site |
| `loaders.nirvana` | `nirvana_json`: the VCF records Illumina's Annotation Engine annotated, one facet per sample |

`loaders.table` reads each file through its spec: the column names from the file (matched against the
spec's names and aliases, case-insensitively), the delimiter, the coordinate convention and the no-call
spellings. Columns the frames do not carry are listed by their file names in `facets.dropped_fields`
(deCODEme's `Variation`, a Final Report's GenCall score and intensities). A file whose columns lack what
its format needs is refused with `TableLayoutError`. MyHeritage exports load as `rsid_result_csv`
(the `ftdna` spec excludes them), with the vendor left unclaimed; no MyHeritage file was available to
check that against. A table with a sample column (a Final Report of several samples, Color) loads one
sample at a time: `load(path, sample=…)`, with `just_vcf.table_samples(path)` listing the ids in the
order they first appear, as `vcf_samples` does for a VCF.

`loaders.cgi` reads both CGI callsets into one per-haplotype shape (a `var` locus's sub-rows joined
in `begin` order, a `masterVarBeta` row's two allele columns) and makes the rows from there, so the
two files of one assembly give the same sites, blocks and allele sets; the haplotype slots can differ,
since neither file's slot order is a phase. Phase comes from `hapLink`, an id per haplotype that is a
phase set of kind `cgi_haplink`, so the two files give the same sets too. The facet's caller is
`CGAPipeline` at the file's `#SOFTWARE_VERSION`, as CGI's own VCF names it, its `callset_type` is
`wgs` where `#TITLE` says whole-genome sequencing (28 of 250 PGP-HMS files have no `#TITLE`), and its
blocks make `absence_means = blocks`. The rows are read one locus range of about
`cgi.CHUNK_ROWS` rows at a time, since every step from a row to a frame row works within one locus: a
whole PGP-HMS `var` (18.9M rows, bz2) loads in about two minutes with a peak of 3.3 GB (2026-10-26).

`loaders.site_tables` reads the tables that write a reference base beside each genotype, as a VCF
does: `source_ref` is the stated base, `alts` the called alleles other than it (sorted by string, since
none of these formats orders them), and `gt` indexes `[source_ref] + alts` at load time; `ref` still
waits for `resolve_ref`. `gt_alleles` has one entry per allele the file names, so a one-letter call
(CASAVA's `Diff`, a plain SOAPsnp base, a GET-Evidence `alleles G`, every Big Y row) is one entry: the
list length is what the file wrote, not a ploidy (SCHEMA Decision 2). Per format:

- CASAVA SNP tables: `ad` from the four base counts over `[source_ref] + alts`, `dp` from `Total Bases`;
  `Call Type` is checked against `Call` and `Ref Base`, not carried. `Sample Name` (free text, a
  person's name in real files) is never read into a frame. The rows hold SNVs only, so the facet's
  `content_filter` is `snv_only` from the rows.
- CASAVA indel tables: `ref/indel` is the reference span over the called allele with `-` padding; the
  dashes are dropped, so a deletion's ALT or an insertion's REF is the empty allele `""`. `pos` is the
  span's first base, or for an insertion the base after the insertion point, which therefore sits at
  `pos - 1` with `end_pos` null (the Genos and CGI rule). Measured on the one PGP-HMS table against
  NCBI36: `ref_upstream` ends before `pos`, and the deleted bases and `ref_downstream` follow, on every
  sampled row. `max_gtype` gives `gt`, `Q(max_gtype)` `gq`, `Q(indel)` `qual`, `ref_reads` and
  `indel_reads` `ad`.
- SOAPsnp: the IUPAC consensus is the base or the pair it names (sorted); the consensus quality is `gq`,
  depth `dp`. Three-base codes stay unread. A depth-0 row keeps `dp = 0`.
- Big Y: every row is chrY (the file names no chromosome, so `contig_map` is `[]`); Known SNP rows have
  no position and are dropped with a count; `?` is `[null]`.
- Color's variant set: `Zygosity` gives `gt` (`[0, 1]` or `[1, 1]` over `[Ref, Alt]`); the table ends at
  its disclaimer; the build comes from `ReferenceVersion` (`GRCh37.p12` reads `GRCh37`).
- GET-Evidence GFF: a `REF` row is a one-site `ref` block, so the facet's `absence_means` is `blocks`
  (reference only where a block says so). A row with `ref_allele` is indexed as above; the build-36
  layout states no REF, so its rows are array-shaped (`gt_alleles` only) until `resolve_ref`.

None of these states a build except GET-Evidence and Color, and `resolve_ref` reads no REF for a facet
that has `source_ref` and no build (a VCF without lengths or a build hint is refused the same way:
inferring the build from the rows would make `ref_check` confirm the comparison that chose it). Nor
does any name its caller; `caller` stays null.

`loaders.nirvana` reads a Nirvana JSON's position objects (the lines between the header line, which
opens `"positions":[`, and the first line starting with `]`; an OMIM or gene section after them holds
no calls) through an NDJSON working copy, as VCF rows: `refAllele` and `altAlleles` are `source_ref` and
`alts` (`[]` where a reference record has none), the chosen sample's `genotype` is the GT text, with
`genotypeQuality`, `totalDepth` and `alleleDepths` as `gq`, `dp` and `ad` (only with one depth per
allele). The trimmed `variants[]` and every annotation key are listed in `dropped_fields`. A schema-v4
header names no assembly and no contig lengths, so the build stays unknown. A file of several samples
loads one at a time (`load(path, sample=…)`).

Not loadable yet, each for a reason the file gives:

- Geno 2.0 (`genographic_geno2`, its headerless variant) and the headerless Final Report name markers
  without positions, so their rows cannot be `calls` rows until a chip manifest places them.
- `promethease_tab_genotype` is a report on the sites SNPedia knows: no chromosome, no position, and
  genotypes oriented as SNPedia orients the site, which may be the minus strand.
- `oxford_gen` (Gencove) states no chromosome in any field defined to hold one: the member name says
  it, and typed rows repeat it in a column the Oxford format defines as the input's SNP id. A decision
  for the person (layer 1); until then `load` refuses it.

Not facets: `cgi_summary` and `decodeme_info` hold no calls. A CGI summary's coverage fractions and
called sex are genome-wide figures about the `var` file beside it (provenance of that facet, layer 2),
not positions a facet could hold; deCODEme's info member states its sibling's build, which the
deCODEme loader already reads.

## Contract gaps the loader step starts from

Found while declaring formats; none changes what `sniff` reports. Where a loader has closed one, the
entry says how.

- **Columns differ by format version.** A spec declares one layout (CGI declares 2.5). Earlier
  versions drop, add or rename columns; CGI 1.x has `totalScore` where 2.x has three score columns.
  The loader should read names from the file's header and look dtypes up by name. *Closed for CGI:*
  the loader reads the `>` line and takes columns by name. The PGP-HMS census (2026-10-23) holds
  `var` in formats 1.6 and 1.7 (one each, CGA 1.11 and 1.12: twelve columns, `totalScore`, no filter
  column), 2.0 (199), 2.2 (1) and 2.5 (178, one with `intervalId`/`intervalAllele`), and
  `masterVarBeta` in 2.0 (2) and 2.5 (123). Loaded whole: a 1.6 `var` and a 2.5 `var` and
  `masterVarBeta`; synthetic files of 1.6, 1.7, 2.0 and 2.5 are tested; 2.2 and the `intervalId`
  variant read by the same by-name rule, not loaded here. A 1.x row has `sample_filters` null.
- **An alias can hide a vocabulary change.** CGI `varQuality` (VQHIGH/VQLOW, ≤2.2) is aliased to
  `varFilter` (PASS/VQLOW/AMBIGUOUS, 2.5); the values need mapping, not just the name. *Kept, not
  mapped:* `sample_filters` holds VQHIGH as written; no CGI document read here equates it with PASS.
- **0-based conversion is more than `+1`.** CGI insertions have `begin == end` and need an anchor base
  from the reference; `reference` is `=` on ref and no-call blocks, so REF comes from the FASTA.
  *Closed in the loader without a FASTA:* a span is `pos = begin + 1`, `end_pos = end`; an insertion
  point sits at `pos = begin` with `end_pos` null and `source_ref = ""`; alleles stay unpadded, `""`
  for an empty one, until RM3 adds the anchor (SCHEMA Decision 2).
- **`no_call` is a flat tuple of spellings.** It cannot express a partial no-call inside an allele
  (CGI `GC?`), `half` zygosity, or 1.x `no-call-rc`. These live in `genotype_encoding` and `notes`.
  *Closed for CGI:* a haplotype written `?` is a null in `gt_alleles`, a `?` inside a sequence is kept
  verbatim, and every no-call varType (`no-call`, `no-call-rc`, `no-call-ri`, `no-ref`,
  `PAR-called-in-X`) reads as no call.
- **Column names outside the table.** Genos names its columns in a `##Columns=` meta line; the spec
  uses `HeaderMode.NONE` with the columns in order. *Closed:* the loader reads the names from that line
  when it is there (the spec's aliases match them) and falls back to the spec's order when it is not.
- **One delimiter per spec.** 2010-era FTDNA Affymetrix exports have a comma header, space-separated
  rows, and quoted comma rows at the end. Quoting has no field either. *Closed in the loader:* an
  `ftdna_space` row is read by its own shape, and every field is stripped of surrounding double quotes.
  Seven PGP-HMS FTDNA uploads are two exports joined, the column line repeated inside; the repeat is
  skipped and counted.
- **Two layouts in one file.** One headerless upload is a 23andMe-shaped export followed by an
  export with one allele per column, the AncestryDNA data layout without its preamble (F3). *Closed in
  the loader, with no new format id:* `rsid_table` reads a run of at least 1,000 rows of one other
  width that ends the file, when its first five rows are in the two-allele layout, as a second facet
  (`relation = appended_export`, SCHEMA.md Decision 8). A run in any other layout is warned about and
  not read.
- **0-based single positions.** Color positions sit one below GRCh37; `Coordinates` has only 1-based
  and 0-based half-open, so the spec says half-open and a note explains. *Closed:* the loader adds 1
  for every 0-based spec. Genos insertion sites (`begin == end`) go to the base before them, where VCF
  places them, with `end_pos` null until RM3 reads the anchor base. Genos's `-` is read as the empty
  allele `""` (SCHEMA Decision 2), `source_gt` keeping the `-`.
- **Build stated elsewhere.** Color states it in a data column; deCODEme states it in a sibling file
  inside the zip, which the archive result does not carry over to the genotype member. *Closed:* Color's
  build is the one every row states (rows that disagree leave it unknown, with a warning), and a
  deCODEme facet takes the build its zip's `decodeme_info` member states, from `Sniffed.members`.
- **Optional columns.** `ColumnSpec` cannot say a column appears only in some versions (Promethease
  `Repute`).
- **A section line before the columns.** Illumina Final Reports put `[Data]` on its own line, and the
  Geno 2.0 export has no column line after it. *Closed for the Final Report:* `HeaderMode.SECTION` with
  `FormatSpec.section = "Data"`. The headerless export keeps `HeaderMode.NONE`.
- **Tab or comma.** Final Reports use either, by export; the spec says `delimiter=None` and the loader
  reads it off the line after `[Data]`. *Closed* that way.
- **Chromosome outside the row.** Oxford `.gen` imputed rows carry it only in the zip member name; Big Y
  is always chrY.
- **Text kinds look at one line.** A headerless BED cannot be named `not_variant_data` without a
  multi-line shape check, so it stays `unknown`.
- **Array or imputed content inside a VCF.** Illumina GenotypingToVCF, 23andMe conversions, Genes for
  Good and Gencove VCFs parse as `vcf`, but an absent site means unknown. `Family` has no axis for it.
  Closed by the loader, not by `sniff`: the `source` detector (`just_vcf/source.py`) reads the full
  header and sets the facet's `callset_type` to `array` or `imputed` where the header says so.
- **Half-calls** (`1|.`, `0/.`) are not expressible in `no_call`. In the array tables a half-call
  is a `-` (FTDNA `-A`) or a `0` (AncestryDNA) beside a letter, and the loader keeps the called allele
  in its slot: `[null, "A"]`. A whole no-call spelling (`spec.no_call`) is a token: in formats that
  write two letters on every called row (FTDNA's `---`, a Final Report's single `-`) it gives two
  nulls.
- **A column line with two `#`.** Two Veritas SV and CNV uploads write `##CHROM` and a third `#` on
  most meta lines (`###contig`, `###FORMAT`). *Closed in the loader:* where no `#CHROM` line exists, a
  meta line that is VCF's fixed column line with one `#` too many is read as the column line, and the
  file is read through a rewritten copy. A `###` line stays what VCF makes of it, a meta line whose key
  starts with `#`: no contig or FORMAT is declared from it (the build hints still read its text).
- **A VCF re-saved by a spreadsheet.** CSV quoting on a tab file (quoted fields, doubled quotes),
  meta lines padded with tabs, bare-CR line ends. *Closed in the loader, under the `vcf` id:* a payload
  whose first 64 KiB has a line starting `"##`, or CRs and no LF, is read through `csv` with universal
  newlines and polars-bio reads the rewritten copy. A separate format id was rejected: the delimiter is
  still a tab, and the quoting is a container's mangling, like CRLF.
- **Variable-width rows.** The key=value table has 12 fixed fields and an annotation tail of any length.
  *Closed in the loader:* only the twelve fixed fields are read, each checked for its key; the tail's
  keys, read from the first 1,000 rows, are listed in `facets.dropped_fields`.
