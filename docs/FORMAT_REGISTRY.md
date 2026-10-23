# Format registry

Every input format just-vcf recognises, precise enough to write a sniffer for in another language, and the broken and odd shapes found in real files. The design behind the registry (why the order is what it is, how a format is declared, which contract gaps the loaders started from) is in [FORMATS.md](FORMATS.md); the frames a loaded file becomes are in [SCHEMA.md](SCHEMA.md).

Generated from the format specs in `src/just_vcf/formats/`, `just_vcf.sniff` and the loaders by `scripts/format_registry.py`: edit those, then run `uv run python scripts/format_registry.py`. A test fails when this file and the code disagree.

**Examples are synthetic.** Ids and positions are real dbSNP sites at their positions on the build the example states (or build 36 where it says so); genotypes, alleles, scores, sample ids, dates and member names are invented and describe no person. Each format's example is a test vector: the suite sniffs it and requires that format back, and no format earlier in the order. Tabs in examples are literal tab characters.

## Contents

1. [How a file is recognised](#how-a-file-is-recognised)
2. [Vocabulary](#vocabulary)
3. [Detection order](#detection-order)
4. [Formats](#formats)
5. [Broken and odd shapes met in real files](#broken-and-odd-shapes-met-in-real-files)

## How a file is recognised

Content decides, never a name, extension or upload label. The steps, in order:

1. **Peel layers by magic bytes**, reading the first 512 bytes of each layer:
   - `1f 8b`: gzip; it is **BGZF** when the FLG byte (offset 3) has FEXTRA (`0x04`) set and bytes 12-13 are `BC`.
   - `BZh`: bz2. `fd 37 7a 58 5a 00`: xz. `28 b5 2f fd`: zstd, reported `unreadable`.
   - `PK 03 04` or `PK 05 06`: zip. `ustar` at offset 257: tar.
   - At most 4 nested layers are followed. A zip inside another layer is buffered, up to 256 MiB; a tar stream is listed for at most 1 GiB.
2. **Archives.** A zip holding `[Content_Types].xml` is an office document (`ooxml`), one holding `mimetype` is `odf`, one holding `index.xml` or `Index/*.iwa` is `iwork`; all three are `not_variant_data`. Members matching `(^|/)(__MACOSX/|\._|\.DS_Store$|Thumbs\.db$)` are skipped. The rest are sniffed largest first, ties by name, at most 64; the answer leads with the largest detected member and lists every member's answer. Zip methods the reader cannot inflate are reported `unreadable` and named: 9 (Deflate64), 98 (PPMd), 99 (AES encryption).
3. **Binary payloads.** On the innermost layer, `42 43 46 02` (`BCF\2`) is the format `bcf`. These prefixes are `not_variant_data` of the named kind: `BAM\x01` bam, `CRAM` cram, `BAI\x01` bai, `TBI\x01` tbi, `CSI\x01` csi, `%PDF` pdf, `\x89PNG` png, `\xff\xd8\xff` jpeg, `GIF8` gif, `\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1` ole2, `Rar!` rar, `7z\xbc\xaf'\x1c` 7z, `bplist00` bplist, `L\x00\x00\x00\x01\x14\x02\x00` lnk, `{\rtf` rtf.
4. **Read a head**: whole lines from the start, stopping once 50 non-empty lines not starting with `#` are complete, or at 1 MiB, cut back to the last line end unless the file ended. An empty or all-whitespace head is `unknown`.
5. **Decode**: a UTF-16 BOM (`ff fe` or `fe ff`) means UTF-16. Otherwise a NUL byte in the first 8 KiB makes the payload binary, which is `unknown`. Otherwise UTF-8 with any BOM dropped (a sequence cut off in the last 4 bytes is trimmed); if that fails, latin-1.
6. **Normalise line ends**: CRLF, then bare CR, become LF. Every pattern below sees LF only.
7. **Match the registry in the order of the table below.** A spec matches when:
   - no `excludes` pattern is found anywhere in the head;
   - every `patterns` regex is found in the head, searched with `^` and `$` as line anchors (multiline) and `\A` as the start of the head;
   - and, when it has a `data_shape`, one of the first 5 data lines matches it in full. Data lines are the non-blank lines not starting with the spec's comment prefix.

   The first spec that matches is the answer; later ones that also match are reported as `also_matched`. A signature (a vendor's own text) is always tried before a fallback (data shape alone).
8. **Text that holds no calls.** If no spec matches, the first non-blank line is tried against these patterns, each naming a `not_variant_data` kind:
   - `^@(HD|SQ|RG|PG)\t`: `sam`
   - `^>`: `fasta`
   - `^@\S`: `fastq`
   - `^<!doctype html|^<html`: `html`
   - `^<\?xml`: `xml`
   - `^0 HEAD\s*$`: `gedcom`
   - `^Seq-submit ::= \{`: `genbank_asn1`
   - `^DYS393,DYS390,DYS19,`: `ftdna_ystr`
   - `^"?SnpName"?,"?Positive"?,"?Test"?\s*$`: `ftdna_ysnp_panel`
   - `^Position,rCRS,RSRS\s*$`: `mtdna_reference_table`
   - `^\s*Type : mtDNA`: `genographic_report`
   - `^Ancestry,Copy,Chromosome,Start Point,End Point`: `ancestry_segments`

   Anything else is `unknown`.
9. **Build**: every build hint of the winning spec is searched in the head (multiline). One build named: that build, with the matched text as evidence. Several builds named: no build, and a warning. None: the build is unknown. No vendor default is ever assumed; inferring a build from contig lengths or rows is the loader's work, not the sniffer's.

The four outcomes are different negatives and must stay apart: `detected`, `unknown` (readable, nothing matched, empty or binary), `unreadable` (a layer could not be opened, with a reason) and `not_variant_data` (recognised, holds no calls, with a kind).

### Regular expressions

Patterns are written in Python's `re` dialect. Besides POSIX extended syntax they use:

- `\A`: start of the head (not of a line)
- `(?i:`: a group matched case-insensitively
- `(?<!`: negative lookbehind
- `(?!`: negative lookahead; `\A(?!)` can never match
- `\b`: word boundary
- `\d`: a digit
- `\s`: whitespace
- `\S`: anything but whitespace
- `\w`: a letter, digit or underscore

## Vocabulary

**Family**, what a format's content is (a closed set):

- `genotype_table`: genotypes at predefined sites (arrays, imputation output); an absent site is unknown.
- `variant_calls`: discovery calls in VCF or a VCF-like layout.
- `callset_native`: a vendor's native whole-genome callset that states its reference and no-call regions.
- `annotated_variants`: a variant list whose purpose is annotation (reports, GFF, annotator JSON).
- `metadata`: describes a callset without holding its calls.

**Coordinates**, the convention of a format's position columns:

- `1-based`: a position names the base itself, as VCF `POS` does.
- `0-based-half-open`: `begin` counts the bases before the span and `end` is exclusive, so a single base at VCF `POS` p is `begin = p - 1, end = p`.

**Column names**, where a delimited format keeps them:

- `none`: the file names no columns; they are the ones listed below, in file order.
- `line`: the first line that is not a comment.
- `commented`: the last comment line before the data (`# rsid\tchromosome…`).
- `vcf`: VCF meta lines and a `#CHROM` line, as the VCF specification defines them.
- `section`: the first line after a `[<section>]` line.

Format ids are lowercase snake_case and form an open set: ids are added, never renamed or reused.

## Detection order

39 formats, 25 of them loadable. Detection tries them in this order.

| # | Format | Family | Vendor | Kind | Loads |
|---|---|---|---|---|---|
| 1 | [`vcf_cgi`](#vcf_cgi) | `variant_calls` | Complete Genomics | signature | yes |
| 2 | [`vcf_freebayes_varscan2`](#vcf_freebayes_varscan2) | `variant_calls` | PGP-HMS | signature | yes |
| 3 | [`vcf`](#vcf) | `variant_calls` | generic | signature | yes |
| 4 | [`vcf_no_meta`](#vcf_no_meta) | `variant_calls` | generic | signature | yes |
| 5 | [`bcf`](#bcf) | `variant_calls` | generic | signature | yes |
| 6 | [`cgi_var`](#cgi_var) | `callset_native` | Complete Genomics | signature | yes |
| 7 | [`cgi_master_var`](#cgi_master_var) | `callset_native` | Complete Genomics | signature | yes |
| 8 | [`cgi_summary`](#cgi_summary) | `metadata` | Complete Genomics | signature | no |
| 9 | [`nirvana_json`](#nirvana_json) | `annotated_variants` | Illumina | signature | no |
| 10 | [`illumina_final_report`](#illumina_final_report) | `genotype_table` | Illumina | signature | yes |
| 11 | [`illumina_final_report_headerless`](#illumina_final_report_headerless) | `genotype_table` | Illumina | signature | no |
| 12 | [`get_evidence_gff`](#get_evidence_gff) | `annotated_variants` | PGP GET-Evidence | signature | no |
| 13 | [`color_variant_set`](#color_variant_set) | `annotated_variants` | Color | signature | no |
| 14 | [`ftdna_bigy`](#ftdna_bigy) | `genotype_table` | FamilyTreeDNA | signature | no |
| 15 | [`casava_snp_table`](#casava_snp_table) | `variant_calls` | Illumina | signature | no |
| 16 | [`casava_indel_table`](#casava_indel_table) | `variant_calls` | Illumina | signature | no |
| 17 | [`ancestrydna`](#ancestrydna) | `genotype_table` | AncestryDNA | signature | yes |
| 18 | [`genos`](#genos) | `genotype_table` | Genos | signature | yes |
| 19 | [`color_genotypes`](#color_genotypes) | `genotype_table` | Color | signature | yes |
| 20 | [`decodeme`](#decodeme) | `genotype_table` | deCODEme | signature | yes |
| 21 | [`decodeme_info`](#decodeme_info) | `metadata` | deCODEme | signature | no |
| 22 | [`genographic_geno2`](#genographic_geno2) | `genotype_table` | National Geographic Genographic Project | signature | no |
| 23 | [`ftdna`](#ftdna) | `genotype_table` | FamilyTreeDNA | signature | yes |
| 24 | [`ftdna_space`](#ftdna_space) | `genotype_table` | FamilyTreeDNA | signature | yes |
| 25 | [`genes_for_good_imputed`](#genes_for_good_imputed) | `genotype_table` | Genes for Good | signature | yes |
| 26 | [`genes_for_good`](#genes_for_good) | `genotype_table` | Genes for Good | signature | yes |
| 27 | [`23andme`](#23andme) | `genotype_table` | 23andMe | signature | yes |
| 28 | [`23andme_csv`](#23andme_csv) | `genotype_table` | 23andMe | signature | yes |
| 29 | [`living_dna`](#living_dna) | `genotype_table` | Living DNA | signature | yes |
| 30 | [`selfdecode`](#selfdecode) | `genotype_table` | SelfDecode | signature | yes |
| 31 | [`plink_23`](#plink_23) | `genotype_table` | PLINK | signature | yes |
| 32 | [`promethease_tab_genotype`](#promethease_tab_genotype) | `annotated_variants` | Promethease | signature | no |
| 33 | [`freebayes_varscan2_keyvalue`](#freebayes_varscan2_keyvalue) | `variant_calls` | PGP-HMS | fallback | yes |
| 34 | [`oxford_gen`](#oxford_gen) | `genotype_table` | generic (IMPUTE2 / Oxford .gen; Gencove) | fallback | no |
| 35 | [`soapsnp_cns`](#soapsnp_cns) | `variant_calls` | BGI (SOAPsnp) | fallback | no |
| 36 | [`rsid_result_csv`](#rsid_result_csv) | `genotype_table` | unknown | fallback | yes |
| 37 | [`genographic_geno2_headerless`](#genographic_geno2_headerless) | `genotype_table` | National Geographic Genographic Project | fallback | no |
| 38 | [`rsid_table`](#rsid_table) | `genotype_table` | unknown | fallback | yes |
| 39 | [`rsid_table_headerless`](#rsid_table_headerless) | `genotype_table` | unknown | fallback | yes |

## Formats

### `vcf_cgi`

Complete Genomics `vcfBeta-*.vcf.bz2` (CGAPipeline 2.x through cgatools mkvcf).

- **Family:** `variant_calls`
- **Vendor:** Complete Genomics
- **Registry position:** 1 of 39, a signature (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.vcf`
- **Delimiter:** VCF: tab-separated, read through polars-bio
- **Comment prefix:** `#`
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** FORMAT/GT allele indices; `|` phased within FORMAT/PS, `/` unphased. Half-calls keep one allele and write `.` for the other (`1|.`, `.|1`, `0|.`, `.|0`, `1/.`, `0/.`). A bare `.` sits on `<CGA_CNVWIN>` records, which are coverage windows, not genotypes
- **No-call spellings:** `./.`, `.|.`, `.`
- **Build:** named only by a header statement matching one of:
  - `^###?reference=.*(?i:grch37|hg19|hs37d5|human_g1k_v37|g1k\.37|(?<![a-z0-9])b37(?![0-9])|ncbi ?37|build ?37)` names GRCh37
  - `^###?reference=.*(?i:grch38|hg38|hs38|GCA_000001405\.15|ncbi ?38|build ?38)` names GRCh38
  - `^###?reference=.*(?i:grch36|hg18|ncbi ?36|build ?36)` names GRCh36
  - `^###?contig=<[^>\n]*assembly=(?i:b37|hg19|grch37)[,>]` names GRCh37
  - `^###?contig=<[^>\n]*assembly=(?i:b38|hg38|grch38)[,>]` names GRCh38
  - `^###?contig=<[^>\n]*assembly=(?i:b36|hg18|ncbi36|grch36)[,>]` names GRCh36
  - `^##source_GENOME_REFERENCE=NCBI build 37$` names GRCh37
  - `^##source_GENOME_REFERENCE=NCBI build 36$` names GRCh36
  - `^##genomemap=.*(?i:ncbi ?37)` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `^##fileformat=VCFv\d`
  - `^##center=Complete Genomics$`
  - `^##source=CGAPipeline_`

**Columns:** the VCF fixed columns, then FORMAT and one column per sample.

**From real files**

- Three record classes share the file. `<CGA_CNVWIN>` ALT rows are 2 kb copy-number windows (INFO CGA_WINEND, FORMAT CGA_GP/CGA_NP/CGA_CP…, GT `.`): drop them from calls. `<CGA_NOCALL>` rows are no-called spans with INFO END. Rows with ALT `.`, GT `./.` and AN=0 are no-called reference spans with no END: their extent is the length of REF, often dozens of bases.
- Half-calls are frequent (a few thousand per 400k rows): one allele known, the other `.`. They are graded evidence, not no-calls and not homozygous.
- FILTER was `.` on every record in the first 400k rows of two files; PASS never appears. A loader that keeps only PASS drops the whole callset.
- Chromosomes have no `chr` prefix; `##contig` lines carry `assembly=B37` and an md5.
- The header states the build twice: `##source_GENOME_REFERENCE=NCBI build 37` and `##reference=…/build37.fa.bz2`.
- INFO carries CGI annotation (CGA_XR dbSNP ids, CGA_FI functional impact, CGA_RPT, CGA_PFAM); the `##source_*` lines record the annotation database versions, not tools.
- PGP-HMS holds one per CGI genome beside the native var/masterVar files of the same assembly; the var file is the fuller callset (it states reference-called spans too).
- FORMAT GL is declared `Type=Integer` and written in phred units, best genotype 0: in 126,604 of 126,604 rows with GQ in one file, GQ was the gap between the two best. The loader reads gl = GL / 10.
- FORMAT AD is declared `Number=2` and counts reads per GT slot (`1/1` writes `3,3`), not per allele: not carried.
- On male X and Y the no-called spans write a haploid `.`; polars-bio 0.36 returns that GT as null, the same as a record with no GT, so the loader reads a null GT on such a record as a no-call.
- No ALT `.` record is called reference (`0/0`) in the first 400k rows: reference-called positions are absent. Compared with the same assembly's `var` file, it covers every base `var` does not call reference, so on a whole genome `absence_means` is `reference` (SCHEMA Decision 4).

**Example**

```text
##fileformat=VCFv4.1
##fileDate=20140505
##center=Complete Genomics
##source=CGAPipeline_2.5.0.37
##source_GENOME_REFERENCE=NCBI build 37
##phasing=partial
##reference=ftp://ftp.completegenomics.com/ReferenceFiles/build37.fa.bz2
##ALT=<ID=CGA_CNVWIN,Description="Copy number analysis window">
##ALT=<ID=CGA_NOCALL,Description="No-called record">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	1	.	N	<CGA_NOCALL>	.	.	END=10000;NS=1;AN=0	GT:PS	./.:.
1	10001	.	T	<CGA_CNVWIN>	.	.	NS=1;CGA_WINEND=12000	GT	.
1	752566	rs3094315	G	A	.	.	NS=1;AN=2;AC=1	GT:PS	0|1:752566
1	752721	rs3131972	A	G	.	.	NS=1;AN=1;AC=1	GT:PS	1|.:752566
```

### `vcf_freebayes_varscan2`

PGP-HMS's hg19 callsets: freeBayes and VarScan2 run per chromosome and joined by vcf-merge.

- **Family:** `variant_calls`
- **Vendor:** PGP-HMS
- **Registry position:** 2 of 39, a signature (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.vcf`
- **Delimiter:** VCF: tab-separated, read through polars-bio
- **Comment prefix:** `#`
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** FORMAT/GT: allele indices separated by / (unphased) or | (phased); one genotype per caller, identified by INFO/SF, never by column name
- **No-call spellings:** `./.`, `.|.`, `.`
- **Build:** named only by a header statement matching one of:
  - `^###?reference=.*(?i:grch37|hg19|hs37d5|human_g1k_v37|g1k\.37|(?<![a-z0-9])b37(?![0-9])|ncbi ?37|build ?37)` names GRCh37
  - `^###?reference=.*(?i:grch38|hg38|hs38|GCA_000001405\.15|ncbi ?38|build ?38)` names GRCh38
  - `^###?reference=.*(?i:grch36|hg18|ncbi ?36|build ?36)` names GRCh36
  - `^###?contig=<[^>\n]*assembly=(?i:b37|hg19|grch37)[,>]` names GRCh37
  - `^###?contig=<[^>\n]*assembly=(?i:b38|hg38|grch38)[,>]` names GRCh38
  - `^###?contig=<[^>\n]*assembly=(?i:b36|hg18|ncbi36|grch36)[,>]` names GRCh36
  - `^##source_GENOME_REFERENCE=NCBI build 37$` names GRCh37
  - `^##source_GENOME_REFERENCE=NCBI build 36$` names GRCh36
  - `^##genomemap=.*(?i:ncbi ?37)` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `^##fileformat=VCFv\d`
  - `^##source=freeBayes`
  - `^##source=VarScan2$`
  - `^##source_\d{8}\.\d+=vcf-merge`

**Columns:** the VCF fixed columns, then FORMAT and one column per sample.

**From real files**

- Two callers, one person. INFO/SF lists the vcf-merge inputs that carry the record, indexing `##sourceFiles_<date>.1=0:…FB.vcf.gz,1:…VS.vcf.gz`: 0 is freeBayes, 1 is VarScan2. Split the two facets by SF.
- Two layouts. Some files have two sample columns, freeBayes's `unknown` then VarScan2's `Sample1`, with `.` in the column of a caller that did not call the site. Most have one data column: the header declares one sample, or declares both (`unknown`, `Sample1`) while every row carries only one. In the one-column layout the column holds freeBayes's values wherever SF has 0 (VarScan2's genotype on a shared site is lost) and VarScan2's where SF is 1 alone.
- polars-bio 0.36 reads the declared-two, carries-one layout and fills the second sample with None; that None is a missing column, not a VarScan2 no-call.
- INFO AC/AN/AF on shared sites sum the two callers as two samples (AN=4, AC=2 for one het). Never derive a genotype from INFO.
- FORMAT key order varies per record (the union of both callers' keys in hash order); look keys up by name. freeBayes keys: RO, AO, QR, QA, GL. VarScan2 keys: FREQ (a percentage string), PVAL, RBQ, ABQ, SDP, RDF, RDR, ADF, ADR.
- `##reference=/vbin/ref/hg19.fa` names the build; chromosomes carry `chr`. Rows are sorted as strings (chr10 before chr1). Some FILTER values are VarScan2's (str10, indelError).
- PGP-HMS arvados collections hold one per sequenced genome. A copy stripped of every meta line detects as `vcf_no_meta`; its rows still carry SF and both callers' keys.
- VarScan2 declares AD as Number=1 (the ALT depth only), so AD is not carried; freeBayes's GL is log10 (Number=G) and is.
- The loader returns one facet per caller (SCHEMA Decision 8, `relation = caller_split`): the caller of an index comes from the `.FB.`/`.VS.` names on `##sourceFiles`; in the one-column layout the VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the shared ones.

**Example**

```text
##fileformat=VCFv4.1
##fileDate=20160919
##source=freeBayes v0.9.21-24-g840b412-dirty
##reference=/vbin/ref/hg19.fa
##phasing=none
##source=VarScan2
##INFO=<ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##source_20160919.1=vcf-merge(r840) chr1.FB.vcf.gz chr1.VS.vcf.gz
##sourceFiles_20160919.1=0:chr1.FB.vcf.gz,1:chr1.VS.vcf.gz
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	unknown	Sample1
chr1	752566	.	G	A	425.3	PASS	AC=2;AN=4;SF=0,1	GT	0/1	0/1
chr1	752721	.	A	G	96.5	.	AC=1;AN=2;SF=0	GT	1/1	.
chr1	776546	.	A	G	.	PASS	AC=1;AN=2;SF=1	GT	.	0/1
```

### `vcf`

VCF with a `##fileformat=VCF` line, first or not.

- **Family:** `variant_calls`
- **Vendor:** generic
- **Registry position:** 3 of 39, a signature (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.vcf`
- **Delimiter:** VCF: tab-separated, read through polars-bio
- **Comment prefix:** `#`
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** FORMAT/GT: allele indices separated by / (unphased) or | (phased)
- **No-call spellings:** `./.`, `.|.`, `.`
- **Build:** named only by a header statement matching one of:
  - `^###?reference=.*(?i:grch37|hg19|hs37d5|human_g1k_v37|g1k\.37|(?<![a-z0-9])b37(?![0-9])|ncbi ?37|build ?37)` names GRCh37
  - `^###?reference=.*(?i:grch38|hg38|hs38|GCA_000001405\.15|ncbi ?38|build ?38)` names GRCh38
  - `^###?reference=.*(?i:grch36|hg18|ncbi ?36|build ?36)` names GRCh36
  - `^###?contig=<[^>\n]*assembly=(?i:b37|hg19|grch37)[,>]` names GRCh37
  - `^###?contig=<[^>\n]*assembly=(?i:b38|hg38|grch38)[,>]` names GRCh38
  - `^###?contig=<[^>\n]*assembly=(?i:b36|hg18|ncbi36|grch36)[,>]` names GRCh36
  - `^##source_GENOME_REFERENCE=NCBI build 37$` names GRCh37
  - `^##source_GENOME_REFERENCE=NCBI build 36$` names GRCh36
  - `^##genomemap=.*(?i:ncbi ?37)` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\A#`
  - `^##fileformat=VCFv\d`

**Columns:** the VCF fixed columns, then FORMAT and one column per sample.

**From real files**

- polars-bio 0.36 refuses a header whose first line is not `##fileformat` ('missing fileformat'). Genos exome exports sort their meta lines, so `##fileformat` can sit hundreds of lines down: move it to line 1 before reading.
- polars-bio 0.36 silently drops INFO and FORMAT keys the header does not define, GT included. A repaired or synthesised header must define every key the loader needs.
- Damaged headers that polars-bio rejects as 'invalid record' until repaired: a spreadsheet round-trip (every meta line tab-padded to the row width, lines containing quotes CSV-quoted as `"##…=""…"""`); SV and CNV files (uploads named as Dante Labs WGZ) with an extra `#` on meta lines (`###FORMAT`, `###contig`) and `##CHROM` for the column line. Their CNV rows write ALT `CNV` without angle brackets.
- gVCF blocks are not a format: infer them from unfiltered rows. Observed block spellings: ALT `<NON_REF>` (GATK, Helix), `<*>` (DeepVariant), ALT `.` with INFO END (Illumina Isaac/starling with BLOCKAVG_min30p3a, Sequencing.com with MinDP, Genos exome). False friends: `##ALT=<ID=NON_REF>` and `##GVCFBlock` headers left on genotyped VCFs (GATK pipelines, and uploads named as Dante Labs and Nebula), END on SV/CNV records, and END=POS on every hom-ref row of the Open Humans 23andMe conversion.
- Genos exome VCFs mark untargeted spans (FILTER and INFO `NT`) and no-called spans (FILTER `NC`) as END blocks with GT `./.`; targeted positions between blocks are absent. FORMAT adds CGIANN_* annotation keys.
- Array and imputed content also arrives as VCF, where an absent site is unknown, never hom-ref: Illumina GenotypingToVCF (`##workflow_type=Illumina_GenotypingToVCF`; unplaceable probes as chrom `NA`, POS 0, REF `N`, ALT both alleles, FILTER GTEX; FILTER NOCALL rows carry GT `./.`), 23andMe conversions (Open Humans, DNA.Land `intersection_dbSNP_23andme`), Genes for Good (ids such as `1:762320C>T` off dbSNP), and imputation (Gencove loimpute and Minimac3 with GP/DS/GL, FILTER LOWCONF).
- A header may declare more sample columns than rows carry; polars-bio 0.36 fills the missing samples with None.
- `##source=aengine` (`##ae_version`; uploads named as FamilyTreeDNA Big Y) writes `DP=…;GTL=…;QUAL=…` into FILTER on failing rows.
- The 23andMe exome VCF names no source and no reference; only its `##contig` lines (`assembly=hg19`) state the build. Its ALT `.` rows are called hom-ref.
- Sequencing.com states `##missingDataClarification=not_sequenced`: coordinates absent from the file were not sequenced, not reference.
- DeepVariant headers name no reference; the build comes from contig lengths, which is the loader's job.
- Header `##reference` paths sometimes name no build (`refData/1055/1055.fa`, `genome.fa`); the build then stays unknown here.

**Example**

```text
##fileformat=VCFv4.2
##reference=file:///ref/hs37d5.fa
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	50	PASS	.	GT	0/1
```

### `vcf_no_meta`

VCF body under a bare `#CHROM` line, with no `##` meta lines at all.

- **Family:** `variant_calls`
- **Vendor:** generic
- **Registry position:** 4 of 39, a signature (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.vcf`
- **Delimiter:** VCF: tab-separated, read through polars-bio
- **Comment prefix:** `#`
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** FORMAT/GT: allele indices separated by / (unphased) or | (phased)
- **No-call spellings:** `./.`, `.|.`, `.`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(\tFORMAT(\t[^\t\n]+)+)?\r?$`
- one of the first 5 data lines matches in full:
  - `[^\t#][^\t]*\t\d+\t[^\t]+\t[ACGTNacgtn]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+(\t[^\t]+)*`

**Columns:** the VCF fixed columns, then FORMAT and one column per sample.

**From real files**

- polars-bio 0.36 rejects the file ('invalid record'). Prepending `##fileformat` makes it readable, but every INFO and FORMAT key without a definition is then dropped silently, GT included: synthesise a definition for each key seen in the rows. The loader does this for the keys the frames carry (GT, GQ, DP, AD, PL, … and INFO END, SVLEN) and lists every other key the first rows use in `dropped_fields`; the facet's `vcf_version` stays null, since the file states none.
- No build statement is possible; the build comes from the data or from a sibling file.
- The one PGP-HMS case is a freeBayes+VarScan2 merge stripped of its meta lines (one `Sample1` column, INFO SF, both callers' FORMAT keys): split it as `vcf_freebayes_varscan2`'s notes describe.

**Example**

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	752566	.	G	A	50.2	.	AB=0.5;SF=0	GT:RO:AO	0/1:7:5
chr1	776546	.	A	G	.	PASS	ADP=28;SF=1	GT:FREQ	0/1:28.57%
```

### `bcf`

Binary VCF. Recognised by its `BCF\2` magic inside BGZF, not by a text head.

- **Family:** `variant_calls`
- **Vendor:** generic
- **Registry position:** 5 of 39, a signature, though its magic bytes are checked before any text spec (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.vcf`
- **Delimiter:** none: binary records (BCF 2), read through polars-bio
- **Comment prefix:** none
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** FORMAT/GT, binary-encoded
- **No-call spellings:** none declared
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- binary magic: after the compression layers are peeled, the payload starts with `42 43 46 02` (`BCF\2`). Its text detection `\A(?!)` never matches.

**Columns:** the VCF fixed columns, then FORMAT and one column per sample.

**From real files**

- polars-bio 0.36 reads BCF only from a path named `.bcf` (any other name fails with 'empty input'), so a BCF under another name or inside an archive is read from a copy named `input.bcf`.
- The header is the VCF header text inside the binary layout; the loader reads it itself, and the first record's contig index and 0-based POS to check polars-bio's coordinates.
- A BCF wrote no genotype text, so `source_gt` stays null; polars-bio's rendering is not the file's.
- No PGP-HMS file on disk is a BCF (2026-10-22); the tests write theirs with pysam.

**Example:** none; the format is binary.

### `cgi_var`

Complete Genomics `var` file (`#TYPE VAR-ANNOTATION`), format versions 1.6 to 2.5.

- **Family:** `callset_native`
- **Vendor:** Complete Genomics
- **Registry position:** 6 of 39, a signature (`just_vcf.formats.cgi`)
- **Loads:** yes, through `just_vcf.loaders.cgi`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment; the marker `>` before the first name is not part of it
- **Coordinates:** 0-based-half-open (`begin` counts the bases before the span and `end` is exclusive, so a single base at VCF `POS` p is `begin = p - 1, end = p`)
- **Genotype:** One row per haplotype at a called or partly called locus: `allele` 1 and 2 (sharing `locus`) each carry that haplotype's varType and alleleSeq, so the genotype is the pair of rows (a het SNP is `1 snp G A` + `2 ref G G`). Uniform blocks (ref, no-call, no-ref) are one row with `allele` `all`. The file has no zygosity column. Allele sequences: `=` means identical to the reference over the locus, `?` an unknown sequence (no call), empty a deletion of the whole span; a `?` inside a sequence (`GC?`, `?GAA?`) marks uncalled bases in a partly called allele. varType seen: ref, no-ref (reference is N, assembly gap), no-call, snp, ins, del, sub, and in format 1.x also no-call-rc (a no-call sub-kind); masterVar adds complex (alleles of different kinds, and every no-call locus outside no-ref gaps).
- **No-call spellings:** `?`
- **Build:** named only by a header statement matching one of:
  - `^#GENOME_REFERENCE\tNCBI build 37\s*$` names GRCh37
  - `^#GENOME_REFERENCE\tNCBI build 36\s*$` names GRCh36

**Recognised by**

- every one of these is found in the head:
  - `^#TYPE\tVAR-ANNOTATION\s*$`
  - `^>locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\t`
- one of the first 5 data lines matches in full:
  - `\d+\t[^\t]*\t(all|\d+)\t[^\t]+\t\d+\t\d+\t[a-z-]+(\t[^\t]*)*`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `locus` | integer |  | locus id; the rows of one locus share it |
| `ploidy` | string |  | ploidy of the locus (2; 1 expected on male X/Y); only digits seen in heads, which cover chr1 only, so read as text and cast by the loader |
| `allele` | string |  | `all` when every haplotype agrees, else `1`, `2` |
| `chromosome` | string |  | `chr1`…`chrX`, `chrY`, `chrM` |
| `begin` | integer |  | 0-based start |
| `end` | integer |  | 0-based exclusive end |
| `varType` | string |  | ref, no-ref, no-call, snp, ins, del, sub (see notes) |
| `reference` | string |  | reference sequence of the span, `=` when not spelled |
| `alleleSeq` | string |  | this haplotype's sequence; `=`, `?`, empty (see notes) |
| `varScoreVAF` | integer |  | variant allele score, VAF model; empty on blocks |
| `varScoreEAF` | integer |  | variant allele score, EAF model; empty on blocks |
| `varFilter` | string | `varQuality` | PASS, VQLOW, AMBIGUOUS, `;`-joined; format ≤2.2 wrote VQHIGH/VQLOW as varQuality |
| `hapLink` | string |  | phasing link: an integer (2.0) or `Phased_a_b_c` (2.5) |
| `xRef` | string |  | `dbsnp.<build>:rs<id>`, `;`-joined |
| `alleleFreq` | string |  | `dbsnp:<freq>`, `;`-aligned with xRef, slots may be empty |
| `alternativeCalls` | string |  | `<seq>:<score>`, `;`-joined |

**From real files**

- The column names are on the first non-comment, non-blank line and start with `>` (`>locus`): strip the `>` from the first name. A blank line separates the `#` metadata from that line.
- Metadata keys come in no fixed order (masterVar sorts them; var and summary do not), so `#FORMAT_VERSION` is not necessarily the first line. `#GENERATED_BY` appears twice in var files (cgatools, then dbsnptool).
- `#SAMPLE`, `#ASSEMBLY_ID` and `#CUSTOMER_SAMPLE_ID` are identifiers; `#CUSTOMER_SAMPLE_ID` can carry the participant's id, sometimes with a leading space.
- Every file on disk states `#GENOME_REFERENCE\tNCBI build 37`. The build 36 hint is the same statement for older deliveries and is not in the PGP-HMS census.
- `begin`/`end` are 0-based half-open. For snp, sub and del the 1-based VCF POS is `begin + 1` and the span is `end - begin` bases. An `ins` has `begin == end` (zero length): VCF needs the anchor base at 1-based position `begin`, read from the reference FASTA.
- `reference` is `=` on ref and no-call blocks (not spelled out) and empty on insertions; REF comes from the reference FASTA, never from an allele column.
- `chromosome` carries the `chr` prefix; map `chrM` to `MT` explicitly, never with a plain `^chr` strip.
- The file tiles the whole genome: ref blocks state reference, no-call and no-ref blocks state uncalled. Absent-means-reference holds only inside ref blocks, and is graded evidence there. A haplotype row with varType no-call beside a called one is a half call, not a hom call.
- Format 1.6/1.7 (12 columns) has `totalScore` where 2.x has varScoreVAF + varScoreEAF + varQuality; `totalScore` is one score and is not either of the 2.x scores.
- Format 2.0/2.2 (14 columns) has `varQuality` (VQHIGH/VQLOW) where 2.5 has `varFilter` (PASS/VQLOW/AMBIGUOUS, `;`-joined); the column is aliased, the vocabulary is not.
- `alleleFreq` and `alternativeCalls` exist from format 2.5. One 2.5 file adds `intervalId` and `intervalAllele` after them. Read the column names from the `>` line, never by position.
- A varType snp can span several reference bases on AMBIGUOUS calls (`TGCC` → `CGCC`).
- Some GET-Evidence copies are named `<assembly>.tsv.bz2` without the `var-` prefix.

**Example**

```text
#GENERATED_BY	cgatools
#SOFTWARE_VERSION	2.5.0.24
#FORMAT_VERSION	2.5
#GENERATED_BY	dbsnptool
#DBSNP_BUILD	dbSNP build 137
#GENOME_REFERENCE	NCBI build 37
#TYPE	VAR-ANNOTATION

>locus	ploidy	allele	chromosome	begin	end	varType	reference	alleleSeq	varScoreVAF	varScoreEAF	varFilter	hapLink	xRef	alleleFreq	alternativeCalls
1	2	all	chr1	0	10000	no-ref	=	?							
2	2	all	chr1	10000	10468	no-call	=	?							
3	2	1	chr1	752565	752566	snp	G	A	120	120	PASS				
3	2	2	chr1	752565	752566	ref	G	G	98	98	PASS				
4	2	all	chr1	752566	752600	ref	=	=							
```

### `cgi_master_var`

Complete Genomics `masterVarBeta` file (`#TYPE VAR-OLPL`), format versions 2.0 and 2.5.

- **Family:** `callset_native`
- **Vendor:** Complete Genomics
- **Registry position:** 7 of 39, a signature (`just_vcf.formats.cgi`)
- **Loads:** yes, through `just_vcf.loaders.cgi`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment; the marker `>` before the first name is not part of it
- **Coordinates:** 0-based-half-open (`begin` counts the bases before the span and `end` is exclusive, so a single base at VCF `POS` p is `begin = p - 1, end = p`)
- **Genotype:** One row per locus with both alleles: allele1Seq and allele2Seq, and zygosity hom, het-ref (one allele is the reference), het-alt (two different non-reference alleles), half (one allele called, the other `?`) or no-call. varType is locus-level. Allele sequences: `=` means identical to the reference over the locus, `?` an unknown sequence (no call), empty a deletion of the whole span; a `?` inside a sequence (`GC?`, `?GAA?`) marks uncalled bases in a partly called allele. varType seen: ref, no-ref (reference is N, assembly gap), no-call, snp, ins, del, sub, and in format 1.x also no-call-rc (a no-call sub-kind); masterVar adds complex (alleles of different kinds, and every no-call locus outside no-ref gaps).
- **No-call spellings:** `?`
- **Build:** named only by a header statement matching one of:
  - `^#GENOME_REFERENCE\tNCBI build 37\s*$` names GRCh37
  - `^#GENOME_REFERENCE\tNCBI build 36\s*$` names GRCh36

**Recognised by**

- every one of these is found in the head:
  - `^#TYPE\tVAR-OLPL\s*$`
  - `^>locus\tploidy\tchromosome\tbegin\tend\tzygosity\tvarType\treference\tallele1Seq\tallele2Seq\t`
- one of the first 5 data lines matches in full:
  - `\d+\t[^\t]*\t[^\t]+\t\d+\t\d+\t[a-z-]+\t[a-z-]+(\t[^\t]*)*`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `locus` | integer |  | locus id, one row each |
| `ploidy` | string |  | ploidy of the locus; text, as in var |
| `chromosome` | string |  | `chr1`…`chrX`, `chrY`, `chrM` |
| `begin` | integer |  | 0-based start |
| `end` | integer |  | 0-based exclusive end |
| `zygosity` | string |  | hom, het-ref, het-alt, half, no-call |
| `varType` | string |  | locus-level: ref, no-ref, snp, ins, del, sub, complex |
| `reference` | string |  | reference sequence of the span, `=` when not spelled |
| `allele1Seq` | string |  | first allele; `=`, `?`, empty (see notes) |
| `allele2Seq` | string |  | second allele; `?` on half calls |
| `allele1VarScoreVAF` | integer |  | allele 1 score, VAF model |
| `allele2VarScoreVAF` | integer |  | allele 2 score, VAF model |
| `allele1VarScoreEAF` | integer |  | allele 1 score, EAF model |
| `allele2VarScoreEAF` | integer |  | allele 2 score, EAF model |
| `allele1VarFilter` | string | `allele1VarQuality` | PASS, VQLOW, AMBIGUOUS, `;`-joined |
| `allele2VarFilter` | string | `allele2VarQuality` | PASS, VQLOW, AMBIGUOUS, `;`-joined |
| `allele1HapLink` | string |  | phasing link: integer or `Phased_a_b_c` |
| `allele2HapLink` | string |  | phasing link: integer or `Phased_a_b_c` |
| `allele1XRef` | string |  | `dbsnp.<build>:rs<id>`, `;`-joined |
| `allele2XRef` | string |  | `dbsnp.<build>:rs<id>`, `;`-joined |
| `allele1Freq` | string |  | `dbsnp:<freq>`, `;`-aligned with allele1XRef |
| `allele2Freq` | string |  | `dbsnp:<freq>`, `;`-aligned with allele2XRef |
| `allele1AlternativeCalls` | string |  | `<seq>:<score>`, `;`-joined |
| `allele2AlternativeCalls` | string |  | `<seq>:<score>`, `;`-joined |
| `evidenceIntervalId` | integer |  | id of the evidence interval behind the call |
| `allele1ReadCount` | integer |  | reads supporting allele 1 |
| `allele2ReadCount` | integer |  | reads supporting allele 2 |
| `referenceAlleleReadCount` | integer |  | reads supporting the reference |
| `totalReadCount` | integer |  | reads over the locus |
| `allele1Gene` | string |  | `geneId:transcript:symbol:component:impact`, `;`-joined |
| `allele2Gene` | string |  | `geneId:transcript:symbol:component:impact`, `;`-joined |
| `pfam` | string |  | `PFAM:<id>:<name>`, `;`-joined |
| `miRBaseId` | string |  | `<name>:<accession>`, `;`-joined |
| `repeatMasker` | string |  | `<repeat>:<family>:<divergence>`, `;`-joined |
| `segDupOverlap` | integer |  | segmental duplications overlapping the locus |
| `relativeCoverageDiploid` | float |  | coverage relative to the diploid baseline |
| `calledPloidy` | string |  | CNV-called ploidy: an integer or `N` |
| `relativeCoverageNondiploid` | float |  | coverage, nondiploid model |
| `calledLevel` | string |  | CNV level: a decimal (2.0) or `N` |
| `bestLAFsingle` | float |  | lesser allele fraction, best estimate |
| `lowLAFsingle` | float |  | lesser allele fraction, lower bound |
| `highLAFsingle` | float |  | lesser allele fraction, upper bound |

**From real files**

- The column names are on the first non-comment, non-blank line and start with `>` (`>locus`): strip the `>` from the first name. A blank line separates the `#` metadata from that line.
- Metadata keys come in no fixed order (masterVar sorts them; var and summary do not), so `#FORMAT_VERSION` is not necessarily the first line. `#GENERATED_BY` appears twice in var files (cgatools, then dbsnptool).
- `#SAMPLE`, `#ASSEMBLY_ID` and `#CUSTOMER_SAMPLE_ID` are identifiers; `#CUSTOMER_SAMPLE_ID` can carry the participant's id, sometimes with a leading space.
- Every file on disk states `#GENOME_REFERENCE\tNCBI build 37`. The build 36 hint is the same statement for older deliveries and is not in the PGP-HMS census.
- `begin`/`end` are 0-based half-open. For snp, sub and del the 1-based VCF POS is `begin + 1` and the span is `end - begin` bases. An `ins` has `begin == end` (zero length): VCF needs the anchor base at 1-based position `begin`, read from the reference FASTA.
- `reference` is `=` on ref and no-call blocks (not spelled out) and empty on insertions; REF comes from the reference FASTA, never from an allele column.
- `chromosome` carries the `chr` prefix; map `chrM` to `MT` explicitly, never with a plain `^chr` strip.
- Like var, the file tiles the genome; a no-call locus carries varType complex (no-ref over assembly gaps), never no-call. zygosity half is a partial call, not a hom call and not a no-call.
- Format 2.0 (35 columns) has allele1/2VarQuality (VQHIGH/VQLOW) where 2.5 has allele1/2VarFilter (aliased, vocabulary differs), and lacks allele1/2Freq, allele1/2AlternativeCalls and the three LAF columns. Read the column names from the `>` line, never by position.
- The coverage columns repeat per CNV window and are filled on no-call and no-ref loci too; calledPloidy and calledLevel mix numbers with `N` (not called).
- It restates the calls of the var file of the same assembly (`#ASSEMBLY_ID`); loading both as separate callsets double-counts them.

**Example**

```text
#CNV_DIPLOID_WINDOW_WIDTH	2000
#CNV_NONDIPLOID_WINDOW_WIDTH	100000
#DBSNP_BUILD	dbSNP build 137
#FORMAT_VERSION	2.5
#GENERATED_BY	cgatools
#GENOME_REFERENCE	NCBI build 37
#SOFTWARE_VERSION	2.5.0.44
#TYPE	VAR-OLPL

>locus	ploidy	chromosome	begin	end	zygosity	varType	reference	allele1Seq	allele2Seq	allele1VarScoreVAF	allele2VarScoreVAF	allele1VarScoreEAF	allele2VarScoreEAF	allele1VarFilter	allele2VarFilter	allele1HapLink	allele2HapLink	allele1XRef	allele2XRef	allele1Freq	allele2Freq	allele1AlternativeCalls	allele2AlternativeCalls	evidenceIntervalId	allele1ReadCount	allele2ReadCount	referenceAlleleReadCount	totalReadCount	allele1Gene	allele2Gene	pfam	miRBaseId	repeatMasker	segDupOverlap	relativeCoverageDiploid	calledPloidy	relativeCoverageNondiploid	calledLevel	bestLAFsingle	lowLAFsingle	highLAFsingle
1	2	chr1	0	10000	no-call	no-ref	=	?	?																																
2	2	chr1	10000	10468	no-call	complex	=	?	?																										0.95	N	1.02	N	0.31	0.27	0.35
3	2	chr1	752565	752566	het-ref	snp	G	A	G	120	98	120	98	PASS	PASS									4	14	12	12	26							1.01	2	0.99	N	0.48	0.44	0.52
```

### `cgi_summary`

Complete Genomics `summary` report (`#TYPE SUMMARY-REPORT`), format versions 2.0 and 2.5.

- **Family:** `metadata`
- **Vendor:** Complete Genomics
- **Registry position:** 8 of 39, a signature (`just_vcf.formats.cgi`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment; the marker `>` before the first name is not part of it
- **Coordinates:** none (the rows are not placed on a genome)
- **Build:** named only by a header statement matching one of:
  - `^#GENOME_REFERENCE\tNCBI build 37\s*$` names GRCh37
  - `^#GENOME_REFERENCE\tNCBI build 36\s*$` names GRCh36

**Recognised by**

- every one of these is found in the head:
  - `^#TYPE\tSUMMARY-REPORT\s*$`
  - `^>Category\tMetric\tValue\tCallConfidence\s*$`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `Category` | string |  | `Genome coverage`, `Genome variations`, `CNV`, `SV`… |
| `Metric` | string |  | metric name, e.g. `Fully called genome fraction` |
| `Value` | string |  | a count, a fraction or a word (`MALE`); text |
| `CallConfidence` | string |  | HIGH, ALL, or NA |

**From real files**

- The column names are on the first non-comment, non-blank line and start with `>` (`>locus`): strip the `>` from the first name. A blank line separates the `#` metadata from that line.
- Metadata keys come in no fixed order (masterVar sorts them; var and summary do not), so `#FORMAT_VERSION` is not necessarily the first line. `#GENERATED_BY` appears twice in var files (cgatools, then dbsnptool).
- `#SAMPLE`, `#ASSEMBLY_ID` and `#CUSTOMER_SAMPLE_ID` are identifiers; `#CUSTOMER_SAMPLE_ID` can carry the participant's id, sometimes with a leading space.
- Every file on disk states `#GENOME_REFERENCE\tNCBI build 37`. The build 36 hint is the same statement for older deliveries and is not in the PGP-HMS census.
- A metric often appears twice, once per CallConfidence (HIGH and ALL); key on (Category, Metric, CallConfidence).
- `Miscellaneous\tGender` is the called sex. Coverage fractions describe how much of the genome the matching var file calls, which grades its absent-means-reference evidence.
- Written uncompressed on disk, unlike var and masterVar (bz2).

**Example**

```text
#SOFTWARE_VERSION	2.5.0.24
#GENERATED_BY	ExportWorkflow&ReportEngine
#FORMAT_VERSION	2.5
#GENOME_REFERENCE	NCBI build 37
#TYPE	SUMMARY-REPORT
#DBSNP_BUILD	dbSNP build 137

>Category	Metric	Value	CallConfidence
Miscellaneous	Gender	FEMALE	NA
Genome coverage	Fully called genome fraction	0.960	HIGH
Genome coverage	Fully called genome fraction	0.960	ALL
Genome variations	SNP total count	3400000	ALL
```

### `nirvana_json`

Illumina Annotation Engine (Nirvana) JSON: a header object, then one position object per line.

- **Family:** `annotated_variants`
- **Vendor:** Illumina
- **Registry position:** 9 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** not delimited text
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** VCF GT string per sample (`0/1`, `1/1`, `1/2`, `0/0`); haploid calls are a single index (`1`, `0`), seen on chrM
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `"genomeAssembly":"(GRCh37|hg19)"` names GRCh37
  - `"genomeAssembly":"(GRCh38|hg38)"` names GRCh38

**Recognised by**

- every one of these is found in the head:
  - `\A\{"header":\{"annotator":"(Illumina Annotation Engine|Nirvana)\b`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `chromosome` | string |  | positions[].chromosome; chr-prefixed, chrM for mtDNA |
| `position` | integer |  | positions[].position; 1-based VCF POS |
| `ref_allele` | string |  | positions[].refAllele; VCF REF with anchor base |
| `alt_alleles` | list of string |  | positions[].altAlleles; VCF ALT |
| `filters` | list of string |  | positions[].filters; PASS, LowGQX… |
| `genotype` | string |  | positions[].samples[].genotype; VCF GT |
| `variant_freq` | float |  | positions[].samples[].variantFreq |
| `total_depth` | integer |  | positions[].samples[].totalDepth |
| `genotype_quality` | integer |  | positions[].samples[].genotypeQuality |
| `allele_depths` | list of integer |  | positions[].samples[].alleleDepths |

**From real files**

- Not NDJSON: line 1 is the header object plus `"positions":[`; each following line is one position object ending in `,`; the file closes with `]}`. Stream it line by line, stripping the trailing comma.
- Schema v4 (Annotation Engine 1.5) headers do not state the assembly. `dataSources[].version` values such as phyloP `hg19` describe the annotation source, not the coordinates; they are not a build hint.
- Each position carries VCF-style REF/ALT; `variants[]` repeats them trimmed (`refAllele` `-` for an insertion, `begin` shifted by one). Load from the position level, not from `variants[]`.
- Filtered records are present (`filters` other than PASS, genotypeQuality 0); keep FILTER, do not drop.

**Example**

```text
{"header":{"annotator":"Illumina Annotation Engine 1.5.1.0","creationTime":"2016-12-02 12:36:47","schemaVersion":4,"dataVersion":"84.22.36","samples":["S1"]},"positions":[
{"chromosome":"chr1","refAllele":"G","position":752566,"quality":40.0,"filters":["PASS"],"altAlleles":["A"],"samples":[{"variantFreq":0.5,"totalDepth":30,"genotypeQuality":40,"alleleDepths":[15,15],"genotype":"0/1"}],"variants":[{"altAllele":"A","refAllele":"G","begin":752566,"chromosome":"chr1","dbsnp":["rs3094315"],"end":752566,"variantType":"SNV","vid":"1:752566:A"}]},
{"chromosome":"chr1","refAllele":"A","position":752721,"quality":40.0,"filters":["PASS"],"altAlleles":["G"],"samples":[{"variantFreq":1,"totalDepth":28,"genotypeQuality":40,"alleleDepths":[0,28],"genotype":"1/1"}],"variants":[{"altAllele":"G","refAllele":"A","begin":752721,"chromosome":"chr1","dbsnp":["rs3131972"],"end":752721,"variantType":"SNV","vid":"1:752721:G"}]}
]}
```

### `illumina_final_report`

GenomeStudio Final Report: `[Header]` key/value block, then `[Data]` with a column line.

- **Family:** `genotype_table`
- **Vendor:** Illumina
- **Registry position:** 10 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** tab or comma, read off the column line
- **Comment prefix:** none
- **Column names:** the first line after `[Data]`
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** Report Wizard chooses the columns: either one two-letter field (`Alleles - Plus`) or one field per allele on a named strand (Forward, Design, Top, Plus). Forward and Design disagree by strand; only Plus is genomic plus strand. I/D on indel probes.
- **No-call spellings:** `-`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A\[Header\][\t,]*\r?$`
  - `^\[Data\][\t,]*\r?\n[^\r\n]*SNP Name[\t,]`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `snp_name` | string | `SNP Name` | probe name: an rsid or a manifest name (`1:103380393`, `CNGB3:c.991-3T>G`, `101SNP8856FG_A`) |
| `sample_id` | string | `Sample ID` | lab sample id |
| `chromosome` | string | `Chr` | 1-22, X, Y, XY (pseudoautosomal), MT; 0 = unmapped |
| `position` | integer | `MapInfo` | manifest position, 1-based |
| `genotype_plus` | string | `Alleles - Plus` | both alleles on the plus strand, one field |
| `allele1_forward` | string | `Allele1 - Forward` | allele 1, dbSNP forward strand |
| `allele2_forward` | string | `Allele2 - Forward` | allele 2, dbSNP forward strand |
| `allele1_design` | string | `Allele1 - Design` | allele 1, probe design strand |
| `allele2_design` | string | `Allele2 - Design` | allele 2, probe design strand |
| `gc_score` | float | `GC Score`, `GenCall Score` | GenCall confidence, 0-1 |
| `log_r_ratio` | float | `LogR` | intensity log R ratio; `nan` on no-calls |
| `b_allele_freq` | float | `BFreq` | B allele frequency; `nan` on no-calls |
| `x` | float | `X` | normalised X intensity |
| `y` | float | `Y` | normalised Y intensity |
| `x_raw` | integer | `Xraw` | raw X intensity |
| `y_raw` | integer | `Yraw` | raw Y intensity |

**From real files**

- Sections, not a header line: `[Header]` holds key/value lines (GSGT Version, Content = the .bpm manifest, Num SNPs…), `[Data]` is followed by the table. Skip to the line after `[Data]`.
- The delimiter is tab in the .txt exports and comma in the .csv ones; in the comma exports the `[Header]` block is padded with commas to the data width. Read the delimiter off the lines after `[Data]`.
- No build is stated. `Content` names the manifest (`GSA-24v2-0_A1.bpm`); mapping a manifest to a build is a vendor convention, not a header statement, so the build stays unknown.
- `Chr 0` marks probes the manifest leaves unmapped; `XY` is the pseudoautosomal region; `MT` is mtDNA.
- Many probe names are not rsids; `-XX` suffixes (`1:153941698-CT`) mark a second probe at one site.
- Indel probes are called `I`/`D` (`II`, `DD`, `ID`). The loader keeps them verbatim: resolving them needs the manifest's alleles and a REF (RM3).
- Files are often CRLF.
- The column set varies per export: seen are {SNP Name, Sample ID, Allele1/2 - Forward, Allele1/2 - Design, GC Score} (Omni 2.5) and {SNP Name, Chr, MapInfo, Alleles - Plus, GenCall Score, LogR, BFreq, X, Y, Xraw, Yraw} (GSA-24v2). `columns` is the union; map by name and keep the strand the file names.
- A no-call is `-` per allele (`-`/`-` in the per-allele layout) and a single `-` in `Alleles - Plus`, where every called row holds two letters, so the loader reads that `-` as two missing alleles.
- The loader needs `Chr`, `MapInfo` and `Alleles - Plus`. An export with only Forward, Design or Top alleles is refused (`TableLayoutError`): turning those into plus-strand letters needs the manifest.
- A report can hold several samples (`Sample ID`); `load(sample=…)` picks one, as for a VCF.
- Exports without Chr/MapInfo carry no coordinates: positions come from the manifest, and that manifest's build must be stated, not assumed.

**Example**

```text
[Header]
Processing Date	08/12/2018 10:10 PM
Content	GSA-24v2-0_A1.bpm
Num SNPs	665608
Num Samples	1
[Data]
SNP Name	Chr	MapInfo	Alleles - Plus	GenCall Score	LogR	BFreq	X	Y	Xraw	Yraw
rs3094315	1	752566	AG	0.9154	0.1658	0.5012	0.4518	0.4692	4302	4589
rs3131972	1	752721	GG	0.8983	0.0588	1.0	0.0235	0.8654	337	6299
rs12124819	1	776546	-	0.0	nan	nan	nan	nan	0	0
```

### `illumina_final_report_headerless`

GenomeStudio Final Report exported without a column line after `[Data]` (Genographic Geno 2.0 `NGS_iSelect` chip): sample, probe, chromosome, two alleles.

- **Family:** `genotype_table`
- **Vendor:** Illumina
- **Registry position:** 11 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** one allele letter per column; I/D on indel probes; strand not stated by the file
- **No-call spellings:** `-`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A\[Header\][\t,]*\r?$`
  - `^\[Data\],*\r?\n[^,\r\n]+,[^,\r\n]+,(0|[1-9]|1\d|2[0-2]|X|Y|XY|MT),[ACGTID-],[ACGTID-]\r?$`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `sample_id` | string |  | lab sample id |
| `snp_name` | string |  | probe name: an rsid or a Y/mt marker name |
| `chromosome` | string |  | 1-22, X, Y; 0 = unmapped (most mt probes) |
| `allele1` | string |  | allele 1; strand not stated |
| `allele2` | string |  | allele 2; strand not stated |

**From real files**

- Sections, not a header line: `[Header]` holds key/value lines (GSGT Version, Content = the .bpm manifest, Num SNPs…), `[Data]` is followed by the table. Skip to the line after `[Data]`.
- The delimiter is tab in the .txt exports and comma in the .csv ones; in the comma exports the `[Header]` block is padded with commas to the data width. Read the delimiter off the lines after `[Data]`.
- No build is stated. `Content` names the manifest (`GSA-24v2-0_A1.bpm`); mapping a manifest to a build is a vendor convention, not a header statement, so the build stays unknown.
- `Chr 0` marks probes the manifest leaves unmapped; `XY` is the pseudoautosomal region; `MT` is mtDNA.
- Many probe names are not rsids; `-XX` suffixes (`1:153941698-CT`) mark a second probe at one site.
- Indel probes are called `I`/`D` (`II`, `DD`, `ID`). The loader keeps them verbatim: resolving them needs the manifest's alleles and a REF (RM3).
- Files are often CRLF.
- The five columns are inferred from content (sample id, probe, chromosome, allele, allele); the file names none of them. There is no position column: positions and strand come from the manifest.
- The chip is a Y/mtDNA-heavy ancestry panel; many probe names are marker names (`CTS100`) that need a marker table to be placed.

**Example**

```text
[Header],,,,
GSGT Version,1.9.4,,,
Content,,NGS_iSelect_v1_15030891_2012_B-wRS-2.bpm,,
Num SNPs,154476,,,
[Data],,,,
S1_A01,rs3094315,1,A,G
S1_A01,rs3131972,1,G,G
S1_A01,rs12124819,1,-,-
```

### `get_evidence_gff`

GET-Evidence genome GFF: a `##genome-build` line, then one site per line with an `alleles` attribute.

- **Family:** `annotated_variants`
- **Vendor:** PGP GET-Evidence
- **Registry position:** 12 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** the `alleles` attribute: `X/Y` is heterozygous, one letter is homozygous for it; `ref_allele`, when present, states REF. A `REF` row (attributes `.`) states hom-ref at that site
- **No-call spellings:** ``
- **Build:** named only by a header statement matching one of:
  - `^##genome-build b?36\b` names GRCh36
  - `^##genome-build b?37\b` names GRCh37
  - `^##genome-build b?38\b` names GRCh38

**Recognised by**

- every one of these is found in the head:
  - `^##genome-build b?3[678]\b`
- one of the first 5 data lines matches in full:
  - `(chr)?([0-9]{1,2}|X|Y|M|MT)\t[^\t]+\t[^\t]+\t\d+\t\d+\t[^\t]+\t[+\-.]\t[^\t]+\t(alleles [ACGTN]*(/[ACGTN]+)*(;[^\t]*)?|\.)`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `chromosome` | string |  | GFF seqid; chr-prefixed |
| `source` | string |  | GFF source: `-`, or `CGI` (not a CGI var file) |
| `type` | string |  | GFF type: `-`, `SNP` or `REF` |
| `start` | integer |  | 1-based, closed |
| `end` | integer |  | 1-based, closed; equals start for a substitution |
| `score` | string |  | `.` |
| `strand` | string |  | `+` |
| `phase` | string |  | `.` |
| `attributes` | string |  | `alleles A/G;db_xref dbsnp:rs…;ref_allele A`, `alleles C`, or `.` on REF rows |

**From real files**

- Two layouts seen. Older (`*_build_36_substitution_variants.gff`): source and type `-`, attributes only `alleles`, no REF. Newer (`##genome-build b37`): source `CGI`, type `SNP` or `REF`, attributes `alleles`, `db_xref dbsnp:rs…` (may list several), `ref_allele`, sometimes `amino_acid` and `ucsc_trans`.
- Without `ref_allele`, REF must come from the stated build's FASTA, never from the `alleles` letters; with it, the FASTA is still the check.
- `REF` rows state hom-ref at a listed site; that is evidence for that site only, not for unlisted ones.
- One row seen with an empty `alleles ` value; treat it as a no-call, not as hom-ref.
- The build line spells the build `36` or `b37`. A source of `CGI` does not make this a CGI var file.
- Rows may be sorted lexicographically (chr10, chr11, …), not karyotypically.

**Example**

```text
##genome-build b37
# File creation date: 2019-01-15 12:00:00.000000
chr1	CGI	REF	752566	752566	.	+	.	.
chr1	CGI	SNP	752721	752721	.	+	.	alleles A/G;db_xref dbsnp:rs3131972;ref_allele A
chr1	CGI	SNP	776546	776546	.	+	.	alleles G;db_xref dbsnp:rs12124819;ref_allele A
```

### `color_variant_set`

Color Health `variant_set.csv`: the panel's reportable variants with HGVS and classification.

- **Family:** `annotated_variants`
- **Vendor:** Color
- **Registry position:** 13 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** Ref and Alt columns plus `Zygosity` (`Heterozygous`, `Homozygous`)
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `,GRCh37(\.p\d+)?(,|\r?$)` names GRCh37
  - `,GRCh38(\.p\d+)?(,|\r?$)` names GRCh38

**Recognised by**

- every one of these is found in the head:
  - `^SampleID,Chromosome,Pos(ition)?,Ref,Alt,Gene,Type,cHGVS,pHGVS,`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `sample_id` | string | `SampleID` | Color sample or kit id |
| `chromosome` | string | `Chromosome` | 1-22, X, Y; no chr prefix |
| `position` | integer | `Position`, `Pos` | 1-based |
| `ref` | string | `Ref` | reference allele as stated by Color |
| `alt` | string | `Alt` | alternate allele |
| `gene` | string | `Gene` | HGNC symbol |
| `type` | string | `Type` | `SNV` or `indel` |
| `c_hgvs` | string | `cHGVS` | cDNA HGVS |
| `p_hgvs` | string | `pHGVS` | protein HGVS; empty off-coding |
| `transcript` | string | `Transcript` | Ensembl transcript |
| `bic` | string | `BIC` | BIC nomenclature (BRCA1/2); usually empty |
| `zygosity` | string | `Zygosity` | `Heterozygous` or `Homozygous` |
| `structural_variant` | string | `StructuralVariant`, `isStructuralVariant` | `no` or empty |
| `classification` | string | `Classification` | ACMG-style class (Benign, Likely Benign…) |
| `comments` | string | `Comments` | free text; usually empty |
| `test_type` | string | `TestType`, `TestOrdered` | test or panel name(s) |
| `reference_version` | string | `ReferenceVersion` | `GRCh37` or `GRCh37.p12` |
| `analysis_completed` | date | `AnalysisCompleted` | ISO date |

**From real files**

- Two layouts. Older (16 columns): `Pos`, an `Exon` column, `isStructuralVariant`, `TestOrdered`, no BIC/Comments/AnalysisCompleted. Newer (18 columns, declared here): `Position`, `BIC`, `StructuralVariant`, `Comments`, `TestType`, `AnalysisCompleted`. Map by name.
- The build is stated per row in `ReferenceVersion`, not in a header line; the build hints read it.
- A footer follows the data: blank line(s), `Raw data disclaimer:` and a quoted paragraph that may span lines. Stop at the first row whose SampleID cell is empty or `Raw data disclaimer:`.
- Only the panel's reportable variants are listed; an absent site is unknown, not reference.
- CRLF line ends.

**Example**

```text
SampleID,Chromosome,Position,Ref,Alt,Gene,Type,cHGVS,pHGVS,Transcript,BIC,Zygosity,StructuralVariant,Classification,Comments,TestType,ReferenceVersion,AnalysisCompleted
S1,1,55523033,A,G,PCSK9,SNV,c.1026A>G,p.Gln342=,ENST00000302118,,Homozygous,,Benign,,Color Hereditary Disease Risk Test 2.0,GRCh37.p12,2025-06-23
S1,1,55524237,G,A,PCSK9,SNV,c.1420G>A,p.Val474Ile,ENST00000302118,,Heterozygous,,Benign,,Color Hereditary Disease Risk Test 2.0,GRCh37.p12,2025-06-23
```

### `ftdna_bigy`

FamilyTreeDNA Big Y results CSV: chrY sites with reference and observed base.

- **Family:** `genotype_table`
- **Vendor:** FamilyTreeDNA
- **Registry position:** 14 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** one base per site (chrY is haploid); `?` where the site was not called
- **No-call spellings:** `?`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^"?Type"?,"?Position"?,"?SNPName"?,"?Derived"?,"?OnTree"?,"?Reference"?,"?Genotype"?\r?$`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `type` | string | `Type` | `Known SNP` or `Novel Variant` |
| `position` | integer | `Position` | chrY position; empty on Known SNP rows |
| `snp_name` | string | `SNPName` | Y-SNP name (`BY1181`); empty on novel rows |
| `derived` | string | `Derived` | `Yes (+)`, `No (-)`, `?`; empty on novel rows |
| `on_tree` | string | `OnTree` | `Yes`/`No`; empty on novel rows |
| `reference` | string | `Reference` | reference base as stated by FTDNA |
| `genotype` | string | `Genotype` | one base (haploid) or `?` |

**From real files**

- There is no chromosome column: every row is chrY, a structural fact of the product.
- Known SNP rows have no position; they are placed only through a Y-SNP name table. Novel Variant rows have a position and mostly a genotype equal to the reference.
- Neither build nor coordinate convention is stated; Big Y moved from hg19 to hg38 over its lifetime, so the build stays unknown.
- Every field is double-quoted in some exports and bare in others.

**Example**

```text
"Type","Position","SNPName","Derived","OnTree","Reference","Genotype"
"Novel Variant","10061836","","","","A","A"
"Known SNP","","A14825","?","Yes","A","?"
"Novel Variant","10061838","","","","A","?"
```

### `casava_snp_table`

Illumina CASAVA-era SNP calls with per-base read counts, one file per chromosome.

- **Family:** `variant_calls`
- **Vendor:** Illumina
- **Registry position:** 15 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** `Call` letters with `Call Type`: `Diff` = one letter, homozygous non-reference; `Het1` = two letters, first equals Ref Base; `Het2` = two letters, second equals Ref Base; `Other` = two letters, neither is the reference
- **No-call spellings:** none declared
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^Index\tSample Name\tChr Id\tPosition\tA Bases\tC Bases\tG Bases\tT Bases\tCall\t`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `index` | integer | `Index` | row index across the export |
| `sample_name` | string | `Sample Name` | free text, often a person's name; drop on load |
| `chromosome` | string | `Chr Id` | 1-22, X, Y, M; no chr prefix |
| `position` | integer | `Position` | 1-based |
| `a_bases` | integer | `A Bases` | reads calling A |
| `c_bases` | integer | `C Bases` | reads calling C |
| `g_bases` | integer | `G Bases` | reads calling G |
| `t_bases` | integer | `T Bases` | reads calling T |
| `call` | string | `Call` | one letter (hom) or two (het) |
| `total_bases` | integer | `Total Bases` | depth |
| `bases_used` | integer | `Bases Used` | depth used for the call |
| `call1_score` | float | `Call 1 Score` | score of the best call |
| `call2_score` | float | `Call 2 Score` | score of the second call |
| `ref_base` | string | `Ref Base` | reference base as stated |
| `call_type` | string | `Call Type` | `Diff`, `Het1`, `Het2`, `Other` |
| `dbsnp` | string | `dbSNP SNP Name` | rsid or empty |
| `genes` | string | `Genes` | UniProt accession(s) or empty |

**From real files**

- Only non-reference calls are listed; there are no hom-ref or no-call rows, so an absent site is unknown.
- `Sample Name` holds free text, in real files a person's name. Drop it on load; never log it.
- `Chr Id` spells mtDNA `M`; map it to `MT` explicitly.
- No build is stated; the producer is inferred (CASAVA read-count calls), the file does not name it.
- CRLF line ends. `Genes` and `dbSNP SNP Name` are often empty, leaving trailing tabs.

**Example**

```text
Index	Sample Name	Chr Id	Position	A Bases	C Bases	G Bases	T Bases	Call	Total Bases	Bases Used	Call 1 Score	Call 2 Score	Ref Base	Call Type	dbSNP SNP Name	Genes
1	S1	1	752566	14	0	12	0	GA	26	26	120.5	98.2	G	Het1	rs3094315	
2	S1	1	752721	0	0	31	0	G	31	31	210.4	0.00	A	Diff	rs3131972	
```

### `casava_indel_table`

Illumina CASAVA `indels.txt` columns, exported with Index and Sample columns in front.

- **Family:** `variant_calls`
- **Vendor:** Illumina
- **Registry position:** 16 of 39, a signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** `max_gtype` `het`/`hom` for the event in `ref/indel`; `-` marks the empty side
- **No-call spellings:** none declared
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^Index\tSample\tChr\tpos\tCIGAR\tref_upstream\tref/indel\tref_downstream\tQ\(indel\)\tmax_gtype\t`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `index` | integer | `Index` | row index |
| `sample_name` | string | `Sample` | free text, often a person's name; drop on load |
| `chromosome` | string | `Chr` | 1-22, X, Y, M; no chr prefix |
| `position` | integer | `pos` | CASAVA indel position |
| `cigar` | string | `CIGAR` | `1D`, `3I`… |
| `ref_upstream` | string | `ref_upstream` | 10 reference bases before the event |
| `ref_indel` | string | `ref/indel` | `G/-` deletion, `-/T` insertion |
| `ref_downstream` | string | `ref_downstream` | 10 reference bases after the event |
| `q_indel` | integer | `Q(indel)` | indel quality |
| `max_gtype` | string | `max_gtype` | `het` or `hom` |
| `q_max_gtype` | integer | `Q(max_gtype)` | genotype quality |
| `max2_gtype` | string | `max2_gtype` | second-best genotype: `het`, `hom`, `ref` |
| `bp1_reads` | integer | `bp1_reads` | reads spanning the left breakpoint |
| `ref_reads` | integer | `ref_reads` | reads supporting the reference |
| `indel_reads` | integer | `indel_reads` | reads supporting the indel |
| `other_reads` | integer | `other_reads` | other reads |
| `repeat_unit` | string | `repeat_unit` | repeat unit at the site |
| `ref_repeat_count` | integer | `ref_repeat_count` | repeat count in the reference |
| `indel_repeat_count` | integer | `indel_repeat_count` | repeat count with the indel |

**From real files**

- Alleles have no anchor base (`G/-`, `-/T`). A VCF record needs the anchor from the reference FASTA; `ref_upstream` gives 10 bases the loader can check the FASTA and the position convention against.
- Only called indels are listed; an absent site is unknown.
- `Sample` holds free text, in real files a person's name. Drop it on load; never log it.
- No build is stated. `Chr` spells mtDNA `M`; map it to `MT` explicitly. CRLF line ends.

**Example**

```text
Index	Sample	Chr	pos	CIGAR	ref_upstream	ref/indel	ref_downstream	Q(indel)	max_gtype	Q(max_gtype)	max2_gtype	bp1_reads	ref_reads	indel_reads	other_reads	repeat_unit	ref_repeat_count	indel_repeat_count
0	S1	20	10000012	1D	GTTTACTACA	G/-	CTTTGTAGTA	3070	hom	76	het	33	0	23	10	G	1	0
1	S1	20	10000058	1I	CTCTAACTTT	-/T	GTTTGTTTTT	2887	het	33	hom	17	9	11	6	T	3	4
```

### `ancestrydna`

AncestryDNA raw data download, array versions V1.0 and V2.0.

- **Family:** `genotype_table`
- **Vendor:** AncestryDNA
- **Registry position:** 17 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two columns, one allele letter each, forward (+) strand per the preamble; I/D on indel probes, which the loader keeps verbatim: resolving them needs the chip manifest's alleles and a REF (RM3); haploid sites (male 23/24, 26) are written as the same letter twice, never as one letter, and the loader keeps both letters: the vendor's convention, not a homozygous call
- **No-call spellings:** `0`
- **Build:** named only by a header statement matching one of:
  - `human reference build 37\.1 coordinates` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\A#AncestryDNA raw data download`
- one of the first 5 data lines matches in full:
  - `\S+\t\d{1,2}\t\d+\t[ACGTDI0]\t[ACGTDI0]`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid where possible, else a vendor probe id |
| `chromosome` | string |  | numeric only: 1-22, then 23, 24, 25 and 26 for the non-autosomes |
| `position` | integer |  | 1-based position on the stated build |
| `allele1` | string |  | first observed allele |
| `allele2` | string |  | second observed allele |

**From real files**

- Chromosomes 23, 24, 25 and 26 are X, Y, the X/Y pseudoautosomal region and MT by AncestryDNA's documentation; the file itself does not say so. Males show homozygous-looking 23 and 24 calls.
- A no-call is `0` in both allele columns; a single `0` beside a letter was not observed.
- Line endings can be mixed within one file (some preamble lines LF, the rest CRLF): strip `\r` per line rather than choosing a line terminator per file.
- The column header `rsid\tchromosome\tposition\tallele1\tallele2` is the first non-comment line.
- Allele order within a heterozygote is not sorted (both `A G` and `G A` occur).

**Example**

```text
#AncestryDNA raw data download
#This file was generated by AncestryDNA at: 01/01/2016 00:00:00 UTC
#Data was collected using AncestryDNA array version: V1.0
#Data is formatted using AncestryDNA converter version: V1.0
#Genetic data is provided below as five TAB delimited columns.  Each line 
#of the SNP using human reference build 37.1 coordinates.  Columns four and five 
rsid	chromosome	position	allele1	allele2
rs4477212	1	82154	A	A
rs3131972	1	752721	A	G
rs12562034	1	768448	0	0
rs11240777	1	798959	G	G
```

### `genos`

Genos exome genotype TSV (`##FileFormat=Genos`) at a request list of sites, including the file Genos labelled `Promethease_format`, which has the same layout.

- **Family:** `genotype_table`
- **Vendor:** Genos
- **Registry position:** 18 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 0-based-half-open (`begin` counts the bases before the span and `end` is exclusive, so a single base at VCF `POS` p is `begin = p - 1, end = p`)
- **Genotype:** alleles joined by `/`, plus strand in reference orientation; an allele may be several bases (indels, MNPs); `-` is an empty allele, used at zero-length insertion sites (begin == end); a single allele with no `/` is a hemizygous call on chrX/chrY
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `^##ReferenceGenome=UCSC\.hg19\r?$` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\A##FileFormat=Genos\r?$`
- one of the first 5 data lines matches in full:
  - `chr([1-9]|1\d|2[0-2]|X|Y|M)\t\d+\t\d+\t\S+\t\S+`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `chromosome` | string | `CHROM` | chr-prefixed: chr1-chr22, chrX, chrY |
| `begin` | integer | `BEGIN` | 0-based start of the reference span |
| `end` | integer | `END` | exclusive end; equals begin at an insertion site |
| `rsid` | string | `ID` | rsid of the requested site |
| `genotype` | string | `GENOTYPE` | alleles joined by `/`, or one allele |

**From real files**

- Column names are stated in the `##Columns=CHROM\tBEGIN\tEND\tID\tGENOTYPE` meta line, not in a header line; `columns` gives them in file order.
- `-/-` is NOT a no-call: it is homozygous for the empty allele (no insertion) at a site whose reference span is empty. No no-call spelling was observed.
- BEGIN is 0-based and END exclusive: for an SNV the 1-based position is END (checked against GRCh37 array positions of the same rsids).
- Rows are only the sites on the `##RequestList=` (SNPedia_jan2017 in every file seen), not the whole exome; the file does not say whether a hom-ref row was actually covered.
- `##Sample=` holds a kit number.

**Example**

```text
##FileFormat=Genos
##DataVersion=1.0.0
##ReferenceGenome=UCSC.hg19
##RequestList=SNPedia_jan2017
##Columns=CHROM	BEGIN	END	ID	GENOTYPE
chr1	935070	935071	rs139739199	C/T
chr1	949607	949608	rs1921	A/G
chr1	949699	949699	rs672601345	-/-
chr1	1147421	1147422	rs17568	T/T
```

### `color_genotypes`

Color `discovery_genotypes.csv`: a handful of trait SNPs with a per-row reference version.

- **Family:** `genotype_table`
- **Vendor:** Color
- **Registry position:** 19 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 0-based-half-open (`begin` counts the bases before the span and `end` is exclusive, so a single base at VCF `POS` p is `begin = p - 1, end = p`)
- **Genotype:** two allele letters; strand not stated
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `^\w+,\w+,\d+,[^,]*,[^,]*,GRCh37\r?$` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\ASampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion\r?$`
- one of the first 5 data lines matches in full:
  - `\w+,([1-9]|1\d|2[0-2]|X|Y|MT?),\d+,rs\d+,[ACGTDI-]{1,2},GRCh\d+`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `sample_id` | string | `SampleID` | Color kit id: all digits or letter-prefixed |
| `chromosome` | string | `Chromosome` | 1-22, no chr prefix |
| `position` | integer | `Position` | 0-based position |
| `rsid` | string | `RSID` | rsid |
| `genotype` | string | `Genotype` | two allele letters |
| `reference_version` | string | `ReferenceVersion` | the build, stated on every row |

**From real files**

- Position is 0-based: every row checked sits exactly one below the GRCh37 1-based position of the same rsid (e.g. rs4988235 at 2:136608645, rs671 at 12:112241765). Add 1 on load.
- The build is stated per row in ReferenceVersion, not in a header.
- The table is followed by trailer rows: `Raw data disclaimer:,,,,,` and a quoted paragraph containing commas. Stop at the first row whose Position is empty.
- A sibling `variant_set.csv` (SampleID,Chromosome,Pos,Ref,Alt,Gene,…) is an annotated clinical variant list, not this format.

**Example**

```text
SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion
10000000000000,2,136608645,rs4988235,AG,GRCh37
10000000000000,7,141673344,rs713598,CG,GRCh37
10000000000000,12,112241765,rs671,GG,GRCh37
10000000000000,16,48258197,rs17822931,CT,GRCh37
```

### `decodeme`

deCODEme genotype CSV, the data member of a deCODEme raw-data zip.

- **Family:** `genotype_table`
- **Vendor:** deCODEme
- **Registry position:** 20 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two allele letters on the strand given in Strand, the same strand as Variation: on `-` rows the loader complements both letters to get the plus strand (rs3094315 is written `C/T`, `-`, where plus-strand arrays report A/G)
- **No-call spellings:** `--`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\AName,Variation,Chromosome,Position,Strand,YourCode\r?$`
- one of the first 5 data lines matches in full:
  - `[^,]+,[ACGTDI-]+(/[ACGTDI-]+)+,\w+,\d+,[+-],(--|[ACGTDI]{1,2})`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string | `Name` | marker name, rsid where possible |
| `variation` | string | `Variation` | the probe's alleles, `A/G` |
| `chromosome` | string | `Chromosome` | no chr prefix |
| `position` | integer | `Position` | 1-based position |
| `strand` | string | `Strand` | `+` or `-` |
| `genotype` | string | `YourCode` | two allele letters |

**From real files**

- The build is not in this member. The zip's small `.txt` member (format `decodeme_info`) states `Reference sequence,NCBI Build 36 …`; every file seen is build 36.
- A third `.txt` member holds only `Signature,` and `Checksum,` lines.
- Sex-chromosome and mitochondrial chromosome spellings were not visible in the heads read.

**Example**

```text
Name,Variation,Chromosome,Position,Strand,YourCode
rs4477212,A/G,1,72017,+,AG
rs2185539,C/T,1,556738,+,TT
rs3094315,C/T,1,742429,-,CC
rs3131972,C/T,1,742584,-,--
```

### `decodeme_info`

The key,value info member of a deCODEme raw-data zip: format version, date, SNP count, build.

- **Family:** `metadata`
- **Vendor:** deCODEme
- **Registry position:** 21 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** none (the rows are not placed on a genome)
- **Build:** named only by a header statement matching one of:
  - `^Reference sequence,NCBI Build 36\b` names GRCh36

**Recognised by**

- every one of these is found in the head:
  - `\AFormat version,\d+\r?$`
  - `^SNP count,\d+\r?$`
  - `^Reference sequence,`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `key` | string |  | `Format version`, `User name`, `Date`, … |
| `value` | string |  | the value, unquoted |

**From real files**

- `User name,` carries the account name, often a real name: never log or copy it.
- Holds no calls; it states the build for the `decodeme` member beside it in the zip.

**Example**

```text
Format version,1
Date,Jan  1 2010
SNP count,1013349
Reference sequence,NCBI Build 36 and Cambridge reference sequence
Public SNP source,dbSNP 128
```

### `genographic_geno2`

Genographic Geno 2.0 results CSV (`SNP,Chr,Allele1,Allele2`), as `.all`, `.autosomal`, `.mtdna` and `.ychromo` files. No position column.

- **Family:** `genotype_table`
- **Vendor:** National Geographic Genographic Project
- **Registry position:** 22 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** none (the rows are not placed on a genome)
- **Genotype:** two columns, one letter each; I/D on indel probes; Y and Mt calls written as the same letter twice; strand not stated
- **No-call spellings:** `-`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\ASNP,Chr,Allele1,Allele2\r?$`
- one of the first 5 data lines matches in full:
  - `[^,]+,([1-9]|1\d|2[0-2]|X|Y|Mt),[ACGTDI-],[ACGTDI-]`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `snp` | string | `SNP` | marker name: rs/kgp ids, Y-SNP names, mt numbers |
| `chromosome` | string | `Chr` | 1-22, X, Y, Mt |
| `allele1` | string | `Allele1` | first allele |
| `allele2` | string | `Allele2` | second allele |

**From real files**

- There is no position column: sites must be placed through the chip manifest by marker name, and the build is whatever that manifest says. The file states none.
- Y rows are named by Y-SNP (CTS…, M…), autosomal rows mostly by kgp ids, and Mt rows by bare numbers that look like mtDNA positions; the same Mt name can repeat on consecutive rows.
- The header line ends LF while the data rows end CRLF.
- Chromosome `Mt`, not `MT`.

**Example**

```text
SNP,Chr,Allele1,Allele2
kgp10004422,12,A,G
CTS100,Y,T,T
CTS10004,Y,-,-
195,Mt,C,C
```

### `ftdna`

FamilyTreeDNA Family Finder raw data CSV (autosomal or X), Affymetrix and Illumina eras.

- **Family:** `genotype_table`
- **Vendor:** FamilyTreeDNA
- **Registry position:** 23 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two plus-strand allele letters in one field; one letter on some hemizygous X sites; `-`+letter where one allele is missing (Affymetrix era)
- **No-call spellings:** `--`, `00`, `---`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\ARSID,CHROMOSOME,POSITION,RESULT\r?$`
- none of these is found:
  - `MyHeritage`
- one of the first 5 data lines matches in full:
  - `"?rs\d+"?,"?([1-9]|1\d|2[0-2]|X|Y|MT)"?,"?\d+"?,"?[ACGTDI0-]{1,3}"?`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string | `RSID` | rsid or a vendor probe id |
| `chromosome` | string | `CHROMOSOME` | 1-22, X, Y, MT; no chr prefix |
| `position` | integer | `POSITION` | 1-based position on an unstated build |
| `genotype` | string | `RESULT` | two allele letters, or one |

**From real files**

- Data rows are usually fully double-quoted (`"rs3094315","1","752566","AA"`) under an unquoted header; some exports are unquoted. Both LF and CRLF occur.
- No build is stated anywhere in the file. FTDNA shipped both build 36 and build 37 exports (rs3094315 at 1:742429 vs 1:752566); only the file name (`o36`/`o37`) or positions tell them apart.
- Autosomes and X come as separate files; the X file holds only X rows.
- Hemizygous X calls are written as one letter in some X files and as the same letter twice in others.
- Affymetrix-era exports (2010-2011) also carry `-A`, `-C`, `-G`, `-T` (one allele missing) and `---`; neither is documented. `---` is treated as a no-call; a `-`+letter call must not be read as a homozygote or a heterozygote.
- No comment lines: a file with a `#` preamble above this header is not an FTDNA export.

**Example**

```text
RSID,CHROMOSOME,POSITION,RESULT
"rs4477212","1","82154","AG"
"rs3094315","1","752566","--"
"rs3131972","1","752721","AA"
"rs12562034","1","768448","AG"
```

### `ftdna_space`

FamilyTreeDNA Affymetrix-era autosomal export whose comma header sits over space-separated rows.

- **Family:** `genotype_table`
- **Vendor:** FamilyTreeDNA
- **Registry position:** 24 of 39, a signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** single space
- **Comment prefix:** none
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two plus-strand allele letters in one field; `-`+letter where one allele is missing
- **No-call spellings:** `--`, `---`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\ARSID,CHROMOSOME,POSITION,RESULT\r?$`
- none of these is found:
  - `MyHeritage`
- one of the first 5 data lines matches in full:
  - `rs\d+ ([1-9]|1\d|2[0-2]|X|Y|MT) \d+ [ACGT-]{1,3}`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string | `RSID` | rsid or a vendor probe id |
| `chromosome` | string | `CHROMOSOME` | 1-22, X, Y, MT; no chr prefix |
| `position` | integer | `POSITION` | 1-based position on an unstated build |
| `genotype` | string | `RESULT` | two allele letters, or one |

**From real files**

- The header is comma-separated but the rows are separated by single spaces: read the header with `,` and the body with ` `.
- The delimiter changes inside the file: the last few hundred rows (end of chromosome 22) are quoted comma rows as in `ftdna`. Parse each row by its own shape.
- CRLF line endings.
- No build is stated anywhere in the file. FTDNA shipped both build 36 and build 37 exports (rs3094315 at 1:742429 vs 1:752566); only the file name (`o36`/`o37`) or positions tell them apart.
- Autosomes and X come as separate files; the X file holds only X rows.
- Hemizygous X calls are written as one letter in some X files and as the same letter twice in others.
- Affymetrix-era exports (2010-2011) also carry `-A`, `-C`, `-G`, `-T` (one allele missing) and `---`; neither is documented. `---` is treated as a no-call; a `-`+letter call must not be read as a homozygote or a heterozygote.
- No comment lines: a file with a `#` preamble above this header is not an FTDNA export.

**Example**

```text
RSID,CHROMOSOME,POSITION,RESULT
rs4477212 1 72017 AG
rs3094315 1 742429 --
rs3131972 1 742584 GG
rs12562034 1 758311 -A
```

### `genes_for_good_imputed`

Genes for Good imputed genotypes in the 23andMe layout (`*_imputed_genotypes_noY_noMT_23andMe`), v1.0 to v1.2.

- **Family:** `genotype_table`
- **Vendor:** Genes for Good
- **Registry position:** 25 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** imputed best-guess calls, not measured ones: two allele letters; one letter on male X; I/D codes on indels and structural sites (II, DD, ID, DI) from v1.2
- **No-call spellings:** `--`, `-`
- **Build:** named only by a header statement matching one of:
  - `^# Human genome reference used: GRCh37` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\A[^\n#]?# Genes for Good v\d+(\.\d+)*: filtered, imputed`
- one of the first 5 data lines matches in full:
  - `\S+\t(\d{1,2}|X|Y|MT)\t\d+\t[ACGTDI-]{1,2}`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid; sites without one are `chrom:posREF>ALT` (v1.0 imputed) or `chrom:pos_REF/ALT` (v1.2), where ALT may be `-` or a symbolic allele such as `<CN0>` or `<INS:ME:ALU>` |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, one on haploid sites, I/D codes, `-` for a missing allele |

**From real files**

- Every export is written by PLINK `--recode 23`, so line 2 is PLINK's own preamble (`generated by PLINK at:` in v1.0 imputed, `on:` elsewhere); the vendor line is line 1.
- One real v1.0 unphased file has a stray byte before the `#` of line 1 (`I# Genes for Good`); a loader that skips lines by comment prefix must drop that line too.
- CRLF line endings in most exports; v1.0 imputed uses LF.
- The column line is `# rsid\tchromosome\tposition\tgenotype`, as in 23andMe.
- Imputed, so every call is an estimate: the loader must grade it below a measured call.
- `noY_noMT`: the file stops after chromosome X; Y and MT are absent, which is unknown, not reference.
- v1.0 imputed headers call the calls phased, but the genotype column carries no phase.
- Around 7 million sites and 200 MB uncompressed.

**Example**

```text
# Genes for Good v1.2: filtered, imputed genotype calls. For more information, go to: https://genesforgood.sph.umich.edu/readme/readme1.2.txt
# This data file generated by PLINK on: 11/29/16
#
# Below is a text version of your data.  Fields are TAB-separated.
# Human genome reference used: GRCh37/Mito:rCRS
#
# rsid	chromosome	position	genotype
rs75333668	1	762320	CT
rs11240777	1	798959	AG
1:861349_C/T	1	861349	CC
rs4475691	1	846808	CC
```

### `genes_for_good`

Genes for Good measured (unphased) genotypes in the 23andMe layout, v1.0 to v1.2.1.

- **Family:** `genotype_table`
- **Vendor:** Genes for Good
- **Registry position:** 26 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two allele letters; one letter on haploid sites (male X/Y, MT); I/D codes (II, DD, ID, DI) on indel probes from v1.2; a half-call keeps its missing allele as `-` (`T-`)
- **No-call spellings:** `--`, `-`
- **Build:** named only by a header statement matching one of:
  - `^# Human genome reference used: GRCh37` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `\A[^\n#]?# Genes for Good v\d+(\.\d+)*: filtered, unphased`
- one of the first 5 data lines matches in full:
  - `\S+\t(\d{1,2}|X|Y|MT)\t\d+\t[ACGTDI-]{1,2}`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid; sites without one are `chrom:posREF>ALT` (v1.0 imputed) or `chrom:pos_REF/ALT` (v1.2), where ALT may be `-` or a symbolic allele such as `<CN0>` or `<INS:ME:ALU>` |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, one on haploid sites, I/D codes, `-` for a missing allele |

**From real files**

- Every export is written by PLINK `--recode 23`, so line 2 is PLINK's own preamble (`generated by PLINK at:` in v1.0 imputed, `on:` elsewhere); the vendor line is line 1.
- One real v1.0 unphased file has a stray byte before the `#` of line 1 (`I# Genes for Good`); a loader that skips lines by comment prefix must drop that line too.
- CRLF line endings in most exports; v1.0 imputed uses LF.
- The column line is `# rsid\tchromosome\tposition\tgenotype`, as in 23andMe.
- v1.0 names sites without an rsid `chrom:posREF>ALT`; v1.1+ uses rsids for the same sites.
- A single `-` is the no-call on a haploid site; `--` on a diploid one.

**Example**

```text
# Genes for Good v1.0: filtered, unphased genotype calls. For more information, go to: https://genesforgood.sph.umich.edu/readme/readme1.0.txt
# This data file generated by PLINK on: 01/06/16
#
# Below is a text version of your data.  Fields are TAB-separated.
# Human genome reference used: GRCh37/Mito:rCRS
#
# rsid	chromosome	position	genotype
1:762320C>T	1	762320	CT
rs11240777	1	798959	AA
rs4475691	1	846808	T-
rs13302982	1	861808	--
```

### `23andme`

23andMe raw data download, all chip versions (v1 to v5).

- **Family:** `genotype_table`
- **Vendor:** 23andMe
- **Registry position:** 27 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** plus-strand allele letters; one letter on haploid sites (male X/Y, MT); D/I on indel probes, which the loader keeps verbatim: resolving them needs the chip manifest's alleles and a REF (RM3)
- **No-call spellings:** `--`
- **Build:** named only by a header statement matching one of:
  - `reference human assembly build 36` names GRCh36
  - `reference human assembly build 37` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `^# This data file (is )?generated by 23andMe`
- one of the first 5 data lines matches in full:
  - `(rs|i)\d+\t(\d{1,2}|X|Y|MT)\t\d+\t([ACGTDI]{1,2}|--)`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid, or a vendor-internal id such as i713426 |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, or one on haploid sites |

**From real files**

- 2011 exports are build 36; the header says so, and the loader must not assume 37.
- Exports from 2026 open with `# file_id:`, `# signature:` and `# timestamp:` lines and say `This data file is generated by 23andMe`.
- Y rows can all be present and all `--` (one real 2026 file has no Y call at all); they are not calls.

**Example**

```text
# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018
#
# We are using reference human assembly build 37 (also known as Annotation Release 104).
#
# rsid	chromosome	position	genotype
rs548049170	1	69869	TT
rs13328684	1	74792	--
i713426	1	726912	AG
rs3094315	1	752566	AA
```

### `23andme_csv`

A 23andMe download re-saved by a spreadsheet as comma-separated text.

- **Family:** `genotype_table`
- **Vendor:** 23andMe
- **Registry position:** 28 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** comma
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** as 23andMe: allele letters, one on haploid sites, D/I on indel probes
- **No-call spellings:** `--`
- **Build:** named only by a header statement matching one of:
  - `reference human assembly build 36` names GRCh36
  - `reference human assembly build 37` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `^\"?# This data file (is )?generated by 23andMe[^\n]*,,,\r?$`
  - `^# rsid,chromosome,position,genotype\r?$`
  - `^(rs|i)\d+,(\d{1,2}|X|Y|MT),\d+,([ACGTDI]{1,2}|--)\r?$`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid, or a vendor-internal id such as i713426 |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, or one on haploid sites |

**From real files**

- Every preamble line is padded to four fields with `,,,`, and preamble lines that contain a comma are quoted, so they start with `"#`: skip both kinds as comments.
- The one real file stops mid chromosome 2 after about 65,500 lines: the spreadsheet cut it. Sites past the cut are unknown, and the loader should warn that coverage is truncated.
- Detection uses patterns only: the quoted comment lines count as data lines to `data_shape`.

**Example**

```text
# This data file generated by 23andMe at: Wed Jan  6 13:54:47 2016,,,
#,,,
# We are using reference human assembly build 37 (also known as Annotation Release 104).,,,
#,,,
# rsid,chromosome,position,genotype
rs3131972,1,752721,AG
rs11240777,1,798959,GG
rs4475691,1,846808,--
```

### `living_dna`

Living DNA raw genotype download (Sirius chip), 23andMe column layout.

- **Family:** `genotype_table`
- **Vendor:** Living DNA
- **Registry position:** 29 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** forward-strand allele letters; an indel is its two allele sequences written back to back (`TGTG`, `CCA`), so a split is ambiguous without the site's alleles
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `^# Human Genome Reference Build 37` names GRCh37

**Recognised by**

- every one of these is found in the head:
  - `^# Living DNA customer genotype data download file`
- one of the first 5 data lines matches in full:
  - `(rs\d+|AX-\d+|i\d+)\t(\d{1,2}|X|Y|MT)\t\d+\t[ACGT]+`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid, or an Affymetrix probe id such as AX-32363003 |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | forward-strand alleles; indels as both allele sequences concatenated |

**From real files**

- Differs from 23andMe in the preamble, the AX- probe ids and the multi-letter indel genotypes; the columns and the commented column line are the same.
- A file named `-autosomal` also holds chromosome X.
- No no-call spelling was seen in the one real file; no-calls may be dropped.

**Example**

```text
# Living DNA customer genotype data download file version: 1.0.2
# Genotype chip: Sirius
# Human Genome Reference Build 37 (GRCh37.p13).
# Genotypes are presented on the forward strand.
#
# rsid	chromosome	position	genotype
rs3131972	1	752721	AG
rs114525117	1	759036	GG
AX-32363003	1	806111	AG
rs4970383	1	838555	CC
```

### `selfdecode`

SelfDecode genotype export (`# Pipeline Version`), 23andMe column layout.

- **Family:** `genotype_table`
- **Vendor:** SelfDecode
- **Registry position:** 30 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two allele letters with heterozygotes in alphabetical order (AG, never GA), so the order carries neither phase nor strand; one letter on male X, Y and MT
- **No-call spellings:** none declared
- **Build:** named only by a header statement matching one of:
  - `^# Assembly: GRCh38` names GRCh38

**Recognised by**

- every one of these is found in the head:
  - `^# Pipeline Version: `
  - `^# Source: SelfDecode`
- one of the first 5 data lines matches in full:
  - `[^\t]+\t(\d{1,2}|X|Y|MT)\t\d+\t([ACGTDI]{1,2}|--)`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid; several ids of one site are joined with `;`, and vendor ids (`seq-…`, `indel…`, `16:28617883_CNV_SULT1A`) carry `:`, `_`, `(` and `>` |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, one on haploid sites |

**From real files**

- The positions really are GRCh38 when the header says so (rs9651229 at chr1:632287).
- `# Sex:` in the header states the sample's sex, which explains the haploid X and Y.
- No no-call spelling was seen in the one real file; no-calls may be dropped.

**Example**

```text
# Pipeline Version: 2.00
# Source: SelfDecode
# Assembly: GRCh38
# Chromosomes: 1-22, X, Y, MT
# rsid	chromosome	position	genotype
rs4477212	1	82154	AG
rs8179455	1	264681	GG
rs9701055	1	630053	CT
rs9651229	1	632287	CC
```

### `plink_23`

PLINK `--recode 23` output: the 23andMe layout written by PLINK from another source.

- **Family:** `genotype_table`
- **Vendor:** PLINK
- **Registry position:** 31 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two allele letters; one on haploid sites; a missing allele is `-` (`T-`)
- **No-call spellings:** `--`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^# This data file generated by PLINK (at|on):`
- none of these is found:
  - `Genes for Good`
- one of the first 5 data lines matches in full:
  - `\S+\t(\d{1,2}|X|Y|MT)\t\d+\t[ACGTDI-]{1,2}`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | the source's variant id: an rsid, or whatever the source used |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, one on haploid sites, `-` for a missing allele |

**From real files**

- The header says to consult the original source for the build, so the build is unknown.
- Which source PLINK converted is not stated; the one real file covers far more sites than an array (positions from 13,896 on chromosome 1), so it may hold imputed or sequenced calls.

**Example**

```text
# This data file generated by PLINK at: Mon Mar 16 03:53:37 2015
#
# Below is a text version of your data.  Fields are TAB-separated.
# For further information (e.g. which reference build was used), consult the
# original source of your data.
#
# rsid	chromosome	position	genotype
rs3094315	1	752566	AG
rs3131972	1	752721	GG
rs11240777	1	798959	--
rs4475691	1	846808	CT
```

### `promethease_tab_genotype`

Promethease `tab-genotype` report: SNPedia genotypes and genosets with magnitude and summary.

- **Family:** `annotated_variants`
- **Vendor:** Promethease
- **Registry position:** 32 of 39, a signature (`just_vcf.formats.rsid_tables`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment
- **Coordinates:** none (the rows are not placed on a genome)
- **Genotype:** inside the name, SNPedia style: `(A;G)`, `(I;I)`, `(D;I)`; oriented as SNPedia orients the site, which may be the minus strand
- **No-call spellings:** `(-;-)`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A#Promethease\ttab-genotype\t`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `name` | string | `Name` | `rs…(A;B)` or `i…(A;B)`: a site with its genotype; `gs…`: a genoset |
| `magnitude` | float | `Magnitude` | SNPedia magnitude; the string `None` when absent |
| `repute` | string | `Repute` | `Good`, `Bad` (case varies) or empty; absent from older reports such as 0.1.110 |
| `population_freq` | string | `PopulationFreq` | a percentage such as `20.4%`, `None`, or empty |
| `description` | string | `Description` | SNPedia summary |

**From real files**

- Not a genotype table: a report on the sites SNPedia knows, with no chromosome or position.
- `#Column-names` sits on its own comment line; the names are on the next, uncommented line.
- `gs…` rows are genosets, conclusions over several sites, not sites.
- `(-;-)` may be a no-call or a real deletion; the report does not say which.
- `None` is written as a string in Magnitude and PopulationFreq.
- The column set changes with the Promethease version: 0.1.110 has no Repute column, 0.1.161 has it. Take the names from the column line, never from position.
- A Description can hold a tab, giving the row an extra field.
- Often found as `report_Tab.txt` inside a zip beside HTML reports and JavaScript.

**Example**

```text
#Promethease	tab-genotype	0.1.161	2014-02-08 23:50
#Column-names	
Name	Magnitude	Repute	PopulationFreq	Description
rs1333049(C;G)	None		None	
rs10490924(G;G)	0.0		None	normal
i3001487(-;-)	None			
```

### `freebayes_varscan2_keyvalue`

PGP-HMS's headerless per-site export of a freeBayes+VarScan2 merge: CHROM, POS, REF, ALT, then both callers' genotype and depth as positional key=value fields, then annotation.

- **Family:** `variant_calls`
- **Vendor:** PGP-HMS
- **Registry position:** 33 of 39, a fallback, tried after every signature (`just_vcf.formats.vcf`)
- **Loads:** yes, through `just_vcf.loaders.keyvalue`
- **Delimiter:** tab
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** VCF-style allele indices into REF and ALT (`0/1`, `1/2`), one per caller, each after its `KEY=` prefix
- **No-call spellings:** `NULL`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- one of the first 5 data lines matches in full:
  - `chr(\d{1,2}|X|Y|M)\t\d+\t[ACGTN]+\t[ACGTN]+(,[ACGTN]+)*\tFBGenoType=[^\t]*\tFBRefAlleleCount=[^\t]*\tFBReferenceAlleleQ=[^\t]*\tFBTotalDepth=[^\t]*\tVSGenoType=[^\t]*\tVSRefAlleleCount=[^\t]*\tVSReferenceAlleleQ=[^\t]*\tVSTotalDepth=[^\t]*(\t.*)?`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `chromosome` | string |  | chr1-chr22, chrX, chrY; `chr` prefix kept |
| `position` | integer |  | 1-based position, the VCF POS of the sibling merge |
| `ref` | string |  | reference allele |
| `alt` | string |  | alternate alleles, comma-separated |
| `fb_genotype` | string |  | `FBGenoType=`: freeBayes GT, NULL when not called |
| `fb_ref_allele_count` | integer |  | `FBRefAlleleCount=`: `.` or NULL when absent |
| `fb_ref_allele_quality` | integer |  | `FBReferenceAlleleQ=`: NULL when absent |
| `fb_total_depth` | integer |  | `FBTotalDepth=` |
| `vs_genotype` | string |  | `VSGenoType=`: VarScan2 GT, NULL when not called |
| `vs_ref_allele_count` | integer |  | `VSRefAlleleCount=`: NULL when absent |
| `vs_ref_allele_quality` | integer |  | `VSReferenceAlleleQ=`: NULL when absent |
| `vs_total_depth` | integer |  | `VSTotalDepth=` |

**From real files**

- It is not a VCF body: no ID, QUAL, FILTER, INFO or FORMAT. Columns 5-12 are `KEY=value` in a fixed order; strip the key before casting. Missing values are `NULL`, and the freeBayes ref count is also `.`.
- Every row carries a variable tail after column 12 (23 to 84 fields in total): EFF_*/VEP_* annotation with bare tokens (EFF_RESULTS, MODIFIER, NOVEL_VAR, CREATE_NEW_VAR) repeated per transcript, and population frequencies (KG, *_MAF, EXAC_FRQ, EVS_FRQ). Drop it or keep it as text; `columns` covers the fixed twelve only.
- Rows are unsorted. Coordinates and alleles match the sibling VCF of the same arvados collection exactly (20,000 of 20,000 sampled rows at offset 0), so POS is 1-based on hg19. The file states no build; the sibling's `##reference=/vbin/ref/hg19.fa` does.
- A NULL genotype means that caller did not report the site, and all four of its fields are NULL. freeBayes also writes `0/0` on some rows (a few percent), mostly where VarScan2 is NULL.

**Example**

```text
chr1	752566	G	A	FBGenoType=0/1	FBRefAlleleCount=12	FBReferenceAlleleQ=475	FBTotalDepth=33	VSGenoType=0/1	VSRefAlleleCount=12	VSReferenceAlleleQ=39	VSTotalDepth=31
chr1	752721	A	G	FBGenoType=1/1	FBRefAlleleCount=.	FBReferenceAlleleQ=NULL	FBTotalDepth=27	VSGenoType=1/1	VSRefAlleleCount=NULL	VSReferenceAlleleQ=NULL	VSTotalDepth=27
chr1	846808	C	T	FBGenoType=0/1	FBRefAlleleCount=23	FBReferenceAlleleQ=869	FBTotalDepth=46	VSGenoType=NULL	VSRefAlleleCount=NULL	VSReferenceAlleleQ=NULL	VSTotalDepth=NULL
```

### `oxford_gen`

Oxford `.gen` genotype probabilities, one site per line, space-separated.

- **Family:** `genotype_table`
- **Vendor:** generic (IMPUTE2 / Oxford .gen; Gencove)
- **Registry position:** 34 of 39, a fallback, tried after every signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** single space
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** three probabilities per sample for AA, AB, BB over the listed alleles A and B; graded evidence, not a hard call. A and B are not REF and ALT: REF comes from the reference
- **No-call spellings:** `0 0 0`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- one of the first 5 data lines matches in full:
  - `\S+ \S+ \d+ ([ACGTN]+|-) ([ACGTN]+|-)( (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+) (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+) (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+))+`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `snp_id` | string |  | `---` on imputed sites; on typed sites the input's first column (the chromosome in Gencove exports, `23` for X) |
| `rsid` | string |  | rsid, or `chr1:13957:D` / `MERGED_DEL_…` for others |
| `position` | integer |  | 1-based |
| `allele_a` | string |  | first allele; `-` for a deleted allele |
| `allele_b` | string |  | second allele; `-` for a deleted allele |
| `p_aa` | float |  | probability of genotype AA |
| `p_ab` | float |  | probability of genotype AB |
| `p_bb` | float |  | probability of genotype BB |

**From real files**

- The chromosome is not in the rows: Gencove ships one member per chromosome (`id_…_chr1.gen`) and only typed rows repeat it in column 1. The loader must take it from the member name or refuse; sniff never trusts names.
- Column 1 `---` marks an imputed site and anything else a typed one; keep that distinction.
- Several samples make 3 probability columns each; `columns` declares one sample.
- No build is stated; the positions seen match GRCh37, but that is not a header statement.
- `0 0 0` is the Oxford no-call spelling; probabilities are rounded and need not sum to 1.

**Example**

```text
--- rs58108140 10583 G A 0.638 0.359 0.004
--- rs189107123 10611 C G 0.981 0.019 0
1 rs3131972 752721 A G 0 0 1
1 rs12124819 776546 A G 0 1 0
```

### `soapsnp_cns`

SOAPsnp consensus output, SNP rows only, 17 SOAPsnp columns plus distance to the nearest SNP.

- **Family:** `variant_calls`
- **Vendor:** BGI (SOAPsnp)
- **Registry position:** 35 of 39, a fallback, tried after every signature (`just_vcf.formats.other_tables`)
- **Loads:** no
- **Delimiter:** tab
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** IUPAC code in `consensus`: a plain letter is homozygous for that base, R/Y/K/M/S/W heterozygous for the pair it stands for
- **No-call spellings:** none declared
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- one of the first 5 data lines matches in full:
  - `chr([0-9]{1,2}|X|Y|M)\t\d+\t[ACGTN]\t[ACGTRYKMSWBDHVN]\t\d+\t[ACGTN]\t\d+\t\d+\t\d+\t[ACGTN]\t\d+\t\d+\t\d+\t\d+\t\d+(\.\d+)?\t\d+(\.\d+)?\t[01]\t\d+`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `chromosome` | string |  | chr-prefixed |
| `position` | integer |  | 1-based |
| `ref` | string |  | reference base as stated |
| `consensus` | string |  | IUPAC genotype code |
| `consensus_quality` | integer |  | Phred quality of the consensus |
| `best_base` | string |  | best-supported base |
| `best_base_quality` | integer |  | mean quality of best base |
| `best_base_unique_reads` | integer |  | uniquely mapped reads, best base |
| `best_base_all_reads` | integer |  | all mapped reads, best base |
| `second_base` | string |  | second-best base |
| `second_base_quality` | integer |  | mean quality of second base |
| `second_base_unique_reads` | integer |  | uniquely mapped reads, second base |
| `second_base_all_reads` | integer |  | all mapped reads, second base |
| `depth` | integer |  | sequencing depth |
| `rank_sum_p` | float |  | rank-sum test p-value |
| `copy_number` | float |  | mean copy number of the region; 255 = none |
| `dbsnp` | integer |  | 1 when the site is a known dbSNP site |
| `nearest_snp_distance` | integer |  | distance to the nearest listed site |

**From real files**

- The producer is inferred from the layout: the file names no tool. Column 18 equals the distance to the neighbouring row on the same chromosome.
- Rows with depth 0 and consensus quality 1 are prior-driven calls at dbSNP sites (`dbsnp` 1) with no read evidence; the loader must carry depth and quality through, never treat them as observed.
- Only SNP rows are present; an absent site is unknown.
- No build is stated. `chrM` must be mapped to `MT` explicitly.

**Example**

```text
chr1	742429	G	R	34	G	36	12	12	A	30	10	10	22	0.412	1.00000	1	15882
chr1	758311	G	A	87	A	38	19	19	G	0	0	0	19	1.00000	1.00000	1	15882
chr1	995669	C	Y	41	C	35	9	9	T	33	8	8	17	0.731	1.00000	1	237358
```

### `rsid_result_csv`

Any other `RSID,CHROMOSOME,POSITION,RESULT` table: the layout FTDNA, MyHeritage and Vitagene share, when the content does not show which vendor wrote it.

- **Family:** `genotype_table`
- **Vendor:** unknown
- **Registry position:** 36 of 39, a fallback, tried after every signature (`just_vcf.formats.vendor_csv`)
- **Loads:** yes, through `just_vcf.loaders.table`
- **Delimiter:** comma
- **Comment prefix:** `#`
- **Column names:** the first line that is not a comment
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** two allele letters in one field, quoted or not; D/I codes (DD, DI, II) on indel probes in GSA-chip exports
- **No-call spellings:** `--`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^RSID,CHROMOSOME,POSITION,RESULT\r?$`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string | `RSID` | rsid or a vendor probe id |
| `chromosome` | string | `CHROMOSOME` | 1-22, X, Y, MT; no chr prefix |
| `position` | integer | `POSITION` | 1-based position on an unstated build |
| `genotype` | string | `RESULT` | two allele letters, or one |

**From real files**

- The vendor is not claimed. A Vitagene GSA export seen with this layout is unquoted, LF, starts with unplaced probes on chromosome `0` at position 0, uses `XY` for the pseudoautosomal region, and has probe ids such as `GSA-rs…`, `exm…`, `kgp…` and gene-named assay ids.
- Rows on chromosome `0` have no position; keep them out of any positional join.
- Lines starting with `#` above the header are comments. No build is stated unless a comment says so, and no hint is declared because no such file was available to check the wording.

**Example**

```text
RSID,CHROMOSOME,POSITION,RESULT
GSA-rs301806,1,8482078,CT
exm2268861,1,6550376,TT
GSA-rs36083022,1,10689678,DI
kgp22779185,0,0,--
```

### `genographic_geno2_headerless`

Genographic Geno 2.0 results with the `SNP,Chr,Allele1,Allele2` line cut off.

- **Family:** `genotype_table`
- **Vendor:** National Geographic Genographic Project
- **Registry position:** 37 of 39, a fallback, tried after every signature (`just_vcf.formats.vendor_csv`)
- **Loads:** no
- **Delimiter:** comma
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** none (the rows are not placed on a genome)
- **Genotype:** two columns, one letter each; I/D on indel probes; Y and Mt calls written as the same letter twice; strand not stated
- **No-call spellings:** `-`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A\n*[^,\n]+,([1-9]|1\d|2[0-2]|X|Y|Mt),[ACGTDI-],[ACGTDI-]$`
- one of the first 5 data lines matches in full:
  - `[^,]+,([1-9]|1\d|2[0-2]|X|Y|Mt),[ACGTDI-],[ACGTDI-]`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `snp` | string | `SNP` | marker name: rs/kgp ids, Y-SNP names, mt numbers |
| `chromosome` | string | `Chr` | 1-22, X, Y, Mt |
| `allele1` | string | `Allele1` | first allele |
| `allele2` | string | `Allele2` | second allele |

**From real files**

- There is no position column: sites must be placed through the chip manifest by marker name, and the build is whatever that manifest says. The file states none.
- Y rows are named by Y-SNP (CTS…, M…), autosomal rows mostly by kgp ids, and Mt rows by bare numbers that look like mtDNA positions; the same Mt name can repeat on consecutive rows.
- The header line ends LF while the data rows end CRLF.
- Chromosome `Mt`, not `MT`.
- Recognised by shape: the first non-blank line must already be a marker row (one export opens with an empty line), and the `Mt` spelling helps keep it apart from other headerless CSVs.

**Example**

```text
kgp10004422,12,A,G
CTS100,Y,T,T
CTS10004,Y,-,-
195,Mt,C,C
```

### `rsid_table`

The 23andMe layout with only its column line: the vendor preamble was removed or never written.

- **Family:** `genotype_table`
- **Vendor:** unknown
- **Registry position:** 38 of 39, a fallback, tried after every signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** `#`
- **Column names:** the last comment line before the data (`# rsid\tchromosome…`)
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** as 23andMe: allele letters, one on haploid sites, D/I on indel probes
- **No-call spellings:** `--`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `^# rsid\tchromosome\tposition\tgenotype\r?$`
- one of the first 5 data lines matches in full:
  - `(rs|i)\d+\t(\d{1,2}|X|Y|MT)\t\d+\t([ACGTDI]{1,2}|--)`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid, or a vendor-internal id such as i713426 |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, or one on haploid sites |

**From real files**

- The build is unknown: real files of this shape are 23andMe exports of both build 36 and build 37 with the preamble cut off.
- Also written by converters: one real file opens `# This file emulates the format of 23andMe` and carries two column lines, the first with an empty genotype name; the last one wins.
- A blank line may precede the column line.

**Example**

```text

# rsid	chromosome	position	genotype
rs4477212	1	82154	AG
rs3131972	1	752721	GG
rs11240777	1	798959	--
rs4970383	1	838555	CC
```

### `rsid_table_headerless`

The 23andMe layout with no comment or column line at all.

- **Family:** `genotype_table`
- **Vendor:** unknown
- **Registry position:** 39 of 39, a fallback, tried after every signature (`just_vcf.formats.rsid_tables`)
- **Loads:** yes, through `just_vcf.loaders.rsid_table`
- **Delimiter:** tab
- **Comment prefix:** none
- **Column names:** the file names no columns; they are the ones listed below, in file order
- **Coordinates:** 1-based (a position names the base itself, as VCF `POS` does)
- **Genotype:** as 23andMe: allele letters, one on haploid sites, D/I on indel probes
- **No-call spellings:** `--`
- **Build:** no header statement is read; the build is unknown until inferred

**Recognised by**

- every one of these is found in the head:
  - `\A(rs|i)\d+\t(\d{1,2}|X|Y|MT)\t\d+\t`
- one of the first 5 data lines matches in full:
  - `(rs|i)\d+\t(\d{1,2}|X|Y|MT)\t\d+\t([ACGTDI]{1,2}|--)`

**Columns**

| Name | Type | Other spellings | Meaning |
|---|---|---|---|
| `rsid` | string |  | rsid, or a vendor-internal id such as i713426 |
| `chromosome` | string |  | 1-22, X, Y, MT; no chr prefix |
| `position` | integer |  | 1-based position on the stated build |
| `genotype` | string |  | two allele letters, or one on haploid sites |

**From real files**

- The build is unknown; one real file of this shape is on build 36 (rs3094315 at 742429).
- Recognised only when line 1 is already a data line, so a TSV with any preamble is not captured.
- One real file is two exports concatenated: a stray first row with an empty fourth field, then 23andMe-shaped 4-column rows, then AncestryDNA-shaped 5-column rows (allele1, allele2; chromosomes 23-26; `0` for no call) for the whole genome again. The loader checks the width of every row and reads that second export as a facet of its own (`relation = appended_export`, F3).

**Example**

```text
rs3094315	1	742429	AG
rs12562034	1	758311	GG
rs3934834	1	995669	--
rs9442372	1	1008567	AA
```

## Broken and odd shapes met in real files

Shapes found in real uploads that a reader has to survive: the dogfooding findings F1-F17 (docs/DOGFOODING.md, docs/PREVIOUS_ISSUES.md) and the dialects behind the contract gaps in docs/FORMATS.md. Each states the shape, gives a synthetic example of it, says what `sniff` answers (the test suite runs every such example) and what just-vcf does with it. *handled* means read; *refused* means refused on purpose with a named error; *open* names a gap.

### A zip member compressed with Deflate64 (F1)

<a id="shape-deflate64-zip"></a>

A zip archive whose member is stored with compression method 9 (Deflate64), which Windows Explorer writes for members over 2 GiB. The archive's directory reads; the member does not inflate with a plain Deflate decoder. PPMd (method 98) and AES encryption (99) fail the same way.

```text
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	50	PASS	.	GT	0/1
```

*The member's text; its local and central zip headers both say compression method 9.*

**`sniff` says:** `unreadable`, reason naming Deflate64 (method 9)

**What just-vcf does** (open): `sniff` answers `unreadable` with the member's name and the reason `zip compression Deflate64 (method 9) is not supported`; `load` refuses anything `sniff` did not detect. Reading it needs a Deflate64 inflater, which is a dependency decision not yet taken.

Code and tests: sniff.py `_ZIP_METHODS`, `test_unsupported_zip_method_is_unreadable_and_named`.

### A 23andMe row with a doubled tab after the id (F2)

<a id="shape-doubled-tab"></a>

Some 2012-2014 23andMe exports write a few dozen rows with an empty field between the id and the chromosome: five fields where the layout has four.

```text
# This data file generated by 23andMe at: Tue Jan  1 00:00:00 2013
#
# We are using reference human assembly build 37 (also known as Annotation Release 104).
#
# rsid	chromosome	position	genotype
rs4477212	1	82154	AA
rs3094315		1	752566	AG
rs3131972	1	752721	GG
rs12562034	1	768448	--
```

**`sniff` says:** `detected`, as `23andme`, build GRCh37

**What just-vcf does** (handled): Detected as `23andme`: one well-formed row among the first data lines is enough. The loader reads a row of exactly five fields whose second field is empty as if that field were not there, and counts those rows in a warning of their own. An empty field anywhere else still drops the row, counted.

Code and tests: loaders/rsid_table.py `read_calls`, `test_a_doubled_tab_after_the_id_is_read_past_and_counted`.

### Two exports in one headerless file (F3)

<a id="shape-appended-export"></a>

A headerless file that holds a 23andMe-shaped export (id, chromosome, position, genotype) and then a second export of the same person in the two-allele layout (id, chromosome, position, allele1, allele2; chromosomes 23-26 for X, Y, XY and MT; `0` for no call), with no line between them. The file can also open with one stray row that has an empty field before the genotype.

```text
rs548049170	1	69869		TT
rs4477212	1	82154	AG
rs3094315	1	752566	AA
rs3131972	1	752721	AG
rs4477212	1	82154	A	G
rs3094315	1	752566	A	A
rs3131972	1	752721	0	0
```

*Shortened: in the file it describes each half runs to hundreds of thousands of rows, and the loader reads the second half only when it is at least 1,000 rows long.*

**`sniff` says:** `detected`, as `rsid_table_headerless`

**What just-vcf does** (handled): Detected as `rsid_table_headerless`. A run of at least 1,000 rows of one other width that ends the file, after other rows, and whose first five rows are in the two-allele layout, is read by that layout's rules as a second facet in the same `Loaded`: `relation = appended_export`, `part = 1` in its facet id, the file's `format_id` kept. Its build is not the one the first export states: it stays null until `resolve_ref` infers it from the rows. A shorter tail, or one in no known layout, is warned about and not read. The stray first row is dropped and counted.

Code and tests: loaders/rsid_table.py `Appended`, `_appended_start`, `_read_appended`, `test_an_appended_export_in_another_layout_is_its_own_facet`, `test_a_tail_that_is_short_or_in_no_known_layout_is_not_read`, `test_the_build_a_preamble_states_is_the_first_exports_only`.

### Indel genotypes written as the two alleles run together (F4, F17)

<a id="shape-run-together-indels"></a>

Living DNA writes an indel call as its two allele sequences back to back, with nothing between them (`CCA` for `C` and `CA`, `GTTGTT` for `GTT` twice). Where one allele ends and the other starts is not written.

```text
# Living DNA customer genotype data download file version: 1.0.2
# Genotype chip: Sirius
# Human Genome Reference Build 37 (GRCh37.p13).
# Genotypes are presented on the forward strand.
#
# rsid	chromosome	position	genotype
rs3131972	1	752721	AG
rs11240777	1	798959	CCA
rs4970383	1	838555	GTTGTT
rs4475691	1	846808	CC
```

**`sniff` says:** `detected`, as `living_dna`, build GRCh37

**What just-vcf does** (open): The loader keeps the text in `source_gt` and leaves `gt_alleles` null. `resolve_ref` then splits a text of three or more bases where exactly one split fits: both pieces start with the same base and one of them equals the reference over its own span from `pos`. REF and ALT follow only where the split says which allele is the reference; a homozygous call (`CA/CA` cannot say whether the site is `CA>C` or `C>CA`) and a het whose two alleles both match the reference keep REF unread. Texts that fit several splits (repeats) or none stay as loaded, counted. A facet whose strand check says `mixed` splits nothing. Settling the rest needs the chip's own site alleles from its manifest, which the library does not read (F17).

Code and tests: reference.py `_split_concatenated`, `CONCATENATED_ALLELES`, `test_living_dna_indel_genotypes_split_against_the_reference`, `test_only_living_dna_texts_are_split`, `test_a_mixed_strand_living_dna_file_splits_nothing`.

### A VCF inside another compression or archive layer (F6)

<a id="shape-archived-vcf"></a>

A VCF compressed with bz2 or xz, or inside a zip or tar, or BGZF under a plain `.vcf` name. The name says nothing reliable about the bytes.

```text
##fileformat=VCFv4.2
##reference=file:///ref/hs37d5.fa
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	50	PASS	.	GT	0/1
```

*The payload; the test compresses it with bz2.*

**`sniff` says:** `detected`, as `vcf`, inside `bz2`, build GRCh37

**What just-vcf does** (handled): `sniff` peels the layers by magic bytes and reports them in `layers` (here `bz2`). polars-bio reads plain, gzip and BGZF VCF directly; any other layer is decompressed into the loader's cache as a BGZF copy (a plain-text copy of a 1.27 GB bz2 gVCF was 14 GB) and read from there.

Code and tests: sniff.py `_sniff_layer`, loaders/vcf.py `readable_path`, `test_compression_is_detected_by_magic_not_extension`, `test_archived_vcf_is_rewritten_and_loaded`.

### A VCF re-saved by a spreadsheet (F7)

<a id="shape-spreadsheet-vcf"></a>

A VCF opened in a spreadsheet and saved again as tab-separated text: every meta line is padded with tabs to the row width, any field holding a comma or a quote is wrapped in quotes with its own quotes doubled (so a meta line with a `Description` starts `"##`), and lines end in a bare CR.

```text
##fileformat=VCFv4.2									
"##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"									
"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=""Allelic depths"">"									
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	50.1	PASS	.	GT:AD	"0/1:12,9"
1	752721	rs3131972	A	G	.	PASS	.	GT	./.
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): Detected as `vcf`: `sniff` reads CR line ends as line ends. The loader looks at the first 64 KiB of the payload: a line starting `"##`, or CRs and no LF, marks a re-save. It then reads the file by universal newlines, cuts each line at tabs, unwraps quoted fields and undoes the doubled quotes, drops the padding from `##` lines, and polars-bio reads the rewritten copy. A warning says so, and that a number whose text the spreadsheet rewrote (`6980.90` as `6980.9`) cannot be told. A separate format id was rejected: the delimiter is still a tab, and the quoting is a container's mangling, like CRLF.

Code and tests: loaders/vcf.py `is_resaved`, `undo_resave`, `_unquote`, `test_a_vcf_resaved_by_a_spreadsheet_loads_as_the_text_it_was_saved_from`, `test_a_resaved_vcf_with_a_field_past_the_csv_limit_loads`, `test_crlf_and_plain_vcfs_are_not_taken_for_a_resave`.

### A VCF column line written with two `#` (F7, F15)

<a id="shape-two-hash-column-line"></a>

Some SV and CNV VCFs write the column line as `##CHROM` and put a third `#` on most meta lines (`###contig`, `###FORMAT`, `###reference`), while `##fileformat` and `##ALT` are written right. Their CNV rows write ALT as the bare letters `CNV`, without angle brackets.

```text
##fileformat=VCFv4.2
##ALT=<ID=DEL,Description="Deletion">
###contig=<ID=1,length=249250621>
###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	.	G	<DEL>	.	PASS	SVTYPE=DEL;END=753566	GT	0/1
1	846808	.	C	CNV	.	PASS	END=860000	GT	0/1
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): Detected as `vcf`. Where no `#CHROM` line exists, the loader takes a meta line that is VCF's column line with one `#` too many as the column line and reads the file through a rewritten copy. A `###` line stays what VCF makes of it, a meta line whose key starts with `#`: no contig or FORMAT is declared from it, so GT and END are declared by the loader as for any undeclared key (the build hints still read its text). An ALT that is not `ACGTN` bases and not a symbolic `<…>` allele, such as `CNV`, is kept verbatim, gets no variant class, is not counted as a small-variant site, and is counted in one warning (F15).

Code and tests: loaders/vcf.py `read_header`, `HASHED_COLUMN_LINE`, callset.py `variant_class_expr`, `test_a_column_line_written_with_two_hashes_is_read`, `test_an_alt_that_is_no_vcf_allele_has_no_variant_class`.

### A sites-only VCF (F8)

<a id="shape-sites-only-vcf"></a>

A VCF with the eight fixed columns and no FORMAT or sample column, as annotation sources write.

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	752566	rs3094315	G	A	.	.	.
1	752721	rs3131972	A	G	.	.	.
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (refused): Detected as `vcf`. `load` refuses it with `NoSampleError` before reading anything: just-vcf holds a person's calls, and a frame of null genotypes would look like a person with nothing called. `vcf_samples` returns `[]`.

Code and tests: loading.py `_load_vcf`, `test_sites_only_vcf_is_refused`.

### An array export that states no build (F9)

<a id="shape-no-build-stated"></a>

FTDNA exports, preamble-less rsid tables, PLINK output and Final Reports name no build, and an array table has no `##contig` lengths to infer one from. FTDNA shipped both build 36 and build 37 exports in one layout.

```text
RSID,CHROMOSOME,POSITION,RESULT
"rs4477212","1","72017","AA"
"rs3094315","1","742429","AG"
"rs3131972","1","742584","GG"
```

*Positions are build 36 ones; nothing in the file says so.*

**`sniff` says:** `detected`, as `ftdna`, build unknown

**What just-vcf does** (handled): `sniff` reports `build = None`: no vendor default is ever assumed. `resolve_ref` infers the build of such a facet from its own rows: for each reference collection it considers, the share of sampled SNV rows whose called alleles include the reference base there. The winner must clear a floor and a margin over the others; it is recorded as `assembly_label` with basis `inferred`, the scores in one warning line. Rows that do not settle it refuse the facet with the scores. A collection whose contig 1 is 247,249,719 bases is GRCh36.

Code and tests: assembly.py `infer_assembly_from_rows`, `MIN_SHARE`, `MIN_MARGIN`, reference.py `_infer_from_rows`, `test_an_array_without_a_build_gets_it_from_its_rows`.

### A GP whose probabilities do not sum to 1 (F10)

<a id="shape-gp-not-summing"></a>

FORMAT GP values all within [0, 1] that do not sum to 1 on some rows (phred GPs rounded to integers and converted back give such rows). From the numbers alone the scale is ambiguous.

```text
##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	.	PASS	.	GT:GP	0/1:0.03,0.9,0.05
1	752721	rs3131972	A	G	.	PASS	.	GT:GP	1/1:0,0.01,0.99
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): Values above 1 are phred and values summing to 1 are probabilities, whatever the header says. Values the numbers cannot settle are read as probabilities under a VCF 4.3 or later header, which defines GP so, and the rows that do not sum to 1 are counted in a warning. Under 4.2 or earlier, or with no version, GP is withheld. The FORMAT line's Description is free text and is not read.

Code and tests: loaders/vcf.py `_gp`, `test_gp_scale_is_read_from_decisive_values`, `test_gp_from_phred_survives_its_own_export`.

### A deletion written as ALT `-` (F12)

<a id="shape-dash-alt"></a>

Some VCF uploads write a deletion as REF the deleted bases and ALT `-`, the dash notation VCF does not have. The anchor base is missing; REF is the deleted span.

```text
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752721	.	AG	-	.	PASS	.	GT	0/0
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): The loader reads `-` as the empty allele `""`, the one spelling of an empty allele (CGI's native deletions use it too), and counts the records in one warning naming them unanchored deletions. `normalize`, once `resolve_ref` has matched REF, anchors it on the base before (`CT>-` at 30 becomes `GCT>G` at 29) and left-aligns it; at a contig's first base it stays as written, counted. `write_vcf` refuses any ALT that is no VCF allele as `alt_not_an_allele`.

Code and tests: loaders/vcf.py `DASH_DELETION`, normalize.py `_align`, writing.py `_unwritable`, `test_a_dash_alt_is_read_as_the_empty_allele_and_counted`, `test_a_dash_deletion_is_anchored_on_the_base_before`, `test_an_alt_that_is_no_vcf_allele_is_never_written`.

### An array export with part of its probes on the minus strand (F13)

<a id="shape-minus-strand-probes"></a>

One FTDNA export writes about one probe in eight on the minus strand, with nothing in the row to say so. It shows only in aggregate: heterozygous rows carrying the complement of the reference base and not the base, which no forward-strand probe can write.

```text
RSID,CHROMOSOME,POSITION,RESULT
"rs4477212","1","82154","TT"
"rs3094315","1","752566","CT"
"rs3131972","1","752721","AG"
```

*On GRCh37 the reference base at 1:752566 is `G`; the het `CT` there carries its complement `C` and not `G`. The `AG` row carries its site's reference base `A`, so it reads the same on either strand; the homozygous `TT` row cannot show its strand at all.*

**`sniff` says:** `detected`, as `ftdna`

**What just-vcf does** (handled): `resolve_ref` checks the strand of every array facet on a sample of its heterozygous two-allele rows. When more than 2% of at least 200 such hets carry the complement and not the base, the facet gets `strand_check = mixed`: every row that is not strand-safe (homozygous calls, half-calls, hets without the base) keeps `ref` and gets `alts` and `gt` null, each reason counted, and the safe hets keep their alleles. Complementing was rejected: an A/T or C/G probe is its own complement, the homozygous rows of a flipped probe cannot be found, and no public strand table names the flipped probes.

Code and tests: strand.py `check_strand`, `withheld_reason`, reference.py `_check_strand`, `test_a_mixed_strand_array_keeps_ref_and_withholds_the_alleles_it_cannot_read`.

### A liftover chain that lists every block twice (F14)

<a id="shape-chain-blocks-twice"></a>

Ensembl's `GRCh38_to_GRCh37` chain repeats every block of its main chains under a second chain id, with the same source span, target and strand. Counted per block, every lifted base looks as if it had two targets. Not an input genotype format, but a file the library reads.

```text
chain 1000 1 248956422 + 10000 10500 1 249250621 + 10000 10500 1
500

chain 1000 1 248956422 + 10000 10500 1 249250621 + 10000 10500 2
500
```

*UCSC chain text: a `chain` header line, then block sizes; the two chains differ only in id.*

**What just-vcf does** (handled): `liftover` cuts the source axis where the set of targets changes, a target being a contig, strand and offset, so two blocks sending every shared base to the same place count once. A position under blocks with different targets is still `multiple_targets`.

Code and tests: liftover.py `_segments`, `test_a_chain_listing_its_blocks_twice_lifts_as_once`, `test_segments_count_targets_not_blocks`.

### A BGZF file with a damaged block (F16)

<a id="shape-damaged-bgzf"></a>

A BGZF file one of whose blocks does not inflate, or fails its CRC32 or its ISIZE, deep inside the file. The head reads, so the file is recognised; the damage shows only when the records are read.

*A byte-level fault, so there is no text example.*

**What just-vcf does** (refused): When a read fails and the input's outer layer is BGZF, `load` walks the file block by block (header, `BC` size, inflate, CRC32, ISIZE) and raises `DamagedBgzfError` with the facet's id, the first bad block's compressed offset and what is wrong with it; it never names the path. A failure whose blocks all read is re-raised unchanged. Damage under another layer (a BGZF inside a zip, a plain gzip) keeps the reader's own error, since only an outer BGZF has block offsets to name.

Code and tests: bgzf.py `find_damage`, loading.py `DamagedBgzfError`, `test_load_names_the_facet_and_the_block_never_the_path`, `test_damage_in_plain_gzip_keeps_the_readers_error`.

### Complete Genomics columns that change with the format version

<a id="shape-cgi-format-versions"></a>

A CGI `var` file's columns depend on its `#FORMAT_VERSION`. 1.6 and 1.7 have one `totalScore` where 2.x has `varScoreVAF` and `varScoreEAF`; 2.0 and 2.2 have `varQuality` (`VQHIGH`, `VQLOW`) where 2.5 has `varFilter` (`PASS`, `VQLOW`, `AMBIGUOUS`, `;`-joined); `alleleFreq` and `alternativeCalls` arrive in 2.5. The name changed and so did the vocabulary.

```text
#FORMAT_VERSION	2.0
#GENERATED_BY	cgatools
#GENOME_REFERENCE	NCBI build 37
#SOFTWARE_VERSION	2.0.0.26
#TYPE	VAR-ANNOTATION

>locus	ploidy	allele	chromosome	begin	end	varType	reference	alleleSeq	varScoreVAF	varScoreEAF	varQuality	hapLink	xRef
1	2	all	chr1	0	10000	no-ref	=	?					
3	2	1	chr1	752565	752566	snp	G	A	120	120	VQHIGH		
3	2	2	chr1	752565	752566	ref	G	G	98	98	VQHIGH		
```

**`sniff` says:** `detected`, as `cgi_var`, build GRCh37

**What just-vcf does** (handled): Detection does not depend on the version: `#TYPE` and the first nine names of the `>` line are the same. The loader reads the names from the `>` line and takes columns by name; `varQuality` is read as `varFilter`, its values kept as written in `sample_filters` (no CGI document read here equates `VQHIGH` with `PASS`). The files loaded so far are format 2.0 and 2.5; `totalScore` is not mapped onto either 2.x score.

Code and tests: loaders/cgi.py `_FILTER_ALIASES`, `test_format_2_0_and_no_title`.

### Complete Genomics loci: 0-based spans, insertions of zero length, partly called alleles

<a id="shape-cgi-unanchored"></a>

CGI writes 0-based half-open spans. An insertion has `begin == end` and an empty `reference`; a deletion has an empty `alleleSeq`; ref and no-call blocks write `reference` as `=`; `?` is an unknown sequence, and a `?` inside a sequence (`GC?`) marks uncalled bases in a partly called allele. A haplotype row of varType `no-call` beside a called one is a half call.

```text
#FORMAT_VERSION	2.5
#GENERATED_BY	cgatools
#GENOME_REFERENCE	NCBI build 37
#SOFTWARE_VERSION	2.5.0.20
#TYPE	VAR-ANNOTATION

>locus	ploidy	allele	chromosome	begin	end	varType	reference	alleleSeq	varScoreVAF	varScoreEAF	varFilter	hapLink	xRef	alleleFreq	alternativeCalls
5	2	1	chr1	752600	752600	ins		T	60	60	PASS				
5	2	2	chr1	752600	752600	ref			60	60	PASS				
6	2	1	chr1	752700	752702	del	AG		55	55	PASS				
6	2	2	chr1	752700	752702	no-call	AG	A?							
```

**`sniff` says:** `detected`, as `cgi_var`, build GRCh37

**What just-vcf does** (handled): A span loads as `pos = begin + 1`, `end_pos = end`; an insertion point sits at `pos = begin` with `end_pos` null and `source_ref = ""`. Alleles stay unpadded, `""` for an empty one, until `normalize` adds the anchor base from the reference. A haplotype written `?` is a null in `gt_alleles`, a `?` inside a sequence is kept verbatim, and every no-call varType (`no-call`, `no-call-rc`, `no-call-ri`, `no-ref`, `PAR-called-in-X`) reads as no call.

Code and tests: `loaders/cgi.py`, `test_variant_loci`, `test_unanchored_indels_count_so_a_genome_does_not_read_snv_only`.

### Half-calls: one allele known, the other not

<a id="shape-half-calls"></a>

A genotype with one allele called and the other missing. VCF writes it `1|.` or `0/.` (Complete Genomics' VCF has thousands); FTDNA's Affymetrix-era exports write `-A`; Genes for Good writes `T-`; CGI's masterVar says zygosity `half`. A half-call is neither a no-call nor a homozygote.

```text
RSID,CHROMOSOME,POSITION,RESULT
rs4477212 1 72017 AG
rs3094315 1 742429 -A
rs3131972 1 742584 ---
rs12562034 1 758311 GG
```

**`sniff` says:** `detected`, as `ftdna_space`

**What just-vcf does** (handled): The called allele stays in its slot and the other is null: `-A` loads as `[null, "A"]`. A whole no-call spelling (`spec.no_call`) is one token: in formats that write two letters on every called row, FTDNA's `---` or a Final Report's single `-` gives two nulls. VCF half-calls keep their phase and are graded evidence.

Code and tests: `loaders/text.py`, `loaders/table.py`, `test_ftdna_half_calls_no_calls_and_hemizygous_letters`, `test_half_calls_phase_and_filters`.

### Column names in a meta line, 0-based positions and the empty allele `-`

<a id="shape-genos-columns-meta"></a>

Genos names its columns in a `##Columns=` meta line instead of a column line, writes 0-based half-open spans, and writes the empty allele as `-`: at an insertion site (`begin == end`) `-/-` is a homozygous call for no insertion, not a no-call. A single allele with no `/` is a hemizygous call on X or Y.

```text
##FileFormat=Genos
##ReferenceGenome=UCSC.hg19
##Columns=CHROM	BEGIN	END	ID	GENOTYPE
chr1	752565	752566	rs3094315	A/G
chr1	752720	752721	rs3131972	G/G
chr1	846808	846808	rs4475691	-/-
```

**`sniff` says:** `detected`, as `genos`, build GRCh37

**What just-vcf does** (handled): The loader reads the column names from `##Columns=` when the line is there (the spec's aliases match them) and falls back to the spec's order when it is not. A span moves to 1-based; an insertion site goes to the base before it, where VCF places it, with `end_pos` null until the anchor base is read. `-` loads as the empty allele `""`, and `source_gt` keeps the `-`.

Code and tests: `loaders/table.py`, `test_genos_spans_insertions_and_hemizygous_calls`, `test_genos_reads_every_field_with_or_without_its_column_line`.

### A column line repeated inside the file, and a delimiter that changes

<a id="shape-joined-exports-and-mixed-delimiters"></a>

Some FTDNA uploads are two exports joined, the column line repeated where the second starts. The 2010-era Affymetrix exports have a comma header over space-separated rows and end in a few hundred quoted comma rows. Rows are double-quoted in some exports and bare in others.

```text
RSID,CHROMOSOME,POSITION,RESULT
"rs4477212","1","72017","AA"
"rs3094315","1","742429","AG"
RSID,CHROMOSOME,POSITION,RESULT
"rs3131972","1","742584","GG"
```

**`sniff` says:** `detected`, as `ftdna`

**What just-vcf does** (handled): The repeated column line is skipped and counted. Every field is stripped of surrounding double quotes. An `ftdna_space` row is read by its own shape, so the quoted comma rows at its end load too.

Code and tests: `loaders/table.py`, `test_ftdna_files_joined_with_their_column_line_repeated`, `test_ftdna_space_export_with_a_quoted_comma_tail`.

### 0-based single positions and a build stated on every row

<a id="shape-build-in-the-rows"></a>

Color's genotype CSV writes each position one below its GRCh37 1-based position, states the build in a `ReferenceVersion` column on every row instead of in a header, and ends in a disclaimer trailer whose quoted paragraph holds commas.

```text
SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion
S1,2,136608645,rs4988235,AG,GRCh37
S1,12,112241765,rs671,GG,GRCh37
Raw data disclaimer:,,,,,
"This paragraph, like the real one, holds commas."
```

**`sniff` says:** `detected`, as `color_genotypes`, build GRCh37

**What just-vcf does** (handled): The loader adds 1 for every 0-based spec. The facet's build is the one every row states; rows that disagree leave it unknown, with a warning. Reading stops at the first row whose Position is empty.

Code and tests: `loaders/table.py`, `test_color_positions_are_shifted_to_one_based_and_the_trailer_is_cut`, `test_color_build_is_read_from_the_rows`.

### A build stated in another member of the same zip

<a id="shape-build-in-a-sibling"></a>

A deCODEme raw-data zip holds the genotype CSV and a small key,value member that states the build; the genotype member itself names none. Rows on the `-` strand are written on that strand.

Member `genotypes.csv`:

```text
Name,Variation,Chromosome,Position,Strand,YourCode
rs4477212,A/G,1,72017,+,AG
rs3094315,C/T,1,742429,-,CT
rs3131972,C/T,1,742584,-,TT
rs12562034,A/G,1,758311,+,GG
```

Member `info.txt`:

```text
Format version,1
SNP count,4
Reference sequence,NCBI Build 36 and Cambridge reference sequence
```

*Two members of one zip; member names are made up.*

**`sniff` says:** `detected`, as `decodeme`, inside `zip`, build unknown

**What just-vcf does** (handled): `sniff` answers for the largest detected member (`decodeme`, build unknown) and carries every member's result in `members`; the info member is `decodeme_info` with build GRCh36. The loader gives the deCODEme facet the build its zip's `decodeme_info` member states, and complements both letters on `-` rows to get the plus strand.

Code and tests: `loading.py`, `loaders/table.py`, `test_decodeme_build_comes_from_the_info_member_of_its_zip`, `test_decodeme_minus_strand_rows_are_complemented`.

### Column names after a section line, tab or comma by export

<a id="shape-section-header"></a>

An Illumina GenomeStudio Final Report opens with a `[Header]` block of key/value lines, then a `[Data]` line, then the column line. The delimiter is tab in `.txt` exports and comma in `.csv` ones, where the header block is padded with commas. Which allele columns appear depends on the export (Plus, Forward, Design or Top strand).

```text
[Header],,,
Content,GSA-24v2-0_A1.bpm,,
Num SNPs,665608,,
[Data],,,
SNP Name,Chr,MapInfo,Alleles - Plus
rs3094315,1,752566,AG
rs3131972,1,752721,GG
rs12124819,1,776546,-
```

**`sniff` says:** `detected`, as `illumina_final_report`

**What just-vcf does** (handled): The column names are read from the line after `[Data]`, and the delimiter off that line. The loader needs `Chr`, `MapInfo` and `Alleles - Plus`; an export with only Forward, Design or Top alleles is refused with `TableLayoutError`, since turning those into plus-strand letters needs the manifest. A report of several samples loads one at a time (`load(path, sample=…)`). A report with no column line and no positions is `illumina_final_report_headerless`, which does not load.

Code and tests: `loaders/table.py`, `test_final_report_tab_and_comma_exports_load_alike`, `test_final_report_without_plus_strand_alleles_is_refused`, `test_final_report_with_two_samples_and_unplaced_probes`.

### A column that only some versions write

<a id="shape-optional-columns"></a>

Promethease's `tab-genotype` report has a `Repute` column from some version on: 0.1.110 does not write it, 0.1.161 does. A column declared for every version is then missing from older files.

```text
#Promethease	tab-genotype	0.1.110	2013-01-01 00:00
#Column-names	
Name	Magnitude	PopulationFreq	Description
rs1333049(C;G)	None	None	
rs10490924(G;G)	0.0	None	normal
```

**`sniff` says:** `detected`, as `promethease_tab_genotype`

**What just-vcf does** (open): Detected as `promethease_tab_genotype` either way. A `ColumnSpec` cannot say a column appears only in some versions, so the spec lists `repute` and its notes say to take names from the column line. The format does not load: it is an annotation report without chromosome or position.

Code and tests: formats/spec.py `ColumnSpec`.

### A chromosome that is not in the rows

<a id="shape-chromosome-outside-the-row"></a>

Gencove's Oxford `.gen` imputation output ships one member per chromosome and repeats the chromosome only on typed rows; imputed rows write `---` in its place. FTDNA Big Y files have no chromosome column at all: every row is chrY.

```text
--- rs4477212 82154 A G 0.9 0.1 0
--- rs3094315 752566 G A 0.2 0.7 0.1
1 rs3131972 752721 A G 0 1 0
```

**`sniff` says:** `detected`, as `oxford_gen`

**What just-vcf does** (open): Detected, not loaded. Taking the chromosome from a member's name would trust a name, which `sniff` never does; the loader for these formats is not built.

Code and tests: formats/other_tables.py `OXFORD_GEN`, `FTDNA_BIGY`.

### A headerless BED file

<a id="shape-headerless-bed"></a>

Three tab-separated columns, chromosome, 0-based start and exclusive end, with no header.

```text
chr1	752565	752566
chr1	752720	752721
```

**`sniff` says:** `unknown`

**What just-vcf does** (open): `unknown`, not `not_variant_data`: the text kinds look at one line, and telling a BED from a headerless table needs a check over several lines that `sniff` does not make.

Code and tests: sniff.py `_TEXT_KINDS`.

### Array or imputed genotypes inside a VCF

<a id="shape-array-content-in-a-vcf"></a>

Illumina GenotypingToVCF, 23andMe conversions, Genes for Good and imputation servers write VCF too. The parse is a VCF's, but an absent site means unknown, never reference. GenotypingToVCF puts probes it cannot place at chromosome `NA`, POS 0.

```text
##fileformat=VCFv4.1
##workflow_type=Illumina_GenotypingToVCF
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	752566	rs3094315	G	A	.	PASS	.	GT	0/1
NA	0	rs3131972	N	A,G	.	GTEX	.	GT	./.
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): `sniff` says `vcf`; `Family` has no axis for array content. The loader's source detector reads the full header and sets the facet's `callset_type` to `array` or `imputed` where the header says so (here `##workflow_type`). polars-bio cannot read a record at POS 0, so those records are left out of the copy it reads and counted.

Code and tests: source.py `detect_source`, `loaders/vcf.py (POS 0)`, `test_load_states_array_content_and_warns_on_a_conflict`, `test_records_at_pos_zero_are_left_out_and_counted`.

### Rows of any width after fixed fields

<a id="shape-variable-width-rows"></a>

PGP-HMS's key=value export of a freeBayes+VarScan2 merge has twelve fixed fields, eight of them `KEY=value`, then an annotation tail of any length (23 to 84 fields per row in the files seen).

```text
chr1	752566	G	A	FBGenoType=0/1	FBRefAlleleCount=10	FBReferenceAlleleQ=400	FBTotalDepth=25	VSGenoType=0/1	VSRefAlleleCount=10	VSReferenceAlleleQ=35	VSTotalDepth=24	EFF_RESULTS	MODIFIER	KG=0.5
chr1	752721	A	G	FBGenoType=1/1	FBRefAlleleCount=.	FBReferenceAlleleQ=NULL	FBTotalDepth=20	VSGenoType=NULL	VSRefAlleleCount=NULL	VSReferenceAlleleQ=NULL	VSTotalDepth=NULL
```

**`sniff` says:** `detected`, as `freebayes_varscan2_keyvalue`

**What just-vcf does** (handled): Only the twelve fixed fields are read, each checked for its key; rows whose fields do not carry the keys in order are counted. The tail's keys, read from the first 1,000 rows, are listed in `facets.dropped_fields`. One facet per caller.

Code and tests: `loaders/keyvalue.py`, `test_one_facet_per_caller`.

### `##fileformat` not on the first line, and keys the header never declares

<a id="shape-fileformat-not-first"></a>

Genos exome VCFs sort their meta lines, so `##fileformat` can sit hundreds of lines down. Other files use INFO or FORMAT keys, GT included, that no meta line declares. polars-bio refuses the first and silently drops the undeclared keys of the second.

```text
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##fileformat=VCFv4.2
##reference=hg19
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	752566	rs3094315	G	A	.	PASS	DP=30	GT	0/1
```

*GT is used by the row and declared nowhere.*

**`sniff` says:** `detected`, as `vcf`, build GRCh37

**What just-vcf does** (handled): Detected as `vcf` (a `##fileformat` line anywhere in the head). The loader moves `##fileformat` to line 1 of the copy it reads and warns that the header was reordered; keys the rows use without a declaration are declared by the loader, named in a warning.

Code and tests: loaders/vcf.py `undeclared_keys`, `declare`, `test_header_repairs`, `test_an_undeclared_end_is_declared_so_a_block_keeps_its_extent`.

### A VCF body under a bare `#CHROM` line

<a id="shape-vcf-without-meta"></a>

A VCF with no `##` meta lines at all, the column line first.

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	752566	.	G	A	50.2	.	DP=12	GT:DP	0/1:12
```

**`sniff` says:** `detected`, as `vcf_no_meta`

**What just-vcf does** (handled): Its own format id, since the parse differs: the loader reads it under a synthesised header (`##fileformat=VCFv4.2` and a definition for each carried key the rows use), lists every other key in `dropped_fields`, leaves `vcf_version` null and warns that the file states no version, contigs or build.

Code and tests: `loading.py`, `loaders/vcf.py`, `test_vcf_no_meta_loads_under_a_synthesised_header`.

### A header that declares more samples than the rows carry

<a id="shape-declared-two-carries-one"></a>

PGP-HMS's freeBayes+VarScan2 merges declare two sample columns (`unknown`, `Sample1`) while most rows carry one. Which caller the one column belongs to is said by INFO/SF, not by its name.

```text
##fileformat=VCFv4.1
##source=freeBayes v0.9.21
##reference=/ref/hg19.fa
##source=VarScan2
##INFO=<ID=SF,Number=.,Type=String,Description="Source File">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##source_20160919.1=vcf-merge chr1.FB.vcf.gz chr1.VS.vcf.gz
##sourceFiles_20160919.1=0:chr1.FB.vcf.gz,1:chr1.VS.vcf.gz
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	unknown	Sample1
chr1	752566	.	G	A	40.1	PASS	SF=0,1	GT	0/1
chr1	776546	.	A	G	.	PASS	SF=1	GT	0/1
```

**`sniff` says:** `detected`, as `vcf_freebayes_varscan2`, build GRCh37

**What just-vcf does** (handled): polars-bio fills the missing second sample with None, which is a missing column, not a no-call. The loader returns one facet per caller (`relation = caller_split`), each caller identified by the `.FB.`/`.VS.` names on `##sourceFiles`; in the one-column layout the VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the shared ones.

Code and tests: `loading.py`, `loaders/vcf.py`, `test_two_declared_columns_one_carried`.

### Reference blocks in a VCF, and headers that only look like a gVCF

<a id="shape-gvcf-blocks"></a>

gVCF blocks are spelled ALT `<NON_REF>` (GATK), `<*>` (DeepVariant), or ALT `.` with INFO END (Illumina, Genos). `##ALT=<ID=NON_REF>` and `##GVCFBlock` survive into genotyped files with no blocks at all, and END also sits on SV records.

```text
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	10001	.	T	<*>	0	.	END=10468	GT	0/0
1	752566	rs3094315	G	A,<*>	50	PASS	.	GT	0/1
```

**`sniff` says:** `detected`, as `vcf`

**What just-vcf does** (handled): gVCF is not a format id: the header marks lie both ways. Whether a callset has blocks is read from its unfiltered rows, and the blocks go to the `blocks` frame, where they decide what an absent site means (`absence_means = blocks`).

Code and tests: `loaders/vcf.py`, `callset.py`, `test_blocks_are_all_star_records_asserting_reference`, `test_illumina_gvcf_blocks`, `test_load_keeps_blocks_deciding_absence_on_a_typed_gvcf`.

### A 23andMe-layout table written by a converter, with two column lines

<a id="shape-converter-column-lines"></a>

Converters emulate the 23andMe layout without its preamble. One opens `# This file emulates the format of 23andMe` and carries two column lines, the first with an empty genotype name; others put a blank line before the column line.

```text
# This file emulates the format of 23andMe
# rsid	chromosome	position	
# rsid	chromosome	position	genotype
rs4477212	1	82154	AG
rs3094315	1	752566	GG
```

**`sniff` says:** `detected`, as `rsid_table`, build unknown

**What just-vcf does** (handled): Detected as `rsid_table`, build unknown. The loader checks the last column line against the four names in order; every comment line is skipped.

Code and tests: loaders/rsid_table.py `column_line_ok`, `test_preamble_quirks_and_line_endings`.

### A stray byte before the first `#`

<a id="shape-stray-byte"></a>

One Genes for Good export opens `I# Genes for Good v1.0: …`: one byte before the comment mark.

```text
I# Genes for Good v1.0: filtered, unphased genotype calls.
# This data file generated by PLINK on: 01/01/16
# Human genome reference used: GRCh37/Mito:rCRS
# rsid	chromosome	position	genotype
rs4477212	1	82154	AG
rs3094315	1	752566	T-
```

**`sniff` says:** `detected`, as `genes_for_good`, build GRCh37

**What just-vcf does** (handled): The detection pattern allows one byte before the `#` of line 1. The loader's comment rule does the same (a line starting with `#`, `"#`, or one byte and then `#`), so the line is skipped.

Code and tests: loaders/rsid_table.py `_COMMENT`, `test_preamble_quirks_and_line_endings`.

### A 23andMe export re-saved by a spreadsheet as CSV, and cut short

<a id="shape-spreadsheet-23andme"></a>

A spreadsheet re-save pads every preamble line to four fields with `,,,` and quotes preamble lines that hold a comma, so they start with `"#`. The one real file stops mid chromosome 2 after about 65,500 lines: the spreadsheet's row limit cut it.

```text
# This data file generated by 23andMe at: Fri Jan  1 00:00:00 2016,,,
#,,,
"# More information on reference human assembly build 37 (a.k.a. Annotation Release 104), see:",,,
# We are using reference human assembly build 37 (also known as Annotation Release 104).,,,
# rsid,chromosome,position,genotype
rs4477212,1,82154,AG
rs3094315,1,752566,GG
```

**`sniff` says:** `detected`, as `23andme_csv`, build GRCh37

**What just-vcf does** (open): Its own format id, since the delimiter changed. Padded and quoted preamble lines are skipped as comments. The cut is not detected: the loader reads the rows that are there, and sites past the cut are absent, which in an array table means unknown, never reference.

Code and tests: formats/rsid_tables.py `TWENTYTHREE_AND_ME_CSV`, `test_a_resaved_23andme_table_is_read_by_its_internal_ids`.

### Byte-order marks, CRLF, bare CR and mixed line ends

<a id="shape-encodings-and-line-ends"></a>

Exports arrive in UTF-16 with a BOM, UTF-8 with a BOM, latin-1, with CRLF or bare-CR line ends, and with line ends mixed inside one file (AncestryDNA writes some preamble lines LF and the rest CRLF).

```text
#AncestryDNA raw data download
#This file was generated by AncestryDNA at: 01/01/2016 00:00:00 UTC
#Data was collected using AncestryDNA array version: V1.0
#of the SNP using human reference build 37.1 coordinates.
rsid	chromosome	position	allele1	allele2
rs4477212	1	82154	A	A
rs3094315	1	752566	A	G
```

*Line 1 ends in LF, the others in CRLF; shown here as plain line breaks.*

**`sniff` says:** `detected`, as `ancestrydna`, build GRCh37

**What just-vcf does** (handled): `sniff` decodes BOM-aware (UTF-16, UTF-8), falls back to latin-1, and turns CRLF and bare CR into LF before matching; the table loaders read lines the same way, so each line's own ending is stripped rather than one terminator chosen per file.

Code and tests: sniff.py `_decode`, loaders/text.py `lines`, `test_utf16_crlf_and_bom_exports_are_read`, `test_cr_only_line_endings_split_into_lines`, `test_ancestrydna_mixed_line_endings_and_bad_rows`.

### Findings that are not file shapes

- **F5**: Memory of the canonical sort on large VCFs: a cost of the load, not a file shape.
- **F11**: PL converted to `gl` one float step off: arithmetic in the loader, not a file shape.
