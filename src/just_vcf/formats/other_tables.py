"""Genotype-bearing files that are neither VCF, CGI, 23andMe-shaped nor a vendor rsid CSV.

Array exports straight out of Illumina GenomeStudio, imputation probabilities, annotator JSON, clinical
variant lists, a build-36 GFF, and the native tables of pre-VCF sequencing pipelines (CASAVA, SOAPsnp).
Every spec here was written against real PGP-HMS files; the notes record what those files do.

Several of these files carry a person's name in a column (CASAVA `Sample Name`) or a lab sample id. A
loader drops those columns; nothing here copies one.
"""

import polars as pl

from just_vcf.formats.spec import (
    BuildHint,
    ColumnSpec,
    Detection,
    Family,
    FormatSpec,
    HeaderMode,
)

# --- Illumina Annotation Engine / Nirvana JSON ------------------------------------------------------

NIRVANA_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "chromosome", "dtype": pl.String, "description": "positions[].chromosome; chr-prefixed, chrM for mtDNA"},
    {"name": "position", "dtype": pl.Int64, "description": "positions[].position; 1-based VCF POS"},
    {"name": "ref_allele", "dtype": pl.String, "description": "positions[].refAllele; VCF REF with anchor base"},
    {"name": "alt_alleles", "dtype": pl.List(pl.String), "description": "positions[].altAlleles; VCF ALT"},
    {"name": "filters", "dtype": pl.List(pl.String), "description": "positions[].filters; PASS, LowGQX…"},
    {"name": "genotype", "dtype": pl.String, "description": "positions[].samples[].genotype; VCF GT"},
    {"name": "variant_freq", "dtype": pl.Float64, "description": "positions[].samples[].variantFreq"},
    {"name": "total_depth", "dtype": pl.Int64, "description": "positions[].samples[].totalDepth"},
    {"name": "genotype_quality", "dtype": pl.Int64, "description": "positions[].samples[].genotypeQuality"},
    {"name": "allele_depths", "dtype": pl.List(pl.Int64), "description": "positions[].samples[].alleleDepths"},
)

NIRVANA_JSON = FormatSpec(
    id="nirvana_json",
    family=Family.ANNOTATED_VARIANTS,
    vendor="Illumina",
    description="Illumina Annotation Engine (Nirvana) JSON: a header object, then one position object per line.",
    detection=Detection(patterns=(r'\A\{"header":\{"annotator":"(Illumina Annotation Engine|Nirvana)\b',)),
    comment_prefix=None,
    delimiter=None,
    header=HeaderMode.NONE,
    columns=NIRVANA_COLUMNS,
    genotype_encoding=(
        "VCF GT string per sample (`0/1`, `1/1`, `1/2`, `0/0`); haploid calls are a single index (`1`, `0`), "
        "seen on chrM"
    ),
    build_hints=(
        BuildHint(r'"genomeAssembly":"(GRCh37|hg19)"', "GRCh37"),
        BuildHint(r'"genomeAssembly":"(GRCh38|hg38)"', "GRCh38"),
    ),
    notes=(
        "Not NDJSON: line 1 is the header object plus `\"positions\":[`; each following line is one position "
        "object ending in `,`; the file closes with `]}`. Stream it line by line, stripping the trailing comma.",
        "Other arrays can follow the positions: the PGP-HMS file closes them with `],\"omim\":[` and 13,343 "
        "OMIM entries before `]}`. The positions end at the first line starting with `]`.",
        "Schema v4 (Annotation Engine 1.5) headers do not state the assembly. `dataSources[].version` values "
        "such as phyloP `hg19` describe the annotation source, not the coordinates; they are not a build hint.",
        "Each position carries VCF-style REF/ALT; `variants[]` repeats them trimmed (`refAllele` `-` for an "
        "insertion, `begin` shifted by one). Load from the position level, not from `variants[]`.",
        "Filtered records are present (`filters` other than PASS, genotypeQuality 0); keep FILTER, do not drop.",
    ),
    example=(
        '{"header":{"annotator":"Illumina Annotation Engine 1.5.1.0","creationTime":"2016-12-02 12:36:47",'
        '"schemaVersion":4,"dataVersion":"84.22.36","samples":["S1"]},"positions":[\n'
        '{"chromosome":"chr1","refAllele":"G","position":752566,"quality":40.0,"filters":["PASS"],'
        '"altAlleles":["A"],"samples":[{"variantFreq":0.5,"totalDepth":30,"genotypeQuality":40,'
        '"alleleDepths":[15,15],"genotype":"0/1"}],"variants":[{"altAllele":"A","refAllele":"G",'
        '"begin":752566,"chromosome":"chr1","dbsnp":["rs3094315"],"end":752566,"variantType":"SNV",'
        '"vid":"1:752566:A"}]},\n'
        '{"chromosome":"chr1","refAllele":"A","position":752721,"quality":40.0,"filters":["PASS"],'
        '"altAlleles":["G"],"samples":[{"variantFreq":1,"totalDepth":28,"genotypeQuality":40,'
        '"alleleDepths":[0,28],"genotype":"1/1"}],"variants":[{"altAllele":"G","refAllele":"A",'
        '"begin":752721,"chromosome":"chr1","dbsnp":["rs3131972"],"end":752721,"variantType":"SNV",'
        '"vid":"1:752721:G"}]}\n'
        "]}\n"
    ),
)

# --- Illumina GenomeStudio FinalReport --------------------------------------------------------------

FINAL_REPORT_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "snp_name", "dtype": pl.String, "description": "probe name: an rsid or a manifest name "
     "(`1:103380393`, `CNGB3:c.991-3T>G`, `101SNP8856FG_A`)", "aliases": ("SNP Name",)},
    {"name": "sample_id", "dtype": pl.String, "description": "lab sample id", "aliases": ("Sample ID",)},
    {"name": "chromosome", "dtype": pl.String, "description": "1-22, X, Y, XY (pseudoautosomal), MT; "
     "0 = unmapped", "aliases": ("Chr",)},
    {"name": "position", "dtype": pl.Int64, "description": "manifest position, 1-based", "aliases": ("MapInfo",)},
    {"name": "genotype_plus", "dtype": pl.String, "description": "both alleles on the plus strand, one field",
     "aliases": ("Alleles - Plus",)},
    {"name": "allele1_forward", "dtype": pl.String, "description": "allele 1, dbSNP forward strand",
     "aliases": ("Allele1 - Forward",)},
    {"name": "allele2_forward", "dtype": pl.String, "description": "allele 2, dbSNP forward strand",
     "aliases": ("Allele2 - Forward",)},
    {"name": "allele1_design", "dtype": pl.String, "description": "allele 1, probe design strand",
     "aliases": ("Allele1 - Design",)},
    {"name": "allele2_design", "dtype": pl.String, "description": "allele 2, probe design strand",
     "aliases": ("Allele2 - Design",)},
    {"name": "gc_score", "dtype": pl.Float64, "description": "GenCall confidence, 0-1",
     "aliases": ("GC Score", "GenCall Score")},
    {"name": "log_r_ratio", "dtype": pl.Float64, "description": "intensity log R ratio; `nan` on no-calls",
     "aliases": ("LogR",)},
    {"name": "b_allele_freq", "dtype": pl.Float64, "description": "B allele frequency; `nan` on no-calls",
     "aliases": ("BFreq",)},
    {"name": "x", "dtype": pl.Float64, "description": "normalised X intensity", "aliases": ("X",)},
    {"name": "y", "dtype": pl.Float64, "description": "normalised Y intensity", "aliases": ("Y",)},
    {"name": "x_raw", "dtype": pl.Int64, "description": "raw X intensity", "aliases": ("Xraw",)},
    {"name": "y_raw", "dtype": pl.Int64, "description": "raw Y intensity", "aliases": ("Yraw",)},
)

_FINAL_REPORT_NOTES = (
    "Sections, not a header line: `[Header]` holds key/value lines (GSGT Version, Content = the .bpm "
    "manifest, Num SNPs…), `[Data]` is followed by the table. Skip to the line after `[Data]`.",
    "The delimiter is tab in the .txt exports and comma in the .csv ones; in the comma exports the "
    "`[Header]` block is padded with commas to the data width. Read the delimiter off the lines after `[Data]`.",
    "No build is stated. `Content` names the manifest (`GSA-24v2-0_A1.bpm`); mapping a manifest to a build "
    "is a vendor convention, not a header statement, so the build stays unknown.",
    "`Chr 0` marks probes the manifest leaves unmapped; `XY` is the pseudoautosomal region; `MT` is mtDNA.",
    "Many probe names are not rsids; `-XX` suffixes (`1:153941698-CT`) mark a second probe at one site.",
    "Indel probes are called `I`/`D` (`II`, `DD`, `ID`). The loader keeps them verbatim: resolving them "
    "needs the manifest's alleles and a REF (RM3).",
    "Files are often CRLF.",
)

ILLUMINA_FINAL_REPORT = FormatSpec(
    id="illumina_final_report",
    family=Family.GENOTYPE_TABLE,
    vendor="Illumina",
    description="GenomeStudio Final Report: `[Header]` key/value block, then `[Data]` with a column line.",
    detection=Detection(
        patterns=(
            r"\A\[Header\][\t,]*\r?$",
            r"^\[Data\][\t,]*\r?\n[^\r\n]*SNP Name[\t,]",
        ),
    ),
    comment_prefix=None,
    delimiter=None,
    header=HeaderMode.SECTION,
    section="Data",
    columns=FINAL_REPORT_COLUMNS,
    genotype_encoding=(
        "Report Wizard chooses the columns: either one two-letter field (`Alleles - Plus`) or one field per "
        "allele on a named strand (Forward, Design, Top, Plus). Forward and Design disagree by strand; only "
        "Plus is genomic plus strand. I/D on indel probes."
    ),
    no_call=("-",),
    notes=(
        *_FINAL_REPORT_NOTES,
        "The column set varies per export: seen are {SNP Name, Sample ID, Allele1/2 - Forward, Allele1/2 - "
        "Design, GC Score} (Omni 2.5) and {SNP Name, Chr, MapInfo, Alleles - Plus, GenCall Score, LogR, BFreq, "
        "X, Y, Xraw, Yraw} (GSA-24v2). `columns` is the union; map by name and keep the strand the file names.",
        "A no-call is `-` per allele (`-`/`-` in the per-allele layout) and a single `-` in `Alleles - Plus`, "
        "where every called row holds two letters, so the loader reads that `-` as two missing alleles.",
        "The loader needs `Chr`, `MapInfo` and `Alleles - Plus`. An export with only Forward, Design or Top "
        "alleles is refused (`TableLayoutError`): turning those into plus-strand letters needs the manifest.",
        "A report can hold several samples (`Sample ID`); `load(sample=…)` picks one, as for a VCF.",
        "Exports without Chr/MapInfo carry no coordinates: positions come from the manifest, and that "
        "manifest's build must be stated, not assumed.",
    ),
    example=(
        "[Header]\n"
        "Processing Date\t08/12/2018 10:10 PM\n"
        "Content\tGSA-24v2-0_A1.bpm\n"
        "Num SNPs\t665608\n"
        "Num Samples\t1\n"
        "[Data]\n"
        "SNP Name\tChr\tMapInfo\tAlleles - Plus\tGenCall Score\tLogR\tBFreq\tX\tY\tXraw\tYraw\n"
        "rs3094315\t1\t752566\tAG\t0.9154\t0.1658\t0.5012\t0.4518\t0.4692\t4302\t4589\n"
        "rs3131972\t1\t752721\tGG\t0.8983\t0.0588\t1.0\t0.0235\t0.8654\t337\t6299\n"
        "rs12124819\t1\t776546\t-\t0.0\tnan\tnan\tnan\tnan\t0\t0\n"
    ),
)

FINAL_REPORT_HEADERLESS_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "sample_id", "dtype": pl.String, "description": "lab sample id"},
    {"name": "snp_name", "dtype": pl.String, "description": "probe name: an rsid or a Y/mt marker name"},
    {"name": "chromosome", "dtype": pl.String, "description": "1-22, X, Y; 0 = unmapped (most mt probes)"},
    {"name": "allele1", "dtype": pl.String, "description": "allele 1; strand not stated"},
    {"name": "allele2", "dtype": pl.String, "description": "allele 2; strand not stated"},
)

ILLUMINA_FINAL_REPORT_HEADERLESS = FormatSpec(
    id="illumina_final_report_headerless",
    family=Family.GENOTYPE_TABLE,
    vendor="Illumina",
    description=(
        "GenomeStudio Final Report exported without a column line after `[Data]` (Genographic Geno 2.0 "
        "`NGS_iSelect` chip): sample, probe, chromosome, two alleles."
    ),
    detection=Detection(
        patterns=(
            r"\A\[Header\][\t,]*\r?$",
            r"^\[Data\],*\r?\n[^,\r\n]+,[^,\r\n]+,(0|[1-9]|1\d|2[0-2]|X|Y|XY|MT),[ACGTID-],[ACGTID-]\r?$",
        ),
    ),
    comment_prefix=None,
    delimiter=",",
    header=HeaderMode.NONE,
    columns=FINAL_REPORT_HEADERLESS_COLUMNS,
    genotype_encoding="one allele letter per column; I/D on indel probes; strand not stated by the file",
    no_call=("-",),
    notes=(
        *_FINAL_REPORT_NOTES,
        "The five columns are inferred from content (sample id, probe, chromosome, allele, allele); the file "
        "names none of them. There is no position column: positions and strand come from the manifest.",
        "The chip is a Y/mtDNA-heavy ancestry panel; many probe names are marker names (`CTS100`) that "
        "need a marker table to be placed.",
    ),
    example=(
        "[Header],,,,\n"
        "GSGT Version,1.9.4,,,\n"
        "Content,,NGS_iSelect_v1_15030891_2012_B-wRS-2.bpm,,\n"
        "Num SNPs,154476,,,\n"
        "[Data],,,,\n"
        "S1_A01,rs3094315,1,A,G\n"
        "S1_A01,rs3131972,1,G,G\n"
        "S1_A01,rs12124819,1,-,-\n"
    ),
)

# --- GET-Evidence genome GFF ------------------------------------------------------------------------

GFF_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "chromosome", "dtype": pl.String, "description": "GFF seqid; chr-prefixed"},
    {"name": "source", "dtype": pl.String, "description": "GFF source: `-`, or `CGI` (not a CGI var file)"},
    {"name": "type", "dtype": pl.String, "description": "GFF type: `-`, `SNP` or `REF`"},
    {"name": "start", "dtype": pl.Int64, "description": "1-based, closed"},
    {"name": "end", "dtype": pl.Int64, "description": "1-based, closed; equals start for a substitution"},
    {"name": "score", "dtype": pl.String, "description": "`.`"},
    {"name": "strand", "dtype": pl.String, "description": "`+`"},
    {"name": "phase", "dtype": pl.String, "description": "`.`"},
    {"name": "attributes", "dtype": pl.String, "description": "`alleles A/G;db_xref dbsnp:rs…;ref_allele A`, "
     "`alleles C`, or `.` on REF rows"},
)

GET_EVIDENCE_GFF = FormatSpec(
    id="get_evidence_gff",
    family=Family.ANNOTATED_VARIANTS,
    vendor="PGP GET-Evidence",
    description="GET-Evidence genome GFF: a `##genome-build` line, then one site per line with an `alleles` attribute.",
    detection=Detection(
        patterns=(r"^##genome-build b?3[678]\b",),
        data_shape=(
            r"(chr)?([0-9]{1,2}|X|Y|M|MT)\t[^\t]+\t[^\t]+\t\d+\t\d+\t[^\t]+\t[+\-.]\t[^\t]+"
            r"\t(alleles [ACGTN]*(/[ACGTN]+)*(;[^\t]*)?|\.)"
        ),
    ),
    comment_prefix="#",
    delimiter="\t",
    header=HeaderMode.NONE,
    columns=GFF_COLUMNS,
    genotype_encoding=(
        "the `alleles` attribute: `X/Y` is heterozygous, one letter is homozygous for it; `ref_allele`, "
        "when present, states REF. A `REF` row (attributes `.`) states hom-ref at that site"
    ),
    no_call=("",),
    build_hints=(
        BuildHint(r"^##genome-build b?36\b", "GRCh36"),
        BuildHint(r"^##genome-build b?37\b", "GRCh37"),
        BuildHint(r"^##genome-build b?38\b", "GRCh38"),
    ),
    notes=(
        "Two layouts seen. Older (`*_build_36_substitution_variants.gff`): source and type `-`, attributes "
        "only `alleles`, no REF. Newer (`##genome-build b37`): source `CGI`, type `SNP` or `REF`, attributes "
        "`alleles`, `db_xref dbsnp:rs…` (may list several), `ref_allele`, sometimes `amino_acid` and "
        "`ucsc_trans`.",
        "Without `ref_allele`, REF must come from the stated build's FASTA, never from the `alleles` letters; "
        "with it, the FASTA is still the check.",
        "`REF` rows state hom-ref at a listed site; that is evidence for that site only, not for unlisted ones.",
        "One row seen with an empty `alleles ` value; treat it as a no-call, not as hom-ref.",
        "The build line spells the build `36` or `b37`. A source of `CGI` does not make this a CGI var file.",
        "Rows may be sorted lexicographically (chr10, chr11, …), not karyotypically.",
    ),
    example=(
        "##genome-build b37\n"
        "# File creation date: 2019-01-15 12:00:00.000000\n"
        "chr1\tCGI\tREF\t752566\t752566\t.\t+\t.\t.\n"
        "chr1\tCGI\tSNP\t752721\t752721\t.\t+\t.\talleles A/G;db_xref dbsnp:rs3131972;ref_allele A\n"
        "chr1\tCGI\tSNP\t776546\t776546\t.\t+\t.\talleles G;db_xref dbsnp:rs12124819;ref_allele A\n"
    ),
)

# --- Color --------------------------------------------------------------------------------------------

COLOR_VARIANT_SET_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "sample_id", "dtype": pl.String, "description": "Color sample or kit id", "aliases": ("SampleID",)},
    {"name": "chromosome", "dtype": pl.String, "description": "1-22, X, Y; no chr prefix", "aliases": ("Chromosome",)},
    {"name": "position", "dtype": pl.Int64, "description": "1-based", "aliases": ("Position", "Pos")},
    {"name": "ref", "dtype": pl.String, "description": "reference allele as stated by Color", "aliases": ("Ref",)},
    {"name": "alt", "dtype": pl.String, "description": "alternate allele", "aliases": ("Alt",)},
    {"name": "gene", "dtype": pl.String, "description": "HGNC symbol", "aliases": ("Gene",)},
    {"name": "type", "dtype": pl.String, "description": "`SNV` or `indel`", "aliases": ("Type",)},
    {"name": "c_hgvs", "dtype": pl.String, "description": "cDNA HGVS", "aliases": ("cHGVS",)},
    {"name": "p_hgvs", "dtype": pl.String, "description": "protein HGVS; empty off-coding", "aliases": ("pHGVS",)},
    {"name": "transcript", "dtype": pl.String, "description": "Ensembl transcript", "aliases": ("Transcript",)},
    {"name": "bic", "dtype": pl.String, "description": "BIC nomenclature (BRCA1/2); usually empty", "aliases": ("BIC",)},
    {"name": "zygosity", "dtype": pl.String, "description": "`Heterozygous` or `Homozygous`", "aliases": ("Zygosity",)},
    {"name": "structural_variant", "dtype": pl.String, "description": "`no` or empty",
     "aliases": ("StructuralVariant", "isStructuralVariant")},
    {"name": "classification", "dtype": pl.String, "description": "ACMG-style class (Benign, Likely Benign…)",
     "aliases": ("Classification",)},
    {"name": "comments", "dtype": pl.String, "description": "free text; usually empty", "aliases": ("Comments",)},
    {"name": "test_type", "dtype": pl.String, "description": "test or panel name(s)", "aliases": ("TestType", "TestOrdered")},
    {"name": "reference_version", "dtype": pl.String, "description": "`GRCh37` or `GRCh37.p12`",
     "aliases": ("ReferenceVersion",)},
    {"name": "analysis_completed", "dtype": pl.Date, "description": "ISO date", "aliases": ("AnalysisCompleted",)},
)

COLOR_VARIANT_SET = FormatSpec(
    id="color_variant_set",
    family=Family.ANNOTATED_VARIANTS,
    vendor="Color",
    description="Color Health `variant_set.csv`: the panel's reportable variants with HGVS and classification.",
    detection=Detection(patterns=(r"^SampleID,Chromosome,Pos(ition)?,Ref,Alt,Gene,Type,cHGVS,pHGVS,",)),
    comment_prefix=None,
    delimiter=",",
    header=HeaderMode.LINE,
    columns=COLOR_VARIANT_SET_COLUMNS,
    genotype_encoding="Ref and Alt columns plus `Zygosity` (`Heterozygous`, `Homozygous`)",
    build_hints=(
        BuildHint(r",GRCh37(\.p\d+)?(,|\r?$)", "GRCh37"),
        BuildHint(r",GRCh38(\.p\d+)?(,|\r?$)", "GRCh38"),
    ),
    notes=(
        "Two layouts. Older (16 columns): `Pos`, an `Exon` column, `isStructuralVariant`, `TestOrdered`, "
        "no BIC/Comments/AnalysisCompleted. Newer (18 columns, declared here): `Position`, `BIC`, "
        "`StructuralVariant`, `Comments`, `TestType`, `AnalysisCompleted`. Map by name.",
        "The build is stated per row in `ReferenceVersion`, not in a header line; the build hints read it.",
        "A footer follows the data: blank line(s), `Raw data disclaimer:` and a quoted paragraph that may "
        "span lines. Stop at the first row whose SampleID cell is empty or `Raw data disclaimer:`.",
        "Only the panel's reportable variants are listed; an absent site is unknown, not reference.",
        "CRLF line ends.",
        "Measured 2026-10-23: `Ref` equals GRCh37's bases at `Position` on all 262 rows of the two PGP-HMS "
        "files, so positions are 1-based, unlike Color's `discovery_genotypes.csv`.",
    ),
    example=(
        "SampleID,Chromosome,Position,Ref,Alt,Gene,Type,cHGVS,pHGVS,Transcript,BIC,Zygosity,StructuralVariant,"
        "Classification,Comments,TestType,ReferenceVersion,AnalysisCompleted\n"
        "S1,1,55523033,A,G,PCSK9,SNV,c.1026A>G,p.Gln342=,ENST00000302118,,Homozygous,,Benign,,"
        "Color Hereditary Disease Risk Test 2.0,GRCh37.p12,2025-06-23\n"
        "S1,1,55524237,G,A,PCSK9,SNV,c.1420G>A,p.Val474Ile,ENST00000302118,,Heterozygous,,Benign,,"
        "Color Hereditary Disease Risk Test 2.0,GRCh37.p12,2025-06-23\n"
    ),
)

# --- FamilyTreeDNA Big Y ------------------------------------------------------------------------------

FTDNA_BIGY_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "type", "dtype": pl.String, "description": "`Known SNP` or `Novel Variant`", "aliases": ("Type",)},
    {"name": "position", "dtype": pl.Int64, "description": "chrY position; empty on Known SNP rows",
     "aliases": ("Position",)},
    {"name": "snp_name", "dtype": pl.String, "description": "Y-SNP name (`BY1181`); empty on novel rows",
     "aliases": ("SNPName",)},
    {"name": "derived", "dtype": pl.String, "description": "`Yes (+)`, `No (-)`, `?`; empty on novel rows",
     "aliases": ("Derived",)},
    {"name": "on_tree", "dtype": pl.String, "description": "`Yes`/`No`; empty on novel rows", "aliases": ("OnTree",)},
    {"name": "reference", "dtype": pl.String, "description": "reference base as stated by FTDNA",
     "aliases": ("Reference",)},
    {"name": "genotype", "dtype": pl.String, "description": "one base (haploid) or `?`", "aliases": ("Genotype",)},
)

FTDNA_BIGY = FormatSpec(
    id="ftdna_bigy",
    family=Family.VARIANT_CALLS,
    vendor="FamilyTreeDNA",
    description="FamilyTreeDNA Big Y results CSV: chrY sites with reference and observed base.",
    detection=Detection(
        patterns=(r'^"?Type"?,"?Position"?,"?SNPName"?,"?Derived"?,"?OnTree"?,"?Reference"?,"?Genotype"?\r?$',),
    ),
    comment_prefix=None,
    delimiter=",",
    header=HeaderMode.LINE,
    columns=FTDNA_BIGY_COLUMNS,
    genotype_encoding="one base per site (chrY is haploid); `?` where the site was not called",
    no_call=("?",),
    notes=(
        "Big Y is targeted sequencing of chrY, not a genotyping array: its rows are sequencing calls with the "
        "reference base stated, as CASAVA's and SOAPsnp's are, so the family is `variant_calls` and "
        "`chip_absence` refuses it as not an array export.",
        "There is no chromosome column: every row is chrY, a structural fact of the product.",
        "Known SNP rows have no position; they are placed only through a Y-SNP name table. Novel Variant "
        "rows have a position and mostly a genotype equal to the reference.",
        "Neither build nor coordinate convention is stated; Big Y moved from hg19 to hg38 over its "
        "lifetime, so the build stays unknown.",
        "Every field is double-quoted in some exports and bare in others.",
        "Measured 2026-10-23: `Reference` equals GRCh38's base at `Position` on every sampled row (about "
        "2,000 each) of the two PGP-HMS exports, and GRCh37's on about 20%, so those exports are 1-based on "
        "GRCh38. That is a measurement, not a statement: the loader still leaves the build unknown.",
    ),
    example=(
        '"Type","Position","SNPName","Derived","OnTree","Reference","Genotype"\n'
        '"Novel Variant","10061836","","","","A","A"\n'
        '"Known SNP","","A14825","?","Yes","A","?"\n'
        '"Novel Variant","10061838","","","","A","?"\n'
    ),
)

# --- Illumina CASAVA tables (exported with Index and Sample columns) -------------------------------

CASAVA_SNP_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "index", "dtype": pl.Int64, "description": "row index across the export", "aliases": ("Index",)},
    {"name": "sample_name", "dtype": pl.String, "description": "free text, often a person's name; drop on load",
     "aliases": ("Sample Name",)},
    {"name": "chromosome", "dtype": pl.String, "description": "1-22, X, Y, M; no chr prefix", "aliases": ("Chr Id",)},
    {"name": "position", "dtype": pl.Int64, "description": "1-based", "aliases": ("Position",)},
    {"name": "a_bases", "dtype": pl.Int64, "description": "reads calling A", "aliases": ("A Bases",)},
    {"name": "c_bases", "dtype": pl.Int64, "description": "reads calling C", "aliases": ("C Bases",)},
    {"name": "g_bases", "dtype": pl.Int64, "description": "reads calling G", "aliases": ("G Bases",)},
    {"name": "t_bases", "dtype": pl.Int64, "description": "reads calling T", "aliases": ("T Bases",)},
    {"name": "call", "dtype": pl.String, "description": "one letter (hom) or two (het)", "aliases": ("Call",)},
    {"name": "total_bases", "dtype": pl.Int64, "description": "depth", "aliases": ("Total Bases",)},
    {"name": "bases_used", "dtype": pl.Int64, "description": "depth used for the call", "aliases": ("Bases Used",)},
    {"name": "call1_score", "dtype": pl.Float64, "description": "score of the best call", "aliases": ("Call 1 Score",)},
    {"name": "call2_score", "dtype": pl.Float64, "description": "score of the second call",
     "aliases": ("Call 2 Score",)},
    {"name": "ref_base", "dtype": pl.String, "description": "reference base as stated", "aliases": ("Ref Base",)},
    {"name": "call_type", "dtype": pl.String, "description": "`Diff`, `Het1`, `Het2`, `Other`",
     "aliases": ("Call Type",)},
    {"name": "dbsnp", "dtype": pl.String, "description": "rsid or empty", "aliases": ("dbSNP SNP Name",)},
    {"name": "genes", "dtype": pl.String, "description": "UniProt accession(s) or empty", "aliases": ("Genes",)},
)

CASAVA_SNP_TABLE = FormatSpec(
    id="casava_snp_table",
    family=Family.VARIANT_CALLS,
    vendor="Illumina",
    description="Illumina CASAVA-era SNP calls with per-base read counts, one file per chromosome.",
    detection=Detection(
        patterns=(r"^Index\tSample Name\tChr Id\tPosition\tA Bases\tC Bases\tG Bases\tT Bases\tCall\t",),
    ),
    comment_prefix=None,
    delimiter="\t",
    header=HeaderMode.LINE,
    columns=CASAVA_SNP_COLUMNS,
    genotype_encoding=(
        "`Call` letters with `Call Type`: `Diff` = one letter, homozygous non-reference; `Het1` = two letters, "
        "first equals Ref Base; `Het2` = two letters, second equals Ref Base; `Other` = two letters, neither "
        "is the reference"
    ),
    notes=(
        "Only non-reference calls are listed; there are no hom-ref or no-call rows, so an absent site is unknown.",
        "`Sample Name` holds free text, in real files a person's name. Drop it on load; never log it.",
        "`Chr Id` spells mtDNA `M`; map it to `MT` explicitly.",
        "No build is stated; the producer is inferred (CASAVA read-count calls), the file does not name it.",
        "Measured 2026-10-23: `Ref Base` equals NCBI36's base at `Position` on every sampled row (2,000 of "
        "one PGP-HMS table), so positions are 1-based; the build is still not stated.",
        "CRLF line ends. `Genes` and `dbSNP SNP Name` are often empty, leaving trailing tabs.",
    ),
    example=(
        "Index\tSample Name\tChr Id\tPosition\tA Bases\tC Bases\tG Bases\tT Bases\tCall\tTotal Bases\tBases Used\t"
        "Call 1 Score\tCall 2 Score\tRef Base\tCall Type\tdbSNP SNP Name\tGenes\n"
        "1\tS1\t1\t752566\t14\t0\t12\t0\tGA\t26\t26\t120.5\t98.2\tG\tHet1\trs3094315\t\n"
        "2\tS1\t1\t752721\t0\t0\t31\t0\tG\t31\t31\t210.4\t0.00\tA\tDiff\trs3131972\t\n"
    ),
)

CASAVA_INDEL_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "index", "dtype": pl.Int64, "description": "row index", "aliases": ("Index",)},
    {"name": "sample_name", "dtype": pl.String, "description": "free text, often a person's name; drop on load",
     "aliases": ("Sample",)},
    {"name": "chromosome", "dtype": pl.String, "description": "1-22, X, Y, M; no chr prefix", "aliases": ("Chr",)},
    {"name": "position", "dtype": pl.Int64, "description": "CASAVA indel position", "aliases": ("pos",)},
    {"name": "cigar", "dtype": pl.String, "description": "`1D`, `3I`…", "aliases": ("CIGAR",)},
    {"name": "ref_upstream", "dtype": pl.String, "description": "10 reference bases before the event",
     "aliases": ("ref_upstream",)},
    {"name": "ref_indel", "dtype": pl.String, "description": "`G/-` deletion, `-/T` insertion",
     "aliases": ("ref/indel",)},
    {"name": "ref_downstream", "dtype": pl.String, "description": "10 reference bases after the event",
     "aliases": ("ref_downstream",)},
    {"name": "q_indel", "dtype": pl.Int64, "description": "indel quality", "aliases": ("Q(indel)",)},
    {"name": "max_gtype", "dtype": pl.String, "description": "`het` or `hom`", "aliases": ("max_gtype",)},
    {"name": "q_max_gtype", "dtype": pl.Int64, "description": "genotype quality", "aliases": ("Q(max_gtype)",)},
    {"name": "max2_gtype", "dtype": pl.String, "description": "second-best genotype: `het`, `hom`, `ref`",
     "aliases": ("max2_gtype",)},
    {"name": "bp1_reads", "dtype": pl.Int64, "description": "reads spanning the left breakpoint",
     "aliases": ("bp1_reads",)},
    {"name": "ref_reads", "dtype": pl.Int64, "description": "reads supporting the reference", "aliases": ("ref_reads",)},
    {"name": "indel_reads", "dtype": pl.Int64, "description": "reads supporting the indel", "aliases": ("indel_reads",)},
    {"name": "other_reads", "dtype": pl.Int64, "description": "other reads", "aliases": ("other_reads",)},
    {"name": "repeat_unit", "dtype": pl.String, "description": "repeat unit at the site", "aliases": ("repeat_unit",)},
    {"name": "ref_repeat_count", "dtype": pl.Int64, "description": "repeat count in the reference",
     "aliases": ("ref_repeat_count",)},
    {"name": "indel_repeat_count", "dtype": pl.Int64, "description": "repeat count with the indel",
     "aliases": ("indel_repeat_count",)},
)

CASAVA_INDEL_TABLE = FormatSpec(
    id="casava_indel_table",
    family=Family.VARIANT_CALLS,
    vendor="Illumina",
    description="Illumina CASAVA `indels.txt` columns, exported with Index and Sample columns in front.",
    detection=Detection(
        patterns=(
            r"^Index\tSample\tChr\tpos\tCIGAR\tref_upstream\tref/indel\tref_downstream\tQ\(indel\)\tmax_gtype\t",
        ),
    ),
    comment_prefix=None,
    delimiter="\t",
    header=HeaderMode.LINE,
    columns=CASAVA_INDEL_COLUMNS,
    genotype_encoding="`max_gtype` `het`/`hom` for the event in `ref/indel`; `-` marks the empty side",
    notes=(
        "Alleles have no anchor base (`G/-`, `-/T`). A VCF record needs the anchor from the reference FASTA; "
        "`ref_upstream` gives 10 bases the loader can check the FASTA and the position convention against.",
        "Only called indels are listed; an absent site is unknown.",
        "`Sample` holds free text, in real files a person's name. Drop it on load; never log it.",
        "No build is stated. `Chr` spells mtDNA `M`; map it to `MT` explicitly. CRLF line ends.",
        "Measured 2026-10-23 against NCBI36 (2,000 sampled rows of the one PGP-HMS table): `ref_upstream` ends "
        "at the base before `pos`; a deletion's bases start at `pos`, with `ref_downstream` after them; an "
        "insertion's `ref_downstream` starts at `pos`. So `pos` is 1-based, the span's first base or the base "
        "after the insertion point. `ref/indel` pads with `-` on both sides of complex events (`C-A/-TA`).",
    ),
    example=(
        "Index\tSample\tChr\tpos\tCIGAR\tref_upstream\tref/indel\tref_downstream\tQ(indel)\tmax_gtype\t"
        "Q(max_gtype)\tmax2_gtype\tbp1_reads\tref_reads\tindel_reads\tother_reads\trepeat_unit\t"
        "ref_repeat_count\tindel_repeat_count\n"
        "0\tS1\t20\t10000012\t1D\tGTTTACTACA\tG/-\tCTTTGTAGTA\t3070\thom\t76\thet\t33\t0\t23\t10\tG\t1\t0\n"
        "1\tS1\t20\t10000058\t1I\tCTCTAACTTT\t-/T\tGTTTGTTTTT\t2887\thet\t33\thom\t17\t9\t11\t6\tT\t3\t4\n"
    ),
)

# --- Fallbacks: recognised by data shape alone -------------------------------------------------------

OXFORD_GEN_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "snp_id", "dtype": pl.String, "description": "`---` on imputed sites; on typed sites the input's "
     "first column (the chromosome in Gencove exports, `23` for X)"},
    {"name": "rsid", "dtype": pl.String, "description": "rsid, or `chr1:13957:D` / `MERGED_DEL_…` for others"},
    {"name": "position", "dtype": pl.Int64, "description": "1-based"},
    {"name": "allele_a", "dtype": pl.String, "description": "first allele; `-` for a deleted allele"},
    {"name": "allele_b", "dtype": pl.String, "description": "second allele; `-` for a deleted allele"},
    {"name": "p_aa", "dtype": pl.Float64, "description": "probability of genotype AA"},
    {"name": "p_ab", "dtype": pl.Float64, "description": "probability of genotype AB"},
    {"name": "p_bb", "dtype": pl.Float64, "description": "probability of genotype BB"},
)

OXFORD_GEN = FormatSpec(
    id="oxford_gen",
    family=Family.GENOTYPE_TABLE,
    vendor="generic (IMPUTE2 / Oxford .gen; Gencove)",
    description="Oxford `.gen` genotype probabilities, one site per line, space-separated.",
    detection=Detection(
        data_shape=(
            r"\S+ \S+ \d+ ([ACGTN]+|-) ([ACGTN]+|-)"
            r"( (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+) (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+)"
            r" (0|1|0?\.\d+|1\.0+|\d(\.\d+)?e-?\d+))+"
        ),
    ),
    comment_prefix=None,
    delimiter=" ",
    header=HeaderMode.NONE,
    columns=OXFORD_GEN_COLUMNS,
    genotype_encoding=(
        "three probabilities per sample for AA, AB, BB over the listed alleles A and B; graded evidence, "
        "not a hard call. A and B are not REF and ALT: REF comes from the reference"
    ),
    no_call=("0 0 0",),
    notes=(
        "The chromosome is not in the rows: Gencove ships one member per chromosome (`id_…_chr1.gen`) and only "
        "typed rows repeat it in column 1. The loader must take it from the member name or refuse; sniff "
        "never trusts names.",
        "Column 1 `---` marks an imputed site and anything else a typed one; keep that distinction.",
        "Several samples make 3 probability columns each; `columns` declares one sample.",
        "No build is stated; the positions seen match GRCh37, but that is not a header statement.",
        "`0 0 0` is the Oxford no-call spelling; probabilities are rounded and need not sum to 1.",
    ),
    example=(
        "--- rs58108140 10583 G A 0.638 0.359 0.004\n"
        "--- rs189107123 10611 C G 0.981 0.019 0\n"
        "1 rs3131972 752721 A G 0 0 1\n"
        "1 rs12124819 776546 A G 0 1 0\n"
    ),
)

SOAPSNP_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "chromosome", "dtype": pl.String, "description": "chr-prefixed"},
    {"name": "position", "dtype": pl.Int64, "description": "1-based"},
    {"name": "ref", "dtype": pl.String, "description": "reference base as stated"},
    {"name": "consensus", "dtype": pl.String, "description": "IUPAC genotype code"},
    {"name": "consensus_quality", "dtype": pl.Int64, "description": "Phred quality of the consensus"},
    {"name": "best_base", "dtype": pl.String, "description": "best-supported base"},
    {"name": "best_base_quality", "dtype": pl.Int64, "description": "mean quality of best base"},
    {"name": "best_base_unique_reads", "dtype": pl.Int64, "description": "uniquely mapped reads, best base"},
    {"name": "best_base_all_reads", "dtype": pl.Int64, "description": "all mapped reads, best base"},
    {"name": "second_base", "dtype": pl.String, "description": "second-best base"},
    {"name": "second_base_quality", "dtype": pl.Int64, "description": "mean quality of second base"},
    {"name": "second_base_unique_reads", "dtype": pl.Int64, "description": "uniquely mapped reads, second base"},
    {"name": "second_base_all_reads", "dtype": pl.Int64, "description": "all mapped reads, second base"},
    {"name": "depth", "dtype": pl.Int64, "description": "sequencing depth"},
    {"name": "rank_sum_p", "dtype": pl.Float64, "description": "rank-sum test p-value"},
    {"name": "copy_number", "dtype": pl.Float64, "description": "mean copy number of the region; 255 = none"},
    {"name": "dbsnp", "dtype": pl.Int64, "description": "1 when the site is a known dbSNP site"},
    {"name": "nearest_snp_distance", "dtype": pl.Int64, "description": "distance to the nearest listed site"},
)

SOAPSNP_CNS = FormatSpec(
    id="soapsnp_cns",
    family=Family.VARIANT_CALLS,
    vendor="BGI (SOAPsnp)",
    description="SOAPsnp consensus output, SNP rows only, 17 SOAPsnp columns plus distance to the nearest SNP.",
    detection=Detection(
        data_shape=(
            r"chr([0-9]{1,2}|X|Y|M)\t\d+\t[ACGTN]\t[ACGTRYKMSWBDHVN]\t\d+\t[ACGTN]\t\d+\t\d+\t\d+"
            r"\t[ACGTN]\t\d+\t\d+\t\d+\t\d+\t\d+(\.\d+)?\t\d+(\.\d+)?\t[01]\t\d+"
        ),
    ),
    comment_prefix=None,
    delimiter="\t",
    header=HeaderMode.NONE,
    columns=SOAPSNP_COLUMNS,
    genotype_encoding=(
        "IUPAC code in `consensus`: a plain letter is homozygous for that base, R/Y/K/M/S/W heterozygous "
        "for the pair it stands for"
    ),
    notes=(
        "The producer is inferred from the layout: the file names no tool. Column 18 equals the distance to "
        "the neighbouring row on the same chromosome.",
        "Rows with depth 0 and consensus quality 1 are prior-driven calls at dbSNP sites (`dbsnp` 1) with no "
        "read evidence; the loader must carry depth and quality through, never treat them as observed.",
        "Only SNP rows are present; an absent site is unknown.",
        "No build is stated. `chrM` must be mapped to `MT` explicitly.",
        "Measured 2026-10-23: `ref` equals NCBI36's base at `position` on every sampled row (2,000 of the one "
        "PGP-HMS file), so positions are 1-based; the build is still not stated.",
    ),
    example=(
        "chr1\t742429\tG\tR\t34\tG\t36\t12\t12\tA\t30\t10\t10\t22\t0.412\t1.00000\t1\t15882\n"
        "chr1\t758311\tG\tA\t87\tA\t38\t19\t19\tG\t0\t0\t0\t19\t1.00000\t1.00000\t1\t15882\n"
        "chr1\t995669\tC\tY\t41\tC\t35\t9\t9\tT\t33\t8\t8\t17\t0.731\t1.00000\t1\t237358\n"
    ),
)

SPECS: tuple[FormatSpec, ...] = (
    NIRVANA_JSON,
    ILLUMINA_FINAL_REPORT,
    ILLUMINA_FINAL_REPORT_HEADERLESS,
    GET_EVIDENCE_GFF,
    COLOR_VARIANT_SET,
    FTDNA_BIGY,
    CASAVA_SNP_TABLE,
    CASAVA_INDEL_TABLE,
)
FALLBACKS: tuple[FormatSpec, ...] = (OXFORD_GEN, SOAPSNP_CNS)
