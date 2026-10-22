"""VCF, gVCF and BCF, and the vendor files that are VCF underneath.

polars-bio reads these; the specs carry detection and the notes a loader needs. A separate id exists
only where the loader must parse differently from a plain VCF read: Complete Genomics' record classes,
PGP's merged freeBayes+VarScan2 calls, a `#CHROM` line with no meta lines above it, and the headerless
key=value export PGP kept beside its merged VCFs. Everything else (vendor, array or imputed content,
gVCF blocks, damaged headers) is a note on `vcf`: the parse is the same, the meaning is inferred later.

gVCF is deliberately not a format id. Its header marks lie both ways: `##ALT=<ID=NON_REF>` and
`##GVCFBlock` survive GenotypeGVCFs into files with no reference blocks, and DeepVariant's `.vcf.gz`
and `.g.vcf.gz` carry identical identifying headers. Blocks are a property of the rows, so whether a
callset has them is callset-type inference on the unfiltered rows, not detection.
"""

import polars as pl

from just_vcf.formats.spec import BuildHint, ColumnSpec, Coordinates, Detection, Family, FormatSpec, HeaderMode

# Header statements only; a file whose header names no build stays unknown (contig-length voting is the
# loader's job). `###` covers hand-edited headers with an extra `#`. `b37` is bounded so it cannot match inside ids.
_VCF_BUILD_HINTS: tuple[BuildHint, ...] = (
    BuildHint(
        r"^###?reference=.*(?i:grch37|hg19|hs37d5|human_g1k_v37|g1k\.37|(?<![a-z0-9])b37(?![0-9])|ncbi ?37|build ?37)",
        "GRCh37",
    ),
    BuildHint(r"^###?reference=.*(?i:grch38|hg38|hs38|GCA_000001405\.15|ncbi ?38|build ?38)", "GRCh38"),
    BuildHint(r"^###?reference=.*(?i:grch36|hg18|ncbi ?36|build ?36)", "GRCh36"),
    BuildHint(r"^###?contig=<[^>\n]*assembly=(?i:b37|hg19|grch37)[,>]", "GRCh37"),
    BuildHint(r"^###?contig=<[^>\n]*assembly=(?i:b38|hg38|grch38)[,>]", "GRCh38"),
    BuildHint(r"^###?contig=<[^>\n]*assembly=(?i:b36|hg18|ncbi36|grch36)[,>]", "GRCh36"),
    BuildHint(r"^##source_GENOME_REFERENCE=NCBI build 37$", "GRCh37"),
    BuildHint(r"^##source_GENOME_REFERENCE=NCBI build 36$", "GRCh36"),
    BuildHint(r"^##genomemap=.*(?i:ncbi ?37)", "GRCh37"),
)

_VCF_NO_CALL: tuple[str, ...] = ("./.", ".|.", ".")
_GT_ENCODING = "FORMAT/GT: allele indices separated by / (unphased) or | (phased)"

VCF_CGI = FormatSpec(
    id="vcf_cgi",
    family=Family.VARIANT_CALLS,
    vendor="Complete Genomics",
    description="Complete Genomics `vcfBeta-*.vcf.bz2` (CGAPipeline 2.x through cgatools mkvcf).",
    detection=Detection(
        patterns=(r"^##fileformat=VCFv\d", r"^##center=Complete Genomics$", r"^##source=CGAPipeline_"),
    ),
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding=(
        "FORMAT/GT allele indices; `|` phased within FORMAT/PS, `/` unphased. Half-calls keep one allele and "
        "write `.` for the other (`1|.`, `.|1`, `0|.`, `.|0`, `1/.`, `0/.`). A bare `.` sits on `<CGA_CNVWIN>` "
        "records, which are coverage windows, not genotypes"
    ),
    no_call=_VCF_NO_CALL,
    build_hints=_VCF_BUILD_HINTS,
    notes=(
        "Three record classes share the file. `<CGA_CNVWIN>` ALT rows are 2 kb copy-number windows "
        "(INFO CGA_WINEND, FORMAT CGA_GP/CGA_NP/CGA_CP…, GT `.`): drop them from calls. `<CGA_NOCALL>` rows "
        "are no-called spans with INFO END. Rows with ALT `.`, GT `./.` and AN=0 are no-called reference spans "
        "with no END: their extent is the length of REF, often dozens of bases.",
        "Half-calls are frequent (a few thousand per 400k rows): one allele known, the other `.`. They are "
        "graded evidence, not no-calls and not homozygous.",
        "FILTER was `.` on every record in the first 400k rows of two files; PASS never appears. A loader that "
        "keeps only PASS drops the whole callset.",
        "Chromosomes have no `chr` prefix; `##contig` lines carry `assembly=B37` and an md5.",
        "The header states the build twice: `##source_GENOME_REFERENCE=NCBI build 37` and "
        "`##reference=…/build37.fa.bz2`.",
        "INFO carries CGI annotation (CGA_XR dbSNP ids, CGA_FI functional impact, CGA_RPT, CGA_PFAM); the "
        "`##source_*` lines record the annotation database versions, not tools.",
        "PGP-HMS holds one per CGI genome beside the native var/masterVar files of the same assembly; the "
        "var file is the fuller callset (it states reference-called spans too).",
        "FORMAT GL is declared `Type=Integer` and written in phred units, best genotype 0: in 126,604 of "
        "126,604 rows with GQ in one file, GQ was the gap between the two best. The loader reads gl = GL / 10.",
        "FORMAT AD is declared `Number=2` and counts reads per GT slot (`1/1` writes `3,3`), not per allele: "
        "not carried.",
        "No ALT `.` record is called reference (`0/0`) in the first 400k rows: reference-called positions "
        "are absent. Whether absence means reference here is not settled, so `absence_means` stays `blocks`.",
    ),
    example=(
        "##fileformat=VCFv4.1\n"
        "##fileDate=20140505\n"
        "##center=Complete Genomics\n"
        "##source=CGAPipeline_2.5.0.37\n"
        "##source_GENOME_REFERENCE=NCBI build 37\n"
        "##phasing=partial\n"
        "##reference=ftp://ftp.completegenomics.com/ReferenceFiles/build37.fa.bz2\n"
        '##ALT=<ID=CGA_CNVWIN,Description="Copy number analysis window">\n'
        '##ALT=<ID=CGA_NOCALL,Description="No-called record">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t1\t.\tN\t<CGA_NOCALL>\t.\t.\tEND=10000;NS=1;AN=0\tGT:PS\t./.:.\n"
        "1\t10001\t.\tT\t<CGA_CNVWIN>\t.\t.\tNS=1;CGA_WINEND=12000\tGT\t.\n"
        "1\t752566\trs3094315\tG\tA\t.\t.\tNS=1;AN=2;AC=1\tGT:PS\t0|1:752566\n"
        "1\t752721\trs3131972\tA\tG\t.\t.\tNS=1;AN=1;AC=1\tGT:PS\t1|.:752566\n"
    ),
)

VCF_FREEBAYES_VARSCAN2 = FormatSpec(
    id="vcf_freebayes_varscan2",
    family=Family.VARIANT_CALLS,
    vendor="PGP-HMS",
    description="PGP-HMS's hg19 callsets: freeBayes and VarScan2 run per chromosome and joined by vcf-merge.",
    detection=Detection(
        patterns=(
            r"^##fileformat=VCFv\d",
            r"^##source=freeBayes",
            r"^##source=VarScan2$",
            r"^##source_\d{8}\.\d+=vcf-merge",
        ),
    ),
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding=_GT_ENCODING + "; one genotype per caller, identified by INFO/SF, never by column name",
    no_call=_VCF_NO_CALL,
    build_hints=_VCF_BUILD_HINTS,
    notes=(
        "Two callers, one person. INFO/SF lists the vcf-merge inputs that carry the record, indexing "
        "`##sourceFiles_<date>.1=0:…FB.vcf.gz,1:…VS.vcf.gz`: 0 is freeBayes, 1 is VarScan2. Split the two "
        "facets by SF.",
        "Two layouts. Some files have two sample columns, freeBayes's `unknown` then VarScan2's `Sample1`, with "
        "`.` in the column of a caller that did not call the site. Most have one data column: the header "
        "declares one sample, or declares both (`unknown`, `Sample1`) while every row carries only one. In "
        "the one-column layout the column holds freeBayes's values wherever SF has 0 (VarScan2's genotype on "
        "a shared site is lost) and VarScan2's where SF is 1 alone.",
        "polars-bio 0.36 reads the declared-two, carries-one layout and fills the second sample with None; "
        "that None is a missing column, not a VarScan2 no-call.",
        "INFO AC/AN/AF on shared sites sum the two callers as two samples (AN=4, AC=2 for one het). Never "
        "derive a genotype from INFO.",
        "FORMAT key order varies per record (the union of both callers' keys in hash order); look keys up by "
        "name. freeBayes keys: RO, AO, QR, QA, GL. VarScan2 keys: FREQ (a percentage string), PVAL, RBQ, ABQ, "
        "SDP, RDF, RDR, ADF, ADR.",
        "`##reference=/vbin/ref/hg19.fa` names the build; chromosomes carry `chr`. Rows are sorted as strings "
        "(chr10 before chr1). Some FILTER values are VarScan2's (str10, indelError).",
        "PGP-HMS arvados collections hold one per sequenced genome. A copy stripped of every meta line "
        "detects as `vcf_no_meta`; its rows still carry SF and both callers' keys.",
        "VarScan2 declares AD as Number=1 (the ALT depth only), so AD is not carried; freeBayes's GL is "
        "log10 (Number=G) and is.",
        "The loader returns one facet per caller (SCHEMA Decision 8, `relation = caller_split`): the caller "
        "of an index comes from the `.FB.`/`.VS.` names on `##sourceFiles`; in the one-column layout the "
        "VarScan2 facet holds only the sites VarScan2 alone found, and a warning counts the shared ones.",
    ),
    example=(
        "##fileformat=VCFv4.1\n"
        "##fileDate=20160919\n"
        "##source=freeBayes v0.9.21-24-g840b412-dirty\n"
        "##reference=/vbin/ref/hg19.fa\n"
        "##phasing=none\n"
        "##source=VarScan2\n"
        '##INFO=<ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "##source_20160919.1=vcf-merge(r840) chr1.FB.vcf.gz chr1.VS.vcf.gz\n"
        "##sourceFiles_20160919.1=0:chr1.FB.vcf.gz,1:chr1.VS.vcf.gz\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\tSample1\n"
        "chr1\t752566\t.\tG\tA\t425.3\tPASS\tAC=2;AN=4;SF=0,1\tGT\t0/1\t0/1\n"
        "chr1\t752721\t.\tA\tG\t96.5\t.\tAC=1;AN=2;SF=0\tGT\t1/1\t.\n"
        "chr1\t776546\t.\tA\tG\t.\tPASS\tAC=1;AN=2;SF=1\tGT\t.\t0/1\n"
    ),
)

VCF = FormatSpec(
    id="vcf",
    family=Family.VARIANT_CALLS,
    vendor="generic",
    description="VCF with a `##fileformat=VCF` line, first or not.",
    detection=Detection(patterns=(r"\A#", r"^##fileformat=VCFv\d")),
    comment_prefix="#",
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding=_GT_ENCODING,
    no_call=_VCF_NO_CALL,
    build_hints=_VCF_BUILD_HINTS,
    notes=(
        "polars-bio 0.36 refuses a header whose first line is not `##fileformat` ('missing fileformat'). "
        "Genos exome exports sort their meta lines, so `##fileformat` can sit hundreds of lines down: move it "
        "to line 1 before reading.",
        "polars-bio 0.36 silently drops INFO and FORMAT keys the header does not define, GT included. A "
        "repaired or synthesised header must define every key the loader needs.",
        "Damaged headers that polars-bio rejects as 'invalid record' until repaired: a spreadsheet round-trip "
        "(every meta line tab-padded to the row width, lines containing quotes CSV-quoted as "
        '`"##…=""…"""`); SV and CNV files (uploads named as Dante Labs WGZ) with an extra `#` on meta '
        "lines (`###FORMAT`, `###contig`) and `##CHROM` for the column line. Their CNV rows write ALT `CNV` "
        "without angle brackets.",
        "gVCF blocks are not a format: infer them from unfiltered rows. Observed block spellings: ALT "
        "`<NON_REF>` (GATK, Helix), `<*>` (DeepVariant), ALT `.` with INFO END (Illumina Isaac/starling with "
        "BLOCKAVG_min30p3a, Sequencing.com with MinDP, Genos exome). False friends: `##ALT=<ID=NON_REF>` and "
        "`##GVCFBlock` headers left on genotyped VCFs (GATK pipelines, and uploads named as Dante Labs and "
        "Nebula), END on SV/CNV records, and END=POS on every hom-ref row of the Open Humans 23andMe "
        "conversion.",
        "Genos exome VCFs mark untargeted spans (FILTER and INFO `NT`) and no-called spans (FILTER `NC`) as "
        "END blocks with GT `./.`; targeted positions between blocks are absent. FORMAT adds CGIANN_* "
        "annotation keys.",
        "Array and imputed content also arrives as VCF, where an absent site is unknown, never hom-ref: "
        "Illumina GenotypingToVCF (`##workflow_type=Illumina_GenotypingToVCF`; unplaceable probes as chrom "
        "`NA`, POS 0, REF `N`, ALT both alleles, FILTER GTEX; FILTER NOCALL rows carry GT `./.`), 23andMe "
        "conversions (Open Humans, DNA.Land `intersection_dbSNP_23andme`), Genes for Good (ids such as "
        "`1:762320C>T` off dbSNP), and imputation (Gencove loimpute and Minimac3 with GP/DS/GL, FILTER "
        "LOWCONF).",
        "A header may declare more sample columns than rows carry; polars-bio 0.36 fills the missing samples "
        "with None.",
        "`##source=aengine` (`##ae_version`; uploads named as FamilyTreeDNA Big Y) writes "
        "`DP=…;GTL=…;QUAL=…` into FILTER on failing rows.",
        "The 23andMe exome VCF names no source and no reference; only its `##contig` lines "
        "(`assembly=hg19`) state the build. Its ALT `.` rows are called hom-ref.",
        "Sequencing.com states `##missingDataClarification=not_sequenced`: coordinates absent from the file "
        "were not sequenced, not reference.",
        "DeepVariant headers name no reference; the build comes from contig lengths, which is the loader's job.",
        "Header `##reference` paths sometimes name no build (`refData/1055/1055.fa`, `genome.fa`); the build "
        "then stays unknown here.",
    ),
    example=(
        "##fileformat=VCFv4.2\n"
        "##reference=file:///ref/hs37d5.fa\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t752566\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
    ),
)

VCF_NO_META = FormatSpec(
    id="vcf_no_meta",
    family=Family.VARIANT_CALLS,
    vendor="generic",
    description="VCF body under a bare `#CHROM` line, with no `##` meta lines at all.",
    detection=Detection(
        patterns=(r"\A#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(\tFORMAT(\t[^\t\n]+)+)?\r?$",),
        data_shape=r"[^\t#][^\t]*\t\d+\t[^\t]+\t[ACGTNacgtn]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+(\t[^\t]+)*",
    ),
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding=_GT_ENCODING,
    no_call=_VCF_NO_CALL,
    notes=(
        "polars-bio 0.36 rejects the file ('invalid record'). Prepending `##fileformat` makes it readable, but "
        "every INFO and FORMAT key without a definition is then dropped silently, GT included: synthesise a "
        "definition for each key seen in the rows. The loader does this for the keys the frames carry "
        "(GT, GQ, DP, AD, PL, … and INFO END, SVLEN) and lists every other key the first rows use in "
        "`dropped_fields`; the facet's `vcf_version` stays null, since the file states none.",
        "No build statement is possible; the build comes from the data or from a sibling file.",
        "The one PGP-HMS case is a freeBayes+VarScan2 merge stripped of its meta lines (one `Sample1` column, "
        "INFO SF, both callers' FORMAT keys): split it as `vcf_freebayes_varscan2`'s notes describe.",
    ),
    example=(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
        "chr1\t752566\t.\tG\tA\t50.2\t.\tAB=0.5;SF=0\tGT:RO:AO\t0/1:7:5\n"
        "chr1\t776546\t.\tA\tG\t.\tPASS\tADP=28;SF=1\tGT:FREQ\t0/1:28.57%\n"
    ),
)

BCF = FormatSpec(
    id="bcf",
    family=Family.VARIANT_CALLS,
    vendor="generic",
    description="Binary VCF. Recognised by its `BCF\\2` magic inside BGZF, not by a text head.",
    detection=Detection(patterns=(r"\A(?!)",)),
    comment_prefix=None,
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding="FORMAT/GT, binary-encoded",
    notes=(
        "polars-bio 0.36 reads BCF only from a path named `.bcf` (any other name fails with 'empty "
        "input'), so a BCF under another name or inside an archive is read from a copy named `input.bcf`.",
        "The header is the VCF header text inside the binary layout; the loader reads it itself, and the "
        "first record's contig index and 0-based POS to check polars-bio's coordinates.",
        "A BCF wrote no genotype text, so `source_gt` stays null; polars-bio's rendering is not the file's.",
        "No PGP-HMS file on disk is a BCF (2026-10-22); the tests write theirs with pysam.",
    ),
)

_FBVS_KV_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "chromosome", "dtype": pl.String, "description": "chr1-chr22, chrX, chrY; `chr` prefix kept"},
    {"name": "position", "dtype": pl.Int64, "description": "1-based position, the VCF POS of the sibling merge"},
    {"name": "ref", "dtype": pl.String, "description": "reference allele"},
    {"name": "alt", "dtype": pl.String, "description": "alternate alleles, comma-separated"},
    {"name": "fb_genotype", "dtype": pl.String, "description": "`FBGenoType=`: freeBayes GT, NULL when not called"},
    {"name": "fb_ref_allele_count", "dtype": pl.Int64, "description": "`FBRefAlleleCount=`: `.` or NULL when absent"},
    {"name": "fb_ref_allele_quality", "dtype": pl.Int64, "description": "`FBReferenceAlleleQ=`: NULL when absent"},
    {"name": "fb_total_depth", "dtype": pl.Int64, "description": "`FBTotalDepth=`"},
    {"name": "vs_genotype", "dtype": pl.String, "description": "`VSGenoType=`: VarScan2 GT, NULL when not called"},
    {"name": "vs_ref_allele_count", "dtype": pl.Int64, "description": "`VSRefAlleleCount=`: NULL when absent"},
    {"name": "vs_ref_allele_quality", "dtype": pl.Int64, "description": "`VSReferenceAlleleQ=`: NULL when absent"},
    {"name": "vs_total_depth", "dtype": pl.Int64, "description": "`VSTotalDepth=`"},
)

FREEBAYES_VARSCAN2_KEYVALUE = FormatSpec(
    id="freebayes_varscan2_keyvalue",
    family=Family.VARIANT_CALLS,
    vendor="PGP-HMS",
    description=(
        "PGP-HMS's headerless per-site export of a freeBayes+VarScan2 merge: CHROM, POS, REF, ALT, then both "
        "callers' genotype and depth as positional key=value fields, then annotation."
    ),
    detection=Detection(
        data_shape=(
            r"chr(\d{1,2}|X|Y|M)\t\d+\t[ACGTN]+\t[ACGTN]+(,[ACGTN]+)*"
            r"\tFBGenoType=[^\t]*\tFBRefAlleleCount=[^\t]*\tFBReferenceAlleleQ=[^\t]*\tFBTotalDepth=[^\t]*"
            r"\tVSGenoType=[^\t]*\tVSRefAlleleCount=[^\t]*\tVSReferenceAlleleQ=[^\t]*\tVSTotalDepth=[^\t]*"
            r"(\t.*)?"
        ),
    ),
    comment_prefix=None,
    delimiter="\t",
    header=HeaderMode.NONE,
    columns=_FBVS_KV_COLUMNS,
    coordinates=Coordinates.ONE_BASED,
    genotype_encoding=(
        "VCF-style allele indices into REF and ALT (`0/1`, `1/2`), one per caller, each after its `KEY=` prefix"
    ),
    no_call=("NULL",),
    notes=(
        "It is not a VCF body: no ID, QUAL, FILTER, INFO or FORMAT. Columns 5-12 are `KEY=value` in a fixed "
        "order; strip the key before casting. Missing values are `NULL`, and the freeBayes ref count is also "
        "`.`.",
        "Every row carries a variable tail after column 12 (23 to 84 fields in total): EFF_*/VEP_* annotation "
        "with bare tokens (EFF_RESULTS, MODIFIER, NOVEL_VAR, CREATE_NEW_VAR) repeated per transcript, and "
        "population frequencies (KG, *_MAF, EXAC_FRQ, EVS_FRQ). Drop it or keep it as text; `columns` covers "
        "the fixed twelve only.",
        "Rows are unsorted. Coordinates and alleles match the sibling VCF of the same arvados collection "
        "exactly (20,000 of 20,000 sampled rows at offset 0), so POS is 1-based on hg19. The file states no "
        "build; the sibling's `##reference=/vbin/ref/hg19.fa` does.",
        "A NULL genotype means that caller did not report the site, and all four of its fields are NULL. "
        "freeBayes also writes `0/0` on some rows (a few percent), mostly where VarScan2 is NULL.",
    ),
    example=(
        "chr1\t752566\tG\tA\tFBGenoType=0/1\tFBRefAlleleCount=12\tFBReferenceAlleleQ=475\tFBTotalDepth=33"
        "\tVSGenoType=0/1\tVSRefAlleleCount=12\tVSReferenceAlleleQ=39\tVSTotalDepth=31\n"
        "chr1\t752721\tA\tG\tFBGenoType=1/1\tFBRefAlleleCount=.\tFBReferenceAlleleQ=NULL\tFBTotalDepth=27"
        "\tVSGenoType=1/1\tVSRefAlleleCount=NULL\tVSReferenceAlleleQ=NULL\tVSTotalDepth=27\n"
        "chr1\t846808\tC\tT\tFBGenoType=0/1\tFBRefAlleleCount=23\tFBReferenceAlleleQ=869\tFBTotalDepth=46"
        "\tVSGenoType=NULL\tVSRefAlleleCount=NULL\tVSReferenceAlleleQ=NULL\tVSTotalDepth=NULL\n"
    ),
)

SPECS: tuple[FormatSpec, ...] = (VCF_CGI, VCF_FREEBAYES_VARSCAN2, VCF, VCF_NO_META, BCF)
FALLBACKS: tuple[FormatSpec, ...] = (FREEBAYES_VARSCAN2_KEYVALUE,)
