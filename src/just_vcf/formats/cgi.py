"""Complete Genomics (CGI) native files: the `var` callset, `masterVarBeta` and the `summary` report.

Every CGI file opens with `#KEY\\tvalue` metadata lines (`#FORMAT_VERSION`, `#GENERATED_BY`, `#TYPE`,
`#GENOME_REFERENCE`…), in no fixed order, then a blank line, then one column-header line that starts
with `>`. `#TYPE` says what the file is; the specs below are told apart by it and by their `>` line.
CGI's own `vcfBeta` export is VCF and belongs to the vcf family, which runs before this one.

Seen on disk (PGP-HMS, 2026-09-29): `#TYPE` VAR-ANNOTATION (`var`, format versions 1.6, 1.7, 2.0, 2.2,
2.5), VAR-OLPL (`masterVarBeta`, 2.0 and 2.5) and SUMMARY-REPORT (`summary`, 2.0 and 2.5). Columns are
declared in their format 2.5 layout; the notes say how earlier versions differ.
"""

import polars as pl

from just_vcf.formats.spec import (
    BuildHint,
    ColumnSpec,
    Coordinates,
    Detection,
    Family,
    FormatSpec,
    HeaderMode,
)

_CGI_BUILD_HINTS = (
    BuildHint(r"^#GENOME_REFERENCE\tNCBI build 37\s*$", "GRCh37"),
    BuildHint(r"^#GENOME_REFERENCE\tNCBI build 36\s*$", "GRCh36"),
)

_CGI_HEADER_NOTES = (
    "The column names are on the first non-comment, non-blank line and start with `>` (`>locus`): strip "
    "the `>` from the first name. A blank line separates the `#` metadata from that line.",
    "Metadata keys come in no fixed order (masterVar sorts them; var and summary do not), so "
    "`#FORMAT_VERSION` is not necessarily the first line. `#GENERATED_BY` appears twice in var files "
    "(cgatools, then dbsnptool).",
    "`#SAMPLE`, `#ASSEMBLY_ID` and `#CUSTOMER_SAMPLE_ID` are identifiers; `#CUSTOMER_SAMPLE_ID` can "
    "carry the participant's id, sometimes with a leading space.",
    "Every file on disk states `#GENOME_REFERENCE\\tNCBI build 37`. The build 36 hint is the same "
    "statement for older deliveries and is not in the PGP-HMS census.",
)

_CGI_COORDINATE_NOTES = (
    "`begin`/`end` are 0-based half-open. For snp, sub and del the 1-based VCF POS is `begin + 1` and "
    "the span is `end - begin` bases. An `ins` has `begin == end` (zero length): VCF needs the anchor "
    "base at 1-based position `begin`, read from the reference FASTA.",
    "`reference` is `=` on ref and no-call blocks (not spelled out) and empty on insertions; REF comes "
    "from the reference FASTA, never from an allele column.",
    "`chromosome` carries the `chr` prefix; map `chrM` to `MT` explicitly, never with a plain `^chr` "
    "strip.",
)

_CGI_ALLELE_SPELLINGS = (
    "Allele sequences: `=` means identical to the reference over the locus, `?` an unknown sequence "
    "(no call), empty a deletion of the whole span; a `?` inside a sequence (`GC?`, `?GAA?`) marks "
    "uncalled bases in a partly called allele."
)

_CGI_VARTYPES = (
    "varType seen: ref, no-ref (reference is N, assembly gap), no-call, snp, ins, del, sub, and in format "
    "1.x also no-call-rc (a no-call sub-kind); masterVar adds complex (alleles of different kinds, and "
    "every no-call locus outside no-ref gaps)."
)


def _tsv(*fields: str, width: int) -> str:
    """One example data line: `fields` padded with empty fields to `width` columns."""
    if len(fields) > width:
        raise ValueError(f"{len(fields)} fields for {width} columns")
    return "\t".join((*fields, *([""] * (width - len(fields))))) + "\n"


# ---------------------------------------------------------------------------------------------------
# var: the callset. One row per haplotype at variant loci, one row per block elsewhere.

CGI_VAR_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "locus", "dtype": pl.Int64, "description": "locus id; the rows of one locus share it"},
    {
        "name": "ploidy",
        "dtype": pl.String,
        "description": "ploidy of the locus (2; 1 expected on male X/Y); only digits seen in heads, which "
        "cover chr1 only, so read as text and cast by the loader",
    },
    {"name": "allele", "dtype": pl.String, "description": "`all` when every haplotype agrees, else `1`, `2`"},
    {"name": "chromosome", "dtype": pl.String, "description": "`chr1`…`chrX`, `chrY`, `chrM`"},
    {"name": "begin", "dtype": pl.Int64, "description": "0-based start"},
    {"name": "end", "dtype": pl.Int64, "description": "0-based exclusive end"},
    {"name": "varType", "dtype": pl.String, "description": "ref, no-ref, no-call, snp, ins, del, sub (see notes)"},
    {"name": "reference", "dtype": pl.String, "description": "reference sequence of the span, `=` when not spelled"},
    {"name": "alleleSeq", "dtype": pl.String, "description": "this haplotype's sequence; `=`, `?`, empty (see notes)"},
    {"name": "varScoreVAF", "dtype": pl.Int64, "description": "variant allele score, VAF model; empty on blocks"},
    {"name": "varScoreEAF", "dtype": pl.Int64, "description": "variant allele score, EAF model; empty on blocks"},
    {
        "name": "varFilter",
        "dtype": pl.String,
        "description": "PASS, VQLOW, AMBIGUOUS, `;`-joined; format ≤2.2 wrote VQHIGH/VQLOW as varQuality",
        "aliases": ("varQuality",),
    },
    {"name": "hapLink", "dtype": pl.String, "description": "phasing link: an integer (2.0) or `Phased_a_b_c` (2.5)"},
    {"name": "xRef", "dtype": pl.String, "description": "`dbsnp.<build>:rs<id>`, `;`-joined"},
    {"name": "alleleFreq", "dtype": pl.String, "description": "`dbsnp:<freq>`, `;`-aligned with xRef, slots may be empty"},
    {"name": "alternativeCalls", "dtype": pl.String, "description": "`<seq>:<score>`, `;`-joined"},
)

CGI_VAR = FormatSpec(
    id="cgi_var",
    family=Family.CALLSET_NATIVE,
    vendor="Complete Genomics",
    description="Complete Genomics `var` file (`#TYPE VAR-ANNOTATION`), format versions 1.6 to 2.5.",
    detection=Detection(
        patterns=(
            r"^#TYPE\tVAR-ANNOTATION\s*$",
            r"^>locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\t",
        ),
        data_shape=r"\d+\t[^\t]*\t(all|\d+)\t[^\t]+\t\d+\t\d+\t[a-z-]+(\t[^\t]*)*",
    ),
    comment_prefix="#",
    delimiter="\t",
    header=HeaderMode.LINE,
    header_prefix=">",
    columns=CGI_VAR_COLUMNS,
    coordinates=Coordinates.ZERO_BASED_HALF_OPEN,
    genotype_encoding=(
        "One row per haplotype at a called or partly called locus: `allele` 1 and 2 (sharing `locus`) "
        "each carry that haplotype's varType and alleleSeq, so the genotype is the pair of rows (a het "
        "SNP is `1 snp G A` + `2 ref G G`). Uniform blocks (ref, no-call, no-ref) are one row with "
        "`allele` `all`. The file has no zygosity column. " + _CGI_ALLELE_SPELLINGS + " " + _CGI_VARTYPES
    ),
    no_call=("?",),
    build_hints=_CGI_BUILD_HINTS,
    notes=(
        *_CGI_HEADER_NOTES,
        *_CGI_COORDINATE_NOTES,
        "The file tiles the whole genome: ref blocks state reference, no-call and no-ref blocks state "
        "uncalled. Absent-means-reference holds only inside ref blocks, and is graded evidence there. "
        "A haplotype row with varType no-call beside a called one is a half call, not a hom call.",
        "Format 1.6/1.7 (12 columns) has `totalScore` where 2.x has varScoreVAF + varScoreEAF + varQuality; "
        "`totalScore` is one score and is not either of the 2.x scores.",
        "Format 2.0/2.2 (14 columns) has `varQuality` (VQHIGH/VQLOW) where 2.5 has `varFilter` "
        "(PASS/VQLOW/AMBIGUOUS, `;`-joined); the column is aliased, the vocabulary is not.",
        "`alleleFreq` and `alternativeCalls` exist from format 2.5. One 2.5 file adds `intervalId` and "
        "`intervalAllele` after them. Read the column names from the `>` line, never by position.",
        "A varType snp can span several reference bases on AMBIGUOUS calls (`TGCC` → `CGCC`).",
        "Some GET-Evidence copies are named `<assembly>.tsv.bz2` without the `var-` prefix.",
    ),
    # A head, not a tiling: loci 3 and 4 sit at rs3094315 (chr1:752566 on GRCh37, reference G), far
    # from locus 2. Alleles are synthetic.
    example=(
        "#GENERATED_BY\tcgatools\n"
        "#SOFTWARE_VERSION\t2.5.0.24\n"
        "#FORMAT_VERSION\t2.5\n"
        "#GENERATED_BY\tdbsnptool\n"
        "#DBSNP_BUILD\tdbSNP build 137\n"
        "#GENOME_REFERENCE\tNCBI build 37\n"
        "#TYPE\tVAR-ANNOTATION\n"
        "\n"
        + _tsv(*(">" + CGI_VAR_COLUMNS[0]["name"], *(c["name"] for c in CGI_VAR_COLUMNS[1:])), width=16)
        + _tsv("1", "2", "all", "chr1", "0", "10000", "no-ref", "=", "?", width=16)
        + _tsv("2", "2", "all", "chr1", "10000", "10468", "no-call", "=", "?", width=16)
        + _tsv("3", "2", "1", "chr1", "752565", "752566", "snp", "G", "A", "120", "120", "PASS", width=16)
        + _tsv("3", "2", "2", "chr1", "752565", "752566", "ref", "G", "G", "98", "98", "PASS", width=16)
        + _tsv("4", "2", "all", "chr1", "752566", "752600", "ref", "=", "=", width=16)
    ),
)


# ---------------------------------------------------------------------------------------------------
# masterVarBeta: one row per locus with both alleles side by side, plus annotation and coverage.

CGI_MASTER_VAR_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "locus", "dtype": pl.Int64, "description": "locus id, one row each"},
    {"name": "ploidy", "dtype": pl.String, "description": "ploidy of the locus; text, as in var"},
    {"name": "chromosome", "dtype": pl.String, "description": "`chr1`…`chrX`, `chrY`, `chrM`"},
    {"name": "begin", "dtype": pl.Int64, "description": "0-based start"},
    {"name": "end", "dtype": pl.Int64, "description": "0-based exclusive end"},
    {"name": "zygosity", "dtype": pl.String, "description": "hom, het-ref, het-alt, half, no-call"},
    {"name": "varType", "dtype": pl.String, "description": "locus-level: ref, no-ref, snp, ins, del, sub, complex"},
    {"name": "reference", "dtype": pl.String, "description": "reference sequence of the span, `=` when not spelled"},
    {"name": "allele1Seq", "dtype": pl.String, "description": "first allele; `=`, `?`, empty (see notes)"},
    {"name": "allele2Seq", "dtype": pl.String, "description": "second allele; `?` on half calls"},
    {"name": "allele1VarScoreVAF", "dtype": pl.Int64, "description": "allele 1 score, VAF model"},
    {"name": "allele2VarScoreVAF", "dtype": pl.Int64, "description": "allele 2 score, VAF model"},
    {"name": "allele1VarScoreEAF", "dtype": pl.Int64, "description": "allele 1 score, EAF model"},
    {"name": "allele2VarScoreEAF", "dtype": pl.Int64, "description": "allele 2 score, EAF model"},
    {
        "name": "allele1VarFilter",
        "dtype": pl.String,
        "description": "PASS, VQLOW, AMBIGUOUS, `;`-joined",
        "aliases": ("allele1VarQuality",),
    },
    {
        "name": "allele2VarFilter",
        "dtype": pl.String,
        "description": "PASS, VQLOW, AMBIGUOUS, `;`-joined",
        "aliases": ("allele2VarQuality",),
    },
    {"name": "allele1HapLink", "dtype": pl.String, "description": "phasing link: integer or `Phased_a_b_c`"},
    {"name": "allele2HapLink", "dtype": pl.String, "description": "phasing link: integer or `Phased_a_b_c`"},
    {"name": "allele1XRef", "dtype": pl.String, "description": "`dbsnp.<build>:rs<id>`, `;`-joined"},
    {"name": "allele2XRef", "dtype": pl.String, "description": "`dbsnp.<build>:rs<id>`, `;`-joined"},
    {"name": "allele1Freq", "dtype": pl.String, "description": "`dbsnp:<freq>`, `;`-aligned with allele1XRef"},
    {"name": "allele2Freq", "dtype": pl.String, "description": "`dbsnp:<freq>`, `;`-aligned with allele2XRef"},
    {"name": "allele1AlternativeCalls", "dtype": pl.String, "description": "`<seq>:<score>`, `;`-joined"},
    {"name": "allele2AlternativeCalls", "dtype": pl.String, "description": "`<seq>:<score>`, `;`-joined"},
    {"name": "evidenceIntervalId", "dtype": pl.Int64, "description": "id of the evidence interval behind the call"},
    {"name": "allele1ReadCount", "dtype": pl.Int64, "description": "reads supporting allele 1"},
    {"name": "allele2ReadCount", "dtype": pl.Int64, "description": "reads supporting allele 2"},
    {"name": "referenceAlleleReadCount", "dtype": pl.Int64, "description": "reads supporting the reference"},
    {"name": "totalReadCount", "dtype": pl.Int64, "description": "reads over the locus"},
    {"name": "allele1Gene", "dtype": pl.String, "description": "`geneId:transcript:symbol:component:impact`, `;`-joined"},
    {"name": "allele2Gene", "dtype": pl.String, "description": "`geneId:transcript:symbol:component:impact`, `;`-joined"},
    {"name": "pfam", "dtype": pl.String, "description": "`PFAM:<id>:<name>`, `;`-joined"},
    {"name": "miRBaseId", "dtype": pl.String, "description": "`<name>:<accession>`, `;`-joined"},
    {"name": "repeatMasker", "dtype": pl.String, "description": "`<repeat>:<family>:<divergence>`, `;`-joined"},
    {"name": "segDupOverlap", "dtype": pl.Int64, "description": "segmental duplications overlapping the locus"},
    {"name": "relativeCoverageDiploid", "dtype": pl.Float64, "description": "coverage relative to the diploid baseline"},
    {"name": "calledPloidy", "dtype": pl.String, "description": "CNV-called ploidy: an integer or `N`"},
    {"name": "relativeCoverageNondiploid", "dtype": pl.Float64, "description": "coverage, nondiploid model"},
    {"name": "calledLevel", "dtype": pl.String, "description": "CNV level: a decimal (2.0) or `N`"},
    {"name": "bestLAFsingle", "dtype": pl.Float64, "description": "lesser allele fraction, best estimate"},
    {"name": "lowLAFsingle", "dtype": pl.Float64, "description": "lesser allele fraction, lower bound"},
    {"name": "highLAFsingle", "dtype": pl.Float64, "description": "lesser allele fraction, upper bound"},
)

CGI_MASTER_VAR = FormatSpec(
    id="cgi_master_var",
    family=Family.CALLSET_NATIVE,
    vendor="Complete Genomics",
    description="Complete Genomics `masterVarBeta` file (`#TYPE VAR-OLPL`), format versions 2.0 and 2.5.",
    detection=Detection(
        patterns=(
            r"^#TYPE\tVAR-OLPL\s*$",
            r"^>locus\tploidy\tchromosome\tbegin\tend\tzygosity\tvarType\treference\tallele1Seq\tallele2Seq\t",
        ),
        data_shape=r"\d+\t[^\t]*\t[^\t]+\t\d+\t\d+\t[a-z-]+\t[a-z-]+(\t[^\t]*)*",
    ),
    comment_prefix="#",
    delimiter="\t",
    header=HeaderMode.LINE,
    header_prefix=">",
    columns=CGI_MASTER_VAR_COLUMNS,
    coordinates=Coordinates.ZERO_BASED_HALF_OPEN,
    genotype_encoding=(
        "One row per locus with both alleles: allele1Seq and allele2Seq, and zygosity hom, het-ref "
        "(one allele is the reference), het-alt (two different non-reference alleles), half (one allele "
        "called, the other `?`) or no-call. varType is locus-level. " + _CGI_ALLELE_SPELLINGS + " "
        + _CGI_VARTYPES
    ),
    no_call=("?",),
    build_hints=_CGI_BUILD_HINTS,
    notes=(
        *_CGI_HEADER_NOTES,
        *_CGI_COORDINATE_NOTES,
        "Like var, the file tiles the genome; a no-call locus carries varType complex (no-ref over "
        "assembly gaps), never no-call. "
        "zygosity half is a partial call, not a hom call and not a no-call.",
        "Format 2.0 (35 columns) has allele1/2VarQuality (VQHIGH/VQLOW) where 2.5 has allele1/2VarFilter "
        "(aliased, vocabulary differs), and lacks allele1/2Freq, allele1/2AlternativeCalls and the three "
        "LAF columns. Read the column names from the `>` line, never by position.",
        "The coverage columns repeat per CNV window and are filled on no-call and no-ref loci too; "
        "calledPloidy and calledLevel mix numbers with `N` (not called).",
        "It restates the calls of the var file of the same assembly (`#ASSEMBLY_ID`); loading both as "
        "separate callsets double-counts them.",
    ),
    # A head, not a tiling: locus 3 is rs3094315 (chr1:752566 on GRCh37, reference G). Alleles,
    # scores, counts and coverage values are synthetic.
    example=(
        "#CNV_DIPLOID_WINDOW_WIDTH\t2000\n"
        "#CNV_NONDIPLOID_WINDOW_WIDTH\t100000\n"
        "#DBSNP_BUILD\tdbSNP build 137\n"
        "#FORMAT_VERSION\t2.5\n"
        "#GENERATED_BY\tcgatools\n"
        "#GENOME_REFERENCE\tNCBI build 37\n"
        "#SOFTWARE_VERSION\t2.5.0.44\n"
        "#TYPE\tVAR-OLPL\n"
        "\n"
        + _tsv(*(">" + CGI_MASTER_VAR_COLUMNS[0]["name"], *(c["name"] for c in CGI_MASTER_VAR_COLUMNS[1:])), width=42)
        + _tsv("1", "2", "chr1", "0", "10000", "no-call", "no-ref", "=", "?", "?", width=42)
        + _tsv(
            "2", "2", "chr1", "10000", "10468", "no-call", "complex", "=", "?", "?",
            *([""] * 25), "0.95", "N", "1.02", "N", "0.31", "0.27", "0.35",
            width=42,
        )
        + _tsv(
            "3", "2", "chr1", "752565", "752566", "het-ref", "snp", "G", "A", "G",
            "120", "98", "120", "98", "PASS", "PASS", "", "", "", "", "", "", "", "",
            "4", "14", "12", "12", "26", "", "", "", "", "", "",
            "1.01", "2", "0.99", "N", "0.48", "0.44", "0.52",
            width=42,
        )
    ),
)


# ---------------------------------------------------------------------------------------------------
# summary: the per-assembly report. Describes a callset, holds no calls.

CGI_SUMMARY_COLUMNS: tuple[ColumnSpec, ...] = (
    {"name": "Category", "dtype": pl.String, "description": "`Genome coverage`, `Genome variations`, `CNV`, `SV`…"},
    {"name": "Metric", "dtype": pl.String, "description": "metric name, e.g. `Fully called genome fraction`"},
    {"name": "Value", "dtype": pl.String, "description": "a count, a fraction or a word (`MALE`); text"},
    {"name": "CallConfidence", "dtype": pl.String, "description": "HIGH, ALL, or NA"},
)

CGI_SUMMARY = FormatSpec(
    id="cgi_summary",
    family=Family.METADATA,
    vendor="Complete Genomics",
    description="Complete Genomics `summary` report (`#TYPE SUMMARY-REPORT`), format versions 2.0 and 2.5.",
    detection=Detection(
        patterns=(
            r"^#TYPE\tSUMMARY-REPORT\s*$",
            r"^>Category\tMetric\tValue\tCallConfidence\s*$",
        ),
    ),
    comment_prefix="#",
    delimiter="\t",
    header=HeaderMode.LINE,
    header_prefix=">",
    columns=CGI_SUMMARY_COLUMNS,
    coordinates=None,
    build_hints=_CGI_BUILD_HINTS,
    notes=(
        *_CGI_HEADER_NOTES,
        "A metric often appears twice, once per CallConfidence (HIGH and ALL); key on "
        "(Category, Metric, CallConfidence).",
        "`Miscellaneous\\tGender` is the called sex. Coverage fractions describe how much of the "
        "genome the matching var file calls, which grades its absent-means-reference evidence.",
        "Written uncompressed on disk, unlike var and masterVar (bz2).",
    ),
    # Values are synthetic.
    example=(
        "#SOFTWARE_VERSION\t2.5.0.24\n"
        "#GENERATED_BY\tExportWorkflow&ReportEngine\n"
        "#FORMAT_VERSION\t2.5\n"
        "#GENOME_REFERENCE\tNCBI build 37\n"
        "#TYPE\tSUMMARY-REPORT\n"
        "#DBSNP_BUILD\tdbSNP build 137\n"
        "\n"
        ">Category\tMetric\tValue\tCallConfidence\n"
        "Miscellaneous\tGender\tFEMALE\tNA\n"
        "Genome coverage\tFully called genome fraction\t0.960\tHIGH\n"
        "Genome coverage\tFully called genome fraction\t0.960\tALL\n"
        "Genome variations\tSNP total count\t3400000\tALL\n"
    ),
)

# #TYPE and the `>` line are both required and differ between the three, so no spec can capture
# another's head; the order is callset, its companion, then the report.
SPECS: tuple[FormatSpec, ...] = (CGI_VAR, CGI_MASTER_VAR, CGI_SUMMARY)
FALLBACKS: tuple[FormatSpec, ...] = ()
