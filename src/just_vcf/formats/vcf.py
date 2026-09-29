"""VCF, gVCF and BCF, and the vendor files that are VCF underneath.

polars-bio reads these; the specs carry detection and the notes a loader needs.
"""

from just_vcf.formats.spec import BuildHint, Detection, Family, FormatSpec, HeaderMode

_VCF_BUILD_HINTS = (
    BuildHint(r"^##reference=.*(GRCh37|hg19|hs37d5|human_g1k_v37)", "GRCh37"),
    BuildHint(r"^##reference=.*(GRCh38|hg38)", "GRCh38"),
)

VCF = FormatSpec(
    id="vcf",
    family=Family.VARIANT_CALLS,
    vendor="generic",
    description="VCF with its `##fileformat=VCF` line.",
    detection=Detection(patterns=(r"\A##fileformat=VCFv\d",)),
    comment_prefix="#",
    delimiter=None,
    header=HeaderMode.VCF,
    genotype_encoding="FORMAT/GT: allele indices separated by / (unphased) or | (phased)",
    no_call=("./.", ".|.", "."),
    build_hints=_VCF_BUILD_HINTS,
    example=(
        "##fileformat=VCFv4.2\n"
        "##reference=file:///ref/hs37d5.fa\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t752566\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
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
)

SPECS: tuple[FormatSpec, ...] = (VCF, BCF)
FALLBACKS: tuple[FormatSpec, ...] = ()

