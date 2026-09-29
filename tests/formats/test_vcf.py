"""The vcf family: each dialect detects as itself, near misses do not, and builds come from the header."""

import bz2
import gzip
import os
import re
from collections.abc import Callable
from pathlib import Path
from typing import IO

import pytest
from dotenv import dotenv_values, find_dotenv

from just_vcf.formats import REGISTRY, vcf
from just_vcf.formats.spec import FormatSpec
from just_vcf.sniff import SniffOutcome, Sniffed, sniff

VCF_FAMILY: frozenset[str] = frozenset(spec.id for spec in (*vcf.SPECS, *vcf.FALLBACKS))
TEXT_SPECS: tuple[FormatSpec, ...] = tuple(s for s in (*vcf.SPECS, *vcf.FALLBACKS) if s.example)

COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
ROW_37 = "1\t752566\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
ROW_38 = "chr1\t817186\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"


def _sniff_text(tmp_path: Path, text: str, name: str = "input") -> Sniffed:
    path = tmp_path / name
    path.write_text(text)
    return sniff(path)


def _vcf(meta: str, row: str = ROW_37) -> str:
    return "##fileformat=VCFv4.2\n" + meta + COLUMNS + row


def test_every_family_spec_is_registered_in_family_order() -> None:
    """SPECS precede every FALLBACK in the registry, in the order the module lists them."""
    positions = [REGISTRY.index(spec) for spec in (*vcf.SPECS, *vcf.FALLBACKS)]
    assert positions == sorted(positions)


@pytest.mark.parametrize("spec", TEXT_SPECS, ids=lambda s: s.id)
def test_example_detects_as_itself_through_gzip(spec: FormatSpec, tmp_path: Path) -> None:
    path = tmp_path / "example.gz"
    path.write_bytes(gzip.compress(spec.example.encode()))
    result = sniff(path)
    assert (result.outcome, result.format, result.layers) == (SniffOutcome.DETECTED, spec.id, ("gzip",))


def test_cgi_vcfbeta_is_its_own_dialect_with_the_stated_build(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, vcf.VCF_CGI.example)
    assert result.format == "vcf_cgi"
    assert result.build == "GRCh37"
    assert "vcf" in result.also_matched


def test_cgi_markers_alone_are_not_enough(tmp_path: Path) -> None:
    """The dialect needs both the center and the pipeline line; either alone is a plain VCF."""
    centre_only = _vcf("##center=Complete Genomics\n")
    pipeline_only = _vcf("##source=CGAPipeline_2.5.0.37\n")
    assert _sniff_text(tmp_path, centre_only, "a").format == "vcf"
    assert _sniff_text(tmp_path, pipeline_only, "b").format == "vcf"


def test_freebayes_varscan2_detects_in_both_column_layouts(tmp_path: Path) -> None:
    """Two sample columns, or two declared and one carried: the same dialect, split later by INFO/SF."""
    two_columns = vcf.VCF_FREEBAYES_VARSCAN2.example
    header, _, _ = two_columns.partition("chr1\t")
    one_column = header + "chr1\t752566\t.\tG\tA\t425.3\tPASS\tAC=2;AN=4;SF=0,1\tGT:RO:AO\t0/1:19:14\n"
    one_declared = one_column.replace("\tunknown\tSample1\n", "\tSample1\n")
    for name, text in {"two": two_columns, "declared_two": one_column, "declared_one": one_declared}.items():
        result = _sniff_text(tmp_path, text, name)
        assert (result.format, result.build) == ("vcf_freebayes_varscan2", "GRCh37"), name


def test_one_caller_or_an_unmerged_pair_is_plain_vcf(tmp_path: Path) -> None:
    example = vcf.VCF_FREEBAYES_VARSCAN2.example
    no_varscan = example.replace("##source=VarScan2\n", "")
    no_merge = re.sub(r"##source_\d+\.\d+=vcf-merge[^\n]*\n", "", example)
    assert _sniff_text(tmp_path, no_varscan, "fb").format == "vcf"
    assert _sniff_text(tmp_path, no_merge, "pair").format == "vcf"


def test_fileformat_need_not_be_the_first_line(tmp_path: Path) -> None:
    """Sorted meta lines put ##fileformat far down; it is still a VCF, and still carries its build."""
    text = (
        '##AnnotationEngineVersion="0.5.4"\n'
        '##INFO=<ID=NT,Number=0,Type=Flag,Description="The region is not targeted">\n'
        "##fileformat=VCFv4.2\n"
        "##reference=file:///ref/ucsc.hg19.fasta\n"
        "##source=SelectVariants\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1\t.\tN\t.\t.\tNT\tEND=65509;NT\tGT\t./.\n"
    )
    result = _sniff_text(tmp_path, text)
    assert (result.format, result.build) == ("vcf", "GRCh37")


def test_damaged_headers_still_detect_as_vcf(tmp_path: Path) -> None:
    """Spreadsheet padding and quoting, and an extra `#` on every meta line, are loader repairs, not formats."""
    padding = "\t" * 9
    spreadsheet = (
        f"##fileformat=VCFv4.1{padding}\n"
        f'"##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"{padding}\n'
        + COLUMNS
        + ROW_37
    )
    hashed = (
        "##fileformat=VCFv4.1\n"
        '###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "###contig=<ID=chr1,length=249250621,assembly=hg19>\n"
        "###reference=hg19.fasta\n"
        "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t869201\t.\tC\tCNV\t.\tPASS\tSVTYPE=CNV;CN=1.2314;END=870400\tGT\t0/1\n"
    )
    assert _sniff_text(tmp_path, spreadsheet, "sheet").format == "vcf"
    dante = _sniff_text(tmp_path, hashed, "dante")
    assert (dante.format, dante.build) == ("vcf", "GRCh37")


@pytest.mark.parametrize(
    "alt, info",
    [("<NON_REF>", "END=10056"), ("<*>", "END=10056"), (".", "END=10056;BLOCKAVG_min30p3a")],
    ids=["gatk", "deepvariant", "isaac"],
)
def test_gvcf_is_not_a_format(alt: str, info: str, tmp_path: Path) -> None:
    """Reference blocks are left to callset-type inference; the file is a plain VCF."""
    row = f"chr1\t10001\t.\tT\t{alt}\t.\t.\t{info}\tGT\t0/0\n"
    assert _sniff_text(tmp_path, _vcf("", row)).format == "vcf"


def test_no_meta_needs_a_vcf_shaped_body(tmp_path: Path) -> None:
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
    good = _sniff_text(tmp_path, vcf.VCF_NO_META.example, "good")
    not_a_body = _sniff_text(tmp_path, columns + "chr1\tnot-a-position\t.\tG\tA\t.\t.\t.\tGT\t0/1\n", "bad")
    assert (good.format, good.build) == ("vcf_no_meta", None)
    assert not_a_body.format not in VCF_FAMILY


def test_keyvalue_body_needs_both_callers_fields(tmp_path: Path) -> None:
    good = _sniff_text(tmp_path, vcf.FREEBAYES_VARSCAN2_KEYVALUE.example, "good")
    fb_only = "chr1\t752566\tG\tA\tFBGenoType=0/1\tFBRefAlleleCount=12\tFBReferenceAlleleQ=475\tFBTotalDepth=33\n"
    assert (good.format, good.build) == ("freebayes_varscan2_keyvalue", None)
    assert _sniff_text(tmp_path, fb_only, "fb_only").format not in VCF_FAMILY


def test_keyvalue_body_accepts_the_annotation_tail(tmp_path: Path) -> None:
    """Real rows carry a variable annotation tail after the twelve fixed fields."""
    tail = "\tKG=0.3196\tEFF_RESULTS\tEFF_EFFECT=intron_variant\tMODIFIER\t\tprotein_coding\tNOVEL_VAR\tNULL\n"
    rows = [line + tail[:-1] for line in vcf.FREEBAYES_VARSCAN2_KEYVALUE.example.splitlines()]
    assert _sniff_text(tmp_path, "\n".join(rows) + "\n").format == "freebayes_varscan2_keyvalue"


def test_genos_and_cgi_native_files_are_not_vcf(tmp_path: Path) -> None:
    genos = (
        "##FileFormat=Genos\n"
        "##DataVersion=1.0\n"
        "##ReferenceGenome=UCSC.hg19\n"
        "chr1\t752565\t752566\trs3094315\tG/A\n"
    )
    cgi_var = (
        "#GENERATED_BY\tcgatools\n"
        "#FORMAT_VERSION\t2.5\n"
        "#GENOME_REFERENCE\tNCBI build 37\n"
        "#TYPE\tVAR-ANNOTATION\n"
        "\n"
        ">locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\n"
        "1\t2\tall\tchr1\t0\t10000\tno-ref\t=\t?\n"
    )
    assert _sniff_text(tmp_path, genos, "genos").format not in VCF_FAMILY
    assert _sniff_text(tmp_path, cgi_var, "var").format not in VCF_FAMILY


@pytest.mark.parametrize(
    "meta, build",
    [
        ("##reference=file:///creph/genomes/G1K.37/G1K.37.nt.fasta\n", "GRCh37"),
        ("##reference=b37_g1k\n", "GRCh37"),
        ("##reference=NCBI37\n", "GRCh37"),
        ("##reference=GRCh37p13\n", "GRCh37"),
        ("##reference=file:///references/grch37/reference.bin\n", "GRCh37"),
        ("##reference=file:///references/grch38/reference.bin\n", "GRCh38"),
        ("##reference=ucsc.hg38.fasta\n", "GRCh38"),
        ("##reference=file:///ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna\n", "GRCh38"),
        ("##contig=<ID=1,length=249250621,assembly=B37>\n", "GRCh37"),
        ("##contig=<ID=chr1,length=248956422,assembly=hg38>\n", "GRCh38"),
        ("##source_GENOME_REFERENCE=NCBI build 36\n", "GRCh36"),
        ("##genomemap=HumanOmni25M-8v1-1.NCBI37.map.txt\n", "GRCh37"),
        ("##reference=file:///mnt/data/refData/1055/1055.fa\n", None),
        ('##DeepVariant_version=1.3.0\n##contig=<ID=1,length=248956422>\n', None),
    ],
)
def test_build_comes_only_from_header_statements(meta: str, build: str | None, tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _vcf(meta))
    assert result.format == "vcf"
    assert result.build == build
    assert result.warnings == ()


def test_conflicting_build_statements_withhold_the_build(tmp_path: Path) -> None:
    meta = "##reference=file:///ref/ucsc.hg19.fasta\n##contig=<ID=chr1,length=248956422,assembly=hg38>\n"
    result = _sniff_text(tmp_path, _vcf(meta, ROW_38))
    assert result.format == "vcf"
    assert result.build is None
    assert any("GRCh37" in w and "GRCh38" in w for w in result.warnings)


def test_b37_is_not_matched_inside_an_identifier(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _vcf("##reference=file:///data/ab37c/genome.fa\n"))
    assert result.build is None


# Opt-in: real PGP-HMS files under JUST_VCF_SAMPLES_DIR (the repo .env sets it). Heads only are read.

def _samples_dir() -> Path:
    """The environment wins, an empty value included; otherwise the repo `.env`, read without exporting it."""
    if "JUST_VCF_SAMPLES_DIR" in os.environ:
        value = os.environ["JUST_VCF_SAMPLES_DIR"]
    else:
        value = dotenv_values(find_dotenv()).get("JUST_VCF_SAMPLES_DIR")
    if not value:
        pytest.skip("JUST_VCF_SAMPLES_DIR is not set")
    participants = Path(value) / "pgp-hms" / "participants"
    if not participants.is_dir():
        pytest.skip(f"{participants} is absent")
    return participants


def _meta_lines(path: Path, opener: Callable[..., IO[str]]) -> list[str]:
    """The `##` lines above `#CHROM`, read independently of the sniffer."""
    lines: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            lines.append(line)
    return lines


def test_real_pgp_arvados_vcfs_detect_in_the_vcf_family() -> None:
    files = sorted(_samples_dir().glob("*/arvados-*/*.vcf.gz"))
    if not files:
        pytest.skip("no arvados VCFs downloaded")
    results = {path: sniff(path) for path in files}
    outside = {path.name: (r.outcome.value, r.format) for path, r in results.items() if r.format not in VCF_FAMILY}
    assert outside == {}

    def is_merge(path: Path) -> bool:
        """Both callers named, and the vcf-merge line that joined them."""
        meta = _meta_lines(path, gzip.open)
        needed = (r"##source=freeBayes", r"##source=VarScan2$", r"##source_\d+\.\d+=vcf-merge")
        return all(any(re.match(p, line) for line in meta) for p in needed)

    merged = {path for path in files if is_merge(path)}
    detected = {path for path, r in results.items() if r.format == "vcf_freebayes_varscan2"}
    assert detected == merged
    for path in merged:
        assert results[path].build == "GRCh37", path.name


def test_real_cgi_vcfbeta_files_detect_as_vcf_cgi() -> None:
    files = sorted(_samples_dir().glob("*/arvados-*/vcfBeta-*.vcf.bz2"))
    if not files:
        pytest.skip("no CGI vcfBeta files downloaded")
    results = {path.name: sniff(path) for path in files}
    assert {name: r.format for name, r in results.items() if r.format != "vcf_cgi"} == {}
    stated = {
        path.name
        for path in files
        if any(line.startswith("##source_GENOME_REFERENCE=NCBI build 37") for line in _meta_lines(path, bz2.open))
    }
    assert {name for name, r in results.items() if r.build == "GRCh37"} == stated
