"""The source detector: vendor, caller, version, stated content and ploidy lines, from VCF headers only.

Each case is one row of RM2's header-evidence table (docs/ROADMAP.md), written the way the survey of
the PGP-HMS headers found it (2026-10-22). Values are header text; nothing is read from a file name.
"""

import collections
import json
from pathlib import Path

import pytest

from just_vcf import load, open_payload, sniff
from just_vcf.loaders import vcf
from just_vcf.source import detect_source, ploidy_statements
from test_load_vcf import vcf_text, write

GATK3_HC = (
    '##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller,Version=3.7-0-gcfedb67,Date="Fri Feb 01 22:42:04 UTC 2019",'
    'Epoch=1549060924000,CommandLineOptions="analysis_type=HaplotypeCaller sample_ploidy=2 genotyping_mode=DISCOVERY">'
)
GATK4_GENOTYPE = (
    '##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs --sample-ploidy 2 --output out.vcf",'
    'Version="4.1.4.1",Date="July 13, 2020 4:33:39 PM HKT">'
)
GATK4_HC = (
    '##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --sample-ploidy 2 --emit-ref-confidence GVCF",'
    'Version="4.1.4.1",Date="July 13, 2020 2:56:16 PM HKT">'
)
GATK_SELECT = '##GATKCommandLine.SelectVariants=<ID=SelectVariants,Version=3.7-0-gcfedb67,Date="x",CommandLineOptions="analysis_type=SelectVariants">'
GATK_VQSR = '##GATKCommandLine=<ID=ApplyVQSR,CommandLine="ApplyVQSR",Version="4.1.4.1",Date="x">'
DRAGEN = (
    '##DRAGENCommandLine=<ID=dragen,Version="SW: 05.021.408.3.4.11, HW: 05.021.408",Date="Wed Apr 01 09:02:15 UTC 2020",'
    'CommandLineOptions="--enable-variant-caller true --sample-sex male">'
)
DRAGEN_HASH = '##DRAGENCommandLine=<ID=HashTableBuild,Version="SW: 01.003.044.3.3.5, HashTableVersion: 7",CommandLineOptions="--ht-build">'


@pytest.mark.parametrize(
    ("meta", "expected"),
    [
        (["##source=freeBayes v0.9.20-4-g42d6a54"], (None, "freeBayes", "0.9.20-4-g42d6a54", None)),
        (["##center=Complete Genomics", "##source=CGAPipeline_2.5.0.37"], ("Complete Genomics", "CGAPipeline", "2.5.0.37", None)),
        (["##source=IsaacVariantCaller", "##source_version=2.0.17", "##content=IsaacVariantCaller small-variant calls"],
         (None, "Isaac Variant Caller", "2.0.17", None)),
        (["##source=starling", "##source_version=2.4.5", "##content=starling small-variant calls"],
         (None, "Isaac Variant Caller", "2.4.5", None)),
        (["##workflow_type=Illumina_WGS_CLIA", "##workflow_version=v2.0.12", "##source=IsaacVariantCaller",
          "##content=starling small-variant calls"], ("Illumina", "Isaac Variant Caller", None, None)),
        (["##source=/opt/CASAVA-VariantCalling-2.12a_gVCF/bin/makeConsensusVcf.pl --chrom=chr18"], (None, "CASAVA", None, None)),
        (["##workflow_type=Illumina_GenotypingToVCF", "##workflow_version=v1.4"], ("Illumina", None, None, "array")),
        ([DRAGEN, DRAGEN_HASH, "##source=DRAGEN_SNV"], (None, "DRAGEN", "05.021.408.3.4.11", None)),
        (['##DRAGENCommandLine=<ID=dragen,Date="Mon Aug 19 11:59:22 UTC 2019",CommandLineOptions="x">',
          '##DRAGENVersion=<ID=dragen,Version="SW: 05.011.281.3.2.8, HW: 05.011.281">'], (None, "DRAGEN", "05.011.281.3.2.8", None)),
        ([DRAGEN, GATK_VQSR, "##source=ApplyVQSR"], (None, "DRAGEN", "05.021.408.3.4.11", None)),
        ([GATK3_HC, GATK_SELECT, "##source=SelectVariants"], (None, "GATK HaplotypeCaller", "3.7-0-gcfedb67", None)),
        ([GATK4_HC, GATK4_GENOTYPE, GATK_VQSR, "##source=HaplotypeCaller", "##source=GenotypeGVCFs"],
         (None, "GATK HaplotypeCaller", "4.1.4.1", None)),
        (['##UnifiedGenotyper="analysis_type=UnifiedGenotyper input_file=[x.bam] sample_ploidy=2"'],
         (None, "GATK UnifiedGenotyper", None, None)),
        (["##DeepVariant_version=1.1.0"], (None, "DeepVariant", "1.1.0", None)),
        (["##source=Platypus_Version_0.8.1"], (None, "Platypus", "0.8.1", None)),
        (["##source=GenerateSVCandidates 1.5.1"], (None, "Manta", "1.5.1", None)),
        (["##source=aengine", "##ae_version=v2.1.3-0-g6c08bf6"], (None, "aengine", "v2.1.3-0-g6c08bf6", None)),
        (["##source=loimpute", "##generatedby=Gencove"], ("Gencove", "loimpute", None, "imputed")),
        (["##source=Minimac3", "##generatedby=Gencove", "##eagleVersion=2.3.5+htslib-1.4"], ("Gencove", "Minimac3", None, "imputed")),
        (["##Genes_for_Good_version=1.2: filtered, imputed, phased genotype calls. For more information, go to: x"],
         ("Genes for Good", None, None, "imputed")),
        (["##Genes_for_Good_version=1.0: filtered, unphased genotype calls. For more information, go to: x"],
         ("Genes for Good", None, None, "array")),
        (["##source=open_humans_data_processing.twenty_three_and_me"], ("23andMe", None, None, "array")),
        (["##source=intersection_dbSNP_23andme"], ("23andMe", None, None, "array")),
        (["##pipeline=helix-v3.0.0"], ("Helix", None, None, None)),
        (["##source=Sequencing.com (30x WGS)", "##dataAnalysisProvider=Sequencing.com", "##bcftools_callVersion=1.11+htslib-1.11",
          "##bcftools_callCommand=call -m --ploidy GRCh38 -Oz", "##bcftools_viewVersion=1.11+htslib-1.11"],
         ("Sequencing.com", "bcftools call", "1.11+htslib-1.11", None)),
        (['##AnnotationEngineVersion="0.5.4"', '##FORMAT=<ID=CGIANN_GENE,Number=1,Type=String,Description="Gene">'],
         ("Genos", None, None, None)),
        (['##AnnotationEngineVersion="0.5.4"'], (None, None, None, None)),
        (["##source=SelectVariants", "##bcftools_concatVersion=1.3.1+htslib-1.3.1", GATK_SELECT], (None, None, None, None)),
        ([], (None, None, None, None)),
    ],
)
def test_header_evidence(meta: list[str], expected: tuple) -> None:
    source = detect_source(meta)
    got = (source.vendor, source.caller, source.caller_version, source.callset_type.value if source.callset_type else None)
    assert got == expected
    assert source.warnings == []


def test_eagle_is_the_stated_phasing_method() -> None:
    assert detect_source(["##source=Minimac3", "##eagleVersion=2.3.5+htslib-1.4"]).phasing_method == "Eagle 2.3.5+htslib-1.4"
    assert detect_source(["##source=Minimac3"]).phasing_method is None


def test_two_callers_conflict_and_both_stay_null() -> None:
    source = detect_source(["##source=freeBayes v0.9.21", "##source=VarScan2", "##DeepVariant_version=1.1.0"])
    assert (source.caller, source.caller_version) == (None, None)
    assert [c.name for c in source.callers] == ["freeBayes", "VarScan2", "DeepVariant"]
    assert len(source.warnings) == 1 and "more than one caller" in source.warnings[0]
    assert all(name in source.warnings[0] for name in ("freeBayes from ##source", "DeepVariant from ##DeepVariant_version"))


def test_a_vcf_merge_names_its_callers_without_a_conflict() -> None:
    source = detect_source(
        [
            "##source=freeBayes v0.9.21-24-g840b412-dirty",
            "##source=VarScan2",
            "##source_20160919.1=vcf-merge(r840) chr1.FB.vcf.gz chr1.VS.vcf.gz",
        ],
        spec_vendor="PGP-HMS",
    )
    assert source.merged and source.warnings == [] and source.caller is None
    assert [(c.name, c.version) for c in source.callers] == [("freeBayes", "0.9.21-24-g840b412-dirty"), ("VarScan2", None)]
    assert source.vendor == "PGP-HMS"


def test_two_versions_of_one_caller_withhold_the_version() -> None:
    source = detect_source([GATK3_HC, GATK4_HC])
    assert (source.caller, source.caller_version) == ("GATK HaplotypeCaller", None)
    assert source.warnings == ["the header gives GATK HaplotypeCaller more than one version (3.7-0-gcfedb67, 4.1.4.1); caller_version left null"]


def test_two_vendors_conflict_without_falling_back_to_the_format() -> None:
    source = detect_source(["##generatedby=Gencove", "##pipeline=helix-v3.0.0"], spec_vendor="Complete Genomics")
    assert source.vendor is None and "more than one vendor" in source.warnings[0]
    assert detect_source([], spec_vendor="Complete Genomics").vendor == "Complete Genomics"
    assert detect_source([], spec_vendor="generic").vendor is None


def test_ploidy_statements_are_verbatim_and_paths_are_withheld() -> None:
    meta = [
        '##INFO=<ID=PLOIDY,Number=1,Type=Integer,Description="ploidy of the site">',
        "##referenceSexKaryotype=XY",
        "##ExpectedSexChromosomeKaryotype=XY",
        "##source_MAX_PLOIDY=10",
        GATK3_HC,
        GATK4_GENOTYPE,
        DRAGEN,
        "##bcftools_callCommand=call -m --ploidy GRCh38 -Oz; Date=Tue Mar 19 03:39:45 2024",
        "##cmdline=/opt/starling --ploidy-region-bed /scratch/run7/ploidy0.bed.gz --min-qscore 3",
        "##fileDate=20200101",
    ]
    assert ploidy_statements(meta) == [
        "##referenceSexKaryotype=XY",
        "##ExpectedSexChromosomeKaryotype=XY",
        "##source_MAX_PLOIDY=10",
        "GATKCommandLine.HaplotypeCaller: sample_ploidy=2",
        "GATKCommandLine GenotypeGVCFs: --sample-ploidy 2",
        "DRAGENCommandLine dragen: --sample-sex male",
        "bcftools_callCommand: --ploidy GRCh38",
        "cmdline: --ploidy-region-bed <path>",
    ]
    assert ploidy_statements(["##DeepVariant_version=1.1.0"]) == []


def test_load_fills_the_facet_from_the_header(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]
    extra = (GATK4_HC, GATK4_GENOTYPE, "##source=HaplotypeCaller")
    facet = load(write(tmp_path, vcf_text(records, extra_meta=extra))).facets.collect().row(0, named=True)
    assert (facet["vendor"], facet["caller"], facet["caller_version"]) == (None, "GATK HaplotypeCaller", "4.1.4.1")
    assert facet["ploidy_statements"] == [
        "GATKCommandLine HaplotypeCaller: --sample-ploidy 2",
        "GATKCommandLine GenotypeGVCFs: --sample-ploidy 2",
    ]
    assert facet["callset_type"] is None

    plain = load(write(tmp_path, vcf_text(records), "plain.vcf")).facets.collect().row(0, named=True)
    assert (plain["vendor"], plain["caller"], plain["caller_version"], plain["ploidy_statements"]) == (None, None, None, [])


def test_load_states_array_content_and_warns_on_a_conflict(tmp_path: Path) -> None:
    records = ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"]
    array = load(write(tmp_path, vcf_text(records, extra_meta=("##workflow_type=Illumina_GenotypingToVCF",)), "array.vcf"))
    facet = array.facets.collect().row(0, named=True)
    assert (facet["vendor"], facet["callset_type"], facet["absence_means"]) == ("Illumina", "array", "unknown")

    both = load(write(tmp_path, vcf_text(records, extra_meta=("##source=freeBayes v1.3.6", "##DeepVariant_version=1.6.0")), "both.vcf"))
    assert both.facets.collect()["caller"][0] is None
    assert sum("more than one caller" in w for w in both.warnings) == 1


def test_rsid_tables_have_no_ploidy_statements(tmp_path: Path) -> None:
    path = write(tmp_path, "# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAG\n", "genome.txt")
    facet = load(path).facets.collect().row(0, named=True)
    assert facet["ploidy_statements"] is None and facet["caller"] is None


def test_pgp_hms_headers(pgp_participants: Path) -> None:
    """Every PGP-HMS VCF header in the manifest: the pipeline files name their caller and nothing
    conflicts, except where the header really names two callers. Aggregates only."""
    manifest = pgp_participants / "manifest.jsonl"
    rows = [json.loads(line) for line in manifest.open()]
    callers: collections.Counter = collections.Counter()
    conflicts = 0
    for row in rows:
        if row["outcome"] not in ("fetched", "already_present") or not row["path"] or row["item"]["source"] != "arvados":
            continue
        path = pgp_participants / Path(row["path"]).relative_to(Path(row["path"]).parents[2])
        if not path.exists() or not any(t in path.name.lower() for t in ("vcf",)):
            continue
        sniffed = sniff(path)
        if sniffed.format not in ("vcf_cgi", "vcf_freebayes_varscan2"):
            continue
        with open_payload(sniffed) as handle:
            header, _ = vcf.read_header(handle)
        source = detect_source(header.meta)
        conflicts += bool(source.warnings)
        callers[(sniffed.format, source.caller, tuple(c.name for c in source.callers))] += 1
    if not callers:
        pytest.skip("no PGP-HMS pipeline VCF on disk")
    assert conflicts == 0
    assert {key for key in callers if key[0] == "vcf_cgi"} == {("vcf_cgi", "CGAPipeline", ("CGAPipeline",))}
    assert {key for key in callers if key[0] == "vcf_freebayes_varscan2"} <= {
        ("vcf_freebayes_varscan2", None, ("freeBayes", "VarScan2")),
    }
