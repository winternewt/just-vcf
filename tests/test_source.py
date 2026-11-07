"""The source detector: vendor, caller, version, stated content and ploidy lines, from VCF headers only.

Each case is one row of RM2's header-evidence table (docs/ROADMAP.md), written the way the survey of
the PGP-HMS headers found it (2026-10-22). Values are header text; nothing is read from a file name.
"""

import collections
import json
import re
from pathlib import Path

import pytest

from just_vcf import load, open_payload, sniff
from just_vcf.loaders import vcf
from just_vcf.schema import STEP_ROLES
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
GATK_SELECT = '##GATKCommandLine.SelectVariants=<ID=SelectVariants,Version=3.7-0-gcfedb67,Date="Fri Feb 01 23:10:00 UTC 2019",CommandLineOptions="analysis_type=SelectVariants">'
GATK_VQSR = '##GATKCommandLine=<ID=ApplyVQSR,CommandLine="ApplyVQSR",Version="4.1.4.1",Date="July 13, 2020 6:01:02 PM HKT">'
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


BEAGLE_INFO = (
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">',
    (
        '##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between '
        'estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">'
    ),
    '##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">',
)
"""The INFO lines Beagle 4.1 and 5.x declare on every run (`vcf/VcfWriter.java`)."""


@pytest.mark.parametrize(
    ("line", "version"),
    [
        # 5.5, 5.4, 5.1 and 5.0: `Main.PROGRAM`, quoted by `VcfWriter.writeMetaLines`.
        ('##source="beagle.27Feb25.75f.jar"', "27Feb25.75f"),
        ('##source="beagle.22Jul22.46e.jar"', "22Jul22.46e"),
        ('##source="beagle.18May20.d20.jar"', "18May20.d20"),
        ('##source="beagle.28Sep18.793.jar"', "28Sep18.793"),
        # 4.1 names its release after the jar; 4.0 names the jar plainly and its revision after it.
        ('##source="beagle.27Jul16.86a.jar (version 4.1)"', "27Jul16.86a"),
        ('##source="beagle.jar (r1399)"', "r1399"),
    ],
)
def test_beagle_is_the_caller_with_its_build(line: str, version: str) -> None:
    for meta in ([line], [line, *BEAGLE_INFO]):
        source = detect_source(meta)
        assert (source.caller, source.caller_version) == ("Beagle", version)
        assert [(s.tool, s.version, s.role, s.line) for s in source.steps] == [("Beagle", version, "call", "##source")]
        # Beagle phases without imputing too: its ##source states no content, with or without DR2 and IMP.
        assert source.callset_type is None
        assert source.warnings == []


def test_a_beagle_vcf_reads_imputed_from_its_info_keys(tmp_path: Path) -> None:
    """The facet's `imputed` comes from Beagle's DR2 and IMP declarations (`callset.imputation_marker`),
    not from its ##source; the caller is Beagle either way."""
    records = ["1\t100\t.\tA\tG\t.\tPASS\tAF=0.1;DR2=0.97\tGT\t0|1"]
    meta = ('##source="beagle.22Jul22.46e.jar"', *BEAGLE_INFO)
    facet = load(write(tmp_path, vcf_text(records, extra_meta=meta), "beagle.vcf")).facets.collect().row(0, named=True)
    assert (facet["caller"], facet["caller_version"], facet["callset_type"], facet["phasing_method"]) == ("Beagle", "22Jul22.46e", "imputed", None)


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


def test_a_caller_version_written_as_a_path_is_withheld() -> None:
    """`##DeepVariant_version=/home/<name>/…` is a path, not a version: `<path>`, as a step's version is."""
    source = detect_source(["##DeepVariant_version=/home/Plantedsmith/bin/deepvariant-1.6.0"])
    assert (source.caller, source.caller_version) == ("DeepVariant", "<path>")
    assert [(c.name, c.version) for c in source.callers] == [("DeepVariant", "<path>")]
    assert source.warnings == [] and "Plantedsmith" not in str([s.as_row() for s in source.steps])


def test_two_caller_versions_written_as_paths_still_conflict() -> None:
    """The versions are compared as written, before either is withheld: two paths are two versions."""
    source = detect_source(["##DeepVariant_version=/home/Plantedsmith/dv", "##DeepVariant_version=C:\\Users\\Plantedjones\\dv"])
    assert (source.caller, source.caller_version) == ("DeepVariant", None)
    assert source.warnings == ["the header gives DeepVariant more than one version (<path>, <path>); caller_version left null"]
    assert source.callers[0].version is None


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
    assert [(s["tool"], s["role"], s["line"]) for s in facet["method_steps"]] == [
        ("GATK HaplotypeCaller", "call", "##GATKCommandLine"),
        ("GATK GenotypeGVCFs", "call", "##GATKCommandLine"),
    ]
    assert facet["method_steps"][0]["options"] == ["--emit-ref-confidence GVCF"]
    assert (facet["reference"], facet["aligner"], facet["aligner_version"]) == (None, None, None)

    plain = load(write(tmp_path, vcf_text(records), "plain.vcf")).facets.collect().row(0, named=True)
    assert (plain["vendor"], plain["caller"], plain["caller_version"], plain["ploidy_statements"]) == (None, None, None, [])
    assert (plain["method_steps"], plain["reference"], plain["phasing_method"]) == ([], None, None)

    dragen = load(write(tmp_path, vcf_text(records, extra_meta=(DRAGEN_MAPPED, DRAGEN_PS)), "dragen.vcf")).facets.collect()
    row = dragen.row(0, named=True)
    assert (row["aligner"], row["aligner_version"], row["reference"], row["phasing_method"]) == (
        "DRAGEN", "07.021.624.3.10.4", "hg38_alt_aware", "DRAGEN 07.021.624.3.10.4"
    )
    assert [s["role"] for s in row["method_steps"]] == ["align", "call", "phase"]


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
    assert (facet["method_steps"], facet["reference"], facet["aligner"]) == (None, None, None)


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
        roles = {(step.role, step.tool) for step in source.steps}
        assert {step.role for step in source.steps} <= STEP_ROLES
        assert not any(re.search(r"[/\\]", re.sub(r"\b(?:INFO|FORMAT)/", "", o)) for step in source.steps for o in step.options)
        if sniffed.format == "vcf_cgi":
            # CGI declares PS in its own words in every pipeline VCF: the caller phased.
            assert ("phase", "CGAPipeline") in roles and source.phasing_method == f"CGAPipeline {source.caller_version}"
        else:
            assert ("merge", "vcf-merge") in roles and {("call", "freeBayes"), ("call", "VarScan2")} <= roles
        assert source.reference is not None
    if not callers:
        pytest.skip("no PGP-HMS pipeline VCF on disk")
    assert conflicts == 0
    assert {key for key in callers if key[0] == "vcf_cgi"} == {("vcf_cgi", "CGAPipeline", ("CGAPipeline",))}
    assert {key for key in callers if key[0] == "vcf_freebayes_varscan2"} <= {
        ("vcf_freebayes_varscan2", None, ("freeBayes", "VarScan2")),
    }


# The method chain (charter 6). Lines are written the way the PGP-HMS and local headers carry them
# (survey of 2026-10-23); the paths in them are invented directories, which every step must withhold.

GATK3_APPLY = (
    '##GATKCommandLine.ApplyRecalibration=<ID=ApplyRecalibration,Version=3.8-1-0-gf15c1c3ef,Date="Tue Mar 02 10:00:00 CST 2021",'
    'Epoch=1614650400000,CommandLineOptions="analysis_type=ApplyRecalibration input_file=[] reference_sequence=/ref/hg38.fa '
    'input=[(RodBinding name=input source=/work/in.vcf)] mode=SNP ts_filter_level=99.5 excludeFiltered=false">'
)
GATK3_GENOTYPE = (
    '##GATKCommandLine.GenotypeGVCFs=<ID=GenotypeGVCFs,Version=3.8-1-0-gf15c1c3ef,Date="Mon Mar 01 10:00:00 CST 2021",'
    'Epoch=1614564000000,CommandLineOptions="analysis_type=GenotypeGVCFs reference_sequence=/ref/hg38.fa out=/work/out.vcf">'
)
GATK3_SELECT = (
    '##GATKCommandLine.SelectVariants=<ID=SelectVariants,Version=3.8-1-0-gf15c1c3ef,Date="Wed Mar 03 10:00:00 CST 2021",'
    'Epoch=1614736800000,CommandLineOptions="analysis_type=SelectVariants sample_name=[NA12878] selectTypeToInclude=[SNP] '
    'selectTypeToExclude=[] excludeFiltered=true">'
)
GATK4_HC_PHASED = (
    '##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --emit-ref-confidence GVCF --output /work/x.g.vcf.gz '
    '--input /work/x.bam --reference /ref/hg38.fa --sample-ploidy 2",Version="4.1.8.1",Date="February 20, 2021 4:33:39 PM CST">'
)
GATK_PS = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">'
NORM_VERSION = "##bcftools_normVersion=1.23.1+htslib-1.23.1"
NORM = (
    "##bcftools_normCommand=norm -m -both -f /ref/hg38.fa --old-rec-tag ORIG -Oz -o /work/norm.vcf.gz /work/dv.vcf; "
    "Date=Tue Mar 19 03:39:45 2024"
)
DRAGEN_MAPPED = (
    '##DRAGENCommandLine=<ID=dragen,Version="SW: 07.021.624.3.10.4, HW: 07.021.624",Date="Wed Apr 01 09:02:15 UTC 2020",'
    'CommandLineOptions="--ref-dir /staging/ref/hg38_alt_aware --fastq-list /work/list.csv --enable-map-align true '
    '--enable-duplicate-marking true --enable-variant-caller true --vc-sample-name NA12878 --RGSM NA12878 --output-directory /work">'
)
DRAGEN_PS = (
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Physical phasing ID information, where each unique ID within a given '
    'sample (but not across samples) connects records within a phasing group">'
)
CLAIR3 = [
    "##clair3_version=1.0.10",
    "##cmdline=/opt/bin/run_clair3.sh --bam_fn=/work/x.bam --ref_fn=/ref/hg38.fa --platform=ilmn --sample_name=NA12878 --output=/work/out",
]


def _chain(meta: list[str]) -> list[tuple[str, str | None, str, tuple[str, ...]]]:
    return [(s.tool, s.version, s.role, s.options) for s in detect_source(meta).steps]


def test_gatk_steps_run_in_date_order_not_in_the_sorted_header_order() -> None:
    """htsjdk writes meta lines sorted by text, so ApplyRecalibration comes first in the header; the
    stated dates put the steps back in the order they ran."""
    meta = [GATK3_APPLY, GATK3_GENOTYPE, GATK3_SELECT, GATK4_HC_PHASED, GATK_PS, "##source=HaplotypeCaller", NORM_VERSION, NORM]
    assert _chain(meta) == [
        ("GATK HaplotypeCaller", "4.1.8.1", "call", ("--emit-ref-confidence GVCF", "--reference hg38.fa")),
        ("GATK HaplotypeCaller", "4.1.8.1", "phase", ()),
        ("GATK GenotypeGVCFs", "3.8-1-0-gf15c1c3ef", "call", ("reference_sequence=hg38.fa",)),
        ("GATK ApplyRecalibration", "3.8-1-0-gf15c1c3ef", "filter", ("mode=SNP", "ts_filter_level=99.5")),
        ("GATK SelectVariants", "3.8-1-0-gf15c1c3ef", "filter", ("selectTypeToInclude=[SNP]", "excludeFiltered=true")),
        ("bcftools norm", "1.23.1+htslib-1.23.1", "normalize", ("-m -both", "-f hg38.fa", "--old-rec-tag ORIG")),
    ]
    # The order is the dates', whatever order the lines come in.
    assert _chain([GATK3_SELECT, GATK4_HC_PHASED, GATK3_APPLY, GATK3_GENOTYPE]) == _chain(
        [GATK3_APPLY, GATK3_GENOTYPE, GATK3_SELECT, GATK4_HC_PHASED]
    )
    source = detect_source(meta)
    assert (source.reference, source.phasing_method, source.aligner) == ("hg38.fa", "GATK HaplotypeCaller 4.1.8.1", None)
    assert source.warnings == []


def test_undated_steps_keep_header_order() -> None:
    meta = [
        "##source=freeBayes v0.9.21-24-g840b412-dirty",
        '##commandline="/opt/bin/freebayes -f /ref/hg19.fa -t /work/targets.bed /work/x.bam"',
        "##source=VarScan2",
        "##reference=/ref/hg19.fa",
        "##source_20160919.1=vcf-merge(r840) /work/chr1.FB.vcf.gz /work/chr1.VS.vcf.gz",
        "##source_20160919.2=/opt/bin/vcf-annotate(r953) -d key=INFO,ID=SF,Number=.,Type=String -a /work/sf.tab.gz -c CHROM,FROM,TO,INFO/SF",
    ]
    assert _chain(meta) == [
        ("freeBayes", "0.9.21-24-g840b412-dirty", "call", ("-f hg19.fa",)),
        ("VarScan2", None, "call", ()),
        ("vcf-merge", "r840", "merge", ()),
        ("vcf-annotate", "r953", "annotate", ("-c CHROM,FROM,TO,INFO/SF",)),
    ]
    assert detect_source(meta).reference == "hg19.fa"


def test_a_dragen_run_that_mapped_the_reads_records_the_aligner() -> None:
    meta = ["##reference=file:///staging/ref/hg38_alt_aware/reference.bin", DRAGEN_MAPPED, DRAGEN_PS, "##source=DRAGEN_SNV"]
    source = detect_source(meta)
    assert _chain(meta) == [
        ("DRAGEN", "07.021.624.3.10.4", "align", ("--ref-dir hg38_alt_aware", "--enable-map-align true", "--enable-duplicate-marking true")),
        ("DRAGEN", "07.021.624.3.10.4", "call", ("--ref-dir hg38_alt_aware", "--enable-variant-caller true")),
        ("DRAGEN", "07.021.624.3.10.4", "phase", ()),
    ]
    assert (source.aligner, source.aligner_version) == ("DRAGEN", "07.021.624.3.10.4")
    assert (source.reference, source.phasing_method, source.warnings) == ("hg38_alt_aware", "DRAGEN 07.021.624.3.10.4", [])
    # From an external BAM DRAGEN maps nothing, so there is no aligner.
    unmapped = detect_source([DRAGEN_MAPPED.replace("--enable-map-align true", "--enable-map-align false")])
    assert unmapped.aligner is None and [s.role for s in unmapped.steps] == ["call"]


def test_options_withhold_paths_samples_inputs_and_outputs() -> None:
    view = (
        "##bcftools_viewCommand=view -s NA12878 -f PASS -R /work/regions.bed -Oz -o /work/out.vcf.gz /work/in.vcf.gz; "
        "Date=Tue Mar 19 04:00:00 2024"
    )
    options = [option for step in detect_source([DRAGEN_MAPPED, GATK3_SELECT, GATK4_HC_PHASED, NORM, view, *CLAIR3]).steps for option in step.options]
    assert not any(word in option for option in options for word in ("/", "\\", "NA12878", "out", "x.bam", "in.vcf", "list.csv"))
    assert "-R <path>" in options and "-f PASS" in options


def test_bcftools_steps_take_their_role_from_the_subcommand_and_its_options() -> None:
    meta = [
        "##bcftools_viewVersion=1.23.1+htslib-1.23.1",
        "##bcftools_viewCommand=view -Oz -o /work/a.vcf.gz /work/b.vcf.gz; Date=Tue Mar 19 03:00:00 2024",
        "##bcftools_viewCommand=view -f PASS,. -Oz -o /work/a.vcf.gz /work/b.vcf.gz; Date=Tue Mar 19 03:10:00 2024",
        "##bcftools_isecVersion=1.23.1+htslib-1.23.1",
        "##bcftools_isecCommand=isec -n=3 -w1 -Oz -o /work/c.vcf.gz /work/a.vcf.gz /work/b.vcf.gz /work/d.vcf.gz; Date=Tue Mar 19 03:20:00 2024",
        "##bcftools_concatCommand=concat -a -f /work/list.txt -Oz; Date=Tue Mar 19 03:30:00 2024",
        "##bcftools_annotateCommand=annotate -a /work/clinvar.vcf.gz -c INFO/CLNSIG -Oz; Date=Tue Mar 19 03:40:00 2024",
    ]
    assert [(s.tool, s.version, s.role, s.options) for s in detect_source(meta).steps] == [
        ("bcftools view", "1.23.1+htslib-1.23.1", "other", ()),
        ("bcftools view", "1.23.1+htslib-1.23.1", "filter", ("-f PASS,.",)),
        ("bcftools isec", "1.23.1+htslib-1.23.1", "merge", ("-n=3", "-w1")),
        ("bcftools concat", None, "merge", ("-a",)),
        ("bcftools annotate", None, "annotate", ("-a <path>", "-c INFO/CLNSIG")),
    ]


def test_a_crossmap_lift_is_a_step_and_explains_a_second_reference() -> None:
    lift = [
        "##liftOverProgram=<CrossMap,version=0.6.3,website=https://sourceforge.net/projects/crossmap>",
        "##liftOverChainFile=/chains/GRCh38_to_GRCh37.chain.gz",
        "##originalFile=/work/someone.hg38.vcf",
        "##targetRefGenome=/ref/Homo_sapiens.GRCh37.75.dna.toplevel.fa",
    ]
    called = [GATK4_HC_PHASED, "##reference=file:///ref/Homo_sapiens.GRCh37.75.dna.toplevel.fa"]
    source = detect_source(called + lift)
    last = source.steps[-1]
    assert (last.tool, last.version, last.role, last.options) == (
        "CrossMap", "0.6.3", "lift", ("liftOverChainFile=GRCh38_to_GRCh37.chain.gz", "targetRefGenome=Homo_sapiens.GRCh37.75.dna.toplevel.fa"),
    )
    assert not any("someone" in option for step in source.steps for option in step.options)
    assert source.reference == "Homo_sapiens.GRCh37.75.dna.toplevel.fa" and source.warnings == []
    # The same pair of references with no lift recorded (Picard LiftoverVcf writes none) is kept, with a warning.
    unexplained = detect_source(called)
    assert unexplained.reference == "Homo_sapiens.GRCh37.75.dna.toplevel.fa"
    assert len(unexplained.warnings) == 1 and "no liftover step is recorded" in unexplained.warnings[0]


def test_a_tool_or_version_written_as_a_path_is_withheld_like_an_option() -> None:
    """A path names a home directory; the frame never holds it, in a step's tool or version either."""
    meta = [
        "##liftOverProgram=</home/plantedname/bin/CrossMap,version=0.6.3>",
        "##bcftools_normVersion=/home/plantedname/opt/bcftools-1.9",
        "##bcftools_normCommand=norm -m -both -Oz -o /work/out.vcf.gz /work/in.vcf.gz; Date=Tue Mar 19 03:39:45 2024",
    ]
    steps = detect_source(meta).steps
    assert [(s.tool, s.version, s.role) for s in steps] == [("<path>", "0.6.3", "lift"), ("bcftools norm", "<path>", "normalize")]
    assert not any("plantedname" in str(value) for step in steps for value in step.as_row().values())


@pytest.mark.parametrize(
    "name",
    ["hg38.fa", "hg19.fa", "ucsc.hg19.fasta", "hs37d5.fa.gz", "hs38DH.fa", "human_g1k_v37_decoy.fasta", "G1K.37.nt.fasta",
     "Homo_sapiens_assembly38.fasta", "Homo_sapiens.GRCh37.75.dna.toplevel.fa", "GRCh38.p13.genome.fa",
     "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna", "GRCh38_full_analysis_set_plus_decoy_hla.fa",
     "GCA_009914755.4_T2T-CHM13v2.0_genomic.fna", "chm13v2.0.fa", "hs1.fa.gz", "b37.fasta", "build37.fa.bz2",
     "hg38_alt_aware", "hg38-alt_masked.cnv.graph.hla.rna-8-r2.0-1", "GRCh38_to_GRCh37.chain.gz", "hg19ToHg38.over.chain.gz"],
)  # fmt: skip
def test_a_known_reference_file_keeps_its_base_name(name: str) -> None:
    source = detect_source([f"##reference=file:///refs/{name}", f'##commandline="/opt/bin/freebayes -f /refs/{name} /work/x.bam"'])
    assert source.reference == name and source.steps[0].options == (f"-f {name}",)


def test_any_other_reference_or_chain_file_name_is_withheld() -> None:
    """A base name can carry a person's name: only a name made of known reference words is kept."""
    meta = [
        "##reference=file:///refs/PLANTEDSMITH-hg38.fa",
        '##commandline="/opt/bin/freebayes -f /refs/hg38_PlantedJones.fa /work/x.bam"',
        "##liftOverProgram=<CrossMap,version=0.6.3>",
        "##liftOverChainFile=/chains/planted_doe_to_hg19.chain.gz",
        "##targetRefGenome=/refs/ref.fa",
        '##DRAGENCommandLine=<ID=dragen,Version="SW: 07.021.624.3.10.4",CommandLineOptions="--ref-dir /staging/PlantedRoe_hg38 --enable-variant-caller true">',
    ]
    source = detect_source(meta)
    row = [step.as_row() for step in source.steps]
    assert source.reference == "<file>"
    assert not any(planted in str(row) for planted in ("PLANTED", "Planted", "planted"))
    assert ("freeBayes", ("-f <file>",)) in [(s.tool, s.options) for s in source.steps]
    assert ("CrossMap", ("liftOverChainFile=<file>", "targetRefGenome=<file>")) in [(s.tool, s.options) for s in source.steps]
    assert ("DRAGEN", ("--ref-dir <file>", "--enable-variant-caller true")) in [(s.tool, s.options) for s in source.steps]
    assert not any(planted in " ".join(source.warnings) for planted in ("PLANTED", "Planted", "planted"))


PLANTED_REFERENCES = ("/refs/PLANTEDSMITH-hg38.fa", "/refs/PLANTEDJONES-hg38.fa")
"""Two reference files whose base names are not a known reference's: both are stored as `<file>`."""


def test_two_different_withheld_reference_files_still_conflict() -> None:
    """The names are compared as written, then withheld: two `<file>`s are not one reference."""
    source = detect_source([f"##reference=file://{path}" for path in PLANTED_REFERENCES])
    assert source.reference is None
    assert source.warnings == ["the header names 2 different ##reference files; reference left null"]
    same = detect_source([f"##reference=file://{PLANTED_REFERENCES[0]}"] * 2)
    assert (same.reference, same.warnings) == ("<file>", [])


def test_two_different_withheld_call_references_still_conflict() -> None:
    meta = [f'##commandline="/opt/bin/freebayes -f {path} /work/x.bam"' for path in PLANTED_REFERENCES]
    source = detect_source(meta)
    assert source.reference is None
    assert source.warnings == ["the call steps name 2 different reference files; reference left null"]


def test_a_withheld_reference_that_differs_from_the_callers_is_still_flagged() -> None:
    """`##reference=` and the caller's `-f` name two different files, both withheld: the warning fires,
    naming neither, and the stored reference is `<file>`."""
    stated, called = PLANTED_REFERENCES
    source = detect_source([f"##reference=file://{stated}", f'##commandline="/opt/bin/freebayes -f {called} /work/x.bam"'])
    assert source.reference == "<file>"
    assert source.warnings == [
        "##reference and the call step name two different reference files, and no liftover step is recorded: the rows "
        "may have been moved to another reference by a tool the header does not name"
    ]
    agree = detect_source([f"##reference=file://{stated}", f'##commandline="/opt/bin/freebayes -f {stated} /work/x.bam"'])
    assert (agree.reference, agree.warnings) == ("<file>", [])


def test_reference_falls_back_to_the_caller_and_two_references_leave_it_null() -> None:
    assert detect_source(CLAIR3).reference == "hg38.fa"
    two = detect_source(["##reference=/ref/hg19.fa", "##reference=/ref/hs37d5.fa"])
    assert two.reference is None and two.warnings == ["the header names 2 different ##reference files; reference left null"]
    assert detect_source(["##reference=GRCh38"]).reference == "GRCh38"
    assert (detect_source([]).reference, detect_source([]).steps) == (None, [])


def test_clair3_is_a_caller_with_its_command_line_options() -> None:
    source = detect_source(CLAIR3)
    assert (source.caller, source.caller_version) == ("Clair3", "1.0.10")
    assert _chain(CLAIR3) == [("Clair3", "1.0.10", "call", ("--ref_fn=hg38.fa", "--platform=ilmn"))]


@pytest.mark.parametrize(
    ("meta", "method"),
    [
        (["##source=Minimac3", "##eagleVersion=2.3.5+htslib-1.4", "##eagleCommand=eagle --vcfTarget /work/x.vcf.gz --chrom 1"],
         "Eagle 2.3.5+htslib-1.4"),
        (["##DeepVariant_version=1.6.0", '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'], "DeepVariant 1.6.0"),
        (["##center=Complete Genomics", "##source=CGAPipeline_2.5.0.37", '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">'],
         "CGAPipeline 2.5.0.37"),
        (['##commandline="(whatshap 1.7) phase -o /work/out.vcf --reference /ref/hg38.fa /work/in.vcf /work/in.bam"'], "WhatsHap 1.7"),
        (["##source=shapeit4.1.3"], "SHAPEIT4 1.3"),
        (['##INFO=<ID=hapcut2,Number=1,Type=Integer,Description="phased by HapCUT2 or not">'], "HapCUT2"),
        # PS in other words than the caller's own: the caller is not shown to have phased.
        (["##DeepVariant_version=1.6.0", '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">'], None),
        (["##DeepVariant_version=1.6.0"], None),
    ],
)
def test_phasing_methods(meta: list[str], method: str | None) -> None:
    source = detect_source(meta)
    assert source.phasing_method == method
    assert source.warnings == []


def test_two_phasers_leave_the_method_null() -> None:
    source = detect_source([
        "##DeepVariant_version=1.6.0",
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        '##commandline="(whatshap 1.7) phase -o /work/out.vcf /work/in.vcf /work/in.bam"',
    ])
    assert source.phasing_method is None
    assert source.warnings == ["the header records more than one phasing step (DeepVariant 1.6.0; WhatsHap 1.7); left null"]
    assert [s.options for s in source.steps if s.tool == "WhatsHap"] == [("phase",)]


def test_gatk2_and_annotation_steps() -> None:
    meta = [
        '##SelectVariants="analysis_type=SelectVariants input_file=[] reference_sequence=/ref/G1K.37.nt.fasta selectTypeToInclude=[INDEL] out=/work/x"',
        '##OriginalSnpEffCmd="SnpEff eff -config /opt/snpEff.config -onlyCoding true -o vcf -i vcf GRCh37.64 /work/x.vcf"',
        '##OriginalSnpEffVersion="2.0.5 (build 2011-12-24), by Pablo Cingolani"',
        "##VEP=v88 cache=/opt/vep/cache db=homo_sapiens_core_88_37 assembly=GRCh37",
    ]
    assert _chain(meta) == [
        ("GATK SelectVariants", None, "filter", ("selectTypeToInclude=[INDEL]",)),
        ("SnpEff", "2.0.5", "annotate", ("-onlyCoding true", "GRCh37.64")),
        ("VEP", "v88", "annotate", ("db=homo_sapiens_core_88_37", "assembly=GRCh37")),
    ]


def test_dates_in_gatk_4_2_form_parse() -> None:
    later = GATK4_HC_PHASED.replace("February 20, 2021 4:33:39 PM CST", "August 13, 2022 at 9:01:02 PM GMT").replace(
        "HaplotypeCaller", "SelectVariants"
    )
    source = detect_source([later, GATK4_HC_PHASED])
    assert [s.tool for s in source.steps] == ["GATK HaplotypeCaller", "GATK SelectVariants"] and source.warnings == []


def test_an_unparsed_date_keeps_header_order_and_says_so() -> None:
    source = detect_source([GATK3_SELECT.replace("Wed Mar 03 10:00:00 CST 2021", "yesterday"), GATK4_HC_PHASED])
    assert [s.tool for s in source.steps] == ["GATK SelectVariants", "GATK HaplotypeCaller"]
    assert source.warnings == ["a step's stated date did not parse; the method chain is in header order"]


def test_every_role_is_in_the_vocabulary() -> None:
    meta = [DRAGEN_MAPPED, DRAGEN_PS, GATK3_APPLY, GATK3_SELECT, NORM, *CLAIR3,
            "##liftOverProgram=<CrossMap,version=0.6.3,website=x>", "##VEP=v88 assembly=GRCh37",
            "##bcftools_concatCommand=concat -a; Date=Tue Mar 19 03:30:00 2024", "##gvcftools_cmdline=/opt/set_haploid_region --ref /ref/hg19.fa"]
    roles = {step.role for step in detect_source(meta).steps}
    assert roles == STEP_ROLES
