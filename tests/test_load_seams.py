"""Merge seams in the loaders: every facet column a load leaves null has a stated reason.

A hand merge once left a per-facet column filled by the VCF loader and by no other, and every test
still passed, because each loader's tests check the columns that loader was written to fill. This
module loads every loadable format's spec `example` through `just_vcf.load` and fails on any `FACETS`
column that is null in a facet unless `NULL_REASONS` says why null is right there. **A new `FACETS`
column therefore fails here for every format until someone states, format by format, why it is null
or makes the loaders fill it.** That is the point: the reasons are the audit, written down.

`NULL_REASONS` maps `(key, column)` to a reason. A key covers facets by, from narrowest to widest:

- `<format_id>/<caller>` for one caller's facet of a caller split, `<format_id>/<relation>` for a
  facet of another kin relation (an appended export), or else `<format_id>`;
- a family value (`genotype_table`, `variant_calls`, `callset_native`);
- `*`, every facet.

A reason must be true of every facet its key covers (`test_every_reason_holds_for_every_facet_it_covers`):
a family-wide or `*` reason for a column that one format of the family fills is refused, so a broad
reason cannot hide the next loader that forgets the column. Where some formats fill a column and
others do not, the others are named one by one. A reason that covers no loaded facet, or whose column
a covered facet fills, is stale and fails too. Rejected: letting the narrowest key win and checking
only it, which would let a family-wide reason excuse a format that fills the column today and drops
it tomorrow.

`KNOWN_GAPS` has the same shape and rule: nulls this guard found that look wrong, excused until the
loaders are fixed, and failing once they are, so the fix removes its entry.

Each example is loaded as the member `example.txt` of a zip archive, so `member` is filled on every
load path. Rejected: loading the plain file and excusing `member` everywhere, which would leave the
one column every loader copies from `sniff` unchecked. BCF has no text example: the `vcf` example is
written as BCF by pysam, with a `##contig=<ID=1>` line added, since htslib refuses a record whose
contig the header does not declare. The appended-export path (F3) has no spec example either; it is
loaded from `test_load_rsid`'s two-export table.

The same loads also check that every load path returns a `links` frame of the `LINKS` schema that
agrees with the facets (`just_vcf.links.link_disagreements`), and that a `subject` passed to `load`
reaches every facet.
"""

import zipfile
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pytest

from just_vcf import LOADABLE_FORMATS, load
from just_vcf.formats import REGISTRY, get_spec
from just_vcf.links import CALLER_SPLIT, link_disagreements
from just_vcf.schema import FACETS, LINKS
from test_load_rsid import _two_exports

try:
    import pysam
except ImportError:  # the [reference] extra; only the BCF case needs it
    pysam = None

SUBJECT = "seam-guard-subject"
MEMBER = "example.txt"

NEVER_EXCUSED = frozenset({"facet_id", "subject"})
"""`facet_id` is required; `subject` is passed to every load here, so a null one is a bug."""

GT, VC, CN = "genotype_table", "variant_calls", "callset_native"
FB_VCF = ("vcf_freebayes_varscan2/freeBayes", "vcf_freebayes_varscan2/VarScan2")
FB_KEYVALUE = ("freebayes_varscan2_keyvalue/freeBayes", "freebayes_varscan2_keyvalue/VarScan2")
APPENDED = "rsid_table_headerless/appended_export"
SINGLE_TABLES = (
    "illumina_final_report", "ancestrydna", "genos", "color_genotypes", "decodeme", "ftdna", "ftdna_space",
    "genes_for_good_imputed", "genes_for_good", "23andme", "23andme_csv", "living_dna", "selfdecode", "plink_23",
    "rsid_result_csv", "rsid_table",
)
"""Genotype-table formats whose load returns one facet (`rsid_table_headerless` can return two)."""


def _each(keys: tuple[str, ...], columns: tuple[str, ...], reason: str) -> dict[tuple[str, str], str]:
    return {(key, column): reason for key in keys for column in columns}


_NOT_A_VCF = "Not a VCF: these columns come from a VCF header, and null means there is none (SCHEMA Decision 5)."
_ONE_FACET = "A load that returns one facet has no kin, and nothing is derived at load."
_CONTENT_FILTER = (
    "The example's rows are far below the 10,000 alleles a rows basis needs (`callset.MIN_FILTER_ALLELES`), no "
    "header line records a filtering command, and the member name `example.txt` carries no class word; the basis "
    "is null with the filter."
)
_TOO_FEW_SITES = (
    "The example's handful of rows falls below every RM4 band (a panel needs 50 sites, `callset.MIN_PANEL_SITES`), "
    "and no header line states array or imputed content."
)
_NO_CHIP_RULE = "The format states no chip and has no probe-set rule (`chip.ChipOutcome.NO_RULE`)."
_FEW_OVERLAP = (
    "Marker overlap names no Illumina family: the example's few rsids are far below the 1,000 the marker "
    "sample needs (`chip.overlap_chip`, a load warning)."
)
_NO_BAND_23 = "The example's row count lies in no 23andMe generation's band (`chip.ChipOutcome.NO_BAND`, a load warning)."
_NO_SAMPLE = "The layout has no sample column; the file holds one person by construction."

NULL_REASONS: dict[tuple[str, str], str] = {
    # Every facet `load` returns.
    ("*", "seqcol_digest"): (
        "No header and no row gives it: the digest is taken over each sequence's refget digest, which needs the "
        "bases (SCHEMA Decision 5, 'Assembly identity from a VCF header'). Null until RM3 reads them; the label "
        "never stands in for it."
    ),
    ("*", "ref_seqcol_digest"): "Null until `resolve_ref` reads REF from a reference collection; `load` reads none.",
    ("*", "strand_check"): "Null until `resolve_ref` compares an array's heterozygous rows with the reference; `load` does not run it.",
    ("*", "derived_from"): (
        "A facet read from a file has no parent. A caller split's and an appended export's facets have none "
        "either: they are kin, linked in `links` (SCHEMA Decision 5, 12)."
    ),
    ("*", "derivation_inputs"): "Null on a facet read from a file: only a derived step writes it (a liftover's chain).",
    **_each(('23andme_csv', '23andme', 'ancestrydna', 'cgi_master_var', 'cgi_var', 'color_genotypes', 'decodeme', 'ftdna', 'ftdna_space', 'genes_for_good_imputed', 'genes_for_good', 'genos', 'illumina_final_report', 'living_dna', 'plink_23', 'rsid_result_csv', 'rsid_table_headerless', 'rsid_table', 'selfdecode'), ("method_steps",), _NOT_A_VCF),
    **_each(("freebayes_varscan2_keyvalue/VarScan2", "freebayes_varscan2_keyvalue/freeBayes"), ("method_steps",), "The key=value table has no header lines to record a method chain."),
    ("rsid_table_headerless/appended_export", "method_steps"): _NOT_A_VCF,
    ("vcf_no_meta", "reference"): "A VCF with no meta lines has no `##reference=` and no recorded command line.",
    # The site tables and Nirvana JSON (layer 1, 2026-10-23): what each states and what it does not.
    **_each(('casava_snp_table', 'casava_indel_table', 'soapsnp_cns', 'ftdna_bigy', 'color_variant_set', 'get_evidence_gff', 'nirvana_json'), ("caller", "caller_version"), "The format names no caller (docs/ROADMAP.md RM2, site tables)."),
    **_each(('casava_snp_table', 'casava_indel_table', 'soapsnp_cns', 'ftdna_bigy', 'color_variant_set', 'get_evidence_gff', 'nirvana_json'), ("ploidy_statements", "header_contigs", "vcf_version", "sample_meta", "method_steps", "reference", "split_marker"), _NOT_A_VCF),
    **_each(('casava_snp_table', 'casava_indel_table', 'soapsnp_cns', 'ftdna_bigy', 'color_variant_set', 'get_evidence_gff', 'nirvana_json'), ("phasing_method",), "The format states no phasing method."),
    **_each(('casava_snp_table', 'casava_indel_table', 'soapsnp_cns', 'ftdna_bigy', 'color_variant_set', 'get_evidence_gff', 'nirvana_json'), ("relation",), _ONE_FACET),
    **_each(("casava_snp_table", "casava_indel_table", "soapsnp_cns", "get_evidence_gff", "ftdna_bigy"), ("sample",), _NO_SAMPLE),
    **_each(("casava_snp_table", "casava_indel_table", "soapsnp_cns", "ftdna_bigy", "nirvana_json"), ("assembly_label", "assembly_label_basis"), (
        "The file states no build; its positions were measured against the reference stores, but a measured build is "
        "not a stated one (docs/LAYERS.md, layer 3 blocker)."
    )),
    **_each(("color_variant_set", "get_evidence_gff", "nirvana_json", "ftdna_bigy"), ("callset_type",), _TOO_FEW_SITES),
    **_each(("color_variant_set", "get_evidence_gff", "nirvana_json"), ("content_filter", "content_filter_basis"), _CONTENT_FILTER),
    **_each(("color_variant_set", "get_evidence_gff", "nirvana_json"), ("sorted_name_length_pairs_digest",), "The format carries no contig lengths."),
    **_each(("color_variant_set", "get_evidence_gff", "nirvana_json", "ftdna_bigy"), ("chip", "chip_basis"), _NO_CHIP_RULE),
    ("rsid_table_headerless/appended_export", "variant_classes"): "Not counted: a genotype table writes no REF or ALT (SCHEMA Decision 5).",
    ("*", "aligner"): "Only a header that records its aligner states one (DRAGEN with `--enable-map-align true`); no example does (SCHEMA Decision 5).",
    ("*", "aligner_version"): "Null with `aligner`: only a recorded aligner has a version.",
    **_each(('23andme_csv', '23andme', 'ancestrydna', 'cgi_master_var', 'cgi_var', 'color_genotypes', 'decodeme', 'ftdna', 'ftdna_space', 'genes_for_good_imputed', 'genes_for_good', 'genos', 'illumina_final_report', 'living_dna', 'plink_23', 'rsid_result_csv', 'rsid_table_headerless', 'rsid_table', 'selfdecode', 'freebayes_varscan2_keyvalue/VarScan2', 'freebayes_varscan2_keyvalue/freeBayes', 'rsid_table_headerless/appended_export'), ("reference",), "Not a VCF: no `##reference=` line or caller reference option to read (SCHEMA Decision 5)."),
    **_each(('23andme', '23andme_csv', 'ancestrydna', 'bcf', 'cgi_master_var', 'cgi_var', 'color_genotypes', 'decodeme', 'ftdna', 'ftdna_space', 'genes_for_good', 'genes_for_good_imputed', 'genos', 'illumina_final_report', 'living_dna', 'plink_23', 'rsid_result_csv', 'rsid_table', 'rsid_table_headerless', 'selfdecode', 'vcf', 'vcf_no_meta', 'freebayes_varscan2_keyvalue/VarScan2', 'freebayes_varscan2_keyvalue/freeBayes', 'rsid_table_headerless/appended_export', 'vcf_freebayes_varscan2/freeBayes', 'vcf_freebayes_varscan2/VarScan2'), ("phasing_method",), (
        "Null unless the header names a phasing step (`source.method_chain`): no other example does. `vcf_cgi`'s "
        "header declares PS in CGAPipeline's wording, so it names CGAPipeline."
    )),
    # Genotype tables.
    **_each((GT,), ("caller", "caller_version"), (
        "The RM2 source detector reads VCF header lines; a genotype table's loader names no calling program, and "
        "an array export states none."
    )),
    **_each((GT,), ("ploidy_statements", "header_contigs", "vcf_version", "sample_meta"), _NOT_A_VCF),
    **_each(('23andme', '23andme_csv', 'ancestrydna', 'color_genotypes', 'decodeme', 'ftdna', 'ftdna_space', 'genes_for_good', 'genes_for_good_imputed', 'genos', 'illumina_final_report', 'living_dna', 'plink_23', 'rsid_result_csv', 'rsid_table', 'rsid_table_headerless', 'selfdecode'), ("variant_classes",), "Not counted: a genotype table writes no REF or ALT (SCHEMA Decision 5)."),
    **_each((GT,), ("content_filter", "content_filter_basis"), (
        "Its variant classes are not counted (no REF or ALT), so no class can be the one the rows were cut to."
    )),
    (GT, "split_marker"): "A genotype table has no INFO field to carry a split marker.",
    (GT, "sorted_name_length_pairs_digest"): "Not a VCF: no `##contig` lines to take names and lengths from.",
    **_each(SINGLE_TABLES, ("relation",), _ONE_FACET),
    ("rsid_table_headerless", "relation"): (
        "A load of one facet, or the first export of an appended pair, which the appended one links to and which "
        "is not itself appended (F3)."
    ),
    **_each(("rsid_result_csv", "rsid_table", "rsid_table_headerless", APPENDED), ("vendor",), (
        "A generic layout: the spec's vendor is `unknown`, and nothing in the file names one."
    )),
    **_each(tuple(f for f in SINGLE_TABLES if f not in ("color_genotypes", "illumina_final_report", "genos")) + ("rsid_table_headerless", APPENDED), ("sample",), _NO_SAMPLE),
    ("illumina_final_report", "sample"): "The example's `[Data]` section has no `Sample ID` column, the one place a Final Report names its sample.",
    ("genos", "sample"): "The example has no `##Sample=` line, the one place a Genos export names its sample.",
    **_each(("illumina_final_report",), ("assembly_label", "assembly_label_basis"), (
        "A Final Report's header names no build (the spec has no build hint); `resolve_ref` infers one from the "
        "rows (F9)."
    )),
    **_each(("decodeme",), ("assembly_label", "assembly_label_basis"), (
        "deCODEme states its build in the info member beside the genotypes (`_sibling_build`); this one-member "
        "archive has none."
    )),
    **_each(("ftdna", "ftdna_space", "rsid_result_csv", "rsid_table", "rsid_table_headerless"), ("assembly_label", "assembly_label_basis"), (
        "The layout states no build (the spec has no build hint); `resolve_ref` infers one from the rows (F9)."
    )),
    **_each(("plink_23",), ("assembly_label", "assembly_label_basis"), (
        "PLINK's preamble names no build: it says to consult the original source of the data."
    )),
    **_each((APPENDED,), ("assembly_label", "assembly_label_basis"), (
        "A preamble that states a build belongs to the first export; the appended one states none of its own (F9)."
    )),
    ("genos", "callset_type"): "Genos reports exome sequencing at a request list of sites, and the file does not say so (`_UNKNOWN_CALLSET_FORMATS`).",
    ("color_genotypes", "callset_type"): "Color does not say which assay its trait SNPs come from (`_UNKNOWN_CALLSET_FORMATS`).",
    ("plink_23", "callset_type"): "PLINK writes this layout from any source; the file does not say which (`_UNKNOWN_CALLSET_FORMATS`).",
    **_each(("genos", "color_genotypes", "decodeme", "genes_for_good_imputed", "genes_for_good", "selfdecode", "plink_23"), ("chip", "chip_basis"), _NO_CHIP_RULE),
    **_each(("rsid_table", "rsid_table_headerless"), ("chip", "chip_basis"), (
        "The layout states no chip, and is read as 23andMe's only when it carries 23andMe's `i…` ids, which the "
        "example does not (`chip._23ANDME_LAYOUT`)."
    )),
    **_each(("ftdna", "ftdna_space", "rsid_result_csv"), ("chip", "chip_basis"), _FEW_OVERLAP),
    **_each(("23andme", "23andme_csv"), ("chip", "chip_basis"), _NO_BAND_23),
    **_each((APPENDED,), ("chip", "chip_basis"), (
        "Not inferred: the 23andMe rule reads the first layout's rows, and no rule reads the appended layout's (F3)."
    )),
    # Complete Genomics' native files.
    **_each((CN,), ("ploidy_statements", "header_contigs", "vcf_version", "sample_meta"), _NOT_A_VCF),
    (CN, "sorted_name_length_pairs_digest"): "Not a VCF: no `##contig` lines to take names and lengths from.",
    (CN, "sample"): "The example has no `#ASSEMBLY_ID` line, which is where a native file names its sample (SCHEMA Decision 5).",
    (CN, "callset_type"): (
        "The example's header has no `#TITLE` stating a whole genome, and its few rows fall below every RM4 band, "
        "so the rows do not say either."
    ),
    **_each((CN,), ("content_filter", "content_filter_basis"), _CONTENT_FILTER),
    (CN, "split_marker"): "CGI's native files have no INFO field to carry a split marker.",
    **_each((CN,), ("chip", "chip_basis"), "Not an array: `chip` names a genotyping chip."),
    (CN, "relation"): _ONE_FACET,
    # VCFs, BCF and the key=value table.
    **_each((VC,), ("chip", "chip_basis"), "Not an array: `chip` names a genotyping chip."),
    (VC, "callset_type"): _TOO_FEW_SITES,
    **_each((VC,), ("content_filter", "content_filter_basis"), _CONTENT_FILTER),
    (VC, "split_marker"): (
        "No example record is tagged as split off a multi-allelic site (vt's `OLD_MULTIALLELIC`, `JUST_VCF_SPLIT`, "
        "a bcftools `--old-rec-tag` key); the key=value table has no INFO field at all."
    ),
    (VC, "sorted_name_length_pairs_digest"): (
        "No example header has a `##contig` line with a length (the one added for the BCF case has none), and the "
        "key=value table has no header."
    ),
    **_each(("vcf", "vcf_cgi", "vcf_no_meta", "bcf"), ("relation",), _ONE_FACET),
    **_each(("vcf", "bcf"), ("vendor", "caller", "caller_version"), (
        "The header names no vendor and no caller (no `##center`, `##source` or command line), and the generic "
        "VCF format implies no vendor."
    )),
    **_each(("vcf_no_meta",), ("vendor", "caller", "caller_version"), "No meta lines: nothing names a vendor or a caller."),
    **_each(("vcf_no_meta",), ("assembly_label", "assembly_label_basis"), "No meta lines: no build hint and no `##contig` lengths."),
    ("vcf_no_meta", "vcf_version"): "No `##fileformat` line: the file states no version (it is read under a synthesised VCFv4.2 header).",
    ("vcf_freebayes_varscan2/VarScan2", "caller_version"): "`##source=VarScan2` carries no version.",
    **_each(FB_KEYVALUE, ("caller_version", "sample"), "The key=value table has no header and no sample column."),
    **_each(FB_KEYVALUE, ("assembly_label", "assembly_label_basis"), "The key=value table has no header to state a build (the spec has no build hint)."),
    **_each(FB_KEYVALUE, ("ploidy_statements", "header_contigs", "vcf_version", "sample_meta"), _NOT_A_VCF),
}
"""(key, column) → why null is right for the facets the key covers."""

KNOWN_GAPS: dict[tuple[str, str], str] = {}
"""Nulls this guard found that look wrong, excused until the loaders are fixed. Held to the same rule as a
reason: an entry fails as soon as every facet it covers fills the column, so the fix removes it."""

EXCUSED = {**NULL_REASONS, **KNOWN_GAPS}


def _zipped(tmp_path: Path, payload: bytes, member: str = MEMBER) -> Path:
    path = tmp_path / "example.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(member, payload)
    return path


def _example_case(format_id: str) -> Callable[[Path], Path]:
    def build(tmp_path: Path) -> Path:
        return _zipped(tmp_path, get_spec(format_id).example.encode())

    return build


def _bcf_case(tmp_path: Path) -> Path:
    text = tmp_path / "example.vcf"
    text.write_text(get_spec("vcf").example.replace("##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1>\n"))
    bcf = tmp_path / "example.bcf"
    with pysam.VariantFile(str(text)) as source, pysam.VariantFile(str(bcf), "wb", header=source.header) as out:
        for record in source:
            out.write(record)
    return _zipped(tmp_path, bcf.read_bytes(), "example.bcf")


def _appended_export_case(tmp_path: Path) -> Path:
    text, _, _ = _two_exports()
    return _zipped(tmp_path, text.encode())


CASES: dict[str, Callable[[Path], Path]] = {
    **{spec.id: _example_case(spec.id) for spec in REGISTRY if spec.id in LOADABLE_FORMATS and spec.example is not None},
    "bcf": _bcf_case,
    "appended_export": _appended_export_case,
}
"""Case id → a builder writing the case's input under a directory and returning its path."""

NEEDS_PYSAM = frozenset({"bcf"})


class Case:
    def __init__(self, facets: pl.DataFrame, links: pl.DataFrame, links_schema: pl.Schema) -> None:
        self.facets = facets
        self.links = links
        self.links_schema = links_schema


@pytest.fixture(scope="module")
def loads(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Case]:
    """Every case loaded once, under a cache of its own (the autouse cache fixture is per test)."""
    out: dict[str, Case] = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("JUST_VCF_CACHE_DIR", str(tmp_path_factory.mktemp("seam-cache")))
        patch.setenv("JUST_VCF_REFERENCE_DIR", "")
        patch.setenv("JUST_VCF_MANIFESTS_DIR", "")
        for case_id, build in CASES.items():
            if case_id in NEEDS_PYSAM and pysam is None:
                continue
            loaded = load(build(tmp_path_factory.mktemp(case_id)), subject=SUBJECT)
            facets, links = pl.collect_all([loaded.facets, loaded.links])
            out[case_id] = Case(facets, links, loaded.links.collect_schema())
    return out


def _case(loads: dict[str, Case], case_id: str) -> Case:
    if case_id not in loads:
        pytest.skip("pysam (the [reference] extra) writes the BCF case")
    return loads[case_id]


def facet_keys(facet: dict) -> list[str]:
    """The keys of `NULL_REASONS` that cover `facet`, narrowest first.

    A kin facet (a caller split's, an appended export) is covered by its own key, not by its format's:
    its columns are filled by the split's rules, and a format-wide reason must not reach it unseen."""
    if facet["relation"] == CALLER_SPLIT and facet["caller"]:
        own = f"{facet['format_id']}/{facet['caller']}"
    elif facet["relation"]:
        own = f"{facet['format_id']}/{facet['relation']}"
    else:
        own = facet["format_id"]
    return [own, facet["family"], "*"]


def _unexplained(facets: pl.DataFrame) -> list[str]:
    missing = []
    for facet in facets.iter_rows(named=True):
        keys = facet_keys(facet)
        for column in FACETS.names():
            if facet[column] is None and not any((key, column) in EXCUSED for key in keys):
                missing.append(f"{keys[0]}: {column}")
    return missing


def test_every_loadable_format_has_a_case() -> None:
    """A format that gains a loader is loaded here, from its example or from a case written for it."""
    assert LOADABLE_FORMATS <= set(CASES), sorted(LOADABLE_FORMATS - set(CASES))


@pytest.mark.parametrize("case_id", list(CASES))
def test_every_null_facet_column_has_a_reason(loads: dict[str, Case], case_id: str) -> None:
    assert _unexplained(_case(loads, case_id).facets) == []


@pytest.mark.parametrize("case_id", list(CASES))
def test_a_subject_reaches_every_facet(loads: dict[str, Case], case_id: str) -> None:
    facets = _case(loads, case_id).facets
    assert facets["subject"].to_list() == [SUBJECT] * facets.height


@pytest.mark.parametrize("case_id", list(CASES))
def test_every_load_path_returns_links_that_agree_with_its_facets(loads: dict[str, Case], case_id: str) -> None:
    case = _case(loads, case_id)
    assert case.links_schema == LINKS and case.links.schema == LINKS
    assert link_disagreements(case.facets, case.links) == []
    ids = set(case.facets["facet_id"])
    endpoints = set(case.links["from_facet_id"]) | set(case.links["to_facet_id"])
    assert endpoints <= ids
    # One facet links to nothing; several facets of one source are kin, every one of them linked.
    assert endpoints == (ids if len(ids) > 1 else set())


def test_every_reason_holds_for_every_facet_it_covers(loads: dict[str, Case]) -> None:
    """A reason is refused when a facet it covers fills the column, or when it covers no facet at all."""
    facets = pl.concat([case.facets for case in loads.values()], how="vertical")
    rows = list(facets.iter_rows(named=True))
    skipped = set(CASES) - set(loads)
    problems = []
    assert not set(NULL_REASONS) & set(KNOWN_GAPS), "a null is either right or a gap, not both"
    for (key, column), reason in EXCUSED.items():
        assert column in FACETS and column not in NEVER_EXCUSED, (key, column)
        assert reason.strip(), (key, column)
        covered = [row for row in rows if key in facet_keys(row)]
        if not covered and not (skipped and key in NEEDS_PYSAM):
            problems.append(f"{key}: {column} covers no facet")
        problems += [f"{key}: {column} is filled in {facet_keys(row)[0]}" for row in covered if row[column] is not None]
    assert problems == []
