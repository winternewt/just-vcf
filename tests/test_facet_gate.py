"""The facet gate: every loadable format, loaded from its spec's `example`, fills every `facets` column
unless this file says why not.

Each loader builds its own `facets` row (VCF, the vendor tables, the 23andMe layouts, CGI's native
files, the key=value table), so a check one builder adds to one loader can silently miss the others,
and each loader's own tests still pass. Here a column left null needs a recorded reason, and a
recorded reason that no longer holds fails too. Adding a facet column, or filling one in a new place,
therefore means touching this table, which is the point: the decision is made for every loader at
once. The examples are synthetic, so this costs no genome and runs in the default suite.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import FACETS, LOADABLE_FORMATS, load
from just_vcf.formats import get_spec

VCF_FAMILY = frozenset({"vcf", "vcf_cgi", "vcf_freebayes_varscan2", "vcf_no_meta", "bcf"})
"""Formats read by the VCF loader: the header columns are theirs alone."""

CALLER_SPLIT = frozenset({"vcf_freebayes_varscan2", "freebayes_varscan2_keyvalue"})
"""Formats loaded as one facet per caller (`relation = caller_split`)."""

ALL = LOADABLE_FORMATS

NULL_REASONS: dict[str, dict[frozenset[str], str]] = {
    "seqcol_digest": {ALL: "needs the bases (RM3); a header cannot give it (SCHEMA Decision 5)"},
    "derived_from": {ALL: "no example is derived from a loaded facet; only rejoin_split sets it"},
    "split_marker": {ALL: "no example tags a record as split off a multi-allelic site"},
    "member": {ALL: "the examples are written as plain files, not archive members"},
    "phasing_method": {ALL: "no example's header names a phasing tool"},
    "sorted_name_length_pairs_digest": {ALL - {"bcf"}: "no text example declares ##contig lines with lengths"},
    "header_contigs": {ALL - VCF_FAMILY: "only a VCF header has ##contig lines"},
    "vcf_version": {ALL - VCF_FAMILY - {"vcf_no_meta"}: "only a VCF header states a version", frozenset({"vcf_no_meta"}): "the file has no meta lines"},
    "ploidy_statements": {ALL - VCF_FAMILY: "null means not a VCF (SCHEMA Decision 5)"},
    "sample_meta": {ALL - VCF_FAMILY: "only a VCF header has ##SAMPLE and ##PEDIGREE lines"},
    "relation": {ALL - CALLER_SPLIT: "a facet of one file and one caller has no relation"},
    "caller": {
        ALL - VCF_FAMILY - {"cgi_var", "cgi_master_var", "freebayes_varscan2_keyvalue"}: "an array or table export names no caller",
        frozenset({"vcf", "vcf_no_meta", "bcf"}): "the examples' headers name no caller",
    },
    "caller_version": {
        ALL - VCF_FAMILY - {"cgi_var", "cgi_master_var"}: "no caller version is stated",
        frozenset({"vcf", "vcf_no_meta", "bcf"}): "the examples' headers state no caller version",
        frozenset({"vcf_freebayes_varscan2"}): "VarScan2's header line states no version",
    },
    "vendor": {
        frozenset({"vcf", "vcf_no_meta", "bcf"}): "a generic VCF names no vendor in the examples",
        frozenset({"rsid_result_csv", "rsid_table", "rsid_table_headerless"}): "a vendor-neutral layout claims no vendor",
    },
    "sample": {
        frozenset({"23andme", "23andme_csv", "ancestrydna", "decodeme", "ftdna", "ftdna_space", "genes_for_good",
                   "genes_for_good_imputed", "living_dna", "plink_23", "rsid_result_csv", "rsid_table",
                   "rsid_table_headerless", "selfdecode", "genos"}): "the format has no sample column or line",
        frozenset({"illumina_final_report"}): "the example names no Sample ID column",
        frozenset({"cgi_var", "cgi_master_var"}): "the examples have no #ASSEMBLY_ID line",
        frozenset({"freebayes_varscan2_keyvalue"}): "the table has no sample column",
    },
    "dropped_fields": {
        frozenset({"23andme", "23andme_csv", "genes_for_good", "genes_for_good_imputed", "living_dna", "plink_23",
                   "rsid_table", "rsid_table_headerless", "selfdecode"}): "the 23andMe layout has no columns beyond the four it reads",
    },
    "callset_type": {
        frozenset({"plink_23", "genos", "color_genotypes"}): "the file does not say which assay (RM4)",
        frozenset({"cgi_var", "cgi_master_var"}): "the examples have no #TITLE",
        frozenset({"vcf", "vcf_cgi", "vcf_freebayes_varscan2", "vcf_no_meta", "bcf", "freebayes_varscan2_keyvalue"}):
            "a few example rows fit no callset class",
    },
    "assembly_label": {
        frozenset({"decodeme", "ftdna", "ftdna_space", "illumina_final_report", "plink_23", "rsid_result_csv",
                   "rsid_table", "rsid_table_headerless", "freebayes_varscan2_keyvalue", "vcf_no_meta"}): "the example states no build",
    },
    "assembly_label_basis": {
        frozenset({"decodeme", "ftdna", "ftdna_space", "illumina_final_report", "plink_23", "rsid_result_csv",
                   "rsid_table", "rsid_table_headerless", "freebayes_varscan2_keyvalue", "vcf_no_meta"}): "no label, no basis",
    },
}
"""Facet column → the formats whose example leaves it null → why. A column absent here must be filled
for every format."""


def _reason(column: str, format_id: str) -> str | None:
    return next((why for formats, why in NULL_REASONS.get(column, {}).items() if format_id in formats), None)


def _write_example(format_id: str, tmp_path: Path) -> Path:
    if format_id == "bcf":
        pysam = pytest.importorskip("pysam", reason="pysam (the [reference] extra) writes the BCF example")
        text = tmp_path / "example.vcf"
        # A BCF indexes contigs by the header, and the text example declares none.
        lines = get_spec("vcf").example.splitlines(keepends=True)
        text.write_text("".join([lines[0], "##contig=<ID=1,length=249250621>\n", *lines[1:]]))
        target = tmp_path / "example.bcf"
        with pysam.VariantFile(str(text)) as source, pysam.VariantFile(str(target), "wb", header=source.header) as out:
            for record in source:
                out.write(record)
        return target
    path = tmp_path / f"{format_id}.example"
    path.write_text(get_spec(format_id).example)
    return path


def test_the_table_names_only_real_columns_and_formats() -> None:
    assert set(NULL_REASONS) <= set(FACETS.names())
    for column, entries in NULL_REASONS.items():
        named = [f for formats in entries for f in formats]
        assert set(named) <= LOADABLE_FORMATS, column
        assert len(named) == len(set(named)), f"{column}: a format has two reasons"


@pytest.mark.parametrize("format_id", sorted(LOADABLE_FORMATS))
def test_every_facet_column_is_filled_or_has_a_reason(format_id: str, tmp_path: Path) -> None:
    facets = load(_write_example(format_id, tmp_path)).facets.collect()
    assert facets.height >= 1
    # A column null on any facet needs a reason; a reason for a column every facet fills is stale. The
    # facets of one file may differ (a merge's VarScan2 facet states no version where freeBayes does).
    null = {name for name in FACETS.names() if facets[name].null_count() > 0}
    unexplained = sorted(name for name in null if _reason(name, format_id) is None)
    stale = sorted(
        column for column in NULL_REASONS if _reason(column, format_id) is not None and column not in null
    )
    assert unexplained == [], f"{format_id}: null without a reason in NULL_REASONS: {unexplained}"
    assert stale == [], f"{format_id}: filled now, so the recorded reason is stale: {stale}"


def test_identity_columns_are_never_null() -> None:
    """The columns every loader must fill whatever the file says are not in the table at all."""
    always = {"facet_id", "source_sha256", "source_path", "format_id", "family", "absence_means", "contig_map"}
    assert always.isdisjoint(NULL_REASONS)
    assert always <= set(FACETS.names())
    assert pl.Series(sorted(always)).is_unique().all()
