"""RM287: an imputed VCF's per-site quality kept on the row (schema 0.0.23).

minimac writes INFO R2 (its estimated r²) on every site and ER2 (the empirical, leave-one-out r²) only
on the sites the chip typed; Beagle writes DR2 (its estimated r², one per ALT). Every fixture is a
synthetic inline VCF; each expected value is read back from the records the test wrote.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import load
from just_vcf.schema import CALLS
from test_load_vcf import write

GT_DS = (
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">',
)
MINIMAC_META = (
    "##fileformat=VCFv4.2",
    "##source=Minimac4.v1.0.2",
    "##contig=<ID=1,length=249250621>",
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">',
    '##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">',
    '##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">',
    '##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">',
    '##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">',
    *GT_DS,
)
BEAGLE_META = (
    "##fileformat=VCFv4.2",
    '##source="beagle.22Jul22.46e.jar"',
    "##contig=<ID=1,length=249250621>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">',
    '##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose and true REF dose">',
    '##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">',
    *GT_DS,
)
COLUMNS = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])

MINIMAC_SITES = [
    # pos, ref, alt, INFO, GT:DS
    (100, "A", "G", "AF=0.1;R2=0.95;IMPUTED", "0|1:1.02"),
    (200, "C", "T", "AF=0.2;R2=0.88;ER2=0.91;TYPED", "0|0:0.05"),
    (300, "G", "A", "AF=0.3;R2=0.21;IMPUTED", "1|1:1.7"),
    (400, "T", "C", "AF=0.3;ER2=0.75;TYPED", "0|1:1"),
    (500, "A", "C", "AF=0.4;IMPUTED", "0|1:1"),
]
"""R2 on every imputed site, ER2 beside R2 on a typed one, ER2 alone on one, and a site with neither."""

BEAGLE_SITES = [
    (100, "A", "G", "AF=0.1;DR2=0.97", "0|1:1"),
    (200, "C", "T,G", "AF=0.1,0.2;DR2=0.5,0.6;IMP", "0|2:0.1,0.9"),
    (300, "G", "A", "AF=0.1;DR2=0.33;IMP", "0|0:0.2"),
    (400, "T", "C", "AF=0.1;DR2=.;IMP", "0|0:0.2"),
]
"""DR2 one per ALT: a multi-allelic site writes two, which no one-value column holds."""


def vcf(tmp_path: Path, meta: tuple[str, ...], sites: list[tuple], name: str) -> Path:
    records = [f"1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:DS\t{sample}" for pos, ref, alt, info, sample in sites]
    return write(tmp_path, "\n".join([*meta, COLUMNS, *records]) + "\n", name)


def info_value(info: str, key: str) -> list[float | None] | None:
    """What the record's INFO writes for `key`: its values (`.` as None), or None when absent."""
    for item in info.split(";"):
        name, _, value = item.partition("=")
        if name == key:
            return [None if part == "." else float(part) for part in value.split(",")]
    return None


def one_value(info: str, key: str) -> float | None:
    """The single value a record writes for `key`, as Float32 holds it; None when it writes none or several."""
    values = info_value(info, key)
    if values is None or len(values) != 1 or values[0] is None:
        return None
    return pl.Series([values[0]], dtype=pl.Float32).item()


def by_pos(path: Path) -> tuple[dict[int, dict], dict, tuple[str, ...]]:
    loaded = load(path)
    calls, facets = pl.collect_all([loaded.calls, loaded.facets])
    return {row["pos"]: row for row in calls.iter_rows(named=True)}, facets.row(0, named=True), loaded.warnings


def test_the_columns_sit_beside_dosage_as_floats() -> None:
    names = CALLS.names()
    at = names.index("dosage")
    assert names[at + 1 : at + 3] == ["imputation_r2", "imputation_empirical_r2"]
    assert CALLS["imputation_r2"] == CALLS["imputation_empirical_r2"] == pl.Float32


def test_minimac_r2_on_every_site_and_er2_on_the_typed_ones(tmp_path: Path) -> None:
    rows, facet, warnings = by_pos(vcf(tmp_path, MINIMAC_META, MINIMAC_SITES, "minimac.dose.vcf"))
    assert set(rows) == {site[0] for site in MINIMAC_SITES}
    for pos, _, _, info, _ in MINIMAC_SITES:
        assert rows[pos]["imputation_r2"] == one_value(info, "R2"), pos
        assert rows[pos]["imputation_empirical_r2"] == one_value(info, "ER2"), pos
    assert rows[200]["imputation_r2"] is not None and rows[200]["imputation_empirical_r2"] is not None, "a typed site has both"
    assert rows[500]["imputation_r2"] is None and rows[500]["imputation_empirical_r2"] is None
    # Carried now, so no longer listed as dropped; AF still is.
    assert {"INFO/R2", "INFO/ER2"}.isdisjoint(facet["dropped_fields"]) and "INFO/AF" in facet["dropped_fields"]
    assert facet["callset_type"] == "imputed"
    assert not any("R2" in line for line in warnings)


def test_beagle_dr2_is_read_where_the_site_has_one_alt(tmp_path: Path) -> None:
    rows, facet, warnings = by_pos(vcf(tmp_path, BEAGLE_META, BEAGLE_SITES, "beagle.vcf"))
    for pos, _, alt, info, _ in BEAGLE_SITES:
        expected = one_value(info, "DR2") if "," not in alt else None
        assert rows[pos]["imputation_r2"] == expected, pos
        assert rows[pos]["imputation_empirical_r2"] is None, "Beagle writes no empirical r²"
    assert rows[100]["imputation_r2"] is not None and rows[400]["imputation_r2"] is None
    several = sum(1 for _, _, _, info, _ in BEAGLE_SITES if len(info_value(info, "DR2") or []) > 1)
    assert several == 1
    assert any(line.startswith(f"{several} INFO DR2 value(s)") and "withheld" in line for line in warnings), warnings
    assert "INFO/DR2" not in facet["dropped_fields"]


def test_a_vcf_without_imputation_keys_has_both_columns_null(tmp_path: Path) -> None:
    meta = ("##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621>", *GT_DS)
    sites = [(pos, ref, alt, ".", gt) for pos, ref, alt, _, gt in MINIMAC_SITES]
    rows, _, warnings = by_pos(vcf(tmp_path, meta, sites, "plain.vcf"))
    assert len(rows) == len(sites)
    assert all(row["imputation_r2"] is None and row["imputation_empirical_r2"] is None for row in rows.values())
    assert not any("R2" in line for line in warnings)


def test_r2_and_dr2_that_disagree_on_a_row_are_withheld_and_counted(tmp_path: Path) -> None:
    """Two tools' estimates of one quantity on one row: kept where they agree, withheld where not."""
    meta = (*MINIMAC_META[:5], BEAGLE_META[4], *GT_DS)
    sites = [
        (100, "A", "G", "R2=0.9;DR2=0.9", "0|1:1"),
        (200, "C", "T", "R2=0.9;DR2=0.4", "0|1:1"),
        (300, "G", "A", "DR2=0.4", "0|1:1"),
    ]
    rows, _, warnings = by_pos(vcf(tmp_path, meta, sites, "both.vcf"))
    assert rows[100]["imputation_r2"] == one_value(sites[0][3], "R2")
    assert rows[200]["imputation_r2"] is None
    assert rows[300]["imputation_r2"] == one_value(sites[2][3], "DR2")
    assert any(line.startswith("1 row(s) write INFO R2 and DR2") for line in warnings), warnings


def test_a_value_that_is_no_number_is_left_out_and_counted(tmp_path: Path) -> None:
    meta = tuple(line.replace("ID=R2,Number=1,Type=Float", "ID=R2,Number=1,Type=String") for line in MINIMAC_META)
    sites = [(100, "A", "G", "R2=0.9", "0|1:1"), (200, "C", "T", "R2=high", "0|1:1")]
    rows, _, warnings = by_pos(vcf(tmp_path, meta, sites, "text.vcf"))
    assert rows[100]["imputation_r2"] == one_value(sites[0][3], "R2")
    assert rows[200]["imputation_r2"] is None
    assert any(line.startswith("1 INFO R2 value(s) not a number") for line in warnings), warnings


@pytest.mark.parametrize("key", ["R2", "ER2", "DR2"])
def test_a_key_the_rows_use_undeclared_is_read(tmp_path: Path, key: str) -> None:
    meta = ("##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621>", *GT_DS)
    sites = [(100, "A", "G", f"{key}=0.8", "0|1:1")]
    rows, _, _ = by_pos(vcf(tmp_path, meta, sites, f"undeclared-{key}.vcf"))
    column = "imputation_empirical_r2" if key == "ER2" else "imputation_r2"
    assert rows[100][column] == one_value(sites[0][3], key)
