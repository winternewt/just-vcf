"""`site_genotypes.read_site_genotypes`, the panel-site reading ancestry and relatedness share, against a
facet that carries no allele of a panel site's class (RM302).

The facets are small VCFs written here and read through `load`, which counts their `variant_classes`
from their ALT alleles; the panel is two SNV sites on contig 1. The finer rules (no-calls, split rows,
absence on the autosomes) are pinned in `tests/test_ancestry.py` through `site_dosages`.
"""

from pathlib import Path

import polars as pl

from just_vcf import FACETS, Loaded, ReferenceStore, load, resolve_ref
from just_vcf.assembly import ASSEMBLY_LENGTHS
from just_vcf.site_genotypes import SiteBasis, facet_of, panel_sites, read_site_genotypes
from test_reference import SEQ, base, homref_table, other, write_fasta

PANEL = panel_sites(pl.DataFrame({"chrom": ["1", "1"], "pos": [100, 800], "ref_allele": ["A", "A"], "effect_allele": ["G", "G"]}))
"""Site 0 at 1:100, where the facets below write a row; site 1 at 1:800, where none does."""

HOM_REF_INSERTION = "1\t100\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0"
"""A hom-ref insertion row at the panel's SNV site: its haplotypes are the base `A` twice, which read as
zero copies of the panel's `G` wherever the row is read."""

INDELS = ["1\t300\t.\tC\tCT\t.\tPASS\t.\tGT\t0/1", "1\t500\t.\tGTT\tG\t.\tPASS\t.\tGT\t1/1"]
SNV_ELSEWHERE = "1\t700\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"


def vcf(path: Path, records: list[str]) -> Loaded:
    length = next(iter(ASSEMBLY_LENGTHS["GRCh37"]["1"]))
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={length}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]),
        *records,
    ]
    path.write_text("\n".join(lines) + "\n")
    return load(path)


def classes(loaded: Loaded) -> set[str] | None:
    counted = loaded.facets.collect()["variant_classes"][0]
    return None if counted is None else {entry["variant_class"] for entry in counted if entry["alleles"]}


def with_columns(loaded: Loaded, **values: object) -> Loaded:
    facets = loaded.facets.with_columns(*[pl.lit(value).cast(FACETS[name]).alias(name) for name, value in values.items()])
    return Loaded(facets=facets, calls=loaded.calls, blocks=loaded.blocks)


def read(loaded: Loaded) -> list[tuple[int | None, str]]:
    got = read_site_genotypes(loaded, facet_of(loaded), PANEL)
    return list(zip(got["count"].to_list(), got["basis"].to_list(), strict=True))


def test_an_indel_only_facet_row_says_nothing_of_a_panel_snv(tmp_path: Path) -> None:
    """RM302. A facet that wrote no SNV (its `variant_classes` hold only indels) has a hom-ref insertion
    row at a panel SNV's position. The row's haplotypes are `A`/`A`, both panel alleles, so it read as
    zero copies of `G`, `observed`; but the facet carries no SNV, so its rows say no more of the panel's
    SNV than its silence does (`callset.class_gate`, RM299's rule): unknown, `class_not_carried`. Read
    as a whole genome too, its silence at 1:800 stays unknown, as before."""
    indels = vcf(tmp_path / "indels.vcf", [HOM_REF_INSERTION, *INDELS])
    assert classes(indels) == {"indel"}
    assert read(indels) == [(None, SiteBasis.CLASS_NOT_CARRIED.value), (None, SiteBasis.ABSENT.value)]
    genome = with_columns(indels, callset_type="wgs", absence_means="weak_reference")
    assert read(genome) == [(None, SiteBasis.CLASS_NOT_CARRIED.value), (None, SiteBasis.ABSENT.value)]


def test_a_class_cut_by_the_content_filter_is_not_read_even_from_a_residue(tmp_path: Path) -> None:
    """A genome cut to its indels (`.indel.` names the cut; the rows are too few to settle it) keeps an
    SNV residue, a het at the panel site, so `snv` is among its `variant_classes`; the filter alone says
    its rows speak for no SNV (`class_cut_by_filter`), as `absence_is_reference` reads its silence."""
    cut = vcf(tmp_path / "cut.indel.vcf", ["1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", *INDELS])
    row = cut.facets.collect().row(0, named=True)
    assert row["content_filter"] == "indel_only"
    assert classes(cut) == {"snv", "indel"}
    assert read(cut) == [(None, SiteBasis.CLASS_NOT_CARRIED.value), (None, SiteBasis.ABSENT.value)]


def test_a_facet_that_carries_snvs_or_was_not_counted_is_read_as_before(tmp_path: Path) -> None:
    """The control: a facet that carries SNVs (one elsewhere, no filter) has its hom-ref insertion row at
    the panel site read as zero copies, `observed`. A facet whose classes were not counted (null
    `variant_classes`) is unknown, not a facet that lacks SNVs, and its row is read the same."""
    mixed = vcf(tmp_path / "mixed.vcf", [HOM_REF_INSERTION, *INDELS, SNV_ELSEWHERE])
    assert classes(mixed) == {"snv", "indel"}
    assert read(mixed) == [(0, SiteBasis.OBSERVED.value), (None, SiteBasis.ABSENT.value)]
    uncounted = with_columns(vcf(tmp_path / "indels.vcf", [HOM_REF_INSERTION, *INDELS]), variant_classes=None)
    assert classes(uncounted) is None
    assert read(uncounted) == [(0, SiteBasis.OBSERVED.value), (None, SiteBasis.ABSENT.value)]


def test_a_facet_that_wrote_no_variant_allele_has_its_rows_read(tmp_path: Path) -> None:
    """Audit 18. A 23andMe-style table whose every row is hom-ref counts no class once `resolve_ref` reads
    REF (`variant_classes` is `[]`). That list says nothing of which classes the table carries, so its
    hom-ref row at a panel SNV is read, zero copies of the panel's ALT, `observed`: it is not a facet
    that carries no SNV. Its silence at the other site is an array's, unknown, as before."""
    store = ReferenceStore.open(tmp_path / "grch37")
    store.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    table = resolve_ref(load(homref_table(tmp_path / "genome.txt")), store)
    assert classes(table) == set()
    ref = {pos: base("1", pos) for pos in (10, 30)}
    panel = panel_sites(pl.DataFrame({"chrom": ["1", "1"], "pos": [10, 30], "ref_allele": [ref[10], ref[30]], "effect_allele": [other(ref[10]), other(ref[30])]}))
    got = read_site_genotypes(table, facet_of(table), panel)
    assert list(zip(got["count"].to_list(), got["basis"].to_list(), strict=True)) == [(0, SiteBasis.OBSERVED.value), (None, SiteBasis.ABSENT.value)]
