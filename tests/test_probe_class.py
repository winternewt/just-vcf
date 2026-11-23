"""An array's row speaks only for its own probe's class (RM303): concordance, consensus and panel reads.

A whole-genome caller's rows are gated by the classes the whole facet carries (`callset.class_gate`,
RM299 to RM302). An array is a set of probes, each for its own allele, so an array facet's row is read
only for alleles of the row's own class (`callset.row_class_expr`). The arrays here are AncestryDNA
exports whose `I`/`D` rows `resolve_ref` reads against `test_probes`' synthetic manifest: `rs100` probes
a deletion at 501, `rs200` an insertion at 1001, both on contig 1 of `test_reference.SEQ`. Every
expected allele is spelled from that contig here; the whole-genome facets are VCFs over the same
sequence, written to spell each indel as the probe's row spells it.
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, resolve_ref
from just_vcf.concordance import concordance
from just_vcf.consensus import consensus
from just_vcf.loading import Loaded
from just_vcf.site_genotypes import facet_of, panel_sites, read_site_genotypes
from test_consensus import E, TRIANGLE, hard, posterior, rows_of
from test_probes import CONTIG, INSERTED, coded_export, grch37_store, manifests  # noqa: F401 (fixtures)
from test_reference import SEQ

DELETION_AT, INSERTION_AT = 501, 1001
"""Where the manifest's two indel probes are spelled (their `MapInfo`, the exports' positions)."""
DELETION_REF, DELETION_ALT = CONTIG[500:504], CONTIG[503]
"""The deletion probe's site as `place_indel` spells it at 501: the deleted bases and the base after."""
INSERTION_REF, INSERTION_ALT = CONTIG[1000], INSERTED + CONTIG[1000]
"""The insertion probe's site at 1001: the base after the gap, with the inserted bases before it."""
ANCHOR = 1300
"""An SNV probe's position, where the whole genome calls an insertion anchored on the probe's base."""


def other(b: str) -> str:
    return next(x for x in "ACGT" if x != b)


def snv_gt(pos: int, copies: int) -> str:
    """An export's two letters at `pos` holding `copies` of `other(REF)`."""
    ref = CONTIG[pos - 1]
    return ref * (2 - copies) + other(ref) * copies


def chip(tmp_path: Path, store: ReferenceStore, rows: list[tuple[str, str, int, str]], name: str = "chip") -> Loaded:
    """An AncestryDNA export (an array, `callset_type = array`) with REF read and its codes resolved."""
    directory = tmp_path / name
    directory.mkdir()
    resolved = resolve_ref(load(coded_export(directory, rows)), store)
    assert resolved.facets.collect()["callset_type"][0] == "array"
    assert resolved.calls.collect()["ref"].null_count() == 0, "every row read"
    return resolved


def write_vcf(path: Path, records: list[tuple[int, str, str, str]]) -> Path:
    """A one-sample VCF on contig 1 of `SEQ`: `(pos, ref, alt, gt)` records, every contig declared."""
    meta = ["##fileformat=VCFv4.3", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items())]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}" for pos, ref, alt, gt in sorted(records)]
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def genome(tmp_path: Path, store: ReferenceStore, name: str, records: list[tuple[int, str, str, str]]) -> Loaded:
    """A VCF facet whose callset type is not an array (too few rows to name one: null), REF read."""
    resolved = resolve_ref(load(write_vcf(tmp_path / f"{name}.vcf", records)), store)
    assert resolved.facets.collect()["callset_type"][0] is None
    return resolved


def snv(pos: int, gt: str) -> tuple[int, str, str, str]:
    ref = CONTIG[pos - 1]
    return (pos, ref, other(ref), gt)


def insertion_at_anchor(gt: str) -> tuple[int, str, str, str]:
    ref = CONTIG[ANCHOR - 1]
    return (ANCHOR, ref, ref + other(CONTIG[ANCHOR]) * 2, gt)


def classes_of(loaded: Loaded) -> set[str]:
    return {entry["variant_class"] for entry in loaded.facets.collect()["variant_classes"][0]}


def one(sites: pl.DataFrame, pos: int, alt: str) -> dict:
    (row,) = sites.filter((pl.col("pos") == pos) & (pl.col("alt") == alt)).iter_rows(named=True)
    return row


def labelled(loaded: Loaded, caller: str) -> Loaded:
    """The facet with its caller named, as the consensus reads a facet's method from the atlas."""
    schema = loaded.facets.collect_schema()
    values = {"caller": caller, "caller_version": "1"}
    return Loaded(
        facets=loaded.facets.with_columns(*[pl.lit(value, dtype=schema[key]).alias(key) for key, value in values.items()]),
        calls=loaded.calls,
        blocks=loaded.blocks,
        links=loaded.links,
    )


def at(result, pos: int, alt: str) -> dict:
    (row,) = [row for row in rows_of(result) if row["pos"] == pos and row["alts"] == [alt]]
    return row


def silent_reasons(result, loaded: Loaded) -> str:
    facet_id = loaded.facets.collect()["facet_id"][0]
    (line,) = [line for line in result.warnings if line.startswith(f"facet {facet_id} silent on")]
    return line


# --- 1. an array facet's row: concordance ----------------------------------------------------------


def test_a_chips_snv_row_says_nothing_of_an_insertion_though_the_chip_carries_an_indel(
    tmp_path: Path, grch37_store: ReferenceStore, manifests: Path
) -> None:
    """A chip with one resolved non-REF indel counts `indel` among its classes, so the facet-level gate
    read its hom-ref SNV probe row at the insertion's anchor as zero copies of the insertion. The row is
    an SNV probe and says nothing of an insertion: the allele stays one-sided (`class_not_carried`)."""
    array = chip(tmp_path, grch37_store, [("rs100", "1", DELETION_AT, "DI"), ("rs1", "1", ANCHOR, snv_gt(ANCHOR, 0)), ("rs2", "1", 1700, snv_gt(1700, 1))])
    assert classes_of(array) == {"snv", "indel"}
    ins = insertion_at_anchor("0/1")
    wgs = genome(tmp_path, grch37_store, "wgs", [ins, snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")])
    result = concordance(wgs, array)
    lone = one(result.sites.collect(), ANCHOR, ins[2])
    assert (lone["variant_class"], lone["key"], lone["agreement"], lone["absence"], lone["a_copies"], lone["b_copies"]) == (
        "indel", None, "only_a", "class_not_carried", 1, None,
    )  # fmt: skip
    assert any("indel alleles of facet a" in line and "class_not_carried" in line for line in result.warnings), result.warnings


def test_a_chips_indel_probe_row_is_read_for_its_indel_and_not_for_an_snv(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    """A chip whose `I`/`D` rows are all hom-ref counts no indel, so the facet-level gate withheld its
    deletion probe row from the deletion it probes. The probe looked and saw none: zero copies, a
    mismatch against the genome's het call (key `site`). Its hom-ref insertion probe row has a one-base
    REF and no ALT, an SNV row's shape, but `source_gt` keeps its codes: it is an indel probe, and says
    nothing of an SNV at that base."""
    array = chip(
        tmp_path, grch37_store,
        [("rs100", "1", DELETION_AT, "II"), ("rs200", "1", INSERTION_AT, "DD"), ("rs1", "1", ANCHOR, snv_gt(ANCHOR, 1))],
    )  # fmt: skip
    assert classes_of(array) == {"snv"}
    calls = {row["pos"]: row for row in array.calls.collect().iter_rows(named=True)}
    assert (calls[DELETION_AT]["ref"], calls[DELETION_AT]["alts"]) == (DELETION_REF, [])
    assert (calls[INSERTION_AT]["ref"], calls[INSERTION_AT]["alts"], calls[INSERTION_AT]["source_gt"]) == (INSERTION_REF, [], "DD")
    snv_there = snv(INSERTION_AT, "0/1")
    wgs = genome(tmp_path, grch37_store, "wgs", [(DELETION_AT, DELETION_REF, DELETION_ALT, "0/1"), snv_there, snv(ANCHOR, "0/1"), insertion_at_anchor("0/1")])
    sites = concordance(wgs, array).sites.collect()
    read = one(sites, DELETION_AT, DELETION_ALT)
    assert (read["key"], read["agreement"], read["absence"], read["a_copies"], read["b_copies"]) == ("site", "mismatch", "other_alleles", 1, 0)
    lone = one(sites, INSERTION_AT, snv_there[2])
    assert (lone["variant_class"], lone["key"], lone["agreement"], lone["absence"], lone["b_copies"]) == ("snv", None, "only_a", "class_not_carried", None)


def test_a_genomes_rows_keep_the_facet_level_gate(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """The control: the same rows written by a facet that is not an array (a VCF of unknown callset
    type) are read as before. It carries an indel, so its hom-ref SNV row at the insertion's anchor is
    read for the insertion (zero copies, a mismatch); one carrying no indel has it withheld."""
    ins = insertion_at_anchor("0/1")
    wgs = genome(tmp_path, grch37_store, "wgs", [ins, snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")])
    carries = genome(tmp_path, grch37_store, "carries", [snv(ANCHOR, "0/0"), snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")])
    read = one(concordance(wgs, carries).sites.collect(), ANCHOR, ins[2])
    assert (read["key"], read["agreement"], read["absence"], read["a_copies"], read["b_copies"]) == ("site", "mismatch", "other_alleles", 1, 0)
    lacks = genome(tmp_path, grch37_store, "lacks", [snv(ANCHOR, "0/0"), snv(1700, "0/1")])
    lone = one(concordance(wgs, lacks).sites.collect(), ANCHOR, ins[2])
    assert (lone["key"], lone["agreement"], lone["absence"]) == (None, "only_a", "class_not_carried")


# --- 1. an array facet's row: consensus ------------------------------------------------------------


def test_the_consensus_reads_a_chips_row_only_for_its_own_class(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    """RM303 in the consensus. A chip with a resolved het deletion carries `indel`; its hom-ref SNV probe
    row at the insertion's anchor is silent on the insertion (`class_not_carried`), never a vote of zero
    copies. Its hom-ref deletion probe row of another chip, which counts no indel, votes zero copies of
    the deletion it probes."""
    ins = insertion_at_anchor("0/1")
    wgs = labelled(genome(tmp_path, grch37_store, "wgs", [ins, snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")]), "X")
    with_indel = labelled(chip(tmp_path, grch37_store, [("rs100", "1", DELETION_AT, "DI"), ("rs1", "1", ANCHOR, snv_gt(ANCHOR, 0)), ("rs2", "1", 1700, snv_gt(1700, 1))]), "Y")
    assert classes_of(with_indel) == {"snv", "indel"}
    result = consensus([wgs, with_indel], atlas=TRIANGLE)
    lone = at(result, ANCHOR, ins[2])
    assert (lone["s_variant_class"], lone["s_supported"], lone["s_contradicted"], lone["s_silent"]) == ("indel", 1, 0, 1)
    assert lone["gp"] == pytest.approx(posterior([hard(1, E["X"])]), rel=1e-5)
    assert "class_not_carried 1" in silent_reasons(result, with_indel)

    hom_ref = labelled(chip(tmp_path, grch37_store, [("rs100", "1", DELETION_AT, "II"), ("rs1", "1", ANCHOR, snv_gt(ANCHOR, 1))], name="hom_ref"), "Y")
    assert classes_of(hom_ref) == {"snv"}
    result = consensus([wgs, hom_ref], atlas=TRIANGLE)
    read = at(result, DELETION_AT, DELETION_ALT)
    assert (read["s_supported"], read["s_contradicted"], read["s_silent"]) == (1, 1, 0)
    assert read["gp"] == pytest.approx(posterior([hard(1, E["X"]), hard(0, E["Y"])]), rel=1e-5)


def test_the_consensus_keeps_the_facet_level_gate_for_a_genome(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """The control: a facet that is not an array and carries an indel votes zero copies of the insertion
    from its hom-ref SNV row at the anchor, as before."""
    ins = insertion_at_anchor("0/1")
    wgs = labelled(genome(tmp_path, grch37_store, "wgs", [ins, snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")]), "X")
    carries = labelled(genome(tmp_path, grch37_store, "carries", [snv(ANCHOR, "0/0"), snv(1700, "0/1"), (DELETION_AT, DELETION_REF, DELETION_ALT, "0/1")]), "Y")
    read = at(consensus([wgs, carries], atlas=TRIANGLE), ANCHOR, ins[2])
    assert (read["s_supported"], read["s_contradicted"], read["s_silent"]) == (1, 1, 0)


# --- 1. an array facet's row: a panel's sites ------------------------------------------------------


def _panel(*positions: int) -> pl.DataFrame:
    return panel_sites(
        pl.DataFrame(
            {
                "chrom": ["1"] * len(positions),
                "pos": list(positions),
                "ref_allele": [CONTIG[p - 1] for p in positions],
                "effect_allele": [other(CONTIG[p - 1]) for p in positions],
            }
        )
    )


def test_a_chips_insertion_probe_row_is_not_read_at_a_panel_snv(tmp_path: Path, grch37_store: ReferenceStore, manifests: Path) -> None:
    """A hom-ref insertion probe row (`DD`) at a panel SNV's position holds REF twice, both panel alleles,
    and was read as zero copies of the panel's ALT. It probes an insertion: `class_not_carried`. The
    chip's SNV probe row at another panel site is read as before."""
    array = chip(tmp_path, grch37_store, [("rs200", "1", INSERTION_AT, "DD"), ("rs1", "1", ANCHOR, snv_gt(ANCHOR, 1))])
    read = read_site_genotypes(array, facet_of(array), _panel(INSERTION_AT, ANCHOR))
    assert read.select("count", "basis").rows() == [(None, "class_not_carried"), (1, "observed")]


def test_a_genomes_hom_ref_insertion_row_is_read_at_a_panel_snv_as_before(tmp_path: Path, grch37_store: ReferenceStore) -> None:
    """The control: a facet that is not an array and carries SNVs has its hom-ref insertion row at a
    panel SNV read as zero copies, as before (the facet-level gate)."""
    wgs = genome(tmp_path, grch37_store, "wgs", [(INSERTION_AT, INSERTION_REF, INSERTION_ALT, "0/0"), snv(ANCHOR, "0/1")])
    read = read_site_genotypes(wgs, facet_of(wgs), _panel(INSERTION_AT, ANCHOR))
    assert read.select("count", "basis").rows() == [(0, "observed"), (1, "observed")]
