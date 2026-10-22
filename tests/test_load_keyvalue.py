"""PGP-HMS's freeBayes+VarScan2 key=value table (`freebayes_varscan2_keyvalue`): one facet per caller.

The synthetic rows copy the shapes counted in the one PGP-HMS table on disk (2026-10-22, counts only):
2,044,829 rows, every one with an annotation tail and the twelve keys in order, no position twice;
1,745,306 with both callers, 226,646 freeBayes only, 72,877 VarScan2 only; 68,901 freeBayes `0/0`;
62,483 multi-allelic ALTs; genotypes such as `2/1` and VarScan2 `0/3`.
"""

import gzip
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import load, sniff
from just_vcf.schema import CALLS, canonical_chrom
from test_load_vcf import assert_gt_consistent

TAIL = "KG=0.3196\tAFR_MAF=0.70\tEFF_RESULTS\tEFF_HGVS=NM_002849:c.471+1175G>T\tVEP_RESULTS"
FIELDS = ("GenoType", "RefAlleleCount", "ReferenceAlleleQ", "TotalDepth")


def row(chrom: str, pos: int, ref: str, alt: str, fb: tuple, vs: tuple, tail: str = TAIL) -> str:
    keyed = [f"{p}{k}={v}" for p, values in (("FB", fb), ("VS", vs)) for k, v in zip(FIELDS, values, strict=True)]
    return "\t".join([chrom, str(pos), ref, alt, *keyed, *([tail] if tail else [])])


NULLS = ("NULL", "NULL", "NULL", "NULL")
ROWS = [
    row("chr1", 752566, "G", "A", ("0/1", "12", "475", "33"), ("0/1", "12", "39", "31")),
    row("chr1", 752721, "A", "G", ("1/1", ".", "NULL", "27"), ("1/1", "NULL", "NULL", "27"), tail="KG=0.1"),
    row("chr1", 846808, "C", "T", ("0/1", "23", "869", "46"), NULLS),
    row("chr2", 1000, "C", "A,T,G", ("2/1", "1", "8", "30"), ("0/3", "NULL", "NULL", "24")),
    row("chr2", 2000, "ACG", "A", NULLS, ("0/1", "NULL", "NULL", "18")),
    row("chr2", 3000, "T", "C", ("0/0", "20", "700", "20"), NULLS, tail=""),
    row("chrX", 5000, "G", "T", ("1/1", ".", "NULL", "9"), ("1/1", "NULL", "NULL", "9")),
    row("chrM", 150, "C", "T", ("1/1", ".", "NULL", "900"), ("1/1", "NULL", "NULL", "880")),
]
MALFORMED = "chr3\t10\tA\tG\tFBGenoType=0/1\tFBTotalDepth=9"


def expected(prefix: str) -> list[dict]:
    """The rows of `ROWS` a caller reported, read back from the text by the test's own split."""
    out = []
    for line in ROWS:
        fields = line.split("\t")
        chrom, pos, ref, alt = fields[:4]
        offset = 4 if prefix == "FB" else 8
        values = [field.split("=", 1)[1] for field in fields[offset : offset + 4]]
        if values[0] == "NULL":
            continue
        alleles = [ref, *alt.split(",")]
        out.append(
            {
                "chrom": canonical_chrom(chrom),
                "pos": int(pos),
                "source_ref": ref,
                "alts": alt.split(","),
                "gt_alleles": [alleles[int(i)] for i in values[0].split("/")],
                "source_gt": values[0],
                "dp": int(values[3]),
            }
        )
    return sorted(out, key=lambda r: ({"1": 0, "2": 1, "X": 22, "MT": 24}[r["chrom"]], r["pos"]))


def write(tmp_path: Path, lines: list[str], name: str = "AB12CD-orig.vcf.gz") -> Path:
    path = tmp_path / name
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    return path


def frames(path: Path):
    loaded = load(path)
    facets, calls, blocks = pl.collect_all([loaded.facets, loaded.calls, loaded.blocks])
    by_caller = {r["caller"]: r for r in facets.iter_rows(named=True)}
    rows = {caller: calls.filter(pl.col("facet_id") == f["facet_id"]) for caller, f in by_caller.items()}
    return loaded, facets, by_caller, rows, calls, blocks


def test_one_facet_per_caller(tmp_path: Path) -> None:
    shuffled = ROWS[:]
    random.Random(7).shuffle(shuffled)
    path = write(tmp_path, [*shuffled, MALFORMED])
    assert sniff(path).format == "freebayes_varscan2_keyvalue"
    loaded, facets, by_caller, rows, calls, blocks = frames(path)
    assert calls.schema == CALLS and blocks.height == 0
    assert facets["facet_id"].to_list() == sorted(facets["facet_id"].to_list())
    assert calls["facet_id"].to_list() == sorted(calls["facet_id"].to_list()), "calls come in facet order"
    for prefix, caller in (("FB", "freeBayes"), ("VS", "VarScan2")):
        got = rows[caller].select("chrom", "pos", "source_ref", "alts", "gt_alleles", "source_gt", "dp").to_dicts()
        assert got == expected(prefix), caller
        facet = by_caller[caller]
        assert (facet["relation"], facet["derived_from"], facet["vendor"]) == ("caller_split", None, "PGP-HMS")
        assert facet["assembly_label"] is None, "the table states no build"
        assert {"FBRefAlleleCount", "VSReferenceAlleleQ", "KG", "EFF_RESULTS", "EFF_HGVS"} <= set(facet["dropped_fields"])
        assert "FBGenoType" not in facet["dropped_fields"] and "VSTotalDepth" not in facet["dropped_fields"]
        assert facet["absence_means"] == "unknown", "eight sites are no genome"
    assert by_caller["freeBayes"]["source_sha256"] == by_caller["VarScan2"]["source_sha256"]
    assert any(w.startswith("1 row(s) whose fields do not carry the twelve fixed keys") for w in loaded.warnings)
    assert_gt_consistent(calls)


def _classes(prefix: str) -> list[dict]:
    """The ALT alleles of each class a caller's rows carry, by the test's own reading of REF and ALT:
    an ALT of another length than its REF is an indel, one as long an SNV (one base) or an MNV."""
    found = {"snv": 0, "mnv": 0, "indel": 0}
    for r in expected(prefix):
        for alt in r["alts"]:
            found["indel" if len(alt) != len(r["source_ref"]) else "snv" if len(alt) == 1 else "mnv"] += 1
    return [{"variant_class": cls, "alleles": n} for cls, n in found.items() if n]


def test_variant_classes_and_content_filter(tmp_path: Path) -> None:
    """Decision 4: a key=value facet counts its ALT alleles by class, as a VCF facet does, so absence
    can be gated per class; too few alleles settle no filter from the rows, and the name's word does."""
    _, _, by_caller, _, _, _ = frames(write(tmp_path, ROWS))
    for prefix, caller in (("FB", "freeBayes"), ("VS", "VarScan2")):
        facet = by_caller[caller]
        assert facet["variant_classes"] == _classes(prefix), caller
        assert (facet["content_filter"], facet["content_filter_basis"]) == (None, None)
    _, _, named, _, _, _ = frames(write(tmp_path, ROWS, "AB12CD-orig.indels.vcf.gz"))
    assert {(f["content_filter"], f["content_filter_basis"]) for f in named.values()} == {("indel_only", "file_name")}


def test_genotype_fields(tmp_path: Path) -> None:
    _, _, _, rows, _, _ = frames(write(tmp_path, ROWS))
    fb = {r["pos"]: r for r in rows["freeBayes"].iter_rows(named=True)}
    vs = {r["pos"]: r for r in rows["VarScan2"].iter_rows(named=True)}
    assert fb[1000]["gt"] == [2, 1] and fb[1000]["gt_alleles"] == ["T", "A"], "kept in the order written"
    assert vs[1000]["gt"] == [0, 3] and vs[1000]["gt_alleles"] == ["C", "G"]
    assert fb[3000]["gt"] == [0, 0] and fb[3000]["alts"] == ["C"], "freeBayes's 0/0 stays a call"
    assert vs[2000]["end_pos"] == 2002 and fb[752566]["end_pos"] == 752566
    assert fb[752566]["phased"] == [False, False] and fb[752566]["ref"] is None
    assert fb[752566]["ref_check"] == "not_run" and fb[752566]["gl"] is None and fb[752566]["ad"] is None
    assert 752566 in vs, "VarScan2 keeps the sites freeBayes also found, unlike a one-column merge"


def test_an_index_past_the_alts_is_counted(tmp_path: Path) -> None:
    lines = [*ROWS, row("chr2", 4000, "A", "G", ("0/2", "1", "1", "5"), NULLS)]
    loaded, _, _, rows, _, _ = frames(write(tmp_path, lines))
    bad = rows["freeBayes"].filter(pl.col("pos") == 4000).row(0, named=True)
    assert bad["gt"] == [0, 2] and bad["gt_alleles"] == ["A", None]
    assert "freeBayes facet: 1 genotype(s) index an allele past the row's ALT list; gt_alleles has a null there" in loaded.warnings


def test_loading_twice_gives_the_same_frames(tmp_path: Path) -> None:
    path = write(tmp_path, ROWS)
    first, second = load(path), load(path)
    assert first.facets.collect().equals(second.facets.collect())
    assert first.calls.collect().equals(second.calls.collect())


def test_real_keyvalue_against_its_merge(pgp_participants: Path) -> None:
    """The PGP-HMS `-orig` table beside its one-column merged VCF, counts only.

    The table is a subset of the merge (on the one pair on disk, 2.04M of 5.30M rows, every chromosome,
    the gene-dense ones most), so it is checked as one: each of its freeBayes rows is a freeBayes row
    of the merge with the same alleles; where both hold a VarScan2 row, the alleles agree; and its
    VarScan2 facet holds sites the merge lost to freeBayes."""
    tables = sorted(pgp_participants.glob("*/arvados-*/*-orig.vcf.gz"), key=lambda p: p.stat().st_size)
    pair = next(
        ((t, t.with_name(t.name.replace("-orig", ""))) for t in tables if t.with_name(t.name.replace("-orig", "")).is_file()),
        None,
    )
    if pair is None:
        pytest.skip("no PGP-HMS key=value table beside its merge on disk")
    orig, merge = load(pair[0]), load(pair[1])
    key = ["chrom", "pos", "source_ref", "gt_alleles"]

    def by_caller(loaded) -> dict[str, pl.LazyFrame]:
        facets = loaded.facets.collect()
        return {
            caller: loaded.calls.filter(pl.col("facet_id") == fid).select(*key)
            for caller, fid in zip(facets["caller"], facets["facet_id"], strict=True)
        }

    table, vcf = by_caller(orig), by_caller(merge)

    def count(frame: pl.LazyFrame) -> int:
        return frame.select(pl.len()).collect().item()

    site = key[:3]
    fb_not_in_merge = count(table["freeBayes"].join(vcf["freeBayes"], on=key, how="anti", nulls_equal=True))
    vs_both = table["VarScan2"].join(vcf["VarScan2"], on=site, how="inner", suffix="_merge")
    vs_disagree = count(vs_both.filter(pl.col("gt_alleles") != pl.col("gt_alleles_merge")))
    vs_recovered = count(table["VarScan2"].join(vcf["VarScan2"], on=site, how="anti"))
    assert fb_not_in_merge == 0 and count(table["freeBayes"]) > 0
    assert count(vs_both) > 0 and vs_disagree == 0
    assert vs_recovered > 0, "the table keeps VarScan2's genotype where both callers found the site"
