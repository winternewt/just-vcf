"""somalier as an optional subprocess (RM274): skipped cleanly without its binary; its pairs.tsv read.

The `pairs.tsv` here is synthetic, in the column layout somalier's `relate.nim` writes; the run itself
is tested only where a `somalier` binary is on PATH (none on this machine: docs/TRUST.md lists no row
for it yet).
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf.relatedness import RELATEDNESS_FRAME, TOO_FEW_SITES, RelationBand, SitePanel
from just_vcf.somalier import SOMALIER_PAIRS_COLUMNS, SomalierNotFoundError, read_somalier_pairs, somalier_binary, somalier_relate

PAIRS = [
    # sample_a, sample_b, relatedness, ibs0, ibs2, concordance, hets_a, hets_b, hets_ab, shared_hets,
    # hom_alts_a, hom_alts_b, shared_hom_alts, n, x_ibs0, x_ibs2, expected_relatedness
    ("f0", "f1", 0.985, 2, 9100, 0.99, 4000, 4010, 4000, 3990, 2000, 2004, 1995, 9200, 0, 120, -1),
    ("f0", "f2", 0.512, 1, 6000, 0.65, 4000, 3900, 3900, 2600, 2000, 1950, 1200, 9150, 3, 90, -1),
    ("f1", "f2", 0.03, 700, 3500, 0.38, 4010, 3900, 3900, 1100, 2004, 1950, 600, 9100, 20, 60, -1),
    ("f0", "f3", 0.9, 0, 300, 0.99, 100, 110, 100, 95, 50, 55, 50, 310, 0, 5, -1),
]


def write_pairs(path: Path) -> Path:
    lines = ["#" + "\t".join(SOMALIER_PAIRS_COLUMNS), *("\t".join(str(value) for value in row) for row in PAIRS)]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_pairs_tsv_reads_into_the_relatedness_frame(tmp_path: Path) -> None:
    frame = read_somalier_pairs(write_pairs(tmp_path / "related.pairs.tsv"), {"f0": "facet-a", "f1": "facet-b"}, min_sites=500)
    assert frame.schema == RELATEDNESS_FRAME
    assert frame.select("facet_a", "facet_b").rows() == [("facet-a", "facet-b"), ("facet-a", "f2"), ("facet-b", "f2"), ("facet-a", "f3")]
    assert frame["somalier_relatedness"].to_list() == [row[2] for row in PAIRS]
    # somalier's relatedness is twice KING's kinship scale.
    assert frame["kinship"].to_list()[:3] == pytest.approx([row[2] / 2 for row in PAIRS[:3]])
    assert frame["band"].to_list()[:3] == [RelationBand.SAME_PERSON.value, RelationBand.FIRST_DEGREE.value, RelationBand.UNRELATED.value]
    assert frame["n_sites"].to_list() == [row[13] for row in PAIRS]
    assert frame["n_het_both"].to_list() == [row[9] for row in PAIRS]
    assert frame["ibs0"].to_list()[:3] == pytest.approx([row[3] / row[13] for row in PAIRS[:3]])
    last = frame.row(3, named=True)
    assert (last["band"], last["reason"], last["kinship"]) == (RelationBand.UNKNOWN.value, TOO_FEW_SITES, None)
    assert set(frame["method"]) == {"somalier"}


def test_without_a_binary_it_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", "")
    assert somalier_binary() is None
    panel = SitePanel(
        assembly_label="GRCh37",
        sites=pl.DataFrame({"chrom": ["1"], "pos": pl.Series([10], dtype=pl.UInt32), "ref_allele": ["A"], "effect_allele": ["G"]}),
        source="synthetic",
    )
    with pytest.raises(SomalierNotFoundError):
        somalier_relate([], panel=panel, sites_vcf=Path("sites.vcf.gz"), fasta=Path("ref.fa"))
