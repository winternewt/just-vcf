"""somalier as an optional subprocess (RM274): skipped cleanly without its binary; its pairs.tsv read.

The `pairs.tsv` here is synthetic, in the column layout somalier's `relate.nim` writes; the run itself
is tested only where a `somalier` binary is on PATH (none on this machine: docs/TRUST.md lists no row
for it yet). What the wrapper hands somalier (its commands, the VCF of each facet) is read where the
binary would have run.
"""

import gzip
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, resolve_ref
from just_vcf.relatedness import RELATEDNESS_FRAME, TOO_FEW_SITES, RelationBand, SitePanel
from just_vcf.site_genotypes import facet_of as site_facet
from just_vcf.site_genotypes import read_site_genotypes
from just_vcf.somalier import SOMALIER_PAIRS_COLUMNS, SomalierNotFoundError, read_somalier_pairs, somalier_binary, somalier_relate
from test_unstated_ref import base, store, vcf, with_facet  # noqa: F401 (store is a fixture)

PAIRS = [
    # sample_a, sample_b, relatedness, ibs0, ibs2, concordance, hets_a, hets_b, hets_ab, shared_hets,
    # hom_alts_a, hom_alts_b, shared_hom_alts, n, x_ibs0, x_ibs2, expected_relatedness
    ("f0", "f1", 0.985, 2, 9100, 0.99, 4000, 4010, 4000, 3990, 2000, 2004, 1995, 9200, 0, 120, -1),
    ("f0", "f2", 0.512, 1, 6000, 0.65, 4000, 3900, 3900, 2600, 2000, 1950, 1200, 9150, 3, 90, -1),
    ("f1", "f2", 0.03, 700, 3500, 0.38, 4010, 3900, 3900, 1100, 2004, 1950, 600, 9100, 20, 60, -1),
    ("f0", "f3", 0.9, 0, 300, 0.99, 100, 110, 100, 95, 50, 55, 50, 310, 0, 5, -1),
    ("f2", "f3", 0.49, 180, 5600, 0.62, 3900, 3950, 3900, 2400, 1950, 1980, 1100, 9000, 4, 80, -1),
]


def write_pairs(path: Path) -> Path:
    lines = ["#" + "\t".join(SOMALIER_PAIRS_COLUMNS), *("\t".join(str(value) for value in row) for row in PAIRS)]
    path.write_text("\n".join(lines) + "\n")
    return path


def transition(b: str) -> str:
    """The transition partner of `b`: a site panel holds no strand-ambiguous (A/T, C/G) pair."""
    return {"A": "G", "G": "A", "C": "T", "T": "C"}[b]


def test_pairs_tsv_reads_into_the_relatedness_frame(tmp_path: Path) -> None:
    frame = read_somalier_pairs(write_pairs(tmp_path / "related.pairs.tsv"), {"f0": "facet-a", "f1": "facet-b"}, min_sites=500)
    assert frame.schema == RELATEDNESS_FRAME
    pairs = [("facet-a", "facet-b"), ("facet-a", "f2"), ("facet-b", "f2"), ("facet-a", "f3"), ("f2", "f3")]
    assert frame.select("facet_a", "facet_b").rows() == pairs
    assert frame["somalier_relatedness"].to_list() == [row[2] for row in PAIRS]
    # somalier's relatedness is twice KING's kinship scale.
    assert frame["kinship"].to_list()[:3] == pytest.approx([row[2] / 2 for row in PAIRS[:3]])
    # First degree is split by IBS0 as KING-robust's is: 1 of 9,150 sites is parent and child, 180 of
    # 9,000 (0.02, past KING's 0.008) full siblings.
    assert frame["band"].to_list()[:3] == ["same_person", "parent_child", "unrelated"]
    assert frame["band"][4] == "full_siblings"
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


def test_somalier_reads_a_withheld_site_as_unknown_and_the_facets_silence_as_reference(
    tmp_path: Path,
    store: ReferenceStore,  # noqa: F811 (the fixture)
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """somalier is given each facet as `site_genotypes` reads it, never a silence to impute.

    `--unknown` is a flag of `somalier relate` (`unk2hr`), not of `extract`, and reads every unknown of
    every sample in the run as homozygous reference. A site whose only rows `write_vcf(resolved_only=True)`
    withheld (`N` → `X`, `ref_unstated`) had no record, so it read as zero copies where `site_genotypes`
    reads it unknown; so did a no-call (`./.`, unusable), under the same flag. The somalier-bound file
    now holds a record only where that reading gives a count: the whole-genome facet's silence (900) as
    `0/0`, no record at the withheld and unusable sites (300, 400, 700), and somalier runs without
    `--unknown`. The binary is not run: the commands and the files handed to it are read where it would
    have run."""
    records = [
        (300, "unstated_hom_alt", "N", transition(base(300)), "1/1"),
        (400, "unstated_het", "N", transition(base(400)), "0/1"),
        (600, "stated_hom_alt", base(600), transition(base(600)), "1/1"),
        (700, "no_call", base(700), transition(base(700)), "./."),
        (800, "beside_stated_unstated", "N", transition(base(800)), "1/1"),
        (800, "beside_stated_stated", base(800), transition(base(800)), "0/1"),
    ]
    positions = [300, 400, 600, 700, 800, 900]
    loaded = with_facet(resolve_ref(load(vcf(tmp_path, records, "sites.vcf")), store), absence_means="weak_reference")
    panel = SitePanel(
        assembly_label="GRCh37",
        sites=pl.DataFrame({"chrom": ["1"] * len(positions), "pos": pl.Series(positions, dtype=pl.UInt32),
                            "ref_allele": [base(p) for p in positions], "effect_allele": [transition(base(p)) for p in positions]}),
        source="synthetic",
    )  # fmt: skip
    commands: list[list[str]] = []
    handed: dict[str, list[tuple[int, str, str, str, str]]] = {}

    def run(command: list[str], **_: object) -> None:
        commands.append(command)
        if command[1] == "extract":
            with gzip.open(command[-1], "rt") as text:
                written = [line.rstrip("\n").split("\t") for line in text if not line.startswith("#")]
            handed[Path(command[-1]).name] = [(int(r[1]), r[3], r[4], r[8], r[9]) for r in written]
        else:
            write_pairs(Path(f"{command[command.index('-o') + 1]}.pairs.tsv"))

    monkeypatch.setattr("just_vcf.somalier.somalier_binary", lambda: Path("/nonexistent/somalier"))
    monkeypatch.setattr("just_vcf.somalier.subprocess.run", run)
    somalier_relate([loaded], panel=panel, sites_vcf=tmp_path / "sites.vcf.gz", fasta=tmp_path / "ref.fa", workdir=tmp_path)

    assert not any(flag in command for command in commands for flag in ("--unknown", "-u"))
    # GT alone: where a header declares AD, somalier reads a record without one as depth zero, unknown.
    assert handed == {
        "f0.vcf.gz": [
            (600, base(600), transition(base(600)), "GT", "1/1"),
            (800, base(800), transition(base(800)), "GT", "0/1"),
            (900, base(900), transition(base(900)), "GT", "0/0"),
        ]
    }
    read = {positions[r["site"]]: r["count"] for r in read_site_genotypes(loaded, site_facet(loaded), panel.sites).iter_rows(named=True)}
    assert {pos for pos, count in read.items() if count is not None} == {pos for pos, *_ in handed["f0.vcf.gz"]}
