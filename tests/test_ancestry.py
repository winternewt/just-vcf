"""Population inference (RM5) on synthetic panels built here with a fixed seed.

Three populations share base allele frequencies and each redraws a third of the sites, so they are
well separated. Every person is drawn from one population's frequencies, written as a real VCF or
23andMe file and read through `load`, so the expected answer is known from the fixture itself.
"""

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

np = pytest.importorskip("numpy")

from just_vcf import FACETS, Loaded, load  # noqa: E402
from just_vcf.ancestry import (  # noqa: E402
    ANCESTRY_REASONS,
    CANONICAL_SUPERPOPS,
    AncestryEstimate,
    build_pca_panel,
    consensus,
    continental_of,
    infer_ancestry,
    build_prive_reference,
    load_pca_panel,
    panel_assemblies,
    prive_proportions,
)
from just_vcf.ancestry._projection import DEFAULT_MISSING, build_reference_pca  # noqa: E402
from just_vcf.ancestry.pca import PANEL_META_KEYS, panel_paths  # noqa: E402
from just_vcf.ancestry.sites import facet_of, panel_sites, site_dosages  # noqa: E402
from just_vcf.assembly import ASSEMBLY_LENGTHS  # noqa: E402
from just_vcf.schema import BLOCKS, conform  # noqa: E402

SEED = 20261022
POPULATIONS = {"AFR": "Yoruba", "EAS": "Han", "EUR": "French"}
N_PER_POPULATION = 150
SITES_PER_CONTIG = 500
CONTIGS = ("1", "2", "3")
ALLELE_PAIRS = (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"), ("A", "C"), ("C", "A"), ("G", "T"), ("T", "G"))
GT_TEXT = {0: "0/0", 1: "0/1", 2: "1/1"}


def _length(assembly: str, contig: str) -> int:
    return next(iter(ASSEMBLY_LENGTHS[assembly][contig]))


def synthetic_panel(seed: int) -> tuple[pl.DataFrame, dict[str, np.ndarray], np.ndarray, pl.DataFrame]:
    """Sites, per-population ALT frequencies, the reference genotype matrix and its labels."""
    rng = np.random.default_rng(seed)
    n = SITES_PER_CONTIG * len(CONTIGS)
    pairs = [ALLELE_PAIRS[i] for i in rng.integers(0, len(ALLELE_PAIRS), n)]
    sites = pl.DataFrame(
        {
            "chrom": [c for c in CONTIGS for _ in range(SITES_PER_CONTIG)],
            "pos": [10_000 + 1_000 * i for _ in CONTIGS for i in range(SITES_PER_CONTIG)],
            "ref": [r for r, _ in pairs],
            "alt": [a for _, a in pairs],
        }
    )
    base = rng.uniform(0.05, 0.5, n)
    freqs = {}
    for group in POPULATIONS:
        freq = base.copy()
        redrawn = rng.choice(n, n // 3, replace=False)
        freq[redrawn] = rng.uniform(0.05, 0.95, redrawn.size)
        freqs[group] = freq
    genotypes = np.hstack(
        [rng.binomial(2, freqs[group][:, None] * np.ones((1, N_PER_POPULATION))) for group in POPULATIONS]
    ).astype(np.int8)
    labels = pl.DataFrame(
        {
            "iid": [f"{group}{i}" for group in POPULATIONS for i in range(N_PER_POPULATION)],
            "superpop": [group for group in POPULATIONS for _ in range(N_PER_POPULATION)],
            "population": [POPULATIONS[group] for group in POPULATIONS for _ in range(N_PER_POPULATION)],
        }
    )
    return sites, freqs, genotypes, labels


@pytest.fixture(scope="module")
def panel(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, pl.DataFrame, dict[str, np.ndarray], float]:
    """Panel `synth` on GRCh37: its directory, sites, frequencies and leave-one-out accuracy."""
    panel_dir = tmp_path_factory.mktemp("panels")
    sites, freqs, genotypes, labels = synthetic_panel(SEED)
    built = build_pca_panel(genotypes, sites, labels, panel="synth", assembly="GRCh37", panel_dir=panel_dir)
    return panel_dir, sites, freqs, built.loo_accuracy


def person(freqs: dict[str, np.ndarray], group: str, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).binomial(2, freqs[group]).astype(np.int8)


def write_vcf(path: Path, sites: pl.DataFrame, dosage: np.ndarray, *, contigs: bool = True, variant_only: bool = False) -> Path:
    meta = ["##fileformat=VCFv4.2"]
    if contigs:
        meta += [f"##contig=<ID={c},length={_length('GRCh37', c)}>" for c in CONTIGS]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines = [*meta, "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])]
    for (chrom, pos, ref, alt), d in zip(sites.iter_rows(), dosage.tolist(), strict=True):
        if variant_only and d == 0:
            continue
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{GT_TEXT[d]}")
    path.write_text("\n".join(lines) + "\n")
    return path


def with_absence(loaded: Loaded, absence: str, blocks: pl.LazyFrame | None = None) -> Loaded:
    """The same facet with another `absence_means` (the loader gives a 1,500-site VCF `unknown`)."""
    facets = loaded.facets.with_columns(pl.lit(absence).cast(FACETS["absence_means"]).alias("absence_means"))
    return Loaded(facets=facets, calls=loaded.calls, blocks=loaded.blocks if blocks is None else blocks)


def test_build_separates_the_three_populations(panel) -> None:
    panel_dir, sites, _, loo = panel
    assert loo >= 0.95
    assert panel_assemblies(panel_dir, "synth") == ["GRCh37"]


@pytest.mark.parametrize(
    ("extra_meta", "match"),
    [({"n_variants": 1}, "would replace"), ({"source": {"files": {"a", "b"}}}, "JSON")],
    ids=["clash", "not_json"],
)
def test_a_refused_rebuild_leaves_the_previous_panel_whole(tmp_path: Path, extra_meta: dict, match: str) -> None:
    """`extra_meta` is checked before any file is written: a rebuild it refuses, from other data,
    leaves the three files of the panel already there as they were, and no temporary file beside them."""
    sites, _, genotypes, labels = synthetic_panel(SEED)
    panel_dir = tmp_path / "panels"
    build_pca_panel(genotypes, sites, labels, panel="synth", assembly="GRCh37", panel_dir=panel_dir)
    before = {path.name: path.read_bytes() for path in sorted(panel_dir.iterdir())}
    half = sites.height // 2
    with pytest.raises(ValueError, match=match):
        build_pca_panel(genotypes[:half], sites.head(half), labels, panel="synth", assembly="GRCh37", panel_dir=panel_dir, extra_meta=extra_meta)
    assert {path.name: path.read_bytes() for path in sorted(panel_dir.iterdir())} == before
    assert load_pca_panel(panel_dir, "synth", "GRCh37").sites.height == sites.height
    written = json.loads(panel_paths(panel_dir, "synth", "GRCh37")["meta"].read_text())
    assert set(written) == PANEL_META_KEYS, "the keys extra_meta is checked against are the keys written"


@pytest.mark.parametrize("group", sorted(POPULATIONS))
def test_each_population_is_placed_from_a_full_vcf(panel, tmp_path: Path, group: str) -> None:
    panel_dir, sites, freqs, _ = panel
    dosage = person(freqs, group, seed=1)
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, dosage))
    result = infer_ancestry(loaded, panel_dir, "synth")
    assert result and result.reason is None
    assert result.estimate == group
    assert result.probabilities[group] >= 0.9 and result.margin >= 0.8
    assert max(result.proportions, key=result.proportions.get) == group
    assert sum(result.proportions.values()) == pytest.approx(1.0, abs=0.01)
    assert result.proportion_method == "pca_nnls"
    assert (result.n_sites_panel, result.n_sites_observed, result.n_sites_used) == (sites.height, sites.height, sites.height)
    assert result.n_sites_absent_as_reference == 0 and result.n_sites_unusable == 0
    assert result.coverage == 1.0 and result.panel_assembly == result.facet_assembly == "GRCh37"
    assert result.group_to_continental == {g: g for g in POPULATIONS}
    assert infer_ancestry(loaded, panel_dir, "synth") == result


def test_population_resolution(panel, tmp_path: Path) -> None:
    panel_dir, sites, freqs, _ = panel
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, person(freqs, "EAS", seed=2)))
    result = infer_ancestry(loaded, panel_dir, "synth", resolution="population")
    assert result.estimate == "Han"
    assert result.group_to_continental == {"French": "EUR", "Han": "EAS", "Yoruba": "AFR"}
    with pytest.raises(ValueError, match="resolution"):
        infer_ancestry(loaded, panel_dir, "synth", resolution="country")


def test_variant_only_vcf_depends_on_what_absence_means(panel, tmp_path: Path) -> None:
    """Hom-ref sites left out: read as reference only when the facet says absence means that."""
    panel_dir, sites, freqs, _ = panel
    dosage = person(freqs, "AFR", seed=3)
    n_carrier = int((dosage > 0).sum())
    loaded = load(write_vcf(tmp_path / "v.vcf", sites, dosage, variant_only=True))
    assert facet_of(loaded).absence_means == "unknown"

    as_loaded = infer_ancestry(loaded, panel_dir, "synth")
    assert (as_loaded.n_sites_observed, as_loaded.n_sites_absent_as_reference) == (n_carrier, 0)
    assert as_loaded.coverage == pytest.approx(n_carrier / sites.height)

    genome = infer_ancestry(with_absence(loaded, "weak_reference"), panel_dir, "synth")
    assert genome.estimate == "AFR"
    assert (genome.n_sites_observed, genome.n_sites_absent_as_reference) == (n_carrier, sites.height - n_carrier)
    assert genome.coverage == 1.0
    full = infer_ancestry(load(write_vcf(tmp_path / "f.vcf", sites, dosage)), panel_dir, "synth")
    assert genome.pc_coords == pytest.approx(full.pc_coords)


def test_too_few_sites_is_unknown_not_a_guess(panel, tmp_path: Path) -> None:
    panel_dir, sites, freqs, _ = panel
    keep = 200
    loaded = load(write_vcf(tmp_path / "few.vcf", sites.head(keep), person(freqs, "EUR", seed=4)[:keep]))
    result = infer_ancestry(loaded, panel_dir, "synth")
    assert not result and result.estimate is None and result.reason == "too_few_sites"
    assert result.n_sites_used == keep and result.coverage == pytest.approx(keep / sites.height)
    assert result.probabilities == {} and result.proportions is None
    assert infer_ancestry(loaded, panel_dir, "synth", coverage_floor=0.1)


def test_array_export_is_read_from_its_letters(panel, tmp_path: Path) -> None:
    """A 23andMe file writes no REF; its letters are matched against the panel's two alleles."""
    panel_dir, sites, freqs, _ = panel
    dosage = person(freqs, "EUR", seed=5)
    no_call = set(range(0, sites.height, 50))
    lines = [
        "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018",
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for i, ((chrom, pos, ref, alt), d) in enumerate(zip(sites.iter_rows(), dosage.tolist(), strict=True)):
        letters = "--" if i in no_call else "".join(sorted(ref * (2 - d) + alt * d))
        lines.append(f"rs{i + 1}\t{chrom}\t{pos}\t{letters}")
    path = tmp_path / "genome.txt"
    path.write_text("\n".join(lines) + "\n")
    loaded = load(path)
    assert facet_of(loaded).absence_means == "unknown"
    result = infer_ancestry(loaded, panel_dir, "synth")
    assert result.estimate == "EUR"
    assert result.n_sites_observed == sites.height - len(no_call)
    assert result.n_sites_unusable == len(no_call) and result.n_sites_absent_as_reference == 0


def _rules_vcf(path: Path) -> Path:
    meta = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=1,length={_length('GRCh37', '1')}>",
        f"##contig=<ID=2,length={_length('GRCh37', '2')}>",
        f"##contig=<ID=X,length={_length('GRCh37', 'X')}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]),
    ]
    records = [
        "1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1",
        "1\t200\t.\tA\tC,G\t.\t.\t.\tGT\t0/2",
        "1\t300\t.\tA\tC,G\t.\t.\t.\tGT\t1/2",
        "1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/.",
        "1\t500\t.\tA\tG\t.\t.\t.\tGT\t./.",
        "1\t600\t.\tA\tG\t.\t.\t.\tGT\t1",
        "1\t700\t.\tA\tG\t.\t.\t.\tGT\t0/1",
        "1\t700\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
        "1\t900\t.\tA\tG\t.\t.\t.\tGT\t1/0",
        "1\t900\t.\tA\tT\t.\t.\t.\tGT\t0/1",
        "1\t1000\t.\tA\tG\t.\t.\t.\tGT\t1/1",
        "1\t1200\t.\tA\t.\t.\t.\t.\tGT\t0/0",
        "1\t1300\t.\tA\tG\t.\t.\t.\tGT\t0/1/.",
    ]
    path.write_text("\n".join([*meta, *records]) + "\n")
    return path


RULE_SITES = pl.DataFrame(
    {
        "chrom": ["1"] * 9 + ["X", "1", "1", "2", "1"],
        "pos": [100, 200, 300, 400, 500, 600, 700, 800, 900, 100, 1000, 1200, 100, 1300],
        "ref_allele": ["A"] * 14,
        "effect_allele": ["G"] * 14,
    }
)
M = DEFAULT_MISSING
# Per site: 0/1; 0/2 over C,G is A/G; 1/2 is C/G (another allele); half-call; no-call; haploid;
# an SNV beside an indel at one position; absent; a split pair (slot pattern); 1/1; absent on X;
# a monomorphic hom-ref record; absent on chromosome 2; triploid.
RULE_DOSAGE = [1, 1, M, M, M, M, 1, M, M, M, 2, 0, M, M]


def test_site_rules_unknown_absence(tmp_path: Path) -> None:
    loaded = load(_rules_vcf(tmp_path / "rules.vcf"))
    facet = facet_of(loaded)
    got = site_dosages(loaded, facet, panel_sites(RULE_SITES))
    assert got.dosage.tolist() == RULE_DOSAGE
    assert (got.n_observed, got.n_absent_as_reference, got.n_unusable) == (5, 0, 6)


def test_site_rules_weak_reference_reads_autosomal_absence_as_reference(tmp_path: Path) -> None:
    loaded = with_absence(load(_rules_vcf(tmp_path / "rules.vcf")), "weak_reference")
    got = site_dosages(loaded, facet_of(loaded), panel_sites(RULE_SITES))
    expected = list(RULE_DOSAGE)
    expected[7] = 0.0  # 1:800, absent
    expected[12] = 0.0  # 2:100, absent; X:100 stays missing (not an autosome)
    assert got.dosage.tolist() == expected
    assert (got.n_observed, got.n_absent_as_reference, got.n_unusable) == (5, 2, 6)


def test_site_rules_blocks_read_absence_inside_reference_blocks_only(tmp_path: Path) -> None:
    loaded = load(_rules_vcf(tmp_path / "rules.vcf"))
    facet_id = facet_of(loaded).facet_id
    blocks = conform(
        pl.LazyFrame(
            {
                "facet_id": [facet_id] * 4,
                "chrom": ["1", "1", "2", "X"],
                "pos": [750, 1500, 1, 50],
                "end_pos": [850, 1600, 99, 150],
                "kind": ["ref", "ref", "ref", "no_call"],
            },
            schema_overrides={"pos": pl.UInt32, "end_pos": pl.UInt32},
        ),
        BLOCKS,
    )
    gvcf = with_absence(loaded, "blocks", blocks)
    got = site_dosages(gvcf, facet_of(gvcf), panel_sites(RULE_SITES))
    expected = list(RULE_DOSAGE)
    expected[7] = 0.0  # 1:800 inside 1:750-850; 2:100 is past 2:1-99, X:100 is in a no-call block
    assert got.dosage.tolist() == expected
    assert got.n_absent_as_reference == 1


def test_site_rules_reference_reads_absence_outside_any_block_as_reference(tmp_path: Path) -> None:
    """`reference` (CGI's VCF): an omitted autosomal site was called reference; a site a no-call block
    covers was not omitted. 1:800 lies in a long block that a shorter one starts inside."""
    loaded = load(_rules_vcf(tmp_path / "rules.vcf"))
    facet_id = facet_of(loaded).facet_id
    blocks = conform(
        pl.LazyFrame(
            {
                "facet_id": [facet_id] * 3,
                "chrom": ["1", "1", "X"],
                "pos": [700, 750, 50],
                "end_pos": [900, 760, 150],
                "kind": ["no_call", "no_call", "no_call"],
            },
            schema_overrides={"pos": pl.UInt32, "end_pos": pl.UInt32},
        ),
        BLOCKS,
    )
    cgi = with_absence(loaded, "reference", blocks)
    got = site_dosages(cgi, facet_of(cgi), panel_sites(RULE_SITES))
    expected = list(RULE_DOSAGE)
    expected[12] = 0.0  # 2:100, omitted; 1:800 is no-called; X:100 is not an autosome
    assert got.dosage.tolist() == expected
    assert got.n_absent_as_reference == 1


def test_assembly_and_panel_gates(panel, tmp_path: Path) -> None:
    panel_dir, sites, freqs, _ = panel
    dosage = person(freqs, "EAS", seed=6)

    unlabelled = load(write_vcf(tmp_path / "u.vcf", sites, dosage, contigs=False))
    assert facet_of(unlabelled).assembly_label is None
    result = infer_ancestry(unlabelled, panel_dir, "synth")
    assert not result and result.reason == "assembly_unknown" and result.n_sites_used is None

    grch38_dir = tmp_path / "grch38"
    small_sites, _, genotypes, labels = synthetic_panel(SEED + 1)
    build_pca_panel(genotypes, small_sites, labels, panel="synth", assembly="GRCh38", panel_dir=grch38_dir)
    labelled = load(write_vcf(tmp_path / "l.vcf", sites, dosage))
    mismatch = infer_ancestry(labelled, grch38_dir, "synth")
    assert not mismatch and mismatch.reason == "assembly_mismatch" and "GRCh38" in mismatch.detail

    for missing_dir, name in ((tmp_path / "empty", "synth"), (panel_dir, "other")):
        missing = infer_ancestry(labelled, missing_dir, name)
        assert not missing and missing.reason == "panel_missing"


def test_estimate_holds_an_answer_or_a_reason() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        AncestryEstimate(method="pca_knn", panel="p")
    with pytest.raises(ValueError, match="exactly one"):
        AncestryEstimate(method="pca_knn", panel="p", estimate="EUR", reason="too_few_sites")
    with pytest.raises(ValueError, match="unknown ancestry reason"):
        AncestryEstimate(method="pca_knn", panel="p", reason="no_idea")
    assert "too_few_sites" in ANCESTRY_REASONS


def test_consensus_of_two_panels(panel, tmp_path: Path) -> None:
    panel_dir, sites, freqs, _ = panel
    other_sites, _, genotypes, labels = synthetic_panel(SEED + 2)
    build_pca_panel(genotypes, other_sites, labels, panel="synth2", assembly="GRCh37", panel_dir=panel_dir)
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, person(freqs, "EAS", seed=7)))
    first = infer_ancestry(loaded, panel_dir, "synth")
    second = infer_ancestry(loaded, panel_dir, "synth2")
    unknown = infer_ancestry(loaded, tmp_path / "empty", "synth")

    fused = consensus([first, second, unknown])
    assert fused.estimate == "EAS" and fused.method == "consensus"
    assert fused.panel == "synth+synth2" and fused.facet_ids == first.facet_ids
    assert set(fused.probabilities) == set(CANONICAL_SUPERPOPS)
    assert sum(fused.probabilities.values()) == pytest.approx(1.0, abs=1e-3)
    assert fused == consensus([first, second])
    assert fused.probabilities["EAS"] >= first.probabilities["EAS"] - 1e-9

    nothing = consensus([unknown])
    assert not nothing and nothing.reason == "no_estimates" and "panel_missing" in nothing.detail


def test_consensus_agreement_sharpens_and_disagreement_flattens() -> None:
    def knn(dist: dict[str, float]) -> AncestryEstimate:
        return AncestryEstimate(method="pca_knn", panel="p", estimate=max(dist, key=dist.get), probabilities=dist)

    agree = consensus([knn({"EUR": 0.9, "AMR": 0.1}), knn({"EUR": 0.85, "SAS": 0.15}), knn({"EUR": 1.0})])
    assert agree.estimate == "EUR" and agree.probabilities["EUR"] > 0.95
    split = consensus([knn({"EUR": 0.9, "EAS": 0.1}), knn({"EAS": 0.9, "EUR": 0.1})])
    assert max(split.probabilities.values()) < 0.95 and split.margin < agree.margin


def test_continental_rollup() -> None:
    prive_groups = {
        "AFR": ("Africa (West)", "Africa (South)", "Africa (East)", "Africa (North)"),
        "EUR": (
            "Ashkenazi", "Italy", "Europe (South East)", "Europe (North East)", "Finland", "Scandinavia",
            "United Kingdom", "Ireland", "Europe (South West)",
        ),
        "SAS": ("Sri Lanka", "Pakistan", "Bangladesh"),
        "EAS": ("Asia (East)", "Japan", "Philippines"),
        "AMR": ("South America",),
    }
    for continental, groups in prive_groups.items():
        assert {continental_of(group) for group in groups} == {continental}
    assert continental_of("Middle East") is None and continental_of("MID") is None
    assert continental_of("CSA") == "SAS" and continental_of("x") is None


def test_core_imports_without_numpy() -> None:
    """`just_vcf` loads with numpy unavailable; `just_vcf.ancestry` says which extra it needs."""
    script = (
        "import sys\n"
        "sys.modules['numpy'] = None\n"  # what an environment without numpy looks like to an import
        "import just_vcf\n"
        "try:\n"
        "    import just_vcf.ancestry\n"
        "except ImportError as error:\n"
        "    print(error)\n"
    )
    done = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True, env=os.environ.copy())
    assert "ancestry" in done.stdout and "just-vcf[ancestry]" in done.stdout


PRIVE_NAMES = {"AFR": "Africa (West)", "EAS": "Asia (East)", "EUR": "Europe (North East)"}


@pytest.fixture(scope="module")
def prive(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, pl.DataFrame, dict[str, np.ndarray]]:
    """A Privé-shaped reference from the synthetic panel: group frequencies and 16 PCA loadings.

    The loadings file is gzip-compressed and the frequency file is not, as both kinds occur.
    """
    root = tmp_path_factory.mktemp("prive")
    sites, freqs, genotypes, _ = synthetic_panel(SEED)
    loadings = build_reference_pca(genotypes, dim_ref=8).U[:, :16]
    key = sites.select(
        pl.col("chrom").cast(pl.Int64).alias("chr"), "pos", pl.col("ref").alias("a0"), pl.col("alt").alias("a1"),
        rsid=pl.format("rs{}", pl.int_range(1, sites.height + 1)),
    )
    key.with_columns(**{PRIVE_NAMES[g]: pl.Series(freqs[g]) for g in POPULATIONS}).write_csv(root / "ref_freqs.csv")
    projection = key.with_columns(**{f"PC{i + 1}": pl.Series(loadings[:, i]) for i in range(16)})
    with gzip.open(root / "projection.csv.gz", "wt") as handle:
        handle.write(projection.write_csv())
    built = build_prive_reference(root / "ref_freqs.csv", root / "projection.csv.gz", root / "reference")
    assert built.n_sites == sites.height and built.groups == tuple(PRIVE_NAMES[g] for g in POPULATIONS)
    return root / "reference", sites, freqs


@pytest.mark.parametrize("group", sorted(POPULATIONS))
def test_prive_proportions_of_each_group(prive, tmp_path: Path, group: str) -> None:
    reference, sites, freqs = prive
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, person(freqs, group, seed=11)))
    result = prive_proportions(loaded, reference)
    assert result and result.method == "prive_qp" and result.proportion_method == "prive_qp"
    assert result.estimate == PRIVE_NAMES[group]
    assert result.proportions[PRIVE_NAMES[group]] >= 0.9
    assert sum(result.proportions.values()) == pytest.approx(1.0, abs=0.01)
    assert set(result.proportion_errors) == set(PRIVE_NAMES.values())
    assert all(0 < e < 0.1 for e in result.proportion_errors.values())
    assert result.n_sites_used == result.n_sites_observed == sites.height
    assert result.group_to_continental == {name: g for g, name in PRIVE_NAMES.items()}
    assert prive_proportions(loaded, reference) == result


def test_a_built_reference_keeps_the_frequency_files_order(tmp_path: Path) -> None:
    """The bootstrap's blocks are runs of the reference's rows, so the rows must come in one order on
    every build: the frequency file's. An order the join chose made `proportion_errors` move between
    builds of the same files, and EUR's error for a group it holds none of sat near 0 (the flaky
    `test_prive_proportions_of_each_group[EUR]`, most likely). Built several times, since one build can fall in
    order by chance."""
    sites, freqs, genotypes, _ = synthetic_panel(SEED)
    loadings = build_reference_pca(genotypes, dim_ref=8).U[:, :16]
    key = sites.select(pl.col("chrom").cast(pl.Int64).alias("chr"), "pos", pl.col("ref").alias("a0"), pl.col("alt").alias("a1"))
    key.with_columns(**{PRIVE_NAMES[g]: pl.Series(freqs[g]) for g in POPULATIONS}).write_csv(tmp_path / "ref_freqs.csv")
    # The loading file in another order than the frequency file, as two files from two sources can be.
    key.with_columns(**{f"PC{i + 1}": pl.Series(loadings[:, i]) for i in range(16)}).reverse().write_csv(tmp_path / "projection.csv")
    expected = key.select(pl.col("chr").cast(pl.String), pl.col("pos").cast(pl.UInt32))
    for build in range(5):
        built = build_prive_reference(tmp_path / "ref_freqs.csv", tmp_path / "projection.csv", tmp_path / f"reference{build}")
        assert pl.read_parquet(built.reference_path, columns=["chr", "pos"]).equals(expected), build


def test_prive_proportions_recover_a_two_way_mixture(prive, tmp_path: Path) -> None:
    """Each allele drawn from an even mix of two groups' frequencies: about half each, within the errors."""
    reference, sites, freqs = prive
    mixed = 0.5 * freqs["AFR"] + 0.5 * freqs["EUR"]
    dosage = np.random.default_rng(12).binomial(2, mixed).astype(np.int8)
    result = prive_proportions(load(write_vcf(tmp_path / "m.vcf", sites, dosage)), reference)
    for g in ("AFR", "EUR"):
        name = PRIVE_NAMES[g]
        assert abs(result.proportions[name] - 0.5) <= max(0.1, 3 * result.proportion_errors[name])
    assert result.proportions[PRIVE_NAMES["EAS"]] <= 0.1
    assert result.margin < 0.3


def test_prive_gates_and_floor(prive, tmp_path: Path) -> None:
    reference, sites, freqs = prive
    dosage = person(freqs, "EAS", seed=13)
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, dosage))
    assert prive_proportions(loaded, tmp_path / "nothing").reason == "panel_missing"
    unlabelled = load(write_vcf(tmp_path / "u.vcf", sites, dosage, contigs=False))
    assert prive_proportions(unlabelled, reference).reason == "assembly_unknown"
    few = prive_proportions(loaded, reference, min_sites=sites.height + 1)
    assert not few and few.reason == "too_few_sites" and few.n_sites_used == sites.height
    plain = prive_proportions(loaded, reference, n_bootstrap=0)
    assert plain.proportion_errors is None
    assert plain.proportions == prive_proportions(loaded, reference).proportions


def test_consensus_reads_prive_groups_through_their_continents(panel, prive, tmp_path: Path) -> None:
    panel_dir, sites, freqs, _ = panel
    reference, _, _ = prive
    loaded = load(write_vcf(tmp_path / "p.vcf", sites, person(freqs, "EUR", seed=14)))
    fused = consensus([infer_ancestry(loaded, panel_dir, "synth"), prive_proportions(loaded, reference)])
    assert fused.estimate == "EUR" and fused.panel == "prive+synth"
