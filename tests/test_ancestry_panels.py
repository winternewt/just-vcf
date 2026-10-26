"""Population inference against real reference panels, opt-in: `JUST_VCF_ANCESTRY_DIR`.

The directory may hold Privé's reference as `build_prive_reference` writes it, and any number of
reference-PCA panels in just-prs's layout. Each test skips when its kind of panel is not there. The
persons are synthetic, drawn from the panel itself, so the expected answer comes from the panel. One
test runs always: the own-sample check's rule, on a synthetic panel and the same panel scrambled.
"""

import json
import shutil
from pathlib import Path

import polars as pl
import pytest

np = pytest.importorskip("numpy")

from just_vcf import load  # noqa: E402
from just_vcf.ancestry import (  # noqa: E402
    AncestryEstimate,
    build_pca_panel,
    continental_of,
    infer_ancestry,
    load_pca_panel,
    prive_proportions,
    to_continental,
)
from just_vcf.ancestry.pca import panel_paths  # noqa: E402
from just_vcf.ancestry.prive import META_FILE, REFERENCE_FILE  # noqa: E402
from just_vcf.assembly import ASSEMBLY_LENGTHS  # noqa: E402
from just_vcf.schema import canonical_chrom  # noqa: E402

SEED = 20261022
GT_TEXT = {0: "0/0", 1: "0/1", 2: "1/1"}


def write_vcf(path: Path, assembly: str, rows: list[tuple[str, int, str, str, int]]) -> Path:
    """A one-sample VCF of `(chrom, pos, ref, alt, dosage)` rows, contig lengths from `assembly`."""
    lengths = ASSEMBLY_LENGTHS[assembly]
    contigs = sorted({chrom for chrom, *_ in rows}, key=lambda c: (len(canonical_chrom(c)), canonical_chrom(c)))
    meta = ["##fileformat=VCFv4.2"]
    meta += [f"##contig=<ID={c},length={next(iter(lengths[canonical_chrom(c)]))}>" for c in contigs if canonical_chrom(c) in lengths]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    body = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{GT_TEXT[d]}" for c, p, r, a, d in rows]
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


def test_prive_places_people_drawn_from_its_groups(ancestry_dir: Path, tmp_path: Path) -> None:
    """A person drawn from one group's frequencies at 50,000 sites lands on that group's continent."""
    if not (ancestry_dir / REFERENCE_FILE).exists():
        pytest.skip(f"no {REFERENCE_FILE} under {ancestry_dir}")
    meta = json.loads((ancestry_dir / META_FILE).read_text())
    groups = ["Africa (West)", "Japan", "United Kingdom"]
    sites = (
        pl.scan_parquet(ancestry_dir / REFERENCE_FILE)
        .select("chr", "pos", "a0", "a1", *groups)
        .filter(pl.col("chr").is_in([str(n) for n in range(1, 23)]))
        .collect()
        .unique(subset=["chr", "pos"], keep="none")
        .sample(n=50_000, seed=SEED)
        .sort(pl.col("chr").cast(pl.Int32), "pos")
    )
    rng = np.random.default_rng(SEED)
    for group in groups:
        dosage = rng.binomial(2, sites[group].to_numpy())
        rows = list(zip(sites["chr"], sites["pos"], sites["a0"], sites["a1"], dosage.tolist(), strict=True))
        loaded = load(write_vcf(tmp_path / f"{group}.vcf", meta["build"], rows))
        result = prive_proportions(loaded, ancestry_dir)
        assert result, result.detail
        assert result.n_sites_used == sites.height
        continental = to_continental(result.proportions, result.group_to_continental)
        assert continental_of(result.estimate) == continental_of(group)
        assert continental[continental_of(group)] >= 0.8, (group, result.proportions)


BORDER_OWN_SHARE = 1 / 3
"""On a border, the sample's own group holds at least this share of its neighbours' votes."""
BORDER_TWO_WAY_SHARE = 0.9
"""On a border, its own group and the one it was placed in together hold at least this share."""


def placed_in_own_region_or_on_its_border(group: str, result: AncestryEstimate) -> bool:
    """Whether a reference sample of `group` is placed where its panel says it belongs.

    In its own group, or on the border between its own group and one other: the vote is two-way (its
    own group and the leader together hold at least `BORDER_TWO_WAY_SHARE` of the neighbours) and its
    own group holds at least `BORDER_OWN_SHARE`. Groups that meet along a cline split a sample's
    neighbours so (aadr_ho, 2026-10-26: a Karakalpak 0.35 Central Asian, 0.60 European; a Xibo 0.40
    Eastern Asian, 0.60 Central Asian). A sample landing among unrelated groups spreads its vote over
    several, or leaves its own group a small share, and fails.
    """
    if result.estimate == group:
        return True
    own = result.probabilities.get(group, 0.0)
    leader = result.probabilities.get(result.estimate, 0.0)
    return own >= BORDER_OWN_SHARE and own + leader >= BORDER_TWO_WAY_SHARE


def place_own_samples(panel_dir: Path, panel: str, assembly: str, tmp_path: Path) -> dict[str, AncestryEstimate]:
    """Each group's first reference sample, rebuilt from the panel's own low-rank factors
    (U·diag(s)·Vᵀ, unstandardised and rounded to 0/1/2), placed on the panel."""
    model = load_pca_panel(panel_dir, panel, assembly)
    raw = pl.read_parquet(panel_paths(panel_dir, panel, assembly)["sites"], columns=["chrom", "pos", "ref", "alt"])
    superpops = model.superpops.tolist()
    placed = {}
    for group in sorted(set(superpops)):
        i = superpops.index(group)
        standardised = model.pca.U @ (model.pca.s * model.pca.V[i])
        dosage = np.clip(np.rint(standardised * model.pca.std.ravel() + model.pca.mean.ravel()), 0, 2).astype(int)
        rows = list(zip(raw["chrom"], raw["pos"], raw["ref"], raw["alt"], dosage.tolist(), strict=True))
        loaded = load(write_vcf(tmp_path / f"{panel_dir.name}-{panel}-{group}.vcf", assembly, rows))
        placed[group] = infer_ancestry(loaded, panel_dir, panel)
    return placed


def scrambled_copy(panel_dir: Path, panel: str, assembly: str, into: Path, seed: int) -> Path:
    """The panel's three files under `into`, with its reference samples' `superpop` labels shuffled."""
    into.mkdir()
    source, target = panel_paths(panel_dir, panel, assembly), panel_paths(into, panel, assembly)
    shutil.copyfile(source["sites"], target["sites"])
    shutil.copyfile(source["meta"], target["meta"])
    refpcs = pl.read_parquet(source["refpcs"])
    shuffled = np.random.default_rng(seed).permutation(refpcs["superpop"].to_numpy())
    refpcs.with_columns(pl.Series("superpop", shuffled, dtype=pl.String)).write_parquet(target["refpcs"])
    return into


def test_the_own_sample_check_passes_a_separated_panel_and_fails_it_scrambled(tmp_path: Path) -> None:
    """The check `test_pca_panels_place_their_own_reference_samples` applies, on a synthetic panel of
    three separated groups: every group's own sample passes, and with the labels shuffled at least one
    fails, since its neighbours then vote by chance."""
    rng = np.random.default_rng(SEED)
    groups, n_per_group, n_sites = ("North", "South", "West"), 40, 600
    pairs = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")]
    picked = [pairs[i] for i in rng.integers(0, len(pairs), n_sites)]
    # Over all 22 autosomes, so the VCFs' contig lines tell their assembly.
    sites = pl.DataFrame(
        {
            "chrom": [str(1 + i % 22) for i in range(n_sites)],
            "pos": [10_000 + 1_000 * (i // 22) for i in range(n_sites)],
            "ref": [r for r, _ in picked],
            "alt": [a for _, a in picked],
        }
    ).sort(pl.col("chrom").cast(pl.Int32), "pos")
    freqs = {group: rng.uniform(0.05, 0.95, n_sites) for group in groups}
    genotypes = np.hstack([rng.binomial(2, freqs[group][:, None] * np.ones((1, n_per_group))) for group in groups]).astype(np.int8)
    labels = pl.DataFrame(
        {
            "iid": [f"{group}{i}" for group in groups for i in range(n_per_group)],
            "superpop": [group for group in groups for _ in range(n_per_group)],
            "population": [group for group in groups for _ in range(n_per_group)],
        }
    )
    panel_dir = tmp_path / "panel"
    build_pca_panel(genotypes, sites, labels, panel="synth", assembly="GRCh37", panel_dir=panel_dir, dim_ref=4)
    placed = place_own_samples(panel_dir, "synth", "GRCh37", tmp_path)
    assert all(placed_in_own_region_or_on_its_border(group, result) for group, result in placed.items()), placed
    scrambled = place_own_samples(scrambled_copy(panel_dir, "synth", "GRCh37", tmp_path / "scrambled", SEED), "synth", "GRCh37", tmp_path)
    assert all(scrambled.values()), "every scrambled placement is an estimate, so the check below judges votes"
    assert not all(placed_in_own_region_or_on_its_border(group, result) for group, result in scrambled.items()), scrambled


def test_pca_panels_place_their_own_reference_samples(ancestry_dir: Path, tmp_path: Path) -> None:
    """Each panel's first reference sample of every group, rebuilt from the panel's own factors, is
    placed in its own group or on its border with one other (`placed_in_own_region_or_on_its_border`:
    a two-way vote, its own group holding at least a third). Never a skip per panel: the same panel
    with its labels shuffled must fail that check for at least one group, so the loosening cannot pass
    a panel that places its samples among unrelated groups."""
    metas = sorted(ancestry_dir.glob("ancestry_model_*_meta.json"))
    if not metas:
        pytest.skip(f"no reference-PCA panel under {ancestry_dir}")
    for meta_path in metas:
        meta = json.loads(meta_path.read_text())
        panel, assembly = meta["panel"], meta["build"]
        placed = place_own_samples(ancestry_dir, panel, assembly, tmp_path)
        for group, result in placed.items():
            assert placed_in_own_region_or_on_its_border(group, result), (panel, group, result.reason, result.probabilities)
        scrambled_dir = scrambled_copy(ancestry_dir, panel, assembly, tmp_path / f"scrambled-{panel}-{assembly}", SEED)
        scrambled = place_own_samples(scrambled_dir, panel, assembly, tmp_path)
        assert all(scrambled.values()), (panel, "every scrambled placement is an estimate, so the check below judges votes")
        assert not all(placed_in_own_region_or_on_its_border(group, result) for group, result in scrambled.items()), (
            panel,
            {group: result.probabilities for group, result in scrambled.items()},
        )
