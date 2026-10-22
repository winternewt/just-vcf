"""Population inference against real reference panels, opt-in: `JUST_VCF_ANCESTRY_DIR`.

The directory may hold Privé's reference as `build_prive_reference` writes it, and any number of
reference-PCA panels in just-prs's layout. Each test skips when its kind of panel is not there. The
persons are synthetic, drawn from the panel itself, so the expected answer comes from the panel.
"""

import json
from pathlib import Path

import polars as pl
import pytest

np = pytest.importorskip("numpy")

from just_vcf import load  # noqa: E402
from just_vcf.ancestry import continental_of, infer_ancestry, load_pca_panel, prive_proportions, to_continental  # noqa: E402
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


def test_pca_panels_place_their_own_reference_samples(ancestry_dir: Path, tmp_path: Path) -> None:
    """Each panel's first reference sample of every group, rebuilt from the panel's own low-rank
    factors (U·diag(s)·Vᵀ, unstandardised and rounded to 0/1/2), is placed in its own group."""
    metas = sorted(ancestry_dir.glob("ancestry_model_*_meta.json"))
    if not metas:
        pytest.skip(f"no reference-PCA panel under {ancestry_dir}")
    for meta_path in metas:
        meta = json.loads(meta_path.read_text())
        panel, assembly = meta["panel"], meta["build"]
        model = load_pca_panel(ancestry_dir, panel, assembly)
        raw = pl.read_parquet(panel_paths(ancestry_dir, panel, assembly)["sites"], columns=["chrom", "pos", "ref", "alt"])
        superpops = model.superpops.tolist()
        for group in sorted(set(superpops)):
            i = superpops.index(group)
            standardised = model.pca.U @ (model.pca.s * model.pca.V[i])
            dosage = np.clip(np.rint(standardised * model.pca.std.ravel() + model.pca.mean.ravel()), 0, 2).astype(int)
            rows = list(zip(raw["chrom"], raw["pos"], raw["ref"], raw["alt"], dosage.tolist(), strict=True))
            loaded = load(write_vcf(tmp_path / f"{panel}-{group}.vcf", assembly, rows))
            result = infer_ancestry(loaded, ancestry_dir, panel)
            assert result.estimate == group, (panel, group, result.reason, result.probabilities)
