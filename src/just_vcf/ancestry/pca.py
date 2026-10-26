"""Reference-PCA panels: build one from a genotype matrix, read one, place a facet on it.

Ported from just-prs `ancestry/__init__.py` (docs/PROVENANCE.md). A panel is three files in one
directory, named as just-prs names them, so a panel just-prs built works here unchanged:

- `ancestry_model_<panel>_<assembly>.parquet`: the panel's sites, `chrom, pos, ref, alt,
  effect_allele, mean, std, u0..u{dim_online-1}` (the variant loadings). `ref` is the assembly's
  reference allele and the effect allele is `alt`.
- `ancestry_refpcs_<panel>_<assembly>.parquet`: the reference samples, `iid, superpop, population,
  v0..v{dim_online-1}` (the sample eigenvectors).
- `ancestry_model_<panel>_<assembly>_meta.json`: `build` (the assembly label), `s` (singular values),
  `dim_ref`, `dim_stu`, `dim_online`, `knn_k`, `loo_accuracy` and counts, plus whatever a builder
  adds (`build_aadr_panel` records its source files and selection rules); readers ignore extra keys.

Which panels exist: just-prs builds them offline, with plink2 for QC and LD pruning, from the PGS
Catalog's 1000 Genomes and HGDP+1kGP reference panels (`pgsc_1000G_v1`, `pgsc_HGDP+1kGP_v1`, GRCh38),
and from AADR Human Origins (GRCh37). That pipeline does not come along; `build_pca_panel` takes the
matrix it would produce.
"""

import json
import logging
from pathlib import Path

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.ancestry._projection import (
    DEFAULT_MISSING,
    ReferencePCA,
    build_reference_pca,
    project_samples,
)
from just_vcf.ancestry.consensus import continental_of
from just_vcf.ancestry.result import AncestryEstimate
from just_vcf.ancestry.sites import facet_of, panel_sites, site_dosages
from just_vcf.loading import Loaded

logger = logging.getLogger(__name__)

DEFAULT_KNN_K = 20
DEFAULT_DIM_REF = 10
DEFAULT_COVERAGE_FLOOR = 0.2
"""Below this fraction of panel sites read, the estimate is unknown (just-prs's value). OADP fills a
missing site with the panel mean, so a thinly covered genome is pulled toward the panel's centre."""

RESOLUTIONS: frozenset[str] = frozenset({"superpop", "population"})
"""`superpop`: the panel's broad groups. `population`: its fine populations, where it has them."""


class PanelBuild(BaseModel):
    """What `build_pca_panel` wrote."""

    model_config = ConfigDict(frozen=True)

    sites_path: Path
    refpcs_path: Path
    meta_path: Path
    loo_accuracy: float
    n_sites: int
    n_reference: int


class PcaPanel(BaseModel):
    """A reference-PCA panel read into memory."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    assembly: str
    pca: ReferencePCA
    sites: pl.DataFrame
    """`chrom, pos, ref_allele, effect_allele` in the panel's order (`sites.panel_sites`)."""
    superpops: np.ndarray
    populations: np.ndarray
    knn_k: int


def panel_paths(panel_dir: Path, panel: str, assembly: str) -> dict[str, Path]:
    """The three files of panel `panel` on `assembly` under `panel_dir`."""
    stem = f"ancestry_model_{panel}_{assembly}"
    return {
        "sites": panel_dir / f"{stem}.parquet",
        "refpcs": panel_dir / f"ancestry_refpcs_{panel}_{assembly}.parquet",
        "meta": panel_dir / f"{stem}_meta.json",
    }


def panel_assemblies(panel_dir: Path, panel: str) -> list[str]:
    """The assemblies panel `panel` has under `panel_dir`, read from its meta files, sorted."""
    if not panel_dir.is_dir():
        return []
    found = []
    for meta in sorted(panel_dir.glob(f"ancestry_model_{panel}_*_meta.json")):
        assembly = json.loads(meta.read_text())["build"]
        if meta == panel_paths(panel_dir, panel, assembly)["meta"]:
            found.append(assembly)
    return sorted(found)


def build_pca_panel(
    genotypes: np.ndarray,
    sites: pl.DataFrame,
    labels: pl.DataFrame,
    *,
    panel: str,
    assembly: str,
    panel_dir: Path,
    dim_ref: int = DEFAULT_DIM_REF,
    knn_k: int = DEFAULT_KNN_K,
    extra_meta: dict | None = None,
) -> PanelBuild:
    """Build a reference-PCA panel and write its three files under `panel_dir`.

    `genotypes` is a (sites x samples) matrix of ALT-allele counts, 0/1/2, -9 where missing; its rows
    follow `sites` (`chrom, pos, ref, alt`, `ref` the assembly's reference allele) and its columns
    follow `labels` (`iid, superpop, population`). `assembly` is the label of the assembly the
    positions are on (`GRCh37`, `GRCh38`, `T2T-CHM13`), the same vocabulary as `facets.assembly_label`.
    `extra_meta` adds keys to the meta file (a builder's record of its source and rules); it may not
    replace one this function writes.
    """
    if genotypes.shape[0] != sites.height:
        raise ValueError(f"genotypes has {genotypes.shape[0]} rows, sites {sites.height}")
    if genotypes.shape[1] != labels.height:
        raise ValueError(f"genotypes has {genotypes.shape[1]} columns, labels {labels.height}")

    model = build_reference_pca(genotypes, dim_ref=dim_ref)
    superpops = labels["superpop"].cast(pl.String).to_numpy()
    loo = loo_accuracy(model.pcs_ref, superpops, knn_k)

    panel_dir.mkdir(parents=True, exist_ok=True)
    paths = panel_paths(panel_dir, panel, assembly)
    sites.select(
        pl.col("chrom").cast(pl.String),
        pl.col("pos").cast(pl.Int64),
        pl.col("ref").cast(pl.String),
        pl.col("alt").cast(pl.String),
    ).with_columns(
        effect_allele=pl.col("alt"),
        mean=pl.Series(model.mean.flatten()),
        std=pl.Series(model.std.flatten()),
        **{f"u{j}": pl.Series(model.U[:, j]) for j in range(model.dim_online)},
    ).write_parquet(paths["sites"])
    labels.select(
        pl.col("iid").cast(pl.String),
        pl.col("superpop").cast(pl.String),
        pl.col("population").cast(pl.String),
    ).with_columns(**{f"v{j}": pl.Series(model.V[:, j]) for j in range(model.dim_online)}).write_parquet(
        paths["refpcs"]
    )
    meta = {
        "panel": panel,
        "build": assembly,
        "dim_ref": model.dim_ref,
        "dim_stu": model.dim_stu,
        "dim_online": model.dim_online,
        "s": model.s.tolist(),
        "knn_k": knn_k,
        "loo_accuracy": loo,
        "n_reference": labels.height,
        "n_variants": sites.height,
        "superpopulations": sorted({str(x) for x in superpops.tolist()}),
    }
    clash = sorted(set(meta) & set(extra_meta or {}))
    if clash:
        raise ValueError(f"extra_meta would replace {clash}")
    meta.update(extra_meta or {})
    paths["meta"].write_text(json.dumps(meta, indent=2))
    logger.info("panel %s on %s: %d sites, %d reference samples, LOO accuracy %.4f", panel, assembly, sites.height, labels.height, loo)
    return PanelBuild(
        sites_path=paths["sites"],
        refpcs_path=paths["refpcs"],
        meta_path=paths["meta"],
        loo_accuracy=loo,
        n_sites=sites.height,
        n_reference=labels.height,
    )


def load_pca_panel(panel_dir: Path, panel: str, assembly: str) -> PcaPanel:
    """Read panel `panel` on `assembly` from `panel_dir`. A missing file raises `FileNotFoundError`."""
    paths = panel_paths(panel_dir, panel, assembly)
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(f"ancestry panel file missing: {path}")
    meta = json.loads(paths["meta"].read_text())
    if meta["build"] != assembly:
        raise ValueError(f"{paths['meta']} says build {meta['build']!r}, its name says {assembly!r}")
    dim_online = int(meta["dim_online"])
    sites = pl.read_parquet(paths["sites"])
    refpcs = pl.read_parquet(paths["refpcs"])
    pca = ReferencePCA(
        mean=sites["mean"].to_numpy().reshape((-1, 1)),
        std=sites["std"].to_numpy().reshape((-1, 1)),
        s=np.asarray(meta["s"], dtype=np.float64),
        V=refpcs.select([f"v{j}" for j in range(dim_online)]).to_numpy(),
        U=sites.select([f"u{j}" for j in range(dim_online)]).to_numpy(),
        dim_ref=int(meta["dim_ref"]),
        dim_stu=int(meta["dim_stu"]),
        dim_online=dim_online,
    )
    return PcaPanel(
        name=panel,
        assembly=assembly,
        pca=pca,
        sites=panel_sites(sites.select("chrom", "pos", pl.col("ref").alias("ref_allele"), "effect_allele")),
        superpops=refpcs["superpop"].to_numpy(),
        populations=refpcs["population"].to_numpy(),
        knn_k=int(meta.get("knn_k", DEFAULT_KNN_K)),
    )


def infer_ancestry(
    loaded: Loaded,
    panel_dir: Path,
    panel: str,
    *,
    facet_id: str | None = None,
    resolution: str = "superpop",
    coverage_floor: float = DEFAULT_COVERAGE_FLOOR,
) -> AncestryEstimate:
    """Place a facet on a reference-PCA panel: the nearest group by a neighbour vote, and proportions.

    The panel is the one under `panel_dir` named `panel` on the facet's own assembly
    (`facets.assembly_label`). A facet with no assembly label, or one the panel does not have, gets an
    unknown estimate saying so; nothing is lifted or assumed. A facet whose sites cover less than
    `coverage_floor` of the panel gets `too_few_sites`.
    """
    if resolution not in RESOLUTIONS:
        raise ValueError(f"resolution {resolution!r}; known: {sorted(RESOLUTIONS)}")
    facet = facet_of(loaded, facet_id)
    unknown = {"method": "pca_knn", "panel": panel, "facet_ids": (facet.facet_id,), "facet_assembly": facet.assembly_label}
    available = panel_assemblies(panel_dir, panel)
    if not available:
        return AncestryEstimate(**unknown, reason="panel_missing", detail=f"no panel {panel!r} under {panel_dir}")
    if facet.assembly_label is None:
        return AncestryEstimate(
            **unknown, reason="assembly_unknown", detail=f"the facet states no assembly; panel {panel!r} has {available}"
        )
    if facet.assembly_label not in available:
        return AncestryEstimate(
            **unknown,
            reason="assembly_mismatch",
            detail=f"the facet is on {facet.assembly_label}; panel {panel!r} has {available}",
        )

    model = load_pca_panel(panel_dir, panel, facet.assembly_label)
    counts = site_dosages(loaded, facet, model.sites)
    n_panel = model.sites.height
    coverage = counts.n_used / n_panel if n_panel else 0.0
    measured = {
        **unknown,
        "panel_assembly": model.assembly,
        "n_sites_panel": n_panel,
        "n_sites_used": counts.n_used,
        "n_sites_observed": counts.n_observed,
        "n_sites_absent_as_reference": counts.n_absent_as_reference,
        "n_sites_unusable": counts.n_unusable,
        "coverage": coverage,
    }
    logger.info(
        "%s on %s/%s: %d of %d sites (%d observed, %d absent read as reference, %d unusable)",
        facet.facet_id, panel, model.assembly, counts.n_used, n_panel,
        counts.n_observed, counts.n_absent_as_reference, counts.n_unusable,
    )
    if counts.n_used == 0 or coverage < coverage_floor:
        return AncestryEstimate(
            **measured,
            reason="too_few_sites",
            detail=f"{counts.n_used} of {n_panel} panel sites read ({coverage:.3f}), below the floor {coverage_floor}",
        )

    pcs = project_samples(model.pca, counts.dosage.reshape((-1, 1)).copy())[0]
    labels = model.populations if resolution == "population" else model.superpops
    leader, votes = knn_vote(model.pca.pcs_ref, labels, pcs, model.knn_k)
    proportions = centroid_proportions(model.pca.pcs_ref, labels, pcs)
    rollup = _population_rollup(model) if resolution == "population" else {}
    return AncestryEstimate(
        **measured,
        estimate=leader,
        probabilities=votes,
        margin=margin(votes),
        proportions=proportions,
        proportion_method="pca_nnls",
        group_to_continental={group: continental_of(rollup.get(group, group)) for group in sorted(set(labels.tolist()))},
        pc_coords=tuple(float(x) for x in pcs),
    )


def _population_rollup(model: PcaPanel) -> dict[str, str]:
    """Each fine population's broad group, as the panel's reference samples pair them."""
    return {str(population): str(superpop) for population, superpop in zip(model.populations, model.superpops, strict=True)}


def knn_vote(ref_pcs: np.ndarray, ref_labels: np.ndarray, query: np.ndarray, k: int) -> tuple[str, dict[str, float]]:
    """Majority vote of the `k` reference samples nearest `query`: (leader, fraction of `k` per group).

    Distance ties keep the reference order (a stable sort) and vote ties go to the group that sorts
    first, so the same inputs give the same answer (charter invariant 12).
    """
    distance = np.linalg.norm(ref_pcs - query, axis=1)
    nearest = np.argsort(distance, kind="stable")[:k]
    groups, counts = np.unique(ref_labels[nearest], return_counts=True)
    votes = {str(group): float(count) / float(len(nearest)) for group, count in zip(groups, counts, strict=True)}
    return str(groups[int(np.argmax(counts))]), votes


def margin(distribution: dict[str, float]) -> float:
    """The leading value minus the runner-up (the leading value when there is no runner-up)."""
    ranked = sorted(distribution.values(), reverse=True)
    return ranked[0] - (ranked[1] if len(ranked) > 1 else 0.0)


def project_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection of `v` onto the probability simplex (non-negative, summing to 1)."""
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    ind = np.arange(1, len(v) + 1)
    cond = u - css / ind > 0
    rho = ind[cond][-1]
    theta = css[cond][-1] / rho
    return np.maximum(v - theta, 0.0)


def simplex_lstsq(M: np.ndarray, y: np.ndarray, iters: int = 800) -> np.ndarray:
    """argmin_w ||M w - y||^2 subject to w >= 0, sum(w) = 1, by projected gradient descent."""
    k = M.shape[1]
    MtM = M.T @ M
    Mty = M.T @ y
    L = float(np.linalg.norm(MtM, 2)) or 1.0
    eta = 1.0 / L
    w = np.full(k, 1.0 / k)
    for _ in range(iters):
        w = project_simplex(w - eta * 2.0 * (MtM @ w - Mty))
    return w


def centroid_proportions(ref_pcs: np.ndarray, ref_labels: np.ndarray, pcs: np.ndarray) -> dict[str, float]:
    """The sample's PC vector as a non-negative, sum-to-one mix of the groups' PC centroids.

    An admixed genome lies between clusters, so this approximates group fractions in PC space. It is
    not haplotype-level admixture.
    """
    groups = sorted({str(x) for x in ref_labels.tolist()})
    centroids = np.column_stack([ref_pcs[ref_labels == group].mean(axis=0) for group in groups])
    weights = simplex_lstsq(centroids, pcs)
    return {group: round(float(w), 4) for group, w in zip(groups, weights, strict=True)}


def loo_accuracy(ref_pcs: np.ndarray, ref_labels: np.ndarray, k: int) -> float:
    """Leave-one-out neighbour-vote accuracy on the reference samples: how well the panel separates."""
    n = ref_pcs.shape[0]
    if n <= k:
        return 0.0
    correct = 0
    for i in range(n):
        distance = np.linalg.norm(ref_pcs - ref_pcs[i], axis=1)
        distance[i] = np.inf
        nearest = np.argsort(distance, kind="stable")[:k]
        groups, counts = np.unique(ref_labels[nearest], return_counts=True)
        correct += groups[int(np.argmax(counts))] == ref_labels[i]
    return float(correct) / float(n)

