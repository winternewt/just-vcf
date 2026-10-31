"""Privé's worldwide reference: group proportions over 21 groups from published PCs and frequencies.

Ported from just-prs `ancestry/prive.py` (docs/PROVENANCE.md), which reimplements the method of the
bigsnpr `ancestry` vignette (Privé 2022, "Using the UK Biobank as a global reference of worldwide
populations", Bioinformatics); only the reference's data files are read, none of bigsnpr's GPL code.

The reference is two public files on figshare, on GRCh37:

- `ref_freqs.csv.gz` (figshare file 31620968): `chr, pos, a0, a1, rsid` and the frequency of `a1` in
  each of 21 groups (`Africa (West)`, …, `Philippines`);
- `projection.csv.gz` (figshare file 31620953): `chr, pos, a0, a1, rsid` and the loadings `PC1..PC16`.

`build_prive_reference` joins them once into `prive_reference.parquet` plus a meta file. A genome is
then projected onto the 16 PCs with the published per-PC shrinkage correction, and its proportions
are the simplex-constrained least squares of that projection against the groups' frequency
centroids in the same PC space, computed over exactly the sites the genome was read at.

`a0` is taken as the reference allele when an absence is read as homozygous reference
(`sites.site_dosages`). Checked 2026-10-22 against a PGP-HMS DRAGEN GRCh37 genome: where its records
met the reference's positions, its REF was `a0` at 2,389,449 of 2,389,956 and `a1` at 252
(docs/ANCESTRY_PANELS.md).
"""

import gzip
import json
import logging
import shutil
from pathlib import Path

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.ancestry._projection import DEFAULT_MISSING
from just_vcf.ancestry.consensus import continental_of
from just_vcf.ancestry.pca import margin, simplex_lstsq
from just_vcf.ancestry.result import AncestryEstimate
from just_vcf.ancestry.sites import facet_of, panel_sites, site_dosages
from just_vcf.loading import Loaded

logger = logging.getLogger(__name__)

PRIVE_ASSEMBLY = "GRCh37"
N_PCS = 16
CORRECTION: tuple[float, ...] = (
    1.0, 1.0, 1.0, 1.008, 1.021, 1.034, 1.052, 1.074,
    1.099, 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443,
)
"""Per-PC shrinkage correction for projecting a new genome onto these published loadings, from the
bigsnpr vignette. It belongs to these loadings; it is not recomputed."""

REFERENCE_FILE = "prive_reference.parquet"
META_FILE = "prive_reference_meta.json"
DEFAULT_MIN_SITES = 1_000
"""A guard against an estimate from almost nothing, not a calibrated threshold: the bootstrap
standard errors are what say how much the sites support the proportions."""
DEFAULT_BOOTSTRAP = 100
DEFAULT_BLOCKS = 200
_KEY = ["chr", "pos", "a0", "a1"]


class PriveBuild(BaseModel):
    """What `build_prive_reference` wrote."""

    model_config = ConfigDict(frozen=True)

    reference_path: Path
    meta_path: Path
    n_sites: int
    groups: tuple[str, ...]


def build_prive_reference(freq_csv: Path, proj_csv: Path, out_dir: Path) -> PriveBuild:
    """Join Privé's frequency and loading files into one parquet and a meta file under `out_dir`.

    Either file may be gzip-compressed. The parquet holds `chr, pos, a0, a1`, `load0..load15` and one
    frequency column per group, in the frequency file's order; a site must be in both files. Rows keep the
    frequency file's order: the bootstrap's blocks are runs of rows, so an order the join chose would make
    `proportion_errors` differ from one build of the same files to the next.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    freq = _scan_csv(freq_csv, work / "freqs.csv")
    proj = _scan_csv(proj_csv, work / "projection.csv")
    groups = [c for c in freq.collect_schema().names() if c not in (*_KEY, "rsid")]
    pcs = [c for c in proj.collect_schema().names() if c not in (*_KEY, "rsid")]
    if len(pcs) != N_PCS:
        raise ValueError(f"expected {N_PCS} loadings in {proj_csv}, found {len(pcs)}: {pcs}")
    reference_path = out_dir / REFERENCE_FILE
    (
        freq.select(*_KEY, *groups)
        .join(proj.select(*_KEY, *pcs).rename({c: f"load{i}" for i, c in enumerate(pcs)}), on=_KEY, how="inner", maintain_order="left")
        .with_columns(pl.col("chr").cast(pl.String), pl.col("pos").cast(pl.UInt32))
        .sink_parquet(reference_path)
    )
    shutil.rmtree(work)
    n_sites = pl.scan_parquet(reference_path).select(pl.len()).collect().item()
    meta_path = out_dir / META_FILE
    meta = {
        "build": PRIVE_ASSEMBLY,
        "n_pcs": N_PCS,
        "correction": list(CORRECTION),
        "groups": groups,
        "n_variants": n_sites,
    }
    meta_path.write_text(json.dumps(meta, indent=2))
    logger.info("Privé reference: %d sites in both files, %d groups", n_sites, len(groups))
    return PriveBuild(reference_path=reference_path, meta_path=meta_path, n_sites=n_sites, groups=tuple(groups))


def _scan_csv(path: Path, scratch: Path) -> pl.LazyFrame:
    """A lazy scan of `path`, decompressed into `scratch` first when it is gzip (by its magic bytes)."""
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        with gzip.open(path, "rb") as source, scratch.open("wb") as target:
            shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
        path = scratch
    return pl.scan_csv(path, schema_overrides={"chr": pl.String, "a0": pl.String, "a1": pl.String})


def prive_proportions(
    loaded: Loaded,
    reference_dir: Path,
    *,
    facet_id: str | None = None,
    min_sites: int = DEFAULT_MIN_SITES,
    n_bootstrap: int = DEFAULT_BOOTSTRAP,
    n_blocks: int = DEFAULT_BLOCKS,
    seed: int = 0,
) -> AncestryEstimate:
    """A facet's proportions over Privé's 21 groups, with block-bootstrap standard errors.

    The reference is on GRCh37; a facet on another assembly, or with none stated, is unknown with a
    reason. Sites are read as `sites.site_dosages` reads them. The projection and the group centroids
    are sums over sites, so the sites (in genome order) are cut into `n_blocks` contiguous blocks and
    resampled `n_bootstrap` times with a generator seeded by `seed`; neighbouring sites move together,
    which keeps linkage between them from shrinking the errors. `n_bootstrap=0` skips it.
    """
    facet = facet_of(loaded, facet_id)
    unknown = {"method": "prive_qp", "panel": "prive", "facet_ids": (facet.facet_id,), "facet_assembly": facet.assembly_label}
    reference_path = reference_dir / REFERENCE_FILE
    meta_path = reference_dir / META_FILE
    if not (reference_path.exists() and meta_path.exists()):
        return AncestryEstimate(
            **unknown, reason="panel_missing", detail=f"no {REFERENCE_FILE} and {META_FILE} under {reference_dir}"
        )
    meta = json.loads(meta_path.read_text())
    if facet.assembly_label is None:
        return AncestryEstimate(**unknown, reason="assembly_unknown", detail=f"the facet states no assembly; the reference is on {meta['build']}")
    if facet.assembly_label != meta["build"]:
        return AncestryEstimate(
            **unknown, reason="assembly_mismatch", detail=f"the facet is on {facet.assembly_label}; the reference is on {meta['build']}"
        )

    groups: list[str] = meta["groups"]
    loads = [f"load{i}" for i in range(int(meta["n_pcs"]))]
    correction = np.asarray(meta["correction"], dtype=np.float64)
    sites = panel_sites(
        pl.read_parquet(reference_path, columns=_KEY).select(
            pl.col("chr").alias("chrom"), "pos", pl.col("a0").alias("ref_allele"), pl.col("a1").alias("effect_allele")
        )
    )
    counts = site_dosages(loaded, facet, sites)
    used = counts.dosage != DEFAULT_MISSING
    n_panel = sites.height
    measured = {
        **unknown,
        "panel_assembly": meta["build"],
        "n_sites_panel": n_panel,
        "n_sites_used": counts.n_used,
        "n_sites_observed": counts.n_observed,
        "n_sites_absent_as_reference": counts.n_absent_as_reference,
        "n_sites_unusable": counts.n_unusable,
        "n_sites_ref_unstated": counts.n_ref_unstated,
        "coverage": counts.n_used / n_panel if n_panel else 0.0,
    }
    logger.info(
        "%s on Privé: %d of %d sites (%d observed, %d absent read as reference, %d unusable, %d with a REF that states no base)",
        facet.facet_id, counts.n_used, n_panel, counts.n_observed, counts.n_absent_as_reference, counts.n_unusable, counts.n_ref_unstated,
    )
    if counts.n_used < max(min_sites, 1):
        return AncestryEstimate(**measured, reason="too_few_sites", detail=f"{counts.n_used} sites read, fewer than {min_sites}")

    table = pl.read_parquet(reference_path, columns=[*loads, *groups]).filter(pl.Series(used))
    loadings = table.select(loads).to_numpy()
    freqs = table.select(groups).to_numpy()
    dose = counts.dosage[used]

    blocks = np.array_split(np.arange(dose.size), min(n_blocks, dose.size))
    y_blocks = np.stack([dose[b] @ loadings[b] for b in blocks])
    x_blocks = np.stack([loadings[b].T @ freqs[b] for b in blocks])
    weights = _solve(y_blocks.sum(axis=0), x_blocks.sum(axis=0), correction)
    proportions = {group: round(float(w), 4) for group, w in zip(groups, weights, strict=True)}

    errors = None
    if n_bootstrap > 0:
        rng = np.random.default_rng(seed)
        draws = rng.multinomial(len(blocks), np.full(len(blocks), 1.0 / len(blocks)), size=n_bootstrap)
        replicates = np.stack(
            [_solve(np.tensordot(c, y_blocks, axes=1), np.tensordot(c, x_blocks, axes=1), correction) for c in draws]
        )
        errors = {group: round(float(e), 4) for group, e in zip(groups, replicates.std(axis=0, ddof=1), strict=True)}

    leader = max(groups, key=lambda group: (proportions[group], -groups.index(group)))
    return AncestryEstimate(
        **measured,
        estimate=leader,
        margin=margin(proportions),
        proportions=proportions,
        proportion_method="prive_qp",
        proportion_errors=errors,
        group_to_continental={group: continental_of(group) for group in groups},
    )


def _solve(y: np.ndarray, centroids: np.ndarray, correction: np.ndarray) -> np.ndarray:
    """Group weights for a projection `y` (uncorrected sums) and the groups' centroids `centroids`.

    The genome's coordinates are `y * correction / 2` (dosage counts two alleles, frequencies one); the
    centroids carry no correction, as in the vignette.
    """
    return simplex_lstsq(centroids, y * correction / 2.0)
