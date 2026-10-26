"""Relatedness by somalier (Brent Pedersen, MIT), run as a subprocess when its binary is on PATH (RM274).

just-vcf's own estimate is `just_vcf.relatedness`; this is the second opinion the person chose beside
it. Nothing here is required: without a `somalier` binary on PATH, `somalier_relate` raises
`SomalierNotFoundError` and `somalier_binary()` is None, so a caller (or a test) skips cleanly.

How it runs. Each facet is written out by `write_vcf`, cut down to the panel's positions, under a
sample name of its own (`f0`, `f1`, …), so two facets of one person that both call their sample `S1`
stay apart. `somalier extract` reads each one against the panel's sites VCF and a FASTA of the facet's
assembly (somalier reads bases from a FASTA; a refget store cannot stand in). A facet whose
`absence_means` is `weak_reference` is extracted with `--unknown`, which reads its silence as
homozygous reference, as `just_vcf.site_genotypes` does; any other facet's silence stays unknown.
`somalier relate` then compares them all, and its `pairs.tsv` is read back by `read_somalier_pairs`.

somalier's `relatedness` column is `2 (shared_hets - 2 ibs0) / het_ab`: on the scale where one person is
1, twice KING's kinship. The frame keeps it verbatim (`somalier_relatedness`) and puts half of it in
`kinship`, banded by `relatedness.band_of`; somalier reports no standard error, so the band is the
point estimate's. A first-degree band is split by somalier's own IBS0 count over its `n` sites, as
just-vcf's estimate is (`relatedness.first_degree_kind`).

**Not run on this machine yet**: `docs/TRUST.md` has no row for somalier's binary, and a source it does
not list is not run here until the person adds its row. The parser is tested on a `pairs.tsv` in the
layout somalier's `relate.nim` writes.
"""

import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from pathlib import Path

import polars as pl

from just_vcf.loading import Loaded
from just_vcf.relatedness import (
    DEFAULT_MIN_SITES,
    RELATEDNESS_FRAME,
    SOMALIER_METHOD,
    TOO_FEW_SITES,
    RelationBand,
    SitePanel,
    band_of,
    split_first_degree,
)
from just_vcf.schema import AbsenceMeans
from just_vcf.site_genotypes import facets_of
from just_vcf.writing import write_vcf

SOMALIER_PAIRS_COLUMNS = (
    "sample_a",
    "sample_b",
    "relatedness",
    "ibs0",
    "ibs2",
    "concordance",
    "hets_a",
    "hets_b",
    "hets_ab",
    "shared_hets",
    "hom_alts_a",
    "hom_alts_b",
    "shared_hom_alts",
    "n",
    "x_ibs0",
    "x_ibs2",
    "expected_relatedness",
)
"""The columns of somalier's `<prefix>pairs.tsv`, in order (its header line starts with `#`)."""


class SomalierNotFoundError(RuntimeError):
    """No `somalier` binary on PATH."""


def somalier_binary() -> Path | None:
    """The `somalier` binary on PATH, or None."""
    found = shutil.which("somalier")
    return None if found is None else Path(found)


def read_somalier_pairs(path: Path, names: dict[str, str] | None = None, *, min_sites: int = DEFAULT_MIN_SITES) -> pl.DataFrame:
    """somalier's `pairs.tsv` as a `RELATEDNESS_FRAME` (method `somalier`), one row per pair.

    `names` maps somalier's sample names to facet ids; a name it lacks is kept as written.
    """
    raw = pl.read_csv(path, separator="\t", has_header=True, infer_schema_length=0)
    raw = raw.rename(dict(zip(raw.columns, SOMALIER_PAIRS_COLUMNS, strict=True)))
    names = names or {}
    rows = []
    for row in raw.iter_rows(named=True):
        n = int(row["n"])
        value = float(row["relatedness"])
        counts = {
            "n_sites": n,
            "n_het_a": int(row["hets_a"]),
            "n_het_b": int(row["hets_b"]),
            "n_het_both": int(row["shared_hets"]),
            "n_ibs0": int(row["ibs0"]),
            "n_ibs2": int(row["ibs2"]),
            "n_vrs_compared": 0,
            "n_vrs_differ": 0,
        }
        known = n >= min_sites
        rows.append(
            {
                "facet_a": names.get(row["sample_a"], row["sample_a"]),
                "facet_b": names.get(row["sample_b"], row["sample_b"]),
                "method": SOMALIER_METHOD,
                "band": split_first_degree(band_of(value / 2), counts["n_ibs0"], n).value if known else RelationBand.UNKNOWN.value,
                "reason": None if known else TOO_FEW_SITES,
                "kinship": value / 2 if known else None,
                "ibs0": counts["n_ibs0"] / n if known and n else None,
                "ibs2": counts["n_ibs2"] / n if known and n else None,
                "somalier_relatedness": value,
                **counts,
            }
        )
    return pl.DataFrame(rows, schema=RELATEDNESS_FRAME) if rows else pl.DataFrame(schema=RELATEDNESS_FRAME)


def somalier_relate(
    loaded: Sequence[Loaded],
    *,
    panel: SitePanel,
    sites_vcf: Path,
    fasta: Path,
    workdir: Path | None = None,
    min_sites: int = DEFAULT_MIN_SITES,
) -> pl.DataFrame:
    """somalier's relatedness of every pair of facets across `loaded`, as a `RELATEDNESS_FRAME`.

    `sites_vcf` is the sites file somalier extracts at (one of its published `sites.*.vcf.gz`), `panel`
    the same sites as a `SitePanel` (the facets are cut down to its positions before writing), `fasta`
    the facets' assembly. Every facet must have its REF read (`resolve_ref`), as `write_vcf` needs.
    Files go under `workdir`, a temporary directory by default.
    """
    binary = somalier_binary()
    if binary is None:
        raise SomalierNotFoundError("no somalier binary on PATH")
    with tempfile.TemporaryDirectory(dir=workdir) as scratch:
        root = Path(scratch)
        extracted = root / "extracted"
        names: dict[str, str] = {}
        index = 0
        for item in loaded:
            for facet in facets_of(item):
                name = f"f{index}"
                index += 1
                names[name] = facet.facet_id
                one = Loaded(
                    facets=item.facets.filter(pl.col("facet_id") == facet.facet_id),
                    calls=item.calls.filter(pl.col("facet_id") == facet.facet_id).join(
                        panel.sites.lazy().select("chrom", "pos"), on=["chrom", "pos"], how="semi"
                    ),
                    blocks=item.blocks.filter(pl.col("facet_id") == facet.facet_id),
                )
                vcf = root / f"{name}.vcf.gz"
                write_vcf(one, vcf, sample=name, resolved_only=True)
                command = [str(binary), "extract", "-d", str(extracted), "--sites", str(sites_vcf), "-f", str(fasta)]
                if facet.absence_means == AbsenceMeans.WEAK_REFERENCE.value:
                    command.append("--unknown")
                subprocess.run([*command, str(vcf)], check=True, capture_output=True)
        prefix = root / "related"
        subprocess.run(
            [str(binary), "relate", "-o", str(prefix), *sorted(str(path) for path in extracted.glob("*.somalier"))],
            check=True,
            capture_output=True,
        )
        # somalier names it `<prefix>.pairs.tsv` (README: `somalier.pairs.tsv` by default).
        return read_somalier_pairs(Path(f"{prefix}.pairs.tsv"), names, min_sites=min_sites)
