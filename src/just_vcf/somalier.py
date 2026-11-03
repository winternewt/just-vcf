"""Relatedness by somalier (Brent Pedersen, MIT), run as a subprocess when its binary is on PATH (RM274).

just-vcf's own estimate is `just_vcf.relatedness`; this is the second opinion the person chose beside
it. Nothing here is required: without a `somalier` binary on PATH, `somalier_relate` raises
`SomalierNotFoundError` and `somalier_binary()` is None, so a caller (or a test) skips cleanly.

How it runs. Each facet is read at the panel's sites by `just_vcf.site_genotypes`, the reading
just-vcf's own relatedness takes, and written out by `write_vcf` as one record per site that reading
gives a count (`somalier_sites`): the panel's REF and effect allele, GT `0/0`, `0/1` or `1/1`, and no
other FORMAT key. A site it leaves unknown gets no record. Each file has a sample name of its own
(`f0`, `f1`, …), so two facets of one person that both call their sample `S1` stay apart. `somalier
extract` reads each one against the panel's sites VCF and a FASTA of the facet's assembly (somalier
reads bases from a FASTA; a refget store cannot stand in), and `somalier relate` compares them all,
without `--unknown`: a site with no record is unknown to somalier, as it is to just-vcf. Its
`pairs.tsv` is read back by `read_somalier_pairs`.

Why not the facet's own rows and `--unknown`, as before 2026-11-03. `--unknown` is a flag of `somalier
relate`, not of `extract` (somalier's `relate.nim`, `unk2hr`; `extract` takes only `--sites`,
`--fasta`, `--out-dir` and `--sample-prefix`), so it cannot be set per facet: it reads every unknown
of every sample in the run as homozygous reference, an array's silence and a no-call included. Nor can
it tell a site the writer withheld (`resolved_only=True` skips a `ref_unstated` or `ref_mismatch` row,
among others) from one the facet never wrote: neither has a record, and both would read as zero
copies, where `site_genotypes` reads the first unknown (`ref_unstated`, `unusable`). So the facet's
absence rule (`absence_means`: `weak_reference`, `blocks`, `reference`) is applied here, by the reading
that already applies it, and written as the `0/0` records it implies. The records carry no AD: where a
header declares AD, somalier reads a record without one as depth zero, unknown, so a `0/0` written
from the facet's silence beside rows with AD would be lost. Without AD somalier reads each GT at fixed
counts (its `fill`, 20 reads a genotype), so it compares the genotypes just-vcf reads rather than
calling its own from allele depths; its relatedness estimate over them is still its own.

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
from just_vcf.schema import BLOCKS, CALLS, EvidenceBasis, canonical_order, conform, empty
from just_vcf.site_genotypes import Facet, SiteBasis, facets_of, read_site_genotypes
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


def somalier_sites(loaded: Loaded, facet: Facet, panel: SitePanel) -> Loaded:
    """`facet` of `loaded` as somalier is given it: one `calls` row per panel site whose genotype
    `site_genotypes.read_site_genotypes` reads, none where it reads unknown.

    Each row is the panel's site: REF its `ref_allele`, the one ALT its `effect_allele`, GT the count
    read there (`0/0`, `0/1`, `1/1`), `evidence_basis` `assay` where the facet called the site and
    `absence` where its silence reads as reference. No other column is filled: somalier reads a
    record by its GT where no AD is declared (module docstring). A site whose rows the writer would
    withhold (`ref_unstated`, `ref_mismatch`, …) is never read as reference here: `read_site_genotypes`
    reads a site with rows by those rows alone, and leaves it unknown when none is usable.
    """
    read = read_site_genotypes(loaded, facet, panel.sites).filter(pl.col("count").is_not_null())
    sites = panel.sites.with_row_index("site").join(read, on="site", how="inner")
    count = pl.col("count")
    ref, effect = pl.col("ref_allele"), pl.col("effect_allele")
    calls = sites.lazy().select(
        pl.lit(facet.facet_id).alias("facet_id"),
        "chrom",
        "pos",
        pl.col("pos").alias("end_pos"),
        ref.alias("ref"),
        pl.concat_list(effect).alias("alts"),
        pl.when(count == 0).then(pl.lit([0, 0])).when(count == 1).then(pl.lit([0, 1])).otherwise(pl.lit([1, 1])).cast(CALLS["gt"]).alias("gt"),
        pl.when(count == 0).then(pl.concat_list(ref, ref)).when(count == 1).then(pl.concat_list(ref, effect)).otherwise(pl.concat_list(effect, effect)).alias("gt_alleles"),
        pl.when(pl.col("basis") == SiteBasis.OBSERVED.value)
        .then(pl.lit(EvidenceBasis.ASSAY.value))
        .otherwise(pl.lit(EvidenceBasis.ABSENCE.value))
        .alias("evidence_basis"),
    )
    return Loaded(
        facets=loaded.facets.filter(pl.col("facet_id") == facet.facet_id),
        calls=canonical_order(conform(calls, CALLS)),
        blocks=empty(BLOCKS),
    )


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
    the same sites as a `SitePanel` (each facet is read at its sites, `somalier_sites`), `fasta` the
    facets' assembly. somalier reads every facet as `site_genotypes` does: a site that reading leaves
    unknown is unknown to somalier too, and `relate` runs without `--unknown`. Files go under
    `workdir`, a temporary directory by default.
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
                vcf = root / f"{name}.vcf.gz"
                write_vcf(somalier_sites(item, facet, panel), vcf, sample=name)
                command = [str(binary), "extract", "-d", str(extracted), "--sites", str(sites_vcf), "-f", str(fasta), str(vcf)]
                subprocess.run(command, check=True, capture_output=True)
        prefix = root / "related"
        subprocess.run(
            [str(binary), "relate", "-o", str(prefix), *sorted(str(path) for path in extracted.glob("*.somalier"))],
            check=True,
            capture_output=True,
        )
        # somalier names it `<prefix>.pairs.tsv` (README: `somalier.pairs.tsv` by default).
        return read_somalier_pairs(Path(f"{prefix}.pairs.tsv"), names, min_sites=min_sites)
