"""A reference-PCA panel built from AADR's Human Origins release (`aadr_ho`, GRCh37).

Ported from just-prs `ancestry/aadr.py` (docs/PROVENANCE.md). AADR, the Allen Ancient DNA Resource
(Reich lab, Harvard Dataverse), publishes ancient and present-day individuals genotyped at the Human
Origins array's sites, in EIGENSTRAT files (`just_vcf.ancestry.eigenstrat`). The panel is its
present-day individuals inside a region box, by default just-prs's West Eurasia, which is what the
1000 Genomes and HGDP panels do not resolve: fine populations within Europe, the Caucasus, the Near
East and Central Asia. docs/ANCESTRY_PANELS.md has the rules in prose and a measured build.

Who goes in, each rule applied in this order, the first one an individual fails being its reason:

1. in the `.anno` file. An `.anno` listing one Genetic ID twice is refused before anything is read
   from the `.geno`: the join would select that person's genotypes twice;
2. present-day: `Date mean in BP` is 0;
3. diploid data: `Data type` is `HumanOrigins` (the array) or `Shotgun.diploid` (a genome with diploid
   calls). Pseudo-haploid shotgun and capture data (one read per site) and the reference genomes
   (chimpanzee, gorilla, the ancestral and human references) are out;
4. `ASSESSMENT` is `Pass`, `PROVISIONAL_PASS` or `MERGE_PASS`;
5. the `.ind` group label does not start with `Ignore` and carries no QC qualifier (`group_label_flag`,
   `QC_QUALIFIERS`). A label is a name, then any number of `-` qualifiers. AADR's README names `Ignore_`
   and outliers as the marks of individuals to filter out; the qualifiers that say so are `-o` and
   every `-o…` (`-oPCA`, `-oRelative`, `-oAfrica`) and the ones in `QC_QUALIFIERS`. A qualifier of
   digits only names a sub-group (`Kumyk-1`, `Hazara_1-1`) and is kept. Any other qualifier is neither,
   so it is counted under a reason of its own rather than taken for either;
6. a latitude and a longitude;
7. inside the box (`WEST_EURASIA` unless the caller passes another);
8. a region for its `Political Entity` (`AADR_REGIONS`, `region_of`);
9. one representation per person: AADR lists some people twice (`Individual ID` shared by two
   `Genetic ID`s, typically a Human Origins genotype and a shotgun genome). The `HumanOrigins` one is
   kept, else the first in `.ind` order.

`population` is the `.ind` group label as written: AADR's numbered sub-groups (`Kumyk-1`, `Kumyk-2`)
are separate groups there and stay separate here. `superpop` is the region: a United Nations M49 geographic
region of the political entity, Europe as one, Sub-Saharan Africa as one, and Russia split at 60°E
(the Urals) into Europe and Siberia, which M49 does not do. It is geography, not genetics: a
Kalmyk sample at 45°E is in Europe. Folding onto the five continental groups goes through
`consensus.continental_of` on the region's name: the European, African and Latin American regions
fold, the Asian ones, Oceania and Siberia do not and are dropped from a consensus, as Middle Eastern
groups are.

Which sites, each rule in this order:

1. on an autosome;
2. one base each, both A, C, G or T;
3. not strand-ambiguous (A/T, C/G). Human Origins has none by design (AADR's README); the rule stays
   for another release's files;
4. one record at the position (both are dropped where there are two);
5. on a contig the GRCh37 reference holds;
6. one allele is the reference base there. The genotypes count `.snp`'s first allele; where the
   reference base is the second, the allele order is swapped and the counts used as ALT dosage as they
   are. REF is never taken from the file or from a genotype;
7. minor allele frequency over the chosen individuals at least `min_maf`.

No LD pruning, as in just-prs: the array's sites were ascertained to be informative already.
"""

import logging
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.ancestry.eigenstrat import MISSING, geno_layout, read_anno, read_genotypes, read_ind, read_snp
from just_vcf.ancestry.pca import DEFAULT_KNN_K, PanelBuild, build_pca_panel
from just_vcf.reference import ReferenceStore
from just_vcf.schema import CONTIG_RANK, canonical_chrom_expr
from just_vcf.site_genotypes import AUTOSOMES

logger = logging.getLogger(__name__)

AADR_ASSEMBLY = "GRCh37"
"""AADR's positions are hg19 (its README)."""

AADR_PANEL = "aadr_ho"
"""just-prs's name for the panel, kept so a panel either built reads in the other."""

DEFAULT_AADR_DIM_REF = 20
"""just-prs's choice: fine structure within Europe needs deeper components than the default 10."""

DEFAULT_MIN_MAF = 0.02

WEST_EURASIA: tuple[float, float, float, float] = (30.0, 75.0, -25.0, 90.0)
"""just-prs's box, `(lat_min, lat_max, lon_min, lon_max)`: Europe, the north of Africa, the Near East,
the Caucasus and western and central Asia."""

DIPLOID_DATA_TYPES: frozenset[str] = frozenset({"HumanOrigins", "Shotgun.diploid"})
PASSING_ASSESSMENTS: frozenset[str] = frozenset({"Pass", "PROVISIONAL_PASS", "MERGE_PASS"})

SIBERIA_FROM_LONGITUDE = 60.0
"""Russia's individuals at or east of this longitude are filed under Siberia, west of it under Europe."""

AADR_REGIONS: dict[str, str | None] = {
    # Europe (M49: Eastern, Northern, Southern and Western Europe)
    **dict.fromkeys(
        [
            "Albania", "Belarus", "Bosnia-Herzegovina", "Bulgaria", "Croatia", "Czechia", "Estonia", "Finland",
            "France", "Greece", "Hungary", "Iceland", "Italy", "Lithuania", "Malta", "Moldova", "Montenegro",
            "North Macedonia", "Norway", "Poland", "Romania", "Russia", "Serbia", "Spain", "Ukraine",
            "United Kingdom",
        ],
        "Europe",
    ),
    **dict.fromkeys(["Algeria", "Egypt", "Libya", "Morocco", "Sudan", "Tunisia"], "Northern Africa"),
    **dict.fromkeys(
        [
            "Angola", "Cameroon", "Central African Republic", "Congo", "Eritrea", "Ethiopia", "Gambia", "Kenya",
            "Lesotho", "Malawi", "Namibia", "Nigeria", "Senegal", "Sierra Leone", "South Africa", "South Sudan",
            "Tanzania",
        ],
        "Sub-Saharan Africa",
    ),
    **dict.fromkeys(
        [
            "Abkhazia", "Armenia", "Azerbaijan", "Cyprus", "Georgia", "Iraq", "Israel", "Jordan", "Lebanon",
            "Saudi Arabia", "Syria", "Turkey", "United Arab Emirates", "Yemen",
        ],
        "Western Asia",
    ),
    **dict.fromkeys(["Kazakhstan", "Kyrgyzstan", "Tajikistan", "Uzbekistan"], "Central Asia"),
    **dict.fromkeys(["Bangladesh", "India", "Iran", "Nepal", "Pakistan", "Sri Lanka"], "Southern Asia"),
    **dict.fromkeys(["China", "Japan", "Mongolia", "Republic of Korea", "Taiwan"], "Eastern Asia"),
    **dict.fromkeys(
        ["Brunei", "Cambodia", "Laos", "Myanmar", "Philippines", "Singapore", "Thailand", "Vietnam"],
        "South-eastern Asia",
    ),
    **dict.fromkeys(["Bolivia", "Brazil", "Colombia", "Mexico", "Peru", "Puerto Rico"], "Latin America and the Caribbean"),
    **dict.fromkeys(["Australia", "Federated States of Micronesia", "New Zealand", "Papua New Guinea"], "Oceania"),
    # Where AADR's present-day samples were collected away from their ancestry's region (1000 Genomes'
    # CEU and ASW in the USA, ACB in Barbados), the place says nothing about the group: no region.
    "USA": None,
    "Barbados": None,
}
"""Each `Political Entity` of AADR v66.p1's present-day diploid individuals to its region (M49, with
Abkhazia under Georgia's and Taiwan under Eastern Asia). An entity not listed, or listed as None, puts
its individuals out of the panel, counted."""

QC_QUALIFIERS: frozenset[str] = frozenset(
    {
        "QCremove",  # failed quality control; v66.p1 `.ind`, e.g. `Kumyk-QCremove`
        "Dup",  # a duplicate of another individual; v66.p1 `.ind`, e.g. `Adygei-Dup`
        "Discovery",  # a Human Origins ascertainment sample (Patterson 2012: the sites were found in its genome); v66.p1 `.ind`, e.g. `French-Discovery`
        "NotPermittedToPublish",  # v66.p1 `.ind`, one ancient individual
        "DontUseHighlyProblematic",  # v66.p1 `.ind`, one ancient individual
    }
)
"""The `.ind` group-label qualifiers, besides `-o` and every `-o…` (outliers, AADR's README), that mark an
individual to leave out, each read off AADR v66.p1's `compatibility_HO` `.ind` file (2026-10-26). The
other qualifiers there are digits (a numbered sub-group, kept) or name an ancient analysis cluster
(`-A`, `-B`, `-C`, `-lowEastEU`, `-highEastEU`, `-Africa`, `-3CentralEuropean`); none of those reaches
rule 5 in v66.p1, since every individual carrying one is ancient."""

LABEL_FLAGGED = "group label flagged (Ignore_ or a QC qualifier)"
LABEL_UNRECOGNISED = "group label qualifier not recognised"

INDIVIDUAL_EXCLUSIONS: tuple[str, ...] = (
    "not in the annotation file",
    "ancient (date mean in BP not 0)",
    "not diploid data (HumanOrigins or Shotgun.diploid)",
    "assessment not a pass",
    LABEL_FLAGGED,
    LABEL_UNRECOGNISED,
    "no latitude or longitude",
    "outside the region box",
    "political entity has no ancestry region",
    "second representation of one individual",
)
"""Why an individual of the `.ind` file is not in the panel, in the order the rules apply."""

SITE_DROPS: tuple[str, ...] = (
    "not on an autosome",
    "not a one-base ACGT SNV",
    "strand-ambiguous (A/T or C/G)",
    "two records at one position",
    "contig not in the reference",
    "neither allele is the reference base",
    "minor allele frequency below the floor",
)
"""Why a site of the `.snp` file is not in the panel, in the order the rules apply."""

_PALINDROMES = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


class AadrBuild(BaseModel):
    """What `build_aadr_panel` read, kept and wrote."""

    model_config = ConfigDict(frozen=True)

    panel: PanelBuild
    n_individuals: int
    """Lines of the `.ind` file."""
    n_selected: int
    excluded: dict[str, int]
    """Per `INDIVIDUAL_EXCLUSIONS` reason, every one listed."""
    regions: dict[str, int]
    """Chosen individuals per region."""
    n_populations: int
    n_snps: int
    """Lines of the `.snp` file."""
    n_sites: int
    dropped: dict[str, int]
    """Per `SITE_DROPS` reason, every one listed."""
    n_swapped: int
    """Kept sites whose `.snp` line lists the reference base second."""


def group_label_flag(group: str) -> str | None:
    """Why rule 5 leaves a `.ind` group label out (`LABEL_FLAGGED` or `LABEL_UNRECOGNISED`), or None
    where it keeps it.

    The label is split at every `-`: the first piece is the name, each later piece a qualifier, and
    every qualifier is read, so `Ulgii-1-Dup` is flagged for its `Dup`. `LABEL_FLAGGED`: the label
    starts with `Ignore`, or a qualifier starts with `o` (an outlier) or is in `QC_QUALIFIERS`.
    Otherwise `LABEL_UNRECOGNISED` where a qualifier is not digits (an empty one included), else kept.
    """
    name, *qualifiers = group.split("-")
    if name.startswith("Ignore") or any(q.startswith("o") or q in QC_QUALIFIERS for q in qualifiers):
        return LABEL_FLAGGED
    if any(not (q.isascii() and q.isdigit()) for q in qualifiers):
        return LABEL_UNRECOGNISED
    return None


def region_of(political_entity: str, longitude: float) -> str | None:
    """The region an individual is filed under (`AADR_REGIONS`, Russia split at the Urals), or None."""
    region = AADR_REGIONS.get(political_entity)
    if political_entity == "Russia" and longitude >= SIBERIA_FROM_LONGITUDE:
        return "Siberia"
    return region


def build_aadr_panel(
    prefix: Path,
    anno: Path,
    panel_dir: Path,
    *,
    reference: ReferenceStore | Path | None = None,
    panel: str = AADR_PANEL,
    box: tuple[float, float, float, float] = WEST_EURASIA,
    min_maf: float = DEFAULT_MIN_MAF,
    dim_ref: int = DEFAULT_AADR_DIM_REF,
    knn_k: int = DEFAULT_KNN_K,
) -> AadrBuild:
    """Build the AADR Human Origins panel and write it under `panel_dir` (`build_pca_panel`'s files).

    `prefix` is the path of the `.geno`, `.snp` and `.ind` files without the suffix (AADR v66.p1:
    `v66.p1_compatibility_HO.aadr.patch.PUB`), `anno` the `.anno` file. `reference` is a refget store
    holding exactly one GRCh37 collection (default: `ReferenceStore.open()`); every kept site's REF is
    read from it, and a store without one is refused. Fewer than `8 * dim_ref` chosen individuals is
    refused, as just-prs did. The module docstring has the rules.
    """
    store = reference if isinstance(reference, ReferenceStore) else ReferenceStore.open(reference)
    digest = store.collection_for(AADR_ASSEMBLY)
    geno, snp_path, ind_path = (Path(f"{prefix}.{suffix}") for suffix in ("geno", "snp", "ind"))
    layout = geno_layout(geno)

    ind = read_ind(ind_path).with_row_index("row")
    annotations = read_anno(anno)
    repeated = annotations["genetic_id"].is_duplicated().sum()
    if repeated:
        # Joined to the .ind, a repeated row would select that person's .geno column twice.
        n_ids = annotations.filter(annotations["genetic_id"].is_duplicated())["genetic_id"].n_unique()
        raise ValueError(f"{anno.name}: {n_ids} Genetic IDs appear more than once ({repeated} rows); each must be listed once")
    chosen, excluded, regions = _choose(ind, annotations, box)
    if chosen.height < 8 * dim_ref:
        raise ValueError(f"{chosen.height} individuals chosen; {dim_ref} components need at least {8 * dim_ref}")

    snp = read_snp(snp_path).with_row_index("row")
    sites, dropped = _sites(snp, store, digest)
    counts = read_genotypes(geno, ind.height, snp.height, chosen["row"].to_numpy())[sites["row"].to_numpy()]
    swapped = sites["swapped"].to_numpy()[:, None]
    dosage = np.where(counts == MISSING, MISSING, np.where(swapped, counts, 2 - counts)).astype(np.int8)
    called = dosage != MISSING
    n_called = called.sum(axis=1)
    alt_alleles = np.where(called, dosage, 0).sum(axis=1, dtype=np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        frequency = alt_alleles / (2.0 * n_called)
    common = np.isfinite(frequency) & (np.minimum(frequency, 1.0 - frequency) >= min_maf)
    dropped["minor allele frequency below the floor"] = int((~common).sum())
    sites = sites.filter(pl.Series(common))
    dosage = dosage[common]

    labels = chosen.select(pl.col("genetic_id").alias("iid"), pl.col("region").alias("superpop"), pl.col("group").alias("population"))
    extra = {
        "source": {
            "geno": geno.name,
            "geno_layout": layout,
            "snp": snp_path.name,
            "ind": ind_path.name,
            "anno": anno.name,
        },
        "selection": {
            "box": list(box),
            "data_types": sorted(DIPLOID_DATA_TYPES),
            "assessments": sorted(PASSING_ASSESSMENTS),
            "excluded": excluded,
            "regions": regions,
        },
        "sites": {"min_maf": min_maf, "dropped": dropped, "swapped": int(sites["swapped"].sum())},
    }
    built = build_pca_panel(
        dosage,
        sites.select("chrom", "pos", "ref", "alt"),
        labels,
        panel=panel,
        assembly=AADR_ASSEMBLY,
        panel_dir=panel_dir,
        dim_ref=dim_ref,
        knn_k=knn_k,
        extra_meta=extra,
    )
    logger.info(
        "%s: %d of %d individuals (%s), %d of %d sites (%s)",
        panel, chosen.height, ind.height, _nonzero(excluded), sites.height, snp.height, _nonzero(dropped),
    )
    return AadrBuild(
        panel=built,
        n_individuals=ind.height,
        n_selected=chosen.height,
        excluded=excluded,
        regions=regions,
        n_populations=chosen["group"].n_unique(),
        n_snps=snp.height,
        n_sites=sites.height,
        dropped=dropped,
        n_swapped=int(sites["swapped"].sum()),
    )


def _nonzero(counts: dict[str, int]) -> str:
    return "; ".join(f"{reason}: {n}" for reason, n in counts.items() if n) or "none dropped"


def _region_of_row(row: dict) -> str | None:
    """`region_of` for a row of the joined frame; no region where the entity or longitude is missing."""
    if row["political_entity"] is None or row["longitude"] is None:
        return None
    return region_of(row["political_entity"], row["longitude"])


def _choose(ind: pl.DataFrame, anno: pl.DataFrame, box: tuple[float, float, float, float]) -> tuple[pl.DataFrame, dict[str, int], dict[str, int]]:
    """The chosen individuals (`row, genetic_id, group, region`, in `.ind` order), the count per
    exclusion reason, and the count per region."""
    lat_min, lat_max, lon_min, lon_max = box
    flags = {group: group_label_flag(group) for group in ind["group"].unique().to_list()}
    joined = ind.join(anno, on="genetic_id", how="left", maintain_order="left").with_columns(
        pl.struct("political_entity", "longitude")
        .map_elements(_region_of_row, return_dtype=pl.String)
        .alias("region"),
        pl.col("group").replace_strict(flags, return_dtype=pl.String).alias("label_flag"),
    )
    rules = [
        pl.col("individual_id").is_not_null(),
        pl.col("date_bp") == 0,
        pl.col("data_type").is_in(DIPLOID_DATA_TYPES),
        pl.col("assessment").is_in(PASSING_ASSESSMENTS),
        pl.col("label_flag").fill_null("") != LABEL_FLAGGED,
        pl.col("label_flag").fill_null("") != LABEL_UNRECOGNISED,
        pl.col("latitude").is_not_null() & pl.col("longitude").is_not_null(),
        pl.col("latitude").is_between(lat_min, lat_max) & pl.col("longitude").is_between(lon_min, lon_max),
        pl.col("region").is_not_null(),
    ]
    why = pl.lit(None, dtype=pl.String)
    for rule, reason in reversed(list(zip(rules, INDIVIDUAL_EXCLUSIONS[:-1], strict=True))):
        why = pl.when(~rule.fill_null(False)).then(pl.lit(reason)).otherwise(why)
    joined = joined.with_columns(why.alias("why"))
    passing = joined.filter(pl.col("why").is_null())
    # One representation per person: the HumanOrigins one, else the first in .ind order.
    kept_rows = (
        passing.sort(pl.col("data_type") != "HumanOrigins", "row")
        .unique("individual_id", keep="first", maintain_order=True)["row"]
    )
    joined = joined.with_columns(
        pl.when(pl.col("why").is_null() & ~pl.col("row").is_in(kept_rows.to_list()))
        .then(pl.lit(INDIVIDUAL_EXCLUSIONS[-1]))
        .otherwise(pl.col("why"))
        .alias("why")
    )
    tally = Counter(joined["why"].drop_nulls().to_list())
    excluded = {reason: tally.get(reason, 0) for reason in INDIVIDUAL_EXCLUSIONS}
    unrecognised = Counter(joined.filter(pl.col("why") == LABEL_UNRECOGNISED)["group"].to_list())
    if unrecognised:
        logger.warning(
            "%d individuals left out for a group-label qualifier that is neither a sub-group number nor a known QC mark: %s",
            sum(unrecognised.values()),
            ", ".join(f"{group} ({n})" for group, n in sorted(unrecognised.items())),
        )
    chosen = joined.filter(pl.col("why").is_null()).select("row", "genetic_id", "group", "region")
    regions = dict(sorted(Counter(chosen["region"].to_list()).items()))
    return chosen, excluded, regions


def _sites(snp: pl.DataFrame, store: ReferenceStore, digest: str) -> tuple[pl.DataFrame, dict[str, int]]:
    """The sites that pass every rule but the frequency floor (`row, chrom, pos, ref, alt, swapped`,
    in canonical contig order then position) and the count per drop reason, the floor's still 0."""
    dropped = dict.fromkeys(SITE_DROPS, 0)
    work = snp.with_columns(canonical_chrom_expr("chrom").alias("chrom"))
    on_autosome = pl.col("chrom").is_in(AUTOSOMES)
    dropped["not on an autosome"] = work.filter(~on_autosome).height
    work = work.filter(on_autosome)
    one_base = pl.col("allele1").is_in(list("ACGT")) & pl.col("allele2").is_in(list("ACGT")) & (pl.col("allele1") != pl.col("allele2"))
    dropped["not a one-base ACGT SNV"] = work.filter(~one_base).height
    work = work.filter(one_base)
    palindrome = pl.concat_str("allele1", "allele2").is_in(["".join(pair) for pair in _PALINDROMES])
    dropped["strand-ambiguous (A/T or C/G)"] = work.filter(palindrome).height
    work = work.filter(~palindrome)
    unique = work.unique(["chrom", "pos"], keep="none", maintain_order=True)
    dropped["two records at one position"] = work.height - unique.height
    work = unique

    contigs = store.contigs_by_chrom(digest)
    parts = []
    for chrom in sorted(work["chrom"].unique().to_list(), key=lambda c: CONTIG_RANK.get(c, len(CONTIG_RANK))):
        rows = work.filter(pl.col("chrom") == chrom)
        serving = contigs.get(chrom, [])
        if len(serving) != 1:
            dropped["contig not in the reference"] += rows.height
            continue
        positions = rows["pos"].to_list()
        bases = store.fetch(serving[0], [(pos - 1, pos) for pos in positions])
        parts.append(rows.with_columns(pl.Series("base", bases, dtype=pl.String)))
    empty = work.clear().with_columns(pl.lit(None, dtype=pl.String).alias("base"))
    checked = pl.concat(parts) if parts else empty
    oriented = checked.with_columns(
        (pl.col("allele2") == pl.col("base")).alias("swapped"),
        ((pl.col("allele1") == pl.col("base")) | (pl.col("allele2") == pl.col("base"))).alias("matches"),
    )
    dropped["neither allele is the reference base"] = oriented.filter(~pl.col("matches")).height
    sites = oriented.filter(pl.col("matches")).select(
        "row",
        "chrom",
        "pos",
        pl.col("base").alias("ref"),
        pl.when(pl.col("swapped")).then(pl.col("allele1")).otherwise(pl.col("allele2")).alias("alt"),
        "swapped",
    )
    return sites.sort(pl.col("chrom").replace_strict(CONTIG_RANK, return_dtype=pl.UInt16), "pos", maintain_order=True), dropped
