"""Relatedness (RM274) on a synthetic family drawn here with a fixed seed.

Founders draw each haplotype from a site's allele frequency; a child takes one haplotype of each
parent at every site. Each person is written as a real VCF or 23andMe file and read through `load`, so
the counts a relatedness reports are checked against the genotypes the fixture drew.
"""

import os
import random
import sys
from pathlib import Path

import polars as pl
import pytest
from dotenv import dotenv_values, find_dotenv

from just_vcf import FACETS, Loaded, load, mint_vrs_ids, resolve_ref
from just_vcf.assembly import ASSEMBLY_LENGTHS
from just_vcf.pedigree import BAND_ORDER, Pedigree, pedigree_from_dict
from just_vcf.person import NO_FACET, PEDIGREE_CHECKS, Person, check_pedigree
from just_vcf.person import people as people_of
from just_vcf.reference import ReferenceStore
from just_vcf.sniff import sniff
from just_vcf.relatedness import (
    ASSEMBLY_MISMATCH,
    ASSEMBLY_UNKNOWN,
    BETWEEN_BANDS,
    NO_PANEL,
    RELATEDNESS_FRAME,
    THIRD_DEGREE_FLOOR,
    TOO_FEW_SITES,
    VARIANT_SITES_ONLY,
    Relatedness,
    RelationBand,
    SitePanel,
    band_of,
    interval_band,
    king_robust,
    packaged_assemblies,
    packaged_panel,
    pairwise_relatedness,
    relatedness,
    same_person,
)

SEED = 20261023
CONTIGS = ("1", "2", "3")
SITES_PER_CONTIG = 3000
N_SITES = SITES_PER_CONTIG * len(CONTIGS)
PAIRS = (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"), ("A", "C"), ("C", "A"), ("G", "T"), ("T", "G"))
GT_TEXT = {(0, 0): "0/0", (0, 1): "0/1", (1, 0): "0/1", (1, 1): "1/1"}


def _length(contig: str) -> int:
    return next(iter(ASSEMBLY_LENGTHS["GRCh37"][contig]))


@pytest.fixture(scope="module")
def family() -> tuple[SitePanel, list[float], dict[str, list[tuple[int, int]]]]:
    """A panel, its ALT frequencies, and the haplotype pairs of eleven people (ALT = 1).

    `a`, `b`, `e`, `u`, `d` and `f` are founders; `c` is the child of `a` and `b`, `s` a full sibling
    of `c`, `g` the child of `c` and `e` (so `a`'s grandchild), `h` the child of `a` and `d`, `k` the
    child of `s` and `f`.
    """
    rng = random.Random(SEED)
    n = SITES_PER_CONTIG * len(CONTIGS)
    pairs = [PAIRS[rng.randrange(len(PAIRS))] for _ in range(n)]
    sites = pl.DataFrame(
        {
            "chrom": [c for c in CONTIGS for _ in range(SITES_PER_CONTIG)],
            "pos": pl.Series([10_000 + 1_000 * i for _ in CONTIGS for i in range(SITES_PER_CONTIG)], dtype=pl.UInt32),
            "ref_allele": [r for r, _ in pairs],
            "effect_allele": [a for _, a in pairs],
        }
    )
    freq = [rng.uniform(0.1, 0.5) for _ in range(n)]

    def founder() -> list[tuple[int, int]]:
        return [(int(rng.random() < p), int(rng.random() < p)) for p in freq]

    def child(mother: list[tuple[int, int]], father: list[tuple[int, int]]) -> list[tuple[int, int]]:
        return [(m[rng.randrange(2)], f[rng.randrange(2)]) for m, f in zip(mother, father)]

    people = {name: founder() for name in ("a", "b", "e", "u")}
    people["c"] = child(people["a"], people["b"])
    people["s"] = child(people["a"], people["b"])
    people["g"] = child(people["c"], people["e"])
    # Drawn after the others, so their genotypes stay as they were: `h` a half sibling of `c` (by `a`
    # and the founder `d`), `k` a first cousin of `g` (by `s` and the founder `f`).
    people |= {name: founder() for name in ("d", "f")}
    people["h"] = child(people["a"], people["d"])
    people["k"] = child(people["s"], people["f"])
    return SitePanel(assembly_label="GRCh37", sites=sites, source="synthetic"), freq, people


def write_vcf(path: Path, panel: SitePanel, haplotypes: list[tuple[int, int]], *, variant_only: bool = True) -> Path:
    """A one-sample VCF on GRCh37; with `variant_only`, hom-ref sites are left out, as a WGS VCF does."""
    lines = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={_length(c)}>" for c in CONTIGS)]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]))
    for (chrom, pos, ref, alt), pair in zip(panel.sites.iter_rows(), haplotypes, strict=True):
        if variant_only and pair == (0, 0):
            continue
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{GT_TEXT[pair]}")
    path.write_text("\n".join(lines) + "\n")
    return path


def wgs(path: Path, subject: str | None = None) -> Loaded:
    """A test VCF loaded as a whole genome: absent means reference (the loader cannot tell a small
    file's callset type and gives it `unknown`)."""
    loaded = load(path, subject=subject)
    facets = loaded.facets.with_columns(pl.lit("weak_reference").cast(FACETS["absence_means"]).alias("absence_means"))
    return Loaded(facets=facets, calls=loaded.calls, blocks=loaded.blocks)


def write_23andme(path: Path, panel: SitePanel, haplotypes: list[tuple[int, int]], *, no_call: set[int] = frozenset(), skip: set[int] = frozenset()) -> Path:
    """A 23andMe build-37 export: every site's letters, `--` at `no_call`, rows at `skip` left out."""
    lines = [
        "# This data file generated by 23andMe at: Wed Apr 11 08:10:59 2018",
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for i, ((chrom, pos, ref, alt), pair) in enumerate(zip(panel.sites.iter_rows(), haplotypes, strict=True)):
        if i in skip:
            continue
        letters = "--" if i in no_call else "".join(sorted((alt if h else ref) for h in pair))
        lines.append(f"rs{i + 1}\t{chrom}\t{pos}\t{letters}")
    path.write_text("\n".join(lines) + "\n")
    return path


def counts_of(first: list[tuple[int, int]], second: list[tuple[int, int]], keep: set[int] | None = None) -> dict[str, int]:
    """The fixture's own counts over the sites both genotypes are known at."""
    indices = range(len(first)) if keep is None else sorted(keep)
    dosages = [(sum(first[i]), sum(second[i])) for i in indices]
    return {
        "n_sites": len(dosages),
        "n_het_a": sum(x == 1 for x, _ in dosages),
        "n_het_b": sum(y == 1 for _, y in dosages),
        "n_het_both": sum(x == 1 and y == 1 for x, y in dosages),
        "n_ibs0": sum(abs(x - y) == 2 for x, y in dosages),
        "n_ibs2": sum(x == y for x, y in dosages),
    }


def test_king_robust_bounds() -> None:
    # Identical genotypes, no IBS0: exactly one half.
    assert king_robust(100, 100, 100, 0) == 0.5
    assert king_robust(0, 10, 0, 0) is None
    assert band_of(0.5) is RelationBand.SAME_PERSON and band_of(0.25) is RelationBand.FIRST_DEGREE
    assert band_of(0.125) is RelationBand.SECOND_DEGREE and band_of(0.0625) is RelationBand.THIRD_DEGREE
    assert band_of(0.0) is RelationBand.UNRELATED and band_of(-0.2) is RelationBand.UNRELATED
    assert interval_band(0.25, 0.01) is RelationBand.FIRST_DEGREE
    assert interval_band(0.18, 0.01) is RelationBand.UNKNOWN


def test_one_person_as_a_vcf_and_as_an_array(family, tmp_path: Path) -> None:
    """A WGS VCF (absent = reference) and a 23andMe export of one person, with no-calls and missing rows."""
    panel, _, people = family
    vcf = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    no_call, skip = set(range(0, N_SITES, 40)), set(range(7, N_SITES, 50))
    array = load(write_23andme(tmp_path / "a.txt", panel, people["a"], no_call=no_call, skip=skip))
    result = relatedness(vcf, array, panel)
    known = set(range(N_SITES)) - no_call - skip
    assert result.model_dump(include=set(counts_of(people["a"], people["a"], known))) == counts_of(people["a"], people["a"], known)
    assert result.kinship == 0.5 and result.n_ibs0 == 0 and result.ibs2 == 1.0
    assert (result.band, result.reason, result.method) == (RelationBand.SAME_PERSON, None, "king_robust")
    assert same_person(vcf, array, panel) is True


def test_genotyping_errors_lower_kinship_but_keep_one_person(family, tmp_path: Path) -> None:
    panel, _, people = family
    rng = random.Random(SEED + 1)
    noisy = [(1 - h1, h2) if rng.random() < 0.02 else (h1, h2) for h1, h2 in people["a"]]
    result = relatedness(wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"])), load(write_23andme(tmp_path / "n.txt", panel, noisy)), panel)
    expected = counts_of(people["a"], noisy)
    assert result.model_dump(include=set(expected)) == expected
    assert result.kinship == pytest.approx(king_robust(expected["n_het_a"], expected["n_het_b"], expected["n_het_both"], expected["n_ibs0"]))
    assert result.kinship < 0.5 and result.band is RelationBand.SAME_PERSON


@pytest.mark.parametrize(
    ("other", "band", "no_ibs0"),
    [
        ("c", "parent_child", True),  # parent and child: IBS0 is impossible
        ("s", "parent_child", True),  # the child's full sibling is `a`'s child too
        ("g", "second_degree", False),
        ("u", "unrelated", False),
    ],
)
def test_relatives_fall_in_their_band(family, tmp_path: Path, other: str, band: str, no_ibs0: bool) -> None:
    panel, _, people = family
    a = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    b = load(write_23andme(tmp_path / f"{other}.txt", panel, people[other]))
    result = relatedness(a, b, panel)
    expected = counts_of(people["a"], people[other])
    assert result.model_dump(include=set(expected)) == expected
    assert result.band.value == band and result.reason is None
    assert result.kinship_se is not None and 0 < result.kinship_se < 0.03
    if no_ibs0:
        assert result.n_ibs0 == 0
    assert same_person(a, b, panel) is False


def test_full_siblings_share_ibs0_where_parent_and_child_cannot(family, tmp_path: Path) -> None:
    panel, _, people = family
    c, s = (wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x])) for x in ("c", "s"))
    result = relatedness(c, s, panel)
    assert result.band.value == "full_siblings"
    assert result.n_ibs0 == counts_of(people["c"], people["s"])["n_ibs0"] > 0


def _with_errors(haplotypes: list[tuple[int, int]], rate: float, seed: int) -> list[tuple[int, int]]:
    """One haplotype's allele flipped at a share `rate` of sites: genotyping error, drawn with `seed`."""
    rng = random.Random(seed)
    return [(1 - h1, h2) if rng.random() < rate else (h1, h2) for h1, h2 in haplotypes]


@pytest.mark.parametrize(
    ("first", "second", "error", "band", "zone"),
    [
        ("a", "c", 0.0, "parent_child", (0.0, 0.0)),  # IBS0 is impossible between parent and child
        ("a", "s", 0.0, "parent_child", (0.0, 0.0)),
        ("c", "g", 0.0, "parent_child", (0.0, 0.0)),
        ("c", "s", 0.0, "full_siblings", (0.008, 1.0)),  # siblings share no allele IBD at about a quarter of sites
        ("a", "c", 0.01, "parent_child", (0.0, 0.005)),  # genotype error below KING's 0.005: still told
        ("a", "c", 0.035, "first_degree", (0.005, 0.008)),  # inside KING's grey zone: not told, not guessed
        # Error past 0.008 reads as siblings: KING's cutoffs assume array-grade genotypes (docs/RELATEDNESS.md).
        ("a", "c", 0.08, "full_siblings", (0.008, 1.0)),
    ],
)
def test_first_degree_is_split_by_ibs0(
    family, tmp_path: Path, first: str, second: str, error: float, band: str, zone: tuple[float, float]
) -> None:
    """KING's split of the first degree (Manichaikul et al. 2010, section 2.3, cut at KING 2.3.2's 0.005
    and 0.008): the band is named only when the IBS0 share's interval lies on one side of both."""
    panel, _, people = family
    genotypes = _with_errors(people[second], error, SEED + 10) if error else people[second]
    a = wgs(write_vcf(tmp_path / f"{first}.vcf", panel, people[first]))
    b = load(write_23andme(tmp_path / f"{second}.txt", panel, genotypes))
    result = relatedness(a, b, panel)
    expected = counts_of(people[first], genotypes)
    assert result.model_dump(include=set(expected)) == expected
    share = expected["n_ibs0"] / expected["n_sites"]
    assert zone[0] <= share <= zone[1]  # the simulation landed where the row says
    # The kinship alone says first degree; the IBS0 share says which.
    assert interval_band(result.kinship, result.kinship_se) is RelationBand.FIRST_DEGREE
    assert (result.band.value, result.reason) == (band, None)
    assert same_person(a, b, panel) is False


def test_too_few_sites_leave_the_first_degree_unsplit(family, tmp_path: Path) -> None:
    """No IBS0 at 700 sites does not rule out KING's 0.005 (two SEs of a share of 0 out of 700 reach
    0.0057), so parent and child stay `first_degree` there; at the panel's 9,000 they read
    `parent_child`. Full siblings' IBS0 at the same 700 sites is already past 0.008."""
    panel, _, people = family
    small = SitePanel(assembly_label="GRCh37", sites=panel.sites.head(700), source="synthetic")
    a, c, s = (wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x])) for x in "acs")
    thin = relatedness(a, c, small)
    assert (thin.n_sites, thin.n_ibs0) == (700, counts_of(people["a"], people["c"], set(range(700)))["n_ibs0"]) == (700, 0)
    assert interval_band(thin.kinship, thin.kinship_se) is RelationBand.FIRST_DEGREE
    assert (thin.band.value, thin.reason) == ("first_degree", None)
    assert relatedness(a, c, panel).band.value == "parent_child"
    siblings = relatedness(c, s, small)
    assert siblings.n_ibs0 == counts_of(people["c"], people["s"], set(range(700)))["n_ibs0"] > 0
    assert siblings.band.value == "full_siblings"


def test_an_array_absence_is_not_reference(family, tmp_path: Path) -> None:
    """Sites an array leaves out are unknown; the same sites left out of a WGS VCF read as reference."""
    panel, _, people = family
    skip = set(range(0, N_SITES, 3))
    array = load(write_23andme(tmp_path / "a.txt", panel, people["a"], skip=skip))
    vcf = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    result = relatedness(vcf, array, panel)
    assert result.n_sites == N_SITES - len(skip)
    assert relatedness(vcf, vcf, panel).n_sites == N_SITES
    # A VCF whose facet does not say it is a whole genome: its silence stays unknown too.
    unknown = load(tmp_path / "a.vcf")
    assert unknown.facets.select("absence_means").collect().item() == "unknown"
    assert relatedness(unknown, vcf, panel).n_sites == sum(pair != (0, 0) for pair in people["a"])


def test_a_facet_known_only_where_it_carries_the_alt_is_not_banded(family, tmp_path: Path) -> None:
    """A variant-only VCF whose silence says nothing (an exome gVCF with only no-call blocks reads the
    same): its genotype chose the shared sites, which moves the kinship (on PGP's Genos exome, up;
    here, with no direction to it), so no band is named."""
    panel, _, people = family
    variants = load(write_vcf(tmp_path / "a.vcf", panel, people["a"]))  # absence_means unknown
    whole = wgs(write_vcf(tmp_path / "w.vcf", panel, people["a"]))
    stranger = wgs(write_vcf(tmp_path / "u.vcf", panel, people["u"]))
    chosen = relatedness(variants, stranger, panel)
    fair = relatedness(whole, stranger, panel)
    assert (chosen.band, chosen.reason) == (RelationBand.UNKNOWN, VARIANT_SITES_ONLY)
    assert chosen.kinship is not None and chosen.kinship_se is not None
    assert chosen.n_sites == sum(pair != (0, 0) for pair in people["a"])
    # The same person's whole genome on the same sites reads the same: the sites, not the genotypes.
    on_chosen = SitePanel(assembly_label="GRCh37", sites=panel.sites.filter(pl.Series([p != (0, 0) for p in people["a"]])), source="chosen")
    assert relatedness(whole, stranger, on_chosen).kinship == pytest.approx(chosen.kinship)
    assert chosen.kinship != fair.kinship and (fair.band, fair.reason) == (RelationBand.UNRELATED, None)
    # One person still reads as one: the offset is far below the same-person floor.
    assert relatedness(variants, whole, panel).reason == VARIANT_SITES_ONLY
    assert same_person(variants, whole, panel) is True


def test_a_pedigree_pair_read_on_variant_sites_only_cannot_tell(family_on_panel, tmp_path: Path) -> None:
    panel, people = family_on_panel
    pedigree = _family_pedigree([("a", KIN_MOTHER, "c")], "acu")
    child = load(write_vcf(tmp_path / "c.vcf", panel, people["c"]), subject="subject-c")  # variant rows only
    others = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x]), subject=f"subject-{x}") for x in "au"]
    rows = {frozenset((r["individual_id"][-1], r["relative_id"][-1])): r for r in check_pedigree(pedigree, child, *others).iter_rows(named=True)}
    for pair in ("ac", "cu"):
        row = rows[frozenset(pair)]
        assert row["kinship"] is not None and (row["agreement"], row["reason"]) == ("cannot_tell", VARIANT_SITES_ONLY)
    assert rows[frozenset("au")]["agreement"] == "agrees"


def test_unknown_reasons(family, tmp_path: Path) -> None:
    panel, _, people = family
    a = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    b = wgs(write_vcf(tmp_path / "b.vcf", panel, people["b"]))
    relabel = lambda loaded, label: Loaded(  # noqa: E731
        facets=loaded.facets.with_columns(pl.lit(label, dtype=FACETS["assembly_label"]).alias("assembly_label")),
        calls=loaded.calls,
        blocks=loaded.blocks,
    )
    assert relatedness(a, relabel(b, "GRCh38"), panel).reason == ASSEMBLY_MISMATCH
    assert relatedness(a, relabel(b, None), panel).reason == ASSEMBLY_UNKNOWN
    on_38 = SitePanel(assembly_label="GRCh38", sites=panel.sites, source="synthetic")
    assert relatedness(a, b, on_38).reason == NO_PANEL
    few = relatedness(a, b, panel, min_sites=10_000)
    assert (few.band, few.reason, few.kinship, few.n_sites) == (RelationBand.UNKNOWN, TOO_FEW_SITES, None, N_SITES)
    assert same_person(a, relabel(b, "GRCh38"), panel) is None
    for result in (relatedness(a, relabel(b, "GRCh38"), panel), few):
        assert result.band is RelationBand.UNKNOWN and result.kinship is None


def test_an_interval_across_a_boundary_is_unknown(family, tmp_path: Path) -> None:
    """Few sites widen the interval until it crosses a band boundary: unknown, and same_person None."""
    panel, _, people = family
    small = SitePanel(assembly_label="GRCh37", sites=panel.sites.head(60), source="synthetic")
    a = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    rng = random.Random(SEED + 2)
    noisy = [(1 - h1, h2) if rng.random() < 0.15 else (h1, h2) for h1, h2 in people["a"]]
    b = wgs(write_vcf(tmp_path / "n.vcf", panel, noisy))
    result = relatedness(a, b, small, min_sites=10)
    assert result.kinship is not None and result.kinship_se is not None
    low, high = result.kinship - 2 * result.kinship_se, result.kinship + 2 * result.kinship_se
    assert band_of(low) is not band_of(high)
    assert (result.band, result.reason) == (RelationBand.UNKNOWN, BETWEEN_BANDS)
    expected = None if low <= 2**-1.5 < high else (low > 2**-1.5)
    assert same_person(a, b, small, min_sites=10) is expected


def test_pairwise_matches_each_pair(family, tmp_path: Path) -> None:
    panel, _, people = family
    names = ("a", "c", "g", "u")
    loaded = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x])) for x in names]
    frame = pairwise_relatedness(loaded, panel)
    assert frame.schema == RELATEDNESS_FRAME
    assert frame.height == len(names) * (len(names) - 1) // 2
    ids = [item.facets.select("facet_id").collect().item() for item in loaded]
    assert frame.select("facet_a", "facet_b").rows() == [(ids[i], ids[j]) for i in range(4) for j in range(i + 1, 4)]
    for row, (i, j) in zip(frame.iter_rows(named=True), [(i, j) for i in range(4) for j in range(i + 1, 4)]):
        assert Relatedness(**row) == relatedness(loaded[i], loaded[j], panel)


def test_a_person_is_checked_by_same_person_by_default(family, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """RM273's `Person.of` runs `same_person` when no check is named: one person's VCF and array read
    as one person; an unrelated array under the same label reads as two, and is kept, warned."""
    panel, _, people = family
    # The default check reads the package's panel; point it at this fixture's sites instead.
    monkeypatch.setattr(sys.modules["just_vcf.relatedness"], "packaged_panel", lambda assembly: panel if assembly == "GRCh37" else None)
    vcf = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]), subject="p")
    array = load(write_23andme(tmp_path / "a.txt", panel, people["a"]), subject="p")
    stranger = load(write_23andme(tmp_path / "u.txt", panel, people["u"]), subject="p")
    one = Person.of(vcf, array)
    assert (one.checks["outcome"].to_list(), one.one_person) == (["same"], True)
    two = Person.of(vcf, stranger)
    assert (two.checks["outcome"].to_list(), two.one_person) == (["different"], False)
    assert any("read as two people" in warning for warning in two.warnings)


KIN_MOTHER, KIN_FATHER, KIN_TWIN = "KIN:027", "KIN:028", "KIN:010"


def _family_pedigree(links: list[tuple[str, str, str]], people: str) -> Pedigree:
    return pedigree_from_dict(
        {
            "id": "SYNTHETIC",
            "individuals": [{"subject": {"id": f"subject-{p}"}} for p in people],
            "relationships": [
                {"individual_id": f"subject-{a}", "relation": {"id": term}, "relative_id": f"subject-{b}"} for a, term, b in links
            ],
        }
    )


@pytest.fixture
def family_on_panel(family, monkeypatch: pytest.MonkeyPatch):
    """The family's panel as the package's GRCh37 panel, so the default reads it."""
    panel, _, people = family
    monkeypatch.setattr(sys.modules["just_vcf.relatedness"], "packaged_panel", lambda assembly: panel if assembly == "GRCh37" else None)
    return panel, people


def test_a_pedigree_is_checked_against_measured_bands(family_on_panel, tmp_path: Path) -> None:
    """Every relationship the family fixture drew, stated as parent links only, read against KING."""
    panel, people = family_on_panel
    parents = {"c": "ab", "s": "ab", "g": "ce", "h": "ad", "k": "sf"}
    links = [(parent, KIN_MOTHER if i == 0 else KIN_FATHER, child) for child, pair in parents.items() for i, parent in enumerate(pair)]
    pedigree = _family_pedigree(links, "abcsegdhfku")
    with_facets = "abcsgdhku"  # `e` and `f` have none
    loaded = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x]), subject=f"subject-{x}") for x in with_facets]
    before = (pedigree.individuals.clone(), pedigree.relationships.clone())
    frame = check_pedigree(pedigree, *loaded)
    assert frame.schema == PEDIGREE_CHECKS
    assert pedigree.individuals.equals(before[0]) and pedigree.relationships.equals(before[1])
    rows = {frozenset((r["individual_id"][-1], r["relative_id"][-1])): r for r in frame.iter_rows(named=True)}
    assert len(rows) == frame.height == 11 * 10 // 2  # one facet each: one row per pair of individuals
    for pair, row in rows.items():
        if not set(pair) <= set(with_facets):
            assert (row["agreement"], row["reason"], row["facet_a"]) == ("cannot_tell", NO_FACET, None)
            continue
        found = relatedness(*(loaded[with_facets.index(x)] for x in (row["individual_id"][-1], row["relative_id"][-1])), panel)
        assert (row["kinship"], row["band"]) == (found.kinship, found.band.value)
        low, high = (BAND_ORDER.index(band_of(found.kinship + side * 2 * found.kinship_se)) for side in (-1, 1))
        # First degree is two bands at one kinship: the measured one where IBS0 told it, both where not.
        told = {found.band.value} & FIRST_DEGREE_KINDS
        reached = set().union(*(_finest(band.value, told) for band in BAND_ORDER[low : high + 1]))
        allowed = set().union(*(_finest(value) for value in row["expected_bands"]))
        if reached <= allowed:
            assert row["agreement"] == "agrees"
        else:
            assert row["agreement"] == ("disagrees" if not reached & allowed else "cannot_tell")
    # Each kind the fixture drew, and what its band must read.
    for pair, kinship, bands, agreement in (
        ("ac", 0.25, ["parent_child"], "agrees"),  # parent and child
        ("gc", 0.25, ["parent_child"], "agrees"),
        ("cs", 0.25, ["full_siblings"], "agrees"),  # full siblings: both parents listed for both
        ("ch", 0.125, ["second_degree"], "agrees"),  # half siblings
        ("ag", 0.125, ["second_degree"], "agrees"),  # grandparent
        ("sg", 0.125, ["second_degree"], "agrees"),  # aunt or uncle
        ("ab", 0.0, ["unrelated"], "agrees"),  # founders
        ("au", 0.0, ["unrelated"], "agrees"),
    ):
        row = rows[frozenset(pair)]
        assert (row["expected_kinship"], row["expected_bands"], row["agreement"]) == (kinship, bands, agreement), pair
    assert (rows[frozenset("ac")]["band"], rows[frozenset("cs")]["band"]) == ("parent_child", "full_siblings")
    cousins = rows[frozenset("gk")]
    assert cousins["expected_kinship"] == 0.0625 and cousins["agreement"] in ("agrees", "cannot_tell")
    assert not frame.filter(pl.col("agreement") == "disagrees").height


FIRST_DEGREE_KINDS = {"parent_child", "full_siblings"}


def _finest(band: str, told: set[str] = FIRST_DEGREE_KINDS) -> set[str]:
    """The bands `band` stands for: first degree is either kind, or the one IBS0 `told`."""
    return set(told) if band == "first_degree" else {band}


def test_a_parent_who_reads_as_a_sibling_disagrees(family_on_panel, tmp_path: Path) -> None:
    """The pedigree's parent-child edge against its sibling edge: `c` stated as `s`'s mother (they are
    full siblings) and `a` stated as `c`'s full sibling (she is the mother) both read `disagrees`, at a
    kinship that agrees; the true mother link `a`-`s` agrees. On 700 sites a mother link reads
    `cannot_tell`, `between_bands`: no IBS0 there cannot yet rule siblings out."""
    panel, people = family_on_panel
    loaded = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x]), subject=f"subject-{x}") for x in "acs"]
    # One pedigree per statement: together they would make `s` a child of two siblings.
    statements = {
        "cs": (("c", KIN_MOTHER, "s"), ["parent_child"], "full_siblings", "disagrees"),
        "ac": (("a", "KIN:008", "c"), ["full_siblings"], "parent_child", "disagrees"),
        "as": (("a", KIN_MOTHER, "s"), ["parent_child"], "parent_child", "agrees"),
    }
    for pair, (link, bands, band, agreement) in statements.items():
        (row,) = check_pedigree(_family_pedigree([link], pair), *loaded).iter_rows(named=True)
        assert (row["expected_kinship"], row["expected_bands"], row["band"], row["agreement"]) == (0.25, bands, band, agreement), pair
    grouped = people_of(*loaded, check=None, pedigree=_family_pedigree([statements["cs"][0]], "acs"))
    (warning,) = [w for w in grouped["subject-s"].warnings if "pedigree" in w]
    assert "expected parent_child, measured full_siblings" in warning
    small = SitePanel(assembly_label="GRCh37", sites=panel.sites.head(700), source="synthetic")
    thin = check_pedigree(_family_pedigree([("a", KIN_MOTHER, "s")], "as"), *loaded, panel=small).row(0, named=True)
    assert (thin["band"], thin["agreement"], thin["reason"]) == ("first_degree", "cannot_tell", BETWEEN_BANDS)


def test_a_stated_relationship_the_genotypes_contradict_is_kept_and_warned(family_on_panel, tmp_path: Path) -> None:
    """`u` stated as `c`'s father, and an MZ twin stated for `a`'s second file: the false father reads
    `disagrees`, the twin `agrees`; `people` warns naming facet ids, never a subject."""
    panel, people = family_on_panel
    pedigree = _family_pedigree([("a", KIN_MOTHER, "c"), ("u", KIN_FATHER, "c"), ("a", KIN_TWIN, "t")], "acut")
    loaded = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x]), subject=f"subject-{x}") for x in "acu"]
    twin = load(write_23andme(tmp_path / "t.txt", panel, people["a"]), subject="subject-t")
    frame = check_pedigree(pedigree, *loaded, twin)
    rows = {frozenset((r["individual_id"][-1], r["relative_id"][-1])): r for r in frame.iter_rows(named=True)}
    assert (rows[frozenset("uc")]["expected_bands"], rows[frozenset("uc")]["band"], rows[frozenset("uc")]["agreement"]) == (
        ["parent_child"],
        "unrelated",
        "disagrees",
    )
    assert (rows[frozenset("at")]["band"], rows[frozenset("at")]["agreement"]) == ("same_person", "agrees")
    assert rows[frozenset("ac")]["agreement"] == "agrees"
    # Grouped with the pedigree: the warning sits on both people of the pair and names their facets.
    grouped = people_of(*loaded, twin, pedigree=pedigree)
    ids = {x: item.facets.select("facet_id").collect().item() for x, item in zip("acu", loaded)}
    for subject in ("subject-c", "subject-u"):
        (warning,) = [w for w in grouped[subject].warnings if "pedigree" in w]
        assert f"{ids['u']}/{ids['c']}" in warning or f"{ids['c']}/{ids['u']}" in warning
    assert not [w for w in grouped["subject-a"].warnings if "pedigree" in w]
    assert all("subject-" not in w for person in grouped.values() for w in person.warnings)
    # The pedigree check is its own request: with no subject check, it still runs.
    unchecked = people_of(*loaded, twin, check=None, pedigree=pedigree)
    assert any("pedigree" in w for w in unchecked["subject-u"].warnings)


def test_a_person_given_as_a_person_reads_as_its_loads(family_on_panel, tmp_path: Path) -> None:
    panel, people = family_on_panel
    pedigree = _family_pedigree([("a", KIN_MOTHER, "c")], "ac")
    loaded = [wgs(write_vcf(tmp_path / f"{x}.vcf", panel, people[x]), subject=f"subject-{x}") for x in "ac"]
    assert check_pedigree(pedigree, *people_of(*loaded, check=None).values()).equals(check_pedigree(pedigree, *loaded))


def test_a_panel_refuses_sites_it_cannot_read_alike() -> None:
    bad = pl.DataFrame(
        {
            "chrom": ["1", "X", "2", "3", "3"],
            "pos": pl.Series([10, 20, 30, 40, 40], dtype=pl.UInt32),
            "ref_allele": ["A", "C", "A", "G", "G"],
            "effect_allele": ["G", "T", "T", "AT", "A"],
        }
    )
    with pytest.raises(ValueError) as refused:
        SitePanel(assembly_label="GRCh37", sites=bad, source="synthetic")
    message = str(refused.value)
    for part in ("1 sites not on an autosome", "1 sites not a one-base ACGT SNV", "1 sites strand-ambiguous", "2 sites two sites at one position"):
        assert part in message


def test_the_packaged_panels_are_valid() -> None:
    assert packaged_assemblies() == ("GRCh37", "GRCh38")
    sizes = {}
    for assembly in packaged_assemblies():
        panel = packaged_panel(assembly)
        assert panel is not None and panel.assembly_label == assembly
        sizes[assembly] = panel.sites.height
        # Panel order is the canonical order, which the jackknife's contiguous blocks rely on.
        ranks = panel.sites.select(pl.col("chrom").cast(pl.Int32), "pos")
        assert ranks.equals(ranks.sort("chrom", "pos"))
    assert min(sizes.values()) > 10_000
    assert packaged_panel("T2T-CHM13") is None


def test_the_packaged_panels_are_derived_from_somaliers_sites(test_reference_dir: Path) -> None:
    """Opt-in: `JUST_VCF_RELATEDNESS_SITES_DIR` holds somalier's two site VCFs (docs/PROVENANCE.md)."""
    if "JUST_VCF_RELATEDNESS_SITES_DIR" in os.environ:
        value = os.environ["JUST_VCF_RELATEDNESS_SITES_DIR"]
    else:
        value = dotenv_values(find_dotenv(usecwd=True)).get("JUST_VCF_RELATEDNESS_SITES_DIR")
    if not value or not Path(value).is_dir():
        pytest.skip("JUST_VCF_RELATEDNESS_SITES_DIR is not set or absent")
    relatedness_sites = pytest.importorskip("relatedness_sites")
    store = ReferenceStore.open(test_reference_dir)
    for assembly, name in (("GRCh37", "sites.GRCh37.vcf.gz"), ("GRCh38", "sites.hg38.vcf.gz")):
        frame, _ = relatedness_sites.derive(relatedness_sites.read_sites(Path(value) / name), store, assembly)
        packaged = packaged_panel(assembly).sites.rename({"ref_allele": "ref", "effect_allele": "alt"})
        assert frame.equals(packaged)


def _first_of(participant: Path, format_id: str) -> Path:
    """The participant's first file of `format_id` on GRCh37."""
    found = next((p for p in sorted(participant.rglob("*")) if p.is_file() and (s := sniff(p)).format == format_id and s.build == "GRCh37"), None)
    if found is None:
        pytest.skip(f"no {format_id} file for {participant.name}")
    return found


def test_pgp_arrays_of_one_person_and_of_two(pgp_participants: Path) -> None:
    """Opt-in, PGP-HMS (public): one participant's 23andMe and AncestryDNA exports read as one person;
    the 23andMe exports of two participants as unrelated (docs/RELATEDNESS.md has the measured ranges)."""
    one = load(_first_of(pgp_participants / "hu832966", "23andme"))
    same = load(_first_of(pgp_participants / "hu832966", "ancestrydna"))
    other = load(_first_of(pgp_participants / "huF31554", "23andme"))
    together = relatedness(one, same)
    assert together.band is RelationBand.SAME_PERSON and together.n_sites > 2_000
    assert together.n_ibs0 <= together.n_sites // 500
    apart = relatedness(one, other)
    assert apart.band is RelationBand.UNRELATED and apart.kinship is not None and apart.kinship < THIRD_DEGREE_FLOOR
    assert same_person(one, same) is True and same_person(one, other) is False


def test_minted_vrs_ids_are_compared_where_both_facets_carry_them(tmp_path: Path) -> None:
    """Two facets of one person on a small labelled reference, VRS ids minted on both: every site where
    both have an ALT row is compared and none differ; an id changed on one facet is counted."""
    rng = random.Random(SEED + 3)
    sequence = "".join(rng.choice("ACGT") for _ in range(6000))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">1\n" + "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
    store = ReferenceStore.open(tmp_path / "store")
    store.add_fasta(fasta, assembly_label="GRCh37")
    positions = list(range(20, 6000, 10))
    ambiguous = ("AT", "TA", "CG", "GC")
    alts = [next(b for b in "ACGT" if b != sequence[p - 1] and sequence[p - 1] + b not in ambiguous) for p in positions]
    sites = pl.DataFrame(
        {
            "chrom": ["1"] * len(positions),
            "pos": pl.Series(positions, dtype=pl.UInt32),
            "ref_allele": [sequence[p - 1] for p in positions],
            "effect_allele": alts,
        }
    )
    panel = SitePanel(assembly_label="GRCh37", sites=sites, source="synthetic")
    haplotypes = [(int(rng.random() < 0.4), int(rng.random() < 0.4)) for _ in positions]
    header = ["##fileformat=VCFv4.2", f"##contig=<ID=1,length={len(sequence)}>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"]))
    body = [f"1\t{p}\t.\t{sequence[p - 1]}\t{a}\t.\tPASS\t.\tGT\t{GT_TEXT[h]}" for p, a, h in zip(positions, alts, haplotypes) if h != (0, 0)]
    minted = []
    for name in ("x.vcf", "y.vcf"):
        sample_line = header[-1].replace("S1", "S2") if name == "y.vcf" else header[-1]
        (tmp_path / name).write_text("\n".join([*header[:-1], sample_line, *body]) + "\n")
        loaded = mint_vrs_ids(resolve_ref(wgs(tmp_path / name), store), store)
        # A contig this short gives no length-read build; the label the store asserts is the build.
        labelled = loaded.facets.with_columns(pl.lit("GRCh37", dtype=FACETS["assembly_label"]).alias("assembly_label"))
        minted.append(Loaded(facets=labelled, calls=loaded.calls, blocks=loaded.blocks))
    first, second = minted
    assert first.facets.select("facet_id").collect().item() != second.facets.select("facet_id").collect().item()
    carriers = [p for p, h in zip(positions, haplotypes) if h != (0, 0)]
    result = relatedness(first, second, panel, min_sites=10)
    assert (result.n_sites, result.n_vrs_compared, result.n_vrs_differ) == (len(positions), len(carriers), 0)
    changed_ids = pl.when(pl.col("pos") == carriers[0]).then(pl.col("vrs_allele_ids").list.eval(pl.element() + "x"))
    changed = Loaded(
        facets=second.facets,
        calls=second.calls.with_columns(changed_ids.otherwise(pl.col("vrs_allele_ids")).alias("vrs_allele_ids")),
        blocks=second.blocks,
    )
    assert relatedness(first, changed, panel, min_sites=10).n_vrs_differ == 1


def test_the_result_is_the_same_bits_every_run(family, tmp_path: Path) -> None:
    """Charter 12: the same frames give the same result to the bit. (It passed before the jackknife's
    blocks were sorted too: group_by's order did not vary here, so it guards, it does not prove.)"""
    panel, _, people = family
    a = wgs(write_vcf(tmp_path / "a.vcf", panel, people["a"]))
    g = load(write_23andme(tmp_path / "g.txt", panel, people["g"]))
    first = relatedness(a, g, panel)
    assert all(relatedness(a, g, panel) == first for _ in range(20))
