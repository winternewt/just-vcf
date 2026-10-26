"""The AADR Human Origins panel builder (RM5), on synthetic EIGENSTRAT files written here.

Every individual, site and genotype is made up with a fixed seed; nothing comes from AADR. The fixture
writes the same genotypes in all three `.geno` layouts (transposed packed, packed, text), with a
`.snp`, `.ind` and `.anno` beside them and a small GRCh37-labelled reference. Each expectation is taken
from how the fixture was made: which individuals it meant to exclude and why, which sites it broke and
how, and the ALT dosages it drew before encoding them as the file's counts of its first allele.

The opt-in test at the end builds from the real AADR files and places one public PGP-HMS genome.
"""

import json
import math
import time
from pathlib import Path

import polars as pl
import pytest

np = pytest.importorskip("numpy")

from just_vcf import load  # noqa: E402
from just_vcf.ancestry import (  # noqa: E402
    AADR_REGIONS,
    INDIVIDUAL_EXCLUSIONS,
    SITE_DROPS,
    WEST_EURASIA,
    build_aadr_panel,
    continental_of,
    infer_ancestry,
    load_pca_panel,
    region_of,
)
from just_vcf.ancestry.eigenstrat import geno_layout, read_anno, read_genotypes, read_ind, read_snp  # noqa: E402
from just_vcf.ancestry.pca import panel_paths  # noqa: E402
from just_vcf.assembly import ASSEMBLY_LENGTHS  # noqa: E402
from just_vcf.reference import ReferenceStore  # noqa: E402
from just_vcf.schema import CONTIG_RANK  # noqa: E402

SEED = 20261026
CONTIG_LENGTH = 40_000
SITES_PER_CONTIG = 300
N_PER_POPULATION = 40
DIM_REF = 4
MISSING_RATE = 0.02
BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
GT_TEXT = {0: "0/0", 1: "0/1", 2: "1/1"}

POPULATIONS = {
    # group: (political entity, latitude, longitude, the region it should be filed under)
    "Alpha": ("France", 46.0, 2.0, "Europe"),
    "Beta": ("Pakistan", 33.0, 70.0, "Southern Asia"),
    "Gamma": ("Algeria", 35.0, 3.0, "Northern Africa"),
}

ANNO_HEADER = [
    "Genetic ID (suffices: made up for a test)",
    "Persistent Genetic ID",
    "Individual ID",
    '"Date mean in BP in years before 1950 CE [made up]"',
    "Group ID",
    "Political Entity",
    "Latitude",
    "Longitude",
    "Data type",
    "ASSESSMENT",
    "ASSESSMENT WARNINGS: the prefix of the column before it",
]


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{name}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n" for name, seq in seqs.items()))
    return path


class Fixture:
    """One synthetic AADR release: what was written, and what a build of it must give."""

    def __init__(self, root: Path) -> None:
        rng = np.random.default_rng(SEED)
        self.root = root
        self.seqs = {c: "".join(rng.choice(list(BASES), CONTIG_LENGTH)) for c in ("1", "2")}

        # Sites: (snp id, chrom, pos, allele1, allele2, kind). `allele1` is the allele the genotypes count.
        sites: list[tuple[str, str, int, str, str, str]] = []
        for chrom in ("1", "2"):
            for k in range(SITES_PER_CONTIG):
                pos = 1_000 + 100 * k
                ref = self.seqs[chrom][pos - 1]
                alt = next(b for b in BASES[rng.integers(0, 4) :] + BASES if b not in (ref, COMPLEMENT[ref]))
                kind = "swapped" if k % 7 == 3 else "kept"
                a1, a2 = (alt, ref) if kind == "swapped" else (ref, alt)
                sites.append((f"s{chrom}_{k}", chrom, pos, a1, a2, kind))
        base = lambda chrom, pos: self.seqs[chrom][pos - 1]  # noqa: E731
        p = 35_123
        sites.append(("x1", "23", 5_000, "A", "G", "not on an autosome"))
        sites.append(("amb1", "1", p, base("1", p), COMPLEMENT[base("1", p)], "strand-ambiguous (A/T or C/G)"))
        sites.append(("indel1", "1", p + 10, "I", "D", "not a one-base ACGT SNV"))
        sites.append(("dup1", "2", p, base("2", p), next(b for b in BASES if b not in (base("2", p), COMPLEMENT[base("2", p)])), "two records at one position"))
        sites.append(("dup2", "2", p, base("2", p), next(b for b in BASES if b not in (base("2", p), COMPLEMENT[base("2", p)])), "two records at one position"))
        q = p + 20
        neither = [b for b in BASES if b != base("1", q)]
        a1, a2 = next((x, y) for x in neither for y in neither if x != y and COMPLEMENT[x] != y)
        sites.append(("neither1", "1", q, a1, a2, "neither allele is the reference base"))
        sites.append(("chr3", "3", 7_000, "C", "T", "contig not in the reference"))
        r = p + 40
        sites.append(("mono1", "1", r, base("1", r), next(b for b in BASES if b not in (base("1", r), COMPLEMENT[base("1", r)])), "minor allele frequency below the floor"))
        self.sites = sites
        n_sites = len(sites)

        # Individuals, in .ind order: (genetic id, individual id, group, entity, lat, lon, data type, assessment, date, why excluded).
        people: list[tuple[str, str, str, str, str, str, str, str, str, str | None]] = []
        people.append(("Alpha0.DG", "IAlpha0", "Alpha", "France", "46.0", "2.0", "Shotgun.diploid", "Pass", "0", "second representation of one individual"))
        assessments = ("Pass", "PROVISIONAL_PASS", "MERGE_PASS")
        for group, (entity, lat, lon, _) in POPULATIONS.items():
            for k in range(N_PER_POPULATION):
                people.append((f"{group}{k}.HO", f"I{group}{k}", group, entity, f"{lat + 0.01 * k}", f"{lon}", "HumanOrigins", assessments[k % 3], "0", None))
        people += [
            ("Anc1.AG", "IAnc1", "Alpha_Neolithic", "France", "46.0", "2.0", "1240k", "Pass", "4500", "ancient (date mean in BP not 0)"),
            ("Sg1.SG", "ISg1", "Alpha", "France", "46.0", "2.0", "Shotgun", "Pass", "0", "not diploid data (HumanOrigins or Shotgun.diploid)"),
            ("Ref1.REF", "IRef1", "Chimp", "..", "..", "..", "Reference Genome", "Pass", "0", "not diploid data (HumanOrigins or Shotgun.diploid)"),
            ("Q1.HO", "IQ1", "Alpha", "France", "46.0", "2.0", "HumanOrigins", "Questionable", "0", "assessment not a pass"),
            ("Fo1.HO", "IFo1", "Alpha-o", "France", "46.0", "2.0", "HumanOrigins", "Pass", "0", "group label flagged (Ignore_ or a -qualifier)"),
            ("Ig1.HO", "IIg1", "Ignore_Alpha", "France", "46.0", "2.0", "HumanOrigins", "Pass", "0", "group label flagged (Ignore_ or a -qualifier)"),
            ("Nl1.HO", "INl1", "Alpha", "France", "..", "..", "HumanOrigins", "Pass", "0", "no latitude or longitude"),
            ("Far1.HO", "IFar1", "Delta", "Nigeria", "10.0", "2.0", "HumanOrigins", "Pass", "0", "outside the region box"),
            ("Us1.HO", "IUs1", "Alpha", "USA", "40.0", "-20.0", "HumanOrigins", "Pass", "0", "political entity has no ancestry region"),
        ]
        self.people = people
        self.not_in_anno = "Missing1.HO"
        ind_rows = [(gid, group) for gid, _, group, *_ in people] + [(self.not_in_anno, "Alpha")]
        n_ind = len(ind_rows)

        # ALT dosage (against the reference base) per site and individual, before encoding.
        frequencies = {}
        shared = rng.uniform(0.15, 0.85, n_sites)
        for group in POPULATIONS:
            redrawn = rng.random(n_sites) < 0.4
            frequencies[group] = np.where(redrawn, rng.uniform(0.15, 0.85, n_sites), shared)
        self.frequencies = frequencies
        alt = rng.integers(0, 3, (n_sites, n_ind)).astype(np.int8)
        for j, (gid, *_rest) in enumerate(people):
            group = people[j][2]
            if group in frequencies:
                alt[:, j] = rng.binomial(2, frequencies[group])
        mono = next(i for i, s in enumerate(sites) if s[0] == "mono1")
        alt[mono, :] = 0
        missing = rng.random((n_sites, n_ind)) < MISSING_RATE
        self.alt = alt
        self.missing = missing

        counted = alt.copy()
        for i, site in enumerate(sites):
            if site[5] != "swapped":
                counted[i] = 2 - alt[i]
        counted[missing] = -9
        self.counted = counted
        self.selected = [j for j, person in enumerate(people) if person[9] is None]

        self.prefix = root / "synthetic.aadr.patch.PUB"
        self.anno = root / "synthetic.aadr.PUB.anno"
        (root / "synthetic.aadr.patch.PUB.ind").write_text("".join(f"{gid:>20} U {group}\n" for gid, group in ind_rows))
        (root / "synthetic.aadr.patch.PUB.snp").write_text(
            "".join(f"{sid:>20} {chrom:>5} {0.0:15.6f} {pos:15d} {a1} {a2}\n" for sid, chrom, pos, a1, a2, _ in sites)
        )
        lines = ["\t".join(ANNO_HEADER)]
        for gid, iid, group, entity, lat, lon, dtype, assess, date, _ in people:
            lines.append("\t".join([gid, gid, iid, date, group, entity, lat, lon, dtype, assess, ".."]))
        self.anno.write_text("\n".join(lines) + "\n")
        self.fasta = write_fasta(root / "grch37.fa", self.seqs)

    def write_geno(self, layout: str, path: Path) -> Path:
        """The genotypes in one `.geno` layout: `tgeno`, `geno` or `text`."""
        n_snp, n_ind = self.counted.shape
        codes = np.where(self.counted == -9, 3, self.counted).astype(np.uint8)
        if layout == "text":
            text = np.where(self.counted == -9, 9, self.counted)
            path.write_text("".join("".join(str(v) for v in row) + "\n" for row in text))
            return path

        def pack(matrix: np.ndarray, width: int) -> np.ndarray:
            records = np.zeros((matrix.shape[0], width), dtype=np.uint8)
            for j in range(matrix.shape[1]):
                records[:, j // 4] |= matrix[:, j] << (6 - 2 * (j % 4))
            return records

        if layout == "tgeno":
            header = f"TGENO {n_ind:7d} {n_snp:7d} {0:x} {0:x}".encode().ljust(48, b"\0")
            body = pack(codes.T, math.ceil(n_snp / 4))
        else:
            rlen = max(48, math.ceil(n_ind / 4))
            header = f"GENO {n_ind:7d} {n_snp:7d} {0:x} {0:x}".encode().ljust(rlen, b"\0")
            body = pack(codes, rlen)
        path.write_bytes(header + body.tobytes())
        return path

    def store(self, root: Path) -> ReferenceStore:
        store = ReferenceStore.open(root / "store")
        store.add_fasta(self.fasta, assembly_label="GRCh37")
        return store

    def kept_sites(self) -> list[tuple[str, int, str, str, int]]:
        """`(chrom, pos, reference base, other allele, site index)` of every site a build keeps."""
        kept = [(s[1], s[2], s[4] if s[5] == "swapped" else s[3], s[3] if s[5] == "swapped" else s[4], i) for i, s in enumerate(self.sites) if s[5] in ("kept", "swapped")]
        return sorted(kept, key=lambda k: (CONTIG_RANK[k[0]], k[1]))


@pytest.fixture(scope="module")
def fixture(tmp_path_factory: pytest.TempPathFactory) -> Fixture:
    return Fixture(tmp_path_factory.mktemp("aadr"))


@pytest.fixture(scope="module")
def built(fixture: Fixture, tmp_path_factory: pytest.TempPathFactory) -> dict:
    """One build per `.geno` layout, from the same fixture and store."""
    store = fixture.store(tmp_path_factory.mktemp("store"))
    out = {}
    for layout in ("tgeno", "geno", "text"):
        fixture.write_geno(layout, Path(f"{fixture.prefix}.geno"))
        panel_dir = tmp_path_factory.mktemp(f"panel-{layout}")
        out[layout] = (build_aadr_panel(fixture.prefix, fixture.anno, panel_dir, reference=store, dim_ref=DIM_REF), panel_dir)
    return out


def test_each_layout_is_recognised_and_reads_the_counts_it_was_written_with(fixture: Fixture, tmp_path: Path) -> None:
    n_snp, n_ind = fixture.counted.shape
    pick = np.array([0, 5, n_ind - 1, 17])
    for layout in ("tgeno", "geno", "text"):
        path = fixture.write_geno(layout, tmp_path / f"{layout}.geno")
        assert geno_layout(path) == layout
        assert np.array_equal(read_genotypes(path, n_ind, n_snp, pick), fixture.counted[:, pick])


def test_a_geno_file_of_the_wrong_size_or_counts_is_refused(fixture: Fixture, tmp_path: Path) -> None:
    n_snp, n_ind = fixture.counted.shape
    path = fixture.write_geno("tgeno", tmp_path / "t.geno")
    with pytest.raises(ValueError, match="individuals"):
        read_genotypes(path, n_ind + 1, n_snp, np.array([0]))
    path.write_bytes(path.read_bytes()[:-1])
    with pytest.raises(ValueError, match="bytes"):
        read_genotypes(path, n_ind, n_snp, np.array([0]))
    text = fixture.write_geno("text", tmp_path / "x.geno")
    with pytest.raises(ValueError, match="line"):
        read_genotypes(text, n_ind + 1, n_snp, np.array([0]))
    bogus = tmp_path / "b.geno"
    bogus.write_bytes(b"\x1f\x8b\x08\x00" + b"\0" * 100)
    with pytest.raises(ValueError, match="layout"):
        geno_layout(bogus)


def test_the_readers_read_what_was_written(fixture: Fixture) -> None:
    ind = read_ind(Path(f"{fixture.prefix}.ind"))
    assert ind["genetic_id"].to_list() == [p[0] for p in fixture.people] + [fixture.not_in_anno]
    snp = read_snp(Path(f"{fixture.prefix}.snp"))
    assert list(zip(snp["chrom"], snp["pos"], snp["allele1"], snp["allele2"], strict=True)) == [(s[1], s[2], s[3], s[4]) for s in fixture.sites]
    anno = read_anno(fixture.anno)
    assert anno["genetic_id"].to_list() == [p[0] for p in fixture.people]
    assert anno["assessment"].to_list() == [p[7] for p in fixture.people], "ASSESSMENT, not the column it prefixes"
    assert anno.filter(pl.col("latitude").is_null())["genetic_id"].to_list() == ["Ref1.REF", "Nl1.HO"]


def test_an_annotation_file_without_a_column_it_needs_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "short.anno"
    path.write_text("\t".join(h for h in ANNO_HEADER if not h.startswith("Political")) + "\n")
    with pytest.raises(ValueError, match="Political Entity"):
        read_anno(path)


def test_the_three_layouts_build_the_same_panel(built: dict) -> None:
    frames = {layout: (pl.read_parquet(build.panel.sites_path), pl.read_parquet(build.panel.refpcs_path)) for layout, (build, _) in built.items()}
    first = frames["tgeno"]
    for layout in ("geno", "text"):
        assert frames[layout][0].equals(first[0]), layout
        assert frames[layout][1].equals(first[1]), layout


def test_individuals_are_kept_and_excluded_by_the_stated_rules(fixture: Fixture, built: dict) -> None:
    build, _ = built["tgeno"]
    expected: dict[str, int] = dict.fromkeys(INDIVIDUAL_EXCLUSIONS, 0)
    for person in fixture.people:
        if person[9] is not None:
            expected[person[9]] += 1
    expected["not in the annotation file"] += 1
    assert build.excluded == expected
    refpcs = pl.read_parquet(build.panel.refpcs_path)
    chosen = [fixture.people[j] for j in fixture.selected]
    assert refpcs["iid"].to_list() == [p[0] for p in chosen], "the HumanOrigins representation, in .ind order"
    assert refpcs["population"].to_list() == [p[2] for p in chosen]
    assert refpcs["superpop"].to_list() == [POPULATIONS[p[2]][3] for p in chosen]
    assert build.n_selected == len(chosen)
    assert build.regions == {region: N_PER_POPULATION for *_, region in POPULATIONS.values()}


def test_sites_are_dropped_by_the_stated_rules_and_oriented_on_the_reference(fixture: Fixture, built: dict) -> None:
    build, _ = built["tgeno"]
    expected: dict[str, int] = dict.fromkeys(SITE_DROPS, 0)
    for site in fixture.sites:
        if site[5] not in ("kept", "swapped"):
            expected[site[5]] += 1
    assert build.dropped == expected
    assert build.n_swapped == sum(1 for s in fixture.sites if s[5] == "swapped")
    sites = pl.read_parquet(build.panel.sites_path)
    kept = fixture.kept_sites()
    assert list(zip(sites["chrom"], sites["pos"], sites["ref"], sites["alt"], strict=True)) == [k[:4] for k in kept]
    assert all(fixture.seqs[c][p - 1] == r for c, p, r in zip(sites["chrom"], sites["pos"], sites["ref"], strict=True))
    assert sites["effect_allele"].equals(sites["alt"], check_names=False)


def test_the_panel_mean_is_the_drawn_alt_dosage_whatever_allele_the_file_counted(fixture: Fixture, built: dict) -> None:
    """The flip from the file's counts to ALT dosage is undone per site where `.snp` lists the
    reference base second: the panel's per-site mean equals the mean of the ALT dosages drawn."""
    build, _ = built["tgeno"]
    sites = pl.read_parquet(build.panel.sites_path)
    rows = [k[4] for k in fixture.kept_sites()]
    alt = fixture.alt[np.ix_(rows, fixture.selected)].astype(float)
    alt[fixture.missing[np.ix_(rows, fixture.selected)]] = np.nan
    assert np.allclose(sites["mean"].to_numpy(), np.nanmean(alt, axis=1))


def test_the_meta_records_the_source_and_the_rules(fixture: Fixture, built: dict) -> None:
    build, _ = built["tgeno"]
    meta = json.loads(build.panel.meta_path.read_text())
    assert meta["build"] == "GRCh37" and meta["panel"] == "aadr_ho"
    assert meta["source"]["geno"] == "synthetic.aadr.patch.PUB.geno"
    assert meta["source"]["geno_layout"] == "tgeno"
    assert meta["source"]["anno"] == fixture.anno.name
    assert meta["selection"]["excluded"] == build.excluded
    assert meta["selection"]["box"] == list(WEST_EURASIA)
    assert meta["sites"]["dropped"] == build.dropped


def test_a_person_drawn_from_one_population_is_placed_in_it(fixture: Fixture, built: dict, tmp_path: Path) -> None:
    _, panel_dir = built["tgeno"]
    rng = np.random.default_rng(SEED + 1)
    lengths = ASSEMBLY_LENGTHS["GRCh37"]
    meta = ["##fileformat=VCFv4.2", *(f"##contig=<ID={c},length={next(iter(lengths[c]))}>" for c in ("1", "2", "3"))]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    for group, (_, _, _, region) in POPULATIONS.items():
        body = []
        for chrom, pos, ref, alt, i in fixture.kept_sites():
            dosage = int(rng.binomial(2, fixture.frequencies[group][i]))
            body.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{GT_TEXT[dosage]}")
        path = tmp_path / f"{group}.vcf"
        path.write_text("\n".join([*meta, header, *body]) + "\n")
        loaded = load(path)
        broad = infer_ancestry(loaded, panel_dir, "aadr_ho")
        assert broad.estimate == region, (group, broad.reason, broad.probabilities)
        assert broad.n_sites_used == len(fixture.kept_sites())
        fine = infer_ancestry(loaded, panel_dir, "aadr_ho", resolution="population")
        assert fine.estimate == group, (group, fine.probabilities)


def test_a_store_without_grch37_is_refused(fixture: Fixture, tmp_path: Path) -> None:
    fixture.write_geno("tgeno", Path(f"{fixture.prefix}.geno"))
    store = ReferenceStore.open(tmp_path / "empty")
    with pytest.raises(ValueError, match="GRCh37"):
        build_aadr_panel(fixture.prefix, fixture.anno, tmp_path / "panel", reference=store, dim_ref=DIM_REF)


def test_too_few_individuals_for_the_requested_components_is_refused(fixture: Fixture, tmp_path: Path) -> None:
    fixture.write_geno("tgeno", Path(f"{fixture.prefix}.geno"))
    with pytest.raises(ValueError, match="individuals"):
        build_aadr_panel(fixture.prefix, fixture.anno, tmp_path / "panel", reference=fixture.store(tmp_path), dim_ref=20)


def test_regions_follow_the_political_entity_and_russia_splits_at_the_urals() -> None:
    assert region_of("France", 2.0) == "Europe"
    assert region_of("Russia", 37.6) == "Europe"
    assert region_of("Russia", 85.0) == "Siberia"
    assert region_of("Iran", 51.0) == "Southern Asia"
    assert region_of("USA", -100.0) is None
    assert region_of("Atlantis", 0.0) is None


def test_only_the_regions_with_a_continental_counterpart_fold_onto_one() -> None:
    """Folding goes through `continental_of` on the region name; the Asian regions, Oceania and
    Siberia withhold rather than be filed under a 1000 Genomes group they only partly match."""
    folds = {region: continental_of(region) for region in sorted(set(AADR_REGIONS.values()) - {None} | {"Siberia"})}
    assert folds == {
        "Central Asia": None,
        "Eastern Asia": None,
        "Europe": "EUR",
        "Latin America and the Caribbean": "AMR",
        "Northern Africa": "AFR",
        "Oceania": None,
        "Siberia": None,
        "South-eastern Asia": None,
        "Southern Asia": None,
        "Sub-Saharan Africa": "AFR",
        "Western Asia": None,
    }


# --- opt-in: the real AADR release ------------------------------------------------------------------


HU6ABACE_SNV = Path("pgp-hms/participants/hu6ABACE/upload-3885/60820188474283.snp.vcf.gz")
"""PGP-HMS hu6ABACE's DRAGEN 3.2.8 SNV genome (GRCh37, hs37d5), consented-public (docs/PGP_HMS.md)."""


def test_the_real_release_builds_a_panel_that_places_a_public_genome(
    ancestry_dir: Path, samples_dir: Path, test_reference_dir: Path, tmp_path: Path
) -> None:
    """Build from AADR's `compatibility_HO` files under `JUST_VCF_ANCESTRY_DIR` (top level or one
    directory down), then place hu6ABACE. Counts and the group only."""
    genos = sorted([*ancestry_dir.glob("*compatibility_HO.aadr.patch.PUB.geno"), *ancestry_dir.glob("*/*compatibility_HO.aadr.patch.PUB.geno")])
    if not genos:
        pytest.skip(f"no AADR compatibility_HO .geno under {ancestry_dir}")
    genome = samples_dir / HU6ABACE_SNV
    if not genome.exists():
        pytest.skip(f"{genome} is absent")
    prefix = genos[-1].with_suffix("")
    anno = prefix.parent / prefix.name.replace(".patch", "")
    anno = anno.with_name(anno.name + ".anno")
    started = time.monotonic()
    build = build_aadr_panel(prefix, anno, tmp_path, reference=ReferenceStore.open(test_reference_dir))
    built_in = time.monotonic() - started
    model = load_pca_panel(tmp_path, "aadr_ho", "GRCh37")
    assert model.sites.height == build.n_sites == build.panel.n_sites
    assert build.n_sites + sum(build.dropped.values()) == build.n_snps
    assert build.n_selected + sum(build.excluded.values()) == build.n_individuals
    assert build.dropped["strand-ambiguous (A/T or C/G)"] == 0, "the README: Human Origins has no A/T or C/G SNPs"
    assert sum(build.regions.values()) == build.n_selected
    started = time.monotonic()
    loaded = load(genome)
    result = infer_ancestry(loaded, tmp_path, "aadr_ho")
    placed_in = time.monotonic() - started
    assert result, (result.reason, result.detail)
    assert result.n_sites_used + result.n_sites_unusable == result.n_sites_panel == build.n_sites
    assert result.estimate in build.regions
    fine = infer_ancestry(loaded, tmp_path, "aadr_ho", resolution="population")
    assert fine.estimate in set(model.populations.tolist())
    print(
        f"\naadr_ho: {build.n_selected} of {build.n_individuals} individuals, {build.n_sites} of {build.n_snps} sites "
        f"({build.n_swapped} swapped), LOO {build.panel.loo_accuracy:.3f}, built in {built_in:.0f} s; "
        f"hu6ABACE {result.n_sites_used} used ({result.n_sites_observed} observed, {result.n_sites_absent_as_reference} absent as reference, "
        f"{result.n_sites_unusable} unusable), {result.estimate} (margin {result.margin:.2f}), population {fine.estimate} "
        f"(margin {fine.margin:.2f}), loaded and placed in {placed_in:.0f} s"
    )
