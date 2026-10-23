"""The method-error atlas (RM279): synthetic VCF pairs on a small refget store, synthetic region sets.

Every expected count comes from the fixture's own case list, which states what each site is; rates
are checked against textbook Wilson intervals and against `concordance(..., region=)` run per class.
"""

import gzip
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceStore, load, resolve_ref
from just_vcf.atlas import (
    ALL,
    ATLAS,
    COUNTS,
    GIAB_V36_FILES,
    RegionSet,
    giab_stratifications,
    method_atlas,
    pair_counts,
    pool_counts,
    read_bed,
    wilson_interval,
)
from just_vcf.concordance import concordance
from just_vcf.loading import Loaded

RNG = random.Random(279)
SEQ = {"1": "".join(RNG.choice("ACGT") for _ in range(2000)), "2": "".join(RNG.choice("ACGT") for _ in range(1000))}


def base(chrom: str, pos: int, length: int = 1) -> str:
    return SEQ[chrom][pos - 1 : pos - 1 + length]


def other(b: str) -> str:
    return next(x for x in "ACGT" if x != b)


def write_vcf(path: Path, source: str, records: list[tuple]) -> Path:
    """`records` are `(chrom, pos, ref, alt, gt, filter, extra FORMAT fields)`."""
    keys = sorted({key for *_, extra in records for key in extra if key != "END"})
    meta = ["##fileformat=VCFv4.3", f"##source={source}", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items())]
    meta += [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
        '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posteriors">',
    ]
    body = []
    for chrom, pos, ref, alt, gt, filt, extra in records:
        info = f"END={extra['END']}" if "END" in extra else "."
        fields = ["GT", *keys]
        values = [gt, *(extra.get(key, ".") for key in keys)]
        body.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\t{':'.join(fields)}\t{':'.join(values)}")
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


@pytest.fixture
def store(tmp_path: Path) -> ReferenceStore:
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in SEQ.items()))
    s = ReferenceStore.open(tmp_path / "store")
    s.add_fasta(fasta)
    return s


def facet(tmp_path: Path, name: str, source: str, records: list[tuple], store: ReferenceStore) -> Loaded:
    """A facet labelled GRCh38, the assembly the synthetic regions say they are on."""
    loaded = resolve_ref(load(write_vcf(tmp_path / f"{name}.vcf", source, records)), store)
    label = pl.lit("GRCh38").cast(loaded.facets.collect_schema()["assembly_label"]).alias("assembly_label")
    return Loaded(facets=loaded.facets.with_columns(label), calls=loaded.calls, blocks=loaded.blocks)


# One site per case: (a's genotype or None for no row, b's, a's FILTER, what the site is).
# What it is: an Agreement value, or "not_carried" for a site neither genotype holds the allele in.
CASES = [
    ("0/1", "0/1", "PASS", "match"),
    ("1/1", "1/1", "PASS", "match"),
    ("0/1", "1/1", "PASS", "mismatch"),
    ("0/0", "0/1", "PASS", "mismatch"),  # carried by b only, both called
    ("1", "1/1", "PASS", "ploidy_differs"),
    ("1/.", "1/1", "PASS", "uncalled"),
    ("0/1", None, "PASS", "only_a"),
    (None, "0/1", "PASS", "only_b"),
    ("0/1", "0/1", "LowQual", "match"),  # left out of carried_pass
    ("0/0", None, "PASS", "not_carried"),  # a stated allele no genotype holds: not a method's error
    ("./.", "./.", "PASS", "not_carried"),
    ("0/1", "1/1", "LowQual", "mismatch"),
]
POSITIONS = [40 + 60 * i for i in range(len(CASES))]
REGIONS = {
    "early": pl.DataFrame({"chrom": ["chr1"], "pos": [1], "end_pos": [POSITIONS[5]]}),
    "late": pl.DataFrame({"chrom": ["1", "1"], "pos": [POSITIONS[4], POSITIONS[8] + 1], "end_pos": [POSITIONS[7], 2000]}),
    "empty": pl.DataFrame({"chrom": ["2"], "pos": [1], "end_pos": [100]}),
}


def records(side: int) -> list[tuple]:
    rows = []
    for pos, case in zip(POSITIONS, CASES):
        gt = case[side]
        if gt is None:
            continue
        ref = base("1", pos)
        rows.append(("1", pos, ref, other(ref), gt, case[2] if side == 0 else "PASS", {}))
    return rows


@pytest.fixture
def pair(tmp_path: Path, store: ReferenceStore) -> tuple[Loaded, Loaded]:
    a = facet(tmp_path, "a", "freeBayes v1.3.6", records(0), store)
    b = facet(tmp_path, "b", "HaplotypeCaller", records(1), store)
    return a, b


def expected(label: str, scope: str) -> dict[str, int]:
    """Counts of each kind of site in region class `label`, read off `CASES`."""
    inside = []
    for pos, case in zip(POSITIONS, CASES):
        frame = REGIONS.get(label)
        held = label == ALL or any(c in ("1", "chr1") and s <= pos <= e for c, s, e in frame.iter_rows())
        if held and case[3] != "not_carried" and (scope == "carried" or case[2] == "PASS"):
            inside.append(case[3])
    return {kind: inside.count(kind) for kind in ("match", "mismatch", "ploidy_differs", "uncalled", "only_a", "only_b")} | {"alleles": len(inside)}


def test_counts_follow_the_cases_per_region_and_scope(pair: tuple[Loaded, Loaded]) -> None:
    regions = RegionSet.from_frames("test set", "GRCh38", REGIONS)
    atlas = method_atlas([pair], regions)
    assert atlas.schema == ATLAS
    assert not {"chrom", "pos", "ref", "alt", "vrs_allele_id"} & set(atlas.columns)  # never a per-site row
    snv = atlas.filter(pl.col("variant_class") == "snv")
    assert set(snv["region_class"].to_list()) == {ALL, "early", "late", "empty"}
    for row in snv.iter_rows(named=True):
        want = expected(row["region_class"], row["scope"])
        assert {kind: row[kind] for kind in want} == want, (row["region_class"], row["scope"])
        assert row["discordant"] == want["mismatch"] + want["ploidy_differs"]
        assert row["called_both"] == want["match"] + want["mismatch"] + want["ploidy_differs"]
    # The other classes have no allele: their rates are unknown, never zero.
    others = atlas.filter(pl.col("variant_class") != "snv")
    assert others["alleles"].sum() == 0
    assert others["discordance_rate"].null_count() == others.height == others["one_sided_ci_high"].null_count()
    empty = snv.filter(pl.col("region_class") == "empty")
    assert empty["alleles"].to_list() == [0, 0] and empty["discordance_rate"].null_count() == 2


def test_the_methods_are_put_in_order_whichever_facet_comes_first(pair: tuple[Loaded, Loaded]) -> None:
    regions = RegionSet.from_frames("test set", "GRCh38", REGIONS)
    a, b = pair
    forward = pair_counts(a, b, regions)
    backward = pair_counts(b, a, regions)
    assert forward.schema == COUNTS
    assert forward.equals(backward)
    # "GATK HaplotypeCaller" sorts before "freeBayes": a is GATK, so a freeBayes-only allele is only_b.
    row = forward.filter((pl.col("region_class") == ALL) & (pl.col("variant_class") == "snv") & (pl.col("scope") == "carried")).row(0, named=True)
    assert (row["a_caller"], row["b_caller"], row["b_caller_version"]) == ("GATK HaplotypeCaller", "freeBayes", "1.3.6")
    assert (row["only_a"], row["only_b"]) == (1, 1)
    cases_only_in_freebayes = sum(1 for case in CASES if case[1] is None and case[3] != "not_carried")
    assert row["only_b"] == cases_only_in_freebayes


def test_one_comparison_labelled_per_class_equals_concordance_cut_per_class(pair: tuple[Loaded, Loaded]) -> None:
    """Without an indel spelled at two positions, placing each allele of one comparison in its classes
    gives the counts `concordance(..., region=)` gives class by class."""
    regions = RegionSet.from_frames("test set", "GRCh38", REGIONS)
    a, b = pair
    counts = pair_counts(a, b, regions).filter((pl.col("scope") == "carried") & (pl.col("variant_class") == "snv"))
    for label, frame in REGIONS.items():
        sites = concordance(a, b, region=frame).sites.collect()
        carried = sites.filter(pl.col("a_gt_alleles").list.contains(pl.col("alt")).fill_null(False) | pl.col("b_gt_alleles").list.contains(pl.col("alt")).fill_null(False))
        by_agreement = dict(carried.group_by(pl.col("agreement").cast(pl.String)).len().iter_rows())
        row = counts.filter(pl.col("region_class") == label).row(0, named=True)
        assert {kind: row[kind] for kind in by_agreement} == by_agreement
        assert row["alleles"] == carried.height


def test_counting_contig_by_contig_gives_the_same_counts(tmp_path: Path, store: ReferenceStore) -> None:
    def snv(chrom: str, pos: int, gt: str) -> tuple:
        ref = base(chrom, pos)
        return (chrom, pos, ref, other(ref), gt, "PASS", {})

    a = facet(tmp_path, "a", "freeBayes v1.3.6", [snv("1", 100, "0/1"), snv("1", 200, "1/1"), snv("2", 50, "0/1"), snv("2", 70, "0/1")], store)
    b = facet(tmp_path, "b", "HaplotypeCaller", [snv("1", 100, "0/1"), snv("1", 300, "0/1"), snv("2", 50, "1/1")], store)
    regions = RegionSet.from_frames("test set", "GRCh38", REGIONS)
    whole = pair_counts(a, b, regions)
    chunked = pair_counts(a, b, regions, by_contig=True)
    assert chunked.equals(whole)
    row = whole.filter((pl.col("region_class") == ALL) & (pl.col("variant_class") == "snv") & (pl.col("scope") == "carried")).row(0, named=True)
    assert (row["alleles"], row["match"], row["mismatch"], row["only_a"] + row["only_b"]) == (5, 1, 1, 3)


def test_wilson_intervals_are_the_textbook_ones() -> None:
    # Newcombe (1998), Statistics in Medicine 17:857, table I, method 3 (the score interval without
    # continuity correction), to the four decimals printed there.
    table = {(81, 263): (0.2553, 0.3662), (15, 148): (0.0624, 0.1605), (0, 20): (0.0, 0.1611), (1, 29): (0.0061, 0.1718)}
    for (k, n), interval in table.items():
        assert wilson_interval(k, n) == pytest.approx(interval, abs=5e-5)
    # The interval is symmetric: k of n failures mirror k of n successes.
    low, high = wilson_interval(0, 20)
    assert wilson_interval(20, 20) == pytest.approx((1 - high, 1 - low))
    assert wilson_interval(0, 0) is None


def test_pooling_sums_counts_and_takes_rates_from_the_sums(pair: tuple[Loaded, Loaded]) -> None:
    regions = RegionSet.from_frames("test set", "GRCh38", REGIONS)
    one = pair_counts(*pair, regions)
    pooled = pool_counts([one, one])
    single = pool_counts([one])
    assert pooled["pairs"].unique().to_list() == [2]
    assert pooled["alleles"].to_list() == [2 * n for n in single["alleles"].to_list()]
    row = pooled.filter((pl.col("region_class") == ALL) & (pl.col("variant_class") == "snv") & (pl.col("scope") == "carried")).row(0, named=True)
    assert row["discordance_rate"] == pytest.approx(row["discordant"] / row["called_both"])
    low, high = wilson_interval(row["discordant"], row["called_both"])
    assert (row["discordance_ci_low"], row["discordance_ci_high"]) == pytest.approx((low, high))
    low, high = wilson_interval(row["only_a"] + row["only_b"], row["alleles"])
    assert (row["one_sided_ci_low"], row["one_sided_ci_high"]) == pytest.approx((low, high))
    # Another method pair pools into its own rows.
    renamed = one.with_columns(pl.lit("DeepVariant").alias("b_caller"))
    assert pool_counts([one, renamed]).height == 2 * single.height
    assert pool_counts([]).schema == ATLAS


def test_p_same_is_pooled_per_basis(tmp_path: Path, store: ReferenceStore) -> None:
    def snv(pos: int, gt: str, **extra: str) -> tuple:
        ref = base("1", pos)
        return ("1", pos, ref, other(ref), gt, "PASS", extra)

    a = facet(tmp_path, "a", "freeBayes v1.3.6", [snv(100, "0/1", PL="20,0,200"), snv(200, "0/1", PL="0,10,90"), snv(300, "0/1", GP="0.1,0.8,0.1"), snv(400, "0/1", PL="5,0,50")], store)
    b = facet(tmp_path, "b", "HaplotypeCaller", [snv(100, "0/1", PL="30,0,50"), snv(200, "1/1", PL="90,10,0"), snv(300, "0/1", GP="0.2,0.6,0.2"), snv(400, "0/1", GP="0.1,0.8,0.1")], store)
    sites = concordance(a, b).sites.collect()
    atlas = method_atlas([(a, b)], RegionSet.from_frames("test set", "GRCh38", {}))
    row = atlas.filter((pl.col("variant_class") == "snv") & (pl.col("scope") == "carried")).row(0, named=True)
    gl = sites.filter(pl.col("basis") == "gl_flat_prior")["p_same"]
    gp = sites.filter(pl.col("basis") == "gp")["p_same"]
    assert (row["graded_gl"], row["graded_gp"]) == (gl.len(), gp.len()) == (2, 1)
    assert row["p_same_gl_mean"] == pytest.approx(gl.mean())
    assert row["p_same_gp_mean"] == pytest.approx(gp.mean())
    assert row["alleles"] == 4  # the likelihood-against-posterior site counts, ungraded


def test_a_one_sided_allele_the_silent_facet_reads_as_reference_is_counted(tmp_path: Path, store: ReferenceStore) -> None:
    ref = base("1", 100)
    a = facet(tmp_path, "a", "freeBayes v1.3.6", [("1", 100, ref, other(ref), "0/1", "PASS", {}), ("1", 500, base("1", 500), other(base("1", 500)), "0/1", "PASS", {})], store)
    blocks = [("1", 90, base("1", 90), "<*>", "0/0", "PASS", {"END": "150"}), ("1", 800, base("1", 800), other(base("1", 800)), "0/1", "PASS", {})]
    b = facet(tmp_path, "b", "HaplotypeCaller", blocks, store)
    atlas = method_atlas([(a, b)], RegionSet.from_frames("test set", "GRCh38", {}))
    row = atlas.filter((pl.col("variant_class") == "snv") & (pl.col("scope") == "carried")).row(0, named=True)
    # 100 lies in b's ref block; 500 lies in nothing b wrote, and b's absence says nothing.
    assert (row["one_sided_reference"], row["only_a"] + row["only_b"]) == (1, 3)


def test_regions_on_another_assembly_are_refused(pair: tuple[Loaded, Loaded]) -> None:
    with pytest.raises(ValueError, match="lift one first"):
        method_atlas([pair], RegionSet.from_frames("test set", "GRCh37", REGIONS))
    with pytest.raises(ValueError, match="no region"):
        RegionSet.from_frames("test set", "GRCh38", {ALL: REGIONS["early"]})


def test_giab_stratifications_are_found_by_name_and_read_one_based(tmp_path: Path) -> None:
    nested = tmp_path / "GRCh38@all" / "Union"
    nested.mkdir(parents=True)
    bed = nested / GIAB_V36_FILES["easy"].format(asm="GRCh38")
    with gzip.open(bed, "wt") as out:
        out.write("#chrom\tstart\tend\nchr1\t0\t10\nchr1\t10\t20\nchr2\t99\t100\textra\n")
    read = read_bed(bed).collect()
    assert read.rows() == [("chr1", 1, 10), ("chr1", 11, 20), ("chr2", 100, 100)]
    regions = giab_stratifications(tmp_path, "GRCh38", ["easy"])
    assert regions.classes == (ALL, "easy")
    assert regions.name == "GIAB v3.6 GRCh38"
    # Touching BED intervals are one closed interval.
    assert regions.intervals["easy"].rows() == [("1", 1, 20), ("2", 100, 100)]
    with pytest.raises(FileNotFoundError, match="GRCh38_alldifficultregions"):
        giab_stratifications(tmp_path, "GRCh38", ["easy", "difficult"])
    with pytest.raises(ValueError, match="unknown GIAB region class"):
        giab_stratifications(tmp_path, "GRCh38", ["easy", "hard"])
    with pytest.raises(ValueError, match="GRCh37, GRCh38"):
        giab_stratifications(tmp_path, "T2T-CHM13")
