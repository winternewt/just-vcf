"""The consensus posterior (RM280): synthetic VCF facets on a small refget store, atlases built from
stated counts.

Every expected rate comes from the counts the test writes into its atlas, and every expected posterior
from the textbook product `prior · ∏ L_f(k)` worked out here in plain Python with the rule the module
docstring states, `L(k) = (1 - e)·q(k) + e·(1 - q(k))/P`.
"""

import math
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import ReferenceMismatchError, ReferenceStore, genotype_evidence, link_disagreements, load, normalize, resolve_ref
from just_vcf.atlas import CI_Z, RegionSet
from just_vcf.consensus import CONSENSUS_SUPPORT, FACET_RATES, METHOD_RATES, Consensus, GpReading, RateBasis, consensus, method_error_rates
from just_vcf.links import CONSENSUS
from just_vcf.loading import Loaded
from just_vcf.person import Person
from just_vcf.schema import CALLS, FACETS

RNG = random.Random(280)
SEQ = {"1": "".join(RNG.choice("ACGT") for _ in range(2000)), "2": "".join(RNG.choice("ACGT") for _ in range(1000))}
N = 10_000_000
"""Alleles behind every atlas row here: large enough that the Wilson floor (about 4e-7) is far below
every rate the tests state."""


def base(chrom: str, pos: int) -> str:
    return SEQ[chrom][pos - 1]


def other(b: str, skip: int = 0) -> str:
    return [x for x in "ACGT" if x != b][skip]


def write_vcf(path: Path, records: list[tuple]) -> Path:
    """`records` are `(chrom, pos, ref, alt, gt, filter, extra FORMAT fields)`."""
    keys = sorted({key for *_, extra in records for key in extra if key != "END"})
    # The file's own name in its header, so facets with the same rows are different files (`facet_id`
    # digests the bytes).
    meta = ["##fileformat=VCFv4.3", f"##test_facet={path.stem}", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items())]
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
        values = [gt, *(extra.get(key, ".") for key in keys)]
        body.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\t{':'.join(['GT', *keys])}\t{':'.join(values)}")
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


def facet(tmp_path: Path, name: str, caller: str | None, version: str | None, records: list[tuple], store: ReferenceStore, **columns: object) -> Loaded:
    """A facet with REF read, labelled GRCh38, its caller and any other `facets` column as given."""
    loaded = resolve_ref(load(write_vcf(tmp_path / f"{name}.vcf", records)), store)
    schema = loaded.facets.collect_schema()
    values = {"assembly_label": "GRCh38", "caller": caller, "caller_version": version, **columns}
    return Loaded(
        facets=loaded.facets.with_columns(*[pl.lit(value, dtype=schema[key]).alias(key) for key, value in values.items()]),
        calls=loaded.calls,
        blocks=loaded.blocks,
        links=loaded.links,
    )


def snv(chrom: str, pos: int, gt: str, filt: str = "PASS", alt: int = 0, **extra: str) -> tuple:
    return (chrom, pos, base(chrom, pos), other(base(chrom, pos), alt), gt, filt, extra)


def atlas_row(a: tuple[str, str | None], b: tuple[str, str | None], *, d: float, region: str = "all", vc: str = "snv",
              only_a: int = 0, only_b: int = 0, n: int = N, ploidy: int = 0) -> dict:  # fmt: skip
    """An atlas row of scope `carried_pass` whose genotype mismatch is `d` over `n` alleles both called
    at one ploidy, beside `ploidy` alleles called at two."""
    return {
        "scope": "carried_pass", "a_caller": a[0], "a_caller_version": a[1], "b_caller": b[0], "b_caller_version": b[1],
        "region_class": region, "variant_class": vc, "match": n - round(d * n), "mismatch": round(d * n),
        "ploidy_differs": ploidy, "called_both": n + ploidy, "discordant": round(d * n) + ploidy,
        "alleles": n + ploidy + only_a + only_b, "only_a": only_a, "only_b": only_b,
    }  # fmt: skip


def atlas(*rows: dict) -> pl.DataFrame:
    counts = ("match", "mismatch", "ploidy_differs", "called_both", "discordant", "alleles", "only_a", "only_b")
    return pl.DataFrame(list(rows)).with_columns(pl.col(name).cast(pl.UInt64) for name in counts)


def floor(n: int) -> float:
    return CI_Z * CI_Z / (n + CI_Z * CI_Z)


def every_class(a: tuple, b: tuple, d: float, **kw: object) -> list[dict]:
    return [atlas_row(a, b, d=d, vc=vc, **kw) for vc in ("snv", "mnv", "indel")]


X, Y, Z = ("X", "1"), ("Y", "1"), ("Z", "1")
D_XY, D_XZ, D_YZ = 0.010, 0.030, 0.034
"""Three methods on the same reads: X errs least, Z most. Triangles give e = 0.003, 0.007, 0.027."""
TRIANGLE = atlas(*every_class(X, Y, D_XY), *every_class(X, Z, D_XZ), *every_class(Y, Z, D_YZ))
E = {"X": (D_XY + D_XZ - D_YZ) / 2, "Y": (D_XY + D_YZ - D_XZ) / 2, "Z": (D_XZ + D_YZ - D_XY) / 2}


def posterior(likelihoods: list[list[float]]) -> list[float]:
    """The flat-prior posterior from per-facet likelihood vectors."""
    logs = [sum(math.log(f[k]) for f in likelihoods) for k in range(len(likelihoods[0]))]
    top = max(logs)
    w = [math.exp(v - top) for v in logs]
    return [x / sum(w) for x in w]


def hard(copies: int, e: float, ploidy: int = 2) -> list[float]:
    return [1 - e if k == copies else e / ploidy for k in range(ploidy + 1)]


def tempered(q: list[float], e: float) -> list[float]:
    p = len(q) - 1
    return [(1 - e) * qk + e * (1 - qk) / p for qk in q]


def flat(pl_values: list[int]) -> list[float]:
    w = [10 ** (-v / 10) for v in pl_values]
    return [x / sum(w) for x in w]


def rows_of(result: Consensus) -> list[dict]:
    calls = result.calls.collect()
    support = result.support.collect()
    assert calls.height == support.height
    return [{**c, **{f"s_{k}": v for k, v in s.items()}} for c, s in zip(calls.iter_rows(named=True), support.iter_rows(named=True))]


# --------------------------------------------------------------------------------------------------
# Rates


def test_method_rates_triangulate_three_methods() -> None:
    rates = method_error_rates(TRIANGLE)
    assert rates.schema == METHOD_RATES
    snv_all = rates.filter((pl.col("region_class") == "all") & (pl.col("variant_class") == "snv"))
    got = {r["caller"]: r for r in snv_all.iter_rows(named=True)}
    assert set(got) == {"X", "Y", "Z"}
    for caller, expected in E.items():
        assert got[caller]["genotype_attribution"] == "triangle"
        assert got[caller]["genotype_error"] == pytest.approx(expected, rel=1e-6)
        assert got[caller]["partners"] == 2


def test_method_rates_halve_a_lone_pair_and_floor_a_zero() -> None:
    lone = atlas(atlas_row(X, Y, d=0.02), atlas_row(X, Z, d=0.0, n=1000))
    rates = {r["caller"]: r for r in method_error_rates(lone).iter_rows(named=True)}
    # X has two partners but no triangle (Y and Z were never compared): the mean of its two halves.
    assert rates["X"]["genotype_attribution"] == "half_pair"
    assert rates["X"]["genotype_error"] == pytest.approx((0.01 + floor(1000)) / 2, rel=1e-9)
    assert rates["Y"]["genotype_error"] == pytest.approx(0.01, rel=1e-9)
    # Zero discordance in 1,000 alleles is not a zero rate.
    assert rates["Z"]["genotype_error"] == pytest.approx(floor(1000), rel=1e-9)


def test_silence_rate_is_directional_and_the_lowest_over_partners() -> None:
    # X is silent on 300 of Y's carried alleles and 50 of Z's; Y and Z on 100 of X's.
    rows = atlas(atlas_row(X, Y, d=0.01, only_a=100, only_b=300), atlas_row(X, Z, d=0.01, only_a=100, only_b=50))
    rates = {r["caller"]: r for r in method_error_rates(rows).iter_rows(named=True)}
    assert rates["X"]["silence_error"] == pytest.approx(min(300 / (N + 300), 50 / (N + 50)), rel=1e-9)
    assert rates["Y"]["silence_error"] == pytest.approx(100 / (N + 100), rel=1e-9)


def test_rates_by_caller_pool_versions_and_drop_self_pairs() -> None:
    rows = atlas(atlas_row(("X", "1"), ("Y", "1"), d=0.02), atlas_row(("X", "2"), ("Y", "1"), d=0.04), atlas_row(("X", "1"), ("X", "2"), d=0.5))
    by_caller = {r["caller"]: r for r in method_error_rates(rows, by="caller").iter_rows(named=True)}
    assert set(by_caller) == {"X", "Y"}
    assert by_caller["X"]["caller_version"] is None
    assert by_caller["X"]["genotype_error"] == pytest.approx(0.03 / 2, rel=1e-6)


def test_a_ploidy_difference_is_no_genotype_error() -> None:
    # A male X: half the alleles both called differ only in ploidy. The atlas's discordance counts
    # them; the genotype rate does not.
    rows = atlas(atlas_row(X, Y, d=0.01, ploidy=N))
    assert rows["discordant"][0] / rows["called_both"][0] > 0.5
    rates = {r["caller"]: r for r in method_error_rates(rows).iter_rows(named=True)}
    assert rates["X"]["genotype_error"] == pytest.approx(0.005, rel=1e-6)
    assert rates["X"]["compared"] == N


def test_rates_keep_only_the_asked_scope() -> None:
    rows = TRIANGLE.with_columns(pl.lit("carried").alias("scope"))
    assert method_error_rates(rows).is_empty()
    assert not method_error_rates(rows, scope="carried").is_empty()


# --------------------------------------------------------------------------------------------------
# The posterior


def three(tmp_path: Path, store: ReferenceStore, gts: tuple[str, str, str], **extra: dict) -> list[Loaded]:
    return [
        facet(tmp_path, name, caller, "1", [snv("1", 100, gt, **extra.get(name, {}))], store)
        for name, caller, gt in zip(("x", "y", "z"), ("X", "Y", "Z"), gts)
    ]


def test_three_agreeing_facets_multiply(tmp_path: Path, store: ReferenceStore) -> None:
    result = consensus(three(tmp_path, store, ("0/1", "0/1", "0/1")), atlas=TRIANGLE)
    (row,) = rows_of(result)
    expected = posterior([hard(1, E[c]) for c in "XYZ"])
    assert row["gp"] == pytest.approx(expected, rel=1e-5)
    assert row["gt"] == [0, 1] and row["s_copies"] == 1
    assert (row["s_supported"], row["s_contradicted"], row["s_silent"], row["s_facets"]) == (3, 0, 0, 3)
    assert row["s_confidence"] == pytest.approx(-10 * math.log10(1 - expected[1]), rel=1e-9)
    assert row["alts"] == [other(base("1", 100))] and row["ref"] == base("1", 100)


def test_two_against_one_where_the_lone_facet_is_better(tmp_path: Path, store: ReferenceStore) -> None:
    # X (error 0.003) says het; Y (0.007) and Z (0.027) say hom-alt. Two noisier facets still win.
    result = consensus(three(tmp_path, store, ("0/1", "1/1", "1/1")), atlas=TRIANGLE)
    (row,) = rows_of(result)
    expected = posterior([hard(1, E["X"]), hard(2, E["Y"]), hard(2, E["Z"])])
    assert row["gp"] == pytest.approx(expected, rel=1e-5)
    assert row["s_copies"] == max(range(3), key=lambda k: expected[k])
    assert (row["s_supported"], row["s_contradicted"]) == (2, 1)


def test_a_graded_facet_is_tempered_by_its_rate(tmp_path: Path, store: ReferenceStore) -> None:
    # X's likelihoods put 30 phred on het over either homozygote; Y and Z say hom-alt.
    pls = [300, 0, 300]
    result = consensus(three(tmp_path, store, ("0/1", "1/1", "1/1"), x={"PL": ",".join(map(str, pls))}), atlas=TRIANGLE)
    (row,) = rows_of(result)
    expected = posterior([tempered(flat(pls), E["X"]), hard(2, E["Y"]), hard(2, E["Z"])])
    assert row["gp"] == pytest.approx(expected, rel=1e-5)
    # Tempered, X can say no more than its rate allows, and the two others carry the site. Untempered,
    # its likelihoods alone would overrule them both.
    untempered = posterior([flat(pls), hard(2, E["Y"]), hard(2, E["Z"])])
    assert row["s_copies"] == 2 == max(range(3), key=lambda k: expected[k])
    assert max(range(3), key=lambda k: untempered[k]) == 1


def test_disagreement_is_kept_in_the_parents(tmp_path: Path, store: ReferenceStore) -> None:
    parents = three(tmp_path, store, ("0/1", "1/1", "1/1"))
    before = [p.calls.collect() for p in parents]
    result = consensus(parents, atlas=TRIANGLE)
    assert [p.calls.collect().equals(b) for p, b in zip(parents, before)] == [True, True, True]
    assert result.facets.collect()["relation"].to_list() == [CONSENSUS]
    assert result.facets.collect()["derived_from"].to_list() == [None]


def test_a_failed_filter_is_silent_unless_asked(tmp_path: Path, store: ReferenceStore) -> None:
    parents = [
        facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store),
        facet(tmp_path, "y", "Y", "1", [snv("1", 100, "1/1", filt="LowQual")], store),
    ]
    (row,) = rows_of(consensus(parents, atlas=TRIANGLE))
    assert (row["s_supported"], row["s_silent"]) == (1, 1)
    assert row["gp"] == pytest.approx(posterior([hard(1, E["X"])]), rel=1e-5)
    carried = TRIANGLE.with_columns(pl.lit("carried").alias("scope"))
    (row,) = rows_of(consensus(parents, atlas=carried, pass_only=False))
    assert row["s_silent"] == 0
    assert row["gp"] == pytest.approx(posterior([hard(1, E["X"]), hard(2, E["Y"])]), rel=1e-5)


def test_silence_reads_as_reference_only_where_absence_says_so(tmp_path: Path, store: ReferenceStore) -> None:
    silence = atlas(*every_class(X, Y, 0.01, only_a=2000, only_b=4000), *every_class(X, Z, 0.01), *every_class(Y, Z, 0.01))
    rates = {r["caller"]: r for r in method_error_rates(silence).filter(pl.col("variant_class") == "snv").iter_rows(named=True)}
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1"), snv("1", 300, "0/1")], store)
    wgs = {"absence_means": "weak_reference", "variant_classes": [{"variant_class": "snv", "alleles": 2}]}
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 300, "0/1"), snv("1", 500, "0/1")], store, **wgs)
    (first, *_rest) = rows_of(consensus([x, y], atlas=silence))
    assert first["pos"] == 100
    assert (first["s_supported"], first["s_contradicted"], first["s_silent"]) == (1, 1, 0)
    expected = posterior([hard(1, rates["X"]["genotype_error"]), hard(0, rates["Y"]["silence_error"])])
    assert first["gp"] == pytest.approx(expected, rel=1e-5)
    # A content filter holds the silence back: the position is then silent on Y's side.
    cut = facet(tmp_path, "y2", "Y", "1", [snv("1", 300, "0/1"), snv("1", 500, "0/1")], store, **wgs, content_filter="snv_only")
    (first, *_rest) = rows_of(consensus([x, cut], atlas=silence))
    assert (first["s_supported"], first["s_silent"]) == (1, 1)
    # An array-like facet whose absence is unknown says nothing.
    unknown = facet(tmp_path, "y3", "Y", "1", [snv("1", 300, "0/1")], store)
    (first, *_rest) = rows_of(consensus([x, unknown], atlas=silence))
    assert first["s_silent"] == 1


def test_a_reference_block_speaks_with_its_likelihoods(tmp_path: Path, store: ReferenceStore) -> None:
    pls = [0, 30, 300]
    block = ("1", 90, base("1", 90), "<*>", "0/0", ".", {"END": "120", "PL": ",".join(map(str, pls))})
    plain = ("1", 150, base("1", 150), "<*>", "0/0", ".", {"END": "160"})
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1"), snv("1", 155, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [block, plain], store)
    rows = rows_of(consensus([x, y], atlas=TRIANGLE))
    rates = {r["caller"]: r for r in method_error_rates(TRIANGLE).filter(pl.col("variant_class") == "snv").iter_rows(named=True)}
    at = {row["pos"]: row for row in rows}
    graded = posterior([hard(1, rates["X"]["genotype_error"]), tempered(flat(pls), rates["Y"]["silence_error"])])
    assert at[100]["gp"] == pytest.approx(graded, rel=1e-5)
    asserted = posterior([hard(1, rates["X"]["genotype_error"]), hard(0, rates["Y"]["silence_error"])])
    assert at[155]["gp"] == pytest.approx(asserted, rel=1e-5)


def test_a_site_row_without_the_allele_reads_as_zero_copies(tmp_path: Path, store: ReferenceStore) -> None:
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1", alt=1)], store)
    rows = rows_of(consensus([x, y], atlas=TRIANGLE))
    assert len(rows) == 2  # two ALTs, one row each
    for row in rows:
        assert (row["s_supported"], row["s_contradicted"]) == (1, 1)


def test_an_unmeasured_method_gets_the_worst_measured_rate_and_is_counted(tmp_path: Path, store: ReferenceStore) -> None:
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    w = facet(tmp_path, "w", "W", "9", [snv("1", 100, "0/1")], store)
    result = consensus([x, w], atlas=TRIANGLE)
    (row,) = rows_of(result)
    assert row["s_defaulted"] == 1
    worst = max(E.values())
    assert row["gp"] == pytest.approx(posterior([hard(1, E["X"]), hard(1, worst)]), rel=1e-5)
    basis = result.rates.filter(pl.col("caller") == "W")["rate_basis"].unique().to_list()
    assert basis == [RateBasis.DEFAULT.value]
    assert result.rates.schema == FACET_RATES
    stated = consensus([x, w], atlas=TRIANGLE, default_rate=0.2)
    (row,) = rows_of(stated)
    assert row["gp"] == pytest.approx(posterior([hard(1, E["X"]), hard(1, 0.2)]), rel=1e-5)
    with pytest.raises(ValueError, match="default_rate"):
        consensus([x, w], atlas=None)


def test_another_version_falls_back_to_the_caller(tmp_path: Path, store: ReferenceStore) -> None:
    x = facet(tmp_path, "x", "X", "7", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1")], store)
    result = consensus([x, y], atlas=TRIANGLE)
    assert set(result.rates.filter(pl.col("caller") == "X")["rate_basis"].to_list()) == {RateBasis.CALLER.value}


def test_regions_weight_a_facet_by_its_worst_class_there(tmp_path: Path, store: ReferenceStore) -> None:
    hard_rows = [atlas_row(a, b, d=d * 10, region="difficult", vc=vc) for a, b, d in ((X, Y, D_XY), (X, Z, D_XZ), (Y, Z, D_YZ)) for vc in ("snv", "mnv", "indel")]
    easy_rows = [atlas_row(a, b, d=d / 10, region="easy", vc=vc) for a, b, d in ((X, Y, D_XY), (X, Z, D_XZ), (Y, Z, D_YZ)) for vc in ("snv", "mnv", "indel")]
    table = pl.concat([TRIANGLE, atlas(*hard_rows, *easy_rows)])
    regions = RegionSet.from_frames("test", "GRCh38", {
        "difficult": pl.DataFrame({"chrom": ["1"], "pos": [50], "end_pos": [150]}),
        "easy": pl.DataFrame({"chrom": ["1", "1"], "pos": [1, 151], "end_pos": [120, 2000]}),
    })  # fmt: skip
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1"), snv("1", 300, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "1/1"), snv("1", 300, "1/1")], store)
    at = {row["pos"]: row for row in rows_of(consensus([x, y], atlas=table, regions=regions))}
    # 100 is in both classes: the difficult rate (the higher) weighs it. 300 is easy only.
    assert at[100]["gp"] == pytest.approx(posterior([hard(1, E["X"] * 10), hard(2, E["Y"] * 10)]), rel=1e-4)
    assert at[300]["gp"] == pytest.approx(posterior([hard(1, E["X"] / 10), hard(2, E["Y"] / 10)]), rel=1e-4)


def test_a_posterior_is_read_as_its_call_unless_asked(tmp_path: Path, store: ReferenceStore) -> None:
    gp = [0.05, 0.9, 0.05]
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1", GP=",".join(map(str, gp)))], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "1/1")], store)
    (row,) = rows_of(consensus([x, y], atlas=TRIANGLE))
    assert row["gp"] == pytest.approx(posterior([hard(1, E["X"]), hard(2, E["Y"])]), rel=1e-5)
    assert row["s_gp_as_likelihood"] == 0
    (row,) = rows_of(consensus([x, y], atlas=TRIANGLE, gp_as=GpReading.LIKELIHOOD))
    assert row["gp"] == pytest.approx(posterior([tempered(gp, E["X"]), hard(2, E["Y"])]), rel=1e-5)
    assert row["s_gp_as_likelihood"] == 1


# --------------------------------------------------------------------------------------------------
# The facet


def test_the_consensus_is_a_facet_linked_to_every_parent(tmp_path: Path, store: ReferenceStore) -> None:
    parents = three(tmp_path, store, ("0/1", "1/1", "1/1"))
    result = consensus(parents, atlas=TRIANGLE)
    assert isinstance(result, Loaded)
    facets = result.facets.collect()
    assert facets.schema == FACETS and result.calls.collect().schema == CALLS
    assert result.support.collect().schema == CONSENSUS_SUPPORT
    links = result.links.collect()
    ids = {p.facets.collect()["facet_id"].item() for p in parents}
    mine = links.filter(pl.col("to_facet_id") == facets["facet_id"].item())
    assert set(mine["from_facet_id"].to_list()) == ids and set(mine["relation"].to_list()) == {CONSENSUS}
    assert link_disagreements(facets, links) == []
    graded = genotype_evidence(result).collect()
    assert graded["grade_basis"].cast(pl.String).unique().to_list() == ["gp"]
    assert facets["ref_seqcol_digest"].item() == parents[0].facets.collect()["ref_seqcol_digest"].item()


def test_the_consensus_is_deterministic_and_order_free(tmp_path: Path, store: ReferenceStore) -> None:
    parents = [
        facet(tmp_path, n, c, "1", [snv("1", 100, g1), snv("1", 300, g2), snv("2", 40, g1)], store)
        for n, c, g1, g2 in (("x", "X", "0/1", "1/1"), ("y", "Y", "1/1", "0/1"), ("z", "Z", "0/1", "0/1"))
    ]
    one = consensus(parents, atlas=TRIANGLE)
    two = consensus(parents, atlas=TRIANGLE)
    flipped = consensus(list(reversed(parents)), atlas=TRIANGLE)
    for other_result in (two, flipped):
        assert one.facets.collect().equals(other_result.facets.collect())
        assert one.calls.collect().equals(other_result.calls.collect())
        assert one.support.collect().equals(other_result.support.collect())
        assert one.links.collect().equals(other_result.links.collect())
    assert one.calls.collect()["chrom"].to_list() == ["1", "1", "2"]


def test_a_person_is_taken_as_its_facets(tmp_path: Path, store: ReferenceStore) -> None:
    parents = [
        facet(tmp_path, n, c, "1", [snv("1", 100, g)], store, subject="p")
        for n, c, g in (("x", "X", "0/1"), ("y", "Y", "0/1"))
    ]
    person = Person.of(*parents, check=None)
    from_person = consensus(person, atlas=TRIANGLE)
    assert from_person.calls.collect().equals(consensus(parents, atlas=TRIANGLE).calls.collect())
    assert from_person.facets.collect()["subject"].item() == "p"


def test_refusals(tmp_path: Path, store: ReferenceStore) -> None:
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1")], store)
    with pytest.raises(ValueError, match="two facets or more"):
        consensus([x], atlas=TRIANGLE)
    with pytest.raises(ValueError, match="more than once"):
        consensus([x, x], atlas=TRIANGLE)
    child = normalize(x, store)
    with pytest.raises(ValueError, match="was made from"):
        consensus([x, child, y], atlas=TRIANGLE)
    elsewhere = facet(tmp_path, "w", "Y", "1", [snv("1", 100, "0/1")], store, ref_seqcol_digest="another", assembly_label="GRCh37")
    with pytest.raises(ReferenceMismatchError):
        consensus([x, elsewhere], atlas=TRIANGLE)
    unread = facet(tmp_path, "u", "Y", "1", [snv("1", 100, "0/1")], store, ref_seqcol_digest=None, assembly_label="GRCh37")
    with pytest.raises(ReferenceMismatchError):
        consensus([x, unread], atlas=TRIANGLE)
    first = consensus([x, y], atlas=TRIANGLE)
    with pytest.raises(ValueError, match="consensus facets"):
        consensus([first, x], atlas=TRIANGLE)
    regions = RegionSet.from_frames("t", "GRCh37", {"easy": pl.DataFrame({"chrom": ["1"], "pos": [1], "end_pos": [10]})})
    with pytest.raises(ValueError, match="regions are on"):
        consensus([x, y], atlas=TRIANGLE, regions=regions)


def test_label_only_assembly_is_warned(tmp_path: Path, store: ReferenceStore) -> None:
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store, ref_seqcol_digest=None)
    y = facet(tmp_path, "y", "Y", "1", [snv("1", 100, "0/1")], store, ref_seqcol_digest=None)
    result = consensus([x, y], atlas=TRIANGLE)
    assert any("not checked by content" in line for line in result.warnings)
    assert result.facets.collect()["ref_seqcol_digest"].item() is None


def test_a_row_over_another_span_says_nothing_about_this_allele(tmp_path: Path, store: ReferenceStore) -> None:
    # A haplotype caller's complex allele over two bases holds no allele of the one-base span.
    ref2 = SEQ["1"][99:101]
    x = facet(tmp_path, "x", "X", "1", [snv("1", 100, "0/1")], store)
    y = facet(tmp_path, "y", "Y", "1", [("1", 100, ref2, other(ref2[0]) + other(ref2[1]), "0/1", "PASS", {})], store)
    result = consensus([x, y], atlas=TRIANGLE)
    rows = rows_of(result)
    assert len(rows) == 2
    assert [(r["s_supported"], r["s_silent"]) for r in rows] == [(1, 1), (1, 1)]
    assert any("other_span 1" in line for line in result.warnings)
