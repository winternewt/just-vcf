"""Assembly label inference and the coordinate-system digest from `##contig` lines (RM4)."""

import json
from itertools import combinations
from pathlib import Path

import polars as pl
import pytest
import refget

from just_vcf import load, sniff
from just_vcf.assembly import (
    ASSEMBLY_LENGTHS,
    MIN_MARGIN,
    MIN_MATCHING_CONTIGS,
    MIN_SCORED_ROWS,
    MIN_SHARE,
    AssemblyOutcome,
    RowsOutcome,
    infer_assembly,
    infer_assembly_from_rows,
    sample_snv_rows,
    sorted_name_length_pairs_digest,
)
from just_vcf.schema import CHM13_ACCESSIONS, canonical_chrom

# Domain constants: the first three chromosomes and the mitochondrion of each assembly, as UCSC's
# chrom.sizes give them (hg19, hg38, hs1), typed here independently of `ASSEMBLY_LENGTHS`.
GRCH37 = {"1": 249250621, "2": 243199373, "3": 198022430}
GRCH38 = {"1": 248956422, "2": 242193529, "3": 198295559}
T2T = {"1": 248387328, "2": 242696752, "3": 201105948}

COMPLIANCE = Path(refget.__file__).parent / "compliance_data" / "test_fasta_digests.json"


def contigs(lengths: dict[str, int], prefix: str = "") -> list[dict]:
    return [{"name": f"{prefix}{name}", "length": length, "md5": None} for name, length in lengths.items()]


@pytest.mark.parametrize(
    ("lengths", "expected"), [(GRCH37, "GRCh37"), (GRCH38, "GRCh38"), (T2T, "T2T-CHM13")]
)
@pytest.mark.parametrize("prefix", ["", "chr"])
def test_three_primary_lengths_name_the_assembly(lengths: dict[str, int], expected: str, prefix: str) -> None:
    inferred = infer_assembly(contigs(lengths, prefix))
    assert inferred and inferred.label == expected
    assert inferred.outcome is AssemblyOutcome.INFERRED
    assert inferred.matched == {expected: len(lengths)}


def test_alt_contigs_and_decoys_neither_help_nor_hurt() -> None:
    header = contigs(GRCH38, "chr") + [
        {"name": "chrUn_KI270302v1", "length": 2274, "md5": None},
        {"name": "chrEBV", "length": 171823, "md5": None},
        {"name": "hs37d5", "length": 35477943, "md5": None},
    ]
    assert infer_assembly(header).label == "GRCh38"


def test_the_mitochondrion_contradicts_but_never_supports() -> None:
    only_mt = infer_assembly([{"name": "chrM", "length": 16569, "md5": None}])
    assert not only_mt and only_mt.outcome is AssemblyOutcome.TOO_FEW
    # hg19 kept the Yoruba chrM (16,571): fine beside GRCh37 lengths, a contradiction beside GRCh38.
    hg19_mt = [{"name": "chrM", "length": 16571, "md5": None}]
    assert infer_assembly(contigs(GRCH37, "chr") + hg19_mt).label == "GRCh37"
    against = infer_assembly(contigs(GRCH38, "chr") + hg19_mt)
    assert not against and against.outcome is AssemblyOutcome.CONFLICTING


def test_too_few_contigs_stay_unknown() -> None:
    two = dict(list(GRCH38.items())[: MIN_MATCHING_CONTIGS - 1])
    inferred = infer_assembly(contigs(two))
    assert not inferred and inferred.label is None
    assert inferred.outcome is AssemblyOutcome.TOO_FEW


def test_mixed_assemblies_are_a_conflict_not_a_vote() -> None:
    # Two GRCh38 lengths outvote one GRCh37 length in just-prs; here the header is not trusted at all.
    mixed = contigs({"1": GRCH37["1"]}) + contigs({k: GRCH38[k] for k in ("2", "3")}) + contigs({"4": 190214555})
    inferred = infer_assembly(mixed)
    assert not inferred and inferred.outcome is AssemblyOutcome.CONFLICTING
    assert inferred.matched == {"GRCh37": 1, "GRCh38": 3}


def test_a_length_no_assembly_has_is_a_conflict() -> None:
    inferred = infer_assembly(contigs(GRCH38) + [{"name": "4", "length": 191154277, "md5": None}])
    assert not inferred and inferred.outcome is AssemblyOutcome.CONFLICTING
    assert inferred.unmatched == ("4",)


def test_no_lengths_is_its_own_answer() -> None:
    assert infer_assembly([]).outcome is AssemblyOutcome.NO_CONTIG_LENGTHS
    no_length = [{"name": name, "length": None, "md5": None} for name in GRCH38]
    assert infer_assembly(no_length).outcome is AssemblyOutcome.NO_CONTIG_LENGTHS


def test_assembly_tables_share_no_length_but_the_mitochondrion() -> None:
    for (a, table_a), (b, table_b) in combinations(ASSEMBLY_LENGTHS.items(), 2):
        shared = {chrom for chrom in table_a if table_a[chrom] & table_b[chrom]}
        assert shared == {"MT"}, (a, b, shared)
    for lengths, name in ((GRCH37, "GRCh37"), (GRCH38, "GRCh38"), (T2T, "T2T-CHM13")):
        assert all(ASSEMBLY_LENGTHS[name][chrom] == {length} for chrom, length in lengths.items())


@pytest.mark.parametrize("case", sorted(json.loads(COMPLIANCE.read_text())))
def test_sorted_name_length_pairs_digest_matches_refgets_compliance_vectors(case: str) -> None:
    vector = json.loads(COMPLIANCE.read_text())[case]
    level2 = vector["level2"]
    header = [{"name": n, "length": length, "md5": None} for n, length in zip(level2["names"], level2["lengths"], strict=True)]
    assert sorted_name_length_pairs_digest(header) == vector["sorted_name_length_pairs_digest"]


def test_sorted_name_length_pairs_digest_is_withheld_without_every_length() -> None:
    assert sorted_name_length_pairs_digest([]) is None
    partial = contigs(GRCH38) + [{"name": "4", "length": None, "md5": None}]
    assert sorted_name_length_pairs_digest(partial) is None
    # Names are verbatim: `chr1` and `1` are two coordinate systems to seqcol.
    assert sorted_name_length_pairs_digest(contigs(GRCH38)) != sorted_name_length_pairs_digest(contigs(GRCH38, "chr"))


def vcf_with_contigs(header_contigs: list[dict], extra_meta: tuple[str, ...] = ()) -> str:
    meta = ["##fileformat=VCFv4.2", *extra_meta]
    meta += [f"##contig=<ID={c['name']}" + (f",length={c['length']}" if c["length"] else "") + ">" for c in header_contigs]
    meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    first = header_contigs[0]["name"] if header_contigs else "1"
    rows = [f"{first}\t{100 * i}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1" for i in range(1, 4)]
    return "\n".join([*meta, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1", *rows]) + "\n"


def facet(tmp_path: Path, text: str) -> tuple[dict, tuple[str, ...]]:
    path = tmp_path / "calls.vcf"
    path.write_text(text)
    loaded = load(path)
    return loaded.facets.collect().row(0, named=True), loaded.warnings


@pytest.mark.parametrize(("lengths", "expected"), [(GRCH37, "GRCh37"), (GRCH38, "GRCh38"), (T2T, "T2T-CHM13")])
def test_load_infers_the_label_from_contig_lengths(tmp_path: Path, lengths: dict[str, int], expected: str) -> None:
    header = contigs(lengths, "chr")
    row, warnings = facet(tmp_path, vcf_with_contigs(header))
    assert (row["assembly_label"], row["assembly_label_basis"]) == (expected, "inferred")
    assert row["sorted_name_length_pairs_digest"] == sorted_name_length_pairs_digest(header)
    assert row["seqcol_digest"] is None
    assert not [w for w in warnings if "assembl" in w or "states" in w]


def test_a_stated_label_is_kept_and_the_disagreement_counted(tmp_path: Path) -> None:
    row, warnings = facet(tmp_path, vcf_with_contigs(contigs(GRCH38, "chr"), ("##reference=file:///refs/GRCh37.fa",)))
    assert (row["assembly_label"], row["assembly_label_basis"]) == ("GRCh37", "stated")
    assert [w for w in warnings if "states" in w] == [
        "the header states GRCh37, but 3 primary chromosome length(s) match GRCh38; the stated label is kept"
    ]


def test_a_stated_label_that_the_lengths_agree_with_warns_nothing(tmp_path: Path) -> None:
    row, warnings = facet(tmp_path, vcf_with_contigs(contigs(GRCH38, "chr"), ("##reference=file:///refs/hg38.fa",)))
    assert (row["assembly_label"], row["assembly_label_basis"]) == ("GRCh38", "stated")
    assert not [w for w in warnings if "assembl" in w or "states" in w]


def test_load_leaves_mixed_missing_or_thin_headers_unknown(tmp_path: Path) -> None:
    for name in ("mixed", "bare", "thin"):
        (tmp_path / name).mkdir()
    mixed = contigs({"1": GRCH37["1"]}, "chr") + contigs({k: GRCH38[k] for k in ("2", "3")}, "chr")
    row, warnings = facet(tmp_path / "mixed", vcf_with_contigs(mixed))
    assert (row["assembly_label"], row["assembly_label_basis"]) == (None, None)
    assert [w for w in warnings if "more than one assembly" in w] == [
        "##contig lengths match more than one assembly (contig lengths disagree: 1 match GRCh37; 2 match GRCh38)"
    ]
    bare = [{"name": name, "length": None, "md5": None} for name in GRCH38]
    row, _ = facet(tmp_path / "bare", vcf_with_contigs(bare))
    assert (row["assembly_label"], row["sorted_name_length_pairs_digest"]) == (None, None)
    row, _ = facet(tmp_path / "thin", vcf_with_contigs(contigs({"1": GRCH38["1"]})))
    assert row["assembly_label"] is None and row["sorted_name_length_pairs_digest"] is not None
    assert row["sorted_name_length_pairs_digest"] == sorted_name_length_pairs_digest(contigs({"1": GRCH38["1"]}))


# T2T-CHM13 v2.0 (RM286): its names as NCBI's assembly report writes them, and its header shapes.

CHM13_REPORT = Path(__file__).resolve().parent.parent / "assets" / "ncbi_assembly_reports" / "GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt"


def chm13_report() -> list[dict[str, str]]:
    """The report's sequence rows, keyed by its own column names."""
    lines = CHM13_REPORT.read_text().splitlines()
    columns = next(line for line in lines if line.startswith("# Sequence-Name")).removeprefix("# ").split("\t")
    return [dict(zip(columns, line.split("\t"), strict=True)) for line in lines if line and not line.startswith("#")]


def test_chm13_accessions_are_the_assembly_reports() -> None:
    rows = chm13_report()
    expected = {row["GenBank-Accn"]: row["Assigned-Molecule"] for row in rows}
    expected |= {row["RefSeq-Accn"]: row["Assigned-Molecule"] for row in rows if row["RefSeq-Accn"] != "na"}
    assert CHM13_ACCESSIONS == expected
    assert len(rows) == 25 and sum(row["RefSeq-Accn"] == "na" for row in rows) == 1
    for row in rows:
        chrom = canonical_chrom(row["UCSC-style-name"])
        assert canonical_chrom(row["GenBank-Accn"]) == chrom
        if chrom != "MT":
            assert ASSEMBLY_LENGTHS["T2T-CHM13"][chrom] == {int(row["Sequence-Length"])}
    assert ASSEMBLY_LENGTHS["T2T-CHM13"]["MT"] == {16569}


@pytest.mark.parametrize("column", ["GenBank-Accn", "RefSeq-Accn"])
def test_an_accession_named_header_names_t2t(tmp_path: Path, column: str) -> None:
    rows = [row for row in chm13_report() if row[column] != "na"]
    header = [{"name": row[column], "length": int(row["Sequence-Length"]), "md5": None} for row in rows]
    inferred = infer_assembly(header)
    assert inferred.label == "T2T-CHM13" and inferred.matched == {"T2T-CHM13": 24}
    row, _ = facet(tmp_path, vcf_with_contigs(header))
    assert (row["assembly_label"], row["assembly_label_basis"]) == ("T2T-CHM13", "inferred")
    first = header[0]["name"]
    assert {entry["source"]: entry["chrom"] for entry in row["contig_map"]} == {first: "1"}


def test_a_t2t_shaped_vcf_header_loads_as_t2t(tmp_path: Path) -> None:
    """The analysis set's names (`chr1`…`chrM`) and lengths, and a `##reference` naming its FASTA."""
    rows = chm13_report()
    header = [{"name": row["UCSC-style-name"], "length": int(row["Sequence-Length"]), "md5": None} for row in rows]
    row, warnings = facet(tmp_path, vcf_with_contigs(header, ("##reference=file:///refs/chm13v2.0.fa",)))
    assert (row["assembly_label"], row["assembly_label_basis"]) == ("T2T-CHM13", "stated")
    assert row["sorted_name_length_pairs_digest"] == sorted_name_length_pairs_digest(header)
    assert {entry["name"]: entry["length"] for entry in row["header_contigs"]} == {c["name"]: c["length"] for c in header}
    assert not [w for w in warnings if "assembl" in w or "states" in w]


@pytest.mark.parametrize(
    ("line", "build"),
    [
        ("##reference=file:///refs/chm13v2.0.fa", "T2T-CHM13"),
        ("##reference=https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz", "T2T-CHM13"),
        ("##reference=GCA_009914755.4_T2T-CHM13v2.0_genomic.fna", "T2T-CHM13"),
        ("##contig=<ID=chr1,length=248387328,assembly=chm13v2.0>", "T2T-CHM13"),
        ("##contig=<ID=chr1,length=248387328,assembly=hs1>", "T2T-CHM13"),
        ("##reference=file:///refs/hs37d5.fa", "GRCh37"),
        ("##reference=file:///refs/hs38DH.fa", "GRCh38"),
        ("##reference=file:///refs/hs10.fa", None),
        ("##reference=file:///refs/GRCh38_to_chm13.fa", None),
    ],
)
def test_vcf_build_hints_name_t2t_and_only_t2t(tmp_path: Path, line: str, build: str | None) -> None:
    lines = ["##fileformat=VCFv4.2", line, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", "chr1\t100\t.\tA\tG\t.\tPASS\t."]
    path = tmp_path / "hint.vcf"
    path.write_text("\n".join(lines) + "\n")
    assert sniff(path).build == build


# F9: an array's build from its own rows. The scorer is a pure function over sampled rows and bases.


def calls_frame(rows: list[tuple[str, int, list[str | None] | None, str | None]]) -> pl.LazyFrame:
    """`(chrom, pos, gt_alleles, source_ref)` rows, the only columns the sampler reads."""
    return pl.LazyFrame(
        rows, schema={"chrom": pl.String, "pos": pl.UInt32, "gt_alleles": pl.List(pl.String), "source_ref": pl.String}, orient="row"
    )


def test_the_sample_keeps_scorable_array_snvs_only() -> None:
    rows = [
        ("1", 10, ["G", "A"], None),  # het, kept as the sorted distinct alleles
        ("1", 11, ["C", None], None),  # half-call: the called base is enough
        ("1", 12, [None, None], None),  # whole no-call
        ("1", 13, ["I", "D"], None),  # indel codes
        ("1", 14, ["A", "A"], "A"),  # a VCF row: it has a source REF
        ("MT", 15, ["A"], None),  # MT is never scored
        ("XY", 16, ["A", "A"], None),  # the PAR code, not a primary chromosome
        ("X", 17, ["T"], None),
        ("1", 18, None, None),  # genotype that could not be split
    ]
    sample = sample_snv_rows(calls_frame(rows))
    assert sample.rows() == [("1", 10, ["A", "G"]), ("1", 11, ["C"]), ("X", 17, ["T"])]


def test_the_sample_is_spread_bounded_and_order_independent() -> None:
    rows = [(chrom, pos, ["A", "C"], None) for chrom in ("1", "2") for pos in range(1, 1001)] + [("3", 5, ["G"], None)]
    shuffled = calls_frame(rows).collect().sample(fraction=1.0, shuffle=True, seed=5).lazy()
    sample = sample_snv_rows(shuffled, per_contig=300)
    assert sample.equals(sample_snv_rows(calls_frame(rows), per_contig=300))
    per_chrom = dict(sample.group_by("chrom").len().iter_rows())
    assert per_chrom == {"1": 250, "2": 250, "3": 1}  # every 4th of 1000 rows, at most 300
    assert sample.filter(pl.col("chrom") == "1")["pos"].to_list() == list(range(1, 1001, 4))


def scored(alleles: list[list[str]], bases: dict[str, list[str | None]], labels: dict[str, str | None]):
    sample = pl.DataFrame({"chrom": ["1"] * len(alleles), "pos": list(range(1, len(alleles) + 1)), "alleles": alleles})
    return infer_assembly_from_rows(sample, {d: pl.Series(b, dtype=pl.String) for d, b in bases.items()}, labels)


def test_rows_name_the_build_that_clears_floor_and_margin() -> None:
    n = MIN_SCORED_ROWS
    alleles = [["A", "G"]] * n
    right = ["A"] * int(n * 0.85) + ["C"] * (n - int(n * 0.85))
    wrong = ["A"] * int(n * 0.40) + ["T"] * (n - int(n * 0.40))
    inferred = scored(alleles, {"d37": right, "d38": wrong}, {"d37": "GRCh37", "d38": "GRCh38"})
    assert inferred and (inferred.label, inferred.digest) == ("GRCh37", "d37")
    assert inferred.outcome is RowsOutcome.INFERRED
    assert inferred.shares == {"d37": pytest.approx(0.85), "d38": pytest.approx(0.40)}
    assert inferred.chance == pytest.approx(0.5)  # two distinct alleles of four bases


def test_a_row_past_the_contig_end_counts_against_and_a_missing_contig_is_not_scored() -> None:
    n = MIN_SCORED_ROWS
    alleles = [["A"]] * (2 * n)
    past_end = ["A"] * n + [""] * n  # half the rows past the end: 50%, below the floor
    missing = ["A"] * n + [None] * n  # half the rows on a contig it lacks: 100% of n
    inferred = scored(alleles, {"short": past_end, "partial": missing}, {"short": "GRCh37", "partial": "GRCh38"})
    assert inferred.scored == {"short": 2 * n, "partial": n}
    assert inferred.shares == {"short": pytest.approx(0.5), "partial": pytest.approx(1.0)}
    assert inferred.label == "GRCh38"


@pytest.mark.parametrize(
    ("right", "rival", "rival_label", "outcome"),
    [
        (MIN_SHARE - 0.05, 0.30, "GRCh38", RowsOutcome.NO_FIT),
        (0.90, 0.90 - MIN_MARGIN + 0.05, "GRCh38", RowsOutcome.NO_MARGIN),
        (0.90, 0.30, None, RowsOutcome.INFERRED),  # an unlabelled rival far below does not block
        (0.30, 0.90, None, RowsOutcome.UNLABELLED),
        (0.90, 0.90 - MIN_MARGIN + 0.05, None, RowsOutcome.NO_MARGIN),  # an unlabelled rival may be another build
    ],
)
def test_the_rows_withhold_a_label_they_do_not_settle(right: float, rival: float, rival_label: str | None, outcome: RowsOutcome) -> None:
    n = 2 * MIN_SCORED_ROWS
    alleles = [["A"]] * n
    bases = {"a": ["A"] * round(n * right) + ["C"] * (n - round(n * right)), "b": ["A"] * round(n * rival) + ["C"] * (n - round(n * rival))}
    inferred = scored(alleles, bases, {"a": "GRCh37", "b": rival_label})
    assert inferred.outcome is outcome
    assert bool(inferred) is (outcome is RowsOutcome.INFERRED)
    assert (inferred.label is None) is (outcome is not RowsOutcome.INFERRED)


def test_chance_is_the_rival_when_one_build_is_scored() -> None:
    n = MIN_SCORED_ROWS
    four = [["A", "C", "G", "T"]] * n  # four alleles called: any base is in them, so chance is 100%
    assert scored(four, {"a": ["A"] * n}, {"a": "GRCh37"}).outcome is RowsOutcome.NO_MARGIN
    assert scored([["A"]] * n, {"a": ["A"] * n}, {"a": "GRCh37"}).label == "GRCh37"


def test_too_few_scored_rows_are_their_own_answer() -> None:
    n = MIN_SCORED_ROWS - 1
    inferred = scored([["A"]] * n, {"a": ["A"] * n}, {"a": "GRCh37"})
    assert inferred.outcome is RowsOutcome.TOO_FEW_ROWS and not inferred
