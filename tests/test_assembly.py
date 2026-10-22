"""Assembly label inference and the coordinate-system digest from `##contig` lines (RM4)."""

import json
from itertools import combinations
from pathlib import Path

import pytest
import refget

from just_vcf import load
from just_vcf.assembly import (
    ASSEMBLY_LENGTHS,
    MIN_MATCHING_CONTIGS,
    AssemblyOutcome,
    infer_assembly,
    sorted_name_length_pairs_digest,
)

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
