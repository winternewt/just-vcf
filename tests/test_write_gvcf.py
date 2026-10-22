"""RM6: `write_vcf(..., blocks=True)` writes the `blocks` frame as gVCF reference blocks, with REF read
from the reference, interleaved with the calls; `load` reads the file back to the same frames.

A synthetic FASTA and gVCF under tmp_path. Expectations come from the sequence and the input frames.
"""

import gzip
import random
from pathlib import Path

import polars as pl
import pytest

from just_vcf import BlockAllele, ReferenceStore, UnwritableRowsError, load, resolve_ref, write_vcf
from just_vcf.loading import Loaded

rng = random.Random(5)
SEQ = {name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in (("1", 3000), ("2", 800), ("X", 1200))}


def at(chrom: str, pos: int) -> str:
    return SEQ[chrom][pos - 1]


def other(base: str) -> str:
    return next(b for b in "ACGT" if b != base)


HEADER = [
    "##fileformat=VCFv4.2",
    *[f"##contig=<ID={name},length={len(seq)}>" for name, seq in SEQ.items()],
    '##ALT=<ID=NON_REF,Description="Any other allele">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Block minimum depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
]


def block(chrom: str, pos: int, end: int, sample: str, fmt: str = "GT:GQ:MIN_DP:PL") -> str:
    return f"{chrom}\t{pos}\t.\t{at(chrom, pos)}\t<NON_REF>\t.\t.\tEND={end}\t{fmt}\t{sample}"


def call(chrom: str, pos: int, alt: str, sample: str, fmt: str = "GT:GQ:DP:AD:PL") -> str:
    return f"{chrom}\t{pos}\t.\t{at(chrom, pos)}\t{alt}\t40\tPASS\t.\t{fmt}\t{sample}"


# A DeepVariant/GATK-shaped gVCF: diploid blocks with PL over [REF, <*>], a no-call block, variant
# records carrying <NON_REF>, a block with nothing but GT, a haploid block on X (PL of two), and a block
# whose END runs over 16 kb windows' worth of nothing (index), on a second contig with no calls at all.
RECORDS = [
    block("1", 1, 99, "0/0:60:23:0,60,900"),
    call("1", 100, f"{other(at('1', 100))},<NON_REF>", "0/1:99:30:14,16,0:51,0,36,93,92,86"),
    block("1", 101, 400, "0/0:30:12:0,30,450"),
    block("1", 401, 450, "./.:.:.:.", fmt="GT:GQ:MIN_DP:PL"),
    block("1", 451, 451, "0/0", fmt="GT"),
    call("1", 452, f"{other(at('1', 452))},<NON_REF>", "1/1:20:8:0,8,0:300,20,0,310,25,320"),
    block("1", 453, 3000, "0/0:50:19:0,50,800"),
    block("2", 1, 800, "0/0:10:3:0,10,100"),
    block("X", 1, 600, "0:45:20:0,45"),
    call("X", 601, f"{other(at('X', 601))},<NON_REF>", "1:50:21:0,21,0:500,0,600", fmt="GT:GQ:DP:AD:PL"),
    block("X", 602, 1200, ".:.:.:.", fmt="GT:GQ:MIN_DP:PL"),
]


def write_fasta(path: Path, seqs: dict[str, str] = SEQ) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in seqs.items()))
    return path


@pytest.fixture
def setup(tmp_path: Path) -> tuple[ReferenceStore, Loaded]:
    store = ReferenceStore.open(tmp_path / "store")
    store.add_fasta(write_fasta(tmp_path / "ref.fa"))
    source = tmp_path / "in.g.vcf"
    source.write_text("\n".join([*HEADER, *RECORDS]) + "\n")
    return store, resolve_ref(load(source), store)


CALL_KEY = ["chrom", "pos", "source_ref", "alts", "gt", "gt_alleles", "end_pos", "gq", "dp", "ad", "gl", "qual", "filters"]
BLOCK_COLUMNS = ["chrom", "pos", "end_pos", "kind", "gq", "dp_min", "gl"]


def test_blocks_round_trip_through_load(tmp_path: Path, setup) -> None:
    store, loaded = setup
    original_blocks = loaded.blocks.collect()
    assert original_blocks.height == 8 and set(original_blocks["kind"]) == {"ref", "no_call"}
    written = write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True, reference=store)
    assert written.blocks_written == 8 and written.blocks_not_written == 0 and written.rows_written == 3
    assert not any("block" in w for w in written.warnings)

    back = load(written.path)
    assert back.blocks.collect().select(BLOCK_COLUMNS).equals(original_blocks.select(BLOCK_COLUMNS))
    assert back.calls.collect().select(CALL_KEY).equals(loaded.calls.collect().select(CALL_KEY))
    assert back.facets.collect()["absence_means"].to_list() == ["blocks"]


def test_block_records_take_ref_from_the_reference_and_interleave_with_calls(tmp_path: Path, setup) -> None:
    store, loaded = setup
    written = write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True, reference=store)
    text = gzip.decompress(written.path.read_bytes()).decode()
    body = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
    assert len(body) == 11
    order = {"1": 0, "2": 1, "X": 2}
    assert [(r[0], int(r[1])) for r in body] == sorted(((r[0], int(r[1])) for r in body), key=lambda k: (order[k[0]], k[1]))
    for r in body:
        assert r[3] == at(r[0], int(r[1])), "REF is the reference base at POS"
        if r[4] == "<*>":
            assert r[7].startswith("END=")
    blocks = {(r[0], int(r[1])): r for r in body if r[4] == "<*>"}
    assert len(blocks) == 8
    gts = {key: r[9].split(":")[0] for key, r in blocks.items()}
    # GT asserts reference (`kind = ref`) or nothing (`no_call`); its ploidy is the block's own GL's
    # (two entries: haploid), else that of the contig's calls (X's call is haploid).
    assert gts[("1", 1)] == gts[("1", 451)] == gts[("2", 1)] == "0/0" and gts[("1", 401)] == "./."
    assert gts[("X", 1)] == "0" and gts[("X", 602)] == "."
    assert blocks[("1", 1)][8] == "GT:GQ:MIN_DP:GL" and blocks[("1", 451)][8] == "GT"
    assert blocks[("1", 453)][7] == "END=3000"
    meta = [line for line in text.splitlines() if line.startswith("##")]
    assert '##ALT=<ID=*,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">' in meta
    assert any(line.startswith("##FORMAT=<ID=MIN_DP,Number=1,Type=Integer") for line in meta)
    assert any(line.startswith("##INFO=<ID=END,") for line in meta)
    # The same frames, the same bytes; and the GATK spelling on request.
    again = write_vcf(loaded, tmp_path / "again.g.vcf.gz", blocks=True, reference=store)
    assert again.sha256 == written.sha256
    gatk = write_vcf(loaded, tmp_path / "gatk.g.vcf.gz", blocks=True, reference=store, block_allele=BlockAllele.NON_REF)
    gatk_text = gzip.decompress(gatk.path.read_bytes()).decode()
    assert "\t<NON_REF>\t" in gatk_text and "\t<*>\t" not in gatk_text
    assert load(gatk.path).blocks.collect().select(BLOCK_COLUMNS).equals(loaded.blocks.collect().select(BLOCK_COLUMNS))


def test_without_a_reference_blocks_are_refused_never_invented(tmp_path: Path, setup) -> None:
    store, loaded = setup
    unresolved = load(tmp_path / "in.g.vcf")
    with pytest.raises(ValueError, match="resolve_ref"):
        write_vcf(unresolved, tmp_path / "out.g.vcf.gz", blocks=True, reference=store)
    with pytest.raises(ValueError, match="reference"):
        write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True)
    assert not (tmp_path / "out.g.vcf.gz").exists()


def test_blocks_on_a_contig_the_reference_cannot_vouch_for_are_counted(tmp_path: Path) -> None:
    store = ReferenceStore.open(tmp_path / "store")
    store.add_fasta(write_fasta(tmp_path / "ref.fa", {name: seq for name, seq in SEQ.items() if name != "2"}))
    source = tmp_path / "in.g.vcf"
    source.write_text("\n".join([*HEADER, *RECORDS]) + "\n")
    loaded = resolve_ref(load(source), store)
    with pytest.raises(UnwritableRowsError) as refused:
        write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True, reference=store)
    assert refused.value.counts == {"block_ref_unread": 1}
    written = write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True, reference=store, resolved_only=True)
    assert written.rows_skipped == {"block_ref_unread": 1} and written.blocks_written == 7 and written.blocks_not_written == 1
    assert load(written.path).blocks.collect().filter(pl.col("chrom") == "2").is_empty()


def test_htslib_reads_the_gvcf_and_finds_a_block_by_its_end(tmp_path: Path, setup) -> None:
    """Opt-in: htslib (pysam) parses every record and the index finds a block from inside its span."""
    pysam = pytest.importorskip("pysam")
    store, loaded = setup
    written = write_vcf(loaded, tmp_path / "out.g.vcf.gz", blocks=True, reference=store)
    with pysam.VariantFile(str(written.path)) as reader:
        assert sum(1 for _ in reader) == 11
        inside = [(r.chrom, r.pos, r.stop) for r in reader.fetch("1", 2000, 2001)]
    assert inside == [("1", 453, 3000)]
