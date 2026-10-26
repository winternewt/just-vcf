"""RM276: one facet's calls exported as VCF Zarr 0.4, read back by polars-bio (an independent reader).

The expected values are computed from the facet's own `calls` frame at runtime; polars-bio's
`scan_vcf_zarr` is the oracle for the bytes. It returns integer sentinels as they are stored (`-1`
missing, `-2` padding), so the comparison translates them.
"""

import gzip
import hashlib
import json
import struct
from pathlib import Path

import polars as pl
import polars_bio as pb
import pytest

from just_vcf import ContigNames, UnsortedRowsError, UnwritableRowsError, load, write_vcf
from just_vcf.loading import Loaded
from just_vcf.vcf_zarr import write_vcf_zarr
from test_load_cgi import LINKED, master_var_text
from test_load_cgi import write as write_cgi
from test_write_breakends import ASSEMBLY, ASSEMBLY_RECORDS, SOURCES, mates, others, record, source_text, written_name
from test_write_vcf import HQ_RECORDS, RICH, RICH_FORMAT, hq_loaded, loaded_from, psl_loaded

GVCF = [
    "1\t100\t.\tG\t<NON_REF>\t.\t.\tEND=113\tGT:GQ:MIN_DP:PL\t0/0:60:23:0,60,900",
    "1\t114\t.\tT\tTC,<NON_REF>\t213\t.\t.\tGT:GQ:MIN_DP:PL\t0/1:99:.:51,0,36,93,92,86",
    "1\t116\t.\tC\t<*>\t.\t.\tEND=120\tGT:GQ:MIN_DP:PL\t./.:.:.:.",
]


def read_back(path: Path) -> pl.DataFrame:
    return pb.scan_vcf_zarr(str(path), use_zero_based=False).collect()


def expected_sites(calls: pl.DataFrame) -> list[tuple]:
    """What the store should say per row, from the frame: the reader's view of VCF's columns."""
    rows = []
    for row in calls.iter_rows(named=True):
        ref = row["ref"] if row["ref"] is not None else row["source_ref"]
        rows.append(
            (
                row["chrom"],
                row["pos"],
                row["end_pos"],
                ";".join(row["ids"]) if row["ids"] else None,
                ref,
                ",".join(row["alts"]) if row["alts"] else ".",
                row["qual"],
            )
        )
    return rows


def expected_calls(calls: pl.DataFrame) -> list[tuple]:
    """A row without a genotype reads as missing at the store's ploidy (the reader cannot tell it
    from `./.`)."""
    ploidy = calls["gt"].list.len().max()
    rows = []
    for row in calls.iter_rows(named=True):
        gt = tuple(row["gt"]) if row["gt"] is not None else (None,) * ploidy
        phased = row["phased"] or []
        rows.append((gt, len(row["gt"] or ()) > 1 and all(bool(p) for p in phased[1:]), row["gq"], row["dp"]))
    return rows


def got_calls(back: pl.DataFrame) -> list[tuple]:
    rows = []
    for row in back.select(pl.col("genotypes").struct.unnest()).iter_rows(named=True):
        gt = tuple(None if allele == -1 else allele for allele in row["GT"][0] if allele != -2)
        gq, dp = row["GQ"][0], row["DP"][0]
        rows.append((gt, row["GT_phased"][0], None if gq == -1 else gq, None if dp == -1 else dp))
    return rows


def tree(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob("*")) if p.is_file()}


def test_polars_bio_reads_back_every_site_and_genotype(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    calls = loaded.calls.collect()
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    assert written.rows_written == calls.height and written.rows_skipped == {} and written.warnings == ()
    back = read_back(written.path)
    sites = list(back.select("chrom", "start", "end", "id", "ref", "alt", "qual").iter_rows())
    assert sites == expected_sites(calls)
    assert got_calls(back) == expected_calls(calls)
    unnested = back.select(pl.col("genotypes").struct.unnest())
    want_ad = [tuple(row) if row else None for row in calls["ad"].to_list()]
    got_ad = [tuple(v for v in row[0] if v != -2) if any(v != -1 for v in row[0]) else None for row in unnested["AD"].to_list()]
    assert got_ad == want_ad
    want_ps = [int(sets[0]) if sets and sets[0] is not None else None for sets in calls["phase_sets"].to_list()]
    assert [None if row[0] == -1 else row[0] for row in unnested["PS"].to_list()] == want_ps
    assert {"call_GL", "call_GP", "call_DS", "call_FT"} <= set(written.arrays)


def test_the_store_is_spec_shaped(tmp_path: Path) -> None:
    """Checked on the bytes, not through a reader: the version, the dimension names, and the float
    sentinels the spec gives for a missing value (`0x7F800001`) and for padding (`0x7F800002`)."""
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    root = write_vcf_zarr(loaded, tmp_path / "out.vcz", variants_chunk_size=8).path
    attrs = json.loads((root / ".zattrs").read_text())
    assert attrs["vcf_zarr_version"] == "0.4"
    facet = loaded.facets.collect().row(0, named=True)
    assert attrs["just_vcf"]["facet_id"] == facet["facet_id"] and "source_path" not in json.dumps(attrs)
    assert json.loads((root / "call_genotype" / ".zattrs").read_text())["_ARRAY_DIMENSIONS"] == ["variants", "samples", "ploidy"]
    qual = loaded.calls.collect()["qual"].to_list()
    bits = struct.unpack("<8I", gzip.decompress((root / "variant_quality" / "1").read_bytes()))
    tail = qual[8:]
    assert bits[: len(tail)] == tuple(0x7F800001 if q is None else struct.unpack("<I", struct.pack("<f", q))[0] for q in tail)
    assert set(bits[len(tail) :]) == {0x7F800002}, "the last chunk is padded past the variants"
    assert qual.count(None) > 0
    # FILTER: `PASS` first, then the others sorted; a row with no FILTER (`.`) sets none of them,
    # which a reader cannot tell from PASS unless it looks at the bytes.
    filters = loaded.calls.collect()["filters"].to_list()
    ids = vlen(gzip.decompress((root / "filter_id" / "0").read_bytes()))
    assert ids == ["PASS", *sorted({name for value in filters if value for name in value} - {"PASS"})]
    flags = b"".join(gzip.decompress((root / "variant_filter" / f"{i}.0").read_bytes()) for i in range(2))
    got = [[name for j, name in enumerate(ids) if flags[row * len(ids) + j]] for row in range(len(filters))]
    assert got == [sorted(value, key=ids.index) if value else [] for value in filters]
    assert [] in got and ["PASS"] in got


def vlen(raw: bytes) -> list[str]:
    """numcodecs' `vlen-utf8`: a count, then each item's length and bytes."""
    (count,), offset, items = struct.unpack_from("<I", raw), 4, []
    for _ in range(count):
        (size,) = struct.unpack_from("<I", raw, offset)
        items.append(raw[offset + 4 : offset + 4 + size].decode())
        offset += 4 + size
    return items


def test_chunking_does_not_change_what_is_read(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    one = read_back(write_vcf_zarr(loaded, tmp_path / "one.vcz").path)
    many = read_back(write_vcf_zarr(loaded, tmp_path / "many.vcz", variants_chunk_size=3).path)
    assert len(list((tmp_path / "many.vcz" / "variant_position").glob("[0-9]*"))) == 4
    assert many.equals(one)


def test_two_writes_of_one_frame_are_byte_identical(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    first = tree(write_vcf_zarr(loaded, tmp_path / "a.vcz").path)
    second = tree(write_vcf_zarr(load(tmp_path / "in.vcf"), tmp_path / "b.vcz").path)
    assert first == second and len(first) > 20


def test_blocks_are_not_representable_and_are_counted(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, GVCF, fmt=("GT", "GQ", "MIN_DP", "PL"), info=("END",))
    blocks = loaded.blocks.collect().height
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    assert written.blocks_not_written == blocks == 2
    assert any(w.startswith(f"{blocks} reference/no-call block(s) not written: VCF Zarr has no blocks") for w in written.warnings)
    assert read_back(written.path).height == written.rows_written == loaded.calls.collect().height == 1
    assert json.loads((written.path / ".zattrs").read_text())["just_vcf"]["blocks_not_written"] == blocks


def test_rows_vcf_cannot_hold_are_refused_as_write_vcf_refuses_them(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1", "1\t200\t.\tCT\t-\t.\t.\t.\tGT\t0/1", "1\t300\t.\tG\tA\t.\t.\t.\tGT\t1/1"]
    loaded = loaded_from(tmp_path, records)
    with pytest.raises(UnwritableRowsError) as as_vcf:
        write_vcf(loaded, tmp_path / "out.vcf.gz")
    with pytest.raises(UnwritableRowsError) as as_zarr:
        write_vcf_zarr(loaded, tmp_path / "out.vcz")
    assert as_zarr.value.counts == as_vcf.value.counts == {"alt_not_an_allele": 1}
    assert not (tmp_path / "out.vcz").exists() and not (tmp_path / "out.vcz.partial").exists()
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz", resolved_only=True)
    assert written.rows_skipped == {"alt_not_an_allele": 1}
    assert read_back(written.path)["start"].to_list() == [100, 300]


def test_split_rows_carry_the_split_flag(tmp_path: Path) -> None:
    records = ["1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/0", "1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1", "1\t300\t.\tC\tG,T\t50\tPASS\t.\tGT\t1/2"]
    loaded = loaded_from(tmp_path, records)
    flagged = loaded.calls.collect()["split_signal"].is_not_null().to_list()
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    back = pb.scan_vcf_zarr(str(written.path), use_zero_based=False, info_fields=["JUST_VCF_SPLIT"]).collect()
    assert back["JUST_VCF_SPLIT"].to_list() == flagged and any(flagged) and not all(flagged)


def test_an_existing_path_is_never_written_over_and_unsorted_rows_leave_nothing(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    (tmp_path / "taken.vcz").mkdir()
    with pytest.raises(FileExistsError):
        write_vcf_zarr(loaded, tmp_path / "taken.vcz")
    shuffled = Loaded(facets=loaded.facets, calls=loaded.calls.collect().reverse().lazy(), blocks=loaded.blocks)
    with pytest.raises(UnsortedRowsError):
        write_vcf_zarr(shuffled, tmp_path / "out.vcz")
    assert not (tmp_path / "out.vcz").exists() and not (tmp_path / "out.vcz.partial").exists()


def test_a_cgi_hap_link_is_written_unphased(tmp_path: Path) -> None:
    """VCF Zarr has no field for a hapLink: as in `write_vcf`, a row it phases is written unphased,
    since `GT_phased` with no PS would put it in the chromosome's implicit set."""
    loaded = load(write_cgi(tmp_path, master_var_text(LINKED), "masterVarBeta-GS000000004-ASM.tsv.bz2"))
    linked = loaded.calls.filter(pl.col("phase_sets").is_not_null()).select(pl.len()).collect().item()
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    back = read_back(written.path)
    phased = back.select(pl.col("genotypes").struct.unnest()).get_column("GT_phased").list.first()
    assert linked > 0 and back.height == loaded.calls.select(pl.len()).collect().item()
    assert not phased.any()
    assert f"{linked} row(s) phased by a CGI hapLink, which VCF Zarr has no field for; written unphased" in written.warnings


def test_a_psl_phase_set_is_written_unphased_under_the_header_write_vcf_writes(tmp_path: Path) -> None:
    """VCF Zarr has one phase flag per call and no PSL: as for a hapLink, a row with a slot in a PSL
    set is written unphased, since `true` would put it in the chromosome's implicit set. A PS row
    beside it keeps its phase and its set. The store's `fileformat` is the one `write_vcf` writes."""
    loaded = psl_loaded(tmp_path)
    calls = loaded.calls.collect()
    in_psl = calls["phase_set_kinds"].list.contains("psl").fill_null(False)
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    back = read_back(written.path)
    phased = back.select(pl.col("genotypes").struct.unnest()).get_column("GT_phased").list.first()
    expected = [
        not psl and len(row["gt"]) > 1 and all(bool(p) for p in row["phased"][1:])
        for psl, row in zip(in_psl.to_list(), calls.iter_rows(named=True))
    ]
    assert in_psl.sum() > 0 and any(expected), "the fixture has PSL rows and a phased PS row"
    assert phased.to_list() == expected
    assert f"{in_psl.sum()} row(s) with a slot in a VCF 4.5 PSL phase set, which VCF Zarr has no field for; written unphased" in written.warnings
    vcf_header = gzip.decompress(write_vcf(loaded, tmp_path / "out.vcf.gz").path.read_bytes()).decode().splitlines()[0]
    attrs = json.loads((written.path / ".zattrs").read_text())
    assert ["fileformat", vcf_header.removeprefix("##fileformat=")] in attrs["vcf_meta_information"]
    assert vcf_header == "##fileformat=VCFv4.5"


def test_a_facet_without_psl_says_vcf_43_as_write_vcf_does(tmp_path: Path) -> None:
    loaded = loaded_from(tmp_path, RICH, fmt=RICH_FORMAT)
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    attrs = json.loads((written.path / ".zattrs").read_text())
    vcf_header = gzip.decompress(write_vcf(loaded, tmp_path / "out.vcf.gz", resolved_only=True).path.read_bytes()).decode().splitlines()[0]
    assert attrs["vcf_meta_information"] == [["fileformat", vcf_header.removeprefix("##fileformat=")]]


# A triploid row beside HQ_RECORDS' diploid and haploid ones: the ploidy dimension holds it, where
# write_vcf's Number=2 HQ cannot.
HQ_ZARR_RECORDS = [*HQ_RECORDS, "2\t100\t.\tA\tG\t50\tPASS\t.\tGT:HQ\t1/1/0:5,.,7"]


def expected_hq(calls: pl.DataFrame, ploidy: int) -> list[list[int]]:
    """call_HQ per row as stored: each GT slot's score, `-1` where it has none (or the row's `hq` is
    null), `-2` past the row's ploidy; a row without a genotype missing throughout."""
    out = []
    for gt, hq in zip(calls["gt"].to_list(), calls["hq"].to_list()):
        if gt is None:
            out.append([-1] * ploidy)
            continue
        values = hq if hq is not None else [None] * len(gt)
        out.append([-1 if v is None else v for v in values] + [-2] * (ploidy - len(gt)))
    return out


def test_haplotype_qualities_are_written_per_ploidy_slot(tmp_path: Path) -> None:
    loaded = hq_loaded(tmp_path, HQ_ZARR_RECORDS, number=".")
    calls = loaded.calls.collect()
    ploidy = calls["gt"].list.len().max()
    assert ploidy == 3 and calls["hq"].null_count() >= 1 and any(None in hq for hq in calls["hq"].drop_nulls().to_list())
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz", variants_chunk_size=4)
    assert "call_HQ" in written.arrays and written.warnings == ()
    root = written.path
    assert json.loads((root / "call_HQ" / ".zattrs").read_text())["_ARRAY_DIMENSIONS"] == ["variants", "samples", "ploidy"]
    meta = json.loads((root / "call_HQ" / ".zarray").read_text())
    assert meta["shape"] == [calls.height, 1, ploidy] and meta["dtype"] == "<i4" and meta["fill_value"] == -2
    # On the bytes: the spec's missing (-1) and fill (-2), the last chunk padded past the variants.
    want = expected_hq(calls, ploidy)
    chunks = sorted(root.joinpath("call_HQ").glob("[0-9]*"), key=lambda p: int(p.name.split(".")[0]))
    stored = [v for chunk in chunks for v in struct.unpack(f"<{4 * ploidy}i", gzip.decompress(chunk.read_bytes()))]
    assert stored[: calls.height * ploidy] == [v for row in want for v in row]
    assert set(stored[calls.height * ploidy :]) == {-2}
    # polars-bio reads back the same values, row by row.
    back = read_back(root).select(pl.col("genotypes").struct.unnest())
    assert [row[0] for row in back["HQ"].to_list()] == want
    # A facet with no haplotype qualities writes no call_HQ.
    assert "call_HQ" not in write_vcf_zarr(loaded_from(tmp_path, RICH, fmt=RICH_FORMAT), tmp_path / "rich.vcz").arrays


def test_haplotype_qualities_not_one_per_gt_slot_are_left_out(tmp_path: Path) -> None:
    """A frame whose `hq` is not parallel to `gt` (no loader writes one) says nowhere which slot each
    value speaks for: written missing and counted, as write_vcf leaves it out of HQ."""
    loaded = hq_loaded(tmp_path)
    calls = loaded.calls.collect()
    odd = pl.col("pos") == 100
    longer = pl.when(odd).then(pl.lit([1, 2, 3], dtype=pl.List(pl.UInt16))).otherwise(pl.col("hq")).alias("hq")
    loaded = loaded.model_copy(update={"calls": loaded.calls.with_columns(longer)})
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    assert "1 row(s) had haplotype qualities that are not one per GT allele; left out of HQ" in written.warnings
    back = read_back(written.path).select(pl.col("genotypes").struct.unnest())
    want = expected_hq(calls.with_columns(pl.when(odd).then(None).otherwise(pl.col("hq")).alias("hq")), calls["gt"].list.len().max())
    assert [row[0] for row in back["HQ"].to_list()] == want


SV_META_45 = (
    '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">',
    '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS">',
    '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END">',
    '##INFO=<ID=CN,Number=A,Type=Float,Description="Copy number of allele">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">',
)
# VCF 4.5 shapes: two offsets per ALT with a missing one, a fractional per-ALT copy number, an ALT
# with no CIEND, a copy-neutral segment's one pair, a sequence ALT with none, and an interval holding
# -1, which VCF Zarr's integer arrays read as missing.
SV_RECORDS_45 = [
    "1\t100\tcnv\tT\t<CNV>,<CNV>\t.\t.\tSVLEN=30,30;CN=1,2.5;CIPOS=-5,5,.,3;CIEND=-10,10,0,4\tGT:CN\t1/2:3.5",
    "1\t200\tdel\tT\t<DEL>\t.\t.\tSVLEN=300;CIPOS=-20,20\tGT\t0/1",
    "1\t600\tins\tT\t<INS>\t.\t.\tSVLEN=40;CIPOS=-1,1;CIEND=3,4\tGT\t0/1",
    "1\t700\tsnv\tA\tG\t.\t.\t.\tGT\t0/1",
    "1\t800\tseg\tN\t.\t.\t.\tEND=900;CIPOS=-50,50\tGT:CN\t0/0:2",
]


def sv_loaded(tmp_path: Path) -> Loaded:
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SYNTHETIC_SAMPLE"])
    path = tmp_path / "sv45.vcf"
    path.write_text("\n".join(["##fileformat=VCFv4.5", "##contig=<ID=1,length=249250621>", *SV_META_45, header, *SV_RECORDS_45]) + "\n")
    return load(path)


def stored_ints(values: list | None, width: int) -> list[int]:
    """A list as an integer array stores it: -1 per missing element (all of them for a null list), -2 past its end."""
    if values is None:
        return [-1] * width
    return [-1 if value is None else value for value in values] + [-2] * (width - len(values))


def test_sv_intervals_and_per_alt_copy_numbers_are_written(tmp_path: Path) -> None:
    loaded = sv_loaded(tmp_path)
    calls = loaded.calls.collect()
    written = write_vcf_zarr(loaded, tmp_path / "out.vcz")
    assert {"variant_CIPOS", "variant_CIEND", "variant_CN"} <= set(written.arrays)
    assert not any("does not write" in warning for warning in written.warnings), written.warnings
    # The interval holding -1 is written missing, counted; the other interval of that row is kept.
    clashing = [any(v in (-1, -2) for v in value if v is not None) if value is not None else False for value in calls["cipos"].to_list()]
    assert sum(clashing) == 1
    assert (
        "1 row(s) with a CIPOS offset of -1 or -2, which VCF Zarr reads as a missing or fill value in an integer "
        "array; their CIPOS written missing"
    ) in written.warnings
    assert not any("CIEND offset" in warning for warning in written.warnings)

    root = written.path
    alleles = calls["alts"].list.len().max() + 1
    back = pb.scan_vcf_zarr(str(root), use_zero_based=False, info_fields=["CIPOS", "CIEND", "CN"]).collect()
    for name, key in (("cipos", "CIPOS"), ("ciend", "CIEND")):
        width = calls[name].list.len().max()
        meta = json.loads((root / f"variant_{key}" / ".zarray").read_text())
        assert meta["shape"] == [calls.height, width] and meta["dtype"] == "<i4" and meta["fill_value"] == -2
        assert json.loads((root / f"variant_{key}" / ".zattrs").read_text())["_ARRAY_DIMENSIONS"] == ["variants", f"{key}_dim"]
        kept = [
            None if value is not None and any(v in (-1, -2) for v in value if v is not None) else value for value in calls[name].to_list()
        ]
        assert back[key].to_list() == [stored_ints(value, width) for value in kept], key

    # CN on alt_alleles, one fewer than variant_allele's width; a float keeps a fractional value, and
    # the spec's NaN sentinels tell a missing value from padding on the bytes.
    width = alleles - 1
    meta = json.loads((root / "variant_CN" / ".zarray").read_text())
    assert meta["shape"] == [calls.height, width] and meta["dtype"] == "<f4"
    assert json.loads((root / "variant_CN" / ".zattrs").read_text())["_ARRAY_DIMENSIONS"] == ["variants", "alt_alleles"]
    raw = gzip.decompress((root / "variant_CN" / "0.0").read_bytes())
    bits = list(struct.unpack(f"<{calls.height * width}I", raw[: 4 * calls.height * width]))
    want_bits = []
    for value in calls["alt_copy_numbers"].to_list():
        row = [0x7F800001] * width if value is None else [struct.unpack("<I", struct.pack("<f", v))[0] for v in value]
        want_bits += row + [0x7F800002] * (width - len(row))
    assert bits == want_bits
    got_cn = [[v for v in row if v == v] or None for row in back["CN"].to_list()]
    assert got_cn == [value or None for value in calls["alt_copy_numbers"].to_list()]
    assert 2.5 in calls["alt_copy_numbers"].explode(empty_as_null=False).to_list()
    # A facet with none of them writes none of these arrays.
    rich = write_vcf_zarr(loaded_from(tmp_path, RICH, fmt=RICH_FORMAT), tmp_path / "rich.vcz")
    assert not {"variant_CIPOS", "variant_CIEND", "variant_CN"} & set(rich.arrays)


def contig_ids(root: Path) -> list[str]:
    return vlen(gzip.decompress((root / "contig_id" / "0").read_bytes()))


@pytest.mark.parametrize("source", sorted(SOURCES))
@pytest.mark.parametrize("style", list(ContigNames))
def test_a_mate_is_named_as_chrom_is_and_listed_in_contig_id(tmp_path: Path, source: str, style: ContigNames) -> None:
    """F32: `write_vcf`'s rule for a breakend mate (test_write_breakends), in the store. The expected
    names come from the source's records and header through `written_name`; the rows are read back by
    polars-bio, and `contig_id` from the bytes, since a contig only an ALT names has no row to show it."""
    contigs, records = SOURCES[source]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(contigs, records))
    root = write_vcf_zarr(load(path), tmp_path / "out.vcz", contig_names=style).path

    back = read_back(root)
    assert back.height == len(records)
    for given, (chrom, alt) in zip(records, back.select("chrom", "alt").iter_rows(), strict=True):
        fields = given.split("\t")
        assert chrom == written_name(style, fields[0], contigs, records)
        given_alts, alts = fields[4].split(","), alt.split(",")
        expected = [(form, bases, written_name(style, mate, contigs, records), pos) for form, bases, mate, pos in mates(given_alts)]
        assert mates(alts) == expected, given
        assert others(alts) == others(given_alts)

    # `contig_id` is the source's contigs, the rows' and the mates', each named once by the style.
    spelled = [*contigs, *(line.split("\t")[0] for line in records)]
    spelled += [mate for line in records for _, _, mate, _ in mates(line.split("\t")[4].split(","))]
    ids = contig_ids(root)
    assert len(ids) == len(set(ids))
    assert set(ids) == {written_name(style, name, contigs, records) for name in spelled}


@pytest.mark.parametrize("style", list(ContigNames))
def test_an_assembly_mate_is_written_verbatim_and_listed_nowhere(tmp_path: Path, style: ContigNames) -> None:
    path = tmp_path / "in.vcf"
    path.write_text("\n".join([*ASSEMBLY, *ASSEMBLY_RECORDS]) + "\n")
    root = write_vcf_zarr(load(path), tmp_path / "out.vcz", contig_names=style).path
    given = [line.split("\t") for line in ASSEMBLY_RECORDS]
    assert read_back(root)["alt"].to_list() == [fields[4] for fields in given]
    assert contig_ids(root) == [written_name(style, "13", ["13"], ASSEMBLY_RECORDS)]


def test_source_refuses_a_contig_only_alts_name_spelled_two_ways(tmp_path: Path) -> None:
    # As `write_vcf` refuses it: contig 3 has no row and no header line, and the ALTs spell it `chr3`
    # and `3`, so which spelling the source meant cannot be told. The other styles write it.
    records = [record("1", 100, "A[chr3:500["), record("1", 200, "A[3:600[")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(["1"], records))
    loaded = load(path)
    with pytest.raises(ValueError, match="spelled 3 more than one way"):
        write_vcf_zarr(loaded, tmp_path / "source.vcz", contig_names=ContigNames.SOURCE)
    assert not (tmp_path / "source.vcz").exists()
    root = write_vcf_zarr(loaded, tmp_path / "chr.vcz", contig_names=ContigNames.CHR).path
    assert read_back(root)["alt"].to_list() == ["A[chr3:500[", "A[chr3:600["]
    assert contig_ids(root) == ["chr1", "chr3"]


def test_a_skipped_rows_mate_is_not_listed(tmp_path: Path) -> None:
    # The second record's GT names an ALT it does not have, so `resolved_only` leaves it out; its mate's
    # contig (7: no row, not in the source header) must not reach `contig_id`.
    records = [record("1", 100, "A[2:500["), record("1", 200, "A[7:500[").replace("0/1", "0/2")]
    path = tmp_path / "in.vcf"
    path.write_text(source_text(["1"], records))
    written = write_vcf_zarr(load(path), tmp_path / "out.vcz", resolved_only=True)
    assert written.rows_skipped == {"genotype_index_past_alts": 1}
    assert read_back(written.path)["alt"].to_list() == ["A[2:500["]
    assert contig_ids(written.path) == ["1", "2"]
