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

from just_vcf import UnsortedRowsError, UnwritableRowsError, load, write_vcf
from just_vcf.loading import Loaded
from just_vcf.vcf_zarr import write_vcf_zarr
from test_load_cgi import LINKED, master_var_text
from test_load_cgi import write as write_cgi
from test_write_vcf import RICH, RICH_FORMAT, loaded_from

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
