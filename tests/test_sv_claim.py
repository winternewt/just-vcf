"""VCF 4.5's INFO SVCLAIM (`calls.sv_claim`, schema 0.0.24): read per ALT, written back under 4.5,
and honoured by `concordance`'s decompose, which turns a `<DEL>` or `<DUP>` into the breakends it
implies only where it claims a junction.

VCF 4.5 "Encoding Structural Variants" (data/interim/hts-specs/VCFv4.5.tex): `D` is an abundance
claim (copy number changed, no adjacency claimed), `J` an adjacency claim, `DJ` both. `<DEL>` and
`<DUP>` must state one; `.` is the missing value. A file before 4.5 states none, and its `<DEL>` and
`<DUP>` imply both, so a null claim there decomposes as before; on a 4.5 file, or a facet whose load
left SVCLAIM out, a null claim on a `<DEL>` or `<DUP>` is unknown and is not decomposed (audit 16).
`<INV>`'s missing value is J. Synthetic records only; the breakends are those
tests/test_concordance_decompose.py pins to the spec's own (`SPEC_DEL_RECORDS`, `SPEC_DUP_RECORDS`,
`SPEC_INV_RECORDS`).
"""

from pathlib import Path

import polars as pl
import pytest

from just_vcf import Store, load, normalize, rejoin_split, resolve_ref, write_vcf
from just_vcf.concordance import concordance
from just_vcf.loading import Loaded
from just_vcf.schema import CALLS, SvClaim
from test_concordance_decompose import (  # noqa: F401  (store is a fixture)
    META,
    SEQ,
    SPEC_DEL,
    SPEC_DEL_RECORDS,
    SPEC_DUP,
    SPEC_DUP_RECORDS,
    SPEC_INV,
    SPEC_INV_RECORDS,
    base,
    bnd,
    facet,
    store,
    sv,
)
from test_store import SUBJECT

SVCLAIM_LINE = (
    '##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call. '
    'Valid values are D, J, DJ for abundance, adjacency and both respectively">'
)
"""VCF 4.5's own declaration of SVCLAIM (section "Encoding Structural Variants")."""

COLUMNS = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1"])


def v45_text(records: list[str], extra: tuple[str, ...] = ()) -> str:
    meta = ["##fileformat=VCFv4.5", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items()), *META, SVCLAIM_LINE, *extra]
    return "\n".join([*meta, COLUMNS, *records]) + "\n"


def v43_text(records: list[str]) -> str:
    """The same records under VCF 4.3, which defines no SVCLAIM."""
    meta = ["##fileformat=VCFv4.3", *(f"##contig=<ID={c},length={len(s)}>" for c, s in SEQ.items()), *META]
    return "\n".join([*meta, COLUMNS, *records]) + "\n"


def v45(tmp_path: Path, name: str, records: list[str], extra: tuple[str, ...] = ()) -> Path:
    path = tmp_path / f"{name}.vcf"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(v45_text(records, extra))
    return path


def claimed(record: str, claim: str | None) -> str:
    """`record` with `SVCLAIM=claim` added to its INFO (unchanged for None)."""
    if claim is None:
        return record
    fields = record.split("\t")
    fields[7] = f"{fields[7]};SVCLAIM={claim}"
    return "\t".join(fields)


# One `<DEL>` and one `<DUP>` per claim, at positions of their own; a multi-ALT record with a missing
# value; a value that is no claim; and a count that is not one per ALT.
CLAIMS = (None, "D", "J", "DJ")
LOADED_RECORDS = [
    *(claimed(sv(("1", 1000 + 1000 * i, 1500 + 1000 * i), "<DEL>"), claim) for i, claim in enumerate(CLAIMS)),
    *(claimed(sv(("1", 10_000 + 1000 * i, 10_500 + 1000 * i), "<DUP>"), claim) for i, claim in enumerate(CLAIMS)),
    f"1\t20000\tmulti\t{base('1', 20000)}\t<DEL>,<DUP>\t.\tPASS\tEND=20500;SVCLAIM=D,.\tGT\t1/2",
    f"1\t21000\tnot_claim\t{base('1', 21000)}\t<DEL>,<DUP>\t.\tPASS\tEND=21500;SVCLAIM=X,J\tGT\t1/2",
    f"1\t22000\tunfit\t{base('1', 22000)}\t<DEL>,<DUP>\t.\tPASS\tEND=22500;SVCLAIM=DJ\tGT\t1/2",
]


def by_pos(loaded: Loaded) -> dict[int, dict]:
    return {row["pos"]: row for row in loaded.calls.collect().iter_rows(named=True)}


def test_sv_claim_is_read_per_alt_as_written(tmp_path: Path) -> None:
    loaded = load(v45(tmp_path, "claims", LOADED_RECORDS))
    assert loaded.calls.collect_schema()["sv_claim"] == CALLS["sv_claim"]
    rows = by_pos(loaded)
    for offset in (1000, 10_000):
        assert [rows[offset + 1000 * i]["sv_claim"] for i in range(len(CLAIMS))] == [None if c is None else [c] for c in CLAIMS]
    # Per ALT: the missing value is a null claim beside the other ALT's; a value that is no claim too.
    assert rows[20000]["sv_claim"] == [SvClaim.ABUNDANCE.value, None]
    assert rows[21000]["sv_claim"] == [None, SvClaim.ADJACENCY.value]
    # Not one per ALT: the record's claims are not known.
    assert rows[22000]["sv_claim"] is None
    assert "1 INFO SVCLAIM value(s) not one per ALT set to null" in loaded.warnings
    assert "1 INFO SVCLAIM claim(s) neither D, J nor DJ (VCF 4.5) read as no claim" in loaded.warnings
    # Read, so no longer among the fields the load left out.
    assert "INFO/SVCLAIM" not in loaded.facets.collect()["dropped_fields"][0]


def test_an_undeclared_sv_claim_is_declared_and_read(tmp_path: Path) -> None:
    text = v45_text([claimed(sv(SPEC_DEL, "<DEL>"), "D")]).replace(SVCLAIM_LINE + "\n", "")
    path = tmp_path / "undeclared.vcf"
    path.write_text(text)
    loaded = load(path)
    assert loaded.calls.collect()["sv_claim"].to_list() == [["D"]]
    assert any("declared by the loader: INFO/SVCLAIM" in w for w in loaded.warnings), loaded.warnings


def test_a_file_without_sv_claim_has_none(tmp_path: Path) -> None:
    loaded = load(v45(tmp_path, "none", [sv(SPEC_DEL, "<DEL>"), f"1\t30000\t.\t{base('1', 30000)}\tG\t.\tPASS\t.\tGT\t0/1"]))
    assert loaded.calls.collect()["sv_claim"].to_list() == [None, None]


def test_sv_claim_round_trips_through_the_store(tmp_path: Path) -> None:
    loaded = load(v45(tmp_path, "claims", LOADED_RECORDS), subject=SUBJECT)
    Store.open(tmp_path / "store", subject=SUBJECT).add(loaded)
    back = Store.open(tmp_path / "store").loaded().calls.collect()
    assert back.schema["sv_claim"] == CALLS["sv_claim"]
    columns = ["pos", "alts", "sv_claim"]
    assert back.select(columns).sort("pos").equals(loaded.calls.collect().select(columns).sort("pos"))


PSL = '##FORMAT=<ID=PSL,Number=P,Type=String,Description="Phase set list">'


def test_sv_claim_round_trips_through_a_vcf_45_export(tmp_path: Path) -> None:
    """`write_vcf` writes VCF 4.5 when a row holds a PSL phase set, and SVCLAIM with it."""
    snv = f"1\t5\tsnv\t{base('1', 5)}\t{'A' if base('1', 5) != 'A' else 'C'}\t.\tPASS\t.\tGT:PSL\t0|1:.,1*5*1"
    first = load(v45(tmp_path, "in", [snv, *LOADED_RECORDS], (PSL,)))
    written = write_vcf(first, tmp_path / "out.vcf", compress=False)
    text = written.path.read_text()
    assert "##fileformat=VCFv4.5\n" in text
    assert "##INFO=<ID=SVCLAIM,Number=A,Type=String," in text
    assert ";SVCLAIM=D,." in text
    back = load(written.path)
    columns = ["pos", "alts", "sv_claim"]
    assert back.calls.collect().select(columns).equals(first.calls.collect().select(columns))
    assert not any("SVCLAIM" in w for w in written.warnings), written.warnings


def test_a_vcf_43_export_leaves_sv_claim_out_and_counts_it(tmp_path: Path) -> None:
    first = load(v45(tmp_path, "in", LOADED_RECORDS))
    held = first.calls.select(pl.col("sv_claim").list.drop_nulls().list.len().fill_null(0).gt(0).sum()).collect().item()
    written = write_vcf(first, tmp_path / "out.vcf", compress=False)
    text = written.path.read_text()
    assert "##fileformat=VCFv4.3\n" in text and "SVCLAIM" not in text
    assert any(w.startswith(f"{held} row(s) with an SV claim (sv_claim) written as VCF 4.3, which has no SVCLAIM") for w in written.warnings), written.warnings
    assert load(written.path).calls.collect()["sv_claim"].is_null().all()


def test_normalize_splits_sv_claim_with_its_alt_and_rejoin_joins_it(tmp_path: Path, store) -> None:  # noqa: F811
    # normalize splits only sequence sites (a symbolic one keeps its record), whose alleles VCF 4.5 says
    # should write `.`; a claim on one is written here only to tell the two slots apart.
    ref = base("1", 30000)
    g, t = [b for b in "ACGT" if b != ref][:2]
    record = f"1\t30000\t.\t{ref}\t{g},{t}\t.\tPASS\tSVCLAIM=J,.\tGT\t1/2"
    loaded = resolve_ref(load(v45(tmp_path, "multi", [record])), store)
    split = normalize(loaded, store)
    rows = {row["alts"][0]: row["sv_claim"] for row in split.calls.collect().iter_rows(named=True)}
    assert rows == {g: ["J"], t: [None]}
    joined = rejoin_split(split).calls.collect()
    assert joined.select("alts", "sv_claim").rows() == [([g, t], ["J", None])]


def sv_sites(tmp_path: Path, store, a: list[str], b: list[str]) -> tuple[pl.DataFrame, tuple[str, ...]]:  # noqa: F811
    left = resolve_ref(load(v45(tmp_path, "a", a)), store)
    right = facet(tmp_path, "b", b, store)
    result = concordance(left, right)
    return result.sites.collect().filter(pl.col("variant_class") == "sv"), result.warnings


NOT_DECOMPOSED = "1 symbolic deletions, duplications or inversions of facet a are not decomposed"
CLAIM_UNKNOWN = "1 claim unknown: VCF 4.4 and later require SVCLAIM on <DEL>/<DUP>"
CLAIM_NOT_READ = "1 claim not read: the load left INFO SVCLAIM out (dropped_fields)"


def _alone(rows: pl.DataFrame) -> bool:
    """Nothing met: the symbolic allele and the breakends each stand alone."""
    return rows["key"].is_null().all() and sorted(set(rows["agreement"].to_list())) == ["only_a", "only_b"]


@pytest.mark.parametrize(("symbol", "call", "records"), [("<DEL>", SPEC_DEL, SPEC_DEL_RECORDS), ("<DUP>", SPEC_DUP, SPEC_DUP_RECORDS)])
@pytest.mark.parametrize("claim", [*CLAIMS, "."])
def test_only_an_adjacency_claim_is_decomposed_in_a_45_file(tmp_path: Path, store, symbol: str, call: tuple, records: list, claim: str | None) -> None:  # noqa: F811
    """VCF 4.5 requires SVCLAIM on `<DEL>` and `<DUP>`: `J` and `DJ` decompose; `D` states no junction;
    `.` or none leaves the claim unknown, so no junction is implied either (audit 16). Each is counted."""
    rows, warnings = sv_sites(tmp_path, store, [claimed(sv(call, symbol), claim)], [bnd(r) for r in records])
    if claim in ("J", "DJ"):
        assert rows["key"].to_list() == ["decomposed"], rows
        assert not any("not decomposed" in w for w in warnings), warnings
    else:
        assert _alone(rows), rows
        why = "abundance claim only" if claim == "D" else CLAIM_UNKNOWN
        assert any(w.startswith(NOT_DECOMPOSED) and why in w for w in warnings), warnings


def test_a_null_claim_is_unknown_from_44_on(tmp_path: Path, store) -> None:  # noqa: F811
    """VCF 4.4 added SVCLAIM ("Changes between VCFv4.4 and VCFv4.3" in VCFv4.5.tex): a 4.4 `<DEL>`
    without one is unknown, as a 4.5 one is."""
    path = tmp_path / "a.vcf"
    path.write_text(v45_text([sv(SPEC_DEL, "<DEL>")]).replace("##fileformat=VCFv4.5\n", "##fileformat=VCFv4.4\n"))
    left = resolve_ref(load(path), store)
    assert left.facets.collect()["vcf_version"].item() == "VCFv4.4"
    result = concordance(left, facet(tmp_path, "b", [bnd(r) for r in SPEC_DEL_RECORDS], store))
    assert _alone(result.sites.collect().filter(pl.col("variant_class") == "sv"))
    assert any(w.startswith(NOT_DECOMPOSED) and CLAIM_UNKNOWN in w for w in result.warnings), result.warnings


@pytest.mark.parametrize(("symbol", "call", "records"), [("<DEL>", SPEC_DEL, SPEC_DEL_RECORDS), ("<DUP>", SPEC_DUP, SPEC_DUP_RECORDS)])
def test_a_null_claim_before_44_implies_both(tmp_path: Path, store, symbol: str, call: tuple, records: list) -> None:  # noqa: F811
    """VCF 4.3 has no SVCLAIM: its `<DEL>` and `<DUP>` claim both, and decompose (truvari does)."""
    path = tmp_path / "a.vcf"
    path.write_text(v43_text([sv(call, symbol)]))
    left = resolve_ref(load(path), store)
    assert left.facets.collect()["vcf_version"].item() == "VCFv4.3"
    result = concordance(left, facet(tmp_path, "b", [bnd(r) for r in records], store))
    rows = result.sites.collect().filter(pl.col("variant_class") == "sv")
    assert rows["key"].to_list() == ["decomposed"], rows
    assert not any("not decomposed" in w for w in result.warnings), result.warnings


@pytest.mark.parametrize("claim", [None, ".", "J"])
def test_an_inversion_with_no_claim_is_a_junction_in_45(tmp_path: Path, store, claim: str | None) -> None:  # noqa: F811
    """VCF 4.5: `<INV>`'s claim is optional, and `.` equals J, so it decomposes whatever it states of it."""
    rows, warnings = sv_sites(tmp_path, store, [claimed(sv(SPEC_INV, "<INV>"), claim)], [bnd(r) for r in SPEC_INV_RECORDS])
    # Taken once: the inversion meets the adjacency its first implied breakend meets (as
    # tests/test_concordance_decompose.py pins it); its other adjacency stands alone.
    assert rows.filter(pl.col("key") == "decomposed").select("alt", "sv_kind").rows() == [("<INV>", "INV")], rows
    assert not any("not decomposed" in w for w in warnings), warnings


def test_a_facet_that_left_sv_claim_out_does_not_decompose_a_null_claim(tmp_path: Path, store) -> None:  # noqa: F811
    """A facet whose `dropped_fields` lists INFO SVCLAIM (a store migrated from 0.0.23) holds claims it
    never read: a null one is unknown, not decomposed, and counted under its own reason. The facet is a
    4.3 file's, where a null claim would otherwise decompose, so only `dropped_fields` decides."""
    path = tmp_path / "a.vcf"
    path.write_text(v43_text([sv(SPEC_DEL, "<DEL>")]))
    loaded = resolve_ref(load(path), store)
    left = loaded.model_copy(
        update={"facets": loaded.facets.with_columns(pl.concat_list("dropped_fields", pl.lit(["INFO/SVCLAIM"])).list.sort().alias("dropped_fields"))}
    )
    result = concordance(left, facet(tmp_path, "b", [bnd(r) for r in SPEC_DEL_RECORDS], store))
    rows = result.sites.collect().filter(pl.col("variant_class") == "sv")
    assert _alone(rows), rows
    assert any(w.startswith(NOT_DECOMPOSED) and CLAIM_NOT_READ in w for w in result.warnings), result.warnings


def test_each_allele_of_a_record_is_decomposed_by_its_own_claim(tmp_path: Path, store) -> None:  # noqa: F811
    # One record, a `<DEL>` claiming abundance only beside a `<DUP>` claiming a junction: the DUP's
    # breakends meet, the DEL's are not implied.
    chrom, pos, end = SPEC_DEL
    record = f"{chrom}\t{pos}\t.\t{base(chrom, pos)}\t<DEL>,<DUP>\t.\tPASS\tEND={end};SVCLAIM=D,J\tGT\t1/2"
    dup = [(chrom, end, "t[p[", chrom, pos + 1), (chrom, pos + 1, "]p]t", chrom, end)]
    rows, _ = sv_sites(tmp_path, store, [record], [bnd(r) for r in [*SPEC_DEL_RECORDS, *dup]])
    met = rows.filter(pl.col("key") == "decomposed")
    assert met.select("alt", "sv_kind").rows() == [("<DUP>", "DUP")]


def test_the_not_decomposed_warning_names_each_reason_apart(tmp_path: Path, store) -> None:  # noqa: F811
    """Three reasons, each counted on its own: an abundance claim, an END the frame does not hold, and an
    SVTYPE naming another type (DRAGEN's CNV) or a duplication subtype other than TANDEM."""
    a = [
        claimed(sv(("1", 1000, 1500), "<DEL>"), "D"),
        claimed(sv(("1", 3000, 3500), "<DUP>"), "D"),
        sv(("1", 5000, 5500), "<DEL>", svtype="CNV"),
        sv(("1", 7000, 7500), "<DUP:DISPERSED>", svtype="DUP"),
        sv(("1", 9000, 9500), "<INV>"),
    ]
    loaded = resolve_ref(load(v45(tmp_path, "a", a)), store)
    # No loader leaves a symbolic allele's `end_pos` null (it falls back to the REF span); a frame built
    # or edited elsewhere can, and such an allele has no END to join.
    no_end = loaded.model_copy(update={"calls": loaded.calls.with_columns(pl.when(pl.col("pos") == 9000).then(None).otherwise(pl.col("end_pos")).alias("end_pos"))})
    right = facet(tmp_path, "b", [bnd(r) for r in SPEC_DEL_RECORDS], store)
    warnings = [w for w in concordance(no_end, right).warnings if "not decomposed" in w]
    assert len(warnings) == 1, warnings
    line = warnings[0]
    assert line.startswith("5 symbolic deletions, duplications or inversions of facet a are not decomposed into breakends: ")
    assert "2 an abundance claim only (VCF 4.5's SVCLAIM=D)" in line
    assert "1 without END, so no junction to join POS to" in line
    assert "2 an SVTYPE naming another type" in line
