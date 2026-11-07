"""`scripts/id_orientation.py`: the per-probe I/D orientation table learned from export and whole-genome
pairs (RM281), on synthetic files only: tiny 23andMe exports, a CGI `var` genome, a small manifest."""

import json
from pathlib import Path

import id_orientation as script
import polars as pl
import pytest
from test_load_cgi import Locus, Seg, var_text, write
from test_probes import deletion, insertion, probe_row, write_manifest
from test_reference import PREAMBLE_23, SEQ, array_text, write_fasta

from just_vcf import ReferenceStore, load
from just_vcf.probes import (
    ORIENTATION_CHIPS,
    ORIENTATION_VARIANT_WINDOW,
    IdOrientation,
    Orientation,
    orientation_tables,
)

CONTIG = SEQ["1"]
INSERTED = "TTGCAG"


@pytest.fixture
def store(tmp_path: Path) -> Path:
    s = ReferenceStore.open(tmp_path / "grch37")
    s.add_fasta(write_fasta(tmp_path / "grch37.fa", SEQ), assembly_label="GRCh37")
    return tmp_path / "grch37"


@pytest.fixture
def manifests(tmp_path: Path) -> Path:
    """rs100 a deletion site, rs200 an insertion site, and three more deletions: rs500 under a no-call,
    rs600 beside a variant, rs800 read both ways."""
    rows = [
        probe_row("rs100", "[I/D]", "1", 501, "PLUS", deletion(500, 3)[1]),
        probe_row("rs200", "[D/I]", "1", 1001, "PLUS", insertion(1000, INSERTED)),
        probe_row("rs500", "[I/D]", "1", 2001, "PLUS", deletion(2000, 2)[1]),
        probe_row("rs600", "[I/D]", "1", 1501, "PLUS", deletion(1500, 2)[1]),
        probe_row("rs800", "[I/D]", "1", 2501, "PLUS", deletion(2500, 2)[1]),
    ]
    write_manifest(tmp_path / "manifests", rows)
    return tmp_path / "manifests"


# The sites as the export's rows resolve them (1-based, inclusive): rs100 [501, 504] (three deleted
# bases and the pad), rs200 [1001, 1001], rs600 [1501, 1503], rs500 [2001, 2003], rs800 [2501, 2503].
RS100_END, RS600_END = 504, 1503
FAR_SNP = RS100_END + ORIENTATION_VARIANT_WINDOW + 1  # just outside rs100's window: rs100 stays evidence
NEAR_SNP = RS600_END + ORIENTATION_VARIANT_WINDOW  # just inside rs600's window: rs600 is no evidence


def snp(position: int) -> Locus:
    """A het SNV at 1-based `position`: a CGI call carrying an ALT."""
    base = CONTIG[position - 1]
    other = next(b for b in "ACGT" if b != base)
    return Locus("chr1", position - 1, position, haps=((Seg(position - 1, position, "snp", base, other),), (Seg(position - 1, position, "ref", base, base),)))


GENOME = [
    Locus("chr1", 0, FAR_SNP - 1, uniform="ref"),
    snp(FAR_SNP),
    Locus("chr1", FAR_SNP, NEAR_SNP - 1, uniform="ref"),
    snp(NEAR_SNP),
    Locus("chr1", NEAR_SNP, 1990, uniform="ref"),
    Locus("chr1", 1990, 2010, uniform="no-call"),
    Locus("chr1", 2010, 3000, uniform="ref"),
]


def export(tmp_path: Path, name: str, codes: dict[str, str]) -> Path:
    positions = {"rs100": 501, "rs200": 1001, "rs500": 2001, "rs600": 1501, "rs800": 2501}
    path = tmp_path / f"{name}.txt"
    path.write_text(array_text([(rsid, "1", positions[rsid], gt) for rsid, gt in codes.items()], PREAMBLE_23))
    return path


def test_the_table_is_learned_from_reference_blocks_only(
    tmp_path: Path, store: Path, manifests: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    genome = write(tmp_path, var_text(GENOME), "var-GS000000001-ASM.tsv.bz2")
    # Over a reference block, a hom-ref reading (I at a deletion site, D at an insertion site) is
    # `same` evidence, a hom-alt reading `flipped` evidence.
    same_and_flipped = {"rs100": "II", "rs200": "II", "rs500": "DD", "rs600": "II", "rs800": "II"}
    people = {
        "hu1": [same_and_flipped],
        "hu2": [same_and_flipped],
        "hu3": [same_and_flipped, same_and_flipped],  # two exports of one person count once
        "hu4": [{"rs100": "DI", "rs800": "DD"}],  # a het is no evidence; rs800 now read both ways
    }
    lines = ["participant\tchip\texport\tgenome"]
    for person, exports in people.items():
        for index, codes in enumerate(exports):
            lines.append(f"{person}\t23andme_v4\t{export(tmp_path, f'{person}-{index}', codes)}\t{genome}")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("\n".join(lines) + "\n")
    out = tmp_path / "id_orientation.json"

    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")  # restored after the test; the script sets it from --manifests
    script.main(pairs=pairs, manifests=manifests, reference=store, out=out, source="synthetic pairs")
    payload = json.loads(out.read_text())
    tables = orientation_tables(payload)
    assert tables["23andme_v4"] == {
        "rs100": IdOrientation(Orientation.SAME, 3, 0),
        "rs200": IdOrientation(Orientation.FLIPPED, 0, 3),
    }, "rs500 has no reference block, rs600 a variant within the window, rs800 evidence both ways"
    assert set(tables) == ORIENTATION_CHIPS and tables["23andme_v3"] == {}
    assert payload["chips"]["23andme_v4"]["pairs"] == 5 and payload["chips"]["23andme_v3"]["pairs"] == 0
    assert payload["provenance"]["learned_from"].startswith("synthetic pairs: 5 pair(s)") and "4 people" in payload["provenance"]["learned_from"]

    logged = capsys.readouterr().err  # the script's own handler writes to stderr
    assert "23andme_v4: hom_ref over a reference block: 8" in logged  # rs100 by 4 exports, rs800 by 4 (hu4's DD is hom-alt)
    for private in ("rs100", "rs200", "501", str(tmp_path), "hu1"):
        assert private not in logged, f"the log names {private!r}"

    first = out.read_text()
    script.main(pairs=pairs, manifests=manifests, reference=store, out=out, source="synthetic pairs")
    assert out.read_text() == first, "the same pairs write the same bytes"


def test_the_ungated_resolution_reads_codes_as_written_whatever_the_chip(tmp_path: Path, store: Path, manifests: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The learner reads a v3 or v4 export's codes though `resolve_ref` withholds them: the probe
    resolution itself, with `I` read as Illumina's `I`."""
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", str(manifests))
    loaded = load(export(tmp_path, "one", {"rs100": "II", "rs200": "II", "rs500": "DI"}))
    relabelled = loaded.model_copy(update={"facets": loaded.facets.with_columns(chip=pl.lit("23andme_v4"))})
    sites, counts = script.homozygous_sites(relabelled, store)
    assert sites.select("rsid", "pos", "end_pos", "kind").rows() == [("rs100", 501, 504, "hom_ref"), ("rs200", 1001, 1001, "hom_alt")]
    assert counts["code_rows"] == 3 and counts["code_resolved"] == 3 and counts["homozygous"] == 2


def test_a_pairs_file_without_its_columns_is_refused(tmp_path: Path) -> None:
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("participant\texport\tgenome\nhu1\ta\tb\n")
    with pytest.raises(ValueError, match="chip"):
        script.read_pairs(pairs)
