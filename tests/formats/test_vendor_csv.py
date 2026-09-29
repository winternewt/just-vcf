"""Vendor genotype tables that are not 23andMe-shaped: detection, near misses, builds, file variants.

Synthetic heads are written to tmp_path and run through `sniff`. The opt-in real-file test reads
PGP-HMS exports under `JUST_VCF_SAMPLES_DIR` and skips when they are absent.
"""

import bz2
import gzip
import io
import zipfile
from pathlib import Path

import pytest

from just_vcf.formats import REGISTRY, get_spec, vendor_csv
from just_vcf.sniff import Sniffed, SniffOutcome, sniff

VENDOR_IDS = frozenset(spec.id for spec in (*vendor_csv.SPECS, *vendor_csv.FALLBACKS))
RSID_RESULT_IDS = frozenset({"ftdna", "ftdna_space", "rsid_result_csv"})


def _example(format_id: str) -> str:
    example = get_spec(format_id).example
    assert example is not None
    return example


def _sniff_text(tmp_path: Path, text: str, name: str = "head.txt") -> Sniffed:
    path = tmp_path / name
    path.write_bytes(text.encode())
    return sniff(path)


def _crlf(text: str) -> str:
    return text.replace("\n", "\r\n")


@pytest.mark.parametrize("spec", [*vendor_csv.SPECS, *vendor_csv.FALLBACKS], ids=lambda s: s.id)
def test_example_detects_as_itself_with_crlf_too(spec, tmp_path: Path) -> None:
    for text in (spec.example, _crlf(spec.example)):
        result = _sniff_text(tmp_path, text)
        assert (result.outcome, result.format, result.family) == (SniffOutcome.DETECTED, spec.id, spec.family)


def test_registry_places_vendor_signatures_before_rsid_tables_and_fallbacks_last() -> None:
    order = [spec.id for spec in REGISTRY]
    for spec in vendor_csv.SPECS:
        assert order.index(spec.id) < order.index("23andme")
    first_fallback = min(order.index(spec.id) for spec in vendor_csv.FALLBACKS)
    assert first_fallback > order.index("23andme")


def test_23andme_files_are_not_captured(tmp_path: Path) -> None:
    for text in (_example("23andme"), _crlf(_example("23andme"))):
        result = _sniff_text(tmp_path, text)
        assert result.format == "23andme"
        assert VENDOR_IDS.isdisjoint(result.also_matched)


def test_23andme_layout_under_an_ancestry_name_is_not_ancestrydna(tmp_path: Path) -> None:
    """Converted files keep the vendor name in the file name and nothing else; content decides."""
    body = "\n".join(line for line in _example("23andme").splitlines() if not line.startswith("# This data file"))
    result = _sniff_text(tmp_path, "# This file emulates the format of 23andMe\n" + body + "\n", "AncestryDNA.txt")
    assert result.format not in VENDOR_IDS


# --- the RSID,CHROMOSOME,POSITION,RESULT layout -------------------------------------------------


def test_ftdna_quoted_and_unquoted_rows_are_both_ftdna(tmp_path: Path) -> None:
    quoted = _example("ftdna")
    unquoted = quoted.replace('"', "")
    assert unquoted != quoted
    for text in (quoted, unquoted, _crlf(quoted), _crlf(unquoted)):
        result = _sniff_text(tmp_path, text)
        assert result.format == "ftdna"
        assert result.build is None and result.build_evidence is None


def test_ftdna_x_file_with_haploid_single_letters(tmp_path: Path) -> None:
    text = 'RSID,CHROMOSOME,POSITION,RESULT\n"rs1000530","X","78124078","T"\n"rs1001874","X","86808168","C"\n'
    assert _sniff_text(tmp_path, text).format == "ftdna"


def test_space_rows_under_a_comma_header_are_ftdna_space_not_ftdna(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _crlf(_example("ftdna_space")))
    assert result.format == "ftdna_space"
    assert "ftdna" not in result.also_matched
    assert get_spec("ftdna_space").delimiter != get_spec("ftdna").delimiter


def test_commented_rsid_result_file_is_not_claimed_as_ftdna(tmp_path: Path) -> None:
    """FTDNA exports carry no comment lines; one with a preamble is left vendor-neutral."""
    text = "# raw data export\n# genotypes on the forward strand\n" + _example("ftdna")
    result = _sniff_text(tmp_path, text)
    assert result.format == "rsid_result_csv"
    assert "ftdna" not in result.also_matched


def test_a_myheritage_mention_blocks_the_ftdna_claim(tmp_path: Path) -> None:
    text = _example("ftdna") + '"rs11240777","1","798959","GG" MyHeritage\n'
    result = _sniff_text(tmp_path, text)
    assert result.format == "rsid_result_csv"


def test_unplaced_gsa_probes_first_fall_back_to_the_neutral_spec(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _example("rsid_result_csv"))
    assert result.format == "rsid_result_csv"
    assert not ({"ftdna", "ftdna_space"} & set(result.also_matched))
    assert get_spec("rsid_result_csv").vendor == "unknown"


def test_ftdna_example_also_matches_the_fallback_but_the_signature_wins(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _example("ftdna"))
    assert result.format == "ftdna"
    assert "rsid_result_csv" in result.also_matched


def test_rsid_result_header_must_be_exact(tmp_path: Path) -> None:
    for header in ("RSID,CHROMOSOME,POSITION,GENOTYPE", "rsid,chromosome,position,result", "RSID\tCHROMOSOME\tPOSITION\tRESULT"):
        text = _example("ftdna").replace("RSID,CHROMOSOME,POSITION,RESULT", header)
        result = _sniff_text(tmp_path, text)
        assert result.format not in RSID_RESULT_IDS


# --- AncestryDNA --------------------------------------------------------------------------------


def test_ancestrydna_build_comes_from_the_preamble(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _example("ancestrydna"))
    assert (result.format, result.build) == ("ancestrydna", "GRCh37")
    assert "37.1" in result.build_evidence


def test_ancestrydna_without_the_build_line_has_unknown_build(tmp_path: Path) -> None:
    text = "\n".join(line for line in _example("ancestrydna").splitlines() if "reference build" not in line) + "\n"
    result = _sniff_text(tmp_path, text)
    assert result.format == "ancestrydna"
    assert result.build is None


def test_ancestrydna_with_mixed_line_endings(tmp_path: Path) -> None:
    lines = _example("ancestrydna").splitlines()
    text = "".join(line + ("\n" if i % 2 else "\r\n") for i, line in enumerate(lines))
    result = _sniff_text(tmp_path, text)
    assert (result.format, result.build) == ("ancestrydna", "GRCh37")


def test_ancestry_preamble_over_a_23andme_body_is_not_ancestrydna(tmp_path: Path) -> None:
    """The data shape guards the vendor line: four columns with a joined genotype are not AncestryDNA."""
    preamble = "".join(line + "\n" for line in _example("ancestrydna").splitlines() if line.startswith("#"))
    body = "rs4477212\t1\t82154\tAA\nrs3131972\t1\t752721\tAG\n"
    result = _sniff_text(tmp_path, preamble + body)
    assert result.format != "ancestrydna"


# --- Genos --------------------------------------------------------------------------------------


def test_genos_build_and_half_open_coordinates(tmp_path: Path) -> None:
    spec = get_spec("genos")
    result = _sniff_text(tmp_path, _example("genos"))
    assert (result.format, result.build) == ("genos", "GRCh37")
    assert "-/-" not in spec.no_call and "-" not in spec.no_call
    rows = [line.split("\t") for line in spec.example.splitlines() if not line.startswith("#")]
    assert [c["name"] for c in spec.columns][1:3] == ["begin", "end"]
    assert all(int(end) - int(begin) in (0, 1) for _, begin, end, _, _ in rows)


def test_genos_without_its_first_line_is_not_genos(tmp_path: Path) -> None:
    text = "\n".join(_example("genos").splitlines()[1:]) + "\n"
    assert _sniff_text(tmp_path, text).format != "genos"


def test_genos_unrecognised_reference_leaves_build_unknown(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _example("genos").replace("UCSC.hg19", "UCSC.hg18"))
    assert result.format == "genos"
    assert result.build is None


# --- Color --------------------------------------------------------------------------------------


def test_color_with_disclaimer_trailer_and_crlf(tmp_path: Path) -> None:
    trailer = (
        "Raw data disclaimer:,,,,,\n"
        '"This data is provided for research use only, and it cannot be relied on for any clinical purpose.",,,,,\n'
    )
    for text in (_example("color_genotypes") + trailer, _crlf(_example("color_genotypes") + trailer)):
        result = _sniff_text(tmp_path, text)
        assert (result.format, result.build) == ("color_genotypes", "GRCh37")


def test_color_letter_prefixed_kit_ids(tmp_path: Path) -> None:
    text = _crlf(_example("color_genotypes").replace("10000000000000,", "L1000000,"))
    result = _sniff_text(tmp_path, text)
    assert (result.format, result.build) == ("color_genotypes", "GRCh37")


def test_color_build_follows_the_reference_version_column(tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, _example("color_genotypes").replace(",GRCh37", ",GRCh38"))
    assert result.format == "color_genotypes"
    assert result.build is None


def test_color_variant_set_is_not_the_genotype_format(tmp_path: Path) -> None:
    text = (
        "SampleID,Chromosome,Pos,Ref,Alt,Gene,Type,cHGVS,pHGVS,Exon,Transcript,Zygosity,"
        "isStructuralVariant,Classification,TestOrdered,ReferenceVersion\n"
        "10000000000000,12,112241766,G,A,ALDH2,SNV,c.1510G>A,p.Glu504Lys,12,ENST00000261733,"
        "Heterozygous,no,Benign,pgx v1,GRCh37\n"
    )
    assert _sniff_text(tmp_path, text).format not in VENDOR_IDS


# --- deCODEme -----------------------------------------------------------------------------------


def test_decodeme_zip_members_carry_data_and_build_separately(tmp_path: Path) -> None:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("genotypes.csv", _crlf(_example("decodeme")))
        zf.writestr("info.txt", _crlf(_example("decodeme_info")))
        zf.writestr("signature.txt", "Signature,0_1\r\nChecksum,0\r\n")
    path = tmp_path / "decodeme.zip"
    path.write_bytes(buffer.getvalue())
    result = sniff(path)
    assert result.format == "decodeme"
    assert result.build is None, "the data member states no build"
    by_format = {m.format: m for m in result.members if m}
    assert by_format["decodeme_info"].build == "GRCh36"
    assert by_format["decodeme_info"].family.value == "metadata"


def test_decodeme_strand_column_is_required(tmp_path: Path) -> None:
    text = _example("decodeme").replace("Position,Strand,YourCode", "Position,YourCode")
    assert _sniff_text(tmp_path, text).format != "decodeme"


# --- Genographic --------------------------------------------------------------------------------


def test_genographic_header_lf_rows_crlf(tmp_path: Path) -> None:
    header, *rows = _example("genographic_geno2").splitlines()
    text = header + "\n" + "".join(row + "\r\n" for row in rows)
    result = _sniff_text(tmp_path, text)
    assert result.format == "genographic_geno2"
    assert result.build is None
    assert get_spec("genographic_geno2").coordinates is None


def test_genographic_needs_allele_columns(tmp_path: Path) -> None:
    text = "SNP,Chr,Allele1,Allele2\nM173,TRUE,DeepSNP\n"
    assert _sniff_text(tmp_path, text).format != "genographic_geno2"


def test_ftdna_y_snp_booleans_are_not_a_genotype_table(tmp_path: Path) -> None:
    text = '"SnpName","Positive","Test"\r\n"M173","True","Deep Clade-R"\r\n"M269","False","Deep Clade-R"\r\n'
    result = _sniff_text(tmp_path, text)
    assert result.format not in VENDOR_IDS


# --- opt-in: real PGP-HMS exports ---------------------------------------------------------------

_REAL_FILE_CAP = 256 << 20
_TWENTYTHREE_FIRST_LINE = "# This data file generated by 23andMe"
_SKIP_SUFFIXES = {".part", ".bam", ".bai", ".cram", ".pdf", ".jpg", ".jpeg", ".png", ".fastq", ".fq"}

_FIRST_LINE_FORMATS: dict[str, frozenset[str]] = {
    "#AncestryDNA raw data download": frozenset({"ancestrydna"}),
    "##FileFormat=Genos": frozenset({"genos"}),
    "SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion": frozenset({"color_genotypes"}),
    "Name,Variation,Chromosome,Position,Strand,YourCode": frozenset({"decodeme"}),
    "SNP,Chr,Allele1,Allele2": frozenset({"genographic_geno2"}),
    "RSID,CHROMOSOME,POSITION,RESULT": RSID_RESULT_IDS,
}



def _first_lines(path: Path) -> list[str]:
    """First line of a plain, gzip or bz2 file, or of each zip member: read independently of `sniff`."""
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic[:2] == b"\x1f\x8b":
        with gzip.open(path) as handle:
            heads = [handle.readline()]
    elif magic[:3] == b"BZh":
        with bz2.open(path) as handle:
            heads = [handle.readline()]
    elif magic == b"PK\x03\x04":
        with zipfile.ZipFile(path) as zf:
            heads = []
            for info in zf.infolist():
                if info.is_dir() or "__MACOSX" in info.filename:
                    continue
                with zf.open(info) as member:
                    heads.append(member.readline())
    else:
        with path.open("rb") as handle:
            heads = [handle.readline()]
    return [head.decode("latin-1").rstrip("\r\n") for head in heads]


def _flatten(result: Sniffed) -> list[Sniffed]:
    return [result, *(m for member in result.members for m in _flatten(member))]


def test_real_pgp_exports(pgp_participants: Path) -> None:
    participants = pgp_participants
    paths = sorted(
        p
        for p in participants.glob("*/upload-*/*")
        if p.is_file() and p.suffix.lower() not in _SKIP_SUFFIXES and p.stat().st_size <= _REAL_FILE_CAP
    )
    checked: dict[str, int] = {}
    checked_23andme = 0
    for path in paths:
        try:
            firsts = _first_lines(path)
        except (OSError, EOFError, ValueError, zipfile.BadZipFile):
            continue
        expected = [_FIRST_LINE_FORMATS[line] for line in firsts if line in _FIRST_LINE_FORMATS]
        is_23andme = [line.startswith(_TWENTYTHREE_FIRST_LINE) for line in firsts]
        if not expected and not any(is_23andme):
            continue
        everything = _flatten(sniff(path))
        if all(is_23andme):
            # Precision: a 23andMe export is never taken by a vendor_csv spec, not even as a runner-up.
            assert all(VENDOR_IDS.isdisjoint({r.format, *r.also_matched}) for r in everything), path.name
            checked_23andme += 1
            continue
        results = [r for r in everything if r.format in VENDOR_IDS]
        detected = {r.format for r in results}
        # Every member whose first line is a vendor header is detected as one of that header's formats.
        for formats in expected:
            assert detected & formats, f"{path.name}: expected one of {sorted(formats)}, got {sorted(detected)}"
        for r in results:
            spec = get_spec(r.format)
            assert r.family is spec.family
            # A build is only ever what a declared hint read from the file.
            hint_builds = {h.build for h in spec.build_hints}
            assert (r.build is None) or (r.build in hint_builds and r.build_evidence), path.name
            if r.format in RSID_RESULT_IDS:
                assert r.build is None
            checked[r.format] = checked.get(r.format, 0) + 1
        if "decodeme" in detected:
            # The data member states no build; the info member beside it does.
            info = [r for r in results if r.format == "decodeme_info"]
            assert info and all(r.build is not None for r in info), path.name
            assert all(r.build is None for r in results if r.format == "decodeme"), path.name
    if not checked:
        pytest.skip("no vendor exports of this family under the samples dir")
    if not checked_23andme:
        pytest.skip("no 23andMe export available as a negative control")


def test_geno2_without_its_column_line_detects_by_shape(tmp_path) -> None:
    from_spec = get_spec("genographic_geno2").example.split("\n", 1)[1]
    path = tmp_path / "geno2.csv"
    path.write_text(from_spec)
    assert sniff(path).format == "genographic_geno2_headerless"
    path.write_text("\n" + from_spec)
    assert sniff(path).format == "genographic_geno2_headerless"
    path.write_text("marker,note,A,B\nfoo,1,A,G\n")
    assert sniff(path).format != "genographic_geno2_headerless"
