"""The other_tables family: Illumina array reports, imputation probabilities, annotator JSON, clinical
lists, GET-Evidence GFF and pre-VCF pipeline tables.

Synthetic heads pin the shapes; the opt-in real-file tests at the bottom check the specs against PGP-HMS
downloads under `JUST_VCF_SAMPLES_DIR` and skip when it is unset.
"""

import gzip
import os
import zipfile
from collections.abc import Iterator
from pathlib import Path

import pytest
from dotenv import dotenv_values, find_dotenv

from just_vcf.formats import other_tables
from just_vcf.formats.spec import FormatSpec
from just_vcf.sniff import SniffOutcome, Sniffed, sniff

OWN_SPECS: tuple[FormatSpec, ...] = (*other_tables.SPECS, *other_tables.FALLBACKS)
OWN_IDS: frozenset[str] = frozenset(spec.id for spec in OWN_SPECS)


def _sniff_text(tmp_path: Path, text: str, name: str = "head") -> Sniffed:
    path = tmp_path / name
    path.write_bytes(text.encode("utf-8"))
    return sniff(path)


def _crlf(text: str) -> str:
    return text.replace("\n", "\r\n")


# --- Every spec detects as itself, and only as itself within the family ---------------------------


@pytest.mark.parametrize("spec", OWN_SPECS, ids=lambda s: s.id)
def test_example_detects_as_itself_and_no_sibling(spec: FormatSpec, tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, spec.example or "")
    assert result.outcome is SniffOutcome.DETECTED, result.reason
    assert result.format == spec.id
    assert result.family is spec.family
    assert not (set(result.also_matched) & OWN_IDS), "two specs of this family claim one head"


@pytest.mark.parametrize(
    "spec",
    [
        other_tables.ILLUMINA_FINAL_REPORT,
        other_tables.ILLUMINA_FINAL_REPORT_HEADERLESS,
        other_tables.COLOR_VARIANT_SET,
        other_tables.FTDNA_BIGY,
        other_tables.CASAVA_SNP_TABLE,
        other_tables.CASAVA_INDEL_TABLE,
    ],
    ids=lambda s: s.id,
)
def test_crlf_heads_detect_the_same(spec: FormatSpec, tmp_path: Path) -> None:
    """Real exports of these formats are CRLF; detection must not depend on the line ending."""
    result = _sniff_text(tmp_path, _crlf(spec.example or ""))
    assert result.format == spec.id


def test_gzipped_final_report_detects_through_the_layer(tmp_path: Path) -> None:
    path = tmp_path / "report.csv.gz"
    path.write_bytes(gzip.compress(_crlf(other_tables.ILLUMINA_FINAL_REPORT_HEADERLESS.example or "").encode()))
    result = sniff(path)
    assert result.format == other_tables.ILLUMINA_FINAL_REPORT_HEADERLESS.id
    assert result.layers == ("gzip",)


def test_comma_final_report_with_a_column_line_is_the_headed_spec(tmp_path: Path) -> None:
    head = (
        "[Header],,,\n"
        "GSGT Version,2.0.4,,\n"
        "Content,,GSA-24v2-0_A1.bpm,\n"
        "[Data],,,\n"
        "SNP Name,Sample ID,Allele1 - Forward,Allele2 - Forward\n"
        "rs3094315,S1,A,G\n"
        "rs3131972,S1,-,-\n"
    )
    assert _sniff_text(tmp_path, head).format == other_tables.ILLUMINA_FINAL_REPORT.id


def test_older_color_layout_detects_as_color(tmp_path: Path) -> None:
    head = (
        "SampleID,Chromosome,Pos,Ref,Alt,Gene,Type,cHGVS,pHGVS,Exon,Transcript,Zygosity,isStructuralVariant,"
        "Classification,TestOrdered,ReferenceVersion\n"
        "S1,1,55523033,A,G,PCSK9,SNV,c.1026A>G,p.Gln342=,7,ENST00000302118,Homozygous,no,Benign,pgx v1,GRCh37\n"
    )
    result = _sniff_text(tmp_path, head)
    assert result.format == other_tables.COLOR_VARIANT_SET.id
    assert result.build == "GRCh37"


def test_multi_member_gen_zip_is_one_format(tmp_path: Path) -> None:
    """Gencove ships one `.gen` per chromosome in a zip; every member is the same format."""
    path = tmp_path / "id.gen.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("part_a.gen", other_tables.OXFORD_GEN.example or "")
        zf.writestr("part_b.gen", "--- rs180734498 13302 C T 0.860 0.139 0\n")
    result = sniff(path)
    assert result.format == other_tables.OXFORD_GEN.id
    assert {member.format for member in result.members} == {other_tables.OXFORD_GEN.id}
    assert not any("different formats" in warning for warning in result.warnings)


# --- Near misses -----------------------------------------------------------------------------------

NEAR_MISSES: dict[str, str] = {
    # `[Header]` block with no `[Data]` section: a truncated or foreign key/value file.
    "final_report_without_data": "[Header]\nGSGT Version\t1.8.4\nContent\t\tHumanOmni2.5-8v1_A.bpm\n",
    # `[Data]` followed by a line that is neither a column line nor a five-field call row.
    "final_report_other_data": "[Header],,\nContent,,x.bpm\n[Data],,\nfoo,bar,baz\n",
    # Two probabilities per site is not an Oxford .gen row.
    "gen_two_probabilities": "--- rs58108140 10583 G A 0.638 0.359\n--- rs189107123 10611 C G 0.981 0.019\n",
    # Three trailing fields that are not probabilities.
    "gen_letters_not_probabilities": "--- rs58108140 10583 G A x y z\n",
    # Color's header with a column out of order.
    "color_reordered": "SampleID,Chromosome,Position,Alt,Ref,Gene,Type,cHGVS,pHGVS\nS1,1,55523033,G,A,PCSK9,SNV,c.1,p.1\n",
    # Color's discovery genotypes: a different (rsid) table from the same vendor.
    "color_discovery_genotypes": (
        "SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion\nS1,1,752566,rs3094315,AG,GRCh37\n"
    ),
    # JSON that is not the annotator's.
    "json_other_producer": '{"header":{"tool":"other"},"positions":[\n{"chromosome":"chr1","position":752566}\n]}\n',
    # GFF3 without a `##genome-build` line or an `alleles` attribute.
    "plain_gff3": "##gff-version 3\nchr1\tensembl\tgene\t65419\t71585\t.\t+\t.\tID=gene1;Name=OR4F5\n",
    # A `##genome-build` line over rows with a different attribute vocabulary.
    "gff_build_line_foreign_rows": "##genome-build 37\nchr1\tsrc\tgene\t65419\t71585\t.\t+\t.\tID=gene1\n",
    # SOAPsnp's documented 17 columns, without the distance column this spec is anchored on.
    "soapsnp_17_columns": "chr1\t742429\tG\tR\t34\tG\t36\t12\t12\tA\t30\t10\t10\t22\t0.412\t1.00000\t1\n",
    # A six-column BED of intervals, the shape of a callable-region file.
    "bed_intervals": "1\t14507\t14512\t.\t.\t.\n1\t14513\t14582\t.\t.\t.\n1\t14583\t14589\t.\t.\t.\n",
    # CASAVA SNP header with a column renamed.
    "casava_renamed_column": "Index\tSample\tChr Id\tPosition\tA Bases\tC Bases\tG Bases\tT Bases\tCall\n1\tS1\t1\t1\t0\t0\t0\t0\tA\n",
    # FTDNA Y-SNP panel: names and a derived flag, no position or allele.
    "ftdna_ysnp_panel": '"SnpName","Positive","Test"\n"L21","False","Deep Clade-R"\n"P312","True","Deep Clade-R"\n',
    # FTDNA Y-STR repeat counts.
    "ftdna_ystr": 'DYS393,DYS390,DYS19,DYS391\n" 13"," 24"," 14"," 11"\n',
    # The rCRS and RSRS reference sequences side by side.
    "mtdna_reference_table": "Position,rCRS,RSRS\n1,G,G\n2,A,A\n73,A,G\n",
    # 23andMe ancestry composition segments.
    "ancestry_segments": "Ancestry,Copy,Chromosome,Start Point,End Point\nBritish & Irish,1,chr1,752566,1005806\n",
}


@pytest.mark.parametrize("head", NEAR_MISSES.values(), ids=NEAR_MISSES.keys())
def test_near_misses_are_not_this_family(head: str, tmp_path: Path) -> None:
    result = _sniff_text(tmp_path, head)
    assert result.format not in OWN_IDS
    assert not (set(result.also_matched) & OWN_IDS)


# --- Build hints -----------------------------------------------------------------------------------

BUILD_CASES: dict[str, tuple[str, str | None, str]] = {
    "gff_build_36": (
        "##genome-build 36\nchr1\t-\t-\t742429\t742429\t.\t+\t.\talleles A/G\n",
        "GRCh36",
        other_tables.GET_EVIDENCE_GFF.id,
    ),
    "gff_build_b37": (other_tables.GET_EVIDENCE_GFF.example or "", "GRCh37", other_tables.GET_EVIDENCE_GFF.id),
    "nirvana_states_assembly": (
        '{"header":{"annotator":"Nirvana 2.0.9","creationTime":"2018-01-01 00:00:00","genomeAssembly":"GRCh38",'
        '"schemaVersion":6},"positions":[\n{"chromosome":"chr1","refAllele":"G","position":817186}\n]}\n',
        "GRCh38",
        other_tables.NIRVANA_JSON.id,
    ),
    # Schema v4 names phyloP `hg19` as a data source; that is not the assembly.
    "nirvana_data_source_is_not_a_build": (
        '{"header":{"annotator":"Illumina Annotation Engine 1.5.1.0","dataSources":[{"name":"phyloP",'
        '"version":"hg19"}]},"positions":[\n',
        None,
        other_tables.NIRVANA_JSON.id,
    ),
    "color_reference_version": (other_tables.COLOR_VARIANT_SET.example or "", "GRCh37", other_tables.COLOR_VARIANT_SET.id),
    # The manifest name is a vendor convention, not a build statement.
    "final_report_manifest_is_not_a_build": (
        other_tables.ILLUMINA_FINAL_REPORT.example or "",
        None,
        other_tables.ILLUMINA_FINAL_REPORT.id,
    ),
    "gen_states_nothing": (other_tables.OXFORD_GEN.example or "", None, other_tables.OXFORD_GEN.id),
    "soapsnp_states_nothing": (other_tables.SOAPSNP_CNS.example or "", None, other_tables.SOAPSNP_CNS.id),
}


@pytest.mark.parametrize("case", BUILD_CASES.values(), ids=BUILD_CASES.keys())
def test_build_is_what_the_file_states(case: tuple[str, str | None, str], tmp_path: Path) -> None:
    head, build, format_id = case
    result = _sniff_text(tmp_path, head)
    assert result.format == format_id
    assert result.build == build


def test_color_rows_naming_two_builds_give_unknown(tmp_path: Path) -> None:
    head = (other_tables.COLOR_VARIANT_SET.example or "") + (
        "S1,1,55524237,G,A,PCSK9,SNV,c.1420G>A,p.Val474Ile,ENST00000302118,,Heterozygous,,Benign,,"
        "Color Hereditary Disease Risk Test 2.0,GRCh38,2025-06-23\n"
    )
    result = _sniff_text(tmp_path, head)
    assert result.format == other_tables.COLOR_VARIANT_SET.id
    assert result.build is None
    assert any("conflicting builds" in warning for warning in result.warnings)


# --- Opt-in: real PGP-HMS files ------------------------------------------------------------------


def _participants_dir() -> Path:
    samples = os.environ.get("JUST_VCF_SAMPLES_DIR") or dotenv_values(find_dotenv(usecwd=True)).get(
        "JUST_VCF_SAMPLES_DIR"
    )
    if not samples:
        pytest.skip("JUST_VCF_SAMPLES_DIR is unset")
    root = Path(samples) / "pgp-hms" / "participants"
    if not root.is_dir():
        pytest.skip(f"{root} is absent")
    return root


def _glob(root: Path, pattern: str) -> list[Path]:
    """Finished downloads matching `pattern` at `<participant>/<collection>/`, sorted."""
    return sorted(p for p in root.glob(f"*/*/{pattern}") if p.is_file() and not p.name.endswith(".part"))


def _who(root: Path, path: Path) -> str:
    """The participant id only: file names can carry a person's name."""
    return path.relative_to(root).parts[0]


def _head_lines(path: Path, limit: int = 64) -> list[str]:
    with path.open("rb") as raw:
        opener = gzip.open if raw.read(2) == b"\x1f\x8b" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        return [line.rstrip("\r\n") for _, line in zip(range(limit), handle)]


def _nonempty(paths: list[Path], what: str) -> Iterator[Path]:
    if not paths:
        pytest.skip(f"no {what} on disk")
    yield from paths


def test_real_final_reports_split_on_the_line_after_data() -> None:
    root = _participants_dir()
    candidates = [*_glob(root, "*.txt"), *_glob(root, "*.csv"), *_glob(root, "*.csv.gz")]
    reports = [p for p in candidates if _head_lines(p, 1)[:1] and _head_lines(p, 1)[0].startswith("[Header]")]
    for path in _nonempty(reports, "[Header] files"):
        lines = _head_lines(path)
        after_data = lines[next(i for i, line in enumerate(lines) if line.startswith("[Data]")) + 1]
        expected = (
            other_tables.ILLUMINA_FINAL_REPORT.id
            if "SNP Name" in after_data
            else other_tables.ILLUMINA_FINAL_REPORT_HEADERLESS.id
        )
        result = sniff(path)
        assert result.format == expected, _who(root, path)
        assert result.build is None, _who(root, path)


def test_real_gen_zips_are_oxford_gen_in_every_member() -> None:
    root = _participants_dir()
    for path in _nonempty(_glob(root, "*.gen.zip"), ".gen.zip files"):
        result = sniff(path)
        assert result.format == other_tables.OXFORD_GEN.id, _who(root, path)
        with zipfile.ZipFile(path) as zf:
            member_count = sum(1 for info in zf.infolist() if not info.is_dir())
        assert len(result.members) == member_count
        assert {member.format for member in result.members} == {other_tables.OXFORD_GEN.id}


def test_real_color_variant_sets_state_their_build() -> None:
    root = _participants_dir()
    for path in _nonempty(_glob(root, "*_variant_set.csv"), "Color variant_set files"):
        lines = _head_lines(path)
        header = lines[0].split(",")
        column = header.index("ReferenceVersion")
        rows = [line.split(",") for line in lines[1:]]
        stated = {row[column].split(".")[0] for row in rows if len(row) == len(header) and row[0]}
        result = sniff(path)
        assert result.format == other_tables.COLOR_VARIANT_SET.id, _who(root, path)
        assert result.build == (stated.pop() if len(stated) == 1 else None)


def test_real_color_zip_leads_with_the_variant_set() -> None:
    root = _participants_dir()
    for path in _nonempty(_glob(root, "Color_Health.zip"), "Color_Health.zip files"):
        result = sniff(path)
        assert result.format == other_tables.COLOR_VARIANT_SET.id
        assert result.member is not None and result.member.endswith("_variant_set.csv")
        others = [member for member in result.members if member.member != result.member]
        assert not any(member.format == other_tables.COLOR_VARIANT_SET.id for member in others)


def test_real_gffs_carry_the_build_line_they_state() -> None:
    root = _participants_dir()
    gffs = [p for p in _glob(root, "*.gz") if _head_lines(p, 1)[:1] and _head_lines(p, 1)[0].startswith("##genome-build")]
    for path in _nonempty(gffs, "GFF files with a ##genome-build line"):
        stated = _head_lines(path, 1)[0].split()[1].removeprefix("b")
        result = sniff(path)
        assert result.format == other_tables.GET_EVIDENCE_GFF.id, _who(root, path)
        assert result.build == f"GRCh{stated}"


def test_real_nirvana_json_detects_and_states_no_build_it_lacks() -> None:
    root = _participants_dir()
    for path in _nonempty(_glob(root, "*.json.gz"), "JSON downloads"):
        first = _head_lines(path, 1)[0]
        if '"annotator":"Illumina Annotation Engine' not in first and '"annotator":"Nirvana' not in first:
            continue
        result = sniff(path)
        assert result.format == other_tables.NIRVANA_JSON.id, _who(root, path)
        assert (result.build is None) == ('"genomeAssembly"' not in first)


def test_real_casava_tables_are_one_format_per_kind() -> None:
    root = _participants_dir()
    snp_formats = {sniff(p).format for p in _nonempty(_glob(root, "*_Chr*_All.txt"), "CASAVA SNP tables")}
    assert snp_formats == {other_tables.CASAVA_SNP_TABLE.id}
    indel_formats = {sniff(p).format for p in _nonempty(_glob(root, "*_All_InDels.txt"), "CASAVA indel tables")}
    assert indel_formats == {other_tables.CASAVA_INDEL_TABLE.id}


def test_real_bed_files_are_not_claimed() -> None:
    root = _participants_dir()
    for path in _nonempty(_glob(root, "*.bed.gz"), "BED files"):
        result = sniff(path)
        assert result.format not in OWN_IDS, _who(root, path)
