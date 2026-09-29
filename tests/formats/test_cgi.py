"""Complete Genomics native files: detection by `#TYPE` and the `>` header, build hints, real files."""

import bz2
import re
from pathlib import Path

import pytest

from just_vcf.formats import cgi
from just_vcf.formats.spec import Coordinates, FormatSpec
from just_vcf.sniff import SniffOutcome, sniff

# The `#TYPE` value CGI writes for each file kind; a domain constant, not read off a dump.
TYPE_OF = {"cgi_var": "VAR-ANNOTATION", "cgi_master_var": "VAR-OLPL", "cgi_summary": "SUMMARY-REPORT"}
SPEC_OF = {spec.id: spec for spec in cgi.SPECS}

_TYPE_LINE = re.compile(r"^#TYPE\t.*$", re.MULTILINE)
_BUILD_LINE = re.compile(r"^#GENOME_REFERENCE\t.*$", re.MULTILINE)


def _write(tmp_path: Path, text: str, name: str = "head.tsv") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def _header_names(spec: FormatSpec) -> list[str]:
    header = next(line for line in spec.example.splitlines() if line.startswith(">"))
    return header[1:].split("\t")


def _rows(spec: FormatSpec) -> list[dict[str, str]]:
    names = _header_names(spec)
    lines = [line for line in spec.example.splitlines() if line.strip() and not line.startswith(("#", ">"))]
    return [dict(zip(names, line.split("\t"), strict=True)) for line in lines]


def test_the_registry_exports_one_spec_per_type() -> None:
    assert set(SPEC_OF) == set(TYPE_OF)
    assert cgi.FALLBACKS == ()


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_example_detects_as_itself_and_nothing_else(spec: FormatSpec, tmp_path: Path) -> None:
    result = sniff(_write(tmp_path, spec.example))
    assert result.outcome is SniffOutcome.DETECTED
    assert result.format == spec.id
    assert result.family is spec.family
    assert result.also_matched == ()


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_example_header_line_names_the_declared_columns(spec: FormatSpec) -> None:
    assert _header_names(spec) == [column["name"] for column in spec.columns]


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_bz2_wrapped_example_is_detected_through_the_layer(spec: FormatSpec, tmp_path: Path) -> None:
    path = tmp_path / "head.tsv"
    path.write_bytes(bz2.compress(spec.example.encode()))
    result = sniff(path)
    assert result.format == spec.id
    assert result.layers == ("bz2",)


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
@pytest.mark.parametrize("other_type", sorted(TYPE_OF.values()))
def test_type_line_decides_the_spec(spec: FormatSpec, other_type: str, tmp_path: Path) -> None:
    """A head carrying another kind's `#TYPE` is not this kind; with a mismatched `>` line it is nothing."""
    text = _TYPE_LINE.sub(f"#TYPE\t{other_type}", spec.example)
    result = sniff(_write(tmp_path, text))
    if other_type == TYPE_OF[spec.id]:
        assert result.format == spec.id
    else:
        assert result.outcome is SniffOutcome.UNKNOWN
        assert result.format is None


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_without_a_type_line_nothing_is_claimed(spec: FormatSpec, tmp_path: Path) -> None:
    text = _TYPE_LINE.sub("", spec.example)
    assert sniff(_write(tmp_path, text)).outcome is SniffOutcome.UNKNOWN


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_type_line_on_a_headerless_body_is_not_enough(spec: FormatSpec, tmp_path: Path) -> None:
    text = "\n".join(line for line in spec.example.splitlines() if not line.startswith(">")) + "\n"
    assert sniff(_write(tmp_path, text)).outcome is SniffOutcome.UNKNOWN


@pytest.mark.parametrize(
    ("header", "type_value", "expected"),
    [
        # Format 1.6/1.7 var: totalScore instead of the 2.x scores, 12 columns.
        (
            ">locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\ttotalScore\thapLink\txRef",
            "VAR-ANNOTATION",
            "cgi_var",
        ),
        # Format 2.0/2.2 var: varQuality, no alleleFreq or alternativeCalls, 14 columns.
        (
            ">locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\tvarScoreVAF"
            "\tvarScoreEAF\tvarQuality\thapLink\txRef",
            "VAR-ANNOTATION",
            "cgi_var",
        ),
        # Format 2.0 masterVar: allele1/2VarQuality, no Freq, AlternativeCalls or LAF columns (35 columns;
        # only the leading ones matter for detection).
        (
            ">locus\tploidy\tchromosome\tbegin\tend\tzygosity\tvarType\treference\tallele1Seq\tallele2Seq"
            "\tallele1VarScoreVAF\tallele2VarScoreVAF\tallele1VarScoreEAF\tallele2VarScoreEAF"
            "\tallele1VarQuality\tallele2VarQuality",
            "VAR-OLPL",
            "cgi_master_var",
        ),
    ],
    ids=["var-1.x", "var-2.0", "mastervar-2.0"],
)
def test_older_format_versions_are_detected(header: str, type_value: str, expected: str, tmp_path: Path) -> None:
    width = header.count("\t") + 1
    if expected == "cgi_var":
        row = ["1", "2", "all", "chr1", "0", "10000", "no-ref", "=", "?"]
    else:
        row = ["1", "2", "chr1", "0", "10000", "no-call", "no-ref", "=", "?", "?"]
    text = (
        "#GENERATED_BY\tcgatools\n"
        "#FORMAT_VERSION\t1.6\n"
        "#GENOME_REFERENCE\tNCBI build 37\n"
        f"#TYPE\t{type_value}\n"
        "\n"
        f"{header}\n" + "\t".join(row + [""] * (width - len(row))) + "\n"
    )
    result = sniff(_write(tmp_path, text))
    assert result.format == expected
    assert result.build == "GRCh37"


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
@pytest.mark.parametrize(("stated", "build"), [("NCBI build 37", "GRCh37"), ("NCBI build 36", "GRCh36")])
def test_build_comes_from_the_genome_reference_line(spec: FormatSpec, stated: str, build: str, tmp_path: Path) -> None:
    text = _BUILD_LINE.sub(f"#GENOME_REFERENCE\t{stated}", spec.example)
    result = sniff(_write(tmp_path, text))
    assert result.format == spec.id
    assert result.build == build
    assert result.build_evidence == f"#GENOME_REFERENCE\t{stated}"


@pytest.mark.parametrize("spec", cgi.SPECS, ids=lambda s: s.id)
def test_build_is_unknown_without_the_line_and_withheld_on_conflict(spec: FormatSpec, tmp_path: Path) -> None:
    silent = sniff(_write(tmp_path, _BUILD_LINE.sub("", spec.example), "silent.tsv"))
    assert silent.format == spec.id
    assert silent.build is None

    both = _BUILD_LINE.sub("#GENOME_REFERENCE\tNCBI build 37\n#GENOME_REFERENCE\tNCBI build 36", spec.example)
    conflicted = sniff(_write(tmp_path, both, "conflicted.tsv"))
    assert conflicted.format == spec.id
    assert conflicted.build is None
    assert conflicted.warnings


def test_build_hints_are_header_statements_only(tmp_path: Path) -> None:
    """`#DBSNP_BUILD` and `#GENE_ANNOTATIONS` mention builds of other things; they are not genome builds."""
    spec = SPEC_OF["cgi_master_var"]
    text = _BUILD_LINE.sub("#GENE_ANNOTATIONS\tNCBI build 37.2", spec.example)
    result = sniff(_write(tmp_path, text))
    assert result.format == spec.id
    assert result.build is None


@pytest.mark.parametrize("spec_id", ["cgi_var", "cgi_master_var"])
def test_example_rows_obey_the_half_open_convention(spec_id: str) -> None:
    """A spelled reference spans `end - begin` bases; an insertion has zero length."""
    spec = SPEC_OF[spec_id]
    assert spec.coordinates is Coordinates.ZERO_BASED_HALF_OPEN
    rows = _rows(spec)
    spelled = [row for row in rows if row["reference"] not in ("=", "")]
    assert spelled, "the example needs at least one row with a spelled reference"
    for row in rows:
        begin, end = int(row["begin"]), int(row["end"])
        assert 0 <= begin <= end
        if row["varType"] == "ins":
            assert begin == end
        if row in spelled:
            assert end - begin == len(row["reference"])


# ---------------------------------------------------------------------------------------------------
# Opt-in: real CGI files under $JUST_VCF_SAMPLES_DIR/pgp-hms.

_REAL_GLOBS = {
    "cgi_var": "var-*-ASM.tsv.bz2",
    "cgi_master_var": "masterVarBeta-*.tsv*",
    "cgi_summary": "summary-*.tsv*",
}



def _stated_header(path: Path) -> tuple[str | None, str | None]:
    """`#TYPE` and `#GENOME_REFERENCE` values, read independently of sniff (bz2 by magic)."""
    with path.open("rb") as raw:
        magic = raw.read(3)
    opener = bz2.open if magic == b"BZh" else open
    stated: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                break
            key, _, value = line.rstrip("\r\n").partition("\t")
            stated.setdefault(key, value)
    return stated.get("#TYPE"), stated.get("#GENOME_REFERENCE")


@pytest.mark.parametrize("spec_id", sorted(_REAL_GLOBS))
def test_real_cgi_files_detect_as_their_type(pgp_participants: Path, spec_id: str) -> None:
    root = pgp_participants
    paths = sorted(p for p in root.glob(f"*/*/{_REAL_GLOBS[spec_id]}") if not p.name.endswith(".part"))
    if not paths:
        pytest.skip(f"no {_REAL_GLOBS[spec_id]} files under {root}")
    builds = {"NCBI build 37": "GRCh37", "NCBI build 36": "GRCh36"}

    wrong: list[str] = []
    for path in paths:
        result = sniff(path)
        stated_type, stated_build = _stated_header(path)
        rel = str(path.relative_to(root))
        if stated_type != TYPE_OF[spec_id]:
            wrong.append(f"{rel}: header says #TYPE {stated_type!r}")
        if result.format != spec_id or result.also_matched:
            wrong.append(f"{rel}: {result.outcome.value} {result.format} also={result.also_matched}")
        if result.build != builds.get(stated_build or ""):
            wrong.append(f"{rel}: build {result.build} but header states {stated_build!r}")
    assert wrong == []
