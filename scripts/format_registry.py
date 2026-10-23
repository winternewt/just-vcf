"""Write docs/FORMAT_REGISTRY.md from the format specs, so the spec never drifts from the code.

Every format section is rendered from its `FormatSpec` in `just_vcf.formats.REGISTRY`, in detection
order; the recognition rules come from `just_vcf.sniff`'s own constants and tables; which formats load
comes from the loaders' `FORMATS` sets. The section on broken and odd shapes is the `SHAPES` table
below: each entry's synthetic example is a test vector (`tests/test_format_registry.py` sniffs it and
requires the outcome the entry states), and each test it cites must exist.

    uv run python scripts/format_registry.py            # rewrite the doc
    uv run python scripts/format_registry.py --check    # exit 1 when the doc is stale

The whole file is generated; there is no hand-edited region. Rejected: generated blocks between
markers with hand-written prose around them, because the prose about broken shapes would then be free
text whose examples nobody runs.
"""

import bz2
import logging
import re
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import typer

from just_vcf import LOADABLE_FORMATS
from just_vcf.formats import BCF_FORMAT, REGISTRY, cgi, other_tables, rsid_tables, vcf, vendor_csv
from just_vcf.formats.spec import Coordinates, Family, FormatSpec, HeaderMode
from just_vcf.loaders import cgi as cgi_loader
from just_vcf.loaders import keyvalue as keyvalue_loader
from just_vcf.loaders import rsid_table as rsid_table_loader
from just_vcf.loaders import table as table_loader
from just_vcf.loaders import vcf as vcf_loader
from just_vcf.sniff import (
    _BINARY_KINDS,
    _TEXT_KINDS,
    _ZIP_JUNK,
    _ZIP_METHODS,
    DATA_SHAPE_LINES,
    HEAD_BYTE_CAP,
    HEAD_DATA_LINES,
    MAX_DEPTH,
    MAX_MEMBERS,
    NESTED_ZIP_CAP,
    TAR_SCAN_CAP,
)

logger = logging.getLogger("format_registry")

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "FORMAT_REGISTRY.md"

FAMILY_MEANING: dict[Family, str] = {
    Family.GENOTYPE_TABLE: "genotypes at predefined sites (arrays, imputation output); an absent site is unknown",
    Family.VARIANT_CALLS: "discovery calls in VCF or a VCF-like layout",
    Family.CALLSET_NATIVE: "a vendor's native whole-genome callset that states its reference and no-call regions",
    Family.ANNOTATED_VARIANTS: "a variant list whose purpose is annotation (reports, GFF, annotator JSON)",
    Family.METADATA: "describes a callset without holding its calls",
}

COORDINATES_MEANING: dict[Coordinates, str] = {
    Coordinates.ONE_BASED: "a position names the base itself, as VCF `POS` does",
    Coordinates.ZERO_BASED_HALF_OPEN: (
        "`begin` counts the bases before the span and `end` is exclusive, so a single base "
        "at VCF `POS` p is `begin = p - 1, end = p`"
    ),
}

HEADER_MEANING: dict[HeaderMode, str] = {
    HeaderMode.NONE: "the file names no columns; they are the ones listed below, in file order",
    HeaderMode.LINE: "the first line that is not a comment",
    HeaderMode.COMMENTED: "the last comment line before the data (`# rsid\\tchromosome…`)",
    HeaderMode.VCF: "VCF meta lines and a `#CHROM` line, as the VCF specification defines them",
    HeaderMode.SECTION: "the first line after a `[<section>]` line",
}

LOADERS: tuple[tuple[str, frozenset[str]], ...] = (
    ("just_vcf.loaders.vcf", vcf_loader.FORMATS),
    ("just_vcf.loaders.rsid_table", rsid_table_loader.FORMATS),
    ("just_vcf.loaders.table", table_loader.FORMATS),
    ("just_vcf.loaders.keyvalue", keyvalue_loader.FORMATS),
    ("just_vcf.loaders.cgi", cgi_loader.FORMATS),
)

MODULES = (vcf, cgi, other_tables, vendor_csv, rsid_tables)
"""The family modules, to name where each spec is declared. Their order here is not the registry's."""

# Constructs beyond POSIX extended regular expressions, with what they mean, so another tool can
# translate the patterns. Only the ones that occur in the registry are listed in the doc.
REGEX_CONSTRUCTS: tuple[tuple[str, str, str], ...] = (
    (r"\A", r"\\A", "start of the head (not of a line)"),
    (r"(?i:", r"\(\?i:", "a group matched case-insensitively"),
    (r"(?<!", r"\(\?<!", "negative lookbehind"),
    (r"(?!", r"\(\?!", "negative lookahead; `\\A(?!)` can never match"),
    (r"\b", r"\\b", "word boundary"),
    (r"\d", r"\\d", "a digit"),
    (r"\s", r"\\s", "whitespace"),
    (r"\S", r"\\S", "anything but whitespace"),
    (r"\w", r"\\w", "a letter, digit or underscore"),
    (r"(?:", r"\(\?:", "a group that captures nothing"),
)


# ---------------------------------------------------------------------------------------------------
# Broken and odd shapes


@dataclass(frozen=True)
class Expect:
    """What `sniff` answers on an entry's example. `build` is checked only when `check_build`."""

    outcome: str
    format: str | None = None
    build: str | None = None
    check_build: bool = False
    layers: tuple[str, ...] = ()
    reason_contains: str | None = None


@dataclass(frozen=True)
class Shape:
    """One broken or odd shape met in real files, with a synthetic example and what just-vcf does."""

    key: str
    title: str
    findings: tuple[str, ...]
    shape: str
    handling: str
    status: str
    """`handled`, `refused` (on purpose, with a named error), or `open` (a gap, said where)."""
    example: str | None = None
    members: tuple[tuple[str, str], ...] = ()
    """For an archive example: (member name, text) pairs, written into one zip."""
    layers: tuple[str, ...] = ()
    """Layers the test wraps the example in: `bz2`, or `zip_deflate64` (a zip whose member says method 9)."""
    example_note: str | None = None
    expect: Expect | None = None
    code: tuple[str, ...] = ()
    tests: tuple[str, ...] = ()


STATUSES = ("handled", "refused", "open")


def _cgi_row(*fields: str, width: int) -> str:
    """One CGI data line padded with empty fields to `width` columns."""
    return "\t".join((*fields, *([""] * (width - len(fields))))) + "\n"


_CGI_VAR_25_HEAD = (
    "#FORMAT_VERSION\t2.5\n"
    "#GENERATED_BY\tcgatools\n"
    "#GENOME_REFERENCE\tNCBI build 37\n"
    "#SOFTWARE_VERSION\t2.5.0.20\n"
    "#TYPE\tVAR-ANNOTATION\n"
    "\n"
    ">locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\tvarScoreVAF\t"
    "varScoreEAF\tvarFilter\thapLink\txRef\talleleFreq\talternativeCalls\n"
)

_GT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

SHAPES: tuple[Shape, ...] = (
    # --- findings F1-F17 ----------------------------------------------------------------------------
    Shape(
        key="deflate64-zip",
        title="A zip member compressed with Deflate64",
        findings=("F1",),
        shape=(
            "A zip archive whose member is stored with compression method 9 (Deflate64), which Windows "
            "Explorer writes for members over 2 GiB. The archive's directory reads; the member does not "
            "inflate with a plain Deflate decoder. PPMd (method 98) and AES encryption (99) fail the same way."
        ),
        example=(
            "##fileformat=VCFv4.2\n"
            + _GT
            + _COLUMNS
            + "1\t752566\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
        ),
        layers=("zip_deflate64",),
        example_note="The member's text; its local and central zip headers both say compression method 9.",
        expect=Expect("unreadable", reason_contains="Deflate64 (method 9)"),
        handling=(
            "`sniff` answers `unreadable` with the member's name and the reason `zip compression Deflate64 "
            "(method 9) is not supported`; `load` refuses anything `sniff` did not detect. Reading it needs a "
            "Deflate64 inflater, which is a dependency decision not yet taken."
        ),
        status="open",
        code=("sniff.py `_ZIP_METHODS`",),
        tests=("test_unsupported_zip_method_is_unreadable_and_named",),
    ),
    Shape(
        key="doubled-tab",
        title="A 23andMe row with a doubled tab after the id",
        findings=("F2",),
        shape=(
            "Some 2012-2014 23andMe exports write a few dozen rows with an empty field between the id and the "
            "chromosome: five fields where the layout has four."
        ),
        example=(
            "# This data file generated by 23andMe at: Tue Jan  1 00:00:00 2013\n"
            "#\n"
            "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
            "#\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs4477212\t1\t82154\tAA\n"
            "rs3094315\t\t1\t752566\tAG\n"
            "rs3131972\t1\t752721\tGG\n"
            "rs12562034\t1\t768448\t--\n"
        ),
        expect=Expect("detected", "23andme", "GRCh37", check_build=True),
        handling=(
            "Detected as `23andme`: one well-formed row among the first data lines is enough. The loader reads "
            "a row of exactly five fields whose second field is empty as if that field were not there, and "
            "counts those rows in a warning of their own. An empty field anywhere else still drops the row, "
            "counted."
        ),
        status="handled",
        code=("loaders/rsid_table.py `read_calls`",),
        tests=("test_a_doubled_tab_after_the_id_is_read_past_and_counted",),
    ),
    Shape(
        key="appended-export",
        title="Two exports in one headerless file",
        findings=("F3",),
        shape=(
            "A headerless file that holds a 23andMe-shaped export (id, chromosome, position, genotype) and "
            "then a second export of the same person in the two-allele layout (id, chromosome, position, "
            "allele1, allele2; chromosomes 23-26 for X, Y, XY and MT; `0` for no call), with no line between "
            "them. The file can also open with one stray row that has an empty field before the genotype."
        ),
        example=(
            "rs548049170\t1\t69869\t\tTT\n"
            "rs4477212\t1\t82154\tAG\n"
            "rs3094315\t1\t752566\tAA\n"
            "rs3131972\t1\t752721\tAG\n"
            "rs4477212\t1\t82154\tA\tG\n"
            "rs3094315\t1\t752566\tA\tA\n"
            "rs3131972\t1\t752721\t0\t0\n"
        ),
        example_note=(
            "Shortened: in the file it describes each half runs to hundreds of thousands of rows, and the "
            "loader reads the second half only when it is at least 1,000 rows long."
        ),
        expect=Expect("detected", "rsid_table_headerless"),
        handling=(
            "Detected as `rsid_table_headerless`. A run of at least 1,000 rows of one other width that ends the "
            "file, after other rows, and whose first five rows are in the two-allele layout, is read by that "
            "layout's rules as a second facet in the same `Loaded`: `relation = appended_export`, `part = 1` in "
            "its facet id, the file's `format_id` kept. Its build is not the one the first export states: it "
            "stays null until `resolve_ref` infers it from the rows. A shorter tail, or one in no known layout, "
            "is warned about and not read. The stray first row is dropped and counted."
        ),
        status="handled",
        code=("loaders/rsid_table.py `Appended`, `_appended_start`, `_read_appended`",),
        tests=(
            "test_an_appended_export_in_another_layout_is_its_own_facet",
            "test_a_tail_that_is_short_or_in_no_known_layout_is_not_read",
            "test_the_build_a_preamble_states_is_the_first_exports_only",
        ),
    ),
    Shape(
        key="run-together-indels",
        title="Indel genotypes written as the two alleles run together",
        findings=("F4", "F17"),
        shape=(
            "Living DNA writes an indel call as its two allele sequences back to back, with nothing between "
            "them (`CCA` for `C` and `CA`, `GTTGTT` for `GTT` twice). Where one allele ends and the other "
            "starts is not written."
        ),
        example=(
            "# Living DNA customer genotype data download file version: 1.0.2\n"
            "# Genotype chip: Sirius\n"
            "# Human Genome Reference Build 37 (GRCh37.p13).\n"
            "# Genotypes are presented on the forward strand.\n"
            "#\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs3131972\t1\t752721\tAG\n"
            "rs11240777\t1\t798959\tCCA\n"
            "rs4970383\t1\t838555\tGTTGTT\n"
            "rs4475691\t1\t846808\tCC\n"
        ),
        expect=Expect("detected", "living_dna", "GRCh37", check_build=True),
        handling=(
            "The loader keeps the text in `source_gt` and leaves `gt_alleles` null. `resolve_ref` then splits a "
            "text of three or more bases where exactly one split fits: both pieces start with the same base and "
            "one of them equals the reference over its own span from `pos`. REF and ALT follow only where the "
            "split says which allele is the reference; a homozygous call (`CA/CA` cannot say whether the site "
            "is `CA>C` or `C>CA`) and a het whose two alleles both match the reference keep REF unread. Texts "
            "that fit several splits (repeats) or none stay as loaded, counted. A facet whose strand check says "
            "`mixed` splits nothing. Settling the rest needs the chip's own site alleles from its manifest, "
            "which the library does not read (F17)."
        ),
        status="open",
        code=("reference.py `_split_concatenated`, `CONCATENATED_ALLELES`",),
        tests=(
            "test_living_dna_indel_genotypes_split_against_the_reference",
            "test_only_living_dna_texts_are_split",
            "test_a_mixed_strand_living_dna_file_splits_nothing",
        ),
    ),
    Shape(
        key="archived-vcf",
        title="A VCF inside another compression or archive layer",
        findings=("F6",),
        shape=(
            "A VCF compressed with bz2 or xz, or inside a zip or tar, or BGZF under a plain `.vcf` name. The "
            "name says nothing reliable about the bytes."
        ),
        example=(
            "##fileformat=VCFv4.2\n"
            "##reference=file:///ref/hs37d5.fa\n"
            + _GT
            + _COLUMNS
            + "1\t752566\trs3094315\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
        ),
        layers=("bz2",),
        example_note="The payload; the test compresses it with bz2.",
        expect=Expect("detected", "vcf", "GRCh37", check_build=True, layers=("bz2",)),
        handling=(
            "`sniff` peels the layers by magic bytes and reports them in `layers` (here `bz2`). polars-bio reads "
            "plain, gzip and BGZF VCF directly; any other layer is decompressed into the loader's cache as a "
            "BGZF copy (a plain-text copy of a 1.27 GB bz2 gVCF was 14 GB) and read from there."
        ),
        status="handled",
        code=("sniff.py `_sniff_layer`", "loaders/vcf.py `readable_path`"),
        tests=("test_compression_is_detected_by_magic_not_extension", "test_archived_vcf_is_rewritten_and_loaded"),
    ),
    Shape(
        key="spreadsheet-vcf",
        title="A VCF re-saved by a spreadsheet",
        findings=("F7",),
        shape=(
            "A VCF opened in a spreadsheet and saved again as tab-separated text: every meta line is padded with "
            "tabs to the row width, any field holding a comma or a quote is wrapped in quotes with its own "
            "quotes doubled (so a meta line with a `Description` starts `\"##`), and lines end in a bare CR."
        ),
        example=(
            "##fileformat=VCFv4.2\t\t\t\t\t\t\t\t\t\r"
            '"##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"\t\t\t\t\t\t\t\t\t\r'
            '"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=""Allelic depths"">"\t\t\t\t\t\t\t\t\t\r'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\r"
            '1\t752566\trs3094315\tG\tA\t50.1\tPASS\t.\tGT:AD\t"0/1:12,9"\r'
            "1\t752721\trs3131972\tA\tG\t.\tPASS\t.\tGT\t./.\r"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "Detected as `vcf`: `sniff` reads CR line ends as line ends. The loader looks at the first 64 KiB "
            "of the payload: a line starting `\"##`, or CRs and no LF, marks a re-save. It then reads the file "
            "by universal newlines, cuts each line at tabs, unwraps quoted fields and undoes the doubled quotes, "
            "drops the padding from `##` lines, and polars-bio reads the rewritten copy. A warning says so, and "
            "that a number whose text the spreadsheet rewrote (`6980.90` as `6980.9`) cannot be told. A separate "
            "format id was rejected: the delimiter is still a tab, and the quoting is a container's mangling, "
            "like CRLF."
        ),
        status="handled",
        code=("loaders/vcf.py `is_resaved`, `undo_resave`, `_unquote`",),
        tests=(
            "test_a_vcf_resaved_by_a_spreadsheet_loads_as_the_text_it_was_saved_from",
            "test_a_resaved_vcf_with_a_field_past_the_csv_limit_loads",
            "test_crlf_and_plain_vcfs_are_not_taken_for_a_resave",
        ),
    ),
    Shape(
        key="two-hash-column-line",
        title="A VCF column line written with two `#`",
        findings=("F7", "F15"),
        shape=(
            "Some SV and CNV VCFs write the column line as `##CHROM` and put a third `#` on most meta lines "
            "(`###contig`, `###FORMAT`, `###reference`), while `##fileformat` and `##ALT` are written right. "
            "Their CNV rows write ALT as the bare letters `CNV`, without angle brackets."
        ),
        example=(
            "##fileformat=VCFv4.2\n"
            '##ALT=<ID=DEL,Description="Deletion">\n'
            "###contig=<ID=1,length=249250621>\n"
            '###FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t752566\t.\tG\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=753566\tGT\t0/1\n"
            "1\t846808\t.\tC\tCNV\t.\tPASS\tEND=860000\tGT\t0/1\n"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "Detected as `vcf`. Where no `#CHROM` line exists, the loader takes a meta line that is VCF's column "
            "line with one `#` too many as the column line and reads the file through a rewritten copy. A `###` "
            "line stays what VCF makes of it, a meta line whose key starts with `#`: no contig or FORMAT is "
            "declared from it, so GT and END are declared by the loader as for any undeclared key (the build "
            "hints still read its text). An ALT that is not `ACGTN` bases and not a symbolic `<…>` allele, such "
            "as `CNV`, is kept verbatim, gets no variant class, is not counted as a small-variant site, and is "
            "counted in one warning (F15)."
        ),
        status="handled",
        code=("loaders/vcf.py `read_header`, `HASHED_COLUMN_LINE`", "callset.py `variant_class_expr`"),
        tests=(
            "test_a_column_line_written_with_two_hashes_is_read",
            "test_an_alt_that_is_no_vcf_allele_has_no_variant_class",
        ),
    ),
    Shape(
        key="sites-only-vcf",
        title="A sites-only VCF",
        findings=("F8",),
        shape="A VCF with the eight fixed columns and no FORMAT or sample column, as annotation sources write.",
        example=(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t752566\trs3094315\tG\tA\t.\t.\t.\n"
            "1\t752721\trs3131972\tA\tG\t.\t.\t.\n"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "Detected as `vcf`. `load` refuses it with `NoSampleError` before reading anything: just-vcf holds a "
            "person's calls, and a frame of null genotypes would look like a person with nothing called. "
            "`vcf_samples` returns `[]`."
        ),
        status="refused",
        code=("loading.py `_load_vcf`",),
        tests=("test_sites_only_vcf_is_refused",),
    ),
    Shape(
        key="no-build-stated",
        title="An array export that states no build",
        findings=("F9",),
        shape=(
            "FTDNA exports, preamble-less rsid tables, PLINK output and Final Reports name no build, and an "
            "array table has no `##contig` lengths to infer one from. FTDNA shipped both build 36 and build 37 "
            "exports in one layout."
        ),
        example=(
            "RSID,CHROMOSOME,POSITION,RESULT\n"
            '"rs4477212","1","72017","AA"\n'
            '"rs3094315","1","742429","AG"\n'
            '"rs3131972","1","742584","GG"\n'
        ),
        example_note="Positions are build 36 ones; nothing in the file says so.",
        expect=Expect("detected", "ftdna", None, check_build=True),
        handling=(
            "`sniff` reports `build = None`: no vendor default is ever assumed. `resolve_ref` infers the build "
            "of such a facet from its own rows: for each reference collection it considers, the share of sampled "
            "SNV rows whose called alleles include the reference base there. The winner must clear a floor and "
            "a margin over the others; it is recorded as `assembly_label` with basis `inferred`, the scores in "
            "one warning line. Rows that do not settle it refuse the facet with the scores. A collection whose "
            "contig 1 is 247,249,719 bases is GRCh36."
        ),
        status="handled",
        code=("assembly.py `infer_assembly_from_rows`, `MIN_SHARE`, `MIN_MARGIN`", "reference.py `_infer_from_rows`"),
        tests=("test_an_array_without_a_build_gets_it_from_its_rows",),
    ),
    Shape(
        key="gp-not-summing",
        title="A GP whose probabilities do not sum to 1",
        findings=("F10",),
        shape=(
            "FORMAT GP values all within [0, 1] that do not sum to 1 on some rows (phred GPs rounded to integers "
            "and converted back give such rows). From the numbers alone the scale is ambiguous."
        ),
        example=(
            "##fileformat=VCFv4.3\n"
            + _GT
            + '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities">\n'
            + _COLUMNS
            + "1\t752566\trs3094315\tG\tA\t.\tPASS\t.\tGT:GP\t0/1:0.03,0.9,0.05\n"
            + "1\t752721\trs3131972\tA\tG\t.\tPASS\t.\tGT:GP\t1/1:0,0.01,0.99\n"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "Values above 1 are phred and values summing to 1 are probabilities, whatever the header says. "
            "Values the numbers cannot settle are read as probabilities under a VCF 4.3 or later header, which "
            "defines GP so, and the rows that do not sum to 1 are counted in a warning. Under 4.2 or earlier, or "
            "with no version, GP is withheld. The FORMAT line's Description is free text and is not read."
        ),
        status="handled",
        code=("loaders/vcf.py `_gp`",),
        tests=("test_gp_scale_is_read_from_decisive_values", "test_gp_from_phred_survives_its_own_export"),
    ),
    Shape(
        key="dash-alt",
        title="A deletion written as ALT `-`",
        findings=("F12",),
        shape=(
            "Some VCF uploads write a deletion as REF the deleted bases and ALT `-`, the dash notation VCF does "
            "not have. The anchor base is missing; REF is the deleted span."
        ),
        example="##fileformat=VCFv4.2\n" + _GT + _COLUMNS + "1\t752721\t.\tAG\t-\t.\tPASS\t.\tGT\t0/0\n",
        expect=Expect("detected", "vcf"),
        handling=(
            "The loader reads `-` as the empty allele `\"\"`, the one spelling of an empty allele (CGI's native "
            "deletions use it too), and counts the records in one warning naming them unanchored deletions. "
            "`normalize`, once `resolve_ref` has matched REF, anchors it on the base before (`CT>-` at 30 "
            "becomes `GCT>G` at 29) and left-aligns it; at a contig's first base it stays as written, counted. "
            "`write_vcf` refuses any ALT that is no VCF allele as `alt_not_an_allele`."
        ),
        status="handled",
        code=("loaders/vcf.py `DASH_DELETION`", "normalize.py `_align`", "writing.py `_unwritable`"),
        tests=(
            "test_a_dash_alt_is_read_as_the_empty_allele_and_counted",
            "test_a_dash_deletion_is_anchored_on_the_base_before",
            "test_an_alt_that_is_no_vcf_allele_is_never_written",
        ),
    ),
    Shape(
        key="minus-strand-probes",
        title="An array export with part of its probes on the minus strand",
        findings=("F13",),
        shape=(
            "One FTDNA export writes about one probe in eight on the minus strand, with nothing in the row to "
            "say so. It shows only in aggregate: heterozygous rows carrying the complement of the reference "
            "base and not the base, which no forward-strand probe can write."
        ),
        example=(
            "RSID,CHROMOSOME,POSITION,RESULT\n"
            '"rs4477212","1","82154","TT"\n'
            '"rs3094315","1","752566","CT"\n'
            '"rs3131972","1","752721","AG"\n'
        ),
        example_note=(
            "On GRCh37 the reference base at 1:752566 is `G`; the het `CT` there carries its complement `C` and "
            "not `G`. The `AG` row carries its site's reference base `A`, so it reads the same on either strand; "
            "the homozygous `TT` row cannot show its strand at all."
        ),
        expect=Expect("detected", "ftdna"),
        handling=(
            "`resolve_ref` checks the strand of every array facet on a sample of its heterozygous two-allele "
            "rows. When more than 2% of at least 200 such hets carry the complement and not the base, the facet "
            "gets `strand_check = mixed`: every row that is not strand-safe (homozygous calls, half-calls, hets "
            "without the base) keeps `ref` and gets `alts` and `gt` null, each reason counted, and the safe hets "
            "keep their alleles. Complementing was rejected: an A/T or C/G probe is its own complement, the "
            "homozygous rows of a flipped probe cannot be found, and no public strand table names the flipped "
            "probes."
        ),
        status="handled",
        code=("strand.py `check_strand`, `withheld_reason`", "reference.py `_check_strand`"),
        tests=("test_a_mixed_strand_array_keeps_ref_and_withholds_the_alleles_it_cannot_read",),
    ),
    Shape(
        key="chain-blocks-twice",
        title="A liftover chain that lists every block twice",
        findings=("F14",),
        shape=(
            "Ensembl's `GRCh38_to_GRCh37` chain repeats every block of its main chains under a second chain id, "
            "with the same source span, target and strand. Counted per block, every lifted base looks as if it "
            "had two targets. Not an input genotype format, but a file the library reads."
        ),
        example=(
            "chain 1000 1 248956422 + 10000 10500 1 249250621 + 10000 10500 1\n"
            "500\n"
            "\n"
            "chain 1000 1 248956422 + 10000 10500 1 249250621 + 10000 10500 2\n"
            "500\n"
        ),
        example_note="UCSC chain text: a `chain` header line, then block sizes; the two chains differ only in id.",
        handling=(
            "`liftover` cuts the source axis where the set of targets changes, a target being a contig, strand "
            "and offset, so two blocks sending every shared base to the same place count once. A position under "
            "blocks with different targets is still `multiple_targets`."
        ),
        status="handled",
        code=("liftover.py `_segments`",),
        tests=("test_a_chain_listing_its_blocks_twice_lifts_as_once", "test_segments_count_targets_not_blocks"),
    ),
    Shape(
        key="damaged-bgzf",
        title="A BGZF file with a damaged block",
        findings=("F16",),
        shape=(
            "A BGZF file one of whose blocks does not inflate, or fails its CRC32 or its ISIZE, deep inside the "
            "file. The head reads, so the file is recognised; the damage shows only when the records are read."
        ),
        example_note="A byte-level fault, so there is no text example.",
        handling=(
            "When a read fails and the input's outer layer is BGZF, `load` walks the file block by block "
            "(header, `BC` size, inflate, CRC32, ISIZE) and raises `DamagedBgzfError` with the facet's id, the "
            "first bad block's compressed offset and what is wrong with it; it never names the path. A failure "
            "whose blocks all read is re-raised unchanged. Damage under another layer (a BGZF inside a zip, a "
            "plain gzip) keeps the reader's own error, since only an outer BGZF has block offsets to name."
        ),
        status="refused",
        code=("bgzf.py `find_damage`", "loading.py `DamagedBgzfError`"),
        tests=(
            "test_load_names_the_facet_and_the_block_never_the_path",
            "test_damage_in_plain_gzip_keeps_the_readers_error",
        ),
    ),
    # --- dialects from the contract gaps in FORMATS.md and the specs' notes ---------------------------
    Shape(
        key="cgi-format-versions",
        title="Complete Genomics columns that change with the format version",
        findings=(),
        shape=(
            "A CGI `var` file's columns depend on its `#FORMAT_VERSION`. 1.6 and 1.7 have one `totalScore` where "
            "2.x has `varScoreVAF` and `varScoreEAF`; 2.0 and 2.2 have `varQuality` (`VQHIGH`, `VQLOW`) where "
            "2.5 has `varFilter` (`PASS`, `VQLOW`, `AMBIGUOUS`, `;`-joined); `alleleFreq` and "
            "`alternativeCalls` arrive in 2.5. The name changed and so did the vocabulary."
        ),
        example=(
            "#FORMAT_VERSION\t2.0\n"
            "#GENERATED_BY\tcgatools\n"
            "#GENOME_REFERENCE\tNCBI build 37\n"
            "#SOFTWARE_VERSION\t2.0.0.26\n"
            "#TYPE\tVAR-ANNOTATION\n"
            "\n"
            ">locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\tvarScoreVAF\t"
            "varScoreEAF\tvarQuality\thapLink\txRef\n"
            + _cgi_row("1", "2", "all", "chr1", "0", "10000", "no-ref", "=", "?", width=14)
            + _cgi_row("3", "2", "1", "chr1", "752565", "752566", "snp", "G", "A", "120", "120", "VQHIGH", width=14)
            + _cgi_row("3", "2", "2", "chr1", "752565", "752566", "ref", "G", "G", "98", "98", "VQHIGH", width=14)
        ),
        expect=Expect("detected", "cgi_var", "GRCh37", check_build=True),
        handling=(
            "Detection does not depend on the version: `#TYPE` and the first nine names of the `>` line are the "
            "same. The loader reads the names from the `>` line and takes columns by name; `varQuality` is read "
            "as `varFilter`, its values kept as written in `sample_filters` (no CGI document read here equates "
            "`VQHIGH` with `PASS`). The files loaded so far are format 2.0 and 2.5; `totalScore` is not mapped "
            "onto either 2.x score."
        ),
        status="handled",
        code=("loaders/cgi.py `_FILTER_ALIASES`",),
        tests=("test_format_2_0_and_no_title",),
    ),
    Shape(
        key="cgi-unanchored",
        title="Complete Genomics loci: 0-based spans, insertions of zero length, partly called alleles",
        findings=(),
        shape=(
            "CGI writes 0-based half-open spans. An insertion has `begin == end` and an empty `reference`; a "
            "deletion has an empty `alleleSeq`; ref and no-call blocks write `reference` as `=`; `?` is an "
            "unknown sequence, and a `?` inside a sequence (`GC?`) marks uncalled bases in a partly called "
            "allele. A haplotype row of varType `no-call` beside a called one is a half call."
        ),
        example=(
            _CGI_VAR_25_HEAD
            + _cgi_row("5", "2", "1", "chr1", "752600", "752600", "ins", "", "T", "60", "60", "PASS", width=16)
            + _cgi_row("5", "2", "2", "chr1", "752600", "752600", "ref", "", "", "60", "60", "PASS", width=16)
            + _cgi_row("6", "2", "1", "chr1", "752700", "752702", "del", "AG", "", "55", "55", "PASS", width=16)
            + _cgi_row("6", "2", "2", "chr1", "752700", "752702", "no-call", "AG", "A?", width=16)
        ),
        expect=Expect("detected", "cgi_var", "GRCh37", check_build=True),
        handling=(
            "A span loads as `pos = begin + 1`, `end_pos = end`; an insertion point sits at `pos = begin` with "
            "`end_pos` null and `source_ref = \"\"`. Alleles stay unpadded, `\"\"` for an empty one, until "
            "`normalize` adds the anchor base from the reference. A haplotype written `?` is a null in "
            "`gt_alleles`, a `?` inside a sequence is kept verbatim, and every no-call varType (`no-call`, "
            "`no-call-rc`, `no-call-ri`, `no-ref`, `PAR-called-in-X`) reads as no call."
        ),
        status="handled",
        code=("loaders/cgi.py",),
        tests=("test_variant_loci", "test_unanchored_indels_count_so_a_genome_does_not_read_snv_only"),
    ),
    Shape(
        key="half-calls",
        title="Half-calls: one allele known, the other not",
        findings=(),
        shape=(
            "A genotype with one allele called and the other missing. VCF writes it `1|.` or `0/.` (Complete "
            "Genomics' VCF has thousands); FTDNA's Affymetrix-era exports write `-A`; Genes for Good writes `T-`; "
            "CGI's masterVar says zygosity `half`. A half-call is neither a no-call nor a homozygote."
        ),
        example=(
            "RSID,CHROMOSOME,POSITION,RESULT\n"
            "rs4477212 1 72017 AG\n"
            "rs3094315 1 742429 -A\n"
            "rs3131972 1 742584 ---\n"
            "rs12562034 1 758311 GG\n"
        ),
        expect=Expect("detected", "ftdna_space"),
        handling=(
            "The called allele stays in its slot and the other is null: `-A` loads as `[null, \"A\"]`. A whole "
            "no-call spelling (`spec.no_call`) is one token: in formats that write two letters on every called "
            "row, FTDNA's `---` or a Final Report's single `-` gives two nulls. VCF half-calls keep their phase "
            "and are graded evidence."
        ),
        status="handled",
        code=("loaders/text.py", "loaders/table.py"),
        tests=("test_ftdna_half_calls_no_calls_and_hemizygous_letters", "test_half_calls_phase_and_filters"),
    ),
    Shape(
        key="genos-columns-meta",
        title="Column names in a meta line, 0-based positions and the empty allele `-`",
        findings=(),
        shape=(
            "Genos names its columns in a `##Columns=` meta line instead of a column line, writes 0-based "
            "half-open spans, and writes the empty allele as `-`: at an insertion site (`begin == end`) `-/-` is "
            "a homozygous call for no insertion, not a no-call. A single allele with no `/` is a hemizygous call "
            "on X or Y."
        ),
        example=(
            "##FileFormat=Genos\n"
            "##ReferenceGenome=UCSC.hg19\n"
            "##Columns=CHROM\tBEGIN\tEND\tID\tGENOTYPE\n"
            "chr1\t752565\t752566\trs3094315\tA/G\n"
            "chr1\t752720\t752721\trs3131972\tG/G\n"
            "chr1\t846808\t846808\trs4475691\t-/-\n"
        ),
        expect=Expect("detected", "genos", "GRCh37", check_build=True),
        handling=(
            "The loader reads the column names from `##Columns=` when the line is there (the spec's aliases "
            "match them) and falls back to the spec's order when it is not. A span moves to 1-based; an "
            "insertion site goes to the base before it, where VCF places it, with `end_pos` null until the "
            "anchor base is read. `-` loads as the empty allele `\"\"`, and `source_gt` keeps the `-`."
        ),
        status="handled",
        code=("loaders/table.py",),
        tests=(
            "test_genos_spans_insertions_and_hemizygous_calls",
            "test_genos_reads_every_field_with_or_without_its_column_line",
        ),
    ),
    Shape(
        key="joined-exports-and-mixed-delimiters",
        title="A column line repeated inside the file, and a delimiter that changes",
        findings=(),
        shape=(
            "Some FTDNA uploads are two exports joined, the column line repeated where the second starts. The "
            "2010-era Affymetrix exports have a comma header over space-separated rows and end in a few hundred "
            "quoted comma rows. Rows are double-quoted in some exports and bare in others."
        ),
        example=(
            "RSID,CHROMOSOME,POSITION,RESULT\n"
            '"rs4477212","1","72017","AA"\n'
            '"rs3094315","1","742429","AG"\n'
            "RSID,CHROMOSOME,POSITION,RESULT\n"
            '"rs3131972","1","742584","GG"\n'
        ),
        expect=Expect("detected", "ftdna"),
        handling=(
            "The repeated column line is skipped and counted. Every field is stripped of surrounding double "
            "quotes. An `ftdna_space` row is read by its own shape, so the quoted comma rows at its end load too."
        ),
        status="handled",
        code=("loaders/table.py",),
        tests=(
            "test_ftdna_files_joined_with_their_column_line_repeated",
            "test_ftdna_space_export_with_a_quoted_comma_tail",
        ),
    ),
    Shape(
        key="build-in-the-rows",
        title="0-based single positions and a build stated on every row",
        findings=(),
        shape=(
            "Color's genotype CSV writes each position one below its GRCh37 1-based position, states the build "
            "in a `ReferenceVersion` column on every row instead of in a header, and ends in a disclaimer "
            "trailer whose quoted paragraph holds commas."
        ),
        example=(
            "SampleID,Chromosome,Position,RSID,Genotype,ReferenceVersion\n"
            "S1,2,136608645,rs4988235,AG,GRCh37\n"
            "S1,12,112241765,rs671,GG,GRCh37\n"
            "Raw data disclaimer:,,,,,\n"
            '"This paragraph, like the real one, holds commas."\n'
        ),
        expect=Expect("detected", "color_genotypes", "GRCh37", check_build=True),
        handling=(
            "The loader adds 1 for every 0-based spec. The facet's build is the one every row states; rows that "
            "disagree leave it unknown, with a warning. Reading stops at the first row whose Position is empty."
        ),
        status="handled",
        code=("loaders/table.py",),
        tests=(
            "test_color_positions_are_shifted_to_one_based_and_the_trailer_is_cut",
            "test_color_build_is_read_from_the_rows",
        ),
    ),
    Shape(
        key="build-in-a-sibling",
        title="A build stated in another member of the same zip",
        findings=(),
        shape=(
            "A deCODEme raw-data zip holds the genotype CSV and a small key,value member that states the build; "
            "the genotype member itself names none. Rows on the `-` strand are written on that strand."
        ),
        members=(
            (
                "genotypes.csv",
                "Name,Variation,Chromosome,Position,Strand,YourCode\n"
                "rs4477212,A/G,1,72017,+,AG\n"
                "rs3094315,C/T,1,742429,-,CT\n"
                "rs3131972,C/T,1,742584,-,TT\n"
                "rs12562034,A/G,1,758311,+,GG\n",
            ),
            (
                "info.txt",
                "Format version,1\n"
                "SNP count,4\n"
                "Reference sequence,NCBI Build 36 and Cambridge reference sequence\n",
            ),
        ),
        example_note="Two members of one zip; member names are made up.",
        expect=Expect("detected", "decodeme", None, check_build=True, layers=("zip",)),
        handling=(
            "`sniff` answers for the largest detected member (`decodeme`, build unknown) and carries every "
            "member's result in `members`; the info member is `decodeme_info` with build GRCh36. The loader "
            "gives the deCODEme facet the build its zip's `decodeme_info` member states, and complements both "
            "letters on `-` rows to get the plus strand."
        ),
        status="handled",
        code=("loading.py", "loaders/table.py"),
        tests=(
            "test_decodeme_build_comes_from_the_info_member_of_its_zip",
            "test_decodeme_minus_strand_rows_are_complemented",
        ),
    ),
    Shape(
        key="section-header",
        title="Column names after a section line, tab or comma by export",
        findings=(),
        shape=(
            "An Illumina GenomeStudio Final Report opens with a `[Header]` block of key/value lines, then a "
            "`[Data]` line, then the column line. The delimiter is tab in `.txt` exports and comma in `.csv` "
            "ones, where the header block is padded with commas. Which allele columns appear depends on the "
            "export (Plus, Forward, Design or Top strand)."
        ),
        example=(
            "[Header],,,\n"
            "Content,GSA-24v2-0_A1.bpm,,\n"
            "Num SNPs,665608,,\n"
            "[Data],,,\n"
            "SNP Name,Chr,MapInfo,Alleles - Plus\n"
            "rs3094315,1,752566,AG\n"
            "rs3131972,1,752721,GG\n"
            "rs12124819,1,776546,-\n"
        ),
        expect=Expect("detected", "illumina_final_report"),
        handling=(
            "The column names are read from the line after `[Data]`, and the delimiter off that line. The "
            "loader needs `Chr`, `MapInfo` and `Alleles - Plus`; an export with only Forward, Design or Top "
            "alleles is refused with `TableLayoutError`, since turning those into plus-strand letters needs the "
            "manifest. A report of several samples loads one at a time (`load(path, sample=…)`). A report with "
            "no column line and no positions is `illumina_final_report_headerless`, which does not load."
        ),
        status="handled",
        code=("loaders/table.py",),
        tests=(
            "test_final_report_tab_and_comma_exports_load_alike",
            "test_final_report_without_plus_strand_alleles_is_refused",
            "test_final_report_with_two_samples_and_unplaced_probes",
        ),
    ),
    Shape(
        key="optional-columns",
        title="A column that only some versions write",
        findings=(),
        shape=(
            "Promethease's `tab-genotype` report has a `Repute` column from some version on: 0.1.110 does not "
            "write it, 0.1.161 does. A column declared for every version is then missing from older files."
        ),
        example=(
            "#Promethease\ttab-genotype\t0.1.110\t2013-01-01 00:00\n"
            "#Column-names\t\n"
            "Name\tMagnitude\tPopulationFreq\tDescription\n"
            "rs1333049(C;G)\tNone\tNone\t\n"
            "rs10490924(G;G)\t0.0\tNone\tnormal\n"
        ),
        expect=Expect("detected", "promethease_tab_genotype"),
        handling=(
            "Detected as `promethease_tab_genotype` either way. A `ColumnSpec` cannot say a column appears only "
            "in some versions, so the spec lists `repute` and its notes say to take names from the column line. "
            "The format does not load: it is an annotation report without chromosome or position."
        ),
        status="open",
        code=("formats/spec.py `ColumnSpec`",),
    ),
    Shape(
        key="chromosome-outside-the-row",
        title="A chromosome that is not in the rows",
        findings=(),
        shape=(
            "Gencove's Oxford `.gen` imputation output ships one member per chromosome and repeats the "
            "chromosome only on typed rows; imputed rows write `---` in its place. FTDNA Big Y files have no "
            "chromosome column at all: every row is chrY."
        ),
        example=(
            "--- rs4477212 82154 A G 0.9 0.1 0\n"
            "--- rs3094315 752566 G A 0.2 0.7 0.1\n"
            "1 rs3131972 752721 A G 0 1 0\n"
        ),
        expect=Expect("detected", "oxford_gen"),
        handling=(
            "Detected, not loaded. Taking the chromosome from a member's name would trust a name, which `sniff` "
            "never does; the loader for these formats is not built."
        ),
        status="open",
        code=("formats/other_tables.py `OXFORD_GEN`, `FTDNA_BIGY`",),
    ),
    Shape(
        key="headerless-bed",
        title="A headerless BED file",
        findings=(),
        shape="Three tab-separated columns, chromosome, 0-based start and exclusive end, with no header.",
        example="chr1\t752565\t752566\nchr1\t752720\t752721\n",
        expect=Expect("unknown"),
        handling=(
            "`unknown`, not `not_variant_data`: the text kinds look at one line, and telling a BED from a "
            "headerless table needs a check over several lines that `sniff` does not make."
        ),
        status="open",
        code=("sniff.py `_TEXT_KINDS`",),
    ),
    Shape(
        key="array-content-in-a-vcf",
        title="Array or imputed genotypes inside a VCF",
        findings=(),
        shape=(
            "Illumina GenotypingToVCF, 23andMe conversions, Genes for Good and imputation servers write VCF too. "
            "The parse is a VCF's, but an absent site means unknown, never reference. GenotypingToVCF puts "
            "probes it cannot place at chromosome `NA`, POS 0."
        ),
        example=(
            "##fileformat=VCFv4.1\n"
            "##workflow_type=Illumina_GenotypingToVCF\n"
            + _GT
            + _COLUMNS
            + "1\t752566\trs3094315\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
            + "NA\t0\trs3131972\tN\tA,G\t.\tGTEX\t.\tGT\t./.\n"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "`sniff` says `vcf`; `Family` has no axis for array content. The loader's source detector reads the "
            "full header and sets the facet's `callset_type` to `array` or `imputed` where the header says so "
            "(here `##workflow_type`). polars-bio cannot read a record at POS 0, so those records are left out "
            "of the copy it reads and counted."
        ),
        status="handled",
        code=("source.py `detect_source`", "loaders/vcf.py (POS 0)"),
        tests=("test_load_states_array_content_and_warns_on_a_conflict", "test_records_at_pos_zero_are_left_out_and_counted"),
    ),
    Shape(
        key="variable-width-rows",
        title="Rows of any width after fixed fields",
        findings=(),
        shape=(
            "PGP-HMS's key=value export of a freeBayes+VarScan2 merge has twelve fixed fields, eight of them "
            "`KEY=value`, then an annotation tail of any length (23 to 84 fields per row in the files seen)."
        ),
        example=(
            "chr1\t752566\tG\tA\tFBGenoType=0/1\tFBRefAlleleCount=10\tFBReferenceAlleleQ=400\tFBTotalDepth=25"
            "\tVSGenoType=0/1\tVSRefAlleleCount=10\tVSReferenceAlleleQ=35\tVSTotalDepth=24"
            "\tEFF_RESULTS\tMODIFIER\tKG=0.5\n"
            "chr1\t752721\tA\tG\tFBGenoType=1/1\tFBRefAlleleCount=.\tFBReferenceAlleleQ=NULL\tFBTotalDepth=20"
            "\tVSGenoType=NULL\tVSRefAlleleCount=NULL\tVSReferenceAlleleQ=NULL\tVSTotalDepth=NULL\n"
        ),
        expect=Expect("detected", "freebayes_varscan2_keyvalue"),
        handling=(
            "Only the twelve fixed fields are read, each checked for its key; rows whose fields do not carry the "
            "keys in order are counted. The tail's keys, read from the first 1,000 rows, are listed in "
            "`facets.dropped_fields`. One facet per caller."
        ),
        status="handled",
        code=("loaders/keyvalue.py",),
        tests=("test_one_facet_per_caller",),
    ),
    Shape(
        key="fileformat-not-first",
        title="`##fileformat` not on the first line, and keys the header never declares",
        findings=(),
        shape=(
            "Genos exome VCFs sort their meta lines, so `##fileformat` can sit hundreds of lines down. Other "
            "files use INFO or FORMAT keys, GT included, that no meta line declares. polars-bio refuses the "
            "first and silently drops the undeclared keys of the second."
        ),
        example=(
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            "##fileformat=VCFv4.2\n"
            "##reference=hg19\n"
            + _COLUMNS
            + "chr1\t752566\trs3094315\tG\tA\t.\tPASS\tDP=30\tGT\t0/1\n"
        ),
        example_note="GT is used by the row and declared nowhere.",
        expect=Expect("detected", "vcf", "GRCh37", check_build=True),
        handling=(
            "Detected as `vcf` (a `##fileformat` line anywhere in the head). The loader moves `##fileformat` to "
            "line 1 of the copy it reads and warns that the header was reordered; keys the rows use without a "
            "declaration are declared by the loader, named in a warning."
        ),
        status="handled",
        code=("loaders/vcf.py `undeclared_keys`, `declare`",),
        tests=("test_header_repairs", "test_an_undeclared_end_is_declared_so_a_block_keeps_its_extent"),
    ),
    Shape(
        key="vcf-without-meta",
        title="A VCF body under a bare `#CHROM` line",
        findings=(),
        shape="A VCF with no `##` meta lines at all, the column line first.",
        example=(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
            "chr1\t752566\t.\tG\tA\t50.2\t.\tDP=12\tGT:DP\t0/1:12\n"
        ),
        expect=Expect("detected", "vcf_no_meta"),
        handling=(
            "Its own format id, since the parse differs: the loader reads it under a synthesised header "
            "(`##fileformat=VCFv4.2` and a definition for each carried key the rows use), lists every other key "
            "in `dropped_fields`, leaves `vcf_version` null and warns that the file states no version, contigs "
            "or build."
        ),
        status="handled",
        code=("loading.py", "loaders/vcf.py"),
        tests=("test_vcf_no_meta_loads_under_a_synthesised_header",),
    ),
    Shape(
        key="declared-two-carries-one",
        title="A header that declares more samples than the rows carry",
        findings=(),
        shape=(
            "PGP-HMS's freeBayes+VarScan2 merges declare two sample columns (`unknown`, `Sample1`) while most "
            "rows carry one. Which caller the one column belongs to is said by INFO/SF, not by its name."
        ),
        example=(
            "##fileformat=VCFv4.1\n"
            "##source=freeBayes v0.9.21\n"
            "##reference=/ref/hg19.fa\n"
            "##source=VarScan2\n"
            '##INFO=<ID=SF,Number=.,Type=String,Description="Source File">\n'
            + _GT
            + "##source_20160919.1=vcf-merge chr1.FB.vcf.gz chr1.VS.vcf.gz\n"
            "##sourceFiles_20160919.1=0:chr1.FB.vcf.gz,1:chr1.VS.vcf.gz\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\tSample1\n"
            "chr1\t752566\t.\tG\tA\t40.1\tPASS\tSF=0,1\tGT\t0/1\n"
            "chr1\t776546\t.\tA\tG\t.\tPASS\tSF=1\tGT\t0/1\n"
        ),
        expect=Expect("detected", "vcf_freebayes_varscan2", "GRCh37", check_build=True),
        handling=(
            "polars-bio fills the missing second sample with None, which is a missing column, not a no-call. "
            "The loader returns one facet per caller (`relation = caller_split`), each caller identified by the "
            "`.FB.`/`.VS.` names on `##sourceFiles`; in the one-column layout the VarScan2 facet holds only the "
            "sites VarScan2 alone found, and a warning counts the shared ones."
        ),
        status="handled",
        code=("loading.py", "loaders/vcf.py"),
        tests=("test_two_declared_columns_one_carried",),
    ),
    Shape(
        key="gvcf-blocks",
        title="Reference blocks in a VCF, and headers that only look like a gVCF",
        findings=(),
        shape=(
            "gVCF blocks are spelled ALT `<NON_REF>` (GATK), `<*>` (DeepVariant), or ALT `.` with INFO END "
            "(Illumina, Genos). `##ALT=<ID=NON_REF>` and `##GVCFBlock` survive into genotyped files with no "
            "blocks at all, and END also sits on SV records."
        ),
        example=(
            "##fileformat=VCFv4.2\n"
            '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">\n'
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
            + _GT
            + _COLUMNS
            + "1\t10001\t.\tT\t<*>\t0\t.\tEND=10468\tGT\t0/0\n"
            + "1\t752566\trs3094315\tG\tA,<*>\t50\tPASS\t.\tGT\t0/1\n"
        ),
        expect=Expect("detected", "vcf"),
        handling=(
            "gVCF is not a format id: the header marks lie both ways. Whether a callset has blocks is read from "
            "its unfiltered rows, and the blocks go to the `blocks` frame, where they decide what an absent "
            "site means (`absence_means = blocks`)."
        ),
        status="handled",
        code=("loaders/vcf.py", "callset.py"),
        tests=(
            "test_blocks_are_all_star_records_asserting_reference",
            "test_illumina_gvcf_blocks",
            "test_load_keeps_blocks_deciding_absence_on_a_typed_gvcf",
        ),
    ),
    Shape(
        key="converter-column-lines",
        title="A 23andMe-layout table written by a converter, with two column lines",
        findings=(),
        shape=(
            "Converters emulate the 23andMe layout without its preamble. One opens `# This file emulates the "
            "format of 23andMe` and carries two column lines, the first with an empty genotype name; others "
            "put a blank line before the column line."
        ),
        example=(
            "# This file emulates the format of 23andMe\n"
            "# rsid\tchromosome\tposition\t\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs4477212\t1\t82154\tAG\n"
            "rs3094315\t1\t752566\tGG\n"
        ),
        expect=Expect("detected", "rsid_table", None, check_build=True),
        handling=(
            "Detected as `rsid_table`, build unknown. The loader checks the last column line against the four "
            "names in order; every comment line is skipped."
        ),
        status="handled",
        code=("loaders/rsid_table.py `column_line_ok`",),
        tests=("test_preamble_quirks_and_line_endings",),
    ),
    Shape(
        key="stray-byte",
        title="A stray byte before the first `#`",
        findings=(),
        shape="One Genes for Good export opens `I# Genes for Good v1.0: …`: one byte before the comment mark.",
        example=(
            "I# Genes for Good v1.0: filtered, unphased genotype calls.\n"
            "# This data file generated by PLINK on: 01/01/16\n"
            "# Human genome reference used: GRCh37/Mito:rCRS\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs4477212\t1\t82154\tAG\n"
            "rs3094315\t1\t752566\tT-\n"
        ),
        expect=Expect("detected", "genes_for_good", "GRCh37", check_build=True),
        handling=(
            "The detection pattern allows one byte before the `#` of line 1. The loader's comment rule does the "
            "same (a line starting with `#`, `\"#`, or one byte and then `#`), so the line is skipped."
        ),
        status="handled",
        code=("loaders/rsid_table.py `_COMMENT`",),
        tests=("test_preamble_quirks_and_line_endings",),
    ),
    Shape(
        key="spreadsheet-23andme",
        title="A 23andMe export re-saved by a spreadsheet as CSV, and cut short",
        findings=(),
        shape=(
            "A spreadsheet re-save pads every preamble line to four fields with `,,,` and quotes preamble lines "
            "that hold a comma, so they start with `\"#`. The one real file stops mid chromosome 2 after about "
            "65,500 lines: the spreadsheet's row limit cut it."
        ),
        example=(
            "# This data file generated by 23andMe at: Fri Jan  1 00:00:00 2016,,,\n"
            "#,,,\n"
            '"# More information on reference human assembly build 37 (a.k.a. Annotation Release 104), see:",,,\n'
            "# We are using reference human assembly build 37 (also known as Annotation Release 104).,,,\n"
            "# rsid,chromosome,position,genotype\n"
            "rs4477212,1,82154,AG\n"
            "rs3094315,1,752566,GG\n"
        ),
        expect=Expect("detected", "23andme_csv", "GRCh37", check_build=True),
        handling=(
            "Its own format id, since the delimiter changed. Padded and quoted preamble lines are skipped as "
            "comments. The cut is not detected: the loader reads the rows that are there, and sites past the "
            "cut are absent, which in an array table means unknown, never reference."
        ),
        status="open",
        code=("formats/rsid_tables.py `TWENTYTHREE_AND_ME_CSV`",),
        tests=("test_a_resaved_23andme_table_is_read_by_its_internal_ids",),
    ),
    Shape(
        key="encodings-and-line-ends",
        title="Byte-order marks, CRLF, bare CR and mixed line ends",
        findings=(),
        shape=(
            "Exports arrive in UTF-16 with a BOM, UTF-8 with a BOM, latin-1, with CRLF or bare-CR line ends, "
            "and with line ends mixed inside one file (AncestryDNA writes some preamble lines LF and the rest "
            "CRLF)."
        ),
        example=(
            "#AncestryDNA raw data download\n"
            "#This file was generated by AncestryDNA at: 01/01/2016 00:00:00 UTC\r\n"
            "#Data was collected using AncestryDNA array version: V1.0\r\n"
            "#of the SNP using human reference build 37.1 coordinates.\r\n"
            "rsid\tchromosome\tposition\tallele1\tallele2\r\n"
            "rs4477212\t1\t82154\tA\tA\r\n"
            "rs3094315\t1\t752566\tA\tG\r\n"
        ),
        example_note="Line 1 ends in LF, the others in CRLF; shown here as plain line breaks.",
        expect=Expect("detected", "ancestrydna", "GRCh37", check_build=True),
        handling=(
            "`sniff` decodes BOM-aware (UTF-16, UTF-8), falls back to latin-1, and turns CRLF and bare CR into "
            "LF before matching; the table loaders read lines the same way, so each line's own ending is "
            "stripped rather than one terminator chosen per file."
        ),
        status="handled",
        code=("sniff.py `_decode`", "loaders/text.py `lines`"),
        tests=(
            "test_utf16_crlf_and_bom_exports_are_read",
            "test_cr_only_line_endings_split_into_lines",
            "test_ancestrydna_mixed_line_endings_and_bad_rows",
        ),
    ),
)

NOT_SHAPES: tuple[tuple[str, str], ...] = (
    ("F5", "Memory of the canonical sort on large VCFs: a cost of the load, not a file shape."),
    ("F11", "PL converted to `gl` one float step off: arithmetic in the loader, not a file shape."),
)
"""Findings in the F1-F17 range that describe no input shape, so a reader checking the range finds them."""


def materialize(shape: Shape, directory: Path) -> Path:
    """Write `shape`'s example as the file it describes, with its layers, and return its path."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / shape.key
    if shape.members:
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            for name, text in shape.members:
                zf.writestr(name, text)
        return path
    if shape.example is None:
        raise ValueError(f"{shape.key} has no example to write")
    data = shape.example.encode()
    for layer in shape.layers:
        if layer == "bz2":
            data = bz2.compress(data)
        elif layer == "zip_deflate64":
            data = _deflate64_zip(path, data)
        else:
            raise ValueError(f"unknown layer {layer!r}")
    path.write_bytes(data)
    return path


def _deflate64_zip(path: Path, payload: bytes) -> bytes:
    """A stored zip whose local and central headers both say method 9, as a Deflate64 archive does."""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
        zf.writestr("calls.vcf", payload)
    data = bytearray(path.read_bytes())
    local, central = data.index(b"PK\x03\x04"), data.index(b"PK\x01\x02")
    data[local + 8 : local + 10] = struct.pack("<H", 9)
    data[central + 10 : central + 12] = struct.pack("<H", 9)
    return bytes(data)


# ---------------------------------------------------------------------------------------------------
# Rendering


def _code(text: str) -> str:
    """`text` as inline code, with a fence longer than any backtick run inside it."""
    longest = max((len(run) for run in re.findall(r"`+", text)), default=0)
    fence = "`" * (longest + 1)
    pad = " " if text.startswith("`") or text.endswith("`") or longest else ""
    return f"{fence}{pad}{text}{pad}{fence}"


def _cell(text: str) -> str:
    """Text safe inside a table cell."""
    return text.replace("|", "\\|").replace("\n", " ")


def _block(text: str) -> list[str]:
    """`text` as a fenced block. Line ends are shown as line breaks; tabs stay literal tabs."""
    shown = text.replace("\r\n", "\n").replace("\r", "\n")
    longest = max((len(run) for run in re.findall(r"`{3,}", shown)), default=2)
    fence = "`" * (longest + 1)
    return [f"{fence}text", *shown.rstrip("\n").split("\n"), fence]


def _anchor(text: str) -> str:
    """The GitHub heading anchor of a heading that is a bare id in backticks."""
    return re.sub(r"[^a-z0-9_ -]", "", text.lower()).replace(" ", "-")


def dtype_name(dtype: object) -> str:
    """A language-neutral name for a column dtype."""
    if isinstance(dtype, pl.List):
        return f"list of {dtype_name(dtype.inner)}"
    for cls, name in ((pl.String, "string"), (pl.Int64, "integer"), (pl.Float64, "float"), (pl.Date, "date")):
        if dtype == cls:
            return name
    raise ValueError(f"no neutral name for dtype {dtype!r}; add one to scripts/format_registry.py")


def _module_of(spec: FormatSpec) -> str:
    for module in MODULES:
        if spec in (*module.SPECS, *module.FALLBACKS):
            return module.__name__
    raise ValueError(f"{spec.id} is in no family module")


def _is_fallback(spec: FormatSpec) -> bool:
    return any(spec in module.FALLBACKS for module in MODULES)


def _loader_of(format_id: str) -> str | None:
    names = [name for name, formats in LOADERS if format_id in formats]
    if len(names) > 1:
        raise ValueError(f"{format_id} is claimed by several loaders: {', '.join(names)}")
    return names[0] if names else None


def _delimiter(spec: FormatSpec) -> str:
    if spec.id == BCF_FORMAT:
        return "none: binary records (BCF 2), read through polars-bio"
    if spec.header is HeaderMode.VCF:
        return "VCF: tab-separated, read through polars-bio"
    if spec.delimiter is None:
        if spec.header is HeaderMode.SECTION:
            return "tab or comma, read off the column line"
        return "not delimited text"
    return {"\t": "tab", ",": "comma", " ": "single space"}.get(spec.delimiter, _code(spec.delimiter))


def _all_regexes() -> list[str]:
    found: list[str] = []
    for spec in REGISTRY:
        found += [*spec.detection.patterns, *spec.detection.excludes, *filter(None, [spec.detection.data_shape])]
        found += [hint.pattern for hint in spec.build_hints]
    return found


def _render_algorithm() -> list[str]:
    binary = ", ".join(f"{_code(_bytes(magic))} {kind}" for magic, kind in _BINARY_KINDS)
    text_kinds = "\n".join(f"   - {_code(regex.pattern)}: `{kind}`" for regex, kind in _TEXT_KINDS)
    methods = ", ".join(f"{number} ({name.split(' (')[0]})" for number, name in _ZIP_METHODS.items())
    regexes = _all_regexes()
    constructs = [
        f"- {_code(token)}: {meaning}"
        for token, pattern, meaning in REGEX_CONSTRUCTS
        if any(re.search(pattern, regex) for regex in regexes)
    ]
    return [
        "## How a file is recognised",
        "",
        "Content decides, never a name, extension or upload label. The steps, in order:",
        "",
        "1. **Peel layers by magic bytes**, reading the first 512 bytes of each layer:",
        "   - `1f 8b`: gzip; it is **BGZF** when the FLG byte (offset 3) has FEXTRA (`0x04`) set and bytes "
        "12-13 are `BC`.",
        "   - `BZh`: bz2. `fd 37 7a 58 5a 00`: xz. `28 b5 2f fd`: zstd, reported `unreadable`.",
        "   - `PK 03 04` or `PK 05 06`: zip. `ustar` at offset 257: tar.",
        f"   - At most {MAX_DEPTH} nested layers are followed. A zip inside another layer is buffered, up to "
        f"{NESTED_ZIP_CAP >> 20} MiB; a tar stream is listed for at most {TAR_SCAN_CAP >> 30} GiB.",
        "2. **Archives.** A zip holding `[Content_Types].xml` is an office document (`ooxml`), one holding "
        "`mimetype` is `odf`, one holding `index.xml` or `Index/*.iwa` is `iwork`; all three are "
        "`not_variant_data`. Members matching "
        f"{_code(_ZIP_JUNK.pattern)} are skipped. The rest are sniffed largest first, ties by name, at most "
        f"{MAX_MEMBERS}; the answer leads with the largest detected member and lists every member's answer. "
        f"Zip methods the reader cannot inflate are reported `unreadable` and named: {methods}.",
        "3. **Binary payloads.** On the innermost layer, `42 43 46 02` (`BCF\\2`) is the format `bcf`. These "
        f"prefixes are `not_variant_data` of the named kind: {binary}.",
        f"4. **Read a head**: whole lines from the start, stopping once {HEAD_DATA_LINES} non-empty lines not "
        f"starting with `#` are complete, or at {HEAD_BYTE_CAP >> 20} MiB, cut back to the last line end "
        "unless the file ended. An empty or all-whitespace head is `unknown`.",
        "5. **Decode**: a UTF-16 BOM (`ff fe` or `fe ff`) means UTF-16. Otherwise a NUL byte in the first "
        "8 KiB makes the payload binary, which is `unknown`. Otherwise UTF-8 with any BOM dropped (a sequence "
        "cut off in the last 4 bytes is trimmed); if that fails, latin-1.",
        "6. **Normalise line ends**: CRLF, then bare CR, become LF. Every pattern below sees LF only.",
        "7. **Match the registry in the order of the table below.** A spec matches when:",
        "   - no `excludes` pattern is found anywhere in the head;",
        "   - every `patterns` regex is found in the head, searched with `^` and `$` as line anchors "
        "(multiline) and `\\A` as the start of the head;",
        f"   - and, when it has a `data_shape`, one of the first {DATA_SHAPE_LINES} data lines matches it in "
        "full. Data lines are the non-blank lines not starting with the spec's comment prefix.",
        "",
        "   The first spec that matches is the answer; later ones that also match are reported as "
        "`also_matched`. A signature (a vendor's own text) is always tried before a fallback (data shape "
        "alone).",
        "8. **Text that holds no calls.** If no spec matches, the first non-blank line is tried against these "
        "patterns, each naming a `not_variant_data` kind:",
        text_kinds,
        "",
        "   Anything else is `unknown`.",
        "9. **Build**: every build hint of the winning spec is searched in the head (multiline). One build "
        "named: that build, with the matched text as evidence. Several builds named: no build, and a warning. "
        "None: the build is unknown. No vendor default is ever assumed; inferring a build from contig lengths "
        "or rows is the loader's work, not the sniffer's.",
        "",
        "The four outcomes are different negatives and must stay apart: `detected`, `unknown` (readable, "
        "nothing matched, empty or binary), `unreadable` (a layer could not be opened, with a reason) and "
        "`not_variant_data` (recognised, holds no calls, with a kind).",
        "",
        "### Regular expressions",
        "",
        "Patterns are written in Python's `re` dialect. Besides POSIX extended syntax they use:",
        "",
        *constructs,
        "",
    ]


def _bytes(magic: bytes) -> str:
    return "".join(chr(b) if 0x20 <= b < 0x7F and chr(b) != "`" else f"\\x{b:02x}" for b in magic)


def _render_vocabulary() -> list[str]:
    lines = ["## Vocabulary", "", "**Family**, what a format's content is (a closed set):", ""]
    lines += [f"- `{family.value}`: {FAMILY_MEANING[family]}." for family in Family]
    lines += ["", "**Coordinates**, the convention of a format's position columns:", ""]
    lines += [f"- `{c.value}`: {COORDINATES_MEANING[c]}." for c in Coordinates]
    lines += ["", "**Column names**, where a delimited format keeps them:", ""]
    lines += [f"- `{mode.value}`: {HEADER_MEANING[mode]}." for mode in HeaderMode]
    lines += [
        "",
        "Format ids are lowercase snake_case and form an open set: ids are added, never renamed or reused.",
        "",
    ]
    return lines


def _render_index() -> list[str]:
    lines = [
        "## Detection order",
        "",
        f"{len(REGISTRY)} formats, {len(LOADABLE_FORMATS)} of them loadable. Detection tries them in this order.",
        "",
        "| # | Format | Family | Vendor | Kind | Loads |",
        "|---|---|---|---|---|---|",
    ]
    for position, spec in enumerate(REGISTRY, start=1):
        kind = "fallback" if _is_fallback(spec) else "signature"
        loads = "yes" if spec.id in LOADABLE_FORMATS else "no"
        lines.append(
            f"| {position} | [`{spec.id}`](#{_anchor(spec.id)}) | `{spec.family.value}` | {_cell(spec.vendor)} | "
            f"{kind} | {loads} |"
        )
    return [*lines, ""]


def _render_spec(position: int, spec: FormatSpec) -> list[str]:
    loader = _loader_of(spec.id)
    kind = "a fallback, tried after every signature" if _is_fallback(spec) else "a signature"
    if spec.id == BCF_FORMAT:
        kind = "a signature, though its magic bytes are checked before any text spec"
    lines = [f"### `{spec.id}`", "", spec.description, ""]
    lines.append(f"- **Family:** `{spec.family.value}`")
    lines.append(f"- **Vendor:** {spec.vendor}")
    lines.append(f"- **Registry position:** {position} of {len(REGISTRY)}, {kind} (`{_module_of(spec)}`)")
    lines.append(f"- **Loads:** yes, through `{loader}`" if loader else "- **Loads:** no")
    lines.append(f"- **Delimiter:** {_delimiter(spec)}")
    prefix = _code(spec.comment_prefix) if spec.comment_prefix else "none"
    lines.append(f"- **Comment prefix:** {prefix}")
    if spec.header is not HeaderMode.VCF and spec.columns:
        names = HEADER_MEANING[spec.header]
        if spec.header is HeaderMode.SECTION:
            names = f"the first line after `[{spec.section}]`"
        if spec.header_prefix:
            names += f"; the marker {_code(spec.header_prefix)} before the first name is not part of it"
        lines.append(f"- **Column names:** {names}")
    if spec.coordinates is None:
        lines.append("- **Coordinates:** none (the rows are not placed on a genome)")
    else:
        lines.append(f"- **Coordinates:** {spec.coordinates.value} ({COORDINATES_MEANING[spec.coordinates]})")
    if spec.genotype_encoding:
        lines.append(f"- **Genotype:** {spec.genotype_encoding}")
    if spec.no_call:
        lines.append(f"- **No-call spellings:** {', '.join(_code(token) for token in spec.no_call)}")
    elif spec.genotype_encoding:
        lines.append("- **No-call spellings:** none declared")
    if spec.build_hints:
        lines.append("- **Build:** named only by a header statement matching one of:")
        lines += [f"  - {_code(hint.pattern)} names {hint.build}" for hint in spec.build_hints]
    else:
        lines.append("- **Build:** no header statement is read; the build is unknown until inferred")
    lines += ["", "**Recognised by**", ""]
    if spec.id == BCF_FORMAT:
        lines.append(
            "- binary magic: after the compression layers are peeled, the payload starts with `42 43 46 02` "
            f"(`BCF\\2`). Its text detection {_code(spec.detection.patterns[0])} never matches."
        )
    else:
        if spec.detection.patterns:
            lines.append("- every one of these is found in the head:")
            lines += [f"  - {_code(pattern)}" for pattern in spec.detection.patterns]
        if spec.detection.excludes:
            lines.append("- none of these is found:")
            lines += [f"  - {_code(pattern)}" for pattern in spec.detection.excludes]
        if spec.detection.data_shape:
            lines.append(f"- one of the first {DATA_SHAPE_LINES} data lines matches in full:")
            lines.append(f"  - {_code(spec.detection.data_shape)}")
    lines.append("")
    if spec.header is HeaderMode.VCF:
        lines += ["**Columns:** the VCF fixed columns, then FORMAT and one column per sample.", ""]
    elif spec.columns:
        lines += ["**Columns**", "", "| Name | Type | Other spellings | Meaning |", "|---|---|---|---|"]
        for column in spec.columns:
            aliases = ", ".join(f"`{alias}`" for alias in column.get("aliases", ())) or ""
            lines.append(
                f"| `{column['name']}` | {dtype_name(column['dtype'])} | {_cell(aliases)} | "
                f"{_cell(column['description'])} |"
            )
        lines.append("")
    if spec.notes:
        lines += ["**From real files**", ""]
        lines += [f"- {note}" for note in spec.notes]
        lines.append("")
    if spec.example is not None:
        lines += ["**Example**", "", *_block(spec.example), ""]
    else:
        lines += ["**Example:** none; the format is binary.", ""]
    return lines


def _render_shape(shape: Shape) -> list[str]:
    found = f" ({', '.join(shape.findings)})" if shape.findings else ""
    lines = [f"### {shape.title}{found}", "", f"<a id=\"shape-{shape.key}\"></a>", "", shape.shape, ""]
    if shape.members:
        for name, text in shape.members:
            lines += [f"Member `{name}`:", "", *_block(text), ""]
    elif shape.example is not None:
        lines += [*_block(shape.example), ""]
    if shape.example_note:
        lines += [f"*{shape.example_note}*", ""]
    if shape.expect is not None:
        lines += [f"**`sniff` says:** {_expect_text(shape.expect)}", ""]
    lines += [f"**What just-vcf does** ({shape.status}): {shape.handling}", ""]
    pointers = [*(f"`{c}`" if "`" not in c else c for c in shape.code), *(f"`{t}`" for t in shape.tests)]
    if pointers:
        lines += [f"Code and tests: {', '.join(pointers)}.", ""]
    return lines


def _expect_text(expect: Expect) -> str:
    parts = [f"`{expect.outcome}`"]
    if expect.format:
        parts.append(f"as `{expect.format}`")
    if expect.layers:
        parts.append(f"inside {', '.join(f'`{layer}`' for layer in expect.layers)}")
    if expect.check_build:
        parts.append(f"build {expect.build}" if expect.build else "build unknown")
    if expect.reason_contains:
        parts.append(f"reason naming {expect.reason_contains}")
    return ", ".join(parts)


def _render_shapes() -> list[str]:
    lines = [
        "## Broken and odd shapes met in real files",
        "",
        "Shapes found in real uploads that a reader has to survive: the dogfooding findings F1-F17 "
        "(docs/DOGFOODING.md, docs/PREVIOUS_ISSUES.md) and the dialects behind the contract gaps in "
        "docs/FORMATS.md. Each states the shape, gives a synthetic example of it, says what `sniff` answers "
        "(the test suite runs every such example) and what just-vcf does with it. *handled* means read; "
        "*refused* means refused on purpose with a named error; *open* names a gap.",
        "",
    ]
    for shape in SHAPES:
        lines += _render_shape(shape)
    lines += ["### Findings that are not file shapes", ""]
    lines += [f"- **{finding}**: {text}" for finding, text in NOT_SHAPES]
    return [*lines, ""]


def render() -> str:
    """The whole of docs/FORMAT_REGISTRY.md."""
    lines = [
        "# Format registry",
        "",
        "Every input format just-vcf recognises, precise enough to write a sniffer for in another language, "
        "and the broken and odd shapes found in real files. The design behind the registry (why the order is "
        "what it is, how a format is declared, which contract gaps the loaders started from) is in "
        "[FORMATS.md](FORMATS.md); the frames a loaded file becomes are in [SCHEMA.md](SCHEMA.md).",
        "",
        "Generated from the format specs in `src/just_vcf/formats/`, `just_vcf.sniff` and the loaders by "
        "`scripts/format_registry.py`: edit those, then run `uv run python scripts/format_registry.py`. A test "
        "fails when this file and the code disagree.",
        "",
        "**Examples are synthetic.** A format's example is its spec's `example`, a head the spec declares "
        "synthetic: real dbSNP ids and positions, never a person's genotypes. The examples of the broken "
        "shapes were written for this document the same way: ids and positions are real sites on the build "
        "the example states (or build 36 where it says so); genotypes, alleles, scores, sample ids, dates and "
        "member names are invented and describe no person. Every example is a test vector: the suite sniffs "
        "a format's example and requires that format back, and no format earlier in the order. Tabs in "
        "examples are literal tab characters.",
        "",
        "## Contents",
        "",
        "1. [How a file is recognised](#how-a-file-is-recognised)",
        "2. [Vocabulary](#vocabulary)",
        "3. [Detection order](#detection-order)",
        "4. [Formats](#formats)",
        "5. [Broken and odd shapes met in real files](#broken-and-odd-shapes-met-in-real-files)",
        "",
        *_render_algorithm(),
        *_render_vocabulary(),
        *_render_index(),
        "## Formats",
        "",
    ]
    for position, spec in enumerate(REGISTRY, start=1):
        lines += _render_spec(position, spec)
    lines += _render_shapes()
    return "\n".join(lines).rstrip("\n") + "\n"


def main(check: bool = typer.Option(False, help="exit 1 when the doc differs from what the code renders")) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    text = render()
    if check:
        if not DOC.is_file() or DOC.read_text(encoding="utf-8") != text:
            logger.error("%s is stale: run uv run python scripts/format_registry.py", DOC.relative_to(ROOT))
            raise typer.Exit(1)
        logger.info("%s is current", DOC.relative_to(ROOT))
        return
    DOC.write_text(text, encoding="utf-8")
    logger.info("wrote %s (%d formats, %d shapes)", DOC.relative_to(ROOT), len(REGISTRY), len(SHAPES))


if __name__ == "__main__":
    typer.run(main)
