"""Write one facet back out as a VCF: the export at the edge (charter, RM6).

`write_vcf(loaded, path)` writes a `Loaded`'s `calls` rows as a single-sample VCF 4.3 (4.5 when a row
holds a PSL phase set), BGZF-compressed
by default. The frames stay the product; the file is one way out of them, written by this module's own
text writer rather than polars-bio's `sink_vcf`, because the writer has to own every byte: the GT text
with its per-haplotype phasing, the float formatting that makes two writes of one frame identical
(charter 12), and the declaration of every key it writes (polars-bio drops an undeclared column
silently on write, P13).

What is written, and the rules it follows:

- **REF** is `ref` (read from the assembly sequence) when there is one, else `source_ref` (what the
  file wrote). It is never taken from the genotype (charter 10). A row whose source REF states no
  base (`ref_check = unstated`, a copy-number caller's `N`, F23) is written with `ref`, the base the
  reference holds there, and its ALTs as the source wrote them when none is a sequence allele (a
  symbolic ALT, or none). A row with neither, or with the
  empty REF `""` (an insertion point whose anchor base `resolve_ref` did not read), a row without
  `alts` (an array row before RM3, or one RM3 withheld in a mixed-strand facet, F13), a row whose
  `ref_check` is `mismatch`, an `unstated` row with a sequence ALT (`ref_unstated`: spelled against the
  `N`, it would read back as a stated `match`), a row with an ALT that is no VCF allele (the empty allele `""` of a
  deletion written without its anchor base, F12's `-` among them, until `normalize` anchors it), a row with a genotype that
  is not indexed (`gt` null beside `gt_alleles`) and a row whose `gt` points past its ALTs cannot be
  written as valid VCF. The writer refuses them with `UnwritableRowsError`, counted per reason, unless
  the caller passes `resolved_only=True`; then the rest is written and the counts are reported.
- **FORMAT** is per row: a key is listed where the row has a value for it, so a field the source
  left out stays out. GT is rebuilt from `gt` and `phased` (`source_gt` is null on derived facets and
  VCF 4.4+ files). A haplotype phased by a CGI hapLink is written `/`: VCF has no field for the
  hapLink, and `|` without its set would read back as the implicit one. `gl` is written as GL: the
  frame holds log10 likelihoods and does not remember
  whether the source wrote PL. `gl`, `gp` and `ad` are written only when `source_alleles`, the order
  they are indexed over, is the written `[REF] + ALTs`; otherwise they are left out and counted.
  `hq` is written as HQ, declared `Number=2` as VCF reserves it: `.` per unscored allele, and a
  haploid row's second value missing (`51,.`, as CGI's vcfBeta writes it), which the loader reads back
  as the one value. A row with more than two GT alleles cannot be written under `Number=2`, and its
  `hq` is left out, counted.
- **VCF 4.3**, not 4.4+: GP is a probability from 4.3 on, which is what `gp` holds, and the writer
  emits no leading phase indicator, so the first haplotype's phase follows 4.3's implicit rule.
  **VCF 4.5 when a row holds a PSL phase set**, which 4.3 has no field for: FORMAT `PSL`, `PSO` and
  `PSQ` (`Number=P`) are written from `phase_sets`, `phase_set_ordinals` and `phase_set_qualities` on
  those rows, the set name on each slot in it and `.` elsewhere, and PS on the others. GT then carries
  a leading `|` or `/` where the first haplotype's known phase differs from the implicit one; just-vcf's
  own load reads that phase back only where PSL names the slot's set (P4), and the rows that lose it so
  are counted. CN is declared Float, as 4.5 reserves it, so a fractional copy number is written.
- **INFO** carries `END` where `end_pos` is not the REF span (symbolic alleles), the SV keys
  `SVTYPE`, `CIPOS`, `CIEND` and `CN` from `sv_type`, `cipos`, `ciend` and `alt_copy_numbers` (FORMAT
  `CN` is `copy_number`; a copy number that is not whole cannot be written as 4.3's Integer CN and
  is left out, counted), `SVLEN` on a length-type SV allele (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`,
  `<CNn>`) as `end_pos - pos`, `.` on the record's other ALTs: the frame keeps no SVLEN of its own,
  so this is END restated in the spelling VCF 4.5 asks for, never the source's SVLEN. VCF 4.3 has it
  `Number=.` and negative for a deletion, as its examples do; 4.5 `Number=A` and unsigned. A record
  with several length-type ALTs gets it only where VCF 4.5 says they share one length (copy-number
  alleles alone, or `<DEL>`/`<DUP>`/`<CNV>` beside FORMAT CN); elsewhere it is left out, counted. And
  `JUST_VCF_SPLIT` on rows flagged as one record of a split multi-allelic site (`split_signal`:
  `normalize`'s split rows, and rows RM270 flagged). Those rows are written as the frame holds them,
  with a warning and the tag, so a reader is told that a `0` there may mean another ALT. `load` reads
  the tag back as a marker (`split_signal = marker`), so a written normalized facet keeps its flag.
  `rejoin_split` joins them first, if that is what the caller wants.
- **Blocks** are written only with `blocks=True`, as gVCF reference blocks; otherwise they are
  counted. A block record needs the reference base at its start, which the `blocks` frame does not
  hold: it is read from the collection `resolve_ref` read the facet's REF from, never made up, and a
  facet whose REF was not read is refused. Each block is `REF <*>` (or `<NON_REF>`), INFO `END`,
  FORMAT GT, GQ, MIN_DP and GL where the frame holds them (it holds no block DP). GT is `0/0` for a
  `ref` block and `./.` for a `no_call` one, at the ploidy the block's own GL implies, else the
  ploidy of the contig's calls (the frame keeps no block GT; `load` reads only whether it called
  anything). Blocks and calls are merged one contig at a time, by position, calls first.
- **Contig names** are written as stored (the canonical `chrom`), with the `chr` prefix (`MT` becomes
  `chrM`), or as the source spelled them (`contig_map`). The header's `##contig` lines follow the same
  names and the rows' order. A breakend ALT's mate (`t[p[` and the other three forms, `p` being
  `contig:position`) is a contig of the file too: its contig is read canonically, as `load` reads
  CHROM, and written under the one mapping CHROM is written under, and it gets a `##contig` line
  whether or not a row stands on it (VCF 4.3 §1.4.7; F30). Under `stored` and `chr`, a mate the style
  has no other name for keeps its canonical name, as such a CHROM does. Under `source`, a mate is
  written as the facet spells its contig (`contig_map`, `header_contigs`), even where the ALT spelled
  it another way, so the file agrees with its own CHROM column; a contig the facet does not know (no
  row, no header line) is written as the ALT spelled it, and refused if ALTs spelled it two ways, as
  a CHROM spelled two ways is. Two canonical names written as one name (a lifted scaffold under the
  chain's name and the FASTA's) get one `##contig` line. Single breakends, symbolic alleles and
  brackets that are none of the four forms are written as held. So is a mate in angle
  brackets (`]<ctg1>:7]T`, VCF 4.3 §5.4.2): it names a contig of the `##assembly` file, not of the
  reference, and gets no `##contig` line (htslib refuses `<ctg1>` as a contig name). Only the written
  rows' mates are declared: a row `resolved_only` leaves out declares none.

- **A tabix index** (`<path>.tbi`) is written beside BGZF output by default, from the byte offsets
  the writer itself produced (`just_vcf.tabix`), so the file is region-queryable by htslib without a
  second pass over it. Each record's interval is the one htslib's VCF preset reads back: from `POS`,
  to INFO `END` where the line writes one, else over REF. A plain-text write has no index, and a
  write without one removes a stale `<path>.tbi`.

The header names the facet (`facet_id`, `source_sha256`, `format_id`, the parent of a derived facet)
and this package's version. It never names `source_path` or the source's `##SAMPLE` lines, which can
carry a real name; the sample column is the facet's `sample`, as given.
"""

import hashlib
import io
import logging
import os
from dataclasses import dataclass
from enum import Enum
from importlib.metadata import version
from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.bgzf import BgzfWriter
from just_vcf.breakends import assembly_mate, breakend_parts, mate_contig, mate_renamed
from just_vcf.loaders.vcf import _SV_PATTERN, FORMAT_KEYS, INFO_KEYS, SPLIT_INFO_KEY, V43_RESERVED_TYPES
from just_vcf.loading import Loaded
from just_vcf.reference import Contig, ReferenceMismatchError, ReferenceStore, _store_for
from just_vcf.schema import CONTIG_RANK, EMPTY_ALLELE, SCHEMA_VERSION, BlockKind, PhaseSetKind, RefCheck, canonical_chrom, unstated_sequence_alt
from just_vcf.tabix import TabixIndexer

logger = logging.getLogger(__name__)

VCF_VERSION = "VCFv4.3"

VCF_VERSION_PSL = "VCFv4.5"
"""The version written when a row holds a PSL phase set, which VCF 4.4 introduced: 4.3 has no PSL."""

_FORMAT_DESCRIPTIONS = {
    "GT": "Genotype",
    "GQ": "Conditional genotype quality",
    "HQ": "Haplotype qualities",
    "DP": "Read depth of this sample",
    "MIN_DP": "Minimum read depth over the block",
    "AD": "Read depth for each allele",
    "GL": "Genotype likelihoods, log10",
    "GP": "Genotype posterior probabilities",
    "DS": "Expected count of the alternate allele",
    "FT": "Sample filter",
    "PS": "Phase set",
    "PSL": "Phase set list",
    "PSO": "Phase set list ordinal",
    "PSQ": "Phase set list quality",
    "CN": "Copy number of this sample",
}
"""The FORMAT keys the writer can emit, in the order they are written. Number and Type come from the
loader's `FORMAT_KEYS`, so a written file declares each key the way the loader reads it back."""

_NO_STRINGS = pl.lit([], dtype=pl.List(pl.String))

_SV_INFO = {
    "SVLEN": (".", "Integer", "Length of the structural variant, END minus POS; under VCF 4.3 negative for a deletion"),
    "SVTYPE": ("1", "String", "Type of structural variant, as the source wrote it"),
    "CIPOS": (".", "Integer", "Confidence interval around POS, two values per ALT"),
    "CIEND": (".", "Integer", "Confidence interval around END, two values per ALT"),
    "CN": ("A", V43_RESERVED_TYPES["INFO/CN"], "Copy number of each ALT allele, as the source wrote it"),
}
"""The structural-variant INFO keys the writer emits (RM282, RM6), in the order they are written: SVLEN
as VCF 4.3 declares it (`Number=.`, signed), CIPOS and CIEND as VCF 4.5 declares them (two values per
ALT, which is how the frame holds them), CN with the Integer Type VCF 4.3 reserves, since polars-bio
refuses a 4.3 CN declared Float. A 4.5 file takes each key's Number and Type from `INFO_KEYS`: SVLEN
`Number=A` and unsigned, CN Float."""

_DELETION_ALLELE = r"^<DEL(:[^>]*)?>$"
_COPY_NUMBER_ALLELE = r"^<(CNV|CN\d+)(:[^>]*)?>$"
_SEGMENT_ALLELE = r"^<(DEL|DUP|CNV|CN\d+)(:[^>]*)?>$"
"""Symbolic alleles by what VCF 4.5 says of their SVLEN: every copy-number allele of a record has the
same SVLEN (the segment), and beside FORMAT CN so has every `<DEL>`, `<DUP>` and `<CNV>`."""


def _whole(value: pl.Expr, v45: bool = False) -> pl.Expr:
    """True where a copy number can be written. VCF 4.3 reserves CN as Integer, and polars-bio
    refuses a 4.3 header declaring it Float, so a fractional one is not written (counted). VCF 4.5
    reserves it as Float, so any value is."""
    return pl.lit(True) if v45 else value == value.round()


def _whole_list(name: str, v45: bool = False) -> pl.Expr:
    return pl.col(name).list.eval(_whole(pl.element(), v45) | pl.element().is_null()).list.all() & (pl.col(name).list.len() > 0)


def _copy_number_text(value: pl.Expr, v45: bool) -> pl.Expr:
    return _number_text(value) if v45 else value.cast(pl.Int64).cast(pl.String)


def _all_alts(pattern: str) -> pl.Expr:
    """True where every length-type SV ALT of the row matches `pattern`."""
    return pl.col("alts").list.eval(~pl.element().str.contains(_SV_PATTERN) | pl.element().str.contains(pattern)).list.all()


def _svlen_count() -> pl.Expr:
    """How many of the row's ALTs are length-type SV alleles (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<CNn>`)."""
    return pl.col("alts").list.eval(pl.element().str.contains(_SV_PATTERN)).list.sum().fill_null(0)


def _svlen_writable() -> pl.Expr:
    """True where `end_pos - pos` is each length-type SV ALT's own length, so SVLEN can be written.

    The frame keeps no SVLEN: `end_pos` is the record's extent, END or `POS + |SVLEN|`, the largest
    over its alleles. With one length-type ALT that is its length. With several it is each one's only
    where VCF 4.5 says they share it: copy-number alleles alone, or `<DEL>`, `<DUP>` and `<CNV>` beside
    FORMAT CN. Elsewhere the ALTs' own lengths are not kept apart, and SVLEN is left out (counted)
    rather than written as a length one of them may not have. REF must be the one padding base a
    symbolic record writes, so that `end_pos` is the event's end and not REF's."""
    count = _svlen_count()
    shared = (count == 1) | _all_alts(_COPY_NUMBER_ALLELE) | (pl.col("copy_number").is_not_null() & _all_alts(_SEGMENT_ALLELE))
    length = pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64)
    return ((count > 0) & shared & (_written_ref().str.len_chars() == 1) & (length > 0)).fill_null(False)


def _svlen_text(v45: bool) -> pl.Expr:
    """SVLEN per ALT: `end_pos - pos` on a length-type SV allele, `.` on any other (`<INS>`, a breakend,
    a sequence allele: VCF 4.5's missing value there, and the frame keeps no inserted length). VCF 4.3
    signs it as its own examples do, negative for a deletion; 4.5 writes it unsigned."""
    deletion = 1 if v45 else -1
    signs = pl.col("alts").list.eval(
        pl.when(pl.element().str.contains(_DELETION_ALLELE))
        .then(pl.lit(deletion, dtype=pl.Int64))
        .when(pl.element().str.contains(_SV_PATTERN))
        .then(pl.lit(1, dtype=pl.Int64))
        .otherwise(None)
    )
    values = signs * (pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64))
    return pl.when(_svlen_writable()).then(values.list.eval(pl.element().cast(pl.String).fill_null(".")).list.join(","))


def _sv_info_values(v45: bool = False) -> dict[str, pl.Expr]:
    """Each SV INFO key's `KEY=value` text per row, null where the row has no value for it."""
    texts = {
        "SVLEN": _svlen_text(v45),
        "SVTYPE": pl.col("sv_type"),
        "CIPOS": _list_text("cipos"),
        "CIEND": _list_text("ciend"),
        "CN": pl.when(_whole_list("alt_copy_numbers", v45)).then(
            pl.col("alt_copy_numbers").list.eval(_copy_number_text(pl.element(), v45).fill_null(".")).list.join(",")
        ),
    }
    return {key: pl.when(text.is_not_null()).then(pl.concat_str(pl.lit(f"{key}="), text)) for key, text in texts.items()}


_UNWRITABLE_REASONS = (
    "ref_unresolved",
    "ref_mismatch",
    "ref_unstated",
    "alts_unresolved",
    "alt_not_an_allele",
    "genotype_unindexed",
    "genotype_index_past_alts",
)

_VCF_ALLELE = r"^([ACGTNacgtn]+|\*|<[^<>]+>|\.[ACGTNacgtn]+|[ACGTNacgtn]+\.)$|\[|\]"
"""An ALT VCF can write: bases, `*`, a symbolic allele, a single breakend (`.A`, `A.`) or a breakend
(`A[2:500[`). Anything else is not one, the empty allele `""` of an unanchored deletion included
(F12's `-` loads as it), and htslib refuses the file."""


class ContigNames(str, Enum):
    """How contig names are spelled in a written VCF."""

    STORED = "stored"
    """The canonical `chrom`: `1`–`22`, `X`, `Y`, `MT`, `XY`; other contigs as the source spelled them."""
    CHR = "chr"
    """`chr1`–`chr22`, `chrX`, `chrY`, `chrM`; every other contig as stored."""
    SOURCE = "source"
    """As the source file spelled them, from the facet's `contig_map` and `header_contigs`; a breakend
    mate's contig the facet does not know, as the ALT spelled it."""


_CHR_NAMES = {**{str(n): f"chr{n}" for n in range(1, 23)}, "X": "chrX", "Y": "chrY", "MT": "chrM"}
"""Canonical → `chr`-prefixed, listed name by name: `MT` is `chrM`, and nothing is prefixed by pattern."""


class BlockAllele(str, Enum):
    """How a gVCF reference block's ALT is spelled. `load` reads both as `<*>`."""

    STAR = "<*>"
    """VCF 4.3's own symbol for any allele other than REF."""
    NON_REF = "<NON_REF>"
    """GATK's spelling, which GATK's own tools (GenotypeGVCFs, CombineGVCFs) expect."""


_BLOCK_ALT_DESCRIPTION = "Represents any possible alternative allele not already represented at this location by REF and ALT"


class UnwritableRowsError(ValueError):
    """Rows that cannot be written as valid VCF, counted per reason. `resolved_only=True` skips them."""

    def __init__(self, facet_id: str, counts: dict[str, int]) -> None:
        self.counts = counts
        detail = ", ".join(f"{count} {reason}" for reason, count in counts.items())
        super().__init__(
            f"{facet_id}: {sum(counts.values())} row(s) cannot be written as VCF ({detail}); "
            "pass resolved_only=True to write the rest"
        )


class UnsortedRowsError(ValueError):
    """The rows are not grouped by contig with positions ascending, which a VCF requires."""


class VcfWritten(BaseModel):
    """What `write_vcf` wrote."""

    model_config = ConfigDict(frozen=True)

    path: Path
    sha256: str
    """Digest of the bytes written, compressed or not."""
    rows_written: int
    rows_skipped: dict[str, int]
    """Rows left out under `resolved_only=True`, per reason. Empty when nothing was left out."""
    blocks_not_written: int
    """Blocks not written: all of them without `blocks=True`, else the ones counted in `rows_skipped`."""
    warnings: tuple[str, ...] = ()
    blocks_written: int = 0
    index_path: Path | None = None
    """The tabix index written beside a BGZF file (`<path>.tbi`); None when none was written."""


class _Hashing:
    """A binary sink that digests what passes through it."""

    def __init__(self, raw: BinaryIO) -> None:
        self._raw = raw
        self.digest = hashlib.sha256()

    def write(self, data: bytes) -> int:
        self.digest.update(data)
        return self._raw.write(data)


def write_vcf(
    loaded: Loaded,
    path: Path | str,
    *,
    contig_names: ContigNames = ContigNames.STORED,
    sample: str | None = None,
    resolved_only: bool = False,
    compress: bool = True,
    index: bool | None = None,
    blocks: bool = False,
    reference: "ReferenceStore | Path | str | None" = None,
    block_allele: BlockAllele = BlockAllele.STAR,
) -> VcfWritten:
    """Write `loaded`'s one facet to `path` as a VCF 4.3, or 4.5 where a row holds a PSL phase set
    (BGZF unless `compress=False`).

    `sample` names the sample column; the default is the facet's `sample`, then its `facet_id` for a
    source that names none (an array export). Rows that cannot be valid VCF raise
    `UnwritableRowsError` unless `resolved_only=True`. The rows must be in the order `load` emits
    them (grouped by contig, positions ascending); otherwise `UnsortedRowsError`. The same frames and
    options give the same bytes.

    `index` writes a tabix index, `<path>.tbi`, beside the file; the default is to write one whenever
    the output is BGZF. A plain-text file cannot be indexed (`index=True` with `compress=False` is
    refused), and a write without an index removes a `<path>.tbi` left from an earlier write, which
    would no longer describe the file.

    `blocks=True` writes the `blocks` frame too, as gVCF reference blocks interleaved with the calls
    (ALT `block_allele`, INFO `END`, FORMAT GT, GQ, MIN_DP and GL as the frame holds them). A block
    record needs the reference base at its start, which the frame does not hold: it is read from
    `reference` (a store or its directory; default the store at `reference_dir()`), from the
    collection `resolve_ref` read the facet's REF from. A facet whose REF was not read is refused;
    a block whose base cannot be read counts as `block_ref_unread`.
    """
    path = Path(path)
    if index and not compress:
        raise ValueError("a tabix index needs BGZF output; pass compress=True or leave index unset")
    indexed = compress if index is None else index
    index_path = path.with_name(path.name + ".tbi")
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"write_vcf takes one facet; this Loaded has {facets.height}")
    facet = facets.row(0, named=True)
    facet_id = facet["facet_id"]
    sample_name = sample if sample is not None else facet["sample"] if facet["sample"] is not None else facet_id
    if not sample_name or any(char in sample_name for char in "\t\n\r"):
        raise ValueError("the sample name must be non-empty and hold no tab or line break")

    calls = loaded.calls.with_columns(_unwritable().alias("_unwritable"))
    stats = _scan(calls)
    skipped = {reason: stats[f"n_{reason}"] for reason in _UNWRITABLE_REASONS if stats[f"n_{reason}"]}
    plan = _block_plan(loaded, facet, calls, reference, block_allele) if blocks else None
    if plan is not None and plan.unreadable:
        skipped["block_ref_unread"] = plan.unreadable
    if skipped and not resolved_only:
        raise UnwritableRowsError(facet_id, skipped)
    block_count = loaded.blocks.select(pl.len()).collect().item()
    blocks_written = block_count - plan.unreadable if plan is not None else 0

    # A breakend's mate is a contig of the file too (VCF 4.3 §1.4.7): it gets a `##contig` line and the
    # name CHROM would get, through the one mapping (F30); under `source`, a contig only ALTs name
    # keeps their spelling.
    mates = stats["mate_chroms"] or []
    chroms = sorted({*stats["chroms"], *(chrom for chrom, _ in mates), *(plan.chroms if plan is not None else [])})
    names = _contig_names(contig_names, facet, chroms, mates)
    warnings: list[str] = []
    if skipped:
        warnings.append(f"{sum(skipped.values())} row(s) not written: " + ", ".join(f"{n} {r}" for r, n in skipped.items()))
    if stats["n_vectors_dropped"]:
        warnings.append(
            f"{stats['n_vectors_dropped']} row(s) had GL/GP/AD indexed over another allele list than the written "
            "REF and ALTs; those vectors were left out"
        )
    if stats["n_hq_polyploid"]:
        warnings.append(
            f"{stats['n_hq_polyploid']} row(s) had haplotype qualities for more than two GT alleles, which HQ "
            "(Number=2) cannot write; left out"
        )
    if stats["n_hq_unfit"]:
        warnings.append(f"{stats['n_hq_unfit']} row(s) had haplotype qualities that are not one per GT allele; left out of HQ")
    if stats["n_ps_dropped"]:
        warnings.append(f"{stats['n_ps_dropped']} row(s) had a PS phase set that is not an integer, which PS cannot write; left out")
    if stats["n_haplink_unwritten"]:
        warnings.append(
            f"{stats['n_haplink_unwritten']} row(s) phased by a CGI hapLink, which VCF has no field for; their "
            "haplotypes written unphased"
        )
    if stats["v45"] and stats["n_first_phase_unread"]:
        warnings.append(
            f"{stats['n_first_phase_unread']} row(s) state their first haplotype's phase, written as VCF 4.5's leading "
            "GT indicator where the implicit one differs; just-vcf's own load reads it back as unknown unless PSL "
            "names the haplotype's set (polars-bio drops the indicator, P4)"
        )
    if stats["n_cn_fractional"]:
        warnings.append(
            f"{stats['n_cn_fractional']} row(s) had a copy number that is not a whole number, which VCF 4.3 "
            "cannot write (it reserves CN as Integer); left out"
        )
    if stats["n_svlen_withheld"]:
        warnings.append(
            f"{stats['n_svlen_withheld']} row(s) with SV alleles whose own lengths the frame does not keep apart; "
            "SVLEN left out, INFO END states the record's extent"
        )
    if stats["n_split"]:
        warnings.append(
            f"{stats['n_split']} row(s) look split one row per ALT (RM270) and were written as held, tagged "
            f"INFO/{SPLIT_INFO_KEY}: a `0` there may mean another ALT, not REF; rejoin_split joins them first"
        )
    if block_count and plan is None:
        warnings.append(
            f"{block_count} reference/no-call block(s) not written: pass blocks=True, with the reference store REF "
            "was read from, to write them as gVCF reference blocks"
        )

    in_blocks = plan.has if plan is not None and blocks_written else {}
    declared = [key for key in _FORMAT_DESCRIPTIONS if stats[f"has_{key}"] or in_blocks.get(key, False)]
    header_stats = {**stats, "has_END": stats["has_END"] or bool(blocks_written)}
    alt_lines = [f'##ALT=<ID={block_allele.value.strip("<>")},Description="{_BLOCK_ALT_DESCRIPTION}">'] if blocks_written else []
    header = _header(facet, sample_name, names, declared, header_stats, alt_lines)
    rows = calls.filter(pl.col("_unwritable").is_null()) if skipped else calls
    lines = rows.select(
        _line(names.mapping, stats["ploidy"], stats["v45"], bool(mates)).alias("line"),
        "chrom",
        "pos",
        *_index_columns(names.mapping),
    )
    if plan is None:
        batches: Iterator[pl.DataFrame] = lines.collect_batches()
    else:
        _check_order(rows.select("chrom", "pos"))
        _check_order(loaded.blocks.select("chrom", "pos"))
        batches = _interleaved(lines, plan, names.mapping)

    # Written beside the target and moved over it at the end, so a refused or failed write leaves no
    # half-written VCF (or index) behind under the requested name.
    partial = path.with_name(path.name + ".partial")
    index_partial = index_path.with_name(index_path.name + ".partial")
    indexer = TabixIndexer() if indexed else None
    try:
        with partial.open("wb") as raw:
            sink = _Hashing(raw)
            if compress:
                with BgzfWriter(sink, level=6) as out:
                    written = _write_body(out, header, batches, indexer)
                if indexer is not None:
                    with index_partial.open("wb") as raw_index, BgzfWriter(raw_index, level=6) as index_out:
                        index_out.write(indexer.build(out.virtual_offset))
            else:
                written = _write_body(sink, header, batches, None)
    except BaseException:
        partial.unlink(missing_ok=True)
        index_partial.unlink(missing_ok=True)
        raise
    os.replace(partial, path)
    if indexer is not None:
        os.replace(index_partial, index_path)
    else:
        index_path.unlink(missing_ok=True)

    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return VcfWritten(
        path=path,
        sha256=sink.digest.hexdigest(),
        rows_written=written - blocks_written,
        rows_skipped=skipped,
        blocks_not_written=block_count - blocks_written,
        blocks_written=blocks_written,
        warnings=tuple(warnings),
        index_path=index_path if indexer is not None else None,
    )


def _write_body(out: BinaryIO | BgzfWriter | _Hashing, header: str, batches: Iterator[pl.DataFrame], indexer: TabixIndexer | None) -> int:
    """The header, then every line batch by batch, checking the order as it goes and handing each
    record's interval and byte offsets to `indexer`. Returns the line count."""
    head = header.encode()
    out.write(head)
    offset = len(head)
    order = _OrderCheck()
    written = 0
    for batch in batches:
        order.check(batch)
        buffer = io.BytesIO()
        batch.select("line").write_csv(buffer, include_header=False, quote_style="never", line_terminator="\n")
        data = buffer.getvalue()
        out.write(data)
        if indexer is not None:
            size = pl.col("line").str.len_bytes().cast(pl.Int64) + 1
            located = batch.select(
                pl.col("_ix_name").alias("name"),
                pl.col("_ix_beg").alias("beg"),
                pl.col("_ix_end").alias("end"),
                (offset + size.cum_sum()).alias("u_end"),
            ).with_columns((pl.col("u_end") - batch.select(size).to_series()).alias("u_start"))
            if batch.height and located["u_end"][-1] != offset + len(data):
                raise RuntimeError("the written lines' sizes do not add up to the bytes written; the index would be wrong")
            indexer.add(located)
        offset += len(data)
        written += batch.height
    return written


def _index_columns(mapping: dict[str, str]) -> list[pl.Expr]:
    """Each record's contig as written and its interval as htslib's VCF preset reads it back: 0-based
    from `POS - 1`, to INFO `END` where the line writes one, else to `POS - 1 + len(REF)`."""
    ref_span_end = pl.col("pos").cast(pl.Int64) + _written_ref().str.len_chars().cast(pl.Int64) - 1
    end = pl.col("end_pos").cast(pl.Int64)
    return [
        pl.col("chrom").replace(mapping).alias("_ix_name"),
        (pl.col("pos").cast(pl.Int64) - 1).alias("_ix_beg"),
        pl.when((end != ref_span_end).fill_null(False)).then(end).otherwise(ref_span_end).alias("_ix_end"),
    ]


def _written_ref() -> pl.Expr:
    """REF as written: `ref` from the sequence when known, else what the file wrote."""
    return pl.coalesce(pl.col("ref"), pl.col("source_ref"))


def _unwritable() -> pl.Expr:
    """Why a row cannot be written as valid VCF, first reason first; null when it can."""
    return (
        # An empty REF is an insertion point whose anchor base `resolve_ref` did not read (CGI's native
        # files, a CASAVA indel table): VCF's REF is at least one base, and that base is the reference's,
        # so the row's REF is unresolved.
        pl.when(_written_ref().is_null() | (_written_ref() == EMPTY_ALLELE))
        .then(pl.lit("ref_unresolved"))
        .when(pl.col("ref_check") == RefCheck.MISMATCH.value)
        .then(pl.lit("ref_mismatch"))
        # A sequence ALT spelled against a source REF of only `N` (F23): written beside the reference's
        # base, `N` → `NAC` would read back as a stated `C` → `NAC`, `ref_check = match`.
        .when(unstated_sequence_alt())
        .then(pl.lit("ref_unstated"))
        .when(pl.col("alts").is_null())
        .then(pl.lit("alts_unresolved"))
        .when(~pl.col("alts").list.eval(pl.element().str.contains(_VCF_ALLELE)).list.all())
        .then(pl.lit("alt_not_an_allele"))
        .when(pl.col("gt").is_null() & pl.col("gt_alleles").is_not_null())
        .then(pl.lit("genotype_unindexed"))
        .when(pl.col("gt").list.max().cast(pl.Int64) > pl.col("alts").list.len().cast(pl.Int64))
        .then(pl.lit("genotype_index_past_alts"))
        .otherwise(None)
    )


def _vectors_indexable() -> pl.Expr:
    """True where `gl`, `gp` and `ad` are indexed over the written `[REF] + ALTs`."""
    return (pl.col("source_alleles") == pl.concat_list(_written_ref(), pl.col("alts"))).fill_null(False)


def _has_vector() -> pl.Expr:
    return pl.any_horizontal(pl.col(name).list.len().fill_null(0) > 0 for name in ("gl", "gp", "ad"))


def _has_haplink() -> pl.Expr:
    """True on a row with a slot phased by a CGI hapLink, which VCF has no field for: a hapLink names a
    molecule, not a GT slot, so it is no PS, and `|` without its set would read back as the
    chromosome's implicit one. Such a slot is written unphased."""
    return pl.col("phase_set_kinds").list.contains(PhaseSetKind.CGI_HAPLINK.value).fill_null(False)


def _phase(slot: int) -> pl.Expr:
    """The phase written for GT slot `slot`: `phased`, null where the slot is phased by a hapLink."""
    kind = pl.col("phase_set_kinds").list.get(slot, null_on_oob=True).cast(pl.String)
    phased = pl.col("phased").list.get(slot, null_on_oob=True)
    return pl.when((kind == PhaseSetKind.CGI_HAPLINK.value).fill_null(False)).then(None).otherwise(phased)


def _has_psl() -> pl.Expr:
    """True on a row with a slot in a PSL phase set."""
    return pl.col("phase_set_kinds").list.contains("psl").fill_null(False)


def _ps_text() -> pl.Expr:
    """PS: the row's PS-kind phase set, as an integer. Null where there is none or the row has a PSL
    set (VCF 4.5: a genotype has PS or PSL, never both; its implicit-set slots are PSL's `.`), or a
    CGI hapLink (`_has_haplink`)."""
    value = pl.col("phase_sets").list.drop_nulls().list.first()
    return pl.when(_has_psl() | _has_haplink()).then(None).otherwise(value.cast(pl.Int64, strict=False)).cast(pl.String)


def _psl_texts(ploidy: int) -> dict[str, pl.Expr]:
    """PSL, PSO and PSQ (`Number=P`, one value per GT allele) on a row with a PSL set: the set name on
    each of its slots and `.` elsewhere (unphased, or the implicit set); PSO and PSQ where the row
    holds any, `.` per slot without one. Null on every other row."""
    names = []
    for slot in range(ploidy):
        kind = pl.col("phase_set_kinds").list.get(slot, null_on_oob=True).cast(pl.String)
        name = pl.col("phase_sets").list.get(slot, null_on_oob=True)
        names.append(pl.when(pl.col("gt").list.len() > slot).then(pl.when(kind == "psl").then(name).otherwise(pl.lit("."))))
    return {
        "PSL": pl.when(_has_psl()).then(pl.concat_str(names, separator=",", ignore_nulls=True)),
        "PSO": pl.when(_has_psl()).then(_list_text("phase_set_ordinals")),
        "PSQ": pl.when(_has_psl()).then(_list_text("phase_set_qualities")),
    }


def _number_text(value: pl.Expr) -> pl.Expr:
    """A number as text: the shortest form that reads back to the same value, and `0.0` for `-0.0`
    (a PL of 0 becomes a GL of -0.0, which is the same likelihood)."""
    return pl.when(value == 0).then(value.abs()).otherwise(value).cast(pl.String)


def _list_text(name: str) -> pl.Expr:
    """A list column as VCF's comma-separated text, `.` per missing element; null when empty or null."""
    column = pl.col(name)
    text = column.list.eval(_number_text(pl.element()).fill_null(".")).list.join(",")
    return pl.when(column.list.len() > 0).then(text).otherwise(None)


def _gt_text(ploidy: int, v45: bool = False) -> pl.Expr:
    """GT from `gt` and `phased`: `.` per uncalled haplotype, `|` before a phased one, else `/`.

    Under VCF 4.5 the first haplotype's phase is its own: a leading `|` or `/` is written where
    `phased[0]` is known and differs from the implicit one (`/` if any later slot is `/`, else `|`),
    and none where it agrees or is unknown. VCF 4.3 has no leading indicator. A slot phased by a CGI
    hapLink is written unphased (`_phase`)."""
    gt = pl.col("gt")
    parts = [gt.list.get(0, null_on_oob=True).cast(pl.String).fill_null(".")]
    later_phased = []
    for slot in range(1, ploidy):
        phased = _phase(slot).fill_null(False)
        separator = pl.when(phased).then(pl.lit("|")).otherwise(pl.lit("/"))
        allele = gt.list.get(slot, null_on_oob=True).cast(pl.String).fill_null(".")
        parts.append(pl.when(gt.list.len() > slot).then(pl.concat_str(separator, allele)))
        later_phased.append(pl.when(gt.list.len() > slot).then(phased).otherwise(pl.lit(True)))
    if v45:
        implicit = pl.all_horizontal(later_phased) if later_phased else pl.lit(True)
        first = _phase(0)
        lead = pl.when(first.is_not_null() & (first != implicit)).then(pl.when(first).then(pl.lit("|")).otherwise(pl.lit("/")))
        parts.insert(0, lead)
    return pl.when(gt.list.len() > 0).then(pl.concat_str(parts, ignore_nulls=True))


def _hq_writable() -> pl.Expr:
    """True where `hq` holds a score and fits HQ's reserved `Number=2`: one value per GT allele, for
    one or two alleles."""
    hq = pl.col("hq")
    return (
        hq.list.drop_nulls().list.len().fill_null(0).gt(0) & (hq.list.len() == pl.col("gt").list.len()) & hq.list.len().is_between(1, 2)
    ).fill_null(False)


def _hq_text() -> pl.Expr:
    """HQ: `hq` with `.` per unscored allele, and a haploid row's second value written missing (`51,.`),
    as `Number=2` asks and as CGI's vcfBeta writes it; the loader reads it back as the one value."""
    text = pl.when(pl.col("hq").list.len() == 1).then(pl.concat_str(_list_text("hq"), pl.lit(",."))).otherwise(_list_text("hq"))
    return pl.when(_hq_writable()).then(text)


def _format_values(ploidy: int, v45: bool = False) -> dict[str, pl.Expr]:
    """Each FORMAT key's text per row, null where the row has no value for it."""
    vectors = _vectors_indexable()
    return {
        "GT": _gt_text(ploidy, v45),
        "GQ": pl.col("gq").cast(pl.String),
        "HQ": _hq_text(),
        "DP": pl.col("dp").cast(pl.String),
        "MIN_DP": pl.lit(None, dtype=pl.String),
        "AD": pl.when(vectors).then(_list_text("ad")),
        "GL": pl.when(vectors).then(_list_text("gl")),
        "GP": pl.when(vectors).then(_list_text("gp")),
        "DS": _number_text(pl.col("dosage")),
        "FT": pl.when(pl.col("sample_filters").list.len() > 0).then(pl.col("sample_filters").list.join(";")),
        "PS": _ps_text(),
        **_psl_texts(ploidy),
        "CN": pl.when(_whole(pl.col("copy_number"), v45)).then(_copy_number_text(pl.col("copy_number"), v45)),
    }


def _scan(calls: pl.LazyFrame) -> dict:
    """One pass for the header and the warnings: which keys occur, the filters, the contigs, the counts."""
    first = (
        calls.select(
            pl.col("gt").list.len().max().alias("ploidy"),
            _has_psl().any().alias("v45"),
            pl.col("alts").list.eval(_bracketed(pl.element())).list.any().any().alias("brackets"),
        )
        .collect()
        .row(0, named=True)
    )
    ploidy, v45 = first["ploidy"] or 1, bool(first["v45"])
    values = _format_values(ploidy, v45)
    ref_span_end = pl.col("pos").cast(pl.Int64) + _written_ref().str.len_chars().cast(pl.Int64) - 1
    row = (
        calls.select(
            *[(pl.col("_unwritable") == reason).sum().alias(f"n_{reason}") for reason in _UNWRITABLE_REASONS],
            *[text.is_not_null().any().alias(f"has_{key}") for key, text in values.items()],
            *[text.is_not_null().any().alias(f"has_info_{key}") for key, text in _sv_info_values(v45).items()],
            (_has_vector() & ~_vectors_indexable() & pl.col("_unwritable").is_null()).sum().alias("n_vectors_dropped"),
            *[
                (pl.col("hq").list.drop_nulls().list.len().fill_null(0).gt(0) & ~_hq_writable() & pl.col("_unwritable").is_null() & fits)
                .sum()
                .alias(name)
                for name, fits in (
                    ("n_hq_polyploid", (pl.col("hq").list.len() == pl.col("gt").list.len()).fill_null(False)),
                    ("n_hq_unfit", (pl.col("hq").list.len() != pl.col("gt").list.len()).fill_null(True)),
                )
            ],
            (pl.col("phase_sets").list.drop_nulls().list.len().fill_null(0).gt(0) & _ps_text().is_null() & ~_has_psl() & ~_has_haplink())
            .sum()
            .alias("n_ps_dropped"),
            _has_haplink().sum().alias("n_haplink_unwritten"),
            (
                _phase(0).is_not_null()
                & (pl.col("phase_set_kinds").list.get(0, null_on_oob=True).cast(pl.String) != "psl").fill_null(True)
            )
            .sum()
            .alias("n_first_phase_unread"),
            (
                (pl.col("copy_number").is_not_null() & ~_whole(pl.col("copy_number"), v45).fill_null(True))
                | (pl.col("alt_copy_numbers").list.len().fill_null(0).gt(0) & ~_whole_list("alt_copy_numbers", v45).fill_null(True))
            )
            .sum()
            .alias("n_cn_fractional"),
            ((_svlen_count() > 0) & ~_svlen_writable() & pl.col("_unwritable").is_null()).sum().alias("n_svlen_withheld"),
            pl.col("split_signal").is_not_null().sum().alias("n_split"),
            (pl.col("end_pos").cast(pl.Int64) != ref_span_end).fill_null(False).any().alias("has_END"),
            pl.concat_list(pl.col("filters").fill_null(_NO_STRINGS), pl.col("sample_filters").fill_null(_NO_STRINGS))
            .explode(empty_as_null=False)
            .drop_nulls()
            .unique()
            .sort()
            .implode()
            .alias("filters"),
            pl.col("chrom").unique().sort().implode().alias("chroms"),
        )
        .collect()
        .row(0, named=True)
    )
    row["ploidy"] = ploidy
    row["v45"] = v45
    row["mate_chroms"] = _mate_chroms(calls) if first["brackets"] else []
    return row


def _bracketed(alt: pl.Expr) -> pl.Expr:
    """True on an ALT holding a bracket: a breakend, or brackets that are none of its forms."""
    return alt.str.contains("[", literal=True) | alt.str.contains("]", literal=True)


def _mate_chroms(calls: pl.LazyFrame) -> list[tuple[str, str]]:
    """Each (canonical contig, spelling) the written rows' breakend ALTs name as mates, the canonical
    as `breakend_parts` reads it and the spelling as the ALT wrote it (`source` writes a contig the
    facet does not know by it). A row left out as unwritable (`_unwritable`) declares nothing, and a
    mate on a contig of the assembly file (`<ctg1>`, `breakends.assembly_mate`) is no reference
    contig and gets no `##contig` line. A pass of its own over the bracketed ALTs alone, run only
    when `_scan`'s first pass found one (on 4.8 million rows without a breakend, the bracket search
    0.1 s, this pass 0.2 s; the regex read inside each row's list, 1.9 s). Not an expression in
    `_scan`'s one select: over `alts` exploded there, polars 1.44's common-subexpression pass mixes
    the exploded length with the rows' and raises a ShapeError."""
    alts = (
        calls.filter(pl.col("_unwritable").is_null())
        .select(pl.col("alts").explode(empty_as_null=False).alias("_alt"))
        .filter(_bracketed(pl.col("_alt")))
        .filter(~assembly_mate(pl.col("_alt")).fill_null(False))
    )
    # `_mate_chrom` is null on an ALT that is none of the four forms, so `drop_nulls` keeps the mates.
    _, canonical, _ = breakend_parts(pl.col("_alt"))
    found = (
        alts.select(canonical, mate_contig(pl.col("_alt")).alias("_spelled"))
        .drop_nulls()
        .unique()
        .sort("_mate_chrom", "_spelled")
        .collect()
    )
    return list(found.iter_rows())


@dataclass(frozen=True)
class _Names:
    """Canonical `chrom` → written name, and the `##contig` lines in the rows' order."""

    mapping: dict[str, str]
    contigs: list[dict]


def _contig_names(style: ContigNames, facet: dict, chroms: list[str], mates: list[tuple[str, str]] | None = None) -> _Names:
    """The written name of every contig in the rows, in the facet's `##contig` lines and among the
    breakend mates (`mates`: each canonical contig with a spelling an ALT gave it, `_mate_chroms`).

    Under `source`, a contig is written as the facet spells it (`contig_map`, then `header_contigs`),
    the one spelling its CHROM column uses, even where an ALT spelled it another way: the file stays
    consistent with its own CHROM column. A contig only ALTs name keeps the ALTs' spelling, not its
    canonical name. A contig spelled two ways (two CHROM or header spellings, or two ALT spellings
    of a contig only ALTs name) is refused: which one the source meant cannot be told."""
    if style is ContigNames.SOURCE:
        # `contig_map` holds the contigs that have rows; a `##contig` with none is spelled by the header.
        by_chrom: dict[str, set[str]] = {}
        for entry in facet["contig_map"] or []:
            by_chrom.setdefault(entry["chrom"], set()).add(entry["source"])
        for contig in facet["header_contigs"] or []:
            by_chrom.setdefault(canonical_chrom(contig["name"]), set()).add(contig["name"])
        # Only after the facet's own spellings, so an ALT's never competes with CHROM's.
        by_alt: dict[str, set[str]] = {}
        for chrom, spelled in mates or []:
            if chrom not in by_chrom:
                by_alt.setdefault(chrom, set()).add(spelled)
        by_chrom.update(by_alt)
        clashes = sorted(chrom for chrom, sources in by_chrom.items() if len(sources) > 1)
        if clashes:
            raise ValueError(f"the source spelled {', '.join(clashes)} more than one way; its spelling cannot be restored")
        source = {chrom: sources.pop() for chrom, sources in by_chrom.items()}

    def written(chrom: str) -> str:
        if style is ContigNames.CHR:
            return _CHR_NAMES.get(chrom, chrom)
        if style is ContigNames.SOURCE:
            return source.get(chrom, chrom)
        return chrom

    declared: dict[str, dict] = {}
    for contig in facet["header_contigs"] or []:
        chrom = canonical_chrom(contig["name"])
        declared.setdefault(chrom, {"length": contig["length"], "md5": contig["md5"]})
    for chrom in chroms:
        declared.setdefault(chrom, {"length": None, "md5": None})
    ordered = sorted(declared, key=lambda chrom: (CONTIG_RANK.get(chrom, len(CONTIG_RANK)), chrom))
    # Two canonical names can be written as one: a lifted row's contig under the chain's name and its
    # mate's under the target FASTA's (`chrUn_KI270742v1` and `KI270742.1`), both `KI270742.1` under
    # `source`. That is one contig, so it gets one `##contig` line, its first stated length and md5.
    contigs: dict[str, dict] = {}
    for chrom in ordered:
        line = contigs.setdefault(written(chrom), {"name": written(chrom), "length": None, "md5": None})
        for key in ("length", "md5"):
            if line[key] is None:
                line[key] = declared[chrom][key]
    return _Names(mapping={chrom: written(chrom) for chrom in ordered}, contigs=list(contigs.values()))


def _header(facet: dict, sample: str, names: _Names, declared: list[str], stats: dict, alt_lines: list[str]) -> str:
    lines = [
        f"##fileformat={VCF_VERSION_PSL if stats['v45'] else VCF_VERSION}",
        f"##source=just-vcf {version('just-vcf')}",
        f"##just_vcf_schema_version={SCHEMA_VERSION}",
        f"##just_vcf_facet_id={facet['facet_id']}",
    ]
    for key in ("source_sha256", "format_id", "derived_from", "relation", "assembly_label"):
        if facet[key] is not None:
            lines.append(f"##just_vcf_{key}={facet[key]}")
    filters = stats["filters"] or []
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
    lines += [
        f'##FILTER=<ID={name},Description="Carried from the source; its description was not kept">'
        for name in filters
        if name != "PASS"
    ]
    if stats["has_END"]:
        lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the record">')
    for key, (number, kind, description) in _SV_INFO.items():
        number, kind = INFO_KEYS[key] if stats["v45"] and key in INFO_KEYS else (number, kind)
        if stats[f"has_info_{key}"]:
            lines.append(f'##INFO=<ID={key},Number={number},Type={kind},Description="{description}">')
    if stats["n_split"]:
        lines.append(
            f'##INFO=<ID={SPLIT_INFO_KEY},Number=1,Type=String,Description="This row looks like one record of a '
            'multi-allelic site written one row per ALT (just-vcf RM270): a 0 in GT may mean another ALT, not REF">'
        )
    for key in declared:
        number, kind = FORMAT_KEYS[key]
        kind = kind if stats["v45"] else V43_RESERVED_TYPES.get(f"FORMAT/{key}", kind)
        lines.append(f'##FORMAT=<ID={key},Number={number},Type={kind},Description="{_FORMAT_DESCRIPTIONS[key]}">')
    lines += alt_lines
    for contig in names.contigs:
        attrs = f"ID={contig['name']}"
        if contig["length"] is not None:
            attrs += f",length={contig['length']}"
        if contig["md5"]:
            attrs += f",md5={contig['md5']}"
        lines.append(f"##contig=<{attrs}>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]))
    return "\n".join(lines) + "\n"


def _line(mapping: dict[str, str], ploidy: int, v45: bool = False, mates: bool = False) -> pl.Expr:
    """One VCF data line per row. With `mates`, each breakend ALT's mate contig is written through
    `mapping`, as CHROM is (F30); without, no ALT names one and the ALTs are written as held."""
    ref = _written_ref()
    alts = pl.col("alts").list.eval(mate_renamed(pl.element(), mapping)) if mates else pl.col("alts")
    ref_span_end = pl.col("pos").cast(pl.Int64) + ref.str.len_chars().cast(pl.Int64) - 1
    end = pl.when(pl.col("end_pos").cast(pl.Int64) != ref_span_end).then(pl.concat_str(pl.lit("END="), pl.col("end_pos").cast(pl.String)))
    split = pl.when(pl.col("split_signal").is_not_null()).then(
        pl.concat_str(pl.lit(f"{SPLIT_INFO_KEY}="), pl.col("split_signal").cast(pl.String))
    )
    info = pl.concat_str([end, *_sv_info_values(v45).values(), split], separator=";", ignore_nulls=True)
    values = _format_values(ploidy, v45)
    keys = pl.concat_str([pl.when(text.is_not_null()).then(pl.lit(key)) for key, text in values.items()], separator=":", ignore_nulls=True)
    sample = pl.concat_str(list(values.values()), separator=":", ignore_nulls=True)
    # A row with no FORMAT value at all still needs both columns: GT, uncalled.
    nothing = keys == ""
    return pl.concat_str(
        [
            pl.col("chrom").replace(mapping),
            pl.col("pos").cast(pl.String),
            pl.when(pl.col("ids").list.len() > 0).then(pl.col("ids").list.join(";")).otherwise(pl.lit(".")),
            ref,
            pl.when(alts.list.len() > 0).then(alts.list.join(",")).otherwise(pl.lit(".")),
            _number_text(pl.col("qual")).fill_null("."),
            pl.when(pl.col("filters").list.len() > 0).then(pl.col("filters").list.join(";")).otherwise(pl.lit(".")),
            pl.when(info == "").then(pl.lit(".")).otherwise(info),
            pl.when(nothing).then(pl.lit("GT")).otherwise(keys),
            pl.when(nothing).then(pl.lit(".")).otherwise(sample),
        ],
        separator="\t",
    )


class _OrderCheck:
    """Raises `UnsortedRowsError` when a contig comes back after another, or a position goes down."""

    def __init__(self) -> None:
        self._seen: set[str] = set()
        self._last: tuple[str, int] | None = None

    def check(self, batch: pl.DataFrame) -> None:
        if batch.is_empty():
            return
        chroms = batch["chrom"]
        positions = batch["pos"].cast(pl.Int64)
        starts = (chroms != chroms.shift(1)).fill_null(True)
        if ((positions.diff().fill_null(0) < 0) & ~starts).any():
            raise UnsortedRowsError("positions go down within a contig; write the frames in the order load emits them")
        for index in starts.arg_true().to_list():
            chrom, pos = chroms[index], positions[index]
            if index == 0 and self._last is not None and chrom == self._last[0]:
                if pos < self._last[1]:
                    raise UnsortedRowsError(f"positions go down within contig {chrom}")
                continue
            if chrom in self._seen:
                raise UnsortedRowsError(f"contig {chrom} comes back after another; the rows are not grouped by contig")
            self._seen.add(chrom)
        self._last = (chroms[-1], positions[-1])


def _check_order(rows: pl.LazyFrame) -> None:
    """`UnsortedRowsError` unless `rows` are grouped by contig with positions ascending, read as a stream."""
    order = _OrderCheck()
    for batch in rows.collect_batches():
        order.check(batch)


@dataclass(frozen=True)
class _BlockPlan:
    """What writing the `blocks` frame needs: the reference contigs to read a base from, the ploidy a
    block's GT is written with when its GL does not say, and which FORMAT keys the frame holds."""

    blocks: pl.LazyFrame
    """The blocks, with `_readable` saying whether the base at `pos` can be read."""
    store: ReferenceStore
    contigs: dict[str, Contig]
    ploidy: dict[str, int]
    fallback_ploidy: int
    allele: BlockAllele
    chroms: list[str]
    unreadable: int
    has: dict[str, bool]


def _block_plan(loaded: Loaded, facet: dict, calls: pl.LazyFrame, reference: "ReferenceStore | Path | str | None", allele: BlockAllele) -> _BlockPlan:
    """Check that each block's REF base can be read, contig by contig, before anything is written."""
    digest = facet["ref_seqcol_digest"]
    if digest is None:
        raise ValueError(
            f"facet {facet['facet_id']}: REF has not been read; a block's REF base comes from the facet's assembly, "
            "so run resolve_ref first, then write_vcf(blocks=True, reference=…)"
        )
    store, _ = _store_for(reference)
    if digest not in store.collections():
        raise ReferenceMismatchError(
            f"facet {facet['facet_id']}: the reference store does not hold collection {digest}, which REF was read "
            "from; pass reference= with the store resolve_ref read"
        )
    # A contig the collection names once, whose length the facet stated the same or not at all: the
    # check resolve_ref made for the calls' contigs, made again for contigs only blocks use.
    stated: dict[str, set[int]] = {}
    for entry in facet["header_contigs"] or []:
        if entry["length"] is not None:
            stated.setdefault(canonical_chrom(entry["name"]), set()).add(entry["length"])
    contigs = {
        chrom: found[0]
        for chrom, found in store.contigs_by_chrom(digest).items()
        if len(found) == 1 and stated.get(chrom, {found[0].length}) == {found[0].length}
    }
    lengths = {chrom: contig.length for chrom, contig in contigs.items()}
    readable = pl.col("pos").cast(pl.Int64) <= pl.col("chrom").replace_strict(lengths, default=0, return_dtype=pl.Int64)
    blocks = loaded.blocks.with_columns(readable.alias("_readable"))
    summary = (
        blocks.select(
            (~pl.col("_readable")).sum().alias("unreadable"),
            pl.col("gq").is_not_null().any().alias("GQ"),
            pl.col("dp_min").is_not_null().any().alias("MIN_DP"),
            (pl.col("gl").list.len() > 0).any().alias("GL"),
            pl.col("chrom").unique().sort().implode().alias("chroms"),
        )
        .collect()
        .row(0, named=True)
    )
    per_contig = calls.group_by("chrom").agg(pl.col("gt").list.len().max().alias("ploidy")).drop_nulls().collect()
    ploidy = dict(per_contig.iter_rows())
    return _BlockPlan(
        blocks=blocks,
        store=store,
        contigs=contigs,
        ploidy=ploidy,
        fallback_ploidy=max(ploidy.values(), default=2),
        allele=allele,
        chroms=summary["chroms"] or [],
        unreadable=summary["unreadable"] or 0,
        has={"GT": True, "GQ": bool(summary["GQ"]), "MIN_DP": bool(summary["MIN_DP"]), "GL": bool(summary["GL"])},
    )


def _interleaved(lines: pl.LazyFrame, plan: _BlockPlan, mapping: dict[str, str]) -> Iterator[pl.DataFrame]:
    """One contig at a time, in the header's order: its call lines and its block lines merged by
    position, calls first where both start at one position. One contig is held in memory at a time."""
    for chrom in mapping:
        calls = lines.filter(pl.col("chrom") == chrom).collect()
        blocks = _block_lines(plan, chrom, mapping)
        if calls.is_empty() and blocks.is_empty():
            continue
        merged = pl.concat([calls.with_columns(pl.lit(0).alias("_source")), blocks.with_columns(pl.lit(1).alias("_source"))], how="vertical")
        yield merged.sort("pos", "_source", maintain_order=True).drop("_source")


_LINE_SCHEMA = {"line": pl.String, "chrom": pl.String, "pos": pl.UInt32, "_ix_name": pl.String, "_ix_beg": pl.Int64, "_ix_end": pl.Int64}


def _block_lines(plan: _BlockPlan, chrom: str, mapping: dict[str, str]) -> pl.DataFrame:
    """One contig's writable blocks as VCF lines, with the REF base read at each block's start."""
    rows = plan.blocks.filter((pl.col("chrom") == chrom) & pl.col("_readable")).collect()
    if rows.is_empty():
        return pl.DataFrame(schema=_LINE_SCHEMA)
    starts = rows["pos"].unique().sort()
    # `pos` is 1-based; the store takes 0-based half-open ranges (charter 13).
    bases = plan.store.fetch(plan.contigs[chrom], [(pos - 1, pos) for pos in starts.to_list()])
    table = pl.DataFrame({"pos": starts, "_ref": [base.upper() for base in bases]})
    rows = rows.join(table, on="pos", how="left", maintain_order="left")
    # GT is not in the frame: a `ref` block asserts reference, a `no_call` block nothing. Its ploidy is
    # the block's own GL's (k + 1 genotypes over [REF, <*>] for ploidy k), else the contig's calls'.
    ploidy = (
        pl.when(pl.col("gl").list.len() > 1)
        .then(pl.col("gl").list.len().cast(pl.Int64) - 1)
        .otherwise(pl.lit(plan.ploidy.get(chrom, plan.fallback_ploidy), dtype=pl.Int64))
    )
    allele = pl.when(pl.col("kind") == BlockKind.REF.value).then(pl.lit("0")).otherwise(pl.lit("."))
    gt = pl.concat_str(allele, pl.concat_str(pl.lit("/"), allele).repeat_by(ploidy - 1).list.join(""))
    values = {
        "GT": gt,
        "GQ": pl.col("gq").cast(pl.String),
        "MIN_DP": pl.col("dp_min").cast(pl.String),
        "GL": _list_text("gl"),
    }
    keys = pl.concat_str([pl.when(text.is_not_null()).then(pl.lit(key)) for key, text in values.items()], separator=":", ignore_nulls=True)
    sample = pl.concat_str(list(values.values()), separator=":", ignore_nulls=True)
    name = pl.col("chrom").replace(mapping)
    return rows.select(
        pl.concat_str(
            [
                name,
                pl.col("pos").cast(pl.String),
                pl.lit("."),
                pl.col("_ref"),
                pl.lit(plan.allele.value),
                pl.lit("."),
                pl.lit("."),
                pl.concat_str(pl.lit("END="), pl.col("end_pos").cast(pl.String)),
                keys,
                sample,
            ],
            separator="\t",
        ).alias("line"),
        "chrom",
        "pos",
        name.alias("_ix_name"),
        (pl.col("pos").cast(pl.Int64) - 1).alias("_ix_beg"),
        pl.col("end_pos").cast(pl.Int64).alias("_ix_end"),
    )
