"""Lift a facet to another assembly through a chain file, as a derived facet (RM6, charter invariant 5).

`liftover(loaded, chain, reference)` returns a `Lifted`: the derived facet (`relation = liftover_to`,
the Cat-VRS term, `derived_from` the parent, which is never changed) and a report with one row per
parent row saying where it went or why it did not. The rules are in docs/SCHEMA.md ("Lifting"); in
short:

1. **The chain is checked against both assemblies first.** A chain header states each contig's size
   on both sides (`tSize`, `qSize`). The source sizes must equal the collection the parent's REF was
   read from, the target sizes the target collection, contig by contig through `canonical_chrom` and,
   for a GRCh37/38 scaffold, the names NCBI's assembly reports give it (`ReferenceStore.contigs_by_chrom`,
   F19: the chain's `chr1_KI270706v1_random` meets the store's `KI270706.1` where the sizes agree). A
   primary contig (1-22, X, Y) that differs refuses the lift: the chain is for another assembly. Any
   other contig that differs leaves its rows unlifted, counted (UCSC hg19's chrM, 16,571 bases,
   against GRCh37's rCRS MT).
2. **Positions.** A row's REF span must lie inside one ungapped block, and inside no block with
   another target (chains that overlap on the source side with different targets make a position's
   target ambiguous; blocks repeating one mapping count once, F14).
   Otherwise it is unlifted: no block at all, a span that crosses a block edge, or more than one
   target. **Across a gap of one chain** (2b) a sequence row still lifts where its first and last
   bases lie in blocks and one chain holds every piece of the span a block covers (whichever other
   chain repeats some of those blocks, F14), every gap between has fewer than `GAP_LIMIT` bases on each side,
   and each ALT's change (the bases it changes, or for an indel in a repeat every place it could
   equally sit) lies inside one block. Its lifted span runs from the image of its first base to the
   image of its last; REF is read there, and each allele is the target's sequence there with that
   allele's change put in, so the change against the source is the same change against the target
   and the gap's own bases are the target's in every allele. Otherwise `split_across_chains`,
   `gap_in_change` or `gap_too_large`.
3. **Alleles.** On a minus-strand block every allele is reverse-complemented. An indel whose alleles
   all start with one base (VCF's padding base) then has that base at their end instead, so it is
   re-anchored: the base is dropped from every allele and the target's base before the lifted span
   is put in front, which moves the row one base left. The indel is not left-aligned on the target
   here; `normalize` on the lifted facet does that.
   REF is then read from the target sequence over the lifted span (charter 10: never carried over).
   Where it equals the parent's REF the row keeps its alleles; where it equals one of the ALTs, REF
   and that ALT swap and the genotype indices follow (`bcftools +liftover`'s and Picard's swap);
   where it equals none, the row is unlifted (`ref_mismatch_after_lift`). "The parent's REF" is
   `ref`, the bases read from the source assembly: a row whose source REF contradicts them
   (`ref_check = mismatch`) is not lifted at all (`source_ref_mismatch`), and one whose source REF
   states no base (`unstated`: a copy-number caller's `N`, F23) lifts as a `match` row does when no
   ALT is a sequence allele (a symbolic SV, a segment, no ALT). One with a sequence ALT stays behind
   (`source_ref_unstated`): that ALT was spelled against the `N`, so comparing the target's base with
   it, swapping and re-indexing the genotype would set an allele the source never placed against a
   base (`N` → `G`, `0/0`, over `A` where the target reads `G`, would lift as `G` → `A`, `1/1`).
4. **Blocks.** A reference or no-call block's span is cut wherever the chain's mapping changes
   (the same axis as step 2): each piece inside one block with one target lifts, with the block's
   `kind`, `gq`, `dp_min` and `gl` unchanged; the pieces in no block, in more than one target, or on
   a target contig the collection lacks stay behind, counted in bases. A block that maps whole to one
   interval lifts whole; one crossing a chain gap or edge becomes several blocks, never one stretched
   over the gap. A piece on a minus-strand block or another contig lifts like a row does (a block has
   no allele to complement). The frame keeps no block REF: an export reads it from the target
   (`write_vcf(..., blocks=True)` reads `ref_seqcol_digest`, which is the target's).
5. **Symbolic SVs** (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<CNn>`, `<INS>`, subtypes included, every
   ALT of the row one of them) and copy-number segments (no ALT, an `end_pos` past REF): POS and
   `end_pos` are lifted as two points, independently, and the row is kept only where both land, on
   one target contig and strand, in their order, with a length (`end_pos - pos`) within
   `SV_LENGTH_TOLERANCE` of the source's (liftoverSV's rule). On a minus-strand block the event's
   image is reversed: a padded SV's new POS is the target's base before the image of END, a segment's
   the image of END itself. REF (one base) is read from the target there, and on a plus-strand block
   must equal the source's. The ALTs stay as written; `cipos` and `ciend` carry, mirrored and swapped
   on a minus-strand block.
6. **Breakends** (`t[p[`, `t]p]`, `]p]t`, `[p[t`, single `t.`/`.t`): POS lifts as one base and the
   mate position inside each ALT through the same chain, on whatever contig it names; the ALT is
   rewritten with the mate as the target collection names it. An end landing on a minus-strand block
   flips the side it keeps (`breakends.turned`): the record's flip moves `t`, reverse-complemented,
   across the brackets, the mate's turns them. REF, the target's base at the new POS, must equal the
   source's (complemented on a minus block). A mate on the record's own source contig must land on
   its target contig, in the order the strand gives, at a distance within `SV_LENGTH_TOLERANCE`
   (liftoverSV's checks). `cipos` mirrors with POS's strand and `ciend` with the mate's; they never
   swap, since POS stays POS. A mate in angle brackets (`C[<ctg1>:7[`, VCF 4.3 §5.4.2) is a contig of
   the `##assembly` file, not of the reference: no chain moves it, so it is kept as written, and the
   form turns by POS's strand alone (`breakends.turned` with the mate on a plus strand), `ciend` as it was.
7. **Withheld in this slice**, each counted: a minus-strand indel or SV whose new anchor would be
   before the target contig's first base, breakend ALTs that are none of the forms or sit beside
   another kind of allele, other symbolic alleles and `*`, and rows flagged as one record of a
   pre-split site.

A lifted row's `source_ref`, `ref_check`, `source_gt` and VRS ids are null, and so are a sequence row's
`cipos` and `ciend`: they described the source's representation on the source assembly (an SV's
intervals are about its two ends, which step 5 places, so they carry). `gl`, `gp` and `ad` stay indexed over
`source_alleles`, never reordered (SCHEMA.md Decision 2); a swapped row's `dosage` is withheld, since it
counted the other allele, and so are its `alt_copy_numbers` (INFO CN, one per ALT): the swapped ALT's
slot now holds the old REF, whose copy number the source did not state. Re-mint ids on the target with `mint_vrs_ids`.

Coordinates: `pos` is 1-based; chains and the store are 0-based half-open. `reference._span` converts
a row to the store's range and `_to_pos` converts a lifted start back, the only two places (charter 13).
"""

import gzip
import hashlib
import json
import logging
from collections import Counter
from collections.abc import Callable
from itertools import pairwise
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf import cache
from just_vcf.breakends import parse_breakend, turned
from just_vcf.links import LIFTOVER_TO, derived_links
from just_vcf.loaders.vcf import _SV_LENGTH_TYPES
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir, canonical_batches, canonical_parts
from just_vcf.reference import PRIMARY, Contig, ReferenceMismatchError, ReferenceStore, _derive_array_alleles, _span, _store_for
from just_vcf.schema import BLOCKS, CALLS, FACETS, SCHEMA_VERSION, AbsenceMeans, AssemblyLabelBasis, RefCheck, canonical_chrom, empty, unstated_sequence_alt

logger = logging.getLogger(__name__)

LIFT_VERSION = 15
"""Part of the derived `facet_id`: bumped whenever the rule changes what `liftover` emits (3: a swapped
row's `alt_copy_numbers` withheld; 4: F19's scaffold aliases on the target; 5: schema 0.0.15's
`phase_set_ordinals` and `phase_set_qualities` carried, null with the phase on another contig;
6: the lifted facet's `reference` is null, since it named the source's file; 7: schema 0.0.19's `hq` carried;
8: the `blocks` frame lifted; 9: symbolic SVs and copy-number segments lifted by both ends; 10: on a
minus strand, a source SV with CIEND and no CIPOS gets `ciend = [0, 0]`, its precise POS, not null;
11: breakends lifted by both ends, the mate rewritten inside ALT and the brackets turned per strand;
12: a breakend whose mate is a contig of the `##assembly` file, `<ctg1>`, lifted by POS alone, the mate
kept as written; 13: a row whose source REF states no base, `ref_check = unstated`, with a sequence ALT
left behind as `source_ref_unstated`; 14: a sequence row whose REF span crosses a gap of one chain lifted,
REF read over the lifted span and each ALT's change carried into the target's sequence there; 15: such a
row judged by the chains that hold each piece of its span, not by the one block per mapping the axis kept,
so a chain repeating only part of another's blocks no longer leaves it `split_across_chains`)."""

RELATION = LIFTOVER_TO
"""The derived facet's `relation`: the Cat-VRS term, verbatim (charter, GA4GH from the start)."""

UNLIFTED_REASONS: frozenset[str] = frozenset(
    {
        # What the parent row is.
        "ref_unread",
        "source_ref_mismatch",
        "source_ref_unstated",
        "allele_not_sequence",
        "split_flagged",
        # What the chain says about its position.
        "contig_not_in_chain",
        "chain_contig_mismatch",
        "no_chain_block",
        "split_across_blocks",
        "multiple_targets",
        "target_contig_unmatched",
        # A span across a gap of one chain (step 2b).
        "split_across_chains",
        "gap_in_change",
        "gap_too_large",
        # What the target says about its alleles.
        "indel_minus_strand",
        "ref_mismatch_after_lift",
        # Symbolic SVs (step 5).
        "sv_end_unlifted",
        "sv_ends_apart",
        "sv_ends_reversed",
        "sv_length_changed",
        # Breakends (step 6).
        "breakend_not_lifted",
        "breakend_mate_not_lifted",
        "breakend_ends_apart",
        "breakend_ends_reversed",
        "breakend_distance_changed",
    }
)
"""Why a parent row is not in the lifted facet; the report's `lift_reason`. An open vocabulary (str)."""

REASON_TEXT: dict[str, str] = {
    "ref_unread": "REF was not read from the source assembly (resolve_ref left it unread)",
    "source_ref_mismatch": "the source REF disagrees with the source assembly (ref_check mismatch)",
    "source_ref_unstated": "the source REF states no base (ref_check unstated) and an ALT is a base sequence, spelled against that N",
    "allele_not_sequence": (
        "an allele that is not a base sequence and no symbolic SV the lift reads (`*`, `<*>`, `<TRA>`, a symbolic SV beside a sequence ALT), "
        "or an SV or segment whose REF is longer than its one padding base"
    ),
    "split_flagged": "flagged as one record of a pre-split site (RM270); run rejoin_split first",
    "contig_not_in_chain": "a contig the chain does not cover",
    "chain_contig_mismatch": "a contig whose size in the chain differs from the source assembly's",
    "no_chain_block": "no chain block covers it (a gap in the alignment)",
    "split_across_blocks": "its REF span runs off a chain block: its first or last base is in no block (a gap's source-only bases, or past the chain's end)",
    "split_across_chains": "no one chain holds every block-covered piece of its REF span: its ends lie in blocks of different chains (another alignment, target contig or strand), or it crosses a block of another chain",
    "gap_in_change": (
        "its REF span crosses a gap of one chain inside the bases an allele changes, or where an indel in a repeat could equally be placed, "
        "so the change has no one place on the target"
    ),
    "gap_too_large": "its REF span crosses a gap of one chain of GAP_LIMIT bases or more on either side",
    "multiple_targets": "it lies in more than one chain block (overlapping chains), so its target is ambiguous",
    "target_contig_unmatched": "the chain's target contig is not in the target collection, or its size differs",
    "indel_minus_strand": "an indel or SV on a minus-strand block lifted to the target contig's first base, with no base before it to anchor on",
    "ref_mismatch_after_lift": "the target's REF equals none of the row's alleles (an SV's: its padding base differs)",
    "breakend_not_lifted": (
        "a breakend ALT that is none of VCF's four forms or a single breakend (brackets that do not pair, bases other than ACGTN), "
        "a breakend beside another kind of allele, or a breakend whose REF is not one base"
    ),
    "breakend_mate_not_lifted": "a breakend whose mate position has no single chain target (no block, several, or a contig the chain or the target lacks)",
    "breakend_ends_apart": "a breakend whose mate is on its own source contig but lands on another target contig",
    "breakend_ends_reversed": "a breakend whose mate is on its own source contig and on one strand with it, but lands on the other side of it",
    "breakend_distance_changed": "a breakend whose mate is on its own source contig and lands at a distance differing by more than SV_LENGTH_TOLERANCE",
    "sv_end_unlifted": "an SV whose end_pos has no single chain target (no block, several, or a target contig the collection lacks)",
    "sv_ends_apart": "an SV whose POS and end_pos land on different target contigs or strands",
    "sv_ends_reversed": "an SV whose POS and end_pos land in reversed order on the target",
    "sv_length_changed": "an SV whose lifted length differs from the source's by more than SV_LENGTH_TOLERANCE",
}

SV_LENGTH_TOLERANCE = 0.05
"""How far a lifted SV's length (`end_pos - pos`) may differ from the source's, as a fraction of the
source's: kept where `|lifted - source| <= SV_LENGTH_TOLERANCE * source`. liftoverSV's rule and default
(`-p/--percent 0.05`, "difference max between both SVLENs < 5%"); Picard LiftoverVcf has none (the
whole POS..END interval must lift, `LIFTOVER_MIN_MATCH` 1.0), CrossMap rewrites END from the REF span.
A point SV (`<INS>`, length 0) must stay a point, which lifting one position always gives."""

GAP_LIMIT = 50
"""A sequence row whose REF span crosses a gap of one chain lifts only where every gap it crosses has
fewer than `GAP_LIMIT` bases on each side (`dt`, the source's own bases, and `dq`, the target's), else
`gap_too_large` (step 2b). The rule never puts a gap inside the change, so no size is forced by the
spelling: the gap's `dq` target bases go into every allele as the target's reference context, and its
`dt` source bases drop out of every allele. What the size bounds is how much sequence a lifted record
states that its source record never covered. No tool or chain fixes that number, so this is the
conservative reading, at the indel/SV boundary of VCF practice (a change of 50 bases or more is an SV):
a gap below it is an indel-sized difference between the assemblies, and the lifted record stays a small
variant. Measured 2026-11-03 on the chains in hand: one-sided gaps under 50 bases are 95% (`dt`) and 96%
(`dq`) of UCSC hg19ToHg38's, 81% and 94% of hg38ToHg19's, 96% and 98% of Ensembl GRCh38_to_GRCh37's,
97% and 98% of hs1ToHg38's; two-sided gaps 20-68%. bcftools +liftover's `--max-indel-inc` (250) bounds
another thing, how far it lets an indel's anchors drift while it realigns the indel through a gap."""

PARENT_COLUMNS_KEPT = ("ids", "phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities", "gl", "gp", "dosage", "imputation_r2", "imputation_empirical_r2", "evidence_basis", "gq", "hq", "dp", "ad", "qual", "filters", "sample_filters", "split_signal")
"""`calls` columns a lifted row carries from its parent as they were (phase sets and dosage aside, below).
An imputation r² is carried on a swapped row too: the squared correlation of a dosage and the truth is
the same whichever allele the dosage counts."""

REPORT = pl.Schema(
    {
        **CALLS,
        "lift_chrom": pl.String,
        "lift_pos": pl.UInt32,
        "lift_strand": pl.String,
        "lift_swapped": pl.Boolean,
        "lift_reason": pl.String,
    }
)
"""`Lifted.report`: the parent's `calls` row, then where it went (`lift_chrom`, 1-based `lift_pos`,
the block's strand `+`/`-`, whether REF and an ALT swapped) or why not (`lift_reason`, null when lifted)."""

BLOCK_REPORT = pl.Schema({"lift_reason": pl.String, "pieces": pl.UInt64, "bases": pl.UInt64})
"""What a lift left behind of the `blocks` frame: per reason, the pieces and the bases in them."""

_SEQUENCE = r"^[ACGTN]+$"
_SV_ALLELE = rf"^<({'|'.join((*_SV_LENGTH_TYPES, 'INS'))})(:[^>]*)?>$"
"""A symbolic SV allele the lift reads: the loader's length types (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`,
`<CNn>`, subtypes included), whose extent is `end_pos`, and `<INS>`, a point at its anchor."""
_BREAKEND = r"[\[\]]|^\.[ACGTN]+$|^[ACGTN]+\.$"
"""A breakend allele (`N[2:123[`, `]2:123]N`) or a single breakend (`.N`, `N.`), as the row's ALT shows it."""
_LIFTABLE_BREAKEND = (
    r"^[ACGTN]+\[[^\[\]]+:\d+\[$|^[ACGTN]+\][^\[\]]+:\d+\]$|^\[[^\[\]]+:\d+\[[ACGTN]+$|^\][^\[\]]+:\d+\][ACGTN]+$|^\.[ACGTN]+$|^[ACGTN]+\.$"
)
"""A breakend ALT the lift reads: one of the four forms with brackets that pair, or a single breakend,
its bases ACGTN (`breakends.parse_breakend` reads the same)."""
_COMPLEMENT = (["A", "C", "G", "T"], ["T", "G", "C", "A"])
_BASE_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_FAR = 1 << 40
"""An end past every contig, closing the last gap segment of a chain's source contig."""


class Lifted(BaseModel):
    """What `liftover` returns: the derived facet, and the fate of every parent row."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    loaded: Loaded
    """The lifted facet (`relation = liftover_to`)."""
    report: pl.LazyFrame
    """One row per parent `calls` row, in the parent's order, with `REPORT`'s columns: CrossMap's
    `.unmap` (the rows with a `lift_reason`) and bcftools' `--write-src` (the rows without) in one frame."""
    reasons: dict[str, int]
    """Unlifted rows per reason, sorted by reason."""
    block_reasons: dict[str, int] = {}
    """Bases of the parent's `blocks` left behind, per reason (the `UNLIFTED_REASONS` a position can
    have), sorted by reason. A block partly lifted counts only its unlifted part here."""


@dataclass(frozen=True)
class Chain:
    """A chain file read into blocks, with each side's contig sizes."""

    name: str
    """The file's name (not its path: a path is a label, and can carry a person's name)."""
    sha256: str
    blocks: pl.DataFrame
    """One row per ungapped block: `t_chrom` (canonical), `t_start`, `t_end` (0-based half-open on the
    source), `q_chrom` (the target's name as the chain writes it), `q_start` (0-based on the chain's
    target strand), `q_size`, `q_minus`, `chain` (the chain's index in the file)."""
    source_sizes: dict[str, set[int]]
    """Source contig sizes by canonical name; more than one size means the chain contradicts itself."""
    target_sizes: dict[str, set[int]]
    """Target contig sizes by the chain's own names."""
    segments: dict[str, pl.DataFrame] = field(default_factory=dict, compare=False, repr=False)
    """`_segments` per source contig, filled as a lift asks for them: breakend mates lie on any contig,
    so one contig's axis is read once per lift rather than once per contig whose breakends point at it."""


def _contig_segments(chain: Chain, chrom: str) -> pl.DataFrame:
    """The chain's axis on source contig `chrom` (`_segments`), computed once per chain."""
    if chrom not in chain.segments:
        chain.segments[chrom] = _segments(chain.blocks.filter(pl.col("t_chrom") == chrom))
    return chain.segments[chrom]


def read_chain(path: Path | str) -> Chain:
    """Read a UCSC chain file (plain or gzipped, by magic bytes) into blocks.

    The format: a header `chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd
    id`, then lines `size dt dq` and a last `size`. `t` is the source, `q` the target; `tStrand` is
    always `+`, and on a `-` chain the q coordinates count on the reverse strand.
    """
    path = Path(path)
    data = path.read_bytes()
    text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
    columns: dict[str, list] = {name: [] for name in ("t_chrom", "t_start", "t_end", "q_chrom", "q_start", "q_size", "q_minus", "chain")}
    source_sizes: dict[str, set[int]] = {}
    target_sizes: dict[str, set[int]] = {}
    chain = -1
    t = q = 0
    t_chrom = q_chrom = ""
    q_size = 0
    q_minus = False
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "chain":
            if len(parts) < 12 or parts[4] != "+":
                raise ValueError(f"{path.name}: not a chain header (or a source strand other than +): {line[:80]}")
            chain += 1
            t_chrom, q_chrom = canonical_chrom(parts[2]), parts[7]
            source_sizes.setdefault(t_chrom, set()).add(int(parts[3]))
            target_sizes.setdefault(q_chrom, set()).add(int(parts[8]))
            q_size, q_minus = int(parts[8]), parts[9] == "-"
            t, q = int(parts[5]), int(parts[10])
            continue
        size = int(parts[0])
        for name, value in zip(columns, (t_chrom, t, t + size, q_chrom, q, q_size, q_minus, chain)):
            columns[name].append(value)
        if len(parts) == 3:
            t += size + int(parts[1])
            q += size + int(parts[2])
    blocks = pl.DataFrame(
        columns,
        schema={"t_chrom": pl.String, "t_start": pl.Int64, "t_end": pl.Int64, "q_chrom": pl.String, "q_start": pl.Int64, "q_size": pl.Int64, "q_minus": pl.Boolean, "chain": pl.Int64},
    ).with_row_index("block")
    return Chain(name=path.name, sha256=hashlib.sha256(data).hexdigest(), blocks=blocks, source_sizes=source_sizes, target_sizes=target_sizes)


def _segments(blocks: pl.DataFrame) -> pl.DataFrame:
    """One source contig's axis cut where the set of targets changes: `seg_start`, how many distinct
    targets cover the segment (`depth`), a block giving it when exactly one does, and running counts of
    covered and multiply-covered segments, so a span over several segments is judged by two lookups.
    That block stands for its mapping only: where several blocks repeat the mapping it is any one of
    them, of any chain, and need not cover the whole segment, so a reader takes its offset, never its
    chain or extent (`_gap_plan` reads the chains from every block of the mapping).

    A target is a mapping, not a block: blocks with the same target contig, strand and offset
    (`q_start - t_start`) send every shared base to the same place, so they count once. Ensembl's
    GRCh38_to_GRCh37 chain lists every block twice, in two chains (F14); counting blocks would make
    the whole genome ambiguous.
    """
    keys = list(zip(blocks["q_chrom"], blocks["q_minus"], (blocks["q_start"] - blocks["t_start"])))
    ids = blocks["block"].to_list()
    events = sorted([(s, 1, i) for i, s in enumerate(blocks["t_start"])] + [(e, -1, i) for i, e in enumerate(blocks["t_end"])])
    starts, depths, singles = [0], [0], [None]
    active: dict[int, tuple] = {}
    index = 0
    while index < len(events):
        at = events[index][0]
        while index < len(events) and events[index][0] == at:
            _, step, i = events[index]
            if step == 1:
                active[i] = keys[i]
            else:
                active.pop(i, None)
            index += 1
        targets = {key: i for i, key in sorted(active.items())}
        depth = len(targets)
        single = next(iter(targets.values())) if depth == 1 else None
        if at == starts[-1]:
            depths[-1], singles[-1] = depth, single
            continue
        # A cut where the one target stays the same (overlapping blocks of one mapping) is no edge.
        previous = singles[-1]
        if depth == depths[-1] and (depth == 0 or (depth == 1 and keys[previous] == keys[single])):
            continue
        starts.append(at)
        depths.append(depth)
        singles.append(single)
    blocks_of = [None if i is None else ids[i] for i in singles]
    frame = pl.DataFrame({"seg_start": starts, "depth": depths, "seg_block": blocks_of}, schema={"seg_start": pl.Int64, "depth": pl.Int64, "seg_block": pl.UInt32})
    return frame.with_row_index("seg").with_columns(
        (pl.col("depth") >= 1).cast(pl.Int64).cum_sum().alias("covered_cum"),
        (pl.col("depth") >= 2).cast(pl.Int64).cum_sum().alias("multi_cum"),
    )


def _to_pos(start: pl.Expr) -> pl.Expr:
    """A lifted 0-based start as a 1-based `pos` (charter 13: the way back from `reference._span`)."""
    return (start + 1).cast(pl.UInt32)


def _complement(allele: pl.Expr, minus: pl.Expr) -> pl.Expr:
    """An allele on the target strand: uppercase, and reverse-complemented on a minus-strand block."""
    upper = allele.str.to_uppercase()
    return pl.when(minus).then(upper.str.reverse().str.replace_many(*_COMPLEMENT)).otherwise(upper)


def _complement_list(alleles: pl.Expr, minus: pl.Expr) -> pl.Expr:
    upper = alleles.list.eval(pl.element().str.to_uppercase())
    flipped = alleles.list.eval(pl.element().str.to_uppercase().str.reverse().str.replace_many(*_COMPLEMENT))
    return pl.when(minus).then(flipped).otherwise(upper)


def liftover(
    loaded: Loaded,
    chain: Path | str | Chain,
    reference: "ReferenceStore | Path | str | None" = None,
    *,
    target: str | None = None,
) -> Lifted:
    """The facet lifted to another assembly through `chain`, as a derived facet, and a report of every row.

    `loaded` holds one facet whose REF `resolve_ref` read; `reference` is a store (or its directory;
    default the store at `reference_dir()`) holding that collection and the target. `target` names the
    target collection's digest; without it, the one collection whose contig sizes match the chain's
    target side on every primary contig is used, and anything else is refused. The input is not changed.
    """
    cache.sweep_stale(cache.LIFT_SPILL, _VERSIONS)
    parent = loaded.facets.collect()
    if parent.height != 1:
        raise ValueError(f"liftover takes one facet; this Loaded has {parent.height}")
    facet = parent.row(0, named=True)
    source_digest = facet["ref_seqcol_digest"]
    if source_digest is None:
        raise ValueError(f"facet {facet['facet_id']}: REF has not been read; run resolve_ref first, then liftover")
    store, _ = _store_for(reference)
    if source_digest not in store.collections():
        raise ReferenceMismatchError(f"facet {facet['facet_id']}: the store does not hold collection {source_digest}, which REF was read from")
    chain = chain if isinstance(chain, Chain) else read_chain(chain)
    source_ok = _check_source(chain, store, source_digest, facet["facet_id"])
    target_digest = _choose_target(chain, store, target, facet["facet_id"])
    target_contigs = _target_contigs(chain, store, target_digest)

    facet_id = _derived_id(facet["facet_id"], chain.sha256, target_digest)
    entry = cache_dir() / "frames" / f"{facet_id}-{SCHEMA_VERSION}"
    with cache.spilling(entry, cache.LIFT_SPILL, _VERSIONS, cache.derived_origin(facet, facet_id), _GROUPS) as spill:
        # Streamed to parquet once, so each group of contigs below is a filtered scan; its bytes key the result,
        # with what else the warnings name, so a run on the same input reuses what an earlier one wrote.
        staged = spill.stage / "input.parquet"
        loaded.calls.sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
        staged_blocks = spill.stage / "input-blocks.parquet"
        loaded.blocks.sink_parquet(staged_blocks, row_group_size=SPILL_ROW_GROUP)
        named = {"chain_name": chain.name, "labels": [store.assembly_label(source_digest), store.assembly_label(target_digest)]}
        key = spill.key([staged, staged_blocks], {**facet, "_lift": named})
        done = spill.find(key)
        if done is None:
            done = _lift_into(spill, key, staged, staged_blocks, chain, source_ok, store, target_contigs, facet_id, source_digest, target_digest)
    calls_parts = done.files("calls")
    report_parts = done.files("report")
    blocks_parts = done.files("blocks")
    block_report_parts = done.files("block_report")
    calls = pl.scan_parquet(calls_parts) if calls_parts else empty(CALLS)
    report = pl.scan_parquet(report_parts) if report_parts else empty(REPORT)
    lifted_blocks = pl.scan_parquet(blocks_parts) if blocks_parts else empty(BLOCKS)
    # Read back from the reports, so a reused generation gives the same counts as the run that wrote it.
    found = dict(report.filter(pl.col("lift_reason").is_not_null()).group_by("lift_reason").len().collect().iter_rows())
    reasons = {reason: found[reason] for reason in sorted(found)}
    block_found = dict(pl.scan_parquet(block_report_parts).select("lift_reason", "bases").collect().iter_rows()) if block_report_parts else {}
    block_reasons = {reason: block_found[reason] for reason in sorted(block_found)}
    warnings = done.warnings
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)

    present = pl.concat([calls.select("chrom"), lifted_blocks.select("chrom")]).unique().collect()["chrom"].to_list()
    names = {chrom: contig.name for chrom, contig in target_contigs.items()}
    # A lifted breakend names its mate's contig inside ALT, as the target collection names it: those
    # contigs are the rows' too, though no row may stand on one. The lift wrote them down beside the rows.
    mate_parts = done.files("mate_contigs")
    mates = pl.read_parquet(mate_parts)["chrom"].to_list() if mate_parts else []
    present = sorted({*present, *mates})
    label = store.assembly_label(target_digest)
    chain_input = [{"role": "chain", "name": chain.name, "sha256": chain.sha256}]
    facets = parent.lazy().with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.lit(facet["facet_id"]).alias("derived_from"),
        pl.lit(RELATION).alias("relation"),
        pl.lit(label, dtype=pl.String).alias("assembly_label"),
        pl.lit(AssemblyLabelBasis.INFERRED.value if label else None, dtype=FACETS["assembly_label_basis"]).alias("assembly_label_basis"),
        pl.lit(None, dtype=pl.String).alias("seqcol_digest"),
        pl.lit(None, dtype=pl.String).alias("sorted_name_length_pairs_digest"),
        pl.lit(target_digest).alias("ref_seqcol_digest"),
        # The parent's header named the source assembly's file; the target is `ref_seqcol_digest`. The
        # aligner stays: the reads were aligned to the source, and that is still what happened.
        pl.lit(None, dtype=FACETS["reference"]).alias("reference"),
        pl.lit(AbsenceMeans.UNKNOWN.value, dtype=FACETS["absence_means"]).alias("absence_means"),
        pl.lit(None, dtype=FACETS["header_contigs"]).alias("header_contigs"),
        pl.lit(
            sorted(({"source": names.get(chrom, chrom), "chrom": chrom} for chrom in present), key=lambda entry: entry["source"]),
            dtype=FACETS["contig_map"],
        ).alias("contig_map"),
        pl.lit(chain_input, dtype=FACETS["derivation_inputs"]).alias("derivation_inputs"),
    )
    derived = Loaded(
        facets=facets.select(FACETS.names()),
        calls=calls,
        blocks=lifted_blocks,
        links=derived_links(
            loaded.links, facet["facet_id"], facet_id, RELATION, step="liftover", step_version=str(LIFT_VERSION), derivation_inputs=chain_input
        ),
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )
    return Lifted(loaded=derived, report=report, reasons=reasons, block_reasons=block_reasons)


LIFT_LAYOUT = 1
"""How a lift spill's files hold the rows, a key in its manifest and in its generation's key: 1 is the
lifted `calls` and `blocks` in canonical order in parts of whole target contigs (`loading.canonical_parts`)
and the report in parts of whole source contigs, both grouped up to `loading.SPILL_PART_ROWS` rows, so
contigs with a few rows share a file (F24). A spill without the key held one file per contig (942
report files for a decoy-heavy GRCh38 Manta file). Its generation key lacks the key, so no run of this
version finds it, and the sweep retires such an entry once it is unused past the grace period. Not a
`LIFT_VERSION` bump: the rows are the same, and that version is part of the derived `facet_id`. Not
`cache.SPILL_LAYOUT` either, which every other step's spill carries too."""

_VERSIONS = {"SCHEMA_VERSION": SCHEMA_VERSION, "LIFT_VERSION": LIFT_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT, "LIFT_LAYOUT": LIFT_LAYOUT}
"""The version keys of a lift spill, in its cache manifest."""

_GROUPS = ("calls", "report", "blocks", "block_report", "mate_contigs")
"""The file groups of a lift's generation, which its key covers with `_VERSIONS` (`cache.Spill.key`)."""


def _lift_into(
    spill: cache.Spill,
    key: str,
    staged: Path,
    staged_blocks: Path,
    chain: Chain,
    source_ok: dict[str, bool],
    store: ReferenceStore,
    target_contigs: dict[str, Contig],
    facet_id: str,
    source_digest: str,
    target_digest: str,
) -> cache.Committed:
    """Lift the staged `calls` and `blocks` into the spill's stage, and commit it as the generation `key`.

    The source rows are read a group of whole contigs at a time (`canonical_batches`), each group one
    filtered scan, its report one file; the lifted rows are then written in canonical order in parts of
    whole target contigs (`canonical_parts`). Within a group a contig the chain covers is lifted on its
    own, through its own chain blocks; the rows of the group's other contigs take their reason together
    (`_lift_group`), with no pass of their own (F24)."""
    source = pl.scan_parquet(staged)
    sizes = dict(source.group_by("chrom").len().drop_nulls("chrom").collect().iter_rows())
    counts: Counter[str] = Counter()
    lifted_parts: list[Path] = []
    report_parts: list[Path] = []
    for index, batch in enumerate(canonical_batches(sizes)):
        report, lifted = _lift_group(_by_contig(source, batch), batch, chain, source_ok, store, target_contigs, facet_id, counts)
        report_parts.append(spill.stage / f"report-{index:05d}.parquet")
        report.write_parquet(report_parts[-1], row_group_size=SPILL_ROW_GROUP)
        if lifted.height:
            lifted_parts.append(spill.stage / f"lifted-{index:05d}.parquet")
            lifted.write_parquet(lifted_parts[-1], row_group_size=SPILL_ROW_GROUP)
    staged.unlink()

    calls_parts = _target_parts(lifted_parts, spill.stage / "lifted.parquet", "calls")
    blocks_parts, block_report = _lift_blocks_into(spill, staged_blocks, chain, source_ok, target_contigs, facet_id, counts)
    block_report_parts: list[Path] = []
    if block_report.height:
        block_report_parts.append(spill.stage / "block-report.parquet")
        block_report.write_parquet(block_report_parts[-1])
    reasons = {reason: counts[f"reason:{reason}"] for reason in sorted(UNLIFTED_REASONS) if counts[f"reason:{reason}"]}
    warnings = _warnings(counts, reasons, block_report, chain, store, source_digest, target_digest)
    # The canonical target contigs lifted breakends name as mates, for the facet's `contig_map`.
    mate_contigs = sorted(key.removeprefix("bnd_mate_chrom:") for key in counts if key.startswith("bnd_mate_chrom:"))
    mate_parts: list[Path] = []
    if mate_contigs:
        mate_parts.append(spill.stage / "mate-contigs.parquet")
        pl.DataFrame({"chrom": mate_contigs}, schema={"chrom": pl.String}).write_parquet(mate_parts[-1])
    groups = {"calls": calls_parts, "report": report_parts, "blocks": blocks_parts, "block_report": block_report_parts, "mate_contigs": mate_parts}
    return spill.commit(key, groups, warnings)


def _by_contig(source: pl.LazyFrame, batch: list[str]) -> dict[str, pl.DataFrame]:
    """The rows of a group of contigs, read in one scan, by contig, each in the source's order. The
    callers take each contig's rows out as they go, so the group is not held twice."""
    rows = source.filter(pl.col("chrom").is_in(batch)).collect()
    return {chrom: frame for (chrom,), frame in rows.partition_by("chrom", as_dict=True, maintain_order=True).items()}


def _lift_group(
    rows: dict[str, pl.DataFrame],
    batch: list[str],
    chain: Chain,
    source_ok: dict[str, bool],
    store: ReferenceStore,
    targets: dict[str, Contig],
    facet_id: str,
    counts: Counter[str],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The rows of a group of whole source contigs (`batch`, in canonical order; `rows` by contig,
    taken out as they are lifted): (the report rows,
    contig by contig in canonical order and in the parent's order within each, as `_lift_contig` gives
    them one contig at a time; the lifted rows, contig by contig likewise).

    A contig the chain covers at the source assembly's size is lifted on its own, through its own chain
    blocks, its breakends' mates placed for the whole group at once (`_mate_places`). The rows of the
    group's other contigs are passed to `_lift_contig` together, once for the contigs the chain does
    not name and once for those it names at another size: each of those rows takes its reason before
    any placement (`contig_not_in_chain`, `chain_contig_mismatch`, or the row's own, which comes
    first), so none is placed and the contig named is only the fitness's.
    """
    reports: dict[str, pl.DataFrame] = {}
    lifted: list[pl.DataFrame] = []
    for fit in (None, False):
        unfit = [chrom for chrom in batch if source_ok.get(chrom) is fit]
        if unfit:
            report = _lift_contig(pl.concat([rows.pop(chrom) for chrom in unfit], how="vertical"), unfit[0], chain, source_ok, store, targets, facet_id, counts, {})[0]
            # Contig by contig as they went in, each in the parent's order.
            reports |= {chrom: frame for (chrom,), frame in report.partition_by("chrom", as_dict=True, maintain_order=True).items()}
    covered = [chrom for chrom in batch if source_ok.get(chrom)]
    mates = _mate_places([rows[chrom] for chrom in covered], chain, source_ok, targets)
    for chrom in covered:
        reports[chrom], rows_lifted = _lift_contig(rows.pop(chrom), chrom, chain, source_ok, store, targets, facet_id, counts, mates)
        lifted.append(rows_lifted)
    report = pl.concat([reports[chrom] for chrom in batch], how="vertical")
    return report, pl.concat(lifted, how="vertical") if lifted else pl.DataFrame()


def _target_parts(pieces: list[Path], staged: Path, name: str) -> list[Path]:
    """The lifted rows in `pieces` (in source order), written as `<name>-<n>.parquet` parts in canonical
    order, of whole target contigs (`canonical_parts`); none when nothing lifted."""
    if not pieces:
        return []
    pl.scan_parquet(pieces).sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
    for piece in pieces:
        piece.unlink()
    parts = canonical_parts(staged, name)
    staged.unlink()
    return parts


def _lift_blocks_into(
    spill: cache.Spill, staged: Path, chain: Chain, source_ok: dict[str, bool], targets: dict[str, Contig], facet_id: str, counts: Counter[str]
) -> tuple[list[Path], pl.DataFrame]:
    """Lift the staged `blocks` a group of whole source contigs at a time, as `_lift_into` lifts the
    calls: the lifted blocks' files, in canonical order in parts of whole target contigs, and the
    `BLOCK_REPORT` of what stayed behind."""
    source = pl.scan_parquet(staged)
    sizes = dict(source.group_by("chrom").len().drop_nulls("chrom").collect().iter_rows())
    pieces: list[Path] = []
    left: list[pl.DataFrame] = []
    for index, batch in enumerate(canonical_batches(sizes)):
        rows = _by_contig(source, batch)
        # The blocks on contigs the chain does not cover, or covers at another size, all stay behind.
        for fit in (None, False):
            unfit = [chrom for chrom in batch if source_ok.get(chrom) is fit]
            if unfit:
                left.append(_lift_block_contig(pl.concat([rows.pop(chrom) for chrom in unfit], how="vertical"), unfit[0], chain, fit, targets, facet_id, counts)[1])
        lifted: list[pl.DataFrame] = []
        for chrom in batch:
            if source_ok.get(chrom):
                pieces_of, behind = _lift_block_contig(rows.pop(chrom), chrom, chain, True, targets, facet_id, counts)
                left.append(behind)
                lifted.append(pieces_of)
        if any(frame.height for frame in lifted):
            pieces.append(spill.stage / f"lifted-blocks-{index:05d}.parquet")
            pl.concat(lifted, how="vertical").write_parquet(pieces[-1], row_group_size=SPILL_ROW_GROUP)
    staged.unlink()
    parts = _target_parts(pieces, spill.stage / "lifted-blocks.parquet", "blocks")
    for part in parts:
        # Two source regions the chain sends to one place give overlapping blocks; both are kept (charter 8).
        overlaps = pl.scan_parquet(part).select((pl.col("pos") <= pl.col("end_pos").cum_max().shift(1).over("chrom")).sum()).collect().item()
        counts["block_overlaps"] += int(overlaps or 0)
    behind = pl.concat(left) if left else pl.DataFrame(schema={"lift_reason": pl.String, "_length": pl.Int64})
    report = (
        behind.group_by("lift_reason")
        .agg(pl.len().cast(pl.UInt64).alias("pieces"), pl.col("_length").sum().cast(pl.UInt64).alias("bases"))
        .sort("lift_reason")
        .select(BLOCK_REPORT.names())
    )
    return parts, report


def _lift_block_contig(
    rows: pl.DataFrame, chrom: str, chain: Chain, contig_ok: bool | None, targets: dict[str, Contig], facet_id: str, counts: Counter[str]
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """One source contig's blocks: (the lifted pieces as `blocks` rows, the pieces left behind as
    `lift_reason` and `_length` in bases).

    Each block is cut at every segment of the chain's axis it overlaps (`_segments`: the axis is cut
    wherever the set of targets changes, so a gap, a block edge or a change of offset is a cut). A
    piece in one segment with one target lifts by that chain block's offset; the rest stays behind.
    """
    counts["blocks"] += rows.height
    start, end = _span(pl.col("pos"), pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64) + 1)
    work = rows.with_row_index("_b").with_columns(start.alias("_s0"), end.alias("_e0"))
    chain_blocks = chain.blocks.filter(pl.col("t_chrom") == chrom)
    if not contig_ok or chain_blocks.is_empty():
        reason = "contig_not_in_chain" if contig_ok is None else "chain_contig_mismatch" if contig_ok is False else "no_chain_block"
        counts["blocks_left"] += rows.height
        behind = work.select(pl.lit(reason).alias("lift_reason"), (pl.col("_e0") - pl.col("_s0")).alias("_length"))
        return pl.DataFrame(schema=BLOCKS), behind
    segments = _contig_segments(chain, chrom).with_columns(pl.col("seg_start").shift(-1).fill_null(_FAR).alias("seg_end"))
    first = (
        work.select("_b", "_s0")
        .sort("_s0")
        .join_asof(segments.select("seg_start", pl.col("seg").alias("_seg_s")), left_on="_s0", right_on="seg_start", strategy="backward")
    )
    last = (
        work.select("_b", (pl.col("_e0") - 1).alias("_last"))
        .sort("_last")
        .join_asof(segments.select("seg_start", pl.col("seg").alias("_seg_e")), left_on="_last", right_on="seg_start", strategy="backward")
    )
    cut = (
        work.join(first.select("_b", "_seg_s"), on="_b", how="left")
        .join(last.select("_b", "_seg_e"), on="_b", how="left")
        .with_columns(pl.int_ranges(pl.col("_seg_s"), pl.col("_seg_e") + 1, dtype=pl.UInt32).alias("seg"))
        .explode("seg", empty_as_null=False)
        .join(segments.select("seg", "seg_start", "seg_end", "depth", "seg_block"), on="seg", how="left")
        .with_columns(pl.max_horizontal("_s0", "seg_start").alias("_p0"), pl.min_horizontal("_e0", "seg_end").alias("_p1"))
        .join(chain_blocks.select(pl.col("block").alias("seg_block"), "t_start", "q_chrom", "q_start", "q_size", "q_minus"), on="seg_block", how="left")
    )
    names = {name: canonical_chrom(name) for name in cut["q_chrom"].drop_nulls().unique().to_list()}
    target = pl.col("q_chrom").replace_strict(names, default=None, return_dtype=pl.String)
    reason = (
        pl.when(pl.col("depth") == 0)
        .then(pl.lit("no_chain_block"))
        .when(pl.col("depth") >= 2)
        .then(pl.lit("multiple_targets"))
        .when(~target.is_in(list(targets)).fill_null(True))
        .then(pl.lit("target_contig_unmatched"))
    )
    length = pl.col("_p1") - pl.col("_p0")
    offset = pl.col("_p0") - pl.col("t_start")
    q0 = pl.when(pl.col("q_minus")).then(pl.col("q_size") - (pl.col("q_start") + offset + length)).otherwise(pl.col("q_start") + offset)
    cut = cut.with_columns(reason.alias("lift_reason"), length.alias("_length"), target.alias("_q_chrom"), q0.alias("_q0"))
    ok = pl.col("lift_reason").is_null()
    per_block = cut.group_by("_b").agg(ok.sum().alias("lifted"), (~ok).sum().alias("left"))
    counts["blocks_whole"] += per_block.filter((pl.col("lifted") == 1) & (pl.col("left") == 0)).height
    counts["blocks_split"] += per_block.filter((pl.col("lifted") >= 1) & ((pl.col("lifted") > 1) | (pl.col("left") > 0))).height
    counts["blocks_left"] += per_block.filter(pl.col("lifted") == 0).height
    lifted = cut.filter(ok)
    counts["block_pieces"] += lifted.height
    counts["block_minus"] += int(lifted["q_minus"].sum() or 0)
    counts["block_other_contig"] += lifted.filter(pl.col("_q_chrom") != chrom).height
    out = lifted.sort("_b", "_p0").select(
        pl.lit(facet_id).alias("facet_id"),
        pl.col("_q_chrom").alias("chrom"),
        _to_pos(pl.col("_q0")).alias("pos"),
        (pl.col("_q0") + pl.col("_length")).cast(pl.UInt32).alias("end_pos"),
        *[pl.col(name) for name in BLOCKS.names() if name not in ("facet_id", "chrom", "pos", "end_pos")],
    )
    return out.cast(dict(BLOCKS)), cut.filter(~ok).select("lift_reason", "_length")


def _derived_id(parent_id: str, chain_sha256: str, target: str) -> str:
    key = json.dumps(
        {"chain": chain_sha256, "derived_from": parent_id, "operation": RELATION, "reference": target, "version": LIFT_VERSION}, sort_keys=True
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]


def _single(store: ReferenceStore, digest: str) -> dict[str, Contig]:
    """The collection's contigs by canonical name; a name it holds twice is left out."""
    return {chrom: contigs[0] for chrom, contigs in store.contigs_by_chrom(digest).items() if len(contigs) == 1}


def _check_source(chain: Chain, store: ReferenceStore, digest: str, facet_id: str) -> dict[str, bool]:
    """Per source contig the chain covers: whether its size there is the source collection's.

    Raises `ReferenceMismatchError` when a primary contig differs: the chain starts from another assembly.
    """
    contigs = _single(store, digest)
    result: dict[str, bool] = {}
    for chrom, sizes in chain.source_sizes.items():
        ok = chrom in contigs and sizes == {contigs[chrom].length}
        if not ok and chrom in PRIMARY and chrom in contigs:
            raise ReferenceMismatchError(
                f"facet {facet_id}: chain {chain.name} gives contig {chrom} {', '.join(map(str, sorted(sizes)))} bases on its source side, "
                f"the facet's assembly {contigs[chrom].length}: the chain starts from another assembly"
            )
        result[chrom] = ok
    return result


def _primary_fit(chain: Chain, store: ReferenceStore, digest: str) -> tuple[int, int]:
    """(primary contigs whose target size the collection matches, primary contigs where it differs)."""
    contigs = _single(store, digest)
    matched = differing = 0
    for name, sizes in chain.target_sizes.items():
        chrom = canonical_chrom(name)
        if chrom not in PRIMARY or chrom not in contigs:
            continue
        if sizes == {contigs[chrom].length}:
            matched += 1
        else:
            differing += 1
    return matched, differing


def _choose_target(chain: Chain, store: ReferenceStore, target: str | None, facet_id: str) -> str:
    if target is not None:
        if target not in store.collections():
            raise ReferenceMismatchError(f"facet {facet_id}: the store does not hold target collection {target}")
        matched, differing = _primary_fit(chain, store, target)
        if differing or not matched:
            raise ReferenceMismatchError(
                f"facet {facet_id}: chain {chain.name}'s target sizes do not fit collection {target} "
                f"({matched} primary contig(s) match, {differing} differ)"
            )
        return target
    fits = []
    for digest in store.collections():
        matched, differing = _primary_fit(chain, store, digest)
        if matched and not differing:
            fits.append(digest)
    if len(fits) != 1:
        found = ", ".join(fits) if fits else "none"
        raise ReferenceMismatchError(
            f"facet {facet_id}: {len(fits)} collection(s) in the store fit chain {chain.name}'s target side ({found}); pass target= to choose"
        )
    return fits[0]


def _target_contigs(chain: Chain, store: ReferenceStore, digest: str) -> dict[str, Contig]:
    """The target collection's contigs by canonical name, kept only where the chain's size for it agrees."""
    contigs = _single(store, digest)
    sizes: dict[str, set[int]] = {}
    for name, found in chain.target_sizes.items():
        sizes.setdefault(canonical_chrom(name), set()).update(found)
    return {chrom: contig for chrom, contig in contigs.items() if sizes.get(chrom) == {contig.length}}


def _lift_contig(
    rows: pl.DataFrame,
    chrom: str,
    chain: Chain,
    source_ok: dict[str, bool],
    store: ReferenceStore,
    targets: dict[str, Contig],
    facet_id: str,
    counts: Counter[str],
    mates: "dict[tuple[str, int], _MateFate]",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """One source contig's rows: (the report rows, the lifted `calls` rows). `source_ok` is every
    source contig's fitness (`_check_source`): a breakend's mate may lie on any of them. `mates` is
    where each breakend mate the rows name lands (`_mate_places`)."""
    contig_ok = source_ok.get(chrom)
    counts["rows"] += rows.height
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    not_sequence = alts.list.eval(~pl.element().str.to_uppercase().str.contains(_SEQUENCE)).list.any()
    breakend = alts.list.eval(pl.element().str.to_uppercase().str.contains(_BREAKEND)).list.any()
    ref_length = pl.col("ref").str.len_chars().cast(pl.Int64)
    liftable_breakend = ((alts.list.len() > 0) & alts.list.eval(pl.element().str.to_uppercase().str.contains(_LIFTABLE_BREAKEND)).list.all() & (ref_length == 1)).fill_null(False)
    symbolic_sv = (alts.list.len() > 0) & alts.list.eval(pl.element().str.contains(_SV_ALLELE)).list.all()
    segment = (alts.list.len() == 0) & (pl.col("end_pos").cast(pl.Int64) > pl.col("pos").cast(pl.Int64) + ref_length - 1)
    extent = (symbolic_sv | segment).fill_null(False)
    reason = (
        pl.when(pl.col("ref").is_null())
        .then(pl.lit("ref_unread"))
        .when((pl.col("ref_check") == RefCheck.MISMATCH.value).fill_null(False))
        .then(pl.lit("source_ref_mismatch"))
        .when(unstated_sequence_alt())
        .then(pl.lit("source_ref_unstated"))
        .when(breakend & ~liftable_breakend)
        .then(pl.lit("breakend_not_lifted"))
        .when((not_sequence & ~symbolic_sv & ~liftable_breakend) | (extent & (ref_length != 1)))
        .then(pl.lit("allele_not_sequence"))
        .when(pl.col("split_signal").is_not_null())
        .then(pl.lit("split_flagged"))
    )
    if contig_ok is None:
        reason = reason.otherwise(pl.lit("contig_not_in_chain"))
    elif not contig_ok:
        reason = reason.otherwise(pl.lit("chain_contig_mismatch"))
    start, end = _span(pl.col("pos"), pl.col("ref").str.len_chars())
    work = rows.with_row_index("_row").with_columns(
        reason.alias("lift_reason"), start.alias("_s0"), end.alias("_e0"), extent.alias("_sv"), liftable_breakend.alias("_bnd")
    )

    blocks = chain.blocks.filter(pl.col("t_chrom") == chrom)
    candidates = work.filter(pl.col("lift_reason").is_null())
    placed_columns = ["_row", "_reason2", "_q_chrom", "_q0", "_minus", "_reanchor", "_t_ref", "_swap", "_new_end", *_BND_PLACED, *_GAP_PLACED]
    if candidates.height and blocks.height:
        segments = _contig_segments(chain, chrom)
        sequence = _place(candidates.filter(~pl.col("_sv") & ~pl.col("_bnd")), blocks, targets, segments)
        # A span crossing a block edge may lie across a gap of one chain (step 2b); the rest read as before.
        across = (pl.col("_reason2") == "split_across_blocks").fill_null(False)
        within = _with_nulls(_read_target(sequence.filter(~across), store, targets), _GAP_PLACED)
        gapped = _place_gap(sequence.filter(across), blocks, segments, store, targets)
        sequence = _with_nulls(
            pl.concat([part.select(placed_columns[:8] + list(_GAP_PLACED)) for part in (within, gapped)], how="vertical").with_columns(
                pl.lit(None, dtype=pl.Int64).alias("_new_end")
            ),
            _BND_PLACED,
        )
        svs = _with_nulls(_with_nulls(_place_sv(candidates.filter(pl.col("_sv")), blocks, segments, store, targets), _BND_PLACED), _GAP_PLACED)
        bnds = _with_nulls(_place_bnd(candidates.filter(pl.col("_bnd")), chrom, blocks, chain, store, targets, counts, mates), _GAP_PLACED)
        placed = pl.concat([part.select(placed_columns) for part in (sequence, svs, bnds)], how="vertical")
        work = work.join(placed, on="_row", how="left", maintain_order="left")
        work = work.with_columns(pl.coalesce("lift_reason", "_reason2").alias("lift_reason"))
    else:
        work = _with_nulls(
            work.with_columns(
                pl.lit(None, dtype=pl.String).alias("_q_chrom"),
                pl.lit(None, dtype=pl.Int64).alias("_q0"),
                pl.lit(None, dtype=pl.Boolean).alias("_minus"),
                pl.lit(None, dtype=pl.Boolean).alias("_reanchor"),
                pl.lit(None, dtype=pl.String).alias("_t_ref"),
                pl.lit(None, dtype=pl.Boolean).alias("_swap"),
                pl.lit(None, dtype=pl.Int64).alias("_new_end"),
            ),
            {**_BND_PLACED, **_GAP_PLACED},
        )
        if candidates.height:
            work = work.with_columns(pl.col("lift_reason").fill_null(pl.lit("no_chain_block")))

    for name, count in work.filter(pl.col("lift_reason").is_not_null()).group_by("lift_reason").len().iter_rows():
        counts[f"reason:{name}"] += count
    ok = pl.col("lift_reason").is_null()
    report = work.select(
        *CALLS.names(),
        pl.when(ok).then(pl.col("_q_chrom")).alias("lift_chrom"),
        pl.when(ok).then(_to_pos(pl.col("_q0"))).alias("lift_pos"),
        pl.when(pl.col("_minus").is_not_null()).then(pl.when(pl.col("_minus")).then(pl.lit("-")).otherwise(pl.lit("+"))).alias("lift_strand"),
        pl.when(ok).then(pl.col("_swap")).alias("lift_swapped"),
        pl.col("lift_reason"),
    )
    lifted = pl.concat(
        [
            _lifted_rows(work.filter(ok & ~pl.col("_sv") & ~pl.col("_bnd")), facet_id, chrom, counts),
            _lifted_sv_rows(work.filter(ok & pl.col("_sv")), facet_id, chrom, counts),
            _lifted_bnd_rows(work.filter(ok & pl.col("_bnd")), facet_id, chrom, counts),
        ],
        how="vertical",
    )
    return report, lifted


_BND_PLACED = {"_bnd_alts": CALLS["alts"], "_bnd_cipos": CALLS["cipos"], "_bnd_ciend": CALLS["ciend"]}
"""What `_place_bnd` adds to a placed breakend row: its rewritten ALTs and intervals (null elsewhere)."""

_GAP_PLACED = {
    "_gap": pl.Boolean,
    "_gap_reanchor": pl.Boolean,
    "_gap_ref": pl.String,
    "_gap_alts": CALLS["alts"],
    "_gap_gt_alleles": CALLS["gt_alleles"],
    "_gap_source_alleles": CALLS["source_alleles"],
}
"""What `_place_gap` adds to a row lifted across a gap of one chain: its alleles already on the target
(the parent's REF carried across, the ALTs, `gt_alleles` and `source_alleles`), and whether it was
re-anchored on a minus strand (null elsewhere)."""


def _with_nulls(frame: pl.DataFrame, columns: dict[str, pl.DataType]) -> pl.DataFrame:
    return frame.with_columns(pl.lit(None, dtype=dtype).alias(name) for name, dtype in columns.items())


def _mirror_pairs(values: list[int | None] | None, flips: list[bool]) -> list[int | None] | None:
    """Confidence-interval pairs (two offsets per ALT), pair `i` mirrored (`[lo, hi]` to `[-hi, -lo]`)
    where `flips[i]` holds (the last flag standing for any pair past the ALTs)."""
    if values is None:
        return None
    out: list[int | None] = []
    for index in range(0, len(values) - 1, 2):
        low, high = values[index], values[index + 1]
        if flips[min(index // 2, len(flips) - 1)]:
            low, high = (None if high is None else -high), (None if low is None else -low)
        out += [low, high]
    return out


def _place_bnd(
    rows: pl.DataFrame,
    chrom: str,
    blocks: pl.DataFrame,
    chain: Chain,
    store: ReferenceStore,
    targets: dict[str, Contig],
    counts: Counter[str],
    mates: "dict[tuple[str, int], _MateFate]",
) -> pl.DataFrame:
    """Each breakend row's target from its own position and every ALT's mate, or why it has none
    (`_reason2`), with the rewritten ALTs (`_bnd_alts`) and intervals.

    POS is placed as a one-base span, as an SNV's is, and the target's base there read: it must be
    the source's REF, complemented on a minus-strand block (the base at POS is part of the joined
    sequence, unlike a symbolic SV's padding base). Each mate is placed the same way through the
    chain blocks of the contig it names, as `mates` holds it (`_mate_places`; `breakend_mate_not_lifted`
    where any has no single target); a mate on a contig of the assembly file (`<ctg1>`) is placed by
    none and stays as written, on a plus strand.
    A mate on the row's own source contig is held to liftoverSV's checks on its two ends: one target
    contig (`breakend_ends_apart`), their order kept on a plus strand and reversed on a minus one when
    both share a strand (`breakend_ends_reversed`), and their distance within `SV_LENGTH_TOLERANCE`
    (`breakend_distance_changed`). The ALT is then `breakends.turned` by both ends' strands, the mate
    named as the target collection names it. Few rows, so the per-ALT work is plain Python.
    """
    if rows.is_empty():
        return _with_nulls(
            rows.select("_row").with_columns(
                pl.lit(None, dtype=pl.String).alias("_reason2"),
                pl.lit(None, dtype=pl.String).alias("_q_chrom"),
                pl.lit(None, dtype=pl.Int64).alias("_q0"),
                pl.lit(None, dtype=pl.Boolean).alias("_minus"),
                pl.lit(False).alias("_reanchor"),
                pl.lit(None, dtype=pl.String).alias("_t_ref"),
                pl.lit(False).alias("_swap"),
                pl.lit(None, dtype=pl.Int64).alias("_new_end"),
            ),
            _BND_PLACED,
        )
    own = _place(rows.select("_row", "_s0", "_e0"), blocks, targets, _contig_segments(chain, chrom))
    bases: dict[int, str] = {}
    readable = own.filter(pl.col("_reason2").is_null())
    for target in sorted(readable["_q_chrom"].unique().to_list()):
        part = readable.filter(pl.col("_q_chrom") == target)
        read = store.fetch(targets[target], [(start, start + 1) for start in part["_q0"].to_list()])
        bases.update(zip(part["_row"].to_list(), read, strict=True))

    # Every ALT's mate on a reference contig, where `_mate_places` put it. A mate on a contig of the
    # assembly file is placed by no chain: it stays where the ALT says, on the strand it was.
    parsed = {row: [parse_breakend(alt) for alt in alts] for row, alts in rows.select("_row", "alts").iter_rows()}
    placed_mates = {
        (row, k): mates[(canonical_chrom(b.mate_contig), b.mate_pos)]
        for row, found in parsed.items()
        for k, b in enumerate(found)
        if b.mate_contig is not None and not b.mate_in_assembly
    }

    out: dict[str, list] = {name: [] for name in ("_row", "_reason2", "_t_ref", "_bnd_alts", "_bnd_cipos", "_bnd_ciend")}
    for (row, pos, ref, cipos, ciend), (_, cause, target, q0, minus) in zip(
        rows.select("_row", "pos", "ref", "cipos", "ciend").iter_rows(), own.select("_row", "_reason2", "_q_chrom", "_q0", "_minus").iter_rows(), strict=True
    ):
        found = parsed[row]
        reason, alts, mate_minus = cause, None, [False] * len(found)
        if reason is None:
            anchor = ref.upper().translate(_BASE_COMPLEMENT) if minus else ref.upper()
            reason = "ref_mismatch_after_lift" if bases.get(row) != anchor else None
        if reason is None:
            failed = [placed_mates[(row, k)][0] for k in range(len(found)) if (row, k) in placed_mates and placed_mates[(row, k)][0] is not None]
            for mate_cause in failed:
                counts[f"bnd_mate:{mate_cause}"] += 1
            reason = "breakend_mate_not_lifted" if failed else None
        if reason is None:
            checks: set[str] = set()
            for k, b in enumerate(found):
                if (row, k) not in placed_mates:
                    continue
                _, mate_target, mate_q0, minus_k = placed_mates[(row, k)]
                mate_minus[k] = bool(minus_k)
                if canonical_chrom(b.mate_contig) != chrom:
                    continue
                source_gap, lifted_gap = b.mate_pos - pos, mate_q0 - q0
                if mate_target != target:
                    checks.add("breakend_ends_apart")
                elif minus_k == minus and (lifted_gap > 0) - (lifted_gap < 0) != ((source_gap > 0) - (source_gap < 0)) * (-1 if minus else 1):
                    checks.add("breakend_ends_reversed")
                elif abs(abs(lifted_gap) - abs(source_gap)) > SV_LENGTH_TOLERANCE * abs(source_gap):
                    checks.add("breakend_distance_changed")
            reason = next((check for check in ("breakend_ends_apart", "breakend_ends_reversed", "breakend_distance_changed") if check in checks), None)
        if reason is None:
            counts["bnd_mate_minus"] += sum(mate_minus)
            alts = []
            for k, b in enumerate(found):
                if b.mate_contig is None:
                    alts.append(turned(b, bool(minus), False, None, None))
                    continue
                if b.mate_in_assembly:
                    counts["bnd_mate_assembly"] += 1
                    alts.append(turned(b, bool(minus), False, b.mate_contig, b.mate_pos))
                    continue
                _, mate_target, mate_q0, _ = placed_mates[(row, k)]
                counts[f"bnd_mate_chrom:{mate_target}"] += 1
                alts.append(turned(b, bool(minus), mate_minus[k], targets[mate_target].name, mate_q0 + 1))
        out["_row"].append(row)
        out["_reason2"].append(reason)
        out["_t_ref"].append(bases.get(row) if reason is None else None)
        out["_bnd_alts"].append(alts)
        out["_bnd_cipos"].append(_mirror_pairs(cipos, [bool(minus)]) if reason is None else None)
        out["_bnd_ciend"].append(_mirror_pairs(ciend, mate_minus) if reason is None else None)
    result = pl.DataFrame(
        out,
        schema={"_row": rows.schema["_row"], "_reason2": pl.String, "_t_ref": pl.String, **_BND_PLACED},
    )
    placed = own.select("_row", "_q_chrom", "_q0", "_minus").join(result, on="_row", how="left", maintain_order="left")
    return placed.with_columns(
        pl.when(pl.col("_reason2").is_null()).then(pl.col(name)).alias(name) for name in ("_q_chrom", "_q0")
    ).with_columns(
        pl.lit(False).alias("_reanchor"),
        pl.lit(False).alias("_swap"),
        (pl.col("_q0") + 1).alias("_new_end"),
    )


_MateFate = tuple[str | None, str | None, int | None, bool | None]
"""Where a breakend mate lands: why it has no single target (None when it has one), the target contig,
the 0-based target position, and whether its chain block is on the minus strand."""


def _mate_places(frames: list[pl.DataFrame], chain: Chain, source_ok: dict[str, bool], targets: dict[str, Contig]) -> dict[tuple[str, int], _MateFate]:
    """Where the mate of every breakend ALT among the rows of `frames` lands, keyed by its canonical
    contig and POS.

    The mates are gathered from every row of a group of source contigs and placed one mate contig at a
    time, each in one pass through its own chain blocks (`_place`, one base), whichever source contigs
    name them: a group costs one placement per mate contig, where it cost one per pair of a source
    contig and a mate contig (F24). A position's fate depends on that position alone, so a mate is
    placed the same in any company; rows that turn out not to be lifted only add positions nobody asks
    for. A mate on a contig the chain does not cover, or covers at another size than the source
    assembly's, has that as its cause.
    """
    bracket = pl.col("alts").list.eval(pl.element().str.contains(r"[\[\]]")).list.any().fill_null(False)
    positions: dict[str, set[int]] = {}
    for alts in (alts for frame in frames for alts in frame.filter(bracket)["alts"]):
        for alt in alts:
            found = parse_breakend(alt) if alt is not None else None
            if found is not None and found.mate_contig is not None and not found.mate_in_assembly:
                positions.setdefault(canonical_chrom(found.mate_contig), set()).add(found.mate_pos)
    fates: dict[tuple[str, int], _MateFate] = {}
    for mate_chrom in sorted(positions):
        at = sorted(positions[mate_chrom])
        ok = source_ok.get(mate_chrom)
        mate_blocks = chain.blocks.filter(pl.col("t_chrom") == mate_chrom)
        if not ok or mate_blocks.is_empty():
            cause = "contig_not_in_chain" if ok is None else "chain_contig_mismatch" if ok is False else "no_chain_block"
            fates.update({(mate_chrom, pos): (cause, None, None, None) for pos in at})
            continue
        s0, e0 = _span(pl.col("_mpos"), pl.lit(1))
        part = pl.DataFrame({"_mpos": at}, schema={"_mpos": pl.Int64}).with_row_index("_row")
        found = _place(part.select("_row", s0.alias("_s0"), e0.alias("_e0")), mate_blocks, targets, _contig_segments(chain, mate_chrom))
        fates.update(zip(((mate_chrom, pos) for pos in at), found.select("_reason2", "_q_chrom", "_q0", "_minus").iter_rows(), strict=True))
    return fates


def _lifted_bnd_rows(work: pl.DataFrame, facet_id: str, source_chrom: str, counts: Counter[str]) -> pl.DataFrame:
    """The lifted `calls` rows of breakends: placed at their own position's image, REF the target's
    base, the ALTs as `_place_bnd` rewrote them, in `gt_alleles` and `source_alleles` too."""
    if work.is_empty():
        return pl.DataFrame(schema=CALLS)
    moved_contig = pl.col("_q_chrom") != source_chrom
    counts["lifted"] += work.height
    counts["bnd_lifted"] += work.height
    counts["bnd_minus"] += int(work["_minus"].sum())
    counts["other_contig"] += int(work.select(moved_contig.sum()).item())
    counts["phase_cleared"] += int(work.select((moved_contig & pl.col("phased").is_not_null()).sum()).item())
    columns: dict[str, list] = {"gt_alleles": [], "source_alleles": []}
    for ref, new_ref, alts, new_alts, gt_alleles, source_alleles in work.select("ref", "_t_ref", "alts", "_bnd_alts", "gt_alleles", "source_alleles").iter_rows():
        moved = {ref.upper(): new_ref, **{old.upper(): new for old, new in zip(alts, new_alts, strict=True)}}
        for name, alleles in (("gt_alleles", gt_alleles), ("source_alleles", source_alleles)):
            columns[name].append(None if alleles is None else [None if allele is None else moved.get(allele.upper(), allele) for allele in alleles])
    out = work.with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.col("_q_chrom").alias("chrom"),
        _to_pos(pl.col("_q0")).alias("pos"),
        pl.col("_new_end").cast(pl.UInt32).alias("end_pos"),
        pl.col("_t_ref").alias("ref"),
        pl.col("_bnd_alts").alias("alts"),
        pl.col("_bnd_cipos").alias("cipos"),
        pl.col("_bnd_ciend").alias("ciend"),
        *[pl.Series(name, values, dtype=CALLS[name]) for name, values in columns.items()],
        pl.lit(None, dtype=CALLS["source_ref"]).alias("source_ref"),
        pl.lit(None, dtype=CALLS["ref_check"]).alias("ref_check"),
        pl.lit(None, dtype=CALLS["source_gt"]).alias("source_gt"),
        pl.lit(None, dtype=CALLS["vrs_allele_ids"]).alias("vrs_allele_ids"),
        pl.lit(None, dtype=CALLS["vrs_errors"]).alias("vrs_errors"),
        *[pl.when(moved_contig).then(None).otherwise(pl.col(name)).alias(name) for name in ("phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities")],
    )
    return out.select(CALLS.names())



def _place_sv(rows: pl.DataFrame, blocks: pl.DataFrame, segments: pl.DataFrame, store: ReferenceStore, targets: dict[str, Contig]) -> pl.DataFrame:
    """Each SV row's target from its two ends, lifted as points, or why it has none (`_reason2`); the
    new 0-based start (`_q0`), the new 1-based `end_pos` (`_new_end`) and the target's REF base.

    A padded SV (a symbolic ALT: POS is the base before the event `[POS+1, END]`) on a minus-strand
    block has its event's image at `[image(END), image(POS) - 1]`, so its new POS is the base before
    that, `image(END) - 1`, and its new END `image(POS) - 1`: the length is kept, and the anchor is the
    target's base before the event, as a re-anchored indel's is. A segment (no ALT, `[POS, END]`
    inclusive, as a block) maps to `[image(END), image(POS)]`.
    """
    if rows.is_empty():
        return rows.with_columns(
            pl.lit(None, dtype=pl.String).alias("_reason2"),
            pl.lit(None, dtype=pl.String).alias("_q_chrom"),
            pl.lit(None, dtype=pl.Int64).alias("_q0"),
            pl.lit(None, dtype=pl.Boolean).alias("_minus"),
            pl.lit(False).alias("_reanchor"),
            pl.lit(None, dtype=pl.String).alias("_t_ref"),
            pl.lit(False).alias("_swap"),
            pl.lit(None, dtype=pl.Int64).alias("_new_end"),
        )
    start_s0, start_e0 = _span(pl.col("pos"), pl.lit(1))
    end_s0, end_e0 = _span(pl.col("end_pos"), pl.lit(1))
    at_pos = _place(rows.select("_row", start_s0.alias("_s0"), start_e0.alias("_e0")), blocks, targets, segments)
    at_end = _place(rows.select("_row", end_s0.alias("_s0"), end_e0.alias("_e0")), blocks, targets, segments).select(
        "_row", pl.col("_reason2").alias("_end_reason"), pl.col("_q_chrom").alias("_end_chrom"), pl.col("_q0").alias("_end_q0"), pl.col("_minus").alias("_end_minus")
    )
    rows = rows.join(at_pos.select("_row", "_reason2", "_q_chrom", pl.col("_q0").alias("_pos_q0"), "_minus"), on="_row", how="left", maintain_order="left")
    rows = rows.join(at_end, on="_row", how="left", maintain_order="left")
    padded = pl.col("alts").list.len() > 0
    minus = pl.col("_minus")
    # 0-based images of POS and END: a, b. New 0-based start and 1-based end per the docstring.
    a, b = pl.col("_pos_q0"), pl.col("_end_q0")
    q0 = pl.when(~minus).then(a).when(padded).then(b - 1).otherwise(b)
    new_end = pl.when(~minus).then(b + 1).when(padded).then(a).otherwise(a + 1)
    source_length = pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64)
    lifted_length = new_end - (q0 + 1)
    reason = (
        pl.when(pl.col("_reason2").is_not_null())
        .then(pl.col("_reason2"))
        .when(pl.col("_end_reason").is_not_null())
        .then(pl.lit("sv_end_unlifted"))
        .when((pl.col("_end_chrom") != pl.col("_q_chrom")) | (pl.col("_end_minus") != minus))
        .then(pl.lit("sv_ends_apart"))
        .when(pl.when(minus).then(b > a).otherwise(b < a))
        .then(pl.lit("sv_ends_reversed"))
        .when((lifted_length - source_length).abs() > SV_LENGTH_TOLERANCE * source_length)
        .then(pl.lit("sv_length_changed"))
        .when(q0 < 0)
        .then(pl.lit("indel_minus_strand"))
    )
    rows = rows.with_columns(reason.alias("_reason2"), q0.alias("_q0"), new_end.alias("_new_end"), pl.lit(False).alias("_reanchor"), pl.lit(False).alias("_swap"))
    rows = rows.with_columns(
        pl.when(pl.col("_reason2").is_null()).then(pl.col(name)).alias(name) for name in ("_q_chrom", "_q0", "_new_end")
    )
    readable = rows.filter(pl.col("_reason2").is_null())
    bases: list[pl.DataFrame] = []
    for chrom in sorted(readable["_q_chrom"].unique().to_list()):
        part = readable.filter(pl.col("_q_chrom") == chrom)
        starts = part["_q0"].to_list()
        read = store.fetch(targets[chrom], [(start, start + 1) for start in starts])
        bases.append(part.select("_row", pl.Series("_t_ref", read, dtype=pl.String)))
    read_frame = pl.concat(bases) if bases else pl.DataFrame(schema={"_row": rows.schema["_row"], "_t_ref": pl.String})
    rows = rows.join(read_frame, on="_row", how="left", maintain_order="left")
    # On a plus-strand block the anchor is the source's own base, so the target's must equal it; on a
    # minus-strand block it is another base (the one before the event's image), read, never compared.
    mismatch = pl.col("_reason2").is_null() & ~minus & (pl.col("_t_ref") != pl.col("ref").str.to_uppercase())
    return rows.with_columns(pl.when(mismatch).then(pl.lit("ref_mismatch_after_lift")).otherwise(pl.col("_reason2")).alias("_reason2"))


def _mirrored(values: list[int | None] | None) -> list[int | None] | None:
    """Confidence-interval pairs (two offsets per ALT) about a position whose image runs the other way:
    each `[lo, hi]` becomes `[-hi, -lo]`."""
    if values is None:
        return None
    out: list[int | None] = []
    for index in range(0, len(values) - 1, 2):
        low, high = values[index], values[index + 1]
        out += [None if high is None else -high, None if low is None else -low]
    return out


def _lifted_sv_rows(work: pl.DataFrame, facet_id: str, source_chrom: str, counts: Counter[str]) -> pl.DataFrame:
    """The lifted `calls` rows of SVs and segments: placed by both ends, REF read at the new POS, the
    ALTs as written. Few rows, so the per-row lists are rebuilt in plain Python."""
    if work.is_empty():
        return pl.DataFrame(schema=CALLS)
    moved_contig = pl.col("_q_chrom") != source_chrom
    counts["lifted"] += work.height
    counts["sv_lifted"] += work.height
    counts["sv_minus"] += int(work["_minus"].sum())
    counts["other_contig"] += int(work.select(moved_contig.sum()).item())
    counts["phase_cleared"] += int(work.select((moved_contig & pl.col("phased").is_not_null()).sum()).item())
    columns: dict[str, list] = {name: [] for name in ("gt_alleles", "source_alleles", "cipos", "ciend")}
    for old, new, minus, gt_alleles, source_alleles, cipos, ciend in work.select(
        pl.col("ref").str.to_uppercase(), "_t_ref", "_minus", "gt_alleles", "source_alleles", "cipos", "ciend"
    ).iter_rows():
        # Every allele that is not one of the ALTs is the padding base REF, now the target's.
        for name, alleles in (("gt_alleles", gt_alleles), ("source_alleles", source_alleles)):
            columns[name].append(None if alleles is None else [new if allele is not None and allele.upper() == old else allele for allele in alleles])
        if minus:
            # The new POS is the image of the old END and the reverse: a missing CIEND is CIPOS (VCF 4.5).
            columns["cipos"].append(_mirrored(ciend if ciend is not None else cipos))
            # A source POS left precise (no CIPOS) beside an imprecise END becomes the new END, and a null
            # there would read as the new CIPOS: it is written precise, `[0, 0]` per ALT.
            columns["ciend"].append(_mirrored(cipos) if cipos is not None or ciend is None else [0] * len(ciend))
        else:
            columns["cipos"].append(cipos)
            columns["ciend"].append(ciend)
    out = work.with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.col("_q_chrom").alias("chrom"),
        _to_pos(pl.col("_q0")).alias("pos"),
        pl.col("_new_end").cast(pl.UInt32).alias("end_pos"),
        pl.col("_t_ref").alias("ref"),
        *[pl.Series(name, values, dtype=CALLS[name]) for name, values in columns.items()],
        pl.lit(None, dtype=CALLS["source_ref"]).alias("source_ref"),
        pl.lit(None, dtype=CALLS["ref_check"]).alias("ref_check"),
        pl.lit(None, dtype=CALLS["source_gt"]).alias("source_gt"),
        pl.lit(None, dtype=CALLS["vrs_allele_ids"]).alias("vrs_allele_ids"),
        pl.lit(None, dtype=CALLS["vrs_errors"]).alias("vrs_errors"),
        *[pl.when(moved_contig).then(None).otherwise(pl.col(name)).alias(name) for name in ("phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities")],
    )
    return out.select(CALLS.names())


def _place(rows: pl.DataFrame, blocks: pl.DataFrame, targets: dict[str, Contig], segments: pl.DataFrame) -> pl.DataFrame:
    """Each row's block, target contig and 0-based target start, or why it has none (`_reason2`).
    `segments` is `_segments(blocks)`, the contig's axis the chain holds once (`_contig_segments`)."""
    by_start = rows.select("_row", "_s0", "_e0").sort("_s0")
    first = by_start.join_asof(
        segments.rename({"seg": "seg_s", "depth": "depth_s", "seg_block": "block_s", "covered_cum": "cov_s", "multi_cum": "multi_s"}),
        left_on="_s0",
        right_on="seg_start",
        strategy="backward",
    )
    last = (
        rows.select("_row", (pl.col("_e0") - 1).alias("_last"))
        .sort("_last")
        .join_asof(segments.select(pl.col("seg_start"), pl.col("seg").alias("seg_e"), pl.col("covered_cum").alias("cov_e"), pl.col("multi_cum").alias("multi_e")), left_on="_last", right_on="seg_start", strategy="backward")
        .select("_row", "seg_e", "cov_e", "multi_e")
    )
    placed = first.join(last, on="_row", how="left").sort("_row")
    # Running counts include the segment itself, so a range [s, e] holds cum[e] - cum[s] + (s's own).
    own_cov = (pl.col("depth_s") >= 1).cast(pl.Int64)
    own_multi = (pl.col("depth_s") >= 2).cast(pl.Int64)
    multi = pl.col("multi_e") - pl.col("multi_s") + own_multi
    covered = pl.col("cov_e") - pl.col("cov_s") + own_cov
    one_segment = pl.col("seg_s") == pl.col("seg_e")
    reason = (
        pl.when(multi > 0)
        .then(pl.lit("multiple_targets"))
        .when(one_segment & (pl.col("depth_s") == 1))
        .then(None)
        .when(covered > 0)
        .then(pl.lit("split_across_blocks"))
        .otherwise(pl.lit("no_chain_block"))
    )
    placed = placed.with_columns(reason.alias("_reason2")).join(
        blocks.select(pl.col("block").alias("block_s"), "t_start", "q_chrom", "q_start", "q_size", "q_minus"), on="block_s", how="left"
    )
    length = pl.col("_e0") - pl.col("_s0")
    offset = pl.col("_s0") - pl.col("t_start")
    q0 = pl.when(pl.col("q_minus")).then(pl.col("q_size") - (pl.col("q_start") + offset + length)).otherwise(pl.col("q_start") + offset)
    target = pl.col("q_chrom").replace_strict({name: canonical_chrom(name) for name in placed["q_chrom"].drop_nulls().unique().to_list()}, default=None, return_dtype=pl.String)
    placed = placed.with_columns(
        pl.when(pl.col("_reason2").is_null()).then(target).alias("_q_chrom"),
        pl.when(pl.col("_reason2").is_null()).then(q0).alias("_q0"),
        pl.when(pl.col("_reason2").is_null()).then(pl.col("q_minus")).alias("_minus"),
    )
    placed = placed.with_columns(
        pl.when(pl.col("_reason2").is_null() & ~pl.col("_q_chrom").is_in(list(targets))).then(pl.lit("target_contig_unmatched")).otherwise(pl.col("_reason2")).alias("_reason2")
    )
    return rows.join(placed.select("_row", "_reason2", "_q_chrom", "_q0", "_minus"), on="_row", how="left", maintain_order="left")


def _read_target(rows: pl.DataFrame, store: ReferenceStore, targets: dict[str, Contig]) -> pl.DataFrame:
    """The target's REF over each placed row's lifted span, and whether the row keeps or swaps its alleles.

    A minus-strand indel whose alleles share their first base is re-anchored (`_reanchor`): its span
    moves one base left on the target, so the base read first is the new padding base.
    """
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    lengths = alts.list.eval(pl.element().str.len_chars())
    ref_length = pl.col("ref").str.len_chars()
    indel = ((lengths.list.min() != ref_length) | (lengths.list.max() != ref_length)).fill_null(False)
    first_bases = pl.concat_list(pl.col("ref"), alts).list.eval(pl.element().str.to_uppercase().str.slice(0, 1))
    padded = first_bases.list.n_unique() == 1
    reanchor = pl.col("_reason2").is_null() & pl.col("_minus") & indel & padded
    rows = rows.with_columns(
        pl.when(reanchor & (pl.col("_q0") == 0)).then(pl.lit("indel_minus_strand")).otherwise(pl.col("_reason2")).alias("_reason2"),
        (reanchor & (pl.col("_q0") > 0)).fill_null(False).alias("_reanchor"),
        pl.lit(None, dtype=pl.String).alias("_t_ref"),
    )
    rows = rows.with_columns(pl.when(pl.col("_reanchor")).then(pl.col("_q0") - 1).otherwise(pl.col("_q0")).alias("_q0"))
    readable = rows.filter(pl.col("_reason2").is_null())
    pieces = []
    for chrom in sorted(readable["_q_chrom"].unique().to_list()):
        part = readable.filter(pl.col("_q_chrom") == chrom).select("_row", "_q0", (pl.col("_q0") + pl.col("_e0") - pl.col("_s0")).alias("_q1"))
        spans = part.select("_q0", "_q1").unique().sort("_q0", "_q1")
        bases = store.fetch(targets[chrom], list(zip(spans["_q0"].to_list(), spans["_q1"].to_list())))
        table = spans.with_columns(pl.Series("_bases", bases, dtype=pl.String))
        pieces.append(part.join(table, on=["_q0", "_q1"], how="left").select("_row", "_bases"))
    if pieces:
        read = pl.concat(pieces)
        rows = rows.join(read, on="_row", how="left", maintain_order="left").with_columns(pl.col("_bases").alias("_t_ref")).drop("_bases")
    lifted_ref = _complement(pl.col("ref"), pl.col("_minus"))
    lifted_alts = _complement_list(alts, pl.col("_minus"))
    # A re-anchored row's alleles are the target's first base, then each complemented allele less its
    # last base (the old padding base). The first base read is that new padding base by construction,
    # so the alleles are compared without it.
    reanchored = pl.col("_reanchor")
    target = pl.when(reanchored).then(pl.col("_t_ref").str.slice(1)).otherwise(pl.col("_t_ref"))
    ref_body = pl.when(reanchored).then(_drop_last(lifted_ref)).otherwise(lifted_ref)
    alts_body = pl.when(reanchored).then(lifted_alts.list.eval(_drop_last(pl.element()))).otherwise(lifted_alts)
    keep = target == ref_body
    swap = alts_body.list.contains(target)
    open_ = pl.col("_reason2").is_null()
    reason = (
        pl.when(open_ & pl.col("_t_ref").is_null())
        .then(pl.lit("target_contig_unmatched"))
        .when(open_ & ~keep & ~swap)
        .then(pl.lit("ref_mismatch_after_lift"))
        .otherwise(pl.col("_reason2"))
    )
    return rows.with_columns(reason.alias("_reason2"), pl.when(open_ & pl.col("_t_ref").is_not_null()).then(~keep).alias("_swap"))


def _drop_last(allele: pl.Expr) -> pl.Expr:
    return allele.str.slice(0, allele.str.len_chars() - 1)


def _revcomp(text: str) -> str:
    return text.translate(_BASE_COMPLEMENT)[::-1]


def _common_prefix(a: str, b: str, limit: int) -> int:
    n = 0
    while n < limit and a[n] == b[n]:
        n += 1
    return n


def _common_suffix(a: str, b: str, limit: int) -> int:
    n = 0
    while n < limit and a[len(a) - 1 - n] == b[len(b) - 1 - n]:
        n += 1
    return n


def _placements(ref: str, allele: str) -> tuple[tuple[int, int, str], tuple[int, int]]:
    """Where `allele` changes `ref`, as offsets into the span: its left-most placement (start, REF
    bases replaced, the bases put in), trimming a common suffix first as `normalize` does, and its
    right-most (start, REF bases replaced), trimming a common prefix first. They differ only for an
    indel in a repeat, which could equally sit anywhere between them."""
    shortest = min(len(ref), len(allele))
    suffix = _common_suffix(ref, allele, shortest)
    prefix = _common_prefix(ref, allele, shortest - suffix)
    left = (prefix, len(ref) - suffix - prefix, allele[prefix : len(allele) - suffix])
    prefix = _common_prefix(ref, allele, shortest)
    suffix = _common_suffix(ref, allele, shortest - prefix)
    return left, (prefix, len(ref) - suffix - prefix)


def _footprint(start: int, replaced: int) -> tuple[int, int]:
    """The source bases a placement needs on one block: those it replaces, or for an insertion
    (nothing replaced) the bases on either side of the point it goes in at."""
    return (start, start + replaced) if replaced else (start - 1, start + 1)


@dataclass
class _GapPlan:
    """One row's lift across gaps of one chain, before the target's bases are read: its blocks in
    order (`t_start`, `t_end`, `q_start`), the chain-strand image of its span (`c0`, `c1`), the forward
    range read (`f0`, `f1`), and each gap's place in the span and on the chain's strand."""

    chrom: str
    minus: bool
    reanchor: bool
    blocks: list[tuple[int, int, int]]
    c0: int
    c1: int
    f0: int
    f1: int
    gaps: list[tuple[int, int, int]]
    """Per gap: the span offset of the last base before it, and the chain-strand range of its target bases."""


def _mapping(block: dict) -> tuple[str, bool, int]:
    """A block's mapping: target contig, strand and offset, what `_segments` counts once."""
    return block["q_chrom"], block["q_minus"], block["q_start"] - block["t_start"]


def _gap_plan(
    s0: int, e0: int, segments: list[tuple[int | None, int, int]], info: dict[int, dict], same: dict[tuple, list[dict]], ref: str, alleles: list[str], targets: dict[str, Contig]
) -> "_GapPlan | str":
    """How a row whose REF span `[s0, e0)` crosses block edges lies across gaps of one chain, or why it
    cannot be lifted so. `segments` are the chain axis segments the span touches, in order, as (the
    block the axis kept for each or None, its start, its end); `info` those blocks by id; `same` every
    block of the contig by mapping; `alleles` every allele but REF.

    The axis keeps one block per mapping, whichever of the blocks repeating it (F14) came last, so
    which chain it names says nothing. The row lifts where one chain holds each piece of the span a
    segment covers, whole, in one block of that segment's mapping; any such chain gives the same
    plan (a piece that is not the span's first or last ends where its segment does, and so does every
    block holding it), and the lowest chain id is taken."""
    if segments[0][0] is None or segments[-1][0] is None:
        return "split_across_blocks"
    holders: list[dict[int, dict]] = []
    for block, start, end in segments:
        if block is None:
            continue
        p0, p1 = max(s0, start), min(e0, end)
        holders.append({b["chain"]: b for b in same[_mapping(info[block])] if b["t_start"] <= p0 and p1 <= b["t_end"]})
    common = set.intersection(*(set(held) for held in holders))
    if not common:
        return "split_across_chains"
    chain = min(common)
    chosen: list[dict] = []
    for held in holders:
        if not chosen or chosen[-1]["block"] != held[chain]["block"]:
            chosen.append(held[chain])
    first = chosen[0]
    chrom = canonical_chrom(first["q_chrom"])
    if chrom not in targets:
        return "target_contig_unmatched"
    blocks = [(block["t_start"], block["t_end"], block["q_start"]) for block in chosen]

    def block_of(x: int) -> int | None:
        return next((index for index, (start, end, _) in enumerate(blocks) if start <= x < end), None)

    for allele in alleles:
        (start, replaced, _), (right, right_replaced) = _placements(ref, allele)
        low, high = _footprint(start, replaced)
        low, high = min(low, _footprint(right, right_replaced)[0]), max(high, _footprint(right, right_replaced)[1])
        home = block_of(s0 + low)
        if home is None or home != block_of(s0 + high - 1):
            return "gap_in_change"
    pairs = list(pairwise(blocks))
    if any(b[0] - a[1] >= GAP_LIMIT or b[2] - (a[2] + a[1] - a[0]) >= GAP_LIMIT for a, b in pairs):
        return "gap_too_large"

    def image(x: int) -> int:
        start, _, q_start = blocks[block_of(x)]
        return q_start + x - start

    c0, c1 = image(s0), image(e0 - 1) + 1
    q_size, minus = first["q_size"], first["q_minus"]
    f0, f1 = (q_size - c1, q_size - c0) if minus else (c0, c1)
    every = [ref, *alleles]
    indel = any(len(allele) != len(ref) for allele in alleles)
    reanchor = minus and indel and len({allele[:1] for allele in every}) == 1
    if reanchor and f0 == 0:
        return "indel_minus_strand"
    gaps = [(a[1] - 1 - s0, a[2] + a[1] - a[0], b[2]) for a, b in pairs]
    return _GapPlan(chrom, minus, reanchor, blocks, c0, c1, f0 - reanchor, f1, gaps)


def _gap_alleles(plan: _GapPlan, s0: int, ref: str, read: str) -> Callable[[str], str]:
    """The function giving each allele on the target, from the row's plan and the target bases read
    over `[plan.f0, plan.f1)`: the parent's REF carried across (each base on a block at its image, the
    gap's own target bases where the chain puts them, the gap's source-only bases dropped) with the
    allele's change put in at its left-most placement; reverse-complemented on a minus strand, and
    re-anchored there when it is a padded indel, as step 3 re-anchors any."""
    main = read[1:] if plan.reanchor else read
    chain_strand = _revcomp(main) if plan.minus else main
    mapped = [any(start <= s0 + i < end for start, end, _ in plan.blocks) for i in range(len(ref))]
    after = {offset: chain_strand[q0 - plan.c0 : q1 - plan.c0] for offset, q0, q1 in plan.gaps}

    def lifted(allele: str) -> str:
        if allele == ref:
            start, replaced, inserted = 0, 0, ""
        else:
            (start, replaced, inserted), _ = _placements(ref, allele)
        out: list[str] = []
        for i in range(len(ref)):
            if i == start:
                out.append(inserted)
            if mapped[i] and not start <= i < start + replaced:
                out.append(ref[i])
            out.append(after.get(i, ""))
        if start == len(ref):
            out.append(inserted)
        text = "".join(out)
        if plan.minus:
            text = _revcomp(text)
        return read[0] + text[:-1] if plan.reanchor else text

    return lifted


def _place_gap(rows: pl.DataFrame, blocks: pl.DataFrame, segments: pl.DataFrame, store: ReferenceStore, targets: dict[str, Contig]) -> pl.DataFrame:
    """Rows whose REF span crosses a block edge (`_place` said `split_across_blocks`), lifted across
    gaps of one chain where that is well defined (step 2b), or why not (`_reason2`). Few rows, so the
    per-row work is plain Python: the target's REF (`_t_ref`), whether it swaps, and every allele on
    the target (`_GAP_PLACED`)."""
    columns = {"_reason2": pl.String, "_q_chrom": pl.String, "_q0": pl.Int64, "_minus": pl.Boolean, "_reanchor": pl.Boolean, "_t_ref": pl.String, "_swap": pl.Boolean, **_GAP_PLACED}
    rows = rows.drop([name for name in columns if name in rows.columns])
    if rows.is_empty():
        return _with_nulls(rows, columns)
    axis = segments.with_columns(pl.col("seg_start").shift(-1).fill_null(_FAR).alias("seg_end"))
    first = rows.select("_row", "_s0").sort("_s0").join_asof(axis.select("seg_start", pl.col("seg").alias("_seg_s")), left_on="_s0", right_on="seg_start", strategy="backward")
    last = (
        rows.select("_row", (pl.col("_e0") - 1).alias("_last"))
        .sort("_last")
        .join_asof(axis.select("seg_start", pl.col("seg").alias("_seg_e")), left_on="_last", right_on="seg_start", strategy="backward")
    )
    touched = (
        first.select("_row", "_seg_s")
        .join(last.select("_row", "_seg_e"), on="_row")
        .with_columns(pl.int_ranges(pl.col("_seg_s"), pl.col("_seg_e") + 1, dtype=pl.UInt32).alias("seg"))
        .explode("seg", empty_as_null=False)
        .join(axis.select("seg", "seg_start", "seg_end", "seg_block"), on="seg", how="left")
        .sort("_row", "seg")
        .group_by("_row", maintain_order=True)
        .agg("seg_block", "seg_start", "seg_end")
    )
    spans = {row: list(zip(block_ids, starts, ends)) for row, block_ids, starts, ends in touched.iter_rows()}
    ids = {block for span in spans.values() for block, _, _ in span if block is not None}
    info = {row["block"]: row for row in blocks.filter(pl.col("block").is_in(list(ids))).iter_rows(named=True)}
    # Every block repeating a mapping the axis kept, so a span is judged by the chains, not by which block was kept.
    kept = pl.DataFrame(sorted({_mapping(block) for block in info.values()}), schema={"q_chrom": pl.String, "q_minus": pl.Boolean, "_offset": pl.Int64}, orient="row")
    repeating = blocks.with_columns((pl.col("q_start") - pl.col("t_start")).alias("_offset")).join(kept, on=["q_chrom", "q_minus", "_offset"], how="semi")
    same: dict[tuple, list[dict]] = {}
    for block in repeating.sort("block").iter_rows(named=True):
        same.setdefault(_mapping(block), []).append(block)

    plans: dict[int, _GapPlan | str] = {}
    for row, s0, e0, ref, alts in rows.select("_row", "_s0", "_e0", pl.col("ref").str.to_uppercase(), "alts").iter_rows():
        alleles = sorted({allele.upper() for allele in alts or []} - {ref})
        plans[row] = _gap_plan(s0, e0, spans[row], info, same, ref, alleles, targets)
    reads: dict[int, str] = {}
    for chrom in sorted({plan.chrom for plan in plans.values() if isinstance(plan, _GapPlan)}):
        wanted = [(row, plan) for row, plan in plans.items() if isinstance(plan, _GapPlan) and plan.chrom == chrom]
        bases = store.fetch(targets[chrom], [(plan.f0, plan.f1) for _, plan in wanted])
        reads |= {row: text for (row, _), text in zip(wanted, bases)}

    out: dict[str, list] = {name: [] for name in ("_row", *columns)}
    for row, s0, ref, alts, gt_alleles, source_alleles in rows.select("_row", "_s0", pl.col("ref").str.to_uppercase(), "alts", "gt_alleles", "source_alleles").iter_rows():
        plan = plans[row]
        out["_row"].append(row)
        if not isinstance(plan, _GapPlan):
            for name in columns:
                out[name].append(plan if name == "_reason2" else None)
            continue
        target = reads[row][:-1] if plan.reanchor else reads[row]
        lifted = _gap_alleles(plan, s0, ref, reads[row])
        old_ref = lifted(ref)
        new_alts = [lifted(allele.upper()) for allele in alts or []]
        gt_lifted, source_lifted = (
            None if alleles is None else [None if allele is None else lifted(allele.upper()) for allele in alleles] for alleles in (gt_alleles, source_alleles)
        )
        swap = target != old_ref
        reason = "ref_mismatch_after_lift" if swap and target not in new_alts else None
        values = {
            "_reason2": reason,
            "_q_chrom": plan.chrom,
            "_q0": plan.f0,
            "_minus": plan.minus,
            "_reanchor": False,
            "_t_ref": target,
            "_swap": swap,
            "_gap": True,
            "_gap_reanchor": plan.reanchor,
            "_gap_ref": old_ref,
            "_gap_alts": new_alts,
            "_gap_gt_alleles": gt_lifted,
            "_gap_source_alleles": source_lifted,
        }
        for name in columns:
            out[name].append(values[name])
    placed = pl.DataFrame(out, schema={"_row": rows.schema["_row"], **columns})
    return rows.join(placed, on="_row", how="left", maintain_order="left")


def _reanchored(rows: pl.DataFrame) -> pl.DataFrame:
    """Re-anchored rows' alleles on the target: the new padding base (REF's first, read from the
    target), then each complemented allele less its old padding base. Few rows, so plain Python."""

    def moved(anchor: str, allele: str | None) -> str | None:
        return None if allele is None else anchor + allele[:-1]

    columns: dict[str, list] = {name: [] for name in ("_lifted_ref", "alts", "gt_alleles", "source_alleles")}
    for ref, lifted_ref, alts, gt_alleles, source_alleles in rows.select("ref", "_lifted_ref", "alts", "gt_alleles", "source_alleles").iter_rows():
        anchor = ref[0]
        columns["_lifted_ref"].append(moved(anchor, lifted_ref))
        columns["alts"].append([moved(anchor, allele) for allele in alts])
        for name, alleles in (("gt_alleles", gt_alleles), ("source_alleles", source_alleles)):
            columns[name].append(None if alleles is None else [moved(anchor, allele) for allele in alleles])
    return rows.with_columns(
        pl.Series("_lifted_ref", columns["_lifted_ref"], dtype=pl.String),
        *[pl.Series(name, columns[name], dtype=CALLS[name]) for name in ("alts", "gt_alleles", "source_alleles")],
    )


def _lifted_rows(work: pl.DataFrame, facet_id: str, source_chrom: str, counts: Counter[str]) -> pl.DataFrame:
    """The lifted `calls` rows from placed, read, kept-or-swapped parent rows."""
    if work.is_empty():
        return pl.DataFrame(schema=CALLS)
    minus = pl.col("_minus")
    alts = pl.col("alts").fill_null(pl.lit([], dtype=CALLS["alts"]))
    # A row lifted across a gap of one chain has its alleles on the target already (`_place_gap`).
    gap = pl.col("_gap").fill_null(False)
    lifted_ref = pl.when(gap).then(pl.col("_gap_ref")).otherwise(_complement(pl.col("ref"), minus))
    moved_contig = pl.col("_q_chrom") != source_chrom
    array = pl.col("source_ref").is_null()
    counts["lifted"] += work.height
    counts["minus"] += int(work["_minus"].sum())
    counts["reanchored"] += int(work["_reanchor"].sum()) + int(work["_gap_reanchor"].fill_null(False).sum())
    counts["gap_lifted"] += int(work.select(gap.sum()).item())
    counts["swapped"] += int(work["_swap"].sum())
    counts["other_contig"] += int(work.select(moved_contig.sum()).item())
    counts["dosage_withheld"] += int(work.select((pl.col("_swap") & pl.col("dosage").is_not_null()).sum()).item())
    counts["copy_numbers_withheld"] += int(work.select((pl.col("_swap") & pl.col("alt_copy_numbers").is_not_null()).sum()).item())
    counts["phase_cleared"] += int(work.select((moved_contig & pl.col("phased").is_not_null()).sum()).item())
    out = work.with_columns(
        pl.lit(facet_id).alias("facet_id"),
        pl.col("_q_chrom").alias("chrom"),
        _to_pos(pl.col("_q0")).alias("pos"),
        (pl.col("_q0") + pl.when(gap).then(pl.col("_t_ref")).otherwise(pl.col("ref")).str.len_chars()).cast(pl.UInt32).alias("end_pos"),
        pl.col("_t_ref").alias("ref"),
        lifted_ref.alias("_lifted_ref"),
        pl.when(gap).then(pl.col("_gap_alts")).otherwise(_complement_list(alts, minus)).alias("alts"),
        pl.when(gap).then(pl.col("_gap_gt_alleles")).when(pl.col("gt_alleles").is_not_null()).then(_complement_list(pl.col("gt_alleles"), minus)).alias("gt_alleles"),
        pl.when(gap)
        .then(pl.col("_gap_source_alleles"))
        .when(pl.col("source_alleles").is_not_null())
        .then(_complement_list(pl.col("source_alleles"), minus))
        .alias("source_alleles"),
        pl.lit(None, dtype=CALLS["source_ref"]).alias("source_ref"),
        pl.lit(None, dtype=CALLS["ref_check"]).alias("ref_check"),
        pl.lit(None, dtype=CALLS["source_gt"]).alias("source_gt"),
        pl.lit(None, dtype=CALLS["vrs_allele_ids"]).alias("vrs_allele_ids"),
        pl.lit(None, dtype=CALLS["vrs_errors"]).alias("vrs_errors"),
        # Offsets from the source's POS and END, which a minus-strand block reverses.
        pl.lit(None, dtype=CALLS["cipos"]).alias("cipos"),
        pl.lit(None, dtype=CALLS["ciend"]).alias("ciend"),
        pl.when(pl.col("_swap")).then(None).otherwise(pl.col("dosage")).alias("dosage"),
        pl.when(pl.col("_swap")).then(None).otherwise(pl.col("alt_copy_numbers")).alias("alt_copy_numbers"),
        *[pl.when(moved_contig).then(None).otherwise(pl.col(name)).alias(name) for name in ("phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities")],
        array.alias("_array"),
    )
    if out["_reanchor"].any():
        out = pl.concat([out.filter(~pl.col("_reanchor")), _reanchored(out.filter(pl.col("_reanchor")))], how="vertical")
    swapped = out.filter(pl.col("_swap") & ~pl.col("_array"))
    if swapped.height:
        out = pl.concat([out.filter(~(pl.col("_swap") & ~pl.col("_array"))), _swap(swapped)], how="vertical")
    # Array rows (no source REF) take their ALTs from the genotype against the new REF, as resolve_ref
    # derives them: the observed alleles other than REF, sorted, with `gt` indexing `[ref] + alts`.
    redo = pl.col("_array") & (pl.col("_swap") | pl.col("_minus"))
    if out.select(redo.any()).item():
        out = _derive_array_alleles(out.with_columns(redo.alias("_array")), Counter())
    return out.select(CALLS.names())


def _swap(rows: pl.DataFrame) -> pl.DataFrame:
    """VCF rows whose target REF is one of their ALTs: that ALT becomes REF and the old REF takes its
    place among the ALTs; the genotype's indices follow. Few rows, so plain Python."""
    alts_out: list[list[str]] = []
    gts_out: list[list[int | None] | None] = []
    for old_ref, alts, gt, new_ref in rows.select("_lifted_ref", "alts", "gt", "ref").iter_rows():
        k = alts.index(new_ref) + 1
        swapped = list(alts)
        swapped[k - 1] = old_ref
        alts_out.append(swapped)
        gts_out.append(None if gt is None else [None if g is None else k if g == 0 else 0 if g == k else g for g in gt])
    return rows.with_columns(pl.Series("alts", alts_out, dtype=CALLS["alts"]), pl.Series("gt", gts_out, dtype=CALLS["gt"]))


def _warnings(
    counts: Counter[str], reasons: dict[str, int], block_report: pl.DataFrame, chain: Chain, store: ReferenceStore, source: str, target: str
) -> list[str]:
    source_label = store.assembly_label(source) or source
    target_label = store.assembly_label(target) or target
    lines = [
        f"{counts['lifted']} of {counts['rows']} row(s) lifted from {source_label} to {target_label} through chain {chain.name} "
        f"(sha256 {chain.sha256}); REF read from collection {target}"
    ]
    if counts["swapped"]:
        lines.append(f"{counts['swapped']} lifted row(s) had the target's REF as an ALT: REF and that ALT swapped, the genotype re-indexed")
    if counts["minus"]:
        lines.append(f"{counts['minus']} lifted row(s) on minus-strand chain blocks: alleles reverse-complemented")
    if counts["sv_lifted"]:
        lines.append(
            f"{counts['sv_lifted']} lifted row(s) were symbolic SVs or copy-number segments, placed by both ends (length within "
            f"{SV_LENGTH_TOLERANCE:.0%}), {counts['sv_minus']} of them on minus-strand blocks (cipos and ciend mirrored and swapped)"
        )
    if counts["bnd_lifted"]:
        lines.append(
            f"{counts['bnd_lifted']} lifted row(s) were breakends, each mate lifted through the chain and named inside ALT as the target "
            f"collection names it; {counts['bnd_minus']} with POS on a minus-strand block (t reverse-complemented to the brackets' other "
            f"side, cipos mirrored), {counts['bnd_mate_minus']} mate(s) on one (brackets turned, ciend mirrored)"
        )
    if counts["bnd_mate_assembly"]:
        lines.append(
            f"{counts['bnd_mate_assembly']} breakend mate(s) on a contig of the ##assembly file kept as written (no chain moves it; "
            "the brackets turned by POS's strand alone)"
        )
    mate_causes = sorted(key for key in counts if key.startswith("bnd_mate:"))
    if mate_causes:
        lines.append(
            "breakend mate(s) not lifted, by the mate position's reason: " + ", ".join(f"{counts[key]} {key.removeprefix('bnd_mate:')}" for key in mate_causes)
        )
    if counts["reanchored"]:
        lines.append(
            f"{counts['reanchored']} of them indels re-anchored on the target's base before the lifted span (one base left); "
            "run normalize on the lifted facet to left-align them on the target"
        )
    if counts["gap_lifted"]:
        lines.append(
            f"{counts['gap_lifted']} lifted row(s) crossed a gap of one chain (under {GAP_LIMIT} bases each side) outside the bases "
            "they change: REF read over the lifted span, each allele's change put into the target's sequence there"
        )
    if counts["dosage_withheld"]:
        lines.append(f"{counts['dosage_withheld']} swapped row(s): dosage withheld (it counted the other allele)")
    if counts["copy_numbers_withheld"]:
        lines.append(f"{counts['copy_numbers_withheld']} swapped row(s): alt_copy_numbers withheld (the swapped ALT's slot now holds the old REF)")
    if counts["other_contig"]:
        lines.append(f"{counts['other_contig']} row(s) lifted to another contig than their source's; phase sets cleared on {counts['phase_cleared']}")
    for reason, count in reasons.items():
        lines.append(f"{count} row(s) not lifted: {REASON_TEXT[reason]} [{reason}]")
    if counts["blocks"]:
        lines.append(
            f"{counts['blocks']} block(s): {counts['blocks_whole']} lifted whole, {counts['blocks_split']} lifted in part or cut at a chain "
            f"gap or edge, {counts['blocks_left']} left behind; {counts['block_pieces']} lifted block(s) written, kind, GQ, MIN_DP and GL "
            "unchanged; the lifted facet's absence_means is unknown"
        )
    if counts["block_minus"] or counts["block_other_contig"]:
        lines.append(
            f"{counts['block_minus']} lifted block(s) on minus-strand chain blocks, {counts['block_other_contig']} on another contig than their source's"
        )
    for reason, pieces, bases in block_report.iter_rows():
        lines.append(f"{bases} block base(s) in {pieces} piece(s) not lifted: {REASON_TEXT[reason]} [{reason}]")
    if counts["block_overlaps"]:
        lines.append(
            f"{counts['block_overlaps']} lifted block(s) overlap another on the target (two source regions the chain sends to one place); both kept"
        )
    return lines
