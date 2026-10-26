"""A facet's assembly identity, from what a VCF header carries (RM4; charter invariants 3 and 9).

Two answers, kept apart because they are different kinds of fact:

- **The label** (`GRCh36`, `GRCh37`, `GRCh38`, `T2T-CHM13`), inferred from the `##contig` lengths of
  the primary chromosomes. The lengths below are domain constants (UCSC `chrom.sizes` for hg19, hg38
  and hs1, checked 2026-10-22; NCBI's NCBI36 assembly report, checked against hg18's, 2026-10-26); no
  two assemblies share a primary chromosome length, on any chromosome. Only the mitochondrion's are
  shared (16,569 by GRCh37, GRCh38 and T2T-CHM13; hg19's and hg18's 16,571), so a matching primary
  length is evidence for one assembly only. A header too thin to tell, or one whose lengths disagree,
  gives no label, never the nearest one.
- **The coordinate system's digest**: the refget Sequence Collection `sorted_name_length_pairs`
  attribute digest of the header's contigs, names verbatim. It is computed exactly as the seqcol
  specification defines it, so it is a conformant value of that attribute. The collection's top-level
  digest is not: it is computed over the `names` and `sequences` attributes, and `sequences` holds the
  refget (sha512t24u) digest of each contig's bases, which no header carries (a `##contig` `md5` is a
  different digest and cannot be converted). `seqcol_digest` therefore stays null until bases are read
  (RM3), and nothing here fills it from a lookup on the label.
"""

from enum import Enum

import polars as pl
from pydantic import BaseModel, ConfigDict
from refget.digests import sha512t24u_digest
from refget.utils import build_sorted_name_length_pairs, canonical_str

from just_vcf.schema import canonical_chrom

_NAMES = [str(n) for n in range(1, 23)] + ["X", "Y"]

ASSEMBLY_LENGTHS: dict[str, dict[str, frozenset[int]]] = {
    # NCBI36 (hg18): NCBI's assembly report for GCF_000001405.12, carried verbatim in
    # assets/ncbi_assembly_reports/; UCSC's hg18 chrom.sizes and Ensembl 54's seq_region give the same.
    "GRCh36": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826,
                    140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254,
                    78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954,
                )),
                strict=True,
            )
        ),
        # NCBI's report lists no mitochondrion; hg18 (assets/ucsc_hg18/) and Ensembl 54 carry the Yoruba
        # chrM (16,571). rCRS (16,569) came with GRCh37, so a build-36 header that has it reads conflicting.
        "MT": frozenset({16571}),
    },
    "GRCh37": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
                    141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                    81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
                )),
                strict=True,
            )
        ),
        # GRCh37 carries the rCRS (16,569); UCSC hg19 kept the older Yoruba chrM (16,571).
        "MT": frozenset({16569, 16571}),
    },
    "GRCh38": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
                    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
                    83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
                )),
                strict=True,
            )
        ),
        "MT": frozenset({16569}),
    },
    "T2T-CHM13": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    248387328, 242696752, 201105948, 193574945, 182045439, 172126628, 160567428, 146259331,
                    150617247, 134758134, 135127769, 133324548, 113566686, 101161492, 99753195, 96330374,
                    84276897, 80542538, 61707364, 66210255, 45090682, 51324926, 154259566, 62460029,
                )),
                strict=True,
            )
        ),
        "MT": frozenset({16569}),
    },
}
"""Primary-chromosome lengths per assembly, by canonical `chrom`. T2T-CHM13 v1.1 has no Y; its other
lengths are v2.0's, so a header of either version matches."""

MIN_MATCHING_CONTIGS = 3
"""Primary chromosomes whose length must match one assembly, and no other, before it is named.
just-prs `detect_genome_build` voted with the same floor. A header with one or two chromosomes (a
region extract, a test file) is too thin: one matching length could be a coincidence of a re-headered
file, three cannot."""


class AssemblyOutcome(str, Enum):
    """Why `infer_assembly` gave the answer it gave. Only `INFERRED` carries a label."""

    INFERRED = "inferred"
    NO_CONTIG_LENGTHS = "no_contig_lengths"
    """The header has no `##contig` line with a length on a primary chromosome."""
    TOO_FEW = "too_few"
    """Fewer than `MIN_MATCHING_CONTIGS` primary chromosomes match an assembly."""
    CONFLICTING = "conflicting"
    """Primary chromosome lengths match different assemblies, or a length matches none."""


class AssemblyInference(BaseModel):
    """The label the header's contig lengths support, or why there is none. Truthy only with a label."""

    model_config = ConfigDict(frozen=True)

    label: str | None
    outcome: AssemblyOutcome
    reason: str
    matched: dict[str, int] = {}
    """Assembly → primary chromosomes whose header length matches it (the mitochondrion excluded)."""
    unmatched: tuple[str, ...] = ()
    """Primary chromosomes whose header length matches no known assembly."""

    def __bool__(self) -> bool:
        return self.label is not None


def infer_assembly(contigs: list[dict]) -> AssemblyInference:
    """The assembly label a header's `##contig` lines support (`Header.contigs`: name, length, md5).

    Only primary chromosomes take part, by canonical name; alt contigs, decoys and unplaced scaffolds
    differ between releases of one assembly and say nothing about the label. The mitochondrion's
    length is shared between assemblies (16,569 by three, 16,571 by hg19 and hg18), so it can
    contradict a label (hg19's 16,571 on GRCh38) but never supports one.
    """
    lengths: dict[str, set[int]] = {}
    for contig in contigs:
        chrom = canonical_chrom(contig["name"])
        if chrom in ASSEMBLY_LENGTHS["GRCh38"] and contig.get("length") is not None:
            lengths.setdefault(chrom, set()).add(int(contig["length"]))
    if not lengths:
        return AssemblyInference(
            label=None, outcome=AssemblyOutcome.NO_CONTIG_LENGTHS, reason="no ##contig line gives a primary chromosome's length"
        )
    matched = {
        name: sum(1 for chrom, seen in lengths.items() if chrom != "MT" and seen <= table[chrom])
        for name, table in ASSEMBLY_LENGTHS.items()
    }
    fitting = sorted(name for name, table in ASSEMBLY_LENGTHS.items() if all(seen <= table[chrom] for chrom, seen in lengths.items()))
    unmatched = tuple(
        sorted(
            (chrom for chrom, seen in lengths.items() if not any(seen <= table[chrom] for table in ASSEMBLY_LENGTHS.values())),
            key=_NAMES_AND_MT.index,
        )
    )
    supported = {name: count for name, count in matched.items() if count}
    if len(fitting) == 1 and matched[fitting[0]] >= MIN_MATCHING_CONTIGS:
        label = fitting[0]
        return AssemblyInference(
            label=label,
            outcome=AssemblyOutcome.INFERRED,
            reason=f"{matched[label]} primary chromosome length(s) match {label} and none contradicts it",
            matched=supported,
        )
    if not fitting:
        parts = [f"{count} match {name}" for name, count in sorted(supported.items())]
        if unmatched:
            parts.append(f"{', '.join(unmatched)} match(es) no known assembly")
        return AssemblyInference(
            label=None,
            outcome=AssemblyOutcome.CONFLICTING,
            reason="contig lengths disagree: " + "; ".join(parts),
            matched=supported,
            unmatched=unmatched,
        )
    best = max(supported.values(), default=0)
    return AssemblyInference(
        label=None,
        outcome=AssemblyOutcome.TOO_FEW,
        reason=f"{best} primary chromosome length(s) match an assembly; {MIN_MATCHING_CONTIGS} are needed",
        matched=supported,
    )


_NAMES_AND_MT = [*_NAMES, "MT"]


ROWS_PER_CONTIG = 2_000
"""SNV rows sampled per primary chromosome when a build is inferred from the rows."""

MIN_SCORED_ROWS = 1_000
"""Sampled rows a collection must be scored on before the rows can name its assembly."""

MIN_SHARE = 0.60
"""The share of scored rows whose genotype carries the reference base, at or above which a collection
fits. On the right build only homozygous non-reference calls miss; on a wrong one a position's base is
unrelated to the probe, so the share falls towards chance. Measured on 153 PGP-HMS array exports
(F9, 2026-10-22): the right build scored 0.665–0.913, every wrong one 0.230–0.334. The low end is
FTDNA's X-chromosome exports, about 2,000 rows of mostly hemizygous calls, where a hom-alt call is
common (a floor of 0.70 refused 27 of the 43), and one FTDNA autosomal export with some probes on the
minus strand (F13)."""

MIN_MARGIN = 0.20
"""How far the best collection's share must stand above the best collection of another label, and above
chance. Two builds that both fit the rows (a wrong label on the same bases) are not told apart."""

_ROW_BASES = ["A", "C", "G", "T"]


class RowsOutcome(str, Enum):
    """Why `infer_assembly_from_rows` gave the answer it gave. Only `INFERRED` carries a label."""

    INFERRED = "inferred"
    TOO_FEW_ROWS = "too_few_rows"
    """Fewer than `MIN_SCORED_ROWS` rows could be scored on any collection."""
    NO_FIT = "no_fit"
    """No collection reaches `MIN_SHARE`."""
    NO_MARGIN = "no_margin"
    """The best share is within `MIN_MARGIN` of another label's, or of chance."""
    UNLABELLED = "unlabelled"
    """The best-scoring collection has no assembly label, so the rows cannot name one."""


class RowsInference(BaseModel):
    """The assembly an array facet's rows support, scored against reference collections. Truthy with a label."""

    model_config = ConfigDict(frozen=True)

    label: str | None
    digest: str | None
    """The collection that scored best, when it named the label."""
    outcome: RowsOutcome
    reason: str
    sampled: int
    """SNV rows in the deterministic sample."""
    scored: dict[str, int] = {}
    """Collection digest → sampled rows on a contig the collection has."""
    shares: dict[str, float] = {}
    """Collection digest → share of its scored rows whose called alleles include its reference base."""
    labels: dict[str, str | None] = {}
    chance: float | None = None
    """The share a random base would get: the mean number of distinct called alleles over four."""

    def __bool__(self) -> bool:
        return self.label is not None

    def summary(self) -> str:
        """The scores as one line, best first: `81.2% on GRCh37 (digest, 24000 rows); …`."""
        order = sorted(self.shares, key=lambda digest: (-self.shares[digest], digest))
        parts = [f"{self.shares[d]:.1%} on {self.labels.get(d) or 'an unlabelled collection'} ({d}, {self.scored[d]} rows)" for d in order]
        if self.chance is not None:
            parts.append(f"chance {self.chance:.1%}")
        return "; ".join(parts)


def sample_snv_rows(calls: pl.LazyFrame, per_contig: int = ROWS_PER_CONTIG) -> pl.DataFrame:
    """A deterministic sample of array SNV rows: `chrom`, `pos`, `alleles` (distinct called bases, sorted).

    Only rows without a source REF (array rows), on 1–22, X and Y, whose called alleles are all single
    bases with at least one called. Per chromosome, rows sorted by position and genotype, every k-th row
    is kept so that at most `per_contig` remain, spread over the whole chromosome. MT is left out: hg19's
    chrM and rCRS differ, and a label covering both cannot be scored there.
    """
    called = pl.col("gt_alleles").list.drop_nulls()
    snv = (
        pl.col("source_ref").is_null()
        & pl.col("chrom").is_in(_NAMES)
        & pl.col("gt_alleles").is_not_null()
        & (called.list.len() > 0)
        & called.list.eval(pl.element().is_in(_ROW_BASES)).list.all()
    )
    rows = (
        calls.filter(snv)
        .select("chrom", "pos", called.list.unique().list.sort().alias("alleles"))
        .sort("chrom", "pos", "alleles")
        .collect()
    )
    index = pl.int_range(pl.len()).over("chrom")
    step = ((pl.len().over("chrom") + per_contig - 1) // per_contig).clip(lower_bound=1)
    return rows.filter(index % step == 0)


def infer_assembly_from_rows(
    sample: pl.DataFrame, bases: dict[str, pl.Series], labels: dict[str, str | None]
) -> RowsInference:
    """The assembly the sampled rows support (F9; charter 9: unknown unless the rows settle it).

    `sample` is `sample_snv_rows`' output; `bases[digest]` is each collection's base at every sampled
    row, parallel to it: null where the collection has no such contig (the row is not scored on it),
    and an empty string where the position is past the contig's end (the row counts against it).
    `labels[digest]` is the collection's assembly label, or None.

    A collection's share is the fraction of its scored rows whose called alleles include its base.
    The best collection names its label when it was scored on at least `MIN_SCORED_ROWS` rows, reaches
    `MIN_SHARE`, and stands `MIN_MARGIN` above both the best collection of any other label and chance
    (`chance`). An unlabelled collection counts as another label: it could be another assembly.
    """
    sampled = sample.height
    chance = sample["alleles"].list.len().mean() / 4 if sampled else None
    scored: dict[str, int] = {}
    shares: dict[str, float] = {}
    for digest, series in sorted(bases.items()):
        on_contig = sample.select("alleles").with_columns(series.alias("_base")).filter(pl.col("_base").is_not_null())
        scored_rows = on_contig.height
        if scored_rows == 0:
            continue
        scored[digest] = scored_rows
        shares[digest] = on_contig.select(pl.col("alleles").list.contains(pl.col("_base")).mean()).item()
    result = {"sampled": sampled, "scored": scored, "shares": shares, "labels": {d: labels.get(d) for d in shares}, "chance": chance}
    if not scored or max(scored.values()) < MIN_SCORED_ROWS:
        most = max(scored.values(), default=0)
        return RowsInference(
            label=None, digest=None, outcome=RowsOutcome.TOO_FEW_ROWS,
            reason=f"{most} SNV row(s) could be scored on a collection; {MIN_SCORED_ROWS} are needed", **result,
        )  # fmt: skip
    eligible = [d for d in shares if scored[d] >= MIN_SCORED_ROWS]
    best = max(eligible, key=lambda d: (shares[d], d))
    label = labels.get(best)
    if shares[best] < MIN_SHARE:
        return RowsInference(
            label=None, digest=None, outcome=RowsOutcome.NO_FIT,
            reason=f"no collection reaches {MIN_SHARE:.0%} of rows carrying its reference base", **result,
        )  # fmt: skip
    if label is None:
        return RowsInference(
            label=None, digest=None, outcome=RowsOutcome.UNLABELLED,
            reason="the collection that fits best has no assembly label", **result,
        )  # fmt: skip
    rivals = [shares[d] for d in shares if labels.get(d) != label]
    second = max([*rivals, chance or 0.0])
    if shares[best] - second < MIN_MARGIN:
        return RowsInference(
            label=None, digest=None, outcome=RowsOutcome.NO_MARGIN,
            reason=f"{label} scores {shares[best]:.1%}, within {MIN_MARGIN:.0%} of {second:.1%}", **result,
        )  # fmt: skip
    return RowsInference(
        label=label, digest=best, outcome=RowsOutcome.INFERRED,
        reason=f"{shares[best]:.1%} of {scored[best]} sampled SNV rows carry {label}'s reference base, against {second:.1%} at best otherwise",
        **result,
    )  # fmt: skip


def sorted_name_length_pairs_digest(contigs: list[dict]) -> str | None:
    """The seqcol `sorted_name_length_pairs` attribute digest of the header's contigs, or None.

    Names are the file's own, verbatim (`chr1` and `1` are different coordinate systems to seqcol).
    None when the header lists no contig, or when any contig lacks a length: a digest over the
    contigs that have one would describe a coordinate system the file never declared. A contig listed
    twice is kept twice, as the header wrote it.
    """
    if not contigs or any(contig.get("length") is None for contig in contigs):
        return None
    collection = {"names": [contig["name"] for contig in contigs], "lengths": [int(contig["length"]) for contig in contigs]}
    return sha512t24u_digest(canonical_str(build_sorted_name_length_pairs(collection)))
