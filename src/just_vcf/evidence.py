"""One graded genotype evidence value across sources, computed on demand (layer 7; charter 7, 9, 11).

Sources write evidence differently: posteriors (GP), likelihoods (GL, PL), a phred confidence of the
called genotype (GQ), an expected allele count (DS), a bare genotype (arrays, many VCFs), reference
blocks, or nothing. The `calls` and `blocks` frames keep each quantity as written (docs/SCHEMA.md
Decision 3). This module reads them into one shape a consumer can rank and compare, per row and per
asked position, and stores nothing: the frames do not change, so `SCHEMA_VERSION` does not either.

Per row it gives a probability vector over the genotypes of the alleles it names, in VCF genotype
order, where the source supports one; a phred-like confidence of one genotype; and the basis the two
rest on. A genotype the source asserts without quantifying gets no number (`confidence_kind =
asserted`), and an unknown one gets none either (`unknown`): null is "not quantified", while 0 is a
quantified "no confidence". Nothing here is a boolean "called" (CLAUDE.md §2). The rules are written
out in docs/SCHEMA.md, "Graded genotype evidence".
"""

import math
from collections.abc import Mapping
from enum import Enum

import polars as pl

from just_vcf.callset import absence_is_reference
from just_vcf.coverage import Coverage, callable_sites
from just_vcf.loading import Loaded
from just_vcf.schema import CONTIG_RANK, BlockKind, EvidenceBasis, GlSource, VariantClass, enum_dtype

__all__ = [
    "GENOTYPE_EVIDENCE",
    "GENOTYPE_EVIDENCE_AT",
    "BlockGrade",
    "ConfidenceKind",
    "GradeBasis",
    "VectorPrior",
    "genotype_evidence",
    "genotype_evidence_at",
]

_LN10 = math.log(10.0)


class GradeBasis(str, Enum):
    """Which quantity a row's graded value rests on. just-vcf's own vocabulary (charter 2)."""

    GP = "gp"
    """The source's posterior probabilities (`gp`), as given."""
    GL = "gl"
    """The source's likelihoods (`gl`, from VCF GL or PL: the frame keeps one column for both, one
    quantity on two scales), normalized under a flat prior. The row's `gl_source` says which key."""
    GQ = "gq"
    """The source's phred confidence of the genotype it called (`gq`), and nothing more."""
    HQ = "hq"
    """The source's phred confidence in each haplotype's call (`hq`, CGI's allele scores), on a fully
    called genotype whose every slot has one: `confidence` is the lowest, the genotype being wrong
    where any haplotype is. No vector: one number per haplotype ranks no other genotype."""
    DOSAGE = "dosage"
    """Only an expected allele count (`dosage`): no genotype vector follows from one number."""
    GT_ONLY = "gt_only"
    """A genotype written with no quantity: asserted, unquantified (arrays, many VCFs)."""
    BLOCK = "block"
    """A reference block with no quantity the block rule reads (`BlockGrade`): reference asserted over the
    block, unquantified. Under the default rule that is a block with neither likelihoods nor, where
    `gq_on_reference` asks for it, a GQ; under another rule, a block lacking what that rule needs."""
    NO_CALL = "no_call"
    """The source wrote that it looked and called nothing (`./.`, a `no_call` block), with no vector."""
    ABSENT_REFERENCE = "absent_reference"
    """No row covers the position, and the facet's silence reads as reference (`absence_is_reference`)."""
    ABSENT_UNKNOWN = "absent_unknown"
    """No row covers the position, and the facet's silence says nothing there."""


class ConfidenceKind(str, Enum):
    """What `confidence` is, so a null is never read as a low number."""

    QUANTIFIED = "quantified"
    """`confidence` is a number the source's own quantity gives (0 included)."""
    ASSERTED = "asserted"
    """The source states a genotype (or reference) without quantifying it; `confidence` is null."""
    UNKNOWN = "unknown"
    """No genotype is stated; `confidence` is null."""


class VectorPrior(str, Enum):
    """How `genotype_probabilities` was made from the source's quantity."""

    SOURCE = "source"
    """GP as the source wrote it: a posterior, under whatever prior the source used."""
    FLAT = "flat"
    """GL or PL normalized with every genotype equally likely beforehand (just-vcf's choice)."""
    EVEN_REMAINDER = "even_remainder"
    """GQ alone: the called genotype gets `1 - 10^(-GQ/10)` and the rest is spread evenly over the
    other genotypes, which the source did not rank."""


class BlockGrade(str, Enum):
    """Which quantity grades a `ref` block (desk #94: the choice is the person's; the default is today's).

    A block's GQ and its likelihoods say different things: GQ is the band's floor, the weakest base the
    block covers, while its PL is one value for the whole band (docs/SCHEMA.md, "Graded genotype
    evidence"). The rule governs `ref` blocks only; a `no_call` block is graded by its likelihoods, or
    is `no_call`, under every rule. A block that lacks what the rule needs is `block` (reference
    asserted, the number withheld): never graded by another rule's quantity."""

    LIKELIHOODS = "likelihoods"
    """By its likelihoods (`gl`), as written; without usable ones, by its GQ where `gq_on_reference`
    asks for it, else `block`. The default."""
    GQ = "gq"
    """By its GQ alone (`gq`, `prior = even_remainder`), whatever `gq_on_reference` says: hom-ref gets
    `1 - 10^(-GQ/10)` and the rest is spread evenly, at the ploidy the block's likelihoods give where
    it has usable ones (no vector otherwise). Without a GQ, `block`."""
    LOWER_OF_GQ_AND_LIKELIHOODS = "lower_of_gq_and_likelihoods"
    """By whichever of the two states less confidence in hom-ref: GQ against the phred its likelihoods
    give genotype 0 (`-10·log10(1 - p_0)`), the likelihoods on a tie. The row's `grade_basis` says which
    (`gq` or `gl`). A block lacking either is `block`."""


_KIND = {
    GradeBasis.GP: ConfidenceKind.QUANTIFIED,
    GradeBasis.GL: ConfidenceKind.QUANTIFIED,
    GradeBasis.GQ: ConfidenceKind.QUANTIFIED,
    GradeBasis.HQ: ConfidenceKind.QUANTIFIED,
    GradeBasis.DOSAGE: ConfidenceKind.ASSERTED,
    GradeBasis.GT_ONLY: ConfidenceKind.ASSERTED,
    GradeBasis.BLOCK: ConfidenceKind.ASSERTED,
    GradeBasis.NO_CALL: ConfidenceKind.UNKNOWN,
    GradeBasis.ABSENT_REFERENCE: ConfidenceKind.ASSERTED,
    GradeBasis.ABSENT_UNKNOWN: ConfidenceKind.UNKNOWN,
}
"""`confidence_kind` follows from the basis alone; it is a column so a consumer never has to map it."""

GENOTYPE_EVIDENCE = pl.Schema(
    {
        "facet_id": pl.String,
        "chrom": pl.String,
        "pos": pl.UInt32,
        "end_pos": pl.UInt32,
        "evidence_basis": enum_dtype(EvidenceBasis),
        "block_kind": enum_dtype(BlockKind),
        "grade_basis": enum_dtype(GradeBasis),
        "gl_source": enum_dtype(GlSource),
        "confidence_kind": enum_dtype(ConfidenceKind),
        "genotype_alleles": pl.List(pl.String),
        "genotype_probabilities": pl.List(pl.Float64),
        "prior": enum_dtype(VectorPrior),
        "genotype_index": pl.UInt32,
        "confidence": pl.Float64,
    }
)
"""What `genotype_evidence` returns: one row per `calls` row and per `blocks` row."""

GENOTYPE_EVIDENCE_AT = pl.Schema(
    {
        "chrom": pl.String,
        "pos": pl.Int64,
        "facet_id": pl.String,
        "record_pos": pl.UInt32,
        "record_end_pos": pl.UInt32,
        **{name: dtype for name, dtype in GENOTYPE_EVIDENCE.items() if name not in ("facet_id", "chrom", "pos", "end_pos")},
        "absence_means": pl.String,
    }
)
"""What `genotype_evidence_at` returns, before the coverage columns `coverage=` adds."""


def _genotype_count(alleles: pl.Expr, ploidy: pl.Expr) -> pl.Expr:
    """VCF's number of genotypes, C(alleles + ploidy - 1, ploidy), for ploidy 1 to 4; null beyond."""
    a = alleles.cast(pl.Int64)
    return (
        pl.when(ploidy == 1).then(a)
        .when(ploidy == 2).then(a * (a + 1) // 2)
        .when(ploidy == 3).then(a * (a + 1) * (a + 2) // 6)
        .when(ploidy == 4).then(a * (a + 1) * (a + 2) * (a + 3) // 24)
        .otherwise(None)
    )  # fmt: skip


def _usable(vector: pl.Expr, genotypes: pl.Expr, *, probabilities: bool) -> pl.Expr:
    """A vector that can be graded: one finite-or-minus-infinite value per genotype, two genotypes or more.

    A probability vector must also be non-negative with a positive sum. Anything else is passed over
    and the row falls to the next basis, never patched."""
    whole = (
        vector.is_not_null()
        & (vector.list.len().cast(pl.Int64) == genotypes)
        & (genotypes >= 2)
        & ~vector.list.eval(pl.element().is_null() | pl.element().is_nan()).list.any()
    )
    if probabilities:
        return whole & ~vector.list.eval(pl.element() < 0).list.any() & (vector.list.sum() > 0)
    return whole & vector.list.max().is_finite()


def _phred_of_top(log10_values: pl.Expr) -> pl.Expr:
    """-10·log10(1 - p_top) from a vector of log10 weights, computed in log space.

    The weights of the genotypes other than the top one are summed relative to the top, so a top
    genotype of probability 1 - 10^-300 still gives 3,000 rather than infinity by rounding. Infinity
    means every other genotype has weight zero as written (a GP of `0,1,0`)."""
    element = pl.element()
    ordered = element.sort(descending=True)
    rest = (ordered - ordered.first()).slice(1)
    second = rest.first()
    others = pl.when(second == float("-inf")).then(float("-inf")).otherwise(second + ((rest - second) * _LN10).exp().sum().log(10))
    log_others = log10_values.list.eval(others).list.first()
    return -10.0 * (log_others - (1.0 + (log_others * _LN10).exp()).log(10))


def _phred_of_first(log10_values: pl.Expr) -> pl.Expr:
    """-10·log10(1 - p_0) of genotype 0 (hom-ref) from a vector of log10 weights, in log space.

    Infinite where every other genotype has weight zero as written."""
    element = pl.element()
    top = element.max()
    log_total = top + ((element - top) * _LN10).exp().sum().log(10)
    rest = element.slice(1)
    rest_top = rest.max()
    log_rest = pl.when(rest_top == float("-inf")).then(float("-inf")).otherwise(rest_top + ((rest - rest_top) * _LN10).exp().sum().log(10))
    return -10.0 * (log10_values.list.eval(log_rest - log_total).list.first())


def _ref_block_grade(gl: pl.Expr, gq: pl.Expr, gl_ok: pl.Expr, *, gq_on_reference: bool, block_grade: BlockGrade) -> pl.Expr:
    """Which quantity grades a `ref` block under `block_grade`: `gl`, `gq`, or null where the block lacks
    what the rule needs (reference asserted, the number withheld). `gl_ok` says its likelihoods are
    usable; `gq_on_reference` speaks only under the default rule. Shared by `genotype_evidence` and the
    consensus's block reading, so the two cannot drift."""
    has_gq = gq.is_not_null()
    by_gl, by_gq = pl.lit(GradeBasis.GL.value), pl.lit(GradeBasis.GQ.value)
    if block_grade is BlockGrade.LIKELIHOODS:
        return pl.when(gl_ok).then(by_gl).when(has_gq & pl.lit(gq_on_reference)).then(by_gq)
    if block_grade is BlockGrade.GQ:
        return pl.when(has_gq).then(by_gq)
    lower = pl.when(gq.cast(pl.Float64) < _phred_of_first(gl)).then(by_gq).otherwise(by_gl)
    return pl.when(gl_ok & has_gq).then(lower)


def _called_index(gt: pl.Expr) -> pl.Expr:
    """The VCF genotype index of a fully called haploid or diploid `gt`; null otherwise."""
    first = gt.list.get(0, null_on_oob=True).cast(pl.Int64)
    second = gt.list.get(1, null_on_oob=True).cast(pl.Int64)
    low, high = pl.min_horizontal(first, second), pl.max_horizontal(first, second)
    return (
        pl.when(gt.list.eval(pl.element().is_null()).list.any()).then(None)
        .when(gt.list.len() == 1).then(first)
        .when(gt.list.len() == 2).then(high * (high + 1) // 2 + low)
        .otherwise(None)
    )  # fmt: skip


def _with_gl_source(rows: pl.LazyFrame, facets: pl.LazyFrame) -> pl.LazyFrame:
    """`rows` with their facet's `gl_source` as `_gl_source`, in their order: the key a row graded by its
    likelihoods names. Joined lazily, so nothing is read until the answer is collected."""
    sources = facets.select("facet_id", pl.col("gl_source").alias("_gl_source"))
    return rows.join(sources, on="facet_id", how="left", maintain_order="left")


def _calls_evidence(calls: pl.LazyFrame, gq_on_reference: bool, facets: pl.LazyFrame, keep: tuple[str, ...] = ()) -> pl.LazyFrame:
    gt, gt_alleles = pl.col("gt"), pl.col("gt_alleles")
    fully_called = gt.is_not_null() & (gt.list.len() > 0) & ~gt.list.eval(pl.element().is_null()).list.any()
    called_any = (
        gt_alleles.list.eval(pl.element().is_not_null()).list.any().fill_null(False)
        | gt.list.eval(pl.element().is_not_null()).list.any().fill_null(False)
        | (gt_alleles.is_null() & gt.is_null() & pl.col("source_gt").is_not_null())
    )
    carries_alt = gt.list.eval(pl.element() > 0).list.any().fill_null(False)
    # The ploidy a vector is indexed at: what GT says, else diploid, as the spec says to assume.
    ploidy = pl.coalesce(gt.list.len(), gt_alleles.list.len(), pl.lit(2, dtype=pl.UInt32)).cast(pl.Int64)
    vector_genotypes = _genotype_count(pl.col("source_alleles").list.len(), ploidy)
    called_alleles = pl.concat_list(pl.coalesce(pl.col("source_ref"), pl.col("ref")), pl.col("alts"))
    called_genotypes = _genotype_count(called_alleles.list.len(), gt.list.len().cast(pl.Int64))
    gq_applies = pl.col("gq").is_not_null() & fully_called & (pl.lit(True) if gq_on_reference else carries_alt)
    # One score per GT slot, none missing. Unlike GQ it grades a reference genotype as well: each
    # haplotype's score is the confidence in that haplotype's call, reference or not.
    hq = pl.col("hq")
    hq_applies = (
        fully_called & hq.is_not_null() & (hq.list.len() == gt.list.len()) & ~hq.list.eval(pl.element().is_null()).list.any()
    ).fill_null(False)

    rows = _with_gl_source(calls, facets).select(
        "facet_id", "chrom", "pos", "end_pos", "evidence_basis", "gt", "gq", *keep,
        pl.col("source_alleles"), called_alleles.alias("_called_alleles"),
        pl.col("gp").cast(pl.List(pl.Float64)).alias("_gp"),
        pl.col("gl").cast(pl.List(pl.Float64)).alias("_gl"),
        "_gl_source",
        pl.col("dosage"),
        vector_genotypes.alias("_vector_genotypes"),
        called_genotypes.alias("_called_genotypes"),
        _called_index(gt).alias("_called_index"),
        called_any.alias("_called_any"),
        gq_applies.alias("_gq_applies"),
        pl.when(hq_applies).then(hq.list.min().cast(pl.Float64)).alias("_hq_min"),
    )  # fmt: skip
    gp_ok = _usable(pl.col("_gp"), pl.col("_vector_genotypes"), probabilities=True)
    gl_ok = _usable(pl.col("_gl"), pl.col("_vector_genotypes"), probabilities=False)
    basis = (
        pl.when(gp_ok).then(pl.lit(GradeBasis.GP.value))
        .when(gl_ok).then(pl.lit(GradeBasis.GL.value))
        .when(pl.col("_gq_applies")).then(pl.lit(GradeBasis.GQ.value))
        .when(pl.col("_hq_min").is_not_null()).then(pl.lit(GradeBasis.HQ.value))
        .when(pl.col("dosage").is_not_null()).then(pl.lit(GradeBasis.DOSAGE.value))
        .when(pl.col("_called_any")).then(pl.lit(GradeBasis.GT_ONLY.value))
        .otherwise(pl.lit(GradeBasis.NO_CALL.value))
    )  # fmt: skip
    rows = rows.with_columns(basis.alias("grade_basis"), pl.lit(None, dtype=enum_dtype(BlockKind)).alias("block_kind"))
    return _graded(rows, called_index=pl.col("_called_index"), called_genotypes=pl.col("_called_genotypes"),
                   vector_alleles=pl.col("source_alleles"), called_alleles=pl.col("_called_alleles"), keep=keep)  # fmt: skip


def _blocks_evidence(blocks: pl.LazyFrame, gq_on_reference: bool, facets: pl.LazyFrame, keep: tuple[str, ...] = (),
                     block_grade: BlockGrade = BlockGrade.LIKELIHOODS) -> pl.LazyFrame:  # fmt: skip
    # A block's likelihoods are over [REF, <*>]: two alleles, so three values diploid and two haploid.
    # The frame keeps no REF base for a block, so that slot of `genotype_alleles` is null.
    gl = pl.col("gl").cast(pl.List(pl.Float64))
    genotypes = gl.list.len().cast(pl.Int64)
    block_alleles = pl.concat_list(pl.lit(None, dtype=pl.String), pl.lit("<*>"))
    is_ref = pl.col("block_kind") == BlockKind.REF.value
    rows = _with_gl_source(blocks, facets).select(
        "facet_id", "chrom", "pos", "end_pos", "gq", *keep,
        pl.lit(EvidenceBasis.REF_BLOCK.value).cast(enum_dtype(EvidenceBasis)).alias("evidence_basis"),
        pl.col("kind").alias("block_kind"),
        pl.lit(None, dtype=pl.List(pl.Float64)).alias("_gp"),
        gl.alias("_gl"),
        "_gl_source",
        block_alleles.alias("_alleles"),
        pl.when(genotypes.is_in([2, 3])).then(genotypes).otherwise(None).alias("_vector_genotypes"),
        pl.lit(None, dtype=pl.Float64).alias("_hq_min"),
    )  # fmt: skip
    gl_ok = _usable(pl.col("_gl"), pl.col("_vector_genotypes"), probabilities=False)
    # A `ref` block by the rule (`BlockGrade`); a `no_call` block by its likelihoods under every rule.
    ref_basis = _ref_block_grade(pl.col("_gl"), pl.col("gq"), gl_ok, gq_on_reference=gq_on_reference, block_grade=block_grade)
    basis = (
        pl.when(is_ref).then(pl.coalesce(ref_basis, pl.lit(GradeBasis.BLOCK.value)))
        .when(gl_ok).then(pl.lit(GradeBasis.GL.value))
        .otherwise(pl.lit(GradeBasis.NO_CALL.value))
    )  # fmt: skip
    rows = rows.with_columns(basis.alias("grade_basis"))
    # Without likelihoods a block's ploidy is not in the frame, so a GQ-graded block gets no vector:
    # the called genotype is hom-ref (index 0) at a ploidy nobody wrote down. With usable ones, their
    # length gives it, and a GQ-graded block's vector is spread over that many genotypes.
    return _graded(rows, called_index=pl.when(is_ref).then(pl.lit(0, dtype=pl.Int64)), called_genotypes=pl.when(gl_ok).then(pl.col("_vector_genotypes")),
                   vector_alleles=pl.col("_alleles"), called_alleles=pl.col("_alleles"), keep=keep)  # fmt: skip


def _graded(rows: pl.LazyFrame, *, called_index: pl.Expr, called_genotypes: pl.Expr, vector_alleles: pl.Expr, called_alleles: pl.Expr,
            keep: tuple[str, ...] = ()) -> pl.LazyFrame:  # fmt: skip
    """The graded columns from `grade_basis` and the source quantities on `rows`."""
    basis = pl.col("grade_basis")
    gp, gl = pl.col("_gp"), pl.col("_gl")
    element = pl.element()
    normalized_gl = gl.list.eval(((element - element.max()) * _LN10).exp() / ((element - element.max()) * _LN10).exp().sum())
    error = (-(pl.col("gq").cast(pl.Float64)) / 10.0 * _LN10).exp()
    spread = error / (called_genotypes - 1).cast(pl.Float64)
    gq_vector = pl.when(called_genotypes >= 2).then(
        pl.concat_list(spread.repeat_by(called_index), 1.0 - error, spread.repeat_by(called_genotypes - called_index - 1))
    )
    is_gp, is_gl, is_gq = basis == GradeBasis.GP.value, basis == GradeBasis.GL.value, basis == GradeBasis.GQ.value
    is_hq = basis == GradeBasis.HQ.value
    # A genotype the row names without a vector: asserted, or quantified by its haplotypes' scores.
    asserted_genotype = basis.is_in([GradeBasis.GT_ONLY.value, GradeBasis.DOSAGE.value, GradeBasis.BLOCK.value, GradeBasis.HQ.value])
    kind = basis.replace_strict({b.value: k.value for b, k in _KIND.items()}, return_dtype=pl.String)
    return rows.with_columns(
        kind.cast(enum_dtype(ConfidenceKind)).alias("confidence_kind"),
        pl.when(is_gp | is_gl).then(vector_alleles).when(is_gq | asserted_genotype).then(called_alleles).alias("genotype_alleles"),
        pl.when(is_gp).then(gp).when(is_gl).then(normalized_gl).when(is_gq).then(gq_vector).alias("genotype_probabilities"),
        pl.when(is_gp).then(pl.lit(VectorPrior.SOURCE.value))
        .when(is_gl).then(pl.lit(VectorPrior.FLAT.value))
        .when(is_gq & gq_vector.is_not_null()).then(pl.lit(VectorPrior.EVEN_REMAINDER.value))
        .cast(enum_dtype(VectorPrior)).alias("prior"),
        pl.when(is_gp).then(gp.list.arg_max().cast(pl.Int64))
        .when(is_gl).then(gl.list.arg_max().cast(pl.Int64))
        .when(is_gq | asserted_genotype).then(called_index)
        .cast(pl.UInt32).alias("genotype_index"),
        pl.when(is_gp).then(_phred_of_top(gp.list.eval(pl.element().log(10))))
        .when(is_gl).then(_phred_of_top(gl))
        .when(is_gq).then(pl.col("gq").cast(pl.Float64))
        .when(is_hq).then(pl.col("_hq_min"))
        .alias("confidence"),
        basis.cast(enum_dtype(GradeBasis)).alias("grade_basis"),
        # The key the likelihoods were read from, on the rows they grade; every other row rests on neither.
        pl.when(is_gl).then(pl.col("_gl_source")).cast(enum_dtype(GlSource)).alias("gl_source"),
    ).select(*GENOTYPE_EVIDENCE.names(), *keep)  # fmt: skip


def genotype_evidence(loaded: Loaded, *, gq_on_reference: bool = False, block_grade: BlockGrade | str = BlockGrade.LIKELIHOODS) -> pl.LazyFrame:
    """One graded value per `calls` row and per `blocks` row of `loaded` (`GENOTYPE_EVIDENCE`).

    Per row, the first basis that applies (docs/SCHEMA.md, "Graded genotype evidence"):

    - `gp`: the posterior vector as written, `prior = source`;
    - `gl`: the likelihoods (GL or PL) normalized with a flat prior, `prior = flat`. `gl_source` names
      the key the facet read them from (`facets.gl_source`), and is null on every row of another basis;
    - `gq`: a fully called genotype carrying an ALT, with a GQ. VCF defines GQ as conditional on the
      site being variant, so by default it grades no reference genotype and no reference block;
      `gq_on_reference=True` lets it, as single-sample callers compute it. The called genotype gets
      `1 - 10^(-GQ/10)`, the rest is spread evenly (`prior = even_remainder`), and `confidence` is GQ;
    - `hq`: a fully called genotype whose every slot has a haplotype score (`hq`, CGI's allele scores),
      reference genotypes included: `confidence` is the lowest score, no vector;
    - `dosage`: a dosage and nothing above. No vector, no number;
    - `gt_only` (a call) or `block` (a `ref` block): asserted, unquantified;
    - `no_call`: nothing was called and no vector given.

    A `ref` block is graded by `block_grade` (`BlockGrade`): by its likelihoods (the default, the order
    above), by its GQ alone, or by the lower of the two; a block lacking what the rule needs is
    `block`, never graded by another rule's quantity. A `no_call` block is graded as above under every
    rule.

    A vector must have one value per genotype at the row's ploidy (GT's, else diploid), at least two;
    otherwise the row falls to the next basis. `genotype_alleles` names what the vector or
    `genotype_index` indexes: `source_alleles` for a vector, `[source_ref or ref] + alts` for a genotype,
    `[null, "<*>"]` for a block (the frame keeps no REF base for one). `confidence` is
    `-10·log10(1 - p)` of the top genotype for a vector (infinite where every other genotype has
    probability zero as written), GQ for `gq`, and null where nothing quantifies it. `genotype_index`
    is the top genotype's index in VCF order, or the called genotype's. Ordered canonically, calls
    before blocks at one position. Lazy: nothing is read until collected.
    """
    facets = loaded.facets
    block_grade = BlockGrade(block_grade)
    frame = pl.concat(
        [_calls_evidence(loaded.calls, gq_on_reference, facets), _blocks_evidence(loaded.blocks, gq_on_reference, facets, block_grade=block_grade)],
        how="vertical",
    )
    rank = pl.col("chrom").replace_strict(CONTIG_RANK, default=len(CONTIG_RANK), return_dtype=pl.UInt16)
    origin = pl.col("evidence_basis").cast(pl.String)
    return frame.sort([pl.col("facet_id"), rank, pl.col("chrom"), pl.col("pos"), origin, pl.col("end_pos")], nulls_last=True, maintain_order=True)


def _one_facet(loaded: Loaded) -> Mapping[str, object]:
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"genotype_evidence_at takes one facet; this Loaded has {facets.height}")
    return facets.row(0, named=True)


def genotype_evidence_at(
    loaded: Loaded,
    sites: pl.DataFrame | pl.LazyFrame,
    *,
    variant_class: VariantClass | str,
    coverage: Coverage | None = None,
    min_depth: int | None = None,
    gq_on_reference: bool = False,
    block_grade: BlockGrade | str = BlockGrade.LIKELIHOODS,
    force: bool = False,
) -> pl.DataFrame:
    """The graded evidence a one-facet `loaded` holds at each of `sites` (`chrom` canonical, `pos` 1-based).

    A site gets one row per record that covers it: a `calls` row from its `pos` to its `end_pos` (its
    `pos` alone where the span is unknown) and a `blocks` row over its span, each graded as
    `genotype_evidence` grades it, with the record's span as `record_pos` and `record_end_pos`. So a
    site inside a called deletion is not absent, and a site with an SNV and a deletion over it gets
    both. A site nothing covers gets one row, `evidence_basis = absence`: `absent_reference` (asserted,
    no number, `genotype_index` 0) where `absence_is_reference(facet, variant_class)` is True, and
    `absent_unknown` otherwise. `variant_class` is the class of variant the question is about; a
    facet's silence speaks only for the classes it called. `force` is `absence_is_reference`'s: it lifts
    a content filter's hold for the class the filter kept (the `.snp.` half of a split delivery). `absence_means` comes along on every row,
    so `reference` and `weak_reference` stay apart. `gq_on_reference` and `block_grade` are
    `genotype_evidence`'s.

    `coverage` (with `min_depth`, which has no default) adds `callable`, `callable_reason` and the
    covering interval's depths per site, from `callable_sites`. They sit beside the answer and change
    no basis: reads reaching a base are not a call (docs/SCHEMA.md, coverage facets). Which coverage
    belongs to which callset is the caller's claim. Sites come back in the order given, a duplicate
    once.
    """
    facet = _one_facet(loaded)
    asked_class = VariantClass(variant_class)
    block_grade = BlockGrade(block_grade)
    if coverage is not None and min_depth is None:
        raise ValueError("coverage needs min_depth: what depth makes a base callable is the caller's choice")
    asked = (
        sites.lazy().select(pl.col("chrom").cast(pl.String), pl.col("pos").cast(pl.Int64)).unique(maintain_order=True).collect()
        .with_columns(pl.int_range(pl.len(), dtype=pl.UInt32).alias("_site"))
    )  # fmt: skip
    contigs = asked["chrom"].unique().to_list()
    # The records are matched first and only the matched rows graded, so a genome's millions of rows
    # are never graded to answer a panel. A call starting at the site is an equality join, and only
    # calls spanning more than their first base go through the range join. A facet's blocks do not
    # overlap one another (a gVCF's and CGI's tile the genome), so the block over a site is the last
    # one starting at or before it, found by an as-of join.
    ask = asked.lazy().select("_site", pl.col("chrom").alias("_site_chrom"), pl.col("pos").alias("_site_pos"))
    calls = loaded.calls.filter(pl.col("chrom").is_in(contigs)).with_columns(
        pl.col("pos").cast(pl.Int64).alias("_start"),
        pl.coalesce(pl.col("end_pos"), pl.col("pos")).cast(pl.Int64).alias("_end"),
    )
    at_start = calls.join(ask, left_on=["chrom", "_start"], right_on=["_site_chrom", "_site_pos"], how="inner")
    # One range join per contig: with two inequalities and no equality, polars runs it as an IEJoin; an
    # equality on `chrom` beside them makes it a filtered cross join (38 GB on a CGI genome's panel).
    spanning = calls.filter(pl.col("_end") > pl.col("_start"))
    inside = [
        spanning.filter(pl.col("chrom") == contig)
        .join_where(ask.filter(pl.col("_site_chrom") == contig), pl.col("_site_pos") > pl.col("_start"), pl.col("_site_pos") <= pl.col("_end"))
        .drop("_site_chrom", "_site_pos")
        for contig in sorted(contigs)
    ]  # fmt: skip
    hit_calls = pl.concat([at_start, *inside], how="vertical_relaxed")
    blocks = loaded.blocks.filter(pl.col("chrom").is_in(contigs)).with_columns(pl.col("pos").cast(pl.Int64).alias("_start")).sort("chrom", "_start")
    hit_blocks = (
        ask.rename({"_site_chrom": "chrom"}).sort("chrom", "_site_pos")
        .join_asof(blocks, left_on="_site_pos", right_on="_start", by="chrom", strategy="backward", check_sortedness=False)
        .filter(pl.col("facet_id").is_not_null() & (pl.col("_site_pos") <= pl.col("end_pos").cast(pl.Int64)))
        .drop("_site_pos", "_start")
    )  # fmt: skip
    graded = pl.concat([
        _calls_evidence(hit_calls, gq_on_reference, loaded.facets, keep=("_site",)),
        _blocks_evidence(hit_blocks, gq_on_reference, loaded.facets, keep=("_site",), block_grade=block_grade),
    ])  # fmt: skip
    covered = (
        graded.rename({"pos": "record_pos", "end_pos": "record_end_pos"}).drop("chrom")
        .join(asked.lazy(), on="_site", how="inner")
        .collect()
    )  # fmt: skip
    absent = asked.join(covered.select("_site").unique(), on="_site", how="anti")
    reads_reference = absence_is_reference(facet, asked_class, force=force) is True
    absence = absent.with_columns(
        pl.lit(facet["facet_id"]).alias("facet_id"),
        pl.lit(EvidenceBasis.ABSENCE.value).cast(enum_dtype(EvidenceBasis)).alias("evidence_basis"),
        pl.lit((GradeBasis.ABSENT_REFERENCE if reads_reference else GradeBasis.ABSENT_UNKNOWN).value).cast(enum_dtype(GradeBasis)).alias("grade_basis"),
        pl.lit((ConfidenceKind.ASSERTED if reads_reference else ConfidenceKind.UNKNOWN).value).cast(enum_dtype(ConfidenceKind)).alias("confidence_kind"),
        pl.lit(0 if reads_reference else None, dtype=pl.UInt32).alias("genotype_index"),
    )  # fmt: skip
    origin = pl.col("evidence_basis").cast(pl.String)
    answer = (
        pl.concat([covered.lazy(), absence.lazy()], how="diagonal_relaxed")
        .with_columns(pl.lit(facet["absence_means"], dtype=pl.String).alias("absence_means"))
        .sort("_site", "record_pos", origin, "record_end_pos", nulls_last=True, maintain_order=True)
    )
    answer = answer.select("_site", *[pl.col(name).cast(dtype) for name, dtype in GENOTYPE_EVIDENCE_AT.items()])
    if coverage is not None:
        depth = callable_sites(coverage, asked.select("chrom", "pos"), min_depth=min_depth).with_columns(pl.int_range(pl.len(), dtype=pl.UInt32).alias("_site"))
        answer = answer.join(depth.drop("chrom", "pos"), on="_site", how="left", maintain_order="left")
    return answer.drop("_site").collect()
