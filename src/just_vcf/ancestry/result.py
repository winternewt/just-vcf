"""The answer every ancestry function gives: `AncestryEstimate` (charter invariant 9).

An estimate places a genome relative to the groups of one reference panel. It is a statistical
estimate with its uncertainty, never a statement of anyone's identity, and it can be "unknown": then
`estimate` is None and `reason` says why. Nothing here defaults to a group.
"""

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

ANCESTRY_METHODS: frozenset[str] = frozenset({"pca_knn", "prive_qp", "consensus"})
"""How an estimate was made. Grows additively, so it is a validated vocabulary, not an enum.

- `pca_knn`: projection onto a reference PCA (FRAPOSA's OADP), then a k-nearest-neighbour vote among
  the panel's reference samples; `proportions` are the PCA-centroid mixture (`pca_nnls`).
- `prive_qp`: Privé's worldwide reference (UK Biobank, 21 groups): projection onto its published PCs,
  then group proportions by simplex-constrained least squares.
- `consensus`: the answered estimates combined over five continental groups, one panel's methods as
  one witness (`consensus.consensus`, `CONSENSUS_COMBINATIONS`).
"""

ANCESTRY_REASONS: frozenset[str] = frozenset(
    {
        "assembly_unknown",
        "assembly_mismatch",
        "panel_missing",
        "too_few_sites",
        "no_estimates",
        "no_continental_groups",
    }
)
"""Why an estimate is unknown. Grows additively.

- `assembly_unknown`: the facet's `assembly_label` is null, so its positions cannot be placed on the
  panel's assembly. Positions are never matched on an assumed build.
- `assembly_mismatch`: the facet's assembly is not the panel's (no liftover happens here).
- `panel_missing`: the panel path holds no panel of the asked name.
- `too_few_sites`: too few panel sites could be read from the facet (`coverage_floor`, `min_sites`).
- `no_estimates`: a consensus was asked of estimates none of which answered.
- `no_continental_groups`: the estimates answered, but every group they put weight on folds onto none of
  the five continental groups (an `aadr_ho` placement in Central Asia or Siberia), so a consensus over
  those groups has nothing to read.
"""

ERROR_METHODS: frozenset[str] = frozenset({"block_bootstrap", "chromosome_jackknife"})
"""How an estimate's standard errors (`proportion_errors`, `pc_errors`) were computed. Grows additively.

- `block_bootstrap`: `prive_qp`'s, a bootstrap over contiguous blocks of the sites read; Privé's
  projection and centroids are sums over sites, so a resample is one weighted sum.
- `chromosome_jackknife`: `pca_knn`'s, a weighted delete-one-block jackknife whose blocks are the
  autosomes the facet and the panel share (Busing, Meijer and van der Leeden 1999, as ADMIXTOOLS uses
  it); each replicate projects the facet again by OADP on the other autosomes' sites
  (`pca.chromosome_jackknife`).
"""

JACKKNIFE_REASONS: frozenset[str] = frozenset({"not_requested", "too_few_blocks", "too_few_sites_per_block"})
"""Why a `pca_knn` estimate that answered carries no jackknife (`jackknife_reason`). Grows additively.

- `not_requested`: the caller passed `jackknife=False`.
- `too_few_blocks`: fewer autosomes with a site read than `MIN_JACKKNIFE_BLOCKS`
  (`n_jackknife_blocks` says how many), too few for a standard error worth reporting.
- `too_few_sites_per_block`: leaving some autosome out would keep no more panel sites than the panel's
  `dim_online`, so that replicate's online loadings already span every kept site and OADP has nothing
  to augment them with (a panel of a few dozen sites). Decided from the sites per autosome, before any
  replicate is projected; `n_jackknife_blocks` still counts the blocks.
"""

PROPORTION_METHODS: frozenset[str] = frozenset({"pca_nnls", "prive_qp"})
"""How `proportions` were computed. Either way they are a projection-space approximation, not
haplotype-level admixture."""

CONSENSUS_COMBINATIONS: frozenset[str] = frozenset({"mean", "product_of_means"})
"""How a consensus combined its witnesses (a witness: the estimates whose panels are not known to have
disjoint reference samples, `ConsensusWitness`). Grows additively.

- `mean`: one witness and no prior. The answer is the mean of the witness's readings, which gives the
  leader no more than the most any one reading gave it: agreement does not sharpen it.
- `product_of_means`: two or more witnesses, or a prior. Each witness's mean, smoothed, is one factor of
  a product, as independent evidence is; the prior is one more.
"""

CONSENSUS_READINGS: frozenset[str] = frozenset({"probabilities", "proportions"})
"""Which distribution of an estimate a witness read: its `probabilities` (`pca_knn`'s neighbour vote) or
its `proportions` (`pca_nnls`, `prive_qp`). One estimate can give a witness both."""


class WitnessReading(BaseModel):
    """One distribution a consensus read, folded onto `CANONICAL_SUPERPOPS`."""

    model_config = ConfigDict(frozen=True)

    panel: str
    method: str
    """The estimate's `method`."""
    reading: str
    """`probabilities` or `proportions` (`CONSENSUS_READINGS`)."""
    facet_ids: tuple[str, ...] = ()
    distribution: dict[str, float] = {}
    """Over the five continental groups, summing to 1; empty for a reading left out because none of its
    groups folds onto one of them."""

    @field_validator("reading")
    @classmethod
    def _known_reading(cls, value: str) -> str:
        if value not in CONSENSUS_READINGS:
            raise ValueError(f"unknown consensus reading {value!r}; known: {sorted(CONSENSUS_READINGS)}")
        return value


class ConsensusWitness(BaseModel):
    """The readings a consensus took as one witness: every estimate on these panels.

    Panels are one witness unless the caller declared their reference samples disjoint, since two
    methods on one set of reference samples and sites are not independent evidence. `mean` is what the
    witness says; `spread` is how far its readings disagree, per group the largest value a reading gave
    it minus the smallest (0 for a witness of one reading).
    """

    model_config = ConfigDict(frozen=True)

    panels: tuple[str, ...]
    readings: tuple[WitnessReading, ...] = ()
    mean: dict[str, float] = {}
    """The readings' mean per continental group; empty when every reading was left out."""
    spread: dict[str, float] = {}
    left_out: tuple[WitnessReading, ...] = ()
    """Readings none of whose groups folds onto a continental group; they enter neither `mean` nor
    `spread`."""


class AncestryEstimate(BaseModel):
    """A genome's placement against one reference panel, or the reason it could not be placed.

    `bool(estimate)` is True only when it answered. The site counts say how much of the panel the
    answer rests on, and how much of that came from absences read as reference rather than from rows
    the source wrote.
    """

    model_config = ConfigDict(frozen=True)

    method: str
    panel: str
    """The panel's name (for a consensus, the names it combined, `+`-joined, sorted)."""
    panel_assembly: str | None = None
    facet_ids: tuple[str, ...] = ()
    facet_assembly: str | None = None
    reason: str | None = None
    """Null exactly when the estimate answered. A value from `ANCESTRY_REASONS` otherwise."""
    detail: str | None = None
    """The reason in words, with the numbers behind it."""

    estimate: str | None = None
    """The panel group the genome lies nearest: the neighbour vote's leader (`pca_knn`), the largest
    proportion (`prive_qp`) or the combined support's leader (`consensus`). None = unknown."""
    probabilities: dict[str, float] = {}
    """`pca_knn`: the fraction of the k nearest reference samples per group (an uncertainty over the
    nearest group, not a genome fraction). `consensus`: the combined support over the continental
    groups, by the rule `combination` names; relative support, not a calibrated probability."""
    margin: float | None = None
    """The leading value of `probabilities` (or of `proportions` for `prive_qp`) minus the runner-up.
    A small margin means the leader barely won."""
    proportions: dict[str, float] | None = None
    """Group proportions summing to 1, an approximation in PC space."""
    proportion_method: str | None = None
    proportion_errors: dict[str, float] | None = None
    """Standard error per proportion, computed as `error_method` says: a block bootstrap (`prive_qp`) or
    a leave-one-chromosome-out jackknife (`pca_knn`). Null where none ran, never 0."""
    error_method: str | None = None
    """How `proportion_errors` and `pc_errors` were computed (`ERROR_METHODS`); null when neither was."""
    group_to_continental: dict[str, str | None] = {}
    """Each group of this panel folded onto `CANONICAL_SUPERPOPS`, None where none fits (Middle
    Eastern groups). The consensus reads it."""
    combination: str | None = None
    """`consensus` only: how its witnesses were combined (`CONSENSUS_COMBINATIONS`). Null otherwise, and
    on a consensus that did not answer."""
    witnesses: tuple[ConsensusWitness, ...] = ()
    """`consensus` only: the witnesses it combined, each with the readings pooled into it."""
    pc_coords: tuple[float, ...] | None = None
    """The projected principal-component coordinates. Their signs depend on the panel's build."""
    pc_errors: tuple[float, ...] | None = None
    """`pca_knn`: the jackknife standard error of each of `pc_coords`, in the same order. Null where the
    jackknife did not run (`jackknife_reason`), never 0."""
    vote_stability: float | None = None
    """`pca_knn`: the share of the jackknife's replicates, one per autosome left out, whose neighbour vote
    names `estimate`. 1.0 means no single chromosome's sites decide the vote. Null where the jackknife did
    not run."""
    n_jackknife_blocks: int | None = None
    """`pca_knn`: the autosomes the facet and the panel share, each with at least one site read: the
    jackknife's blocks. Set when they were counted (the jackknife ran, or was withheld as
    `too_few_blocks` or `too_few_sites_per_block`); null when it was not requested or the estimate did
    not answer."""
    jackknife_reason: str | None = None
    """`pca_knn` only: why an estimate that answered carries no jackknife (`JACKKNIFE_REASONS`). Null when
    the jackknife ran, and on an estimate that did not answer (its `reason` says why)."""

    n_sites_panel: int | None = None
    n_sites_used: int | None = None
    """Panel sites that entered the estimate: observed plus absent-as-reference."""
    n_sites_observed: int | None = None
    """Sites the source wrote a usable diploid call for."""
    n_sites_absent_as_reference: int | None = None
    """Sites the facet did not write, read as homozygous reference because its `absence_means` says
    absence is (weak) evidence of reference, or because a reference block covers them."""
    n_sites_unusable: int | None = None
    """Sites the source wrote but that could not be read against the panel's two alleles: a no-call
    or half-call, a haploid call, another allele, a row flagged by `split_signal`, or two rows that
    both fit."""
    n_sites_ref_unstated: int | None = None
    """Sites the source wrote only rows for whose REF states no base (`ref_check = unstated`, a
    copy-number caller's `N`) and whose ALT is a sequence allele spelled against it: unknown, never
    read as a genotype and never as an absence."""
    coverage: float | None = None
    """`n_sites_used / n_sites_panel`."""

    @field_validator("method")
    @classmethod
    def _known_method(cls, value: str) -> str:
        if value not in ANCESTRY_METHODS:
            raise ValueError(f"unknown ancestry method {value!r}; known: {sorted(ANCESTRY_METHODS)}")
        return value

    @field_validator("reason")
    @classmethod
    def _known_reason(cls, value: str | None) -> str | None:
        if value is not None and value not in ANCESTRY_REASONS:
            raise ValueError(f"unknown ancestry reason {value!r}; known: {sorted(ANCESTRY_REASONS)}")
        return value

    @field_validator("proportion_method")
    @classmethod
    def _known_proportion_method(cls, value: str | None) -> str | None:
        if value is not None and value not in PROPORTION_METHODS:
            raise ValueError(f"unknown proportion method {value!r}; known: {sorted(PROPORTION_METHODS)}")
        return value

    @field_validator("error_method")
    @classmethod
    def _known_error_method(cls, value: str | None) -> str | None:
        if value is not None and value not in ERROR_METHODS:
            raise ValueError(f"unknown error method {value!r}; known: {sorted(ERROR_METHODS)}")
        return value

    @field_validator("jackknife_reason")
    @classmethod
    def _known_jackknife_reason(cls, value: str | None) -> str | None:
        if value is not None and value not in JACKKNIFE_REASONS:
            raise ValueError(f"unknown jackknife reason {value!r}; known: {sorted(JACKKNIFE_REASONS)}")
        return value

    @field_validator("combination")
    @classmethod
    def _known_combination(cls, value: str | None) -> str | None:
        if value is not None and value not in CONSENSUS_COMBINATIONS:
            raise ValueError(f"unknown consensus combination {value!r}; known: {sorted(CONSENSUS_COMBINATIONS)}")
        return value

    @model_validator(mode="after")
    def _answer_or_reason(self) -> "AncestryEstimate":
        if (self.estimate is None) == (self.reason is None):
            raise ValueError("an estimate carries exactly one of `estimate` and `reason`")
        if self.jackknife_reason is not None and (
            self.estimate is None or self.pc_errors is not None or self.vote_stability is not None
        ):
            raise ValueError("`jackknife_reason` says why an answered estimate has no jackknife; it cannot sit beside one")
        return self

    def __bool__(self) -> bool:
        return self.estimate is not None
