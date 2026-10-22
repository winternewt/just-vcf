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
- `consensus`: a product of the answered estimates' distributions over five continental groups.
"""

ANCESTRY_REASONS: frozenset[str] = frozenset(
    {
        "assembly_unknown",
        "assembly_mismatch",
        "panel_missing",
        "too_few_sites",
        "no_estimates",
    }
)
"""Why an estimate is unknown. Grows additively.

- `assembly_unknown`: the facet's `assembly_label` is null, so its positions cannot be placed on the
  panel's assembly. Positions are never matched on an assumed build.
- `assembly_mismatch`: the facet's assembly is not the panel's (no liftover happens here).
- `panel_missing`: the panel path holds no panel of the asked name.
- `too_few_sites`: too few panel sites could be read from the facet (`coverage_floor`, `min_sites`).
- `no_estimates`: a consensus was asked of estimates none of which answered.
"""

PROPORTION_METHODS: frozenset[str] = frozenset({"pca_nnls", "prive_qp"})
"""How `proportions` were computed. Either way they are a projection-space approximation, not
haplotype-level admixture."""


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
    proportion (`prive_qp`) or the posterior's leader (`consensus`). None = unknown."""
    probabilities: dict[str, float] = {}
    """`pca_knn`: the fraction of the k nearest reference samples per group (an uncertainty over the
    nearest group, not a genome fraction). `consensus`: the posterior over the continental groups."""
    margin: float | None = None
    """The leading value of `probabilities` (or of `proportions` for `prive_qp`) minus the runner-up.
    A small margin means the leader barely won."""
    proportions: dict[str, float] | None = None
    """Group proportions summing to 1, an approximation in PC space."""
    proportion_method: str | None = None
    proportion_errors: dict[str, float] | None = None
    """Block-bootstrap standard error per proportion, where the method allows one (`prive_qp`)."""
    group_to_continental: dict[str, str | None] = {}
    """Each group of this panel folded onto `CANONICAL_SUPERPOPS`, None where none fits (Middle
    Eastern groups). The consensus reads it."""
    pc_coords: tuple[float, ...] | None = None
    """The projected principal-component coordinates. Their signs depend on the panel's build."""

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

    @model_validator(mode="after")
    def _answer_or_reason(self) -> "AncestryEstimate":
        if (self.estimate is None) == (self.reason is None):
            raise ValueError("an estimate carries exactly one of `estimate` and `reason`")
        return self

    def __bool__(self) -> bool:
        return self.estimate is not None
