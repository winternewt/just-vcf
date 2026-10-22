"""Population inference: where a facet's genotypes lie against a reference panel (RM5).

Behind the `[ancestry]` extra (numpy). `just_vcf` itself never imports this package; import it as
`just_vcf.ancestry`. Every function takes a `Loaded` and answers with an `AncestryEstimate`, which is
an estimate with its uncertainty and its site counts, or unknown with a reason. A facet is matched to a
panel on its own `assembly_label`; a facet without one, or on another assembly, is never matched on an
assumed build.

Ported from just-prs `ancestry/` (docs/PROVENANCE.md), without its `PRSCatalog` wrappers and without
its Hugging Face download: a panel is a local path the caller passes.
"""

try:
    import numpy  # noqa: F401
except ImportError as error:
    raise ImportError(
        "just_vcf.ancestry needs numpy, which comes with the `ancestry` extra: "
        "`uv add 'just-vcf[ancestry]'` or `pip install 'just-vcf[ancestry]'`"
    ) from error

from just_vcf.ancestry.consensus import CANONICAL_SUPERPOPS, consensus, continental_of, to_continental
from just_vcf.ancestry.pca import (
    DEFAULT_COVERAGE_FLOOR,
    DEFAULT_DIM_REF,
    DEFAULT_KNN_K,
    RESOLUTIONS,
    PanelBuild,
    PcaPanel,
    build_pca_panel,
    infer_ancestry,
    load_pca_panel,
    panel_assemblies,
)
from just_vcf.ancestry.prive import PRIVE_ASSEMBLY, PriveBuild, build_prive_reference, prive_proportions
from just_vcf.ancestry.result import ANCESTRY_METHODS, ANCESTRY_REASONS, PROPORTION_METHODS, AncestryEstimate

__all__ = [
    "ANCESTRY_METHODS",
    "ANCESTRY_REASONS",
    "CANONICAL_SUPERPOPS",
    "DEFAULT_COVERAGE_FLOOR",
    "DEFAULT_DIM_REF",
    "DEFAULT_KNN_K",
    "PRIVE_ASSEMBLY",
    "PROPORTION_METHODS",
    "RESOLUTIONS",
    "AncestryEstimate",
    "PanelBuild",
    "PcaPanel",
    "PriveBuild",
    "build_pca_panel",
    "build_prive_reference",
    "consensus",
    "continental_of",
    "infer_ancestry",
    "load_pca_panel",
    "panel_assemblies",
    "prive_proportions",
    "to_continental",
]
