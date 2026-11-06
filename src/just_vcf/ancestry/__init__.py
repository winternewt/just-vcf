"""Population inference: where a facet's genotypes lie against a reference panel (RM5).

Behind the `[ancestry]` extra (numpy). `just_vcf` itself never imports this package; import it as
`just_vcf.ancestry`. Every function takes a `Loaded` and answers with an `AncestryEstimate`, which is
an estimate with its uncertainty and its site counts, or unknown with a reason. A facet is matched to a
panel on its own `assembly_label`; a facet without one, or on another assembly, is never matched on an
assumed build.

Ported from just-prs `ancestry/` (docs/PROVENANCE.md), without its `PRSCatalog` wrappers and without
its Hugging Face download: a panel is a local path the caller passes. `build_aadr_panel` builds one
from AADR's Human Origins release (`just_vcf.ancestry.aadr`).
"""

try:
    import numpy  # noqa: F401
except ImportError as error:
    raise ImportError(
        "just_vcf.ancestry needs numpy, which comes with the `ancestry` extra: "
        "`uv add 'just-vcf[ancestry]'` or `pip install 'just-vcf[ancestry]'`"
    ) from error

from just_vcf.ancestry.aadr import (
    AADR_REGIONS,
    INDIVIDUAL_EXCLUSIONS,
    QC_QUALIFIERS,
    SITE_DROPS,
    WEST_EURASIA,
    AadrBuild,
    build_aadr_panel,
    group_label_flag,
    region_of,
)
from just_vcf.ancestry.consensus import CANONICAL_SUPERPOPS, consensus, continental_of, to_continental
from just_vcf.ancestry.pca import (
    DEFAULT_COVERAGE_FLOOR,
    DEFAULT_DIM_REF,
    DEFAULT_KNN_K,
    MIN_JACKKNIFE_BLOCKS,
    RESOLUTIONS,
    PanelBuild,
    PcaPanel,
    build_pca_panel,
    infer_ancestry,
    load_pca_panel,
    panel_assemblies,
)
from just_vcf.ancestry.prive import PRIVE_ASSEMBLY, PriveBuild, build_prive_reference, prive_proportions
from just_vcf.ancestry.result import (
    ANCESTRY_METHODS,
    ANCESTRY_REASONS,
    CONSENSUS_COMBINATIONS,
    CONSENSUS_READINGS,
    ERROR_METHODS,
    JACKKNIFE_REASONS,
    PROPORTION_METHODS,
    AncestryEstimate,
    ConsensusWitness,
    WitnessReading,
)

__all__ = [
    "AADR_REGIONS",
    "ANCESTRY_METHODS",
    "ANCESTRY_REASONS",
    "CANONICAL_SUPERPOPS",
    "CONSENSUS_COMBINATIONS",
    "CONSENSUS_READINGS",
    "DEFAULT_COVERAGE_FLOOR",
    "DEFAULT_DIM_REF",
    "DEFAULT_KNN_K",
    "ERROR_METHODS",
    "INDIVIDUAL_EXCLUSIONS",
    "JACKKNIFE_REASONS",
    "MIN_JACKKNIFE_BLOCKS",
    "PRIVE_ASSEMBLY",
    "PROPORTION_METHODS",
    "QC_QUALIFIERS",
    "RESOLUTIONS",
    "SITE_DROPS",
    "WEST_EURASIA",
    "AadrBuild",
    "AncestryEstimate",
    "ConsensusWitness",
    "PanelBuild",
    "PcaPanel",
    "PriveBuild",
    "WitnessReading",
    "build_aadr_panel",
    "build_pca_panel",
    "build_prive_reference",
    "consensus",
    "continental_of",
    "group_label_flag",
    "infer_ancestry",
    "load_pca_panel",
    "panel_assemblies",
    "prive_proportions",
    "region_of",
    "to_continental",
]
