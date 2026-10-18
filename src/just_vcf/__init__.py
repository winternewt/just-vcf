"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

from just_vcf.loading import LOADABLE_FORMATS, Loaded, MultiSampleError, UnsupportedFormatError, load, vcf_samples
from just_vcf.schema import (
    BLOCKS,
    CALLS,
    FACETS,
    SCHEMA_VERSION,
    AbsenceMeans,
    AssemblyLabelBasis,
    BlockKind,
    CallsetType,
    EvidenceBasis,
    PhaseSetKind,
    RefCheck,
    SplitSignal,
)
from just_vcf.sniff import Sniffed, SniffOutcome, open_payload, sniff

__version__ = version("just-vcf")

__all__ = [
    "BLOCKS",
    "CALLS",
    "FACETS",
    "LOADABLE_FORMATS",
    "SCHEMA_VERSION",
    "AbsenceMeans",
    "AssemblyLabelBasis",
    "BlockKind",
    "CallsetType",
    "Loaded",
    "MultiSampleError",
    "EvidenceBasis",
    "PhaseSetKind",
    "RefCheck",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "UnsupportedFormatError",
    "__version__",
    "load",
    "open_payload",
    "sniff",
    "vcf_samples",
]
