"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

from just_vcf.loading import (
    LOADABLE_FORMATS,
    Loaded,
    MultiSampleError,
    NoSampleError,
    UnsupportedFormatError,
    load,
    vcf_samples,
)
from just_vcf.loaders.table import TableLayoutError
from just_vcf.rejoin import rejoin_split
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
    "NoSampleError",
    "EvidenceBasis",
    "PhaseSetKind",
    "RefCheck",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "TableLayoutError",
    "UnsupportedFormatError",
    "__version__",
    "load",
    "open_payload",
    "rejoin_split",
    "sniff",
    "vcf_samples",
]
