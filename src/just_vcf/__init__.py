"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

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
)
from just_vcf.sniff import Sniffed, SniffOutcome, sniff

__version__ = version("just-vcf")

__all__ = [
    "BLOCKS",
    "CALLS",
    "FACETS",
    "SCHEMA_VERSION",
    "AbsenceMeans",
    "AssemblyLabelBasis",
    "BlockKind",
    "CallsetType",
    "EvidenceBasis",
    "PhaseSetKind",
    "RefCheck",
    "Sniffed",
    "SniffOutcome",
    "__version__",
    "sniff",
]
