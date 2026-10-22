"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

from just_vcf.karyotype import KaryotypicSex, SexChromosomeEvidence, XCopies, sex_chromosome_evidence
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
    "KaryotypicSex",
    "PhaseSetKind",
    "RefCheck",
    "SexChromosomeEvidence",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "TableLayoutError",
    "UnsupportedFormatError",
    "XCopies",
    "__version__",
    "load",
    "open_payload",
    "rejoin_split",
    "sex_chromosome_evidence",
    "sniff",
    "vcf_samples",
]
