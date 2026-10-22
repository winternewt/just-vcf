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
from just_vcf.reference import ReferenceMismatchError, ReferenceStore, reference_dir, resolve_ref
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
from just_vcf.writing import ContigNames, UnsortedRowsError, UnwritableRowsError, VcfWritten, write_vcf

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
    "ContigNames",
    "Loaded",
    "MultiSampleError",
    "NoSampleError",
    "EvidenceBasis",
    "PhaseSetKind",
    "RefCheck",
    "ReferenceMismatchError",
    "ReferenceStore",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "TableLayoutError",
    "UnsortedRowsError",
    "UnsupportedFormatError",
    "UnwritableRowsError",
    "VcfWritten",
    "__version__",
    "load",
    "open_payload",
    "reference_dir",
    "rejoin_split",
    "resolve_ref",
    "sniff",
    "vcf_samples",
    "write_vcf",
]
