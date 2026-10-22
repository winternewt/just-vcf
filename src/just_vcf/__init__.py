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
from just_vcf.normalize import normalize
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
    ChipBasis,
    EvidenceBasis,
    PhaseSetKind,
    RefCheck,
    SplitSignal,
    StrandCheck,
)
from just_vcf.sniff import Sniffed, SniffOutcome, open_payload, sniff
from just_vcf.writing import ContigNames, UnsortedRowsError, UnwritableRowsError, VcfWritten, write_vcf
from just_vcf.vrs import mint_vrs_ids, vrs_withheld

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
    "ChipBasis",
    "ContigNames",
    "Loaded",
    "MultiSampleError",
    "NoSampleError",
    "EvidenceBasis",
    "KaryotypicSex",
    "PhaseSetKind",
    "RefCheck",
    "ReferenceMismatchError",
    "ReferenceStore",
    "SexChromosomeEvidence",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "StrandCheck",
    "TableLayoutError",
    "UnsortedRowsError",
    "UnsupportedFormatError",
    "UnwritableRowsError",
    "VcfWritten",
    "XCopies",
    "__version__",
    "load",
    "mint_vrs_ids",
    "normalize",
    "open_payload",
    "reference_dir",
    "rejoin_split",
    "resolve_ref",
    "sex_chromosome_evidence",
    "sniff",
    "vcf_samples",
    "vrs_withheld",
    "write_vcf",
]
