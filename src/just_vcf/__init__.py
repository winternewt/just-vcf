"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

from just_vcf.cache import CacheClearReport, CacheEntryInfo, clear_cache
from just_vcf.callset import absence_is_reference
from just_vcf.karyotype import KaryotypicSex, SexChromosomeEvidence, XCopies, sex_chromosome_evidence
from just_vcf.liftover import Chain, Lifted, liftover, read_chain
from just_vcf.links import link_disagreements
from just_vcf.loading import (
    LOADABLE_FORMATS,
    DamagedBgzfError,
    Loaded,
    MultiSampleError,
    NoSampleError,
    UnsupportedFormatError,
    load,
    table_samples,
    vcf_samples,
)
from just_vcf.loaders.table import TableLayoutError
from just_vcf.normalize import normalize
from just_vcf.person import SUBJECT_CHECKS, Person, SameSubjectCheck, SubjectCheckOutcome, SubjectError, people
from just_vcf.query import Dimension, GenotypeClass, Presence, PresenceBasis, Query, query
from just_vcf.reference import ReferenceMismatchError, ReferenceStore, reference_dir, resolve_ref
from just_vcf.rejoin import rejoin_split
from just_vcf.schema import (
    BLOCKS,
    CALLS,
    FACETS,
    LINKS,
    SCHEMA_VERSION,
    AbsenceMeans,
    AssemblyLabelBasis,
    BlockKind,
    CallsetType,
    ChipBasis,
    ContentFilter,
    ContentFilterBasis,
    EvidenceBasis,
    PhaseSetKind,
    RefCheck,
    SplitSignal,
    StrandCheck,
    VariantClass,
)
from just_vcf.sniff import Sniffed, SniffOutcome, open_payload, sniff
from just_vcf.writing import BlockAllele, ContigNames, UnsortedRowsError, UnwritableRowsError, VcfWritten, write_vcf
from just_vcf.vrs import mint_vrs_ids, vrs_withheld

__version__ = version("just-vcf")

__all__ = [
    "BLOCKS",
    "CALLS",
    "FACETS",
    "LINKS",
    "LOADABLE_FORMATS",
    "SCHEMA_VERSION",
    "SUBJECT_CHECKS",
    "AbsenceMeans",
    "AssemblyLabelBasis",
    "BlockAllele",
    "BlockKind",
    "CacheClearReport",
    "CacheEntryInfo",
    "CallsetType",
    "ChipBasis",
    "Chain",
    "ContentFilter",
    "ContentFilterBasis",
    "ContigNames",
    "DamagedBgzfError",
    "Dimension",
    "Lifted",
    "Loaded",
    "MultiSampleError",
    "NoSampleError",
    "EvidenceBasis",
    "GenotypeClass",
    "KaryotypicSex",
    "PhaseSetKind",
    "Person",
    "Presence",
    "PresenceBasis",
    "Query",
    "RefCheck",
    "ReferenceMismatchError",
    "ReferenceStore",
    "SameSubjectCheck",
    "SexChromosomeEvidence",
    "Sniffed",
    "SniffOutcome",
    "SplitSignal",
    "StrandCheck",
    "SubjectCheckOutcome",
    "SubjectError",
    "TableLayoutError",
    "UnsortedRowsError",
    "UnsupportedFormatError",
    "UnwritableRowsError",
    "VariantClass",
    "VcfWritten",
    "XCopies",
    "__version__",
    "absence_is_reference",
    "clear_cache",
    "liftover",
    "link_disagreements",
    "load",
    "mint_vrs_ids",
    "normalize",
    "people",
    "query",
    "open_payload",
    "read_chain",
    "reference_dir",
    "rejoin_split",
    "resolve_ref",
    "sex_chromosome_evidence",
    "sniff",
    "table_samples",
    "vcf_samples",
    "vrs_withheld",
    "write_vcf",
]
