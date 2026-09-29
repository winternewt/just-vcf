"""Conversion of consumer and personal genomics files into VCF."""

from importlib.metadata import version

from just_vcf.sniff import Sniffed, SniffOutcome, sniff

__version__ = version("just-vcf")

__all__ = ["Sniffed", "SniffOutcome", "__version__", "sniff"]
