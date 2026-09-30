"""Load a genotype file into the `facets`, `calls` and `blocks` frames (docs/SCHEMA.md).

`load(path)` sniffs the file, reopens the payload `sniff` recognised and hands it to the loader for
that format. A format without a loader is refused with `UnsupportedFormatError`, never read by a
loader for a look-alike format.
"""

import hashlib
import json
import logging
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.loaders import rsid_table
from just_vcf.schema import (
    BLOCKS,
    FACETS,
    SCHEMA_VERSION,
    AbsenceMeans,
    AssemblyLabelBasis,
    CallsetType,
    canonical_order,
    contig_map,
    empty,
)
from just_vcf.sniff import Sniffed, open_payload, sniff

logger = logging.getLogger(__name__)

LOADABLE_FORMATS: frozenset[str] = rsid_table.FORMATS
"""Format ids `load` can read. Grows as loaders land."""

_IMPUTED_FORMATS = frozenset({"genes_for_good_imputed"})
_UNKNOWN_CALLSET_FORMATS = frozenset({"plink_23"})
"""PLINK's 23andMe layout is written from any source (the one real file looks imputed or sequenced)."""


class UnsupportedFormatError(ValueError):
    """The file was recognised, but no loader reads its format yet, or it was not recognised at all."""

    def __init__(self, sniffed: Sniffed) -> None:
        self.sniffed = sniffed
        what = sniffed.format or sniffed.kind or sniffed.reason or sniffed.outcome.value
        super().__init__(f"{sniffed.path}: no loader for {what} ({sniffed.outcome.value})")


class Loaded(BaseModel):
    """One loaded input: its facet and its rows. The frames are lazy."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    facets: pl.LazyFrame
    calls: pl.LazyFrame
    blocks: pl.LazyFrame
    warnings: tuple[str, ...] = ()
    """Rows the loader could not take as written, one line per reason with a count."""
    schema_version: str = SCHEMA_VERSION


def load(source: Path | str | Sniffed) -> Loaded:
    """Load `source` (a path, or a `sniff` result for it) into the three frames."""
    sniffed = source if isinstance(source, Sniffed) else sniff(Path(source))
    if not sniffed or sniffed.format not in LOADABLE_FORMATS:
        raise UnsupportedFormatError(sniffed)
    spec = rsid_table.spec_for(sniffed.format)
    with open_payload(sniffed) as handle:
        data = handle.read()

    calls, warnings = rsid_table.read_calls(data, spec)
    if spec.columns and rsid_table.column_line_ok(data, spec) is False:
        warnings.append("the column line does not read rsid, chromosome, position, genotype; columns taken by position")
    source_sha256 = _sha256(sniffed.path)
    facet_id = _facet_id(source_sha256, sniffed.member, spec.id, options={})

    callset = (
        CallsetType.IMPUTED
        if spec.id in _IMPUTED_FORMATS
        else None
        if spec.id in _UNKNOWN_CALLSET_FORMATS
        else CallsetType.ARRAY
    )
    facet = {
        "facet_id": facet_id,
        "source_sha256": source_sha256,
        "source_path": str(sniffed.path),
        "member": sniffed.member,
        "format_id": spec.id,
        "family": spec.family.value,
        "vendor": None if spec.vendor == "unknown" else spec.vendor,
        "assembly_label": sniffed.build,
        "assembly_label_basis": AssemblyLabelBasis.STATED.value if sniffed.build else None,
        "callset_type": callset.value if callset else None,
        "absence_means": AbsenceMeans.UNKNOWN.value,
        "contig_map": contig_map(calls["source_chrom"].unique().to_list()),
    }
    facets = pl.LazyFrame([facet], schema={name: FACETS[name] for name in facet}).select(
        [pl.col(name) if name in facet else pl.lit(None, dtype=dtype).alias(name) for name, dtype in FACETS.items()]
    )
    for warning in warnings:
        logger.warning("%s: %s", facet_id, warning)
    return Loaded(
        facets=facets,
        calls=canonical_order(rsid_table.to_calls_frame(calls, facet_id)),
        blocks=empty(BLOCKS),
        warnings=tuple(warnings),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _facet_id(source_sha256: str, member: str | None, format_id: str, options: dict) -> str:
    """Deterministic facet id: the source bytes, the member, the format and the load options."""
    key = json.dumps(
        {"sha256": source_sha256, "member": member, "format": format_id, "options": options},
        sort_keys=True,
    )
    return hashlib.sha256(key.encode()).hexdigest()[:32]
