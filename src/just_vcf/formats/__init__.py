"""The format registry: every recognisable input format, in detection order.

Each family module exports `SPECS` (formats with a vendor signature) and `FALLBACKS` (formats
recognised by data shape alone, such as a headerless rsid table). The registry is every module's
`SPECS` in `_MODULES` order, then every module's `FALLBACKS` in the same order, so a signature always
beats a shape. Within a module, the more specific format comes first. `sniff` takes the first match.
"""

from just_vcf.formats import cgi, other_tables, rsid_tables, vendor_csv, vcf
from just_vcf.formats.spec import FormatSpec

# Signatures are specific enough that module order rarely matters; where two could match, the
# stricter module comes first. Fallbacks follow the same order after every signature.
_MODULES = (vcf, cgi, other_tables, vendor_csv, rsid_tables)

REGISTRY: tuple[FormatSpec, ...] = (
    *(spec for module in _MODULES for spec in module.SPECS),
    *(spec for module in _MODULES for spec in module.FALLBACKS),
)

FORMATS: frozenset[str] = frozenset(spec.id for spec in REGISTRY)
"""Every registered format id. Grows additively; validate with `validate_format`."""

BCF_FORMAT = "bcf"
"""The one format recognised by binary magic (`BCF\\2`) rather than by a text head."""

_duplicates = sorted({s.id for s in REGISTRY if sum(t.id == s.id for t in REGISTRY) > 1})
if _duplicates:
    raise RuntimeError(f"format ids registered twice: {', '.join(_duplicates)}")
if BCF_FORMAT not in FORMATS:
    raise RuntimeError(f"the registry must define {BCF_FORMAT!r}")


def validate_format(format_id: str) -> str:
    """`format_id` if it is registered; ValueError otherwise."""
    if format_id not in FORMATS:
        raise ValueError(f"unknown format {format_id!r}; registered: {', '.join(sorted(FORMATS))}")
    return format_id


def get_spec(format_id: str) -> FormatSpec:
    """The spec registered under `format_id`."""
    validate_format(format_id)
    return next(spec for spec in REGISTRY if spec.id == format_id)
