"""Who made a VCF, read from its header: vendor, caller, caller version, and the ploidy statements.

Evidence comes from header lines only, never from a file name or a sample name, and a field the
header does not settle stays null (charter 9). Each rule in `_RULES` reads one kind of meta line;
the table it implements is RM2's header-evidence table in docs/ROADMAP.md.

Three kinds of header facts are kept apart, because real headers carry all three:

- **A caller** wrote the genotypes: HaplotypeCaller, DeepVariant, DRAGEN, freeBayes, an imputation
  tool. Two callers in one header is a conflict, and both fields stay null, with a warning naming
  both lines. A `vcf-merge` of several callers is not a conflict: the header says it is a merge.
- **A tool that touched the file** afterwards (GATK SelectVariants, ApplyVQSR, VariantFiltration,
  bcftools view or concat) says nothing about the caller and is not evidence for one.
- **A vendor** is the party that delivered the file, named by the header (`##center=Complete
  Genomics`, `##generatedby=Gencove`). A caller's maker is not a vendor: DRAGEN is Illumina's, but a
  DRAGEN file from Sequencing.com was delivered by Sequencing.com.

`callset_type` is set only where the header says what the content is: array genotypes or imputed
genotypes inside a VCF, where an absent site is unknown. Whole-genome and exome statements are left
to RM4's inference.

`ploidy_statements` are the header's ploidy lines verbatim, with no vocabulary (RM269 half 1): a meta
line whose key names a ploidy or karyotype, and each ploidy option inside a recorded command line,
prefixed with the line it came from. An option whose value is a path keeps its name and writes
`<path>` for the value: the file's content is not in the header, and paths can carry lab ids.
"""

import re
from collections.abc import Callable
from dataclasses import dataclass, field

from just_vcf.schema import CallsetType


@dataclass(frozen=True)
class Caller:
    """One caller the header names, with its version when the header states one."""

    name: str
    version: str | None
    line: str
    """The meta key the evidence came from, for warnings (`##source`, `##GATKCommandLine`)."""


@dataclass
class Source:
    """What the header says about who made the file."""

    vendor: str | None = None
    caller: str | None = None
    caller_version: str | None = None
    callset_type: CallsetType | None = None
    phasing_method: str | None = None
    ploidy_statements: list[str] = field(default_factory=list)
    callers: list[Caller] = field(default_factory=list)
    """Every caller named, one entry per distinct name, in header order."""
    merged: bool = False
    """The header records a `vcf-merge` of the callers' outputs."""
    warnings: list[str] = field(default_factory=list)


@dataclass
class _Evidence:
    callers: list[Caller] = field(default_factory=list)
    vendors: list[tuple[str, str]] = field(default_factory=list)
    """(vendor, meta key)."""
    content: list[tuple[CallsetType, str]] = field(default_factory=list)
    versions: dict[str, list[str]] = field(default_factory=dict)
    """`##source_version` values, applied to the one caller of a `##source` line."""
    phasing: list[str] = field(default_factory=list)
    merged: bool = False
    annotation_engine: bool = False
    """`##AnnotationEngineVersion`: Genos, when FORMAT also declares its `CGIANN_*` keys."""


_STRUCTURAL = ("INFO", "FORMAT", "FILTER", "ALT", "contig", "SAMPLE", "PEDIGREE", "META")

_GATK_CALLERS = frozenset({"HaplotypeCaller", "UnifiedGenotyper", "Mutect2", "MuTect2"})
"""GATK tools that write genotypes. Every other GATK tool (GenotypeGVCFs, SelectVariants, ApplyVQSR,
VariantFiltration, CombineVariants…) works on calls another tool made. GenotypeGVCFs genotypes a
gVCF HaplotypeCaller wrote, so it is the same caller's pipeline, not a second caller."""

_ISAAC = "Isaac Variant Caller"
"""Illumina's small-variant caller, also called starling: headers name it either way."""


def _attr(body: str, name: str) -> str | None:
    match = re.search(rf'(?:^|,|<){name}=("(?:[^"\\]|\\.)*"|[^,>]*)', body)
    return match.group(1).strip('"') if match else None


def _source(value: str, ev: _Evidence) -> None:
    """`##source=`: the many tools that write it, callers and not."""
    patterns: tuple[tuple[str, Callable[[re.Match], Caller | None]], ...] = (
        (r"^freeBayes\s+v?(\S+)", lambda m: Caller("freeBayes", m.group(1), "##source")),
        (r"^freeBayes$", lambda m: Caller("freeBayes", None, "##source")),
        (r"^VarScan2?$", lambda m: Caller("VarScan2", None, "##source")),
        (r"^CGAPipeline_(\S+)", lambda m: Caller("CGAPipeline", m.group(1), "##source")),
        (r"^(IsaacVariantCaller|starling)\b", lambda m: Caller(_ISAAC, None, "##source")),
        (r"^DRAGEN_\w+", lambda m: Caller("DRAGEN", None, "##source")),
        (r"^(HaplotypeCaller|UnifiedGenotyper|Mutect2)\b", lambda m: Caller(f"GATK {m.group(1)}", None, "##source")),
        (r"^Platypus_Version_(\S+)", lambda m: Caller("Platypus", m.group(1), "##source")),
        (r"^GenerateSVCandidates(?:\s+(\S+))?", lambda m: Caller("Manta", m.group(1), "##source")),
        (r"^aengine\b", lambda m: Caller("aengine", None, "##source")),
        (r"^loimpute\b", lambda m: Caller("loimpute", None, "##source")),
        (r"^Minimac(\d)\b", lambda m: Caller(f"Minimac{m.group(1)}", None, "##source")),
    )
    for pattern, make in patterns:
        match = re.search(pattern, value)
        if match:
            caller = make(match)
            if caller:
                ev.callers.append(caller)
            break
    else:
        if re.search(r"CASAVA", value) and re.search(r"makeConsensusVcf|VariantCalling", value):
            ev.callers.append(Caller("CASAVA", None, "##source"))
    if value.startswith("open_humans_data_processing.twenty_three_and_me"):
        ev.vendors.append(("23andMe", "##source"))
        ev.content.append((CallsetType.ARRAY, "##source"))
    elif value.startswith("intersection_dbSNP_23andme"):
        ev.vendors.append(("23andMe", "##source"))
        ev.content.append((CallsetType.ARRAY, "##source"))
    elif value.startswith("loimpute") or re.match(r"Minimac\d", value):
        ev.content.append((CallsetType.IMPUTED, "##source"))
    elif value.startswith("Sequencing.com"):
        ev.vendors.append(("Sequencing.com", "##source"))


def _gatk_command_line(key: str, body: str, ev: _Evidence) -> None:
    tool = _attr(body, "ID") or key.removeprefix("GATKCommandLine.").split(".")[0]
    if tool in _GATK_CALLERS:
        ev.callers.append(Caller(f"GATK {tool.replace('MuTect2', 'Mutect2')}", _attr(body, "Version"), f"##{key}"))


def _dragen(key: str, body: str, ev: _Evidence) -> None:
    if (_attr(body, "ID") or "dragen") != "dragen":
        return  # `ID=HashTableBuild` records the hash table's build, not the caller run.
    version = _attr(body, "Version")
    match = re.match(r"SW:\s*([^,\s]+)", version or "")
    ev.callers.append(Caller("DRAGEN", match.group(1) if match else version, f"##{key}"))


def _genes_for_good(value: str, ev: _Evidence) -> None:
    ev.vendors.append(("Genes for Good", "##Genes_for_Good_version"))
    # The vendor genotypes on an array; its header says which of its two deliveries this is.
    if re.search(r"\bimputed\b", value):
        ev.content.append((CallsetType.IMPUTED, "##Genes_for_Good_version"))
    elif re.search(r"\bgenotype calls\b", value):
        ev.content.append((CallsetType.ARRAY, "##Genes_for_Good_version"))


def _workflow_type(value: str, ev: _Evidence) -> None:
    if value.startswith("Illumina_"):
        ev.vendors.append(("Illumina", "##workflow_type"))
    if value == "Illumina_GenotypingToVCF":
        ev.content.append((CallsetType.ARRAY, "##workflow_type"))


_PLAIN_RULES: dict[str, Callable[[str, _Evidence], None]] = {
    "source": _source,
    "source_version": lambda v, ev: ev.versions.setdefault("source", []).append(v),
    "content": lambda v, ev: (
        ev.callers.append(Caller(_ISAAC, None, "##content")) if re.match(r"(IsaacVariantCaller|starling)\b", v) else None
    ),
    "center": lambda v, ev: ev.vendors.append((v.strip(), "##center")) if v.strip() else None,
    "workflow_type": _workflow_type,
    "DeepVariant_version": lambda v, ev: ev.callers.append(Caller("DeepVariant", v, "##DeepVariant_version")),
    "bcftools_callVersion": lambda v, ev: ev.callers.append(Caller("bcftools call", v, "##bcftools_callVersion")),
    "UnifiedGenotyper": lambda v, ev: ev.callers.append(Caller("GATK UnifiedGenotyper", None, "##UnifiedGenotyper")),
    "generatedby": lambda v, ev: ev.vendors.append(("Gencove", "##generatedby")) if v.startswith("Gencove") else None,
    "eagleVersion": lambda v, ev: ev.phasing.append(f"Eagle {v}"),
    "Genes_for_Good_version": _genes_for_good,
    "pipeline": lambda v, ev: ev.vendors.append(("Helix", "##pipeline")) if v.startswith("helix-") else None,
    "dataAnalysisProvider": lambda v, ev: ev.vendors.append((v.strip(), "##dataAnalysisProvider")) if v.strip() else None,
    "AnnotationEngineVersion": lambda v, ev: setattr(ev, "annotation_engine", True),
    "ae_version": lambda v, ev: ev.versions.setdefault("aengine", []).append(v),
}
"""Rules keyed by the exact meta key, for unstructured `##key=value` lines."""

_PLOIDY_KEY = re.compile(r"(?i)ploidy|karyotype")
_PLOIDY_OPTION = re.compile(
    r"(?<![\w-])("
    r"sample_ploidy=|--sample-ploidy[= ]|-ploidy[= ]|--ploidy[= ]|--ploidy-file[= ]|--ploidy-region-bed[= ]"
    r"|--sample-sex[= ]|--haploid[-_]contigs[= ]"
    r")\s*(\"[^\"]*\"|[^\s\",>]+)"
)


def ploidy_statements(meta: list[str]) -> list[str]:
    """The header's ploidy lines verbatim, and each ploidy option of a recorded command line."""
    out: list[str] = []
    for line in meta:
        key = line[2:].split("=", 1)[0]
        if key.split(".")[0] in _STRUCTURAL:
            continue
        if _PLOIDY_KEY.search(key):
            out.append(line)
            continue
        body = line[2 + len(key) + 1 :]
        tool = _attr(body, "ID") if body.startswith("<") else None
        where = f"{key} {tool}" if tool and tool not in key else key
        for match in _PLOIDY_OPTION.finditer(body):
            option, value = match.group(1), match.group(2)
            shown = "<path>" if re.search(r"[/\\]", value) else value
            out.append(f"{where}: {option}{shown}" if option.endswith("=") else f"{where}: {option.rstrip('= ')} {shown}")
    return out


def detect_source(meta: list[str], spec_vendor: str | None = None) -> Source:
    """What `meta` (a VCF's `##` lines, as `read_header` returns them) says about who made the file.

    `spec_vendor` is the vendor the detected format implies (`vcf_cgi`: Complete Genomics), used when
    no line names one. `generic` and `unknown` are not vendors.
    """
    ev = _Evidence()
    for line in meta:
        if not line.startswith("##") or "=" not in line:
            continue
        key, value = line[2:].split("=", 1)
        if re.match(r"source_\d{8}\.\d+$", key) and value.startswith("vcf-merge"):
            ev.merged = True
        elif key.startswith("GATKCommandLine"):
            _gatk_command_line(key, value, ev)
        elif key in ("DRAGENCommandLine", "DRAGENVersion"):
            _dragen(key, value, ev)
        elif key in _PLAIN_RULES:
            _PLAIN_RULES[key](value, ev)

    if ev.annotation_engine and any(line.startswith("##FORMAT=<ID=CGIANN_") for line in meta):
        ev.vendors.append(("Genos", "##AnnotationEngineVersion"))

    source = Source(ploidy_statements=ploidy_statements(meta), merged=ev.merged)
    source.callers = _distinct_callers(ev, source.warnings)
    names = [caller.name for caller in source.callers]
    if len(names) == 1:
        source.caller, source.caller_version = source.callers[0].name, source.callers[0].version
    elif len(names) > 1 and not ev.merged:
        lines = "; ".join(f"{c.name} from {c.line}" for c in source.callers)
        source.warnings.append(f"the header names more than one caller ({lines}); caller left null")
    fallback = spec_vendor if spec_vendor not in (None, "generic", "unknown") else None
    source.vendor = _one("vendor", ev.vendors, source.warnings) if ev.vendors else fallback
    content = _one("callset type", [(kind.value, key) for kind, key in ev.content], source.warnings)
    source.callset_type = CallsetType(content) if content else None
    source.phasing_method = ev.phasing[0] if len(set(ev.phasing)) == 1 else None
    return source


def _distinct_callers(ev: _Evidence, warnings: list[str]) -> list[Caller]:
    """One entry per caller name. Versions stated twice must agree, or the version is withheld."""
    by_name: dict[str, list[Caller]] = {}
    for caller in ev.callers:
        by_name.setdefault(caller.name, []).append(caller)
    out = []
    for name, seen in by_name.items():
        extra = ev.versions.get("source", []) if any(c.line == "##source" for c in seen) and len(by_name) == 1 else []
        if name == "aengine":
            extra = extra + ev.versions.get("aengine", [])
        versions = sorted({c.version for c in seen if c.version} | set(extra))
        if len(versions) > 1:
            warnings.append(f"the header gives {name} more than one version ({', '.join(versions)}); caller_version left null")
        out.append(Caller(name, versions[0] if len(versions) == 1 else None, seen[0].line))
    return out


def _one(what: str, found: list[tuple[str, str]], warnings: list[str]) -> str | None:
    values = sorted({value for value, _ in found})
    if len(values) > 1:
        lines = "; ".join(f"{value} from {key}" for value, key in found)
        warnings.append(f"the header names more than one {what} ({lines}); left null")
        return None
    return values[0] if values else None
