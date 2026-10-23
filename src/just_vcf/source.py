"""Who made a VCF and how, read from its header: vendor, caller, ploidy statements and the method chain.

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

The **method chain** (`method_chain`, charter 6) is every step the header records, in the order the
steps ran: alignment where a DRAGEN line says it mapped the reads, the call, and each tool that
worked on the calls afterwards (filtering, normalization, merging, liftover, annotation, phasing),
each with its stated version, a role from `schema.STEP_ROLES` and its key options. The reference
the rows are against comes from `##reference=`, or from the caller's command line. A step is read
from the line that names it, never from a file name. just-vcf's own shape (charter 2).
"""

import re
import shlex
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime

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
    """The one phasing step's tool and version (`Eagle 2.3.5`, `DRAGEN 05.021.408.3.4.11`); null when
    the header records none, or more than one, with a warning."""
    reference: str | None = None
    """The reference the header says the rows are against, path-redacted (`method_chain`)."""
    aligner: str | None = None
    aligner_version: str | None = None
    steps: list["Step"] = field(default_factory=list)
    """The method chain: every step the header records, in the order they ran (`method_chain`)."""
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
    merged: bool = False
    annotation_engine: bool = False
    """`##AnnotationEngineVersion`: Genos, when FORMAT also declares its `CGIANN_*` keys."""
    caller_at: dict[str, int] = field(default_factory=dict)
    """The index of the first meta line naming each caller, where its call step sits in the chain."""


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
    "clair3_version": lambda v, ev: ev.callers.append(Caller("Clair3", v, "##clair3_version")),
    "bcftools_callVersion": lambda v, ev: ev.callers.append(Caller("bcftools call", v, "##bcftools_callVersion")),
    "UnifiedGenotyper": lambda v, ev: ev.callers.append(Caller("GATK UnifiedGenotyper", None, "##UnifiedGenotyper")),
    "generatedby": lambda v, ev: ev.vendors.append(("Gencove", "##generatedby")) if v.startswith("Gencove") else None,
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
    for index, line in enumerate(meta):
        if not line.startswith("##") or "=" not in line:
            continue
        key, value = line[2:].split("=", 1)
        named = len(ev.callers)
        _read_line(key, value, ev)
        for caller in ev.callers[named:]:
            ev.caller_at.setdefault(caller.name, index)

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
    chain = method_chain(meta, source.callers, ev.caller_at)
    source.steps, source.reference = chain.steps, chain.reference
    source.warnings += chain.warnings
    source.phasing_method = _one_step("phasing", [s for s in chain.steps if s.role == "phase"], source.warnings)
    aligners = [s for s in chain.steps if s.role == "align"]
    if _one_step("aligner", aligners, source.warnings):
        source.aligner, source.aligner_version = aligners[0].tool, aligners[0].version
    return source


def _read_line(key: str, value: str, ev: _Evidence) -> None:
    """Apply the caller, vendor and content rules to one `##key=value` line."""
    if re.match(r"source_\d{8}\.\d+$", key) and value.startswith("vcf-merge"):
        ev.merged = True
    elif key.startswith("GATKCommandLine"):
        _gatk_command_line(key, value, ev)
    elif key in ("DRAGENCommandLine", "DRAGENVersion"):
        _dragen(key, value, ev)
    elif key in _PLAIN_RULES:
        _PLAIN_RULES[key](value, ev)


def _one_step(what: str, steps: list["Step"], warnings: list[str]) -> str | None:
    """`tool version` of the one distinct step, or None; more than one distinct step is a warning."""
    named = list(dict.fromkeys(" ".join(filter(None, (s.tool, s.version))) for s in steps))
    if len(named) > 1:
        warnings.append(f"the header records more than one {what} step ({'; '.join(named)}); left null")
        return None
    return named[0] if named else None


def _distinct_callers(ev: _Evidence, warnings: list[str]) -> list[Caller]:
    """One entry per caller name. Versions stated twice must agree, or the version is withheld.

    The versions are compared as written; only then is one written as a path or a file name kept as
    `<path>` (`_shown`), as a step's version is: `##DeepVariant_version=/home/<name>/…` names a home
    directory. Comparing after would make two different paths agree."""
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
            shown = ", ".join(_shown(version, "v") for version in versions)
            warnings.append(f"the header gives {name} more than one version ({shown}); caller_version left null")
        out.append(Caller(name, _shown(versions[0], "v") if len(versions) == 1 else None, seen[0].line))
    return out


def _one(what: str, found: list[tuple[str, str]], warnings: list[str]) -> str | None:
    values = sorted({value for value, _ in found})
    if len(values) > 1:
        lines = "; ".join(f"{value} from {key}" for value, key in found)
        warnings.append(f"the header names more than one {what} ({lines}); left null")
        return None
    return values[0] if values else None


# The method chain (charter 6: reference, aligner, caller, phasing method, post-processing).


@dataclass(frozen=True)
class Step:
    """One step the header records: a tool, its stated version, what it did, and its key options."""

    tool: str
    """The tool's name; a path or file name is written `<path>`, as an option's value is."""
    version: str | None
    """As stated; a path or file name is written `<path>`."""
    role: str
    """One of `schema.STEP_ROLES`."""
    options: tuple[str, ...] = ()
    """The options that say what the step did, as written (`-m -both`, `selectTypeToInclude=[SNP]`).
    A value that is a path or a file name is written `<path>`, except a reference or a chain, which
    keeps its base name where it is a known reference's (`_known_reference`), else `<file>`. Inputs,
    outputs and sample-naming options are never kept."""
    line: str = ""
    """The meta key the step was read from (`##bcftools_normCommand`, `##GATKCommandLine.SelectVariants.2`)."""

    def as_row(self) -> dict:
        return {"tool": self.tool, "version": self.version, "role": self.role, "options": list(self.options), "line": self.line}


@dataclass
class MethodChain:
    """What `method_chain` read: the ordered steps, the stated reference, and a warning per conflict."""

    steps: list[Step] = field(default_factory=list)
    reference: str | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class _Found:
    """A step as read: where it sits in the header, when it ran if the line says so, the references it names."""

    position: float
    step: Step
    dated: bool = False
    """The line states a run date, parsed or not."""
    when: datetime | None = None
    references: list[str] = field(default_factory=list)


_PATHLIKE = re.compile(
    r"[/\\]|\.(fa|fasta|fna|2bit|gz|bgz|zip|vcf|bcf|bam|cram|sam|bed|chain|txt|tsv|csv|json|dict|fai|list|intervals|config)$",
    re.IGNORECASE,
)
_EMPTY_VALUES = frozenset({"", "[]", "null", "false", "NA", "None"})
"""GATK lists every option with its default; these say the option was not used."""

_Spec = dict[str, str]
"""An option's name and how its value is kept: `v` as written (a path or file name becomes `<path>`),
`ref` the base name of a reference or chain, `flag` no value."""

_GATK_ROLES: dict[str, str] = {
    "GenotypeGVCFs": "call",
    "SelectVariants": "filter",
    "ApplyRecalibration": "filter",
    "ApplyVQSR": "filter",
    "VariantFiltration": "filter",
    "FilterVariantTranches": "filter",
    "CombineVariants": "merge",
    "CombineGVCFs": "merge",
    "GatherVcfs": "merge",
    "MergeVcfs": "merge",
    "VariantAnnotator": "annotate",
    "CNNScoreVariants": "annotate",
    "LeftAlignAndTrimVariants": "normalize",
    "ReadBackedPhasing": "phase",
}
"""GATK tools that work on calls another step made; any GATK tool not here or in `_GATK_CALLERS` is `other`."""

_GATK_REFERENCE: _Spec = {"--reference": "ref", "-R": "ref", "reference_sequence": "ref"}
_GATK_OPTIONS: dict[str, _Spec] = {
    "HaplotypeCaller": {**_GATK_REFERENCE, "--emit-ref-confidence": "v", "-ERC": "v", "emitRefConfidence": "v", "--dragen-mode": "v"},
    "UnifiedGenotyper": {**_GATK_REFERENCE, "genotype_likelihoods_model": "v", "-glm": "v", "output_mode": "v", "genotyping_mode": "v"},
    "Mutect2": dict(_GATK_REFERENCE),
    "MuTect2": dict(_GATK_REFERENCE),
    "GenotypeGVCFs": dict(_GATK_REFERENCE),
    "SelectVariants": {
        "selectTypeToInclude": "v", "selectTypeToExclude": "v", "excludeFiltered": "v", "excludeNonVariants": "v",
        "select_expressions": "v", "selectexpressions": "v", "selectExpressions": "v",
        "--select-type-to-include": "v", "-select-type": "v", "--select-type-to-exclude": "v", "-xl-select-type": "v",
        "--exclude-filtered": "v", "--exclude-non-variants": "v", "--selectExpressions": "v", "-select": "v",
    },
    "ApplyRecalibration": {"mode": "v", "ts_filter_level": "v"},
    "ApplyVQSR": {"--mode": "v", "-mode": "v", "--truth-sensitivity-filter-level": "v", "-ts-filter-level": "v"},
    "VariantFiltration": {
        "filterName": "v", "filterExpression": "v", "genotypeFilterName": "v", "genotypeFilterExpression": "v",
        "--filter-name": "v", "--filter-expression": "v", "--genotype-filter-name": "v", "--genotype-filter-expression": "v",
    },
    "CombineVariants": {"genotypemergeoption": "v", "filteredrecordsmergetype": "v"},
    "LeftAlignAndTrimVariants": {"--split-multi-allelics": "v"},
}

_REGIONS: _Spec = {"-r": "v", "--regions": "v", "-R": "v", "--regions-file": "v", "-t": "v", "--targets": "v", "-T": "v", "--targets-file": "v"}
_EXPRESSIONS: _Spec = {"-i": "v", "--include": "v", "-e": "v", "--exclude": "v"}
_BCFTOOLS: dict[str, tuple[str, _Spec]] = {
    "norm": ("normalize", {
        "-m": "v", "--multiallelics": "v", "-f": "ref", "--fasta-ref": "ref", "--old-rec-tag": "v", "-a": "flag",
        "--atomize": "flag", "-c": "v", "--check-ref": "v", "-d": "v", "--rm-dup": "v", "-N": "flag", "--do-not-normalize": "flag",
    }),
    "view": ("filter", {
        **_REGIONS, **_EXPRESSIONS, "-f": "v", "--apply-filters": "v", "-v": "v", "--types": "v", "-V": "v", "--exclude-types": "v",
        "-m": "v", "--min-alleles": "v", "-M": "v", "--max-alleles": "v", "-c": "v", "--min-ac": "v", "-C": "v", "--max-ac": "v",
        "-g": "v", "--genotype": "v", "-k": "flag", "--known": "flag", "-n": "flag", "--novel": "flag", "-p": "flag",
        "--phased": "flag", "-P": "flag", "--exclude-phased": "flag", "-a": "flag", "--trim-alt-alleles": "flag",
        "-U": "flag", "--exclude-uncalled": "flag",
    }),
    "filter": ("filter", {
        **_REGIONS, **_EXPRESSIONS, "-s": "v", "--soft-filter": "v", "-m": "v", "--mode": "v", "-g": "v", "--SnpGap": "v",
        "-G": "v", "--IndelGap": "v", "-S": "v", "--set-GTs": "v",
    }),
    "annotate": ("annotate", {
        **_REGIONS, **_EXPRESSIONS, "-a": "v", "--annotations": "v", "-c": "v", "--columns": "v", "-x": "v", "--remove": "v",
        "-I": "v", "--set-id": "v",
    }),
    "concat": ("merge", {
        **_REGIONS, "-a": "flag", "--allow-overlaps": "flag", "-D": "flag", "--remove-duplicates": "flag", "-d": "v", "--rm-dups": "v",
        "-n": "flag", "--naive": "flag",
    }),
    "merge": ("merge", {
        **_REGIONS, "-m": "v", "--merge": "v", "-i": "v", "--info-rules": "v", "-f": "v", "--apply-filters": "v", "-0": "flag",
        "--missing-to-ref": "flag",
    }),
    "isec": ("merge", {
        **_REGIONS, **_EXPRESSIONS, "-n": "v", "--nfiles": "v", "-w": "v", "--write": "v", "-c": "v", "--collapse": "v",
        "-f": "v", "--apply-filters": "v", "-C": "flag", "--complement": "flag",
    }),
    "call": ("call", {
        "-m": "flag", "--multiallelic-caller": "flag", "-c": "flag", "--consensus-caller": "flag", "-v": "flag",
        "--variants-only": "flag", "-g": "v", "--gvcf": "v", "-A": "flag", "--keep-alts": "flag", "-P": "v", "--prior": "v",
    }),
    "mpileup": ("call", {
        **_REGIONS, "-f": "ref", "--fasta-ref": "ref", "-q": "v", "--min-MQ": "v", "-Q": "v", "--min-BQ": "v", "-a": "v",
        "--annotate": "v", "-g": "v", "--gvcf": "v", "-B": "flag", "--no-BAQ": "flag", "-C": "v", "--adjust-MQ": "v",
    }),
    "csq": ("annotate", {"-f": "ref", "--fasta-ref": "ref", "-g": "v", "--gff-annot": "v"}),
    "convert": ("other", {"--tsv2vcf": "flag", "--gvcf2vcf": "flag", "--vcf2gvcf": "flag", "-f": "ref", "--fasta-ref": "ref"}),
}
"""bcftools subcommands: the role and the options kept. `view` is a filter only when a filtering option
is present; otherwise it only rewrote the file (`other`). Any other subcommand is `other`. Ploidy
options are in `ploidy_statements`, not here."""

_VCFTOOLS_ROLES = {"vcf-merge": "merge", "vcf-concat": "merge", "vcf-isec": "merge", "vcf-annotate": "annotate"}
_CALLER_COMMANDS: dict[str, tuple[str, _Spec]] = {
    "freebayes": ("freeBayes", {"-f": "ref", "--fasta-reference": "ref"}),
    "starling": (_ISAAC, {"-samtools-reference": "ref", "--do-short-range-phasing": "flag"}),
    "starling2": (_ISAAC, {"-samtools-reference": "ref", "--do-short-range-phasing": "flag"}),
    "run_clair3.sh": ("Clair3", {"--ref_fn": "ref", "--platform": "v", "--gvcf": "flag"}),
}
"""`##commandline=`/`##cmdline=` lines whose program is a caller: their options go to its call step."""

_DRAGEN_CALL: _Spec = {
    "--ref-dir": "ref", "--enable-variant-caller": "v", "--vc-emit-ref-confidence": "v", "--vc-enable-phasing": "v",
    "--enable-cnv": "v", "--enable-sv": "v",
}
_DRAGEN_ALIGN: _Spec = {"--ref-dir": "ref", "--enable-map-align": "v", "--enable-duplicate-marking": "v"}

_CALLER_PHASE_SETS: dict[str, str] = {
    "DRAGEN": "Physical phasing ID information",
    "GATK HaplotypeCaller": "Phasing set (typically the position of the first variant in the set)",
    "DeepVariant": "Phase set",
    "CGAPipeline": "Phase Set",
}
"""The FORMAT `PS` description each caller writes when it phases its own calls (read from the PGP-HMS
and local headers, 2026-10-23). DRAGEN's is a prefix, the rest are whole. A `PS` declared in the
caller's own words means the caller wrote the phase sets."""

_WHATSHAP = re.compile(r'^"?\(whatshap ([^)\s]+)\)\s+(\S+)(.*?)"?$')
"""WhatsHap's `##commandline=`: `(whatshap <version>) <subcommand> <args>` (whatshap/cli/phase.py)."""
_SHAPEIT = re.compile(r"^shapeit(\d)\.(\S+)$")
"""SHAPEIT4's `##source=shapeit4.1.3` (shapeit4 src/io/haplotype_writer.cpp)."""
_HAPCUT2 = "##INFO=<ID=hapcut2,"
"""The INFO key HapCUT2 tags the records it phased with (HapCUT2 hapcut2-src/output_phasedvcf.c)."""
_HEADER_VERSIONS = ("SnpEffVersion", "OriginalSnpEffVersion", "gvcftools_version", "eagleVersion")


def method_chain(meta: list[str], callers: list[Caller], caller_at: dict[str, int] | None = None) -> MethodChain:
    """The steps `meta` records, in the order they ran, and the reference the rows are against.

    `callers` are `detect_source`'s, one call step each, at its first line (`caller_at`), with the
    options of its command lines. A DRAGEN line that says it mapped the reads adds an `align` step; a
    caller's FORMAT `PS` in its own words adds a `phase` step after its call. Every later tool the
    header names (GATK, bcftools, vcftools, SnpEff, VEP, gvcftools, CrossMap, a phaser) is a step.

    Order: header order, except that the steps whose line states a run date (GATK, DRAGEN, bcftools)
    are sorted by that date among the places they hold: htsjdk writes meta lines sorted by their text,
    which puts `ApplyRecalibration` before `HaplotypeCaller`. Dates are compared without their zone.
    When a stated date does not parse, the order is the header's, with a warning.

    `reference` is `##reference=`'s base name (DRAGEN's `reference.bin` gives its directory's); with
    none, the one reference the call steps name. Two of them leave it null with a warning; a
    `##reference=` that differs from the caller's with no `lift` step recorded is kept, with a warning.
    """
    caller_at = caller_at or {}
    chain = MethodChain()
    found: list[_Found] = []
    calls = {caller.name: _Found(float(caller_at.get(caller.name, 0)), Step(caller.name, caller.version, "call", (), caller.line)) for caller in callers}
    stated: list[str] = []
    gatk_named: set[str] = set()
    gatk_sources: list[tuple[int, str]] = []
    phases_by_caller: list[str] = []
    versions = {key: value for key, value in (_meta_pair(line) for line in meta) if key in _HEADER_VERSIONS}

    def add(item: _Found, caller: str | None = None) -> None:
        """A step, or, for a caller's own command line, its options added to the caller's step."""
        call = calls.get(caller) if caller else None
        if call is None:
            found.append(item)
            return
        options = call.step.options + tuple(o for o in item.step.options if o not in call.step.options)
        call.step = Step(call.step.tool, call.step.version, "call", options, call.step.line)
        call.position = min(call.position, item.position)
        if item.dated and (not call.dated or (item.when and call.when and item.when < call.when)):
            call.dated, call.when = True, item.when
        call.references += [r for r in item.references if r not in call.references]

    for index, line in enumerate(meta):
        if line.startswith(_HAPCUT2):
            add(_Found(index, Step("HapCUT2", None, "phase", (), "##INFO=<ID=hapcut2>")))
            continue
        if line.startswith("##FORMAT=<ID=PS,"):
            described = _attr(line[len("##FORMAT=") :], "Description") or ""
            phases_by_caller += [
                name for name, words in _CALLER_PHASE_SETS.items()
                if name in calls and (described == words or (name == "DRAGEN" and described.startswith(words)))
            ]
            continue
        key, value = _meta_pair(line)
        if not key or key.split(".")[0] in _STRUCTURAL:
            continue
        if key.startswith("bcftools") and key.endswith("Version"):
            versions[key.removesuffix("Version")] = value.strip()
        elif key == "reference":
            stated.append(_reference_name(value))
        elif key.startswith("GATKCommandLine") or _gatk2(key, value):
            tool, version, text = _gatk_parts(key, value)
            gatk_named.add(tool)
            spec = _GATK_OPTIONS.get(tool, {})
            options = _options(text, spec)
            dated, when = _dated(value)
            role = "call" if tool in _GATK_CALLERS else _GATK_ROLES.get(tool, "other")
            step = Step(f"GATK {tool.replace('MuTect2', 'Mutect2')}", version, role, tuple(options), f"##{key}")
            add(_Found(index, step, dated, when, _references(options, spec)), step.tool if tool in _GATK_CALLERS else None)
        elif key == "source" and value.strip() in _GATK_ROLES:
            gatk_sources.append((index, value.strip()))
        elif key == "source" and (match := _SHAPEIT.match(value.strip())):
            add(_Found(index, Step(f"SHAPEIT{match.group(1)}", match.group(2), "phase", (), "##source")))
        elif key == "DRAGENCommandLine" and (_attr(value, "ID") or "dragen") == "dragen":
            version, text = _dragen_version(value), _attr(value, "CommandLineOptions") or ""
            dated, when = _dated(value)
            align, call = _options(text, _DRAGEN_ALIGN), _options(text, _DRAGEN_CALL)
            add(_Found(index, Step("DRAGEN", version, "call", tuple(call), f"##{key}"), dated, when, _references(call, _DRAGEN_CALL)), "DRAGEN")
            if "--enable-map-align true" in align:
                add(_Found(index - 0.5, Step("DRAGEN", version, "align", tuple(align), f"##{key}"), dated, when))
        elif re.fullmatch(r"bcftools_?\w*Command", key):
            item = _bcftools(index, key, value, versions)
            add(item, "bcftools call" if item.step.tool == "bcftools call" else None)
        elif re.fullmatch(r"source_\d{8}\.\d+", key):
            program, version, rest = _program(value)
            if program:
                options = _options(rest, {"-c": "v", "--columns": "v"}) if program == "vcf-annotate" else []
                add(_Found(index, Step(program, version, _VCFTOOLS_ROLES.get(program, "other"), tuple(options), f"##{key}")))
        elif key == "commandline" and (match := _WHATSHAP.match(value.strip())):
            sub = match.group(2)
            spec = {"--reference": "ref", "-r": "ref", "--indels": "flag", "--distrust-genotypes": "flag"}
            role = "phase" if sub in ("phase", "polyphase") else "call" if sub == "genotype" else "other"
            add(_Found(index, Step("WhatsHap", match.group(1), role, (sub, *_options(match.group(3), spec)), "##commandline")))
        elif key in ("commandline", "cmdline"):
            program, _, rest = _program(value)
            if program in _CALLER_COMMANDS:
                name, spec = _CALLER_COMMANDS[program]
                options = _options(rest, spec)
                add(_Found(index, Step(name, None, "call", tuple(options), f"##{key}"), references=_references(options, spec)), name)
        elif key in ("SnpEffCmd", "OriginalSnpEffCmd"):
            version = versions.get(key.removesuffix("Cmd") + "Version")
            add(_Found(index, Step("SnpEff", version, "annotate", tuple(_snpeff_options(value)), f"##{key}")))
        elif key == "VEP":
            words = value.strip().split()
            options = [w for w in words[1:] if re.match(r"(assembly|db|ensembl)=", w) and not _PATHLIKE.search(w.split("=", 1)[1])]
            add(_Found(index, Step("VEP", words[0].strip('"') if words else None, "annotate", tuple(options), "##VEP")))
        elif key == "gvcftools_cmdline":
            program, _, rest = _program(value)
            options = _options(rest, {"--ref": "ref"})
            tool = f"gvcftools {program}" if program else "gvcftools"
            add(_Found(index, Step(tool, versions.get("gvcftools_version"), "other", tuple(options), f"##{key}")))
        elif key in ("eagleVersion", "eagleCommand") and not any(item.step.tool == "Eagle" for item in found):
            add(_Found(index, Step("Eagle", versions.get("eagleVersion"), "phase", (), f"##{key}")))
        elif key == "liftOverProgram":
            program = re.match(r"<?([^,>]+)", value)
            options = [f"{k}={_known_reference(_basename(v))}" for k, v in map(_meta_pair, meta) if k in ("liftOverChainFile", "targetRefGenome")]
            add(_Found(index, Step(program.group(1) if program else value, _attr(value, "version"), "lift", tuple(options), f"##{key}")))

    for index, tool in gatk_sources:
        if tool not in gatk_named:
            add(_Found(index, Step(f"GATK {tool}", None, _GATK_ROLES[tool], (), "##source")))
    found += calls.values()
    for name in dict.fromkeys(phases_by_caller):
        # The caller's own phase sets: right after its call step, dated like it.
        call = calls[name]
        found.append(_Found(call.position + 0.25, Step(name, call.step.version, "phase", (), "##FORMAT=<ID=PS>"), call.dated, call.when))

    chain.steps = [_without_paths(item.step) for item in _ordered(found, chain.warnings)]
    called = list(dict.fromkeys(ref for item in found if item.step.role == "call" for ref in item.references))
    chain.reference = _reference(stated, called, chain.steps, chain.warnings)
    return chain


def _without_paths(step: Step) -> Step:
    """`step` with a tool or version that is a path or a file name written `<path>`, as an option's
    value is: `##liftOverProgram=</home/<name>/bin/CrossMap,…>` names a home directory."""
    version = None if step.version is None else _shown(step.version, "v")
    return Step(_shown(step.tool, "v"), version, step.role, step.options, step.line)


def _meta_pair(line: str) -> tuple[str, str]:
    if not line.startswith("##") or "=" not in line:
        return "", ""
    key, value = line[2:].split("=", 1)
    if key in _HEADER_VERSIONS:
        words = value.strip().strip('"').split()
        value = words[0] if words else ""
    return key, value


def _ordered(found: list[_Found], warnings: list[str]) -> list[_Found]:
    """Header order, with the dated steps sorted by date among the places they hold."""
    by_line = sorted(found, key=lambda item: item.position)
    dated = [item for item in by_line if item.dated]
    if any(item.when is None for item in dated):
        warnings.append("a step's stated date did not parse; the method chain is in header order")
        return by_line
    by_date = iter(sorted(dated, key=lambda item: item.when or datetime.min))
    return [next(by_date) if item.dated else item for item in by_line]


def _reference(stated: list[str], called: list[str], steps: list[Step], warnings: list[str]) -> str | None:
    distinct = list(dict.fromkeys(name for name in stated if name))
    if len(distinct) > 1:
        warnings.append(f"the header names more than one ##reference ({', '.join(distinct)}); reference left null")
        return None
    if not distinct:
        if len(called) > 1:
            warnings.append(f"the call steps name more than one reference ({', '.join(called)}); reference left null")
            return None
        return called[0] if called else None
    if called and distinct[0] not in called and not any(step.role == "lift" for step in steps):
        warnings.append(
            f"##reference names {distinct[0]} but the call step names {', '.join(called)}, and no liftover step is "
            "recorded: the rows may have been moved to another reference by a tool the header does not name"
        )
    return distinct[0]


def _reference_name(value: str) -> str:
    """`##reference=`'s base name (`_known_reference`); DRAGEN names the `reference.bin` inside its
    reference directory."""
    text = _strip_url(value)
    parts = [part for part in re.split(r"[/\\]", text) if part]
    if len(parts) > 1 and parts[-1] == "reference.bin":
        return _known_reference(parts[-2])
    return _known_reference(parts[-1]) if parts else ""


FILE_PLACEHOLDER = "<file>"
"""What a reference or chain file's base name is written as when it is not a known reference's."""

_ASSEMBLY_WORD = re.compile(
    r"hg(?:18|19|38)|grch3[678]|hs37d5|hs38(?:dh)?|hs1|b3[678]|build3[678]|ncbi3[678]|g1k|1kg|v37|assembly(?:19|38)"
    r"|chm13(?:v\d+)?|t2t|gca|gcf|(?:hg\d+|grch3\d)to(?:hg\d+|grch3\d|chm13)",
    re.IGNORECASE,
)
"""A word that names an assembly in a reference or chain file's name: UCSC's, GRC's, the 1000 Genomes
and GATK bundles', NCBI's accessions (`GCA_000001405`, T2T's `GCA_009914755`) and UCSC's chains."""
_REFERENCE_WORD = re.compile(
    r"\d+|[prv]\d+|ucsc|human|homo|sapiens|assembly|primary|toplevel|dna|sm|rm|genome|genomic|nt|alt|aware|masked"
    r"|decoy|hla|nohla|noalt|no|analysis|set|full|plus|ebv|mt|chrm|rcrs|cnv|graph|rna|to|over|chain|ref|reference"
    r"|1000genomes|fa|fasta|fna|2bit|gz|bgz|bz2|bin|dict|fai",
    re.IGNORECASE,
)
"""The other words such a name is made of: versions, the bundles' descriptors, extensions."""


def _known_reference(name: str) -> str:
    """`name` (a reference or chain file's base name) if it is a known reference's, else `<file>`.

    Known: every word between `.`, `_`, `-` and `+` is an assembly word or a descriptor
    (`_ASSEMBLY_WORD`, `_REFERENCE_WORD`), and one names an assembly: `hg38.fa`, `hs37d5.fa.gz`,
    `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, `GRCh38_to_GRCh37.chain.gz`. Any other base name
    may carry a person's or a lab's name (`ALICESMITH-hg38.fa`); one made of known words cannot."""
    words = [word for word in re.split(r"[._+\-]+", name) if word]
    if not words:
        return name
    known = all(_ASSEMBLY_WORD.fullmatch(word) or _REFERENCE_WORD.fullmatch(word) for word in words)
    return name if known and any(_ASSEMBLY_WORD.fullmatch(word) for word in words) else FILE_PLACEHOLDER


def _strip_url(value: str) -> str:
    return re.sub(r"^file:/*", "/", value.strip().strip('"').strip("<>"))


def _basename(value: str) -> str:
    parts = [part for part in re.split(r"[/\\]", _strip_url(value)) if part]
    return parts[-1] if parts else ""


def _gatk2(key: str, value: str) -> bool:
    """A GATK 1 or 2 command line, written under the tool's own key (`##SelectVariants="analysis_type=SelectVariants …"`)."""
    return value.lstrip('"').startswith(f"analysis_type={key} ")


def _gatk_parts(key: str, value: str) -> tuple[str, str | None, str]:
    """The tool, its version and its options text, in each of GATK's three header styles."""
    if not key.startswith("GATKCommandLine"):
        return key, None, value.strip('"')
    tool = _attr(value, "ID") or key.removeprefix("GATKCommandLine.").split(".")[0]
    text = _attr(value, "CommandLine") or _attr(value, "CommandLineOptions") or ""
    return tool, _attr(value, "Version"), text.replace('\\"', '"')


def _dragen_version(value: str) -> str | None:
    version = _attr(value, "Version")
    match = re.match(r"SW:\s*([^,\s]+)", version or "")
    return match.group(1) if match else version


def _bcftools(index: int, key: str, value: str, versions: dict[str, str]) -> _Found:
    """A `##bcftools_<sub>Command` (or `##bcftoolsCommand`) line, with the version line before it."""
    command, _, date = value.partition("; Date=")
    words = command.split()
    sub = words[0] if words else key.removeprefix("bcftools_").removesuffix("Command")
    role, spec = _BCFTOOLS.get(sub, ("other", {}))
    options = _options(command.strip()[len(sub) :], spec)
    if sub == "view" and not options:
        role = "other"
    version = versions.get(f"bcftools_{sub}") or versions.get("bcftools")
    step = Step(f"bcftools {sub}", version, role, tuple(options), f"##{key}")
    return _Found(index, step, bool(date.strip()), _parse_date(date) if date.strip() else None, _references(options, spec))


def _program(value: str) -> tuple[str | None, str | None, str]:
    """A command line's program base name, the `(rNNN)` revision vcftools appends to it, and the rest."""
    first, _, rest = value.strip().strip('"').partition(" ")
    match = re.fullmatch(r"([\w.+-]+?)(?:\((r\d+)\))?", re.split(r"[/\\]", first)[-1])
    if not match or re.search(r"\.(vcf|gz|txt|bam|cram)$", match.group(1)):
        return None, None, rest
    return match.group(1), match.group(2), rest


_SNPEFF_KEPT: _Spec = {
    "-onlyCoding": "v", "-ud": "v", "-upDownStreamLen": "v", "-noStats": "flag", "-no-intergenic": "flag",
    "-no-downstream": "flag", "-no-upstream": "flag", "-no-intron": "flag", "-canon": "flag", "-lof": "flag", "-hgvs": "flag",
}
_SNPEFF_SKIPPED = frozenset({"-config", "-c", "-o", "-i", "-s", "-stats", "-dataDir"})
"""SnpEff options whose value is a file or a format, skipped with their value."""


def _snpeff_options(value: str) -> list[str]:
    """SnpEff's options that change what it annotates, then its genome database (the first bare word)."""
    words = value.strip().strip('"').split()[1:]
    if words and words[0] in ("eff", "ann"):
        words = words[1:]
    out: list[str] = []
    i = 0
    while i < len(words):
        word = words[i]
        i += 1
        if word in _SNPEFF_SKIPPED:
            i += 1
        elif _SNPEFF_KEPT.get(word) == "flag":
            out.append(word)
        elif word in _SNPEFF_KEPT and i < len(words):
            out.append(f"{word} {_shown(words[i], 'v')}")
            i += 1
        elif not word.startswith("-") and not _PATHLIKE.search(word):
            out.append(word)
            break
    return out


def _options(text: str, spec: _Spec) -> list[str]:
    """The options `spec` names, in command order, written as they were with their values redacted.

    Reads `--name value`, `--name=value`, `-xVALUE` and GATK 3's `name=value` and `name=[a, b]`."""
    if not spec or not text:
        return []
    out: list[str] = []
    for match in re.finditer(r"(?:^|\s)(\w+)=(\[[^\]]*\]|\"[^\"]*\"|\S+)", text):
        name, value = match.group(1), match.group(2)
        if name in spec and value.strip('"') not in _EMPTY_VALUES:
            out.append(f"{name}={_shown(value, spec[name])}")
    tokens = _tokens(text)
    i = 0
    while i < len(tokens):
        token = tokens[i]
        i += 1
        name, sep, value = token, " ", None
        if token.startswith("--") and "=" in token:
            name, value = token.split("=", 1)
            sep = "="
        elif token.startswith("-") and not token.startswith("--") and token not in spec and spec.get(token[:2], "flag") != "flag":
            name, value, sep = token[:2], token[2:], ""
        kind = spec.get(name) if name.startswith("-") else None
        if kind == "flag":
            out.append(name)
            continue
        if kind is None:
            continue
        if value is None:
            if i >= len(tokens):
                continue
            value = tokens[i]
            i += 1
        if value not in _EMPTY_VALUES:
            out.append(f"{name}{sep}{_shown(value, kind)}")
    return list(dict.fromkeys(out))


def _tokens(text: str) -> list[str]:
    try:
        return shlex.split(text)
    except ValueError:
        return text.split()


def _shown(value: str, kind: str) -> str:
    """A value as a step keeps it: a reference's base name, `<path>` for any other path or file."""
    if kind == "ref":
        return _known_reference(_basename(value))
    # `INFO/DP` and `FORMAT/GQ` in a column list or an expression are field names, not paths.
    return "<path>" if _PATHLIKE.search(re.sub(r"\b(?:INFO|FORMAT|FMT)/", "", value.strip('"'))) else value


def _references(options: list[str], spec: _Spec) -> list[str]:
    """The reference base names among a step's options."""
    out = []
    names = [name for name, kind in spec.items() if kind == "ref"]
    for option in options:
        for name in names:
            rest = option[len(name) :]
            if option.startswith(name) and rest and (rest[0] in " =" or len(name) == 2):
                out.append(rest[1:] if rest[0] in " =" else rest)
                break
    return out


def _dated(value: str) -> tuple[bool, datetime | None]:
    date = _attr(value, "Date")
    return (True, _parse_date(date)) if date else (False, None)


def _parse_date(text: str) -> datetime | None:
    """A run date as GATK 3, GATK 4 (`July 13, 2020 4:33:39 PM HKT`, from 4.2 `… 2020 at 4:33:39 PM`),
    DRAGEN and bcftools write it, without its zone."""
    text = re.sub(r"(?<=\d{4}) at (?=\d)", " ", text.strip())
    text = re.sub(r"\s+[A-Z]{2,5}(?=\s+\d{4}$)", "", text)
    text = re.sub(r"(?<=[AP]M)\s+[A-Z]{2,5}$", "", text)
    for fmt in ("%a %b %d %H:%M:%S %Y", "%B %d, %Y %I:%M:%S %p", "%b %d, %Y %I:%M:%S %p"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None
