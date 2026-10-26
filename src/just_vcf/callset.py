"""What kind of callset a VCF facet is, and so what an absent position means (RM4; charter 9, 11).

The scope test is ported from just-dna-lite `annotation/restoration.py` `detect_callset_scope`
(docs/PROVENANCE.md): a whole-genome callset carries millions of sites, and nearly all of its span
lies within one flank of a call. That file separated whole-genome from "targeted" only. Here the
targeted side is split into exome and panel, and every band between two classes answers unknown, so
weak evidence never lands in the nearest class.

The classes, by what the calls look like:

- `wgs`: at least `MIN_WGS_SITES` sites on the primary assembly and at least `MIN_WGS_BREADTH` of the
  span within `FLANK_BP` of a call (the ported rule, thresholds unchanged).
- `wes`: an exome's spread. Between `MIN_WES_SITES` and `MAX_WES_SITES` sites, on every autosome, and
  clustered (breadth at most `MAX_TARGETED_BREADTH`).
- `panel`: calls confined to a small region set. Between `MIN_PANEL_SITES` and `MAX_PANEL_SITES` sites,
  all within `MAX_PANEL_REACH_BP` of genome. A region extract of a larger callset (`bcftools view -r`)
  looks the same and cannot be told apart from the calls: either way absence outside the regions means
  nothing was looked at, which is what the facet then says.
- `imputed`: the header names an imputation tool or declares an imputation INFO key. Imputed sites are
  a panel's sites, whatever their spread.

Each sequencing class also needs the header to declare a read-evidence FORMAT key (AD, DP, GQ, PL, GL,
MIN_DP): a genome-wide file of genotypes with no read evidence at all is what imputation and array
conversions write, so its spread alone is not taken as a whole-genome callset.
"""

import re
from collections.abc import Mapping
from enum import Enum
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.schema import (
    CONTENT_FILTER_KEEPS,
    AbsenceMeans,
    CallsetType,
    ContentFilter,
    ContentFilterBasis,
    VariantClass,
)

FLANK_BP = 10_000
"""just-dna-lite's `restoration_max_flank_bp` default, the flank the breadth is measured at."""

MIN_WGS_SITES = 1_000_000
"""just-dna-lite `MIN_WGS_SITES`: about 4x below the WGS floor it measured (4.3–4.7M sites)."""

MIN_WGS_BREADTH = 0.75
"""just-dna-lite `MIN_WGS_BREADTH`: WGS measured 0.942–0.950 there, a clustered callset 0.21."""

MIN_WES_SITES = 10_000
MAX_WES_SITES = 400_000
"""An exome carries roughly 20,000–150,000 sites, more with off-target calls; the bounds sit below and
above that with room on both sides. 400,000 to 1,000,000 is unknown."""

MAX_TARGETED_BREADTH = 0.6
"""An exome's breadth stays well under a genome's 0.94; 0.6 to 0.75 is unknown."""

MIN_PANEL_SITES = 50
MAX_PANEL_SITES = 5_000
"""Fewer than 50 sites (a test file, a single-gene extract) say too little to name a class; 5,000 to
10,000 is the unknown band between a large panel and a small exome."""

MAX_PANEL_REACH_BP = 50_000_000
"""Genome within `FLANK_BP` of a panel's calls. A few hundred genes, each with its flanks, stay far
below it; an exome's calls reach hundreds of megabases."""

READ_EVIDENCE_KEYS = frozenset({"AD", "DP", "GQ", "PL", "GL", "MIN_DP"})
"""FORMAT keys a caller that read the sample's reads writes; imputation and array conversions do not."""

IMPUTATION_INFO_KEYS = frozenset({"R2", "ER2", "DR2", "IMP", "IMPUTED", "TYPED"})
"""INFO keys minimac (R2, ER2, IMPUTED, TYPED) and Beagle (DR2, IMP) declare."""

_IMPUTATION_SOURCE = re.compile(r"^##source=.*(?i:minimac|beagle|impute[245]?|glimpse|loimpute)", re.MULTILINE)
"""A `##source` line naming an imputation tool (Gencove writes `loimpute`)."""

_PRIMARY = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]
_AUTOSOMES = frozenset(str(n) for n in range(1, 23))


class CallsetInference(BaseModel):
    """The callset type the evidence supports, or why there is none. Truthy only with a type."""

    model_config = ConfigDict(frozen=True)

    callset_type: CallsetType | None
    reason: str
    sites: int | None = None
    """Distinct `(chrom, pos)` on the primary assembly; None when the spread was not measured."""
    breadth: float | None = None
    reach_bp: int | None = None
    autosomes: int | None = None

    def __bool__(self) -> bool:
        return self.callset_type is not None


def imputation_marker(meta: list[str], info_keys: list[str]) -> str | None:
    """What in the header says the file is imputation output, or None."""
    declared = sorted(IMPUTATION_INFO_KEYS & set(info_keys))
    if declared:
        return f"the header declares imputation INFO key(s) {', '.join(declared)}"
    for line in meta:
        if _IMPUTATION_SOURCE.match(line):
            return f"the header names an imputation tool ({line[:60]})"
    return None


def spread(calls: pl.LazyFrame, flank_bp: int = FLANK_BP) -> dict:
    """Site count, breadth, reach and autosome count of a `calls` frame's small-variant sites.

    A site is a distinct `(chrom, pos)` on the primary assembly (canonical `chrom`) of a row that
    writes a sequence allele, or none (a monomorphic record). Rows whose ALTs are all symbolic
    (`<DEL>`, `<*>`, breakends, `*`) do not count: a genome-wide structural-variant callset says
    nothing about small variants, and must not read as a whole-genome one. Alt contigs and decoys do
    not count either; they carry a handful of calls each.

    Breadth is just-dna-lite's `_coverage_breadth`: each gap between neighbouring sites counts up to
    `2 * flank_bp`, summed (the reach), divided by the span the sites spread over per contig.
    """
    # The empty allele is a sequence allele (a deletion written without its anchor base, F12).
    sequence_allele = pl.col("alts").list.eval(pl.element().str.contains(r"^[ACGTNacgtn]*$")).list.any()
    primary = (
        calls.filter(pl.col("chrom").is_in(_PRIMARY) & ((pl.col("alts").list.len() == 0) | sequence_allele))
        .select("chrom", "pos")
        .unique()
    )
    per_contig = (
        primary.group_by("chrom")
        .agg(
            pl.len().alias("n"),
            (pl.col("pos").max().cast(pl.Int64) - pl.col("pos").min().cast(pl.Int64)).alias("span"),
            pl.col("pos").sort().cast(pl.Int64).diff().clip(upper_bound=2 * flank_bp).sum().alias("reach"),
        )
        .collect()
    )
    total = int(per_contig["n"].sum() or 0)
    span = int(per_contig["span"].sum() or 0)
    reach = int(per_contig["reach"].sum() or 0)
    return {
        "sites": total,
        "breadth": reach / span if span else 0.0,
        "reach_bp": reach,
        "autosomes": int(per_contig.filter(pl.col("chrom").is_in(_AUTOSOMES)).height),
    }


def classify(calls: pl.LazyFrame, format_keys: list[str], imputed_by: str | None) -> CallsetInference:
    """The callset type of a VCF facet from its `calls` (`chrom`, `pos`, `alts`) and its header."""
    if imputed_by:
        return CallsetInference(callset_type=CallsetType.IMPUTED, reason=imputed_by)
    measured = spread(calls)
    n, breadth, reach, autosomes = measured["sites"], measured["breadth"], measured["reach_bp"], measured["autosomes"]
    where = f"{n:,} site(s) on the primary assembly, {breadth:.1%} of the span within {FLANK_BP:,} bp of a call, {reach:,} bp reached, {autosomes} autosome(s)"

    def answer(callset_type: CallsetType | None, why: str) -> CallsetInference:
        return CallsetInference(callset_type=callset_type, reason=f"{why}: {where}", **measured)

    evidence = sorted(READ_EVIDENCE_KEYS & set(format_keys))
    if not evidence:
        return answer(None, "unknown, the header declares no read-evidence FORMAT key (AD, DP, GQ, PL, GL, MIN_DP), so the calls may not come from reads")
    if n >= MIN_WGS_SITES and breadth >= MIN_WGS_BREADTH:
        return answer(CallsetType.WGS, "whole genome")
    if MIN_WES_SITES <= n <= MAX_WES_SITES and breadth <= MAX_TARGETED_BREADTH and autosomes == len(_AUTOSOMES):
        return answer(CallsetType.WES, "exome")
    if MIN_PANEL_SITES <= n <= MAX_PANEL_SITES and reach <= MAX_PANEL_REACH_BP:
        return answer(CallsetType.PANEL, "panel or region extract")
    return answer(None, "unknown, the spread fits no class with room to spare")


NOT_A_VARIANT = frozenset({"<*>", "<NON_REF>", "<CGA_NOCALL>", "<CGA_CNVWIN>", "*"})
"""ALT values that name no variant: the gVCF reference allele, Complete Genomics' no-call and
copy-number window, the spanning deletion `*` (the variant is another record's)."""


def variant_class_expr(ref: pl.Expr, alt: pl.Expr) -> pl.Expr:
    """The `VariantClass` of one ALT against its REF, as a string, or null for no variant.

    A sequence ALT is VCF's: bases `ACGTN` only, or the empty allele (`schema.EMPTY_ALLELE`: a deletion
    written without its anchor base, beside a REF of the deleted bases; an insertion point's REF). An
    ALT of other letters is no VCF allele (Veritas's `CNV`, F15) and has no class: counting it as an
    indel would claim content the file does not hold.
    """
    sequence = alt.str.contains(r"^[ACGTNacgtn]*$") & ref.str.contains(r"^[A-Za-z]*$") & ~((alt == "") & (ref == ""))
    symbolic = alt.str.starts_with("<") | alt.str.contains(r"[\[\]]") | alt.str.contains(r"^\.[A-Za-z]+$|^[A-Za-z]+\.$")
    return (
        pl.when(alt.is_in(list(NOT_A_VARIANT)) | alt.is_null())
        .then(None)
        .when(symbolic)
        .then(pl.lit(VariantClass.SV.value))
        .when(sequence & (alt.str.len_bytes() != ref.str.len_bytes()))
        .then(pl.lit(VariantClass.INDEL.value))
        .when(sequence & (alt.str.len_bytes() == 1))
        .then(pl.lit(VariantClass.SNV.value))
        .when(sequence)
        .then(pl.lit(VariantClass.MNV.value))
        .otherwise(None)
    )


def variant_class_counts(calls: pl.LazyFrame) -> dict[VariantClass, int]:
    """ALT alleles per variant class over a `calls` frame (`source_ref`, `alts`), every class listed."""
    counted = (
        calls.select("source_ref", "alts")
        .explode("alts", empty_as_null=False)
        .select(variant_class_expr(pl.col("source_ref"), pl.col("alts")).alias("variant_class"))
        .drop_nulls()
        .group_by("variant_class")
        .len()
        .collect()
    )
    found = dict(zip(counted["variant_class"].to_list(), counted["len"].to_list(), strict=True))
    return {cls: int(found.get(cls.value, 0)) for cls in VariantClass}


MIN_FILTER_ALLELES = 10_000
"""ALT alleles of one class the rows need before they can show a content filter: an exome's worth.
Fewer, and a file with no indels may just be small."""

MAX_FILTER_RESIDUE = 0.01
"""Share of ALT alleles outside the kept class at or under which the rows show a content filter.
`bcftools view -v snps` keeps a multi-allelic site whose other ALT is an indel, so a filtered file is
not always pure: measured 2026-10-22 on PGP-HMS split DRAGEN genomes, 0.16% in a `-v snps` half and
0.48% in a `-v indels` half; eleven Isaac, CASAVA, GATK SelectVariants and `bcftools filter -i
'TYPE=…'` halves 0%."""

MIN_UNFILTERED_RESIDUE = 0.05
"""Share at or above which the rows show the file was not cut to one class: an unsplit genome's
indels are 21–26% of its ALT alleles (measured on a DRAGEN and a DeepVariant genome). Between the two
the rows do not say."""

_TYPE_WORDS: dict[str, ContentFilter] = {
    "snp": ContentFilter.SNV_ONLY,
    "snps": ContentFilter.SNV_ONLY,
    "indel": ContentFilter.INDEL_ONLY,
    "indels": ContentFilter.INDEL_ONLY,
    "symbolic": ContentFilter.SV_ONLY,
}
"""bcftools (`snps`, `TYPE="snp"`) and GATK (`SNP`, `INDEL`, `SYMBOLIC`) type names of one class."""

_BCFTOOLS_TYPES = frozenset({"snps", "mnps", "indels", "other"})
"""bcftools `-v`/`-V` types that hold variants (`ref` and `bnd` aside)."""

_NAME_WORDS: dict[str, ContentFilter] = {
    "snv": ContentFilter.SNV_ONLY,
    "snvs": ContentFilter.SNV_ONLY,
    "snp": ContentFilter.SNV_ONLY,
    "snps": ContentFilter.SNV_ONLY,
    "indel": ContentFilter.INDEL_ONLY,
    "indels": ContentFilter.INDEL_ONLY,
    "sv": ContentFilter.SV_ONLY,
    "svs": ContentFilter.SV_ONLY,
    "cnv": ContentFilter.SV_ONLY,
}
"""File-name words that name one variant class (DRAGEN's `.snv.`/`.sv.`/`.cnv.`, Illumina's
`.snps.`/`.indels.`). A copy-number file holds symbolic alleles only, so `cnv` names `sv_only`."""

_NAME_WORD = re.compile(r"(?<![a-z0-9])(" + "|".join(sorted(_NAME_WORDS, key=len, reverse=True)) + r")(?=[._-])")


def _quoted(text: str) -> str:
    return text[1:-1] if len(text) > 1 and text[0] == text[-1] and text[0] in "'\"" else text


def _bcftools_filter(line: str) -> ContentFilter | None:
    """The one class a bcftools `view`/`filter` command line keeps, or None."""
    command = line.split("=", 1)[1].split("; Date=", 1)[0]
    if re.search(r"(?:^|\s)(?:-s|--soft-filter)[\s=]", command):
        return None  # a soft filter flags rows in FILTER and removes none
    kept: set[ContentFilter] = set()
    for words in re.findall(r"(?:^|\s)(?:-v|--types)[\s=](\S+)", command):
        types = set(words.split(","))
        if len(types) == 1 and words in _TYPE_WORDS:
            kept.add(_TYPE_WORDS[words])
        else:
            return None
    for words in re.findall(r"(?:^|\s)(?:-V|--exclude-types)[\s=](\S+)", command):
        left = _BCFTOOLS_TYPES - set(words.split(","))
        if len(left) == 1 and next(iter(left)) in _TYPE_WORDS:
            kept.add(_TYPE_WORDS[next(iter(left))])
        else:
            return None
    for flag, expression in re.findall(r"(?:^|\s)(-i|--include|-e|--exclude)[\s=]('[^']*'|\"[^\"]*\"|\S+)", command):
        include = flag in ("-i", "--include")
        match = re.fullmatch(r'\s*TYPE\s*(==?|!=)\s*"(\w+)"\s*', _quoted(expression))
        if not match or (match.group(1) == "!=") is include or match.group(2) not in _TYPE_WORDS:
            return None
        kept.add(_TYPE_WORDS[match.group(2)])
    return kept.pop() if len(kept) == 1 else None


def _gatk_select(line: str) -> ContentFilter | None:
    """The one class a GATK SelectVariants command line keeps (`selectTypeToInclude=[SNP]` in GATK 3,
    `--select-type-to-include SNP` in GATK 4), or None."""
    included: list[str] = []
    for listed in re.findall(r"selectTypeToInclude=\[([^\]]*)\]", line):
        included += [word.strip() for word in listed.split(",") if word.strip()]
    included += re.findall(r"(?:--select-type-to-include|-select-type|-selectType)[\s=](\w+)", line)
    words = {word.lower() for word in included}
    if len(words) == 1 and next(iter(words)) in _TYPE_WORDS:
        return _TYPE_WORDS[next(iter(words))]
    return None


def stated_content_filters(meta: list[str]) -> list[tuple[ContentFilter, str]]:
    """Each header line that records cutting the rows down to one variant class, with its meta key.

    Read: bcftools `view`/`filter` with `-v`/`--types`, `-V`/`--exclude-types`, or an
    `-i`/`-e` expression that is a bare `TYPE` comparison; GATK SelectVariants with one select type.
    `##source=DRAGEN_SNV` is not one: it names DRAGEN's small-variant caller, and its unsplit
    `hard-filtered` output carries indels under it.
    """
    found: list[tuple[ContentFilter, str]] = []
    for line in meta:
        key = line[2:].split("=", 1)[0]
        if key.startswith(("bcftools_viewCommand", "bcftools_filterCommand")):
            kept = _bcftools_filter(line)
        elif key.startswith("GATKCommandLine") and "SelectVariants" in line:
            kept = _gatk_select(line)
        else:
            continue
        if kept:
            found.append((kept, f"##{key}"))
    return found


def named_content_filter(name: str) -> tuple[ContentFilter, str] | None:
    """The class a file name's words name (`x.snp.vcf.gz`, `x.indels.vcf.gz`), with the word, or None.

    Only one class may be named; a name with two (`snps_indels`) names none."""
    words = _NAME_WORD.findall(Path(name).name.lower())
    named = {_NAME_WORDS[word] for word in words}
    if len(named) != 1:
        return None
    return named.pop(), f".{words[0]}."


class ContentFilterInference(BaseModel):
    """The content filter the evidence supports and its basis, with a warning per conflict."""

    model_config = ConfigDict(frozen=True)

    content_filter: ContentFilter | None
    basis: ContentFilterBasis | None
    reason: str
    warnings: tuple[str, ...] = ()


def rows_content_filter(counts: dict[VariantClass, int]) -> tuple[ContentFilter | None, bool, str]:
    """What the allele counts show: (the filter, whether the rows settle it, why).

    `(filter, True, …)`: one class with at most `MAX_FILTER_RESIDUE` of the rest. `(None, True, …)`:
    the other classes are at least `MIN_UNFILTERED_RESIDUE`, so nothing was cut to one class.
    `(None, False, …)`: too few alleles, or a residue between the two, or the main class is MNV.
    """
    total = sum(counts.values())
    main = max(VariantClass, key=lambda cls: (counts.get(cls, 0), -list(VariantClass).index(cls)))
    kept = counts.get(main, 0)
    residue = (total - kept) / total if total else 0.0
    where = f"{kept:,} {main.value} of {total:,} ALT allele(s), {residue:.3%} other"
    if kept < MIN_FILTER_ALLELES:
        return None, False, f"too few alleles to tell ({where})"
    if residue >= MIN_UNFILTERED_RESIDUE:
        return None, True, f"more than one class ({where})"
    filters = {cls: f for f, cls in CONTENT_FILTER_KEEPS.items()}
    if residue <= MAX_FILTER_RESIDUE and main in filters:
        return filters[main], True, f"one class ({where})"
    return None, False, f"a residue the rows do not settle ({where})"


def infer_content_filter(
    counts: dict[VariantClass, int],
    stated: list[tuple[ContentFilter, str]],
    named: tuple[ContentFilter, str] | None,
) -> ContentFilterInference:
    """The facet's content filter: the rows first, then the header, then the file name.

    The strongest evidence that settles it wins, and every weaker piece that says otherwise is a
    warning. A header that records a filter its rows contradict (a `bcftools concat` of the `.snp.` and
    `.indel.` halves keeps the first half's `view -v snps` line) is reported, not believed.
    """
    from_rows, settled, rows_why = rows_content_filter(counts)
    stated_filters = sorted({kept for kept, _ in stated}, key=list(ContentFilter).index)
    stated_why = ", ".join(f"{key} keeps {kept.value}" for kept, key in stated)
    warnings: list[str] = []
    if len(stated_filters) > 1:
        warnings.append(f"content filter: header lines keep different classes ({stated_why}); the header settles none")
    header = stated_filters[0] if len(stated_filters) == 1 else None

    def say(what: ContentFilter | None) -> str:
        return what.value if what else "no filter"

    if settled:
        chosen, basis, reason = from_rows, ContentFilterBasis.ROWS if from_rows else None, f"the rows show {rows_why}"
        if header is not None and header is not from_rows:
            warnings.append(f"content filter: the header says {header.value} ({stated_why}), the rows show {say(from_rows)}: {rows_why}; the rows are kept")
        if named and named[0] is not from_rows:
            warnings.append(f"content filter: the file name ({named[1]}) says {named[0].value}, the rows show {say(from_rows)}: {rows_why}; the rows are kept")
    elif header is not None:
        chosen, basis, reason = header, ContentFilterBasis.STATED, f"the header says so ({stated_why}); the rows do not settle it: {rows_why}"
        if named and named[0] is not header:
            warnings.append(f"content filter: the file name ({named[1]}) says {named[0].value}, the header says {header.value} ({stated_why}); the header is kept")
    elif named and not stated_filters:
        chosen, basis, reason = named[0], ContentFilterBasis.FILE_NAME, f"only the file name says so ({named[1]}); the rows do not settle it: {rows_why}"
    else:
        chosen, basis, reason = None, None, f"no evidence settles one; the rows: {rows_why}"
    return ContentFilterInference(content_filter=chosen, basis=basis, reason=reason, warnings=tuple(warnings))


REFERENCE_BY_OMISSION_CALLERS: frozenset[str] = frozenset({"CGAPipeline"})
"""Callers whose VCF leaves out exactly the positions they called reference, so an omitted position
is a reference call (`AbsenceMeans.REFERENCE`). Complete Genomics' CGAPipeline: on a PGP-HMS genome
its vcfBeta covers every base its native `var` file does not call reference (0 of 317.7M missed,
2026-10-22, `test_real_var_converges_with_its_vcf_beta`)."""


def absence_means(callset_type: CallsetType | None, has_blocks: bool, *, omitted_called_reference: bool = False) -> AbsenceMeans:
    """What a position absent from the facet means (SCHEMA.md Decision 4).

    A producer that leaves out exactly what it called reference (`omitted_called_reference`, a caller
    in `REFERENCE_BY_OMISSION_CALLERS`) makes absence a reference call, but only where the rows show a
    whole genome: a region cut from it, which RM4 cannot always name, would otherwise claim reference
    for the regions cut away. Otherwise reference blocks decide it wherever they exist. Without them
    only a whole-genome callset makes absence weak evidence of reference; an exome, a panel, an array
    or imputed content says nothing outside what it looked at, and neither does a callset of unknown
    type.
    """
    if omitted_called_reference and callset_type is CallsetType.WGS:
        return AbsenceMeans.REFERENCE
    if has_blocks:
        return AbsenceMeans.BLOCKS
    return AbsenceMeans.WEAK_REFERENCE if callset_type is CallsetType.WGS else AbsenceMeans.UNKNOWN


class AbsenceReason(str, Enum):
    """Why a position absent from a facet reads as reference, or why it says nothing (`absence_support`).

    just-vcf's own vocabulary (charter 2). The first three read as reference; every other one is
    unknown, never a negative."""

    WHOLE_GENOME_SILENCE = "whole_genome_silence"
    """`absence_means = weak_reference`: a whole-genome caller reported no variant there."""
    OMITTED_CALLED_REFERENCE = "omitted_called_reference"
    """`absence_means = reference`: the producer called every position it left out reference (CGI's VCF)."""
    CONTENT_FILTER_FORCED = "content_filter_forced"
    """Either of the two above on a whole genome cut to one class, read as reference for the kept class
    only because the caller passed `force=True` (the person's 2026-10-22 ruling: restoration on a
    filtered genome is additionally gated). `absence_means` still says which of the two it was."""
    OUTSIDE_BLOCKS = "outside_blocks"
    """`absence_means = blocks`: the facet states its reference in blocks, and none covers the position."""
    ARRAY = "array"
    """An array: it probes chosen sites, and an absent one is unknown (RM281's `chip_absence` says why)."""
    TARGETED = "targeted"
    """An exome or a panel: nothing outside its targets was looked at."""
    IMPUTED = "imputed"
    """Imputed content: its sites are a reference panel's, so its silence says nothing about others."""
    CALLSET_TYPE_UNKNOWN = "callset_type_unknown"
    """`absence_means = unknown` and RM4 could not name the callset type."""
    CLASS_NOT_COUNTED = "class_not_counted"
    """The facet's `variant_classes` were not counted (a genotype table writes no REF/ALT to count)."""
    CLASS_NOT_CALLED = "class_not_called"
    """The facet wrote no allele of the asked class, so its silence says nothing about that class."""
    CONTENT_FILTER_HELD = "content_filter_held"
    """A whole genome cut to the asked class: held back unless `force=True`."""
    CLASS_CUT_BY_FILTER = "class_cut_by_filter"
    """A whole genome cut to another class: the asked class was cut away, and `force` does not lift it."""


READS_REFERENCE = frozenset({AbsenceReason.WHOLE_GENOME_SILENCE, AbsenceReason.OMITTED_CALLED_REFERENCE, AbsenceReason.CONTENT_FILTER_FORCED})
"""The reasons under which an absent position reads as reference."""

_UNKNOWN_BY_CALLSET = {
    CallsetType.ARRAY.value: AbsenceReason.ARRAY,
    CallsetType.WES.value: AbsenceReason.TARGETED,
    CallsetType.PANEL.value: AbsenceReason.TARGETED,
    CallsetType.IMPUTED.value: AbsenceReason.IMPUTED,
}


class AbsenceSupport(BaseModel):
    """What a position absent from a facet says about reference, and why. Truthy only when it reads as
    reference; `reads_reference` is `True` or `None` (unknown), never `False`."""

    model_config = ConfigDict(frozen=True)

    reads_reference: bool | None
    reason: AbsenceReason

    def __bool__(self) -> bool:
        return self.reads_reference is True


def absence_support(facet: Mapping[str, object], variant_class: VariantClass | str, *, force: bool = False) -> AbsenceSupport:
    """`absence_is_reference` with its reason: the gates below, in that order, each with its own reason.

    A facet whose `absence_means` is `unknown` (or missing) is answered by its `callset_type`: `array`,
    `targeted` (an exome or a panel), `imputed`, else `callset_type_unknown`.
    """
    asked = VariantClass(variant_class)
    means = facet.get("absence_means")
    if means == AbsenceMeans.BLOCKS.value:
        return AbsenceSupport(reads_reference=None, reason=AbsenceReason.OUTSIDE_BLOCKS)
    if means not in (AbsenceMeans.WEAK_REFERENCE.value, AbsenceMeans.REFERENCE.value):
        unknown = _UNKNOWN_BY_CALLSET.get(str(facet.get("callset_type")), AbsenceReason.CALLSET_TYPE_UNKNOWN)
        return AbsenceSupport(reads_reference=None, reason=unknown)
    counted = facet.get("variant_classes")
    if not isinstance(counted, list):
        return AbsenceSupport(reads_reference=None, reason=AbsenceReason.CLASS_NOT_COUNTED)
    if asked.value not in {entry["variant_class"] for entry in counted if entry["alleles"]}:
        return AbsenceSupport(reads_reference=None, reason=AbsenceReason.CLASS_NOT_CALLED)
    content_filter = facet.get("content_filter")
    if content_filter is not None:
        if CONTENT_FILTER_KEEPS[ContentFilter(content_filter)] is not asked:
            return AbsenceSupport(reads_reference=None, reason=AbsenceReason.CLASS_CUT_BY_FILTER)
        if not force:
            return AbsenceSupport(reads_reference=None, reason=AbsenceReason.CONTENT_FILTER_HELD)
        return AbsenceSupport(reads_reference=True, reason=AbsenceReason.CONTENT_FILTER_FORCED)
    if means == AbsenceMeans.REFERENCE.value:
        return AbsenceSupport(reads_reference=True, reason=AbsenceReason.OMITTED_CALLED_REFERENCE)
    return AbsenceSupport(reads_reference=True, reason=AbsenceReason.WHOLE_GENOME_SILENCE)


def absence_is_reference(facet: Mapping[str, object], variant_class: VariantClass | str, *, force: bool = False) -> bool | None:
    """Whether a position absent from `facet` reads as evidence of reference for `variant_class`.

    `facet` is one `facets` row (`load(...).facets.collect().row(0, named=True)`). A position is
    absent when no `calls` row and no `blocks` row covers it. The answer is Kleene with no negative:
    `True` or `None` (unknown), never `False`, because a facet's silence never shows that a variant is
    there. `absence_support` gives the same answer with its reason. `True` needs all of:

    - `absence_means` is `weak_reference` (a whole-genome callset without blocks) or `reference` (a
      producer that called every omitted position reference; how strong the evidence is, the member
      says). A gVCF's absence is decided per position by its `blocks` frame, which this facet-level
      answer cannot see: `None`.
    - `variant_class` is among the facet's `variant_classes`. A callset that wrote no SVs says
      nothing about SVs, and a facet whose classes were not counted says nothing at all.
    - the facet has no `content_filter`. A whole genome cut to one class (the `.snp.` half of a split
      delivery) still reads as `wgs`, but its absences say nothing about the classes cut away, and the
      cut is a reason to hold back for the kept class too. `force=True` lifts that hold for the kept
      class only; a class the filter cut stays `None`.

    Hom-ref restoration from absence, in just-vcf or in a consumer, gates on this and refuses a
    position where the answer is not `True` (docs/SCHEMA.md Decision 4).
    """
    return absence_support(facet, variant_class, force=force).reads_reference
