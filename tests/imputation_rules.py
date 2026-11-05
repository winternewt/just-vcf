"""What an imputation server checks in an uploaded VCF, read off a file `write_vcf` wrote (RM287).

A test helper, not library API: it reads the VCF as text, the way the servers' own Java reader does
(`VCFLineParser` keeps ALT and FILTER verbatim), never through `load` or htslib, and reports each rule
as `met` True, False, or None with counts. None is "not checked": the rule needs something this file
does not hold (the reference panel's legend, the job's build), and is never reported as passed.

The rules, their sources and the outcome just-vcf's output gets are in docs/SCHEMA.md, "What an
imputation server checks". Sources, read 2026-11-05:

- MIS2 docs: https://genepi.github.io/michigan-imputationserver/prepare-your-data,
  .../getting-started, .../pipeline, .../qc-allele-swaps
- TOPMed (TIS v2) docs: https://statgen.github.io/tis-v2-docs/imputation/data-preparation,
  .../pipeline-overview
- the servers' QC code, genepi/imputationserver-utils at 621293a1 (`VcfFileUtil`,
  `InputValidationCommand`, `QualityControlCommand`, `fastqc/StatisticsTask`, `util/GenomicTools`,
  `vcf/VCFLineParser`), and genepi/imputationserver2 at 5506abe3 (`nextflow.config`)
- checkVCF, which the MIS docs recommend before upload: zhanxw/checkVCF at f26d497c (`checkVCF.py`)
"""

import gzip
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

BGZF_MAGIC = b"\x1f\x8b\x08\x04"
"""A gzip member with FEXTRA set; BGZF puts its `BC` subfield there (SAM spec §4.1)."""

VALID_CHROMOSOMES = frozenset([*(str(n) for n in range(1, 24)), "X", "MT"])
"""`VcfFileUtil.validChromosomes`, bare; the set also holds each name with a `chr` prefix. No Y, no
`M` (so `chrM` is refused, `chrMT` accepted), no `XY`."""

NON_PAR = {"GRCh37": (2_699_520, 154_931_044), "GRCh38": (2_781_479, 155_701_383)}
"""`StatisticsTask.prepareChrX`: a chrX record with `nonParStart <= POS <= nonParEnd` is nonPAR
(hg19 the default, hg38 its own bounds)."""

MIN_SAMPLES = 20
"""imputationserver2 `nextflow.config`, `min_samples = 20`; TOPMed's docs, "a minimum of 20 samples"."""

MAX_MIXED_CHRX = 0.10
"""The pipeline docs: "the proportion of mixed genotypes (e.g., 1/.) is < 10%" on chrX nonPAR."""

MIN_SITE_CALL_RATE = 0.90
"""`StatisticsTask`: a site whose no-call share is above 0.10 is excluded ("SNP call rate < 90%")."""

ACGT = frozenset("ACGT")
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

SITE_EXCLUSIONS = ("multiallelic", "invalid_alleles", "filtered", "indel", "monomorphic")
"""`StatisticsTask.processLine`'s exclusions before the panel lookup, in its order: a record is
excluded by the first that applies, so each is counted once."""

PANEL_RULES = {
    "panel_overlap": "at least 50% of a chunk's variants in the reference panel",
    "chunk_min_variants": "at least 3 variants in the panel per chunk (HRC: 20)",
    "allele_switches": "allele switches against the panel, the job failing past the panel's limit (100 in the docs)",
    "strand_flips": "strand flips against the panel, the job failing past the panel's limit (100 in the docs)",
    "allele_mismatch": "alleles that match the panel's neither as written, switched nor flipped",
    "chunk_sample_call_rate": "each sample called at 50% of a chunk's variants that passed QC",
    "allele_frequency": "the study's allele frequencies against the panel's (chi-square)",
}
"""Rules that read the reference panel's legend, which no local file holds: always not checked."""

RefBases = Callable[[str, list[tuple[int, int]]], list[str | None]]
"""The reference's bases for a written contig name and 1-based `(pos, length)` spans, None where it
holds none."""


@dataclass(frozen=True)
class Outcome:
    met: bool | None
    """True: the file passes the rule. False: it does not. None: not checked (see `why`)."""
    counts: dict[str, int] = field(default_factory=dict)
    why: str = ""


@dataclass(frozen=True)
class Record:
    chrom: str
    pos: int
    ref: str
    alt: str
    filters: str
    format_keys: list[str]
    samples: list[str]

    def gt(self, sample: int) -> str | None:
        """The sample's GT text, None when FORMAT has no GT."""
        if "GT" not in self.format_keys:
            return None
        values = self.samples[sample].split(":")
        index = self.format_keys.index("GT")
        return values[index] if index < len(values) else "."


@dataclass(frozen=True)
class VcfText:
    bgzf: bool
    gzip: bool
    meta: list[str]
    samples: list[str]
    records: list[Record]

    @property
    def fileformat(self) -> str | None:
        found = [line.split("=", 1)[1] for line in self.meta if line.startswith("##fileformat=")]
        return found[0] if found else None


def read_vcf_text(path: Path) -> VcfText:
    """The file as the server's reader sees it: compression by magic bytes, then tab-split lines."""
    data = path.read_bytes()
    is_gzip = data[:2] == b"\x1f\x8b"
    text = (gzip.decompress(data) if is_gzip else data).decode()
    meta: list[str] = []
    samples: list[str] = []
    records: list[Record] = []
    for line in text.splitlines():
        if line.startswith("##"):
            meta.append(line)
        elif line.startswith("#"):
            samples = line.split("\t")[9:]
        elif line:
            tiles = line.split("\t")
            records.append(Record(tiles[0], int(tiles[1]), tiles[3], tiles[4], tiles[6], tiles[8].split(":") if len(tiles) > 8 else [], tiles[9:]))
    return VcfText(bgzf=data[:4] == BGZF_MAGIC and data[12:14] == b"BC", gzip=is_gzip, meta=meta, samples=samples, records=records)


def _alleles(gt: str | None) -> list[str]:
    return [] if gt is None else gt.replace("|", "/").split("/")


def _is_no_call(gt: str | None) -> bool:
    """`VCFLineParser`'s no-call: no allele called (`./.`, `.`); `0/.` is a called hom-ref there."""
    return gt is None or all(a == "." for a in _alleles(gt))


def _is_mixed(gt: str | None) -> bool:
    """htsjdk's `Genotype.isMixed`: some alleles called and some not (`1/.`)."""
    alleles = _alleles(gt)
    return any(a == "." for a in alleles) and any(a != "." for a in alleles)


def _strip_chr(name: str) -> str:
    """`VcfFileUtil.load`: every `chr` in the name removed before chromosomes are compared."""
    return name.replace("chr", "")


def _valid_base(allele: str) -> bool:
    """`GenomicTools.isValid`: the allele, uppercased, is one of A, C, G and T."""
    return allele.upper() in ACGT


def _site_exclusion(record: Record, n_samples: int) -> str | None:
    if "," in record.alt:
        return "multiallelic"
    if not (_valid_base(record.ref) and _valid_base(record.alt)):
        return "invalid_alleles"
    if record.filters not in ("PASS", "."):
        return "filtered"
    if len(record.ref) != len(record.alt):
        return "indel"
    if n_samples > 1 and all(set(_alleles(record.gt(i))) == {"0"} for i in range(n_samples)):
        return "monomorphic"
    return None


def check_vcf(path: Path, build: str | None = None, reference: RefBases | None = None) -> dict[str, Outcome]:
    """Every rule's outcome for the VCF at `path`.

    `build` is the job's build, `GRCh37` or `GRCh38` (the servers' hg19 and hg38); without it the
    build rules are not checked. `reference` gives the reference's bases, as checkVCF reads its FASTA;
    without it REF is not checked.
    """
    vcf = read_vcf_text(path)
    records, n = vcf.records, len(vcf.samples)
    out: dict[str, Outcome] = {}

    out["bgzf"] = Outcome(vcf.bgzf, {"bgzf": int(vcf.bgzf), "gzip": int(vcf.gzip)})
    version = vcf.fileformat
    if version is None or not version.startswith("VCFv"):
        out["vcf_version_at_most_4_2"] = Outcome(None, why="no ##fileformat line")
    else:
        major, minor = (int(x) for x in version.removeprefix("VCFv").split(".")[:2])
        out["vcf_version_at_most_4_2"] = Outcome((major, minor) <= (4, 2), {f"v{major}_{minor}": 1})

    chromosomes = Counter(_strip_chr(r.chrom) for r in records)
    out["one_chromosome"] = Outcome(len(chromosomes) <= 1, {"chromosomes": len(chromosomes)})
    invalid = Counter(r.chrom for r in records if _strip_chr(r.chrom) not in VALID_CHROMOSOMES)
    out["valid_chromosome"] = Outcome(not invalid, {f"records_on_{name}": count for name, count in sorted(invalid.items())})

    prefixed = sum(r.chrom.startswith("chr") for r in records)
    counts = {"records_with_chr": prefixed, "records_without_chr": len(records) - prefixed}
    if build not in NON_PAR:
        out["chr_prefix_for_build"] = Outcome(None, counts, "the job's build was not given")
    elif build == "GRCh37":
        out["chr_prefix_for_build"] = Outcome(prefixed == 0, counts)
    else:
        out["chr_prefix_for_build"] = Outcome(prefixed == len(records), counts)

    out_of_order = 0
    seen: list[str] = []
    for before, record in zip(records, records[1:]):
        if record.chrom != before.chrom:
            if record.chrom in seen:
                out_of_order += 1
            seen.append(before.chrom)
        elif record.pos < before.pos:
            out_of_order += 1
    out["sorted"] = Outcome(out_of_order == 0, {"records_out_of_order": out_of_order})

    duplicates = sum(1 for before, record in zip(records, records[1:]) if (record.chrom, record.pos) == (before.chrom, before.pos))
    dup_filter = sum("DUP" in r.filters for r in records)
    out["no_duplicate_positions"] = Outcome(duplicates == 0 and dup_filter == 0, {"records_at_a_previous_position": duplicates, "dup_filter": dup_filter})

    missing_gt = sum("GT" not in r.format_keys for r in records)
    out["gt_in_format"] = Outcome(missing_gt == 0, {"records_without_gt": missing_gt})
    same = sum(r.ref == r.alt for r in records)
    out["ref_differs_from_alt"] = Outcome(same == 0, {"records_ref_equals_alt": same})

    out["samples_at_least_20"] = Outcome(n >= MIN_SAMPLES, {"samples": n})
    out["unique_sample_ids"] = Outcome(len(set(vcf.samples)) == n, {"samples": n, "distinct": len(set(vcf.samples))})

    excluded = Counter(reason for r in records if (reason := _site_exclusion(r, n)) is not None)
    for reason in SITE_EXCLUSIONS:
        out[reason] = Outcome(excluded[reason] == 0, {"records_excluded": excluded[reason]})
    alt_unstated = sum(r.alt == "." for r in records)
    out["invalid_alleles"].counts["alt_unstated"] = alt_unstated

    low_call = sum(1 for r in records if n and sum(_is_no_call(r.gt(i)) for i in range(n)) / n > 1 - MIN_SITE_CALL_RATE)
    out["site_call_rate"] = Outcome(low_call == 0, {"records_excluded": low_call})

    if reference is None:
        out["ref_matches_reference"] = Outcome(None, why="no reference given")
    else:
        mismatch = unread = 0
        by_chrom: dict[str, list[Record]] = {}
        for r in records:
            by_chrom.setdefault(r.chrom, []).append(r)
        for chrom, rows in by_chrom.items():
            for r, b in zip(rows, reference(chrom, [(r.pos, len(r.ref)) for r in rows])):
                if b is None:
                    unread += 1
                elif b.upper() != r.ref.upper():
                    mismatch += 1
        out["ref_matches_reference"] = Outcome(None if unread else mismatch == 0, {"records_mismatched": mismatch, "records_unread": unread})

    out.update(_chrx_rules(records, n, build))
    out["palindromic_sites"] = Outcome(None, palindromes(records), "the server counts A/T and C/G sites against the panel and keeps them (`complicatedGenotypes`); which strand each is on needs the panel")
    for rule, text in PANEL_RULES.items():
        out[rule] = Outcome(None, why=f"needs the reference panel: {text}")
    return out


def _chrx_rules(records: list[Record], n: int, build: str | None) -> dict[str, Outcome]:
    if build not in NON_PAR:
        why = "the job's build was not given"
        return {"chrx_ploidy": Outcome(None, why=why), "chrx_mixed_genotypes": Outcome(None, why=why)}
    start, end = NON_PAR[build]
    nonpar = [r for r in records if _strip_chr(r.chrom) == "X" and start <= r.pos <= end]
    if not nonpar or not n:
        why = "no chrX nonPAR record"
        return {"chrx_ploidy": Outcome(None, why=why), "chrx_mixed_genotypes": Outcome(None, why=why)}
    # `StatisticsTask.checkPloidy`: a sample seen haploid, then with a ploidy other than one, fails the job.
    haploid_seen: set[int] = set()
    after_haploid = haploid = 0
    mixed = Counter[int]()
    for r in nonpar:
        for i in range(n):
            gt = r.gt(i)
            ploidy = len(_alleles(gt))
            if i in haploid_seen and ploidy != 1:
                after_haploid += 1
            if ploidy == 1:
                haploid_seen.add(i)
                haploid += 1
            mixed[i] += _is_mixed(gt)
    worst = max(mixed[i] / len(nonpar) for i in range(n))
    return {
        "chrx_ploidy": Outcome(after_haploid == 0, {"nonpar_records": len(nonpar), "haploid_calls": haploid, "non_haploid_after_haploid": after_haploid}),
        "chrx_mixed_genotypes": Outcome(worst <= MAX_MIXED_CHRX, {"nonpar_records": len(nonpar), "mixed_calls": sum(mixed.values())}),
    }


def palindromes(records: Iterable[Record]) -> dict[str, int]:
    """A/T and C/G sites as written: `het` reads the same on either strand, `hom_alt` does not (`TT` at
    a reference `A` is a minus-strand hom-ref read forward); `alt_unstated` is a record whose ALT is `.`,
    which cannot be told apart from the record."""
    counts = Counter({"het": 0, "hom_alt": 0, "alt_unstated": 0})
    for r in records:
        if r.alt == ".":
            counts["alt_unstated"] += 1
        elif _valid_base(r.ref) and _valid_base(r.alt) and COMPLEMENT[r.ref.upper()] == r.alt.upper():
            alleles = {a for a in _alleles(r.gt(0)) if a != "."}
            counts["het" if alleles == {"0", "1"} else "hom_alt" if alleles == {"1"} else "other"] += 1
    return dict(counts)


def strand_withheld(calls: pl.DataFrame) -> dict[str, int]:
    """Array rows `resolve_ref` read REF for and withheld the alleles of (F13: a `mixed` facet), by what
    the export called: the reference base only (a forward hom-ref, or a minus hom-alt of an A/T or C/G
    probe), its complement only (a forward hom-alt of such a probe, or a minus hom-ref), or anything
    else. These never reach the VCF: `write_vcf` counts them as `alts_unresolved`."""
    withheld = calls.filter(pl.col("ref").is_not_null() & pl.col("alts").is_null() & pl.col("source_ref").is_null())
    counts = Counter({"ref_base_only": 0, "complement_only": 0, "other": 0})
    for ref, alleles in withheld.select("ref", "gt_alleles").iter_rows():
        called = {a for a in alleles if a is not None}
        counts["ref_base_only" if called == {ref} else "complement_only" if called == {COMPLEMENT.get(ref)} else "other"] += 1
    return dict(counts)
