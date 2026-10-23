"""REF read from the reference sequence (RM3, charter invariant 10).

`resolve_ref(loaded, reference)` reads each call's REF from the facet's assembly sequence, held in a
gtars `RefgetStore`, and returns the same facet with the slots SCHEMA.md leaves open until then filled:

- a VCF row (it has `source_ref`): `ref` is the reference's bases over the same span, and `ref_check`
  says whether the source agreed (`match`, `mismatch`). `source_ref` is never touched;
- an array row (no `source_ref`, single-base alleles): `ref` is the reference base, `alts` are the
  observed alleles other than REF, sorted, and `gt` indexes `[ref] + alts`. REF is never taken from
  the genotype; the genotype only decides which of its alleles are not REF;
- a Living DNA row whose genotype is two allele sequences written back to back (`CCA`, F4): split
  where exactly one split fits the reference (`_split_concatenated`), and given REF and ALT only where
  the two alleles say which of them is the reference.

Before any base is read the reference has to be shown to be the facet's assembly, contig by contig
(`ContigCheck`). A facet that stated contig lengths (VCF `##contig`) is checked by length and md5;
one that did not (arrays) by its assembly label against the label the reference's own lengths give.
An array that states no label either has one inferred from its rows against each collection first
(F9, `just_vcf.assembly.infer_assembly_from_rows`), or is refused with the scores.
A primary contig of another length means another assembly, and the whole facet is refused
(`ReferenceMismatchError`); any other contig that cannot be shown to match is left unresolved, with
its rows counted in the warnings. Nothing falls back to an unchecked read.

Coordinates: `pos` is 1-based; the store takes 0-based half-open ranges. `_span` is the one place
that converts between them (charter 13).
"""

import hashlib
import logging
import os
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import polars as pl
from gtars.refget import RefgetStore

from just_vcf.assembly import RowsInference, infer_assembly_from_rows, sample_snv_rows
from just_vcf import cache
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.schema import CALLS, CONTIG_RANK, FACETS, SCHEMA_VERSION, AssemblyLabelBasis, RefCheck, StrandCheck, canonical_chrom, canonical_order, empty
from just_vcf.strand import StrandResult, check_strand, withheld_reason

logger = logging.getLogger(__name__)

RESOLVE_VERSION = 4
"""Part of the spill key: bumped whenever the rule changes what `resolve_ref` writes."""

CONCATENATED_ALLELES: frozenset[str] = frozenset({"living_dna"})
"""Formats whose multi-letter genotype text is the two allele sequences written back to back, as the
vendor states it ("both allele sequences concatenated"). Only these are split against the reference:
another format's long genotype text may mean something else."""

_CONCATENATED = r"^[ACGT]{3,}$"
"""A genotype text that can be two allele sequences: bases only, at least one of them two long."""

ASSEMBLY_NAMESPACE = "just_vcf_assembly"
"""Collection-alias namespace for an assembly label a caller asserted when adding a FASTA."""

_FASTA_NAMESPACE = "just_vcf_fasta"
"""Collection-alias namespace mapping a FASTA file's fingerprint to its collection, so a FASTA already
added is not digested again."""

CHR1_LENGTHS: dict[int, str] = {247_249_719: "GRCh36", 249_250_621: "GRCh37", 248_956_422: "GRCh38", 248_387_328: "T2T-CHM13"}
"""Contig 1's length per assembly, the fingerprint that names a reference's assembly. Domain constants."""

PRIMARY: frozenset[str] = frozenset([str(n) for n in range(1, 23)] + ["X", "Y"])
"""Contigs whose length differing from the reference means another assembly, not another contig."""

LABEL_UNVERIFIABLE: dict[str, frozenset[str]] = {"GRCh36": frozenset({"MT"}), "GRCh37": frozenset({"MT"})}
"""Contigs a label cannot vouch for: `GRCh37` covers both hg19, whose chrM is the Yoruba sequence
(16,571 bases), and GRCh37, whose MT is rCRS (16,569). Only a stated length tells them apart.
`GRCh36` likewise: NCBI36 and hg18 carry the Yoruba sequence, and an array labelled build 36 does not
say whether its MT positions are on it or on rCRS."""

_REF_SPAN_FOR = {"XY": "X"}
"""Canonical names read from another reference contig: array code 25 (`XY`, the X pseudo-autosomal
region) is positioned on X."""

_IUPAC: dict[str, frozenset[str]] = {
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}  # fmt: skip
"""Ambiguity codes and the bases each stands for. VCF 4.5 writes such a reference base as a concrete
one, so a mismatch against one of these is expected and is counted apart."""

_BASES = ["A", "C", "G", "T"]


class ReferenceMismatchError(ValueError):
    """The reference cannot serve this facet: another assembly, or nothing could be shown to match."""


class ContigCheck(str, Enum):
    """How a facet's contig was matched against a reference contig. The first two allow reading."""

    LENGTH = "length"
    """The facet stated the contig's length (and md5, if any), and the reference agrees."""
    LABEL = "label"
    """The facet stated no length; its assembly label is the one the reference's lengths give."""
    NOT_IN_REFERENCE = "not_in_reference"
    AMBIGUOUS = "ambiguous"
    """Two reference contigs have the same canonical name (`1` and `chr1`)."""
    LENGTH_MISMATCH = "length_mismatch"
    MD5_MISMATCH = "md5_mismatch"
    LABEL_UNVERIFIABLE = "label_unverifiable"
    """Verified by label, but the label covers two sequences for this contig (`LABEL_UNVERIFIABLE`)."""
    UNVERIFIED = "unverified"
    """No stated length, and no assembly label on both sides that agrees."""


_READABLE = frozenset({ContigCheck.LENGTH, ContigCheck.LABEL})


@dataclass(frozen=True)
class Contig:
    """One sequence of a reference collection."""

    name: str
    length: int
    md5: str
    digest: str
    """The refget sequence digest (`SQ.` prefix stripped)."""


def reference_dir() -> Path:
    """Where the reference store lives: `JUST_VCF_REFERENCE_DIR`, else `reference/` in the cache dir."""
    configured = os.environ.get("JUST_VCF_REFERENCE_DIR")
    return Path(configured) if configured else cache_dir() / "reference"


class ReferenceStore:
    """A gtars `RefgetStore` on disk holding one or more reference collections (assemblies)."""

    def __init__(self, store: RefgetStore, path: Path | None) -> None:
        store.set_quiet(True)
        self.store = store
        self.path = path

    @classmethod
    def open(cls, path: Path | str | None = None) -> "ReferenceStore":
        """The store at `path` (default `reference_dir()`), created empty when there is none yet."""
        where = Path(path) if path is not None else reference_dir()
        if RefgetStore.store_exists(str(where)):
            return cls(RefgetStore.open_local(str(where)), where)
        where.mkdir(parents=True, exist_ok=True)
        return cls(RefgetStore.on_disk(str(where)), where)

    @classmethod
    def in_memory(cls) -> "ReferenceStore":
        """A store held in memory only, for a FASTA read once."""
        return cls(RefgetStore.in_memory(), None)

    def add_fasta(self, fasta: Path | str, assembly_label: str | None = None) -> str:
        """Add a FASTA (plain or gzipped) as a collection and return its seqcol digest.

        `assembly_label` names the assembly when the sequence's own lengths cannot (a reference
        without contig 1). A label that the lengths contradict is refused. A FASTA added before, same
        path, size and modification time, is not read again.
        """
        fasta = Path(fasta).resolve()
        stat = fasta.stat()
        fingerprint = hashlib.sha256(f"{fasta}\0{stat.st_size}\0{stat.st_mtime_ns}".encode()).hexdigest()[:32]
        known = _FASTA_NAMESPACE in self._namespaces() and fingerprint in self.store.list_collection_aliases(_FASTA_NAMESPACE)
        digest = self.store.get_collection_metadata_by_alias(_FASTA_NAMESPACE, fingerprint).digest if known else None
        if digest is None:
            meta, _ = self.store.add_sequence_collection_from_fasta(str(fasta))
            digest = meta.digest
            self.store.add_collection_alias(_FASTA_NAMESPACE, fingerprint, digest)
        if assembly_label is not None:
            measured = self._label_from_lengths(digest)
            if measured is not None and measured != assembly_label:
                raise ReferenceMismatchError(
                    f"{fasta.name}: its contig 1 is {measured}'s length, not {assembly_label}'s"
                )
            if assembly_label not in self.asserted_labels(digest):
                self.store.add_collection_alias(ASSEMBLY_NAMESPACE, assembly_label, digest)
        return digest

    def collections(self) -> list[str]:
        """The seqcol digests of the collections in the store, sorted."""
        digests: list[str] = []
        page = 0
        while True:
            listing = self.store.list_collections(page=page, page_size=100)
            digests += [meta.digest for meta in listing["results"]]
            if len(digests) >= listing["pagination"]["total"] or not listing["results"]:
                return sorted(digests)
            page += 1

    def contigs(self, digest: str) -> list[Contig]:
        """The collection's sequences, in the collection's order."""
        level2 = self.store.get_collection_level2(digest)
        found = []
        for name, length, sequence in zip(level2["names"], level2["lengths"], level2["sequences"]):
            bare = sequence.removeprefix("SQ.")
            meta = self.store.get_sequence_metadata(bare)
            found.append(Contig(name=name, length=int(length), md5=meta.md5 if meta else "", digest=bare))
        return found

    def contigs_by_chrom(self, digest: str) -> dict[str, list[Contig]]:
        """The collection's sequences under the canonical `chrom` each serves, in collection order.

        A canonical name served by more than one sequence (`1` and `chr1`) lists them all; a caller
        reads none of them. `XY` (array code 25) is served by X (`_REF_SPAN_FOR`).
        """
        by_canonical: dict[str, list[Contig]] = {}
        for contig in self.contigs(digest):
            by_canonical.setdefault(canonical_chrom(contig.name), []).append(contig)
        for chrom, served_by in _REF_SPAN_FOR.items():
            by_canonical[chrom] = by_canonical.get(served_by, [])
        return by_canonical

    def assembly_label(self, digest: str) -> str | None:
        """The collection's assembly: from contig 1's length, else the one label asserted for it."""
        measured = self._label_from_lengths(digest)
        if measured is not None:
            return measured
        asserted = self.asserted_labels(digest)
        return asserted[0] if len(asserted) == 1 else None

    def asserted_labels(self, digest: str) -> list[str]:
        """Labels a caller asserted for the collection with `add_fasta(…, assembly_label=…)`."""
        return sorted(alias for namespace, alias in self.store.get_aliases_for_collection(digest) if namespace == ASSEMBLY_NAMESPACE)

    def fetch(self, contig: Contig, ranges: list[tuple[int, int]]) -> list[str]:
        """The bases over 0-based half-open `ranges` of `contig`, in order. Uppercase."""
        if not ranges:
            return []
        return self.store.get_substrings(contig.digest, ranges)

    def _label_from_lengths(self, digest: str) -> str | None:
        level2 = self.store.get_collection_level2(digest)
        lengths = {canonical_chrom(name): int(length) for name, length in zip(level2["names"], level2["lengths"])}
        return CHR1_LENGTHS.get(lengths.get("1", -1))

    def _namespaces(self) -> list[str]:
        return list(self.store.list_collection_alias_namespaces())


@dataclass(frozen=True)
class _Plan:
    """What one collection can do for one facet."""

    digest: str
    label: str | None
    checks: dict[str, ContigCheck]
    contigs: dict[str, Contig]
    detail: dict[str, str]
    refusal: str | None
    readable_rows: int


def resolve_ref(
    loaded: Loaded,
    reference: "ReferenceStore | Path | str | None" = None,
    *,
    collection: str | None = None,
) -> Loaded:
    """The facet with `ref` read from its assembly sequence, `ref_check` compared, and array alleles derived.

    `reference` is a `ReferenceStore`, a store directory, or a FASTA file (added to the store at
    `reference_dir()` on first use, and then the only collection considered). Without one, the store
    at `reference_dir()` is used. The store's collection that matches the facet is chosen by the
    contig checks; `collection` names one instead.

    The result is the same facet (same `facet_id`) with the slots SCHEMA.md leaves open until the
    sequence is read now filled, and `facets.ref_seqcol_digest` naming the collection that was read.
    The input `Loaded` is not changed. Rows that could not be read keep `ref` null (and `ref_check`
    `not_run`), each reason counted once in the warnings.

    Raises `ReferenceMismatchError` when no collection can be shown to be the facet's assembly, when
    a primary contig's stated length or md5 disagrees with the reference, or when more than one
    collection fits equally and `collection` does not choose.
    """
    versions = {"SCHEMA_VERSION": SCHEMA_VERSION, "RESOLVE_VERSION": RESOLVE_VERSION, "SPILL_LAYOUT": cache.SPILL_LAYOUT}
    cache.sweep_stale(cache.RESOLVE_SPILL, versions)
    facets = loaded.facets.collect()
    if facets.height != 1:
        raise ValueError(f"resolve_ref takes one facet; this Loaded has {facets.height}")
    facet = facets.row(0, named=True)
    if facet["ref_seqcol_digest"] is not None:
        raise ValueError(f"facet {facet['facet_id']}: REF was already read (from {facet['ref_seqcol_digest']}); resolve the facet as loaded")
    store, given = _store_for(reference)
    collection = collection if collection is not None else given

    rows_per_chrom = dict(loaded.calls.group_by("chrom").agg(pl.len()).collect().iter_rows())
    candidates = [collection] if collection is not None else store.collections()
    if not candidates:
        raise ReferenceMismatchError("the reference store holds no collection; add a FASTA with ReferenceStore.add_fasta")
    inferred = _infer_from_rows(store, loaded, facet, candidates)
    if inferred is not None:
        facet = {**facet, "assembly_label": inferred.label, "assembly_label_basis": AssemblyLabelBasis.INFERRED.value}
    plans = [_plan(store, digest, facet, rows_per_chrom) for digest in candidates]
    plan = _choose(plans, facet["facet_id"], explicit=collection is not None)

    warnings = [f"assembly inferred from the rows as {inferred.label}: {inferred.summary()}"] if inferred is not None else []
    warnings += _contig_warnings(plan, rows_per_chrom)
    strand = _check_strand(store, loaded, facet, plan.digest)
    if strand is not None:
        warnings.append(strand.summary())
    withhold = strand is not None and strand.check is StrandCheck.MIXED
    concatenated = facet["format_id"] in CONCATENATED_ALLELES
    entry = cache_dir() / "frames" / f"{facet['facet_id']}-{SCHEMA_VERSION}-ref-{plan.digest}-{RESOLVE_VERSION}"
    with cache.spilling(entry, cache.RESOLVE_SPILL, versions, cache.derived_origin(facet, facet["facet_id"])) as spill:
        # The input is streamed to parquet once, in its canonical (contig-grouped) order, so each contig
        # below is a cheap filtered scan rather than a replay of whatever plan built the input. Its
        # bytes also key the result: a run on the same input reuses what an earlier one wrote.
        staged = spill.stage / "input.parquet"
        loaded.calls.sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
        key = cache.content_key([staged], facet)
        done = spill.find(key)
        if done is None:
            counts: Counter[str] = Counter()
            source = pl.scan_parquet(staged)
            parts: list[Path] = []
            for index, chrom in enumerate(sorted(rows_per_chrom, key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))):
                rows = source.filter(pl.col("chrom") == chrom).collect()
                contig = plan.contigs.get(chrom) if plan.checks.get(chrom) in _READABLE else None
                resolved = _resolve_rows(rows, store, contig, counts, withhold, concatenated) if contig is not None else rows
                part = spill.stage / f"calls-{index:05d}.parquet"
                canonical_order(resolved.lazy()).collect().write_parquet(part, row_group_size=SPILL_ROW_GROUP)
                parts.append(part)
            staged.unlink()
            done = spill.commit(key, {"calls": parts}, _row_warnings(counts, plan.digest))
    warnings += done.warnings
    parts = done.files("calls")
    for warning in warnings:
        logger.warning("%s: %s", facet["facet_id"], warning)

    labelled = (
        [pl.lit(facet["assembly_label"]).alias("assembly_label"), pl.lit(facet["assembly_label_basis"]).cast(FACETS["assembly_label_basis"]).alias("assembly_label_basis")]
        if inferred is not None
        else []
    )
    return Loaded(
        facets=loaded.facets.with_columns(
            pl.lit(plan.digest).alias("ref_seqcol_digest"),
            pl.lit(strand.check.value if strand is not None and strand.check is not None else None)
            .cast(FACETS["strand_check"])
            .alias("strand_check"),
            *labelled,
        ),
        calls=pl.scan_parquet(parts) if parts else empty(CALLS),
        blocks=loaded.blocks,
        links=loaded.links,
        warnings=loaded.warnings + tuple(warnings),
        schema_version=loaded.schema_version,
    )


def _infer_from_rows(store: ReferenceStore, loaded: Loaded, facet: dict, candidates: list[str]) -> RowsInference | None:
    """The facet's assembly from its own rows (F9), for a facet that states neither a label nor contig lengths.

    None when the facet needs no inference (it states one or the other). Raises `ReferenceMismatchError`
    when it needs one and its rows do not settle it, with the scores, rather than reading an unchecked
    reference. Only array rows (no source REF) are scored: for a VCF, `ref_check` would then confirm
    the very comparison that chose its build.
    """
    stated_lengths = any(entry["length"] is not None for entry in facet["header_contigs"] or [])
    if facet["assembly_label"] is not None or stated_lengths or facet["vcf_version"] is not None:
        return None
    sample = sample_snv_rows(loaded.calls)
    bases = {digest: _bases_at(store, digest, sample) for digest in candidates}
    inference = infer_assembly_from_rows(sample, bases, {digest: store.assembly_label(digest) for digest in candidates})
    if not inference:
        scores = f": {inference.summary()}" if inference.shares else ""
        raise ReferenceMismatchError(
            f"facet {facet['facet_id']}: the facet states neither contig lengths nor an assembly, and its rows do not "
            f"settle one ({inference.reason}){scores}"
        )
    return inference


def _check_strand(store: ReferenceStore, loaded: Loaded, facet: dict, digest: str) -> StrandResult | None:
    """The strand an array facet's genotypes are on (F13), from the same row sample F9 scores, against
    the collection REF is read from. None for a VCF, whose REF the source wrote."""
    if facet["vcf_version"] is not None:
        return None
    sample = sample_snv_rows(loaded.calls)
    if sample.is_empty():
        return None
    return check_strand(sample, _bases_at(store, digest, sample))


def _bases_at(store: ReferenceStore, digest: str, sample: pl.DataFrame) -> pl.Series:
    """The collection's base at each sampled row: null on a contig it lacks or names twice, `""` past its end."""
    by_canonical: dict[str, list[Contig]] = {}
    for contig in store.contigs(digest):
        by_canonical.setdefault(canonical_chrom(contig.name), []).append(contig)
    keyed = sample.select("chrom", "pos").with_row_index("_row")
    parts = []
    for chrom in sorted(keyed["chrom"].unique().to_list()):
        matches = by_canonical.get(chrom, [])
        if len(matches) != 1:
            continue
        rows = keyed.filter(pl.col("chrom") == chrom)
        start, end = _span(pl.col("pos"), pl.lit(1))
        spans = rows.with_columns(start.alias("_start"), end.alias("_end"))
        read = _read(spans, store, matches[0]).fill_null("")
        parts.append(rows.select("_row").with_columns(read.alias("_base")))
    found = pl.concat(parts) if parts else pl.DataFrame(schema={"_row": pl.UInt32, "_base": pl.String})
    return keyed.select("_row").join(found, on="_row", how="left", maintain_order="left")["_base"]


def _store_for(reference: "ReferenceStore | Path | str | None") -> tuple[ReferenceStore, str | None]:
    """The store to read, and the collection to use when `reference` was one FASTA file."""
    if isinstance(reference, ReferenceStore):
        return reference, None
    if reference is None:
        return ReferenceStore.open(), None
    path = Path(reference)
    if RefgetStore.store_exists(str(path)):
        return ReferenceStore.open(path), None
    if path.is_file():
        store = ReferenceStore.open()
        return store, store.add_fasta(path)
    raise FileNotFoundError(f"{path}: neither a reference store nor a FASTA file")


def _plan(store: ReferenceStore, digest: str, facet: dict, rows_per_chrom: dict[str, int]) -> _Plan:
    """Check every contig the facet's calls use against one collection."""
    label = store.assembly_label(digest)
    by_canonical = store.contigs_by_chrom(digest)
    stated: dict[str, list[dict]] = {}
    sources = {entry["source"]: entry["chrom"] for entry in facet["contig_map"] or []}
    for entry in facet["header_contigs"] or []:
        if entry["length"] is not None:
            stated.setdefault(sources.get(entry["name"], canonical_chrom(entry["name"])), []).append(entry)

    facet_label = facet["assembly_label"]
    checks: dict[str, ContigCheck] = {}
    chosen: dict[str, Contig] = {}
    detail: dict[str, str] = {}
    refusal = None
    for chrom in rows_per_chrom:
        matches = by_canonical.get(chrom, [])
        if not matches:
            checks[chrom] = ContigCheck.NOT_IN_REFERENCE
            continue
        if len(matches) > 1:
            checks[chrom] = ContigCheck.AMBIGUOUS
            detail[chrom] = " and ".join(sorted(contig.name for contig in matches))
            continue
        contig = matches[0]
        chosen[chrom] = contig
        if chrom in stated:
            lengths = {entry["length"] for entry in stated[chrom]}
            md5s = {entry["md5"].lower() for entry in stated[chrom] if entry["md5"]}
            if lengths != {contig.length}:
                checks[chrom] = ContigCheck.LENGTH_MISMATCH
                detail[chrom] = f"stated {', '.join(map(str, sorted(lengths)))}, reference {contig.length}"
            elif md5s and md5s != {contig.md5}:
                checks[chrom] = ContigCheck.MD5_MISMATCH
            else:
                checks[chrom] = ContigCheck.LENGTH
            if checks[chrom] is not ContigCheck.LENGTH and chrom in PRIMARY and refusal is None:
                refusal = f"contig {chrom}: {checks[chrom].value} ({detail.get(chrom, 'md5 differs')}), so the reference is another assembly"
        elif facet_label is not None and label is not None and facet_label == label:
            unverifiable = chrom in LABEL_UNVERIFIABLE.get(label, frozenset())
            checks[chrom] = ContigCheck.LABEL_UNVERIFIABLE if unverifiable else ContigCheck.LABEL
        else:
            checks[chrom] = ContigCheck.UNVERIFIED
    readable = sum(rows_per_chrom[chrom] for chrom, check in checks.items() if check in _READABLE)
    if refusal is None and rows_per_chrom and readable == 0:
        if facet_label is not None and label is not None and facet_label != label:
            refusal = f"the facet's assembly is {facet_label} and the reference is {label}"
        elif not stated and facet_label is None:
            refusal = "the facet states neither contig lengths nor an assembly, so no reference can be shown to be its own"
        elif not stated and label is None:
            refusal = "the facet states no contig lengths and the reference's assembly is unknown (add_fasta(…, assembly_label=…))"
        else:
            refusal = "no contig of the facet's calls could be matched with the reference"
    return _Plan(digest, label, checks, chosen, detail, refusal, readable)


def _choose(plans: list[_Plan], facet_id: str, explicit: bool) -> _Plan:
    usable = [plan for plan in plans if plan.refusal is None]
    if not usable:
        reasons = "; ".join(f"{plan.digest}: {plan.refusal}" for plan in plans)
        raise ReferenceMismatchError(f"facet {facet_id}: no reference collection fits ({reasons})")
    if explicit:
        return usable[0]
    best = max(plan.readable_rows for plan in usable)
    top = [plan for plan in usable if plan.readable_rows == best]
    if len(top) > 1:
        raise ReferenceMismatchError(
            f"facet {facet_id}: {len(top)} collections fit equally ({', '.join(plan.digest for plan in top)}); pass collection= to choose"
        )
    return top[0]


def _contig_warnings(plan: _Plan, rows_per_chrom: dict[str, int]) -> list[str]:
    """One line per kind of contig that will not be read, with its rows and contig names."""
    grouped: dict[ContigCheck, list[str]] = {}
    for chrom, check in plan.checks.items():
        if check not in _READABLE:
            grouped.setdefault(check, []).append(chrom)
    explain = {
        ContigCheck.NOT_IN_REFERENCE: "contig(s) the reference does not have",
        ContigCheck.AMBIGUOUS: "contig(s) the reference names twice",
        ContigCheck.LENGTH_MISMATCH: "contig(s) whose stated length differs from the reference",
        ContigCheck.MD5_MISMATCH: "contig(s) whose stated md5 differs from the reference",
        ContigCheck.LABEL_UNVERIFIABLE: f"contig(s) the {plan.label} label cannot vouch for without a stated length",
        ContigCheck.UNVERIFIED: "contig(s) with no stated length and no agreeing assembly label",
    }
    lines = []
    for check in ContigCheck:
        names = sorted(grouped.get(check, []), key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))
        if not names:
            continue
        rows = sum(rows_per_chrom[name] for name in names)
        shown = ", ".join(f"{name} ({plan.detail[name]})" if name in plan.detail else name for name in names[:10])
        more = f" and {len(names) - 10} more" if len(names) > 10 else ""
        lines.append(f"{rows} row(s) left with REF unread, on {explain[check]}: {shown}{more}")
    return lines


def _span(pos: pl.Expr, length: pl.Expr) -> tuple[pl.Expr, pl.Expr]:
    """1-based `pos` and a REF length as the store's 0-based half-open `[start, end)`."""
    start = pos.cast(pl.Int64) - 1
    return start, start + length.cast(pl.Int64)


def _read(rows: pl.DataFrame, store: ReferenceStore, contig: Contig) -> pl.Series:
    """The reference bases over each row's `_start`/`_end`, null where the span leaves the contig."""
    inside = (pl.col("_start") >= 0) & (pl.col("_end") <= contig.length) & (pl.col("_end") > pl.col("_start"))
    spans = rows.filter(inside).select("_start", "_end").unique().sort("_start", "_end")
    bases = store.fetch(contig, list(zip(spans["_start"].to_list(), spans["_end"].to_list())))
    table = spans.with_columns(pl.Series("_bases", bases, dtype=pl.String))
    return rows.select("_start", "_end").join(table, on=["_start", "_end"], how="left", maintain_order="left")["_bases"]


def _resolve_rows(
    rows: pl.DataFrame, store: ReferenceStore, contig: Contig, counts: Counter[str], withhold: bool = False, concatenated: bool = False
) -> pl.DataFrame:
    """One contig's rows with REF read; `counts` gathers the outcome per reason. `withhold`: the facet's
    strand is mixed (F13), so array rows whose genotype is not strand-safe get no `alts` or `gt`.
    `concatenated`: the format writes an indel genotype as its two allele sequences back to back (F4)."""
    if rows.is_empty():
        return rows
    is_vcf = pl.col("source_ref").is_not_null()
    called = pl.col("gt_alleles").list.drop_nulls()
    single_bases = pl.col("gt_alleles").list.eval(pl.element().is_null() | pl.element().is_in(_BASES)).list.all()
    array_ok = pl.col("source_ref").is_null() & pl.col("gt_alleles").is_not_null() & (called.list.len() > 0) & single_bases
    length = pl.when(is_vcf).then(pl.col("source_ref").str.len_chars()).when(array_ok).then(pl.lit(1)).otherwise(None)
    start, end = _span(pl.col("pos"), length)
    work = rows.with_row_index("_row").with_columns(start.alias("_start"), end.alias("_end"), is_vcf.alias("_vcf"), array_ok.alias("_array"))
    work = work.with_columns(pl.when(pl.col("_start").is_not_null()).then(_read(work, store, contig)).otherwise(None).alias("_bases"))
    work = _split_concatenated(work, store, contig, counts, withhold) if concatenated else work.with_columns(
        *[pl.lit(None, dtype=CALLS[name]).alias(f"_split_{name}") for name in ("ref", "alts", "gt", "end_pos")]
    )

    vcf = work.filter(pl.col("_vcf"))
    if vcf.height:
        _count_vcf(vcf, counts)
    array_rows = work.filter(~pl.col("_vcf"))
    if array_rows.height:
        _count_array_unreadable(array_rows, counts)

    agree = pl.col("source_ref").str.to_uppercase() == pl.col("_bases")
    work = work.with_columns(
        pl.coalesce("_bases", "_split_ref").alias("ref"),
        pl.when(pl.col("_vcf") & pl.col("_bases").is_not_null())
        .then(pl.when(agree).then(pl.lit(RefCheck.MATCH.value)).otherwise(pl.lit(RefCheck.MISMATCH.value)))
        .otherwise(pl.col("ref_check").cast(pl.String))
        .cast(CALLS["ref_check"])
        .alias("ref_check"),
    )
    work = _derive_array_alleles(work, counts, withhold)
    split = pl.col("_split_ref").is_not_null()
    work = work.with_columns(
        *[pl.when(split).then(pl.col(f"_split_{name}")).otherwise(pl.col(name)).alias(name) for name in ("alts", "gt", "end_pos")]
    )
    return work.select(CALLS.names())


def _split_concatenated(work: pl.DataFrame, store: ReferenceStore, contig: Contig, counts: Counter[str], withhold: bool) -> pl.DataFrame:
    """Living DNA's indel genotypes (F4) split against the reference.

    The text is two allele sequences back to back (`CCA`), split nowhere by the source. A split is
    taken when it is the only one, of all the points where both pieces start with the same base (VCF's
    padding base), at which one piece equals the reference over its own span from `pos`. Several such
    splits (a repeat: `AAAA` over `AAA…` fits `A/AAA`, `AA/AA` and `AAA/A`) or none leave the row as
    loaded, counted. A split row gets `gt_alleles`, in the order written. It gets `ref`, `alts`, `gt`
    and `end_pos` only when exactly one of its two alleles equals the reference over its own span: a
    biallelic site's REF is one of its alleles, so the other is ALT. Where both do (`C/CA` over
    `CA…`: an insertion of A and a deletion of A fit alike), or the call is homozygous (`CA/CA` shows
    one of the site's two alleles, so REF's span, and whether the call is reference, is unknown), REF
    stays unread. In a facet of mixed strand (F13) nothing is split. Returns `work` with `gt_alleles`
    filled where split, and `_split_ref`, `_split_alts`, `_split_gt`, `_split_end_pos` set where REF is known.
    """
    target = (
        pl.col("source_ref").is_null() & pl.col("gt_alleles").is_null() & pl.col("source_gt").str.contains(_CONCATENATED)
    ).fill_null(False)
    rows = work.filter(target).select("_row", "pos", "source_gt")
    split_columns: dict[str, list] = {"_row": [], "_gt_alleles": [], "_split_ref": [], "_split_alts": [], "_split_gt": [], "_split_end_pos": []}
    if rows.height and withhold:
        counts["concat_mixed_strand"] += rows.height
    elif rows.height:
        texts = rows["source_gt"].to_list()
        starts = [pos - 1 for pos in rows["pos"].to_list()]
        spans = [(start, min(start + len(text), contig.length)) for start, text in zip(starts, texts)]
        inside = [index for index, (start, stop) in enumerate(spans) if stop > start]
        windows = [""] * len(texts)
        for index, bases in zip(inside, store.fetch(contig, [spans[index] for index in inside])):
            windows[index] = bases.upper()
        for row, pos, text, window in zip(rows["_row"].to_list(), rows["pos"].to_list(), texts, windows):
            fits = [
                (text[:k], text[k:])
                for k in range(1, len(text))
                if text[0] == text[k] and window.startswith((text[:k], text[k:]))
            ]
            if len(fits) != 1:
                counts["concat_several" if fits else "concat_no_fit"] += 1
                continue
            first, second = fits[0]
            matches = [allele for allele in (first, second) if window.startswith(allele)]
            ref = None
            if first == second:
                counts["concat_homozygous"] += 1
            elif len(matches) == 2:
                counts["concat_repeat"] += 1
            else:
                counts["concat_resolved"] += 1
                ref = matches[0]
            alt = next((allele for allele in (first, second) if allele != ref), None)
            for name, value in zip(
                split_columns,
                (row, [first, second], ref, None if ref is None else [alt], None if ref is None else [int(a != ref) for a in (first, second)], None if ref is None else pos + len(ref) - 1),
            ):
                split_columns[name].append(value)
    found = pl.DataFrame(
        split_columns,
        schema={
            "_row": pl.UInt32,
            "_gt_alleles": CALLS["gt_alleles"],
            "_split_ref": CALLS["ref"],
            "_split_alts": CALLS["alts"],
            "_split_gt": CALLS["gt"],
            "_split_end_pos": CALLS["end_pos"],
        },
    )
    return work.join(found, on="_row", how="left", maintain_order="left").with_columns(
        pl.coalesce("gt_alleles", "_gt_alleles").alias("gt_alleles")
    ).drop("_gt_alleles")


def _count_vcf(vcf: pl.DataFrame, counts: Counter[str]) -> None:
    read = vcf.filter(pl.col("_bases").is_not_null())
    counts["vcf_rows"] += vcf.height
    counts["vcf_outside"] += vcf.height - read.height
    agree = read.filter(pl.col("source_ref").str.to_uppercase() == pl.col("_bases"))
    counts["vcf_match"] += agree.height
    mismatched = read.filter(pl.col("source_ref").str.to_uppercase() != pl.col("_bases"))
    counts["vcf_mismatch"] += mismatched.height
    for source, ref in mismatched.select("source_ref", "_bases").iter_rows():
        if _iupac_explains(source.upper(), ref):
            counts["vcf_mismatch_iupac"] += 1


def _iupac_explains(source: str, ref: str) -> bool:
    """Whether every base where `source` differs from `ref` is an ambiguity code covering the source base."""
    if len(source) != len(ref):
        return False
    differing = [(s, r) for s, r in zip(source, ref) if s != r]
    return bool(differing) and all(s in _IUPAC.get(r, frozenset()) for s, r in differing)


def _count_array_unreadable(rows: pl.DataFrame, counts: Counter[str]) -> None:
    counts["array_rows"] += rows.height
    if rows.filter(pl.col("_array")).height:
        counts["array_outside"] += rows.filter(pl.col("_array") & pl.col("_bases").is_null()).height
    rest = rows.filter(~pl.col("_array"))
    if not rest.height:
        return
    codes = rest["gt_alleles"].list.eval(pl.element().is_in(["I", "D"])).list.any().fill_null(False)
    no_call = rest["gt_alleles"].is_not_null() & (rest["gt_alleles"].list.drop_nulls().list.len() == 0)
    counts["array_code"] += int(codes.sum())
    counts["array_no_call"] += int((no_call & ~codes).sum())
    counts["array_unreadable"] += int(rest["gt_alleles"].is_null().sum())


def _derive_array_alleles(work: pl.DataFrame, counts: Counter[str], withhold: bool = False) -> pl.DataFrame:
    """For array rows whose base was read: `alts` the observed non-REF alleles, sorted; `gt` into `[ref] + alts`.

    With `withhold` (a mixed-strand facet, F13) a row whose genotype is not strand-safe keeps `ref` and
    gets neither, each reason counted.
    """
    target = pl.col("_array") & pl.col("ref").is_not_null()
    if withhold:
        reason = withheld_reason(pl.col("gt_alleles"), pl.col("ref"))
        work = work.with_columns(pl.when(target).then(reason).otherwise(None).alias("_withheld"))
        for name, count in work.filter(pl.col("_withheld").is_not_null()).group_by("_withheld").len().iter_rows():
            counts[f"strand_{name}"] += count
        target = target & pl.col("_withheld").is_null()
    if not work.select(target.any()).item():
        return work
    alleles = pl.col("gt_alleles")
    work = work.with_columns(
        pl.when(target)
        .then(alleles.list.drop_nulls().list.set_difference(pl.concat_list("ref")).list.unique().list.sort())
        .otherwise(pl.col("alts"))
        .alias("alts")
    )
    # Vectorised over haplotype slots and alt positions: each allele is a single base, so a slot's
    # index is 0 when it is REF and j + 1 when it is the j-th sorted alt.
    widths = work.filter(target).select(alleles.list.len().max().alias("slots"), pl.col("alts").list.len().max().alias("alts")).row(0)
    indexes = []
    for slot in range(widths[0]):
        allele = alleles.list.get(slot, null_on_oob=True)
        index = pl.when(allele.is_null()).then(None).when(allele == pl.col("ref")).then(0)
        for position in range(widths[1]):
            index = index.when(allele == pl.col("alts").list.get(position, null_on_oob=True)).then(position + 1)
        indexes.append(index.otherwise(None).cast(pl.UInt8))
    gt = pl.concat_list(indexes).list.head(alleles.list.len())
    work = work.with_columns(pl.when(target).then(gt).otherwise(pl.col("gt")).alias("gt"))
    derived = work.filter(target)
    counts["array_read"] += derived.height
    counts["array_ref_not_acgt"] += derived.filter(~pl.col("ref").is_in(_BASES)).height
    counts["array_no_ref_allele"] += derived.filter(~pl.col("gt").list.contains(0)).height
    return work


def _row_warnings(counts: Counter[str], digest: str) -> list[str]:
    lines = []
    if counts["vcf_rows"]:
        line = (
            f"REF read from collection {digest} for {counts['vcf_match'] + counts['vcf_mismatch']} of "
            f"{counts['vcf_rows']} row(s) with a source REF on readable contigs: {counts['vcf_match']} match, "
            f"{counts['vcf_mismatch']} mismatch"
        )
        if counts["vcf_mismatch_iupac"]:
            line += f" ({counts['vcf_mismatch_iupac']} of them against an IUPAC ambiguity code, which VCF writes as a concrete base)"
        lines.append(line)
    if counts["vcf_outside"]:
        lines.append(f"{counts['vcf_outside']} row(s) with a source REF reaching past their contig, or at POS 0: REF left unread")
    if counts["array_rows"]:
        lines.append(
            f"REF read from collection {digest} for {counts['array_read']} of {counts['array_rows']} array row(s) "
            "on readable contigs; ALT derived from the genotype against it"
        )
    reasons = {
        "array_code": "with an I/D indel code (resolving it needs the chip manifest)",
        "array_no_call": "with no allele called (the probe's REF span is unknown)",
        "array_unreadable": "whose genotype could not be split into alleles",
        "array_outside": "past the end of their contig, or at position 0",
    }
    for key, reason in reasons.items():
        if counts[key]:
            lines.append(f"{counts[key]} array row(s) {reason}: REF left unread")
    split = counts["concat_resolved"] + counts["concat_homozygous"] + counts["concat_repeat"]
    if split:
        lines.append(
            f"{split} genotype(s) written as two allele sequences back to back split where exactly one split fits the "
            f"reference (both alleles start with the same base, one equals the reference over its own span): "
            f"{counts['concat_resolved']} with REF and ALT read; REF left unread on {counts['concat_homozygous']} "
            f"homozygous call(s) (one of the site's two alleles shown, so REF's span is unknown) and on "
            f"{counts['concat_repeat']} where both alleles equal the reference over their own spans (an insertion and a "
            "deletion fit alike)"
        )
    concat = {
        "concat_several": "fit the reference at more than one split (a repeat)",
        "concat_no_fit": "fit the reference at no split",
        "concat_mixed_strand": "were not split, the facet's strand being mixed",
    }
    for key, reason in concat.items():
        if counts[key]:
            lines.append(f"{counts[key]} genotype(s) written as two allele sequences back to back {reason}: left unsplit")
    if counts["array_ref_not_acgt"]:
        lines.append(f"{counts['array_ref_not_acgt']} array row(s) where the reference base is not A, C, G or T")
    withheld = {
        "strand_minus": "heterozygous rows carrying the reference base's complement and not the base (minus-strand probes)",
        "strand_neither": "heterozygous rows carrying neither the reference base nor its complement",
        "strand_undecidable": "homozygous rows and half-calls, which read differently on the two strands",
    }
    for key, reason in withheld.items():
        if counts[key]:
            lines.append(f"{counts[key]} array row(s) with alts and gt withheld, the facet's strand being mixed: {reason}")
    if counts["array_no_ref_allele"]:
        share = counts["array_no_ref_allele"] / counts["array_read"]
        lines.append(
            f"{counts['array_no_ref_allele']} array row(s) ({share:.1%} of those read) carry no REF allele in their "
            "genotype: homozygous non-reference, or reported on the other strand or another build"
        )
    return lines

