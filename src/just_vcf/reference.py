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
  the two alleles say which of them is the reference;
- an insertion point (an empty `source_ref`: CGI's native files, a CASAVA indel table; a Genos row of
  zero length): the anchor base VCF writes, the reference base at `pos` (the base after, at a contig's
  start), is `ref` and is written into every allele (`_anchor_points`); a point that calls no inserted
  allele (every called haplotype the empty allele: "no insertion here") is left unread;
- a Genos row (`SPANNED_GENOTYPES`): REF over the span the row states, whatever its length, and its
  alleles derived against it as an array row's are.

Before any base is read the reference has to be shown to be the facet's assembly, contig by contig
(`ContigCheck`). A facet that stated contig lengths (VCF `##contig`) is checked by length and md5;
one that did not (arrays) by its assembly label against the label the reference's own lengths give.
An array that states no label either has one inferred from its rows against each collection first
(F9, `just_vcf.assembly.infer_assembly_from_rows`), or is refused with the scores.
An autosome of another length or md5 means another assembly, and the whole facet is refused
(`ReferenceMismatchError`). Any other contig that cannot be shown to match is left unresolved, with
its rows counted in the warnings: X, Y and MT among them, where assemblies of one build legitimately
differ (F20). Nothing falls back to an unchecked read.

Coordinates: `pos` is 1-based; the store takes 0-based half-open ranges. `_span` is the one place
that converts between them (charter 13).
"""

import hashlib
import json
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
from just_vcf.contig_names import scaffold_names
from just_vcf.loading import SPILL_ROW_GROUP, Loaded, cache_dir
from just_vcf.probes import IndelSite, Manifest, catalog_digest, id_calibration, indel_catalog, locate_indel, place_indel, snv_catalog, window_for
from just_vcf.schema import CALLS, CONTIG_RANK, EMPTY_ALLELE, FACETS, SCHEMA_VERSION, AssemblyLabelBasis, RefCheck, StrandCheck, canonical_chrom, canonical_order, empty
from just_vcf.strand import COMPLEMENT, StrandResult, check_strand, withheld_reason

logger = logging.getLogger(__name__)

RESOLVE_VERSION = 12
"""Part of the spill key: bumped whenever the rule changes what `resolve_ref` writes (10: an insertion
point's anchor base, and REF over a Genos row's stated span; 11: a point whose every called haplotype is
the empty allele is left unread, not anchored as a hom-ref call at the base before it; 12: the row sample
the build inference and the strand check score leaves out rows not over one base, a Genos insertion site
among them, so a facet's inferred label and `strand_check`, and the rows a mixed strand withholds, can
change; the strand is not in the spill key)."""

CONCATENATED_ALLELES: frozenset[str] = frozenset({"living_dna"})
"""Formats whose multi-letter genotype text is the two allele sequences written back to back, as the
vendor states it ("both allele sequences concatenated"). Only these are split against the reference:
another format's long genotype text may mean something else."""

SPANNED_GENOTYPES: frozenset[str] = frozenset({"genos"})
"""Genotype-table formats whose rows state their reference span (Genos's 0-based half-open BEGIN/END,
loaded as `pos`..`end_pos`) and spell each allele as a plus-strand sequence, `""` the empty allele. REF
is read over that span, whatever its length, and a row of zero length (an insertion site: `end_pos`
null) is anchored. Only these: in an array export a null `end_pos` means the span is unknown (`I`/`D`
codes), not that it is empty."""

_SEQUENCE_ALLELE = r"^[ACGTN]*$"
"""A spanned row's allele that REF and ALT can be read against: bases, or the empty allele."""

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
"""The primary chromosomes. `concordance` and `liftover` take one of these whose sequence differs as
another assembly; `resolve_ref` takes only an autosome so (`AUTOSOMES`)."""

AUTOSOMES: frozenset[str] = frozenset(str(n) for n in range(1, 23))
"""Contigs whose stated length or md5 differing from the reference means another assembly, so
`resolve_ref` refuses the facet. X, Y and MT are not among them: assemblies of one build differ there
legitimately (a PAR N-masked in one release and not in another; Complete Genomics' GRCh37 chrY has
another md5 than Ensembl's, F20; rCRS against the Yoruba MT), so a mismatch there leaves that contig
unread and the matching contigs are read."""

LABEL_UNVERIFIABLE: dict[str, frozenset[str]] = {
    "GRCh36": frozenset({"MT"}),
    "GRCh37": frozenset({"MT"}),
    "T2T-CHM13": frozenset({"MT"}),
}
"""Contigs a label cannot vouch for: `GRCh37` covers both hg19, whose chrM is the Yoruba sequence
(16,571 bases), and GRCh37, whose MT is rCRS (16,569). Only a stated length tells them apart.
`GRCh36` likewise: NCBI36 and hg18 carry the Yoruba sequence, and an array labelled build 36 does not
say whether its MT positions are on it or on rCRS. `T2T-CHM13` covers `chm13v2.0.fa`, whose chrM is
the CHM13 cell line's own mitochondrion (CP068254.1), and the T2T analysis set's
`chm13v2.0_maskedY_rCRS.fa.gz`, which puts rCRS (NC_012920.1) there: both 16,569 bases, so even a
stated length does not tell them apart, only a stated md5 or the REF check."""

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

_MISMATCHES = frozenset({ContigCheck.LENGTH_MISMATCH, ContigCheck.MD5_MISMATCH})
"""Checks whose warning names each contig with its own row count and what was stated against the reference."""


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

        A GRCh37 or GRCh38 scaffold is also listed under the other names NCBI's assembly reports give
        it (F19, `contig_names.scaffold_names`): Ensembl's `KI270706.1` under UCSC's
        `chr1_KI270706v1_random` and RefSeq's `NT_187361.1`, and the other way round. Only where the
        sequence's length is the report's for that name, and only under a name the collection does not
        hold a sequence of its own by: a name is never taken from a pattern, and never over a contig the
        collection itself calls so.
        """
        by_canonical: dict[str, list[Contig]] = {}
        contigs = self.contigs(digest)
        for contig in contigs:
            by_canonical.setdefault(canonical_chrom(contig.name), []).append(contig)
        aliased: dict[str, list[Contig]] = {}
        for contig in contigs:
            for name in sorted(scaffold_names(contig.name, contig.length) - set(by_canonical)):
                aliased.setdefault(name, []).append(contig)
        by_canonical.update(aliased)
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
    sequence is read now filled, `facets.ref_seqcol_digest` naming the collection that was read,
    `facets.resolve_version` the `RESOLVE_VERSION` that read it, and `facets.resolve_inputs` a digest of the
    chip manifests it read an array's `I`/`D` codes or a mixed strand's allele pairs with (null where it
    read none).
    The input `Loaded` is not changed. Rows that could not be read keep `ref` null (and `ref_check`
    `not_run`), each reason counted once in the warnings.

    Raises `ReferenceMismatchError` when no collection can be shown to be the facet's assembly, when
    an autosome's stated length or md5 disagrees with the reference, or when more than one collection
    fits equally and `collection` does not choose. A stated X, Y or MT that disagrees leaves its rows
    unread, counted per contig in the warnings (F20).
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
    spanned = facet["format_id"] in SPANNED_GENOTYPES
    codes = _code_catalog(loaded, facet, plan.label)
    pairs = _pair_catalog() if withhold else None
    entry = cache_dir() / "frames" / f"{facet['facet_id']}-{SCHEMA_VERSION}-ref-{plan.digest}-{RESOLVE_VERSION}"
    with cache.spilling(entry, cache.RESOLVE_SPILL, versions, cache.derived_origin(facet, facet["facet_id"])) as spill:
        # The input is streamed to parquet once, in its canonical (contig-grouped) order, so each contig
        # below is a cheap filtered scan rather than a replay of whatever plan built the input. Its
        # bytes also key the result: a run on the same input reuses what an earlier one wrote.
        staged = spill.stage / "input.parquet"
        loaded.calls.sink_parquet(staged, row_group_size=SPILL_ROW_GROUP)
        key = cache.content_key([staged], {**facet, "probe_catalog": codes.digest, "pair_catalog": pairs.digest if pairs else None})
        done = spill.find(key)
        if done is None:
            counts: Counter[str] = Counter()
            source = pl.scan_parquet(staged)
            parts: list[Path] = []
            for index, chrom in enumerate(sorted(rows_per_chrom, key=lambda name: (CONTIG_RANK.get(name, len(CONTIG_RANK)), name))):
                rows = source.filter(pl.col("chrom") == chrom).collect()
                contig = plan.contigs.get(chrom) if plan.checks.get(chrom) in _READABLE else None
                resolved = _resolve_rows(rows, store, contig, counts, withhold, concatenated, codes, pairs, spanned) if contig is not None else rows
                part = spill.stage / f"calls-{index:05d}.parquet"
                canonical_order(resolved.lazy()).collect().write_parquet(part, row_group_size=SPILL_ROW_GROUP)
                parts.append(part)
            staged.unlink()
            done = spill.commit(key, {"calls": parts}, _row_warnings(counts, plan.digest, codes, pairs))
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
            # The rule that read the rows, since the facet keeps its id: a store tells two REF states apart by it.
            pl.lit(RESOLVE_VERSION, dtype=FACETS["resolve_version"]).alias("resolve_version"),
            # What else the rows were read with: the same facet resolved without the manifests keeps its
            # I/D rows unread, so a store must tell that state from this one.
            pl.lit(resolve_inputs(codes.digest, pairs.digest if pairs else None), dtype=FACETS["resolve_inputs"]).alias("resolve_inputs"),
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


def resolve_inputs(probe_catalog: str | None, pair_catalog: str | None) -> str | None:
    """`facets.resolve_inputs`: one digest of the manifests `resolve_ref` read rows with, the indel probes
    (`probe_catalog`) and the SNV allele pairs (`pair_catalog`), each `probes.catalog_digest` or None
    where none applied. None when neither did: the rows then depend on the collection and the
    `RESOLVE_VERSION` alone."""
    if probe_catalog is None and pair_catalog is None:
        return None
    text = json.dumps({"pair_catalog": pair_catalog, "probe_catalog": probe_catalog}, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest()[:16]


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
                detail[chrom] = f"length stated {', '.join(map(str, sorted(lengths)))}, reference {contig.length}"
            elif md5s and md5s != {contig.md5}:
                checks[chrom] = ContigCheck.MD5_MISMATCH
                detail[chrom] = f"md5 stated {', '.join(sorted(md5s))}, reference {contig.md5}"
            else:
                checks[chrom] = ContigCheck.LENGTH
            if checks[chrom] is not ContigCheck.LENGTH and chrom in AUTOSOMES and refusal is None:
                refusal = f"contig {chrom}: {checks[chrom].value} ({detail[chrom]}), so the reference is another assembly"
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
    """One line per kind of contig that will not be read, with its rows and contig names, and one for the
    scaffolds read under the reference's other name for them (F19)."""
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
        per_contig = check in _MISMATCHES
        shown = ", ".join(
            f"{name} ({rows_per_chrom[name]} rows; {plan.detail[name]})" if per_contig else f"{name} ({plan.detail[name]})" if name in plan.detail else name
            for name in names[:10]
        )
        more = f" and {len(names) - 10} more" if len(names) > 10 else ""
        lines.append(f"{rows} row(s) left with REF unread, on {explain[check]}: {shown}{more}")
    aliased = sorted(
        (chrom, contig.name)
        for chrom, contig in plan.contigs.items()
        if plan.checks.get(chrom) in _READABLE and chrom not in _REF_SPAN_FOR and canonical_chrom(contig.name) != chrom
    )
    if aliased:
        rows = sum(rows_per_chrom[chrom] for chrom, _ in aliased)
        shown = ", ".join(f"{chrom} as {name}" for chrom, name in aliased[:10])
        more = f" and {len(aliased) - 10} more" if len(aliased) > 10 else ""
        lines.append(f"{rows} row(s) read on scaffold(s) the reference names otherwise, as NCBI's assembly report pairs them: {shown}{more}")
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
    rows: pl.DataFrame,
    store: ReferenceStore,
    contig: Contig,
    counts: Counter[str],
    withhold: bool = False,
    concatenated: bool = False,
    codes: "_Codes | None" = None,
    pairs: "_Codes | None" = None,
    spanned: bool = False,
) -> pl.DataFrame:
    """One contig's rows with REF read; `counts` gathers the outcome per reason. `withhold`: the facet's
    strand is mixed (F13), so array rows whose genotype is not strand-safe get no `alts` or `gt`.
    `concatenated`: the format writes an indel genotype as its two allele sequences back to back (F4).
    `codes`: the manifests' indel probes an `I`/`D` row is resolved against (RM281). `pairs`: the
    manifests' SNV allele pairs a mixed facet's withheld rows are read against (RM281). `spanned`: the
    format states each row's REF span and spells its alleles as sequences (`SPANNED_GENOTYPES`)."""
    if rows.is_empty():
        return rows
    rows = _anchor_points(rows, store, contig, counts, spanned)
    # An empty source REF is an insertion point (`_anchor_points` read its anchor): no REF bases to compare.
    is_vcf = pl.col("source_ref").is_not_null() & (pl.col("source_ref") != EMPTY_ALLELE)
    called = pl.col("gt_alleles").list.drop_nulls()
    single_bases = pl.col("gt_alleles").list.eval(pl.element().is_null() | pl.element().is_in(_BASES)).list.all()
    # A spanned row (Genos) is read over the span it states, its alleles any sequence; an insertion site
    # has its span by now, the anchor base.
    stated_span = pl.lit(spanned) & pl.col("end_pos").is_not_null()
    sequences = pl.col("gt_alleles").list.eval(pl.element().is_null() | pl.element().str.contains(_SEQUENCE_ALLELE)).list.all()
    alleles_ok = pl.when(stated_span).then(sequences).otherwise(single_bases)
    array_ok = pl.col("source_ref").is_null() & pl.col("gt_alleles").is_not_null() & (called.list.len() > 0) & alleles_ok
    array_length = pl.when(stated_span).then(pl.col("end_pos").cast(pl.Int64) - pl.col("pos").cast(pl.Int64) + 1).otherwise(pl.lit(1))
    length = pl.when(is_vcf).then(pl.col("source_ref").str.len_chars()).when(array_ok).then(array_length).otherwise(None)
    start, end = _span(pl.col("pos"), length)
    work = rows.with_row_index("_row").with_columns(start.alias("_start"), end.alias("_end"), is_vcf.alias("_vcf"), array_ok.alias("_array"))
    work = work.with_columns(pl.when(pl.col("_start").is_not_null()).then(_read(work, store, contig)).otherwise(None).alias("_bases"))
    work = _split_concatenated(work, store, contig, counts, withhold) if concatenated else work.with_columns(
        *[pl.lit(None, dtype=CALLS[name]).alias(f"_split_{name}") for name in ("ref", "alts", "gt", "end_pos")]
    )
    work = _resolve_codes(work, store, contig, counts, codes)

    vcf = work.filter(pl.col("_vcf"))
    if vcf.height:
        _count_vcf(vcf, counts)
    # A spanned insertion site left unanchored (nothing called, no insertion, past the contig) was
    # counted by `_anchor_points`, under its own reason.
    unanchored_point = pl.lit(spanned) & pl.col("end_pos").is_null()
    array_rows = work.filter(~pl.col("_vcf") & pl.col("source_ref").is_null() & ~unanchored_point)
    if array_rows.height:
        _count_array_unreadable(array_rows, counts)

    agree = pl.col("source_ref").str.to_uppercase() == pl.col("_bases")
    work = work.with_columns(
        # `ref` is null on input (a facet already read is refused) but on an anchored insertion point.
        pl.coalesce("_bases", "_split_ref", "ref").alias("ref"),
        pl.when(pl.col("_vcf") & pl.col("_bases").is_not_null())
        .then(pl.when(agree).then(pl.lit(RefCheck.MATCH.value)).otherwise(pl.lit(RefCheck.MISMATCH.value)))
        .otherwise(pl.col("ref_check").cast(pl.String))
        .cast(CALLS["ref_check"])
        .alias("ref_check"),
    )
    work = _derive_array_alleles(work, counts, withhold, pairs)
    split = pl.col("_split_ref").is_not_null()
    work = work.with_columns(
        *[pl.when(split).then(pl.col(f"_split_{name}")).otherwise(pl.col(name)).alias(name) for name in ("alts", "gt", "end_pos")]
    )
    return work.select(CALLS.names())


_ANCHORED = ("pos", "end_pos", "ref", "ref_check", "alts", "gt_alleles", "source_alleles")
"""The columns `_anchor_points` rewrites on an insertion point it anchors."""


def _anchor_points(rows: pl.DataFrame, store: ReferenceStore, contig: Contig, counts: Counter[str], spanned: bool) -> pl.DataFrame:
    """Insertion points given VCF's anchor base, read from the reference.

    An insertion point is a row over no reference base: `source_ref` empty (CGI's native files, a
    CASAVA indel table), or, in a spanned format (`SPANNED_GENOTYPES`: Genos), a row with no `end_pos`.
    Every loader places it at the base before the point (1-based `pos = begin` of a zero-length 0-based
    locus), where VCF puts the anchor, so the anchor is the reference base at `pos`, prefixed onto every
    allele. A point before a contig's first base (`pos` 0) has no base before: VCF's rule takes the base
    after, appended to every allele, at `pos` 1. The row gets `end_pos` (the anchor's position) and its
    alleles in `gt_alleles`, `alts` and `source_alleles` anchored in place, so `gt` and the vectors stay
    indexed; a partly called allele (`?A`) gets the anchor too and stays one. A row with a source REF
    gets `ref` (the anchor) and `ref_check` null: the source wrote no base to compare. A Genos row gets
    REF from the span like any spanned row, after this. A point where no haplotype was called (`gt_alleles`
    all null: CGI's zero-length no-calls) has no allele to anchor, and REF with no ALT would state a
    no-call at the anchor base; it is left as loaded, counted. A point whose every called haplotype is
    the empty allele (CGI `""/""` and `""/?`, Genos `-/-`) says nothing is inserted there: evidence about
    the gap between two bases, not about the base before it. Anchored it would read `REF=X`, no ALT,
    `0/0`, a hom-ref call at that base the source never made, beside whatever the source did call
    there. No frame carries a statement about a gap without a genotype at a base (a `calls` row at `pos`
    is one, and a block spans bases), so it is left as loaded, counted apart. A point with an inserted
    allele on any haplotype is anchored: there the `0` is the haplotype without the insertion.
    """
    point = (pl.col("source_ref") == EMPTY_ALLELE).fill_null(False)
    if spanned:
        point = point | (pl.col("source_ref").is_null() & pl.col("end_pos").is_null() & pl.col("gt_alleles").is_not_null())
    indexed = rows.with_row_index("_p")
    targets = indexed.filter(point)
    if targets.is_empty():
        return rows
    called = pl.col("gt_alleles").list.drop_nulls()
    any_called = (called.list.len() > 0).fill_null(False)
    # Every called haplotype the empty allele: "no insertion here", which says nothing about the base before.
    no_insertion = any_called & called.list.eval(pl.element() == EMPTY_ALLELE).list.all().fill_null(False)
    counts["point_no_call"] += targets.filter(~any_called).height
    counts["point_no_insertion"] += targets.filter(no_insertion).height
    targets = targets.filter(any_called & ~no_insertion)
    inside = targets.filter(pl.col("pos").cast(pl.Int64) <= contig.length)
    counts["point_outside"] += targets.height - inside.height
    if inside.is_empty():
        return rows
    # The anchor: the base at 1-based `pos`, 0-based [pos - 1, pos); before the first base, base 1.
    starts = [max(pos, 1) - 1 for pos in inside["pos"].to_list()]
    spans = sorted(set(starts))
    bases = dict(zip(spans, (text.upper() for text in store.fetch(contig, [(start, start + 1) for start in spans]))))
    found: dict[str, list] = {"_p": [], **{f"_a_{name}": [] for name in _ANCHORED}}
    for (index, pos, source_ref, alts, gt_alleles, source_alleles), start in zip(
        inside.select("_p", "pos", "source_ref", "alts", "gt_alleles", "source_alleles").iter_rows(), starts
    ):
        anchor = bases[start]
        before = pos >= 1
        at = max(pos, 1)
        counts["point_anchored"] += 1
        counts["point_at_start"] += not before
        values = (
            index, at, at, anchor if source_ref is not None else None, None,
            _with_anchor(alts, anchor, before), _with_anchor(gt_alleles, anchor, before), _with_anchor(source_alleles, anchor, before),
        )  # fmt: skip
        for name, value in zip(found, values):
            found[name].append(value)
    fixed = pl.DataFrame(found, schema={"_p": pl.UInt32, **{f"_a_{name}": CALLS[name] for name in _ANCHORED}})
    hit = pl.col("_hit")
    joined = indexed.join(fixed.with_columns(pl.lit(True).alias("_hit")), on="_p", how="left", maintain_order="left")
    return joined.with_columns(
        *[pl.when(hit.fill_null(False)).then(pl.col(f"_a_{name}")).otherwise(pl.col(name)).alias(name) for name in _ANCHORED]
    ).select(rows.columns)


def _with_anchor(alleles: list[str | None] | None, anchor: str, before: bool) -> list[str | None] | None:
    """Each allele with the anchor base written before it (after it, at a contig's start); null stays null."""
    if alleles is None:
        return None
    return [None if allele is None else (anchor + allele if before else allele + anchor) for allele in alleles]


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


CODES = ["I", "D"]
"""An array's indel genotype letters: `I` the allele with the inserted bases, `D` the one without."""

_CODE_NAMES = ("ref", "alts", "gt", "end_pos", "gt_alleles")


@dataclass(frozen=True)
class _Codes:
    """The indel probes a facet's `I`/`D` rows are resolved against, or why there are none."""

    catalog: pl.DataFrame | None
    manifests: tuple[Manifest, ...]
    absent: str | None
    """Why `catalog` is None, for the warning; None also when the facet has no code rows at all."""
    unverified: str | None = None
    """Why the facet's chip has its codes withheld (`ID_CALIBRATION`): its I/D orientation is not
    verified. None when the chip's calibration passed, or there are no code rows."""

    @property
    def digest(self) -> str | None:
        return catalog_digest(list(self.manifests)) if self.manifests else None


def _code_catalog(loaded: Loaded, facet: dict, label: str | None) -> _Codes:
    """The configured manifests' indel probes on the facet's assembly, read only when it has `I`/`D` rows."""
    if facet["vcf_version"] is not None:
        return _Codes(None, (), None)
    code = pl.col("source_ref").is_null() & pl.col("gt_alleles").list.eval(pl.element().is_in(CODES)).list.any()
    if not loaded.calls.select(code.fill_null(False).any()).collect().item():
        return _Codes(None, (), None)
    calibration = id_calibration(facet["chip"])
    if not calibration.resolve:
        return _Codes(None, (), None, unverified=calibration.evidence if facet["chip"] is None else f"{facet['chip']}: {calibration.evidence}")
    if label is None:
        return _Codes(None, (), "the reference's assembly is unknown, so no manifest's build can be matched with it")
    catalog, manifests = indel_catalog(label)
    if not manifests:
        return _Codes(None, (), f"no manifest with indel probes on {label} is configured (JUST_VCF_MANIFESTS_DIR)")
    return _Codes(catalog, tuple(manifests), None)


def _pair_catalog() -> _Codes:
    """The configured manifests' SNV allele pairs on the plus strand, for a mixed-strand facet (RM281)."""
    catalog, manifests = snv_catalog()
    if not manifests:
        return _Codes(None, (), "no manifest is configured (JUST_VCF_MANIFESTS_DIR)")
    return _Codes(catalog, tuple(manifests), None)


_PALINDROMIC = [["A", "T"], ["C", "G"]]


def _read_by_pair(work: pl.DataFrame, catalog: pl.DataFrame, counts: Counter[str]) -> pl.DataFrame:
    """A mixed facet's withheld rows read by their probe's two alleles on the plus strand (RM281).

    A row is matched with a manifest SNV probe by an rsid among its `ids`, on the same contig, and read
    only when the probe's pair holds the row's reference base (so the pair is on the facet's own plus
    strand) and is not its own complement (`A/T`, `C/G`, which read alike on both strands). Called
    alleles all in the pair: the row is on the forward strand, and its alleles are derived as written.
    All in the pair's complement: it is on the minus strand, and `gt_alleles` is complemented before
    they are derived; `source_gt` keeps what the export wrote. Anything else stays withheld, counted.
    Returns `work` with `_withheld` cleared on the rows read.
    """
    rows = work.filter(pl.col("_withheld").is_not_null()).select("_row", "chrom", "ids", "gt_alleles", "ref")
    if rows.is_empty():
        return work
    pairs = catalog.filter(pl.col("chrom") == rows["chrom"][0]).select("rsid", "plus")
    matched = (
        rows.select("_row", pl.col("ids").list.unique().alias("rsid")).explode("rsid", empty_as_null=False)
        .join(pairs, on="rsid", how="inner")
        .group_by("_row").agg(pl.col("plus").unique())
        .filter(pl.col("plus").list.len() == 1)
        .with_columns(pl.col("plus").list.first())
    )  # fmt: skip
    called = pl.col("gt_alleles").list.drop_nulls()
    complement = pl.col("plus").list.eval(pl.element().replace_strict(COMPLEMENT, default=None))
    decision = (
        pl.when(pl.col("plus").is_null()).then(pl.lit("no_probe"))
        .when(pl.col("plus").is_in(_PALINDROMIC)).then(pl.lit("palindromic"))
        .when(~pl.col("plus").list.contains(pl.col("ref"))).then(pl.lit("ref_outside"))
        .when(called.list.set_difference(pl.col("plus")).list.len() == 0).then(pl.lit("forward"))
        .when(called.list.set_difference(complement).list.len() == 0).then(pl.lit("minus"))
        .otherwise(pl.lit("neither"))
    )  # fmt: skip
    decided = rows.join(matched, on="_row", how="left").select("_row", decision.alias("_pair"))
    for name, count in decided.group_by("_pair").len().iter_rows():
        counts[f"pair_{name}"] += count
    flipped = pl.col("gt_alleles").list.eval(pl.element().replace_strict(COMPLEMENT, default=None))
    work = work.join(decided, on="_row", how="left", maintain_order="left")
    return work.with_columns(
        pl.when(pl.col("_pair") == "minus").then(flipped).otherwise(pl.col("gt_alleles")).alias("gt_alleles"),
        pl.when(pl.col("_pair").is_in(["forward", "minus"])).then(None).otherwise(pl.col("_withheld")).alias("_withheld"),
    ).drop("_pair")


def _pair_warnings(counts: Counter[str], pairs: _Codes | None) -> list[str]:
    """What the manifests' allele pairs read of a mixed facet's withheld rows."""
    if pairs is None:
        return []
    asked = sum(counts[f"pair_{name}"] for name in ("forward", "minus", "palindromic", "ref_outside", "neither", "no_probe"))
    if not pairs.manifests:
        return [f"strand of withheld rows not read by probe: {pairs.absent}"]
    if not asked:
        return []
    used = ", ".join(manifest.descriptor for manifest in pairs.manifests)
    return [
        f"of {asked} array row(s) withheld for the facet's mixed strand, {counts['pair_forward'] + counts['pair_minus']} read by their "
        f"probe's two alleles on the plus strand ({used}, matched by rsid): {counts['pair_forward']} on the forward strand, "
        f"{counts['pair_minus']} on the minus strand and complemented; still withheld: {counts['pair_palindromic']} on A/T or C/G "
        f"probes, {counts['pair_ref_outside']} whose reference base is not one of the probe's alleles, {counts['pair_neither']} "
        f"carrying neither the pair nor its complement, {counts['pair_no_probe']} with no probe (or an rsid with two pairs)"
    ]


def _resolve_codes(work: pl.DataFrame, store: ReferenceStore, contig: Contig, counts: Counter[str], codes: _Codes | None) -> pl.DataFrame:
    """Array rows whose genotype is `I`/`D` codes, resolved against a manifest's indel probe (RM281).

    A row is matched with a probe by an rsid among its `ids`, on the same contig. The probe's site is
    located in the reference by its flanks (`locate_indel`) and spelled with REF starting at the row's
    `pos` (`place_indel`), REF read from the reference: `I` is the longer of the site's two alleles and
    `D` the shorter. The row gets `ref`, `alts` (the called alleles other than REF, as for any array
    row), `gt`, `end_pos`, and `gt_alleles` spelled as the sequences; `source_gt` keeps the codes. A row
    whose rsids name no probe or more than one site, or whose site is not located or cannot be spelled
    at its position, is left as loaded, each reason counted. A facet whose chip's I/D orientation is
    not verified (`probes.ID_CALIBRATION`) has every code row left as loaded, counted
    `id_orientation_unverified`. Returns `work` with the `_split_*` slots filled for the rows resolved
    here, and their `gt_alleles` replaced.
    """
    is_code = pl.col("source_ref").is_null() & pl.col("gt_alleles").list.eval(pl.element().is_in(CODES)).list.any().fill_null(False)
    only_codes = pl.col("gt_alleles").list.eval(pl.element().is_null() | pl.element().is_in(CODES)).list.all()
    rows = work.filter(is_code).select("_row", "chrom", "pos", "ids", "gt_alleles", only_codes.alias("_only"))
    found: dict[str, list] = {"_row": [], **{f"_code_{name}": [] for name in _CODE_NAMES}}
    if rows.height:
        counts["code_rows"] += rows.height
        counts["code_mixed"] += rows.filter(~pl.col("_only")).height
        rows = rows.filter(pl.col("_only"))
        if codes is not None and codes.unverified is not None:
            counts["code_id_orientation_unverified"] += rows.height
        elif codes is None or codes.catalog is None:
            counts["code_no_manifest"] += rows.height
        elif rows.height:
            probes = codes.catalog.filter(pl.col("chrom") == rows["chrom"][0])
            _resolve_code_rows(rows, probes, store, contig, counts, found)
    resolved = pl.DataFrame(found, schema={"_row": pl.UInt32, **{f"_code_{name}": CALLS[name] for name in _CODE_NAMES}})
    work = work.join(resolved, on="_row", how="left", maintain_order="left")
    return work.with_columns(
        *[pl.coalesce(f"_split_{name}", f"_code_{name}").alias(f"_split_{name}") for name in ("ref", "alts", "gt", "end_pos")],
        pl.when(pl.col("_code_gt_alleles").is_not_null()).then(pl.col("_code_gt_alleles")).otherwise(pl.col("gt_alleles")).alias("gt_alleles"),
    ).drop([f"_code_{name}" for name in _CODE_NAMES])


def _resolve_code_rows(
    rows: pl.DataFrame, probes: pl.DataFrame, store: ReferenceStore, contig: Contig, counts: Counter[str], found: dict[str, list]
) -> None:
    """`_resolve_codes` for one contig's code rows against that contig's probes; appends to `found`."""
    matched = (
        rows.select("_row", pl.col("ids").list.unique().alias("rsid")).explode("rsid", empty_as_null=False)
        .join(probes, on="rsid", how="inner")
        .group_by("_row", maintain_order=True)
        .agg(pl.struct("map_info", "left_flank", "indel_seq", "right_flank", "source_strand").unique(maintain_order=True).alias("probes"))
    )  # fmt: skip
    by_row = dict(zip(matched["_row"].to_list(), matched["probes"].to_list()))
    located: dict[tuple, IndelSite | str] = {}
    wanted = {tuple(probe.values()) for probe_list in by_row.values() for probe in probe_list}
    spans = {key: window_for(key[0], key[1], key[2], key[3], contig.length) for key in wanted}
    keys = sorted(spans)
    windows = dict(zip(keys, store.fetch(contig, [spans[key] for key in keys]))) if keys else {}
    for key in keys:
        located[key] = locate_indel(windows[key], spans[key][0], key[1], key[2], key[3], key[4])
    for row, pos, alleles in rows.select("_row", "pos", "gt_alleles").iter_rows():
        probe_list = by_row.get(row)
        if not probe_list:
            counts["code_no_probe"] += 1
            continue
        sites = {(key, located[key]) for key in (tuple(probe.values()) for probe in probe_list)}
        good = {site for _, site in sites if isinstance(site, IndelSite)}
        if len(good) > 1:
            counts["code_ambiguous"] += 1
            continue
        if not good:
            counts[f"code_{sorted(str(site) for _, site in sites)[0]}"] += 1
            continue
        site = good.pop()
        key = min(key for key, candidate in sites if candidate == site)  # deterministic among probes of one site
        placed = place_indel(windows[key], spans[key][0], site, pos)
        if isinstance(placed, str):
            counts[f"code_{placed}"] += 1
            continue
        ref, alt = placed
        longer, shorter = (ref, alt) if site.reference_carries else (alt, ref)
        spelled = [None if code is None else (longer if code == "I" else shorter) for code in alleles]
        called = {allele for allele in spelled if allele is not None}
        counts["code_resolved"] += 1
        for name, value in zip(
            found,
            (row, ref, sorted(called - {ref}), [None if a is None else int(a != ref) for a in spelled], pos + len(ref) - 1, spelled),
        ):
            found[name].append(value)


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
    codes = rest["gt_alleles"].list.eval(pl.element().is_in(CODES)).list.any().fill_null(False)
    no_call = rest["gt_alleles"].is_not_null() & (rest["gt_alleles"].list.drop_nulls().list.len() == 0)
    # Rows with an I/D code are counted by `_resolve_codes`, by outcome.
    counts["array_no_call"] += int((no_call & ~codes).sum())
    counts["array_unreadable"] += int(rest["gt_alleles"].is_null().sum())


def _derive_array_alleles(work: pl.DataFrame, counts: Counter[str], withhold: bool = False, pairs: "_Codes | None" = None) -> pl.DataFrame:
    """For array rows whose base was read: `alts` the observed non-REF alleles, sorted; `gt` into `[ref] + alts`.

    With `withhold` (a mixed-strand facet, F13) a row whose genotype is not strand-safe keeps `ref` and
    gets neither, each reason counted, unless a manifest's allele pair for its probe decides its strand
    (`pairs`, `_read_by_pair`).
    """
    target = pl.col("_array") & pl.col("ref").is_not_null()
    if withhold:
        reason = withheld_reason(pl.col("gt_alleles"), pl.col("ref"))
        work = work.with_columns(pl.when(target).then(reason).otherwise(None).alias("_withheld"))
        if pairs is not None and pairs.catalog is not None:
            work = _read_by_pair(work, pairs.catalog, counts)
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
    # Vectorised over haplotype slots and alt positions: a slot's index is 0 when its allele is REF and
    # j + 1 when it is the j-th sorted alt (a single base, or a spanned row's sequence).
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
    counts["array_ref_not_acgt"] += derived.filter(pl.col("ref").str.contains("[^ACGT]")).height
    counts["array_no_ref_allele"] += derived.filter(~pl.col("gt").list.contains(0)).height
    return work


_CODE_REASONS = {
    "code_mixed": "an I/D code beside a base letter",
    "code_no_probe": "no indel probe of the manifests carries one of their rsids on their contig",
    "code_ambiguous": "their rsids name more than one indel site",
    "code_not_found": "their probe's flanks are not found in the reference near its position",
    "code_several": "their probe's flanks are found at more than one place",
    "code_strand": "their probe's flanks are found only on the strand its manifest does not give its source sequence on",
    "code_before_pos": "the site lies before their position, and no shift along a repeat brings it there",
    "code_far": "the site lies too far after their position",
    "code_window": "the site reaches past the reference read for it",
}


def _code_warnings(counts: Counter[str], codes: _Codes | None) -> list[str]:
    """One line for the I/D rows resolved, one per reason the others were left unread."""
    if not counts["code_rows"]:
        return []
    lines = []
    if codes is not None and codes.manifests:
        used = ", ".join(f"{manifest.descriptor} (sha256 {manifest.sha256[:12]}…)" for manifest in codes.manifests)
        lines.append(
            f"{counts['code_resolved']} of {counts['code_rows']} array row(s) with an I/D indel code resolved against the indel "
            f"probes of {used}: matched by rsid, the site located by the probe's flanks in the reference, REF read from it"
        )
    if counts["code_id_orientation_unverified"]:
        unverified = codes.unverified if codes is not None else None
        lines.append(
            f"{counts['code_id_orientation_unverified']} array row(s) with an I/D indel code left unread "
            f"(id_orientation_unverified): the chip's I/D orientation against the manifests' probes is not verified; {unverified}"
        )
    if counts["code_no_manifest"]:
        absent = codes.absent if codes is not None and codes.absent else "no manifest was given"
        lines.append(f"{counts['code_no_manifest']} array row(s) with an I/D indel code left unread: {absent}")
    for key, reason in _CODE_REASONS.items():
        if counts[key]:
            lines.append(f"{counts[key]} array row(s) with an I/D indel code left unread: {reason}")
    return lines


def _row_warnings(counts: Counter[str], digest: str, codes: _Codes | None = None, pairs: _Codes | None = None) -> list[str]:
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
    if counts["point_anchored"]:
        at_start = counts["point_at_start"]
        start = f" ({at_start} at a contig's start, anchored on the base after)" if at_start else ""
        lines.append(
            f"{counts['point_anchored']} insertion point(s) given their anchor base from the reference{start}: REF is that "
            "base and every allele carries it, as VCF writes an insertion; the source wrote no REF base to compare it with"
        )
    if counts["point_no_call"]:
        lines.append(f"{counts['point_no_call']} insertion point(s) with no haplotype called: REF left unread (no allele to anchor)")
    if counts["point_no_insertion"]:
        lines.append(
            f"{counts['point_no_insertion']} insertion point(s) where every called haplotype is the empty allele: REF left "
            "unread. They say nothing is inserted there, which is not a reference call at the base before; anchored, "
            "they would read as one"
        )
    if counts["point_outside"]:
        lines.append(f"{counts['point_outside']} insertion point(s) past the end of their contig: REF left unread")
    if counts["array_rows"]:
        lines.append(
            f"REF read from collection {digest} for {counts['array_read']} of {counts['array_rows']} array row(s) "
            "on readable contigs; ALT derived from the genotype against it"
        )
    reasons = {
        "array_no_call": "with no allele called (the probe's REF span is unknown)",
        "array_unreadable": "whose genotype could not be split into alleles",
        "array_outside": "past the end of their contig, or at position 0",
    }
    for key, reason in reasons.items():
        if counts[key]:
            lines.append(f"{counts[key]} array row(s) {reason}: REF left unread")
    lines += _code_warnings(counts, codes)
    lines += _pair_warnings(counts, pairs)
    split =counts["concat_resolved"] + counts["concat_homozygous"] + counts["concat_repeat"]
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

