"""A facet's assembly identity, from what a VCF header carries (RM4; charter invariants 3 and 9).

Two answers, kept apart because they are different kinds of fact:

- **The label** (`GRCh37`, `GRCh38`, `T2T-CHM13`), inferred from the `##contig` lengths of the
  primary chromosomes. The lengths below are domain constants (UCSC `chrom.sizes` for hg19, hg38 and
  hs1, checked 2026-10-22); no two assemblies share a primary chromosome length except the
  mitochondrion's 16,569, so a matching length is evidence for one assembly only. A header too thin to
  tell, or one whose lengths disagree, gives no label, never the nearest one.
- **The coordinate system's digest**: the refget Sequence Collection `sorted_name_length_pairs`
  attribute digest of the header's contigs, names verbatim. It is computed exactly as the seqcol
  specification defines it, so it is a conformant value of that attribute. The collection's top-level
  digest is not: it is computed over the `names` and `sequences` attributes, and `sequences` holds the
  refget (sha512t24u) digest of each contig's bases, which no header carries (a `##contig` `md5` is a
  different digest and cannot be converted). `seqcol_digest` therefore stays null until bases are read
  (RM3), and nothing here fills it from a lookup on the label.
"""

from enum import Enum

from pydantic import BaseModel, ConfigDict
from refget.digests import sha512t24u_digest
from refget.utils import build_sorted_name_length_pairs, canonical_str

from just_vcf.schema import canonical_chrom

_NAMES = [str(n) for n in range(1, 23)] + ["X", "Y"]

ASSEMBLY_LENGTHS: dict[str, dict[str, frozenset[int]]] = {
    "GRCh37": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
                    141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                    81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
                )),
                strict=True,
            )
        ),
        # GRCh37 carries the rCRS (16,569); UCSC hg19 kept the older Yoruba chrM (16,571).
        "MT": frozenset({16569, 16571}),
    },
    "GRCh38": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
                    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
                    83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
                )),
                strict=True,
            )
        ),
        "MT": frozenset({16569}),
    },
    "T2T-CHM13": {
        **dict(
            zip(
                _NAMES,
                (frozenset({n}) for n in (
                    248387328, 242696752, 201105948, 193574945, 182045439, 172126628, 160567428, 146259331,
                    150617247, 134758134, 135127769, 133324548, 113566686, 101161492, 99753195, 96330374,
                    84276897, 80542538, 61707364, 66210255, 45090682, 51324926, 154259566, 62460029,
                )),
                strict=True,
            )
        ),
        "MT": frozenset({16569}),
    },
}
"""Primary-chromosome lengths per assembly, by canonical `chrom`. T2T-CHM13 v1.1 has no Y; its other
lengths are v2.0's, so a header of either version matches."""

MIN_MATCHING_CONTIGS = 3
"""Primary chromosomes whose length must match one assembly, and no other, before it is named.
just-prs `detect_genome_build` voted with the same floor. A header with one or two chromosomes (a
region extract, a test file) is too thin: one matching length could be a coincidence of a re-headered
file, three cannot."""


class AssemblyOutcome(str, Enum):
    """Why `infer_assembly` gave the answer it gave. Only `INFERRED` carries a label."""

    INFERRED = "inferred"
    NO_CONTIG_LENGTHS = "no_contig_lengths"
    """The header has no `##contig` line with a length on a primary chromosome."""
    TOO_FEW = "too_few"
    """Fewer than `MIN_MATCHING_CONTIGS` primary chromosomes match an assembly."""
    CONFLICTING = "conflicting"
    """Primary chromosome lengths match different assemblies, or a length matches none."""


class AssemblyInference(BaseModel):
    """The label the header's contig lengths support, or why there is none. Truthy only with a label."""

    model_config = ConfigDict(frozen=True)

    label: str | None
    outcome: AssemblyOutcome
    reason: str
    matched: dict[str, int] = {}
    """Assembly → primary chromosomes whose header length matches it (the mitochondrion excluded)."""
    unmatched: tuple[str, ...] = ()
    """Primary chromosomes whose header length matches no known assembly."""

    def __bool__(self) -> bool:
        return self.label is not None


def infer_assembly(contigs: list[dict]) -> AssemblyInference:
    """The assembly label a header's `##contig` lines support (`Header.contigs`: name, length, md5).

    Only primary chromosomes take part, by canonical name; alt contigs, decoys and unplaced scaffolds
    differ between releases of one assembly and say nothing about the label. The mitochondrion's
    length is shared by all three assemblies, so it can contradict a label (hg19's 16,571 on GRCh38)
    but never supports one.
    """
    lengths: dict[str, set[int]] = {}
    for contig in contigs:
        chrom = canonical_chrom(contig["name"])
        if chrom in ASSEMBLY_LENGTHS["GRCh38"] and contig.get("length") is not None:
            lengths.setdefault(chrom, set()).add(int(contig["length"]))
    if not lengths:
        return AssemblyInference(
            label=None, outcome=AssemblyOutcome.NO_CONTIG_LENGTHS, reason="no ##contig line gives a primary chromosome's length"
        )
    matched = {
        name: sum(1 for chrom, seen in lengths.items() if chrom != "MT" and seen <= table[chrom])
        for name, table in ASSEMBLY_LENGTHS.items()
    }
    fitting = sorted(name for name, table in ASSEMBLY_LENGTHS.items() if all(seen <= table[chrom] for chrom, seen in lengths.items()))
    unmatched = tuple(
        sorted(
            (chrom for chrom, seen in lengths.items() if not any(seen <= table[chrom] for table in ASSEMBLY_LENGTHS.values())),
            key=_NAMES_AND_MT.index,
        )
    )
    supported = {name: count for name, count in matched.items() if count}
    if len(fitting) == 1 and matched[fitting[0]] >= MIN_MATCHING_CONTIGS:
        label = fitting[0]
        return AssemblyInference(
            label=label,
            outcome=AssemblyOutcome.INFERRED,
            reason=f"{matched[label]} primary chromosome length(s) match {label} and none contradicts it",
            matched=supported,
        )
    if not fitting:
        parts = [f"{count} match {name}" for name, count in sorted(supported.items())]
        if unmatched:
            parts.append(f"{', '.join(unmatched)} match(es) no known assembly")
        return AssemblyInference(
            label=None,
            outcome=AssemblyOutcome.CONFLICTING,
            reason="contig lengths disagree: " + "; ".join(parts),
            matched=supported,
            unmatched=unmatched,
        )
    best = max(supported.values(), default=0)
    return AssemblyInference(
        label=None,
        outcome=AssemblyOutcome.TOO_FEW,
        reason=f"{best} primary chromosome length(s) match an assembly; {MIN_MATCHING_CONTIGS} are needed",
        matched=supported,
    )


_NAMES_AND_MT = [*_NAMES, "MT"]


def sorted_name_length_pairs_digest(contigs: list[dict]) -> str | None:
    """The seqcol `sorted_name_length_pairs` attribute digest of the header's contigs, or None.

    Names are the file's own, verbatim (`chr1` and `1` are different coordinate systems to seqcol).
    None when the header lists no contig, or when any contig lacks a length: a digest over the
    contigs that have one would describe a coordinate system the file never declared. A contig listed
    twice is kept twice, as the header wrote it.
    """
    if not contigs or any(contig.get("length") is None for contig in contigs):
        return None
    collection = {"names": [contig["name"] for contig in contigs], "lengths": [int(contig["length"]) for contig in contigs]}
    return sha512t24u_digest(canonical_str(build_sorted_name_length_pairs(collection)))
