"""vrs-python as the oracle for just-vcf's VRS ids (RM271): a data proxy over a gtars refget store.

Dev-only (vrs-python is a dev dependency, never imported from `src/`). The proxy reads bases from the
same `ReferenceStore` `resolve_ref` read REF from, so an id the oracle mints names the same sequence
as one just-vcf mints. It reads nothing over the network.

`oracle_id` is vrs-python's own path for a VCF-style allele (as its gnomAD translator builds it): a
`LiteralSequenceExpression` of ALT over REF's span, normalized, then identified.
"""

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrs.dataproxy import _DataProxy

from just_vcf.reference import Contig, ReferenceStore


class StoreProxy(_DataProxy):
    """A vrs-python data proxy serving the sequences of one or more collections of a store."""

    def __init__(self, store: ReferenceStore, collections: list[str] | None = None) -> None:
        self.store = store
        self.contigs: dict[str, Contig] = {}
        for digest in collections if collections is not None else store.collections():
            for contig in store.contigs(digest):
                self.contigs[f"SQ.{contig.digest}"] = contig

    def _contig(self, identifier: str) -> Contig:
        return self.contigs[identifier.removeprefix("ga4gh:")]

    def get_sequence(self, identifier: str, start: int | None = None, end: int | None = None) -> str:
        contig = self._contig(identifier)
        start = 0 if start is None else max(start, 0)
        end = contig.length if end is None else min(end, contig.length)
        if end <= start:
            return ""
        return self.store.fetch(contig, [(start, end)])[0]

    def get_metadata(self, identifier: str) -> dict:
        contig = self._contig(identifier)
        return {"length": contig.length, "aliases": [f"ga4gh:SQ.{contig.digest}"], "alphabet": "ACGTN"}


def oracle_id(proxy: _DataProxy, accession: str, pos: int, ref: str, alt: str) -> str:
    """vrs-python's id for the VCF allele `ref` → `alt` at 1-based `pos` on `accession` (`SQ.…`)."""
    start = pos - 1
    allele = models.Allele(
        location=models.SequenceLocation(
            sequenceReference=models.SequenceReference(refgetAccession=accession), start=start, end=start + len(ref)
        ),
        state=models.LiteralSequenceExpression(sequence=alt.upper()),
    )
    return ga4gh_identify(normalize(allele, proxy))
