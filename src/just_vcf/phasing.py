"""VRS CisPhasedBlocks, minted on demand from the phase columns (RM284, charter: phasing is VRS
CisPhasedBlock).

Phase sets stay columns of `calls` (`phased`, `phase_sets`, `phase_set_kinds`, docs/SCHEMA.md
Decision 2). `cis_phased_blocks(loaded)` reads them and returns one row per haplotype of each phase
set: the VRS 2.x `CisPhasedBlock` of the alleles that haplotype carries, with its computed id. Nothing
is stored back into `calls`.

What a block is:

- **A PS set** is keyed `(facet_id, chrom, PS)`, the implicit set a null PS. A PS does not say which
  molecule is which, the GT order does: slot `h` of every record in the set is one haplotype, so a
  PS set gives one block per slot, `haplotype = h`.
- **A PSL set** is keyed `(facet_id, name)`, `chrom` null: VCF 4.5 ties a PSL name to one haplotype
  and lets it span chromosomes, so the set itself is one block, `haplotype` null, whatever slot each
  allele sits in.
- **The members** are the ids in `vrs_allele_ids` of the alleles the haplotype carries: by default
  the non-reference ones only. A `0` on a record split one row per ALT, or under an overlapping
  deletion, means "not this ALT", not the reference base, so a REF member would state something the
  file did not. `reference_alleles=True` takes REF slots too, for a caller who knows the facet holds
  no such rows. `*` is never a member: it stands for a deletion that is a member at its own record.
  An allele carried twice (a duplicated record) is one member: a block is a set of alleles in cis.
- **Haploid calls** have no separator, so VCF reads them phased, and they sit in their chromosome's
  implicit set even in a file nobody phased. On a chromosome with one copy they are in cis by ploidy.

When a block is withheld, and why (`withheld`; `cpb_id` and `members` are then null): a block that
leaves out an allele the haplotype carries would carry an id that looks like the haplotype's and is
not (charter 4), so any slot of the haplotype that cannot be stated withholds the whole block. The
first such slot in row order gives the reason: its phase is unknown (`phase_unknown`, a VCF 4.4+ slot
0 polars-bio's rewrite hides, P4, when its allele would be a member), it is uncalled (`uncalled`),
its allele has no index into the site's alleles (`not_indexed`), or the allele has no VRS id (the
reason `vrs_withheld` gives, `not_run` when `mint_vrs_ids` has not run). A haplotype with fewer than two members is
`too_few_members`: the spec requires two, and one allele is an Allele, not a block.

The id is the VRS 2.x computed identifier: `ga4gh:CPB.` + sha512t24u of
`{"members":[…],"type":"CisPhasedBlock"}`, the members being the Allele digests (the id without
`ga4gh:VA.`), sorted. That is what the spec's digest serialization makes of inline Allele members, and
what vrs-python 2.3 computes for them. vrs-python digests an IRI member as written, prefix and all,
so a block re-identified from `members` as IRIs gets another id: resolve them to Alleles first.
`sequenceReference` is not an inherent property and is not emitted.
"""

import polars as pl
from gtars.refget import sha512t24u_digest

from just_vcf.loading import Loaded
from just_vcf.schema import PhaseSetKind, enum_dtype
from just_vcf.vrs import WITHHELD_REASONS, vrs_withheld

CPB_WITHHELD_REASONS: frozenset[str] = frozenset({"phase_unknown", "uncalled", "not_indexed", "too_few_members"}) | WITHHELD_REASONS
"""Every reason `cis_phased_blocks` can give for a withheld block: its own, and the member allele's
`vrs_withheld` reason. Open, like `vrs.WITHHELD_REASONS`."""

CIS_PHASED_BLOCKS = pl.Schema(
    {
        "facet_id": pl.String,
        "chrom": pl.String,
        "phase_set_kind": enum_dtype(PhaseSetKind),
        "phase_set": pl.String,
        "haplotype": pl.UInt32,
        "sites": pl.UInt32,
        "members": pl.List(pl.String),
        "cpb_id": pl.String,
        "withheld": pl.String,
    }
)
"""The frame `cis_phased_blocks` returns (docs/SCHEMA.md, "CisPhasedBlock export")."""

_VA = "ga4gh:VA."
_CPB = "ga4gh:CPB."


def cis_phased_blocks(loaded: Loaded, *, reference_alleles: bool = False) -> pl.LazyFrame:
    """One row per haplotype of each phase set in `loaded.calls`, with its VRS CisPhasedBlock.

    `members` are the haplotype's allele ids, sorted, and `cpb_id` the block's id; both are null and
    `withheld` says why when the block cannot be stated whole (module docstring). `sites` counts the
    records whose slot is in the set on this haplotype. Rows come in the order their set first
    appears in `calls`, ties by kind, set and haplotype. Unphased rows are in no set; a facet with
    none gives an empty frame.

    Ids come from `vrs_allele_ids`, so run `mint_vrs_ids` first; until then every block that needs an
    allele's id is withheld `not_run`.
    """
    phased_rows = loaded.calls.with_row_index("_order").filter(pl.col("phased").list.eval(pl.element().fill_null(False)).list.any().fill_null(False))
    rows = (
        vrs_withheld(phased_rows)
        .select("_order", "facet_id", "chrom", "gt", "gt_alleles", "phased", "phase_sets", "phase_set_kinds", "vrs_allele_ids", "vrs_withheld")
    )
    is_ps = pl.col("phase_set_kinds").list.eval(pl.element() == PhaseSetKind.PS.value).list.any().fill_null(False)
    # A PS is one value per record, on each of its phased slots; null for the implicit set.
    row_ps = pl.col("phase_sets").list.drop_nulls().list.first()
    slots = (
        rows.with_columns(is_ps.alias("_is_ps"), row_ps.alias("_row_ps"), pl.int_ranges(0, pl.col("phased").list.len(), dtype=pl.UInt32).alias("_h"))
        .explode("_h", empty_as_null=False)
        .with_columns(
            pl.col("phased").list.get(pl.col("_h")).alias("_phased"),
            pl.col("phase_sets").list.get(pl.col("_h"), null_on_oob=True).alias("_set"),
            pl.col("phase_set_kinds").list.get(pl.col("_h"), null_on_oob=True).cast(pl.String).alias("_kind"),
            pl.col("gt_alleles").list.get(pl.col("_h"), null_on_oob=True).alias("_allele"),
            pl.col("gt").list.get(pl.col("_h"), null_on_oob=True).alias("_index"),
        )
        .with_columns(
            pl.col("vrs_allele_ids").list.get(pl.col("_index").cast(pl.Int64), null_on_oob=True).alias("_id"),
            pl.col("vrs_withheld").list.get(pl.col("_index").cast(pl.Int64), null_on_oob=True).alias("_why"),
        )
    )
    in_set = pl.col("_phased").fill_null(False)
    skipped = (pl.col("_allele") == "*") | ((pl.col("_index") == 0) & pl.lit(not reference_alleles))
    member = in_set & pl.col("_id").is_not_null() & ~skipped.fill_null(False)
    # A slot that would not be a member whatever its phase (REF by default, `*`) withholds nothing.
    problem = (
        pl.when(pl.col("_phased") == False)  # noqa: E712 - null is a different answer
        .then(None)
        .when(pl.col("_allele").is_not_null() & pl.col("_index").is_not_null() & skipped.fill_null(False))
        .then(None)
        .when(pl.col("_phased").is_null())
        .then(pl.lit("phase_unknown"))
        .when(pl.col("_allele").is_null())
        .then(pl.lit("uncalled"))
        .when(pl.col("_index").is_null())
        .then(pl.lit("not_indexed"))
        .when(pl.col("_id").is_null())
        .then(pl.coalesce(pl.col("_why"), pl.lit("not_run")))
        .otherwise(None)
    )
    slots = slots.with_columns(member.alias("_member"), problem.alias("_problem"), in_set.alias("_in_set"))

    # A PS set: every slot of its records but the unphased ones, one block per slot.
    ps = slots.filter(pl.col("_is_ps") & (pl.col("_phased").is_null() | pl.col("_phased"))).select(
        "_order", "facet_id", "chrom", pl.lit(PhaseSetKind.PS.value).alias("_kind"), pl.col("_row_ps").alias("phase_set"), pl.col("_h").alias("haplotype"),
        "_member", "_problem", "_in_set", "_id",
    )
    # A PSL set: the phased slots that name it, wherever they sit, one block per set.
    psl = slots.filter(in_set & (pl.col("_kind") == PhaseSetKind.PSL.value)).select(
        "_order", "facet_id", pl.lit(None, dtype=pl.String).alias("chrom"), "_kind", pl.col("_set").alias("phase_set"), pl.lit(None, dtype=pl.UInt32).alias("haplotype"),
        "_member", "_problem", "_in_set", "_id",
    )
    key = ["facet_id", "chrom", "_kind", "phase_set", "haplotype"]
    grouped = pl.concat([ps, psl]).group_by(key).agg(
        pl.col("_order").min().alias("_first"),
        pl.col("_in_set").sum().cast(pl.UInt32).alias("sites"),
        pl.col("_problem").sort_by("_order").drop_nulls().first().alias("_problem"),
        pl.col("_id").filter(pl.col("_member")).unique().sort().alias("_members"),
    )
    withheld = pl.coalesce(pl.col("_problem"), pl.when(pl.col("_members").list.len() < 2).then(pl.lit("too_few_members")))
    canonical = pl.concat_str(
        pl.lit('{"members":["'), pl.col("members").list.eval(pl.element().str.strip_prefix(_VA)).list.join('","'), pl.lit('"],"type":"CisPhasedBlock"}')
    )
    return (
        grouped.with_columns(withheld.alias("withheld"))
        .with_columns(pl.when(pl.col("withheld").is_null()).then(pl.col("_members")).alias("members"))
        .with_columns((pl.lit(_CPB) + canonical.map_batches(_digest, return_dtype=pl.String)).alias("cpb_id"))
        .sort("_first", "_kind", "phase_set", "haplotype", nulls_last=True)
        .select(
            "facet_id",
            "chrom",
            pl.col("_kind").cast(enum_dtype(PhaseSetKind)).alias("phase_set_kind"),
            "phase_set",
            "haplotype",
            "sites",
            "members",
            "cpb_id",
            "withheld",
        )
    )


def _digest(texts: pl.Series) -> pl.Series:
    """gtars' sha512t24u of each string; null stays null."""
    return pl.Series([None if text is None else sha512t24u_digest(text) for text in texts.to_list()], dtype=pl.String)
