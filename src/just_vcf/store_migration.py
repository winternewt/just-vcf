"""The steps that migrate a durable store's frames from one schema version to the next (RM276).

`store.migrate_store` rewrites a store written under an older `SCHEMA_VERSION` to this one by applying
the steps here in order, one per version bump from the store's version to `SCHEMA_VERSION`. Each step
is a pure function per frame, `LazyFrame -> LazyFrame`, that turns the frame as that version's loaders
wrote it into the frame as the next version's do; a frame the bump did not change passes through. The
shapes come from docs/SCHEMA.md, "Changes since sign-off"; "Store migration" there has the table.

A step never invents a value. A column the bump added is filled null. Where that null, or a value an
older loader wrote, reads differently from what a load of the same source gives today (a VCF header's
method chain was not read, so `method_steps` null there is "not read", not "not a VCF"), the step
names the facets it applies to and why, as a `FilledUnknown` in the migration's report and the store's
manifest; `Store.loaded()` repeats it as a warning. Loading the source again is the only way to get
those values.

Each step states its dtypes as they were at its own version, not as `schema` has them now, so a later
bump that grows an enum again does not change what an earlier step writes. The suite checks the chain
is contiguous from the store's first version to `SCHEMA_VERSION` and ends in exactly `schema`'s
frames, so a bump without a step fails it.
"""

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import datetime

import polars as pl
from pydantic import BaseModel, ConfigDict

from just_vcf.loaders.cgi import FORMATS as CGI_FORMATS
from just_vcf.loaders.cgi import MASTER_VAR
from just_vcf.loaders.vcf import FORMATS as VCF_FORMATS
from just_vcf.schema import SCHEMA_VERSION

STORE_FIRST_VERSION = "0.0.14"
"""The schema version the store shipped under (merge 6830fb7, 2026-10-23): no store is older."""

FrameStep = Callable[[pl.LazyFrame], pl.LazyFrame]


class FilledUnknown(BaseModel):
    """Values a migration could not carry honestly: what they are, on which facets, and why."""

    model_config = ConfigDict(frozen=True)

    version: str
    """The schema version whose step this is (its target)."""
    frame: str
    columns: tuple[str, ...]
    facet_ids: tuple[str, ...]
    """The facets it applies to, in `facet_id` order. Never empty in a report."""
    reason: str


class MigrationApplied(BaseModel):
    """One step a migration applied."""

    model_config = ConfigDict(frozen=True)

    from_version: str
    to_version: str
    summary: str
    added: tuple[str, ...] = ()
    """`frame.column` the step added, filled null."""
    retyped: tuple[str, ...] = ()
    """`frame.column` the step cast to a wider dtype, values unchanged."""


class StoreMigration(BaseModel):
    """What `migrate_store` did: the steps applied, the facets rewritten, what was filled unknown.

    From a store already at this version nothing is applied and `migrated` is null.
    """

    model_config = ConfigDict(frozen=True)

    from_version: str
    to_version: str
    migrated: datetime | None = None
    """When the migrated store was committed, ISO-8601 UTC."""
    steps: tuple[MigrationApplied, ...] = ()
    facets: tuple[str, ...] = ()
    """Every facet rewritten, in `facet_id` order."""
    filled_unknown: tuple[FilledUnknown, ...] = ()


@dataclass(frozen=True)
class Caveat:
    """A null a step fills, or an older value it keeps, that reads differently from a load today, on
    the facets `applies` selects (an expression over that version's `facets` frame) and, when `rows` is
    given, only those of them holding a row it selects (an expression over that version's `frame`)."""

    frame: str
    columns: tuple[str, ...]
    applies: pl.Expr
    reason: str
    rows: pl.Expr | None = None


@dataclass(frozen=True)
class MigrationStep:
    """One schema version bump, as a pure function per frame (a frame not named passes through)."""

    from_version: str
    to_version: str
    summary: str
    frames: Mapping[str, FrameStep] = field(default_factory=dict)
    added: tuple[str, ...] = ()
    retyped: tuple[str, ...] = ()
    caveats: tuple[Caveat, ...] = ()

    def apply(self, name: str, frame: pl.LazyFrame) -> pl.LazyFrame:
        step = self.frames.get(name)
        return frame if step is None else step(frame)

    def applied(self) -> MigrationApplied:
        return MigrationApplied(
            from_version=self.from_version,
            to_version=self.to_version,
            summary=self.summary,
            added=self.added,
            retyped=self.retyped,
        )


def _add_after(after: str, columns: Mapping[str, pl.DataType]) -> FrameStep:
    """A step adding `columns`, null, right after the column `after`, as the bump placed them."""

    def step(frame: pl.LazyFrame) -> pl.LazyFrame:
        names = frame.collect_schema().names()
        present = [name for name in columns if name in names]
        if present:
            raise ValueError(f"column(s) {', '.join(present)} already present: not the shape this step migrates")
        at = names.index(after) + 1
        return frame.select(
            *[pl.col(name) for name in names[:at]],
            *[pl.lit(None, dtype=dtype).alias(name) for name, dtype in columns.items()],
            *[pl.col(name) for name in names[at:]],
        )

    return step


def _cast(column: str, dtype: pl.DataType) -> FrameStep:
    """A step casting `column` strictly to a wider `dtype`: a value that does not fit raises."""

    def step(frame: pl.LazyFrame) -> pl.LazyFrame:
        return frame.with_columns(pl.col(column).cast(dtype, strict=True))

    return step


def _dropped(*fields: str) -> pl.Expr:
    """A facet whose `dropped_fields` lists any of `fields`."""
    return pl.col("dropped_fields").list.eval(pl.element().is_in(list(fields))).list.any().fill_null(False)


_VCF = pl.col("format_id").is_in(sorted(VCF_FORMATS))
_CGI = pl.col("format_id").is_in(sorted(CGI_FORMATS))

_METHOD_STEP_0016 = pl.Struct(
    {"tool": pl.String, "version": pl.String, "role": pl.String, "options": pl.List(pl.String), "line": pl.String}
)
_PHASE_SET_KINDS_0017 = pl.List(pl.Enum(["ps", "psl", "cgi_haplink"]))
_GL_SOURCE_0018 = pl.Enum(["gl", "pl", "gl_over_pl"])
_REF_CHECK_0022 = pl.Enum(["not_run", "match", "mismatch", "unstated"])

_RELOAD = "load the source again to read them"
_READ = pl.col("ref_seqcol_digest").is_not_null()

MIGRATIONS: tuple[MigrationStep, ...] = (
    MigrationStep(
        "0.0.14",
        "0.0.15",
        "PSL phase sets: calls.phase_set_ordinals (PSO) and phase_set_qualities (PSQ)",
        frames={"calls": _add_after("phase_set_kinds", {"phase_set_ordinals": pl.List(pl.UInt32), "phase_set_qualities": pl.List(pl.UInt16)})},
        added=("calls.phase_set_ordinals", "calls.phase_set_qualities"),
        caveats=(
            Caveat(
                "calls",
                ("phased", "phase_sets", "phase_set_kinds", "phase_set_ordinals", "phase_set_qualities"),
                _dropped("FORMAT/PSL", "FORMAT/PSO", "FORMAT/PSQ"),
                "the load left out FORMAT PSL, PSO and PSQ (in dropped_fields), which 0.0.15 reads as phase sets: "
                f"these rows carry no PSL set and their ordinals and qualities are null for not read; {_RELOAD}",
            ),
        ),
    ),
    MigrationStep(
        "0.0.15",
        "0.0.16",
        "layer 2's method chain: facets.reference, aligner, aligner_version, method_steps",
        frames={
            "facets": _add_after(
                "phasing_method",
                {"reference": pl.String, "aligner": pl.String, "aligner_version": pl.String, "method_steps": pl.List(_METHOD_STEP_0016)},
            )
        },
        added=("facets.reference", "facets.aligner", "facets.aligner_version", "facets.method_steps"),
        caveats=(
            Caveat(
                "facets",
                ("reference", "aligner", "aligner_version", "method_steps", "phasing_method"),
                _VCF,
                "a VCF header's method chain was not read before 0.0.16: method_steps null here means not read, not "
                "'not a VCF', and reference, aligner and aligner_version null mean not read, not 'not recorded'; "
                f"phasing_method keeps the value the earlier rule gave; {_RELOAD}",
            ),
        ),
    ),
    MigrationStep(
        "0.0.16",
        "0.0.17",
        "CGI hapLink: PhaseSetKind gains cgi_haplink, widening calls.phase_set_kinds",
        frames={"calls": _cast("phase_set_kinds", _PHASE_SET_KINDS_0017)},
        retyped=("calls.phase_set_kinds",),
        caveats=(
            Caveat(
                "calls",
                ("phased", "phase_sets", "phase_set_kinds"),
                _CGI,
                "CGI's hapLink was not read before 0.0.17: these rows' phase is null for not read, not 'unphased'; "
                f"{_RELOAD}",
            ),
            Caveat(
                "calls",
                ("ad", "source_alleles"),
                pl.col("format_id") == MASTER_VAR,
                f"masterVarBeta's per-allele read counts were not read before 0.0.17; {_RELOAD}",
            ),
        ),
    ),
    MigrationStep(
        "0.0.17",
        "0.0.18",
        "the PL/GL source: facets.gl_source",
        frames={"facets": _add_after("split_marker", {"gl_source": _GL_SOURCE_0018})},
        added=("facets.gl_source",),
        caveats=(
            Caveat(
                "facets",
                ("gl_source",),
                _VCF,
                "which FORMAT key gl was read from was not recorded before 0.0.18: gl_source null here means not "
                f"recorded, not 'no likelihoods'; {_RELOAD}",
            ),
        ),
    ),
    MigrationStep(
        "0.0.18",
        "0.0.19",
        "CGI's allele scores: calls.hq",
        frames={"calls": _add_after("gq", {"hq": pl.List(pl.UInt16)})},
        added=("calls.hq",),
        caveats=(
            Caveat(
                "calls",
                ("hq",),
                _CGI | _dropped("FORMAT/HQ"),
                "CGI's varScoreEAF, and a VCF's FORMAT HQ (in dropped_fields), were not read into hq before 0.0.19 "
                f"and its first update: hq null here means not read, so these rows grade gt_only; {_RELOAD}",
            ),
        ),
    ),
    MigrationStep(
        "0.0.19",
        "0.0.20",
        "the rule that read REF: facets.resolve_version",
        frames={"facets": _add_after("ref_seqcol_digest", {"resolve_version": pl.UInt32})},
        added=("facets.resolve_version",),
        caveats=(
            Caveat(
                "facets",
                ("resolve_version",),
                _READ,
                "which RESOLVE_VERSION read these facets' REF was not recorded before 0.0.20: resolve_version null "
                "here means not recorded, not 'REF not read', so the store refuses the same facet resolved today as "
                f"another state of it; {_RELOAD} and resolve them into a new store",
            ),
        ),
    ),
    MigrationStep(
        "0.0.20",
        "0.0.21",
        "the manifests that read REF: facets.resolve_inputs",
        frames={"facets": _add_after("resolve_version", {"resolve_inputs": pl.String})},
        added=("facets.resolve_inputs",),
        caveats=(
            Caveat(
                "facets",
                ("resolve_inputs",),
                _READ & (pl.col("family") == "genotype_table"),
                "which chip manifests read these array facets' I/D codes and mixed-strand rows was not recorded before "
                "0.0.21: resolve_inputs null here means not recorded, not 'no manifest', so the store takes the same "
                f"facet resolved today with manifests for another state of it; {_RELOAD} and resolve them into a new store",
            ),
            Caveat(
                "facets",
                ("resolve_inputs",),
                pl.col("relation") == "consensus",
                "which REF states these consensus facets' parents were in was not recorded before 0.0.21: resolve_inputs "
                "null here means not recorded, so the store takes the same consensus built today for another state of "
                "it; build it again into a new store",
            ),
        ),
    ),
    MigrationStep(
        "0.0.21",
        "0.0.22",
        "a source REF that states no base: RefCheck gains unstated, widening calls.ref_check",
        frames={"calls": _cast("ref_check", _REF_CHECK_0022)},
        retyped=("calls.ref_check",),
        caveats=(
            Caveat(
                "calls",
                ("ref_check",),
                _READ,
                "a source REF of only N was compared as bases before 0.0.22 (RESOLVE_VERSION 12 and earlier): these rows keep "
                "ref_check mismatch where today's rule reads unstated (the source stated no base), so liftover leaves "
                "them behind and the store refuses the same facet resolved today as another state of it; "
                f"{_RELOAD} and resolve them into a new store",
                rows=pl.col("source_ref").str.contains(r"^[Nn]+$").fill_null(False) & (pl.col("ref_check") == "mismatch"),
            ),
        ),
    ),
    MigrationStep(
        "0.0.22",
        "0.0.23",
        "per-site imputation quality: calls.imputation_r2 and imputation_empirical_r2",
        frames={"calls": _add_after("dosage", {"imputation_r2": pl.Float32, "imputation_empirical_r2": pl.Float32})},
        added=("calls.imputation_r2", "calls.imputation_empirical_r2"),
        caveats=(
            Caveat(
                "calls",
                ("imputation_r2", "imputation_empirical_r2"),
                _dropped("INFO/R2", "INFO/ER2", "INFO/DR2"),
                "the load left out INFO R2, ER2 and DR2 (in dropped_fields), which 0.0.23 reads as each site's "
                f"imputation quality: imputation_r2 and imputation_empirical_r2 null here mean not read; {_RELOAD}",
            ),
        ),
    ),
)
"""One step per schema version bump since the store shipped, in order."""


def version_key(version: str) -> tuple[int, ...] | None:
    """`0.0.19` as `(0, 0, 19)`, so `0.0.9` sorts before it; None for anything not dotted numbers."""
    parts = version.split(".")
    if not all(part.isdigit() for part in parts):
        return None
    return tuple(int(part) for part in parts)


def is_newer(version: str) -> bool:
    """A schema version later than this just-vcf's: never migrated (no downgrade)."""
    key, current = version_key(version), version_key(SCHEMA_VERSION)
    return key is not None and current is not None and key > current


def migration_chain(version: str) -> tuple[MigrationStep, ...] | None:
    """The steps from `version` to `SCHEMA_VERSION`, in order: `()` at this version, None when there is
    no path (a later version, or one no store was written under)."""
    by_source = {step.from_version: step for step in MIGRATIONS}
    chain: list[MigrationStep] = []
    at = version
    while at != SCHEMA_VERSION:
        step = by_source.get(at)
        if step is None:
            return None
        chain.append(step)
        at = step.to_version
    return tuple(chain)
