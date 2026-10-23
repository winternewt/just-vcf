"""The `just-vcf` command line: `sniff`, `load`, `convert` and `check` one genotype file.

The CLI is a privacy boundary (RM272, CLAUDE.md §2): it prints counts and aggregates, never a row, a
sample name, a path or an archive member name. `--sample` takes a name and the output never echoes it.
Library warnings and refusal messages pass through a redaction of the strings this run knows to be
labels (the input's path and name, the sample, the archive member); a library error's own message is
never printed, since several embed the input's path. Every command takes `--json` for a machine-
readable report on stdout; refusals go to stderr.

Exit codes tell a refusal (the library declining a file it understood, `ExitCode`) from a crash (1)
and from a usage error (2, Typer's own).
"""

import json
import logging
import re
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from enum import IntEnum
from pathlib import Path
from typing import Annotated, Any

import polars as pl
import typer
from dotenv import find_dotenv, load_dotenv
from gtars.refget import RefgetStore

from just_vcf.liftover import liftover
from just_vcf.loaders.table import TableLayoutError
from just_vcf.loading import LOADABLE_FORMATS, DamagedBgzfError, Loaded, MultiSampleError, NoSampleError, UnsupportedFormatError, load
from just_vcf.normalize import normalize
from just_vcf.reference import ReferenceMismatchError, ReferenceStore, resolve_ref
from just_vcf.sniff import Sniffed, SniffOutcome, sniff
from just_vcf.vrs import mint_vrs_ids, vrs_withheld
from just_vcf.writing import UnsortedRowsError, UnwritableRowsError, write_vcf


class ExitCode(IntEnum):
    """What the process exit status means. 0 is success, 1 a crash, 2 a usage error (Typer's)."""

    OK = 0
    CRASH = 1
    """An error the library did not mean as a refusal. Its message is withheld (it can carry a path)."""
    USAGE = 2
    UNSUPPORTED = 3
    """Not recognised, not variant data, or a format no loader reads (`UnsupportedFormatError`)."""
    MULTI_SAMPLE = 4
    """Several samples and no `--sample`, or the `--sample` given is not in the file."""
    NO_SAMPLE = 5
    """A sites-only VCF: no sample column to load."""
    REFERENCE_MISMATCH = 6
    """The reference cannot serve the facet (`ReferenceMismatchError`), an empty store included."""
    DAMAGED_BGZF = 7
    TABLE_LAYOUT = 8
    """A table whose columns lack what its format needs (`TableLayoutError`)."""
    UNWRITABLE_ROWS = 9
    """Rows that cannot be valid VCF; `--resolved-only` writes the rest."""
    UNSORTED_ROWS = 10
    SEVERAL_FACETS = 11
    """The file loads as several facets and `--facet` did not choose one for a one-facet step."""


REDACTED = "<redacted>"

app = typer.Typer(
    name="just-vcf",
    help="Sniff, load, convert and check a personal genomics file. Prints counts and aggregates only.",
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_enable=False,
)

_JSON = Annotated[bool, typer.Option("--json", help="Print the report as JSON.")]
_FILE = Annotated[Path, typer.Argument(exists=True, dir_okay=False, help="The genotype file.")]
_SAMPLE = Annotated[
    str | None,
    typer.Option(help="The sample to load from a multi-sample VCF or table. Never echoed in the output."),
]
_FACET = Annotated[
    str | None,
    typer.Option(help="The facet to use when the file loads as several: its index in `load`'s list, or its facet_id."),
]
_REFERENCE = Annotated[
    Path | None,
    typer.Option(
        help="A reference store directory, or a FASTA added to the store at JUST_VCF_REFERENCE_DIR. "
        "Default: that store."
    ),
]


class _Redactor:
    """Replaces the labels this run knows (paths, file names, sample names, members) in library text."""

    def __init__(self) -> None:
        self.secrets: set[str] = set()

    def add(self, value: str | None) -> None:
        if value:
            self.secrets.add(value)

    def add_path(self, path: Path) -> None:
        self.add(str(path))
        self.add(str(path.absolute()))
        self.add(str(path.resolve()))
        self.add(path.name)
        stem = path.name.split(".", 1)[0]
        # A short stem ("a", "in") would cut ordinary words out of a warning; a name is longer.
        if len(stem) >= 3:
            self.add(stem)

    def add_member(self, member: str | None) -> None:
        if member:
            self.add(member)
            for part in re.split(r"[!/\\]", member):
                if part:
                    self.add_path(Path(part))

    def __call__(self, text: str) -> str:
        for secret in sorted(self.secrets, key=len, reverse=True):
            text = text.replace(secret, REDACTED)
        return text


class _Refusal(Exception):
    def __init__(self, code: ExitCode, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(message)


@app.callback()
def main(
    ctx: typer.Context,
    traceback: Annotated[
        bool,
        typer.Option("--traceback", help="Re-raise a crash with its traceback. It can print paths and sample names."),
    ] = False,
) -> None:
    """Sniff, load, convert and check a personal genomics file. Prints counts and aggregates only."""
    load_dotenv(find_dotenv(usecwd=True))
    # The library logs every warning it also returns, and polars-bio's import configures the root logger
    # (`logging.basicConfig()`), so each would reach stderr unredacted. The CLI prints them itself,
    # redacted and once; logging stays off while the command runs, and comes back when it closes.
    logging.disable(logging.CRITICAL)
    ctx.call_on_close(lambda: logging.disable(logging.NOTSET))
    ctx.obj = {"traceback": traceback}


@app.command("sniff")
def sniff_command(ctx: typer.Context, file: _FILE, as_json: _JSON = False) -> None:
    """What the file is: format, family, container layers and the build its header names.

    Exits 0 when a format is detected, 3 otherwise (the report still says what was found)."""
    redact = _Redactor()
    redact.add_path(file)
    with _guarded(ctx, redact):
        found = sniff(file)
        redact.add_member(found.member)
        for member in found.members:
            redact.add_member(member.member)
        report = _sniffed(found, redact)
        if found.members:
            report["members"] = [_sniffed(member, redact) for member in found.members]
        _emit(report, as_json)
    if not found:
        raise typer.Exit(ExitCode.UNSUPPORTED)


@app.command("load")
def load_command(ctx: typer.Context, file: _FILE, sample: _SAMPLE = None, as_json: _JSON = False) -> None:
    """Load the file and summarise each facet: what it is, how sure, and how many rows."""
    redact = _redactor(file, sample)
    with _guarded(ctx, redact, sample):
        loaded = load(file, sample=sample)
        facets = _learn(loaded, redact)
        report = {"facets": _facet_summaries(loaded, facets), "warnings": _warnings(loaded.warnings, redact)}
        _emit(report, as_json)


@app.command("check")
def check_command(
    ctx: typer.Context,
    file: _FILE,
    sample: _SAMPLE = None,
    reference: _REFERENCE = None,
    vrs: Annotated[bool, typer.Option("--vrs", help="Mint VRS ids too, so the withheld counts say why each missing id is missing.")] = False,
    as_json: _JSON = False,
) -> None:
    """Load the file, read REF from the reference, and count ref_check outcomes and withheld VRS ids per facet."""
    redact = _redactor(file, sample)
    with _guarded(ctx, redact, sample):
        store, collection = _reference(reference)
        loaded = load(file, sample=sample)
        facets = _learn(loaded, redact)
        checked = []
        for facet_id in facets["facet_id"].to_list():
            one = _one_facet(loaded, facet_id)
            resolved = resolve_ref(one, store, collection=collection)
            if vrs:
                resolved = mint_vrs_ids(resolved, store)
            checked.append(_check_summary(resolved, len(one.warnings), redact))
        _emit({"facets": checked, "load_warnings": _warnings(loaded.warnings, redact)}, as_json)


@app.command("convert")
def convert_command(
    ctx: typer.Context,
    file: _FILE,
    out: Annotated[Path, typer.Argument(dir_okay=False, help="The VCF to write (BGZF with a tabix index, unless --plain).")],
    sample: _SAMPLE = None,
    facet: _FACET = None,
    reference: _REFERENCE = None,
    normalize_rows: Annotated[bool, typer.Option("--normalize", help="Split multi-allelic sites, left-align and trim indels.")] = False,
    lift: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="A UCSC chain file to lift the facet through.")] = None,
    target_reference: Annotated[
        Path | None,
        typer.Option(help="The lift's target assembly: a FASTA, added to the reference store, or that same store."),
    ] = None,
    vrs: Annotated[bool, typer.Option("--vrs", help="Mint VRS ids and count them (the VCF written does not carry them yet).")] = False,
    gvcf: Annotated[bool, typer.Option("--gvcf", help="Write the reference blocks too, as gVCF.")] = False,
    resolved_only: Annotated[bool, typer.Option("--resolved-only", help="Write the rows that can be VCF and count the rest.")] = False,
    plain: Annotated[bool, typer.Option("--plain", help="Write plain text, without BGZF or an index.")] = False,
    as_json: _JSON = False,
) -> None:
    """Convert the file to VCF (or gVCF), optionally reading REF, lifting, normalizing and minting VRS ids.

    The steps run in this order: load, resolve_ref, liftover, normalize, mint_vrs_ids, write_vcf. REF is
    read whenever a reference is given or a later step needs it."""
    if target_reference is not None and lift is None:
        raise typer.BadParameter("--target-reference names the lift's target; it needs --lift", param_hint="--target-reference")
    redact = _redactor(file, sample)
    redact.add_path(out)
    with _guarded(ctx, redact, sample):
        needs_ref = reference is not None or normalize_rows or lift is not None or vrs or gvcf
        store, collection = _reference(reference) if needs_ref else (None, None)
        loaded = load(file, sample=sample)
        facets = _learn(loaded, redact)
        loaded = _one_facet(loaded, _choose_facet(facets, facet))
        steps = ["load"]
        report: dict[str, Any] = {}
        if store is not None:
            loaded = resolve_ref(loaded, store, collection=collection)
            steps.append("resolve_ref")
        if lift is not None:
            target = _target(store, target_reference)
            lifted = liftover(loaded, lift, store, target=target)
            loaded = lifted.loaded
            report["liftover"] = {"unlifted_by_reason": lifted.reasons, "unlifted": sum(lifted.reasons.values())}
            steps.append("liftover")
        if normalize_rows:
            loaded = normalize(loaded, store)
            steps.append("normalize")
        if vrs:
            loaded = mint_vrs_ids(loaded, store)
            report["vrs"] = _vrs_counts(loaded.calls)
            report["vrs"]["note"] = "the VCF written does not carry VRS ids yet"
            steps.append("mint_vrs_ids")
        written = write_vcf(loaded, out, resolved_only=resolved_only, compress=not plain, blocks=gvcf, reference=store if gvcf else None)
        steps.append("write_vcf")
        final = loaded.facets.collect()
        report = {
            "steps": steps,
            "facet": _facet_summaries(loaded, final)[0],
            **report,
            "written": {
                "rows_written": written.rows_written,
                "rows_skipped": dict(written.rows_skipped),
                "blocks_written": written.blocks_written,
                "blocks_not_written": written.blocks_not_written,
                "sha256": written.sha256,
                "indexed": written.index_path is not None,
            },
            "warnings": _warnings(loaded.warnings + written.warnings, redact),
        }
        _emit(report, as_json)


# --- what a command prints ------------------------------------------------------------------------


def _sniffed(found: Sniffed, redact: _Redactor) -> dict[str, Any]:
    return {
        "outcome": found.outcome.value,
        "format": found.format,
        "family": found.family.value if found.family is not None else None,
        "loadable": found.format in LOADABLE_FORMATS,
        "kind": found.kind,
        "layers": list(found.layers),
        "build": found.build,
        "build_evidence": found.build_evidence,
        "also_matched": list(found.also_matched),
        "reason": redact(found.reason) if found.reason else None,
        "warnings": _warnings(found.warnings, redact),
    }


_FACET_FIELDS = (
    "facet_id",
    "format_id",
    "family",
    "vendor",
    "caller",
    "caller_version",
    "relation",
    "derived_from",
    "vcf_version",
    "assembly_label",
    "assembly_label_basis",
    "ref_seqcol_digest",
    "callset_type",
    "absence_means",
    "chip",
    "chip_basis",
    "content_filter",
    "content_filter_basis",
    "phasing_method",
    "strand_check",
)
"""The facet columns a summary prints. Never `source_path`, `member`, `sample` or `sample_meta`: they
are labels, and can carry a person's name."""


def _facet_summaries(loaded: Loaded, facets: pl.DataFrame) -> list[dict[str, Any]]:
    calls, blocks = pl.collect_all(
        [
            loaded.calls.group_by("facet_id").agg(pl.len().alias("rows"), pl.col("chrom").n_unique().alias("contigs")),
            loaded.blocks.group_by("facet_id").agg(pl.len().alias("blocks")),
        ]
    )
    rows = {row["facet_id"]: row for row in calls.iter_rows(named=True)}
    block_counts = dict(blocks.iter_rows())
    summaries = []
    for index, facet in enumerate(facets.iter_rows(named=True)):
        summary: dict[str, Any] = {"index": index, **{name: facet[name] for name in _FACET_FIELDS}}
        summary["variant_classes"] = {entry["variant_class"]: entry["alleles"] for entry in facet["variant_classes"] or []} or None
        summary["dropped_fields"] = facet["dropped_fields"]
        counted = rows.get(facet["facet_id"], {})
        summary["rows"] = counted.get("rows", 0)
        summary["contigs"] = counted.get("contigs", 0)
        summary["blocks"] = block_counts.get(facet["facet_id"], 0)
        summaries.append(summary)
    return summaries


def _check_summary(resolved: Loaded, inherited: int, redact: _Redactor) -> dict[str, Any]:
    facet = resolved.facets.collect().row(0, named=True)
    ref_check = dict(resolved.calls.group_by("ref_check").len().collect().iter_rows())
    return {
        "facet_id": facet["facet_id"],
        "format_id": facet["format_id"],
        "assembly_label": facet["assembly_label"],
        "assembly_label_basis": facet["assembly_label_basis"],
        "ref_seqcol_digest": facet["ref_seqcol_digest"],
        "strand_check": facet["strand_check"],
        "rows": sum(ref_check.values()),
        "ref_check": {str(key): ref_check[key] for key in sorted(ref_check, key=str)},
        "vrs": _vrs_counts(resolved.calls),
        "warnings": _warnings(resolved.warnings[inherited:], redact),
    }


def _vrs_counts(calls: pl.LazyFrame) -> dict[str, Any]:
    """Allele slots (`[ref] + alts`) with a minted id, and those withheld per reason (`vrs_withheld`)."""
    minted, withheld = pl.collect_all(
        [
            calls.select(pl.col("vrs_allele_ids").explode(empty_as_null=False).is_not_null().sum()),
            vrs_withheld(calls).select(pl.col("vrs_withheld").explode(empty_as_null=False).alias("reason")).drop_nulls().group_by("reason").len(),
        ]
    )
    reasons = dict(withheld.iter_rows())
    return {"minted": minted.item(), "withheld_by_reason": {reason: reasons[reason] for reason in sorted(reasons)}}


def _warnings(warnings: tuple[str, ...], redact: _Redactor) -> list[str]:
    """One line per distinct warning, redacted, with a count where the same line came more than once."""
    counted = Counter(redact(warning) for warning in warnings)
    return [line if n == 1 else f"{line} (x{n})" for line, n in counted.items()]


def _emit(report: dict[str, Any], as_json: bool) -> None:
    if as_json:
        typer.echo(json.dumps(report, indent=2, default=str))
        return
    for line in _render(report, 0):
        typer.echo(line)


def _render(value: dict[str, Any], indent: int) -> Iterator[str]:
    pad = " " * indent
    for key, item in value.items():
        if isinstance(item, dict):
            if not item:
                yield f"{pad}{key}: none"
                continue
            yield f"{pad}{key}:"
            yield from _render(item, indent + 2)
        elif isinstance(item, list) and any(isinstance(entry, dict) for entry in item):
            yield f"{pad}{key}:"
            for entry in item:
                lines = list(_render(entry, indent + 4))
                if lines:
                    yield f"{pad}  - {lines[0].lstrip()}"
                    yield from lines[1:]
        elif isinstance(item, list):
            if key == "warnings":
                yield f"{pad}{key}:" if item else f"{pad}{key}: none"
                for entry in item:
                    yield f"{pad}  - {entry}"
            else:
                yield f"{pad}{key}: {', '.join(str(entry) for entry in item) if item else 'none'}"
        else:
            yield f"{pad}{key}: {_scalar(item)}"


def _scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


# --- plumbing ---------------------------------------------------------------------------------------


def _redactor(file: Path, sample: str | None) -> _Redactor:
    redact = _Redactor()
    redact.add_path(file)
    redact.add(sample)
    return redact


def _learn(loaded: Loaded, redact: _Redactor) -> pl.DataFrame:
    """The facets frame, with every label it holds added to the redaction."""
    facets = loaded.facets.collect()
    for facet in facets.iter_rows(named=True):
        redact.add(facet["sample"])
        redact.add_member(facet["member"])
        if facet["source_path"]:
            redact.add_path(Path(facet["source_path"]))
    return facets


def _choose_facet(facets: pl.DataFrame, choice: str | None) -> str:
    ids = facets["facet_id"].to_list()
    if choice is None:
        if len(ids) == 1:
            return ids[0]
        raise _Refusal(ExitCode.SEVERAL_FACETS, f"the file loads as {len(ids)} facets; choose one with --facet (an index or a facet_id from `load`)")
    if choice in ids:
        return choice
    if choice.isdigit() and int(choice) < len(ids):
        return ids[int(choice)]
    raise typer.BadParameter(f"no facet {choice!r}: the file loads as {len(ids)} facet(s)", param_hint="--facet")


def _one_facet(loaded: Loaded, facet_id: str) -> Loaded:
    """`loaded` narrowed to one facet. The warnings stay whole: a load's warnings are not split by facet."""
    keep = pl.col("facet_id") == facet_id
    return Loaded(
        facets=loaded.facets.filter(keep),
        calls=loaded.calls.filter(keep),
        blocks=loaded.blocks.filter(keep),
        warnings=loaded.warnings,
        schema_version=loaded.schema_version,
    )


def _reference(reference: Path | None) -> tuple[ReferenceStore, str | None]:
    """The store to read and, for a FASTA, the one collection to use: resolved once, then shared by
    every step, so `normalize`, `liftover`, `mint_vrs_ids` and `write_vcf` read the collection
    `resolve_ref` read from."""
    if reference is None:
        return ReferenceStore.open(), None
    if reference.is_dir():
        if not RefgetStore.store_exists(str(reference)):
            raise typer.BadParameter("not a reference store (and a directory is no FASTA)", param_hint="--reference")
        return ReferenceStore.open(reference), None
    if reference.is_file():
        store = ReferenceStore.open()
        return store, store.add_fasta(reference)
    raise typer.BadParameter("neither a reference store nor a FASTA file", param_hint="--reference")


def _target(store: ReferenceStore, target_reference: Path | None) -> str | None:
    """The lift's target collection: a FASTA's digest once added to `store`, or None to let `liftover`
    choose among the store's collections (the store itself was named)."""
    if target_reference is None:
        return None
    if target_reference.is_file():
        return store.add_fasta(target_reference)
    if target_reference.is_dir() and store.path is not None and target_reference.resolve() == store.path.resolve():
        return None
    raise typer.BadParameter(
        "liftover reads the source and the target from one store: pass the target's FASTA, or the store given as --reference",
        param_hint="--target-reference",
    )


@contextmanager
def _guarded(ctx: typer.Context, redact: _Redactor, sample: str | None = None) -> Iterator[None]:
    """Turn the library's refusals into their exit codes, and any other error into a crash whose
    message is withheld. Every message printed is the CLI's own or redacted."""
    try:
        yield
    except (typer.Exit, typer.BadParameter, typer.Abort):
        raise
    except _Refusal as err:
        _refuse(err.code, err.message)
    except UnsupportedFormatError as err:
        found = err.sniffed
        what = found.format or found.kind or found.outcome.value
        why = "no loader reads this format yet" if found.outcome is SniffOutcome.DETECTED else "not a format just-vcf recognises"
        _refuse(ExitCode.UNSUPPORTED, f"refused: {why} ({what})")
    except MultiSampleError:
        if sample is None:
            _refuse(ExitCode.MULTI_SAMPLE, "refused: the file holds more than one sample; choose one with --sample")
        _refuse(ExitCode.MULTI_SAMPLE, "refused: the sample given with --sample is not in the file")
    except NoSampleError:
        _refuse(ExitCode.NO_SAMPLE, "refused: a sites-only VCF, with no sample column to load")
    except TableLayoutError as err:
        _refuse(ExitCode.TABLE_LAYOUT, f"refused: {redact(str(err))}")
    except DamagedBgzfError as err:
        where = f"facet {err.facet_id}" if err.facet_id is not None else f"source sha256 {err.source_sha256}"
        _refuse(ExitCode.DAMAGED_BGZF, f"refused: {where}: the BGZF input is damaged at compressed byte {err.offset} ({err.reason})")
    except ReferenceMismatchError as err:
        _refuse(ExitCode.REFERENCE_MISMATCH, f"refused: {redact(str(err))}")
    except UnwritableRowsError as err:
        counts = ", ".join(f"{n} {reason}" for reason, n in err.counts.items())
        _refuse(ExitCode.UNWRITABLE_ROWS, f"refused: {sum(err.counts.values())} row(s) cannot be written as VCF ({counts}); --resolved-only writes the rest")
    except UnsortedRowsError:
        _refuse(ExitCode.UNSORTED_ROWS, "refused: the rows are not in the order load emits them")
    except Exception as err:
        if ctx.obj and ctx.obj.get("traceback"):
            raise
        typer.echo(
            f"crash: {type(err).__name__}; its message is withheld, since it can carry paths or sample names "
            "(rerun with --traceback before the command to see it)",
            err=True,
        )
        raise typer.Exit(ExitCode.CRASH) from None


def _refuse(code: ExitCode, message: str) -> None:
    typer.echo(message, err=True)
    raise typer.Exit(code)
