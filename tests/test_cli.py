"""RM272: the `just-vcf` command line, through Typer's CliRunner on synthetic files.

Every expected count is recomputed from the library's own frames or from the synthetic reference, and
every test that names a sample or a person-like file checks the name never reaches the output.
"""

import json
import logging
import gzip
import zipfile
from collections import Counter
from importlib.metadata import entry_points
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from just_vcf import ReferenceStore, load, normalize, resolve_ref
from just_vcf import cli
from just_vcf.cli import ExitCode, app
from just_vcf.formats import get_spec
from just_vcf.normalize import RELATION as NORMALIZED
from just_vcf.sniff import sniff
from test_bgzf import _damaged, _vcf_bgzf
import test_liftover
import test_write_gvcf
from test_load_rsid import _two_exports
from test_write_vcf import as_set
from test_reference import LENGTHS, SEQ, base, contig_lines, other, write_fasta

runner = CliRunner()

PERSON = "Jane_Roe"
"""A person-like label used as a file name, a sample name and an archive member name."""


def run(*args: str | Path) -> object:
    return runner.invoke(app, [str(arg) for arg in args])


def report(result) -> dict:
    assert result.exit_code == 0, (result.exit_code, result.stdout, result.stderr)
    return json.loads(result.stdout)


def assert_private(result, *labels: str) -> None:
    for label in labels:
        assert label not in result.stdout
        assert label not in result.stderr


def vcf(path: Path, records: list[tuple[str, int, str, str, str]], samples: tuple[str, ...] = (PERSON,), contigs: list[str] | None = None) -> Path:
    """A VCF of `(chrom, pos, ref, alt, gt)` records, the same genotype for every sample."""
    meta = ["##fileformat=VCFv4.2", *(contigs if contigs is not None else contig_lines(LENGTHS)), '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", *(["FORMAT", *samples] if samples else [])])
    body = ["\t".join([c, str(p), ".", r, a, ".", "PASS", ".", *(["GT", *([g] * len(samples))] if samples else [])]) for c, p, r, a, g in records]
    path.write_text("\n".join([*meta, header, *body]) + "\n")
    return path


# Positions on the synthetic reference: REF as the sequence has it, and two rows whose REF disagrees.
MATCHING = [("1", 100, "A", "0/1"), ("1", 250, "C", "1/1"), ("2", 40, "G", "0/1"), ("X", 77, "T", "0/1")]
MISMATCHING = [("1", 500, "1/1"), ("2", 900, "0/1")]


def checked_vcf(path: Path) -> Path:
    records = [(c, p, base(c, p), other(base(c, p)), g) for c, p, _, g in MATCHING]
    # A REF the sequence contradicts, with the sequence's own base as the ALT.
    records += [(c, p, other(base(c, p)), base(c, p), g) for c, p, g in MISMATCHING]
    records.sort(key=lambda r: (list(LENGTHS).index(r[0]), r[1]))
    return vcf(path, records)


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    return write_fasta(tmp_path / "ref.fa", SEQ)


# --- sniff --------------------------------------------------------------------------------------------


def test_sniff_reports_the_format_and_exits_zero(tmp_path: Path) -> None:
    path = vcf(tmp_path / f"{PERSON}.vcf", [("1", 100, "A", "G", "0/1")])
    found = report(run("sniff", path, "--json"))
    sniffed = sniff(path)
    assert (found["outcome"], found["format"], found["family"], found["loadable"]) == ("detected", sniffed.format, sniffed.family.value, True)
    text = run("sniff", path)
    assert text.exit_code == 0 and "format: vcf" in text.stdout
    assert_private(text, PERSON, str(tmp_path))


def test_sniff_of_something_not_recognised_exits_with_the_unsupported_code(tmp_path: Path) -> None:
    path = tmp_path / f"{PERSON}.txt"
    path.write_text("just some notes\n")
    result = run("sniff", path, "--json")
    assert result.exit_code == ExitCode.UNSUPPORTED
    assert json.loads(result.stdout)["outcome"] == "unknown"
    assert_private(result, PERSON)


def test_sniff_does_not_name_an_archive_member(tmp_path: Path) -> None:
    path = tmp_path / "export.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{PERSON}_genome.txt", get_spec("23andme").example)
    result = run("sniff", path, "--json")
    assert report(result)["format"] == "23andme"
    assert_private(result, PERSON)


# --- load ---------------------------------------------------------------------------------------------


def test_load_summarises_each_facet_with_the_library_s_counts(tmp_path: Path) -> None:
    path = vcf(tmp_path / f"{PERSON}.vcf", [("1", 100, "A", "G", "0/1"), ("1", 200, "C", "T", "1/1"), ("2", 300, "G", "GA", "0/1")])
    result = run("load", path, "--json")
    summary = report(result)
    loaded = load(path)
    facets, calls = pl.collect_all([loaded.facets, loaded.calls])
    [facet] = summary["facets"]
    assert facet["facet_id"] == facets["facet_id"].item()
    assert facet["rows"] == calls.height
    assert facet["contigs"] == calls["chrom"].n_unique()
    assert facet["blocks"] == 0
    assert facet["absence_means"] == facets["absence_means"].item()
    assert facet["variant_classes"] == {e["variant_class"]: e["alleles"] for e in facets["variant_classes"].item()}
    assert not {"source_path", "sample", "member", "sample_meta"} & facet.keys()
    assert_private(result, PERSON, str(tmp_path))


def test_load_of_several_facets_lists_each(tmp_path: Path) -> None:
    text, _, _ = _two_exports()
    path = tmp_path / "genome.txt"
    path.write_text(text)
    summary = report(run("load", path, "--json"))
    calls = load(path).calls.collect()
    expected = dict(calls.group_by("facet_id").len().iter_rows())
    assert {f["facet_id"]: f["rows"] for f in summary["facets"]} == expected
    assert [f["index"] for f in summary["facets"]] == list(range(len(expected)))


def test_sample_choice_is_taken_and_never_echoed(tmp_path: Path) -> None:
    other_person = "John_Doe"
    path = vcf(tmp_path / "family.vcf", [("1", 100, "A", "G", "0/1")], samples=(PERSON, other_person))

    several = run("load", path)
    assert several.exit_code == ExitCode.MULTI_SAMPLE and "--sample" in several.stderr
    assert_private(several, PERSON, other_person)

    chosen = run("load", path, "--sample", other_person, "--json")
    assert report(chosen)["facets"][0]["rows"] == 1
    assert_private(chosen, PERSON, other_person)

    absent = run("load", path, "--sample", "Nobody_Here")
    assert absent.exit_code == ExitCode.MULTI_SAMPLE
    assert_private(absent, PERSON, other_person, "Nobody_Here")


def test_refusals_have_their_own_exit_codes(tmp_path: Path) -> None:
    sites_only = vcf(tmp_path / "sites.vcf", [("1", 100, "A", "G", "0/1")], samples=())
    assert run("load", sites_only).exit_code == ExitCode.NO_SAMPLE

    no_loader = tmp_path / f"{PERSON}_summary.tsv"
    no_loader.write_text(get_spec("cgi_summary").example)
    refused = run("load", no_loader)
    assert refused.exit_code == ExitCode.UNSUPPORTED and "no loader" in refused.stderr
    assert_private(refused, PERSON, str(tmp_path))

    damaged = tmp_path / f"{PERSON}.vcf.gz"
    starts = _vcf_bgzf(damaged)
    offset = _damaged(damaged, starts, len(starts) // 2, "crc32")
    result = run("load", damaged)
    assert result.exit_code == ExitCode.DAMAGED_BGZF and str(offset) in result.stderr
    assert_private(result, PERSON, str(tmp_path))

    assert run("load", tmp_path / "absent.vcf").exit_code == ExitCode.USAGE


def test_a_crash_is_told_from_a_refusal_and_its_message_withheld(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = vcf(tmp_path / "in.vcf", [("1", 100, "A", "G", "0/1")])

    def broken(*args, **kwargs):
        raise RuntimeError(f"cannot read {tmp_path}/{PERSON}.vcf")

    monkeypatch.setattr(cli, "load", broken)
    result = run("load", path)
    assert result.exit_code == ExitCode.CRASH and "RuntimeError" in result.stderr
    assert_private(result, PERSON, str(tmp_path))
    with pytest.raises(RuntimeError):
        runner.invoke(app, ["--traceback", "load", str(path)], catch_exceptions=False)


def test_library_logging_is_off_while_a_command_runs(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """polars-bio's import configures the root logger, so a library warning would reach stderr as logged,
    past the redaction; the CLI prints it itself, once."""
    text, _, _ = _two_exports()
    path = tmp_path / "genome.txt"
    path.write_text(text)
    with caplog.at_level(logging.WARNING):
        summary = report(run("load", path, "--json"))
    assert summary["warnings"] and not [r for r in caplog.records if r.name.startswith("just_vcf")]
    assert logging.root.manager.disable == logging.NOTSET


def test_a_library_warning_that_names_a_label_is_redacted() -> None:
    redact = cli._Redactor()
    redact.add_path(Path(f"/data/genomes/{PERSON}.vcf.gz"))
    redact.add("SAMPLE_7")
    warnings = (f"/data/genomes/{PERSON}.vcf.gz: 3 rows dropped", "SAMPLE_7 has no GT", "SAMPLE_7 has no GT", "2 rows past the end")
    assert cli._warnings(warnings, redact) == [f"{cli.REDACTED}: 3 rows dropped", f"{cli.REDACTED} has no GT (x2)", "2 rows past the end"]


# --- check --------------------------------------------------------------------------------------------


def test_check_counts_ref_check_against_the_reference(tmp_path: Path, fasta: Path) -> None:
    path = checked_vcf(tmp_path / f"{PERSON}.vcf")
    result = run("check", path, "--reference", fasta, "--json")
    [facet] = report(result)["facets"]
    assert facet["ref_check"] == {"match": len(MATCHING), "mismatch": len(MISMATCHING)}
    assert facet["rows"] == len(MATCHING) + len(MISMATCHING)
    assert facet["ref_seqcol_digest"] == ReferenceStore.open().add_fasta(fasta)
    # Not minted: every slot is withheld, the mismatching rows for their REF, the rest as not run.
    withheld = facet["vrs"]["withheld_by_reason"]
    assert facet["vrs"]["minted"] == 0
    assert withheld == {"not_run": 2 * len(MATCHING), "ref_mismatch": 2 * len(MISMATCHING)}
    assert_private(result, PERSON, str(tmp_path))


def test_check_with_vrs_mints_the_matching_rows(tmp_path: Path, fasta: Path) -> None:
    path = checked_vcf(tmp_path / "in.vcf")
    [facet] = report(run("check", path, "--reference", fasta, "--vrs", "--json"))["facets"]
    assert facet["vrs"]["minted"] == 2 * len(MATCHING)
    assert facet["vrs"]["withheld_by_reason"] == {"ref_mismatch": 2 * len(MISMATCHING)}


def test_check_reads_the_reference_store_a_dotenv_names(tmp_path: Path, fasta: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    store = ReferenceStore.open(tmp_path / "store")
    store.add_fasta(fasta)
    path = checked_vcf(tmp_path / "in.vcf")
    monkeypatch.delenv("JUST_VCF_REFERENCE_DIR")
    monkeypatch.chdir(tmp_path)
    (tmp_path / ".env").write_text(f"JUST_VCF_REFERENCE_DIR={tmp_path / 'store'}\n")
    [facet] = report(run("check", path, "--json"))["facets"]
    assert facet["ref_check"]["match"] == len(MATCHING)


def test_check_refuses_a_reference_of_another_assembly(tmp_path: Path, fasta: Path) -> None:
    longer = {**LENGTHS, "1": LENGTHS["1"] + 10}
    path = vcf(tmp_path / "in.vcf", [("1", 100, base("1", 100), other(base("1", 100)), "0/1")], contigs=contig_lines(longer))
    assert run("check", path, "--reference", fasta).exit_code == ExitCode.REFERENCE_MISMATCH
    # An empty store is a refusal too, not a crash.
    assert run("check", path).exit_code == ExitCode.REFERENCE_MISMATCH


# --- convert ------------------------------------------------------------------------------------------


def test_convert_writes_a_vcf_that_loads_back_with_the_same_rows(tmp_path: Path) -> None:
    path = vcf(tmp_path / f"{PERSON}.vcf", [("1", 100, "A", "G", "0/1"), ("1", 200, "C", "T", "1/1"), ("2", 300, "G", "GA", "0/1")])
    out = tmp_path / "out.vcf.gz"
    result = run("convert", path, out, "--json")
    written = report(result)
    source, back = load(path).calls.collect(), load(out).calls.collect()
    columns = ["chrom", "pos", "source_ref", "alts", "gt"]
    assert written["written"]["rows_written"] == source.height
    assert as_set(back, columns) == as_set(source, columns)
    assert written["written"]["indexed"] and (tmp_path / "out.vcf.gz.tbi").is_file()
    assert written["steps"] == ["load", "write_vcf"]
    assert_private(result, PERSON, str(tmp_path))


def test_convert_of_an_array_needs_its_ref_read(tmp_path: Path, fasta: Path) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(get_spec("23andme").example)
    unread = run("convert", path, tmp_path / "out.vcf.gz")
    assert unread.exit_code == ExitCode.UNWRITABLE_ROWS and "--resolved-only" in unread.stderr
    assert not (tmp_path / "out.vcf.gz").exists()


def test_convert_normalizes_through_the_reference(tmp_path: Path, fasta: Path) -> None:
    path = checked_vcf(tmp_path / "in.vcf")
    out = tmp_path / "out.vcf"
    assert run("convert", path, out, "--reference", fasta, "--normalize").exit_code == ExitCode.UNWRITABLE_ROWS
    written = report(run("convert", path, out, "--reference", fasta, "--normalize", "--vrs", "--plain", "--resolved-only", "--json"))
    assert written["steps"] == ["load", "resolve_ref", "normalize", "mint_vrs_ids", "write_vcf"]
    expected = normalize(resolve_ref(load(path), fasta)).calls.collect()
    mismatched = expected.filter(pl.col("ref_check") == "mismatch").height
    assert written["written"]["rows_skipped"] == {"ref_mismatch": mismatched}
    assert written["written"]["rows_written"] == expected.height - mismatched
    assert written["facet"]["relation"] == NORMALIZED
    assert not written["written"]["indexed"] and out.read_text().startswith("##fileformat=VCFv4.3")
    assert written["vrs"]["minted"] == 2 * len(MATCHING)


def test_convert_writes_the_blocks_as_gvcf(tmp_path: Path) -> None:
    fasta = test_write_gvcf.write_fasta(tmp_path / "ref.fa")
    path = tmp_path / "in.g.vcf"
    path.write_text("\n".join([*test_write_gvcf.HEADER, *test_write_gvcf.RECORDS]) + "\n")
    out = tmp_path / "out.g.vcf.gz"
    written = report(run("convert", path, out, "--reference", fasta, "--gvcf", "--json"))
    loaded = load(path)
    blocks, calls = pl.collect_all([loaded.blocks, loaded.calls])
    assert written["written"]["blocks_written"] == blocks.height and written["written"]["blocks_not_written"] == 0
    assert written["written"]["rows_written"] == calls.height
    back = load(out).blocks.collect()
    assert back.select(test_write_gvcf.BLOCK_COLUMNS).equals(blocks.select(test_write_gvcf.BLOCK_COLUMNS))


def test_convert_lifts_to_a_target_fasta_added_to_the_store(tmp_path: Path) -> None:
    store_dir = tmp_path / "store"
    ReferenceStore.open(store_dir).add_fasta(test_liftover.write_fasta(tmp_path / "source.fa", test_liftover.SOURCE), assembly_label="GRCh37")
    target = test_liftover.write_fasta(tmp_path / "target.fa", test_liftover.TARGET)
    chain = tmp_path / "source_to_target.chain.gz"
    chain.write_bytes(gzip.compress(test_liftover.chain_text().encode()))
    order = sorted(test_liftover.CASES.values(), key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    path = tmp_path / f"{PERSON}.vcf"
    path.write_text("\n".join([*test_liftover.HEADER, *order]) + "\n")

    out = tmp_path / "out.vcf"
    args = ("convert", path, out, "--reference", store_dir, "--lift", chain, "--target-reference", target, "--resolved-only", "--plain", "--json")
    result = run(*args)
    written = report(result)
    expected = Counter(reason for reason in test_liftover.EXPECTED_REASON.values() if reason is not None)
    assert written["liftover"]["unlifted_by_reason"] == dict(sorted(expected.items()))
    assert written["steps"] == ["load", "resolve_ref", "liftover", "write_vcf"]
    assert written["facet"]["relation"] == "liftover_to"
    assert written["facet"]["ref_seqcol_digest"] == ReferenceStore.open(store_dir).add_fasta(target)
    lifted = len(test_liftover.CASES) - sum(expected.values())
    assert written["written"]["rows_written"] + sum(written["written"]["rows_skipped"].values()) == lifted
    assert_private(result, PERSON, str(tmp_path))


def test_convert_of_several_facets_takes_one_chosen(tmp_path: Path) -> None:
    text, _, _ = _two_exports()
    path = tmp_path / "genome.txt"
    path.write_text(text)
    out = tmp_path / "out.vcf"
    assert run("convert", path, out).exit_code == ExitCode.SEVERAL_FACETS
    facet_ids = load(path).facets.collect()["facet_id"].to_list()
    by_index = report(run("convert", path, out, "--facet", "1", "--resolved-only", "--plain", "--json"))
    assert by_index["facet"]["facet_id"] == facet_ids[1]
    by_id = report(run("convert", path, out, "--facet", facet_ids[0], "--resolved-only", "--plain", "--json"))
    assert by_id["facet"]["facet_id"] == facet_ids[0]
    assert run("convert", path, out, "--facet", "7").exit_code == ExitCode.USAGE


def test_convert_options_that_do_not_fit_are_usage_errors(tmp_path: Path, fasta: Path) -> None:
    path = vcf(tmp_path / "in.vcf", [("1", 100, "A", "G", "0/1")])
    out = tmp_path / "out.vcf.gz"
    assert run("convert", path, out, "--target-reference", fasta).exit_code == ExitCode.USAGE
    assert run("convert", path, out, "--reference", tmp_path).exit_code == ExitCode.USAGE


def test_the_entry_point_is_registered() -> None:
    [script] = [ep for ep in entry_points(group="console_scripts") if ep.name == "just-vcf"]
    assert script.value == "just_vcf.cli:app"
