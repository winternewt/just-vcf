"""The suite's own harness: a repo `.env` never opts a plain run into private data."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

CONFTEST = Path(__file__).with_name("conftest.py")

_PROBE = '''
import os
from dotenv import find_dotenv, load_dotenv

def test_samples(samples_dir):
    raise AssertionError("samples_dir was handed to a test from .env")

def test_the_cli_dotenv_load_keeps_it_empty():
    load_dotenv(find_dotenv(usecwd=True))
    assert os.environ["JUST_VCF_SAMPLES_DIR"] == ""
    assert os.environ["JUST_VCF_MANIFESTS_DIR"] == ""
'''


def test_a_dotenv_does_not_opt_a_plain_run_in(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    run = tmp_path / "run"
    run.mkdir()
    (run / ".env").write_text(f"JUST_VCF_SAMPLES_DIR={genomes}\nJUST_VCF_MANIFESTS_DIR={genomes}\n")
    shutil.copy(CONFTEST, run / "conftest.py")
    (run / "test_probe.py").write_text(_PROBE)
    env = {k: v for k, v in os.environ.items() if not k.startswith("JUST_VCF_") and not k.startswith("PYTEST")}
    done = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-rs", str(run)],
        cwd=run, env=env, capture_output=True, text=True,
    )  # fmt: skip
    assert done.returncode == 0, done.stdout + done.stderr
    assert "1 passed, 1 skipped" in done.stdout
    assert "JUST_VCF_SAMPLES_DIR is not set" in done.stdout


def test_the_shell_still_opts_a_run_in(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    run = tmp_path / "run"
    run.mkdir()
    shutil.copy(CONFTEST, run / "conftest.py")
    (run / "test_probe.py").write_text("def test_samples(samples_dir):\n    assert samples_dir.name == 'genomes'\n")
    env = {k: v for k, v in os.environ.items() if not k.startswith("JUST_VCF_") and not k.startswith("PYTEST")}
    env["JUST_VCF_SAMPLES_DIR"] = str(genomes)
    done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(run)], cwd=run, env=env, capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    assert "1 passed" in done.stdout
