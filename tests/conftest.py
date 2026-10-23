import os
import sys
from pathlib import Path

import pytest
from dotenv import dotenv_values, find_dotenv

# scripts/ holds dev tools that import each other by module name (pgp_participants imports pgp_leech).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Loader cache under a temporary directory, never the user's real cache."""
    path = tmp_path_factory.mktemp("cache")
    monkeypatch.setenv("JUST_VCF_CACHE_DIR", str(path))
    # No configured reference store: a test that wants one passes it or sets the variable itself.
    monkeypatch.setenv("JUST_VCF_REFERENCE_DIR", "")
    # No configured chip manifests either (RM281): a test that reads the real ones asks for `manifests_dir`.
    monkeypatch.setenv("JUST_VCF_MANIFESTS_DIR", "")
    return path


@pytest.fixture(scope="session")
def samples_dir() -> Path:
    """`JUST_VCF_SAMPLES_DIR`: the environment wins, an empty value included; else the repo `.env`.

    Read without exporting, so a test that neutralises the variable with `setenv(…, "")` still skips.
    """
    if "JUST_VCF_SAMPLES_DIR" in os.environ:
        value = os.environ["JUST_VCF_SAMPLES_DIR"]
    else:
        value = dotenv_values(find_dotenv(usecwd=True)).get("JUST_VCF_SAMPLES_DIR")
    if not value:
        pytest.skip("JUST_VCF_SAMPLES_DIR is not set")
    path = Path(value)
    if not path.is_dir():
        pytest.skip(f"{path} is absent")
    return path


@pytest.fixture(scope="session")
def ancestry_dir() -> Path:
    """`JUST_VCF_ANCESTRY_DIR`, the directory of reference panels, read like `samples_dir`."""
    if "JUST_VCF_ANCESTRY_DIR" in os.environ:
        value = os.environ["JUST_VCF_ANCESTRY_DIR"]
    else:
        value = dotenv_values(find_dotenv(usecwd=True)).get("JUST_VCF_ANCESTRY_DIR")
    if not value:
        pytest.skip("JUST_VCF_ANCESTRY_DIR is not set")
    path = Path(value)
    if not path.is_dir():
        pytest.skip(f"{path} is absent")
    return path


_MANIFESTS_DIR = os.environ.get("JUST_VCF_MANIFESTS_DIR")
"""`JUST_VCF_MANIFESTS_DIR` as the session started, before `isolated_cache` empties it for every test."""


@pytest.fixture(scope="session")
def manifests_dir() -> Path:
    """`JUST_VCF_MANIFESTS_DIR`, Illumina's manifest zips (docs/PROVENANCE.md), read like `samples_dir`.

    The library reads the variable too (RM281), so every test sees it empty; a test that reads the real
    manifests asks for this fixture and sets it itself.
    """
    value = _MANIFESTS_DIR if _MANIFESTS_DIR is not None else dotenv_values(find_dotenv(usecwd=True)).get("JUST_VCF_MANIFESTS_DIR")
    if not value:
        pytest.skip("JUST_VCF_MANIFESTS_DIR is not set")
    path = Path(value)
    if not path.is_dir():
        pytest.skip(f"{path} is absent")
    return path


@pytest.fixture(scope="session")
def test_reference_dir() -> Path:
    """`JUST_VCF_TEST_REFERENCE_DIR`, a real reference store (GRCh37 and GRCh38), read like `samples_dir`.

    Apart from `JUST_VCF_REFERENCE_DIR`, which every test sees empty: a test that reads a real assembly
    asks for this fixture, and only reads the store.
    """
    if "JUST_VCF_TEST_REFERENCE_DIR" in os.environ:
        value = os.environ["JUST_VCF_TEST_REFERENCE_DIR"]
    else:
        value = dotenv_values(find_dotenv(usecwd=True)).get("JUST_VCF_TEST_REFERENCE_DIR")
    if not value:
        pytest.skip("JUST_VCF_TEST_REFERENCE_DIR is not set")
    path = Path(value)
    if not path.is_dir():
        pytest.skip(f"{path} is absent")
    return path


@pytest.fixture(scope="session")
def pgp_root(samples_dir: Path) -> Path:
    """The PGP-HMS download tree under the samples directory (docs/PGP_HMS.md)."""
    root = samples_dir / "pgp-hms"
    if not root.is_dir():
        pytest.skip(f"{root} is absent")
    return root


@pytest.fixture(scope="session")
def pgp_participants(pgp_root: Path) -> Path:
    """`pgp-hms/participants/<hu>/<source>-<id>/<file>`, as written by scripts/pgp_participants.py."""
    root = pgp_root / "participants"
    if not root.is_dir():
        pytest.skip(f"{root} is absent")
    return root
