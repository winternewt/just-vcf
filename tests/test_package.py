import shutil
import subprocess
import tarfile
import zipfile
from importlib.metadata import version
from pathlib import Path

import pytest

import just_vcf

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "just_vcf"

# What a source distribution may carry besides the package itself. The repository's docs/ hold
# documents derived from private genomes, and its root holds agent and team files: an sdist built from
# hatchling's default selection (every file git does not ignore) would publish them with the package.
SDIST_EXTRA = frozenset({"pyproject.toml", "README.md", "LICENSE", "PKG-INFO", ".gitignore"})


def symlinks(root: Path) -> list[str]:
    """Every symlink under `root`, file or directory, relative to it; a linked directory is not entered."""
    found = []
    for directory, dirnames, filenames in root.walk(follow_symlinks=False):
        found += [str((directory / name).relative_to(root)) for name in (*dirnames, *filenames) if (directory / name).is_symlink()]
    return sorted(found)


def test_symlinks_are_found_whether_they_name_a_file_or_a_directory(tmp_path: Path) -> None:
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "NOTES.md").write_text("private")
    package = tmp_path / "src" / "pkg"
    (package / "sub").mkdir(parents=True)
    (package / "internal_docs").symlink_to(docs, target_is_directory=True)
    (package / "sub" / "notes.md").symlink_to(docs / "NOTES.md")
    assert symlinks(package) == ["internal_docs", "sub/notes.md"]


def test_the_package_holds_no_symlink() -> None:
    """hatchling follows a symlink under the package: `internal_docs -> ../../docs` would ship docs/ in
    the wheel and the sdist under the package's own prefix, past the allow-lists the tests below read."""
    assert symlinks(PACKAGE) == []


def test_version_comes_from_distribution_metadata() -> None:
    assert just_vcf.__version__ == version("just-vcf")


@pytest.fixture(scope="module")
def dists(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("uv is not on PATH")
    out = tmp_path_factory.mktemp("dist")
    done = subprocess.run([uv, "build", "--out-dir", str(out), str(ROOT)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-2000:]
    name = f"just_vcf-{version('just-vcf')}"
    return out / f"{name}.tar.gz", out / f"{name}-py3-none-any.whl"


def test_sdist_carries_only_the_package_and_its_metadata(dists: tuple[Path, Path]) -> None:
    sdist, _ = dists
    prefix = f"just_vcf-{version('just-vcf')}/"
    with tarfile.open(sdist) as tar:
        members = [m.name for m in tar.getmembers() if m.isfile()]
    assert members and all(m.startswith(prefix) for m in members)
    paths = {m.removeprefix(prefix) for m in members}
    stray = sorted(p for p in paths if not p.startswith("src/just_vcf/") and p not in SDIST_EXTRA)
    assert stray == []
    assert {"pyproject.toml", "README.md", "LICENSE", "PKG-INFO", "src/just_vcf/__init__.py"} <= paths


def test_wheel_carries_only_the_package_and_its_dist_info(dists: tuple[Path, Path]) -> None:
    _, wheel = dists
    info = f"just_vcf-{version('just-vcf')}.dist-info/"
    with zipfile.ZipFile(wheel) as whl:
        names = [n for n in whl.namelist() if not n.endswith("/")]
    stray = sorted(n for n in names if not n.startswith(("just_vcf/", info)))
    assert stray == []
    assert {"just_vcf/__init__.py", "just_vcf/py.typed", f"{info}licenses/LICENSE"} <= set(names)


def test_sdist_package_files_match_the_wheel(dists: tuple[Path, Path]) -> None:
    sdist, wheel = dists
    prefix = f"just_vcf-{version('just-vcf')}/src/"
    with tarfile.open(sdist) as tar:
        source = {m.name.removeprefix(prefix) for m in tar.getmembers() if m.isfile() and m.name.startswith(prefix)}
    with zipfile.ZipFile(wheel) as whl:
        built = {n for n in whl.namelist() if n.startswith("just_vcf/") and not n.endswith("/")}
    assert source == built
