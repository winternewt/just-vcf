from importlib.metadata import version

import just_vcf


def test_version_comes_from_distribution_metadata() -> None:
    assert just_vcf.__version__ == version("just-vcf")
