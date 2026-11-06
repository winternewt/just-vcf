"""The weighted delete-one-block jackknife (`just_vcf.jackknife`), against cases with a known answer."""

import math

import pytest

from just_vcf.jackknife import weighted_block_jackknife

# Unequal blocks, so that a formula that ignored the weights would give another answer.
BLOCK_MEANS = (0.31, -0.12, 0.05, 0.44, 0.18, -0.27, 0.09)
BLOCK_SIZES = (120, 35, 410, 18, 260, 77, 150)


def _weighted_mean(values, sizes) -> float:
    return sum(v * m for v, m in zip(values, sizes, strict=True)) / sum(sizes)


def test_a_weighted_mean_has_the_block_means_as_its_pseudo_values() -> None:
    """For the site-weighted mean of block means, block j's pseudo-value is its own mean exactly, so the
    jackknife estimate is the full mean and the variance is `(1/g) sum (x_j - mean)^2 / (h_j - 1)`,
    computed here from the block means alone, without the replicates."""
    n = sum(BLOCK_SIZES)
    g = len(BLOCK_SIZES)
    full = _weighted_mean(BLOCK_MEANS, BLOCK_SIZES)
    replicates = [
        (full * n - x * m) / (n - m) for x, m in zip(BLOCK_MEANS, BLOCK_SIZES, strict=True)
    ]
    estimate, se = weighted_block_jackknife(full, replicates, BLOCK_SIZES)
    expected = math.sqrt(
        sum((x - full) ** 2 / (n / m - 1.0) for x, m in zip(BLOCK_MEANS, BLOCK_SIZES, strict=True)) / g
    )
    assert estimate == pytest.approx(full, rel=1e-12)
    assert se == pytest.approx(expected, rel=1e-12)
    unweighted = weighted_block_jackknife(full, replicates, [1] * g)[1]
    assert unweighted != pytest.approx(se, rel=1e-3), "the block sizes change the answer"


def test_equal_blocks_give_the_ordinary_delete_one_jackknife() -> None:
    replicates = [0.52, 0.49, 0.55, 0.47, 0.51, 0.50]
    g = len(replicates)
    mean = sum(replicates) / g
    ordinary = math.sqrt((g - 1) / g * sum((r - mean) ** 2 for r in replicates))
    for size in (1, 40):
        assert weighted_block_jackknife(0.5, replicates, [size] * g)[1] == pytest.approx(ordinary, rel=1e-12)


@pytest.mark.parametrize(
    ("replicates", "sizes", "match"),
    [([0.1], [10], "two blocks"), ([0.1, 0.2], [10], "replicates but"), ([0.1, 0.2], [10, 0], "positive size")],
    ids=["one_block", "lengths", "empty_block"],
)
def test_refusals(replicates: list, sizes: list, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        weighted_block_jackknife(0.1, replicates, sizes)
