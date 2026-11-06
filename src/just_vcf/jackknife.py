"""The weighted delete-one-block jackknife (Busing, Meijer and van der Leeden 1999), stdlib only.

Busing, F. M. T. A., Meijer, E. and van der Leeden, R. (1999), "Delete-m jackknife for unequal m",
Statistics and Computing 9, 3-8: the jackknife for blocks of unequal size, as ADMIXTOOLS uses it with
chromosomes or runs of sites as blocks (Patterson et al. 2012, Genetics 192, appendix). Kept out of
`just_vcf.ancestry` so that it imports, and is tested, without numpy.
"""

import math
from collections.abc import Sequence

__all__ = ["weighted_block_jackknife"]


def weighted_block_jackknife(full: float, replicates: Sequence[float], sizes: Sequence[int]) -> tuple[float, float]:
    """The bias-corrected estimate and the standard error of a statistic, from its delete-one-block values.

    `full` is the statistic on every block; `replicates[j]` is it with block `j` left out; `sizes[j]` is
    block `j`'s weight, its count of sites (`m_j`). With `n = sum(m_j)`, `g` blocks and `h_j = n / m_j`:

    - estimate `theta_J = g * full - sum_j (1 - m_j / n) * replicates[j]`;
    - pseudo-values `tau_j = h_j * full - (h_j - 1) * replicates[j]`;
    - variance `(1 / g) * sum_j (tau_j - theta_J)^2 / (h_j - 1)`.

    With equal blocks this is the ordinary delete-one jackknife, `(g - 1) / g * sum (replicates[j] -
    mean)^2`. Needs at least two blocks, each with a positive size; anything else is refused, since a
    single block has no variance to estimate (`h_j - 1 = 0`).
    """
    g = len(replicates)
    if g != len(sizes):
        raise ValueError(f"{g} replicates but {len(sizes)} block sizes")
    if g < 2:
        raise ValueError(f"a jackknife needs at least two blocks, got {g}")
    if any(m <= 0 for m in sizes):
        raise ValueError(f"every block needs a positive size, got {list(sizes)}")
    n = float(sum(sizes))
    estimate = g * full - math.fsum((1.0 - m / n) * r for r, m in zip(replicates, sizes, strict=True))
    variance = 0.0
    for r, m in zip(replicates, sizes, strict=True):
        h = n / m
        tau = h * full - (h - 1.0) * r
        variance += (tau - estimate) ** 2 / (h - 1.0)
    return estimate, math.sqrt(variance / g)
