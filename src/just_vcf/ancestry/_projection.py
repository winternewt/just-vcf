"""Reference-PCA build and shrinkage-corrected projection (vendored from FRAPOSA, via just-prs).

The functions ``standardize``, ``svd_eigcov``, ``svd_online``, ``procrustes``,
``procrustes_diffdim``, ``ref_aug_procrustes`` and ``oadp`` are adapted (numpy-only,
array-first, logging and file I/O removed) from FRAPOSA / fraposa_pgsc:

    https://github.com/PGScatalog/fraposa_pgsc  (MIT License)
    Copyright (c) Daiwei Zhang (2023) and the PGS Catalog team.

The full upstream MIT notice ships alongside this module as ``FRAPOSA_LICENSE``. just-prs vendored
them first (``just_prs/ancestry/_projection.py``, MIT); this copy is that file with this docstring
rewritten and just-vcf's changes (docs/PROVENANCE.md): ``svd_online`` reads a residual below
``RESIDUAL_TOLERANCE`` as none instead of dividing it by its norm, and sums it with ``np.sum``; its update
proper is ``_augment``, which ``svd_online_from_projection`` (just-vcf's) reaches from ``U1^T b`` and
``|b|^2`` without U1.

Why vendor instead of depend: ``fraposa-pgsc`` pins ``numpy<2.0`` and ``pandas<2.0``.

OADP (online augmentation, decomposition and Procrustes) corrects the shrinkage bias of naive PCA
projection when the variants far outnumber the reference samples, which is why FRAPOSA exists. It is
pure numpy, with no compiled extension.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "standardize",
    "svd_eigcov",
    "svd_online",
    "svd_online_from_projection",
    "procrustes",
    "procrustes_diffdim",
    "ref_aug_procrustes",
    "oadp",
    "build_reference_pca",
    "ReferencePCA",
    "project_samples",
    "DEFAULT_MISSING",
    "RESIDUAL_TOLERANCE",
]

# Missing-genotype sentinel used in standardisation. A panel matrix and a sample vector carry -9 for
# a missing genotype; standardisation treats it as missing (set to the column mean,
# i.e. 0 after centering).
DEFAULT_MISSING: float = -9.0

RESIDUAL_TOLERANCE: float = float(np.sqrt(np.finfo(np.float64).eps))
"""just-vcf's, not FRAPOSA's: the largest residual `b - U1 U1^T b`, relative to `|b|`, that
`svd_online` and `svd_online_from_projection` read as `b` lying in U1's span (about 1.5e-8). With
U1's columns orthonormal the residual is orthogonal to them, so its only trace in the update is
`b_tilde^T b = |residual|`, one entry of the small matrix whose SVD is taken: zeroing it moves that
matrix by at most this share of `|b|`, below any digit a PC coordinate is read to. And it is far above
what rounding leaves of a vector in the span (about `p` float64 epsilons, 1e-14 at the few dozen sites
where the span can be all of R^p), while a target over a real panel (sites far more than U1's columns)
keeps most of its norm outside the span."""


def standardize(
    X: np.ndarray,
    mean: np.ndarray | None = None,
    std: np.ndarray | None = None,
    miss: float = DEFAULT_MISSING,
) -> tuple[np.ndarray, np.ndarray]:
    """In-place per-variant standardisation of a (variants x samples) float matrix.

    Missing entries (== ``miss``) are excluded from the mean/std and set to 0 after
    centering. When ``mean``/``std`` are given they are reused (study-sample path);
    otherwise they are computed (reference-build path). Returns ``(mean, std)`` as
    column vectors.
    """
    assert np.issubdtype(X.dtype, np.floating)
    p, _ = X.shape
    is_miss = X == miss
    if (mean is None) or (std is None):
        mean = np.zeros(p)
        std = np.zeros(p)
        for i in range(p):
            row_nomiss = X[i, :][~is_miss[i, :]]
            if row_nomiss.size:
                mean[i] = np.mean(row_nomiss)
                std[i] = np.std(row_nomiss)
        std[std == 0] = 1
    mean = mean.reshape((-1, 1))
    std = std.reshape((-1, 1))
    X -= mean
    X /= std
    X[is_miss] = 0
    return mean, std


def svd_eigcov(XTX: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """SVD of a covariance matrix via eigendecomposition. Returns ``(s, V)`` desc."""
    ssq, V = np.linalg.eigh(XTX)
    V = np.squeeze(V)
    ssq = np.squeeze(ssq)
    s = np.sqrt(abs(ssq))
    V = V.T[::-1].T
    s = s[::-1]
    return s, V


def svd_online(
    U1: np.ndarray, d1: np.ndarray, V1: np.ndarray, b: np.ndarray, l: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Rank-1 online SVD update augmenting the reference SVD with one new sample ``b``.

    The update reads ``b`` only through its coordinates ``U1^T b`` and one number, ``b_tilde^T b``
    with ``b_tilde`` the unit residual off U1's span (``_augment``); here both come from U1 itself.
    ``svd_online_from_projection`` is the same update from the coordinates and ``|b|^2`` alone.
    """
    assert U1.shape[1] == V1.shape[1]
    p = U1.shape[0]
    b = b.reshape((p, 1))
    coordinates = U1.T @ b
    b_tilde = b - U1 @ coordinates
    residual = np.sqrt(np.sum(np.square(b_tilde)))
    # just-vcf (audit 16): a residual of rounding size means b lies in U1's span (always so when U1 has
    # p <= k orthonormal columns, a chromosome jackknife replicate keeping few sites); the update then
    # adds no direction, so the residual's entry is zero instead of 0 / 0 or a unit vector of noise.
    # Otherwise the normalisation is FRAPOSA's, with `np.sum` for its builtin `sum` (RM306).
    if residual > RESIDUAL_TOLERANCE * np.sqrt(np.sum(np.square(b))):
        b_tilde = b_tilde / residual
        tail = float((b_tilde.transpose() @ b)[0, 0])
    else:
        tail = 0.0
    return _augment(d1, V1, coordinates, tail, l)


def svd_online_from_projection(
    d1: np.ndarray, V1: np.ndarray, coordinates: np.ndarray, b_norm_sq: float, l: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """``svd_online``'s update from ``b``'s coordinates ``c = U1^T b`` and ``|b|^2``, U1 never formed.

    just-vcf's, not FRAPOSA's (RM306). With U1's columns orthonormal the residual ``b - U1 c`` is
    orthogonal to them, so ``b_tilde^T b`` is the residual's norm ``r`` and ``r^2 = |b|^2 - |c|^2``.
    A caller that has ``c`` and ``|b|^2`` from sums it already holds (a jackknife replicate, from the
    panel's sums less one chromosome's) never touches a row of U1.

    The same tolerance as ``svd_online``: ``r^2 <= (RESIDUAL_TOLERANCE |b|)^2`` reads as ``b`` in
    U1's span and the residual's entry is zero; a negative ``r^2`` or ``|b|^2`` (rounding) falls under
    it. The subtraction is safe there because the residual enters the update only as that one entry of
    the small matrix whose SVD is taken, never as a direction, and that SVD moves continuously with
    it. Its rounding error is a few float64 epsilons of ``|b|^2`` (and of the sums it came from): far
    from the tolerance, as every replicate the jackknife runs is (more kept sites than loadings), ``r``
    keeps every digit a coordinate is read to; near it, ``r`` is uncertain by about
    ``RESIDUAL_TOLERANCE |b|``, which moves the update no more than the tolerance's own zeroing does.
    """
    coordinates = np.asarray(coordinates, dtype=np.float64).reshape((-1, 1))
    residual_sq = float(b_norm_sq) - float(np.sum(np.square(coordinates)))
    if residual_sq > RESIDUAL_TOLERANCE**2 * max(float(b_norm_sq), 0.0):
        tail = float(np.sqrt(residual_sq))
    else:
        tail = 0.0
    return _augment(d1, V1, coordinates, tail, l)


def _augment(
    d1: np.ndarray, V1: np.ndarray, coordinates: np.ndarray, tail: float, l: int | None
) -> tuple[np.ndarray, np.ndarray]:
    """FRAPOSA's update proper: the SVD of ``[[diag(d1), U1^T b], [0, b_tilde^T b]]`` rotating
    ``[[V1^T, 0], [0, 1]]``, with ``coordinates = U1^T b`` and ``tail = b_tilde^T b``."""
    n, k = V1.shape
    if l is None:
        l = k
    assert len(d1) == k
    assert coordinates.shape == (k, 1)
    assert l <= k
    R = np.concatenate((np.diag(d1), coordinates), axis=1)
    R_tail = np.concatenate((np.zeros((1, k)), np.array([[tail]])), axis=1)
    R = np.concatenate((R, R_tail), axis=0)
    d2, R_Vt = np.linalg.svd(R, full_matrices=False)[1:]
    V_new = np.zeros((k + 1, n + 1))
    V_new[:k, :n] = V1.transpose()
    V_new[k, n] = 1
    V2 = (R_Vt @ V_new).transpose()[:, :l]
    return d2, V2


def procrustes(
    Y_mat: np.ndarray, X_mat: np.ndarray, return_transformed: bool = False
):
    """Best similarity transform (rotation R, scale rho, translation c) from X to Y."""
    X = np.array(X_mat, dtype=np.double, copy=True)
    Y = np.array(Y_mat, dtype=np.double, copy=True)
    X_mean = np.mean(X, 0)
    Y_mean = np.mean(Y, 0)
    X -= X_mean
    Y -= Y_mean
    C = Y.T @ X
    U, s, VT = np.linalg.svd(C, full_matrices=False)
    trXX = np.sum(X**2)
    trS = np.sum(s)
    R = VT.T @ U.T
    rho = trS / trXX
    c = Y_mean - rho * X_mean @ R
    if return_transformed:
        X_new = X_mat @ R * rho + c
        return R, rho, c, X_new
    return R, rho, c


def procrustes_diffdim(
    Y_mat: np.ndarray,
    X_mat: np.ndarray,
    n_iter_max: int = 10000,
    epsilon_min: float = 1e-6,
    return_transformed: bool = False,
):
    """Procrustes when X has more dimensions than Y (iterative padding)."""
    X = np.array(X_mat, dtype=np.double, copy=True)
    Y = np.array(Y_mat, dtype=np.double, copy=True)
    n_X, p_X = X.shape
    n_Y, p_Y = Y.shape
    assert n_X == n_Y
    assert p_X >= p_Y
    if p_X == p_Y:
        return procrustes(Y, X, return_transformed)
    Z = np.zeros((n_X, p_X - p_Y))
    R = rho = c = X_new = None
    for _ in range(n_iter_max):
        W = np.hstack((Y, Z))
        R, rho, c = procrustes(W, X)
        X_new = X @ R * rho + c
        Z_new = X_new[:, p_Y:]
        Z_new_centered = Z_new - np.mean(Z_new, 0)
        Z_diff = Z_new - Z
        epsilon = np.sum(Z_diff**2) / np.sum(Z_new_centered**2)
        if epsilon < epsilon_min:
            break
        Z = Z_new
    if return_transformed:
        return R, rho, c, X_new
    return R, rho, c


def ref_aug_procrustes(pcs_ref: np.ndarray, pcs_aug: np.ndarray) -> np.ndarray:
    """Map the augmented (ref + 1 study) PCs back onto the reference PC frame."""
    n_ref, _ = pcs_ref.shape
    n_aug, _ = pcs_aug.shape
    assert n_aug == n_ref + 1
    pcs_aug_head = pcs_aug[:-1, :]
    pcs_aug_tail = pcs_aug[-1, :].reshape((1, -1))
    R, rho, c = procrustes_diffdim(pcs_ref, pcs_aug_head)
    pcs_aug_tail_trsfed = pcs_aug_tail @ R * rho + c
    return pcs_aug_tail_trsfed.flatten()


def oadp(
    U: np.ndarray,
    s: np.ndarray,
    V: np.ndarray,
    b: np.ndarray,
    dim_ref: int = 4,
    dim_stu: int | None = None,
    dim_online: int | None = None,
) -> np.ndarray:
    """Project a single standardised study sample ``b`` with shrinkage correction."""
    if dim_stu is None:
        dim_stu = dim_ref * 2
    if dim_online is None:
        dim_online = dim_stu * 2
    pcs_ref = V[:, :dim_ref] * s[:dim_ref]
    s_aug, V_aug = svd_online(U[:, :dim_online], s[:dim_online], V[:, :dim_online], b)
    s_aug, V_aug = s_aug[:dim_stu], V_aug[:, :dim_stu]
    pcs_aug = V_aug * s_aug
    pcs_stu = ref_aug_procrustes(pcs_ref, pcs_aug)
    return pcs_stu[:dim_ref]


# ---------------------------------------------------------------------------
# just-prs wrappers (not from FRAPOSA): reference build + batched projection.
# ---------------------------------------------------------------------------


class ReferencePCA:
    """Reference PCA model produced by :func:`build_reference_pca`.

    Holds everything OADP projection + KNN classification need:
    ``mean``/``std`` (per variant), ``s`` (singular values), ``V`` (n_samples x
    dim_online sample eigenvectors), ``U`` (n_variants x dim_online loadings), and the
    derived ``pcs_ref`` (n_samples x dim_ref reference PC scores).
    """

    __slots__ = ("mean", "std", "s", "V", "U", "dim_ref", "dim_stu", "dim_online")

    def __init__(self, mean, std, s, V, U, dim_ref, dim_stu, dim_online):
        self.mean = mean
        self.std = std
        self.s = s
        self.V = V
        self.U = U
        self.dim_ref = dim_ref
        self.dim_stu = dim_stu
        self.dim_online = dim_online

    @property
    def pcs_ref(self) -> np.ndarray:
        """Reference-sample PC scores (n_samples x dim_ref)."""
        return self.V[:, : self.dim_ref] * self.s[: self.dim_ref]


def build_reference_pca(
    X: np.ndarray, dim_ref: int = 10, missing: float = DEFAULT_MISSING
) -> ReferencePCA:
    """Build a reference PCA from a (variants x samples) genotype dosage matrix.

    Mirrors FRAPOSA's ``pca`` reference-build path: standardise per variant, eigen-
    decompose the sample covariance ``X.T @ X``, then form loadings ``U = X (V/s)``.
    ``dim_stu = 2*dim_ref``, ``dim_online = 2*dim_stu`` (FRAPOSA defaults).
    """
    dim_stu = dim_ref * 2
    dim_online = dim_stu * 2
    Xf = X.astype(np.float64, copy=True)
    mean, std = standardize(Xf, miss=missing)
    n_samples = Xf.shape[1]
    dim_online = min(dim_online, n_samples)
    XTX = Xf.T @ Xf
    s, V = svd_eigcov(XTX)
    V = V[:, :dim_online]
    s = s[:dim_online]
    # Guard against zero singular values when forming loadings.
    s_safe = np.where(s == 0, 1.0, s)
    U = Xf @ (V / s_safe)
    return ReferencePCA(mean, std, s, V, U, dim_ref, dim_stu, dim_online)


def project_samples(model: ReferencePCA, W: np.ndarray, missing: float = DEFAULT_MISSING) -> np.ndarray:
    """OADP-project standardised study genotypes ``W`` (variants x samples) -> PC scores.

    Returns an (n_samples x dim_ref) array. ``W`` is standardised in place against the
    reference ``mean``/``std`` (missing -> column mean) before projection.
    """
    p_stu, n_stu = W.shape
    pcs = np.zeros((n_stu, model.dim_ref))
    for i in range(n_stu):
        w = W[:, i].astype(np.float64).reshape((-1, 1))
        standardize(w, model.mean, model.std, miss=missing)
        pcs[i, :] = oadp(
            model.U, model.s, model.V, w,
            dim_ref=model.dim_ref, dim_stu=model.dim_stu, dim_online=model.dim_online,
        )
    return pcs
