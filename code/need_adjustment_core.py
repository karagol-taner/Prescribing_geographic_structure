"""Shared functions for the need-adjusted spatial analyses (residual Moran's I with Freedman-Lane
permutation). W is passed explicitly; areas are in the analysis order of codes_af.pkl."""
import numpy as np, pandas as pd
TOL = 1e-12


def load_covariates(codes, path='icb_covariates.csv'):
    return pd.read_csv(path).set_index('ods').loc[codes]


def pcs(C, cols, k):
    """First k principal components of the standardised columns (population SD)."""
    Z = (C[cols] - C[cols].mean()) / C[cols].std(ddof=0)
    u, s, vt = np.linalg.svd(Z.values, full_matrices=False)
    return Z.values @ vt[:k].T, (s ** 2 / (s ** 2).sum())[:k]


def primary_covariates(codes, path='icb_covariates.csv'):
    """Primary need model: % aged 65+, % aged under 18, IMD 2025 score, first two PCs of 21 QOF prevalences."""
    C = load_covariates(codes, path)
    q, _ = pcs(C, [c for c in C.columns if c.startswith('qof_')], 2)
    return np.column_stack([C.age_65plus, C.age_u18, C.imd2025, q])


def resid_maker(Xc):
    """Residual-maker matrix M = I - X(X'X)^-1 X' for covariates Xc (intercept added) and the number of columns of X."""
    n = len(Xc); X = np.column_stack([np.ones(n), Xc])
    return np.eye(n) - X @ np.linalg.pinv(X.T @ X) @ X.T, X.shape[1]


def moran_rows(R, Wm):
    """Moran's I (row-standardised W) of each row of R, which must already be centred or residualised."""
    return (R * (R @ Wm.T)).sum(-1) / (R ** 2).sum(-1)


def fl_perm(E, Mres, Wm, I_obs, N, seed, B=100):
    """One-sided Freedman-Lane permutation P-values: rows of E are residual vectors; each permutation
    permutes the residuals and re-projects them through M; P=(b+1)/(N+1), ties within 1e-12 counted."""
    n = E.shape[1]; g = np.random.default_rng(seed); c = np.zeros(len(E)); done = 0
    while done < N:
        b = min(B, N - done); P = np.array([g.permutation(n) for _ in range(b)])
        R = E[:, P] @ Mres
        c += (moran_rows(R, Wm) >= I_obs[:, None] - TOL).sum(1); done += b
    return (c + 1) / (N + 1)


def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1)
