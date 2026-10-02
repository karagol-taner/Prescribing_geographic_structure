"""Jensen-Shannon distance between area profiles and the per-item Mantel test (phylogeographic signal).
JSD: smooth profiles by +1e-10, renormalise rows, JS divergence with natural-log
KL (scipy.stats.entropy default base), distance = sqrt(JS divergence).
Mantel: per drug, pairwise |x_i-x_j| on the area proportion vector, Pearson
correlation against the off-diagonal of the JSD matrix; significance by
area-relabelling permutation, p = max(#(r_perm>=r_obs)/n_perm, 1/n_perm).
"""
import numpy as np
from scipy.stats import entropy


def jsd_matrix(A, eps=1e-10):
    """A: (n_areas, n_drugs) non-negative; rows are renormalised internally."""
    P = A.astype(float) + eps
    P = P / P.sum(1, keepdims=True)
    n = P.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p, q = P[i], P[j]
            m = 0.5 * (p + q)
            js = 0.5 * entropy(p, m) + 0.5 * entropy(q, m)   # natural log
            D[i, j] = D[j, i] = np.sqrt(js)
    return D


def _rvec(X, iu, dc, dss):
    M = np.abs(X[:, iu[0]] - X[:, iu[1]])
    Mc = M - M.mean(1, keepdims=True)
    num = Mc @ dc
    den = np.sqrt((Mc ** 2).sum(1) * dss)
    with np.errstate(invalid='ignore', divide='ignore'):
        r = num / den
    return np.where(np.isfinite(r), r, 0.0)


def mantel_all(A, D, n_perm=1000, seed=0):
    """Per-drug Mantel r and permutation p. A: (n_areas, n_drugs) proportions."""
    n = A.shape[0]
    iu = np.triu_indices(n, k=1)
    dvec = D[iu]
    dc = dvec - dvec.mean()
    dss = (dc ** 2).sum()
    X = A.T.astype(float)                      # (n_drugs, n_areas)
    r_obs = _rvec(X, iu, dc, dss)
    rng = np.random.default_rng(seed)
    count = np.zeros(len(r_obs))
    for _ in range(n_perm):
        perm = rng.permutation(n)
        count += (_rvec(X[:, perm], iu, dc, dss) >= r_obs)
    p = np.maximum(count / n_perm, 1.0 / n_perm)
    return r_obs, p
