import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Need-adjusted spatial clustering.
Each item's within-area share is regressed on ICB-level need covariates (OLS, 42 areas) and global
Moran's I is computed on the residuals. Significance uses the Freedman-Lane permutation scheme for
regression residuals: residuals of the covariate-only model are permuted and re-projected through
the residual-maker matrix M = I - X(X'X)^-1 X', so every permuted statistic is a valid residual
Moran's I. One-sided test for positive autocorrelation; 9,999 permutations, refined to 99,999 for
P<0.05 in every model (the k-nearest-neighbour sensitivity analysis reports nominal counts from 9,999); P=(b+1)/(N+1) with ties within 1e-12 counted; Benjamini-Hochberg FDR across all items;
clustered = positive residual I and q<0.05.
Primary model: % aged 65+, % aged under 18, IMD 2025 score, and the first two principal components
of 21 QOF disease prevalences (morbidity). Sensitivity: (A) age and deprivation only; (B) first six
principal components of all 29 covariates; (C) primary model with k-nearest-neighbour (k=4) weights."""
import numpy as np, pandas as pd, pickle, json
from need_adjustment_core import load_covariates, pcs, primary_covariates, resid_maker, moran_rows, fl_perm, bh
POS = 1e-12
A = np.load('A_af.npy'); W = np.load('W_contig.npy'); Wk = np.load('W_knn4.npy')
codes = pickle.load(open('codes_af.pkl', 'rb')); M0 = pd.read_csv('moran_af.csv')
C = load_covariates(codes)
n = len(codes); Y = A.T.astype(float)                       # items x areas
qcols = [c for c in C.columns if c.startswith('qof_')]
allc = ['age_0_4', 'age_5_14', 'age_u18', 'age_15plus', 'age_65plus', 'age_75plus', 'age_85plus', 'imd2025'] + qcols
_, qev = pcs(C, qcols, 2); apc, aev = pcs(C, allc, 6)
models = {
    'primary': primary_covariates(codes),
    'age_deprivation': np.column_stack([C.age_65plus, C.age_u18, C.imd2025]),
    'pca6': apc,
}
out = {'qof_pc_explained': qev.tolist(), 'all_pc_explained': aev.tolist(), 'n_items': int(len(Y)),
       'corr_imd2019_imd2025_29_icbs': float(C[['imd2019', 'imd2025']].dropna().corr().iloc[0, 1])}
res = pd.DataFrame({'Drug': M0.Drug, 'moran_I_raw': M0.moran_I, 'q_raw': M0.q_moran, 'clustered_raw': M0.clustered})
for name, Xc in models.items():
    Mres, k = resid_maker(Xc)
    E = Y @ Mres; sst = ((Y - Y.mean(1, keepdims=True)) ** 2).sum(1); r2 = 1 - (E ** 2).sum(1) / sst
    I = moran_rows(E, W); EI = np.trace(Mres @ W) / (n - k)
    p = fl_perm(E, Mres, W, I, 9999, seed=11)
    cand = np.where(p < 0.05)[0]                                # refine every model's candidates
    p[cand] = fl_perm(E[cand], Mres, W, I[cand], 99999, seed={'primary': 12, 'age_deprivation': 14, 'pca6': 15}[name], B=200)
    q = bh(p); clu = (q < 0.05) & (I > POS); nom = (p < 0.05) & (I > POS)
    res[f'R2_{name}'] = r2; res[f'I_{name}'] = I; res[f'p_{name}'] = p; res[f'q_{name}'] = q
    res[f'nominal_{name}'] = nom; res[f'clustered_{name}'] = clu
    out[name] = dict(k_covariates=k - 1, expected_I=float(EI), median_R2=float(np.median(r2)), median_I=float(np.median(I)),
                     n_nominal=int(nom.sum()), n_clustered=int(clu.sum()))
    if name == 'primary':
        Ik = moran_rows(E, Wk); pk = fl_perm(E, Mres, Wk, Ik, 9999, seed=13)
        res['I_primary_knn4'] = Ik; res['p_primary_knn4'] = pk
        out[name]['n_nominal_knn4'] = int(((pk < 0.05) & (Ik > POS)).sum())
        lam = 0.5; out[name]['storey_pi0'] = float(min(1, (p > lam).mean() / (1 - lam)))
    print(name, out[name])
res.to_csv('need_adjusted_moran.csv', index=False)
json.dump(out, open('need_summary.json', 'w'), indent=2)
print('saved need_adjusted_moran.csv, need_summary.json')
