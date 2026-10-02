import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Access versus choice with need adjustment. For each class, total use (the class's within-area share)
and choice of agent (one agent's share of the class) are tested with Moran's I by calendar year, both
unadjusted and on residuals from a class-specific need model (Freedman-Lane permutation, 9,999,
fixed seeds). Need models: SGLT2 inhibitors - QOF diabetes, heart failure and CKD prevalence, % aged
65+, IMD 2025; GLP-1 and GIP/GLP-1 agonists - QOF diabetes and obesity prevalence, % aged 65+, IMD 2025;
thiazide-type diuretics - QOF hypertension prevalence, % aged 65+, IMD 2025."""
import numpy as np, pandas as pd, pickle, json
W = np.load('W_contig.npy'); shares = np.load('shares_af.npy'); years = np.load('years_af.npy')
areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb')); codes = pickle.load(open('codes_af.pkl', 'rb'))
assert list(areas_td) == list(codes)
C = pd.read_csv('icb_covariates.csv').set_index('ods').loc[codes]
n = len(codes); yi = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
S = ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin']
GL = ['semaglutide', 'dulaglutide', 'liraglutide', 'exenatide', 'lixisenatide', 'tirzepatide']
TH = ['bendroflumethiazide', 'indapamide', 'chlortalidone', 'hydrochlorothiazide', 'metolazone']
need = {'SGLT2': ['qof_dm', 'qof_hf', 'qof_ckd', 'age_65plus', 'imd2025'],
        'GLP1': ['qof_dm', 'qof_obes', 'age_65plus', 'imd2025'],
        'THZ': ['qof_hyp', 'age_65plus', 'imd2025']}


def resid_maker(cols):
    X = np.column_stack([np.ones(n), C[cols].values]); return np.eye(n) - X @ np.linalg.pinv(X.T @ X) @ X.T


def moran_test(x, Mres=None, N=9999, seed=0):
    z = x - x.mean() if Mres is None else Mres @ x
    I = z @ W @ z / (z @ z)
    g = np.random.default_rng(seed); P = np.array([g.permutation(n) for _ in range(N)])
    Zp = z[P] if Mres is None else z[P] @ Mres            # Freedman-Lane: permute residuals, re-project
    Ip = ((Zp @ W.T) * Zp).sum(1) / (Zp ** 2).sum(1)
    return float(I), float(((Ip >= I - 1e-12).sum() + 1) / (N + 1))


def series(k):
    s = shares[k]; col = lambda m: s[:, yi[m]]
    tS = sum(col(m) for m in S); tG = sum(col(m) for m in GL); tT = sum(col(m) for m in TH)
    return {'SGLT2 inhibitors: class total': (tS, 'SGLT2'), 'Dapagliflozin share of SGLT2 inhibitors': (col('dapagliflozin') / tS, 'SGLT2'),
            'GLP-1 and GIP/GLP-1 agonists: class total': (tG, 'GLP1'), 'Semaglutide share of GLP-1 class': (col('semaglutide') / tG, 'GLP1'),
            'Thiazide-type diuretics: class total': (tT, 'THZ'),
            'Indapamide share of indapamide plus bendroflumethiazide': (col('indapamide') / (col('indapamide') + col('bendroflumethiazide')), 'THZ')}


Ms = {k: resid_maker(v) for k, v in need.items()}
out = {}
for k, y in enumerate(years):
    for name, (v, cls) in series(k).items():
        I0, p0 = moran_test(v, None, seed=500 + k)              # same seeds as the unadjusted analysis in make_geography_figures.py
        Ia, pa = moran_test(v, Ms[cls], seed=700 + k)
        X = np.column_stack([np.ones(n), C[need[cls]].values]); r2 = 1 - ((Ms[cls] @ v) ** 2).sum() / ((v - v.mean()) ** 2).sum()
        out.setdefault(name, {})[int(y)] = dict(I=I0, p=p0, I_adj=Ia, p_adj=pa, R2_need=float(r2))
json.dump(dict(need_models=need, results=out), open('class_need.json', 'w'), indent=1)
for name, d in out.items():
    print(name)
    for y, r in d.items():
        print(f"  {y}: I {r['I']:+.3f} (P={r['p']:.4f})  need-adjusted I {r['I_adj']:+.3f} (P={r['p_adj']:.4f})  R2 {r['R2_need']:.2f}")
