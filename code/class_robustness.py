import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Robustness checks for the class-level (access versus choice) results of class_need.py.
1. Dapagliflozin's share of SGLT2 inhibitors, 2023-2025: need-adjusted Moran's I (Freedman-Lane, 9,999
   permutations) under alternative need models (class-specific, primary, age and deprivation, six
   components of all covariates), k-nearest-neighbour (k=4) weights, and leaving out each area in turn
   (2025, class-specific model, contiguity weights re-standardised after the area is removed).
2. What the class-specific need models rest on: standardised coefficients for SGLT2 total use (2025)
   and for indapamide's share (2021, 2025), and indapamide's share adjusted for age alone or for
   recorded hypertension alone.
3. Multiplicity: Benjamini-Hochberg q-values across the 24 tests of Table 2 (six measures, 2021 and
   2025, unadjusted and adjusted) and across all 30 need-adjusted class-year tests."""
import numpy as np, pandas as pd, pickle, json
from need_adjustment_core import load_covariates, primary_covariates, pcs, resid_maker, bh
W = np.load('W_contig.npy'); Wk = np.load('W_knn4.npy'); Wb = np.load('Wb_contig.npy')
shares = np.load('shares_af.npy'); years = [int(y) for y in np.load('years_af.npy')]
areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb')); codes = pickle.load(open('codes_af.pkl', 'rb'))
assert list(areas_td) == list(codes)
C = load_covariates(codes); n = len(codes)
nm = pd.read_csv('icb_order.csv').set_index('code').ICB23NM.str.replace('NHS ', '').str.replace(' Integrated Care Board', '')
col = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
S = ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin']
qcols = [c for c in C.columns if c.startswith('qof_')]
allc = ['age_0_4', 'age_5_14', 'age_u18', 'age_15plus', 'age_65plus', 'age_75plus', 'age_85plus', 'imd2025'] + qcols
models = {'class-specific': C[['qof_dm', 'qof_hf', 'qof_ckd', 'age_65plus', 'imd2025']].values,
          'primary': primary_covariates(codes),
          'age and deprivation': C[['age_65plus', 'age_u18', 'imd2025']].values,
          'six components of all covariates': pcs(C, allc, 6)[0]}


def fl(x, Xc, Wm, N=9999, seed=0):
    """Need-adjusted Moran's I with Freedman-Lane permutation; I scaled by n/S0 (S0 = sum of weights)."""
    m = len(x); Mres, _ = resid_maker(Xc); e = Mres @ x; f = m / Wm.sum()
    I = f * (e @ Wm @ e) / (e @ e); g = np.random.default_rng(seed)
    U = e[np.array([g.permutation(m) for _ in range(N)])] @ Mres
    Ip = f * ((U @ Wm.T) * U).sum(1) / (U ** 2).sum(1)
    return float(I), float(((Ip >= I - 1e-12).sum() + 1) / (N + 1))


def dshare(y):
    s = shares[years.index(y)]; return s[:, col['dapagliflozin']] / sum(s[:, col[d]] for d in S)


out = {'dapagliflozin_share': {}}
for y in [2023, 2024, 2025]:
    x = dshare(y); r = {}
    for k, (nmod, Xc) in enumerate(models.items()):
        r[f'{nmod}, contiguity'] = fl(x, Xc, W, seed=10 * y + k)
    r['class-specific, k-nearest neighbours'] = fl(x, models['class-specific'], Wk, seed=10 * y + 9)
    out['dapagliflozin_share'][y] = r
    print(y, {k: (round(v[0], 3), v[1]) for k, v in r.items()}, flush=True)
loo = []
x = dshare(2025)
for i in range(n):
    keep = np.delete(np.arange(n), i); B = Wb[np.ix_(keep, keep)]; d = B.sum(1)
    Wr = np.divide(B, d[:, None], out=np.zeros_like(B), where=d[:, None] > 0)
    I_, p_ = fl(x[keep], models['class-specific'][keep], Wr, N=9999, seed=500 + i)
    loo.append(dict(dropped=nm[codes[i]], I=I_, p=p_, areas_without_neighbours=int((d == 0).sum())))
L = pd.DataFrame(loo)
out['dapagliflozin_share']['leave_one_out_2025'] = dict(I_min=float(L.I.min()), I_max=float(L.I.max()), p_max=float(L.p.max()),
                                                        n_p_ge_005=int((L.p >= 0.05).sum()), weakest=L.sort_values('p').tail(3).to_dict('records'))
print('LOO', json.dumps(out['dapagliflozin_share']['leave_one_out_2025'], default=float), flush=True)


def std_coefs(yv, cols):
    """Standardised OLS coefficients with their t statistics and two-sided P-values."""
    from scipy.stats import t as tdist
    Z = (C[cols] - C[cols].mean()) / C[cols].std(ddof=0); X = np.column_stack([np.ones(n), Z.values])
    yz = (yv - yv.mean()) / yv.std(ddof=0); b = np.linalg.lstsq(X, yz, rcond=None)[0]
    res = yz - X @ b; df = n - X.shape[1]; s2 = (res @ res) / df
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return {c: dict(beta=round(float(b[i + 1]), 3), t=round(float(b[i + 1] / se[i + 1]), 2),
                    p=round(float(2 * tdist.sf(abs(b[i + 1] / se[i + 1]), df)), 4)) for i, c in enumerate(cols)}


s25 = shares[years.index(2025)]; tS = sum(s25[:, col[d]] for d in S)
out['sglt2_total_2025'] = dict(standardised_coefficients=std_coefs(tS, ['qof_dm', 'qof_hf', 'qof_ckd', 'age_65plus', 'imd2025']),
                               correlations={c: round(float(np.corrcoef(tS, C[c])[0, 1]), 3) for c in ['qof_dm', 'qof_hf', 'qof_ckd', 'age_65plus', 'imd2025']})
out['indapamide_share'] = {}
for y in [2021, 2025]:
    s = shares[years.index(y)]; ind = s[:, col['indapamide']] / (s[:, col['indapamide']] + s[:, col['bendroflumethiazide']])
    out['indapamide_share'][y] = dict(standardised_coefficients=std_coefs(ind, ['qof_hyp', 'age_65plus', 'imd2025']),
                                      adjusted_for_age_only=fl(ind, C[['age_65plus']].values, W, seed=y + 1),
                                      adjusted_for_hypertension_only=fl(ind, C[['qof_hyp']].values, W, seed=y + 2),
                                      class_model=fl(ind, C[['qof_hyp', 'age_65plus', 'imd2025']].values, W, seed=y + 3),
                                      r_hypertension_age65=round(float(np.corrcoef(C.qof_hyp, C.age_65plus)[0, 1]), 3))
    print(y, 'indapamide', json.dumps(out['indapamide_share'][y]), flush=True)
cn = json.load(open('class_need.json'))['results']
t24 = [(m, y, a) for m in cn for y in ['2021', '2025'] for a in ['p', 'p_adj']]
q24 = bh(np.array([cn[m][y][a] for m, y, a in t24]))
t30 = [(m, y) for m in cn for y in cn[m]]
q30 = bh(np.array([cn[m][y]['p_adj'] for m, y in t30]))
out['multiplicity'] = dict(table2_q={f'{m} | {y} | {"adjusted" if a == "p_adj" else "unadjusted"}': round(float(q), 4) for (m, y, a), q in zip(t24, q24)},
                           adjusted_class_year_q={f'{m} | {y}': round(float(q), 4) for (m, y), q in zip(t30, q30)})
json.dump(out, open('class_robustness.json', 'w'), indent=1, default=float)
print(json.dumps(out['multiplicity'], indent=1))
print(json.dumps(out['sglt2_total_2025'], indent=1))
