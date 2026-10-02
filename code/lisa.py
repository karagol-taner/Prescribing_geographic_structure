import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Local indicators of spatial association (local Moran's I, Anselin 1995) for SGLT2 inhibitor use and
dapagliflozin's share of SGLT2 inhibitors in 2025, before and after need adjustment (class-specific
need model of class_need.py: QOF diabetes, heart failure and CKD prevalence, % aged 65+, IMD 2025).
I_i = z_i * sum_j w_ij z_j with z standardised (residuals standardised for the adjusted versions) and
row-standardised contiguity weights. Conditional permutation: z_i is held fixed and the other 41 values
are permuted (9,999 permutations); folded one-sided pseudo P = (b+1)/(N+1) in the direction of the
observed I_i. Locations with P<0.05 are classed High-High, Low-Low, High-Low or Low-High; Benjamini-Hochberg
q-values across the 42 areas are reported to flag locations that survive control of the false discovery rate."""
import numpy as np, pandas as pd, pickle, json
W = np.load('W_contig.npy'); shares = np.load('shares_af.npy'); years = np.load('years_af.npy')
areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb')); codes = pickle.load(open('codes_af.pkl', 'rb'))
C = pd.read_csv('icb_covariates.csv').set_index('ods').loc[codes]
n = len(codes); yi = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
k = list(years).index(2025); s = shares[k]
S = ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin']
tot = sum(s[:, yi[m]] for m in S); dshare = s[:, yi['dapagliflozin']] / tot
X = np.column_stack([np.ones(n), C[['qof_dm', 'qof_hf', 'qof_ckd', 'age_65plus', 'imd2025']].values])
Mres = np.eye(n) - X @ np.linalg.pinv(X.T @ X) @ X.T


def lisa(x, N=9999, seed=0):
    z = (x - x.mean()) / x.std(); lag = W @ z; I = z * lag
    g = np.random.default_rng(seed); p = np.empty(n)
    for i in range(n):
        others = np.delete(np.arange(n), i); wi = np.delete(W[i], i)
        Pm = np.array([g.permutation(others) for _ in range(N)])
        Ii = z[i] * (z[Pm] @ wi)
        b = (Ii >= I[i] - 1e-12).sum() if I[i] >= 0 else (Ii <= I[i] + 1e-12).sum()
        p[i] = (b + 1) / (N + 1)
    cat = np.where(p >= 0.05, 'Not significant', np.where(z > 0, np.where(lag > 0, 'High-High', 'High-Low'), np.where(lag < 0, 'Low-Low', 'Low-High')))
    o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return z, lag, I, p, np.minimum(q, 1), cat


rows = []; summary = {}
for name, x in [('sglt2_use', tot), ('sglt2_use_need_adjusted', Mres @ tot), ('dapagliflozin_share', dshare), ('dapagliflozin_share_need_adjusted', Mres @ dshare)]:
    z, lag, I, p, q, cat = lisa(x, seed=sum(map(ord, name)))
    summary[name] = {c: int((cat == c).sum()) for c in ['High-High', 'Low-Low', 'High-Low', 'Low-High']}
    for i, code in enumerate(codes):
        rows.append(dict(variable=name, ods=code, value=float(x[i]), z=float(z[i]), lag=float(lag[i]), local_I=float(I[i]), p=float(p[i]), q=float(q[i]), category=cat[i]))
L = pd.DataFrame(rows); L.to_csv('lisa_2025.csv', index=False)
order = pd.read_csv('icb_order.csv').set_index('code'); nm = order.ICB23NM.str.replace('NHS ', '').str.replace(' Integrated Care Board', '')
for v, d in summary.items():
    print(v, d)
    sub = L[(L.variable == v) & (L.category != 'Not significant')]
    for _, r in sub.iterrows(): print(f'   {r.category:10s} {nm[r.ods]} (value {r.value:.4f}, P={r.p:.4f}, q={r.q:.3f})')
json.dump(summary, open('lisa_summary.json', 'w'), indent=1)
