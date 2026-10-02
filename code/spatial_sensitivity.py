import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Step 8: sensitivity analyses of the unadjusted spatial clustering.
(1) Year-by-year spatial clustering (Moran's I per calendar year): formulary-wide
    and for named agents, to test whether new drugs spread as regional waves.
(2) Size of variation: P90/P10 ratio of within-area share across the 42 areas.
(3) Coding-stable sensitivity: restrict to items dispensed in every year 2021-2025."""
import numpy as np, pandas as pd, pickle, json
from method import jsd_matrix, mantel_all
W = np.load('W_contig.npy'); A = np.load('A_af.npy'); R = pd.read_csv('moran_af.csv')
shares = np.load('shares_af.npy'); years = np.load('years_af.npy')
areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb')); codes = pickle.load(open('codes_af.pkl', 'rb'))
assert list(areas_td) == list(codes), 'area order mismatch'
n = len(codes); out = {}
def moran_all(Xd, N, seed):
    Z = Xd - Xd.mean(1, keepdims=True); ss = (Z ** 2).sum(1); ok = ss > 0
    I = np.full(len(Xd), np.nan); p = np.full(len(Xd), np.nan)
    Zo = Z[ok]; sso = ss[ok]; I[ok] = (Zo * (Zo @ W.T)).sum(1) / sso
    g = np.random.default_rng(seed); c = np.zeros(ok.sum()); done = 0
    while done < N:
        b = min(100, N - done); P = np.array([g.permutation(n) for _ in range(b)])
        Zp = Zo[:, P]; c += ((Zp * (Zp @ W.T)).sum(2) / sso[:, None] >= I[ok][:, None] - 1e-12).sum(1); done += b
    p[ok] = (c + 1) / (N + 1); return I, p
# (1) per-year
tdidx = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
yr_rows = []; per_year = {}
named = ['Dapagliflozin', 'Empagliflozin', 'Semaglutide', 'Tirzepatide', 'Liraglutide', 'Chlortalidone', 'Indapamide',
         'Bendroflumethiazide', 'Hydrocortisone butyrate', 'Lisinopril', 'Ketoprofen', 'Loperamide hydrochloride', 'Mebeverine hydrochloride']
for k, y in enumerate(years):
    Xy = shares[k].T                                   # (drugs, areas)
    present = Xy.sum(1) > 0
    I, p = moran_all(Xy, 999, seed=100 + k)
    per_year[int(y)] = dict(n_items=int(present.sum()), median_I=float(np.nanmedian(I[present])),
                            n_nominal=int(np.nansum((p[present] < 0.05) & (I[present] > 1e-12))), pct_nominal=float(100 * np.nanmean((p[present] < 0.05) & (I[present] > 1e-12))))
    for d in named:
        j = tdidx.get(d.lower())
        if j is None or not present[j]: yr_rows.append(dict(Drug=d, year=int(y), I=np.nan, p=np.nan)); continue
        Ij, pj = moran_all(Xy[[j]], 9999, seed=1000 + k)
        yr_rows.append(dict(Drug=d, year=int(y), I=float(Ij[0]), p=float(pj[0])))
YR = pd.DataFrame(yr_rows); YR.to_csv('moran_by_year.csv', index=False)
out['per_year'] = per_year
print('formulary-wide Moran by year:')
for y, v in per_year.items(): print(f"  {y}: items {v['n_items']}, median I {v['median_I']:.3f}, nominal P<0.05 {v['n_nominal']} ({v['pct_nominal']:.1f}%)")
print('\nnamed agents, Moran I by year (P):')
piv = YR.pivot(index='Drug', columns='year', values='I').round(3); pp = YR.pivot(index='Drug', columns='year', values='p').round(4)
for d in named:
    print(f'  {d:26s} ' + '  '.join(f"{int(y)}: {piv.loc[d, y]:+.3f} ({pp.loc[d, y]:.4f})" if not np.isnan(piv.loc[d, y]) else f'{int(y)}: n/a' for y in years))
# (2) size of variation
X = A.T
p90 = np.quantile(X, .9, axis=1); p10 = np.quantile(X, .1, axis=1)
ratio = np.where(p10 > 0, p90 / np.where(p10 > 0, p10, 1), np.nan)
vol = R.Drug.map(dict(zip(R.Drug, range(len(R)))))
okr = ~np.isnan(ratio)
out['p90p10'] = dict(n_items_with_p10_gt0=int(okr.sum()), median=float(np.nanmedian(ratio)), iqr=[float(np.nanquantile(ratio, .25)), float(np.nanquantile(ratio, .75))],
                     median_clustered=float(np.nanmedian(ratio[R.clustered.values & okr])), median_not_clustered=float(np.nanmedian(ratio[~R.clustered.values & okr])))
print('\nP90/P10 ratio of within-area share: items with P10>0: %d; median %.2f (IQR %.2f to %.2f); clustered items %.2f vs not clustered %.2f' % (
    okr.sum(), out['p90p10']['median'], *out['p90p10']['iqr'], out['p90p10']['median_clustered'], out['p90p10']['median_not_clustered']))
for d in named:
    i = np.where(R.Drug.str.lower() == d.lower())[0]
    if len(i): i = i[0]; print(f'  {d:26s} P90/P10 {ratio[i]:.2f}  (max/min {X[i].max() / X[i].min() if X[i].min() > 0 else np.inf:.1f})')
R['p90_p10'] = ratio; R.to_csv('moran_af.csv', index=False)
# (3) coding-stable sensitivity
yt = shares.sum(1)                                     # (years, drugs) summed shares > 0 means present
stable = (yt > 0).all(0)
stable_names = set(np.array(drugs_td)[stable])
sel = np.array([d in stable_names for d in R.Drug])
As = A[:, sel]; As = As / As.sum(1, keepdims=True)
Ds = jsd_matrix(As); rs, ps = mantel_all(As, Ds, n_perm=1000, seed=0)
Is, pis = moran_all(As.T, 9999, seed=7)
out['coding_stable'] = dict(n_items=int(sel.sum()), mantel_nominal=int((ps < .05).sum()), pct_mantel=float(100 * (ps < .05).mean()),
                            moran_nominal=int(((pis < .05) & (Is > 1e-12)).sum()), pct_moran=float(100 * ((pis < .05) & (Is > 1e-12)).mean()),
                            corr_mantel_r_with_main=float(np.corrcoef(rs, R.Mantel_r.values[sel])[0, 1]),
                            corr_moran_I_with_main=float(np.corrcoef(Is, R.moran_I.values[sel])[0, 1]))
print('\ncoding-stable items (dispensed in every year): %d; Mantel nominal %d (%.1f%%); Moran nominal %d (%.1f%%); corr with main: Mantel r %.3f, Moran I %.3f' % (
    sel.sum(), out['coding_stable']['mantel_nominal'], out['coding_stable']['pct_mantel'], out['coding_stable']['moran_nominal'],
    out['coding_stable']['pct_moran'], out['coding_stable']['corr_mantel_r_with_main'], out['coding_stable']['corr_moran_I_with_main']))
json.dump(out, open('extras_summary.json', 'w'), indent=2, default=float)
print('\nsaved moran_by_year.csv, extras_summary.json, moran_af.csv (with p90_p10)')
