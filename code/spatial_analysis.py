import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Step 7: unadjusted spatial analyses.
(1) Does overall prescribing similarity (Jensen-Shannon) follow geography?
    Mantel test of JSD against centroid distance; neighbour vs non-neighbour JSD.
(2) Global Moran's I for every item (within-area share), contiguity weights,
    one-sided permutation test; nominal P from 9,999 permutations, refined to
    99,999 for nominally significant items, Benjamini-Hochberg FDR.
    k-nearest-neighbour (k=4) weights as a sensitivity analysis.
(3) Concordance with the Mantel (phylogeographic) signal."""
import numpy as np, pandas as pd, json, time
from scipy.stats import spearmanr
A = np.load('A_af.npy'); D = np.load('D_af.npy'); R = pd.read_csv('mantel_af_fdr.csv')
Wb = np.load('Wb_contig.npy'); W = np.load('W_contig.npy'); Wk = np.load('W_knn4.npy'); G = np.load('Dgeo_km.npy')
n = A.shape[0]; iu = np.triu_indices(n, 1); rng = np.random.default_rng(42)
TOL = 1e-12   # permuted statistics within 1e-12 of the observed value count as ties (at least as large)
out = {}

# ---------- (1) overall prescribing similarity vs geography ----------
dv, gv = D[iu], G[iu]
r_obs = np.corrcoef(dv, gv)[0, 1]; rho_obs = spearmanr(dv, gv).correlation
NP = 9999; cnt = 0; cnt_s = 0
for _ in range(NP):
    p = rng.permutation(n); gp = G[np.ix_(p, p)][iu]
    cnt += np.corrcoef(dv, gp)[0, 1] >= r_obs - TOL
    cnt_s += spearmanr(dv, gp).correlation >= rho_obs - TOL
out['mantel_jsd_vs_distance'] = dict(r=r_obs, p=(cnt + 1) / (NP + 1), spearman=rho_obs, p_spearman=(cnt_s + 1) / (NP + 1))
nb = Wb[iu] == 1
stat = dv[~nb].mean() - dv[nb].mean(); c2 = 0
for _ in range(NP):
    p = rng.permutation(n); nbp = Wb[np.ix_(p, p)][iu] == 1
    c2 += (dv[~nbp].mean() - dv[nbp].mean()) >= stat - TOL
out['neighbours'] = dict(mean_jsd_neighbours=dv[nb].mean(), mean_jsd_non_neighbours=dv[~nb].mean(),
                         pct_more_similar=100 * stat / dv[~nb].mean(), p=(c2 + 1) / (NP + 1), n_neighbour_pairs=int(nb.sum()))
q = np.quantile(gv, [0, .2, .4, .6, .8, 1]); cls = np.clip(np.searchsorted(q, gv, side='right') - 1, 0, 4)
out['distance_classes'] = [dict(km_from=float(q[k]), km_to=float(q[k + 1]), mean_jsd=float(dv[cls == k].mean()),
                                se=float(dv[cls == k].std(ddof=1) / np.sqrt((cls == k).sum())), n_pairs=int((cls == k).sum())) for k in range(5)]
print('JSD vs distance: Mantel r = %.3f (P = %.4f); Spearman %.3f (P = %.4f)' % (r_obs, out['mantel_jsd_vs_distance']['p'], rho_obs, out['mantel_jsd_vs_distance']['p_spearman']))
print('neighbours JSD %.4f vs non-neighbours %.4f (%.1f%% more similar), P = %.4f, %d neighbour pairs' % (
    dv[nb].mean(), dv[~nb].mean(), out['neighbours']['pct_more_similar'], out['neighbours']['p'], nb.sum()))
for c in out['distance_classes']: print('  %4.0f-%4.0f km: mean JSD %.4f (SE %.4f, %d pairs)' % (c['km_from'], c['km_to'], c['mean_jsd'], c['se'], c['n_pairs']))

# ---------- (2) Moran's I for every item ----------
X = A.T.astype(float); Z = X - X.mean(1, keepdims=True); ss = (Z ** 2).sum(1)
def moran(Zm, Wm): return (Zm * (Zm @ Wm.T)).sum(-1) / ss_sel
def perm_p(Zsel, Wm, I_obs, N, B=100, seed=1):
    g_ = np.random.default_rng(seed); c = np.zeros(len(I_obs)); done = 0
    while done < N:
        b = min(B, N - done); P = np.array([g_.permutation(n) for _ in range(b)])
        Zp = Zsel[:, P]; Ip = (Zp * (Zp @ Wm.T)).sum(2) / ss_sel[:, None]
        c += (Ip >= I_obs[:, None] - TOL).sum(1); done += b
    return (c + 1) / (N + 1)
ss_sel = ss; I = moran(Z, W); Ik = moran(Z, Wk)
t = time.time(); p9 = perm_p(Z, W, I, 9999); pk = perm_p(Z, Wk, Ik, 9999, seed=2)
print(f'\nMoran permutations (all items, contiguity + kNN) in {time.time() - t:.0f}s')
cand = np.where(p9 < 0.05)[0]; ss_sel = ss[cand]
t = time.time(); pref = perm_p(Z[cand], W, I[cand], 99999, B=200, seed=3)
print(f'refined {len(cand)} candidates with 99,999 permutations in {time.time() - t:.0f}s')
ss_sel = ss
p_all = p9.copy(); p_all[cand] = pref; m = len(p_all); o = np.argsort(p_all)
qv = np.empty(m); qv[o] = np.minimum.accumulate((p_all[o] * m / np.arange(1, m + 1))[::-1])[::-1]; qv = np.minimum(qv, 1)
M = pd.DataFrame({'Drug': R.Drug, 'n_areas_used': (A > 0).sum(0), 'moran_I': I, 'p_moran': p_all, 'q_moran': qv, 'moran_I_knn4': Ik, 'p_moran_knn4': pk,
                  'Mantel_r': R.Mantel_r, 'p_mantel_1000': R.P_Value, 'p_mantel_refined': R.p_refined, 'q_mantel': R.q_BH, 'appliance': R.appliance})
# an item is classed as spatially clustered when Moran's I is positive and significant: for items used in
# only one or a few areas the permutation distribution is discrete, and a negative I can reach P<0.05
POS = 1e-12   # 'positive' means above floating-point zero (a few single-area items have I of order 1e-18 with kNN weights)
M['nominal'] = (M.p_moran < 0.05) & (M.moran_I > POS); M['clustered'] = (M.q_moran < 0.05) & (M.moran_I > POS)
M.to_csv('moran_af.csv', index=False)
EI = -1 / (n - 1)
out['moran'] = dict(expected_I=EI, median_I=float(np.median(I)), n_nominal=int(M.nominal.sum()), n_fdr05=int(M.clustered.sum()),
                    n_fdr10=int(((qv < 0.10) & (I > POS)).sum()), n_nominal_knn=int(((pk < 0.05) & (Ik > POS)).sum()),
                    n_negative_I_with_p_below_05=int(((p_all < 0.05) & (I <= POS)).sum()), n_negative_I_with_q_below_05=int(((qv < 0.05) & (I <= POS)).sum()),
                    spearman_contig_knn=float(spearmanr(I, Ik).correlation),
                    spearman_moran_mantel=float(spearmanr(I, R.Mantel_r).correlation))
print("Moran's I: E[I] = %.4f, median I = %.4f; nominal P<0.05: %d (expected by chance %.0f); FDR q<0.05: %d; q<0.10: %d" % (
    EI, np.median(I), out['moran']['n_nominal'], 0.05 * m, out['moran']['n_fdr05'], out['moran']['n_fdr10']))
print("kNN-4 sensitivity: nominal %d; Spearman(I contiguity, I kNN) = %.3f" % (out['moran']['n_nominal_knn'], out['moran']['spearman_contig_knn']))
print("Concordance with Mantel: Spearman(Moran I, Mantel r) = %.3f" % out['moran']['spearman_moran_mantel'])
nm_ = (M.p_mantel_1000 < .05); mo_ = M.nominal
print(f'cross-tab nominal: both {int((nm_ & mo_).sum())}, Mantel only {int((nm_ & ~mo_).sum())}, Moran only {int((~nm_ & mo_).sum())}')
print('\nFDR-significant Moran items (q<0.10):')
print(M[M.q_moran < 0.10].sort_values('q_moran')[['Drug', 'moran_I', 'p_moran', 'q_moran', 'Mantel_r', 'appliance']].to_string(index=False))
print('\nnamed agents:')
for d in ['Chlortalidone', 'Ketoprofen', 'Hydrocortisone butyrate', 'Lisinopril', 'Padimate O', 'Cladribine (Immunomodulating)',
          'Thyrotropin alfa', 'Sucrase', 'Catheters', 'Night Drainage Bags', 'Incontinence Sheaths', 'Skin Fillers And Protectives',
          'Dapagliflozin', 'Empagliflozin', 'Semaglutide', 'Tirzepatide', 'Liraglutide', 'Indapamide', 'Bendroflumethiazide']:
    r = M[M.Drug.str.lower() == d.lower()]
    if len(r): x = r.iloc[0]; print(f'  {d:30s} I={x.moran_I:+.3f} p={x.p_moran:.5f} q={x.q_moran:.3f} | I_knn={x.moran_I_knn4:+.3f} p_knn={x.p_moran_knn4:.4f} | Mantel r={x.Mantel_r:.3f}')
json.dump(out, open('geo_summary.json', 'w'), indent=2, default=float)
print('\nsaved moran_af.csv, geo_summary.json')
