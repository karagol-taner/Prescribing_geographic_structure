import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Do neighbouring areas prescribe more alike once need is accounted for?
Overall dissimilarity is the Jensen-Shannon distance between area profiles (as in method.jsd_matrix:
+1e-10 smoothing, rows renormalised, natural log, square root). The primary analysis excludes the 49
dispensing-appliance items, whose dispensing is attributed to the few areas that host dispensing
appliance contractors rather than to patients' areas (profiles are renormalised over the remaining
2,111 items); the analysis of all 2,160 items is a sensitivity analysis.
Need-adjusted profiles: each item's share in an area is replaced by the mean share across areas plus
the area's OLS residual from the need model, i.e. the profile the area would have if its need
covariates were average; negative values are set to zero and profiles renormalised. With an
intercept-only model this reproduces the unadjusted profiles exactly.
Statistics: percentage by which bordering pairs are more similar than other pairs, and the rank
(Spearman) correlation of dissimilarity with centroid distance. Unadjusted: permutation of area labels.
Need-adjusted: Freedman-Lane-type permutation, in which the rows (areas) of the residual matrix are
permuted and re-projected through the residual-maker matrix before the profiles and distances are
recomputed, so that the null distribution reflects the dependence induced by the regression.
P=(b+1)/(N+1)."""
import numpy as np, pandas as pd, pickle, json, time
from need_adjustment_core import primary_covariates, resid_maker, load_covariates
from method import jsd_matrix
A = np.load('A_af.npy'); Wb = np.load('Wb_contig.npy'); G = np.load('Dgeo_km.npy')
codes = pickle.load(open('codes_af.pkl', 'rb')); n = len(codes); iu = np.triu_indices(n, 1)
gv = G[iu]; nb = Wb[iu] == 1
app = pd.read_csv('moran_af.csv').appliance.astype(bool).values
C = load_covariates(codes)
M1 = np.eye(n) - np.ones((n, n)) / n
Mp, _ = resid_maker(primary_covariates(codes))
Ma, _ = resid_maker(np.column_stack([C.age_65plus, C.age_u18, C.imd2025]))
q = np.quantile(gv, [0, .2, .4, .6, .8, 1]); cls = np.clip(np.searchsorted(q, gv, side='right') - 1, 0, 4)


def jsd_pairs(P, eps=1e-10):
    """Vectorised Jensen-Shannon distance for all pairs (upper triangle) of the rows of P (areas x items)."""
    P = P + eps; P = P / P.sum(1, keepdims=True)
    p, r = P[iu[0]], P[iu[1]]; m = 0.5 * (p + r)
    return np.sqrt(np.maximum(0.5 * (p * np.log(p / m)).sum(1) + 0.5 * (r * np.log(r / m)).sum(1), 0))


def rank(v):
    o = np.argsort(v); r = np.empty(len(v)); r[o] = np.arange(len(v)); return r


rg = rank(gv); rg = (rg - rg.mean()) / np.sqrt(((rg - rg.mean()) ** 2).sum())


def spearman(dv):
    rd = rank(dv); rd = rd - rd.mean(); return float(rd @ rg / np.sqrt((rd ** 2).sum()))


def contrast(dv, mask=nb):
    return float(100 * (dv[~mask].mean() - dv[mask].mean()) / dv[~mask].mean())


def profiles(Y, Mr):
    m = Y.mean(1); P = m[:, None] + Y @ Mr          # items x areas
    return np.clip(P, 0, None).T, int((P < 0).sum())


def analyse(Y, Mr, N, seed, adjusted):
    P, nneg = profiles(Y, Mr); dv = jsd_pairs(P); s = contrast(dv); rho = spearman(dv)
    rng = np.random.default_rng(seed); cs = cr = 0; null_s = []
    E = Y @ Mr; m = Y.mean(1)
    for _ in range(N):
        pi = rng.permutation(n)
        if adjusted:                                    # Freedman-Lane-type: permute residual rows, re-project
            d_ = jsd_pairs(np.clip(m[:, None] + E[:, pi] @ Mr, 0, None).T); s_ = contrast(d_); r_ = spearman(d_)
        else:                                           # permutation of area labels
            Pm = np.zeros((n, n)); Pm[iu] = dv; Pm = Pm + Pm.T; d_ = Pm[np.ix_(pi, pi)][iu]; s_ = contrast(d_); r_ = spearman(d_)
        cs += s_ >= s - 1e-12; cr += r_ >= rho - 1e-12; null_s.append(s_)
    return dict(mean_neighbours=float(dv[nb].mean()), mean_non_neighbours=float(dv[~nb].mean()), pct_more_similar=s,
                p_neighbours=float((cs + 1) / (N + 1)), null_mean_pct=float(np.mean(null_s)), spearman_distance=rho,
                p_spearman=float((cr + 1) / (N + 1)), n_perm=N, test='Freedman-Lane-type' if adjusted else 'area labels',
                relative_distance_bordering=float(dv[nb].mean() / dv.mean()),
                relative_distance_by_quintile=[float(dv[cls == k].mean() / dv.mean()) for k in range(5)],
                distance_quintile_km=[float(x) for x in q], clipped_cells=nneg), dv


t0 = time.time(); out = {}
Yall = A.T.astype(float); Yna = A[:, ~app].T.astype(float)
out['n_items_all'] = int(len(Yall)); out['n_appliances'] = int(app.sum()); out['n_items_primary'] = int(len(Yna))
# check: intercept-only construction reproduces the unadjusted distances of method.jsd_matrix
D = np.load('D_af.npy'); out['max_abs_diff_intercept_only_vs_raw_jsd'] = float(np.abs(jsd_pairs(profiles(Yall, M1)[0]) - D[iu]).max())
out['max_abs_diff_vectorised_vs_method_jsd'] = float(np.abs(jsd_pairs(A[:, ~app]) - jsd_matrix(A[:, ~app])[iu]).max())
# contribution of the appliances to overall dissimilarity
dv_all = jsd_pairs(Yall.T); dv_na = jsd_pairs(Yna.T)
out['appliance_share_of_mean_shares'] = float(Yall[app].sum() / Yall.sum())
res = {}
for label, Y, N0, N1, s0 in [('primary_excluding_appliances', Yna, 9999, 9999, 41), ('all_items', Yall, 9999, 4999, 51)]:
    r = {}
    r['raw'], d0 = analyse(Y, M1, N0, s0, False)
    r['age_deprivation_adjusted'], dd = analyse(Y, Ma, N1, s0 + 1, True)
    r['need_adjusted'], da = analyse(Y, Mp, N1, s0 + 2, True)
    res[label] = r
    if label == 'primary_excluding_appliances':
        for nm_, dvv in [('D_unadjusted_excl_appliances.npy', d0), ('D_age_deprivation_adjusted.npy', dd), ('D_need_adjusted.npy', da)]:
            Dm = np.zeros((n, n)); Dm[iu] = dvv; np.save(nm_, Dm + Dm.T)
    print(label, json.dumps({k: (round(v['pct_more_similar'], 2), v['p_neighbours'], round(v['spearman_distance'], 3), v['p_spearman']) for k, v in r.items()}),
          f'{time.time() - t0:.0f}s', flush=True)
out.update(res)
# items dispensed in every area (excluding appliances): clipping is then rare
every = (Yna > 0).all(1)
full = pd.read_csv('icb_area_substance_items_full.csv', usecols=['DRUG', 'ITEMS']); full['DRUG'] = full.DRUG.astype(str).str.strip()
vol = full.groupby('DRUG').ITEMS.sum()
drugs = pd.read_csv('moran_af.csv').Drug.str.strip().values
v_all = vol.reindex(drugs).fillna(0).values
v_na = v_all[~app]
out['items_used_in_every_area'] = dict(n_items=int(every.sum()), share_of_items_dispensed_nationally=float(v_na[every].sum() / v_all.sum()))
Pn = Yna.mean(1)[:, None] + Yna @ Mp
out['items_used_in_every_area']['share_of_clipped_cells_in_items_not_used_everywhere'] = float((Pn[~every] < 0).sum() / (Pn < 0).sum())
out['items_used_in_every_area']['raw'], _ = analyse(Yna[every], M1, 4999, 61, False)
out['items_used_in_every_area']['need_adjusted'], _ = analyse(Yna[every], Mp, 4999, 62, True)
out['appliance_share_of_items_dispensed_nationally'] = float(v_all[app].sum() / v_all.sum())
json.dump(out, open('need_similarity.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if not isinstance(v, dict)}, indent=1))
print('every-area:', json.dumps({k: (round(v['pct_more_similar'], 2), v['p_neighbours']) for k, v in out['items_used_in_every_area'].items() if isinstance(v, dict)}))
