import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Numbers of the unadjusted spatial and phylogeographic analyses quoted in the paper (Table S2, its note and
the main text), printed to standard output: python code/summary_numbers.py > data/summary_numbers.txt"""
import pandas as pd, numpy as np, json
from scipy.stats import spearmanr
M = pd.read_csv('moran_af.csv'); F = pd.read_csv('mantel_af_fdr.csv'); A = np.load('A_af.npy')
full = pd.read_csv('icb_area_substance_items_full.csv', usecols=['DRUG', 'SECTION', 'ITEMS'])
full['DRUG'] = full.DRUG.astype(str).str.strip()
vol = full.groupby('DRUG').ITEMS.sum()
sec = full.groupby(['DRUG', 'SECTION']).ITEMS.sum().reset_index().sort_values('ITEMS').groupby('DRUG').last().SECTION   # section with most items, as in the Methods
M['vol'] = M.Drug.map(vol); M['sec'] = M.Drug.map(sec)
print('volume missing:', M.vol.isna().sum(), ' total items', M.vol.sum())
cl = M.clustered.astype(bool); nomM = M.nominal.astype(bool)
print("Moran I quantiles: median %.3f IQR %.3f to %.3f; mean %.3f" % (M.moran_I.median(), M.moran_I.quantile(.25), M.moran_I.quantile(.75), M.moran_I.mean()))
print('FDR-clustered items: %d; share of all items dispensed %.1f%%' % (cl.sum(), 100 * M.vol[cl].sum() / M.vol.sum()))
print('nominal-clustered share of volume %.1f%%' % (100 * M.vol[nomM].sum() / M.vol.sum()))
mf = F.q_BH < 0.05; print('Mantel FDR items share of volume %.3f%%' % (100 * M.vol[mf.values].sum() / M.vol.sum()))
# cross-tabulation with the 1,000-permutation Mantel P
mn = (F.P_Value < 0.05).values
print('non-positive-I items with P<0.05 (excluded from nominal):', int(((M.p_moran < .05) & (M.moran_I <= 1e-12)).sum()), '; with q<0.05:', int(((M.q_moran < .05) & (M.moran_I <= 1e-12)).sum()))
print('kNN-4 nominal (positive I):', int(((M.p_moran_knn4 < .05) & (M.moran_I_knn4 > 1e-12)).sum()))
print('Mantel nominal (1000 perm):', mn.sum(), '| refined <0.05:', (F.p_refined < .05).sum())
print('cross-tab (Mantel 1000-perm nominal vs Moran nominal): both %d, Mantel only %d, Moran only %d, neither %d' % ((mn & nomM).sum(), (mn & ~nomM).sum(), (~mn & nomM).sum(), (~mn & ~nomM).sum()))
print('of Mantel nominal, Moran FDR clustered: %d' % (mn & cl).sum())
print('Spearman Moran I vs Mantel r: %.3f' % spearmanr(M.moran_I, M.Mantel_r).correlation)
ap = M.appliance.astype(bool)
print('appliances flagged %d; Moran nominal %d; Moran FDR %d; median I %.3f; Mantel nominal appliances %d (mean r %.2f)' % (ap.sum(), (ap & nomM).sum(), (ap & cl).sum(), M.moran_I[ap].median(), (ap & mn).sum(), M.Mantel_r[ap & mn].mean()))
print('strong Mantel (r>0.6 & nominal):', ((M.Mantel_r > 0.6) & mn).sum(), 'of which appliances', ((M.Mantel_r > 0.6) & mn & ap).sum(), '; their median Moran I %.3f; any clustered: %d' % (M.moran_I[(M.Mantel_r > 0.6) & mn].median(), ((M.Mantel_r > 0.6) & mn & cl).sum()))
print('\nMantel FDR survivors:'); print(M[mf.values][['Drug', 'Mantel_r', 'moran_I', 'p_moran']].to_string(index=False))
# P90/P10
r = M.p90_p10
print('\nP90/P10: n %d, median %.2f IQR %.2f-%.2f; clustered %.2f; not %.2f' % (r.notna().sum(), r.median(), r.quantile(.25), r.quantile(.75), r[cl].median(), r[~cl].median()))
t = pd.qcut(M.vol.rank(method='first'), 3, labels=['low', 'mid', 'high'])
for k in ['low', 'mid', 'high']:
    s = t == k; print('  tertile %s: n %d, items with P10>0 %d, median P90/P10 %.2f, clustered %.1f%%, clustered median ratio %.2f' % (k, s.sum(), r[s].notna().sum(), r[s].median(), 100 * cl[s].mean(), r[s & cl].median()))
# volume-weighted median ratio
ok = r.notna(); o = np.argsort(r[ok].values); w = M.vol[ok].values[o]; cw = np.cumsum(w) / w.sum()
print('  volume-weighted median P90/P10 %.2f' % r[ok].values[o][np.searchsorted(cw, 0.5)])
top100 = M.nlargest(100, 'vol'); print('top-100 most dispensed: clustered %d, nominal %d, median ratio %.2f, median I %.3f' % (top100.clustered.sum(), top100.nominal.sum(), top100.p90_p10.median(), top100.moran_I.median()))
top20 = M.nlargest(20, 'vol'); print('top-20:'); print(top20[['Drug', 'moran_I', 'q_moran', 'p90_p10']].round(3).to_string(index=False))
print('\nTop 12 by Moran I:'); print(M.nlargest(12, 'moran_I')[['Drug', 'sec', 'moran_I', 'q_moran', 'moran_I_knn4', 'p90_p10', 'vol']].round(3).to_string(index=False))
for d in ['Loperamide hydrochloride', 'Mebeverine hydrochloride', 'Dapagliflozin', 'Empagliflozin', 'Chlortalidone', 'Ketoprofen', 'Catheters', 'Lisinopril', 'Hydrocortisone butyrate', 'Semaglutide', 'Indapamide', 'Bendroflumethiazide', 'Magnesium lactate']:
    x = M[M.Drug == d].iloc[0]; i = M.index[M.Drug == d][0]; a = A[:, i]
    print('  %-26s I %.3f p %.5f q %.4f | knn %.3f p %.4f | P90/P10 %.2f max/min %.1f' % (d, x.moran_I, x.p_moran, x.q_moran, x.moran_I_knn4, x.p_moran_knn4, x.p90_p10, a.max() / a.min() if a.min() > 0 else np.inf))
