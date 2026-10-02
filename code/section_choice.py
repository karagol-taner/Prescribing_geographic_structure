import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Formulary-wide access versus choice. Each item is assigned to its BNF section (the section under
which most of its items were dispensed). Access = a section's total within-area share; choice = an
item's share of its section's dispensing in each area (sections with at least two items). Both are
tested with need-adjusted Moran's I (choice only for sections used in all 42 areas) (primary need model of need_adjustment.py: % aged 65+, % aged under
18, IMD 2025, first two principal components of 21 QOF prevalences; Freedman-Lane permutation, 9,999,
refined to 99,999 for P<0.05; BH FDR within each family; clustered = positive I and q<0.05).
Unadjusted Moran's I is reported alongside."""
import numpy as np, pandas as pd, pickle, json
from need_adjustment_core import resid_maker, moran_rows, fl_perm, bh, primary_covariates
TOL, POS = 1e-12, 1e-12
A = np.load('A_af.npy'); W = np.load('W_contig.npy'); codes = pickle.load(open('codes_af.pkl', 'rb'))
M0 = pd.read_csv('moran_af.csv'); n = len(codes)
full = pd.read_csv('icb_area_substance_items_full.csv', usecols=['DRUG', 'SECTION', 'ITEMS'])
full['DRUG'] = full.DRUG.astype(str).str.strip()
sec = full.groupby(['DRUG', 'SECTION']).ITEMS.sum().reset_index().sort_values('ITEMS').groupby('DRUG').last().SECTION
vol = full.groupby('DRUG').ITEMS.sum()
drugs = M0.Drug.str.strip(); s_of = drugs.map(sec).values; v_of = drugs.map(vol).values
assert not pd.isna(s_of).any()
Mres, k = resid_maker(primary_covariates(codes))
secs = pd.Series(s_of).value_counts(); secs = secs[secs >= 2].index.tolist()
T = np.array([A[:, s_of == s].sum(1) for s in secs])                       # section totals (sections x areas)
rows, ch = [], []
for s, t in zip(secs, T):
    if not (t > 0).all():                                                      # choice undefined where a section is unused
        continue
    for j in np.where(s_of == s)[0]:
        rows.append((s, M0.Drug[j], v_of[j])); ch.append(A[:, j] / t)
Ch = np.array(ch)
out = {'n_sections': len(secs), 'n_sections_used_in_all_areas': int((T > 0).all(1).sum()), 'n_choice_items': len(Ch)}
res = {}
for fam, Y in [('access', T), ('choice', Ch)]:
    E0 = Y - Y.mean(1, keepdims=True); I0 = moran_rows(E0, W)
    p0 = fl_perm(E0, np.eye(n) - np.ones((n, n)) / n, W, I0, 9999, seed=21)       # unadjusted (centring only)
    E = Y @ Mres; I = moran_rows(E, W); p = fl_perm(E, Mres, W, I, 9999, seed=22)
    cand = np.where(p < 0.05)[0]; p[cand] = fl_perm(E[cand], Mres, W, I[cand], 99999, seed=23, B=200)
    cand0 = np.where(p0 < 0.05)[0]; p0[cand0] = fl_perm(E0[cand0], np.eye(n) - np.ones((n, n)) / n, W, I0[cand0], 99999, seed=24, B=200)
    q0, q = bh(p0), bh(p)
    res[fam] = pd.DataFrame({'I_raw': I0, 'p_raw': p0, 'q_raw': q0, 'I_adj': I, 'p_adj': p, 'q_adj': q})
    out[fam] = dict(n=len(Y), clustered_raw=int(((q0 < 0.05) & (I0 > POS)).sum()), clustered_adj=int(((q < 0.05) & (I > POS)).sum()),
                    nominal_raw=int(((p0 < 0.05) & (I0 > POS)).sum()), nominal_adj=int(((p < 0.05) & (I > POS)).sum()),
                    median_I_raw=float(np.median(I0)), median_I_adj=float(np.median(I)))
acc = res['access']; acc.insert(0, 'section', secs); acc.insert(1, 'n_items', [int((s_of == s).sum()) for s in secs])
acc.insert(2, 'items_dispensed', [float(v_of[s_of == s].sum()) for s in secs]); acc.to_csv('section_access.csv', index=False)
cho = res['choice']; cho.insert(0, 'section', [r[0] for r in rows]); cho.insert(1, 'Drug', [r[1] for r in rows]); cho.insert(2, 'items_dispensed', [r[2] for r in rows])
cho.to_csv('section_choice.csv', index=False)
# volume-restricted comparison: choices for items in the most dispensed half of choice items
big = cho.items_dispensed >= cho.items_dispensed.median()
out['choice_upper_half_volume'] = dict(n=int(big.sum()), clustered_adj=int(((cho.q_adj < 0.05) & (cho.I_adj > POS) & big).sum()),
                                       nominal_adj=int(((cho.p_adj < 0.05) & (cho.I_adj > POS) & big).sum()))
for fam in ['access', 'choice']:
    d = res[fam]; w = acc.items_dispensed.values if fam == 'access' else cho.items_dispensed.values
    out[fam]['volume_share_clustered_adj'] = float(w[((d.q_adj < 0.05) & (d.I_adj > POS)).values].sum() / w.sum())
json.dump(out, open('section_choice_summary.json', 'w'), indent=2)
print(json.dumps(out, indent=1))
