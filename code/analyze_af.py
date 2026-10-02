import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Step 1: 42-area profiles and the per-item phylogeographic (Mantel) signal.
Sums the aggregated monthly totals to the 42 ICB area codes (over the STP and ICB names of each area and
over BNF sections; item names stripped of leading and trailing spaces, giving 2,160 items), computes the
Jensen-Shannon distance between area profiles and each item's Mantel r with 1,000 permutations (method.py).
"""
import sys, numpy as np, pandas as pd, pickle, json
from method import jsd_matrix, mantel_all

FULL = sys.argv[1] if len(sys.argv) > 1 else 'icb_area_substance_items_full.csv'
df = pd.read_csv(FULL)
cols = {c.strip().upper(): c for c in df.columns}
ac, dc, it = cols['AREA_CODE'], cols['DRUG'], cols['ITEMS']
df[dc] = df[dc].astype(str).str.strip()
df[ac] = df[ac].astype(str).str.strip()

# --- collapse to area code x item (sum over the two area names and over sections) ---
ad = df.groupby([ac, dc], as_index=False)[it].sum()
codes = sorted(ad[ac].unique())
mat = ad.pivot(index=ac, columns=dc, values=it).reindex(codes).fillna(0.0)
counts = mat.values.astype(float)              # (n_areas, n_drugs) raw items
drugs = list(mat.columns)
A = counts / counts.sum(1, keepdims=True)      # row proportions
print(f'areas: {len(codes)}  drugs: {len(drugs)}  total items: {counts.sum():,.0f}')

# --- Jensen-Shannon distance and Mantel test (method.py) ---
D = jsd_matrix(A)
r, p = mantel_all(A, D, n_perm=1000, seed=0)
res = pd.DataFrame({'Drug': drugs, 'Mantel_r': r, 'P_Value': p,
                    'items': counts.sum(0), 'top_share': counts.max(0) / counts.sum(0)})
res.to_csv('mantel_af.csv', index=False)
np.save('A_af.npy', A); np.save('D_af.npy', D)
pickle.dump(codes, open('codes_af.pkl', 'wb'))

nsig = int((p < 0.05).sum())
strong = int(((p < 0.05) & (r > 0.6)).sum())
print(f'\n42 areas: significant {nsig}/{len(drugs)} ({100*nsig/len(drugs):.1f}%) | strong r>0.6 {strong} | median r {np.median(r):.4f}')
for name in ['Chlortalidone', 'Dapagliflozin', 'Catheters', 'Thyrotropin Alfa', 'Ketoprofen', 'Lisinopril', 'Tirzepatide', 'Semaglutide']:
    row = res[res.Drug.str.lower() == name.lower()]
    print(f'  {name:18s} r {row.Mantel_r.iloc[0]:.4f}  P {row.P_Value.iloc[0]:.3f}')

summary = dict(areas=len(codes), drugs=len(drugs), total_items=float(counts.sum()),
               n_sig=nsig, pct_sig=100*nsig/len(drugs), strong=strong, median_r=float(np.median(r)))
json.dump(summary, open('af_summary.json', 'w'), indent=2)
print('\nsaved mantel_af.csv, A_af.npy, D_af.npy, af_summary.json')
