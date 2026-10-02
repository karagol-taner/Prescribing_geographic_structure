import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Numeric columns of Supplementary Tables S3 and S4 (the phylogeographic comparison), from the outputs of
analyze_af.py, temporal_af.py, fdr_refine.py, mantel_84.py, spatial_analysis.py and spatial_sensitivity.py.
Table S3: the 16 non-appliance items with the strongest nominally significant Mantel signals (1,000-permutation
screen, P<0.05; appliance flag of fdr_refine.py), with national items, Mantel r and P, unadjusted Moran's I and
its nominal P (spatial_analysis.py), and two concentration measures computed on the volume-adjusted shares
(each area's within-area share, so that large areas do not dominate): the share of the total in the single
highest area and the number of areas needed to reach 80% of the total.
Table S4: for nine named agents, Mantel r on the 84 area-period profiles (mantel_84.py) and on the 42 areas,
pattern stability (correlation of the 2021 and 2025 area-share vectors) and national change (items 2025 / 2021)
from temporal_af.py, and unadjusted Moran's I by calendar year from spatial_sensitivity.py.
Writes data/tableS3_items.csv and data/tableS4_agents.csv."""
import numpy as np, pandas as pd

R = pd.read_csv('mantel_af.csv'); F = pd.read_csv('mantel_af_fdr.csv'); M = pd.read_csv('moran_af.csv'); A = np.load('A_af.npy')
assert list(R.Drug) == list(F.Drug) == list(M.Drug)


def areas_for_80(v):
    v = np.sort(v)[::-1]; c = np.cumsum(v) / v.sum()
    return int(np.searchsorted(c, 0.80) + 1)


sel = F[(F.P_Value < 0.05) & ~F.appliance.astype(bool)].sort_values('Mantel_r', ascending=False).head(16)
rows = []
for rank, (i, r) in enumerate(sel.iterrows(), 1):
    v = A[:, i]
    rows.append(dict(rank=rank, item=r.Drug, items_2021_2025=int(R['items'].iloc[i]), mantel_r=r.Mantel_r, mantel_p=r.P_Value,
                     moran_I=M.moran_I.iloc[i], moran_p=M.p_moran.iloc[i], moran_q=M.q_moran.iloc[i],
                     top_area_share_pct=100 * v.max() / v.sum(), areas_for_80pct=areas_for_80(v)))
S3 = pd.DataFrame(rows)
S3.to_csv('tableS3_items.csv', index=False)
print(S3.round(3).to_string(index=False))

m84 = pd.read_csv('mantel_84.csv'); m84['k'] = m84.Drug_Name.str.strip().str.lower()
T = pd.read_csv('temporal_af.csv'); T['k'] = T.Drug.str.strip().str.lower()
Y = pd.read_csv('moran_by_year.csv')
rows = []
for d in ['Dapagliflozin', 'Tirzepatide', 'Liraglutide', 'Semaglutide', 'Empagliflozin', 'Chlortalidone',
          'Hydrocortisone butyrate', 'Ketoprofen', 'Lisinopril']:
    k = d.lower(); f = F[F.Drug.str.strip().str.lower() == k].iloc[0]; t = T[T.k == k].iloc[0]
    row = dict(agent=d, mantel_r_84_profiles=float(m84[m84.k == k].Mantel_r.iloc[0]), mantel_r_42_areas=f.Mantel_r,
               mantel_p_42_areas=f.P_Value, stability_2021_2025=t.stability, national_change_2025_2021=t.growth)
    for y in range(2021, 2026):
        yy = Y[(Y.Drug == d) & (Y.year == y)]
        row[f'moran_I_{y}'] = float(yy.I.iloc[0]) if len(yy) else np.nan
    rows.append(row)
S4 = pd.DataFrame(rows)
S4.to_csv('tableS4_agents.csv', index=False)
print(S4.round(3).to_string(index=False))
