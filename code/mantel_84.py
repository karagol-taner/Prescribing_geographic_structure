import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Phylogeographic signal on the naive 84 area-period profiles (Table S2, Table S4 and Figure S6A).
Each ICB appears in the source data under its former Sustainability and Transformation Partnership (STP)
name and under its ICB name; treating the two commissioning-era profiles of each area as separate units
gives 84 profiles. Item names are used exactly as they appear in the source data (six names carry a leading
or trailing space and are kept as separate items, giving 2,161 items). Mantel r is the Pearson correlation
between each item's pairwise absolute differences in within-profile share and the Jensen-Shannon distance
between profiles (method.py); P from 1,000 area-label permutations (seed 0). Writes mantel_84.csv."""
import numpy as np, pandas as pd
from method import jsd_matrix, mantel_all

F = pd.read_csv('icb_area_substance_items_full.csv')
prof = (F.groupby(['AREA_CODE', 'AREA_NAME', 'DRUG'], as_index=False).ITEMS.sum()
        .pivot_table(index=['AREA_CODE', 'AREA_NAME'], columns='DRUG', values='ITEMS', fill_value=0))
assert prof.shape[0] == 84, prof.shape
A = prof.values / prof.values.sum(1, keepdims=True)
D = jsd_matrix(A)
r, p = mantel_all(A, D, n_perm=1000, seed=0)
out = pd.DataFrame({'Drug_Name': prof.columns, 'Mantel_r': r, 'P_Value': p})
out.to_csv('mantel_84.csv', index=False)
nom = int((p < 0.05).sum())
print(f'{len(out):,} items on {prof.shape[0]} profiles; P<0.05: {nom} ({100 * nom / len(out):.1f}%)')
for d in ['Dapagliflozin', 'Tirzepatide', 'Liraglutide', 'Semaglutide', 'Empagliflozin', 'Chlortalidone',
          'Hydrocortisone butyrate', 'Ketoprofen', 'Lisinopril']:
    print(f'  {d:24s} r = {float(out.loc[out.Drug_Name == d, "Mantel_r"].iloc[0]):.2f}')
print('saved mantel_84.csv')
