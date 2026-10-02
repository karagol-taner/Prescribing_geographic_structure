import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Step 2: year-resolved shares, 2021-2025, to separate change over time from differences between places.
For each drug, the within-area share vector across 42 areas is computed per year.
- pattern stability  = Pearson r between the 2021 and 2025 area-share vectors
- dispersion (CV)    = between-area CV of the within-area share, per year
- national growth    = national items 2025 / 2021
Genuine spatial signals: high stability, ~flat CV. Temporal adoption signals:
low stability, CV that collapses as the drug diffuses nationally.
"""
import numpy as np, pandas as pd
BY = 'icb_substance_items_by_year.csv'
df = pd.read_csv(BY)
df['DRUG'] = df['DRUG'].astype(str).str.strip()
df['AREA_CODE'] = df['AREA_CODE'].astype(str).str.strip()
years = sorted(df['YEAR'].unique())
areas = sorted(df['AREA_CODE'].unique())
drugs = sorted(df['DRUG'].unique())
print('years', years, '| areas', len(areas), '| drugs', len(drugs))
# national items per year (completeness check)
nat = df.groupby('YEAR')['ITEMS'].sum()
print('national items/yr:', {int(y): f'{v/1e6:.0f}M' for y, v in nat.items()})

ai = {a: i for i, a in enumerate(areas)}
di = {d: j for j, d in enumerate(drugs)}
yi = {y: k for k, y in enumerate(years)}
# shares[year, area, drug] = within-area share
items = np.zeros((len(years), len(areas), len(drugs)))
for r in df.itertuples(index=False):
    items[yi[r.YEAR], ai[r.AREA_CODE], di[r.DRUG]] = r.ITEMS
area_tot = items.sum(2, keepdims=True)                 # (Y, A, 1) area basket per year
shares = np.divide(items, area_tot, out=np.zeros_like(items), where=area_tot > 0)

def cv(x):
    m = x.mean()
    return x.std() / m if m > 0 else np.nan

nat_drug = items.sum(1)                                 # (Y, D) national items per drug/year
rows = []
e, l = yi[years[0]], yi[years[-1]]                      # early=2021, late=2025
for j, d in enumerate(drugs):
    se, sl = shares[e, :, j], shares[l, :, j]
    stab = np.corrcoef(se, sl)[0, 1] if se.std() > 0 and sl.std() > 0 else np.nan
    cvs = [cv(shares[k, :, j]) for k in range(len(years))]
    growth = (nat_drug[l, j] / nat_drug[e, j]) if nat_drug[e, j] > 0 else np.inf
    rows.append(dict(Drug=d, stability=stab, cv_2021=cvs[0], cv_2025=cvs[-1],
                     cv_ratio=(cvs[-1] / cvs[0] if cvs[0] and cvs[0] == cvs[0] else np.nan),
                     items_2021=nat_drug[e, j], items_2025=nat_drug[l, j], growth=growth))
T = pd.DataFrame(rows)
man = pd.read_csv('mantel_af.csv'); man['k'] = man.Drug.str.lower()
T['k'] = T.Drug.str.lower()
T = T.merge(man[['k', 'Mantel_r', 'P_Value']], on='k', how='left').drop(columns='k')
T.to_csv('temporal_af.csv', index=False)
np.save('shares_af.npy', shares); np.save('years_af.npy', np.array(years))
import pickle; pickle.dump((areas, drugs), open('td_index.pkl', 'wb'))

print('\n{:20s} {:>6s} {:>7s} {:>7s} {:>7s} {:>7s}'.format('drug', 'r', 'stab', 'cv21', 'cv25', 'growth'))
for name in ['Chlortalidone', 'Ketoprofen', 'Lisinopril', 'Hydrocortisone butyrate',
             'Dapagliflozin', 'Empagliflozin', 'Tirzepatide', 'Semaglutide', 'Liraglutide']:
    r = T[T.Drug.str.lower() == name.lower()]
    if len(r):
        x = r.iloc[0]
        print('{:20s} {:6.3f} {:7.3f} {:7.3f} {:7.3f} {:7.1f}'.format(
            name[:20], x.Mantel_r, x.stability, x.cv_2021, x.cv_2025, x.growth))
# aggregate contrast: significant drugs split by national growth
sig = T[T.P_Value < 0.05].dropna(subset=['stability', 'growth'])
fast = sig[sig.growth > 2]; slow = sig[sig.growth < 1.5]
print(f'\nsignificant & fast-growing (>2x): n={len(fast)}  median stability={fast.stability.median():.3f}  median cv_ratio={fast.cv_ratio.median():.3f}')
print(f'significant & stable-volume (<1.5x): n={len(slow)}  median stability={slow.stability.median():.3f}  median cv_ratio={slow.cv_ratio.median():.3f}')
print('saved temporal_af.csv, shares_af.npy')
