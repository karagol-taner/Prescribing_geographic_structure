import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""National dispensing trends, 2021-2025 (Table S5 and the numbers behind Figure S7).
Quarterly national totals of items and net ingredient cost, mean unique chemicals per
area, and the between-area coefficient of variation (CV, sample SD / mean) of items.
Monotonic trend: Mann-Kendall test (normal approximation with continuity correction and
tie-corrected variance); slope: Theil-Sen; seasonality: seasonal Kendall test by quarter."""
import json
import numpy as np, pandas as pd
from scipy.stats import norm, theilslopes

D = pd.read_csv('temporal_data.csv')          # Period, STP_CODE, Total_Items, Total_Cost, Unique_Chemicals
D['Year'] = D.Period.str[:4].astype(int)
q = D.groupby('Period').agg(items=('Total_Items', 'sum'), cost=('Total_Cost', 'sum'), chem=('Unique_Chemicals', 'mean'),
                            cv=('Total_Items', lambda v: v.std(ddof=1) / v.mean()))
q['Year'] = q.index.str[:4].astype(int)
season = np.array([int(p[-1]) for p in q.index])


def mann_kendall(x):
    x = np.asarray(x, float); n = len(x)
    s = sum(np.sign(x[j] - x[i]) for i in range(n) for j in range(i + 1, n))
    _, c = np.unique(x, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - sum(t * (t - 1) * (2 * t + 5) for t in c)) / 18
    z = (s - np.sign(s)) / np.sqrt(var)
    return float(s), float(2 * (1 - norm.cdf(abs(z))))


def seasonal_kendall(x, seasons):
    x = np.asarray(x, float); S = 0.0; V = 0.0
    for s_ in np.unique(seasons):
        xs = x[seasons == s_]; n = len(xs)
        S += sum(np.sign(xs[j] - xs[i]) for i in range(n) for j in range(i + 1, n))
        V += n * (n - 1) * (2 * n + 5) / 18
    z = (S - np.sign(S)) / np.sqrt(V)
    return float(S), float(2 * (1 - norm.cdf(abs(z))))


yearly = pd.DataFrame({'items_bn': D.groupby('Year').Total_Items.sum() / 1e9,
                       'cost_gbp_bn': D.groupby('Year').Total_Cost.sum() / 1e9,
                       'mean_unique_chemicals': D.groupby('Year').Unique_Chemicals.mean(),
                       'between_area_cv_mean_of_quarters': q.groupby('Year').cv.mean()})
print('Table S5 (yearly):'); print(yearly.round(3).to_string())
print('totals: items %.2f bn, cost GBP %.1f bn' % (D.Total_Items.sum() / 1e9, D.Total_Cost.sum() / 1e9))
out = {'yearly': yearly.round(4).to_dict(orient='index')}
for col, lab in [('items', 'items'), ('cost', 'net ingredient cost'), ('cv', 'between-area CV'), ('chem', 'mean unique chemicals')]:
    s, p = mann_kendall(q[col].values); slope = theilslopes(q[col].values, np.arange(len(q)))[0]
    S2, p2 = seasonal_kendall(q[col].values, season)
    out[col] = dict(mann_kendall_S=s, mann_kendall_p=p, theil_sen_slope_per_quarter=float(slope), seasonal_kendall_S=S2, seasonal_kendall_p=p2)
    print(f'{lab:24s} Mann-Kendall P={p:.1e}  Theil-Sen slope per quarter={slope:.4g}  seasonal Kendall P={p2:.1e}')
json.dump(out, open('trends_summary.json', 'w'), indent=2)
print('saved trends_summary.json')
