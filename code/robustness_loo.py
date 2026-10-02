import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Leave-one-out robustness check.
Each item contributes to the overall Jensen-Shannon matrix against which it is
tested. Here each tested item's Mantel r is recomputed against a distance matrix
built without that item, for ten named agents and the 20 items with the highest mean within-area share."""
import numpy as np, pandas as pd
from method import jsd_matrix

A = np.load('A_af.npy'); R = pd.read_csv('mantel_af.csv')
n = A.shape[0]; iu = np.triu_indices(n, 1); share = A.mean(0)

def mantel_r(x, Dm):
    return np.corrcoef(np.abs(x[iu[0]] - x[iu[1]]), Dm[iu])[0, 1]

idx = {d.strip().lower(): i for i, d in enumerate(R.Drug)}
named = ['chlortalidone', 'ketoprofen', 'hydrocortisone butyrate', 'lisinopril', 'thyrotropin alfa',
         'padimate o', 'catheters', 'dapagliflozin', 'empagliflozin', 'semaglutide']
targets = sorted(set([idx[x] for x in named if x in idx] + list(np.argsort(-share)[:20])))
rows = []
for j in targets:
    r_loo = mantel_r(A[:, j], jsd_matrix(np.delete(A, j, axis=1)))
    rows.append(dict(Drug=R.Drug.iloc[j], mean_within_area_share=share[j],
                     r=R.Mantel_r.iloc[j], r_leave_one_out=r_loo, change=r_loo - R.Mantel_r.iloc[j]))
out = pd.DataFrame(rows).sort_values('change')
out.to_csv('loo_robustness.csv', index=False)
print(out.to_string(index=False))
print('max |change| = %.4f' % out.change.abs().max())
