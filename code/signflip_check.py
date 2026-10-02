import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Robustness of the need-adjusted item-level clustering to unequal error variances.
The Freedman-Lane permutation test assumes exchangeable errors. A sign-flip test, in which the signs of
the residuals are randomised (Rademacher weights) and the result re-projected through the residual-maker
matrix, remains valid when the error variance differs between areas. The sign-flip test is applied to every
item under the primary need model with 9,999 random draws and compared with the Freedman-Lane P-values of
need_adjusted_moran.csv (9,999 permutations, refined to 99,999 for P<0.05); the script reports how many of the
items classed as clustered there also reach P<0.05 and P<0.01 with the sign-flip test."""
import numpy as np, pandas as pd, pickle, json
from need_adjustment_core import primary_covariates, resid_maker, moran_rows
TOL = 1e-12; POS = 1e-12
A = np.load('A_af.npy'); W = np.load('W_contig.npy'); codes = pickle.load(open('codes_af.pkl', 'rb')); n = len(codes)
R = pd.read_csv('need_adjusted_moran.csv')
Mres, k = resid_maker(primary_covariates(codes))
Y = A.T.astype(float); E = Y @ Mres; I = moran_rows(E, W)
g = np.random.default_rng(77); N = 9999; cnt = np.zeros(len(E)); done = 0
while done < N:
    b = min(200, N - done); S = g.choice([-1.0, 1.0], size=(b, n))
    for s in S:
        cnt += moran_rows((E * s) @ Mres, W) >= I - TOL
    done += b
p = (cnt + 1) / (N + 1)
clu = R.clustered_primary.values.astype(bool)
assert np.allclose(I, R.I_primary.values)
out = dict(n_clustered=int(clu.sum()), clustered_with_signflip_p_lt_005=int(((p < 0.05) & clu).sum()),
           clustered_with_signflip_p_lt_001=int(((p < 0.01) & clu).sum()),
           nominal_signflip_all_items=int(((p < 0.05) & (I > POS)).sum()), nominal_freedman_lane_all_items=int(R.nominal_primary.sum()),
           draws=N)
pd.DataFrame(dict(Drug=R.Drug, I_primary=I, p_freedman_lane=R.p_primary, p_signflip=p)).to_csv('signflip_check.csv', index=False)
json.dump(out, open('signflip_check.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
