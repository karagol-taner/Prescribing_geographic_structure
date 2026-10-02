import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""FDR sensitivity: refine permutation p-values to 99,999 permutations for the
nominally significant items, then Benjamini-Hochberg across all 2,160 items.
Permutation p uses (count+1)/(N+1) (Phipson & Smyth). Permuted coefficients are compared
with a 1e-12 tolerance, so that permutations reproducing the observed arrangement (exact ties,
which batched floating-point arithmetic can place a few ulps below the observed value) are
counted as at least as large; without this, single-area items obtain spuriously small P."""
import numpy as np, pandas as pd, time
A = np.load('A_af.npy'); D = np.load('D_af.npy'); R = pd.read_csv('mantel_af.csv')
dk = ['catheter','bag','ostomy','stoma','appliance','sheath','irrigation','faecal','drainage','adhesive','belt','tubing','protector','filler','plug','valve','filter','shield','plate','deodorant','lubricant','discharge solidifying','sensor','suspensory','collar','pouch','wipe','glove','dressing','bottle','teat','cup','syringe','needle','notification','urinal']
n = A.shape[0]; iu = np.triu_indices(n, 1); dvec = D[iu]; dc = dvec - dvec.mean(); dss = (dc**2).sum()
cand = np.where(R.P_Value.values < 0.05)[0]
X = A.T[cand].astype(float)                                  # (k, 42)
def rbatch(Xp):                                              # Xp (k, B, 42)
    M = np.abs(Xp[:, :, iu[0]] - Xp[:, :, iu[1]])
    Mc = M - M.mean(2, keepdims=True)
    num = Mc @ dc; den = np.sqrt((Mc**2).sum(2) * dss)
    with np.errstate(invalid='ignore', divide='ignore'): r = num / den
    return np.where(np.isfinite(r), r, 0.0)
r_obs = rbatch(X[:, None, :])[:, 0]
TOL = 1e-12; NP, B = 99999, 50; rng = np.random.default_rng(20261001); cnt = np.zeros(len(cand)); done = 0; t = time.time()
while done < NP:
    b = min(B, NP - done); P = np.array([rng.permutation(n) for _ in range(b)])
    cnt += (rbatch(X[:, P]) >= r_obs[:, None] - TOL).sum(1); done += b
p_ref = (cnt + 1) / (NP + 1)
print(f'{len(cand)} candidates refined with {NP} permutations in {time.time()-t:.0f}s')
p_all = R.P_Value.values.copy(); p_all[cand] = p_ref
m = len(p_all); o = np.argsort(p_all); q = np.empty(m)
q[o] = np.minimum.accumulate((p_all[o] * m / np.arange(1, m+1))[::-1])[::-1]; q = np.minimum(q, 1)
R['p_refined'] = p_all; R['q_BH'] = q
EXCLUDE = ['baguette', 'eye tear']                         # keyword false positives (food rolls, eye drops)
R['appliance'] = R.Drug.str.lower().apply(lambda s: any(k in s for k in dk) and not any(x in s for x in EXCLUDE))
R.to_csv('mantel_af_fdr.csv', index=False)
print('nominal P<0.05 with refined p:', int((p_all < 0.05).sum()))
for qt in [0.05, 0.10]:
    s = R[R.q_BH < qt]
    print(f'BH q<{qt}: {len(s)} items | appliances {int(s.appliance.sum())} | other {int((~s.appliance).sum())}')
print('\nnon-appliance items with q<0.10:')
print(R[(R.q_BH < 0.10) & (~R.appliance)].sort_values('q_BH')[['Drug','Mantel_r','p_refined','q_BH']].to_string(index=False))
print('\nnamed agents:')
for d in ['Chlortalidone','Ketoprofen','Hydrocortisone butyrate','Lisinopril','Thyrotropin alfa','Padimate O','Cladribine (Immunomodulating)','Catheters','Empagliflozin']:
    r = R[R.Drug.str.lower() == d.lower()]
    if len(r): print(f'  {d:30s} r={r.Mantel_r.iloc[0]:.3f} p_ref={r.p_refined.iloc[0]:.5f} q={r.q_BH.iloc[0]:.3f}')
lam = 0.5; pi0 = min(1, (p_all > lam).sum() / (m * (1 - lam)))
print(f'\nStorey pi0 (lambda=0.5) = {pi0:.3f} -> estimated share of items with a true Mantel signal {100*(1-pi0):.1f}%')
