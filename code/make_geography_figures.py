import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Supplementary Figure S1 (maps), plus the means, ranges and unadjusted Moran's I of the access-versus-choice
measures by year (key by_year in access_choice_unadjusted.json) and the Moran values of Table S3 (key tableS3_moran).
Permutation settings fixed: 9,999 permutations, fixed seeds."""
import numpy as np, pandas as pd, geopandas as gpd, pickle, json, os, shutil
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlesize': 10.5, 'axes.titleweight': 'bold'})
CT = '#2c3e50'; GENUINE = '#1a7a5a'; POINT = '#c0392b'; GREY = '#95a5a6'; BLUE = '#2f6db3'; ACC = '#e08a1e'
OUT = '../figures'; os.makedirs(OUT, exist_ok=True)
W = np.load('W_contig.npy'); A = np.load('A_af.npy'); M = pd.read_csv('moran_af.csv'); D = np.load('D_af.npy'); G = np.load('Dgeo_km.npy')
Wb = np.load('Wb_contig.npy'); shares = np.load('shares_af.npy'); years = np.load('years_af.npy')
areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb')); codes = pickle.load(open('codes_af.pkl', 'rb')); order = pd.read_csv('icb_order.csv')
n = W.shape[0]; pi = {d.strip().lower(): i for i, d in enumerate(M.Drug)}; yi = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
def moran(x, N=9999, seed=0):
    z = x - x.mean(); ss = (z ** 2).sum(); I = z @ W @ z / ss
    g = np.random.default_rng(seed); P = np.array([g.permutation(n) for _ in range(N)]); Zp = z[P]
    return float(I), float(((((Zp @ W.T) * Zp).sum(1) / ss) >= I - 1e-12).sum() + 1) / (N + 1)
cv = lambda v: float(v.std() / v.mean())
S = ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin']
GL = ['semaglutide', 'dulaglutide', 'liraglutide', 'exenatide', 'lixisenatide', 'tirzepatide']
TH = ['bendroflumethiazide', 'indapamide', 'chlortalidone', 'hydrochlorothiazide', 'metolazone']
def yseries(k):
    s = shares[k]; col = lambda m: s[:, yi[m]]
    tS = sum(col(m) for m in S if m in yi); tG = sum(col(m) for m in GL if m in yi); tT = sum(col(m) for m in TH if m in yi)
    sema = col('semaglutide') / tG
    return {'SGLT2 inhibitors: class total': tS, 'Dapagliflozin share of SGLT2 inhibitors': col('dapagliflozin') / tS,
            'GLP-1 and GIP/GLP-1 agonists: class total': tG, 'Semaglutide share of GLP-1 class': sema,
            'Thiazide-type diuretics: class total': tT,
            'Indapamide share of indapamide plus bendroflumethiazide': col('indapamide') / (col('indapamide') + col('bendroflumethiazide'))}
T5 = {}
for k, y in enumerate(years):
    for name, v in yseries(k).items():
        I, p = moran(v, seed=500 + k); T5.setdefault(name, {})[int(y)] = dict(I=I, p=p, cv=cv(v), mean=float(v.mean()), min=float(v.min()), max=float(v.max()))
T2 = {}
for d in ['Padimate O', 'Cladribine (Immunomodulating)', 'Chlortalidone', 'Thyrotropin alfa', 'Sucrase', 'Busulfan', 'Hydrocortisone butyrate',
          'Magnesium lactate', 'Haemophilus influenzae B/meningococcal C', 'Ketoprofen', 'Deferasirox', 'Paricalcitol', 'Butobarbital',
          'Amobarbital sodium', 'Methylprednisolone sodium succinate', 'Lisinopril']:
    r = M[M.Drug.str.lower() == d.lower()].iloc[0]; T2[d] = dict(I=float(r.moran_I), p=float(r.p_moran), q=float(r.q_moran))
json.dump(dict(by_year=T5, tableS3_moran=T2), open('access_choice_unadjusted.json', 'w'), indent=1)

# ---------------- Figure S1: maps ----------------
gdf = gpd.read_file('icb_boundaries_small.geojson').to_crs(27700)
gdf = gdf.set_index('ICB23CD').loc[order.ICB23CD].reset_index()           # now in analysis (codes) order
dapa25 = yseries(4)['Dapagliflozin share of SGLT2 inhibitors']
I_d25 = T5['Dapagliflozin share of SGLT2 inhibitors'][2025]['I']
panels = [(A[:, pi['loperamide hydrochloride']], f"Loperamide\nMoran's I = {M.loc[pi['loperamide hydrochloride'], 'moran_I']:.2f}: smooth regional gradient"),
          (dapa25, f"Dapagliflozin share of SGLT2 inhibitors, 2025\nMoran's I = {I_d25:.2f}: regional choice of agent"),
          (A[:, pi['chlortalidone']], f"Chlortalidone\nMoran's I = {M.loc[pi['chlortalidone'], 'moran_I']:.2f}: weak autocorrelation"),
          (A[:, pi['catheters']], f"Catheters\nMoran's I = {M.loc[pi['catheters'], 'moran_I']:.2f}: isolated dispensing hubs")]
nrm = TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=1.5)
fig, axes = plt.subplots(2, 2, figsize=(8.6, 10.4))
for ax, (v, ttl) in zip(axes.flat, panels):
    g = gdf.copy(); g['val'] = np.clip(v / v.mean(), 0.5, 1.5)
    g.plot(column='val', cmap='RdBu_r', norm=nrm, ax=ax, edgecolor='0.55', linewidth=0.25)
    ax.set_title(ttl, fontsize=9.5, fontweight='bold'); ax.axis('off')
fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.09, wspace=0.02, hspace=0.12)
cax = fig.add_axes([0.25, 0.055, 0.5, 0.018]); sm = ScalarMappable(norm=nrm, cmap='RdBu_r'); sm.set_array([])
cb = fig.colorbar(sm, cax=cax, orientation='horizontal', extend='both'); cb.set_ticks([0.5, 0.75, 1.0, 1.25, 1.5])
cb.set_label('Relative to the mean across areas (1.0 = average area)', fontsize=8.5); cb.ax.tick_params(labelsize=8)
plt.savefig(f'{OUT}/figS1_maps.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()

