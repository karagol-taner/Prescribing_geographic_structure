import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Figures 1 to 3 of the manuscript (need adjustment, access versus choice, LISA maps) and
Supplementary Figures S2 (spatial clustering versus phylogeographic signal) and S8 (placebo covariates).
Figure 1A uses the distance-band summaries that need_similarity.py writes to need_similarity.json (dispensing appliances excluded)."""
import numpy as np, pandas as pd, geopandas as gpd, pickle, json, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlesize': 10.5, 'axes.titleweight': 'bold'})
INK = '#2c3e50'; BLUE = '#2f6db3'; ORANGE = '#e08a1e'; GREY = '#95a5a6'; LIGHT = '#e6e6e6'; RED = '#c0392b'
OUT = os.environ.get('FIG_OUT', '../figures'); os.makedirs(OUT, exist_ok=True)
GEO = os.environ.get('GEO', 'icb_boundaries_small.geojson')
FULL = os.environ.get('FULL', 'icb_area_substance_items_full.csv')
codes = pickle.load(open('codes_af.pkl', 'rb')); n = len(codes)
R = pd.read_csv('need_adjusted_moran.csv'); M = pd.read_csv('moran_af.csv'); years = np.load('years_af.npy')

# ---------------- Figure 1: need ----------------
fig, ax = plt.subplots(2, 2, figsize=(11, 8.8))
NSJ = json.load(open('need_similarity.json'))['primary_excluding_appliances']   # written by need_similarity.py
q = NSJ['raw']['distance_quintile_km']
xt = ['Bordering'] + [f'{q[k]:.0f}-{q[k + 1]:.0f}' for k in range(5)]
for key, lab, col in [('raw', 'Unadjusted', INK), ('age_deprivation_adjusted', 'Adjusted for age and deprivation', BLUE),
                      ('need_adjusted', 'Adjusted for age, deprivation and morbidity', ORANGE)]:
    yv = [NSJ[key]['relative_distance_bordering']] + NSJ[key]['relative_distance_by_quintile']
    ax[0, 0].plot(range(6), yv, '-', marker='o', color=col, lw=2, ms=6, label=lab)
ax[0, 0].axhline(1, color=INK, lw=.7, ls=':'); ax[0, 0].set_xticks(range(6)); ax[0, 0].set_xticklabels(xt, fontsize=8.5)
ax[0, 0].set_xlabel('Area pairs: bordering, then by centroid distance (km)'); ax[0, 0].set_ylabel('Prescribing dissimilarity\n(relative to the mean of all pairs)')
ax[0, 0].set_title('A  Need adjustment removes most neighbour similarity', loc='left'); ax[0, 0].legend(frameon=False, fontsize=8, loc='lower right')
sc = ax[0, 1].scatter(R.moran_I_raw, R.I_primary, c=R.R2_primary.clip(0, 1), cmap='Blues', vmin=0, vmax=1, s=9, edgecolor='none', alpha=.85)
lim = [-0.3, 0.75]; ax[0, 1].plot(lim, lim, color=INK, lw=.8, ls=':'); ax[0, 1].set_xlim(lim); ax[0, 1].set_ylim(-0.45, 0.65)
for d, dx, dy in [('Metformin hydrochloride', -40, -18), ('Loperamide hydrochloride', -40, 40), ('Amitriptyline hydrochloride', -70, 20),
                  ('Insulin glargine', -95, 12), ('Mebeverine hydrochloride', -54, -34), ('Dapagliflozin', -78, -18), ('Lansoprazole', -86, 16)]:
    r = R[R.Drug == d].iloc[0]
    ax[0, 1].annotate(d.replace(' hydrochloride', ''), (r.moran_I_raw, r.I_primary), textcoords='offset points', xytext=(dx, dy), fontsize=7.5, color=INK,
                      arrowprops=dict(arrowstyle='-', color=GREY, lw=.6))
cb = fig.colorbar(sc, ax=ax[0, 1], fraction=.045, pad=.02); cb.set_label('Share of between-area variance\nexplained by need (R²)', fontsize=8); cb.ax.tick_params(labelsize=8)
ax[0, 1].set_xlabel("Moran's I, unadjusted"); ax[0, 1].set_ylabel("Moran's I, need-adjusted")
ax[0, 1].set_title('B  Clustering shrinks once need is removed', loc='left')
full = pd.read_csv(FULL, usecols=['DRUG', 'ITEMS']); full['DRUG'] = full.DRUG.astype(str).str.strip()
vol = R.Drug.str.strip().map(full.groupby('DRUG').ITEMS.sum())
tert = pd.qcut(vol.rank(method='first'), 3, labels=['Least used third', 'Middle third', 'Most used third'])
raw = [100 * R.clustered_raw[tert == t].mean() for t in tert.cat.categories]; adj = [100 * R.clustered_primary[tert == t].mean() for t in tert.cat.categories]
x = np.arange(3); w = .38
b1 = ax[1, 0].bar(x - w / 2, raw, w, color=GREY, edgecolor='white', label='Unadjusted'); b2 = ax[1, 0].bar(x + w / 2, adj, w, color=BLUE, edgecolor='white', label='Need-adjusted')
for bars in (b1, b2):
    for rr in bars: ax[1, 0].text(rr.get_x() + rr.get_width() / 2, rr.get_height() + 1, f'{rr.get_height():.0f}%', ha='center', fontsize=8.5)
ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels(tert.cat.categories); ax[1, 0].set_ylim(0, 60); ax[1, 0].set_ylabel('Items spatially clustered (%, q<0.05)')
ax[1, 0].legend(frameon=False, fontsize=8.5, loc='upper left'); ax[1, 0].set_title('C  Clustering before and after adjustment', loc='left')
top = vol.rank(ascending=False) <= 100
groups = [R.R2_primary[tert == t].clip(0, 1) for t in tert.cat.categories] + [R.R2_primary[top].clip(0, 1)]
bp = ax[1, 1].boxplot(groups, widths=.55, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.5))
for patch in bp['boxes']: patch.set_facecolor('#c9dcf2'); patch.set_edgecolor(BLUE)
ax[1, 1].set_xticks(range(1, 5)); ax[1, 1].set_xticklabels(['Least used\nthird', 'Middle\nthird', 'Most used\nthird', '100 most\ndispensed'], fontsize=8.5)
ax[1, 1].set_ylabel('Variance explained by need (R²)'); ax[1, 1].set_ylim(0, 1)
for k, g in enumerate(groups): ax[1, 1].text(k + 1, g.median() + .02, f'{g.median():.2f}', ha='center', va='bottom', fontsize=8, color=INK)
chance = 5 / (n - 1)                       # expected R2 of an OLS fit with five covariates and 42 areas under no association
ax[1, 1].axhline(chance, color=INK, lw=.8, ls=':'); ax[1, 1].text(1.5, chance + .012, 'Chance', ha='center', fontsize=7.5, color=INK)
ax[1, 1].set_title('D  Need explains more for common medicines', loc='left')
plt.tight_layout(); plt.savefig(f'{OUT}/fig1_need.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()

# ---------------- Figure 2: access versus choice, unadjusted and need-adjusted ----------------
cn = json.load(open('class_need.json'))['results']
panels = [('SGLT2 inhibitors', 'SGLT2 inhibitors: class total', 'Dapagliflozin share of SGLT2 inhibitors', 'dapagliflozin share'),
          ('GLP-1 and GIP/GLP-1 receptor agonists', 'GLP-1 and GIP/GLP-1 agonists: class total', 'Semaglutide share of GLP-1 class', 'semaglutide share'),
          ('Thiazide-type diuretics', 'Thiazide-type diuretics: class total', 'Indapamide share of indapamide plus bendroflumethiazide', 'indapamide vs bendroflumethiazide')]
fig, ax = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
yrs = [int(y) for y in years]
for a, (ttl, tot, cho, chlab) in zip(ax, panels):
    for key, col, lab in [(tot, BLUE, 'Total use'), (cho, ORANGE, chlab.capitalize())]:
        for adjk, ls in [('', '-'), ('_adj', '--')]:
            I = [cn[key][str(y)]['I' + adjk] for y in yrs]; p = [cn[key][str(y)]['p' + adjk] for y in yrs]
            a.plot(yrs, I, ls, color=col, lw=1.8)
            a.scatter(yrs, I, s=36, color=[col if pp < 0.05 else 'white' for pp in p], edgecolor=col, zorder=5, lw=1.4)
    a.axhline(0, color=INK, lw=.7, ls=':'); a.set_xticks(yrs); a.set_title(ttl, loc='left', fontsize=10)
    a.text(0.02, 0.97, f'Choice: {chlab}', transform=a.transAxes, fontsize=8, color=INK, va='top')
ax[0].set_ylabel("Moran's I (filled: P<0.05)")
handles = [Line2D([], [], color=BLUE, lw=2, label='Total use of the class'), Line2D([], [], color=ORANGE, lw=2, label='Choice of agent within the class'),
           Line2D([], [], color=INK, lw=1.6, ls='-', label='Unadjusted'), Line2D([], [], color=INK, lw=1.6, ls='--', label='Need-adjusted')]
fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.04))
plt.tight_layout(rect=(0, 0.06, 1, 1)); plt.savefig(f'{OUT}/fig2_access_choice.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()

# ---------------- Figure 3: LISA maps ----------------
order = pd.read_csv('icb_order.csv'); gdf = gpd.read_file(GEO).to_crs(27700)
gdf = gdf.set_index('ICB23CD').loc[order.ICB23CD].reset_index(); gdf['ods'] = order.code.values
L = pd.read_csv('lisa_2025.csv')
CAT = {'High-High': '#a50f15', 'Low-Low': '#08519c', 'High-Low': '#ee7d4f', 'Low-High': '#4a9fd8', 'Not significant': LIGHT}
mp = [('sglt2_use', 'A  SGLT2 inhibitor use, 2025'), ('sglt2_use_need_adjusted', 'B  SGLT2 inhibitor use, need-adjusted'),
      ('dapagliflozin_share', 'C  Dapagliflozin share of SGLT2 inhibitors, 2025'), ('dapagliflozin_share_need_adjusted', 'D  Dapagliflozin share, need-adjusted')]
fig, axes = plt.subplots(2, 2, figsize=(9.6, 11))
for a, (v, ttl) in zip(axes.flat, mp):
    d = L[L.variable == v].set_index('ods').loc[gdf.ods]
    g = gdf.copy(); g['col'] = d.category.map(CAT).values
    g.plot(color=g['col'], ax=a, edgecolor='0.6', linewidth=0.3)
    sig = g[(d.q < 0.05).values]
    if len(sig): sig.boundary.plot(ax=a, color='black', linewidth=1.6)
    a.set_title(ttl, fontsize=10, fontweight='bold', loc='left'); a.axis('off')
leg = [Patch(facecolor=c, edgecolor='0.6', label=k) for k, c in CAT.items()] + [Patch(facecolor='white', edgecolor='black', lw=1.6, label='Survives FDR control (q<0.05)')]
fig.legend(handles=leg, loc='lower center', ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.0))
plt.tight_layout(rect=(0, 0.06, 1, 1)); plt.savefig(f'{OUT}/fig3_lisa.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()

# ---------------- Supplementary Figure S2: spatial clustering vs phylogeographic signal ----------------
cl = M.clustered.astype(bool); app = M.appliance.astype(bool); mfdr = M.q_mantel < 0.05
fig, a = plt.subplots(figsize=(6.6, 5.2))
a.scatter(M.Mantel_r[~cl & ~app], M.moran_I[~cl & ~app], s=6, color=GREY, alpha=.35, edgecolor='none', label='Not clustered')
a.scatter(M.Mantel_r[cl & ~app], M.moran_I[cl & ~app], s=7, color=BLUE, alpha=.45, edgecolor='none', label="Clustered (Moran's I, q<0.05)")
a.scatter(M.Mantel_r[app], M.moran_I[app], s=18, marker='^', color=ORANGE, alpha=.85, edgecolor='none', label='Dispensing appliances')
a.scatter(M.Mantel_r[mfdr], M.moran_I[mfdr], s=55, facecolor='none', edgecolor=RED, lw=1.2, label='Mantel signal, q<0.05')
HALO = [pe.withStroke(linewidth=2.5, foreground='white')]      # keeps labels legible over points
for d, dx, dy, ha in [('Loperamide hydrochloride', 6, 0, 'left'), ('Mebeverine hydrochloride', 0, 6, 'center'), ('Catheters', 0, -12, 'center'),
                      ('Chlortalidone', 6, 2, 'left'), ('Dapagliflozin', 6, 2, 'left'), ('Ketoprofen', 6, -3, 'left')]:
    r = M[M.Drug == d].iloc[0]
    a.annotate(d.replace(' hydrochloride', ''), (r.Mantel_r, r.moran_I), textcoords='offset points', xytext=(dx, dy), ha=ha, fontsize=7.5, color=INK,
               path_effects=HALO)
r = M[M.Drug == 'Padimate O'].iloc[0]      # labelled above the dotted line, with a leader
a.annotate('Padimate O', (r.Mantel_r, r.moran_I), xytext=(0.745, 0.015), textcoords='data', ha='right', va='bottom', fontsize=7.5, color=INK,
           path_effects=HALO, arrowprops=dict(arrowstyle='-', color=GREY, lw=.6, shrinkA=1, shrinkB=3))
a.axhline(-1 / (n - 1), color=INK, lw=.7, ls=':'); a.set_xlabel('Phylogeographic signal (Mantel r)'); a.set_ylabel("Spatial clustering (Moran's I, unadjusted)")
a.legend(frameon=True, framealpha=0.95, edgecolor='none', facecolor='white', fontsize=7, loc='upper right', bbox_to_anchor=(1.0, 0.92))
plt.tight_layout(); plt.savefig(f'{OUT}/figS2_moran_mantel.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()
# ---------------- Supplementary Figure S8: placebo covariates ----------------
PS = json.load(open('placebo_summary.json')); PR = pd.read_csv('placebo_runs.csv'); PRn = PR[PR.variant == 'need']; real = PS['real']
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
for a, colk, lab, ttl in [(ax[0], 'clustered', 'Items clustered after FDR control', 'A  Item-level clustering'),
                          (ax[1], 'neighbour_pct', 'Neighbouring pairs more similar than other pairs (%)', 'B  Similarity of neighbouring areas')]:
    a.hist(PRn[colk], bins=20, color='#c9dcf2', edgecolor=BLUE, lw=.6)
    for key, c, txt in [('unadjusted', GREY, 'Unadjusted'), ('primary', ORANGE, 'Need-adjusted')]:
        a.axvline(real[key][colk], color=c, lw=2); a.text(real[key][colk], a.get_ylim()[1] * .97, ' ' + txt, color=c, fontsize=8, va='top',
                                                           ha='left' if key == 'primary' else 'right')
    a.set_xlabel(lab); a.set_ylabel('Placebo covariate sets'); a.set_title(ttl, loc='left')
plt.tight_layout(); plt.savefig(f'{OUT}/figS8_placebo.png', dpi=200, bbox_inches='tight', facecolor='white'); plt.close()
