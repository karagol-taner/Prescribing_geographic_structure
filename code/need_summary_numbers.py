import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Every number quoted in the paper that is not printed by another script: need-adjusted summaries by
volume third, covariate summary (Table S1), the most clustered items and sections after adjustment
(Tables S6 and S7), correlates of the class-level results, dispensing appliances, items used in every area,
BNF sections, dispensing volume per registered patient, and the remaining numbers of the Methods and
Results. Prints to standard output (saved as data/need_summary_numbers.txt) and writes the numbers used
in the text to need_summary_numbers.json."""
import numpy as np, pandas as pd, pickle, json, os
from decimal import Decimal, ROUND_HALF_UP
from fractions import Fraction


def r1(x, d=1):
    """Round half up, as in the paper (369/720 = 51.25% is reported as 51.3%); x may be an exact Fraction."""
    v = Decimal(x.numerator) / Decimal(x.denominator) if isinstance(x, Fraction) else Decimal(repr(float(x)))
    return float(v.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))


def pct(mask):
    """Exact percentage of True values (a Fraction, so that halves are rounded correctly by r1)."""
    mask = np.asarray(mask, dtype=bool); return Fraction(int(mask.sum()), len(mask)) * 100


FULL = os.environ.get('FULL', 'icb_area_substance_items_full.csv')
codes = pickle.load(open('codes_af.pkl', 'rb'))
R = pd.read_csv('need_adjusted_moran.csv'); M = pd.read_csv('moran_af.csv')
C = pd.read_csv('icb_covariates.csv').set_index('ods').loc[codes]
full = pd.read_csv(FULL, usecols=['DRUG', 'SECTION', 'ITEMS']); full['DRUG'] = full.DRUG.astype(str).str.strip()
vol = R.Drug.str.strip().map(full.groupby('DRUG').ITEMS.sum())
sec = full.groupby(['DRUG', 'SECTION']).ITEMS.sum().reset_index().sort_values('ITEMS').groupby('DRUG').last().SECTION
R['vol'] = vol; R['sec'] = R.Drug.str.strip().map(sec)
tert = pd.qcut(vol.rank(method='first'), 3, labels=['least', 'middle', 'most']); top = vol.rank(ascending=False) <= 100
print('== clustered by model: n, % of items, % of volume, top-100, by third (least/middle/most)')
for k in ['raw', 'age_deprivation', 'pca6', 'primary']:
    cl = R['clustered_raw'] if k == 'raw' else R[f'clustered_{k}']
    nom = M['nominal'] if k == 'raw' else R[f'nominal_{k}']
    print(f'{k:16s} clustered {cl.sum():4d} ({100 * cl.mean():.1f}%) | nominal {nom.sum():4d} ({100 * nom.mean():.1f}%) | volume {100 * vol[cl].sum() / vol.sum():.1f}% | top100 {int(cl[top].sum())} | thirds',
          [r1(pct(cl[tert == t])) for t in ['least', 'middle', 'most']])
for k in ['age_deprivation', 'pca6', 'primary']:
    print(f'median R2 {k}: all {R[f"R2_{k}"].median():.2f}; thirds', [r1(R[f"R2_{k}"][tert == t].median(), 2) for t in ['least', 'middle', 'most']], f'; top100 {R[f"R2_{k}"][top].median():.2f}')
print('primary: nominal by third', [r1(pct(R.nominal_primary[tert == t])) for t in ['least', 'middle', 'most']], '; top100 nominal', int(R.nominal_primary[top].sum()))
print('Spearman raw vs adjusted I: %.3f' % R[['moran_I_raw', 'I_primary']].corr(method='spearman').iloc[0, 1])
print('\n== Table S1: covariates (42 ICBs)')
lab = {'age_65plus': '% aged 65+', 'age_u18': '% aged under 18', 'age_85plus': '% aged 85+', 'imd2025': 'IMD 2025 score', 'imd2019': 'IMD 2019 score (29 ICBs)'}
for c in ['age_65plus', 'age_u18', 'age_85plus', 'imd2025', 'imd2019'] + [c for c in C.columns if c.startswith('qof_')]:
    v = C[c].dropna(); print(f'  {lab.get(c, c):26s} mean {v.mean():7.2f}  SD {v.std():6.2f}  min {v.min():7.2f}  max {v.max():7.2f}  r(age65) {np.corrcoef(C.loc[v.index, "age_65plus"], v)[0, 1]:+.2f}  r(IMD25) {np.corrcoef(C.loc[v.index, "imd2025"], v)[0, 1]:+.2f}')
print('\n== Table S6: 25 items most strongly clustered after need adjustment (q<0.05)')
S6 = R[R.clustered_primary].sort_values('I_primary', ascending=False).head(25)
for _, r in S6.iterrows():
    print(f'  {r.Drug[:60]:60s} | {r.sec[:45]:45s} | raw {r.moran_I_raw:.2f} | adj {r.I_primary:.2f} | q {r.q_primary:.3f} | R2 {r.R2_primary:.2f} | items {r.vol:,.0f}')
S6.to_csv('tableS6_items.csv', index=False)
print('\n== Table S7: BNF sections clustered after need adjustment')
A = pd.read_csv('section_access.csv'); A['clu'] = (A.q_adj < 0.05) & (A.I_adj > 1e-12)
for _, r in A[A.clu].sort_values('I_adj', ascending=False).iterrows():
    print(f'  {r.section[:60]:60s} | n {r.n_items:3d} | items {r.items_dispensed:14,.0f} | raw {r.I_raw:.2f} | adj {r.I_adj:.2f} | q {r.q_adj:.3f}')
A[A.clu].sort_values('I_adj', ascending=False).to_csv('tableS7_sections.csv', index=False)
big = A.sort_values('items_dispensed', ascending=False).head(10)
print('10 largest sections: raw clustered', int(((big.q_raw < .05) & (big.I_raw > 1e-12)).sum()), '; adjusted clustered', int(big.clu.sum()))
print(big[['section', 'items_dispensed', 'I_raw', 'q_raw', 'I_adj', 'q_adj']].round(3).to_string(index=False))
print('\n== thiazide choice correlates (2021, 2025)')
shares = np.load('shares_af.npy'); years = list(np.load('years_af.npy')); areas_td, drugs_td = pickle.load(open('td_index.pkl', 'rb'))
yi = {d.strip().lower(): j for j, d in enumerate(drugs_td)}
for y in (2021, 2025):
    s = shares[years.index(y)]; ind = s[:, yi['indapamide']] / (s[:, yi['indapamide']] + s[:, yi['bendroflumethiazide']])
    print(f'  {y}: r(indapamide share, % aged 65+) {np.corrcoef(ind, C.age_65plus)[0, 1]:+.2f}; r(., IMD 2025) {np.corrcoef(ind, C.imd2025)[0, 1]:+.2f}; r(., QOF hypertension) {np.corrcoef(ind, C.qof_hyp)[0, 1]:+.2f}')
    if y == 2025:
        tS = sum(s[:, yi[m]] for m in ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin'])
        print(f'  2025: r(SGLT2 total, QOF diabetes) {np.corrcoef(tS, C.qof_dm)[0, 1]:+.2f}; r(., IMD) {np.corrcoef(tS, C.imd2025)[0, 1]:+.2f}; r(., % 65+) {np.corrcoef(tS, C.age_65plus)[0, 1]:+.2f}')
        d = s[:, yi['dapagliflozin']] / tS
        print(f'  2025: r(dapagliflozin share, QOF HF/DM) {np.corrcoef(d, C.qof_hf / C.qof_dm)[0, 1]:+.2f}; r(., QOF CKD/DM) {np.corrcoef(d, C.qof_ckd / C.qof_dm)[0, 1]:+.2f}')

print('\n== appliances, items used everywhere, sections, dispensing volume per registered patient')
BY = os.environ.get('BY', 'icb_substance_items_by_year.csv')
app = M.appliance.astype(bool).values; A = np.load('A_af.npy')
v_all = M.Drug.str.strip().map(full.groupby('DRUG').ITEMS.sum()).fillna(0).values
print(f'appliances: {app.sum()} items, {100 * v_all[app].sum() / v_all.sum():.2f}% of items dispensed nationally')
# share of the overall Jensen-Shannon divergence (summed over all pairs) contributed by appliance items, unadjusted profiles
P = A + 1e-10; P = P / P.sum(1, keepdims=True); iu = np.triu_indices(len(P), 1); p, r = P[iu[0]], P[iu[1]]; mm = 0.5 * (p + r)
contrib = 0.5 * p * np.log(p / mm) + 0.5 * r * np.log(r / mm)
print(f'appliance share of summed JS divergence across pairs: {100 * contrib[:, app].sum() / contrib.sum():.1f}%')
every = (A > 0).all(0)
print(f'items dispensed in every area: {every.sum()} ({(every & ~app).sum()} excluding appliances); national share of items dispensed {100 * v_all[every].sum() / v_all.sum():.3f}%')
fsec = full.groupby(['DRUG', 'SECTION']).ITEMS.sum().reset_index().sort_values('ITEMS').groupby('DRUG').last().SECTION
ns = fsec.value_counts(); multi = ns[ns >= 2].index
print(f'BNF sections: {len(ns)} in all; {len(multi)} with at least two items, containing {int(ns[multi].sum())} items and {100 * vol[R.sec.isin(multi)].sum() / vol.sum():.1f}% of items dispensed')
pop = pd.read_csv('icb_populations.csv').set_index('ods').loc[codes][['pop2023', 'pop2024', 'pop2025']].mean(axis=1)
by = pd.read_csv(BY); tot25 = by[by.YEAR == 2025].groupby('AREA_CODE').ITEMS.sum().reindex(codes)
pp = tot25 / pop; nmz = pd.read_csv('icb_order.csv').set_index('code').ICB23NM.str.replace('NHS ', '').str.replace(' Integrated Care Board', '')
print(f'items per registered patient, 2025: mean across areas {pp.mean():.1f} (median {pp.median():.1f}, range {pp.min():.1f} to {pp.max():.1f}); highest: ' +
      '; '.join(f'{nmz[c]} {v:.1f}' for c, v in pp.sort_values(ascending=False).head(3).items()))
top3 = pp.sort_values(ascending=False).head(3)
Xv = np.column_stack([np.ones(len(codes)), C.age_65plus.values, C.imd2025.values]); bv = np.linalg.lstsq(Xv, pp.values, rcond=None)[0]
resid_pp = pd.Series(pp.values - Xv @ bv, index=codes)
print(f'items per registered patient, 2025, relative to expectation from % aged 65+ and IMD: West Yorkshire {resid_pp["QWO"]:+.1f}; largest residuals ' +
      '; '.join(f'{nmz[c]} {v:+.1f}' for c, v in resid_pp.sort_values(ascending=False).head(3).items()))
json.dump(dict(appliance_items=int(app.sum()), appliance_pct_items=float(100 * v_all[app].sum() / v_all.sum()),
               appliance_pct_divergence=float(100 * contrib[:, app].sum() / contrib.sum()),
               items_every_area=int(every.sum()), items_every_area_excl_appliances=int((every & ~app).sum()),
               pct_items_dispensed_in_items_every_area=float(100 * v_all[every].sum() / v_all.sum()),
               sections_all=int(len(ns)), sections_multi=int(len(multi)), section_items=int(ns[multi].sum()),
               section_pct_items_dispensed=float(100 * vol[R.sec.isin(multi)].sum() / vol.sum()),
               items_per_patient_2025=dict(mean=float(pp.mean()), median=float(pp.median()), min=float(pp.min()), max=float(pp.max()),
                                           top3=[[nmz[c], float(v)] for c, v in top3.items()],
                                           west_yorkshire_excess_over_age_imd_expectation=float(resid_pp["QWO"]),
                                           largest_excess_area=nmz[resid_pp.idxmax()])),
          open('need_summary_numbers.json', 'w'), indent=1)
print('\n== numbers quoted in the paper that no other script prints')
from need_adjustment_core import pcs
qpc, _ = pcs(C, [c for c in C.columns if c.startswith('qof_')], 2)
print(f'|r(QOF PC1, % aged 65+)| {abs(np.corrcoef(qpc[:, 0], C.age_65plus)[0, 1]):.2f}; |r(QOF PC2, IMD 2025)| {abs(np.corrcoef(qpc[:, 1], C.imd2025)[0, 1]):.2f}')
L = pd.read_csv('lisa_2025.csv'); blk = L[(L.variable == 'dapagliflozin_share') & (L.category == 'Low-Low')].ods
s25 = shares[years.index(2025)]; SG = ['dapagliflozin', 'empagliflozin', 'canagliflozin', 'ertugliflozin']
cls = pd.DataFrame({m: s25[:, yi[m]] for m in SG}, index=codes); cls = cls.div(cls.sum(axis=1), axis=0); inb = cls.index.isin(blk)
print(f'empagliflozin share of SGLT2 items, 2025: low-share block {cls.empagliflozin[inb].mean():.2f}, elsewhere {cls.empagliflozin[~inb].mean():.2f}')
byd = pd.read_csv(BY); byd['DRUG'] = byd.DRUG.astype(str).str.strip()
first = byd[byd.ITEMS > 0].groupby('DRUG').YEAR.min(); new23 = first.index[first == 2023]
print(f'items first dispensed in 2023: {len(new23)}, of which oral or enteral nutrition: {int(sec.reindex(new23).isin(["Oral nutrition", "Enteral nutrition"]).sum())}')
print(f'Mantel r: median {M.Mantel_r.median():.3f}, IQR {M.Mantel_r.quantile(.25):.3f} to {M.Mantel_r.quantile(.75):.3f}')
Pp = pd.read_csv('icb_populations.csv').set_index('ods')
print(f'registered population below 2023 by up to {100 * (1 - Pp.pop2021 / Pp.pop2023).max():.1f}% (2021) and {100 * (1 - Pp.pop2022 / Pp.pop2023).max():.1f}% (2022)')
print(f'non-appliance items used in every area: {(every & ~app).sum()}, {100 * v_all[every & ~app].sum() / v_all[~app].sum():.2f}% of non-appliance items dispensed')
PS = json.load(open('placebo_summary.json')); r_, pn = PS['real'], PS['placebo']['need']
for k in ['clustered', 'neighbour_pct']:
    d = r_['unadjusted'][k] - r_['primary'][k]; f = lambda x: 100 * (r_['unadjusted'][k] - x) / d
    print(f'placebo share of the real reduction ({k}): mean {f(pn[k]["mean"]):.0f}%, 90% range {f(pn[k]["p95"]):.0f}% to {f(pn[k]["p05"]):.0f}%')
strong = (M.Mantel_r > 0.6) & (M.p_mantel_1000 < 0.05); print(f'strongest Mantel signals: {strong.sum()}, Moran nominal among them: {int((strong & M.nominal).sum())}')
NSJ = json.load(open('need_similarity.json'))['primary_excluding_appliances']['need_adjusted']
print(f'clipped item-area cells, primary model: {100 * NSJ["clipped_cells"] / (int((~app).sum()) * len(codes)):.1f}%')
