import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Markdown for Supplementary Tables S1 (need covariates), S6 (items clustered after need adjustment),
S7 (BNF sections clustered after need adjustment) and S8 (robustness of the choice result), generated from
the analysis outputs so that the supplement is not hand-transcribed. Writes tableS1.md, tableS6.md, tableS7.md
and tableS8.md."""
import pandas as pd, numpy as np, pickle, os
FULL = os.environ.get('FULL', 'icb_area_substance_items_full.csv')
codes = pickle.load(open('codes_af.pkl', 'rb'))
C = pd.read_csv('icb_covariates.csv').set_index('ods').loc[codes]
qf = lambda q: '<0.001' if q < 0.001 else f'{q:.3f}'
rows = [('% of registered population aged 65 or over', 'age_65plus', '93468', 'Mean of 2023 to 2025'),
        ('% of registered population aged under 18', 'age_u18', '93468', 'Mean of 2023 to 2025'),
        ('% of registered population aged 0 to 4 (component analysis only)', 'age_0_4', '93468', 'Mean of 2023 to 2025'),
        ('% of registered population aged 5 to 14 (component analysis only)', 'age_5_14', '93468', 'Mean of 2023 to 2025'),
        ('% of registered population aged 15 or over (component analysis only)', 'age_15plus', '93468', 'Mean of 2023 and 2025'),
        ('% of registered population aged 75 or over (component analysis only)', 'age_75plus', '93468', 'Mean of 2023 to 2025'),
        ('% of registered population aged 85 or over (component analysis only)', 'age_85plus', '93468', 'Mean of 2023 to 2025'),
        ('IMD 2025 score', 'imd2025', '94240', '2025'),
        ('IMD 2019 score (29 ICBs only; not used in models)', 'imd2019', '93553', '2019')]
qnames = {'qof_dm': ('Diabetes (17+)', '241'), 'qof_hyp': ('Hypertension', '219'), 'qof_ast': ('Asthma (6+)', '90933'), 'qof_copd': ('COPD', '253'),
          'qof_chd': ('Coronary heart disease', '273'), 'qof_hf': ('Heart failure', '262'), 'qof_af': ('Atrial fibrillation', '280'), 'qof_ckd': ('Chronic kidney disease (18+)', '258'),
          'qof_dep': ('Depression (18+)', '848'), 'qof_stroke': ('Stroke and transient ischaemic attack', '212'), 'qof_dem': ('Dementia', '247'),
          'qof_smi': ('Serious mental illness', '90581'), 'qof_epi': ('Epilepsy (18+)', '224'), 'qof_cancer': ('Cancer', '276'), 'qof_osteo': ('Osteoporosis (50+)', '90443'),
          'qof_ra': ('Rheumatoid arthritis (16+)', '91269'), 'qof_obes': ('Obesity (18+)', '94136'), 'qof_pad': ('Peripheral arterial disease', '92590'),
          'qof_ld': ('Learning disability', '200'), 'qof_ndh': ('Non-diabetic hyperglycaemia (18+)', '93797'), 'qof_smok': ('Smoking (15+)', '91280')}
per = {'qof_dep': 'Mean of 2021/22, 2022/23 and 2024/25', 'qof_obes': '2023/24'}
for k, (nm, iid) in qnames.items():
    rows.append((f'QOF prevalence: {nm}', k, iid, per.get(k, 'Mean of 2021/22 to 2024/25')))
out = ['| Covariate | Fingertips indicator | Period | Mean (SD) | Range |', '|---|---|---|---|---|']
for lab, col, iid, period in rows:
    v = C[col].dropna()
    out.append(f'| {lab} | {iid} | {period} | {v.mean():.2f} ({v.std():.2f}) | {v.min():.2f} to {v.max():.2f} |')
open('tableS1.md', 'w').write('\n'.join(out) + '\n')
S6 = pd.read_csv('tableS6_items.csv')
out = ['| Item (BNF name) | BNF section | Items dispensed, 2021-2025 | Moran\'s I, unadjusted | Moran\'s I, need-adjusted (q) | R² of need model |', '|---|---|---|---|---|---|']
for _, r in S6.iterrows():
    out.append(f'| {r.Drug} | {r.sec} | {r.vol:,.0f} | {r.moran_I_raw:.2f} | {r.I_primary:.2f} ({qf(r.q_primary)}) | {r.R2_primary:.2f} |')
open('tableS6.md', 'w').write('\n'.join(out) + '\n')
S7 = pd.read_csv('tableS7_sections.csv')
# largest items of each section (items assigned to the section in which most of their items were dispensed)
full = pd.read_csv(FULL, usecols=['DRUG', 'SECTION', 'ITEMS']); full['DRUG'] = full.DRUG.astype(str).str.strip()
sec = full.groupby(['DRUG', 'SECTION']).ITEMS.sum().reset_index().sort_values('ITEMS').groupby('DRUG').last().SECTION
vol = full.groupby('DRUG').ITEMS.sum()
# readable names for item names that are abbreviated or misspelt in the source data
label = {'Fluorouracil (Sunscreen)': 'Fluorouracil', 'Simple': 'Simple linctus', 'Ciprofloxain/dexameth': 'Ciprofloxacin with dexamethasone',
       'Co-careldopa (Carbidopa/levodopa)': 'Co-careldopa', 'Co-beneldopa (Benserazide/levodopa)': 'Co-beneldopa', 'Influenza': 'Influenza vaccine'}


def largest(section):
    v = vol[sec[sec == section].index].sort_values(ascending=False); sh = v / v.sum()
    parts = [f'{label.get(d, d)} ({100 * x:.0f}%)' for d, x in sh.head(2).items() if x >= 0.05]
    return '; '.join(parts[:1] + [p[0].lower() + p[1:] if not p.startswith('RtS') else p for p in parts[1:]])


out = ['| BNF section | Items in section | Largest items (% of section) | Items dispensed, 2021-2025 | Moran\'s I, unadjusted (q) | Moran\'s I, need-adjusted (q) |', '|---|---|---|---|---|---|']
for _, r in S7.iterrows():
    assert abs(vol[sec[sec == r.section].index].sum() - r.items_dispensed) < 1
    out.append(f'| {r.section} | {int(r.n_items)} | {largest(r.section)} | {r.items_dispensed:,.0f} | {r.I_raw:.2f} ({qf(r.q_raw)}) | {r.I_adj:.2f} ({qf(r.q_adj)}) |')
open('tableS7.md', 'w').write('\n'.join(out) + '\n')
import json
R8 = json.load(open('class_robustness.json'))['dapagliflozin_share']
CN = json.load(open('class_need.json'))['results']['Dapagliflozin share of SGLT2 inhibitors']
for y in ['2023', '2024', '2025']:      # the main-analysis row is taken from class_need.py so that it matches Table 2
    R8[y]['class-specific, contiguity'] = [CN[y]['I_adj'], CN[y]['p_adj']]
labs = ['class-specific, contiguity', 'primary, contiguity', 'age and deprivation, contiguity', 'six components of all covariates, contiguity', 'class-specific, k-nearest neighbours']
names = {'class-specific, contiguity': 'Class-specific need model (main analysis)', 'primary, contiguity': 'Primary need model of the item-level analysis',
         'age and deprivation, contiguity': 'Age and deprivation only', 'six components of all covariates, contiguity': 'Six principal components of all 29 covariates',
         'class-specific, k-nearest neighbours': 'Class-specific model, nearest-neighbour weights (k=4)'}
pf = lambda p: '<0.001' if p < 0.001 else f'{p:.3f}'
out = ['| Need model and spatial weights | 2023 | 2024 | 2025 |', '|---|---|---|---|']
for l in labs:
    out.append(f'| {names[l]} | ' + ' | '.join(f'{R8[y][l][0]:.2f} ({pf(R8[y][l][1])})' for y in ['2023', '2024', '2025']) + ' |')
L = R8['leave_one_out_2025']
out.append(f'| Class-specific model, each area omitted in turn (range) | not assessed | not assessed | {L["I_min"]:.2f} to {L["I_max"]:.2f} (largest P {pf(L["p_max"])}) |')
open('tableS8.md', 'w').write('\n'.join(out) + '\n')
print(open('tableS8.md').read())
print(open('tableS1.md').read()); print(open('tableS6.md').read()[:600]); print(open('tableS7.md').read()[:600])
