import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Spatial weights for the 42 ICBs, aligned to the analysis order (codes_af.pkl).
Fuzzy queen contiguity (polygons buffered by 100 m; links are identical from 0 to 250 m, so no
gaps were introduced by simplification), row-standardised; k-nearest-neighbour (k=4) weights as a
sensitivity; centroid distances in km (British National Grid). Also writes the contiguity and k-nearest-neighbour
weights as tables (spatial_weights_contiguity.csv, spatial_weights_knn4.csv) and the ICB order (icb_order.csv)."""
import numpy as np, pandas as pd, geopandas as gpd, re, pickle
GEO = 'icb_boundaries_small.geojson'
FULL = 'icb_area_substance_items_full.csv'
gdf = gpd.read_file(GEO).to_crs(27700)
def norm(s): return re.sub('[^a-z]', '', str(s).lower().replace('nhs', '').replace('integrated care board', ''))
gdf['key'] = gdf['ICB23NM'].map(norm)
df = pd.read_csv(FULL, usecols=['AREA_CODE', 'AREA_NAME'])
icbname = df[df['AREA_NAME'].str.contains('INTEGRATED CARE', case=False, na=False)].groupby('AREA_CODE')['AREA_NAME'].first()
key2code = {norm(n): c for c, n in icbname.items()}
gdf['code'] = gdf['key'].map(key2code)
assert gdf['code'].notna().all()
codes = pickle.load(open('codes_af.pkl', 'rb'))
g = gdf.set_index('code').loc[codes].reset_index()
n = len(g)
buf = g.geometry.buffer(100)
Wb = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        if buf.iloc[i].intersects(buf.iloc[j]):
            Wb[i, j] = Wb[j, i] = 1
nn = Wb.sum(1)
print(f'contiguity: neighbours per ICB min {nn.min():.0f}, median {np.median(nn):.0f}, max {nn.max():.0f}; isolates {int((nn == 0).sum())}; links {int(Wb.sum() / 2)}')
cent = g.geometry.centroid
xy = np.c_[cent.x, cent.y]
Dgeo = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)) / 1000
k = 4; Wk = np.zeros((n, n))
for i in range(n):
    Wk[i, np.argsort(Dgeo[i])[1:k + 1]] = 1
np.save('Wb_contig.npy', Wb); np.save('W_contig.npy', Wb / Wb.sum(1, keepdims=True))
np.save('W_knn4.npy', Wk / Wk.sum(1, keepdims=True)); np.save('Dgeo_km.npy', Dgeo)
g[['code', 'ICB23CD', 'ICB23NM']].to_csv('icb_order.csv', index=False)
print(f'centroid distances: min {Dgeo[Dgeo > 0].min():.0f} km, median {np.median(Dgeo[np.triu_indices(n, 1)]):.0f} km, max {Dgeo.max():.0f} km')
nm = dict(zip(g.index, g['ICB23NM'].str.replace('NHS ', '').str.replace(' Integrated Care Board', '')))
for probe in ['Cornwall', 'Kent', 'North East and North Cumbria', 'North Central London']:
    i = [k for k, v in nm.items() if v.startswith(probe)][0]
    print(f'  {nm[i]}: neighbours = {[nm[j] for j in np.where(Wb[i])[0]]}')
# the neighbour sets are the same for every tolerance from 0 to 250 m (one more link appears at 300 m)
for tol in (0, 50, 150, 200, 250, 300):
    bt = g.geometry.buffer(tol) if tol else g.geometry
    i_, j_ = bt.sindex.query(bt, predicate='intersects'); Wt = np.zeros((n, n)); Wt[i_, j_] = 1; np.fill_diagonal(Wt, 0)
    print(f'tolerance {tol} m: {int(Wt.sum() / 2)} links; identical to 100 m: {bool((Wt == Wb).all())}')
# the weights as tables (binary; ODS codes; row-standardise before use)
pd.DataFrame(Wb.astype(int), index=codes, columns=codes).to_csv('spatial_weights_contiguity.csv')
pd.DataFrame(Wk.astype(int), index=codes, columns=codes).to_csv('spatial_weights_knn4.csv')
