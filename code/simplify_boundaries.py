import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""How data/icb_boundaries_small.geojson was made from the full-resolution ONS boundaries.
Input: the ONS 'Integrated Care Boards (April 2023) EN BFE' boundaries as GeoJSON in WGS84, as downloaded
by code/02_fetch_boundaries_and_maps.ipynb (icb_boundaries.geojson, about 52 MB; not included in the
repository). Coordinates are rounded to 3 decimal places (about 70 to 110 m), consecutive duplicate
vertices are removed, and polygons whose outer ring has fewer than four points after rounding (a few
very small islands) are dropped; the result is written as compact JSON. Spatial weights are not affected
by the simplification (spatial_weights.py: contiguity links are identical with 0 to 250 m tolerance).
With the input file used for the paper (SHA-256 86b8ee43371839c2af75be45bf2418d867ef10831d1a855c79aa8de3877c422a)
the output is byte-identical to data/icb_boundaries_small.geojson (SHA-256 8d4e11b4703cca3b100006f07a1c4f0ee0eb7771dab00197de04b2407c4c691d);
ONS may revise the boundary file. The output is written to icb_boundaries_small_rebuilt.geojson and
compared with data/icb_boundaries_small.geojson, which is not overwritten.
Usage: python code/simplify_boundaries.py [path to the full boundary file, absolute or relative to the
repository root; default data/icb_boundaries.geojson]"""
import hashlib, json, os, sys

SRC = os.path.abspath(os.path.join(_ROOT, sys.argv[1])) if len(sys.argv) > 1 else 'icb_boundaries.geojson'


def ring(r, nd=3):
    out = []
    for x, y in r:
        q = [round(x, nd), round(y, nd)]
        if not out or out[-1] != q:
            out.append(q)
    return out


def polygon(rings):
    rr = [ring(r) for r in rings]
    return rr if len(rr[0]) >= 4 else None


raw = open(SRC, 'rb').read()
print('input SHA-256', hashlib.sha256(raw).hexdigest())
gj = json.loads(raw)
for f in gj['features']:
    g = f['geometry']
    if g['type'] == 'Polygon':
        g['coordinates'] = polygon(g['coordinates'])
    else:
        g['coordinates'] = [p for p in (polygon(x) for x in g['coordinates']) if p is not None]
s = json.dumps(gj, separators=(',', ':')).encode()
open('icb_boundaries_small_rebuilt.geojson', 'wb').write(s)
print('output SHA-256', hashlib.sha256(s).hexdigest(), f'({len(gj["features"])} areas)')
print('identical to data/icb_boundaries_small.geojson:', s == open('icb_boundaries_small.geojson', 'rb').read())
