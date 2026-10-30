"""Pack slice results into compact binary blobs for the HTML viewer.
usage: python3 prep_viewer_data.py slices_flow.npz slice_stats.csv DISTRICTS.geojson OUT.json"""
import sys, json, base64, numpy as np, pandas as pd
from shapely.geometry import shape
SF, ST, DG, OUT = sys.argv[1:5]
d = np.load(SF, allow_pickle=True); FL, TT = d['flow'], d['t']; NS = FL.shape[0]
hw = d['highway']; way = d['way_id']; ul, ua, vl, va = d['u_lon'], d['u_lat'], d['v_lon'], d['v_lat']; a, b = d['a'], d['b']; ln = d['length']; cap = d['cap']; t0 = d['t0']
ART = {'trunk','primary','secondary','tertiary','trunk_link','primary_link','secondary_link','tertiary_link'}
keep = np.isin(hw, list(ART)) | (FL.max(0) >= 40)
idx = np.where(keep)[0]
# chain edges of one way (same direction) in travel order, cut into chunks of <= 350 m
by = {}
for i in idx: by.setdefault(int(way[i]), []).append(int(i))
chunks = []
for w, es in by.items():
    nxt = {int(a[i]): i for i in es}; heads = set(int(b[i]) for i in es); seen = set()
    starts = [i for i in es if int(a[i]) not in heads] or [es[0]]
    for s0 in starts + es:
        if s0 in seen: continue
        cur = s0; run = []
        while cur is not None and cur not in seen:
            seen.add(cur); run.append(cur); cur = nxt.get(int(b[cur]))
        cl, acc = [], 0.0
        for i in run:
            cl.append(i); acc += ln[i]
            if acc >= 350: chunks.append(cl); cl, acc = [], 0.0
        if cl: chunks.append(cl)
nch = len(chunks); print('edges kept', len(idx), 'chunks', nch)
lon0, lat0 = 69.04, 41.08; Q = 1e-5
pts, off, flow_q, ratio_q, capc, lenc, t0c = [], [0], np.zeros((NS, nch), np.uint16), np.zeros((NS, nch), np.uint8), np.zeros(nch, np.uint16), np.zeros(nch, np.uint16), np.zeros(nch, np.uint16)
for k, c in enumerate(chunks):
    c = np.array(c); L = ln[c]; pl = [(ul[c[0]], ua[c[0]])] + [(vl[i], va[i]) for i in c]
    for x, y in pl: pts += [int(round((x-lon0)/Q)), int(round((y-lat0)/Q))]
    off.append(len(pts)//2)
    capc[k] = min(65535, int(round(np.average(cap[c], weights=L)))); lenc[k] = min(65535, int(L.sum())); t0c[k] = min(65535, int(round(t0[c].sum()*10)))
    for s in range(NS):
        flow_q[s, k] = min(65535, int(round(np.average(FL[s, c], weights=L))))
        ratio_q[s, k] = min(255, int(round(TT[s, c].sum()/max(t0[c].sum(), 1e-6)*20)))
pts = np.array(pts, np.uint16); off = np.array(off, np.uint32)
# district outlines (simplified)
g = json.load(open(DG)); rings = []
for f in g['features']:
    geom = shape(f['geometry']).simplify(0.0006)
    polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    for p in polys:
        rings.append(dict(d=int(f['properties']['d_id']), p=[[int(round((x-lon0)/Q)), int(round((y-lat0)/Q))] for x, y in p.exterior.coords]))
b64 = lambda arr: base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode()
stats = pd.read_csv(ST)
out = dict(nch=nch, ns=NS, lon0=lon0, lat0=lat0, q=Q, pts=b64(pts), off=b64(off), flow=b64(flow_q), ratio=b64(ratio_q), cap=b64(capc), len=b64(lenc), t0=b64(t0c), rings=rings,
           stats=stats.to_dict('records'))
json.dump(out, open(OUT, 'w'), separators=(',', ':'))
import os; print('json MB', round(os.path.getsize(OUT)/1e6, 2), '; points', len(pts)//2)
