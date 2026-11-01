"""Zone graph data: how many people are in each mahalla / district over time, and the flows between them.
usage: python3 build_od_graph.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz slices_flow.npz DISTRICTS.geojson OUT.json
Run with the same CAR_WORK / BETA_* / LANES_FILE environment as congestion_slices.py (only the gravity part and the congested times are used).
Who is counted: every person who makes a morning trip in the model (all modes, not only drivers). At 06:00 they are at their home mahalla.
A trip leaves in its 15-minute slice (spread evenly over three 5-minute bins) and arrives after the congested car time between the two zones
(ASSUMPTION: walkers and public transport are assumed to take the same time as a car; they are really slower, so arrivals are early).
Trips inside one mahalla never change a node and are ignored. The timeline runs to 11:00 so late arrivals are shown."""
import sys, os, json, base64, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
from shapely.geometry import shape
EF, ZN, PA, SF, DG, OUT = sys.argv[1:7]
E = lambda k, d: float(os.environ.get(k, d))
BETA = {'work': E('BETA_WORK', 0.08), 'school': E('BETA_SCHOOL', 0.30), 'college': E('BETA_COLLEGE', 0.10), 'university': E('BETA_UNIVERSITY', 0.03)}
INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0
PUR = ['work', 'school', 'college', 'university']
PROFILE = {   # keep identical to congestion_slices.py (WORK_PROFILE env is honoured below)
 'work':       [.01,.01,.03,.05, .02,.03,.04,.06, .08,.11,.12,.14, .07,.06,.04,.03],
 'school':     [.005,.005,.01,.02, .10,.17,.23,.25, .08,.04,.02,.01, .01,.005,.003,.002],
 'college':    [.01,.01,.03,.04, .09,.14,.18,.19, .10,.07,.05,.03, .02,.01,.01,.01],
 'university': [.005,.005,.01,.03, .08,.11,.14,.17, .16,.12,.08,.04, .02,.01,.01,.01]}
def _shift(p, k): q = p[k:] + [0.0]*k; f = sum(p)/sum(q); return [v*f for v in q]
_wp = os.environ.get('WORK_PROFILE', 'base'); W0 = PROFILE['work']
if _wp == 'early30': PROFILE['work'] = _shift(W0, 2)
elif _wp == 'mix': PROFILE['work'] = [0.5*a + 0.5*b for a, b in zip(W0, _shift(W0, 4))]
elif _wp in ('start9', 'stagger'):     # same rule as congestion_slices.py: leave 45..5 min before the start; shift workers (slices 0-7 of the base profile) are kept
    def _kern(S): return [max(0.0, min(360+15*j+15, S-5) - max(360+15*j, S-45))/40.0 for j in range(16)]
    _mass = sum(W0[8:]); _early = W0[:8] + [0.0]*8
    _ws = [(540, 1.0)] if _wp == 'start9' else [(510, .25), (540, .25), (570, .25), (600, .25)]
    PROFILE['work'] = [a + _mass*sum(w*_kern(S)[j] for S, w in _ws) for j, a in enumerate(_early)]
# ---- graph exactly as in congestion_slices.py, so edge order matches the saved slice times
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids)
e['a'] = ix[e.u].values; e['b'] = ix[e.v].values; e['t0'] = e.length_m / (e.speed_kmh / 3.6)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
S = np.load(SF, allow_pickle=True); TT = S['t']; NS = TT.shape[0]
assert len(e) == TT.shape[1] and (e.a.values == S['a']).all() and (e.b.values == S['b']).all(), 'graph does not match slices_flow.npz (different inputs or LANES_FILE?)'
N = len(ids); indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr); indices = e.b.values.astype(np.int32)
G = lambda t: sp.csr_matrix((t, indices, indptr), shape=(N, N))
nc, lab = connected_components(G(e.t0.values), directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; kc = np.cos(np.radians(41.3))
_, j = cKDTree(np.c_[xy.lon.values[cand]*kc, xy.lat.values[cand]]).query(np.c_[z.centroid_lon.values*kc, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60
pa = np.load(PA); P = {'work': z.population.values/z.population.sum()*z.employment.sum()}; A = {'work': z.employment.values.astype(float)}
for p in PUR[1:]: P[p] = pa[f'{p}_P']; A[p] = pa[f'{p}_A']
def ipf(Kk, Pv, Av, iters=400, tol=1e-6):
    r = np.ones(n); c = np.ones(n)
    for _ in range(iters):
        r = np.divide(Pv, Kk@c, out=np.zeros(n), where=(Kk@c) > 0); c = np.divide(Av, Kk.T@r, out=np.zeros(n), where=(Kk.T@r) > 0)
        M = r[:, None]*Kk*c[None]
        if max(np.abs(M.sum(1)-Pv).max()/Pv.max(), np.abs(M.sum(0)-Av).max()/Av.max()) < tol: break
    return M
def zone_times(t):
    Gt = G(t); T = np.zeros((n, n))
    for s in range(0, n, 80):
        ch = np.arange(s, min(s+80, n)); T[ch] = dijkstra(Gt, directed=True, indices=zn[ch])[:, zn]/60.0
    return T
# ---- people in each zone at 5-minute boundaries, 06:00..11:00 (61 boundaries)
NB = 60; BIN = 5.0
leave = np.zeros((NB, n)); arrive = np.zeros((NB, n)); base = np.zeros(n)
dep_slice = []          # per slice: person departures o->d (n x n), for the arcs
for s_ in range(NS):
    T = zone_times(TT[s_].astype(np.float64)); Tm = T.copy(); Tm[np.diag_indices(n)] = intra
    pers = np.zeros((n, n))
    for p in PUR:
        M = ipf(np.exp(-BETA[p]*Tm), P[p], A[p]); pers += M*PROFILE[p][s_]
        if s_ == 0: base += M.sum(1)              # everybody who will travel, at home at 06:00 (row sums = production, fixed by the gravity constraints)
    np.fill_diagonal(pers, 0); dep_slice.append(pers)
    off = np.rint(T/BIN).astype(int)
    for jb in range(3):
        b = 3*s_ + jb; part = pers/3.0; leave[b] += part.sum(1)
        for k in np.unique(off[pers > 0]):
            m = pers*(off == k)/3.0; a = min(NB-1, b + int(k)); arrive[a] += m.sum(0)
    print(f'slice {s_:2d}: {pers.sum():,.0f} people on the move', flush=True)
present = np.zeros((NB+1, n)); present[0] = base
for b in range(NB): present[b+1] = present[b] - leave[b] + arrive[b]
transit = np.r_[0, np.cumsum(leave.sum(1) - arrive.sum(1))]; departed = np.r_[0, np.cumsum(leave.sum(1))]
# ---- district level
dids = sorted(z.district_id.unique()); dmap = {d: i for i, d in enumerate(dids)}; zd = z.district_id.map(dmap).values
agg = np.zeros((n, len(dids))); agg[np.arange(n), zd] = 1
dpresent = present@agg; dbase = base@agg
dflow = np.array([agg.T@d@agg for d in dep_slice]); dflow[:, np.arange(len(dids)), np.arange(len(dids))] = 0    # slice x from x to (persons leaving in that 15 min)
w = z.population.values.astype(float); dcx = (agg.T@(z.centroid_lon.values*w))/(agg.T@w); dcy = (agg.T@(z.centroid_lat.values*w))/(agg.T@w)
# ---- strongest mahalla-to-mahalla links per slice
TOPN = int(E('TOPN', 300)); links = []
for pers in dep_slice:
    f = pers.ravel(); idx = np.argpartition(-f, TOPN)[:TOPN]; idx = idx[np.argsort(-f[idx])]
    links.append([[int(i // n), int(i % n), round(float(f[i]), 1)] for i in idx if f[i] > 0])
lon0, lat0 = 69.04, 41.08
g = json.load(open(DG)); rings = []
for f_ in g['features']:
    geom = shape(f_['geometry']).simplify(0.0006)
    for pg in ([geom] if geom.geom_type == 'Polygon' else list(geom.geoms)):
        rings.append(dict(d=int(f_['properties']['d_id']), p=[[round(x, 4), round(y, 4)] for x, y in pg.exterior.coords]))
r1 = lambda a: np.round(a, 1).tolist()
out = dict(nb=NB, bin=BIN, rings=rings, ns=NS,
           transit=r1(transit), departed=r1(departed), people=float(base.sum()),
           mah=dict(name=z.mahalla_name.tolist(), d=z.district_id.tolist(), x=np.round(z.centroid_lon.values, 5).tolist(), y=np.round(z.centroid_lat.values, 5).tolist(), base=r1(base), present=[r1(r) for r in present], links=links),
           dist=dict(id=[int(d) for d in dids], x=np.round(dcx, 5).tolist(), y=np.round(dcy, 5).tolist(), base=r1(dbase), present=[r1(r) for r in dpresent], flow=[[[round(float(v), 1) for v in row] for row in m] for m in dflow]))
json.dump(out, open(OUT, 'w'), ensure_ascii=False, separators=(',', ':'))
print(f'people who travel: {base.sum():,.0f}; json {os.path.getsize(OUT)/1e6:.2f} MB')
print('biggest end-of-window gains (district):', [(int(d), int(round(a - b))) for d, a, b in sorted(zip(dids, dpresent[-1], dbase), key=lambda r: -(r[1]-r[2]))[:3]])
