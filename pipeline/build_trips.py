"""Sample individual vehicle trips (origin zone -> destination zone, with their route) for the viewer's trip layer.
usage: python3 build_trips.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz slices_flow.npz OUT.json
Run it with the SAME environment variables (CAR_WORK, OCC_WORK, BETA_*, LANES_FILE) that were used for congestion_slices.py,
otherwise the sampled trips will not match the slice results. Extra knobs: SAMPLE (1 drawn car stands for SAMPLE real cars, default 80), SEED.
How it works: for every 15-minute slice the gravity OD matrices are rebuilt from that slice's congested travel times (saved in slices_flow.npz),
the cars departing in the slice are drawn from them (Poisson per OD cell), each gets a random departure minute inside the slice and
its Dijkstra shortest path at that slice's congested edge times. Intrazonal trips never load the graph, so they are not drawn."""
import sys, os, json, base64, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
EF, ZN, PA, SF, OUT = sys.argv[1:6]
E = lambda k, d: float(os.environ.get(k, d))
BETA = {'work': E('BETA_WORK', 0.08), 'school': E('BETA_SCHOOL', 0.30), 'college': E('BETA_COLLEGE', 0.10), 'university': E('BETA_UNIVERSITY', 0.03)}
CAR = {'work': E('CAR_WORK', 0.45), 'school': 0.10, 'college': 0.02, 'university': 0.02}
OCC = {'work': E('OCC_WORK', 1.15), 'school': 1.3, 'college': 1.3, 'university': 1.3}
K = E('SAMPLE', 80); rng = np.random.default_rng(int(E('SEED', 7)))
INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0
PUR = ['work', 'school', 'college', 'university']
PROFILE = {   # keep identical to congestion_slices.py
 'work':       [.01,.01,.03,.05, .02,.03,.04,.06, .08,.11,.12,.14, .07,.06,.04,.03],
 'school':     [.005,.005,.01,.02, .10,.17,.23,.25, .08,.04,.02,.01, .01,.005,.003,.002],
 'college':    [.01,.01,.03,.04, .09,.14,.18,.19, .10,.07,.05,.03, .02,.01,.01,.01],
 'university': [.005,.005,.01,.03, .08,.11,.14,.17, .16,.12,.08,.04, .02,.01,.01,.01]}
# ---- graph, built exactly as in congestion_slices.py so edge order matches the saved slices
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids)
e['a'] = ix[e.u].values; e['b'] = ix[e.v].values; e['t0'] = e.length_m / (e.speed_kmh / 3.6)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
S = np.load(SF, allow_pickle=True); TT = S['t']; NS = TT.shape[0]
assert len(e) == TT.shape[1] and (e.a.values == S['a']).all() and (e.b.values == S['b']).all(), 'graph does not match slices_flow.npz (different inputs or LANES_FILE?)'
N = len(ids); ne = len(e)
indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr); indices = e.b.values.astype(np.int32)
ekey = e.a.values.astype(np.int64)*N + e.b.values            # sorted ascending: edges are ordered by (a, b)
G = lambda t: sp.csr_matrix((t, indices, indptr), shape=(N, N))
nc, lab = connected_components(G(e.t0.values), directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
NLON, NLAT = xy.lon.values, xy.lat.values
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; kc = np.cos(np.radians(41.3))
_, j = cKDTree(np.c_[NLON[cand]*kc, NLAT[cand]]).query(np.c_[z.centroid_lon.values*kc, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60
# ---- gravity (same as congestion_slices.py)
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
        ch = np.arange(s, min(s+80, n)); dist = dijkstra(Gt, directed=True, indices=zn[ch]); T[ch] = dist[:, zn]/60.0
    return T
# ---- RDP simplification of a path (metres), keeps first/last vertex
def rdp(px, py, tol):
    m = len(px); keep = np.zeros(m, bool); keep[0] = keep[-1] = True; stack = [(0, m-1)]
    while stack:
        i, k = stack.pop()
        if k <= i+1: continue
        dx, dy = px[k]-px[i], py[k]-py[i]; L = np.hypot(dx, dy)
        d = np.hypot(px[i+1:k]-px[i], py[i+1:k]-py[i]) if L < 1e-9 else np.abs(dx*(py[i+1:k]-py[i]) - dy*(px[i+1:k]-px[i]))/L
        m_ = int(np.argmax(d))
        if d[m_] > tol: keep[i+1+m_] = True; stack += [(i, i+1+m_), (i+1+m_, k)]
    return np.where(keep)[0]
lon0, lat0, Q = 69.04, 41.08, 1e-5
X, Y, TQ, OFF, DEP, PUR_I, ORI, DST = [], [], [], [0], [], [], [], []
for s_ in range(NS):
    t = TT[s_].astype(np.float64); T = zone_times(t); Tm = T.copy(); Tm[np.diag_indices(n)] = intra
    trips = []   # (o, d, purpose)
    for pi, p in enumerate(PUR):
        M = ipf(np.exp(-BETA[p]*Tm), P[p], A[p]); v = M*CAR[p]/OCC[p]*PROFILE[p][s_]     # cars leaving in this 15-minute slice
        if p == 'school': v = v + v.T                                                   # empty return leg after drop-off
        np.fill_diagonal(v, 0); c = rng.poisson(v/K); oo, dd = np.nonzero(c)
        for o_, d_ in zip(np.repeat(oo, c[oo, dd]), np.repeat(dd, c[oo, dd])): trips.append((int(o_), int(d_), pi))
    by = {}
    for tr in trips: by.setdefault(tr[0], []).append(tr)
    Gt = G(t); origins = np.array(sorted(by)); made = 0
    for c0 in range(0, len(origins), 80):
        ch = origins[c0:c0+80]; _, pred = dijkstra(Gt, directed=True, indices=zn[ch], return_predecessors=True)
        for r, o_ in enumerate(ch):
            for (_, d_, pi) in by[o_]:
                s_node, d_node = zn[o_], zn[d_]
                if s_node == d_node or pred[r][d_node] < 0: continue
                nodes = [d_node]
                while nodes[-1] != s_node: nodes.append(int(pred[r][nodes[-1]]))
                nodes = np.array(nodes[::-1]); eids = np.searchsorted(ekey, nodes[:-1].astype(np.int64)*N + nodes[1:])
                cum = np.r_[0.0, np.cumsum(t[eids])]                                  # seconds since departure at each node
                lon, lat = NLON[nodes], NLAT[nodes]; keep = rdp((lon-lon[0])*kc*111320, (lat-lat[0])*111320, 4.0)
                X.append(np.round((lon[keep]-lon0)/Q).astype(np.int64)); Y.append(np.round((lat[keep]-lat0)/Q).astype(np.int64)); TQ.append(np.minimum(65535, np.round(cum[keep]*10)).astype(np.int64))
                OFF.append(OFF[-1] + len(keep)); DEP.append(int(round(s_*900 + rng.uniform(0, 900)))); PUR_I.append(pi); ORI.append(int(o_)); DST.append(int(d_)); made += 1
    print(f'slice {s_:2d} ({6+s_//4:02d}:{15*(s_%4):02d}): {len(trips)} cars drawn, {made} routed', flush=True)
x = np.concatenate(X).astype(np.uint16); y = np.concatenate(Y).astype(np.uint16); tq = np.concatenate(TQ).astype(np.uint16)
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()
out = dict(k=K, n=len(DEP), x=b64(x), y=b64(y), tq=b64(tq), off=b64(np.array(OFF, np.uint32)), dep=b64(np.array(DEP, np.uint16)), pur=b64(np.array(PUR_I, np.uint8)),
           o=b64(np.array(ORI, np.uint16)), d=b64(np.array(DST, np.uint16)),
           zones=[dict(n=r.mahalla_name, d=int(r.district_id), x=int(round((r.centroid_lon-lon0)/Q)), y=int(round((r.centroid_lat-lat0)/Q)), pop=int(r.population)) for r in z.itertuples()])
json.dump(out, open(OUT, 'w'), ensure_ascii=False, separators=(',', ':'))
print(f'{len(DEP)} sampled cars (1 drawn = {K:.0f} real), {len(x)} vertices, json {os.path.getsize(OUT)/1e6:.2f} MB')
