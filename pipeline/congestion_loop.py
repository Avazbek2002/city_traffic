"""Equilibrium loop: gravity demand <-> shortest-path (Dijkstra) assignment with BPR congestion, MSA averaging.
usage: python3 congestion_loop.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz OUTDIR [scenario]
 scenario: am0809 (default, work-dominated hour; infographic congestion window) | am0708 (student-dominated hour)
Every driver takes the least-time path (Dijkstra) at the current congested times. Time units: minutes, flows: vehicles/hour.
All ASSUMPTION values can be overridden by env vars (BETA_WORK=0.08 ...)."""
import sys, os, time, numpy as np, pandas as pd, scipy.sparse as sp, numba as nb
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
EF, ZN, PA, OUT = sys.argv[1:5]; SC = sys.argv[5] if len(sys.argv) > 5 else 'am0809'
os.makedirs(OUT, exist_ok=True); E = lambda k, d: float(os.environ.get(k, d))
# ------------------------------------------------------------------ ASSUMPTIONS
BETA = {'work': E('BETA_WORK', 0.08), 'school': E('BETA_SCHOOL', 0.30), 'college': E('BETA_COLLEGE', 0.10), 'university': E('BETA_UNIVERSITY', 0.03)}  # 1/min decay on congested time
CAR = {'work': E('CAR_WORK', 0.45), 'school': 0.10, 'college': 0.02, 'university': 0.02}   # car share of person trips
OCC = {'work': E('OCC_WORK', 1.15), 'school': 1.3, 'college': 1.3, 'university': 1.3}      # persons per car
HOUR = {'am0809': {'work': 0.45, 'school': 0.15, 'college': 0.25, 'university': 0.40},     # share of the purpose's morning trips in the hour
        'am0708': {'work': 0.15, 'school': 0.75, 'college': 0.60, 'university': 0.50}}[SC]
CAP_LANE = {'motorway': 1800, 'trunk': 1100, 'primary': 1000, 'secondary': 900, 'tertiary': 800, 'unclassified': 600, 'residential': 500}  # veh/h/lane effective (signals included)
BPR_A, BPR_B = 0.15, 4.0
ITERS = int(E('ITERS', 15)); DROPOFF_RETURN = True      # school car trips: return leg (empty) also loads the network
INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0
t0 = time.time()
# ------------------------------------------------------------------ routing graph (no alleys)
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids)
e['a'] = ix[e.u].values; e['b'] = ix[e.v].values
e['t0'] = e.length_m / (e.speed_kmh / 3.6)                                  # s
LF = os.environ.get('LANES_FILE')          # xlsx from build_lane_checklist.py: column 'Typical lanes per direction (main part)' by street name
if LF:
    ov = pd.read_excel(LF, sheet_name='Streets'); ov = ov[ov['Typical lanes per direction (main part)'].notna()]
    m = dict(zip(ov['Street name (as in OpenStreetMap)'], ov['Typical lanes per direction (main part)'].astype(float)))
    st_ = e.names.fillna('').str.split('|').str[0]; sel = st_.isin(m) & (e.lanes_src == 'default')
    e.loc[sel, 'lanes'] = st_[sel].map(m); print(f'lane overrides: {len(m)} streets, {int(sel.sum())} road sections changed')
base = e.highway.str.replace('_link', '', regex=False)
e['cap'] = e.lanes * base.map(CAP_LANE).fillna(500) * np.where(e.highway.str.endswith('_link'), 0.8, 1.0)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
N = len(ids); ne = len(e)
indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr)
indices = e.b.values.astype(np.int32); eid = np.arange(ne, dtype=np.int32)
t0s = e.t0.values; cap = e.cap.values; length = e.length_m.values
G0 = sp.csr_matrix((t0s, indices, indptr), shape=(N, N))
nc, lab = connected_components(G0, directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; k = np.cos(np.radians(41.3))
d, j = cKDTree(np.c_[xy.lon.values[cand]*k, xy.lat.values[cand]]).query(np.c_[z.centroid_lon.values*k, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
print(f'routing graph: {N} nodes, {ne} edges ({e.length_m.sum()/1000:.0f} km), SCC {big.mean():.1%}; centroid snap median {np.median(d)*111320:.0f} m, max {d.max()*111320:.0f} m', flush=True)

@nb.njit(cache=True)
def accumulate(order, pred, indptr, indices, eid, acc, flow):
    for q in range(len(order)):
        x = order[q]; p = pred[x]
        if p < 0: continue
        a = acc[x]
        if a == 0.0: continue
        acc[p] += a
        for jj in range(indptr[p], indptr[p+1]):
            if indices[jj] == x:
                flow[eid[jj]] += a; break

def tt_from(t):                       # edge time (s) -> graph
    return sp.csr_matrix((t, indices, indptr), shape=(N, N))

def zone_times_and_assign(t, od=None):
    """Dijkstra from every zone node. Returns zone-zone times (min); if od given also AON edge flows."""
    G = tt_from(t); T = np.zeros((n, n)); flow = np.zeros(ne)
    for s in range(0, n, 80):
        ch = np.arange(s, min(s+80, n))
        dist, pred = dijkstra(G, directed=True, indices=zn[ch], return_predecessors=True)
        T[ch] = dist[:, zn]/60.0
        if od is not None:
            for r, o in enumerate(ch):
                acc = np.zeros(N); np.add.at(acc, zn, od[o])
                order = np.argsort(-np.where(np.isfinite(dist[r]), dist[r], -1.0)).astype(np.int32)
                accumulate(order, pred[r].astype(np.int32), indptr, indices, eid, acc, flow)
    return T, flow
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60

# ------------------------------------------------------------------ gravity per purpose
pa = np.load(PA); PUR = ['work', 'school', 'college', 'university']
P = {'work': z.population.values/z.population.sum()*z.employment.sum()}; A = {'work': z.employment.values.astype(float)}
for p in PUR[1:]: P[p] = pa[f'{p}_P']; A[p] = pa[f'{p}_A']
def ipf(K, Pv, Av, iters=400, tol=1e-6):
    r = np.ones(n); c = np.ones(n)
    for _ in range(iters):
        r = np.divide(Pv, K@c, out=np.zeros(n), where=(K@c) > 0); c = np.divide(Av, K.T@r, out=np.zeros(n), where=(K.T@r) > 0)
        M = r[:, None]*K*c[None]
        if max(np.abs(M.sum(1)-Pv).max()/Pv.max(), np.abs(M.sum(0)-Av).max()/Av.max()) < tol: break
    return M
def demand(T):
    Tm = T.copy(); Tm[np.diag_indices(n)] = intra
    M = {p: ipf(np.exp(-BETA[p]*Tm), P[p], A[p]) for p in PUR}
    veh = np.zeros((n, n))
    for p in PUR:
        v = M[p]*CAR[p]/OCC[p]*HOUR[p]
        veh += v + (v.T if (p == 'school' and DROPOFF_RETURN) else 0)
    veh[np.diag_indices(n)] = 0                       # intrazonal trips do not load the graph
    return M, veh, Tm
def bpr(x): return t0s*(1 + BPR_A*(x/cap)**BPR_B)

# ------------------------------------------------------------------ MSA equilibrium
x = np.zeros(ne); t = t0s.copy(); log = []
for it in range(1, ITERS+1):
    T, _ = zone_times_and_assign(t)
    M, veh, Tm = demand(T)
    _, y = zone_times_and_assign(t, veh)
    xn = y if it == 1 else x + (y-x)/it
    gap = float((t*(x-y)).sum()/max((t*x).sum(), 1e-9)) if it > 1 else np.nan   # relative gap before the update
    x = xn; t = bpr(x)
    vht = float((x*t).sum()/3600); vkt = float((x*length).sum()/1000)
    mean_t = {p: float((M[p]*Tm).sum()/M[p].sum()) for p in PUR}
    row = dict(it=it, rel_gap=round(gap, 4) if gap == gap else None, veh_total=round(float(veh.sum())), VKT=round(vkt), VHT=round(vht), avg_speed_kmh=round(vkt/vht, 1),
               edges_vc_gt1_pct=round(100*float((x/cap > 1).mean()), 2), **{f'mean_min_{p}': round(v, 1) for p, v in mean_t.items()})
    log.append(row); print(row, f'[{time.time()-t0:.0f}s]', flush=True)
pd.DataFrame(log).to_csv(f'{OUT}/loop_log_{SC}.csv', index=False)
T, _ = zone_times_and_assign(t); T[np.diag_indices(n)] = intra
np.save(f'{OUT}/travel_time_congested_{SC}.npy', T.astype(np.float32))
np.savez_compressed(f'{OUT}/edge_flows_{SC}.npz', a=e.a.values, b=e.b.values, u=e.u.values, v=e.v.values, flow=x, cap=cap, t0=t0s, t=t, length=length, highway=e.highway.values, names=e.names.fillna('').values,
                    u_lon=e.u_lon.values, u_lat=e.u_lat.values, v_lon=e.v_lon.values, v_lat=e.v_lat.values)
print('free-flow mean zone-zone', round(float(np.mean(zone_times_and_assign(t0s)[0][~np.eye(n, dtype=bool)])), 1), 'min; congested', round(float(T[~np.eye(n, dtype=bool)].mean()), 1), 'min')
