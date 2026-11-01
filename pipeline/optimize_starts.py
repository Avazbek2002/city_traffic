"""Assign each WORKPLACE mahalla one work start time (08:30 / 09:00 / 09:30 / 10:00) so that the morning load on the roads is spread out as well as possible.
usage: python3 optimize_starts.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz BASE_slices_flow.npz OUT_classes.csv
Run with the same environment as the base run that produced BASE_slices_flow.npz (SIGNALS_FILE, SIGNAL_DELAY, CROSSING_DELAY, CAR_WORK, BETA_*, LANES_FILE).
How it works (a fast proxy, the real test is a full congestion_slices.py run with WORK_PROFILE=zones WORK_CLASS_FILE=OUT_classes.csv):
 1. Footprint: at the congested times of the 08:45 slice, route all work cars going TO each workplace mahalla over the shortest paths and record how much of each road they use (F[d, edge]).
 2. A start time moves a workplace's cars to different 15-minute slices (same departure rule as congestion_slices.py). The load of a road in a slice = base-run load - base work cars + re-timed work cars.
 3. Local search over workplaces: move one workplace to another start time when it lowers the BPR vehicle-hours of the whole morning; every start time keeps 20-30% of the jobs.
ASSUMPTION: routes of a workplace's cars stay as at 08:45 (no rerouting inside the search), so the proxy is only used to pick the assignment; the numbers come from the full run."""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, numba as nb
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
EF, ZN, PA, SF, OUT = sys.argv[1:6]
E = lambda k, d: float(os.environ.get(k, d))
BETA = {'work': E('BETA_WORK', 0.08)}; CAR, OCC = E('CAR_WORK', 0.45), E('OCC_WORK', 1.15)
CAP_LANE = {'motorway': 1800, 'trunk': 1100, 'primary': 1000, 'secondary': 900, 'tertiary': 800, 'unclassified': 600, 'residential': 500}
BPR_A, BPR_B = 0.15, 4.0; INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0; SEED = int(E('SEED', 1)); LO, HI = E('SHARE_LO', 0.20), E('SHARE_HI', 0.30)
STARTS = [510, 540, 570, 600]
W0 = [.01,.01,.03,.05, .02,.03,.04,.06, .08,.11,.12,.14, .07,.06,.04,.03]
# ---- graph exactly as in congestion_slices.py
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids)
e['a'] = ix[e.u].values; e['b'] = ix[e.v].values; e['t0'] = e.length_m / (e.speed_kmh / 3.6)
SG = os.environ.get('SIGNALS_FILE')
if SG:
    sg = pd.read_csv(SG); dj, dc = E('SIGNAL_DELAY', 20), E('CROSSING_DELAY', 5)
    isj = e.v.isin(set(sg.node[sg.type == 'junction'])); isc = e.v.isin(set(sg.node[sg.type == 'crossing'])) & ~isj
    e['t0'] = e.t0 + np.where(isj, dj, np.where(isc, dc, 0.0))
LF = os.environ.get('LANES_FILE')
if LF: raise SystemExit('LANES_FILE not supported here: copy the override block from congestion_slices.py')
base = e.highway.str.replace('_link', '', regex=False)
e['cap'] = e.lanes * base.map(CAP_LANE).fillna(500) * np.where(e.highway.str.endswith('_link'), 0.8, 1.0)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
N = len(ids); ne = len(e); S = np.load(SF, allow_pickle=True); FL, TT = S['flow'].astype(np.float64), S['t'].astype(np.float64); NS = FL.shape[0]
assert ne == FL.shape[1] and (e.a.values == S['a']).all() and (e.b.values == S['b']).all(), 'graph does not match the base run (different environment?)'
t0s, cap = e.t0.values, e.cap.values
indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr); indices = e.b.values.astype(np.int32)
G0 = sp.csr_matrix((t0s, indices, indptr), shape=(N, N)); nc, lab = connected_components(G0, directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; kc = np.cos(np.radians(41.3))
_, j = cKDTree(np.c_[xy.lon.values[cand]*kc, xy.lat.values[cand]]).query(np.c_[z.centroid_lon.values*kc, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60
Pw = z.population.values/z.population.sum()*z.employment.sum(); Aw = z.employment.values.astype(float)
def ipf(K, Pv, Av, iters=400, tol=1e-6):
    r = np.ones(n); c = np.ones(n)
    for _ in range(iters):
        r = np.divide(Pv, K@c, out=np.zeros(n), where=(K@c) > 0); c = np.divide(Av, K.T@r, out=np.zeros(n), where=(K.T@r) > 0)
        M = r[:, None]*K*c[None]
        if max(np.abs(M.sum(1)-Pv).max()/Pv.max(), np.abs(M.sum(0)-Av).max()/Av.max()) < tol: break
    return M
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
# ---- 1. footprints at the 08:45 slice (index 11)
SL = 11; t = TT[SL]
Gf = sp.csr_matrix((t, indices, indptr), shape=(N, N)); T = np.zeros((n, n))
for s in range(0, n, 80):
    ch = np.arange(s, min(s+80, n)); T[ch] = dijkstra(Gf, directed=True, indices=zn[ch])[:, zn]/60.0
Tm = T.copy(); Tm[np.diag_indices(n)] = intra
Wv = ipf(np.exp(-BETA['work']*Tm), Pw, Aw)*CAR/OCC; Wv[np.diag_indices(n)] = 0       # cars per whole morning, origin x destination (intrazonal do not load the graph)
order_r = np.lexsort((e.a.values, e.b.values)); ip_r = np.zeros(N+1, np.int64); np.add.at(ip_r, e.b.values+1, 1); ip_r = np.cumsum(ip_r)
idx_r = e.a.values[order_r].astype(np.int32); eid_r = order_r.astype(np.int32)
Gr = sp.csr_matrix((t[order_r], idx_r, ip_r), shape=(N, N))
F = np.zeros((n, ne), np.float32)
for s in range(0, n, 80):
    ch = np.arange(s, min(s+80, n)); dist, pred = dijkstra(Gr, directed=True, indices=zn[ch], return_predecessors=True)
    for r, d_ in enumerate(ch):
        acc = np.zeros(N); np.add.at(acc, zn, Wv[:, d_])
        order = np.argsort(-np.where(np.isfinite(dist[r]), dist[r], -1.0)).astype(np.int32); fl = np.zeros(ne)
        accumulate(order, pred[r].astype(np.int32), ip_r, idx_r, eid_r, acc, fl); F[d_] = fl
    print(f'footprints {min(s+80, n)}/{n}', flush=True)
Ltot = F.sum(0).astype(np.float64)
# ---- 2. loads
def kern(S_): return np.array([max(0.0, min(360+15*j+15, S_-5) - max(360+15*j, S_-45))/40.0 for j in range(16)])
MASS = sum(W0[8:]); PRE = np.array(W0[:8] + [0.0]*8); W0a = np.array(W0)
WK = np.array([MASS*kern(S_) for S_ in STARTS])                                       # class x slice: share of a workplace's morning cars leaving in the slice
REST = FL - 4*W0a[:, None]*Ltot[None] + 4*PRE[:, None]*Ltot[None]                      # slice x edge: everything except the re-timed work cars (the early shift workers stay)
tt0 = t0s/3600.0
def cost(x, ee): return (tt0[ee]*(x + BPR_A*x**(BPR_B+1)/cap[ee]**BPR_B)).sum()
allE = np.arange(ne)
def total_cost(L):          # L: class x edge (sum of F over the workplaces in the class)
    X = np.maximum(REST + 4*(WK.T@L), 0); return sum(cost(X[s_], allE) for s_ in range(NS))
# ---- 3. local search
rng = np.random.default_rng(SEED); jobs = Aw/Aw.sum(); supp = [np.nonzero(F[d_])[0] for d_ in range(n)]
def random_assign():
    cls = np.zeros(n, int); perm = rng.permutation(n); cum = np.cumsum(jobs[perm]); cls[perm] = np.minimum(3, (cum/ (1+1e-9) * 4).astype(int)); return cls
def Lmat(cls): return np.array([F[cls == k].sum(0) for k in range(4)], dtype=np.float64)
uniform = np.tile(0.25*Ltot, (4, 1))
base_cost = sum(cost(np.maximum(FL[s_], 0), allE) for s_ in range(NS))
c_start9 = total_cost(np.array([0*Ltot, Ltot, 0*Ltot, 0*Ltot]))
c_stag = total_cost(uniform)
rand_costs = []; best = None
for rep in range(3):
    cls = random_assign(); L = Lmat(cls); rand_costs.append(total_cost(L))
    X = np.maximum(REST + 4*(WK.T@L), 0); share = np.array([jobs[cls == k].sum() for k in range(4)])
    for it in range(4):
        moved = 0
        for d_ in rng.permutation(n):
            a = cls[d_]; ee = supp[d_]
            if len(ee) == 0: continue
            f = F[d_, ee].astype(np.float64); cur = sum(cost(X[s_, ee], ee) for s_ in range(NS)); bestk, bestc = a, cur
            for k in range(4):
                if k == a or share[a]-jobs[d_] < LO or share[k]+jobs[d_] > HI: continue
                xn = np.maximum(X[:, ee] + 4*(WK[k]-WK[a])[:, None]*f[None], 0); c = sum(cost(xn[s_], ee) for s_ in range(NS))
                if c < bestc - 1e-9: bestk, bestc = k, c
            if bestk != a:
                X[:, ee] = np.maximum(X[:, ee] + 4*(WK[bestk]-WK[a])[:, None]*f[None], 0); share[a] -= jobs[d_]; share[bestk] += jobs[d_]; cls[d_] = bestk; moved += 1
        print(f'rep {rep} pass {it}: moved {moved}, shares {np.round(share, 3)}', flush=True)
        if moved == 0: break
    cfin = total_cost(Lmat(cls))
    if best is None or cfin < best[0]: best = (cfin, cls.copy())
rel = lambda c: f'{(c/base_cost-1)*100:+.1f}% vs base profile'
print(f'proxy vehicle-hours (relative to the current assumed profile): all 09:00 {rel(c_start9)} | 25% each at every workplace {rel(c_stag)} | random one-time-per-workplace {rel(np.mean(rand_costs))} | optimised by location {rel(best[0])}')
pd.DataFrame(dict(zone_row=np.arange(n), mahalla_id=z.mahalla_id.values, start_min=np.array(STARTS)[best[1]])).to_csv(OUT, index=False)
print('written', OUT, 'jobs share by start', {int(S_): round(float(jobs[best[1] == k].sum()), 3) for k, S_ in enumerate(STARTS)})
