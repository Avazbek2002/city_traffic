"""How many cars are on one street at the peak, which roads they come from and which roads they go to.
usage: python3 street_demand.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz slices_flow.npz OUT.json [--name REGEX] [--centre LAT,LON] [--radius M] [--slice N]
Run with the same environment as the congestion_slices.py run that made slices_flow.npz (CAR_WORK, OCC_WORK, BETA_*, WORK_PROFILE ...).
Defaults: the Mahtumquli street (name matched by 'mahtum|makhtum|махтум', without the minor "проезд" side passages), the stretch within 1500 m of 41.307206,69.302682.
How it works
 1. The street's road segments in the macro model are found by name. The "stretch" is those within --radius metres of --centre (the part the micro model draws).
 2. Cars on the stretch in a 15-minute slice = sum over its segments of (flow veh/h x congested travel time / 3600). The slice with the most cars is the peak (or use --slice).
 3. For that slice the gravity demand between every pair of mahallas (cars per hour, all purposes) is routed along its shortest path at the slice's congested times,
    exactly as in build_trips.py. For each path through the stretch we record the road it came from (the segment before the first stretch segment) and the road it goes to
    (the segment after the last one), with the junction where it joins / leaves. Summed over all pairs this gives "from road -> to road" flows in cars/hour.
 4. These all-or-nothing flows are scaled so that their total on the stretch equals the macro model's converged (averaged) flow, and the factor is printed.
ASSUMPTION: the macro model has no through traffic from outside the city, buses, taxis or trucks, so the numbers are cars of the modelled morning purposes only."""
import sys, os, json, re, argparse, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
ap = argparse.ArgumentParser(); ap.add_argument('edges'); ap.add_argument('zones'); ap.add_argument('pa'); ap.add_argument('slices'); ap.add_argument('out')
ap.add_argument('--name', default='mahtum|makhtum|махтум'); ap.add_argument('--exclude', default='проезд'); ap.add_argument('--centre', default='41.307206,69.302682')
ap.add_argument('--radius', type=float, default=1500); ap.add_argument('--slice', type=int, default=-1); ap.add_argument('--min-vph', type=float, default=0.5)
a = ap.parse_args(); EF, ZN, PA, SF, OUT = a.edges, a.zones, a.pa, a.slices, a.out
E = lambda k, d: float(os.environ.get(k, d))
BETA = {'work': E('BETA_WORK', 0.08), 'school': E('BETA_SCHOOL', 0.30), 'college': E('BETA_COLLEGE', 0.10), 'university': E('BETA_UNIVERSITY', 0.03)}
CAR = {'work': E('CAR_WORK', 0.45), 'school': 0.10, 'college': 0.02, 'university': 0.02}; OCC = {'work': E('OCC_WORK', 1.15), 'school': 1.3, 'college': 1.3, 'university': 1.3}
INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0; PUR = ['work', 'school', 'college', 'university']
PROFILE = {   # keep identical to congestion_slices.py (WORK_PROFILE env is honoured below)
 'work':       [.01,.01,.03,.05, .02,.03,.04,.06, .08,.11,.12,.14, .07,.06,.04,.03],
 'school':     [.005,.005,.01,.02, .10,.17,.23,.25, .08,.04,.02,.01, .01,.005,.003,.002],
 'college':    [.01,.01,.03,.04, .09,.14,.18,.19, .10,.07,.05,.03, .02,.01,.01,.01],
 'university': [.005,.005,.01,.03, .08,.11,.14,.17, .16,.12,.08,.04, .02,.01,.01,.01]}
def _shift(p, k): q = p[k:] + [0.0]*k; f = sum(p)/sum(q); return [v*f for v in q]
_wp = os.environ.get('WORK_PROFILE', 'base'); W0 = PROFILE['work']
if _wp == 'early30': PROFILE['work'] = _shift(W0, 2)
elif _wp == 'mix': PROFILE['work'] = [0.5*x + 0.5*y for x, y in zip(W0, _shift(W0, 4))]
elif _wp in ('start9', 'stagger'):
    def _kern(S): return [max(0.0, min(360+15*j+15, S-5) - max(360+15*j, S-45))/40.0 for j in range(16)]
    _mass = sum(W0[8:]); _early = W0[:8] + [0.0]*8; _ws = [(540, 1.0)] if _wp == 'start9' else [(510, .25), (540, .25), (570, .25), (600, .25)]
    PROFILE['work'] = [x + _mass*sum(w*_kern(S)[j] for S, w in _ws) for j, x in enumerate(_early)]
elif _wp != 'base': raise SystemExit('WORK_PROFILE zones is not supported here')
# ---- graph exactly as in congestion_slices.py / build_trips.py
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids)
e['a'] = ix[e.u].values; e['b'] = ix[e.v].values; e['t0'] = e.length_m / (e.speed_kmh / 3.6)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
S = np.load(SF, allow_pickle=True); TT = S['t'].astype(np.float64); FL = S['flow'].astype(np.float64); NS = TT.shape[0]
assert len(e) == TT.shape[1] and (e.a.values == S['a']).all() and (e.b.values == S['b']).all(), 'graph does not match slices_flow.npz (different inputs?)'
N = len(ids); ne = len(e)
indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr); indices = e.b.values.astype(np.int32)
ekey = e.a.values.astype(np.int64)*N + e.b.values; G = lambda t: sp.csr_matrix((t, indices, indptr), shape=(N, N))
nc, lab = connected_components(G(e.t0.values), directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
NLON, NLAT = xy.lon.values, xy.lat.values
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; kc = np.cos(np.radians(41.3))
_, j = cKDTree(np.c_[NLON[cand]*kc, NLAT[cand]]).query(np.c_[z.centroid_lon.values*kc, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60
pa = np.load(PA); P = {'work': z.population.values/z.population.sum()*z.employment.sum()}; Aw = {'work': z.employment.values.astype(float)}
for p in PUR[1:]: P[p] = pa[f'{p}_P']; Aw[p] = pa[f'{p}_A']
def ipf(Kk, Pv, Av, iters=400, tol=1e-6):
    r = np.ones(n); c = np.ones(n)
    for _ in range(iters):
        r = np.divide(Pv, Kk@c, out=np.zeros(n), where=(Kk@c) > 0); c = np.divide(Av, Kk.T@r, out=np.zeros(n), where=(Kk.T@r) > 0); M = r[:, None]*Kk*c[None]
        if max(np.abs(M.sum(1)-Pv).max()/Pv.max(), np.abs(M.sum(0)-Av).max()/Av.max()) < tol: break
    return M
def zone_times(t):
    Gt = G(t); T = np.zeros((n, n))
    for s in range(0, n, 80): ch = np.arange(s, min(s+80, n)); T[ch] = dijkstra(Gt, directed=True, indices=zn[ch])[:, zn]/60.0
    return T
# ---- 1. the street
lat0, lon0 = (float(v) for v in a.centre.split(','))
nm = e.names.fillna('').astype(str); isst = nm.str.contains(a.name, case=False, regex=True) & ~nm.str.contains(a.exclude, case=False, regex=True)
mx, my = (e.u_lon+e.v_lon)/2, (e.u_lat+e.v_lat)/2; dist = np.hypot((mx-lon0)*111320*kc, (my-lat0)*111320).values
stretch = (isst & (dist <= a.radius)).values; allst = isst.values
print(f'street segments: {allst.sum()} in the macro graph, {stretch.sum()} within {a.radius:.0f} m of the centre ({e.length_m[stretch].sum()/1000:.2f} km)')
du = np.c_[(e.v_lon-e.u_lon)*111320*kc, (e.v_lat-e.u_lat)*111320]; pts = np.c_[(mx[stretch]-lon0)*111320*kc, (my[stretch]-lat0)*111320]
axis = np.linalg.svd(pts - pts.mean(0))[2][0]; axis = axis if axis[1] >= 0 else -axis; fwd = (du@axis > 0)
bearing = (np.degrees(np.arctan2(axis[0], axis[1])) + 360) % 360
# ---- 2. cars on the stretch in every slice, and the peak
cars = (FL*TT/3600.0)[:, stretch]; carsF = cars[:, fwd[stretch]].sum(1); carsR = cars[:, ~fwd[stretch]].sum(1); tot = cars.sum(1)
cen = np.argmin(dist + np.where(stretch, 0, 1e9)); series = []
for s_ in range(NS):
    fl = lambda m_: FL[s_, stretch & m_]
    series.append(dict(slice=s_, clock=f'{6+s_//4:02d}:{15*(s_ % 4):02d}', cars_on_stretch=round(float(tot[s_]), 1), cars_forward=round(float(carsF[s_]), 1), cars_reverse=round(float(carsR[s_]), 1),
                       flow_forward_vph=round(float(np.median(fl(fwd))), 0) if (stretch & fwd).any() else 0, flow_reverse_vph=round(float(np.median(fl(~fwd))), 0) if (stretch & ~fwd).any() else 0,
                       speed_kmh=round(float(e.length_m[stretch].sum()/max(TT[s_, stretch].sum(), 1e-9)*3.6), 1)))
pk = int(np.argmax(tot)) if a.slice < 0 else a.slice
print('slice  time   cars on stretch (fwd / rev)   median flow veh/h fwd / rev   avg speed'); [print(f"  {r['slice']:2d}  {r['clock']}   {r['cars_on_stretch']:7.0f} ({r['cars_forward']:.0f} / {r['cars_reverse']:.0f})        {r['flow_forward_vph']:5.0f} / {r['flow_reverse_vph']:5.0f}      {r['speed_kmh']:.1f} km/h" + ('   <- peak' if r['slice'] == pk else '')) for r in series[3:14]]
# ---- 3. route-level demand in the peak slice
t = TT[pk]; T = zone_times(t); Tm = T.copy(); Tm[np.diag_indices(n)] = intra; D = np.zeros((n, n))
for p in PUR:
    M = ipf(np.exp(-BETA[p]*Tm), P[p], Aw[p]); v = M*CAR[p]/OCC[p]*PROFILE[p][pk]*4.0
    if p == 'school': v = v + v.T
    D += v
np.fill_diagonal(D, 0)
def road_name(i):
    names = [s.strip() for s in str(e.names.iat[i]).split('|') if s.strip() and s.strip() != 'nan']; lat = [s for s in names if re.search('[A-Za-z]', s)]
    return (lat[0] if lat else (names[0] if names else '(unnamed)')), names
flows = {}; aon = np.zeros(ne); multi = 0; used = 0; Gt = G(t)
for s0 in range(0, n, 80):
    ch = np.arange(s0, min(s0+80, n)); _, pred = dijkstra(Gt, directed=True, indices=zn[ch], return_predecessors=True)
    for r, o_ in enumerate(ch):
        for d_ in np.nonzero(D[o_] >= 0.02)[0]:
            sn, dn = zn[o_], zn[d_]
            if sn == dn or pred[r][dn] < 0: continue
            nodes = [dn]
            while nodes[-1] != sn: nodes.append(int(pred[r][nodes[-1]]))
            nodes = np.array(nodes[::-1]); eids = np.searchsorted(ekey, nodes[:-1].astype(np.int64)*N + nodes[1:]); w = D[o_, d_]; np.add.at(aon, eids, w)
            ins = stretch[eids]
            if not ins.any(): continue
            i0 = int(np.argmax(ins)); i1 = i0
            while i1 + 1 < len(eids) and ins[i1 + 1]: i1 += 1
            if ins[i1+1:].any(): multi += 1
            used += 1
            if i0 > 0: pe = int(eids[i0-1]); kin = 'street_end' if allst[pe] else 'side_road'; ink = (kin, int(nodes[i0]), road_name(pe)[0]); inf = (pe, nodes[i0])
            else: ink = ('local_origin', int(nodes[0]), '(trips that start on the street)'); inf = (None, nodes[0])
            if i1 < len(eids) - 1: qe = int(eids[i1+1]); kout = 'street_end' if allst[qe] else 'side_road'; outk = (kout, int(nodes[i1+1]), road_name(qe)[0]); outf = (qe, nodes[i1+1])
            else: outk = ('local_destination', int(nodes[-1]), '(trips that end on the street)'); outf = (None, nodes[-1])
            k = (ink, outk); fl = flows.setdefault(k, [0.0, inf, outf]); fl[0] += w
print(f'routed {used} OD pairs through the stretch ({multi} leave and re-join it)')
ratio = FL[pk, stretch].sum() / max(aon[stretch].sum(), 1e-9); print(f'all-or-nothing flow on the stretch is {aon[stretch].sum():,.0f} veh/h segment-sum vs {FL[pk, stretch].sum():,.0f} in the converged macro result; scaling factor {ratio:.2f}')
# ---- 4. export
def end(kind, node, name, edge, outgoing):
    d = dict(kind=kind, node=int(node), lon=round(float(NLON[node]), 6), lat=round(float(NLAT[node]), 6), road=name)
    if edge is not None:
        d['road_names'] = road_name(edge)[1]; d['highway'] = str(e.highway.iat[edge]); d['lanes'] = float(e.lanes.iat[edge])
        dx, dy = (e.v_lon.iat[edge]-e.u_lon.iat[edge])*kc, (e.v_lat.iat[edge]-e.u_lat.iat[edge]); d['bearing'] = round(float((np.degrees(np.arctan2(dx, dy)) + 360) % 360), 1)
    return d
od = []
for (ink, outk), (w, inf, outf) in sorted(flows.items(), key=lambda kv: -kv[1][0]):
    if w * ratio < a.min_vph: continue
    od.append(dict(vph=round(float(w*ratio), 1), src=end(*ink[:1], inf[1], ink[2], inf[0], False), dst=end(*outk[:1], outf[1], outk[2], outf[0], True)))
sumv = sum(o['vph'] for o in od); print(f'{len(od)} flows kept, {sumv:,.0f} cars/h in total entering the stretch at the peak')
agg = lambda key: sorted(((k, round(sum(o['vph'] for o in od if key(o) == k))) for k in {key(o) for o in od}), key=lambda x: -x[1])
print('\nWhere the cars come from (cars/h at the peak):')
for k, v in agg(lambda o: (o['src']['kind'], o['src']['road']))[:14]: print(f'  {v:5d}  {k[1]}   [{k[0]}]')
print('\nWhere they go:'); [print(f'  {v:5d}  {k[1]}   [{k[0]}]') for k, v in agg(lambda o: (o['dst']['kind'], o['dst']['road']))[:14]]
json.dump(dict(street=a.name, centre=[lat0, lon0], radius_m=a.radius, bearing_deg=round(float(bearing), 1), slice=pk, clock=f'{6+pk//4:02d}:{15*(pk % 4):02d}', scenario=os.environ.get('WORK_PROFILE', 'base'), car_share_work=CAR['work'],
                cars_on_stretch=round(float(tot[pk]), 1), scale=round(float(ratio), 3), series=series, od=od), open(OUT, 'w'), ensure_ascii=False, indent=1)
print('written', OUT)
