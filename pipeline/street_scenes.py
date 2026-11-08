"""Cars on EVERY road of a 1500 m x 1500 m window around the most congested part of each congested street, from the city (macro) model.
usage: python3 street_scenes.py road_edges.csv.gz mahalla_zones.csv pa_vectors.npz slices_flow.npz OUTDIR [--half 750] [--only id,id]
Run with the same environment as the congestion_slices.py run that made slices_flow.npz (CAR_WORK, OCC_WORK, BETA_* ...).
Writes OUTDIR/scene_<id>.demand.json for every scene in SCENES and OUTDIR/scenes_index.json.
How it works (same idea as street_demand.py, but for the whole window instead of one named street)
 1. Scene centre = the point of the street where the most vehicle-hours of delay fall inside a window of 2*half metres (all macro road segments of that street name).
 2. Peak slice = the 15-minute slice with the most cars on that street inside the window.
 3. For the peak slice, the gravity demand of every mahalla pair is routed along its shortest path at the slice's congested times (as in build_trips.py).
    Every stretch of a path inside the window becomes one trip: it enters on a road ~120 m before the window edge and leaves on a road ~120 m after it
    (trips that start or end inside the window start / end there). Trips are summed per (entry road, exit road).
 4. These all-or-nothing flows are scaled so that their total over the roads of the window equals the macro model's converged (averaged) flow on those roads.
ASSUMPTION: no through traffic from outside the city, buses, taxis or trucks; only the modelled morning purposes."""
import sys, os, json, re, argparse, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
ap = argparse.ArgumentParser(); [ap.add_argument(x) for x in ('edges', 'zones', 'pa', 'slices', 'outdir')]
ap.add_argument('--half', type=float, default=750); ap.add_argument('--only', default=''); ap.add_argument('--min-vph', type=float, default=0.15); ap.add_argument('--gate-m', type=float, default=120); ap.add_argument('--factors-only', action='store_true', help='only (re)write edge_factors into existing scene_*.demand.json')
a = ap.parse_args(); EF, ZN, PA, SF, OUT = a.edges, a.zones, a.pa, a.slices, a.outdir; os.makedirs(OUT, exist_ok=True); H = a.half
# id, label, regex on the street-name variants (case-insensitive)
SCENES = [('ahmad_donish', 'Ahmad Donish', r'ahmad donish|ахмад дониш'), ('bunyodkor', 'Bunyodkor shoh', r'bunyodkor shoh|бунёдкор шох'), ('qodiriy', 'Abdulla Qodiriy', r'abdulla qodiriy|абдулла кодирий|абдулла қодирий'),
          ('miskin', 'Miskin', r'miskin|мискин'), ('sebzor', 'Sebzor / Sabzor', r'sebzor|sabzor|себзор|сабзор'), ('furqat', 'Furqat', r'furqat|фуркат|фурқат'),
          ('ipakchi', 'Ipakchi', r'ipakchi|ипакчи'), ('mahtumquli', 'Mahtumquli (for comparison)', r'mahtum|makhtum|махтум')]
if a.only: SCENES = [s for s in SCENES if s[0] in a.only.split(',')]
E = lambda k, d: float(os.environ.get(k, d))
BETA = {'work': E('BETA_WORK', 0.08), 'school': E('BETA_SCHOOL', 0.30), 'college': E('BETA_COLLEGE', 0.10), 'university': E('BETA_UNIVERSITY', 0.03)}
CAR = {'work': E('CAR_WORK', 0.45), 'school': 0.10, 'college': 0.02, 'university': 0.02}; OCC = {'work': E('OCC_WORK', 1.15), 'school': 1.3, 'college': 1.3, 'university': 1.3}
INTRA_F, CIRC, SPEED = 0.5, 1.35, 25.0; PUR = ['work', 'school', 'college', 'university']
PROFILE = {'work': [.01,.01,.03,.05, .02,.03,.04,.06, .08,.11,.12,.14, .07,.06,.04,.03], 'school': [.005,.005,.01,.02, .10,.17,.23,.25, .08,.04,.02,.01, .01,.005,.003,.002],
           'college': [.01,.01,.03,.04, .09,.14,.18,.19, .10,.07,.05,.03, .02,.01,.01,.01], 'university': [.005,.005,.01,.03, .08,.11,.14,.17, .16,.12,.08,.04, .02,.01,.01,.01]}
assert os.environ.get('WORK_PROFILE', 'base') == 'base', 'only the base profile is supported here'
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids); e['a'] = ix[e.u].values; e['b'] = ix[e.v].values; e['t0'] = e.length_m / (e.speed_kmh / 3.6)
e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
S = np.load(SF, allow_pickle=True); TT = S['t'].astype(np.float64); FL = S['flow'].astype(np.float64); NS = TT.shape[0]
assert len(e) == TT.shape[1] and (e.a.values == S['a']).all() and (e.b.values == S['b']).all(), 'graph does not match slices_flow.npz'
N = len(ids); ne = len(e); LEN = e.length_m.values
indptr = np.zeros(N+1, np.int64); np.add.at(indptr, e.a.values+1, 1); indptr = np.cumsum(indptr); indices = e.b.values.astype(np.int32)
ekey = e.a.values.astype(np.int64)*N + e.b.values; G = lambda t: sp.csr_matrix((t, indices, indptr), shape=(N, N))
nc, lab = connected_components(G(e.t0.values), directed=True, connection='strong'); big = lab == np.argmax(np.bincount(lab))
xy = pd.concat([e[['a', 'u_lon', 'u_lat']].set_axis(['n', 'lon', 'lat'], axis=1), e[['b', 'v_lon', 'v_lat']].set_axis(['n', 'lon', 'lat'], axis=1)]).drop_duplicates('n').set_index('n').sort_index()
z = pd.read_csv(ZN); n = len(z); cand = np.where(big)[0]; kc = np.cos(np.radians(41.3))
_, j = cKDTree(np.c_[xy.lon.values[cand]*kc, xy.lat.values[cand]]).query(np.c_[z.centroid_lon.values*kc, z.centroid_lat.values]); zn = cand[j].astype(np.int32)
intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60
pa = np.load(PA); P = {'work': z.population.values/z.population.sum()*z.employment.sum()}; Aw = {'work': z.employment.values.astype(float)}
for p in PUR[1:]: P[p] = pa[f'{p}_P']; Aw[p] = pa[f'{p}_A']
def ipf(Kk, Pv, Av, iters=400, tol=1e-6):
    r = np.ones(n); c = np.ones(n)
    for _ in range(iters):
        r = np.divide(Pv, Kk@c, out=np.zeros(n), where=(Kk@c) > 0); c = np.divide(Av, Kk.T@r, out=np.zeros(n), where=(Kk.T@r) > 0); M = r[:, None]*Kk*c[None]
        if max(np.abs(M.sum(1)-Pv).max()/Pv.max(), np.abs(M.sum(0)-Av).max()/Av.max()) < tol: break
    return M
def demand(t):
    Gt = G(t); T = np.zeros((n, n))
    for s in range(0, n, 80): ch = np.arange(s, min(s+80, n)); T[ch] = dijkstra(Gt, directed=True, indices=zn[ch])[:, zn]/60.0
    Tm = T.copy(); Tm[np.diag_indices(n)] = intra; parts = {}
    for p in PUR:
        M = ipf(np.exp(-BETA[p]*Tm), P[p], Aw[p]); D_p = M*CAR[p]/OCC[p]*4.0         # cars/h per unit of profile share
        parts[p] = D_p + D_p.T if p == 'school' else D_p
    return parts
def road_names(i):
    names = [s.strip() for s in str(e.names.iat[i]).split('|') if s.strip() and s.strip() != 'nan']; lat = [s for s in names if re.search('[A-Za-z]', s)]
    return (lat[0] if lat else (names[0] if names else '(unnamed)')), names
STREET = np.array([road_names(i)[0] for i in range(ne)])
mx, my = (e.u_lon.values+e.v_lon.values)/2, (e.u_lat.values+e.v_lat.values)/2; LON0, LAT0 = 69.2, 41.3; MX, MY = (mx-LON0)*111320*kc, (my-LAT0)*111320
nm = e.names.fillna('').astype(str); tfree = TT.min(0); DEL = ((TT - tfree[None]) * FL * 0.25 / 3600).sum(0)
# ---- 1. scene centres and windows
sc = []
for sid, label, rx in SCENES:
    isst = nm.str.contains(rx, case=False, regex=True).values & ~nm.str.contains('проезд', case=False).values; idx = np.nonzero(isst)[0]
    if len(idx) == 0: print('no street found for', sid); continue
    best = (-1, None)
    for i in idx[np.argsort(-DEL[idx])][:400]:
        m = (np.abs(MX[idx]-MX[i]) <= H) & (np.abs(MY[idx]-MY[i]) <= H); s_ = DEL[idx][m].sum()
        if s_ > best[0]: best = (s_, i)
    c = best[1]; cx, cy = MX[c], MY[c]; win = (np.abs(MX-cx) <= H) & (np.abs(MY-cy) <= H); cars = (FL*TT/3600.0)[:, win & isst].sum(1); pk = int(np.argmax(cars))
    clon, clat = float(LON0 + cx/(111320*kc)), float(LAT0 + cy/111320)
    print(f'{sid:13s} street edges {len(idx):4d}  centre {clat:.5f},{clon:.5f}  delay in window {best[1] and best[0]:.0f} veh-h  peak slice {pk} ({6+pk//4:02d}:{15*(pk % 4):02d})  cars on street in window {cars[pk]:.0f}  window edges {win.sum()}')
    sc.append(dict(id=sid, label=label, rx=rx, isst=isst, win=win, pk=pk, centre=(clat, clon), cars=cars, flows={}, cx=cx, cy=cy, delay=float(best[0])))
anywin = np.any([s['win'] for s in sc], axis=0)
# ---- 2/3. route the peak-slice demand and cut it at the window edges
for pk in ([] if a.factors_only else sorted({s['pk'] for s in sc})):
    group = [s for s in sc if s['pk'] == pk]; t = TT[pk]; parts = demand(t)
    D = sum(parts[p] * PROFILE[p][pk] for p in PUR); np.fill_diagonal(D, 0)
    aon = np.zeros(ne); Gt = G(t); used = 0
    for s0 in range(0, n, 80):
        ch = np.arange(s0, min(s0+80, n)); _, pred = dijkstra(Gt, directed=True, indices=zn[ch], return_predecessors=True)
        for r, o_ in enumerate(ch):
            for d_ in np.nonzero(D[o_] >= 0.02)[0]:
                sn, dn = zn[o_], zn[d_]
                if sn == dn or pred[r][dn] < 0: continue
                nodes = [dn]
                while nodes[-1] != sn: nodes.append(int(pred[r][nodes[-1]]))
                nodes = np.array(nodes[::-1]); eids = np.searchsorted(ekey, nodes[:-1].astype(np.int64)*N + nodes[1:]); w = D[o_, d_]; np.add.at(aon, eids, w)
                if not anywin[eids].any(): continue
                used += 1
                for s in group:
                    ins = s['win'][eids]
                    if not ins.any(): continue
                    ii = np.flatnonzero(ins); runs = np.split(ii, np.flatnonzero(np.diff(ii) > 1)+1)
                    for run in runs:
                        g0 = int(run[0]); acc = 0.0
                        while g0 > 0 and acc < a.gate_m: g0 -= 1; acc += LEN[eids[g0]]
                        g1 = int(run[-1]); acc = 0.0
                        while g1 < len(eids)-1 and acc < a.gate_m: g1 += 1; acc += LEN[eids[g1]]
                        k = (int(eids[g0]), int(eids[g1])); s['flows'][k] = s['flows'].get(k, 0.0) + w
    print(f'slice {pk}: routed {used} OD pairs through at least one window')
    for s in group:
        ratio = FL[pk, s['win']].sum() / max(aon[s['win']].sum(), 1e-9); s['ratio'] = float(ratio); print(f"  {s['id']}: all-or-nothing {aon[s['win']].sum():,.0f} vs converged {FL[pk, s['win']].sum():,.0f} segment-sum veh/h; scale {ratio:.2f}; {len(s['flows'])} entry/exit pairs")
# ---- 4. export
def edge_info(i):
    dx, dy = (e.v_lon.iat[i]-e.u_lon.iat[i])*kc, (e.v_lat.iat[i]-e.u_lat.iat[i]); nmz, names = road_names(i)
    return dict(lon=round(float((e.u_lon.iat[i]+e.v_lon.iat[i])/2), 6), lat=round(float((e.u_lat.iat[i]+e.v_lat.iat[i])/2), 6), bearing=round(float((np.degrees(np.arctan2(dx, dy))+360) % 360), 1), road=nmz, road_names=names, highway=str(e.highway.iat[i]), lanes=float(e.lanes.iat[i]))
# congestion factor of every road segment: congested time (incl. signal delay) / driving time at the limit; the SUMO router uses it so that cars choose the macro model's routes
fac = TT / (LEN / (e.speed_kmh.values / 3.6))[None]
index = []
for s in sc:
    pk = s['pk']; win = s['win']
    if a.factors_only:
        f_ = f"{OUT}/scene_{s['id']}.demand.json"; out = json.load(open(f_)); near = np.nonzero((np.abs(MX-s['cx']) <= H+450) & (np.abs(MY-s['cy']) <= H+450))[0]
        dx_, dy_ = (e.v_lon.values[near]-e.u_lon.values[near])*kc, (e.v_lat.values[near]-e.u_lat.values[near])
        out['edge_factors'] = [[round(float(mx[i]), 6), round(float(my[i]), 6), round(float((np.degrees(np.arctan2(dx_[k], dy_[k]))+360) % 360)), round(float(fac[out['slice'], i]), 2)] for k, i in enumerate(near)]
        json.dump(out, open(f_, 'w'), ensure_ascii=False); print('factors', s['id'], len(near)); index.append({k: out[k] for k in ('id', 'label', 'centre', 'clock', 'cars_on_street', 'delay_vehicle_hours', 'entering_vph', 'window_total_cars')}); continue
    wi = np.nonzero(win)[0]; od = []
    for (ei, eo), w in sorted(s['flows'].items(), key=lambda kv: -kv[1]):
        if w*s['ratio'] < a.min_vph or ei == eo: continue
        od.append(dict(vph=round(float(w*s['ratio']), 2), src=edge_info(ei), dst=edge_info(eo)))
    # per street inside the window: what the macro model says at the peak
    df = pd.DataFrame(dict(street=STREET[wi], len=LEN[wi], flow=FL[pk, wi], cars=(FL*TT/3600.0)[pk, wi], spd=LEN[wi]/np.maximum(TT[pk, wi], 1e-9)*3.6))
    streets = []
    for nmx, g in df.groupby('street'):
        if g.len.sum() < 80: continue
        streets.append(dict(street=nmx, km=round(float(g.len.sum()/1000), 2), flow_vph=round(float(np.average(g.flow, weights=g.len))), cars=round(float(g.cars.sum()), 1), speed_kmh=round(float(g.len.sum()/max((g.len/np.maximum(g.spd/3.6, 1e-9)).sum(), 1e-9)*3.6), 1)))
    streets.sort(key=lambda x: -x['cars']); sm = s['isst'] & win
    series = [dict(clock=f'{6+k//4:02d}:{15*(k % 4):02d}', cars_on_street=round(float(s['cars'][k]), 1), flow_vph=round(float(np.average(FL[k, sm], weights=LEN[sm]))) if sm.any() else 0,
                   speed_kmh=round(float(LEN[sm].sum()/max(TT[k, sm].sum(), 1e-9)*3.6), 1)) for k in range(NS)]
    out = dict(id=s['id'], label=s['label'], name_regex=s['rx'], centre=list(s['centre']), half_m=H, slice=pk, clock=f'{6+pk//4:02d}:{15*(pk % 4):02d}', scale=round(s['ratio'], 3), cars_on_street=round(float(s['cars'][pk]), 1),
               delay_vehicle_hours=round(s['delay'], 1), entering_vph=round(sum(o['vph'] for o in od)), window_total_cars=round(float((FL*TT/3600.0)[pk, wi].sum())), street_km_in_window=round(float(LEN[sm].sum()/1000), 2), series=series, streets=streets[:14], od=od, edge_factors=[[round(float(mx[i]), 6), round(float(my[i]), 6), round(float((np.degrees(np.arctan2((e.v_lon.iat[i]-e.u_lon.iat[i])*kc, e.v_lat.iat[i]-e.u_lat.iat[i]))+360) % 360)), round(float(fac[pk, i]), 2)] for i in np.nonzero((np.abs(MX-s['cx']) <= H+450) & (np.abs(MY-s['cy']) <= H+450))[0]])
    json.dump(out, open(f"{OUT}/scene_{s['id']}.demand.json", 'w'), ensure_ascii=False)
    index.append({k: out[k] for k in ('id', 'label', 'centre', 'clock', 'cars_on_street', 'delay_vehicle_hours', 'entering_vph', 'window_total_cars')}); print(f"written {s['id']}: {len(od)} flows, {out['entering_vph']:,} cars/h entering, {out['window_total_cars']} cars in the window (macro)")
json.dump(index, open(f'{OUT}/scenes_index.json', 'w'), ensure_ascii=False, indent=1)
