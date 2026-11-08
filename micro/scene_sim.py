"""Lane-level SUMO simulation of one 1500 m window, fed with the city-model (macro) flows on EVERY road, exported for the dashboard.
usage: python3 scene_sim.py scene_ID.demand.json tashkent.osm WORKDIR OUT.json [--warmup 600] [--record 300] [--seed 1]
 1. netconvert cuts the window (+300 m margin) out of the OSM roads (OSM signals kept, fixed 90 s cycle: an assumption).
 2. Each macro entry -> exit flow (cars/h at the peak 15 minutes) is placed on the SUMO roads by position, direction and street name; duarouter finds the SUMO routes.
 3. SUMO runs --warmup s so the roads fill up, then --record s are saved: every car's position each second (lane level), plus lane geometry.
 4. The result is compared with the macro model: cars/h per street in the window, macro vs SUMO."""
import sys, os, re, json, math, html, base64, argparse, subprocess, shutil, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('demand'); ap.add_argument('osm'); ap.add_argument('workdir'); ap.add_argument('out')
ap.add_argument('--warmup', type=int, default=600); ap.add_argument('--record', type=int, default=300); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--cycle', type=int, default=90)
ap.add_argument('--drain', default='', help='test: comma-separated road ids whose cars are removed as they arrive (an idealised junction that passes everything), to check whether jams upstream disappear'); ap.add_argument('--relieve', default='', help='test: comma-separated junction ids whose traffic light is switched off (cars pass by priority rules), to check whether a jam upstream disappears'); ap.add_argument('--edge-stats', help='write per-road queue time series (for micro/cascade.py) to this json'); ap.add_argument('--tls', default='actuated', choices=['static', 'actuated']); ap.add_argument('--sigma', type=float, default=0.5); ap.add_argument('--margin', type=float, default=650); ap.add_argument('--tol', type=float, default=30); ap.add_argument('--teleport', type=int, default=120)
a = ap.parse_args(); os.makedirs(a.workdir, exist_ok=True)
try:
    import pyproj; os.environ.setdefault('PROJ_DATA', pyproj.datadir.get_data_dir())
except ImportError: pass
import sumolib, traci
D = json.load(open(a.demand)); clat, clon = D['centre']; H = D['half_m']; kc = math.cos(math.radians(clat)); sid = D['id']
dlat, dlon = (H + a.margin) / 111320, (H + a.margin) / (111320 * kc)
net_f = f'{a.workdir}/{sid}.{a.tls}.net.xml'
if not os.path.isfile(net_f):
    cmd = ['netconvert', '--osm-files', a.osm, '-o', net_f, '--keep-edges.in-geo-boundary', f'{clon-dlon:.6f},{clat-dlat:.6f},{clon+dlon:.6f},{clat+dlat:.6f}', '--geometry.remove', '--roundabouts.guess',
           '--junctions.join', '--tls.join', '--tls.default-type', a.tls, '--tls.cycle.time', str(a.cycle), '--keep-edges.by-vclass', 'passenger', '--remove-edges.isolated', '--output.street-names',
           '--ignore-errors', '--no-warnings']
    tm = '/usr/local/share/sumo/data/typemap/osmNetconvert.typ.xml'
    for base in (os.environ.get('SUMO_HOME'), '/usr/local/share/sumo', '/usr/share/sumo'):
        if base and os.path.isfile(os.path.join(base, 'data/typemap/osmNetconvert.typ.xml')): cmd += ['--type-files', os.path.join(base, 'data/typemap/osmNetconvert.typ.xml')]; break
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or not os.path.isfile(net_f): sys.exit('netconvert failed:\n' + r.stderr[-1500:])
net = sumolib.net.readNet(net_f, withInternal=True); C = np.array(net.convertLonLat2XY(clon, clat))
cand = [e for e in net.getEdges(withInternal=False) if e.allows('passenger') and e.getLength() > 3]
norm = lambda s: re.sub(r'[^a-z0-9]', '', html.unescape(s or '').lower())
bearing = lambda p, q: (math.degrees(math.atan2(q[0]-p[0], q[1]-p[1])) + 360) % 360
angd = lambda x, y: min(abs(x-y) % 360, 360 - abs(x-y) % 360)
def near(P, sh):
    best = (1e9, 0.0)
    for q, r in zip(sh[:-1], sh[1:]):
        q, r = np.array(q), np.array(r); v = r-q; L2 = float(v@v) or 1e-9; t = min(1.0, max(0.0, float((P-q)@v)/L2)); dd = float(np.linalg.norm(q+t*v-P))
        if dd < best[0]: best = (dd, bearing(q, r))
    return best
cache = {}
def match(end):
    key = (end['lon'], end['lat'], end['bearing'])
    if key in cache: return cache[key]
    P = np.array(net.convertLonLat2XY(end['lon'], end['lat'])); names = [norm(x) for x in end.get('road_names', []) if x]; res = None
    for tol, atol in ((a.tol, 50), (a.tol * 2, 70)):
        best = None
        for e, d0 in net.getNeighboringEdges(P[0], P[1], tol):
            if e.getFunction() != '' or not e.allows('passenger') or e.getLength() <= 3: continue
            dd, brg = near(P, e.getShape()); ang = angd(brg, end['bearing'])
            if dd > tol or ang > atol: continue
            en = norm(e.getName()); hit = bool(en) and any(en in n_ or n_ in en for n_ in names if n_)
            cost = dd + 0.5 * ang - (25 if hit else 0)
            if best is None or cost < best[0]: best = (cost, e)
        if best: res = best[1]; break
    cache[key] = res; return res
pairs = {}; lost = 0.0; unm = []
for o in D['od']:
    s_, d_ = match(o['src']), match(o['dst'])
    if not s_ or not d_ or s_.getID() == d_.getID(): lost += o['vph']; unm.append((o['vph'], o['src']['road'] if not s_ else '-', o['src']['highway'] if not s_ else '', o['dst']['road'] if not d_ else '-', o['dst']['highway'] if not d_ else '')); continue
    pairs[(s_.getID(), d_.getID())] = pairs.get((s_.getID(), d_.getID()), 0.0) + o['vph']
tot = sum(o['vph'] for o in D['od']); placed = sum(pairs.values()); print(f"{sid}: {len(D['od'])} macro flows, {tot:,.0f} cars/h; placed on the SUMO roads {placed:,.0f} cars/h ({100*placed/max(tot,1):.0f}%) as {len(pairs)} edge pairs")
for u_ in sorted(unm, reverse=True)[:6]: print('   not placed: %6.0f veh/h  src %s [%s]  dst %s [%s]' % u_)
T_END = a.warmup + a.record + 1; trips = f'{a.workdir}/{sid}.trips.xml'; routes = f'{a.workdir}/{sid}.rou.xml'
with open(trips, 'w') as f:
    f.write(f'<routes>\n<vType id="car" accel="2.6" decel="4.5" sigma="{a.sigma}" length="4.5" speedFactor="normc(1,0.1,0.8,1.2)"/>\n')
    for k, ((fe, te), v) in enumerate(sorted(pairs.items(), key=lambda kv: -kv[1])):
        f.write(f'<flow id="f{k}" type="car" from="{fe}" to="{te}" begin="0" end="{T_END}" vehsPerHour="{v:.2f}" departLane="best" departSpeed="max"/>\n')
    f.write('</routes>\n')
# travel-time weights for the router: driving time on every SUMO road x the congestion factor of the nearest macro road heading the same way (macro routes are shortest paths at congested times)
from scipy.spatial import cKDTree
EF_ = np.array(D['edge_factors']); MXY = np.c_[(EF_[:, 0] - clon) * 111320 * kc, (EF_[:, 1] - clat) * 111320]; tree = cKDTree(MXY); wfile = f'{a.workdir}/{sid}.weights.xml'; nf_ = 0
with open(wfile, 'w') as f:
    f.write('<meandata><interval begin="0" end="999999">\n')
    for e in net.getEdges(withInternal=False):
        sh = e.getShape(); i = len(sh) // 2; p0, p1 = np.array(sh[max(i-1, 0)]), np.array(sh[min(i, len(sh)-1)]); mid = (p0 + p1) / 2 - C; brg = bearing(p0, p1) if np.any(p0 != p1) else 0.0; fac = 1.0
        dd, ii = tree.query(mid, k=6, distance_upper_bound=60)
        for d_, i_ in zip(np.atleast_1d(dd), np.atleast_1d(ii)):
            if np.isfinite(d_) and i_ < len(EF_) and angd(brg, EF_[i_, 2]) < 70: fac = max(1.0, EF_[i_, 3]); nf_ += 1; break
        f.write(f'<edge id="{e.getID()}" traveltime="{e.getLength() / max(e.getSpeed(), 1.0) * fac:.2f}"/>\n')
    f.write('</interval></meandata>\n')
r = subprocess.run(['duarouter', '-n', net_f, '--route-files', trips, '-o', routes, '-w', wfile, '--ignore-errors', '--no-warnings', '--seed', str(a.seed)], capture_output=True, text=True)
if r.returncode or not os.path.isfile(routes): sys.exit('duarouter failed:\n' + r.stderr[-1500:])
# ---- simulate
traci.start(['sumo', '-n', net_f, '-r', routes, '--step-length', '0.5', '--seed', str(a.seed), '--no-warnings', '--no-step-log', '--time-to-teleport', str(a.teleport), '--collision.action', 'warn',
             '--end', str(T_END + 5), '--lateral-resolution', '0'], label=sid)
tc = traci.constants; subs = [tc.VAR_POSITION, tc.VAR_SPEED, tc.VAR_ROAD_ID]
tracks = {}; segs = []; cross = {}; series = []; teleports = 0; pend = []
step = 0.0
for j_ in [x for x in a.relieve.split(',') if x]:
    tl_ids = traci.trafficlight.getIDList(); tid = j_ if j_ in tl_ids else None
    if tid is None and j_ in {n.getID() for n in net.getNodes()}:
        jc = np.array(net.getNode(j_).getCoord()); d_ = [(float(np.linalg.norm(np.array(net.getNode(t_).getCoord()) - jc)), t_) for t_ in tl_ids if net.hasNode(t_)]; tid = min(d_)[1] if d_ and min(d_)[0] < 40 else None
    if tid: traci.trafficlight.setProgram(tid, 'off'); print('traffic light', tid, 'switched off (test)')
    else: print('no traffic light found for junction', j_)
EDG = []; ES = {}; DRAIN = [x for x in a.drain.split(',') if x]
if a.edge_stats:       # queue (halting cars) and car count of every road in the window, every 2 s
    for e in net.getEdges(withInternal=False):
        sh = np.array(e.getShape()); m = sh[len(sh)//2] - C
        if abs(m[0]) <= H + 100 and abs(m[1]) <= H + 100 and e.allows('passenger'):
            fn, tn = e.getFromNode(), e.getToNode(); EDG.append(dict(id=e.getID(), frm=fn.getID(), to=tn.getID(), len=round(e.getLength(), 1), lanes=e.getLaneNumber(), name=html.unescape(e.getName() or ''), tls=tn.getType() == 'traffic_light', x=round(float(m[0]), 1), y=round(float(m[1]), 1),
                                                                                                       fx=round(float(fn.getCoord()[0] - C[0]), 1), fy=round(float(fn.getCoord()[1] - C[1]), 1), tx=round(float(tn.getCoord()[0] - C[0]), 1), ty=round(float(tn.getCoord()[1] - C[1]), 1)))
            traci.edge.subscribe(e.getID(), [tc.LAST_STEP_VEHICLE_HALTING_NUMBER, tc.LAST_STEP_VEHICLE_NUMBER]); ES[e.getID()] = ([], [])
def inwin(x, y, m=0): return abs(x) <= H + m and abs(y) <= H + m
while step < a.warmup + a.record:
    traci.simulationStep(); step += 0.5; teleports += traci.simulation.getStartingTeleportNumber()
    for eid_ in DRAIN:
        for v_ in traci.edge.getLastStepVehicleIDs(eid_): traci.vehicle.remove(v_)
    for v in traci.simulation.getDepartedIDList(): traci.vehicle.subscribe(v, subs)
    if abs(step - round(step)) > 1e-6: continue
    sec = int(round(step))
    if sec % 30 == 0 and sec <= a.warmup + a.record: pend.append(len(traci.simulation.getPendingVehicles()))
    if sec < a.warmup: continue
    res = traci.vehicle.getAllSubscriptionResults(); n_in = 0; vs = []
    for v, d in res.items():
        x, y = d[tc.VAR_POSITION][0] - C[0], d[tc.VAR_POSITION][1] - C[1]; spd = d[tc.VAR_SPEED]; road = d[tc.VAR_ROAD_ID]
        t = tracks.setdefault(v, dict(road=None, seg=None))
        if t['seg'] is not None and math.hypot(x - t['seg']['x'][-1], y - t['seg']['y'][-1]) > 40: t['seg'] = None      # SUMO teleported this car out of a jam: draw it as a new car
        if inwin(x, y, 60):
            if t['seg'] is None: t['seg'] = dict(t0=sec - a.warmup, x=[], y=[], v=[]); segs.append(t['seg'])
            t['seg']['x'].append(x); t['seg']['y'].append(y); t['seg']['v'].append(spd)
        else: t['seg'] = None
        if road != t['road'] and road and not road.startswith(':') and inwin(x, y): cross[road] = cross.get(road, 0) + 1
        if road: t['road'] = road
        if inwin(x, y): n_in += 1; vs.append(spd)
    series.append((sec - a.warmup, n_in, float(np.mean(vs)) * 3.6 if vs else 0.0))
    if ES and sec % 2 == 0:
        for eid, r in traci.edge.getAllSubscriptionResults().items(): ES[eid][0].append(r[tc.LAST_STEP_VEHICLE_HALTING_NUMBER]); ES[eid][1].append(r[tc.LAST_STEP_VEHICLE_NUMBER])
backlog = [len(traci.simulation.getPendingVehicles())]; traci.close()
if a.edge_stats: json.dump(dict(id=sid, record=a.record, dt=2, teleports=teleports, edges=EDG, halt={k: v[0] for k, v in ES.items()}, nveh={k: v[1] for k, v in ES.items()}), open(a.edge_stats, 'w'), separators=(',', ':')); print('edge stats:', len(EDG), 'roads')
# ---- street flows, SUMO vs macro
out_streets = {}
for e in net.getEdges(withInternal=False):
    nm = html.unescape(e.getName() or '')
    if not nm or e.getID() not in cross and False: continue
    x, y = np.array(e.getShape()[len(e.getShape())//2]) - C
    if not inwin(x, y): continue
    k = norm(nm); o = out_streets.setdefault(k, dict(name=nm, len=0.0, cross=0.0)); o['len'] += e.getLength(); o['cross'] += cross.get(e.getID(), 0) * e.getLength()
micro_streets = {k: dict(name=o['name'], flow_vph=round(o['cross'] / max(o['len'], 1) * 3600 / a.record), km=round(o['len'] / 1000, 2)) for k, o in out_streets.items() if o['len'] > 80}
cmp_rows = []
for s in D['streets']:
    m = micro_streets.get(norm(s['street']))
    cmp_rows.append(dict(street=s['street'], macro_vph=s['flow_vph'], sumo_vph=m['flow_vph'] if m else None, macro_cars=s['cars'], macro_kmh=s['speed_kmh']))
# ---- export geometry and cars
q = lambda v: np.clip(np.round(np.asarray(v) * 10), -32000, 32000).astype(np.int16); B = lambda arr: base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode()
lanes_pts, lanes_off, lanes_w, lanes_cls, labels = [], [0], [], [], {}
for e in net.getEdges(withInternal=True):
    internal = e.getFunction() == 'internal'; nm = html.unescape(e.getName() or '')
    for ln in e.getLanes():
        sh = np.array(ln.getShape()) - C
        if not (np.abs(sh).max(1) <= H + 150).any(): continue
        spd = ln.getSpeed(); cls = 3 if internal else 2 if spd >= 15 else 1 if spd >= 9 else 0
        lanes_pts.append(q(sh).reshape(-1)); lanes_off.append(lanes_off[-1] + len(sh)); lanes_w.append(min(60, int(round(ln.getWidth() * 10)))); lanes_cls.append(cls)
    if not internal and nm and e.getLength() > 45:
        sh = np.array(e.getLanes()[0].getShape()) - C; L = np.r_[0, np.cumsum(np.hypot(*np.diff(sh, axis=0).T))]; h = np.searchsorted(L, L[-1] / 2); h = max(1, min(h, len(sh) - 1))
        p = sh[h-1] + (sh[h] - sh[h-1]) * ((L[-1] / 2 - L[h-1]) / max(L[h] - L[h-1], 1e-9)); ang = math.degrees(math.atan2(sh[h][1] - sh[h-1][1], sh[h][0] - sh[h-1][0]))
        if inwin(p[0], p[1], -20):
            prev = labels.get(nm)
            if prev is None or e.getLength() > prev[3]: labels[nm] = (round(float(p[0]), 1), round(float(p[1]), 1), round(ang, 1), e.getLength())
tls = [[round(float(n.getCoord()[0] - C[0]), 1), round(float(n.getCoord()[1] - C[1]), 1)] for n in net.getNodes() if n.getType() == 'traffic_light' and inwin(n.getCoord()[0] - C[0], n.getCoord()[1] - C[1], 60)]
tr = [(i, t) for i, t in enumerate(segs) if len(t['x']) >= 2]
off = np.cumsum([0] + [len(t['x']) for _, t in tr]).astype(np.uint32)
def dlt(key):      # first sample absolute, then differences (decimetres): stationary cars give zeros, which compress very well
    out = []
    for _, t in tr:
        z = q(t[key]).astype(np.int32); d = z.copy(); d[1:] = np.diff(z); out.append(np.clip(d, -32000, 32000).astype(np.int16))
    return np.concatenate(out)
X = dlt('x'); Y = dlt('y')
V = np.clip(np.round(np.concatenate([t['v'] for _, t in tr]) * 3), 0, 255).astype(np.uint8); T0 = np.array([t['t0'] for _, t in tr], np.uint16)
mean_in = float(np.mean([s[1] for s in series])); mean_v = float(np.mean([s[2] for s in series if s[1]])) if series else 0
res = dict(id=sid, label=D['label'], centre=D['centre'], half=H, record=a.record, clock=D['clock'],
           lanes=dict(n=len(lanes_w), off=B(np.array(lanes_off, np.uint32)), pts=B(np.concatenate(lanes_pts)), w=B(np.array(lanes_w, np.uint8)), cls=B(np.array(lanes_cls, np.uint8))),
           labels=[dict(t=k, x=v[0], y=v[1], a=v[2]) for k, v in labels.items()], tls=tls, cars=dict(n=len(tr), t0=B(T0), off=B(off), x=B(X), y=B(Y), v=B(V)),
           macro=dict(cars_on_street=D['cars_on_street'], entering_vph=D['entering_vph'], window_total_cars=D['window_total_cars'], street_km=D['street_km_in_window'], scale=D['scale'], series=D['series'], streets=D['streets']),
           micro=dict(placed_vph=round(placed), placed_share=round(placed / max(tot, 1), 3), mean_cars_in_window=round(mean_in), mean_kmh=round(mean_v, 1), teleports=teleports, waiting_to_enter=pend[-1] if pend else 0, compare=cmp_rows,
                      series=[[s[0], s[1], round(s[2], 1)] for s in series[::5]]))
json.dump(res, open(a.out, 'w'), ensure_ascii=False, separators=(',', ':'))
print(f"{sid}: {len(tr)} cars recorded, {mean_in:.0f} cars in the window on average (macro: {D['window_total_cars']:.0f}), mean speed {mean_v:.1f} km/h, teleports {teleports}, waiting to enter {res['micro']['waiting_to_enter']}; {os.path.getsize(a.out)/1e3:.0f} kB")
for r_ in cmp_rows[:8]: print(f"   {r_['street'][:34]:34s} macro {r_['macro_vph']:>6} veh/h   SUMO {r_['sumo_vph']}")
