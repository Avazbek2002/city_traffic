"""
Live traffic animation of Makhtumkuli street (Tashkent) on REAL OpenStreetMap geometry.

What it does
  1. Gets the street + 300 m around it from OpenStreetMap
        - from uzbekistan-latest.osm.pbf (Geofabrik) if that file is in the folder, else
        - from the Overpass API (needs internet; tries several mirrors)
     If several streets share the name, it keeps the one with the most traffic lights.
  2. Builds a SUMO road network from it (real lanes, real turns, real traffic lights).
  3. Adds traffic and opens a LIVE animation: cars coloured by speed, every signal
     coloured by its real state, junctions numbered along the street.

By default the window shows two copies of the street side by side, same traffic:
     left  - lights NOT coordinated (random timing offsets)
     right - lights as a GREEN WAVE for 55 km/h
with live counters (stopped cars, average through-trip time).

Setup (once):
    SUMO installed (macOS: official .pkg from https://sumo.dlr.de/docs/Installing/)
    pip install traci sumolib pyproj numpy matplotlib osmium      (osmium only for the .pbf route; sumolib + pyproj only for --macro-demand)
    ffmpeg optional (only for --save video.mp4)

Run (in the folder that has uzbekistan-latest.osm.pbf, or with internet):
    python makhtumkuli_live.py                    # live, one version (green wave)
    python makhtumkuli_live.py --mode uncoordinated   # the other version
    python makhtumkuli_live.py --mode compare     # both, stacked
    python makhtumkuli_live.py --save mq.mp4 --no-show --end 900   # video file
    python makhtumkuli_live.py --osm my.osm       # use your own extract

Traffic from the macro (city) simulation instead of assumed volumes:
    python street_demand.py ... macro_demand.json        (in the city repo: pipeline/street_demand.py)
    python makhtumkuli_live.py --macro-demand macro_demand.json --stats-only --end 900
  Every flow in the file (cars/h at the peak 15 minutes, from one road onto the street and off to another road) becomes a SUMO flow between the
  matching SUMO edges, so the junction approaches carry the volumes the city model computed. --demand-scale multiplies all of them (e.g. 0.8 for a quieter day).

IMPORTANT - what is real and what is assumed
  Real: street layout, lanes, turn restrictions and traffic-light positions (as mapped
  in OpenStreetMap, which can be wrong - check the printed junction table).
  Assumed: signal phase lengths (90 s cycle), traffic volumes (--through-vph,
  --side-vph) and the 55 km/h design speed. Replace them with counts you measured.
"""
import argparse, glob, itertools, math, os, random, re, shutil, subprocess, sys
import urllib.parse, urllib.request, json, html
import xml.etree.ElementTree as ET
import numpy as np
try:                                                       # let SUMO's PROJ library find its data (avoids "Cannot find proj.db")
    import pyproj
    os.environ.setdefault("PROJ_DATA", pyproj.datadir.get_data_dir())
    os.environ.setdefault("PROJ_LIB", pyproj.datadir.get_data_dir())
except Exception:
    pass

NAME_PATTERN = r"mahtum|makhtum|махтум"
TASHKENT = (41.15, 69.05, 41.45, 69.45)            # S, W, N, E search window
MARGIN_M = 300
ROADS = ("motorway|trunk|primary|secondary|tertiary|unclassified|residential|living_street|"
         "service|motorway_link|trunk_link|primary_link|secondary_link|tertiary_link").split("|")
OVERPASS = ["https://overpass-api.de/api/interpreter",
            "https://overpass.kumi.systems/api/interpreter",
            "https://overpass.private.coffee/api/interpreter"]


# ============================ SUMO discovery ================================
def sumo_bin(tool):
    cands = ([os.path.join(os.environ["SUMO_HOME"], "bin")] if os.environ.get("SUMO_HOME") else [])
    cands += sorted(glob.glob("/Library/Frameworks/EclipseSUMO.framework/Versions/*/EclipseSUMO/bin"))
    cands += [os.path.join(sys.prefix, "bin"), "/opt/homebrew/bin", "/usr/local/bin"]
    seen, tried = set(), []
    for exe in [os.path.join(d, tool) for d in cands] + [shutil.which(tool) or ""]:
        if not os.path.isfile(exe) or os.path.realpath(exe) in seen:
            continue
        seen.add(os.path.realpath(exe))
        r = subprocess.run([exe, "--version"], capture_output=True, text=True)
        if r.returncode == 0:
            return exe
        tried.append(f"  {exe}: exit {r.returncode}")
    sys.exit(f"Could not find a working '{tool}'.\n" + "\n".join(tried) +
             "\nInstall SUMO (macOS: official .pkg from https://sumo.dlr.de/docs/Installing/) "
             "and open a new terminal.")


def import_traci():
    try:
        import traci
        return traci
    except ImportError:
        pass
    homes = [os.environ.get("SUMO_HOME", "")] + glob.glob(
        "/Library/Frameworks/EclipseSUMO.framework/Versions/*/EclipseSUMO")
    for h in homes:
        if h and os.path.isdir(os.path.join(h, "tools")):
            sys.path.append(os.path.join(h, "tools"))
            try:
                import traci
                return traci
            except ImportError:
                pass
    sys.exit("Python package 'traci' not found. Run:  pip install traci")


# ============================ 1. OpenStreetMap data =========================
def choose_street(ways, signals, name_re):
    """ways: [{'nodes':[ids], 'pts':[(lat,lon)]}]. Keep the connected group of ways
    with the most traffic signals (ties: longest). Returns (ways, bbox S,W,N,E)."""
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for w in ways:
        for n in w["nodes"][1:]:
            parent[find(n)] = find(w["nodes"][0])
    groups = {}
    for w in ways:
        groups.setdefault(find(w["nodes"][0]), []).append(w)

    def score(g):
        sig = len({n for w in g for n in w["nodes"] if n in signals})
        ln = sum(111_320 * math.hypot(b[0] - a[0], (b[1] - a[1]) * 0.75)
                 for w in g for a, b in zip(w["pts"], w["pts"][1:]))
        return sig, ln
    best = max(groups.values(), key=score)
    print(f"Found {len(groups)} separate street(s) with this name; using the one with "
          f"{score(best)[0]} traffic signals, {score(best)[1] / 1000:.1f} km.")
    lats = [p[0] for w in best for p in w["pts"]]
    lons = [p[1] for w in best for p in w["pts"]]
    dlat = MARGIN_M / 111_320
    dlon = MARGIN_M / (111_320 * math.cos(math.radians(sum(lats) / len(lats))))
    return best, (min(lats) - dlat, min(lons) - dlon, max(lats) + dlat, max(lons) + dlon)


def street_bbox_from_pbf(pbf, name_re):
    import osmium
    print(f"Reading {pbf} (first pass: find the street) ...")
    ways, signals = [], {}
    for o in osmium.FileProcessor(pbf, osmium.osm.NODE | osmium.osm.WAY).with_locations():
        if o.is_node():
            if o.tags.get("highway") == "traffic_signals":
                signals[o.id] = (o.location.lat, o.location.lon)
        elif o.tags.get("highway") in ROADS and any(
                name_re.search(v) for k, v in o.tags if k.startswith("name") or k in ("alt_name", "int_name")):
            ways.append({"nodes": [n.ref for n in o.nodes if n.location.valid()],
                         "pts": [(n.location.lat, n.location.lon) for n in o.nodes if n.location.valid()]})
    ways = [w for w in ways if len(w["nodes"]) > 1]
    if not ways:
        sys.exit("No street with 'Mahtumquli' in its name was found in the file.")
    return choose_street(ways, signals, name_re)[1]


def get_osm_from_pbf(pbf, out, name_re, bbox=None):
    try:
        import osmium
    except ImportError:
        sys.exit("pip install osmium   (needed to read the .pbf file)")
    S, W, N, E = bbox or street_bbox_from_pbf(pbf, name_re)
    print(f"Cutting S={S:.5f} W={W:.5f} N={N:.5f} E={E:.5f} out of {pbf} ...")
    inside = lambda l: l.valid() and S <= l.lat <= N and W <= l.lon <= E
    keep_w, keep_n = set(), set()
    for o in osmium.FileProcessor(pbf, osmium.osm.NODE | osmium.osm.WAY).with_locations():
        if o.is_way() and o.tags.get("highway") in ROADS and any(inside(n.location) for n in o.nodes):
            keep_w.add(o.id)
            keep_n.update(n.ref for n in o.nodes)
    with osmium.SimpleWriter(out, overwrite=True) as wr:
        for o in osmium.FileProcessor(pbf, osmium.osm.NODE | osmium.osm.WAY):
            if (o.is_node() and o.id in keep_n) or (o.is_way() and o.id in keep_w):
                wr.add(o)
    print(f"Saved {out} ({len(keep_w):,} road segments)")


def overpass(query, as_json):
    last = None
    for url in OVERPASS:
        try:
            req = urllib.request.Request(url, urllib.parse.urlencode({"data": query}).encode(),
                                         headers={"User-Agent": "makhtumkuli-sim/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                raw = r.read()
            return json.loads(raw) if as_json else raw
        except Exception as e:
            last = e
            print(f"  {url} failed: {e}")
    sys.exit(f"All Overpass mirrors failed ({last}). Download uzbekistan-latest.osm.pbf from "
             "https://download.geofabrik.de/asia/uzbekistan.html into this folder and rerun.")


def street_bbox_from_overpass(name_re):
    S, W, N, E = TASHKENT
    print("Asking Overpass for the street ...")
    d = overpass(f'[out:json][timeout:120];(way["highway"]["name"~"Mahtum|Makhtum|Махтум",i]'
                 f'({S},{W},{N},{E});way["highway"]["name:en"~"Makhtum|Mahtum",i]({S},{W},{N},{E}););'
                 f'out geom;', True)
    ways = [{"nodes": e["nodes"], "pts": [(g["lat"], g["lon"]) for g in e["geometry"]]}
            for e in d["elements"] if e["type"] == "way" and len(e.get("nodes", [])) > 1]
    if not ways:
        sys.exit("Overpass found no street with 'Mahtumquli' in its name.")
    sig = overpass(f'[out:json][timeout:120];node["highway"="traffic_signals"]({S},{W},{N},{E});out;', True)
    signals = {e["id"]: (e["lat"], e["lon"]) for e in sig["elements"]}
    return choose_street(ways, signals, name_re)[1]


def get_osm_from_overpass(out, name_re, bbox=None):
    s, w, n, e = bbox or street_bbox_from_overpass(name_re)
    print("Downloading all roads in the box ...")
    raw = overpass(f'[out:xml][timeout:180];(way["highway"~"^({"|".join(ROADS)})$"]({s},{w},{n},{e});>;);out;',
                   False)
    open(out, "wb").write(raw)
    print(f"Saved {out}")


# ============================ 2. SUMO network ===============================
def build_net(osm, net, cycle):
    exe = sumo_bin("netconvert")
    cmd = [exe, "--osm-files", osm, "-o", net, "--geometry.remove", "--roundabouts.guess",
           "--junctions.join", "--tls.join", "--tls.default-type", "static",
           "--tls.cycle.time", str(cycle), "--keep-edges.by-vclass", "passenger",
           "--remove-edges.isolated", "--output.street-names", "--ignore-errors", "--no-warnings"]
    for base in (os.environ.get("SUMO_HOME"), os.path.dirname(os.path.dirname(os.path.realpath(exe)))):
        tm = base and os.path.join(base, "data", "typemap", "osmNetconvert.typ.xml")
        if tm and os.path.isfile(tm):
            cmd += ["--type-files", tm]
            break
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or not os.path.isfile(net):
        sys.exit("netconvert failed:\n" + r.stderr[-1500:])


def shape(s):
    return np.array([[float(v) for v in p.split(",")[:2]] for p in s.split()])


def load_net(path):
    root = ET.parse(path).getroot()
    edges = {}
    for e in root.iter("edge"):
        if e.get("function") == "internal":
            continue
        edges[e.get("id")] = dict(name=e.get("name", ""), frm=e.get("from"), to=e.get("to"),
                                  lanes={l.get("id"): shape(l.get("shape")) for l in e.findall("lane")},
                                  w={l.get("id"): float(l.get("width", 3.2)) for l in e.findall("lane")})
    polys = [shape(j.get("shape")) for j in root.iter("junction")
             if j.get("shape") and j.get("type") != "internal"]
    return edges, polys


# ============================ 3. traffic ====================================
def corridor_axis(edges, name_re):
    cor = {i: e for i, e in edges.items() if name_re.search(e["name"])}
    if not cor:
        sys.exit("The built network has no edges named like Mahtumquli (names lost in conversion?).")
    pts = np.vstack([s for e in cor.values() for s in e["lanes"].values()])
    c = pts.mean(0)
    axis = np.linalg.svd(pts - c)[2][0]
    if (axis[1] if abs(axis[1]) >= abs(axis[0]) else axis[0]) < 0:
        axis = -axis                                    # point north (or east)
    return cor, c, axis


def ang_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def bearing_of(p, q):                                      # compass bearing of the step p -> q (net coordinates, y = north)
    return (math.degrees(math.atan2(q[0] - p[0], q[1] - p[1])) + 360) % 360


def make_macro_demand(edges, cor, c, axis, args, files):
    """Flows from the city (macro) model: each entry (road, junction) -> exit (road, junction) with cars/h.
    Entries / exits are matched to SUMO edges by position (junction within --map-tol-m), road name and direction of travel."""
    try:
        import sumolib
    except ImportError:
        sys.exit("pip install sumolib pyproj   (needed to place the macro-model flows on the SUMO network)")
    d = json.load(open(args.macro_demand))
    if abs(d["centre"][0] - float(args.center.split(",")[0])) > 1e-4 or abs(d["centre"][1] - float(args.center.split(",")[1])) > 1e-4 or abs(d["radius_m"] - args.radius) > 1:
        print(f"WARNING: the demand file was made for centre {d['centre']} radius {d['radius_m']:.0f} m, this run uses {args.center} / {args.radius:.0f} m.")
    net = sumolib.net.readNet(files["net"], withInternal=False)
    cand = [e for e in net.getEdges() if e.allows("passenger") and e.getLength() > 3]
    name_hit = lambda sumo_name, variants: any(v.lower() in sumo_name.lower() or sumo_name.lower() in v.lower() for v in variants if v and v != "(unnamed)") if sumo_name else False
    is_street = lambda e: bool(re.search(args.name_pattern, html.unescape(e.getName()), re.I))

    def seg_dist(P, sh):                                   # distance from P to a polyline, and the bearing of the nearest piece
        best = (1e9, 0.0)
        for q, r in zip(sh[:-1], sh[1:]):
            q, r = np.array(q), np.array(r); v = r - q; L2 = float(v @ v) or 1e-9
            t = min(1.0, max(0.0, float((P - q) @ v) / L2)); dd = float(np.linalg.norm(q + t * v - P))
            if dd < best[0]:
                best = (dd, bearing_of(q, r))
        return best

    def match(end, incoming):
        P = np.array(net.convertLonLat2XY(end["lon"], end["lat"]))
        want = end.get("bearing")
        for tol, loose in ((args.map_tol_m, False), (args.map_tol_m * 1.7, True)):
            best = None
            for e in cand:
                sh = e.getShape(); on_street = is_street(e)
                if end["kind"] == "street_end":
                    # the macro node is on the street but usually in the middle of a SUMO edge: take the edge of the right direction that passes through it
                    if not on_street:
                        continue
                    dist, brg = seg_dist(P, sh)
                    if dist > min(tol, 40):
                        continue
                    ang = ang_diff(brg, want) if want is not None else 0.0
                    if ang > 60:
                        continue
                    cost = dist + 0.5 * ang
                else:
                    if on_street:
                        continue
                    q, r = (np.array(sh[-1]), np.array(sh[-2])) if incoming else (np.array(sh[0]), np.array(sh[1]))
                    dist = float(np.linalg.norm(q - P))
                    if dist > tol:
                        continue
                    brg = bearing_of(r, q) if incoming else bearing_of(q, r)
                    ang = ang_diff(brg, want) if want is not None else 0.0
                    nm = name_hit(html.unescape(e.getName()), end.get("road_names", []))
                    cost = dist + 0.5 * ang - (25 if nm else 0)
                if best is None or cost < best[0]:
                    best = (cost, e, dist, ang, loose)
            if best:
                return best
        return None

    cache, rows = {}, []
    for k, o in enumerate(d["od"]):
        ends = []
        for key, inc in (("src", True), ("dst", False)):
            ck = (key, o[key]["node"], o[key]["road"])
            if ck not in cache:
                cache[ck] = match(o[key], inc)
            ends.append(cache[ck])
        rows.append((o, ends[0], ends[1]))
    ok = [(o, a, b) for o, a, b in rows if a and b and a[1].getID() != b[1].getID()]
    tot, got = sum(o["vph"] for o in d["od"]), sum(o["vph"] for o, _, _ in ok)
    print(f"Macro demand: {len(d['od'])} flows, {tot:,.0f} cars/h at {d['clock']} ({d['scenario']} scenario, {d['cars_on_stretch']:.0f} cars on the stretch in the city model);"
          f" {len(ok)} flows / {got:,.0f} cars/h could be placed on the SUMO network ({100 * got / max(tot, 1):.0f}%).")
    print("  entries (cars/h  SUMO edge  road  distance m / angle deg):")
    for (key, node, road), m in sorted(cache.items(), key=lambda kv: -sum(o["vph"] for o in d["od"] if (key, o[key]["node"], o[key]["road"]) == kv[0])):
        v = sum(o["vph"] for o in d["od"] if o[key]["node"] == node and o[key]["road"] == road)
        if v < 20:
            continue
        print(f"   {key} {v:6.0f}  " + (f"{m[1].getID():24s} {html.unescape(m[1].getName())[:24]:24s} {m[2]:4.0f} m / {m[3]:3.0f}°{'  (loose match)' if m[4] else ''}" if m else "NOT FOUND within tolerance"))
    first = lambda e: np.array(e.getShape()[0])
    out = ['<routes>', '<vType id="car" accel="2.6" decel="4.5" sigma="0.5" length="4.5" speedFactor="normc(1,0.1,0.8,1.2)"/>']
    total = args.warmup + args.end + 600
    for k, (o, a, b) in enumerate(ok):
        ea, eb = a[1], b[1]
        if o["src"]["kind"] == "street_end" and o["dst"]["kind"] == "street_end":
            tag = "thru_fwd" if float((first(eb) - first(ea)) @ axis) > 0 else "thru_rev"
        else:
            tag = "side"
        out.append(f'<flow id="{tag}_{k}" type="car" from="{ea.getID()}" to="{eb.getID()}" begin="0" end="{total}" probability="{o["vph"] * args.demand_scale / 3600:.6f}" '
                   f'departLane="best" departSpeed="max"/>')
    out.append("</routes>")
    open(files["trips"], "w").write("\n".join(out))
    r = subprocess.run([sumo_bin("duarouter"), "-n", files["net"], "--route-files", files["trips"], "-o", files["routes"], "--ignore-errors", "--no-warnings",
                        "--seed", str(args.seed)], capture_output=True, text=True)
    if r.returncode or not os.path.isfile(files["routes"]):
        sys.exit("duarouter failed:\n" + r.stderr[-1500:])
    print(f"Traffic: {got * args.demand_scale:,.0f} cars/h in total entering the street from {len({(a[1].getID()) for _, a, _ in ok})} entry edges (macro model, scale {args.demand_scale}).")


def make_demand(edges, cor, c, axis, args, files):
    if args.macro_demand:
        return make_macro_demand(edges, cor, c, axis, args, files)
    proj = lambda p: float((np.asarray(p) - c) @ axis)
    first = lambda e: next(iter(e["lanes"].values()))
    fwd = [i for i, e in cor.items() if proj(first(e)[-1]) > proj(first(e)[0])]
    rev = [i for i in cor if i not in fwd]
    out = ['<routes>', '<vType id="car" accel="2.6" decel="4.5" sigma="0.5" length="4.5" '
           'speedFactor="normc(1,0.1,0.8,1.2)"/>']
    total = args.warmup + args.end + 600
    for tag, group, key in (("thru_fwd", fwd, lambda i: proj(first(cor[i])[0])),
                            ("thru_rev", rev, lambda i: -proj(first(cor[i])[-1]))):
        if group:
            a, b = min(group, key=key), max(group, key=key)
            out.append(f'<flow id="{tag}" type="car" from="{a}" to="{b}" begin="0" end="{total}" '
                       f'vehsPerHour="{args.through_vph}" departLane="best" departSpeed="max"/>')
    nin, nout = {}, {}
    for i, e in edges.items():
        nout.setdefault(e["frm"], []).append(i)
        nin.setdefault(e["to"], []).append(i)
    ln = lambda i: np.linalg.norm(np.diff(first(edges[i]), axis=0), axis=1).sum()
    # fringe = nodes at the edge of the map: touch at most one neighbouring node
    nbr = lambda n: {edges[i]["to"] for i in nout.get(n, [])} | {edges[i]["frm"] for i in nin.get(n, [])}
    fringe = {n for n in set(nin) | set(nout) if len(nbr(n)) <= 1}
    src = [i for i, e in edges.items() if e["frm"] in fringe and ln(i) > 20]
    dst = [i for i, e in edges.items() if e["to"] in fringe and ln(i) > 20]
    rng = random.Random(args.seed)
    n = int(args.side_vph * total / 3600)
    trips = []
    for _ in range(n * 3):
        if len(trips) >= n or not src or not dst:
            break
        a, b = rng.choice(src), rng.choice(dst)
        if a != b and np.linalg.norm(first(edges[a])[0] - first(edges[b])[-1]) > 500:
            trips.append((rng.uniform(0, total), a, b))
    for k, (t, a, b) in enumerate(sorted(trips)):
        out.append(f'<trip id="side_{k}" type="car" depart="{t:.1f}" from="{a}" to="{b}" '
                   f'departLane="best" departSpeed="max"/>')
    out.append("</routes>")
    open(files["trips"], "w").write("\n".join(out))
    r = subprocess.run([sumo_bin("duarouter"), "-n", files["net"], "--route-files", files["trips"],
                        "-o", files["routes"], "--ignore-errors", "--no-warnings",
                        "--seed", str(args.seed)], capture_output=True, text=True)
    if r.returncode or not os.path.isfile(files["routes"]):
        sys.exit("duarouter failed:\n" + r.stderr[-1500:])
    print(f"Traffic: {args.through_vph}/h each way through the street ({len(fwd)}+{len(rev)} "
          f"street edges), {len(trips)} side trips over {total / 3600:.1f} h (assumed volumes).")


# ============================ 4. simulation =================================
class Sim:
    def __init__(self, traci, label, plan, files, args, tls_info=None):
        self.tc = traci.constants
        traci.start([sumo_bin("sumo"), "-n", files["net"], "-r", files["routes"],
                     "--step-length", str(args.step), "--no-step-log", "--no-warnings",
                     "--time-to-teleport", "120", "--seed", str(args.seed)], label=label)
        self.conn, self.plan, self.args = traci.getConnection(label), plan, args
        self.depart, self.tt = {}, {"thru_fwd": [], "thru_rev": [], "other": []}

    def step(self):
        c = self.conn
        c.simulationStep()
        now = c.simulation.getTime()
        for v in c.simulation.getDepartedIDList():
            c.vehicle.subscribe(v, [self.tc.VAR_POSITION, self.tc.VAR_SPEED, self.tc.VAR_ANGLE])
            self.depart[v] = now
        for v in c.simulation.getArrivedIDList():
            d = self.depart.pop(v, None)
            if d is not None and d >= self.args.warmup:
                key = v.rsplit(".", 1)[0] if v.startswith("thru") else "other"
                self.tt["thru_fwd" if key.startswith("thru_fwd") else "thru_rev" if key.startswith("thru_rev") else "other"].append(now - d)

    def snapshot(self):
        res = self.conn.vehicle.getAllSubscriptionResults()
        xy = np.array([r[self.tc.VAR_POSITION] for r in res.values()]).reshape(-1, 2)
        sp = np.array([r[self.tc.VAR_SPEED] for r in res.values()])
        an = np.array([r[self.tc.VAR_ANGLE] for r in res.values()])
        return xy, sp, an

    def close(self):
        self.conn.close()


def find_junctions(conn, edges, name_re, merge_m):
    """Signals on the street, merged into physical junctions (signals closer than merge_m
    metres to the cluster's first signal are one junction: dual carriageways, crossings)."""
    raw = []
    for tid in conn.trafficlight.getIDList():
        links = conn.trafficlight.getControlledLinks(tid)
        segs, cor_idx, cross = [], [], set()
        for i, grp in enumerate(links):
            e = edges.get(grp[0][0].rsplit("_", 1)[0]) if grp else None
            if e is None:
                segs.append(None)
                continue
            lane = grp[0][0]
            sh = e["lanes"][lane]
            d = sh[-1] - sh[-2]
            d = d / (np.linalg.norm(d) or 1)
            nrm = np.array([-d[1], d[0]]) * e["w"][lane] / 2 * 0.95
            segs.append(np.array([sh[-1] + nrm, sh[-1] - nrm]))
            if name_re.search(e["name"]):
                cor_idx.append(i)
            elif e["name"]:
                cross.add(html.unescape(e["name"]))
        if cor_idx:
            pos = np.mean([segs[i].mean(0) for i in cor_idx], axis=0)
            raw.append(dict(id=tid, pos=pos, segs=segs, cor=cor_idx, cross=cross))
    return raw


def merge_junctions(raw, c, axis, merge_m):
    raw.sort(key=lambda j: float((j["pos"] - c) @ axis))
    out = []
    for j in raw:
        if out and np.linalg.norm(j["pos"] - out[-1]["first"]) <= merge_m:
            g = out[-1]
            g["members"].append(j)
            g["cross"] |= j["cross"]
            g["pos"] = np.mean([m["pos"] for m in g["members"]], axis=0)
        else:
            out.append(dict(first=j["pos"], pos=j["pos"], members=[j], cross=set(j["cross"])))
    return out


def apply_plan(conn, junctions, plan, axis, c, args):
    if plan == "default":
        return
    rng = random.Random(args.seed + 7)
    proj = [float((j["pos"] - c) @ axis) for j in junctions]
    v = args.wave_kmh / 3.6
    for j, p in zip(junctions, proj):
        if plan == "uncoordinated":
            s = None
        else:                                             # green wave: same offset for the whole junction
            s = ((p - min(proj)) if args.wave_direction == "forward" else (max(proj) - p)) / v
        for m in j["members"]:
            ph = conn.trafficlight.getAllProgramLogics(m["id"])[0].phases
            dur = [x.duration for x in ph]
            cyc = sum(dur)
            score = [sum(1 for i in m["cor"] if i < len(x.state) and x.state[i] in "Gg") for x in ph]
            k0 = int(np.argmax(score))                    # phase that serves the street
            sm = rng.uniform(0, cyc) if s is None else s
            elapsed = (-sm) % cyc                         # time since street-green started, at t=0
            k, acc = k0, 0.0
            while acc + dur[k] <= elapsed:
                acc += dur[k]
                k = (k + 1) % len(ph)
            conn.trafficlight.setPhase(m["id"], k)
            conn.trafficlight.setPhaseDuration(m["id"], acc + dur[k] - elapsed)


# ============================ 4b. scheme (road-map) background ==============
TILES = {   # name: (url, attribution)
    "voyager": ("https://basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}@2x.png",
                "© OpenStreetMap contributors © CARTO"),
    "light": ("https://basemaps.cartocdn.com/light_all/{z}/{x}/{y}@2x.png",
              "© OpenStreetMap contributors © CARTO"),
    "dark": ("https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}@2x.png",
             "© OpenStreetMap contributors © CARTO"),
    "osm": ("https://tile.openstreetmap.org/{z}/{x}/{y}.png", "© OpenStreetMap contributors"),
}   # CARTO tiles (voyager / light / dark) need an API key since 2025-26; osm needs none (its usage policy asks for a few cached tiles only, which is what this does)


def fetch_tile(url, path):
    """Download one map tile (cached on disk). Returns False if it can't be fetched."""
    if os.path.isfile(path):
        return True
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "makhtumkuli-sim/1.0 (personal project)"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = r.read()
        with open(path, "wb") as fh:
            fh.write(data)
        return True
    except Exception as e:
        print(f"  tile failed ({e})")
        return False


def tile_xy(lon, lat, z):                                  # lon/lat -> fractional tile index
    n = 2 ** z
    return ((lon + 180) / 360 * n,
            (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n)


def tile_lonlat(tx, ty, z):                                # tile corner -> lon/lat
    n = 2 ** z
    return tx / n * 360 - 180, math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * ty / n))))


def build_basemap(kind, bbox, conn, T, full, cor_y, cache):
    """Road-map tiles, warped into the street-aligned frame the animation is drawn in.
    Returns (RGBA array, extent) or None."""
    from PIL import Image
    url = TILES[kind][0]
    S, W, N, E = bbox
    z = 17 if (full[1] - full[0]) < 6000 else 16
    tx0, ty0 = (int(v) for v in tile_xy(W, N, z))
    tx1, ty1 = (int(v) for v in tile_xy(E, S, z))
    X0, X1, Y0, Y1 = full[0] - 100, full[1] + 100, cor_y[0] - 300, cor_y[1] + 300
    res = max(0.5, (X1 - X0) / 12000)                      # metres per output pixel
    wpx, hpx = int((X1 - X0) / res), int((Y1 - Y0) / res)
    out = Image.new("RGBA", (wpx, hpx), (0, 0, 0, 0))
    os.makedirs(cache, exist_ok=True)
    corner = {}

    def xy(tx, ty):                                        # tile corner -> street-aligned metres
        if (tx, ty) not in corner:
            lon, lat = tile_lonlat(tx, ty, z)
            corner[tx, ty] = T(np.array(conn.simulation.convertGeo(lon, lat, True)))
        return corner[tx, ty]
    todo = []
    for ty in range(ty0, ty1 + 1):
        for tx in range(tx0, tx1 + 1):
            cs = np.array([xy(tx, ty), xy(tx + 1, ty), xy(tx, ty + 1), xy(tx + 1, ty + 1)])
            if cs[:, 0].max() >= X0 and cs[:, 0].min() <= X1 and cs[:, 1].max() >= Y0 and cs[:, 1].min() <= Y1:
                todo.append((tx, ty, cs))
    print(f"Map background '{kind}': {len(todo)} tiles (first run downloads them; cached after).")
    ok = 0
    for tx, ty, cs in todo:
        path = os.path.join(cache, f"{kind}_{z}_{tx}_{ty}.png")
        if not fetch_tile(url.format(z=z, x=tx, y=ty), path):
            continue
        tile = Image.open(path).convert("RGBA")
        sz = tile.width
        p00, p10, p01 = cs[0], cs[1], cs[2]
        M = np.array([[(p10[0] - p00[0]) / sz, (p01[0] - p00[0]) / sz, p00[0]],
                      [(p10[1] - p00[1]) / sz, (p01[1] - p00[1]) / sz, p00[1]], [0, 0, 1]])
        Q = np.linalg.inv(M) @ np.array([[res, 0, X0], [0, -res, Y1], [0, 0, 1]])   # pixel -> tile pixel
        i0 = max(0, int((cs[:, 0].min() - X0) / res) - 1)
        i1 = min(wpx, int((cs[:, 0].max() - X0) / res) + 2)
        j0 = max(0, int((Y1 - cs[:, 1].max()) / res) - 1)
        j1 = min(hpx, int((Y1 - cs[:, 1].min()) / res) + 2)
        if i1 <= i0 or j1 <= j0:
            continue
        a, b, c0 = Q[0]
        d, e, f0 = Q[1]
        sub = tile.transform((i1 - i0, j1 - j0), Image.AFFINE,
                             (a, b, a * i0 + b * j0 + c0, d, e, d * i0 + e * j0 + f0), Image.BICUBIC)
        out.alpha_composite(sub, (i0, j0))
        ok += 1
    if not ok:
        print("Could not download any map tiles - continuing without the background.")
        return None
    return np.asarray(out), (X0, X0 + wpx * res, Y1 - hpx * res, Y1)


# ============================ 5. drawing ====================================
BG, ROAD, STREET, TXT = "#0e1217", "#2c3541", "#56708f", "#e8edf3"
MUTED, BADGE_FC, BADGE_TX = "#8b97a5", "#dfe6ee", "#0e1217"
from matplotlib.colors import LinearSegmentedColormap
SPEED_CMAP = LinearSegmentedColormap.from_list("spd", ["#ff3b30", "#ff9f0a", "#ffe14a", "#5ee6a8"])


def set_theme(light):
    """Light theme for road-map (scheme) backgrounds, dark theme otherwise."""
    global BG, ROAD, STREET, TXT, MUTED, BADGE_FC, BADGE_TX, SPEED_CMAP
    if light:
        BG, ROAD, STREET, TXT = "#f4f3ef", "#9aa4b0", "#5c7596", "#1b1f24"
        MUTED, BADGE_FC, BADGE_TX = "#5f6b77", "#1b1f24", "#ffffff"
        SPEED_CMAP = LinearSegmentedColormap.from_list("spd", ["#d7263d", "#f28c28", "#e6b800", "#1f9d6b"])
SIG = {"G": "#2ecc71", "g": "#2ecc71", "y": "#f1c40f", "Y": "#f1c40f", "r": "#e74c3c", "R": "#e74c3c"}


class Live:
    """State of one simulation, read once per frame and shared by all its map panels."""
    def __init__(self, sim, junctions, R):
        self.sim, self.junctions, self.R = sim, junctions, R
        self.sig_map = [(ji, mi, k) for ji, j in enumerate(junctions) for mi, m in enumerate(j["members"])
                        for k, sg in enumerate(m["segs"]) if sg is not None]
        self.segs = [junctions[ji]["members"][mi]["segs"][k] for ji, mi, k in self.sig_map]
        self.pt, self.h, self.sp = np.empty((0, 2)), np.empty((0, 2)), np.empty(0)
        self.cols = ["#888"] * len(self.sig_map)

    def refresh(self):
        s = self.sim
        xy, self.sp, an = s.snapshot()
        if len(xy):
            self.pt = (xy - s.c) @ self.R
            a = np.radians(an)
            self.h = np.stack([np.sin(a), np.cos(a)], -1) @ self.R
        else:
            self.pt, self.h = np.empty((0, 2)), np.empty((0, 2))
        st = {(ji, mi): s.conn.trafficlight.getRedYellowGreenState(m["id"])
              for ji, j in enumerate(self.junctions) for mi, m in enumerate(j["members"])}
        self.cols = [SIG.get(st[ji, mi][k], "#888") if k < len(st[ji, mi]) else "#888"
                     for ji, mi, k in self.sig_map]

    def stats(self):
        fw, rv = self.sim.tt["thru_fwd"], self.sim.tt["thru_rev"]
        return (np.mean(fw) if len(fw) > 3 else None, np.mean(rv) if len(rv) > 3 else None)

    def hud_text(self):
        s, sp = self.sim, self.sp
        t = s.conn.simulation.getTime() - s.args.warmup
        mean = lambda a: f"{np.mean(a) / 60:4.1f} min" if len(a) else "   -    "
        fw, rv = s.tt["thru_fwd"], s.tt["thru_rev"]
        return (f"time {int(t // 60):02d}:{int(t % 60):02d}   cars {len(sp):3d}\n"
                f"stopped {int((sp < 0.2).sum()):3d}   avg speed {sp.mean() * 3.6 if len(sp) else 0:4.1f} km/h\n"
                f"end-to-end trip, forward {mean(fw)} (n={len(fw)})\n"
                f"end-to-end trip, reverse {mean(rv)} (n={len(rv)})")


class Panel:
    """One map (or one strip of the street). Lanes are drawn at true width (re-scaled on
    zoom), cars as true-size rectangles, signals as a stop-line across each lane."""
    def __init__(self, ax, live, edges, polys, name_re, xr, title="", hud=False, tsize=13, bg=None):
        from matplotlib.collections import LineCollection, PolyCollection
        self.ax, self.live = ax, live
        m = 150
        near = lambda sh: sh[:, 0].max() >= xr[0] - m and sh[:, 0].min() <= xr[1] + m
        ax.set_facecolor(BG)
        ax.axis("off")
        if bg is not None:                                 # scheme (road map) underneath
            ax.imshow(bg[0], extent=bg[1], aspect="auto", zorder=0, interpolation="bilinear")
        a_side, a_main = (0.55, 0.9) if bg is not None else (1.0, 1.0)
        ax.add_collection(PolyCollection([p for p in polys if near(p)], facecolors=ROAD,
                                         edgecolors="none", alpha=a_side, zorder=1))
        lanes = {0: ([], []), 1: ([], [])}
        for e in edges.values():
            k = 1 if name_re.search(e["name"]) else 0
            for lid, sh in e["lanes"].items():
                if near(sh):
                    lanes[k][0].append(sh)
                    lanes[k][1].append(e["w"][lid])
        self.cor_pts = np.vstack([sh for e in edges.values() if name_re.search(e["name"])
                                  for sh in e["lanes"].values()])
        self.lane_w = [np.array(lanes[k][1]) for k in (0, 1)]
        self.lanes = [LineCollection(lanes[k][0], colors=[ROAD, STREET][k], alpha=[a_side, a_main][k],
                                  zorder=2 + k) for k in (0, 1)]
        for lc in self.lanes:
            ax.add_collection(lc)
        self.sig_idx = [n for n, sg in enumerate(live.segs)
                        if xr[0] - m <= sg[:, 0].mean() <= xr[1] + m]
        self.stop = LineCollection([live.segs[n] for n in self.sig_idx], colors="#888", zorder=6,
                                   capstyle="butt")
        ax.add_collection(self.stop)
        self.cars = PolyCollection([], cmap=SPEED_CMAP, edgecolors="none", zorder=5)
        self.cars.set_clim(0, 14)
        ax.add_collection(self.cars)
        for n, j in enumerate(live.junctions, 1):
            if xr[0] <= j["pos"][0] <= xr[1]:
                ax.annotate(str(n), j["pos"], xytext=(0, 14), textcoords="offset points", color=BADGE_TX,
                            fontsize=8, weight="bold", ha="center", va="center", zorder=8,
                            bbox=dict(boxstyle="circle,pad=0.25", fc=BADGE_FC, ec="none"))
        ax.text(0.008, 0.95, title, transform=ax.transAxes, color=TXT, fontsize=tsize,
                weight="bold", va="top", zorder=9)
        self.hud = ax.text(0.012, 0.86, "", transform=ax.transAxes, color=TXT, fontsize=9.5, va="top",
                           linespacing=1.5, family="monospace", zorder=9) if hud else None
        self.x = xr

    def ppm(self):                                         # screen points per metre
        bb = self.ax.get_window_extent()
        x0, x1 = self.ax.get_xlim()
        return bb.width * 72 / self.ax.figure.dpi / (x1 - x0)

    def set_view(self, x0, x1, shift=0.25):
        ax = self.ax
        bb = ax.get_window_extent()
        ax.set_xlim(x0, x1)
        sel = self.cor_pts[(self.cor_pts[:, 0] >= x0) & (self.cor_pts[:, 0] <= x1)]
        ym = (sel if len(sel) else self.cor_pts)[:, 1].mean()
        half = (x1 - x0) * bb.height / bb.width / 2          # equal scale in x and y
        ax.set_ylim(ym - half * (1 - shift), ym + half * (1 + shift))
        k = self.ppm()
        for lc, w in zip(self.lanes, self.lane_w):
            lc.set_linewidths(np.maximum(w * k * 0.80, 0.5))
        self.stop.set_linewidths(max(1.4, 1.3 * k))

    def update(self):
        lv = self.live
        x0, x1 = self.ax.get_xlim()
        if len(lv.pt):
            sel = (lv.pt[:, 0] > x0 - 20) & (lv.pt[:, 0] < x1 + 20)
            pt, h = lv.pt[sel], lv.h[sel]
            k = self.ppm()
            L, Wd = max(4.5, 3.0 / k), max(1.9, 1.6 / k)
            nrm = np.stack([-h[:, 1], h[:, 0]], -1) * Wd / 2
            back = pt - h * L
            self.cars.set_verts(np.stack([back + nrm, pt + nrm, pt - nrm, back - nrm], 1))
            self.cars.set_array(lv.sp[sel])
        else:
            self.cars.set_verts([])
        self.stop.set_color([lv.cols[n] for n in self.sig_idx])
        if self.hud is not None:
            self.hud.set_text(lv.hud_text())


# ============================ main ==========================================
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", default="wave", choices=["compare", "default", "uncoordinated", "wave"])
    ap.add_argument("--osm", help="use this .osm file instead of fetching")
    ap.add_argument("--pbf", default="uzbekistan-latest.osm.pbf")
    ap.add_argument("--workdir", default="makhtumkuli_sim")
    ap.add_argument("--refresh", action="store_true", help="re-fetch the map data")
    ap.add_argument("--name-pattern", default=NAME_PATTERN)
    ap.add_argument("--cycle", type=int, default=90, help="signal cycle (s), assumed")
    ap.add_argument("--wave-kmh", type=float, default=55)
    ap.add_argument("--wave-direction", default="forward", choices=["forward", "reverse"])
    ap.add_argument("--through-vph", type=int, default=1000, help="cars/h each way, end to end")
    ap.add_argument("--side-vph", type=int, default=900, help="cars/h on other trips (all roads)")
    ap.add_argument("--warmup", type=int, default=300)
    ap.add_argument("--end", type=int, default=1800, help="animated sim seconds")
    ap.add_argument("--step", type=float, default=0.5)
    ap.add_argument("--speed", type=float, default=10, help="x real time")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", help="write a video (.mp4 or .gif) instead of showing a window")
    ap.add_argument("--window", type=float, default=600, help="view width in metres (0 = whole street)")
    ap.add_argument("--focus", type=int, default=0, help="zoom layout: junction number to start on (0 = centre)")
    ap.add_argument("--center", default="41.307206,69.302682", help="lat,lon of the centre of the simulation")
    ap.add_argument("--radius", type=float, default=1500, help="metres shown each side of the centre along the street")
    ap.add_argument("--basemap", default="osm", choices=["osm", "custom", "voyager", "light", "dark", "none"],
                    help="road-map ('scheme') background: osm = OpenStreetMap standard tiles (no key, default), custom = your own --tile-url, none = no background. "
                         "voyager / light / dark are CARTO tiles, which now need an API key (use --tile-url instead)")
    ap.add_argument("--tile-url", help="tile server template with {z}/{x}/{y}, e.g. https://api.maptiler.com/maps/streets-v2/256/{z}/{x}/{y}.png?key=YOURKEY  (use with --basemap custom)")
    ap.add_argument("--tile-attribution", default="© map provider", help="credit printed under the video when --basemap custom")
    ap.add_argument("--tile-dark", action="store_true", help="the custom tiles are dark (switches to the light-on-dark colours)")
    ap.add_argument("--dpi", type=int, default=110, help="video resolution (110 -> 1760 px wide; 240 -> 3840 px)")
    ap.add_argument("--layout", default="strips", choices=["strips", "zoom"],
                    help="strips: whole street at once as stacked rows; zoom: one window, pan with keys")
    ap.add_argument("--row-m", type=float, default=1200, help="strips: street metres per row")
    ap.add_argument("--rows", type=int, default=0, help="strips: number of rows (0 = automatic)")
    ap.add_argument("--tour", action="store_true", help="auto-pan along the street (good for videos)")
    ap.add_argument("--merge-m", type=float, default=45, help="signals closer than this = one junction")
    ap.add_argument("--macro-demand", help="traffic from the city model (json made by street_demand.py) instead of --through-vph / --side-vph")
    ap.add_argument("--demand-scale", type=float, default=1.0, help="multiply every macro-model flow (1 = the peak 15 minutes as computed)")
    ap.add_argument("--map-tol-m", type=float, default=60, help="how far (m) an entry/exit junction may be from the matching SUMO edge end")
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--stats-only", action="store_true", help="print trip-time numbers, no animation")
    args = ap.parse_args()

    name_re = re.compile(args.name_pattern, re.I)
    if args.tile_url and args.basemap == "osm":
        args.basemap = "custom"
    if args.basemap == "custom":
        if not args.tile_url:
            sys.exit("--basemap custom needs --tile-url")
        TILES["custom"] = (args.tile_url, args.tile_attribution)
    set_theme(args.basemap in ("voyager", "light", "osm") or (args.basemap == "custom" and not args.tile_dark))
    os.makedirs(args.workdir, exist_ok=True)
    lat0, lon0 = (float(v) for v in args.center.split(","))
    tag = f"{lat0:.5f}_{lon0:.5f}_{int(args.radius)}"
    f = {k: os.path.join(args.workdir, v) for k, v in dict(
        osm=f"map_{tag}.osm", net=f"map_{tag}.net.xml", trips="trips.xml", routes="routes.rou.xml").items()}
    reach = args.radius + 350                              # map kept around the centre (m)
    dlat = reach / 111_320
    dlon = reach / (111_320 * math.cos(math.radians(lat0)))
    bbox = (lat0 - dlat, lon0 - dlon, lat0 + dlat, lon0 + dlon)
    if args.osm:
        f["osm"] = args.osm
    elif args.refresh or not os.path.isfile(f["osm"]):
        if not os.path.isfile(args.pbf):                   # look in the usual places
            home = os.path.expanduser("~")
            for d in (".", "..", "~/Desktop/post", "~/Desktop", "~/Downloads", "~/Documents", "~"):
                hit = os.path.join(os.path.expanduser(d), os.path.basename(args.pbf))
                if os.path.isfile(hit):
                    print(f"Found {hit}")
                    args.pbf = hit
                    break
        if os.path.isfile(args.pbf):
            get_osm_from_pbf(args.pbf, f["osm"], name_re, bbox)
        else:
            print(f"{args.pbf} not found here -> using Overpass.")
            get_osm_from_overpass(f["osm"], name_re, bbox)
    print("Building SUMO network ...")
    build_net(f["osm"], f["net"], args.cycle)
    edges, polys = load_net(f["net"])
    cor, c, axis = corridor_axis(edges, name_re)
    make_demand(edges, cor, c, axis, args, f)

    traci = import_traci()
    plans = ["uncoordinated", "wave"] if args.mode == "compare" else [args.mode]
    labels = {"default": "Netconvert default (all in sync)", "uncoordinated": "Lights NOT coordinated",
              "wave": f"GREEN WAVE {args.wave_kmh:.0f} km/h"}
    sims = [Sim(traci, f"s{i}", p, f, args) for i, p in enumerate(plans)]
    raw = find_junctions(sims[0].conn, edges, name_re, args.merge_m)
    if not raw:
        sys.exit("No traffic lights found on the street in the network.")
    junctions = merge_junctions(raw, c, axis, args.merge_m)
    cxy = np.array(sims[0].conn.simulation.convertGeo(lon0, lat0, True))      # centre in net coordinates
    xc = float((cxy - c) @ axis)                                              # its position along the street
    shown = [j for j in junctions if abs(float((j["pos"] - c) @ axis) - xc) <= args.radius]
    if not shown:
        sys.exit("No signalised junction within the radius of the centre point - raise --radius.")
    bearing = (math.degrees(math.atan2(axis[0], axis[1])) + 360) % 360
    print(f"\n{len(raw)} signal nodes on the street -> {len(junctions)} junctions, "
          f"{len(shown)} within {args.radius:.0f} m of the centre "
          f"(numbered along the street; 'forward' = toward {bearing:.0f}° compass):")
    for n, j in enumerate(shown, 1):
        gap = "" if n == 1 else f"  +{np.linalg.norm(j['pos'] - shown[n - 2]['pos']):.0f} m"
        names = ", ".join(sorted(j["cross"])) or "(cross street unnamed in OSM)"
        print(f"  {n:2d}  {names}{gap}   [{len(j['members'])} signal node(s)]")
    for s in sims:
        apply_plan(s.conn, junctions, s.plan, axis, c, args)
    print(f"Warm-up {args.warmup}s ...")
    for _ in range(int(args.warmup / args.step)):
        for s in sims:
            s.step()

    if args.stats_only:                      # numbers only, no window / video
        for _ in range(int(args.end / args.step)):
            for s in sims:
                s.step()
        print(f"\nEnd-to-end trip times after {args.end}s (mean, n):")
        for s in sims:
            f_, r_, o_ = s.tt["thru_fwd"], s.tt["thru_rev"], s.tt["other"]
            m_ = lambda a: f"{np.mean(a):5.0f}s" if len(a) else "    - "
            print(f"  {labels[s.plan]:28s} end-to-end forward {m_(f_)} (n={len(f_)})   reverse {m_(r_)} (n={len(r_)})   other trips {m_(o_)} (n={len(o_)})")
        for s in sims:
            s.close()
        return
    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    perp = np.array([-axis[1], axis[0]])
    R = np.stack([axis, perp], 1)
    T = lambda a: (np.asarray(a) - c) @ R
    edges_t = {i: dict(e, lanes={k: T(v) for k, v in e["lanes"].items()}) for i, e in edges.items()}
    polys_t = [T(p) for p in polys]
    junc_t = [dict(j, pos=T(j["pos"]), members=[dict(m, segs=[None if q is None else T(q) for q in m["segs"]])
                                                for m in j["members"]]) for j in shown]
    for s in sims:
        s.c = c
    n = len(sims)
    lives = [Live(s, junc_t, R) for s in sims]
    jx = [float(j["pos"][0]) for j in junc_t]
    xs = np.vstack([p for e in edges_t.values() if name_re.search(e["name"]) for p in e["lanes"].values()])[:, 0]
    full = (xc - args.radius, xc + args.radius)
    bg = None
    if args.basemap != "none":
        ys = np.vstack([p for e in edges_t.values() if name_re.search(e["name"]) for p in e["lanes"].values()])[:, 1]
        bg = build_basemap(args.basemap, bbox, sims[0].conn, T, full, (ys.min(), ys.max()),
                           os.path.join(args.workdir, "tiles"))
    foot = fig = None
    panels, headers, view = [], [], {}
    if args.layout == "strips":
        total = full[1] - full[0]
        rows = args.rows or max(1, round(total / args.row_m))
        seg = total / rows
        ratios = ([1.0] * n + [0.3]) * rows
        fig = plt.figure(figsize=(16, 10), facecolor=BG)
        gs = fig.add_gridspec(len(ratios) - 1, 1, height_ratios=ratios[:-1], left=0.005, right=0.995,
                              top=0.84, bottom=0.07, hspace=0.0)
        for r in range(rows):
            x0 = full[0] + r * seg
            for i, lv in enumerate(lives):
                ax = fig.add_subplot(gs[r * (n + 1) + i])
                pn = Panel(ax, lv, edges_t, polys_t, name_re, (x0, x0 + seg),
                           title=(labels[lv.sim.plan] if n > 1 else "") +
                                 (f"   {(x0 - full[0]) / 1000:.1f}-{(x0 + seg - full[0]) / 1000:.1f} km"
                                  if i == 0 or n == 1 else ""), tsize=9, bg=bg)
                pn.set_view(x0, x0 + seg, shift=0)
                panels.append(pn)
        for i, lv in enumerate(lives):
            xh = 0.012 + i * (0.5 if n == 2 else 0)
            fig.text(xh, 0.985, labels[lv.sim.plan], color=TXT, fontsize=15, weight="bold", va="top")
            headers.append((lv, fig.text(xh, 0.945, "", color=TXT, fontsize=10, va="top",
                                         linespacing=1.5, family="monospace")))
        print(f"Layout: {rows} strips of {seg:.0f} m, {len(jx)} junctions all visible.")
    else:
        fig, axs = plt.subplots(n, 1, figsize=(16, 9 if n == 2 else 5.6), facecolor=BG)
        axs = np.atleast_1d(axs)
        fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.075, hspace=0.01)
        panels = [Panel(a, lv, edges_t, polys_t, name_re, full, title=labels[lv.sim.plan], hud=True, bg=bg)
                  for a, lv in zip(axs, lives)]
        view = {"cx": xc if args.focus <= 0 else jx[min(args.focus, len(jx)) - 1],
                "w": args.window if 0 < args.window < full[1] - full[0] else full[1] - full[0]}

        def show_view():
            if view["w"] >= full[1] - full[0] - 1:
                lo, hi = full
            else:
                lo, hi = view["cx"] - view["w"] / 2, view["cx"] + view["w"] / 2
            for p in panels:
                p.set_view(lo, hi)
        show_view()

        def on_key(ev):
            k = ev.key
            if k in ("right", "left"):
                view["cx"] += (1 if k == "right" else -1) * view["w"] / 2
            elif k in ("+", "=", "up"):
                view["w"] = max(120, view["w"] / 1.5)
            elif k in ("-", "down"):
                view["w"] = min(full[1] - full[0], view["w"] * 1.5)
            elif k in ("n", "p"):
                nxt = [x for x in jx if (x > view["cx"] + 5 if k == "n" else x < view["cx"] - 5)]
                if nxt:
                    view["cx"] = min(nxt) if k == "n" else max(nxt)
            elif k == "a":
                view["w"] = full[1] - full[0]
            else:
                return
            show_view()
            fig.canvas.draw_idle()
        fig.canvas.mpl_connect("key_press_event", on_key)
        if not args.no_show and not args.save:
            fig.text(0.995, 0.004, "← →  pan    + −  zoom    n / p  next / previous junction    a  whole street",
                     color=MUTED, ha="right", fontsize=8)
    foot = fig.text(0.5, 0.036, "", color=TXT, ha="center", fontsize=11, weight="bold")
    fig.text(0.5, 0.006, "Simulation (SUMO) on OpenStreetMap geometry. Signal timing, volumes and 55 km/h "
             "are assumptions, not measurements." + (f"   Map: {TILES[args.basemap][1]}" if bg is not None else ""), color=MUTED, ha="center", fontsize=7.5)
    spf = max(1, round(args.speed * 0.05 / args.step))
    nframes = int(args.end / (spf * args.step))

    def frame(i):
        if args.tour and args.layout == "zoom":
            view["cx"] = jx[0] + (jx[-1] - jx[0]) * i / max(1, nframes - 1)
            show_view()
        for s in sims:
            for _ in range(spf):
                s.step()
        for lv in lives:
            lv.refresh()
        for p in panels:
            p.update()
        for lv, t in headers:
            t.set_text(lv.hud_text())
        res = [lv.stats() for lv in lives]
        if n == 2 and all(x is not None for r in res for x in r):
            f_, r_ = res[0][0] - res[1][0], res[0][1] - res[1][1]
            word = lambda d: f"{'saves' if d > 0 else 'costs'} {abs(d):.0f} s"
            foot.set_text(f"Green wave vs uncoordinated, per car:  forward direction {word(f_)}"
                          f"   |   reverse direction {word(r_)}")
        return []
    anim = FuncAnimation(fig, frame, frames=nframes, interval=50, repeat=False, blit=False)
    if args.save:
        from matplotlib import animation
        writer = (animation.FFMpegWriter(fps=20, bitrate=4000) if args.save.endswith(".mp4")
                  and animation.writers.is_available("ffmpeg") else animation.PillowWriter(fps=20))
        print(f"Rendering {nframes} frames to {args.save} ...")
        anim.save(args.save, writer=writer, dpi=args.dpi, savefig_kwargs={"facecolor": BG})
        print("Saved", args.save)
    if not args.no_show and not args.save:
        plt.show()
    for s in sims:
        try:
            s.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()