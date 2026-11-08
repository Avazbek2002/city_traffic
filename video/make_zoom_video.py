"""Render the "city to a single lane" clip (1080x1350, 30 fps, ~37 s): the whole-city load map, then zoom in to one congested junction at lane level, zoom out, zoom in to the next ...
usage: python3 make_zoom_video.py viewer_data.json SCENES_DIR OUT.mp4 [first_frame last_frame] [--order bunyodkor,ahmad_donish,qodiriy,furqat]
 viewer_data.json : prep_viewer_data.py output of the base run (city road loads); SCENES_DIR: site/scenes (lane-level street files made by micro/scene_sim.py)
The junction shown for each street = the traffic signal inside the 1.5 km window with the most stopped cars (speed < 5.4 km/h) during the recorded 5 minutes.
Needs: playwright (chromium) and ffmpeg. The optional frame range renders only part of the film (for tests)."""
import sys, os, json, base64, math, subprocess, tempfile, shutil, threading, http.server, functools, numpy as np
from playwright.sync_api import sync_playwright
args = [a for a in sys.argv[1:] if not a.startswith('--')]; order = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--order=')), 'bunyodkor,ahmad_donish,qodiriy,furqat').split(',')
VD, SD, OUT = args[:3]; rng = [int(v) for v in args[3:5]] if len(args) > 4 else None
here = os.path.dirname(os.path.abspath(__file__)); tmp = tempfile.mkdtemp(prefix='zoomvid_'); kc = math.cos(math.radians(41.3)); LON0, LAT0 = 69.2, 41.3
b64 = lambda s, t: np.frombuffer(base64.b64decode(s), t)
# ---- city layer: road chunks (global metres, y north) and load class at the slices that are needed
D = json.load(open(VD)); Q = D['q']; pts = b64(D['pts'], np.uint16).astype(np.float64); off = b64(D['off'], np.uint32); N = D['nch']; NS = D['ns']
flow = b64(D['flow'], np.uint16).reshape(NS, N); cap = b64(D['cap'], np.uint16)
gx = ((D['lon0'] + pts[0::2] * Q) - LON0) * 111320 * kc; gy = ((D['lat0'] + pts[1::2] * Q) - LAT0) * 111320
scenes = {sid: json.load(open(f'{SD}/{sid}.json')) for sid in order}; pks = sorted({int(round((sum(int(x) * m for x, m in zip(s['clock'].split(':'), (60, 1))) - 360) / 15)) for s in scenes.values()})
def cls(sl):
    f = flow[sl].astype(float); vc = f / np.maximum(cap, 1); return np.where(f < 15, 5, np.where(vc < .3, 0, np.where(vc < .6, 1, np.where(vc < .85, 2, np.where(vc < 1, 3, 4))))).astype(int).tolist()
rings = [[[round(((D['lon0'] + x * Q) - LON0) * 111320 * kc), round(((D['lat0'] + y * Q) - LAT0) * 111320)] for x, y in r['p']] for r in D['rings']]
city = dict(off=off.tolist(), x=np.round(gx).astype(int).tolist(), y=np.round(gy).astype(int).tolist(), cls={str(p): cls(p) for p in pks}, rings=rings)
json.dump(city, open(f'{tmp}/city.json', 'w'), separators=(',', ':'))
# ---- the junction of every street: the signal with the most stopped cars
tg = []
for sid in order:
    s = scenes[sid]; C = s['cars']; co = b64(C['off'], np.uint32); dx = b64(C['x'], np.int16).astype(np.int64); dy = b64(C['y'], np.int16).astype(np.int64); v = b64(C['v'], np.uint8) / 3.0
    X = np.zeros(len(dx)); Y = np.zeros(len(dy))
    for k in range(C['n']): a, e = co[k], co[k+1]; X[a:e] = np.cumsum(dx[a:e]) / 10; Y[a:e] = np.cumsum(dy[a:e]) / 10
    st = v < 1.5; best = (-1, None)
    for tx, ty in s['tls']:
        if abs(tx) > 600 or abs(ty) > 600: continue
        n_ = int((st & (np.hypot(X - tx, Y - ty) < 60)).sum())
        if n_ > best[0]: best = (n_, (tx, ty))
    cl, cn = s['centre']; ser = s['macro']['series']; pk = int(round((int(s['clock'][:2]) * 60 + int(s['clock'][3:]) - 360) / 15)); name = s['label'].replace(' (for comparison)', '')
    first = name.split()[0].lower()[:6]; fl = max((r for r in s['macro']['streets'] if first in r['street'].lower()), key=lambda r: r['flow_vph'], default=dict(flow_vph=ser[pk]['flow_vph']))['flow_vph']
    tg.append(dict(id=sid, label=name, clock=s['clock'], pk=pk, gx=(cn - LON0) * 111320 * kc, gy=(cl - LAT0) * 111320, jx=best[1][0], jy=best[1][1], stopped=best[0], flow_vph=fl, rec=s['record'], window_cars=s['macro']['window_total_cars']))
    print(f"{sid}: junction at ({best[1][0]:.0f}, {best[1][1]:.0f}) m from the window centre, {best[0]} stopped car-seconds; city model {fl} cars/h")
    shutil.copy(f'{SD}/{sid}.json', f'{tmp}/{sid}.json')
json.dump(tg, open(f'{tmp}/targets.json', 'w')); shutil.copy(f'{here}/zoom_template.html', f'{tmp}/stage.html')
srv = http.server.ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(http.server.SimpleHTTPRequestHandler, directory=tmp)); srv.RequestHandlerClass.log_message = lambda *a, **k: None
threading.Thread(target=srv.serve_forever, daemon=True).start(); port = srv.server_address[1]
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1080, 'height': 1350}); errs = []; pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto(f'http://127.0.0.1:{port}/stage.html'); pg.wait_for_function('window.READY === true', timeout=120000); n = pg.evaluate('TOTALF'); a, z = rng if rng else (0, n)
    if os.environ.get('FRAMES'):                 # test: FRAMES=0,90,150 FRAMES_DIR=/path -> only those frames as png
        for f in [int(v) for v in os.environ['FRAMES'].split(',')]: pg.evaluate(f'renderFrame({f})'); pg.screenshot(path=f"{os.environ['FRAMES_DIR']}/z{f:05d}.png", type='png')
        b.close(); sys.exit(0)
    for f in range(a, z):
        pg.evaluate(f'renderFrame({f})'); pg.screenshot(path=f'{tmp}/f{f - a:05d}.png', type='png')
        if f % 100 == 0: print('frame', f, '/', z, flush=True)
    b.close()
if errs: print('page errors:', errs); sys.exit(1)
subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '30', '-i', f'{tmp}/f%05d.png', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '17', '-preset', 'medium', '-movflags', '+faststart', OUT], check=True)
shutil.rmtree(tmp); print('wrote', OUT)
