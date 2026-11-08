"""Pack one work-start-time scenario for the dashboard's scenario switch (lazy-loaded file scenarios/ID.json).
usage: python3 prep_scenario.py ID "NOTE" viewer_data_scn.json od_graph_scn.json trips_scn.json MAHALLAS.geojson CLASSES.csv OUT.json
 viewer_data_scn.json: prep_viewer_data.py output of the scenario run made with KEEP_FROM=<base slices_flow.npz> (same road chunks as the base data); od_graph / trips: made with
 WORK_PROFILE=zones WORK_CLASS_FILE=CLASSES.csv; CLASSES.csv: zone_row, mahalla_id, start_min (pipeline/make_ring_starts.py or optimize_starts.py).
The file also carries the mahalla outlines coloured by start time (the 'areas')."""
import sys, json, pandas as pd
from shapely.geometry import shape
sid, note, VD, GR, TRF, MG, CL, OUT = sys.argv[1:9]
d = json.load(open(VD)); b = json.load(open('/dev/null')) if False else None
g = json.load(open(GR)); t = json.load(open(TRF)); cl = pd.read_csv(CL); st = dict(zip(cl.mahalla_id, cl.start_min)); mg = json.load(open(MG)); lon0, lat0, Q = d['lon0'], d['lat0'], d['q']
areas = []
for f in mg['features']:
    s = st.get(f['properties']['mahalla_id'])
    if s is None: continue
    geom = shape(f['geometry']).simplify(0.0003); polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
    for p in polys:
        if p.is_empty: continue
        areas.append(dict(s=int(s), p=[[int(round((x-lon0)/Q)), int(round((y-lat0)/Q))] for x, y in p.exterior.coords]))
out = dict(id=sid, note=note, flow=d['flow'], ratio=d['ratio'], stats=d['stats'], graph=dict(dist=dict(present=g['dist']['present'], flow=g['dist']['flow']), mah=dict(present=g['mah']['present'], links=g['mah']['links']), transit=g['transit'], departed=g['departed']), trips=t, areas=areas)
json.dump(out, open(OUT, 'w'), separators=(',', ':')); import os; print('areas', len(areas), 'MB', round(os.path.getsize(OUT)/1e6, 2))
