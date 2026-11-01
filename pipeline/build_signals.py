"""Traffic-signal nodes from an OSM .pbf (inside the mahalla bounding box).
usage: python3 build_signals.py country.pbf MAHALLAS.geojson OUT_signals.csv
type 'junction'  = highway=traffic_signals (signalised intersection)
type 'crossing' = pedestrian crossing with signals (crossing=traffic_signals, no junction tag on the node)"""
import sys, json, osmium, pandas as pd
pbf, mah, out = sys.argv[1:4]
xs, ys = [], []
for f in json.load(open(mah))['features']:
    g = f['geometry']; st = [g['coordinates']] if g['type'] == 'Polygon' else g['coordinates']
    for poly in st:
        for x, y in poly[0]: xs.append(x); ys.append(y)
B = 0.03; x0, x1, y0, y1 = min(xs)-B, max(xs)+B, min(ys)-B, max(ys)+B
rows = []
class H(osmium.SimpleHandler):
    def node(self, n):
        t = n.tags
        if t.get('highway') == 'traffic_signals': typ = 'junction'
        elif t.get('crossing') == 'traffic_signals' or t.get('crossing:signals') == 'yes': typ = 'crossing'
        else: return
        if x0 <= n.location.lon <= x1 and y0 <= n.location.lat <= y1: rows.append((n.id, n.location.lon, n.location.lat, typ))
H().apply_file(pbf)
d = pd.DataFrame(rows, columns=['node', 'lon', 'lat', 'type']); d.to_csv(out, index=False)
print(len(d), 'signal nodes:', d.type.value_counts().to_dict())
