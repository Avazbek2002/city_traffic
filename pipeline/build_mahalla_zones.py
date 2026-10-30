"""Join place assignments (from build_attraction.py) to mahalla polygons -> one row per mahalla (OD zone).
Usage: python3 build_mahalla_zones.py places_assigned.csv MAHALLAS.geojson OUT.csv
"""
import sys, json, numpy as np, pandas as pd
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
PL, MAH, OUT = sys.argv[1:4]
m = json.load(open(MAH))
props = pd.DataFrame([f['properties'] for f in m['features']])
geoms = [shape(f['geometry']) for f in m['features']]
tree = STRtree(geoms)
props['area_km2'] = [g.area * 111.32 * 111.32 * np.cos(np.radians(41.3)) for g in geoms]
props['centroid_lon'] = [g.centroid.x for g in geoms]; props['centroid_lat'] = [g.centroid.y for g in geoms]

df = pd.read_csv(PL)
mid = []
for x, y in zip(df.lon, df.lat):
    h = tree.query(Point(x, y), predicate='within')
    mid.append(props.mahalla_id[h[0]] if len(h) else np.nan)
df['mahalla_id'] = mid
a = df[df.mahalla_id.notna() & (df.weight > 0)]
cnt = a.pivot_table(index='mahalla_id', columns='purpose', values='placeId', aggfunc='count', fill_value=0).add_suffix('_count')
wt = a.pivot_table(index='mahalla_id', columns='purpose', values='weight', aggfunc='sum', fill_value=0).add_suffix('_weighted')
bz = a[a.subtype == 'bazaar'].groupby('mahalla_id').size().rename('bazaar_count')
z = props.set_index('mahalla_id').join([cnt, wt, bz]).fillna(0).reset_index()
for c in [c for c in z.columns if c.endswith(('_count', '_weighted'))]:
    z[c] = z[c].astype(float if c.endswith('_weighted') else int)

# zone-level inputs for the gravity model (each column sums to 1 over the city)
z['production_share'] = z.population / z.population.sum()                 # trip production: residents
z['attr_work_share'] = z.employment / z.employment.sum()                  # OBSERVED jobs (mahalla file)
z['attr_work_2gis_share'] = z.work_weighted / z.work_weighted.sum()       # 2GIS-derived, for comparison
for p in ['education', 'health', 'retail', 'industry']:
    z[f'attr_{p}_share'] = z[f'{p}_weighted'] / z[f'{p}_weighted'].sum()
z['pop_density_per_km2'] = (z.population / z.area_km2).round(0)
z.round(5).to_csv(OUT, index=False)
df.to_csv(OUT.replace('.csv', '_places.csv'), index=False)   # places with mahalla_id, read by the student/university scripts
print(len(z), 'mahallas; places inside a mahalla:', int(a.shape[0]), '/', int((df.weight > 0).sum()), 'with weight>0')
print(z[['population', 'employment']].sum().to_dict())
print(z.filter(like='_share').sum().round(4).to_dict())
