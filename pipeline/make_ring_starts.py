"""Ring-based work start times. usage: python3 make_ring_starts.py mahalla_zones.csv OUT_prefix [bandwidth_km]
Centre = the point where jobs are densest (job-weighted kernel density over the mahalla centroids, bandwidth 2 km by default).
Rings = four bands of distance from that centre, each holding 25% of all jobs. Writes OUT_prefix_outer_early.csv (centre 10:00 ... edge 08:30) and OUT_prefix_centre_early.csv (centre 08:30 ... edge 10:00)."""
import sys, numpy as np, pandas as pd
z = pd.read_csv(sys.argv[1]); out = sys.argv[2]; bw = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
kc = np.cos(np.radians(41.3)); x = z.centroid_lon.values*111.32*kc; y = z.centroid_lat.values*111.32; w = z.employment.values.astype(float)
dens = np.array([(w*np.exp(-0.5*(((x-xi)**2+(y-yi)**2)/bw**2))).sum() for xi, yi in zip(x, y)]); c = int(np.argmax(dens)); cx, cy = x[c], y[c]
r = np.hypot(x-cx, y-cy); order = np.argsort(r); cum = np.cumsum(w[order])/w.sum(); ring = np.empty(len(z), int); ring[order] = np.minimum(3, (cum*4-1e-9).astype(int))
print(f'job-density centre: mahalla {z.mahalla_name[c]} (district {z.district_id[c]}), lon {z.centroid_lon[c]:.4f} lat {z.centroid_lat[c]:.4f}; job-weighted mean point lon {np.average(z.centroid_lon, weights=w):.4f} lat {np.average(z.centroid_lat, weights=w):.4f}')
print('ring outer radius (km):', [round(float(r[ring == k].max()), 1) for k in range(4)], '| zones', [int((ring == k).sum()) for k in range(4)], '| jobs share', [round(float(w[ring == k].sum()/w.sum()), 2) for k in range(4)])
for name, st in (('outer_early', [600, 570, 540, 510]), ('centre_early', [510, 540, 570, 600])):
    pd.DataFrame(dict(zone_row=np.arange(len(z)), mahalla_id=z.mahalla_id.values, start_min=np.array(st)[ring])).to_csv(f'{out}_{name}.csv', index=False)
