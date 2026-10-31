"""Compare sensitivity runs. usage: python3 summarize_sensitivity.py data/work/sens  (writes sensitivity_summary.csv there)"""
import sys, os, pandas as pd
D = sys.argv[1]; rows = []
for r in sorted(os.listdir(D)):
    f = f'{D}/{r}/slice_stats.csv'
    if not os.path.isfile(f): continue
    s = pd.read_csv(f); h = lambda c, a, b: s[c].iloc[a:b].mean()
    rows.append(dict(run=r, trips_0708=h('veh_per_hour_rate', 4, 8), trips_0809=h('veh_per_hour_rate', 8, 12), ratio_0809_0708=h('veh_per_hour_rate', 8, 12)/h('veh_per_hour_rate', 4, 8),
                     peak_veh_on_road=s.veh_on_network.max(), peak_clock=s.clock[s.veh_on_network.idxmax()], min_speed_kmh=s.avg_speed_kmh.min(), max_edges_over_cap=s.edges_vc_gt1.max(),
                     mean_min_work=s.mean_min_work.mean()))
t = pd.DataFrame(rows).round(2)
if 'base' in set(t.run):
    b = t[t.run == 'base'].iloc[0]
    for c in ['trips_0708', 'trips_0809', 'peak_veh_on_road', 'max_edges_over_cap']: t[c + '_vs_base_%'] = ((t[c]/b[c]-1)*100).round(0)
t.to_csv(f'{D}/sensitivity_summary.csv', index=False); print(t.T.to_string())
