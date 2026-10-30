"""Street-level lane checklist: the N busiest streets at the peak slice, lanes used in the model, check links, input column.
usage: python3 build_lane_checklist.py road_edges.csv.gz slices_flow.npz OUT.xlsx [N=100] [slice=11]
Fill the yellow column, then run the model with LANES_FILE=OUT.xlsx (congestion_loop.py / congestion_slices.py).
The override replaces ASSUMED lane counts only (road sections whose lane count is tagged in OpenStreetMap keep their tag)."""
import sys, numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
EF, SF, OUT = sys.argv[1:4]; N = int(sys.argv[4]) if len(sys.argv) > 4 else 100; SL = int(sys.argv[5]) if len(sys.argv) > 5 else 11
e = pd.read_csv(EF); e = e[~e.highway.isin(['service', 'living_street'])].reset_index(drop=True)      # same filtering and ordering as congestion_loop.py
ids = pd.unique(np.r_[e.u.values, e.v.values]); ix = pd.Series(np.arange(len(ids)), index=ids); e['a'] = ix[e.u].values; e['b'] = ix[e.v].values
e['t0'] = e.length_m / (e.speed_kmh / 3.6); e = e.sort_values(['a', 'b', 't0']).drop_duplicates(['a', 'b']).reset_index(drop=True)
d = np.load(SF, allow_pickle=True); assert len(e) == len(d['way_id']) and (e.way_id.values == d['way_id']).all() and np.allclose(e.length_m.values, d['length']), 'edge order mismatch'
e['flow'] = d['flow'][SL]; e['cap'] = d['cap']; e['vc'] = e.flow / e.cap; e['vkt'] = e.flow * e.length_m / 1000; e['tag'] = (e.lanes_src == 'tag')
e['street'] = e.names.fillna('').str.split('|').str[0]; tot = e.vkt.sum()
st = e[e.street != '']; rows = []
for s, g in st.groupby('street'):
    w = g.groupby('way_id').agg(vkt=('vkt', 'sum'), L=('length_m', 'sum'), fl=('flow', 'mean'), vc=('vc', 'mean'), lon=('u_lon', 'mean'), lat=('u_lat', 'mean'))
    b = w.vkt.idxmax(); gb = g[g.way_id == b]; j = ((gb.u_lon - w.lon[b]) ** 2 + (gb.u_lat - w.lat[b]) ** 2).idxmin()
    nm = [n for n in dict.fromkeys(g.names.iloc[0].split('|')) if n and n != s][:3]
    lo, hi = int(g.lanes.min()), int(g.lanes.max())
    rows.append(dict(street=s, other=' | '.join(nm), hw=', '.join(g.groupby('highway').length_m.sum().sort_values(ascending=False).index[:2]), km=g.length_m.sum() / 1000, share=100 * g.vkt.sum() / tot, vkt=g.vkt.sum(),
                     flow=float(w.fl[b]), vc=float(w.vc[b]), lanes=f'{lo}' if lo == hi else f'{lo}–{hi}', tagged=100 * g.length_m[g.tag].sum() / g.length_m.sum(), nways=len(w), way=int(b), lat=float(g.u_lat[j]), lon=float(g.u_lon[j])))
r = pd.DataFrame(rows).sort_values('vkt', ascending=False).head(N).reset_index(drop=True); r['cum'] = r.share.cumsum()
print(f'top {N} streets = {r.share.sum():.1f}% of vehicle-km at the peak; lane count mostly tagged (>=50%) on {(r.tagged>=50).sum()} of them; ways {int(r.nways.sum())}')
F = lambda **k: Font(name='Arial', **{'size': 10, **k}); wb = Workbook(); ws = wb.active; ws.title = 'Streets'
H = ['Rank', 'Street name (as in OpenStreetMap)', 'Other names', 'Road class', 'Length, all directions (km)', 'Share of city vehicle-km (%)', 'Cumulative share (%)', 'Peak flow, busiest section (veh/h)', 'Load vs capacity, busiest section',
     'Lanes per direction in model', 'Length with lane tag in OSM (%)', 'Road sections', 'OSM way ID (busiest section)', 'Google satellite', 'Yandex satellite', 'Mapillary', 'Latitude', 'Longitude', 'Typical lanes per direction (main part)', 'Notes (lanes vary? tunnel? reversible?)']
ws.append(H); yel = PatternFill('solid', fgColor='FFFF99'); thin = Side(style='thin', color='D0D0D0')
for i, x in r.iterrows():
    ws.append([i + 1, x.street, x.other, x.hw, round(x.km, 1), round(x.share, 2), round(x.cum, 1), round(x.flow), round(x.vc, 2), x.lanes, round(x.tagged), x.nways, x.way, 'open', 'open', 'open', round(x.lat, 6), round(x.lon, 6), None, None]); n = ws.max_row
    for col, url in ((13, f'https://www.openstreetmap.org/way/{x.way}'), (14, f'https://www.google.com/maps/@?api=1&map_action=map&center={x.lat:.6f},{x.lon:.6f}&zoom=19&basemap=satellite'),
                     (15, f'https://yandex.com/maps/?ll={x.lon:.6f},{x.lat:.6f}&z=19&l=sat'), (16, f'https://www.mapillary.com/app/?lat={x.lat:.6f}&lng={x.lon:.6f}&z=18')):
        c = ws.cell(n, col); c.hyperlink = url; c.font = F(color='0563C1', underline='single')
    ws.cell(n, 19).fill = yel; ws.cell(n, 20).fill = yel
for c in ws[1]: c.font = F(bold=True, color='FFFFFF'); c.fill = PatternFill('solid', fgColor='1F3A5F'); c.alignment = Alignment(wrap_text=True, vertical='center')
for row in ws.iter_rows(min_row=2):
    for c in row:
        if c.column not in (13, 14, 15, 16): c.font = F()
        c.border = Border(bottom=thin)
for k, wd in enumerate([6, 34, 40, 18, 12, 12, 12, 13, 13, 12, 13, 10, 14, 14, 14, 11, 11, 11, 18, 40], 1): ws.column_dimensions[get_column_letter(k)].width = wd
ws.row_dimensions[1].height = 62; ws.freeze_panes = 'C2'; ws.auto_filter.ref = f'A1:T{ws.max_row}'
L = wb.create_sheet('How to use', 0)
lines = [('How to use this checklist', True), ('', False),
 ('Purpose: the traffic model guesses lane counts where OpenStreetMap has none. These are the streets that carry the most vehicle-km at the 08:45 peak, so a wrong lane count here changes the results most.', False),
 (f'The {N} streets below carry {r.share.sum():.0f}% of all vehicle-km in that 15-minute step. Each street is made of many short road sections (the "Road sections" column).', False), ('', False),
 ('1. Open the street with the Google, Yandex or Mapillary link. The point is on its busiest section.', False),
 ('2. Count the lanes in ONE direction of travel along the main part of the street. Ignore parking, bus-only lanes and short turn pockets.', False),
 ('3. Type that number in the yellow column "Typical lanes per direction". Leave it empty if you are not sure: the model then keeps its own value.', False),
 ('4. For a one-way street, enter its total number of lanes.', False),
 ('5. If the lane count changes along the street (for example a ring road), enter the most common value and describe the change in Notes. Tell me the section and I can add a section-specific override.', False),
 ('6. The model applies your number to every section of that street whose lane count is not already tagged in OpenStreetMap. Tagged sections keep their tag.', False), ('', False),
 ('Example (not a real entry): row "Amir Temur shoh ko\'chasi", lanes in model 2, typical lanes per direction: 3, notes: "bus lane on the right, 2 lanes near the metro".', False), ('', False),
 ('Where the model values come from: a lane count is "tagged" when it is in OpenStreetMap; otherwise it is my default by road class (trunk 3, primary 2, secondary 2, all others 1 per direction).', False),
 ('Flow and load are model results for the 08:45-09:00 step, not measurements. Capacity per lane is also an assumption.', False)]
for t, b in lines: L.append([t]); c = L.cell(L.max_row, 1); c.font = F(bold=b, size=12 if b else 10); c.alignment = Alignment(wrap_text=True, vertical='top')
L.column_dimensions['A'].width = 120; wb.save(OUT); print('saved', OUT)
