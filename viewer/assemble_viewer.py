"""Build a viewer page from a template and the packed data.
usage: python3 assemble_viewer.py TEMPLATE.html viewer_data.json OUT.html [trips.json [od_graph.json [scenes_index.json]]]
 roads page      : viewer_template.html  viewer_data.json  index.html  [trips.json]
 zone graph page : od_graph_template.html od_graph.json    od_graph.html
 dashboard (all) : dashboard_template.html viewer_data.json dashboard.html trips.json od_graph.json
(trips.json from pipeline/build_trips.py and od_graph.json from pipeline/build_od_graph.py; a missing trips file hides the trip layer)"""
import sys
t, d, o = sys.argv[1:4]; tr = sys.argv[4] if len(sys.argv) > 4 else None; gr = sys.argv[5] if len(sys.argv) > 5 else None; sc = sys.argv[6] if len(sys.argv) > 6 else None
rd = lambda f: open(f, encoding='utf-8').read()
html = rd(t).replace('__DATA__', rd(d)).replace('__TRIPS__', rd(tr) if tr else 'null').replace('__GRAPH__', rd(gr) if gr else 'null').replace('__SCENES__', rd(sc) if sc else 'null')
open(o, 'w', encoding='utf-8').write(html)
print('wrote', o)
