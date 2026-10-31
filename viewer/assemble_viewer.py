"""Build viewer/index.html from the template and the packed data.
usage: python3 assemble_viewer.py viewer_template.html viewer_data.json index.html [trips.json]
(trips.json from pipeline/build_trips.py is optional; without it the trip layer is hidden)"""
import sys
t, d, o = sys.argv[1:4]; tr = sys.argv[4] if len(sys.argv) > 4 else None
open(o, 'w', encoding='utf-8').write(open(t, encoding='utf-8').read().replace('__DATA__', open(d, encoding='utf-8').read()).replace('__TRIPS__', open(tr, encoding='utf-8').read() if tr else 'null'))
print('wrote', o)
