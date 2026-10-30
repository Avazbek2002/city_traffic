"""Build viewer/index.html from the template and the packed data.
usage: python3 assemble_viewer.py viewer_template.html viewer_data.json index.html"""
import sys
t, d, o = sys.argv[1:4]
open(o, 'w').write(open(t).read().replace('__DATA__', open(d).read()))
print('wrote', o)
