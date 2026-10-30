"""Drivable road graph from an OSM .pbf, clipped to the mahalla bounding box.
usage: python3 build_graph.py country.pbf MAHALLAS.geojson OUT_EDGES.csv.gz
Edges are directed way segments (consecutive OSM nodes). speed_kmh / lanes are tag values or ASSUMED class defaults (flagged by *_src)."""
import sys, json, math, numpy as np, pandas as pd, osmium
pbf, mah, out = sys.argv[1:4]
xs=[];ys=[]
for f in json.load(open(mah))['features']:
    g=f['geometry']; st=[g['coordinates']] if g['type']=='Polygon' else g['coordinates']
    for poly in st:
        for x,y in poly[0]: xs.append(x); ys.append(y)
B=0.03; x0,x1,y0,y1=min(xs)-B,max(xs)+B,min(ys)-B,max(ys)+B; print('bbox',x0,y0,x1,y1)
DRIVE={'motorway','trunk','primary','secondary','tertiary','unclassified','residential','living_street','service',
       'motorway_link','trunk_link','primary_link','secondary_link','tertiary_link'}
SPEED={'motorway':80,'trunk':60,'primary':50,'secondary':45,'tertiary':40,'unclassified':30,'residential':30,'living_street':15,'service':20}  # ASSUMPTION km/h
LANES={'motorway':3,'trunk':3,'primary':2,'secondary':2,'tertiary':1,'unclassified':1,'residential':1,'living_street':1,'service':1}              # ASSUMPTION per direction
def hav(a,b,c,d):
    p=math.pi/180; h=math.sin((c-a)*p/2)**2+math.cos(a*p)*math.cos(c*p)*math.sin((d-b)*p/2)**2
    return 12742000*math.asin(math.sqrt(h))
def num(s):
    try: return float(str(s).split(';')[0].split()[0])
    except: return None
rows=[]
class H(osmium.SimpleHandler):
    def way(self,w):
        hw=w.tags.get('highway')
        if hw not in DRIVE: return
        if w.tags.get('access') in ('no','private') and w.tags.get('motor_vehicle') not in ('yes','designated'): return
        if w.tags.get('motor_vehicle')=='no': return
        try: pts=[(n.ref,n.lon,n.lat) for n in w.nodes]
        except osmium.InvalidLocationError: return
        if not pts or not (x0<=pts[0][1]<=x1 and y0<=pts[0][2]<=y1): return
        base=hw.replace('_link',''); ms=num(w.tags.get('maxspeed'))
        sp=ms if ms and 5<=ms<=130 else SPEED[base]*(0.7 if hw.endswith('_link') else 1); ssrc='tag' if ms and 5<=ms<=130 else 'default'
        ln=num(w.tags.get('lanes')); lsrc='tag' if ln else 'default'
        ow=w.tags.get('oneway'); rev = ow=='-1'; one = ow in ('yes','true','1','-1') or w.tags.get('junction')=='roundabout' or hw in('motorway','motorway_link')
        tot=ln if ln else None
        nm='|'.join(x for x in [w.tags.get('name',''),w.tags.get('name:ru',''),w.tags.get('name:uz',''),w.tags.get('name:en',''),w.tags.get('ref','')] if x)
        for (a,ax,ay),(b,bx,by) in zip(pts[:-1],pts[1:]):
            L=hav(ay,ax,by,bx)
            if L==0: continue
            fwd=not rev
            def add(u,ux,uy,v,vx,vy,dirn):
                lanes = (tot if one else max(1,round(tot/2))) if tot else LANES[base]
                rows.append((u,v,ux,uy,vx,vy,L,hw,sp,ssrc,lanes,lsrc,w.id,nm))
            if one:
                if fwd: add(a,ax,ay,b,bx,by,1)
                else: add(b,bx,by,a,ax,ay,1)
            else: add(a,ax,ay,b,bx,by,1); add(b,bx,by,a,ax,ay,1)
H().apply_file(pbf,locations=True)
df=pd.DataFrame(rows,columns=['u','v','u_lon','u_lat','v_lon','v_lat','length_m','highway','speed_kmh','speed_src','lanes','lanes_src','way_id','names'])
df.to_csv(out,index=False,compression='gzip')
print(len(df),'directed edges;',df.u.nunique(),'nodes;',round(df.length_m.sum()/1000),'km')
print(df.highway.value_counts().to_string()); print(df.speed_src.value_counts().to_dict(), df.lanes_src.value_counts().to_dict())
