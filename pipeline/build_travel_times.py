"""Free-flow mahalla-to-mahalla travel times (minutes) by Dijkstra on the road graph.
usage: python3 build_travel_times.py road_edges.csv.gz mahalla_zones.csv OUT_DIR
Centroids are snapped to the nearest node of the largest strongly connected component that touches a non-alley road (not service/living_street).
Intrazonal time = 0.5*sqrt(area)*1.35 / 25 km/h (as before). Free-flow only: congestion comes from the BPR loop."""
import sys, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.spatial import cKDTree
E, ZN, OUT = sys.argv[1:4]
e=pd.read_csv(E); z=pd.read_csv(ZN); n=len(z)
ids=pd.unique(np.r_[e.u.values,e.v.values]); ix=pd.Series(np.arange(len(ids)),index=ids)
u=ix[e.u].values; v=ix[e.v].values; t=e.length_m.values/(e.speed_kmh.values/3.6)   # seconds
G=sp.csr_matrix((t,(u,v)),shape=(len(ids),)*2)
# keep the min-time edge if parallel edges exist (csr sums duplicates), so rebuild with minimum
df=pd.DataFrame({'u':u,'v':v,'t':t}).groupby(['u','v']).t.min().reset_index()
G=sp.csr_matrix((df.t,(df.u,df.v)),shape=(len(ids),)*2)
nc,lab=connected_components(G,directed=True,connection='strong'); big=np.argmax(np.bincount(lab)); inbig=lab==big
print('nodes',len(ids),'largest SCC',int(inbig.sum()),f'({inbig.mean():.1%})')
main=e[~e.highway.isin(['service','living_street'])]
okn=np.zeros(len(ids),bool); okn[ix[main.u].values]=True; okn[ix[main.v].values]=True; okn&=inbig
xy=pd.concat([e[['u','u_lon','u_lat']].set_axis(['n','lon','lat'],axis=1),e[['v','v_lon','v_lat']].set_axis(['n','lon','lat'],axis=1)]).drop_duplicates('n').set_index('n').loc[ids]
cand=np.where(okn)[0]; k=np.cos(np.radians(41.3))
tree=cKDTree(np.c_[xy.lon.values[cand]*k,xy.lat.values[cand]])
d,j=tree.query(np.c_[z.centroid_lon.values*k,z.centroid_lat.values]); src=cand[j]
print('snap distance m: median',round(float(np.median(d)*111320)),'max',round(float(d.max()*111320)))
D=dijkstra(G,directed=True,indices=src)[:,src]/60.0          # minutes, origin x destination
intra=0.5*np.sqrt(z.area_km2.values)*1.35/25*60; D[np.diag_indices(n)]=intra
print('unreachable pairs',int(np.isinf(D).sum()))
np.save(f'{OUT}/travel_time_freeflow_min.npy',D.astype(np.float32))
off=~np.eye(n,dtype=bool)
import os
P=np.load(f'{OUT}/travel_time_min.npy') if os.path.exists(f'{OUT}/travel_time_min.npy') else None
print('mean time free-flow',round(float(D[off].mean()),1),'min', ('vs proxy '+str(round(float(P[off].mean()),1))) if P is not None else '')
print('percentiles of free-flow (min):',np.percentile(D[off],[5,25,50,75,95]).round(1))
