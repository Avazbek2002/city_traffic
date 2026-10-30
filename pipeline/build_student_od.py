"""Student gravity models (school, college/lyceum, university). Doubly-constrained, IPF, exp(-beta*T), person trips per morning (arrive ~08:00).
usage: python3 build_student_od.py mahalla_zones.csv mahalla_student_production.csv mahalla_zones_places.csv mahalla_university_students.csv OUTDIR [travel_times.npy]
Everything marked ASSUMPTION is a placeholder; override with env vars (e.g. TARGET_SCHOOL=8)."""
import sys, os, numpy as np, pandas as pd
ZN, SP, PL, UN, OUT = sys.argv[1:6]; TT = sys.argv[6] if len(sys.argv) > 6 else None
os.makedirs(OUT, exist_ok=True)
E = lambda k, d: float(os.environ.get(k, d))
# ---- ASSUMPTIONS ----
TARGET = {'school': E('TARGET_SCHOOL', 10), 'college': E('TARGET_COLLEGE', 20), 'university': E('TARGET_UNIVERSITY', 30)}  # mean trip time, min
PUPILS_D = {2603:64411,2602:58107,2609:57335,2610:54856,2606:50799,2605:46512,2608:46192,2607:36751,2601:30043,2612:22845,2604:19579,2611:11987}  # stat.uz, 1 Jan 2023
PUPIL_GROWTH = E('PUPIL_GROWTH', 1.0)      # scale 2023 pupils to the model year
SHARE_1617_SCHOOL = E('SHARE_1617_SCHOOL', 0.5)   # of ages 16-17, share in grades 10-11 (rest: college/lyceum)
COLLEGE_1819 = E('COLLEGE_1819', 0.15)     # share of ages 18-19 in college/lyceum
UNIV_TOTAL = E('UNIV_TOTAL', 354000)       # university students making a daily trip (stat.uz 2022/23 = 354k; OTM sheet suggests up to ~783k)
UNIV_RES_SHARE = E('UNIV_RES_SHARE', 179675/783465)  # share of students with permanent residence in Tashkent (OTM sheet)
NONRES_DORM = E('NONRES_DORM', 0.5)        # non-resident students living near their campus (dorm/rent) vs spread like residents
MODE = {  # walk / public transport / car share (user: students almost all use public transport)
  'school': dict(walk=0.10, pt=0.80, car=0.10), 'college': dict(walk=0.07, pt=0.91, car=0.02), 'university': dict(walk=0.05, pt=0.93, car=0.02)}
CAR_OCC = E('CAR_OCC', 1.3)
CIRC, SPEED, INTRA_F = 1.35, 25.0, 0.5     # proxy travel time (REPLACE with Dijkstra/A* times)

z = pd.read_csv(ZN); sp = pd.read_csv(SP); n = len(z)
assert (z.mahalla_id.values == sp.mahalla_id.values).all()
lat = np.radians(z.centroid_lat.values); lon = np.radians(z.centroid_lon.values)
if TT: T = np.load(TT); src = 'supplied'
else:
    a = np.sin((lat[:,None]-lat[None])/2)**2 + np.cos(lat[:,None])*np.cos(lat[None])*np.sin((lon[:,None]-lon[None])/2)**2
    T = 2*6371.0*np.arcsin(np.sqrt(a))*CIRC/SPEED*60; T[np.diag_indices(n)] = INTRA_F*np.sqrt(z.area_km2.values)*CIRC/SPEED*60; src = 'PROXY'
idx = {int(m): i for i, m in enumerate(z.mahalla_id)}; did = z.district_id.values
pl = (pd.read_csv(PL) if str(PL).endswith('.csv') else pd.read_pickle(PL)); pl = pl[pl.mahalla_id.notna()]
def count(sub):
    v = np.zeros(n)
    for m in pl[pl.subtype == sub].mahalla_id.astype(int): v[idx[m]] += 1
    return v

# ---- production (persons) and attraction per level ----
p67, p815, p1617, p1819, p2024 = (sp[c].values for c in ['age6_7','age8_15','age16_17','age18_19','age20_24'])
# SCHOOL: attraction = district pupils spread over that district's schools (2GIS school count)
sch = count('school'); A_s = np.zeros(n)
for d, pup in PUPILS_D.items():
    m = did == d
    A_s[m] = pup*PUPIL_GROWTH*sch[m]/sch[m].sum()
P_s = p67 + p815 + SHARE_1617_SCHOOL*p1617
# COLLEGE / LYCEUM: attraction = 2GIS college listings (only ~15: LOW RELIABILITY)
P_c = (1-SHARE_1617_SCHOOL)*p1617 + COLLEGE_1819*p1819; A_c = count('college')
# UNIVERSITY: attraction = students by campus (OTM sheet x 2GIS listings); production = residents + non-residents
un = pd.read_csv(UN); A_u = np.zeros(n)
for m, s in zip(un.mahalla_id.astype(int), un.univ_students): A_u[idx[m]] += s
res_w = p1819 + p2024
P_u = UNIV_RES_SHARE*res_w/res_w.sum() + (1-UNIV_RES_SHARE)*((1-NONRES_DORM)*res_w/res_w.sum() + NONRES_DORM*A_u/A_u.sum())
levels = {'school': (P_s, A_s, A_s.sum()), 'college': (P_c, A_c, P_c.sum()), 'university': (P_u, A_u, UNIV_TOTAL)}

def ipf(K, P, A, iters=3000, tol=1e-10):
    r = np.ones(len(P)); c = np.ones(len(A))
    for it in range(iters):
        r = np.divide(P, K@c, out=np.zeros_like(P), where=(K@c) > 0); c = np.divide(A, K.T@r, out=np.zeros_like(A), where=(K.T@r) > 0)
        M = r[:,None]*K*c[None]; err = max(np.abs(M.sum(1)-P).max()/P.max(), np.abs(M.sum(0)-A).max()/A.max())
        if err < tol: break
    return M, it, err
mt = lambda M: float((M*T).sum()/M.sum())
res, rows, mode_rows = {}, [], []
for k, (P0, A0, tot) in levels.items():
    P = P0/P0.sum()*tot; A = A0/A0.sum()*tot; lo, hi = 1e-4, 1.0
    for _ in range(40):
        b = np.sqrt(lo*hi); M, it, err = ipf(np.exp(-b*T), P, A)
        if mt(M) > TARGET[k]: lo = b
        else: hi = b
    b = np.sqrt(lo*hi); M, it, err = ipf(np.exp(-b*T), P, A); res[k] = M
    rows.append(dict(level=k, persons=round(M.sum()), beta=round(b,4), mean_min=round(mt(M),1), target=TARGET[k], intrazonal_share=round(float(np.trace(M)/M.sum()),3),
                     origin_zones_zero_attr=int((A==0).sum()), ipf_err=float(f'{err:.0e}')))
    for mode, s in MODE[k].items(): res[f'{k}_{mode}'] = M*s
    mode_rows.append(dict(level=k, pt_persons=round(M.sum()*MODE[k]['pt']), car_vehicles=round(M.sum()*MODE[k]['car']/CAR_OCC), walk_persons=round(M.sum()*MODE[k]['walk'])))
np.savez_compressed(f'{OUT}/od_student_matrices.npz', mahalla_id=z.mahalla_id.values, **res)
s = pd.DataFrame(rows).merge(pd.DataFrame(mode_rows)); s.to_csv(f'{OUT}/od_student_summary.csv', index=False); print('travel times:', src); print(s.to_string(index=False))
# district checks: school pupils attracted vs resident school-age
dm = lambda M: pd.DataFrame(M, index=did, columns=did).groupby(level=0).sum().T.groupby(level=0).sum().T
d = dm(res['school']); print('\nschool: share of pupils attending school in own district:', round(float(np.trace(d.values)/d.values.sum()),3))
chk = pd.DataFrame({'pupils_in_district_schools': d.sum(0), 'pupils_living_in_district': d.sum(1)}).round(0); chk['net_inflow'] = chk.iloc[:,0]-chk.iloc[:,1]; print(chk.astype(int).to_string())
for k in ['college','university']:
    d = dm(res[k]); print(k,'own-district share:', round(float(np.trace(d.values)/d.values.sum()),3))
d.round(0).astype(int).to_csv(f'{OUT}/od_university_district_matrix.csv')
dm(res['school']).round(0).astype(int).to_csv(f'{OUT}/od_school_district_matrix.csv')

np.savez(f'{OUT}/pa_vectors.npz', mahalla_id=z.mahalla_id.values, **{f'{k}_{w}': v for k,(P0,A0,tot) in levels.items() for w,v in (('P',P0/P0.sum()*tot),('A',A0/A0.sum()*tot))})
