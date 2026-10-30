"""Doubly-constrained gravity model on mahalla zones, one matrix per purpose.

Usage: python3 build_od.py mahalla_zones.csv OUTDIR [travel_times.npy]
 - travel_times.npy (optional): 585x585 minutes, rows/cols in mahalla_zones.csv order. If omitted, a PROXY is used:
   haversine centroid distance * CIRCUITY / SPEED_KMH (intrazonal from zone area).
All parameters marked ASSUMPTION are placeholders to calibrate against the infographic / real data.
"""
import sys, os, json, numpy as np, pandas as pd
ZONES, OUT = sys.argv[1:3]
TT = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(OUT, exist_ok=True)

CIRCUITY  = 1.35   # ASSUMPTION: road distance / straight-line distance
SPEED_KMH = 25.0   # ASSUMPTION: average urban speed (peak)
INTRA_F   = 0.5    # ASSUMPTION: intrazonal distance = INTRA_F * sqrt(area) * CIRCUITY
# ASSUMPTION: target mean trip time (minutes) per purpose - sets the distance-decay beta
TARGET_MIN = {'work': 22.0, 'education': 12.0, 'health': 17.0, 'retail': 10.0}
for _k in TARGET_MIN:   # optional override, e.g. TARGET_WORK=20
    if os.environ.get('TARGET_'+_k.upper()): TARGET_MIN[_k] = float(os.environ['TARGET_'+_k.upper()])

z = pd.read_csv(ZONES)
n = len(z)
lat = np.radians(z.centroid_lat.values); lon = np.radians(z.centroid_lon.values)

# ---------------------------------------------------------------- impedance (minutes)
if TT:
    T = np.load(TT); assert T.shape == (n, n)
    src = 'supplied travel times'
else:
    dlat = lat[:, None] - lat[None, :]; dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat/2)**2 + np.cos(lat[:, None])*np.cos(lat[None, :])*np.sin(dlon/2)**2
    D = 2*6371.0*np.arcsin(np.sqrt(a))                       # km, straight line
    T = D*CIRCUITY/SPEED_KMH*60.0
    intra = INTRA_F*np.sqrt(z.area_km2.values)*CIRCUITY/SPEED_KMH*60.0
    T[np.diag_indices(n)] = intra
    src = f'PROXY (haversine x {CIRCUITY} / {SPEED_KMH} km/h)'
np.save(f'{OUT}/travel_time_min.npy', T.astype(np.float32))

# ---------------------------------------------------------------- production / attraction per purpose
pop = z.population.values.astype(float)
purposes = {
  'work':      dict(P=pop, A=z.employment.values.astype(float)),
  'education': dict(P=pop, A=z.education_weighted.values.astype(float)),
  'health':    dict(P=pop, A=z.health_weighted.values.astype(float)),
  'retail':    dict(P=pop, A=z.retail_weighted.values.astype(float)),
}
TOTAL = {'work': float(z.employment.sum())}      # work: one commuter per job (persons, one-way, morning)
for k in purposes:                                 # other purposes: shares (sum = 1) until volumes are calibrated
    TOTAL.setdefault(k, 1.0)

def ipf(K, P, A, iters=2000, tol=1e-10):
    """Scale K so row sums = P and col sums = A (P and A must have equal totals)."""
    r = np.ones(len(P)); c = np.ones(len(A))
    for it in range(iters):
        r = np.divide(P, K@c, out=np.zeros_like(P), where=(K@c) > 0)
        c = np.divide(A, K.T@r, out=np.zeros_like(A), where=(K.T@r) > 0)
        M = (r[:, None]*K)*c[None, :]
        err = max(np.abs(M.sum(1)-P).max()/max(P.max(), 1e-12), np.abs(M.sum(0)-A).max()/max(A.max(), 1e-12))
        if err < tol: break
    return M, it, err

def mean_time(M): return float((M*T).sum()/M.sum())

def fit(purpose):
    P0, A0 = purposes[purpose]['P'], purposes[purpose]['A']
    tot = TOTAL[purpose]
    P = P0/P0.sum()*tot; A = A0/A0.sum()*tot
    lo, hi = 1e-4, 1.0                                # beta bracket (1/min)
    for _ in range(40):                               # bisection on log beta
        b = np.sqrt(lo*hi); M, it, err = ipf(np.exp(-b*T), P, A)
        if mean_time(M) > TARGET_MIN[purpose]: lo = b
        else: hi = b
    b = np.sqrt(lo*hi); M, it, err = ipf(np.exp(-b*T), P, A)
    return M, b, it, err, P, A

res = {}; summary = []
for p in purposes:
    M, b, it, err, P, A = fit(p)
    res[p] = M
    summary.append(dict(purpose=p, beta_per_min=round(b, 4), mean_time_min=round(mean_time(M), 2), target_min=TARGET_MIN[p],
                        total=round(M.sum(), 3), ipf_iters=it, ipf_rel_err=float(f'{err:.1e}'),
                        intrazonal_share=round(float(np.trace(M)/M.sum()), 4),
                        row_err=float(np.abs(M.sum(1)-P).max()), col_err=float(np.abs(M.sum(0)-A).max()),
                        zones_zero_attraction=int((A == 0).sum())))
np.savez_compressed(f'{OUT}/od_matrices.npz', mahalla_id=z.mahalla_id.values, **res)
pd.DataFrame(summary).to_csv(f'{OUT}/od_fit_summary.csv', index=False)
print('travel times:', src); print(pd.DataFrame(summary).to_string(index=False))

# ---------------------------------------------------------------- long-form top flows and district summary
name = z.mahalla_name.values; did = z.district_id.values; mid = z.mahalla_id.values
rows = []
for p, M in res.items():
    flat = np.argsort(M, axis=None)[::-1][:300]
    for f in flat:
        i, j = divmod(int(f), n)
        rows.append(dict(purpose=p, origin_id=mid[i], origin=name[i], origin_district=did[i], dest_id=mid[j], dest=name[j],
                         dest_district=did[j], flow=round(float(M[i, j]), 5 if p != 'work' else 1), time_min=round(float(T[i, j]), 1)))
pd.DataFrame(rows).to_csv(f'{OUT}/od_top_flows.csv', index=False)
dd = pd.DataFrame(res['work'], index=did, columns=did).groupby(level=0).sum().T.groupby(level=0).sum().T
dd.round(0).astype(int).to_csv(f'{OUT}/od_work_district_matrix.csv')
print('\nwork trips district->district (rows = origin):'); print(dd.round(0).astype(int).to_string())
print('within-district share of work trips:', round(float(np.trace(dd.values)/dd.values.sum()), 3))
