"""Congestion cascades between nearby junctions in one SUMO street window: which junctions start a jam and which ones just inherit it.
usage: python3 cascade.py NAME.edges.json OUT_PREFIX [--spill 0.8] [--min-share 0.10] [--jam 0.5]
NAME.edges.json comes from  micro/scene_sim.py ... --edge-stats NAME.edges.json  (halting cars and cars on every road, every 2 s).
Definitions
 * occupancy of a road = halting cars x 7.5 m / (length x lanes): 1 = the queue fills the whole road.
 * a road A->B is SPILLED while its occupancy >= --spill: the queue of junction B reaches back to junction A and blocks A's exit.
 * cascade graph: a directed edge B -> A (cause -> victim) when road A->B is spilled for at least --min-share of the recorded time.
 * root bottleneck: a junction that blocks others but is not blocked itself (none of its outgoing roads is spilled). Its score = halting car-seconds on the roads feeding all junctions it blocks, directly or through a chain.
 * jam share of a junction = share of time its fullest incoming road has occupancy >= --jam.
Writes OUT_PREFIX.json (numbers) and OUT_PREFIX.png (map)."""
import sys, json, argparse, numpy as np, networkx as nx
ap = argparse.ArgumentParser(); ap.add_argument('edges'); ap.add_argument('out'); ap.add_argument('--spill', type=float, default=0.8); ap.add_argument('--min-share', type=float, default=0.10); ap.add_argument('--jam', type=float, default=0.5)
ap.add_argument('--min-len', type=float, default=25); ap.add_argument('--no-plot', action='store_true'); a = ap.parse_args()
R = json.load(open(a.edges)); DT = R['dt']; E = {e['id']: e for e in R['edges']}; HALT = {k: np.array(v, float) for k, v in R['halt'].items()}; T = len(next(iter(HALT.values())))
cap = {k: max(E[k]['len'] * E[k]['lanes'] / 7.5, 1e-9) for k in E}; occ = {k: np.minimum(HALT[k] / cap[k], 1.0) for k in E}
tot_halt = sum(HALT[k].sum() for k in E) * DT
# ---- spill events and cascade graph
H = nx.DiGraph(); spill = {}
for k, e in E.items():
    if e['len'] < a.min_len or cap[k] < 3 or e['frm'] == e['to']: continue
    sh = float((occ[k] >= a.spill).mean()); spill[k] = sh
    if sh >= a.min_share:
        u, v = e['to'], e['frm']                                   # cause = the junction the road leads to, victim = the junction it starts at
        if H.has_edge(u, v): H[u][v]['share'] = max(H[u][v]['share'], sh)
        else: H.add_edge(u, v, share=sh, road=k, name=e['name'], len=e['len'])
# ---- node jam share (fullest incoming road)
inc = {}
for k, e in E.items():
    if e['len'] >= a.min_len: inc.setdefault(e['to'], []).append(k)
node_xy = {}
for e in E.values(): node_xy[e['frm']] = (e['fx'], e['fy']); node_xy[e['to']] = (e['tx'], e['ty'])
jam = {n: float((np.max([occ[k] for k in ks], axis=0) >= a.jam).mean()) for n, ks in inc.items()}
feed = {n: sum(HALT[k].sum() for k in ks) * DT for n, ks in inc.items()}      # halting car-seconds on the roads that feed the junction
# ---- roots and cascades
sccs = [c for c in nx.strongly_connected_components(H) if len(c) > 1]
C = nx.condensation(H); members = C.graph['mapping']; roots = []
for cn in C.nodes:
    if C.in_degree(cn) == 0 and C.out_degree(cn) > 0:
        nodes = [n for n, m in members.items() if m == cn]; desc = set(); [desc.update(nx.descendants(H, n)) for n in nodes]; blocked = desc | set(nodes)
        roots.append(dict(junctions=nodes, ids=[str(n) for n in nodes], x=float(np.mean([node_xy[n][0] for n in nodes])), y=float(np.mean([node_xy[n][1] for n in nodes])), blocked=len(blocked) - len(nodes), gridlock_loop=len(nodes) > 1,
                          score_car_s=float(sum(feed.get(n, 0) for n in blocked)), own_jam=float(max(jam.get(n, 0) for n in nodes))))
roots.sort(key=lambda r: -r['score_car_s'])
inCasc = set(H.nodes); inCascFeed = sum(feed.get(n, 0) for n in inCasc); allFeed = sum(feed.values())
top1 = roots[0]['score_car_s'] / max(allFeed, 1e-9) if roots else 0; top2 = sum(r['score_car_s'] for r in roots[:2]) / max(allFeed, 1e-9) if roots else 0
# ---- lagged correlation along cascade edges (does the victim's queue follow the cause's?)
lags = []
for u, v, d in H.edges(data=True):
    if u not in inc or v not in inc: continue
    x = np.max([occ[k] for k in inc[u]], axis=0); y = np.max([occ[k] for k in inc[v]], axis=0); best = (-2, 0)
    for lag in range(-30, 31):                                   # +lag: the victim follows the cause by lag*DT seconds
        xs, ys = (x[:T-lag], y[lag:]) if lag >= 0 else (x[-lag:], y[:T+lag])
        if len(xs) > 20 and xs.std() > 1e-6 and ys.std() > 1e-6:
            c = float(np.corrcoef(xs, ys)[0, 1]);  best = max(best, (c, lag))
    lags.append(dict(cause=u, victim=v, road=d['name'], corr=round(best[0], 2), lag_s=best[1] * DT))
follow = float(np.mean([l['lag_s'] > 0 for l in lags if l['corr'] > 0.3])) if any(l['corr'] > 0.3 for l in lags) else None
out = dict(scene=R['id'], record_s=R['record'], roads=len(E), junctions_with_signal_or_node=len(inc), teleports=R.get('teleports'), spilled_roads=sum(1 for v in spill.values() if v >= a.min_share),
           cascade_junctions=len(H.nodes), cascade_links=len(H.edges), gridlock_loops=[len(c) for c in sccs], share_of_queueing_in_cascade_junctions=round(inCascFeed / max(allFeed, 1e-9), 2),
           roots=[dict(ids=r['ids'], x=round(r['x']), y=round(r['y']), blocked_junctions=r['blocked'], gridlock_loop=r['gridlock_loop'], share_of_all_queueing=round(r['score_car_s'] / max(allFeed, 1e-9), 2), own_jam_share=round(r['own_jam'], 2)) for r in roots[:6]],
           top1_root_share=round(top1, 2), top2_roots_share=round(top2, 2), victims_follow_cause_share=None if follow is None else round(follow, 2), lags=lags[:40])
json.dump(out, open(a.out + '.json', 'w'), indent=1)
print(f"{R['id']}: {len(E)} roads, {len(H.nodes)} junctions in cascades ({len(H.edges)} spill links), {len(roots)} root bottlenecks; top root explains {100*top1:.0f}% of all queueing (top 2: {100*top2:.0f}%); cascade junctions hold {100*out['share_of_queueing_in_cascade_junctions']:.0f}% of it;"
      f" gridlock loops {out['gridlock_loops']}; victim follows cause in {out['victims_follow_cause_share']} of correlated links")
for r in out['roots'][:4]: print('   root at', (r['x'], r['y']), 'm: blocks', r['blocked_junctions'], 'junctions, share', r['share_of_all_queueing'], 'loop' if r['gridlock_loop'] else '')
# ---- map
if not a.no_plot:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    fig, ax = plt.subplots(figsize=(10, 10), dpi=130); fig.patch.set_facecolor('#0a1017'); ax.set_facecolor('#0a1017')
    for k, e in E.items(): ax.plot([e['fx'], e['tx']], [e['fy'], e['ty']], color=plt.cm.YlOrRd(min(1, float(occ[k].mean()) * 1.6)) if occ[k].mean() > 0.05 else '#1f2f40', lw=0.6 + 2.2 * min(1, float(occ[k].mean())), solid_capstyle='round', zorder=1)
    for u, v, d in H.edges(data=True):
        (x0, y0), (x1, y1) = node_xy[u], node_xy[v]; ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=14, lw=1.4 + 5 * d['share'], color='#3fd0b0', alpha=0.9, zorder=3, shrinkA=4, shrinkB=4))
    for n, (x, y) in node_xy.items():
        if n in jam and jam[n] > 0.05: ax.scatter([x], [y], s=20 + 260 * jam[n], facecolor='none', edgecolor='#ff6a5c', lw=1.2, zorder=4)
    for r in roots[:3]: ax.scatter([r['x']], [r['y']], s=260, marker='*', color='#f2b64c', edgecolor='#0a1017', zorder=6)
    ax.set_xlim(-800, 800); ax.set_ylim(-800, 800); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f"{R['id']}: roads coloured by average queue · teal arrows = cause → victim (queue spills back) · red rings = jam share · star = root bottleneck", color='#d9e4ee', fontsize=9)
    fig.savefig(a.out + '.png', facecolor=fig.get_facecolor(), bbox_inches='tight'); print('wrote', a.out + '.png')
