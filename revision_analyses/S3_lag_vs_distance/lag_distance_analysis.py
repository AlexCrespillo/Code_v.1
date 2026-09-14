"""
================================================================================
Propagation lag vs. along-channel distance  (AUDITED, tie-break included)
--------------------------------------------------------------------------------
Response analysis for Reviewer #4 comment 2 (physical support / flow routing)
and the reframing of the routing limitation (comment 16).

Question: does the drought propagation lag between connected upstream-origin
station pairs scale with the hydraulic (along-channel) distance separating them?
If it does not, lag is governed by catchment storage/memory rather than routing,
which justifies the deliberately non-routing design of the framework.

Inputs (edit PATHS):
  Locations.csv        station_id, xutm, yutm, NEAR_FID  (';' sep, decimal comma)
  Grafo.csv            river network node-arc list (arcid, from_node, to_node,
                       ORIG_FID, Shape_Length)  -> exported from
                       river_connectivity_singlepart
  Propagation_raw.csv   candidate upstream-origin pair-events (PRE tie-break),
                       with origin_event_id, upstream_event_id, overlap_days,
                       upstream_severity, lag_days, ...

Outputs: fig_lag_vs_distance.(png|pdf), pairwise_lag_distance.csv, printed stats.
================================================================================
"""
import pandas as pd, numpy as np, networkx as nx
from scipy import stats
import matplotlib.pyplot as plt

# ---------------------------- CONFIG ---------------------------------------
UP  = "./"
OUT = "./"
LOC_CSV, GRAPH_CSV, PAIRS_CSV = UP+"Locations.csv", UP+"Grafo.csv", UP+"Propagation_raw.csv"
REGULATED = [9101]     # stations downstream of major reservoirs (Yesa -> 9101)

# ------------------ 1. River network as undirected weighted graph ----------
g = pd.read_csv(GRAPH_CSV, sep=";")
g["length_m"] = g["Shape_Length"].astype(str).str.replace(",", ".").astype(float)
G = nx.Graph()
for r in g.itertuples():
    a, b = int(r.from_node), int(r.to_node)
    if G.has_edge(a, b):
        G[a][b]["weight"] = min(G[a][b]["weight"], r.length_m)
    else:
        G.add_edge(a, b, weight=r.length_m)
print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
      f"{nx.number_connected_components(G)} components, is_forest={nx.is_forest(G)}")
# is_forest=True  => no cycles => shortest path is the UNIQUE along-channel path.

# ------------------ 2. Station -> node (validated mapping) ------------------
loc = pd.read_csv(LOC_CSV, sep=";", decimal=",", encoding="utf-8-sig")
loc["NEAR_FID"] = loc["NEAR_FID"].astype(int)
seg = {int(o): (int(f), int(t)) for o, f, t in zip(g.ORIG_FID, g.from_node, g.to_node)}
node_of = {int(r.station_id): seg[int(r.NEAR_FID)][1] for r in loc.itertuples()}  # to_node
coord   = {int(r.station_id): (r.xutm, r.yutm) for r in loc.itertuples()}
def euclid_km(u, o):
    (x1,y1),(x2,y2) = coord[u], coord[o]; return np.hypot(x1-x2, y1-y2)/1000.0
def network_km(u, o):
    try: return nx.shortest_path_length(G, node_of[u], node_of[o], weight="weight")/1000.0
    except (nx.NetworkXNoPath, KeyError): return np.nan

# ------------------ 3. Candidate pairs -> apply deterministic tie-break -----
# Propagation_raw.csv is the PRE-tie-break candidate set: an origin event may have
# several candidate upstream EVENTS at the same upstream station. The published
# method (Sec. 2.1.3) retains ONE per station via the ranking:
#   (i) max overlap, (ii) max severity, (iii) min |onset difference|=min|lag|, (iv) event id.
# We reproduce that here so the analysis is consistent with the 525 chains.
p = pd.read_csv(PAIRS_CSV)
p["abs_lag"] = p["lag_days"].abs()
p = (p.sort_values(["origin_event_id", "upstream_station",
                    "overlap_days", "upstream_severity", "abs_lag", "upstream_event_id"],
                   ascending=[True, True, False, False, True, True])
       .drop_duplicates(["origin_event_id", "upstream_station"], keep="first")
       .reset_index(drop=True))
print(f"Matched pair-events after tie-break: {len(p)}  "
      f"(station-pairs: {p.groupby(['upstream_station','origin_station']).ngroups})")

# distances (one value per station pair, broadcast to events)
keys = p[["upstream_station","origin_station"]].drop_duplicates()
net = {(int(u),int(o)): network_km(int(u),int(o)) for u,o in zip(keys.upstream_station,keys.origin_station)}
euc = {(int(u),int(o)): euclid_km(int(u),int(o)) for u,o in zip(keys.upstream_station,keys.origin_station)}
p["net_km"] = [net[(int(u),int(o))] for u,o in zip(p.upstream_station,p.origin_station)]
p["euc_km"] = [euc[(int(u),int(o))] for u,o in zip(p.upstream_station,p.origin_station)]
assert p["net_km"].notna().all(), "unreachable pair(s) - check mapping/graph"

pair_tab = (p.groupby(["upstream_station","origin_station"])
              .agg(net_km=("net_km","first"), euc_km=("euc_km","first"),
                   med_abs_lag=("abs_lag","median"), n_events=("abs_lag","size")).reset_index())
print(f"Mapping check: corr(network, euclidean) = "
      f"{np.corrcoef(pair_tab.net_km, pair_tab.euc_km)[0,1]:.2f}; "
      f"median network/euclidean = {np.median(pair_tab.net_km/pair_tab.euc_km):.2f}")

# ------------------ 4. Statistics ------------------------------------------
def sp(x, y): r, pv = stats.spearmanr(x, y); return r, pv
re, pe = sp(p.net_km, p.abs_lag)
rp, pp = sp(pair_tab.net_km, pair_tab.med_abs_lag)
subp = pair_tab[~pair_tab.upstream_station.isin(REGULATED)]
rs, ps = sp(subp.net_km, subp.med_abs_lag)
p2 = p.copy()
p2["al_res"] = p2.abs_lag - p2.groupby("origin_station").abs_lag.transform("mean")
p2["nk_res"] = p2.net_km  - p2.groupby("origin_station").net_km.transform("mean")
rw, pw = sp(p2.nk_res, p2.al_res)
reg = p[p.upstream_station.isin(REGULATED)].abs_lag.mean()
oth = p[~p.upstream_station.isin(REGULATED)].abs_lag.mean()
print("\n--- |lag| vs along-channel distance ---")
print(f"event level          : n={len(p)}   rho={re:+.3f}  p={pe:.2e}")
print(f"station-pair level    : n={len(pair_tab)}   rho={rp:+.3f}  p={pp:.2f}")
print(f"pair, excl. {REGULATED} : n={len(subp)}   rho={rs:+.3f}  p={ps:.2f}")
print(f"within-origin partial : rho={rw:+.3f}  p={pw:.2e}")
print(f"regulation signature  : mean |lag| {REGULATED}={reg:.0f} d vs {oth:.0f} d ({reg/oth:.1f}x)")

# ------------------ 5. Figure ----------------------------------------------
plt.rcParams.update({"font.size":13,"axes.titlesize":13,"axes.labelsize":13,
    "xtick.labelsize":12,"ytick.labelsize":12,"font.family":"DejaVu Sans",
    "axes.spines.top":False,"axes.spines.right":False})
BLUE, ORANGE, PURPLE = "#2E75B6", "#E1812C", "#7030A0"
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
ax[0].scatter(p.net_km, p.abs_lag, s=9, alpha=0.20, color=BLUE, edgecolors="none")
bins = np.arange(0, p.net_km.max()+60, 60); mid, med = [], []
for i in range(len(bins)-1):
    m = (p.net_km >= bins[i]) & (p.net_km < bins[i+1])
    if m.sum() >= 10: mid.append(0.5*(bins[i]+bins[i+1])); med.append(p.abs_lag[m].median())
ax[0].plot(mid, med, "o-", color="k", lw=1.6, ms=5, label=r"binned median ($\geq$10 pairs)")
ax[0].set_xlabel("Along-channel distance (km)"); ax[0].set_ylabel("|lag| (days)")
ax[0].set_title(f"(a) Event level (n={len(p)})\nSpearman "+r"$\rho$"+f"={re:.2f}, p={pe:.1e}")
ax[0].legend(frameon=False, fontsize=11)
reg_m = pair_tab.upstream_station.isin(REGULATED)
ax[1].scatter(pair_tab.net_km[~reg_m], pair_tab.med_abs_lag[~reg_m], s=pair_tab.n_events[~reg_m]*3,
              alpha=0.65, color=ORANGE, edgecolors="k", linewidths=0.3, label="unregulated pairs")
ax[1].scatter(pair_tab.net_km[reg_m], pair_tab.med_abs_lag[reg_m], s=pair_tab.n_events[reg_m]*3,
              alpha=0.9, color=PURPLE, edgecolors="k", linewidths=0.5, label="regulated (Yesa, 9101)")
ax[1].set_xlabel("Along-channel distance (km)"); ax[1].set_ylabel("median |lag| (days)")
ax[1].set_title(f"(b) Station-pair level (n={len(pair_tab)})\nSpearman "+r"$\rho$"+f"={rp:.2f}, p={pp:.2f}")
ax[1].legend(frameon=False, fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT+"fig_lag_vs_distance.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT+"fig_lag_vs_distance.pdf", bbox_inches="tight")
pair_tab.sort_values("net_km").to_csv(OUT+"pairwise_lag_distance.csv", index=False)
print("\nSaved figure, pdf and pairwise table.")
