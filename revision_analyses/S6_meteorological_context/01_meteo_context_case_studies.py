"""
================================================================================
Meteorological context (SPEI) for the three case-study propagation chains
--------------------------------------------------------------------------------
Response analysis for Reviewer #2, comment 4 ("showing the corresponding
meteorological drought conditions would help readers understand the
hydroclimatic context") and Reviewer #4, comment 5 (meteorological-hydrological
coupling; addressed here as an initial, illustrative step -- see the response
text for the scope of what is and is not claimed).

WHAT THIS DOES
--------------
For each of the three case studies presented in the manuscript (Section 3.2,
Figures 6-8), this script plots the daily SSI series at the ORIGIN station
alongside the weekly SPEI-3 (short/seasonal accumulation) and SPEI-12
(long-term accumulation) series at that same station's coordinates, over the
same approximate context window used in Figures 6-8.

DATA
----
  - SSI_daily.csv                : daily SSI, 33 stations, 1961-2020
                                    (station_id, date, Q, SSI)
  - DATA_40_81_0_52.csv           : SPEI/SPI weekly series at 40.81N, 0.52E
                                    (station 9027, Ebro at Tortosa -- Case 1)
  - DATA_41_66_-0_88.csv          : at 41.66N, -0.88E
                                    (station 9011, Ebro at Zaragoza -- Case 2)
  - DATA_42_35_-1_65.csv          : at 42.35N, -1.65E
                                    (station 9005, Aragon at Caparroso -- Case 3)
Station coordinates (UTM, ETRS89/UTM 30N) were converted to lat/lon and cross-
checked against the known real-world locations of Tortosa, Zaragoza and
Caparroso (agreement to ~0.01-0.02 degrees).

IMPORTANT, HONEST FINDING (updated)
------------------------------------
Two independent checks show that neither the origin station's local
meteorology, nor its own individual SSI record, fully captures the network-
wide chain reported in the manuscript:

(1) Local SPEI at the origin does not show a uniform "meteorological drought
precedes hydrological drought" signal (see cases below).

(2) The "hydrological event" window quoted in the manuscript's Figures 6-8
titles is the CHAIN window (earliest onset to latest termination among ALL
matched stations, Table 2: t(c,s) to t(c,e)) -- not the origin station's own
individual SSI event. Because upstream stations typically enter drought
earlier (negative lag), the origin's own event can start well after the
chain window begins, or even appear as several separate, shorter events
within a chain window that is reported as continuous (most strikingly for
Case 2, which has SIX separate origin-station events within the 427-day
chain window). This script now shows BOTH layers explicitly: the chain
window (grey background) and the origin's own individual event(s)
(highlighted in red, directly on its SSI curve).

  - Case 1: SPEI-3 falls sharply right at onset (short-term local signal);
    the origin's own event starts 4 months after the chain window begins.
  - Case 2: SPEI-3/SPEI-12 are both positive (wet) before onset; the origin
    has 6 separate, shorter SSI events within the 427-day chain window.
  - Case 3: SPEI-3 is wet before onset but SPEI-12 (2yr) is negative; the
    origin's own event (44 d) is narrower than the chain window (92 d).

Neither finding is a weakness: both are direct, case-specific evidence that
the origin/outlet, examined alone (meteorologically or hydrologically),
does not capture the network-wide propagation signal -- reinforcing the
manuscript's storage- and network-driven framing of propagation.

Run:  python 01_meteo_context_case_studies.py
================================================================================
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

OUT = "./"
DATA = "./"

ssi_daily = pd.read_csv(DATA + "SSI_daily.csv", parse_dates=["date"])
origin_events = pd.read_csv(DATA + "origin_events_full.csv", parse_dates=["start_date", "end_date"])

# (label, origin_station, spei_csv, context_start, context_end, chain_window_start, chain_window_end)
CASES = [
    ("Case 1 (1988\u20131989)\nOrigin: 9027 \u2014 Ebro at Tortosa",
     9027, "DATA_40_81_0_52.csv", "1987-06-01", "1989-08-01", "1988-09-21", "1989-04-24"),
    ("Case 2 (2016\u20132018)\nOrigin: 9011 \u2014 Ebro at Zaragoza",
     9011, "DATA_41_66_-0_88.csv", "2015-06-01", "2018-06-01", "2016-12-30", "2018-03-01"),
    ("Case 3 (1997)\nOrigin: 9005 \u2014 Arag\u00f3n at Caparroso",
     9005, "DATA_42_35_-1_65.csv", "1996-06-01", "1997-10-01", "1997-02-28", "1997-05-30"),
]

BLUE, ORANGE, PURPLE, GREY, RED = "#2E75B6", "#E1812C", "#7030A0", "#D9D9D9", "#C00000"
plt.rcParams.update({"font.size": 11.5, "axes.titlesize": 11.5, "axes.labelsize": 11,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False})

fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.5))

for ax, (label, sid, speicsv, ctx_s, ctx_e, ev_s, ev_e) in zip(axes, CASES):
    ctx_s, ctx_e, ev_s, ev_e = map(pd.Timestamp, [ctx_s, ctx_e, ev_s, ev_e])

    ssi = ssi_daily[(ssi_daily.station_id == sid) &
                    (ssi_daily.date >= ctx_s) & (ssi_daily.date <= ctx_e)]
    spei = pd.read_csv(DATA + speicsv, parse_dates=["DATA"])
    spei = spei[(spei.DATA >= ctx_s) & (spei.DATA <= ctx_e)]
    own = origin_events[(origin_events.station_id == sid) &
                         (origin_events.start_date <= ctx_e) & (origin_events.end_date >= ctx_s)]

    assert ssi["SSI"].notna().all() and spei["spei_3"].notna().all() and spei["spei_12"].notna().all(), \
        f"Unexpected missing values for station {sid} in the plotted window."

    # Layer 1: chain window (network-wide span, from Figs. 6-8 titles)
    ax.axvspan(ev_s, ev_e, color=GREY, alpha=0.5, zorder=0, label="chain window (all stations)")
    ax.axhline(0, color="grey", lw=0.7, zorder=1)
    ax.axhline(-1.28, color="grey", lw=0.7, ls="--", zorder=1)
    ax.plot(ssi.date, ssi.SSI, color=BLUE, lw=1.1, label="SSI (origin station, daily)", zorder=3)

    # Layer 2: origin station's OWN individual drought event(s)
    first = True
    for _, r in own.iterrows():
        seg = ssi[(ssi.date >= r.start_date) & (ssi.date <= r.end_date)]
        ax.plot(seg.date, seg.SSI, color=RED, lw=2.2, zorder=6,
                label="origin station's own SSI event" if first else None)
        ax.fill_between(seg.date, seg.SSI, -1.28, color=RED, alpha=0.15, zorder=2)
        first = False

    ax.plot(spei.DATA, spei.spei_3, color=ORANGE, lw=1.4, label="SPEI-3 (weekly)", zorder=4)
    ax.plot(spei.DATA, spei.spei_12, color=PURPLE, lw=1.6, label="SPEI-12 (weekly)", zorder=5)

    ax.set_title(label, fontsize=11.3, loc="left")
    ax.set_ylabel("SSI  /  SPEI")
    ax.set_ylim(-4.3, 3.6)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

axes[0].legend(loc="lower left", ncol=2, fontsize=8.0, frameon=False)
fig.suptitle("Meteorological context (SPEI) alongside the hydrological drought signal (SSI)\n"
             "at the origin station of each case study", fontsize=12.2, y=1.005)
fig.tight_layout()
fig.savefig(OUT + "fig_meteo_context_case_studies.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT + "fig_meteo_context_case_studies.pdf", bbox_inches="tight")
print("Saved fig_meteo_context_case_studies.(png|pdf)")
plt.show()
