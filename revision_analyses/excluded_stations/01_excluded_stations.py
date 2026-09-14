"""
================================================================================
Excluded stations: metadata and reservoir-proximity test
--------------------------------------------------------------------------------
Response analysis for the Academic Editor, comment 4:
  "the final manuscript should contain a supplementary table listing the
  excluded stations and the reason for exclusion... This is important as the
  exclusion of a station may change the network topology and the
  reconstruction of propagation chains."

The 12 stations removed at the second screening stage (45 -> 33) were
identified by visual inspection of their SSI distributions and temporal
evolution, which showed abrupt regime shifts and inhomogeneous high/low
sequences inconsistent with a stationary index over 1961-2020.

This script does two things:
  1. Attaches official metadata (name, catchment area, elevation) to each
     excluded station, from the CEDEX Anuario de Aforos ESTAF table.
  2. Tests, rather than assumes, the working hypothesis that these
     irregularities are regulation-related, by comparing each station's
     distance to the nearest reservoir (EMBALSE table) between the excluded
     and retained sets.

RESULT: excluded stations are significantly closer to reservoirs than
retained ones (median 5.0 vs 12.6 km; Mann-Whitney U, p = 0.005), with 50%
within 5 km versus 18% of retained stations. This supports — but, being
observational and based on proximity alone, does not prove — a regulation
origin for the observed inhomogeneities.

Run:   python 01_excluded_stations.py
Needs: estaf.csv, embalse.csv (CEDEX Anuario de Aforos), Locations.csv
================================================================================
"""
import numpy as np
import pandas as pd
from scipy import stats

EXCLUDED = [9010, 9034, 9035, 9039, 9041, 9042, 9052, 9058, 9073, 9086, 9110, 9125]

est = pd.read_csv('estaf.csv', sep=';', encoding='latin-1')
emb = pd.read_csv('embalse.csv', sep=';', encoding='latin-1')
loc = pd.read_csv('Locations.csv', sep=';')
retained = loc.iloc[:, 0].astype(int).tolist()

res_x = emb['xutm30'].values.astype(float)
res_y = emb['yutm30'].values.astype(float)


def nearest_reservoir(station_id):
    """Euclidean distance (km) to the nearest reservoir, in UTM30 coordinates."""
    m = est[est.indroea == station_id]
    if len(m) == 0:
        return np.nan, None
    r = m.iloc[0]
    d = np.sqrt((res_x - float(r.xutm30)) ** 2 + (res_y - float(r.yutm30)) ** 2) / 1000.0
    i = int(np.nanargmin(d))
    return float(d[i]), emb.iloc[i].nom_embalse.strip()


# ============================================================================
# STEP 1 — Build the supplementary table
# ============================================================================
rows = []
for sid in EXCLUDED:
    r = est[est.indroea == sid].iloc[0]
    dist, res_name = nearest_reservoir(sid)
    rows.append({
        'Station code': sid,
        'Gauging station': r.lugar.strip().title(),
        'Catchment area (km2)': int(r.suprest),
        'Elevation (m a.s.l.)': int(r.alti),
        'Nearest reservoir': res_name.title(),
        'Distance to reservoir (km)': round(dist, 1),
        'Reservoir-proximate': 'Yes (<10 km)' if dist < 10 else 'No',
        'Exclusion stage': 'Stage 2 (SSI homogeneity screening)',
        'Reason for exclusion': ('Structural irregularity in the SSI distribution and '
                                  'temporal evolution (abrupt regime shifts, inhomogeneous '
                                  'high/low sequences)'),
    })

table = pd.DataFrame(rows).sort_values('Distance to reservoir (km)')
table.to_csv('table_S4_excluded_stations.csv', index=False)
print(table[['Station code', 'Gauging station', 'Catchment area (km2)',
             'Nearest reservoir', 'Distance to reservoir (km)']].to_string(index=False))

# ============================================================================
# STEP 2 — Test the regulation hypothesis
# ============================================================================
d_excl = np.array([nearest_reservoir(s)[0] for s in EXCLUDED])
d_kept = np.array([nearest_reservoir(s)[0] for s in retained])
d_excl = d_excl[~np.isnan(d_excl)]
d_kept = d_kept[~np.isnan(d_kept)]

u, p = stats.mannwhitneyu(d_excl, d_kept, alternative='less')

print(f"\n{'=' * 62}\nRESERVOIR-PROXIMITY TEST\n{'=' * 62}")
print(f"Excluded (n={len(d_excl)}): median {np.median(d_excl):.1f} km, mean {d_excl.mean():.1f} km")
print(f"Retained (n={len(d_kept)}): median {np.median(d_kept):.1f} km, mean {d_kept.mean():.1f} km")
print(f"\nMann-Whitney U (one-sided, excluded < retained): U={u:.1f}, p={p:.4f}")
print(f"Within 5 km of a reservoir: excluded {100*(d_excl<5).mean():.0f}%, "
      f"retained {100*(d_kept<5).mean():.0f}%")
print("\nInterpretation: excluded stations are significantly closer to reservoirs,")
print("supporting a regulation origin for the observed SSI inhomogeneities. This is")
print("proximity-based and observational: it does not establish that reservoir")
print("operation caused each individual irregularity.")
