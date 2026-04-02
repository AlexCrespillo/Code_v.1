"""
Publication-quality figure – Drought propagation chain metric distributions.
2×3 grid of KDE density plots, one per chain descriptor.

Assumes chains_df already exists, or loads from CSV below.
Output: PNG (300 dpi) + PDF to the project directory.
"""

import os
import matplotlib
matplotlib.use("Agg")          # raster backend for dual PNG+PDF export
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

# ─────────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────────
CSV_PATH = os.path.join("data", "Chains_final.csv")
chains_df = pd.read_csv(CSV_PATH)

# ─────────────────────────────────────────────────────────────────────────────
# Variable definitions
#   col     : column name in chains_df
#   label   : x-axis label
#   log10   : use log10 x-axis (right-skewed variables)
#   abs_val : take absolute value before plotting (lag_mean ≤ 0 by construction)
# ─────────────────────────────────────────────────────────────────────────────
VARS = [
    dict(col="chain_duration",       label="Duration (days)",               log10=False, abs_val=False),
    dict(col="chain_size",           label="Chain size (stations)",          log10=False, abs_val=False),
    dict(col="propagation_fraction", label="Propagation fraction",           log10=False, abs_val=False),
    dict(col="lag_mean",             label="Lag (days)",                     log10=False, abs_val=True),
    dict(col="chain_severity",       label="Severity (–)",                   log10=True,  abs_val=False),
    dict(col="chain_severity_hm3",   label=r"Volumetric severity (hm$^3$)",  log10=True,  abs_val=False),
]

PANEL_LABELS = list("abcdef")

# ─────────────────────────────────────────────────────────────────────────────
# Colour palette
# ─────────────────────────────────────────────────────────────────────────────
C_KDE   = "#1a4f72"   # deep navy       – KDE line
C_FILL  = "#c5dff3"   # very light blue – KDE fill
C_ANNOT = "#1b2631"   # near-black navy – annotation text

# ─────────────────────────────────────────────────────────────────────────────
# Global rcParams
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.labelsize":   8.0,
    "ytick.labelsize":   8.0,
    "pdf.fonttype":      42,     # TrueType fonts in PDF (editable)
    "ps.fonttype":       42,
})

FS_XLABEL = 8.5
FS_YLABEL = 8.0
FS_PANEL  = 9.5
FS_STATS  = 6.2    # stats box text
LINE_W    = 1.5

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    2, 3,
    figsize=(7.0, 4.4),      # two-column journal width ≈ 7 in
    constrained_layout=True,
)
axes_flat = axes.flatten()

for ax, vdef, plabel in zip(axes_flat, VARS, PANEL_LABELS):

    # ── Extract and clean data ────────────────────────────────────────────────
    raw = chains_df[vdef["col"]].dropna().values.astype(float)

    if vdef["abs_val"]:
        raw = np.abs(raw)

    if vdef["log10"]:
        raw = raw[raw > 0]          # drop zeros/negatives before log transform
        data = np.log10(raw)
    else:
        data = raw

    if len(data) < 5:
        ax.set_visible(False)
        continue

    # ── KDE ───────────────────────────────────────────────────────────────────
    kde = gaussian_kde(data, bw_method="scott")
    pad = (data.max() - data.min()) * 0.08
    xs  = np.linspace(data.min() - pad, data.max() + pad, 500)
    ys  = kde(xs)

    ax.fill_between(xs, ys, alpha=0.55, color=C_FILL, linewidth=0)
    ax.plot(xs, ys, color=C_KDE, linewidth=LINE_W)

    # ── Summary stats box (top-right) ─────────────────────────────────────────
    # Compute stats on original (non-log) scale for readability
    orig = raw   # raw already has abs_val applied; log10 not yet applied
    q1, q3 = np.percentile(orig, [25, 75])

    def fmt(v):
        """Format a value compactly without scientific notation."""
        av = abs(v)
        if av >= 10000:
            return f"{v:,.0f}"
        elif av >= 1000:
            return f"{v:,.0f}"
        elif av >= 100:
            return f"{v:.1f}"
        elif av >= 10:
            return f"{v:.2f}"
        else:
            return f"{v:.3f}"

    stats_text = (
        f"Mean   {fmt(orig.mean())}\n"
        f"Median {fmt(np.median(orig))}\n"
        f"Min    {fmt(orig.min())}\n"
        f"Max    {fmt(orig.max())}\n"
        f"IQR    {fmt(q3 - q1)}"
    )
    ax.text(
        0.97, 0.97, stats_text,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=FS_STATS,
        fontfamily="monospace",
        color="black",
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white",
            edgecolor="black",
            linewidth=0.7,
        ),
        zorder=5,
    )

    # ── x-axis ticks (log10: label as 10^n) ──────────────────────────────────
    if vdef["log10"]:
        ticks = np.arange(np.floor(data.min() - pad), np.ceil(data.max() + pad) + 1)
        ax.set_xticks(ticks)
        ax.set_xticklabels(
            [f"$10^{{{int(t)}}}$" for t in ticks],
            fontsize=7.5,
        )

    # ── Axis formatting ───────────────────────────────────────────────────────
    ax.set_xlabel(vdef["label"], fontsize=FS_XLABEL, color=C_ANNOT, labelpad=3)
    ax.set_ylabel("Density", fontsize=FS_YLABEL, color=C_ANNOT, labelpad=3)
    # Show 3 evenly-spaced tick values on y-axis (no scientific notation).
    # y_max is captured as default arg to avoid the closure-in-loop pitfall:
    # without it, all formatters would use the last panel's y_max at render time.
    y_max = ys.max()
    yticks = np.linspace(0, y_max, 4)[1:]   # skip 0; show ~3 ticks
    ax.set_yticks(yticks)
    def _yfmt(v, _, _ymax=y_max):
        if _ymax < 0.01:
            return f"{v:.4f}"
        elif _ymax < 0.1:
            return f"{v:.3f}"
        elif _ymax < 1:
            return f"{v:.2f}"
        else:
            return f"{v:.1f}"
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(_yfmt))
    ax.set_ylim(bottom=0)
    ax.set_xlim(xs[0], xs[-1])
    ax.tick_params(colors=C_ANNOT)
    for spine in ax.spines.values():
        spine.set_edgecolor(C_ANNOT)

    # ── Panel label (a)–(f) ───────────────────────────────────────────────────
    ax.text(
        -0.15, 1.02, f"({plabel})",
        transform=ax.transAxes,
        fontsize=FS_PANEL, fontweight="bold",
        color=C_ANNOT, va="bottom", ha="left",
    )

# ─────────────────────────────────────────────────────────────────────────────
# Export
# ─────────────────────────────────────────────────────────────────────────────
OUT_BASE = os.path.join("output", "figures", "chain_distributions")

fig.savefig(OUT_BASE + ".png", format="png", dpi=300, bbox_inches="tight")

# Re-save as PDF using the pdf backend
matplotlib.use("pdf")
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(OUT_BASE + ".pdf") as pdf:
    pdf.savefig(fig, bbox_inches="tight",
                metadata={"Title": "Drought chain metric distributions"})

plt.close(fig)
print(f"Saved:\n  {OUT_BASE}.png\n  {OUT_BASE}.pdf")
