import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "short-range.csv"   # <-- point this at your file

# ----------------------------------------------------------------------
# Read the two header rows:
#   row 1 = category name (merged cells -> blank for repeated columns)
#   row 2 = column name, e.g. "dimension", "y2,1", "y3,1", ...
# ----------------------------------------------------------------------
df = pd.read_csv(CSV_PATH, header=[0, 1])

# pandas turns blank merged cells into "Unnamed: N_level_0" -> treat as NaN,
# then forward-fill so every y-column knows which category it belongs to.
top = pd.Series(df.columns.get_level_values(0))
top = top.where(~top.str.startswith("Unnamed", na=False), np.nan).ffill()
sub = df.columns.get_level_values(1)
df.columns = pd.MultiIndex.from_arrays([top, sub])

# Locate the dimension column (its sub-label is "dimension")
dim_col = [c for c in df.columns if c[1].strip().lower() == "dimension"][0]
dimension = df[dim_col].astype(float).values

# ----------------------------------------------------------------------
# Style: colors cycle per category, markers cycle per k-index.
# Add more entries here if you have more than 7 categories / 7 k-values.
# ----------------------------------------------------------------------
color_cycle = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"]
marker_cycle = ["o", "s", "^", "D", "v", "*", "X"]

color_map = {}

fig, ax = plt.subplots(figsize=(8, 6))

for category, col in df.columns:
    if col.strip().lower() == "dimension":
        continue

    m = re.match(r"\s*y\s*(\d+)\s*,\s*(\d+)\s*", col)
    if not m:
        continue  # skip any column that isn't of the form y<n>,<k>
    n, k = int(m.group(1)), int(m.group(2))

    if category not in color_map:
        color_map[category] = color_cycle[len(color_map) % len(color_cycle)]

    y = pd.to_numeric(df[(category, col)], errors="coerce").values

    ax.plot(
        dimension, y,
        color=color_map[category],
        marker=marker_cycle[(k - 1) % len(marker_cycle)],
        markersize=6,
        linewidth=1.6,
        label=fr"$y_{{{n},{k}}}$ ({category})",
    )

ax.set_xlabel("dimension")
ax.set_ylabel(r"$y_{n,k}$")
ax.set_title("RG eigenvalues as a function of effective dimension")
ax.invert_xaxis()  # data runs from high dimension down to ~2
ax.legend(fontsize=8, ncol=2, loc="best")
ax.grid(alpha=1, linestyle="--")

fig.tight_layout()
fig.savefig("ynk_vs_dimension.png", dpi=200)
plt.show()