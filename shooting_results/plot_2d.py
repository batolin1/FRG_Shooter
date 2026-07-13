"""
2D plot: Anomalous dimension vs Dimension, by multicritical order.
Requires: pandas, matplotlib
    pip install pandas matplotlib
"""

import pandas as pd
import matplotlib.pyplot as plt

# --- Load data -----------------------------------------------------------
# Expects a CSV with columns: Category,Dimension,Anomalous,S,Sigma,Rho
df = pd.read_csv("data.csv")

categories = ["Bicritical", "Tricritical", "Quadricritical", "Pentacritical", "Hexacritical"]
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"]
markers = ["o", "s", "^", "D", "v", "*", "X"]
theoretical_pts = [1, 1.33, 1.5, 1.6, 1.66, 1.714, 1.75, 1.77, 1.80]

# --- Plot ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
i = 0
for cat, color, marker in zip(categories, colors, markers):
    sub = df[df["Category"] == cat].sort_values("Dimension", ascending=False)
    ax.plot(
        sub["Dimension"], sub["Sigma"],
        marker=marker, color=color, label=cat,
        linewidth=1.8, markersize=6,
    )
    # ax.scatter (theoretical_pts[i], 0.0, marker="D", color=color)
    # i+=1


ax.set_xlabel("Long Range Exponent s", fontsize=12)
ax.set_ylabel("Initial Condition $\sigma$", fontsize=12)
ax.set_title("Emergence of Multicritical Solutions as function of Long Range Exponent", fontsize=13)
ax.invert_xaxis()
ax.legend(title="Critical order", frameon=False)
ax.grid(True, which="both", linestyle="--", alpha=1)

fig.tight_layout()
fig.savefig("anomalous_vs_dimension.png", dpi=150)
plt.show()
