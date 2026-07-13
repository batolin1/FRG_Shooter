"""
Interactive 3D plot: Dimension vs Anomalous dimension vs Sigma, by multicritical order.
Requires: pandas, plotly
    pip install pandas plotly

Produces an HTML file you can open in a browser and rotate/zoom freely.
"""

import pandas as pd
import plotly.graph_objects as go

# --- Load data -----------------------------------------------------------
# Expects a CSV with columns: Category,Dimension,Anomalous,S,Sigma,Rho
df = pd.read_csv("data.csv")

categories = ["Bicritical", "Tricritical", "Quadricritical", "Pentacritical", "Hexacritical"]
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

# --- Build one 3D line+marker trace per category --------------------------
traces = []
for cat, color in zip(categories, colors):
    sub = df[df["Category"] == cat].sort_values("Dimension", ascending=False)
    traces.append(
        go.Scatter3d(
            x=sub["Dimension"],
            y=sub["Anomalous"],
            z=sub["Sigma"],
            mode="lines+markers",
            name=cat,
            line=dict(color=color, width=4),
            marker=dict(color=color, size=4),
        )
    )

# --- Layout ----------------------------------------------------------------
layout = go.Layout(
    title=dict(text="Dimension vs Anomalous dimension vs Sigma, by critical order", font=dict(size=24)),
    scene=dict(
        xaxis=dict(title=dict(text="Dimension (d)", font=dict(size=18)),
                    autorange="reversed", tickfont=dict(size=13)),
        yaxis=dict(title=dict(text="Anomalous dimension", font=dict(size=18)),
                    tickfont=dict(size=13)),
        zaxis=dict(title=dict(text="Sigma", font=dict(size=18)),
                    tickfont=dict(size=13)),
    ),
    margin=dict(l=0, r=0, b=0, t=60),
    legend=dict(title=dict(text="Critical order", font=dict(size=16)), font=dict(size=14)),
)

fig = go.Figure(data=traces, layout=layout)

# --- Export ------------------------------------------------------------------
fig.write_html("3d_interactive.html")
fig.show()
