"""
Analysis of OG6_121348 — DNA polymerase II in Bacteria.

Selected because of negative ncd1000_host values.
"""
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from src import analysis

OGCODE = "OG6_121348"  # polymerase II
KINGDOMS = ["Archaea"]

# %% codecell
topGroups = {k: v for k, v in analysis.orthoGroups.items() if k == OGCODE}
colors = {k: v for k, v in analysis.COLORS.items() if k in KINGDOMS}
data = analysis.load_interpro()
data = data[
    (data.group == OGCODE) &
    (data.kingdom.isin(KINGDOMS))
]
regressor = LinearRegression()
title = f"{topGroups[data.group.iloc[0]]}"


# %% codecell
# Correlations
for k, subset in data.groupby("kingdom"):
    if subset.shape[0] < 100:
        continue
    correlations = {}
    for varname in analysis.physcols:
        correlations[varname] = subset[
            [f"{varname}_intein", f"{varname}_host"]
        ].corr().iloc[0, 1].round(2)

    print(
        k, f"(n={subset.shape[0]})",
        {k: v for k, v in correlations.items() if abs(v) > 0.5}
    )


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
for k, subset in data.groupby("kingdom"):
    analysis.pairkde(ax, subset, "ncd1000", "intein", "host")
    ax.set_title(f"{title} — {k} (n={subset.shape[0]})")


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
for k, subset in data.groupby("kingdom"):
    subset.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
    regressor.fit(
        subset.ncd1000_host.values.reshape(-1, 1),
        subset.ncd1000_intein.values.reshape(-1, 1)
    )
    x = np.linspace(subset.ncd1000_host.min(), subset.ncd1000_host.max(), 100)
    y = regressor.predict(x.reshape(-1, 1))
    ax.plot(x, y, "--r")
    r = subset[['ncd1000_host', 'ncd1000_intein']].corr().iloc[0, 1]
    ax.set_title(f"{title}; {k}; R={r.round(2)}")


# %% codecell
g = sns.jointplot(
    data=data,
    x="ncd1000_host",
    y="ncd1000_intein",
    hue="kingdom",
    kind="kde"
)
ax = g.figure.axes[0]
for k, subset in data.groupby("kingdom"):
    ax.axvline(
        subset.ncd1000_host.median(),
        color='k', linestyle="--", zorder=-1
    )
    ax.axhline(
        subset.ncd1000_intein.median(),
        color='k', linestyle="--", zorder=-1
    )


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
for k, subset in data.groupby("kingdom"):
    subset.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
    ax.set_title(f"{title}; {k} (n={subset.shape[0]})")
    ax.axvline(subset.ncd1000_host.median(), color="k", linestyle="--")
    ax.axhline(subset.ncd1000_intein.median(), color="k", linestyle="--")
