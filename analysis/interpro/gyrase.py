import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from src import analysis


OGCODE = "OG6_105177"


# %% codecell
cellular = ["Bacteria", "Archaea", "Eukaryota"]
topGroups = {k: v for k, v in analysis.orthoGroups.items() if k == OGCODE}
colors = {k: v for k, v in analysis.COLORS.items() if k in cellular}
data = analysis.load_interpro()
data = data[data.group.isin(topGroups)]
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
fig = plt.figure(figsize=(15, 10))
axs = fig.subplots(2, 3).flatten()
for ax, varname in zip(axs, analysis.physcols):
    analysis.pairkde(ax, data, varname, "intein", "host")
    ax.set_title(f"{title} -- {varname}")


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
data.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
r2 = data[['ncd1000_host', 'ncd1000_intein']].corr().iloc[0, 1]**2
ax.set_title(f"{title}; R={r2.round(2)}")


# %% codecell
fig = plt.figure(figsize=(15, 5))
axs = fig.subplots(1, 3).flatten()
for ax, (label, subset) in zip(axs, data.groupby("kingdom")):
    subset.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
    r2 = subset[['ncd1000_host', 'ncd1000_intein']].corr().iloc[0, 1]**2
    ax.set_title(f"{title}; {label}; R={r2.round(2)}")
    regressor.fit(
        subset.ncd1000_host.values.reshape(-1, 1),
        subset.ncd1000_intein.values.reshape(-1, 1)
    )
    x = np.linspace(subset.ncd1000_host.min(), subset.ncd1000_host.max(), 100)
    y = regressor.predict(x.reshape(-1, 1))
    ax.plot(x, y, "--r")


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
sns.kdeplot(
    data=data[data.kingdom == "Bacteria"],
    x="ncd1000_host",
    y="ncd1000_intein",
    hue="kingdom",
)
ax.set_aspect('equal', 'datalim')


# %% codecell
grid = sns.jointplot(
    data=data[data.kingdom == "Bacteria"],
    x="ncd1000_host",
    y="ncd1000_intein",
    kind="reg",
    color="#4CB391",
    height=5,
    ratio=2,
)
grid.figure.axes[0].set(xlim=(-2, 2), ylim=(-2, 2))
grid.plot_joint(sns.kdeplot)

# %% codecell
grid = sns.jointplot(
    data=data[(data.kingdom == "Bacteria") & (data.ncd1000_host >= 0)],
    x="ncd1000_host",
    y="ncd1000_intein",
    kind="reg",
    color="#4CB391",
    height=5,
    ratio=2,
)
# grid.plot_joint(sns.kdeplot)


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
subset = data[data.kingdom == "Bacteria"]
subset.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
r = subset[['ncd1000_host', 'ncd1000_intein']].corr().iloc[0, 1]
ax.set_title(f"{title}; R={r.round(2)}")
ax.axvline(0, color="k", linestyle="--")
ax.axhline(0, color="k", linestyle="--")
