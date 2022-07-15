import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

csv = Path("csv")
COLORS = {
    "Bacteria": ("lime", "tab:green"),
    "Archaea": ("lightcoral", "tab:red"),
    "Eukaryota": ("cyan", "tab:blue"),
    # "Viruses": ("violet", "purple")
}

topGroups = {
    # "OG6_104977",
    # "OG6_107432": "replicative DNA helicase",  # replicative DNA helicase
    # "OG6_118012",
    "OG6_100448": "ribonucleotide reductase",  # rnr
    # "OG6_186238",
    # "OG6_109432",
    # "OG6_107988",
    # "OG6_100284",
    # "OG6_105177",
    # "OG6_101445",
    # "OG6_121348",
}

varnames = [
    # 'A',
    # 'C',
    # 'D',
    # 'E',
    # 'F',
    # 'G',
    # 'H',
    # 'I',
    # 'K',
    # 'L',
    # 'M',
    # 'N',
    # 'P',
    # 'Q',
    # 'R',
    # 'S',
    # 'T',
    # 'V',
    # 'W',
    # 'Y',
    'fCharged',
    'fFatty',
    'fNeg',
    'fPos',
    'length',
    'nc',
    'ncd1000',
]


# %% codecell
def pairplot(ax, df, x, y, i, j, color="w", label=""):
    props = dict(s=30, edgecolor="k", ax=ax, alpha=.5)
    df.plot.scatter(
        f"{x}_{i}", f"{y}_{i}",
        marker="X",
        label=label,
        color=color,
        **props
    )
    df.plot.scatter(
        f"{x}_{j}", f"{y}_{j}",
        marker="o",
        color="w",
        **props
    )
    for _, series in df.iterrows():
        ax.plot(
            [series[f"{x}_{i}"], series[f"{x}_{j}"]],
            [series[f"{y}_{i}"], series[f"{y}_{j}"]],
            linewidth=0.5,
            marker="",
            zorder=0,
            color="k",
        )


def pairkde(
    ax, df, x, i, j, colors=(
        "tab:orange", "tab:cyan"), labels=(
            "", "")):
    df[f"{x}_{i}"].plot.kde(label=i, color=colors[0], ax=ax)
    df[f"{x}_{j}"].plot.kde(label=j, color=colors[1], ax=ax)
    ax.legend()
    return ax


def load():
    hh = pd.read_csv(csv / "IPR036844_host_phys.csv")
    ii = pd.read_csv(csv / "IPR036844_intein_phys.csv")

    phys = ii.merge(hh, on="refseq_accno", suffixes=("_intein", "_host"))
    groups = pd.read_csv("csv/ipr036844_groups_map.csv")[["accession", "group"]]
    taxonomy = pd.read_csv("csv/ipr036844_taxonomy.csv")

    grouped = phys.merge(
        groups,
        left_on="refseq_accno",
        right_on="accession",
    )

    return grouped.merge(
        taxonomy,
        on="accession",
    )


# %% codecell
data = load()
data = data[
    data.group.isin(topGroups) &
    data.kingdom.isin(["Bacteria", "Archaea"])
]
regressor = LinearRegression()
title = f"{topGroups[data.group.iloc[0]]}"


# %% codecell
# Correlations
for k, subset in data.groupby(["kingdom"]):
    if subset.shape[0] < 100:
        continue
    correlations = {}
    for varname in varnames:
        correlations[varname] = subset[
            [f"{varname}_intein", f"{varname}_host"]
        ].corr().iloc[0, 1].round(2)

    print(
        k, f"(n={subset.shape[0]})",
        {k: v for k, v in correlations.items() if abs(v) > 0.5}
    )


# %% codecell
fig = plt.figure()
ax = fig.subplots(1, 1)
pairkde(ax, data, "ncd1000", "intein", "host")
ax.set_title(title)


# %% codecell
fig = plt.figure(figsize=(5, 5))
ax = fig.subplots(1, 1)
data.plot.scatter(x="ncd1000_host", y="ncd1000_intein", ax=ax)
r2 = data[['ncd1000_host', 'ncd1000_intein']].corr().iloc[0, 1]**2
ax.set_title(f"{title}; R={r2.round(2)}")


# %% codecell
fig = plt.figure(figsize=(10, 5))
axs = fig.subplots(1, 2).flatten()
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
ax = sns.kdeplot(
    data=data,
    x="ncd1000_host",
    y="ncd1000_intein",
    hue="kingdom",
    levels=15,
)
ax.set_aspect('equal', 'box')
ax.set(xlim=(-6, 4), ylim=(-6, 4))
