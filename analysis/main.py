import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

csv = Path("csv")
subsets = {
    "COG0209": "RNR (Archae and Bacteria)",
    "COG0417": "Helicase (Archaea)",
    "COG0305": "Helicase (Bacteria)",
    "COG0086": "RNA polymerase (Archaea and Bacteria)",
    "COG0587": "DNA polymerase subunit alpha (Archaea and Bacteria)",
    "COG0468": "Recombinase (Archaea and Bacteria)",
    # "COG5362": "Terminase 6"
}
colors = {
    'COG0209': "tab:orange",
    "COG0305": "tab:green",
    "COG0417": "tab:red",
    "COG0086": "tab:blue",
    "COG0587": "tab:pink",
    "COG0468": "tab:purple",
    # "COG5362": "tab:olive"
}


# %% codecell
def pairplot(ax, df, x, y, i, j, color="w", label=""):
    props = dict(color=color, s=40, edgecolor="k", ax=ax)
    df.plot.scatter(f"{x}_{i}", f"{y}_{i}", marker="X", label=label, **props)
    df.plot.scatter(f"{x}_{j}", f"{y}_{j}", marker="o", **props)
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


# %% codecell
hosts = []
inteins = []
for cogid, name in subsets.items():
    h = pd.read_csv(csv / f"{cogid}_host_phys.csv")
    i = pd.read_csv(csv / f"{cogid}_intein_phys.csv")
    h["subset"] = cogid
    hosts.append(h)
    inteins.append(i)

hh = pd.concat(hosts)
ii = pd.concat(inteins)

data = ii.merge(hh, on="refseq_accno", suffixes=("_intein", "_host"))


# %% codecell
data[data["refseq_accno"] == 'WP_011011558.1'][[
    "refseq_accno", "ncd1000_host", "ncd1000_intein"
]]


# %% codecell
for cogid, subdata in data.groupby("subset"):
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    pairplot(
        ax,
        subdata,
        "ncd1000",
        "mwkda",
        "intein",
        "host",
        color=colors[cogid],
        label=subsets[cogid]
    )


# %% codecell
for cogid, subdata in data.groupby("subset"):
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    pairkde(
        ax,
        subdata,
        "ncd1000",
        "intein",
        "host",
    )
    ax.set_title(subsets[cogid] + f" (n = {subdata.shape[0]})")

# %% codecell
# def phyloplot(ax, df, x, y, rank):
#     names = df[rank].value_counts().index
#     cmap = sns.color_palette("tab10")
#     for name, color in zip(names, cmap):
#         data = df[df[rank] == name]
#         pairplot(ax, data, x, y, "intein", "rnr", color = color, label=name)
#     ax.set_xlabel(x)
#     ax.set_ylabel(y)
#     ax.legend(bbox_to_anchor=(1,1), markerscale=1.5)
