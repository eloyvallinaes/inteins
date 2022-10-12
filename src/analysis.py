"""
Abstracted functions for intein analysis and plots.
"""

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
    "Viruses": ("violet", "purple")
}

orthoGroups = {
    "OG6_104977": "hedgehog protein",
    "OG6_107432": "replicative DNA helicase",
    "OG6_118012": "N/A",
    "OG6_100448": "ribonucleotide reductase",
    "OG6_186238": "hintN domain containing protein",
    "OG6_109432": "hintN domain containing protein",
    "OG6_107988": "DNA polymerase III",
    "OG6_100284": "DNA topoisomerase",
    "OG6_105177": "DNA gyrase",
    "OG6_101445": "tRNA-splicigng ligase rtcB",
    "OG6_121348": "DNA polymerase II",
    "OG6_104530": "protein recA",
    "OG6_128365": "DNA-methyltransferase"
}

AACOLS = [
    'A',
    'C',
    'D',
    'E',
    'F',
    'G',
    'H',
    'I',
    'K',
    'L',
    'M',
    'N',
    'P',
    'Q',
    'R',
    'S',
    'T',
    'V',
    'W',
    'Y',
]

PHYSCOLS = [
    'fCharged',
    'fFatty',
    'fNeg',
    'fPos',
    'nc',
    'ncd1000',
]

linReg = LinearRegression()


def load_interpro():
    hh = pd.read_csv(csv / "IPR036844_host_phys.csv", dtype={"taxid": str})
    ii = pd.read_csv(csv / "IPR036844_intein_phys.csv", dtype={"taxid": str})

    phys = ii.merge(hh, on="refseq_accno", suffixes=("_intein", "_host"))
    groups = pd.read_csv("csv/ipr036844_groups_map.csv")[["accession", "group"]]
    taxonomy = pd.read_csv("csv/ipr036844_taxonomy.csv", dtype={"taxid": str})

    grouped = phys.merge(
        groups,
        left_on="refseq_accno",
        right_on="accession",
    )

    return grouped.merge(
        taxonomy,
        on="accession",
    )


def pairkde(
    ax, df, x, i, j, colors=(
        "tab:orange", "tab:cyan"), labels=(
            "", "")):
    """
    Kernel Density Estimate plots for property x in groups {i, j} for a
    DataFrame df such that there are columns x_i and x_j.
    """
    df[f"{x}_{i}"].plot.kde(label=i, color=colors[0], ax=ax)
    df[f"{x}_{j}"].plot.kde(label=j, color=colors[1], ax=ax)
    ax.legend()
    return ax


def pairplot(ax, df, x, y, i="intein", j="host", color="w", label=""):
    """
    Scatter plot connecting properties (x,y) between groups {i, j} assuming
    that df has columns names x_i, x_j, y_i and y_j.
    """
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
    # connecting lines
    for _, series in df.iterrows():
        ax.plot(
            [series[f"{x}_{i}"], series[f"{x}_{j}"]],
            [series[f"{y}_{i}"], series[f"{y}_{j}"]],
            linewidth=0.5,
            marker="",
            zorder=0,
            color="k",
        )


def load_proteomes(db="ncbi", metric="avg"):
    if db not in ["ncbi", "uniprot"]:
        raise ValueError(f"db must be 'ncbi' or 'uniprot', not {db}")

    if metric not in ["avg", "median"]:
        raise ValueError(f"metric must be 'avg' or 'median', not {metric}")

    datasets = Path("analysis/proteomeData")
    pattern = f"*_{db}_{metric}.csv"

    df = pd.concat(pd.read_csv(dataset) for dataset in datasets.glob(pattern))
    df.columns = [
        'length',
        'mass',
        'SASAmiller',
        'nc',
        'ncd',
        'fFatty',
        'fPos',
        'fNeg',
        'Species',
        'taxid',
        'seqnum',
        'fCharged',
        'ncd1000',
        'mwkda',
        'I', 'L', 'V', 'A', 'G', 'P', 'F', 'M', 'W', 'Y', 'H', 'T', 'S', 'N',
        'Q', 'D', 'E', 'K', 'R', 'C', 'X', 'U', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species', 'GC', 'AssemblyID'
    ]
    df.taxid = df.taxid.astype(str)
    return df.drop(["Species", "mass"], axis="columns")


def taxkde(data, x, rank, name, into, n=10) -> plt.figure:
    fig = plt.figure(figsize=(10, 5))
    axs = fig.subplots(1, 2, sharex=True, sharey=True)
    top10 = data[data[rank] == name][into].value_counts().index[:n]
    subset = data[
        (data[into].isin(top10)) &
        (data[rank] == name)
    ]

    for var, ax in zip(["host", "intein"], axs):
        g = sns.kdeplot(
            data=subset,
            x=f"{x}_proteome",
            y=f"{x}_{var}",
            hue=into,
            hue_order=top10,
            ax=ax,
        )
        ax.set_title(var)
        ax.set_xlabel(f"{x} organism avg.")
        ax.set_ylabel(f"{x} protein")

    fig.get_axes()[0].get_legend().remove()
    handles = fig.get_axes()[1].get_legend().get_lines()
    texts = [t.get_text() for t in fig.get_axes()[1].get_legend().get_texts()]
    labels = [
        f"{name} (n={subset[subset[into] == name].shape[0]})"
        for name in texts
    ]
    fig.get_axes()[1].legend(
        handles,
        labels,
        bbox_to_anchor=(1.5, 0.9),
        framealpha=1
    )
    return fig


def load_cogs(*cogids):
    """
    Load COG  data into dataframes and merge with suffixes for the intein and
    the host protein.
    """
    suffs = ("_intein", "_host")
    frames = []
    for cogid in cogids:
        h = pd.read_csv(csv / f"{cogid}_host_phys.csv")
        i = pd.read_csv(csv / f"{cogid}_intein_phys.csv")
        h["subset"] = cogid
        frames.append(
            i.merge(h, on="refseq_accno", suffixes=suffs, how='inner')
        )
    return pd.concat(frames)


def corrScatter(data, groupcol, x, y, n=5):
    """
    Make a figure with n subgroups and plot linear regression fits above each
    data scatter.
    """
    palette = sns.color_palette("muted", n)
    groups = data[groupcol].value_counts().index[:n]
    fig = plt.figure(figsize=(6, 6))
    ax = fig.subplots(1, 1)
    for group, color in zip(groups, palette):
        subset = data[data[groupcol] == group]
        r2 = subset[[x, y]].corr().iloc[0, 1]**2
        subset.plot.scatter(
            x=x,
            y=y,
            ax=ax,
            color=color,
            edgecolor="w",
            linewidth=0.2,
            label=f"{group} (n={subset.shape[0]}); R2={r2.round(2)}"
        )
        ax.legend(bbox_to_anchor=(1, 1))
        # regression lines
        linReg.fit(
            subset[x].values.reshape(-1, 1),
            subset[y].values.reshape(-1, 1)
        )
        xx = np.linspace(subset[x].min(), subset[x].max(), 100)
        yy = linReg.predict(xx.reshape(-1, 1))
        ax.plot(xx, yy, "--", color=color)
    return fig


def corrKde(data, groupcol, x, y, n=5):
    """
    Make a figure with n subgroups as 2D kde contours.
    """
    groups = data[groupcol].value_counts().index[:n]
    fig = plt.figure(figsize=(6, 6))
    ax = fig.subplots(1, 1)
    sns.kdeplot(
        data=data,
        x=x,
        y=y,
        ax=ax,
        hue=groupcol,
        hue_order=groups,
    )
    return fig


def corrCoefficients(data, varnames=PHYSCOLS, groupcol="kingdom"):
    """
    Find correlation coefficients over subgroups and variables.
    """
    for k, subset in data.groupby(groupcol):
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


def reframeToEntities(
    data: pd.DataFrame,
    columns: list[str],
    categories: list[str],
    suffixes=("intein", "host"),
    newcolname='entity'
):
    """
    Reorganise datraframe. Column suffixes are turned into a
    new categorical column. Example:

    >>>df = pd.DataFrame({"a_x": 0, "a_y": 1, "name": "rhino"})
    >>>print(df)
        a_x a_y name
    0   0   1   rhino
    >>>df = reframeToEntities(
        df,
        columns=["a"],
        categories=["rhino"],
        suffixes=("x", "y")
    )
        a    name   entity
    0   0   rhino   x
    1   1   rhino   y
    """

    dff = []
    for suff in suffixes:
        labels = [f"{c}_{suff}" for c in columns] + categories
        df = data.reindex(labels, axis='columns')
        df.columns = columns + categories
        df[newcolname] = suff
        dff.append(df)

    return pd.concat(dff, axis=0).reset_index(drop=True)
