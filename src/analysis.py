"""
Abstracted functions for intein analysis and plots.
"""

import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

csv = Path("csv/interpro")


COLORS = {
    "Bacteria": ("lime", "tab:green"),
    "Archaea": ("lightcoral", "tab:red"),
    "Eukaryota": ("cyan", "tab:blue"),
    "Viruses": ("violet", "purple")
}

orthoGroups = {
    "OG6_104977": "hedgehog protein",
    "OG6_107432": "replicative DNA helicase",
    "OG6_118012": "hintN domain containing protein (OG6_118012)",
    "OG6_100448": "ribonucleotide reductase",
    "OG6_186238": "OG6_186238",
    "OG6_109432": "hintN domain containing protein (OG6_109432)",
    "OG6_107988": "DNA polymerase III-alpha",
    "OG6_100284": "DNA topoisomerase",
    "OG6_105177": "DNA gyrase",
    "OG6_101445": "tRNA-splicing ligase rtcB",
    "OG6_121348": "DNA polymerase II-large",
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

TAXCOLS = [
    'superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species'
]

linReg = LinearRegression()


def load_phys():
    # separate datasets
    hh = pd.read_csv(csv / "hosts.phys", dtype={"taxid": str})
    ii = pd.read_csv(csv / "inteins.phys")  # taxid column omitted
    # join properties from intein and intein
    return ii.merge(
        hh,
        on="accession",
        suffixes=("_intein", "_host")
    )


def load_mini():
    phys = load_phys()
    mini = pd.read_csv(csv / "mini.phys")
    endo = pd.read_csv(csv / "endo.phys")
    return phys.merge(
        mini.add_suffix("_mini"),
        left_on="intein_id",
        right_on="intein_id_mini"
    ).merge(
        endo.add_suffix("_endo"),
        left_on="intein_id",
        right_on="intein_id_endo"
    )


def load_interpro(mini=False):
    phys = load_mini() if mini is True else load_phys()
    # collect groups
    groups = pd.read_csv(csv / "groups.csv")[["accession", "group"]]
    # collect taxonomy
    taxonomy = pd.read_csv(csv / "taxonomy.csv", dtype={"taxid": str})
    # merge
    df = phys.merge(
        groups,
        on="accession",
        how='left'
    ).merge(
        taxonomy,
        on="taxid",
        how="left"
    )
    # add groupnames
    df["groupname"] = df.group.map(orthoGroups)
    return df


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


def load_merged_data(mini=True):
    """
    Return the interpro data merged with the proteome average data. Columns
    labeled with _proteome suffix.
    """
    # interpro sequences
    interpro = load_interpro(mini=mini)
    # proteome data
    proteomes = load_proteomes().drop_duplicates(
        subset="taxid",
        keep='first'
    )
    # merged dataset
    return interpro.merge(
            proteomes[PHYSCOLS + ['taxid']].add_suffix("_proteome"),
            left_on='taxid',
            right_on="taxid_proteome",
            how='left'
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


def entities_kde_ax(ax, data, varname, entities=[]):
    """
    Kernel density estimate plot for a variable across entities.

    """
    entlist = entities if entities else ['host', 'intein', 'endo', 'mini']
    for entity in entlist:
        sns.kdeplot(
            ax=ax,
            data=data,
            x=f"{varname}_{entity}",
            label=entity
        )
    ax.set_xlabel(f"segment {varname}")


def entities_proteome_ax(ax, data, varname, entities=[], **kwargs):
    """
    Scatter plot of variable across entities
    """
    for entity in ['host', 'intein', 'endo', 'mini']:
        sns.scatterplot(
            data=data,
            x=f"{varname}_proteome",
            y=f"{varname}_{entity}",
            label=entity,
            ax=ax,
            **kwargs
        )
    ax.axline([0, 0], slope=1, linestyle="--", color="k")
    ax.set_ylabel(f"segment {varname}")
    return ax


def get_accessions(data, short=True):
    """
    A handy function to print accessions such that you can copy/paste them in
    uniprot ID mapping.

    """
    lim = 6 if short else np.Inf
    print("\n".join(
        data[data.accession.str.len() <= lim].accession.values
    ))
