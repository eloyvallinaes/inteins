"""
Comparisons across orthologous groups.
"""

import numpy as np
import seaborn as sns
from src import analysis
import matplotlib.pyplot as plt

# %% codecell
data = analysis.load_interpro()

# %% codecell
data["groupname"] = data.group.map(analysis.orthoGroups)
fig = plt.figure(figsize=(15, 5))
ax = fig.subplots(1, 1)
sns.violinplot(
    data=data[
        (data.kingdom == "Bacteria") &
        (data.group.isin(analysis.orthoGroups))
    ],
    x="groupname",
    y="ncd1000_host",
    ax=ax,
    common_norm=True,
)
ax.tick_params(axis='x', labelrotation=30)


# %% codecell
data.groupby("group")[
    "ncd1000_host"
].describe()[
    ["count", "50%"]
].sort_values(
    ["50%", "count"],
    ascending=(False, False)
).head(12)


# %% codecell
piv = data.pivot_table(
    index="kingdom",
    columns="group",
    values="ncd1000_host",
    aggfunc=(np.median, 'count'),
)

# %% codecell
piv.loc[
    slice(None),
    ("count", list(analysis.orthoGroups.keys()))
]

# %% codecell
piv.loc[
    slice(None),
    ("median", list(analysis.orthoGroups.keys()))
]
