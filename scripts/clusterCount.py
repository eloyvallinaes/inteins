#!/bin/python

import re
import pandas as pd
import matplotlib.pyplot as plt


def countSorted(filename):
    counts = {}
    with open(filename, "r") as clusterFile:
        cc = 0
        clusterId = ""
        while True:
            line = clusterFile.readline()
            if not line:
                break

            elif line.startswith(">"):
                if clusterId and cc > 0:
                    counts[clusterId] = cc
                    clusterId = ""

                clusterId = int(re.search("[0-9]{1,}", line).group(0))
                cc = 0
            else:
                cc += 1

    return counts


counts = countSorted("cdhit/1663922554.fas.1.clstr.sorted")
cc = pd.DataFrame.from_dict(counts, orient="index", columns=["clustSize"])


# %% codecell
fig = plt.figure()
ax = fig.subplots(1, 1)
ax.plot(cc.clustSize.cumsum(), marker="o")
ax.plot([0, 20000], [10000, 30000], "-r")
# ax.set_xlim([0, 400])
# ax.set_ylim([0, 3500])


# %% codecell
cc[cc.clustSize > 10]
