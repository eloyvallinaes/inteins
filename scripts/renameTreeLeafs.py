#!/bin/python
"""
"""

import re
import pandas as pd

# interpro = analysis.load_interpro()
taxonomy = pd.read_csv("csv/IPR036844_taxonomy.csv")\
             .set_index("accession")["taxid"]\
             .to_dict()

pattern = re.compile("[A-Z][A-Z0-9]{5,}")


def acc2taxid(match):
    acc = match.group(0)
    return taxonomy[acc]


with open("trees/replicativeDNA/trees/helicase_host.tree", "r") as treefile:
    tree = treefile.read()
    newTree = re.sub(pattern, acc2taxid, tree)
