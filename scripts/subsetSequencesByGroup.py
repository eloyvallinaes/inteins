#!/bin/python

from pathlib import Path
from src import analysis, parse

"""
Prepare subsets of sequences in FASTA format for COG categories.
"""


def prepareData():
    interpro = analysis.load_interpro()
    interpro["groupname"] = interpro.group.map(analysis.orthoGroups)
    jj = interpro.dropna(
        subset=["groupname"]
        ).drop_duplicates(
            ["taxid"],
            keep="first"
        )
    # unique accession and taxid
    assert jj.taxid.is_unique and jj.accession.is_unique
    return jj


def fastaSubset(df, entity, groupname):
    """
    Load FASTA file as dictionary and write subset by accession.
    """
    fasta = Path("fasta")
    sequences = fasta / f"IPR036844_{entity}.fasta"
    accessions = df[df.groupname == groupname].accession.unique().tolist()
    parser = parse.Fasta2Dict(sequences, subset=accessions)
    label = "_".join(groupname.split())
    parser.write(fasta / f"IPR036844_{entity}_{label}.fasta")


def main():
    df = prepareData()
    for entity in ["host", "intein"]:
        for groupname in ["hedgehog protein", "replicative DNA helicase"]:
            fastaSubset(df, entity, groupname)


if __name__ == '__main__':
    main()
