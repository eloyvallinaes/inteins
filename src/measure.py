"""
A Measure class for measuring physicochemical properties on proteins and their
inteins
"""

import pandas as pd
from pathlib import Path

from src import parse

AA = 'ACDEFGHIKLMNPQRSTVWY'
COLNAMES = [
    "length",
    "sasa",
    "mwkda",
    "fCharged",
    "fFatty",
    "fPos",
    "fNeg",
    "nc",
    "ncd1000",
]


def physcoprops(seq: str):
    """
    Measure amino acid composition and physicochemical properties of sequences.
    """
    # start props dict with aa composition
    L = len(seq)
    props = {f"{letter}": seq.count(letter) / L for letter in AA}
    # calculate physicochemical properties
    fPos = sum([seq.count(letter) for letter in "KR"]) / L
    fNeg = sum([seq.count(letter) for letter in "DE"]) / L
    fCharged = fPos + fNeg
    fFatty = sum([seq.count(letter) for letter in "FLIV"]) / L
    nc = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
    mw = L * 110
    sasa = 6.3 * mw**0.73
    ncd1000 = nc / sasa * 1000
    mwkda = mw / 1000
    # update props dict with physicochemical properties
    values = [L, sasa, mwkda, fCharged, fFatty, fPos, fNeg, nc, ncd1000]
    props.update({
        f"{key}": val
        for key, val in zip(COLNAMES, values)
    })
    return props


if __name__ == '__main__':
    fasta = Path("fasta/interpro")
    parser = parse.Fasta2Dict(fasta / 'IPR036844.fasta')
    # data = []
    # for acc, inteins in parser.inteins.items():
    #     for i, record in enumerate(inteins):
    #         m = physcoprops(record['seq'])
    #         m.update({'accession': acc, 'intein_id': f"{acc}_{str(i)}"})
    #         data.append(m)
    #
    # df = pd.DataFrame(data)
    # df.to_csv("inteins.phys", index=False)

    data = []
    for acc, hosts in parser.hosts.items():
        record = hosts[0]
        m = physcoprops(record['seq'])
        m.update({'accession': acc, 'taxid': parser.taxids[acc]})
        data.append(m)

    df = pd.DataFrame(data)
    df.to_csv("hosts.phys", index=False)
