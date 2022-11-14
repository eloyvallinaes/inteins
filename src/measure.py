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
    if len(seq) == 0:
        return dict()

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


def measure_segments(segments):
    """
    Measure over a segments records dictionary.

    """
    data = []
    for acc, annotations in segments.items():
        for i, record in enumerate(annotations):
            m = physcoprops(record['seq'])
            m.update({
                'accession': acc,
                'segment_id': f"{acc}_{str(i)}",
            })
            data.append(m)

    return data


def measure(parser, entity, outname):
    """
    Common measure function for sequences annotated with one interpro code
    and entity either 'annots' or 'hosts'.

    """
    if entity not in ['annots', 'hosts']:
        raise ValueError(
            f"{entity} not understood; must be 'annots' or 'hosts'"
        )
    segments = parser.annots if entity == "annots" else parser.hosts
    # measure the segment collection
    data = measure_segments(segments)
    # add the taxid column
    for m in data:
        acc = m['accession']
        m.update({
            'taxid': parser.taxids[acc]
        })

    df = pd.DataFrame(data)
    df.to_csv(outname, index=False)


if __name__ == '__main__':
    fasta = Path("fasta/interpro")
    nterm = parse.Fasta2Dict(fasta / "IPR003587_head.fasta")
    cterm = parse.Fasta2Dict(fasta / "IPR003586.fasta")
    splicing = nterm + cterm
    measure(splicing, 'annots', 'splicing.phys')
