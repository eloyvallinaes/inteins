"""
A Measure class for measuring physicochemical properties on proteins and their
inteins
"""

import pandas as pd
from parse import InterproSegments, COGSegments


class Measure:
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

    @classmethod
    def physcoprops(cls, seq, acc):
        # aa composition
        L = len(seq)
        props = {f"{letter}": seq.count(letter) / L for letter in Measure.AA}
        # accession
        props['refseq_accno'] = acc
        # properties
        fPos = sum([seq.count(letter) for letter in "KR"]) / L
        fNeg = sum([seq.count(letter) for letter in "DE"]) / L
        fCharged = fPos + fNeg
        fFatty = sum([seq.count(letter) for letter in "FLIV"]) / L
        nc = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
        mw = L * 110
        sasa = 6.3 * mw**0.73
        ncd1000 = nc / sasa * 1000
        mwkda = mw / 1000

        values = [L, sasa, mwkda, fCharged, fFatty, fPos, fNeg, nc, ncd1000]
        props.update({
            f"{key}": val
            for key, val in zip(Measure.COLNAMES, values)
        })
        return props

    def measure_inteins(self):
        data = []
        for acc, regions in self.regions.items():
            for intein in regions['inteins']:
                data.append(Measure.physcoprops(intein['seq'], acc))
        return data

    def measure_host(self):
        data = []
        for acc, regions in self.regions.items():
            data.append(Measure.physcoprops(regions['host'], acc))
        return data

    def write_hostdata(self, outname):
        pd.DataFrame(self.hostdata).to_csv(outname, index=False)

    def write_inteindata(self, outname):
        pd.DataFrame(self.inteindata).to_csv(outname, index=False)


class MeasureCOG(Measure):
    def __init__(self, seqfile, ntermfile, ctermfile):
        self.regions = COGSegments.parse(seqfile, ntermfile, ctermfile)
        self.hostdata = self.measure_host()
        self.inteindata = self.measure_inteins()


class MeasureInterpro(Measure):
    def __init__(self, fastafilename):
        self.regions = InterproSegments.parse(fastafilename)
        self.hostdata = self.measure_host()
        self.inteindata = self.measure_inteins()


if __name__ == '__main__':
    from pathlib import Path
    fasta = Path("fasta")
    cogs = Path("cogs")
    terms = Path("terms")
    for subset in ["IPR036844"]:
        seqfile = fasta / f"{subset}.fasta"
        m = MeasureInterpro(seqfile)


    for subset in ["COG0417"]:
        seqfile = cogs / (subset + ".tsv")
        ntermfile = terms / f"TIGR01445.1_{subset}.csv"
        ctermfile = terms / f"TIGR01443.1_{subset}.csv"
        m = MeasureCOG(seqfile, ntermfile, ctermfile)
