"""
A Measure class for measuring physicochemical properties on proteins and their
inteins
"""

import pandas as pd


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

    def __init__(self, seqfile, ntermfile, ctermfile):
        self.sequences = pd.read_csv(seqfile, sep="\t")\
            .to_dict(orient="records")
        self.segments = self.extract_segments(ntermfile, ctermfile)\
                            .to_dict(orient="records")
        self.seqmap = {
            row['refseq_accno']: row['sequence']
            for row in self.sequences
        }
        self.regions = self.mark_regions()
        self.hostdata = self.measure_host()
        self.inteindata = self.measure_inteins()

    def extract_segments(self, ntermfile, ctermfile):
        xcolumns = ["end", "motifname", "motifacc"]
        ycolumns = ["start", "motifname", "motifacc"]
        nterm = pd.read_csv(ntermfile).drop(xcolumns, axis="columns")
        cterm = pd.read_csv(ctermfile).drop(ycolumns, axis="columns")
        ss = nterm.merge(cterm, on="refseq_accno")
        ss['span'] = ss.end - ss.start
        # select index combinations giving positive, shortest span for each acc
        # and output sorted by starting positions -- proper sequence slicing
        return ss[ss.span > 0].sort_values("span")\
                              .drop_duplicates(
            ["refseq_accno", "start"],
            keep="first")\
            .sort_values("start")

    def mark_regions(self):
        data = {}
        for acc, sequence in self.seqmap.items():
            segments = [segment
                        for segment in self.segments
                        if segment['refseq_accno'] == acc
                        ]
            inteins = Measure.extract_inteins(sequence, segments)
            hostseq = Measure.extract_host(sequence, inteins)
            data[acc] = {'inteins': inteins, 'host': hostseq}
        return data

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

    @classmethod
    def extract_host(cls, sequence, inteins):
        indeces = [0]
        for i in inteins:
            indeces += [i['start']] + [i['end']]
        indeces += [-1]

        host = ''
        for e in range(0, len(indeces), 2):
            host += sequence[indeces[e]:indeces[e + 1]]

        return host

    @classmethod
    def extract_inteins(cls, sequence, segments):
        inteins = []
        for segment in segments:
            start = segment['start']
            end = segment['end']
            inteins.append({
                'start': start,
                'end': end,
                'span': end - start,
                'seq': sequence[start:end]
            })
        return inteins


if __name__ == '__main__':
    from pathlib import Path
    cogs = Path("cogs")
    csv = Path("csv")
    terms = Path("terms")
    results = {}
    for subset in ["COG0209"]:
        seqfile = cogs / f"{subset}.tsv"
        nterms = terms / f"TIGR01445.1_{subset}.csv"
        cterms = terms / f"TIGR01443.1_{subset}.csv"

        m = Measure(seqfile, nterms, cterms)
        segments = m.extract_segments(nterms, cterms)
