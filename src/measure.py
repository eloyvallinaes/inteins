import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
        self.sequences = pd.read_csv(seqfile, sep="\t")
        self.segments = self.extract_segments(ntermfile, ctermfile)
        self.accs = self.segments.refseq_accno.values
        self.inteindata = self.mark_inteins()
        self.records = self.measure()

    def extract_segments(self, ntermfile, ctermfile):
        nsegs = pd.read_csv(ntermfile)
        csegs = pd.read_csv(ctermfile)
        segments = pd.concat([nsegs, csegs])
        # Extract proteins with exactly two intein hits
        segcounts = segments.refseq_accno.value_counts()
        accs = segcounts[segcounts == 2].index
        return segments[segments.refseq_accno.isin(accs)]


    def mark_inteins(self):
        seqs =  self.sequences
        segments = self.segments
        accs = self.accs
        data = segments[segments.refseq_accno.isin(accs)]
        Ndata = data[data.motifname == "intein_Nterm"][["refseq_accno", "start"]]
        Cdata = data[data.motifname == "intein_Cterm"][["refseq_accno", "end"]]
        inteins =  Ndata.merge(Cdata, on = "refseq_accno", how = "inner")
        inteins["span"] = inteins.end - inteins.start
        inteins = inteins[inteins.span > 0]
        return inteins.merge(seqs[["refseq_accno","sequence"]], on="refseq_accno")

    @classmethod
    def measure_row(cls, row):
        acc = row.refseq_accno
        data = {}
        for name, seq in Measure.splitseq(row).items():
            L = len(seq)
            data.update( {f"{letter}_{name}": seq.count(letter)/L for letter in Measure.AA} )
            fPos = sum([seq.count(letter) for letter in "KR"]) / L
            fNeg = sum([seq.count(letter) for letter in "DE"]) / L
            fCharged = fPos + fNeg
            fFatty = sum([seq.count(letter) for letter in "FLIV"]) / L
            nc = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
            mw = L * 110
            sasa = 6.3 * mw**0.73
            ncd1000 =  nc / sasa * 1000
            mwkda = mw / 1000

            values = [L, sasa, mwkda, fCharged, fFatty, fPos, fNeg, nc, ncd1000]
            data.update( {f"{key}_{name}": val for key, val in zip(Measure.COLNAMES, values)} )
            data.update( {'refseq_accno': acc} )
        return data

    @classmethod
    def splitseq(cls, row):
        seq = row.sequence
        start = row.start
        end = row.end
        intein = seq[start:end]
        rnr = seq[0:start] + seq[end:-1]
        return {'rnr': rnr, 'intein': intein}


    def measure(self):
        records = []
        for ind, row in self.inteindata.iterrows():
            records.append( Measure.measure_row(row) )
        return pd.DataFrame(records)


    def writeInteins(self, outfilename):
        with open(outfilename, "a") as fastafile:
            for ind, row in self.inteindata.iterrows():
                seq = Measure.splitseq(row)["intein"]
                fastafile.write(f">{row.refseq_accno}|intein\n{seq}\n")

    def writeMeasurements(self, outfilename):
        self.records.to_csv(outfilename, index=False)



if __name__=='__main__':
    m = Measure("csv/rnr_refseq_proteins.tsv", "csv/rnr_nterm_segments.csv", "csv/rnr_cterm_segments.csv")
