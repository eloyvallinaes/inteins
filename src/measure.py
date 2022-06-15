import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def load_data():
    sequences = pd.read_csv("csv/rnr_refseq_proteins.tsv", sep="\t")
    nsegs = pd.read_csv("csv/nterm_segments.csv")
    csegs = pd.read_csv("csv/cterm_segments.csv")
    segments = pd.concat([nsegs, csegs])

    # Extract proteins with exactly two intein hits
    segcounts = segments.refseq_accno.value_counts()
    accs = segcounts[segcounts == 2].index

    return sequences, segments, accs


def mark_inteins(seqs, segments, accs):
    data = segments[segments.refseq_accno.isin(accs)]
    Ndata = data[data.motifname == "intein_Nterm"][["refseq_accno", "start"]]
    Cdata = data[data.motifname == "intein_Cterm"][["refseq_accno", "end"]]
    inteins =  Ndata.merge(Cdata, on = "refseq_accno", how = "inner")
    inteins["span"] = inteins.end - inteins.start
    return inteins.merge(seqs[["refseq_accno","sequence"]], on="refseq_accno")


def measure(row):
    acc = row.refseq_accno
    data = {}
    for name, seq in splitseq(row).items():
        L = len(seq)
        data.update( {f"{letter}_{name}": seq.count(letter)/L for letter in AA} )
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
        data.update( {f"{key}_{name}": val for key, val in zip(COLNAMES, values)} )
        data.update( {'refseq_accno': acc} )
    return data


def splitseq(row):
    seq = row.sequence
    start = row.start
    end = row.end
    intein = seq[start:end]
    rnr = seq[0:start] + seq[end:-1]
    return {'rnr': rnr, 'intein': intein}


def main():
    sequences, segments, accs = load_data()
    inteins = mark_inteins(sequences, segments, accs)
    records = []
    for ind, row in inteins[inteins.span > 0].iterrows():
        records.append( measure(row) )

    return pd.DataFrame(records)


def writeInteinRow(row, outfilename):
    with open(outfilename, "a") as fastafile:
        intein = splitseq(row)['intein']
        fastafile.write(f">{row.refseq_accno}|intein\n{intein}\n")


def writeInteins(outfilename):
    sequences, segments, accs = load_data()
    inteins = mark_inteins(sequences, segments, accs)
    records = []
    for ind, row in inteins[inteins.span > 0].iterrows():
        writeInteinRow(row, outfilename)


if __name__=='__main__':
    writeInteins("inteins.fasta")
