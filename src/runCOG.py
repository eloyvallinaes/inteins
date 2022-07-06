import os
import math
import argparse
from pathlib import Path

from fetch import COGFetch
from signatureScan import scan
from parse import HMMSearchParser, CSV2Fasta
from measure import Measure

fasta = Path("fasta")
cogs = Path("cogs")
csv = Path("csv")
scans = Path("scans")
terms = Path("terms")


def run(cogid, domE):
    cogfile = cogs / f"{cogid}.tsv"
    # fetch sequences from Conserved Orthologous Groups
    if os.path.exists(cogfile):
        print("Re-starting from sequences already in disk")
    else:
        print(f"Fetching sequences from COG and RefSeq for {cogid}")
        fetcher = COGFetch(cogid, pagefrom=1, pageto=math.inf)
        fetcher.writeout()

    # parse COG file to fasta
    fastafile = fasta / f"{cogid}.fasta"
    parser = CSV2Fasta()
    parser.load(cogfile)
    parser.write(fastafile)

    # scan with hmmsearch
    print("Scanning")
    scan(cogid, domE)

    # parse hmmer results
    parser = HMMSearchParser()
    parser.convert(
        scans / f"TIGR01445.1_{cogid}.tab",
        terms / f"TIGR01445.1_{cogid}.csv"
    )
    parser.convert(
        scans / f"TIGR01443.1_{cogid}.tab",
        terms / f"TIGR01443.1_{cogid}.csv"
    )

    # measure properties
    print("Measuring")
    seqfile = cogs / f"{cogid}.tsv"
    nterms = terms / f"TIGR01445.1_{cogid}.csv"
    cterms = terms / f"TIGR01443.1_{cogid}.csv"
    m = Measure(seqfile, nterms, cterms)
    m.write_hostdata(csv / f"{cogid}_host_phys.csv")
    m.write_inteindata(csv / f"{cogid}_intein_phys.csv")

    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
                "cogids",
                type=str,
                nargs='+',
                help="One or serveral Conserved Orthologous Group id(s)"
    )
    parser.add_argument(
                "-domE",
                type=float,
                help="Conditional evalue for detection of intein signatures",
                default=0.01
    )

    args = parser.parse_args()
    for cogid in args.cogids:
        run(cogid, args.domE)
