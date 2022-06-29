import os
import math
import argparse
from pathlib import Path

from fetch import COGFetch
from signatureScan import scan
from measure import Measure

cogs = Path("cogs")
csv = Path("csv")
fasta = Path("fasta")

def run(cogid):
    # fetch sequences from Conserved Orthologous Groups
    fastafile = fasta / f"{cogid}.fasta"
    cogfile = cogs / f"{cogid}.tsv"
    if os.path.exists(fastafile) and os.path.exists(cogfile):
        print("Re-starting from sequences already in disk")
    else:
        print(f"Fetching sequences from COG and RefSeq for {cogid}")
        fetcher = COGFetch(cogid, pagefrom=1, pageto=math.inf)
        fetcher.writeout()
        fetcher.writefasta()

    # scan with hmmsearch
    print("Scanning")
    scan(cogid)

    # measure properties
    print("Measuring")
    seqfile = cogs / f"{cogid}.tsv"
    nterms = csv / f"TIGR01445.1_{cogid}.csv"
    cterms = csv / f"TIGR01443.1_{cogid}.csv"
    m = Measure(seqfile, nterms, cterms)
    m.writeMeasurements(csv / f"{cogid}_phys.csv")

    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
                        "cogids",
                        type = str,
                        nargs='+',
                        help = "One or serveral Conserved Orthologous Group id(s)"
                    )

    args = parser.parse_args()
    for cogid in args.cogids:
        run(cogid)
