import os
import math
import argparse
from pathlib import Path

from measure import MeasureInterpro

fasta = Path("fasta")
interpro = Path("interpro")
csv = Path("csv")


def runInterpro(ipcode):
    fastafile = fasta / f"{ipcode}.fasta"
    m = MeasureInterpro(fastafile)
    m.write_hostdata(csv / f"{ipcode}_host_phys.csv")
    m.write_inteindata(csv / f"{ipcode}_intein_phys.csv")

    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
                "ipcodes",
                type=str,
                nargs='+',
                help="One or serveral Interpro ids"
    )

    args = parser.parse_args()
    for ipcode in args.ipcodes:
        runInterpro(ipcode)
