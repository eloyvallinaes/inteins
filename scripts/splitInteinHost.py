#!/bin/python

"""
Split sequences into inteins and host proteins and write to FASTA file.
"""

from src import parse
from pathlib import Path

fasta = Path("fasta")

parser = parse.InterproSegments(fasta / "IPR036844.fasta")
parser.to_fasta("host", fasta / "IPR036884_host.fasta")
parser.to_fasta("intein", fasta / "IPR036884_intein.fasta")
