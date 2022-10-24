#!/bin/python

"""
Split sequences into inteins and host proteins and write to FASTA file.
"""

from src import parse
from pathlib import Path

fasta = Path("fasta/interpro")

parser = parse.Fasta2Dict(fasta / "IPR036844.fasta")
parser.write_segments(fasta / "inteins.fasta", "inteins")
parser.write_segments(fasta / "hosts.fasta", "hosts")
