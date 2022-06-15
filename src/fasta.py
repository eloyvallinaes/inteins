"""
Dealigned MSA in FASTA format and obtaine ungapped FASTA with whole sequence
in one line.
"""

from Bio import SeqIO
import csv

def csv2fasta(inname, outname):
    with open(outname, "w") as outfile:
        with open(inname, "r") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                    acc = row['refseq_accno']
                    seq = row['sequence']
                    outfile.write(f">{acc}\n{seq}\n")

def aln2plain(inname, outname):
    with open(outname, "w") as outfile:
        with open(inname) as infile:
            for record in SeqIO.parse(infile, "fasta"):
                sequence = record.seq.replace("-","")
                outfile.write(f">{record.id}\n{sequence}\n")


if __name__ == '__main__':
    csv2fasta("rnr_refseq_proteins.tsv", "rnr_refseq.fasta")
