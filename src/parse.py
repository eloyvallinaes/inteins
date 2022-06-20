"""
Parsers for different file formats used across this project.
"""

import csv
from pathlib import Path
from Bio import SearchIO
from Bio import SeqIO



class HMMSearchParser:
    KEYS = [
        'motifname',
        'motifacc',
        'refseq_accno',
        'start',
        'end',
        'evalue'
    ]
    def __init__(self, infilename):
        self.records = self.parse(infilename)

    def parse(self, filename):
        """
        Extract fields of interest and return a list of dictionaries
        """
        parser = SearchIO.parse(filename, format="hmmsearch3-domtab")
        records = []
        for results in parser:
            queryacc = results.accession
            queryname = results.id
            for hit in results:
                hitname = hit.id
                for domain in hit:
                    values = [
                        queryname,
                        queryacc,
                        hitname,
                        domain.hit_start,
                        domain.hit_end,
                        domain.evalue
                    ]
                    records.append({
                        key: val
                        for key, val in zip(HMMSearchParser.KEYS, values)
                    })
        return records


    def write(self, outfilename):
        """
        Write fields of interest in HMMSearchParser.KEYS to CSV file.
        """
        with open(outfilename, "w") as outfile:
            writer = csv.DictWriter(outfile, fieldnames = HMMSearchParser.KEYS)
            writer.writeheader()
            for row in self.records:
                writer.writerow(row)


class CSV2Fasta:
    KEYS = [
        "refseq_accno",
        "sequence"
    ]
    def __init__(self, infilename):
        self.records = self.parse(infilename)

    def parse(self, infilename):
        records = []
        with open(infilename, "r") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                acc = row['refseq_accno']
                seq = row['sequence']
                records.append({
                    key: val
                    for key, val in row.items()
                    if key in CSV2Fasta.KEYS
                })
        return records

    def write(self, outfilename):
        with open(outfilename, "w") as outfile:
            for record in self.records:
                outfile.write(">{}\n{}\n".format(*record.values()))



class Aln2Fasta:
    KEYS = [
        "id",
        "sequence"
    ]
    def __init__(self, infilename):
        self.records = self.parse(infilename)

    def parse(self, infilename):
        records = []
        with open(infilename) as infile:
            for record in SeqIO.parse(infile, "fasta"):
                records.append({
                        "id": record.id,
                        "sequence": record.seq.replace("-","").upper()
                })
        return records

    def write(self, outfilename):
        with open(outfilename, "w") as outfile:
            for record in self.records:
                outfile.write(
                        ">{}\n{}\n".format(*record.values)
                )


if __name__ == '__main__':
    CSV2Fasta('csv/COG0209.csv').write("test.fasta")
