"""
Parsers for different file formats used across this project.
"""

import re
import csv
import json
from pathlib import Path
from Bio import SeqIO
from Bio import SearchIO


class HMMSearchParser:
    KEYS = [
        'motifname',
        'motifacc',
        'refseq_accno',
        'start',
        'end',
        'evalue'
    ]

    def parse(self, filename):
        """
        Extract fields of interest and add a list of dictionaries to
        the records attribute.
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
        self.records = records

    def convert(self, infilename, outfilename):
        """
        Write parsed records HMMSearchParser.KEYS to CSV file.
        """
        self.parse(infilename)
        with open(outfilename, "w") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=HMMSearchParser.KEYS)
            writer.writeheader()
            for row in self.records:
                writer.writerow(row)


class CSV2Fasta:
    def load(self, infilename):
        self.records = self.parse(infilename)

    def parse(self, infilename):
        with open(infilename, "r") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            records = {row['refseq_accno']: row['sequence'] for row in reader}
        return records

    def write(self, outfilename):
        with open(outfilename, "w") as outfile:
            for acc, seq in self.records.items():
                outfile.write(f">{acc}\n{seq}\n")


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
                        "sequence": record.seq.replace("-", "").upper()
                })
        return records

    def write(self, outfilename):
        with open(outfilename, "w") as outfile:
            for record in self.records:
                outfile.write(
                        ">{}\n{}\n".format(*record.values)
                )


class Fasta2Dict:
    @staticmethod
    def parse(filename):
        records = {}
        with open(filename) as fastafile:
            for row in SeqIO.parse(fastafile, "fasta"):
                records[row.id] = row.seq.upper()
        return records


class InterproSegments:
    def parse(self, filename):
        records = {}
        with open(filename) as fastafile:
            for row in SeqIO.parse(fastafile, "fasta"):
                acc, limits, name = row.id.split("|")
                seq = str(row.seq.upper())
                inteins = InterproSegments.extractInteins(seq, limits)
                host = InterproSegments.extractHost(seq, inteins)
                records[acc] = {
                    'inteins': inteins,
                    'host': host,
                    }
        return records

    @staticmethod
    def parselimits(limits: str):
        pattern = re.compile(r"([0-9]+)\.\.\.([0-9]+)")
        for start, end in pattern.findall(limits):
            yield int(start), int(end)

    @staticmethod
    def extractInteins(sequence, limits: str):
        inteins = []
        for start, end in InterproSegments.parselimits(limits):
            inteins.append(
                {
                    'start': start,
                    'end': end,
                    'span': end-start,
                    'seq': sequence[start:end]
                }
            )
        return inteins

    @staticmethod
    def extractHost(sequence, inteins):
        indeces = [0]
        for i in inteins:
            indeces += [i['start']] + [i['end']]
        indeces += [-1]

        host = ''
        for e in range(0, len(indeces), 2):
            host += sequence[indeces[e]:indeces[e + 1]]

        return host


if __name__ == '__main__':
    segments = InterproSegments().parse("fasta/IPR036844.fasta")
