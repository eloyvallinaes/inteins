"""
Scan sequences for prosite motifs
"""

import csv
import json
from Bio.ExPASy import ScanProsite

SIGNATURES = ["PS50817", "PS50819", "PS50818"]

class Hit():
    def __init__(
        self,
        sequence_ac,
        sequence_id,
        sequence_db,
        start,
        stop,
        signature_ac,
        signature_id,
        score,
        level
    ):
        self.acc = sequence_ac
        self.signature = signature_ac
        self.signame = signature_id
        self.start = start
        self.stop = stop
        self.level = level


def scan(tsvfile):
    with open(tsvfile, "r") as tsvfile:
        out = []
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            id = row["refseq_accno"]
            sequence = row['sequence']
            seq = f">{id}:{sequence}:"
            xmlhandle = ScanProsite.scan(
                seq = seq,
                sig = " and ".join(SIGNATURES)
            )
            out.append(ScanProsite.read(xmlhandle))
    return out


def load(filename):
    with open(filename, "r") as jsonfile:
        results =  json.load(jsonfile)

    return [Hit(**result) for result in results]

def edges(acc, hits):
    records = [hit for hit in hits if hit.acc == acc]
    starts = []
    stops = []
    for sig in SIGNATURES:
        starts.extend([hit.start for hit in records if hit.signature == sig])
        stops.extend([hit.stop for hit in records if hit.signature == sig])

    return min(starts), max(stops), max(stops) - min(starts)


if __name__ == "__main__":
    results = scan("csv/rnr_refseq_proteins.tsv")
