# /usr/bin/python

import csv
import requests
from pathlib import Path
from requests.adapters import HTTPAdapter
from time import sleep

from src import parse

url = "https://www.ebi.ac.uk/proteins/api/taxonomy/lineage/{}"

RANKS = [
    'superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species'
]

FIELDNAMES = RANKS.copy()
FIELDNAMES.append('taxid')


def taxonomy(taxid):
    s = requests.Session()
    a = requests.adapters.HTTPAdapter(max_retries=3)
    s.mount('http://', a)
    r = s.get(url.format(taxid))
    if r.ok:
        data = r.json()
        row = {
            pair["rank"]: pair["scientificName"]
            for pair in data['taxonomies'] if pair["rank"] in RANKS
        }

        return row


if __name__ == '__main__':
    fasta = Path("fasta/interpro")
    try:
        known = set([
            line.split(",")[0]
            for line in open('taxonomy.csv', "r").readlines()
        ])
    except FileNotFoundError:
        known = set()

    rows = []
    i = 0
    with open("taxonomy.csv", "a") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
        if len(known) == 0:
            writer.writeheader()

        parser = parse.Fasta2Dict(fasta / "IPR036844.fasta")
        taxids = set(parser.taxids.values())
        for taxid in taxids:
            if taxid in known:
                continue

            try:
                row = taxonomy(taxid)
                row.update({"taxid": taxid})
                writer.writerow(row)

            except KeyError:
                print(f"{acc}: lineage not found")

            known.add(taxid)
            sleep(0.5)
