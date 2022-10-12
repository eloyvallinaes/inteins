# /usr/bin/python

import csv
import json
import requests
from requests.adapters import HTTPAdapter
from pathlib import Path

url = "https://rest.uniprot.org/uniprotkb/search?query={}&fields=lineage"

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
FIELDNAMES.extend(['taxid', 'accession'])


def load(filename):
    with open(filename, "r") as jsonfile:
        j = json.load(jsonfile)
        for result in j['results']:
            yield result


def taxonomy(result):
    acc = result['to']['primaryAccession']
    taxid = result['to']['organism']['taxonId']
    s = requests.Session()
    a = requests.adapters.HTTPAdapter(max_retries=3)
    s.mount('http://', a)
    r = s.get(url.format(acc))
    if r.ok:
        data = r.json()
        lineages = data["results"][0]["lineages"]
        row = {
            pair["rank"]: pair["scientificName"]
            for pair in lineages if pair["rank"] in RANKS
        }
        row.update({
            "accession": acc, "taxid": taxid
        })

        return row


def main(jsonname, outname):
    with open(outname, "w") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=FIELDNAMES)
        writer.writeheader()
        for result in load(jsonname):
            writer.writerow(taxonomy(result))


if __name__ == '__main__':
    jsonPath = Path("json")
    try:
        known = set([
            line.split(",")[0]
            for line in open('ipr036844_taxonomy.csv', "r").readlines()
        ])
    except FileNotFoundError:
        known = {}

    rows = []
    with open("ipr036844_taxonomy.csv", "a") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
        if len(known) > 0:
            writer.writeheader()

        for result in load(jsonPath / "ipr036844_uniprot.json"):
            acc = result['to']['primaryAccession']
            if acc in known:
                continue

            known.add(acc)

            try:
                writer.writerow(taxonomy(result))

            except KeyError:
                print(f"{acc}: lineage not found")
