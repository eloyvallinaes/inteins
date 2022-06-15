"""
Read results in domain hits table (--domtblout)
"""

import csv
from Bio import SearchIO

KEYS = [
    'motifname',
    'motifacc',
    'refseq_accno',
    'start',
    'end',
    'evalue'
]

def parse(filename):
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

                records.append( {key: val for key, val in zip(KEYS, values)} )
    return records


def writeout(records, filename):
    with open(filename, "w") as outfile:
        writer = csv.DictWriter(outfile, fieldnames = KEYS)
        writer.writeheader()
        for row in records:
            writer.writerow(row)

def main(infilename, outfilename):
    records = parse(infilename)
    writeout(records, outfilename)


if __name__ == '__main__':
    main("csv/nterm_results.tab", "csv/nterm_segments.csv")
